#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 19:17:27 2020

@author: aero_user
"""

import pandas as pd
import math
from typing import List, Dict, Tuple


def cluster_junctions(junc_counts, ref_length, threshold, width=6) -> Dict[int, int]:

    clusters = {}

    kmer_count = {}
    for pos in junc_counts['pos']:
        kmer_count[pos] = junc_counts[(junc_counts['pos'] >= pos) &
                                      (junc_counts['pos'] < pos + width)]['count'].sum()

    while True:

        if len(kmer_count) == 0:
            break

        max_count = max(kmer_count.values())

        if max_count < threshold:
            break

        sorted_keys = sorted(kmer_count.keys())

        for pos in sorted_keys:
            if kmer_count[pos] >= threshold:
                curr_pos = pos
                curr_count = kmer_count[curr_pos]
                break
        # max_pos = max(kmer_count, key=kmer_count.get)

        # clusters[max_pos] = max_count

        assert(curr_pos not in clusters.keys())
        clusters[curr_pos] = curr_count

        # junc_counts.drop(junc_counts[(junc_counts['pos'] >= max_pos) &
        #                   (junc_counts['pos'] < max_pos + width)].index, inplace=True)
        junc_counts.drop(junc_counts[(junc_counts['pos'] >= curr_pos) &
                                     (junc_counts['pos'] < curr_pos + width)].index, inplace=True)
        for idx in range(sorted_keys.index(curr_pos) + 1):
            del kmer_count[sorted_keys[idx]]
        for pos in range(curr_pos + 1, curr_pos + width):
            if pos in kmer_count.keys():
                del kmer_count[pos]

        # for pos in range(max_pos-width, max_pos+width):
        #     kmer_count[pos] = junc_counts[(junc_counts['pos'] >= pos) &
        #                                 (junc_counts['pos'] < pos + width)]['count'].sum()
        for pos in range(curr_pos-width, curr_pos+width):
            if pos in kmer_count.keys():
                kmer_count[pos] = junc_counts[(junc_counts['pos'] >= pos) &
                                              (junc_counts['pos'] < pos + width)]['count'].sum()

    # print(min(clusters.keys()))
    return clusters


def readInputPhasing(inputPhasing) -> List[List[Tuple[int, int, str]]]:

    phasing = []
    phase = []
    with open(inputPhasing, 'r') as input:
        for line in input:
            if line.startswith('>'):
                if len(phase) > 1:
                    phasing.append(phase)
                    phase = []
            else:
                data = line[1:-1].split(',')
                phase.append(
                    (int(data[0]), int(data[1]), data[2].split('\'')[1]))
        phasing.append(phase)

    return phasing


def readInputSJreadsPhasing(inputCSV, inputPhasing) -> Tuple[pd.DataFrame, List[List[Tuple[int, int, str]]]]:
    df_sj_reads = pd.read_csv(inputCSV, sep='\t')

    phasing = []
    phase = []
    with open(inputPhasing, 'r') as input:
        for line in input:
            if line.startswith('>'):
                if len(phase) > 1:
                    phasing.append(phase)
                    phase = []
            else:
                data = line[1:-1].split(',')
                phase.append(
                    (int(data[0]), int(data[1]), data[2].split('\'')[1]))
        phasing.append(phase)

    return df_sj_reads, phasing


def generateGraphInput(SJcounter, bam, df_sj_reads, phasing, mutList, df_depth, ref, width, sj_threshold):

    graph_edges = {}
    #breakpoints = set()

    start_junction_clusters = cluster_junctions(df_sj_reads.groupby(['5prime']).sum()['count'].reset_index()
                                                .set_axis(['pos', 'count'], axis=1, inplace=False),
                                                ref.lengths[0], sj_threshold, width=width)

    end_junction_clusters = cluster_junctions(df_sj_reads.groupby(['3prime']).sum()['count'].reset_index()
                                              .set_axis(['pos', 'count'], axis=1, inplace=False),
                                              ref.lengths[0], sj_threshold, width=width)

    # splice edges
    for _, row in df_sj_reads.iterrows():
        start_dist = [abs(x + float(width)/2 - row['5prime'])
                      for x in list(start_junction_clusters.keys())]
        min_start = min(start_dist)

        end_dist = [abs(x + float(width)/2 - row['3prime'])
                    for x in list(end_junction_clusters.keys())]
        min_end = min(end_dist)

        if min_start <= float(width)/2 and min_end <= float(width)/2:
            start_pos = list(start_junction_clusters.keys())[
                start_dist.index(min_start)]
            end_pos = list(end_junction_clusters.keys())[
                end_dist.index(min_end)]

            # breakpoints.add(start_pos)
            # breakpoints.add(end_pos)

            if (start_pos, end_pos) not in graph_edges.keys():
                graph_edges[(start_pos, end_pos)] = {'splice': row['count']}
            else:
                graph_edges[(start_pos, end_pos)]['splice'] += row['count']

    # prune the splice edges and build the breakpoint locations
    breakpoints = set()
    # pruned_edges = [(left_pos, right_pos) for (left_pos, right_pos), count in graph_edges.items()
    #                 if count < sj_threshold]
    pruned_edges = []
    for (left_pos, right_pos) in graph_edges.keys():
        if graph_edges[(left_pos, right_pos)]['splice'] < sj_threshold:
            pruned_edges.append((left_pos, right_pos))
        else:
            breakpoints = breakpoints.union({left_pos, right_pos})
    for entry in pruned_edges:
        del graph_edges[entry]

    # add mutation breakpoints
    breakpoints = breakpoints.union(set(mutList))

    # mutation and adjacency edges
    for pos in breakpoints:
        if pos not in mutList:
            # adjacency edge
            if str(pos) not in list(df_depth['pos']):
                graph_edges[(pos, pos)] = {}
            else:
                graph_edges[(pos, pos)] = {ref.fetch(SJcounter.contig, pos-1, pos):
                                           int(df_depth[df_depth['pos'] == str(pos)]['count'].values[0])}
        else:
            # mutation edges
            graph_edges[(pos, pos)] = SJcounter.get_allele_counts(bam, pos)

    # build clustered phasing reads
    clustered_phasing = []
    for phase in phasing:
        newphase = []
        for left_pos, right_pos, edgeType in phase:
            if edgeType == 'splice':
                left_dist = [abs(x - left_pos)
                             for x in list(start_junction_clusters.keys())]
                min_left_dist = min(left_dist)

                right_dist = [abs(x - right_pos)
                              for x in list(end_junction_clusters.keys())]
                min_right_dist = min(right_dist)

                if min_left_dist > width or min_right_dist > width:
                    break
                else:
                    newphase.append((list(start_junction_clusters.keys())[left_dist.index(min_left_dist)],
                                     list(end_junction_clusters.keys())[
                        right_dist.index(min_right_dist)],
                        edgeType))
            else:
                newphase.append((left_pos, right_pos, edgeType))
        if len(newphase) == len(phase):
            clustered_phasing.append(newphase)

    return breakpoints, graph_edges, clustered_phasing


def readGraphInput(breakpointFile, graphEdgeFile, phasingFile):

    breakpoints = set()
    if breakpointFile:
        with open(breakpointFile, 'r') as input:
            for line in input:
                data = int(line.rstrip('\n'))
                breakpoints.add(data)

    graph_edges = {}
    if graphEdgeFile:
        with open(graphEdgeFile, 'r') as input:
            for line in input:
                data = line.rstrip('\n').split('\t')
                left_pos = int(data[0])
                right_pos = int(data[1])
                edgetype = data[2]
                count = float(data[3])
                if (left_pos, right_pos) not in graph_edges.keys():
                    graph_edges[(left_pos, right_pos)] = {edgetype: count}
                else:
                    graph_edges[(left_pos, right_pos)][edgetype] = count

    clustered_phasing = []
    if phasingFile:
        phase = []
        with open(phasingFile, 'r') as input:
            for line in input:
                if line.startswith('>'):
                    if len(phase) > 0:
                        clustered_phasing.append(phase)
                    phase = []
                else:
                    data = line.rstrip('\n').split('\t')
                    left_pos = int(data[0])
                    right_pos = int(data[1])
                    edgetype = data[2]
                    phase.append((left_pos, right_pos, edgetype))
        if len(phase) > 0:
            clustered_phasing.append(phase)

    return breakpoints, graph_edges, clustered_phasing


def getGraphSkeleton(df_sj_reads, ref, width, sj_threshold):

    splice_edges = {}
    adjacency_edges = {}
    #breakpoints = set()

    start_junction_clusters = cluster_junctions(df_sj_reads.groupby(['5prime']).sum()['count'].reset_index()
                                                .set_axis(['pos', 'count'], axis=1, inplace=False),
                                                ref.lengths[0], sj_threshold, width=width)

    end_junction_clusters = cluster_junctions(df_sj_reads.groupby(['3prime']).sum()['count'].reset_index()
                                              .set_axis(['pos', 'count'], axis=1, inplace=False),
                                              ref.lengths[0], sj_threshold, width=width)

    # splice edges
    for _, row in df_sj_reads.iterrows():
        start_dist = [abs(x + float(width)/2 - row['5prime'])
                      for x in list(start_junction_clusters.keys())]
        min_start = min(start_dist)

        end_dist = [abs(x + float(width)/2 - row['3prime'])
                    for x in list(end_junction_clusters.keys())]
        min_end = min(end_dist)

        if min_start <= float(width)/2 and min_end <= float(width)/2:
            start_pos = list(start_junction_clusters.keys())[
                start_dist.index(min_start)]
            end_pos = list(end_junction_clusters.keys())[
                end_dist.index(min_end)]

            # breakpoints.add(start_pos)
            # breakpoints.add(end_pos)

            if (start_pos, end_pos) not in splice_edges.keys():
                splice_edges[(start_pos, end_pos)] = {'splice': row['count']}
            else:
                splice_edges[(start_pos, end_pos)]['splice'] += row['count']

    # prune the splice edges and build the breakpoint locations
    breakpoints = set()
    # pruned_edges = [(left_pos, right_pos) for (left_pos, right_pos), count in graph_edges.items()
    #                 if count < sj_threshold]
    pruned_edges = []
    for (left_pos, right_pos) in splice_edges.keys():
        if splice_edges[(left_pos, right_pos)]['splice'] < sj_threshold:
            pruned_edges.append((left_pos, right_pos))
        else:
            breakpoints = breakpoints.union({left_pos, right_pos})
    for entry in pruned_edges:
        del splice_edges[entry]

    # mutation and adjacency edges
    for pos in breakpoints:
        adjacency_edges[(pos, pos)] = {ref.fetch(
            ref.references[0], pos-1, pos): 0}

    return breakpoints, splice_edges, adjacency_edges


def readGraphSkeleton(breakpointFile, graphEdgeFile):

    breakpoints = set()
    if breakpointFile:
        with open(breakpointFile, 'r') as input:
            for line in input:
                data = int(line.rstrip('\n'))
                breakpoints.add(data)

    splice_edges = {}
    adjacency_edges = {}
    if graphEdgeFile:
        with open(graphEdgeFile, 'r') as input:
            for line in input:
                data = line.rstrip('\n').split('\t')
                left_pos = int(data[0])
                right_pos = int(data[1])
                edgetype = data[2]
                count = float(data[3])
                if edgetype == 'splice':
                    if (left_pos, right_pos) not in splice_edges.keys():
                        splice_edges[(left_pos, right_pos)] = {edgetype: count}
                    else:
                        splice_edges[(left_pos, right_pos)][edgetype] = count
                else:
                    if (left_pos, right_pos) not in adjacency_edges.keys():
                        adjacency_edges[(left_pos, right_pos)] = {
                            edgetype: count}
                    else:
                        adjacency_edges[(left_pos, right_pos)
                                        ][edgetype] = count

    return breakpoints, splice_edges, adjacency_edges


def writeBreakpoints(fname, sorted_breakpoints):
    with open(fname, 'w') as output:
        for entry in sorted_breakpoints:
            output.write(f"{entry}\n")


def writePhasing(fname, phasing_reads):
    with open(fname, 'w') as output:
        for idx, (splus, sminus) in enumerate(phasing_reads.keys()):
            output.write(f">phasing{idx}\t{phasing_reads[splus, sminus]}\n")
            for left, right in splus:
                output.write(f"+\t{left}\t{right}\n")
            for left, right in sminus:
                output.write(f"-\t{left}\t{right}\n")
            #output.write(f">phasing {idx} -- {phasing_reads[(splus, sminus)]}\n>>{splus}\n>>{sminus}\n")


def readPhasing(fname):
    phasing = {}
    print(f"reading {fname}")
    idx = 0
    with open(fname, 'r') as input:
        flag = False
        for line in input:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if not flag:
                    splus = set()
                    sminus = set()
                    count = int(line.split('\t')[-1])
                    flag = True
                else:
                    if tuple([tuple(splus), tuple(sminus)]) in phasing.keys():
                        print(
                            f"phasing repeat -- {idx} -- {list(phasing.keys()).index(tuple([tuple(splus),tuple(sminus)]))}")
                        print(f"{splus}")
                        print(f"{sminus}")
                    phasing[tuple(
                        [tuple(sorted(splus)), tuple(sorted(sminus))])] = count
                    splus = set()
                    sminus = set()
                    count = int(line.split('\t')[-1])
                idx += 1
            else:
                if line.startswith('+'):
                    data = line.split('\t')
                    splus.add((int(data[1]), int(data[2])))
                elif line.startswith('-'):
                    data = line.split('\t')
                    sminus.add((int(data[1]), int(data[2])))
        phasing[tuple([tuple(sorted(splus)), tuple(sorted(sminus))])] = count

    return phasing


def writeEdges(fname, edgeList):
    with open(fname, 'w') as output:
        for edges in edgeList:
            for (left_pos, right_pos), edge_dict in edges.items():
                for edgeType, count in edge_dict.items():
                    output.write(
                        f"{left_pos}\t{right_pos}\t{edgeType}\t{count}\n")


def isoverlap(edge1, edge2):
    if edge1.left.right <= edge2.right.left and edge1.right.left >= edge2.left.right:
        return True
    else:
        return False
