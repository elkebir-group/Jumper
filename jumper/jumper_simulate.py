#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 2020

@author: Palash Sashittal
"""

import pysam
import pandas as pd
import numpy as np
import sys
import os
import argparse
from collections import Counter
import math
import random

src = os.path.join(os.path.dirname(os.path.abspath(__file__)))
if not os.path.isdir(src):
    raise ValueError("jumper directory cannot be found i.e. {}, is anything been moved?".format(src))
sys.path.append(src)

from jumper_main import segmentGraph
from segment_graph_aux import *


def getPathProb(path):
    prob = 1
    for edge in path:
        prob *= edge.weight
    return prob


def getVectorBernoulli(size, p=0.5):
    result = []
    for i in range(size):
        result.append(np.random.binomial(1, p))
    return result


def chooseFromProb(probs, cum=False):
    if not cum:
        assert(sum(probs) == 1)
        cumprobs = [sum(probs[:i+1])/sum(probs) for i in range(len(probs))]
    else:
        assert(probs[-1] == 1)
        cumprobs = probs

    rand_num = np.random.random()
    for idx in range(len(cumprobs)):
        if rand_num <= cumprobs[idx]:
            return idx
    raise Exception(
        f"rand_num ({rand_num}) more than cumprobs[-1] ({cumprobs[-1]})")


def getTerminalNodes(fragment, path_nodes):
    assert(len(fragment) == 2)

    # print(f"node = 0 left pos is {path_nodes[0].left}")
    left_node = None
    right_node = None
    shift = 0
    prev = 1
    for node in path_nodes:
        curr = node.left
        shift += curr - prev
        # print(f"{node.left} -- {node.right} -- {fragment[0] + shift} -- {shift}")
        if (node.left <= fragment[0] + shift and
                node.right > fragment[0] + shift):
            # print("here")
            left_node = node
            break
        prev = node.right

    shift = 0
    prev = 1
    for node in path_nodes:
        curr = node.left
        shift += curr - prev
        # print(f"{node.left} -- {node.right} -- {fragment[1] + shift} -- {shift}")
        if (node.left < fragment[1] + shift and
                node.right >= fragment[1] + shift):
            # print("here")
            right_node = node
            break
        prev = node.right

    # print(f"{left_node.name} \t {right_node.name}")

    # print(f"path nodes are --")
    # for node in path_nodes:
    #     print(f"{node.left} \t {node.right}")

    if not left_node:
        raise Exception(
            f"left node not found for fragment ({fragment[0]}, {fragment[1]})")

    if not right_node:
        raise Exception(
            f"right node not found for fragment ({fragment[0]}, {fragment[1]})")

    return left_node, right_node


def getSplusSminus(fragment, curr_path_nodes, curr_path, spliceEdges):

    splus = set()
    sminus = set()

    left_node, right_node = getTerminalNodes(fragment, curr_path_nodes)

    left_index = curr_path_nodes.index(left_node)
    right_index = curr_path_nodes.index(right_node)

    for edge in curr_path[left_index:right_index]:
        if edge.type == 'splice':
            # splus.add((edge.left.right, edge.right.left))
            splus.add(edge)

    for node in curr_path_nodes[left_index:right_index+1]:
        for edge in spliceEdges:
            if edge.left.right < node.right and edge.right.left > node.left:
                # sminus.add((edge.left.right, edge.right.left))
                sminus.add(edge)

    return splus, sminus


class simulator():

    def __init__(self, ref, breakpoints, edges, sense=True, verbose=True, canonicalprior=None, canonicalFlag=False,
                 nreads=None, npaths=None, readlength=None, insertlength=None):
        self.ref = ref
        self.verbose = verbose
        self.graph = segmentGraph(breakpoints=breakpoints, edges=edges, verbose=self.verbose,
                                  ref=self.ref)
        self.graph.initCanonicalTranscripts()
        self.sense = sense
        if nreads:
            self.nreads = nreads
        else:
            self.nreads = 1000
        if npaths:
            self.npaths = npaths
        else:
            self.npaths = 1
        if readlength:
            self.readlength = readlength
        else:
            self.readlength = 100
        if insertlength:
            self.insertlength = insertlength
        else:
            self.insertlength = 250
        if canonicalprior:
            self.canonicalprior = canonicalprior
        else:
            self.canonicalprior = 10

        self.canonicalFlag = canonicalFlag

        self.initProbabilityEdgeWeights()
        self.paths, self.abundances = self.drawRandomPaths()
        self.path_lengths = None
        self.path_probs = [getPathProb(path) for path in self.paths]
        self.total_log_prob_path = 0
        for idx, prob in enumerate(self.path_probs):
            self.total_log_prob_path += self.abundances[idx] * math.log(prob)
        self.phasing_reads = self.generatePhasingReads()

    def initProbabilityEdgeWeights(self):

        for node in self.graph.nodes:
            if self.sense:

                alpha = []
                for edge in self.graph.outArcs(node):
                    if edge not in self.graph.nonTranscriptSpliceEdges:
                        alpha.append(self.canonicalprior)
                    else:
                        alpha.append(1)

                prob = np.random.dirichlet(alpha)

                idx = 0
                for edge in self.graph.outArcs(node):
                    edge.weight = prob[idx]
                    idx += 1

            else:

                alpha = []
                for edge in self.graph.inArcs(node):
                    if edge not in self.graph.nonTranscriptSpliceEdges:
                        alpha.append(self.canonicalprior)
                    else:
                        alpha.append(1)

                prob = np.random.dirichlet(alpha)

                idx = 0
                for edge in self.graph.inArcs(node):
                    edge.weight = prob[idx]
                    idx += 1

    def drawRandomPaths(self):

        drawn_paths = Counter()
        for i in range(self.npaths):

            if self.sense:
                start_node = self.graph.nodes[0]
                end_node = self.graph.nodes[-1]

                curr_node = start_node
                path = []
                while curr_node != end_node:
                    rand_num = np.random.random()
                    cumsum = 0
                    for edge in self.graph.outArcs(curr_node):
                        cumsum += edge.weight
                        if rand_num <= cumsum:
                            break
                    path.append(edge)
                    curr_node = edge.right

            else:
                start_node = self.graph.nodes[-1]
                end_node = self.graph.nodes[0]

                curr_node = start_node
                path = []
                while curr_node != end_node:
                    rand_num = np.random.random()
                    cumsum = 0
                    for edge in self.graph.inArcs(curr_node):
                        cumsum += edge.weight
                        if rand_num <= cumsum:
                            path.append(edge)
                            break
                    curr_node = edge.left

                path.reverse()

            if not self.canonicalFlag:
                drawn_paths[tuple(path)] += 1
            else:
                if path in self.graph.transcripts:
                    drawn_paths[tuple(path)] += 1

        # print(f"{drawn_paths.values()}")

        # print('_'*50)
        # for path, count in drawn_paths.items():
        #     print(f"{self.graph.getPathIndex(path)}\t{count}")

        # print('_'*50)
        return list(drawn_paths.keys()), list(drawn_paths.values())

    def generatePhasingReads(self):

        spliceEdges = self.graph.transcriptSpliceEdges + \
            self.graph.nonTranscriptSpliceEdges
        # splice_lengths = [edge.right.legt - edge.left.right for edge in spliceEdges]

        numPaths = len(self.paths)

        self.path_lengths = [29903] * numPaths
        for idx, path in enumerate(self.paths):
            for edge in path:
                if edge.type == 'splice':
                    self.path_lengths[idx] -= (edge.right.left -
                                               edge.left.right)

        probs = [self.abundances[idx] * self.path_lengths[idx] /
                 self.npaths for idx in range(numPaths)]
        cumprobs = [sum(probs[:i+1])/sum(probs) for i in range(numPaths)]

        self.total_log_prob_reads = 0
        phasing_reads = Counter()
        for idx in range(self.nreads):

            pathIndex = chooseFromProb(cumprobs, cum=True)
            curr_path = self.paths[pathIndex]
            self.total_log_prob_reads += math.log(
                probs[pathIndex] / sum(probs) / self.path_lengths[pathIndex])

            curr_length = self.path_lengths[pathIndex]

            # read1_start = 0
            # read2_start = 0
            # while abs(read2_start - read1_start) <= self.readlength + self.insertlength:
            #     read1_start = np.random.random() * curr_length
            #     read2_start = np.random.random() * curr_length

            # read1_start = np.random.random() * curr_length
            # read2_start = np.random.random() * curr_length

            read1_start = np.random.random() * (curr_length - 1 - 2 *
                                                self.readlength - self.insertlength) + 1
            read2_start = read1_start + self.readlength + self.insertlength

            # print(f"read start locations are {read1_start}, {read2_start}")

            curr_path_nodes = [curr_path[0].left] + \
                [edge.right for edge in curr_path]

            # read 1 analysis
            read1_splus, read1_sminus = getSplusSminus((read1_start, min(read1_start + self.readlength, curr_length)),
                                                       curr_path_nodes, curr_path, spliceEdges)

            read2_splus, read2_sminus = getSplusSminus((read2_start, min(read2_start + self.readlength, curr_length)),
                                                       curr_path_nodes, curr_path, spliceEdges)

            splus = read1_splus | read2_splus
            sminus = read1_sminus | read2_sminus

            for plus_edge in splus:
                removal = []
                for minus_edge in sminus:
                    if isoverlap(plus_edge, minus_edge):
                        removal.append(minus_edge)
                for edge in removal:
                    sminus.remove(edge)

            splus_tuples = [(edge.left.right, edge.right.left)
                            for edge in splus]
            sminus_tuples = [(edge.left.right, edge.right.left)
                             for edge in sminus]

            phasing_reads[tuple(
                [tuple(sorted(splus_tuples)), tuple(sorted(sminus_tuples))])] += 1

        return dict(phasing_reads)

    def writeReadCounts(self, fname):

        pathTuples = [(self.graph.getPathIndex(path),
                       self.abundances[idx] * self.path_lengths[idx] / self.npaths)
                      for idx, path in enumerate(self.paths)]
        # pathTuples = [(self.graph.getPathIndex(path),
        #               self.abundances[idx]) for idx, path in enumerate(self.paths)]
        pathTuples.sort(key=lambda x: x[0])

        # minAbundance = min(self.abundances)
        totalEffectiveLength = sum([x[1] for x in pathTuples])
        with open(fname, 'w') as output:
            for _, val in pathTuples:
                output.write(
                    f"{math.floor(self.nreads * val / totalEffectiveLength)}\n")
                # output.write(f"{val}\n")


def main(args):

    ref = pysam.FastaFile(args.fasta)

    breakpoints, splice_edges, adjacency_edges = readGraphSkeleton(
        args.inputBreakpoints, args.inputEdges)

    instance = simulator(ref=ref, breakpoints=breakpoints, edges={**splice_edges, **adjacency_edges},
                         sense=args.sense, verbose=args.verbose, canonicalprior=args.canonicalprior,
                         npaths=args.npaths, nreads=args.nreads, readlength=args.readlength,
                         insertlength=args.insertlength, canonicalFlag=args.canonical)

    if args.outputCanonical:
        instance.graph.writeTranscripts(args.outputCanonical)

    if args.outputEdges:
        writeEdges(args.outputEdges, [instance.graph.getEdges()])

    if args.outputDOT:
        instance.graph.writeDOT(args.outputDOT)

    if args.verbose:
        print(f"num of ST paths is {instance.graph.getNumSTpaths()}")
        print(
            f"abudances are -- {[x / instance.npaths for x in instance.abundances]}")
        print(
            f"total log probability paths is -- {instance.total_log_prob_path}")
        print(
            f"total log probability reads is -- {instance.total_log_prob_reads}")
        print(f"sense is {args.sense}")

    if args.outputPaths:
        instance.graph.writePaths(
            instance.paths, args.outputPaths, instance.abundances)

    if args.outputFasta:
        # instance.graph.writeFasta(instance.paths, args.outputFasta, instance.abundances)
        instance.graph.writeFasta(instance.paths, args.outputFasta)

    if args.outputReadCounts:
        instance.writeReadCounts(args.outputReadCounts)

    if args.outputPhasing:
        writePhasing(args.outputPhasing, instance.phasing_reads)


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1', '+', 'pos'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0', '-', 'neg'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", type=str,
                        help="fasta file", required=True)
    parser.add_argument("--inputBreakpoints", type=str,
                        help="input file contiaining breakpoints")
    parser.add_argument("--inputEdges", type=str,
                        help="input file contiaining graph edges")
    parser.add_argument("--outputEdges", type=str,
                        help="output file containing graph edges")
    parser.add_argument("--outputDOT", type=str,
                        help="output DOT file for splice graph")
    parser.add_argument("--outputCanonical", type=str,
                        help="output file for transcripts")
    parser.add_argument("--noverbose", dest="verbose",
                        help="do not output statements from internal solvers [default is true]", action='store_false')
    parser.add_argument("--npaths", type=int,
                        help='number of paths to draw [1]', default=1)
    parser.add_argument("--outputPaths", type=str,
                        help="output file randomly chosen paths")
    parser.add_argument("--sense", default=True, const=True, type=str2bool,
                        nargs='?', help="set pos strand flag for transcript generation")
    parser.add_argument(
        "--nreads", type=int, help="number of reads to be generated [100]", default=100)
    parser.add_argument("--outputPhasing", type=str,
                        help="output file for phasing reads")
    parser.add_argument("--readlength", type=int,
                        help="length of the reads (one end)", default=100)
    parser.add_argument("--insertlength", type=int,
                        help="length of the gap between the paired ends", default=50)
    parser.add_argument("--canonicalprior", type=int,
                        help="alpha value canonical discontinuous edge", default=10)
    parser.add_argument(
        "--seed", type=int, help="seed for the random number generators", default=0)
    parser.add_argument("--outputFasta", type=str,
                        help="output fasta file for drawn paths")
    parser.add_argument("--outputReadCounts", type=str,
                        help="output file with number of reads drawn per transcript")
    parser.add_argument("--canonical", type=str2bool,
                        help="flag to output only canonical transcripts [False]", default=False)

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    random.seed(args.seed)
    np.random.seed(args.seed)
    return args


def main_cli():
    arguments = get_options()
    main(arguments)


if __name__ == "__main__":
    arguments = get_options()
    main(arguments)
