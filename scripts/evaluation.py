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
import argparse
from typing import List, Dict, Tuple, Optional
from collections import Counter
import math
import itertools

from jumper.segment_graph_aux import *


class solution():
    
    def __init__(self, graph, fname, matching_fname = None):
        
        self.graph = graph
        self.paths = []
        self.path_abundances = []
        self.path_indices = []
        flag = False
        with open(fname, 'r') as inp:
            for line in inp:
                if line.startswith('>transcript'):
                    if flag:
                        assert(self.graph.getPathIndex(path) == pathIndex)
                        self.paths.append(tuple(path))
                        self.path_abundances.append(pathAbundance)
                        self.path_indices.append(pathIndex)
                        
                        path = []
                    pathIndex = int(line.rstrip('\n').split('\t')[0].split('_')[1])
                    pathAbundance = float(line.rstrip('\n').split(':')[1])
                    # print(f"transcript {pathIndex} is drawn")
                    self.paths.append(tuple(self.graph.transcripts[pathIndex]))
                    self.path_abundances.append(pathAbundance)
                    self.path_indices.append(pathIndex)
                    flag = False
                elif line.startswith('>path') or line.startswith('>Path'):
                    if flag:
                        assert(self.graph.getPathIndex(path) == pathIndex)
                        self.paths.append(tuple(path))
                        self.path_abundances.append(pathAbundance)
                        self.path_indices.append(pathIndex)
                    else:
                        flag = True
                    path = []
                    pathIndex = int(line.rstrip('\n').split(':')[0].split('_')[1])
                    pathAbundance = float(line.rstrip('\n').split(':')[1])
                elif flag:
                    data = line.rstrip('\n').split('\t')

                    left_node_name = data[0]
                    right_node_name = data[1]
                    edge_type = data[2]
                    
                    curr_edge = [ edge for edge in self.graph.edges if (edge.left.name == left_node_name and
                                                                        edge.right.name == right_node_name and
                                                                        edge.type == edge_type)]
                    assert(len(curr_edge) == 1)
                    path.append(curr_edge[0])

            if flag:
                assert(self.graph.getPathIndex(path) == pathIndex)
                self.paths.append(tuple(path))
                self.path_abundances.append(pathAbundance)
                self.path_indices.append(pathIndex)
        self.npaths = len(self.paths)

        self.pathCounter = Counter(dict(zip(self.path_indices, self.path_abundances)))

        if matching_fname:
            self.matching = {}
            with open(matching_fname, 'r') as inp:
                for line in inp:
                    if line.startswith('>'):
                        pathIndex = int(line.rstrip('\n').split('_')[1])
                        self.matching[pathIndex] = {}
                    else:
                        data = line.rstrip('\n').split('\t')
                        phaseID = int(data[0])
                        count = int(data[1])
                        self.matching[pathIndex][phaseID] = count


class alternate_solution():
    
    def __init__(self, fname, gtf = False):
       
        if gtf:
            self.read_gtf_file(fname)
        else:
            #self.read_gtf_file(fname)
            self.read_jumper_file(fname)
        

    def read_jumper_file(self, fname):

        self.path_jump_dict = {}
        self.path_abundances_dict = {}
        self.path_exon_dict = {}
        key = None
        with open(fname, 'r') as inp:
            for line in inp:
                if line.startswith('>'):                    
                    key = line.rstrip('\n').split('>')[-1].split('\t')[0]
                    if key in self.path_jump_dict.keys():
                        raise Exception(f"repeated entries {key}")
                    else:
                        self.path_jump_dict[key] = []
                        self.path_abundances_dict[key] = float(line.rstrip('\n').split(':')[1])
                else:
                    if not key:
                        raise Exception(f"no identifier found at line: {line}")
                    edgetype = line.rstrip('\n').split('\t')[-1]
                    if edgetype == 'splice':
                        left = int(line.split('\t')[0].split(',')[-1].rstrip(']'))
                        right = int(line.split('\t')[1].split(',')[0].lstrip('['))
                        self.path_jump_dict[key].append((left, right))

        # remove 0 weight paths
        remove_keys = [k for k, v in self.path_abundances_dict.items() if v == 0]
        for key in remove_keys:
            del self.path_abundances_dict[key]
            del self.path_jump_dict[key]

        for key, val in self.path_jump_dict.items():
            self.path_exon_dict[key] = []
            if len(val) == 0:
                self.path_exon_dict[key].append((1, 29903))
            else:
                for idx in range(len(val)):
                    if idx == 0:
                        self.path_exon_dict[key].append((1,val[idx][0]))
                    else:
                        self.path_exon_dict[key].append((val[idx-1][-1], val[idx][0]))
                self.path_exon_dict[key].append((val[-1][-1], 29903))
        self.path_length_dict = {}

        for key, exons in self.path_exon_dict.items():
            self.path_length_dict[key] = 0
            for exon in exons:
                self.path_length_dict[key] += exon[1] - exon[0] + 1

    def read_gtf_file(self, fname):
        
        self.path_jump_dict = {}
        self.path_abundances_dict = {}
        self.path_exon_dict = {}

        df = pd.read_csv(fname, sep='\t', header=None, comment='#')

        path_rpkm_dict = {}
        transcript_id = 0
        for _, row in df.iterrows():
            if row[2] == 'transcript':
                key = f"transcript_{transcript_id}"
                transcript_id += 1

                entry_dict = {entry.lstrip(' ').split(' ')[0]: entry.lstrip(' ').split(' ')[-1][1:-1] 
                              for entry in row[8].split(';')}

                if 'RPKM' in entry_dict.keys():     
                    path_rpkm_dict[key] = float(entry_dict['RPKM'])
                elif 'TPM' in entry_dict.keys():
                    path_rpkm_dict[key] = float(entry_dict['TPM'])
                else:
                    raise Exception(f"TPM or RKPM not found in gtf file in the following row\n {row}")
                self.path_exon_dict[key] = []
            else:
                self.path_exon_dict[key].append((int(row[3]), int(row[4])))
                
        total_rpkm = sum(path_rpkm_dict.values())
        self.path_abundances_dict = {key: val / total_rpkm for key, val in path_rpkm_dict.items()}
        
        for key, val in self.path_exon_dict.items():
            self.path_jump_dict[key] = []
            for idx in range(len(val)-1):
                self.path_jump_dict[key].append((val[idx][-1], val[idx+1][0]))         

        self.path_length_dict = {}
        for key, exons in self.path_exon_dict.items():
            self.path_length_dict[key] = 0
            for exon in exons:
                self.path_length_dict[key] += exon[1] - exon[0] + 1

    def top_k_jump_dict(self, kdx, length = 29903):

        chosen_keys = [key for key, path_jumps in sorted(self.path_jump_dict.items(), 
            key = lambda x: (-self.path_abundances_dict[x[0]] * self.path_length_dict[x[0]], -self.path_length_dict[x[0]]))][:kdx]

        return {k: self.path_jump_dict[k] for k in chosen_keys}
    
    def top_k_jump_list(self, kdx, length = 29903):

        chosen_keys = [key for key, path_jumps in sorted(self.path_jump_dict.items(), 
            key = lambda x: (-self.path_abundances_dict[x[0]] * self.path_length_dict[x[0]], -self.path_length_dict[x[0]]))][:kdx]

        return [ self.path_jump_dict[k] for k in chosen_keys]

    def top_k_perc_jump_dict(self, k, length = 29903):

        kdx = max(1, math.ceil(len(self.path_jump_dict) * k / 100))
        chosen_keys = [key for key, path_jumps in sorted(self.path_jump_dict.items(), 
            key = lambda x: (-self.path_abundances_dict[x[0]] * self.path_length_dict[x[0]], -self.path_length_dict[x[0]]))][:kdx]

        return {k: self.path_jump_dict[k] for k in chosen_keys}

class evaluator():
    
    def __init__(self, ground_fname):
        
        self.read_ground_truth_file(ground_fname)
        self.solutions = {}
        self.alternate_solutions = {}

    def read_ground_truth_file(self, fname):

        self.path_jump_dict = {}
        path_counts = {}
        self.path_exon_dict = {}
        self.path_length_dict = {}
        key = None
        with open(fname, 'r') as inp:
            for line in inp:
                if line.startswith('>'):                    
                    key = line.rstrip('\n').split('>')[-1].split(':')[0]
                    if key in self.path_jump_dict.keys():
                        raise Exception(f"repeated entries {key}")
                    else:
                        self.path_jump_dict[key] = []
                        path_counts[key] = float(line.rstrip('\n').split(':')[1])
                else:
                    if not key:
                        raise Exception(f"no identifier found at line: {line}")
                    edgetype = line.rstrip('\n').split('\t')[-1]
                    if edgetype == 'splice':
                        left = int(line.split('\t')[0].split(',')[-1].rstrip(']'))
                        right = int(line.split('\t')[1].split(',')[0].lstrip('['))
                        self.path_jump_dict[key].append((left, right))

        for key, val in self.path_jump_dict.items():
            self.path_exon_dict[key] = []
            for idx in range(len(val)):
                if idx == 0:
                    self.path_exon_dict[key].append((1,val[idx][0]))
                else:
                    self.path_exon_dict[key].append((val[idx-1][-1], val[idx][0]))
            if len(val) > 0:
                self.path_exon_dict[key].append((val[-1][-1], 29903))
            else:
                self.path_exon_dict[key].append((1, 29903))

        total_counts = sum(path_counts.values())
        self.path_abundances_dict = {key: val / total_counts for key, val in path_counts.items()}
        
        self.path_length_dict = {}
        for key, exons in self.path_exon_dict.items():
            self.path_length_dict[key] = 0
            for exon in exons:
                self.path_length_dict[key] += exon[1] - exon[0] + 1

    def top_k_jump_dict(self, kdx, length = 29903):

        chosen_keys = [key for key, path_jumps in sorted(self.path_jump_dict.items(), 
            key = lambda x: (-self.path_abundances_dict[x[0]] * self.path_length_dict[x[0]], -self.path_length_dict[x[0]]))][:kdx]

        return {k: self.path_jump_dict[k] for k in chosen_keys}

    def top_k_jump_list(self, kdx, length = 29903):

        chosen_keys = [key for key, path_jumps in sorted(self.path_jump_dict.items(), 
            key = lambda x: (-self.path_abundances_dict[x[0]] * self.path_length_dict[x[0]], -self.path_length_dict[x[0]]))][:kdx]

        return [ self.path_jump_dict[k] for k in chosen_keys]

    def getCaonicalCount(self, ref, tolerance = 6):
        canonical = 0
        noncanonical = 0
        for key, jumps in self.path_jump_dict.items():
            if evaluator.getORF(ref, jumps, tolerance) == 'non-canonical':
                noncanonical += 1
            else:
                canonical += 1

        return canonical, noncanonical 

    @staticmethod
    def getORF(ref, jumps, tolerance = 6):
        
        sgRNAs_from_pos = {
            21562: "S",
            25392: "3a",
            26244: "E",
            26522: "M",
            27201: "6",
            27393: "7a",
            27755: "7b",
            27893: "8",
            28273: "N",
            29557: "10"
            }        
        
        if len(jumps) == 1:
            left_pos = jumps[0][0]
            right_pos = jumps[0][1]
            
            codon_pos = ref.fetch(ref.references[0], right_pos, ref.lengths[0]).find('ATG') + right_pos
            
            if 55 - tolerance <= left_pos <= 85 + tolerance and codon_pos in sgRNAs_from_pos.keys():
                return sgRNAs_from_pos[codon_pos]
            else:
                return 'non-canonical'
        else:
            return 'non-canonical'        


    def add_solution(self, sol_id, fname, matching_fname = None):
        if sol_id not in self.solutions.keys():
            self.solutions[sol_id] = solution(self.graph, fname, matching_fname)
        else:
            raise Exception(f"sol_id {sol_id} already exists!")

    def add_alternate_solution(self, sol_id, fname, gtf = False):
        if sol_id not in self.alternate_solutions.keys():
            self.alternate_solutions[sol_id] = alternate_solution(fname, gtf)
        else:
            raise Exception(f"alternate_sol_id {sol_id} already exists!")

    @staticmethod
    def check_match(source, target, tolerance = 10):
        
        if len(source) == len(target):
            flag = True
            for idx in range(len(source)):
                source_left = source[idx][0]
                source_right = source[idx][1]
                
                target_left = target[idx][0]
                target_right = target[idx][1]
                
                if abs(source_left - target_left) > tolerance or abs(source_right - target_right) > tolerance:
                    flag = False
                    break
            
            return flag
        else:
            return False

    @staticmethod
    def comparedtranscript_abundance(source_dict, source_abundance, 
                                     target_dict, target_abundance, tolerance = 10):
        
        match_error = 0
        unmatch_error = 0
        unmatched_key_set = set(target_dict.keys())
        #print(f"{len(unmatched_key_set)}")
        error_dict = target_abundance.copy()
        for source_key, source_transcript in source_dict.items():
            match_key_list = []
            
            for target_key, target_transcript in target_dict.items():
                if evaluator.check_match(source_transcript, target_transcript, tolerance):
                    match_key_list.append(target_key)
            
            nmatches = len(match_key_list)
            if nmatches >= 1:
                for target_key in match_key_list:
                    #match_error += (target_abundance[target_key] - source_abundance[source_key] / nmatches)**2
                    error_dict[target_key] -= source_abundance[source_key] / nmatches
            unmatched_key_set = unmatched_key_set.difference(match_key_list)

        #print('-'*50)
        #print(f"{len(unmatched_key_set)}")
        #print('-'*50)

        for target_key in unmatched_key_set:
            unmatch_error += target_abundance[target_key]**2

        return match_error, unmatch_error, sum([val**2 for val in error_dict.values()])

    @staticmethod
    def unique_transcripts(source_dict, target_dict, tolerance = 10):

        unmatched_key_set = set(target_dict.keys())
        for source_key, source_transcript in source_dict.items():
            for target_key, target_transcript in target_dict.items():
                if evaluator.check_match(source_transcript, target_transcript, tolerance):
                    unmatched_key_set = unmatched_key_set.difference([target_key])
                    break

        return len(target_dict.keys()) - len(unmatched_key_set)

    @staticmethod
    def unique_canonical_noncanonical_from_jumps(source_dict, target_dict, ref, tolerance = 10):

        unmatched_key_set = set(target_dict.keys())
        for source_key, source_transcript in source_dict.items():
            for target_key, target_transcript in target_dict.items():
                if evaluator.check_match(source_transcript, target_transcript, tolerance):
                    unmatched_key_set = unmatched_key_set.difference([target_key])
                    break

        canonical_count = 0
        noncanonical_count = 0
        for target_key, target_transcript in target_dict.items():
            if target_key not in unmatched_key_set:
                curr_orf = evaluator.getORF(ref, target_transcript)
                if curr_orf == 'non-canonical':
                    noncanonical_count += 1
                else:
                    canonical_count += 1

        return canonical_count, noncanonical_count

    @staticmethod
    def compare_transcripts(source, target, tolerance = 10):
        
        true_pos = 0
        false_pos = 0
        
        for source_transcript in source:
            
            flag = False
            
            for target_transcript in target:                
                if evaluator.check_match(source_transcript, target_transcript, tolerance):
                    true_pos += 1
                    flag = True
                    break
            if not flag:
                false_pos += 1
        
        return true_pos, false_pos

    def get_false_negative_exons(self, source, tolerance = 10):
        
        false_neg = 0
        for target_transcript in self.path_exon_dict.values():
            flag = False
            for source_transcript in source:
                if check_match(target_transcript, source_transcript, tolerance):
                    flag = True
                    break
            if not flag:
                false_neg += 1

        return false_neg
