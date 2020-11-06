#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  2 00:58:51 2020

@author: Palash Sashittal
"""

import pysam
import pandas as pd
import sys
import argparse
import statistics
import math

import bam_io
from solveTranscription import solveTranscription
from segment_graph_aux import *

class mynode():
            
    def __init__(self, name = None, left = None, right = None, length = 29903):
        if left is None and right is None:
            self.name = None
            self.left = None
            self.right = None
        elif left is None and right is not None:
            self.left = 1
            self.right = right
            self.name = f"[start,{right}]"
        elif right is None and left is not None:
            self.left = left
            self.right = length
            self.name = f"[{left},end]"
        else:
            assert(left <= right), "left must be <= right in mynode()"
            self.name = f"[{left},{right}]"
            self.left = left
            self.right = right

        if name is not None:
            self.name = name
        
class myedge():
    
    def __init__(self, left_node = None, right_node = None, weight = None, name = None, edgetype = None):
        
        if left_node is None and right_node is None:
            self.left = None
            self.right = None
            self.weight = None
            self.name = None
        elif left_node is not None and right_node is not None:
            assert(left_node.right <= right_node.left), "nodes must be ordered [myedge()]"
            self.left = left_node
            self.right = right_node
            self.weight = None
            if edgetype is not None:
                self.name = f"{left_node.name}->{edgetype}->{right_node.name}"
            else:
                self.name = f"{left_node.name}->{right_node.name}"
        else:
            raise Exception("either both nodes must be defined or neither! [myedge()]")
        
        if name is not None:
            self.name = name
        if weight is not None:
            self.weight = weight
        if type is not None:
            self.type = edgetype

class segmentGraph():
    
    def __init__(self, breakpoints = None, edges = None, phasing = None, ref = None, verbose = True, edge_threshold = 0, phasing_threshold = 0, width=0, threads = 1):

        breakpoints = sorted(list(breakpoints))
        
        self.nodes = []
        if len(breakpoints) > 0:
            self.nodes.append(mynode(right=breakpoints[0], length = ref.lengths[0]))
            for idx in range(len(breakpoints)-1):
                self.nodes.append(mynode(left=breakpoints[idx],
                                         right=breakpoints[idx+1],
                                         length = ref.lengths[0]))
            if len(breakpoints) > 1:
                self.nodes.append(mynode(left=breakpoints[-1], length = ref.lengths[0]))
        
        self.edges = []
        key2edge = {}
        self.edge2Index = {}
        edgeIndex = 0
        if len(edges) > 0:
            for key, value in edges.items():
                left_idx = breakpoints.index(key[0])
                right_idx = breakpoints.index(key[1])
                for edgetype, count in value.items():
                    edge = myedge(left_node = self.nodes[left_idx], 
                                  right_node = self.nodes[right_idx + 1],
                                  weight = count, edgetype = edgetype)
                    self.edges.append(edge)

                    if edgetype == 'splice':
                        key2edge[(key[0], key[1])] = edge
                    self.edge2Index[edge] = edgeIndex
                    edgeIndex += 1

        if phasing is not None:
            self.phasing = []
            idx = 0
            for key, count in phasing.items():
                # print(f"{key} -- {count}")
                if count > phasing_threshold:
                    splus = set()
                    sminus = set()

                    for left, right in key[0]:
                        splus.add(key2edge[(left, right)])
                    for left, right in key[1]:
                        sminus.add(key2edge[(left, right)])

                    if len(splus.intersection(sminus)) > 0:
                        print(f"error phase {idx} {count}")
                        print(f"splus edges --")
                        for edge in splus:
                            print(f"{edge.name}")
                        print(f"sminus edges --")
                        for edge in sminus:
                            print(f"{edge.name}")
                    assert(len(splus.intersection(sminus)) == 0)
                    
                    # remove unncessary sminus edges
                    for splus_edge in splus:
                        removeEdges = []
                        for sminus_edge in sminus:
                            if isoverlap(sminus_edge, splus_edge):
                                removeEdges.append(sminus_edge)
                        for edge in removeEdges:
                            sminus.remove(edge)
                    
                    #if len(splus) >= 1 or len(sminus) >= 1:
                    self.phasing.append((list(splus), list(sminus), count))
                    idx += 1
        else:
            self.phasing = None
                        
        self.ref = ref
        self.verbose = verbose
        self.edge_threshold = edge_threshold 
        self.width = width

        self.numNodes = len(self.nodes)
        self.numEdges = len(self.edges)
        self.numSpliceEdges = len([edge for edge in self.edges if edge.type == 'splice'])
        self.transcripts = None
        if self.phasing:
            self.numPhasing = len(self.phasing)
        else:
            self.numPhasing = 0
        self.deviation = None
        self.total_weight = None
        self.nonTranscriptSpliceEdges = None
        self.threads = threads
        self.transcriptSpliceEdges = None
        self.numTranscripts = None
        
    def getBreakpoints(self) -> List[int]:
        
        breakpoints = []
        for node in self.nodes[:-1]:
            if node.right:
                breakpoints.append(node.right)

        return breakpoints
    
    def getEdges(self) -> Dict[Tuple[int, int], Dict[str, int]]:
        
        edgeDict = {}
        for edge in self.edges:
            key = (edge.left.right, edge.right.left)
            if key not in edgeDict.keys():
                edgeDict[key] = {}
            edgeDict[key][edge.type] = edge.weight            
        
        return edgeDict
    
    def getTotalFlow(self) -> float:
        
        return sum([edge.weight for edge in self.outArcs(self.nodes[0])])
    
    def getPhasing(self) -> List[Tuple[int, int, str]]:
        
        phasingList = []
        for phase in self.phasing:
            curr_phase = []
            for edge in phase:
                curr_phase.append((edge.left.right, edge.right.left, edge.type))
            phasingList.append(curr_phase)
            
        return phasingList
        
    def __sub__(self, other: 'segmentGraph') -> 'segmentGraph':
                
        breakpoints = self.getBreakpoints()
        if not breakpoints == other.getBreakpoints():
            raise Exception('the breakpoints do not match')
                
        phasing = self.getPhasing()
        # if not phasing == other.getPhasing():
        #     raise Exception('the phasing do not match')
        
        edgeDict = self.getEdges()
        keyList = [(edge.left.right, edge.right.left, edge.type) for edge in self.edges]
        # edges2 = other.getEdges()
        
        # for key, val in edgeDict.items():
        #     for mytype, weight in val.items():
        #         print(f"{key[0]}->{key[1]} -- {mytype} : {weight} with {type(weight)}")
        
        for edge in other.edges:
            edgeKey = (edge.left.right, edge.right.left)
            key = (edge.left.right, edge.right.left, edge.type)
            if key in keyList:
                # print(f"checking {key} with {edge.type} and {edge.weight}")
                edgeDict[edgeKey][edge.type] -= edge.weight
            else:
                if edgeKey in edgeDict.keys():
                    edgeDict[edgeKey][edge.type] = -edge.weight
                else:
                    edgeDict[edgeKey] = {edge.type: -edge.weight}            
        
        return segmentGraph(breakpoints = breakpoints, edges = edgeDict, phasing = phasing, ref = self.ref,
                           verbose = self.verbose, threads = self.threads)

    def copy(self) -> 'segmentGraph':        
        return segmentGraph(breakpoints = self.getBreakpoints(), edges = self.getEdges(), phasing = self.getPhasing(), 
                           ref = self.ref, verbose = self.verbose, threads = self.threads)

    def outArcs(self, node, index = False) -> myedge:
        for idx, edge in enumerate(self.edges):
            if edge.left == node:
                if index:
                    yield idx, edge
                else:   
                    yield edge

    def inArcs(self, node, index = False) -> myedge:
        for idx, edge in enumerate(self.edges):
            if edge.right == node:
                if index:
                    yield idx, edge
                else:
                    yield edge

    def allArcs(self, node, index = False) -> myedge:
        for idx, edge in enumerate(self.edges):
            if edge.left == node or edge.right == node:
                if index:
                    yield idx, edge
                else:
                    yield edge

    def countOutArcs(self, node):
        count = 0
        for edge in self.edges:
            if edge.left == node:
                count += 1
        return count
    
    def countInArcs(self, node):
        count = 0
        for edge in self.edges:
            if edge.right == node:
                count += 1
        return count    

    def writeGraph(self, fname) -> None:
        
        with open(fname, 'w') as output:
            output.write(f">nodes\t{self.numNodes}\n")
            for node in self.nodes:
                output.write(f"{node.name}\n")
            output.write(f">edges\t{self.numEdges}\n")
            for edge in self.edges:
                output.write(f"{edge.left.name}\t{edge.right.name}\t{edge.type}\t{edge.weight}\n")

    def writeDOT(self, fname) -> None:
        
        edgeColor = {'splice': 'black',
                     'A': 'blue',
                     'C': 'red',
                     'T': 'yellow',
                     'G': 'green',
                     'adjacency': 'black'}
        
        with open(fname, 'w') as output:
        
            output.write(f"digraph N {{\n")
            output.write(f"\toverlap=\"false\"\n")
            output.write(f"\trankdir=\"LR\"\n")
            
            idx_dict = {}
            idx = 0
            for node in self.nodes:
                idx_dict[node.name] = idx
                output.write(f"\t{idx} [label=\"{node.name}\", style=\"bold\"];\n")
                idx += 1
            
            for edge in self.edges:
                left_idx = idx_dict[edge.left.name]
                right_idx = idx_dict[edge.right.name]
                #if edge.weight >0 :
                output.write(f"\t{left_idx} -> {right_idx} [label=\"{edge.weight}\", fillcolor={edgeColor[edge.type]}, style=\"bold\"];\n")
                # else:
                #     output.write(f"\t{left_idx} -> {right_idx} [fillcolor={edgeColor[edge.type]}, style=\"bold\"];\n")

            output.write(f"}}")
            
    def getNumSTpaths(self) -> int:
        
        nodeCount = {}
        for idx, node in enumerate(self.nodes):
            if idx == 0:
                nodeCount[node] = 1
            else:
                nodeCount[node] = 0
                for edge in self.inArcs(node):
                    nodeCount[node] += nodeCount[edge.left]

        return nodeCount[self.nodes[-1]]

    def stitchPath(self, path_spliceEdges):
        path_spliceEdges.sort(key = lambda x: x.left.right)

        path = []
        # piece together the path
        left_id = 0
        for edge in path_spliceEdges:
            right_id = self.nodes.index(edge.left)
            subpath = self.getCanonicalPaths(left_id, right_id)
            if len(subpath) > 1:
                raise Exception(f"more than 1 canonical subpath for {left_id} -> {right_id}")
            if len(subpath) == 1:
                for subedge in subpath[0]:
                    #path.append((subedge.left.name, subedge.right.name, subedge.type))
                    path.append(subedge)
            #path.append((edge.left.name, edge.right.name, edge.type))
            path.append(edge)
            left_id = self.nodes.index(edge.right)

        if left_id < self.numNodes - 1:
            subpath = self.getCanonicalPaths(left_id, self.numNodes - 1)
            if len(subpath) > 1:
                raise Exception(f"more than 1 canonical subpath for {left_id} -> {right_id}")                    
            for subedge in subpath[0]:
                #path.append((subedge.left.name, subedge.right.name, subedge.type))
                path.append(subedge)
        return path

    def getPathIndex(self, path):
        pathIndex = 0
        count = 0
        for idx, edge in enumerate(self.transcriptSpliceEdges):
            if edge in path:
                pathIndex += idx + 1
                count += 1
        assert(count <= 1)
        nTranscriptSpliceEdges = len(self.transcriptSpliceEdges)
        for idx, edge in enumerate(self.nonTranscriptSpliceEdges):
            if edge in path:
                pathIndex += (nTranscriptSpliceEdges + 1) * (2**(idx + 1))

        return pathIndex
    
    def getPathLength(self, path):
        
        length = 29903
        for edge in path:
            if edge.type == 'splice':
                length -= (edge.right.left - edge.left.right)
        
        return length
    
    def writePaths(self, paths, fname, info = [None]):
        if all(info):
            assert(len(paths) == len(info))
            pathTuples = [(self.getPathIndex(path), path, info[idx]) for idx, path in enumerate(paths)]
        else:
            pathTuples = [(self.getPathIndex(path), path) for path in paths]
        pathTuples.sort(key = lambda x: x[0])
        with open(fname, 'w') as output:
            for element in pathTuples:
                pathIndex = element[0]
                path = element[1]
                if all(info):
                    if pathIndex <= self.numTranscripts:
                        output.write(f">transcript_{pathIndex}: {element[2]}\n")
                    else:
                        output.write(f">path_{pathIndex}: {element[2]}\n")
                else:
                    if pathIndex <= self.numTranscripts:
                        output.write(f">transcript_{pathIndex}:\n")
                    else:   
                        output.write(f">path_{pathIndex}:\n")
                for edge in path:
                    output.write(f"{edge.left.name}\t{edge.right.name}\t{edge.type}\n")

    def writeFasta(self, paths, fname, info = [None]):
        if all(info):
            assert(len(paths) == len(info))
            pathTuples = [(self.getPathIndex(path), path, info[idx]) for idx, path in enumerate(paths)]
        else:
            pathTuples = [(self.getPathIndex(path), path) for path in paths]
        pathTuples.sort(key = lambda x: x[0])
        with open(fname, 'w') as output:
            for element in pathTuples:
                pathIndex = element[0]
                path = element[1]
                if all(info):
                    if pathIndex <= self.numTranscripts:
                        output.write(f">transcript_{pathIndex}: {element[2]}\n")
                    else:
                        output.write(f">path_{pathIndex}: {element[2]}\n")
                else:
                    if pathIndex <= self.numTranscripts:
                        output.write(f">transcript_{pathIndex}:\n")
                    else:   
                        output.write(f">path_{pathIndex}:\n")
                genome = ''
                for edge in path:
                    if edge.left.left == 1:
                        genome += self.ref.fetch(self.ref.references[0], 0, edge.left.right)
                    else:
                        genome += self.ref.fetch(self.ref.references[0], edge.left.left, edge.left.right)
                genome += self.ref.fetch(self.ref.references[0], path[-1].right.left, path[-1].right.right)
                output.write(f"{genome}\n")


    def getTotalWeight(self) -> float:
        return sum([edge.weight for edge in self.edges])
        
    def writeReport(self, fname) -> None:
        with open(fname, 'w') as output:
            output.write(f"nodes\t{self.numNodes}\n")
            output.write(f"edges\t{self.numEdges}\n")
            output.write(f"splices\t{self.numSpliceEdges}\n")
            output.write(f"paths\t{self.getNumSTpaths()}\n")
            output.write(f"deviation\t{self.deviation}\n")
            output.write(f"total_flow\t{self.getTotalFlow()}\n")
            output.write(f"total_weigth\t{self.total_weight}\n")
            if not self.total_weight or not self.deviation:
                output.write(f"fraction\t{None}\n")
            else:
                output.write(f"fraction\t{self.deviation/self.total_weight}\n")
            if not self.transcripts:
                self.initCanonicalTranscripts()
                output.write(f"transcripts\t{len(self.transcripts)}\n")
            
            inbyout = []
            for node in self.nodes[1:-1]:
                total_out = sum([edge.weight for edge in self.outArcs(node)])
                total_in = sum([edge.weight for edge in self.inArcs(node)])
                inbyout.append(total_out / total_in)
            output.write(f"median_in_out\t{statistics.median(inbyout)}\n")
            output.write(f"mean_in_out\t{statistics.mean(inbyout)}\n")
        
    def edgesBetween(self, leftNode, rightNode) -> myedge:
        
        for edge in self.outArcs(leftNode):
            if edge.right == rightNode:
                yield edge
        
    def getCanonicalPaths(self, startIndex, endIndex):
        
        if startIndex > endIndex:
            raise Exception(f"start index {startIndex} must be less than or equal to index {endIndex}!")

        paths = []
        if startIndex == endIndex:
            return paths
        
        lastNode = self.nodes[endIndex]
        penultimateNode = self.nodes[endIndex-1]                
        if startIndex == endIndex - 1:
            for edge in self.edgesBetween(penultimateNode, lastNode):
                paths.append([edge])
        else:
            for path in self.getCanonicalPaths(startIndex, endIndex - 1):
                for edge in self.edgesBetween(penultimateNode, lastNode):
                    paths.append(path + [edge])
        return paths
        
    def initCanonicalTranscripts(self, ref = None):
        
        if not ref:
            if not self.ref:
                raise Exception('cannot initialize canonical transcripts without reference!')
            else:
                ref = self.ref
        else:
            self.ref = ref
        
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
        
        spliceTypeList = []
        self.nonTranscriptSpliceEdges = []
        self.transcriptSpliceEdges = []
        for edge in self.edges:
            left_pos = edge.left.right
            right_pos = edge.right.left
            
            codon_pos = ref.fetch(ref.references[0], right_pos, ref.lengths[0]).find('ATG') + right_pos
            
            if 55 - self.width <= left_pos <= 85 + self.width and codon_pos in sgRNAs_from_pos.keys():
                self.transcriptSpliceEdges.append(edge)
                spliceTypeList.append(sgRNAs_from_pos[codon_pos])
            else:
                if edge.type == 'splice':
                    self.nonTranscriptSpliceEdges.append(edge)
            # else:
            #     print(f"failed ({left_pos} -> {right_pos}) with codon_pos as {codon_pos}")
        #self.transcriptSpliceEdges = spliceList.copy()

        # self.transcripts = []
        # self.transcript_type = []
        self.transcripts = self.getCanonicalPaths(0, self.numNodes - 1)
        assert(len(self.transcripts) == 1)
        self.transcript_type = ["1ab"]
        
        for idx, edge in enumerate(self.transcriptSpliceEdges):
            _5primeIndex = self.nodes.index(edge.left)
            _3primeIndex = self.nodes.index(edge.right)
            
            leftpaths = self.getCanonicalPaths(0, _5primeIndex)
            rightpaths = self.getCanonicalPaths(_3primeIndex, self.numNodes - 1)
            
            if len(leftpaths) > 0 and len(rightpaths) > 0:
                for leftpath in leftpaths:
                    for rightpath in rightpaths:        
                        self.transcripts.append(leftpath + [edge] + rightpath)
                        # self.transcript_type.append(spliceTypeList[idx])
            elif len(leftpaths) > 0 and len(rightpaths) == 0:
                for leftpath in leftpaths:
                    self.transcripts.append(leftpath + [edge])
                    # self.transcript_type.append(spliceTypeList[idx])
            elif len(leftpaths) == 0 and len(rightpaths) > 0:
                for rightpath in rightpaths:
                    self.transcripts.append([edge] + rightpath)
                    # self.transcript_type.append(spliceTypeList[idx])     
            else:
                self.transcripts.append([edge])
                # self.transcript_type.append(spliceTypeList[idx])
            self.transcript_type.append(spliceTypeList[idx])

        # self.transcripts += self.getCanonicalPaths(0, self.numNodes - 1)
        # self.transcript_type.append("1ab")
        self.numTranscripts = len(self.transcripts)
        
        print(f"number of transcripts is {len(self.transcripts)}")
        print(f"number of transcript splice edges is {len(self.transcriptSpliceEdges)}")
        print('-'*50)

    def coarsenEdgeWeights(self, resolution):
        for edge in self.edges:
            edge.weight = int(math.floor(float(edge.weight) / float(resolution)))

    def getMaxNonCanonicalFlow(self) -> float:
        if not self.nonTranscriptSpliceEdges:
            self.initCanonicalTranscripts()
        
        return sum([edge.weight for edge in self.nonTranscriptSpliceEdges])

    def writeTranscripts(self, fname):
        with open(fname, 'w') as output:
            for idx, transcript in enumerate(self.transcripts):
                output.write(f">transcript_{idx}\t{self.transcript_type[idx]}\n")
                for edge in transcript:
                    output.write(f"{edge.left.name}\t{edge.right.name}\t{edge.type}\n")

    def getTranscriptsFromPaths(self, paths):        
        return None
    
    def prune_edges(self):        
        self.edges = [edge for edge in self.edges if edge.weight >= self.edge_threshold or edge.type != 'splice']
        
    def prune_nodes(self):
        return None
    
    def prune(self):
        self.prune_edges()
        self.prune_nodes()

def main(args):
    
    ref = pysam.FastaFile(args.fasta)    

    SJcounter = bam_io.SpliceJunction(args.fasta, min_base_qual=args.min_base_qual, 
                                      min_mapping_qual=args.min_mapping_qual,
                                      contig = ref.references[0])
    
    if args.inputBreakpoints and args.inputEdges:

        breakpoints, splice_edges, adjacency_edges = readGraphSkeleton(args.inputBreakpoints, args.inputEdges)        
    else:
        
        if args.inputCSV:
            df_sj_reads = pd.read_csv(args.inputCSV, sep='\t')
        else:
            df_sj_reads = SJcounter.get_sj_reads(args.bam, paired = args.paired)
    
        breakpoints, splice_edges, adjacency_edges = getGraphSkeleton(df_sj_reads, ref, args.width, args.sj_threshold)

    # for idx, (left, right) in enumerate(splice_edges.keys()):
    #     print(f"{idx} -- {left} -- {right}")

    if args.inputPhasing:
        phasing_reads = readPhasing(args.inputPhasing)
    else:
        phasing_reads = SJcounter.get_phasing_reads(args.bam, splice_edges.keys(), args.width, args.paired)
    
    if args.verbose:
        print(f"number of splice edges -- {len(splice_edges)}")
        print(f"number of phasing reads -- {len(phasing_reads)}")
        print(f"total number of reads -- {sum([val for val in phasing_reads.values()])}")
    
    if args.outputCSV:
        df_sj_reads.to_csv(args.outputCSV, sep='\t', index=False)
    
    
    sorted_breakpoints = sorted(list(breakpoints))
    if args.outputBreakpoints:
        writeBreakpoints(args.outputBreakpoints, sorted_breakpoints)

    if args.outputPhasing:
        writePhasing(args.outputPhasing, phasing_reads)

    if args.outputEdges:
        writeEdges(args.outputEdges, [splice_edges, adjacency_edges])

    sampleGraph = segmentGraph(breakpoints = breakpoints, edges = {**splice_edges, **adjacency_edges}, 
                              phasing = phasing_reads, verbose = args.verbose, phasing_threshold=args.phasing_threshold,
                              edge_threshold = args.sj_threshold, ref = ref, width=args.width)

    if args.verbose:
        print(f"number of ST paths -- {sampleGraph.getNumSTpaths()}")

    if args.outputDOT:
        sampleGraph.writeDOT(args.outputDOT)
        
    if args.outputGraph:
        sampleGraph.writeGraph(args.outputGraph)
    
    if args.outputTranscripts:
        sampleGraph.initCanonicalTranscripts(ref = ref)
        if args.outputTranscripts:
            sampleGraph.writeTranscripts(args.outputTranscripts)

    if args.outputDecomposition:
            
        solver = solveTranscription(graph = sampleGraph, numPaths = args.numPaths, verbose = args.verbose, 
                                    threads = args.threads, timelimit = args.timelimit, max_iter = args.maxIter,
                                    samplingFrequency = args.samplingFrequency, ref_length = ref.lengths[0])

        if solver.solve(args.greedy):

            if args.cov_threshold > 0:
                solver.filterSolution(args.cov_threshold)
            
            solver.writeSolution(args.outputDecomposition)
            
            if args.outputMatching:
                solver.writeMatching(args.outputMatching)
            
            if args.outputGTF:
                solver.writeGTF(args.outputGTF, ref.references[0])
        else:
            with open(args.outputDecomposition, 'w') as output:
                output.write('no solution!')

def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bam", type=str, help="aligned bam file")
    parser.add_argument("-f", "--fasta", type=str, help="fasta file", required=True)
    parser.add_argument("--min-base-qual", type=int, help="minimum base quality [20]", default=20)
    parser.add_argument("--min-mapping-qual", type=int, help="minimum mapping quality [20]", default=20)
    parser.add_argument("--min-snv-qual", type=int, help = "minimum snv quality [30]", default = 30)
    parser.add_argument("--sj_threshold", type=int, help="minimum support for splicing junction [20]", default=20)
    parser.add_argument("--outputCSV", type=str, help="output csv file for sj reads")
    parser.add_argument("--outputPhasing", type=str, help="output file containing phasing reads")
    parser.add_argument("--inputCSV", type=str, help="input csv file with sj reads")
    parser.add_argument("--inputPhasing", type=str, help="input phasing file")
    parser.add_argument("--outputGraph", type=str, help="output graph file")
    parser.add_argument("--outputDOT", type=str, help="output DOT file for splice graph")
    parser.add_argument("--report", type=str, help="output file for report on the splice graph")
    parser.add_argument("-w", "--width", type=int, help="spliced junction width parameter [6]", default=6)
    parser.add_argument("--outputTranscripts", type = str, help="output file for transcripts")
    parser.add_argument("-k", type=int, dest='numPaths', help="number of paths for the flow decomposition", default = 1)
    parser.add_argument("--outputBreakpoints", type=str, help="output file containing breakpoints")
    parser.add_argument("--inputBreakpoints", type=str, help="input file containing breakpoints")
    parser.add_argument("--outputEdges", type=str, help="output file containing graph edges")
    parser.add_argument("--inputEdges", type=str, help="input file containing graph edges")
    parser.add_argument("--outputDecomposition", type=str, help="output file for the decomposed non-canonical transcripts")
    parser.add_argument("--outputMatching", type=str, help="output file for the matching of phasing reads to inferred transcripts")    
    parser.add_argument("--noverbose", dest="verbose", help="do not output statements from internal solvers [default is false]", action='store_false')
    parser.add_argument("--threads", type=int, default = 1, help = 'number of threads allowed to be used')
    parser.set_defaults(verbose=True)
    parser.add_argument("--samplingFrequency", type = int, help="number of sampling points for the likelihood function", default = 4)
    parser.add_argument("--timelimit", type=int, help="time limt for the gurobi solvers in seconds [None]", default = None)
    parser.add_argument("--greedy", type=str2bool, default=False, help="set greedy flag to TRUE")
    parser.add_argument("--maxIter", type=int, help="maximum iterations for the greedy algorithm [100]", default = 100)
    parser.add_argument("--outputGTF", type=str, help="output file in GTF format")
    parser.add_argument("--phasing_threshold", type=int, help='coverage threshold for transcripts [0]', default = 0)

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    args.vcf = None
    args.paired = True
    main(args)
