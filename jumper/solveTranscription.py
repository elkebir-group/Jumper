#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 20:13:46 2020

@author: Palash Sashittal
"""

import gurobipy as gp
import math

class solveTranscription():
    
    def __init__(self, graph = None, numPaths = 1, verbose = True, ref = None,
                 samplingFrequency = 10, threads = 1, timelimit = None, max_iter = 100,
                 ref_length = 29903):

        self.graph = graph
        self.numPaths = numPaths
        self.verbose = verbose
        self.ref = ref

        self.transcriptWeights = None
        self.solution = None
        self.chosen_paths = None
        self.supported_phasing_read_indices = None
        self.threads = threads
        self.timelimit = timelimit
        self.max_iter = max_iter
        self.ref_length = ref_length
        if samplingFrequency <= 1:
            raise Exception(f"sampling frequecny {samplingFrequency} too low for solveTranscription")
        else:
            self.samplingFrequency = samplingFrequency
    
    @staticmethod
    def find_shortest_transcript_length(spliceEdges, spliceLengths, ref_length = 29903):
        
        model = gp.Model('shortestPathProblem')
        
        nsplice = len(spliceEdges)
        assert(nsplice == len(spliceLengths))
        
        w = model.addVars(nsplice, lb = 0, vtype = gp.GRB.BINARY, name = 'w')
        
        for idx1, edge1 in enumerate(spliceEdges):
            for idx2, edge2 in enumerate(spliceEdges):
                u1 = edge1.left.right
                v1 = edge1.right.left
                u2 = edge2.left.right
                v2 = edge2.right.left
                #if (u1 <= u2 and v1 >= u2) or (u2 <= u1 and v2 >= u1):
                if (u1 <= u2 and v1 > u2) and (idx1 != idx2):
                    model.addConstr(w[idx1] + w[idx2] <= 1)
                
        expr_obj = gp.LinExpr()
        for idx in range(nsplice):
            expr_obj += spliceLengths[idx] * w[idx]

        model.setObjective(expr_obj, gp.GRB.MAXIMIZE)
        
        model.setParam(gp.GRB.Param.OutputFlag, False)
        
        model.optimize()
        
        return ref_length - model.getObjective().getValue()
    
    def solve_greedy_likelihood(self):
        
        # initialize canonical and noncanonical splice edges if not done already
        if not self.graph.nonTranscriptSpliceEdges or not self.graph.transcriptSpliceEdges:
            self.graph.initCanonicalTranscripts()

        # list of splice edges
        spliceEdges = self.graph.transcriptSpliceEdges + self.graph.nonTranscriptSpliceEdges
    
        # splice edge to index
        spliceEdge2Index = {}
        for idx, edge in enumerate(spliceEdges):
            spliceEdge2Index[edge] = idx
    
        # splice lengths
        spliceLengths = []
        for edge in spliceEdges:
            spliceLengths.append(edge.right.left - edge.left.right)


        lmin = solveTranscription.find_shortest_transcript_length(spliceEdges, spliceLengths, self.ref_length)
        self.chosen_paths = []
        
        
        if self.verbose:
            print('-'*50)
            print(f"max phasing count is {max([count for _, _, count in self.graph.phasing])}")
            print(f"total reads is {sum([count for _, _, count in self.graph.phasing])}")
            print('-'*50)
        
        
        # function values and sampling points
        F = [0]*self.samplingFrequency 
        xh = [0]*self.samplingFrequency
        delta_x = 1.0/2**(self.samplingFrequency- 1)
        for i in range(self.samplingFrequency):
            if i == 0:
                xh[i] = 0
            else:
                xh[i] = 2**(i - 1) * delta_x
            if i == 0:
                F[i] = math.log(delta_x / 100)
            else:
                F[i] = math.log(xh[i])

        obj_scaling_factor = 10**4 / max([count for _, _, count in self.graph.phasing]) / abs(F[0])

        flag = 2
        iter_count = 0
        kdx = 0
        npaths = 1
        while iter_count <= self.max_iter or flag in [1,2,3]:

            model = gp.Model(f"solveGreedyTranscription_{kdx}")
            if self.threads > 1:
                model.setParam(gp.GRB.Param.Threads, self.threads)
            if self.timelimit:
                model.setParam(gp.GRB.Param.TimeLimit, self.timelimit)

            c = model.addVars(len(self.chosen_paths), lb = 0, ub = 1, 
                              vtype = gp.GRB.CONTINUOUS, name="c")
            g = model.addVars(npaths, lb = 0, ub = 1, vtype = gp.GRB.CONTINUOUS, name = 'g')
            w = model.addVars(self.graph.numSpliceEdges, npaths, lb = 0, 
                              vtype = gp.GRB.BINARY, name = 'w')
            y = model.addVars(self.graph.numPhasing, npaths, lb = 0, ub = 1,
                              vtype = gp.GRB.CONTINUOUS, name = 'y')
            ysplice = model.addVars(self.graph.numSpliceEdges, npaths, lb = 0, ub = 1,
                                    vtype = gp.GRB.CONTINUOUS, name = 'ysplice')
            ld = model.addVars(self.graph.numPhasing, self.samplingFrequency, lb = 0, ub = 1, 
                               vtype = gp.GRB.CONTINUOUS, name = 'lambda')

            observed_splice_edges = set()
            for phase in self.graph.phasing:
                observed_splice_edges |= set(phase[0])

            # remove unseen nontranscript splice edges
            for edge in spliceEdges:
                if edge not in observed_splice_edges:
                    for pid in range(npaths):
                        model.addConstr(w[spliceEdge2Index[edge], pid] == 0)

            # constraints
            # non-overlapping splice edges
            for idx1, edge1 in enumerate(spliceEdges):
                for idx2, edge2 in enumerate(spliceEdges):
                    u1 = edge1.left.right
                    v1 = edge1.right.left
                    u2 = edge2.left.right
                    v2 = edge2.right.left
                    #if (u1 <= u2 and v1 >= u2) or (u2 <= u1 and v2 >= u1):
                    if (u1 <= u2 and v1 > u2) and (idx1 != idx2):
                        for pid in range(npaths):
                            model.addConstr(w[idx1, pid] + w[idx2, pid] <= 1)

            # transcript contribution to phasing reads
            for idx in range(self.graph.numPhasing):
                for pid in range(npaths):
                    model.addConstr(y[idx, pid] <= g[pid])
                    
                    splus = self.graph.phasing[idx][0]
                    sminus = self.graph.phasing[idx][1]
                    expr_plus = gp.LinExpr()
                    expr_minus = gp.LinExpr()

                    for edge in splus:
                        # spliceIndex = spliceEdges.index(self.graph.edges[edgeIndex])
                        spliceIndex = spliceEdge2Index[edge]
                        model.addConstr(y[idx, pid] <= w[spliceIndex, pid])
                        expr_plus += w[spliceIndex, pid]

                    for edge in sminus:
                        # spliceIndex = spliceEdges.index(self.graph.edges[edgeIndex])
                        spliceIndex = spliceEdge2Index[edge]
                        model.addConstr(y[idx, pid] <= 1 - w[spliceIndex, pid])
                        expr_minus += 1 - w[spliceIndex, pid]
                    
                    nplus = len(splus)
                    nminus = len(sminus)
                    model.addConstr(y[idx, pid] >= expr_plus + expr_minus + 
                                    g[pid] - nplus - nminus)

            # transcript contribution to splice edges
            for idx in range(self.graph.numSpliceEdges):
                for pid in range(npaths):
                    model.addConstr(ysplice[idx, pid] <= g[pid])                    
                    model.addConstr(ysplice[idx, pid] <= w[idx, pid])
                    model.addConstr(ysplice[idx, pid] >= g[pid] + w[idx, pid] - 1)

            # normalization
            ## new path
            expr_norm = gp.LinExpr()
            for pid in range(npaths):
                expr_norm +=  self.ref_length * g[pid]
                for idx in range(self.graph.numSpliceEdges):
                    expr_norm -= (spliceLengths[idx] * ysplice[idx, pid])
                        
            
            # chosen paths
            for idx, path in enumerate(self.chosen_paths):
                path_length = self.ref_length
                for edge in path:
                    if edge.type == 'splice':
                        path_length -= spliceLengths[spliceEdge2Index[edge]]
                expr_norm += c[idx] * path_length

            model.addConstr(expr_norm == lmin)


            # lambda sum to 1
            for i in range(self.graph.numPhasing):
                expr_sum = gp.LinExpr()
                for j in range(self.samplingFrequency):
                    expr_sum += ld[i, j]
                model.addConstr(expr_sum == 1)
    
            # values at which we want to evalute the objective function
            # expr_q = [gp.LinExpr()] * len(self.chosen_paths)
            for idx in range(self.graph.numPhasing):
                expr_x = gp.LinExpr()
                for j in range(self.samplingFrequency):
                    expr_x += xh[j] * ld[idx, j]
                
                expr_y = gp.LinExpr()
                for pid in range(npaths):
                    expr_y += y[idx, pid]
                
                splus = self.graph.phasing[idx][0]
                sminus = self.graph.phasing[idx][1]
                
                for pathIndex, path in enumerate(self.chosen_paths):
                    pathSpliceEdges = [edge for edge in path if edge.type == 'splice']
                    
                    if (set(splus).issubset(set(pathSpliceEdges)) and
                        len(set(pathSpliceEdges).intersection(set(sminus))) == 0):                                                
                        expr_y += c[pathIndex]

                model.addConstr(expr_x == expr_y)
                # each phasing read must be supported
                # model.addConstr(expr_x >= delta_x)


            # objective
            expr_obj = gp.LinExpr()
            total_reads = 0
            for idx in range(self.graph.numPhasing):
                nreads = self.graph.phasing[idx][2]
                total_reads += nreads
                for i in range(self.samplingFrequency):
                    expr_obj += nreads * F[i] * ld[idx, i] * obj_scaling_factor

            model.setObjective(expr_obj, gp.GRB.MAXIMIZE)
            # else:
            #     model.modelSense = gp.GRB.MAXIMIZE
            #     model.setObjectiveN(expr_obj, 0, 1)

            #     expr_obj_splice_count = gp.LinExpr()
            #     for idx, path in enumerate(self.chosen_paths):
            #         spliceCount = len([edge for edge in path if edge.type == 'splice'])
            #         expr_obj_splice_count += -1 * c[idx] * spliceCount
            #     model.setObjectiveN(expr_obj_splice_count, 1, 0)

            model.setParam(gp.GRB.Param.OutputFlag, self.verbose)
            model.setParam(gp.GRB.Param.ConcurrentMIP, 1)
            model.setParam(gp.GRB.Param.Method, 4)

            model.setParam(gp.GRB.Param.IntFeasTol, 1e-9)
            model.setParam(gp.GRB.Param.FeasibilityTol, 1e-9)
            model.setParam(gp.GRB.Param.OptimalityTol, 1e-9)
            
            if self.verbose:
                print('-'*50)
                print(f"greedy solver iteration {kdx + 1}")
                print(f"number of splice edges is {self.graph.numSpliceEdges}")
                print(f"number of phasing reads is {self.graph.numPhasing}")
                print(f"length of shorted possible transcript is {lmin}")
                print('-'*50)

            model.optimize()
            
            if model.status == gp.GRB.Status.OPTIMAL or model.status == gp.GRB.Status.TIME_LIMIT:
    
                solC = model.getAttr('x', c)
                solG = model.getAttr('x', g)
                solW = model.getAttr('x', w)

                if self.verbose:
                    print("solution G:")
                    print(f"{solG}")
                    print("solution C:")
                    print(f"{solC}")

                if flag == 1:
                    kdx += 1
                    iter_count += 1
                    self.choose_best_k_paths(solC, kdx)
                    if kdx == self.numPaths:
                        flag = 3
                        npaths = 0
                    else:
                        flag = 2
                        npaths = 1
                elif flag == 2:
                    npaths = 0
                    npaths_added = self.add_path_to_chosen_with_expansion(solW, spliceEdges)
                    if npaths_added == 0:
                        if self.verbose:
                            print('+'*50)
                            print(f"no paths were added")
                            print('+'*50)
                        flag = 3
                    else:
                        #self.expand_all_paths()
                        flag = 1
                elif flag == 3:
                    self.decodeSolutionGreedyLikelihood(solC)
                    break
                    flag = 4
                else:
                    raise Exception(f"unknown value of flag {flag}")

            else:
                
                return False
                
        return True    
    
    
    
    def choose_best_k_paths(self, solC, kdx):
        
        if self.verbose:
            print("*"*50)
            print(f"iteration {kdx}")
            print(f"{solC}")
            print(f"{[self.graph.getPathIndex(path) for path in self.chosen_paths]}")
            print("*"*50)
        
        assert(kdx <= len(solC))
        # chosen_indices = [i for i,v in sorted(solC.items(), key=lambda x: (-x[1], len([edge for edge in self.chosen_paths[x[0]] if edge.type == 'splice'])))][:kdx]
        chosen_indices = [i for i,v in sorted(solC.items(), key=lambda x: (-x[1] * self.graph.getPathLength(self.chosen_paths[x[0]]), self.graph.getPathIndex(self.chosen_paths[x[0]])))][:kdx]

        self.chosen_paths = [self.chosen_paths[i] for i in chosen_indices]
       
        if self.verbose:
            print("*"*50)
            print(f"number of chosen paths = {len(self.chosen_paths)}")
            print(f"{chosen_indices}")
            print(f"{[self.graph.getPathIndex(path) for path in self.chosen_paths]}")
            print("*"*50)
        
    def add_path_to_chosen_with_expansion(self, solW, spliceEdges):
        
        new_path_spliceEdges = []
        nspliceEdges = 0
        for idx, edge in enumerate(spliceEdges):
            if solW[idx, 0] > 0.5:
                new_path_spliceEdges.append(edge)
                nspliceEdges += 1

        # print('+'*50)
        # print(f"adding {2**nspliceEdges} new paths")
        # print('+'*50)
        
        npaths_added = 0
        for case in range(2**nspliceEdges):
            expansionIndex = bin(case).replace("0b","").zfill(nspliceEdges)
            
            curr_path_spliceEdges = []
            for idx, edge in enumerate(new_path_spliceEdges):
                if expansionIndex[idx] == '1':
                    curr_path_spliceEdges.append(edge)
        
            curr_path = self.graph.stitchPath(curr_path_spliceEdges)
            if curr_path not in self.chosen_paths:
                self.chosen_paths.append(curr_path)
                npaths_added += 1
        
        return npaths_added
    
    def expand_all_paths(self):
        
        original_paths = self.chosen_paths[:]
        npaths_added = 0
        for path in original_paths:
            path_spliceEdges = [edge for edge in path if edge.type == 'splice']
            npath_spliceEdges = len(path_spliceEdges)
            
            for case in range(2**npath_spliceEdges - 1):
                expansionIndex = bin(case).replace("0b", "").zfill(npath_spliceEdges)

                case_path_spliceEdges = []
                for idx, edge in enumerate(path_spliceEdges):
                    if expansionIndex[idx] == '1':
                        case_path_spliceEdges.append(edge)
                
                case_path = self.graph.stitchPath(case_path_spliceEdges)
                if case_path not in self.chosen_paths:
                    self.chosen_paths.append(case_path)
                    npaths_added += 1
              
    def add_path_to_chosen(self, solW, spliceEdges = None):
        
        if not spliceEdges:
            spliceEdges = self.graph.transcriptSpliceEdges + self.graph.nonTranscriptSpliceEdges

        # get the new path splice edges
        new_path_spliceEdges = []
        for idx, edge in enumerate(spliceEdges):
            if solW[idx, 0] > 0.5:
                new_path_spliceEdges.append(edge)

        # piece together path
        self.chosen_paths.append(self.graph.stitchPath(new_path_spliceEdges))

    def decodeSolutionGreedyLikelihood(self, solCWeight):
    
        total_abundance = sum(solCWeight.values())
        
        solPaths = [(self.graph.getPathIndex(path), path, solCWeight[pathIndex] / total_abundance)
                    for pathIndex, path in enumerate(self.chosen_paths) if solCWeight[pathIndex] > 0]

        self.solution = sorted(solPaths, key = lambda x: x[0])
        
        self.supported_phasing_read_indices = {}
        for pathIndex, path in enumerate(self.chosen_paths):

            if solCWeight[pathIndex] == 0:
                continue

            path_id = self.graph.getPathIndex(path)
            path_spliceEdges = [edge for edge in path if edge.type == 'splice']
            
            supported_reads = []
            # find which phasing reads are supported
            for phase_id, phase in enumerate(self.graph.phasing):
                splus = phase[0]
                sminus = phase[1]
                if (set(splus).issubset(set(path_spliceEdges)) and 
                    len(set(sminus).intersection(set(path_spliceEdges))) == 0):                        
                    supported_reads.append((phase_id, phase[2]))
            
            self.supported_phasing_read_indices[f"path_{path_id}"] = supported_reads.copy()
        
    def solve_likelihood(self):
                
        model = gp.Model("solveTranscription")
        
        model.setParam(gp.GRB.Param.Threads, self.threads)
        if self.timelimit:
            model.setParam(gp.GRB.Param.TimeLimit, self.timelimit)

        if not self.graph.nonTranscriptSpliceEdges or not self.graph.transcriptSpliceEdges:
            self.graph.initCanonicalTranscripts()

        spliceEdges = self.graph.transcriptSpliceEdges + self.graph.nonTranscriptSpliceEdges
        numTranscriptSpliceEdges = len(self.graph.transcriptSpliceEdges)
        
        # splice edge to index
        spliceEdge2Index = {}
        for idx, edge in enumerate(spliceEdges):
            spliceEdge2Index[edge] = idx
        
        # splice lengths
        spliceLengths = []
        for edge in spliceEdges:
            spliceLengths.append(edge.right.left - edge.left.right)

        # function values and sampling points
        F = [0]*self.samplingFrequency 
        xh = [0]*self.samplingFrequency
        delta_x = 1.0/2**(self.samplingFrequency- 1)
        for i in range(self.samplingFrequency):
            if i == 0:
                xh[i] = 0
            else:
                xh[i] = 2**(i - 1) * delta_x
            if i == 0:
                F[i] = math.log(delta_x / 100)
            else:
                F[i] = math.log(xh[i])

        obj_scaling_factor = 10**4 / max([count for _, _, count in self.graph.phasing]) / abs(F[0])

        # finding lower limit of shortest possible transcript
        lmin = solveTranscription.find_shortest_transcript_length(spliceEdges, spliceLengths)
        lend = self.ref_length
        
        g = model.addVars(self.numPaths, lb = 0, ub = 1, vtype = gp.GRB.CONTINUOUS, name="g")
        w = model.addVars(self.graph.numSpliceEdges, self.numPaths, lb = 0, 
                          vtype = gp.GRB.BINARY, name = 'w')
        y = model.addVars(self.graph.numPhasing, self.numPaths, lb = 0, ub = 1,
                          vtype = gp.GRB.CONTINUOUS, name = 'y')
        ysplice = model.addVars(self.graph.numSpliceEdges, self.numPaths, lb = 0, ub = 1,
                                vtype = gp.GRB.CONTINUOUS, name = 'ysplice')
        ld = model.addVars(self.graph.numPhasing, self.samplingFrequency, lb = 0, ub = 1, 
                           vtype = gp.GRB.CONTINUOUS, name = 'lambda')

        observed_splice_edges = set()
        for phase in self.graph.phasing:
            observed_splice_edges |= set(phase[0])

        # remove unseen nontranscript splice edges
        for edge in spliceEdges:
            if edge not in observed_splice_edges:
                for pathIndex in range(self.numPaths):
                    model.addConstr(w[spliceEdge2Index[edge], pathIndex] == 0)

        # constraints
        # non-overlapping splice edges
        for idx1, edge1 in enumerate(spliceEdges):
            for idx2, edge2 in enumerate(spliceEdges):
                u1 = edge1.left.right
                v1 = edge1.right.left
                u2 = edge2.left.right
                v2 = edge2.right.left
                #if (u1 <= u2 and v1 >= u2) or (u2 <= u1 and v2 >= u1):
                if (u1 <= u2 and v1 > u2) and (idx1 != idx2):
                    for pathIndex in range(self.numPaths):
                        model.addConstr(w[idx1, pathIndex] + w[idx2, pathIndex] <= 1)

        # transcript contribution to phasing reads
        for idx in range(self.graph.numPhasing):
            for pathIndex in range(self.numPaths):
                model.addConstr(y[idx, pathIndex] <= g[pathIndex])
                
                splus = self.graph.phasing[idx][0]
                sminus = self.graph.phasing[idx][1]
                expr_plus = gp.LinExpr()
                expr_minus = gp.LinExpr()

                for edge in splus:
                    # spliceIndex = spliceEdges.index(self.graph.edges[edgeIndex])
                    spliceIndex = spliceEdge2Index[edge]
                    model.addConstr(y[idx, pathIndex] <= w[spliceIndex, pathIndex])
                    expr_plus += w[spliceIndex, pathIndex]

                for edge in sminus:
                    # spliceIndex = spliceEdges.index(self.graph.edges[edgeIndex])
                    spliceIndex = spliceEdge2Index[edge]
                    model.addConstr(y[idx, pathIndex] <= 1 - w[spliceIndex, pathIndex])
                    expr_minus += 1 - w[spliceIndex, pathIndex]
                
                nplus = len(splus)
                nminus = len(sminus)
                model.addConstr(y[idx, pathIndex] >= expr_plus + expr_minus + 
                                g[pathIndex] - nplus - nminus)

        # transcript contribution to splice edges
        for idx in range(self.graph.numSpliceEdges):
            for pathIndex in range(self.numPaths):
                model.addConstr(ysplice[idx, pathIndex] <= g[pathIndex])                    
                model.addConstr(ysplice[idx, pathIndex] <= w[idx, pathIndex])
                model.addConstr(ysplice[idx, pathIndex] >= g[pathIndex] + 
                                w[idx, pathIndex] - 1)

        # normalization
        ## noncanonical
        expr_norm = gp.LinExpr()        
        for pathIndex in range(self.numPaths):
            expr_norm += self.ref_length * g[pathIndex]
            for idx in range(self.graph.numSpliceEdges):
                expr_norm -= (spliceLengths[idx] * ysplice[idx, pathIndex])
        
        model.addConstr(expr_norm == lmin)

        # lambda sum to 1
        for i in range(self.graph.numPhasing):
            expr_sum = gp.LinExpr()
            for j in range(self.samplingFrequency):
                expr_sum += ld[i, j]
            model.addConstr(expr_sum == 1)
    
        # values at which we want to evalute the objective function
        for idx in range(self.graph.numPhasing):
            expr_x = gp.LinExpr()
            for j in range(self.samplingFrequency):
                expr_x += xh[j] * ld[idx, j]
            
            expr_y = gp.LinExpr()
            for pathIndex in range(self.numPaths):
                expr_y += y[idx, pathIndex]
            
            model.addConstr(expr_x == expr_y)

        # objective
        expr_obj = gp.LinExpr()
        total_reads = 0
        for idx in range(self.graph.numPhasing):
            nreads = self.graph.phasing[idx][2]
            total_reads += nreads
            for i in range(self.samplingFrequency):
                expr_obj += nreads * F[i] * ld[idx, i] * obj_scaling_factor

        model.setObjective(expr_obj, gp.GRB.MAXIMIZE)

        model.setParam(gp.GRB.Param.OutputFlag, self.verbose)

        if self.verbose:
            print('-'*50)
            print('decomposition solver')
            print(f"number of splice edges is {self.graph.numSpliceEdges}")
            print(f"number of transcript splice edges is {numTranscriptSpliceEdges}")
            print(f"number of phasing reads is {self.graph.numPhasing}")
            print('-'*50)

        model.optimize()
                
        if model.status == gp.GRB.Status.OPTIMAL:

            solG = model.getAttr('x', g)
            solW = model.getAttr('x', w)

            print(f"{solG}")
            
            self.decodeContinuousSolutionLikelihood(solW, solG, spliceEdges)
            
            return True
        else:
            return False

    def decodeContinuousSolutionLikelihood(self, solW, solWeight, spliceEdges = None):
        
        if not spliceEdges:
            spliceEdges = self.graph.transcriptSpliceEdges + self.graph.nonTranscriptSpliceEdges
        
        total_abundance = sum(solWeight.values())
        
        self.transcriptWeights = []
        self.supported_phasing_read_indices = {}
                    
        solPaths = []
        for pathIndex in range(self.numPaths):

            if solWeight[pathIndex] == 0:
                continue

            path = []
            # get the path splice edges
            path_spliceEdges = []
            for idx, edge in enumerate(spliceEdges):
                if solW[idx, pathIndex] > 0.5:
                    path_spliceEdges.append(edge)
                        
            # piece together path
            path = self.graph.stitchPath(path_spliceEdges)            
            path_id = self.graph.getPathIndex(path)
            
            supported_reads = []
            # find which phasing reads are supported
            for phase_id, phase in enumerate(self.graph.phasing):
                splus = phase[0]
                sminus = phase[1]
                if (set(splus).issubset(set(path_spliceEdges)) and 
                    len(set(sminus).intersection(set(path_spliceEdges))) == 0):                        
                    supported_reads.append((phase_id, phase[2]))
            
            self.supported_phasing_read_indices[f"path_{path_id}"] = supported_reads.copy()
            
            solPaths.append((path_id, path, solWeight[pathIndex] / total_abundance))

        self.solution = sorted(solPaths, key = lambda x: x[0])


    def solve(self, greedy = False):

        if not greedy:
            if self.solve_likelihood():
                return True
            else:
                return False
        else:
            if self.solve_greedy_likelihood():
                return True
            else:
                return False

    def filterSolution(self, threshold, read_length = 100):
        
        total_reads = sum([x[2] for x in self.graph.phasing])
        total_length = sum([weight * self.graph.getPathLength(path) 
                            for _, path, weight in self.solution])

        prefilter_number = len(self.solution)
        self.solution = [(index, path, weight)
                         for index, path , weight in self.solution 
                         if weight * total_reads * self.graph.getPathLength(path) / total_length >= threshold]

        # renormalize the weights
        total_weight = sum([weight for _, _, weight in self.solution])
        self.solution = [(index, path, weight/total_weight)
                         for index, path , weight in self.solution] 

        if self.verbose:
            print(f"filtered out {prefilter_number - len(self.solution)} transcripts")
    
    def writeSolution(self, filename):
        
        with open(filename, 'w') as output:
            if self.transcriptWeights:
                for idx, weight in self.transcriptWeights:
                    output.write(f">transcript_{idx}\t{self.graph.transcript_type[idx]}:{weight}\n")
                    for edge in self.graph.transcripts[idx]:
                        output.write(f"{edge.left.name}\t{edge.right.name}\t{edge.type}\n")

            if self.solution:
                for pathIndex, path, weight in self.solution:
                    output.write(f">Path_{pathIndex}:{weight}\n")
                    for edge in path:
                        output.write(f"{edge.left.name}\t{edge.right.name}\t{edge.type}\n")
        
    def writeGTF(self, filename, contig='NC_045512.2'):
        
        with open(filename, 'w') as output:
            if self.transcriptWeights:
                for idx, weight in self.transcriptWeights:
                    path = self.graph.transcripts[idx]
                    spliceEdges = [edge for edge in path if edge.type == 'splice']
                    assert(len(spliceEdges) <= 1)
                    if len(spliceEdges) == 0:
                        output.write(f"{contig}\tjumper\ttranscript\t{path[0].left.left}\t{path[-1].right.right}\t1000\t.\t.\tRPKM \"{weight}\"\n")
                    else:
                        output.write(f"{contig}\tjumper\ttranscript\t{path[0].left.left}\t{path[-1].right.right}\t1000\t.\t.\tRPKM \"{weight}\"\n")
                        output.write(f"{contig}\tjumper\texon\t{path[0].left.left}\t{spliceEdges[0].left.right}\t1000\t.\t.\t.\n")
                        output.write(f"{contig}\tjumper\texon\t{spliceEdges[0].right.left}\t{path[-1].right.right}\t1000\t.\t.\t.\n")
            
            if self.solution:
                for pathIndex, path, weight in self.solution:
                    output.write(f"{contig}\tjumper\ttranscript\t{path[0].left.left}\t{path[-1].right.right}\t1000\t.\t.\tRPKM \"{weight}\"\n")
                    spliceEdges = [edge for edge in path if edge.type == 'splice']
                    if len(spliceEdges) == 0:
                        output.write(f"{contig}\tjumper\ttranscript\t{path[0].left.left}\t{path[-1].right.right}\t1000\t.\t.\tRPKM \"{weight}\"\n")
                    else:
                        for idx, edge in enumerate(spliceEdges):
                            if idx == 0:
                                output.write(f"{contig}\tjumper\texon\t{path[0].left.left}\t{edge.left.right}\t1000\t.\t.\t.\n")
                            else:
                                output.write(f"{contig}\tjumper\texon\t{spliceEdges[idx-1].right.left}\t{edge.left.right}\t1000\t.\t.\t.\n")
                        output.write(f"{contig}\tjumper\texon\t{edge.right.left}\t{path[-1].right.right}\t1000\t.\t.\t.\n")
        
    def writeMatching(self, filename):
        
        with open(filename, 'w') as output:
            for key, val in self.supported_phasing_read_indices.items():
                output.write(f">{key}\n")
                for idx, count in val:
                    output.write(f"{idx}\t{count}\n")
