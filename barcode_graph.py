#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2024 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import math
import numpy as np
import editdistance
import networkx as nx
import matplotlib.pyplot as plt
import igraph as ig
from collections import defaultdict
import pandas as pd
import logging
import edlib
import itertools
from Levenshtein import distance
from concurrent.futures import ProcessPoolExecutor
from statistics import median, mean

from index import QGramIndex
from common import get_score, dfs, dfs_without_recursion, rank, unrank

logger = logging.getLogger("BarcodeGraph")

READ_CHUNK_SIZE = 100000
BC_CHUNK_SIZE = 10000

class Chunk:
    
    index = QGramIndex(1, 16, 6)
    counts = dict()
    
class Edgedist:
    
    edges = dict()
    dists = dict()

class BarcodeGraph:
    
    def __init__(self, threshold):
        
        self.threshold = threshold
        self.counts = defaultdict(int) #list of my nodes with their counts #changing this to ranked instead of sequence as the key, kills all stats functions probably, need to look at them before I use them again
        self.edges = defaultdict(list) # idea: change list to set
        self.dists = defaultdict(int)
        #self.barcodes = defaultdict(str)
        #self.numbered = defaultdict(int)
        self.clusters = defaultdict(list)
        self.clustering = dict()
        self.clustered = defaultdict(bool)
        self.index = None
    
    #not sure I love this here, think about moving that to a different place
    def index_chunk(self, barcode_chunk, bc_len, num):
        
        #print("Processing chunk" + str(num))
        index = QGramIndex(self.threshold, bc_len, 6)
        counts = dict() #wäre hier defaultdict sinnvoller?
        # kann ich nicht theoretisch direkt mit self.counts arbeiten? Oder funktioniert nur daraus lesen aber nicht es verändern
        # könnte mir vorstellen dass wenn ich es verändere nicht alle veränderungen übernommen werden
        # aber lesen müsste eigentlich gehen und dann hat eh alles alle self.counts
        
        for sequence in barcode_chunk:
            if len(sequence) == bc_len + 1:
                sequence = sequence[:-1]
            if len(sequence) == bc_len:
                bc_rank = rank(sequence, bc_len)
                if bc_rank in counts.keys():
                    counts[bc_rank] += 1
                else:
                    # bc_rank = rank(sequence, bc_len) #TODO add unrank function
                    counts[bc_rank] = 1
                    index.add_to_index(sequence, bc_rank)
        
        chunk = Chunk()
        chunk.index = index
        chunk.counts = counts
        
        logger.info("Finished indexing chunk" + str(num))
        return chunk
        
    def compare_chunk(self, bc_chunk, bc_len, num): # qär vielleicht sinnvoll zu self.index zu machen damit ich es nicht jedes mal weitergeben muss
        # zu bc_chunk: gibt es irgendeine Option, dass ich das nicht machen muss? Die sind ja alle in self.counts, da kann ich drauf zugreifen
        
        edges = defaultdict(list)
        dists = defaultdict(int)
        
        for bc_rank in bc_chunk:
            
            barcode = unrank(bc_rank, bc_len)
        
            closest = self.index.get_close(barcode, bc_rank)
                
            for seq_rank in closest: 
                
                if seq_rank > bc_rank:
                    
                    sequence = unrank(seq_rank, bc_len)
                    
                    if sequence == "" or sequence == barcode: 
                        continue
                    else: 
                        
                        dist = min(editdistance.eval(barcode, sequence), editdistance.eval(barcode[:-1],sequence), editdistance.eval(barcode,sequence[:-1]))
                        
                        if dist <= self.threshold:
                            edges[bc_rank].append(seq_rank)
                            edges[seq_rank].append(bc_rank)
                            dists[(bc_rank, seq_rank)] = dist
                            dists[(seq_rank, bc_rank)] = dist 
                        
        output = Edgedist()
        output.edges = edges
        output.dists = dists
        
        logger.info("Finished comparing chunk" + str(num))
        
        return output 
        
    def get_read_chunks(self, barcodes):
        current_chunk = []
        for r in barcodes:
            current_chunk.append(r)
            if len(current_chunk) >= READ_CHUNK_SIZE:
                yield current_chunk
                current_chunk = []
        yield current_chunk
        
    def get_barcode_chunks(self, barcodes):
        current_chunk = []
        for bc in barcodes:
            current_chunk.append(bc)
            if len(current_chunk) >= BC_CHUNK_SIZE:
                yield current_chunk
                current_chunk = []
        yield current_chunk
        
    def index_bc_in_parallel(self, barcodes, bc_len, threads):
        
        #self.index = QGramIndex(self.threshold, bc_len, 6)
        barcode_chunks = self.get_read_chunks(barcodes)
        
        indexing_gen = (
            self.index_chunk,
            barcode_chunks,
            itertools.repeat(bc_len),
            itertools.count(start=0, step=1),
        )
        
        with ProcessPoolExecutor(max_workers = threads) as proc:
            chunks = proc.map(*indexing_gen, chunksize = 1)
            
        merged_index = self.index.index
        
        for chunk in chunks:
            ix = chunk.index
            i = ix.index
            c = chunk.counts
            
            for j in range(len(i)):
                kmer_dict = i[j]
                for key in kmer_dict.keys():
                    if key not in merged_index[j].keys():
                        merged_index[j][key] = kmer_dict[key]
            
            for key in c.keys():
                self.counts[key] += c[key]
        
        self.index.index = merged_index
        
        #return index
    
    def compare_in_parallel(self, bc_len, threads):
        
        barcodes = sorted(self.counts.keys())
        
        bc_chunks = self.get_barcode_chunks(barcodes)
        
        comparing_gen = (
            self.compare_chunk,
            bc_chunks,
            itertools.repeat(bc_len),
            itertools.count(start=0, step=1),
        ) 
        
        with ProcessPoolExecutor(max_workers = threads) as proc:
            chunks = proc.map(*comparing_gen, chunksize = 1)
            
        for chunk in chunks:
            
            edges = chunk.edges
            dists = chunk.dists
            
            for key in edges.keys():
                self.edges[key].extend(edges[key])
            
            for key in dists.keys():
                self.dists[key] = dists[key]
                
        
    def index_bc_single_thread(self, barcodes, bc_len):
        
        for sequence in barcodes:
            if len(sequence) == bc_len + 1:
                sequence = sequence[:-1] 
            if len(sequence) == bc_len:
                num = rank(sequence, bc_len)
                if num in self.counts.keys():
                    self.counts[num] +=1
                else:
                    #self.barcodes[num] = sequence
                    #self.numbered[sequence] = num
                    num = rank(sequence, bc_len)
                    self.counts[num] = 1
                    self.index.add_to_index(sequence, num)
                    
        #return index
    
    def graph_construction(self, barcodes, bc_len, threads): #could bc_len be a self.len or something like that?
        self.index = QGramIndex(self.threshold, bc_len, 6)
        num = 0
        
        if threads > 1:
            
            self.index_bc_in_parallel(barcodes, bc_len, threads)
            
            #bc_chunks # don't knoe how I get them yet
            self.compare_in_parallel(bc_len, threads)
            
            
            
  
        else: 
            logger.info("Using edlib")
            
            self.index_bc_single_thread(barcodes, bc_len)
            
            for bc_rank in self.counts.keys():
                # info: keys of counts are still the sequences but everywhere else I use the rank
                # that means that I need to rank once in the beginning
                # and then for everything that is in the closest set I need to unrank
                # is this really more efficient than keeping them in a dictionary and looking them up?
                # for multithreaded it probably is, and much more space efficient, but I am not sure about single threaded
                
                barcode = unrank(bc_rank, bc_len)
                
                # this is comparing to everything in the index so I do some comparisons multiple times! I might need to rewrite the "closest" function to only compare to lexicographically higher
                # however that is only partly what I want for the parallel right? no I think maybe that still makes sense
                # If I rewrite already closest might also make 
                # I might actually have to use the other ranking function because then it should be lexicographically sorted in counts
                # and then I just do only compare to things that are higher in the counts 
                # but does that work with closest or do I then need to do it entirely different
                # because in index it is first sorted by k-mer rank and then a dict
                # and dicts are not sorted
                closest = self.index.get_close(barcode, bc_rank)
                
                for seq_rank in closest:
                    
                    if seq_rank > bc_rank:
                    
                        sequence = unrank(seq_rank, bc_len)
                        
                        if sequence == "" or sequence == barcode: 
                            continue
                        else: 
                            
                            dist = min(editdistance.eval(barcode, sequence), editdistance.eval(barcode[:-1],sequence), editdistance.eval(barcode,sequence[:-1]))
                            
                            if dist <= self.threshold:
                                self.edges[bc_rank].append(seq_rank)
                                self.edges[seq_rank].append(bc_rank)
                                self.dists[(bc_rank, seq_rank)] = dist
                                self.dists[(seq_rank, bc_rank)] = dist 
                        #for s in self.barcodes.keys():
                        #    seq = self.barcodes[s]
                        # closest = index.get_close(sequence, num)
                        # for s in closest:
                            # seq = self.barcodes[s]
                            # if seq == "" or seq == sequence:
                                # continue
                            # else:
                                # # test: min of three distances to possibly ignore last position for indels
                                # #dist = editdistance.eval(sequence, seq)
                                # #dist = min(editdistance.eval(sequence, seq), editdistance.eval(sequence[:-1],seq), editdistance.eval(sequence,seq[:-1]))
                                # #r1 = edlib.align(sequence, seq, mode = "SHW")
                                # #r2 = edlib.align(seq, sequence, mode = "SHW")
                                # #r = edlib.align(sequence, seq, mode = "NW", task = "path")
                                # #dist = r["editDistance"]
                                # #path = r["cigar"]
                                # #if path[-1] == 'I' or path[-1] == 'D':
                                # #    dist = dist -1
                                # #dist = min(r1["editDistance"], r2["editDistance"])
                                # dist = min(distance(sequence, seq, score_cutoff = 1), distance(sequence[:-1],seq, score_cutoff = 1), distance(sequence,seq[:-1], score_cutoff = 1))
                                # #score = get_score(sequence, seq)
                                # if dist <= self.threshold:
                                # #if score >= 16*3 - 4: ## fully equal is 3*len, 2 indels would be -4
                                    # self.edges[num].append(s)
                                    # self.edges[s].append(num)
                                    # # #self.dists[(sequence,seq)] = dist # think if I should add both sides or one is enough
                                    # self.dists[(num,s)] = dist 
                                    # self.dists[(s,num)] = dist
                        # num += 1
                        # if num%500000 == 0:
                            # logger.info(f"processed {num} distinct barcodes")
                            
        
    def get_cluster_centers(self, true_barcodes, bc_len, barcode_list, n_cells, interval):
        sorted_counts = dict(sorted(self.counts.items(), key=lambda item: item[1],reverse = True))
        bc_by_counts = list(sorted_counts.keys())
        cutoff = mean(list(self.counts.values())[:(n_cells)])
        cutoff = max(cutoff/5.0, 5)
        tbcs = []
        n = 0
        i = 0
        if true_barcodes:
            tbcs = [rank(bc, bc_len) for bc in true_barcodes]
        elif barcode_list:
            while i < len(bc_by_counts) and self.counts[bc_by_counts[i]] > cutoff and n <= n_cells + n_cells*interval*0.01:
                #print(bc_by_counts[i])
                if unrank(bc_by_counts[i], bc_len) in barcode_list:
                    tbcs.append(bc_by_counts[i])
                    n+=1
                i+=1
        else:
            while self.counts[bc_by_counts[i]] > cutoff and n <= n_cells + n_cells*interval*0.01:
                tbcs.append(bc_by_counts[i])
                i += 1
                n += 1
            #tbcs = bc_by_counts[:n_cells]
        while n < n_cells - n_cells*interval*0.01:
            tbcs.append(bc_by_counts[i])
            i += 1
            n += 1
        return tbcs
    
    def cluster(self, true_barcodes, barcode_list, n_cells, bc_len, interval):
        #self.clustered = [False for node in self.counts.keys()]
        
        tbcs = self.get_cluster_centers(true_barcodes, bc_len, barcode_list, n_cells, interval)
        
        for tbc in tbcs:
            # tbc = rank(tbc, bc_len)
            self.clusters[tbc] = [tbc]
            self.clustering[tbc] = (tbc, 0)
            self.clustered[tbc] = True
            
        for i in range(1,3):
            print(i)
            for center in self.clusters.keys():
                for n in range(len(self.clusters[center])):
                    node = self.clusters[center][n]
                    for neighbor in self.edges[node]:
                        if not self.clustered[neighbor]:
                            self.clusters[center].append(neighbor)
                            self.clustering[neighbor] = (center,i)
                            self.clustered[neighbor] = True 
                        elif self.clustering[neighbor][0] != center and self.clustering[neighbor][0] != -1:
                            if self.clustering[neighbor][1] == i:
                                self.clusters[self.clustering[neighbor][0]].remove(neighbor)
                                self.clustering[neighbor] = (-1, -1)
        
    def get_assignments(self, true_barcodes, components):
        observed_assignments = defaultdict(str)
        for component in components:
            correct_nodes = []
            for node in component:
                bc = self.barcodes[node]
                if bc in true_barcodes:
                    correct_nodes.append(bc)
            for node in component:
                bc = self.barcodes[node]
                min_dist = 32
                for tbc in correct_nodes:
                    dist = editdistance.eval(tbc, bc)
                    #dist = min(editdistance.eval(tbc, bc),editdistance.eval(tbc[:-1], bc),editdistance.eval(tbc, bc[:-1]))
                    if dist < min_dist:
                        min_dist = dist
                        observed_assignments[bc] = tbc
        return observed_assignments
        
    def assign_by_cluster(self, bc_len):
        observed_assignments = defaultdict(str)
        for node in self.counts.keys():
            if self.clustered[node] and self.clustering[node][0] != -1:
                bc = unrank(node, bc_len)
                tbc = unrank(self.clustering[node][0], bc_len)
                observed_assignments[bc] = tbc
        return observed_assignments
                
    def isoquant_output(self, read_assignment, true_barcodes):
        
        components = []
        visited = [False for node in self.barcodes.keys()]
        for node in self.barcodes.keys():
            #print(node)
            if not visited[node]:
                component, visited = dfs_without_recursion(visited, node, [], self.edges)
                components.append(component)
                
        print(len(read_assignment))
        
        assignments = self.assign_by_cluster()
        #assignments = self.get_assignments(true_barcodes, components)
        read_ids = []
        results = []
        real = []
        
        for read in read_assignment:
            
            #print(read)
            
            observed_bc = read[2]
            if observed_bc != "*":  ## don't know yet if I am keeping those in?
                assigned_bc = assignments[observed_bc]
                real_bc = read[1]
                if assigned_bc != "":
                    #print(read[0])
                    #print(real_bc)
                    #print(assigned_bc)
                    read_ids.append(read[0])
                    real.append(real_bc)
                    results.append(assigned_bc)
                    if real_bc != assigned_bc:
                        print("test")
            
        print(len(read_ids))
        print(len(real))
        print(len(results))
        correct = pd.DataFrame({"readID": read_ids,
                                "barcode": real})
        correct.to_csv('Human_R9_V5.3_tt_read_groups_tbcs.tsv', sep = '\t', index = False)
        res = pd.DataFrame({"readID": read_ids,
                            "barcode": results})
        res.to_csv('Human_R9_V5.3_tt_read_groups_tool.tsv', sep = '\t', index = False)
        
    def postprocessing(self, assignments, bc_len):
        bcs = self.counts.keys()
        cluster_centers = set(assignments.values())
        for bc_rank in bcs:
            bc = unrank(bc_rank, bc_len)
            if assignments[bc] == "" or assignments[bc] == "*":
                min_dist = 16
                min_bc = "*"
                for tbc in cluster_centers:
                    dist = editdistance.eval(bc, tbc)
                    if dist < min_dist:
                        min_dist = dist
                        min_bc = tbc
                if min_dist < 3:
                    assignments[bc] = min_bc
        return assignments
            
                
    def output_file(self, read_assignment, out, true_barcodes, bc_len, post):
        
        # components = []
        # visited = [False for node in self.barcodes.keys()]
        # for node in self.counts.keys():
            # if not visited[node]:
                # component, visited = dfs_without_recursion(visited, node, [], self.edges)
                # components.append(component)
        
        assignments = self.assign_by_cluster(bc_len)
        #assignments = self.get_assignments(true_barcodes, components)
        read_ids = []
        results = []
        if post:
            assignments = self.postprocessing(assignments, bc_len)
        for read in read_assignment:
            observed_bc = read[1]
            assigned_bc = "*"
            if observed_bc != "*":
                assigned_bc = assignments[observed_bc]
                if assigned_bc == "":
                    assigned_bc = "*"
                    
            read_ids.append(read[0])
            results.append(assigned_bc)
        
        out_file = out + "_output_file.tsv"

        res = pd.DataFrame({"readID": read_ids,
                            "barcode": results})
        res.to_csv(out_file, sep = '\t', index = False)

