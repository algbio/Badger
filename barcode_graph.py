#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2024 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import editdistance
from collections import defaultdict
import pandas as pd
import logging
import itertools
from concurrent.futures import ProcessPoolExecutor
from statistics import median, mean

from index import QGramIndex
from common import dfs_without_recursion, rank, unrank

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
        self.counts = defaultdict(int) #list of my nodes with their counts #changing this to ranked instead of sequence as the key
        self.edges = defaultdict(list) 
        self.dists = defaultdict(int)
        self.clusters = defaultdict(list)
        self.clustering = dict()
        self.clustered = defaultdict(bool)
        self.index = None
    
    def index_chunk(self, barcode_chunk, bc_len, num):
        
        index = QGramIndex(self.threshold, bc_len, 6)
        counts = dict() 
        
        for sequence in barcode_chunk:
            if len(sequence) == bc_len + 1:
                sequence = sequence[:-1]
            if len(sequence) == bc_len:
                bc_rank = rank(sequence, bc_len)
                if bc_rank in counts.keys():
                    counts[bc_rank] += 1
                else:
                    counts[bc_rank] = 1
                    index.add_to_index(sequence, bc_rank)
        
        chunk = Chunk()
        chunk.index = index
        chunk.counts = counts
        
        if num%10 == 0:
            logger.info("Finished indexing chunk" + str(num))
        return chunk
        
    def compare_chunk(self, bc_chunk, bc_len, num): 
        
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
        
        if num%10 == 0:
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
                    num = rank(sequence, bc_len)
                    self.counts[num] = 1
                    self.index.add_to_index(sequence, num)
                    
    
    def graph_construction(self, barcodes, bc_len, threads): 
        self.index = QGramIndex(self.threshold, bc_len, 6)
        num = 0
        
        if threads > 1:
            
            self.index_bc_in_parallel(barcodes, bc_len, threads)
            
            self.compare_in_parallel(bc_len, threads)

        else: 
            
            self.index_bc_single_thread(barcodes, bc_len)
            
            for bc_rank in self.counts.keys():
                # info: keys of counts are still the sequences but everywhere else I use the rank
                # that means that I need to rank once in the beginning
                # and then for everything that is in the closest set I need to unrank
                
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
                                self.edges[bc_rank].append(seq_rank)
                                self.edges[seq_rank].append(bc_rank)
                                self.dists[(bc_rank, seq_rank)] = dist
                                self.dists[(seq_rank, bc_rank)] = dist 
                            
        
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
                if unrank(bc_by_counts[i], bc_len) in barcode_list:
                    tbcs.append(bc_by_counts[i])
                    n+=1
                i+=1
        else:
            while self.counts[bc_by_counts[i]] > cutoff and n <= n_cells + n_cells*interval*0.01:
                tbcs.append(bc_by_counts[i])
                i += 1
                n += 1
        while n < n_cells - n_cells*interval*0.01:
            tbcs.append(bc_by_counts[i])
            i += 1
            n += 1
        return tbcs
    
    def cluster(self, true_barcodes, barcode_list, n_cells, bc_len, interval):
        
        tbcs = self.get_cluster_centers(true_barcodes, bc_len, barcode_list, n_cells, interval)
        
        for tbc in tbcs:
            self.clusters[tbc] = [tbc]
            self.clustering[tbc] = (tbc, 0)
            self.clustered[tbc] = True
            
        for i in range(1,3):
            #print(i)
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
            if not visited[node]:
                component, visited = dfs_without_recursion(visited, node, [], self.edges)
                components.append(component)
                
        #print(len(read_assignment))
        
        assignments = self.assign_by_cluster()
        read_ids = []
        results = []
        real = []
        
        for read in read_assignment:
            observed_bc = read[2]
            if observed_bc != "*":  
                assigned_bc = assignments[observed_bc]
                real_bc = read[1]
                if assigned_bc != "":
                    read_ids.append(read[0])
                    real.append(real_bc)
                    results.append(assigned_bc)
                    if real_bc != assigned_bc:
                        print("test")
            
        #print(len(read_ids))
        #print(len(real))
        #print(len(results))
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

    def output_file(self, read_assignments, out_file, barcode_detector, post):
        bc_len = barcode_detector.barcode_length
        assignments = self.assign_by_cluster(bc_len)
        if post:
            assignments = self.postprocessing(assignments, bc_len)

        with open(out_file, "w") as outf:
            outf.write(barcode_detector.header() + "\n")
            for read in read_assignments:
                observed_bc = read.detected_results["Barcode"].seq
                assigned_bc = "*"
                if observed_bc != "*":
                    assigned_bc = assignments[observed_bc]
                    if assigned_bc == "":
                        assigned_bc = "*"
                read.detected_results["Barcode"].seq = assigned_bc
                outf.write(barcode_detector.format_result(read) + "\n")


