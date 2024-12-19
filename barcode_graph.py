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
        self.assignment = None
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
                            
    
                
    def visualize_graph(self):
        #G = nx.Graph()
        #G.add_nodes_from(self.counts.keys())
        edges = []
        for edge in self.dists.keys():
            #G.add_edge(edge[0], edge[1], weight = self.dists[edge])
            edges.append([edge[0], edge[1]])
        #G = nx.complete_graph(5)
        #nx.draw(G)
        #plt.show()
        g = ig.Graph(n = len(self.counts.keys()), edges = edges)
        unconnected = []
        for i in range(0, len(self.counts.keys())):
            if g.degree(i) == 0:
                unconnected.append(i)
        g.delete_vertices(unconnected) 
        layout = g.layout("fr")
        ig.plot(g, target = 'graph_dist_2.pdf', vertex_size = 1, layout = layout, edge_color = ['red', 'black'])
        
    def cluster(self, true_barcodes, barcode_list, n_cells, bc_len):
        #self.clustered = [False for node in self.counts.keys()]
        
        sorted_counts = dict(sorted(self.counts.items(), key=lambda item: item[1],reverse = True))
        bc_by_counts = list(sorted_counts.keys())
        tbcs = []
        n = 0
        i = 0
        if true_barcodes:
            tbcs = [rank(bc, bc_len) for bc in true_barcodes]
        elif barcode_list:
            while n < n_cells and i < len(bc_by_counts):
                #print(bc_by_counts[i])
                if unrank(bc_by_counts[i], bc_len) in barcode_list:
                    tbcs.append(bc_by_counts[i])
                    n+=1
                i+=1
        else:
            tbcs = bc_by_counts[:n_cells]
        
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
                                
    #move everything statistics related to a different file, instead of self then just give it the graph
    #doesn't have to be in a class, more like common
    def graph_statistics(self, true_barcodes):
        components = []
        singletons = []
        lengths = []
        false_components = []
        visited = [False for node in self.barcodes.keys()]
        both = 0
        degree_better = 0
        count_better = 0
        min_dist = 32
        dists = []
        dists_d = []
        dists_c = []
        dists_n = []
        best_is_max = 0
        n = 0
        correct_in_component = []
        for node in self.barcodes.keys():
            if not visited[node]:
                component, visited = dfs_without_recursion(visited, node, [], self.edges)
                if len(component) == 1:
                    singletons.append(component[0])
                    #print(component)
                else:
                    lengths.append(len(component))
                    #print("Component length:", len(component))
                    correct = 0
                    correct_nodes = []
                    # if len(component) > 80:
                        # G = nx.Graph()
                        # G.add_nodes_from(component)
                        # for node in component:
                            # for neighbour in self.edges[node]:
                                # G.add_edge(node, neighbour, weight = self.dists[(node, neighbour)] * 10)
                        # figname = "lserved_assignments[bc] = tbcarge_component_" + str(n) + ".png"
                        # nx.draw(G, node_size = 50)
                        # plt.savefig(figname)
                        # n += 1
                    min_dist_n = 32
                    min_bc_n = -1
                    min_node = -1
                    max_degree = 0
                    max_degree_node = -1
                    max_count = 0
                    max_count_node = -1
                    for node in component:
                        node_dist = 32
                        if len(self.edges[node]) == 0:
                            print(component)
                            print("Singleton?")
                        if len(self.edges[node]) > max_degree:
                            max_degree = len(self.edges[node])
                            max_degree_node = node
                        if self.counts[self.barcodes[node]] > max_count:
                            max_count = self.counts[self.barcodes[node]]
                            max_count_node = node
                        for bc in true_barcodes:
                            dist = editdistance.eval(self.barcodes[node], bc)
                            if dist < min_dist_n:
                                min_dist_n = dist
                                min_node = node
                                min_bc_n = bc
                            if dist < node_dist:
                                node_dist = dist
                        if node_dist == 0:
                            correct_nodes.append(node)
                            correct +=1
                        dists_n.append(min_dist_n)
                    if len(self.edges[max_count_node]) == max_degree:
                        max_degree_node = max_count_node
                    # #best_is_max += (min_node == max_degree_node or min_node == max_count_node)
                    #print("New component")
                    #print("Number of nodes:", len(component))
                    #print("Node with highest degree:", max_degree_node, "degree:", max_degree)
                    #print("Count of the node with the highest degree:", self.counts[self.barcodes[max_degree_node]])
                    #print("Node with highest count:", max_count_node, "count:", max_count)
                    #print("Degree of the node with the highest count:", len(self.edges[max_count_node]))
                    min_bc_c = -1
                    min_bc_d = -1
                    min_dist_c = 32
                    min_dist_d = 32
                    for bc in true_barcodes:
                        dist_c = editdistance.eval(self.barcodes[max_count_node], bc)
                        dist_d = editdistance.eval(self.barcodes[max_degree_node], bc)
                        if dist_c < min_dist_c:
                            min_dist_c = dist_c
                            min_bc_c = bc
                        if dist_d < min_dist_d:
                            min_dist_d = dist_d
                            min_bc_d = bc
                    #print("Closest true barcode to max degree:", min_bc_d, "Distance:", min_dist_d)
                    #print("Closest true barcode to max count:", min_bc_c, "Distance:", min_dist_c)
                    #print("Closest true barcode to any node:", min_bc_n, "Distance:", min_dist_n)
                    both += (max_degree_node == max_count_node)
                    count_better += (min_dist_c < min_dist_d)
                    degree_better += (min_dist_d < min_dist_c)
                    if min_dist_c < min_dist: 
                        min_dist = min_dist_c
                    if min_dist_d < min_dist:
                        min_dist = min_dist_d
                    dists.append(min_dist_c)
                    dists.append(min_dist_d)
                    dists_c.append(min_dist_c)
                    dists_d.append(min_dist_d)
                    best_is_max += (min_dist_n == min_dist_d or min_dist_n == min_dist_c)
                    if min(min_dist_d, min_dist_c, min_dist_n) > 1:
                        #print(len(component))
                        false_components.append(len(component))
                    #print("Number of correct barcodes in the component:", correct)
                    if correct > 1:
                        bc_in_component = []
                        G = nx.Graph()
                        color_map = []
                        G.add_nodes_from(component)
                        for node in component:
                            for neighbour in self.edges[node]:
                                G.add_edge(node, neighbour, weight = self.dists[(node, neighbour)] * 10)
                        for node in G:
                            if self.barcodes[node] in true_barcodes:
                                bc_in_component.append(self.barcodes[node])
                                color_map.append('red')
                            else:
                                color_map.append('blue')
                        #nx.draw(G, node_size = 50, node_color = color_map)
                        #plt.show()
                        for bc in bc_in_component:
                            for bc2 in bc_in_component:
                                if bc!= bc2:
                                    dist = editdistance.eval(bc,bc2)
                                    if dist <= 1:
                                        print(bc, bc2, "distance:", dist)
                    correct_in_component.append(correct)
                components.append(component)
        single_counts = []
        for node in singletons:
            single_counts.append(self.counts[self.barcodes[node]])
        print("number of components:", len(components))
        print("number of singletons", len(singletons))
        print("maximal component size", max(lengths))
        print("Number of components with equal max degree and max count node:", both)
        print("Number of times max count has closer match than max degree:", count_better)
        print("Number of times max degree has closer match than max count:", degree_better)
        print("Minimum distance of any max node to a true barcode:", min_dist)
        print("Number of times the node with minimum distance to a true barcode is a max node:", best_is_max)
        # plt.hist(dists)
        # plt.title("Distances of all max nodes")
        # plt.show()
        # plt.hist(dists_c)
        # plt.title("Distances of maximum count nodes")
        # plt.show()
        # plt.hist(dists_d)
        # plt.title("Distances of maximum degree nodes")
        # plt.show()
        # plt.hist(dists_n)
        # plt.title("Minimum distance in each component")
        # plt.show()
        plt.hist(lengths, bins = 100)
        plt.title("Size of the components")
        plt.yscale("log")
        plt.show()
        plt.hist(lengths, bins = 500)
        plt.title("Size of the components")
        plt.yscale("log")
        plt.show()
        lengths.remove(max(lengths))
        plt.hist(lengths, bins = 100)
        plt.title("Size of the components without the largest component")
        plt.yscale("log")
        plt.show()
        plt.hist(correct_in_component, bins = 1000)
        plt.title("Number of correct barcodes in components")
        plt.yscale("log")
        plt.show()
        plt.hist(correct_in_component, bins = 100)
        plt.title("Number of correct barcodes in components")
        plt.yscale("log")
        plt.show()
        correct_in_component.remove(max(correct_in_component))
        plt.hist(correct_in_component, bins = 5)
        plt.title("Number of correct barcodes in components")
        plt.yscale("log")
        plt.show()
        plt.hist(false_components, bins = 100)
        plt.title("Size of components with no true barcodes")
        plt.yscale("log")
        plt.show()
        plt.hist(single_counts, bins = 20)
        plt.title("Counts of singleton nodes")
        plt.yscale("log")
        plt.show()
        # self.closest_true(singletons, true_barcodes)
        
    def closest_true(self, singletons, true_barcodes):
        true_barcodes = set(true_barcodes)
        closest = []
        dists = []
        seqs = []
        for node in singletons:
            min_bc = -1
            min_dist = 32
            for bc in true_barcodes:
                dist = editdistance.eval(self.barcodes[node], bc)
                if dist < min_dist:
                    min_bc = bc
                    min_dist = dist
            closest.append(min_bc)
            dists.append(min_dist)
            seqs.append(self.barcodes[node])
            if min_bc == -1:
                print("nothing remotely close")
        df = pd.DataFrame(list(zip(seqs, closest, dists)), columns = ['singletons', 'closest true barcode', 'distance'])
        df.to_csv('singletons.tsv', sep = '\t')
        print(min(dists))
        plt.hist(dists)
        plt.title("Minimum distance of each singleton")
        plt.show()
        
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
        
    def assigned_stats(self, barcodes, title):
        counts = []
        degree = []
        for bc in barcodes:
            node = self.numbered[bc]
            counts.append(self.counts[bc])
            degree.append(len(self.edges[node]))
        plt.hist(counts, bins = 30)
        plt.title(("Counts of" +title +"barcodes"))
        plt.show()
        plt.hist(degree, bins = 30)
        plt.title(("Degree of" +title +"barcodes"))
        plt.show()
        
    def assign_by_cluster(self, bc_len):
        observed_assignments = defaultdict(str)
        for node in self.counts.keys():
            if self.clustered[node] and self.clustering[node][0] != -1:
                bc = unrank(node, bc_len)
                tbc = unrank(self.clustering[node][0], bc_len)
                observed_assignments[bc] = tbc
        return observed_assignments
    
    def compare_results(self, true_assignment, true_barcodes):
        components = []
        visited = [False for node in self.barcodes.keys()]
        for node in self.barcodes.keys():
            #print(node)
            if not visited[node]:
                component, visited = dfs_without_recursion(visited, node, [], self.edges)
                components.append(component)
        observed_assignments = self.get_assignments(true_barcodes, components)
        #observed_assignments = self.assign_by_cluster()
        n_correct_assignments = 0
        n_correct_in_component = 0
        n_incorrect = 0
        n_unobserved = 0
        n_unassigned = 0
        unobserved_barcodes = set()
        distances_correct = []
        distances_in_component = []
        correctly_assigned = []
        incorrectly_assigned = set()
        correct_in_component = []
        incorrect_component_size = []
        unassigned = set()
        correctly_tbcs = []
        component_tbcs = []
        
        for component in components:
            #if len(component) > 10000:
            observed_true = 0
            for node in component:
                if self.barcodes[node] in true_barcodes:
                    observed_true += 1
            for node in component:
                bc = self.barcodes[node]
                for tbc in true_assignment[bc].keys():
                #for tbc in true_assignment[bc]:
                    if tbc in self.numbered.keys():
                        correct_bc = self.numbered[tbc]
                        #if correct_bc in self.edges[node]:
                            #n_correct_assignments += 1
                            #n_correct_in_component += 1
                        if tbc == observed_assignments[bc]:
                            n_correct_assignments += true_assignment[bc][tbc]
                            correctly_assigned.append(bc)
                            correctly_tbcs.append(tbc)
                            #n_correct_assignments += 1
                            for i in range(true_assignment[bc][tbc]):
                                distances_correct.append(editdistance.eval(bc, tbc))
                        elif correct_bc in component:
                            if observed_assignments[bc] == "":
                                n_unassigned += true_assignment[bc][tbc]
                                unassigned.add(bc)
                            else:
                                n_correct_in_component += true_assignment[bc][tbc]
                                # if n_correct_in_component <= 100:
                                # print("correct in component")
                                # print("barcode:", bc)
                                # print("assigned barcode:", observed_assignments[bc], "count", self.counts[bc], "dist", editdistance.eval(bc, observed_assignments[bc]))
                                # print("true barcode:", tbc, "count", true_assignment[bc][tbc], "dist", editdistance.eval(bc,tbc))
                                #n_correct_in_component += 1
                                for i in range(true_assignment[bc][tbc]):
                                    distances_in_component.append(editdistance.eval(bc, observed_assignments[bc]))
                                correct_in_component.append(bc)
                                component_tbcs.append(observed_assignments[bc])
                                #print("distance of bc and true bc:" , editdistance.eval(bc, tbc))
                                #print("distance of bc and observed assignment:" , editdistance.eval(bc, observed_assignments[bc]))
                        else:
                            #print("Correct barcode in different component")
                            #print("Barcode:", self.barcodes[node], "count:", self.counts[self.barcodes[node]])
                            #print("True assignment:", true_assignment[self.barcodes[node]])
                            #print("Size of the component:", len(component))
                            #print("Number of true barcodes in the component:", observed_true)
                            if observed_assignments[bc] == "":
                                n_unassigned += true_assignment[bc][tbc]
                                unassigned.add(bc)
                            else:
                                n_incorrect += true_assignment[bc][tbc]
                                # if n_incorrect <= 100:
                                # print("incorrect")
                                # print("barcode:", bc)
                                # print("assigned barcode:", observed_assignments[bc], "count", self.counts[bc], "dist", editdistance.eval(bc, observed_assignments[bc]))
                                # print("true barcode:", tbc, "count", true_assignment[bc][tbc], "dist", editdistance.eval(bc,tbc))
                            #n_incorrect += 1
                                incorrectly_assigned.add(bc)
                                incorrect_component_size.append(len(component))
                    else:
                        #print("correct barcode never observed")
                        unobserved_barcodes.add(tbc)
                        #print(tbc)
                        n_incorrect += 1
                        n_unobserved += 1
        print("adding count for each distinct barcode")
        #print("only for largest component")
        print("n_correct_in_component:", n_correct_in_component)
        print("n_correct_assignments:", n_correct_assignments)
        print("n_incorrect:", n_incorrect)
        print("n_unassigned:", n_unassigned)
        print("Number of never observed barcodes:", len(unobserved_barcodes))
        print("Number of times a barcode is unobserved:", n_unobserved)
        
        # correct = pd.DataFrame({"barcode": correctly_assigned,
                                # "cluster": correctly_tbcs})
        # correct.to_csv('correctly_3it.tsv', sep = '\t')
        # comp = pd.DataFrame({"barcode": correct_in_component,
                             # "cluster": component_tbcs})
        # comp.to_csv('component_3it.tsv', sep = '\t')

            
        #plt.hist(incorrect_component_size, bins = 50)
        #plt.title("Size of components containing bcs with true bc in another component")
        #plt.show()
        #plt.hist(distances_correct)
        #plt.title("Distance to closest true barcode in the component, assignment correct, with counts")
        #plt.show()
        #plt.hist(distances_in_component)
        #plt.title("Distance to closest true barcode in the component, assignment incorrect, true barcode in component, with counts")
        #plt.show()
        #self.assigned_stats(unassigned, " not assigned ")
                # correct_barcode = self.numbered[true_assignment[self.barcodes[node]]] # figure out what to do if not in the list, maybe start my numbering with 1 so 0 means it doesn't exist or can I change defaultdict to default = -1?
                # if correct_barcode in self.edges[node]:
                    # n_correct_assignments += 1
                    # n_correct_in_component += 1
                # elif correct_barcode in component:
                    # n_correct_in_component += 1
                # else:
                    # print("Correct barcode in different component")
                    # print("Barcode:", node, "count:", self.counts[self.barcodes[node]])
                    # print("True assignment:", correct_barcode, true_assignment[self.barcodes[node]])
                    
    def compare_to_cluster(self, true_barcodes, true_assignment):
        
        components = []
        visited = [False for node in self.barcodes.keys()]
        for node in self.barcodes.keys():
            #print(node)
            if not visited[node]:
                component, visited = dfs_without_recursion(visited, node, [], self.edges)
                components.append(component)
        cluster_assignment = self.assign_by_cluster()
        assignment = self.get_assignments(true_barcodes, components)
        for bc in self.numbered.keys():
            if cluster_assignment[bc] != assignment[bc]:
                print("barcode:", bc)
                print("cluster assignment:", cluster_assignment[bc], "dist:", editdistance.eval(bc, cluster_assignment[bc]))
                print("original assignment:", assignment[bc], "dist:", editdistance.eval(bc, assignment[bc]))
                if cluster_assignment[bc] in true_assignment[bc].keys():
                    print("Cluster assignment correct")
                if assignment[bc] in true_assignment[bc].keys():
                    print("Original assignment correct")
                    
    def true_barcode_stats(self, true_barcodes):
        counts = []
        f_counts = []
        degree = []
        f_degree = []
        component_size = []
        f_component_size = []
        components = []
        visited = [False for node in self.barcodes.keys()]
        for node in self.barcodes.keys():
            if not visited[node]:
                component, visited = dfs_without_recursion(visited, node, [], self.edges)
                components.append(component)
        print("Before")
        for component in components:
            #if len(component) < 30000:
            for node in component:
                if self.barcodes[node] in true_barcodes:
                    component_size.append(len(component))
                    counts.append(self.counts[self.barcodes[node]])
                    degree.append(len(self.edges[node]))
                else:
                    f_component_size.append(len(component))
                    f_counts.append(self.counts[self.barcodes[node]])
                    f_degree.append(len(self.edges[node]))
                    if self.counts[self.barcodes[node]] > 50:
                        print(self.barcodes[node], self.counts[self.barcodes[node]])
        #for bc in true_barcodes:
        #    counts.append(self.counts[bc])
        #    tbc = self.numbered[bc]
        #    degree.append(len(self.edges[tbc]))
        print("After")
        print("Minumum count:", min(counts))
        print("Minimum degree:", min(degree))
        plt.figure()
        plt.hist(counts, bins = 75, log = True)
        plt.hist(f_counts, bins = 14, color = "red", log = True)
        plt.title("Counts of all barcodes, true blue, other red")
        plt.savefig("counts_all_P7.10p.png")
        #plt.show()
        plt.figure()
        plt.hist(degree, bins = 25)
        plt.hist(f_degree, bins = 25, color = "red", log = True)
        plt.title("Degrees of all barcodes, true blue, other red")
        plt.savefig("degrees_all_P7.10p.png")
        #plt.show()
        plt.figure()
        plt.scatter(component_size, counts, alpha = 0.2)
        plt.scatter(f_component_size, f_counts, c = "red", alpha = 0.2)
        plt.title("Counts of barcodes by component size, true blue, other red")
        plt.savefig("counts_all_by_size_P7.10p.png")
        #plt.show()
        plt.figure()
        plt.scatter(component_size, degree, alpha = 0.2)
        plt.scatter(f_component_size, f_degree, c = "red", alpha = 0.2)
        plt.title("Degrees of barcodes by component size, true blue, other red")
        plt.savefig("degrees_all_by_size_P7.10p.png")
        #plt.show()
        plt.figure()
        plt.scatter(f_component_size, f_degree, c = "red")
        plt.title("Degrees of not true barcodes by component size")
        plt.savefig("degrees_nt_by_size_P7.10p.png")
        #plt.show()
        plt.figure()
        plt.scatter(f_component_size, f_counts, c = "red")
        plt.title("Counts of not true barcodes by component size")
        plt.savefig("counts_nt_by_size_P7.10p.png")
        #plt.show()
        plt.figure()
        plt.hist(f_counts, bins = 80, color = "red", log = True)
        plt.title("Counts of not true barcodes")
        plt.savefig("counts_nt_P7.10p.png")
        #plt.show()
        plt.figure()
        plt.hist(f_degree, bins = 25, color = "red")
        plt.title("Degrees of not true barcodes")
        plt.savefig("degrees_nt_P7.10p.png")
        #plt.show()
        plt.figure()
        plt.scatter(component_size, counts, c = "blue")
        plt.title("Degrees of true barcodes by component size")
        plt.savefig("degrees_t_by_size_P7.10p.png")
        #plt.show()
        plt.figure()
        plt.scatter(component_size, degree, c = "blue")
        plt.title("Counts of true barcodes by component size")
        plt.savefig("counts_t_by_size_P7.10p.png")
        #plt.show()
        plt.figure()
        plt.hist(counts, bins = 80, color = "blue", log = True)
        plt.title("Counts of true barcodes")
        plt.savefig("counts_t_P7.10p.png")
        #plt.show()
        plt.figure()
        plt.hist(degree, bins = 25, color = "blue")
        plt.title("Degrees of true barcodes")
        plt.savefig("degrees_t_P7.10p.png")
        #plt.show()
        
    def components_without_true(self, true_barcodes, true_assignment):
        sizes = []
        dists = []
        components = []
        num = 0
        comp = 0
        bc_to_comp = {}
        visited = [False for node in self.barcodes.keys()]
        for node in self.barcodes.keys():
            if not visited[node]:
                component, visited = dfs_without_recursion(visited, node, [], self.edges)
                components.append(component)
                min_dist = 32
                min_bc = -1
                for node in component:
                    bc = self.barcodes[node]
                    for tbc in true_barcodes:
                        dist = editdistance.eval(bc,tbc)
                        if dist < min_dist:
                            min_dist = dist 
                            min_bc = tbc
                        if dist == 0:
                            bc_to_comp[tbc] = component
                comp += 1
        for component in components:
            min_dist = 32
            min_bc = -1
            actual_bc = []
            for node in component:
                bc = self.barcodes[node]
                for tbc in true_assignment[bc].keys():
                    actual_bc.append(tbc)
                for tbc in true_barcodes:
                    dist = editdistance.eval(bc,tbc)
                    if dist < min_dist:
                        min_dist = dist 
                        min_bc = tbc
            if min_dist > 0:
                sizes.append(len(component))
                dists.append(min_dist)
                if len(component) > 10:
                    print("Component without true barcode")
                    print([self.barcodes[x] for x in component])
                    print("Closest true barcode:", min_bc)
                    print("Distance to closest true barcode:", min_dist)
                    print("Actual true barcodes:")
                    print(actual_bc)
                    #print("Component containing closest true")
                    #print([self.barcodes[x] for x in bc_to_comp[min_bc]])
                    ## investigate looking into what is actually the barcode something is connected to in these => is it the closest true for any of them? If not, extraction error?
                    num += 1
                    # G = nx.Graph()
                    # G.add_nodes_from(component)
                    # for node in component:
                        # for neighbour in self.edges[node]:
                            # G.add_edge(node, neighbour, weight = self.dists[(node, neighbour)] * 10)
                    # figname = "no_true_" + str(num) + ".png"
                    # plt.figure()
                    # nx.draw(G, node_size = 20)
                    # plt.savefig(figname)
        print("Number of components without true barcodes > 10:", num)
        # plt.hist(sizes, bins = 100, log = True)
        # plt.title("Sizes of components without true barcodes")
        # plt.show()
        # plt.hist(dists)
        # plt.title("Distances to closest true barcode from components without")
        # plt.show()
             
    def large_component(self, true_barcodes):
        components = []
        visited = [False for node in self.barcodes.keys()]
        subgraphs = []
        l_component = []
        for node in self.barcodes.keys():
            #print(node)
            if not visited[node]:
                component, visited = dfs_without_recursion(visited, node, [], self.edges)
                components.append(component)
                if len(component) > 10000:
                    l_component = component
                    break
        visited = [False for node in self.barcodes.keys()]
        #print(l_component)
        for tbc in true_barcodes: 
            tbc = self.numbered[tbc]
            if tbc in l_component:
                print("here")
                subgraph = []
                c = 100
                if not visited[tbc]:
                    stack = []
                    stack.append(tbc)
                    while stack and c > 0:
                        node = stack.pop()
                        if not visited[node]:
                            subgraph.append(node)
                            print("subgraph")
                            visited[node] = True
                            c = c - 1
                            for i in self.edges[node]:
                                stack.append(i)
                    bc_in_component = []
                    G = nx.Graph()
                    color_map = []
                    G.add_nodes_from(subgraph)
                    for node in subgraph:
                        for neighbour in self.edges[node]:
                            G.add_edge(node, neighbour, weight = self.dists[(node, neighbour)] * 10)
                    for node in G:
                        edge_node = False
                        for bc in self.edges[node]:
                            if not visited[bc]:
                                edge_node = True
                        if self.barcodes[node] in true_barcodes:
                            bc_in_component.append(self.barcodes[node])
                            color_map.append('red')
                        elif edge_node:
                            color_map.append('green')
                        else:
                            color_map.append('blue')
                    plt.figure()
                    nx.draw(G, node_size = 50, node_color = color_map)
                    plt.show()
            
            
    def choose_true(self, true_barcodes, barcode_list, n_cells):
        sorted_counts = dict(sorted(self.counts.items(), key=lambda item: item[1],reverse = True))
        bc_by_counts = list(sorted_counts.keys())
        maybe_true = []
        n = 0
        i = 0
        if barcode_list:
            while n < n_cells and i < len(bc_by_counts):
                if bc_by_counts[i] in barcode_list:
                    maybe_true.append(bc_by_counts[i])
                    n+=1
                i+=1
        else:
            maybe_true = bc_by_counts[:n_cells]
        wrong = 0
        for bc in true_barcodes:
            if not bc in barcode_list:
                print("True barcode not in barcode list")
                print(bc)
            if bc not in maybe_true:
                print("True barcode not included")
                print(bc, self.counts[bc])
                wrong += 1
        for bc in maybe_true:
            if bc not in true_barcodes:
                print("Barcode included but not true")
                print(bc, self.counts[bc])
                wrong += 1
        print(wrong)
                
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
                
    def output_file(self, read_assignment, out, true_barcodes, bc_len):
        
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
        
        for read in read_assignment.keys():
            observed_bcs = read_assignment[read]
            assigned_bcs = []
            for observed_bc in observed_bcs:
                assigned_bc = "*"
                if observed_bc != "*":
                    assigned_bc = assignments[observed_bc]
                    if assigned_bc == "":
                        assigned_bc = "*"
                assigned_bcs.append(assigned_bc)
                    
            read_ids.append(read)
            results.append(set(assigned_bcs))
        
        out_file = out + "_output_file.tsv"

        res = pd.DataFrame({"readID": read_ids,
                            "barcode": results})
        res.to_csv(out_file, sep = '\t', index = False)
        
        self.assignment = dict(zip(read_ids, results))
                
## /abga/work/rebpfeil/single_cell/read_group_files/


# later maybe switch dictionaries 
# edges maybe string to set

    def print_components(self, true_barcodes):
        components = []
        visited = [False for node in self.barcodes.keys()]
        for node in self.barcodes.keys():
            #print(node)
            if not visited[node]:
                component, visited = dfs_without_recursion(visited, node, [], self.edges)
                components.append(component)
        for component in components:
            if len(component) > 10 and len(component) < 1000:
                G = nx.Graph()
                color_map = []
                color_map_1 = []
                color_map_2 = []
                color_map_3 = []
                G.add_nodes_from(component)
                for node in component:
                    for neighbour in self.edges[node]:
                        G.add_edge(node, neighbour, weight = self.dists[(node, neighbour)] * 10)
                for node in G:
                    if self.barcodes[node] in true_barcodes:
                        color_map.append('red')
                        color_map_1.append('red')
                        color_map_2.append('red')
                        color_map_3.append('red')
                    elif self.clustered[node] and self.clustering[node][1] == 1:
                        color_map.append('blue')
                        color_map_1.append('limegreen')
                        color_map_2.append('limegreen')
                        color_map_3.append('limegreen')
                    elif self.clustered[node] and self.clustering[node][1] == 2:
                        color_map.append('blue')
                        color_map_1.append('blue')
                        color_map_2.append('darkorchid')
                        color_map_3.append('darkorchid')
                    elif self.clustered[node] and self.clustering[node][1] == -1:
                        color_map.append('blue')
                        color_map_1.append('blue')
                        color_map_2.append('blue')
                        color_map_3.append('darkorange')
                    else:
                        color_map.append('blue')
                        color_map_1.append('blue')
                        color_map_2.append('blue')
                        color_map_3.append('blue')
                    
                plt.figure()
                nx.draw(G, node_size = 50, node_color = color_map)
                plt.show()
                plt.figure()
                nx.draw(G, node_size = 50, node_color = color_map_1)
                plt.show()
                plt.figure()
                nx.draw(G, node_size = 50, node_color = color_map_2)
                plt.show()
                plt.figure()
                nx.draw(G, node_size = 50, node_color = color_map_3)
                plt.show()
                
    def result_statistics(self, true_assignment):
        
        lengths = []
        wrong = 0
        correct = 0
        unassigned = 0
        
        for read in true_assignment.keys():
            assigned = list(self.assignment[read])
            lengths.append(len(assigned))
            t = true_assignment[read]
            if len(assigned) == 1:
                if assigned[0] == "*":
                    unassigned += 1
                elif t == assigned[0]:
                    correct += 1
                else:
                    wrong += 1
            else:
                unass = True
                corr = False
                for a in assigned:
                    if a != "*":
                        unass = False
                    if a == t:
                        corr = True
                if corr:
                    correct += 1
                elif unass:
                    unassigned += 1
                else:
                    wrong += 1
        
        plt.hist(lengths)
        plt.title("Number of barcodes (including *) in the assignment")
        plt.yscale("log")
        plt.show()
                        
        print("assignment statistics:")
        print("correctly assigned:", correct)
        print("incorrectly assigned:", wrong)
        print("unassigned:", unassigned)
