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

from index import QGramIndex
from common import get_score, dfs, dfs_without_recursion

logger = logging.getLogger("BarcodeGraph")

class BarcodeGraph:
    
    def __init__(self, threshold):
        
        self.threshold = threshold
        self.counts = dict() # basically also a list of my nodes
        self.edges = defaultdict(list) # idea: change list to set
        self.dists = defaultdict(int) # I feel like I should have one dict that for every node contains all other nodes it is connected to and one dict having edges as keys and the distance of the endpoints as values
        # might change as I move on but right now this feels correct
        self.barcodes = defaultdict(str)
        self.numbered = defaultdict(int)
        self.clusters = defaultdict(list)
        self.clustering = dict()
        self.clustered = []
    
    def graph_construction(self, barcodes, bc_len):
        #think about cutting one base off the right
        index = QGramIndex(self.threshold, bc_len, 6)
        num = 0
        
        logger.info("Using edit distance")
        
        for sequence in barcodes:
            if len(sequence) == bc_len + 1:
                sequence = sequence[:-1] 
            if len(sequence) == bc_len:
                if sequence in self.counts.keys():
                    self.counts[sequence] +=1
                else:
                    self.barcodes[num] = sequence
                    self.numbered[sequence] = num
                    self.counts[sequence] = 1
                    index.add_to_index(sequence, num)
                    #for s in self.barcodes.keys():
                    #    seq = self.barcodes[s]
                    closest = index.get_close(sequence, num)
                    for s in closest:
                        seq = self.barcodes[s]
                        if seq == "" or seq == sequence:
                            continue
                        else:
                            # test: min of three distances to possibly ignore last position for indels
                            #dist = editdistance.eval(sequence, seq)
                            dist = min(editdistance.eval(sequence, seq), editdistance.eval(sequence[:-1],seq), editdistance.eval(sequence,seq[:-1]))
                            #score = get_score(sequence, seq)
                            if dist <= self.threshold:
                            #if score >= 16*3 - 4: ## fully equal is 3*len, 2 indels would be -4
                                self.edges[num].append(s)
                                self.edges[s].append(num)
                                # #self.dists[(sequence,seq)] = dist # think if I should add both sides or one is enough
                                self.dists[(num,s)] = dist 
                                self.dists[(s,num)] = dist
                    num += 1
                    if num%500000 == 0:
                        logger.info(f"processed {num} distinct barcodes")
                
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
        
    def cluster(self, true_barcodes, barcode_list, n_cells):
        # might need another dict keeping cluster centre for every node
        # what if in this dict I keep bc: (tbc, round in which it was added to the cluster)
        # then when I want to add to a cluster I check whether it is already in a cluster and if so if it was added in an earlier round
        # if it was added in a previous round leave it where it is, if it is the same round discard bc from cluster
        self.clustered = [False for node in self.barcodes.keys()]
        
        sorted_counts = dict(sorted(self.counts.items(), key=lambda item: item[1],reverse = True))
        bc_by_counts = list(sorted_counts.keys())
        tbcs = []
        n = 0
        i = 0
        if true_barcodes:
            tbcs = true_barcodes
        elif barcode_list:
            while n < n_cells and i < len(bc_by_counts):
                if bc_by_counts[i] in barcode_list:
                    tbcs.append(bc_by_counts[i])
                    n+=1
                i+=1
        else:
            tbcs = bc_by_counts[:n_cells]
        
        # do I iterate through true barcodes or components? Also how do I decide number of rounds for clustering? Do I just choose a cutoff or make it depend on the graph?
        for tbc in tbcs:
            tbc = self.numbered[tbc]
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
                            #interesting: check that tbc is not the same before discarding any nodes, can be connected to the same tbc over two different other barcodes
                            self.clusters[center].append(neighbor)
                            self.clustering[neighbor] = (center,i)
                            self.clustered[neighbor] = True 
                        elif self.clustering[neighbor][0] != center and self.clustering[neighbor][0] != -1:
                            if self.clustering[neighbor][1] == i:
                                self.clusters[self.clustering[neighbor][0]].remove(neighbor)
                                self.clustering[neighbor] = (-1, -1)
        
    def graph_statistics(self, true_barcodes):
        components = []
        singletons = []
        lengths = []
        # I want to try only visualizing large components
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
                    # I want to make statistics if they have one node that is connected to most others
                    # maybe also check on the counts if one is seen more often
                    # and when I find one that is either connected to most others or seen much more often compare it to the true barcodes
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
        
    def assign_by_cluster(self):
        observed_assignments = defaultdict(str)
        for node in self.barcodes.keys():
            if self.clustered[node] and self.clustering[node][0] != -1:
                bc = self.barcodes[node]
                tbc = self.barcodes[self.clustering[node][0]]
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
                        # careful: I think the list might contain every correct barcode as many times as it is in the reads. This is technically what we want, just remember that
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
                                # think about the case of equal edit distance
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
        #plt.savefig("counts_all_P7.10p.png")
        plt.show()
        plt.figure()
        plt.hist(degree, bins = 25)
        plt.hist(f_degree, bins = 25, color = "red", log = True)
        plt.title("Degrees of all barcodes, true blue, other red")
        #plt.savefig("degrees_all_P7.10p.png")
        plt.show()
        plt.figure()
        plt.scatter(component_size, counts, alpha = 0.2)
        plt.scatter(f_component_size, f_counts, c = "red", alpha = 0.2)
        plt.title("Counts of barcodes by component size, true blue, other red")
        #plt.savefig("counts_all_by_size_P7.10p.png")
        plt.show()
        plt.figure()
        plt.scatter(component_size, degree, alpha = 0.2)
        plt.scatter(f_component_size, f_degree, c = "red", alpha = 0.2)
        plt.title("Degrees of barcodes by component size, true blue, other red")
        #plt.savefig("degrees_all_by_size_P7.10p.png")
        plt.show()
        plt.figure()
        plt.scatter(f_component_size, f_degree, c = "red")
        plt.title("Degrees of not true barcodes by component size")
        #plt.savefig("degrees_nt_by_size_P7.10p.png")
        plt.show()
        plt.figure()
        plt.scatter(f_component_size, f_counts, c = "red")
        plt.title("Counts of not true barcodes by component size")
        #plt.savefig("counts_nt_by_size_P7.10p.png")
        plt.show()
        plt.figure()
        plt.hist(f_counts, bins = 80, color = "red", log = True)
        plt.title("Counts of not true barcodes")
        #plt.savefig("counts_nt_P7.10p.png")
        plt.show()
        plt.figure()
        plt.hist(f_degree, bins = 25, color = "red")
        plt.title("Degrees of not true barcodes")
        #plt.savefig("degrees_nt_P7.10p.png")
        plt.show()
        plt.figure()
        plt.scatter(component_size, counts, c = "blue")
        plt.title("Degrees of true barcodes by component size")
        #plt.savefig("degrees_t_by_size_P7.10p.png")
        plt.show()
        plt.figure()
        plt.scatter(component_size, degree, c = "blue")
        plt.title("Counts of true barcodes by component size")
        #plt.savefig("counts_t_by_size_P7.10p.png")
        plt.show()
        plt.figure()
        plt.hist(counts, bins = 80, color = "blue", log = True)
        plt.title("Counts of true barcodes")
        #plt.savefig("counts_t_P7.10p.png")
        plt.show()
        plt.figure()
        plt.hist(degree, bins = 25, color = "blue")
        plt.title("Degrees of true barcodes")
        #plt.savefig("degrees_t_P7.10p.png")
        plt.show()
        
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
                
    def output_file(self, read_assignment, out, true_barcodes):
        
        components = []
        visited = [False for node in self.barcodes.keys()]
        for node in self.barcodes.keys():
            if not visited[node]:
                component, visited = dfs_without_recursion(visited, node, [], self.edges)
                components.append(component)
        
        #assignments = self.assign_by_cluster()
        assignments = self.get_assignments(true_barcodes, components)
        read_ids = []
        results = []
        
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
    
