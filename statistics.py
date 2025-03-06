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

def visualize_graph(graph):
    edges = []
    for edge in graph.dists.keys():
        edges.append([edge[0], edge[1]])
    g = ig.Graph(n = len(graph.counts.keys()), edges = edges)
    unconnected = []
    for i in range(0, len(graph.counts.keys())):
        if g.degree(i) == 0:
            unconnected.append(i)
    g.delete_vertices(unconnected) 
    layout = g.layout("fr")
    ig.plot(g, target = 'graph_dist_2.pdf', vertex_size = 1, layout = layout, edge_color = ['red', 'black'])
    
def graph_statistics(graph, true_barcodes, bc_len):
    components = []
    singletons = []
    lengths = []
    false_components = []
    visited = [False for node in graph.counts.keys()]
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
    for node in graph.counts.keys():
        if not visited[node]:
            component, visited = dfs_without_recursion(visited, node, [], graph.edges)
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
                        # for neighbour in graph.edges[node]:
                            # G.add_edge(node, neighbour, weight = graph.dists[(node, neighbour)] * 10)
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
                    if len(graph.edges[node]) == 0:
                        print(component)
                        print("Singleton?")
                    if len(graph.edges[node]) > max_degree:
                        max_degree = len(graph.edges[node])
                        max_degree_node = node
                    if graph.counts[node] > max_count:
                        max_count = graph.counts[node]
                        max_count_node = node
                    for bc in true_barcodes:
                        dist = editdistance.eval(unrank(node, bc_len), bc)
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
                if len(graph.edges[max_count_node]) == max_degree:
                    max_degree_node = max_count_node
                # #best_is_max += (min_node == max_degree_node or min_node == max_count_node)
                #print("New component")
                #print("Number of nodes:", len(component))
                #print("Node with highest degree:", max_degree_node, "degree:", max_degree)
                #print("Count of the node with the highest degree:", graph.counts[graph.barcodes[max_degree_node]])
                #print("Node with highest count:", max_count_node, "count:", max_count)
                #print("Degree of the node with the highest count:", len(graph.edges[max_count_node]))
                min_bc_c = -1
                min_bc_d = -1
                min_dist_c = 32
                min_dist_d = 32
                for bc in true_barcodes:
                    dist_c = editdistance.eval(unrank(max_count_node, bc_len), bc)
                    dist_d = editdistance.eval(unrank(max_degree_node, bc_len), bc)
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
                        for neighbour in graph.edges[node]:
                            G.add_edge(node, neighbour, weight = graph.dists[(node, neighbour)] * 10)
                    for node in G:
                        if graph.barcodes[node] in true_barcodes:
                            bc_in_component.append(graph.barcodes[node])
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
        single_counts.append(graph.counts[node])
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
    # graph.closest_true(singletons, true_barcodes)
    
def closest_true(graph, singletons, true_barcodes, bc_len):
    true_barcodes = set(true_barcodes)
    closest = []
    dists = []
    seqs = []
    for node in singletons:
        min_bc = -1
        min_dist = 32
        for bc in true_barcodes:
            dist = editdistance.eval(unrank(node, bc_len), bc)
            if dist < min_dist:
                min_bc = bc
                min_dist = dist
        closest.append(min_bc)
        dists.append(min_dist)
        seqs.append(unrank(node, bc_len))
        if min_bc == -1:
            print("nothing remotely close")
    df = pd.DataFrame(list(zip(seqs, closest, dists)), columns = ['singletons', 'closest true barcode', 'distance'])
    df.to_csv('singletons.tsv', sep = '\t')
    print(min(dists))
    plt.hist(dists)
    plt.title("Minimum distance of each singleton")
    plt.show()

def assigned_stats(graph, barcodes, title, bc_len):
    counts = []
    degree = []
    for bc in barcodes:
        node = rank(bc, bc_len)
        counts.append(graph.counts[node])
        degree.append(len(graph.edges[node]))
    plt.hist(counts, bins = 30)
    plt.title(("Counts of" +title +"barcodes"))
    plt.show()
    plt.hist(degree, bins = 30)
    plt.title(("Degree of" +title +"barcodes"))
    plt.show()
    
def compare_results(graph, true_assignment, true_barcodes, bc_len):
    components = []
    visited = [False for node in graph.counts.keys()]
    for node in graph.counts.keys():
        #print(node)
        if not visited[node]:
            component, visited = dfs_without_recursion(visited, node, [], graph.edges)
            components.append(component)
    #observed_assignments = graph.get_assignments(true_barcodes, components)
    observed_assignments = graph.assign_by_cluster()
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
            if unrank(node, bc_len) in true_barcodes:
                observed_true += 1
        for node in component:
            bc = unrank(node, bc_len)
            for tbc in true_assignment[bc].keys():
            #for tbc in true_assignment[bc]:
                if rank(tbc, bc_len) in graph.counts.keys():
                    correct_bc = rank(tbc, bc_len)
                    #if correct_bc in graph.edges[node]:
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
                            # print("assigned barcode:", observed_assignments[bc], "count", graph.counts[bc], "dist", editdistance.eval(bc, observed_assignments[bc]))
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
                        #print("Barcode:", unrank(node, bc_len), "count:", graph.counts[graph.barcodes[node]])
                        #print("True assignment:", true_assignment[unrank(node, bc_len)])
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
                            # print("assigned barcode:", observed_assignments[bc], "count", graph.counts[bc], "dist", editdistance.eval(bc, observed_assignments[bc]))
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
    #graph.assigned_stats(unassigned, " not assigned ")
            # correct_barcode = graph.numbered[true_assignment[graph.barcodes[node]]] # figure out what to do if not in the list, maybe start my numbering with 1 so 0 means it doesn't exist or can I change defaultdict to default = -1?
            # if correct_barcode in graph.edges[node]:
                # n_correct_assignments += 1
                # n_correct_in_component += 1
            # elif correct_barcode in component:
                # n_correct_in_component += 1
            # else:
                # print("Correct barcode in different component")
                # print("Barcode:", node, "count:", graph.counts[graph.barcodes[node]])
                # print("True assignment:", correct_barcode, true_assignment[graph.barcodes[node]])
                
def compare_to_cluster(graph, true_barcodes, true_assignment, bc_len):
    
    components = []
    visited = [False for node in graph.counts.keys()]
    for node in graph.counts.keys():
        #print(node)
        if not visited[node]:
            component, visited = dfs_without_recursion(visited, node, [], graph.edges)
            components.append(component)
    cluster_assignment = graph.assign_by_cluster()
    assignment = graph.get_assignments(true_barcodes, components)
    for node in graph.counts.keys():
        bc = unrank(node, bc_len)
        if cluster_assignment[bc] != assignment[bc]:
            print("barcode:", bc)
            print("cluster assignment:", cluster_assignment[bc], "dist:", editdistance.eval(bc, cluster_assignment[bc]))
            print("original assignment:", assignment[bc], "dist:", editdistance.eval(bc, assignment[bc]))
            if cluster_assignment[bc] in true_assignment[bc].keys():
                print("Cluster assignment correct")
            if assignment[bc] in true_assignment[bc].keys():
                print("Original assignment correct")
                
def true_barcode_stats(graph, true_barcodes, bc_len):
    counts = []
    f_counts = []
    degree = []
    f_degree = []
    component_size = []
    f_component_size = []
    components = []
    visited = [False for node in graph.counts.keys()]
    for node in graph.counts.keys():
        if not visited[node]:
            component, visited = dfs_without_recursion(visited, node, [], graph.edges)
            components.append(component)
    print("Before")
    for component in components:
        #if len(component) < 30000:
        for node in component:
            if unrank(node, bc_len) in true_barcodes:
                component_size.append(len(component))
                counts.append(graph.counts[node])
                degree.append(len(graph.edges[node]))
            else:
                f_component_size.append(len(component))
                f_counts.append(graph.counts[node])
                f_degree.append(len(graph.edges[node]))
                if graph.counts[node] > 50:
                    print(unrank(node, bc_len), graph.counts[node])
    #for bc in true_barcodes:
    #    tbc = rank(bc, bc_len)
    #    counts.append(graph.counts[tbc])
    #    degree.append(len(graph.edges[tbc]))
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
    
def components_without_true(graph, true_barcodes, true_assignment, bc_len):
    sizes = []
    dists = []
    components = []
    num = 0
    comp = 0
    bc_to_comp = {}
    visited = [False for node in graph.counts.keys()]
    for node in graph.counts.keys():
        if not visited[node]:
            component, visited = dfs_without_recursion(visited, node, [], graph.edges)
            components.append(component)
            min_dist = 32
            min_bc = -1
            for node in component:
                bc = unrank(node, bc_len)
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
            bc = unrank(node, bc_len)
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
                print([graph.barcodes[x] for x in component])
                print("Closest true barcode:", min_bc)
                print("Distance to closest true barcode:", min_dist)
                print("Actual true barcodes:")
                print(actual_bc)
                #print("Component containing closest true")
                #print([graph.barcodes[x] for x in bc_to_comp[min_bc]])
                ## investigate looking into what is actually the barcode something is connected to in these => is it the closest true for any of them? If not, extraction error?
                num += 1
                # G = nx.Graph()
                # G.add_nodes_from(component)
                # for node in component:
                    # for neighbour in graph.edges[node]:
                        # G.add_edge(node, neighbour, weight = graph.dists[(node, neighbour)] * 10)
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
    
def large_component(graph, true_barcodes, bc_len):
    components = []
    visited = [False for node in graph.counts.keys()]
    subgraphs = []
    l_component = []
    for node in graph.counts.keys():
        #print(node)
        if not visited[node]:
            component, visited = dfs_without_recursion(visited, node, [], graph.edges)
            components.append(component)
            if len(component) > 10000:
                l_component = component
                break
    visited = [False for node in graph.barcodes.keys()]
    #print(l_component)
    for tbc in true_barcodes: 
        tbc = rank(node, bc_len)
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
                        for i in graph.edges[node]:
                            stack.append(i)
                bc_in_component = []
                G = nx.Graph()
                color_map = []
                G.add_nodes_from(subgraph)
                for node in subgraph:
                    for neighbour in graph.edges[node]:
                        G.add_edge(node, neighbour, weight = graph.dists[(node, neighbour)] * 10)
                for node in G:
                    edge_node = False
                    for bc in graph.edges[node]:
                        if not visited[bc]:
                            edge_node = True
                    if graph.barcodes[node] in true_barcodes:
                        bc_in_component.append(graph.barcodes[node])
                        color_map.append('red')
                    elif edge_node:
                        color_map.append('green')
                    else:
                        color_map.append('blue')
                plt.figure()
                nx.draw(G, node_size = 50, node_color = color_map)
                plt.show()
                
def choose_true(graph, true_barcodes, barcode_list, n_cells, bc_len):
    sorted_counts = dict(sorted(graph.counts.items(), key=lambda item: item[1],reverse = True))
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
    for tbc in true_barcodes:
        bc = unrank(tbc, bc_len)
        if not bc in barcode_list:
            print("True barcode not in barcode list")
            print(bc)
        if bc not in maybe_true:
            print("True barcode not included")
            print(bc, graph.counts[tbc])
            wrong += 1
    for bc in maybe_true:
        if bc not in true_barcodes:
            print("Barcode included but not true")
            print(bc, graph.counts[tbc])
            wrong += 1
    print(wrong)
    
def print_components(graph, true_barcodes):
    components = []
    visited = [False for node in graph.counts.keys()]
    for node in graph.counts.keys():
        #print(node)
        if not visited[node]:
            component, visited = dfs_without_recursion(visited, node, [], graph.edges)
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
                for neighbour in graph.edges[node]:
                    G.add_edge(node, neighbour, weight = graph.dists[(node, neighbour)] * 10)
            for node in G:
                if graph.barcodes[node] in true_barcodes:
                    color_map.append('red')
                    color_map_1.append('red')
                    color_map_2.append('red')
                    color_map_3.append('red')
                elif graph.clustered[node] and graph.clustering[node][1] == 1:
                    color_map.append('blue')
                    color_map_1.append('limegreen')
                    color_map_2.append('limegreen')
                    color_map_3.append('limegreen')
                elif graph.clustered[node] and graph.clustering[node][1] == 2:
                    color_map.append('blue')
                    color_map_1.append('blue')
                    color_map_2.append('darkorchid')
                    color_map_3.append('darkorchid')
                elif graph.clustered[node] and graph.clustering[node][1] == -1:
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
