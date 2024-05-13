#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2024 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import pandas as pd
import sys
from collections import defaultdict

def compare_results(true_assignment, graph_assignment, scTagger_assignment):
        correct_graph = 0
        correct_1_scTagger = 0
        correct_list_scTagger = 0
        wrong_graph = 0
        wrong_scTagger = 0
        unassigned_scTagger = 0
        unassigned_graph = 0
        
        
        reads = true_assignment.keys()
        for read in reads:
            graph = graph_assignment[read]
            scTagger = scTagger_assignment[read]
            t = true_assignment[read]
            if graph == "" or graph == "*":
                unassigned_graph += 1
            elif t == graph:
                correct_graph += 1
            else:
                wrong_graph += 1
            if scTagger == []:
                unassigned_scTagger += 1
            elif t in scTagger:
                if len(scTagger) == 1:
                    correct_1_scTagger += 1
                else:
                    correct_list_scTagger += 1
            else:
                wrong_scTagger += 1
            
            
        print("graph assignment statistics:")
        print("correctly assigned:", correct_graph)
        print("incorrectly assigned:", wrong_graph)
        print("unassigned:", unassigned_graph)
        print("---------------------------------------------------------")
        print("scTagger statistics:")
        print("correctly assigned:", correct_1_scTagger)
        print("correct barcode in list:", correct_list_scTagger)
        print("incorrectly assigned:", wrong_scTagger)
        print("unassigned:", unassigned_scTagger)
        
def dict_from_scTagger(scTagger_df):
    scTagger_dict = defaultdict(list)
    for i in scTagger_df.index:
        read = scTagger_df.iloc[i,0]
        assignments = scTagger_df.iloc[i,4]
        assignments = assignments.split(',')
        scTagger_dict[read] = assignments
    return scTagger_dict

def dict_from_graph(graph_df):
    graph_dict = defaultdict(str)
    for i in graph_df.index:
        read = graph_df.loc[i,'readID']
        assignment = graph_df.loc[i,'barcode']
        graph_dict[read] = assignment
    return graph_dict

scTagger_df = pd.read_csv(sys.argv[1], sep = "\t")
graph_df = pd.read_csv(sys.argv[2], sep = "\t")

scTagger_dict = dict_from_scTagger(scTagger_df)
graph_dict = dict_from_graph(graph_df)

reads = pd.read_csv(sys.argv[3], sep = "\t", header = None)
ids = reads.iloc[1:,0].tolist()
true_assignment = defaultdict(str)
seen = defaultdict(set)
for i in range(len(ids)): 
    if ids[i] != "#read_id":
        true_bc = ids[i].split('_')[3]
        if true_bc == "PAR":
            true_bc = ids[i].split('_')[5]
        true_assignment[ids[i]] = true_bc

compare_results(true_assignment, graph_dict, scTagger_dict)
