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
from Bio import SeqIO

def compare_results(true_assignment, graph_assignment, blaze_assignment):
        correct_graph = 0
        correct_blaze = 0
        wrong_graph = 0
        wrong_blaze = 0
        unassigned_blaze = 0
        unassigned_graph = 0
        
        
        reads = true_assignment.keys()
        for read in reads:
            graph = graph_assignment[read]
            blaze = blaze_assignment[read]
            t = true_assignment[read]
            if graph == "" or graph == "*":
                unassigned_graph += 1
            elif t == graph:
                correct_graph += 1
            else:
                wrong_graph += 1
                #if t == blaze: 
                    #print("Badger wrong, Blaze correct")
                    #print("Read:", read)
                    #print("Badger:", graph)
                    #print("Blaze:", blaze)
            if blaze == "":
                unassigned_blaze += 1
            elif t == blaze:
                correct_blaze += 1
            else:
                wrong_blaze += 1
            
            
        print("graph assignment statistics:")
        print("correctly assigned:", correct_graph)
        print("incorrectly assigned:", wrong_graph)
        print("unassigned:", unassigned_graph)
        print("---------------------------------------------------------")
        print("Blaze statistics:")
        print("correctly assigned:", correct_blaze)
        print("incorrectly assigned:", wrong_blaze)
        print("unassigned:", unassigned_blaze)
        
def dict_from_blaze(blaze_data):
    blaze_dict = defaultdict(str)
    
    for read in blaze_data:
        rid = read.id
        r = rid.split('_')
        bc = r[0]
        #print(bc)
        rid = rid.split('#')
        readid = rid[1][:-2]
        #print(readid)
        
        blaze_dict[readid] = bc
    
    return blaze_dict

def dict_from_graph(graph_df):
    graph_dict = defaultdict(str)
    for i in graph_df.index:
        read = graph_df.loc[i,'readID']
        assignment = graph_df.loc[i,'barcode']
        graph_dict[read] = assignment
    return graph_dict
    
blaze_data = SeqIO.parse(sys.argv[1], "fastq")

graph_df = pd.read_csv(sys.argv[2], sep = "\t")

graph_dict = dict_from_graph(graph_df)

blaze_dict = dict_from_blaze(blaze_data)

reads = pd.read_csv(sys.argv[3], sep = "\t", header = None)
ids = reads.iloc[1:,0].tolist()
true_assignment = defaultdict(str)
for i in range(len(ids)): 
    if ids[i] != "#read_id":
        true_bc = ids[i].split('_')[3]
        if true_bc == "PAR":
            true_bc = ids[i].split('_')[5]
        true_assignment[ids[i]] = true_bc

compare_results(true_assignment, graph_dict, blaze_dict)
