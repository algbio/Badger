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
import editdistance
import matplotlib.pyplot as plt
from itertools import islice

def compare_results(true_assignment, graph_assignment, blaze_assignment, read_assignment):
        correct_graph = 0
        correct_blaze = 0
        wrong_graph = 0
        wrong_blaze = 0
        unassigned_blaze = 0
        unassigned_graph = 0
        bl2e = []
        ba2t = []
        ba2e = []
        extracted2true = []
        unassignedDists = []
        unextracted = 0
        
        
        reads = read_assignment.keys()
        for read in reads:
            graph = graph_assignment[read]
            blaze = blaze_assignment[read]
            extracted = read_assignment[read]#[:-1]
            t = true_assignment[read]
            if extracted =="*":
                unextracted += 1
            if graph == "" or graph == "*":
                unassigned_graph += 1
                if extracted != "*":
                    #print("extracted barcode: ", extracted)
                    #print("blaze barcode: ", blaze)
                    #print("correct barcode: ", t)
                    dist = editdistance.eval(extracted, t)
                    #print("distance to true: ", dist)
                    unassignedDists.append(dist)
                    
            elif t == graph:
                correct_graph += 1
            else:
                wrong_graph += 1
                #if t == blaze: 
                    # print("Badger wrong, Blaze correct")
                    # print("Read:", read)
                    # print("Badger:", graph)
                    # print("Blaze:", blaze)
                    # print("True:", t)
                    # print("extracted:", extracted)
                    # badger2true = editdistance.eval(graph, t)
                    # blaze2extracted = editdistance.eval(blaze, extracted)
                    # badger2extracted = editdistance.eval(graph, extracted)
                    # bl2e.append(blaze2extracted)
                    # ba2t.append(badger2true)
                    # ba2e.append(badger2extracted)
                    # extracted2true.append(editdistance.eval(extracted, t))
                    # print("Distance of Badger to true:", badger2true)
                    # print("Distance of Badger to extracted:", badger2extracted)
                    # print("Distance of Blaze to extracted:", blaze2extracted)
                    # print("Distance of true to extracted:", editdistance.eval(extracted, t))
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
        print("unextracted reads:", unextracted)
        
        # plt.hist(ba2t)
        # plt.title("Distances of Badger assignment to truth, blaze is correct")
        # plt.show()
        # plt.hist(bl2e)
        # plt.title("Distances of Blaze assignment to extracted, blaze is correct")
        # plt.show()
        # plt.hist(ba2e)
        # plt.title("Distances of Badger assignment to extracted, blaze is correct")
        # plt.show()
        # plt.hist(extracted2true)
        # plt.title("Distances of truth to extracted, blaze is correct")
        # plt.show()
        plt.hist(unassignedDists)
        plt.title("Distances of extracted to truth, Badger is unassigned")
        plt.show()
        
def dict_from_blaze(blaze_data):
    blaze_dict = defaultdict(str)
    
    for read in blaze_data:
    #for read in islice(blaze_data, 20):
        rid = read.id
        r = rid.split('_')
        bc = r[0]
        #print(bc)
        rid = rid.split('#')
        readid = rid[1][:-2]
        #print(rid[1][:-2])
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

reads = pd.read_csv(sys.argv[3], sep = "\t")
ids = reads.iloc[:,0].tolist()
observed = reads["barcode"]
observed = observed.fillna('*')
observed = observed.tolist()
read_assignment = defaultdict(str)
true_assignment = defaultdict(str)
unextracted = 0
for i in range(len(ids)): 
    if ids[i] != "#read_id":
        true_bc = ids[i].split('_')[3]
        if true_bc == "PAR":
            true_bc = ids[i].split('_')[5]
        true_assignment[ids[i]] = true_bc
        read_assignment[ids[i]] = observed[i]
        #print(ids[i])
        if observed[i] == '*':
            unextracted += 1

print("unextracted:", unextracted)

compare_results(true_assignment, graph_dict, blaze_dict, read_assignment)
