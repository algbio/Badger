#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2024 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import pandas as pd
import argparse
import logging
import sys
from collections import defaultdict
from io import StringIO
from traceback import print_exc

from barcode_graph import BarcodeGraph

logger = logging.getLogger('BarcodeGraph')

def parse_args(args):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    #parser.add_argument("--barcodes", "-b", help="tsv file containing the observed cell barcodes", ##this at some point I can just get from the barcode extraction output directly, but do that after thesis
    #                    type = str, dest = "bar_file", required = True)
    parser.add_argument("--threshold", "-t", help = "Maximal accepted difference between barcodes",
                        type = int, dest = "threshold", default = 1) 
    parser.add_argument("--reads", "-r", help = "output file of barcode extraction algorithm",
                        type = str, dest = "reads", required = True)
    parser.add_argument("--ground_truth", help = "File connecting each observed barcode to its read ID containing true barcode, only used for statistics",
                        type = str, default = None)
    parser.add_argument("--barcode_list", "-l", help = "List of all possible barcodes for the used method, helps identify correct barcodes",
                        type = str, dest = "barcode_list", default = None)
    parser.add_argument("--data_type", "-d", help = "Type of single cell sequencing data in the input, options are 10x and Double", 
                        choices = ["10x", "Visium"], type = str)
    parser.add_argument("--true_barcodes", help = "List of all true barcodes of the input data, for example obtained from short read data",
                        type = str, default = None)
    parser.add_argument("--n_cells", "-c", help = "expected number of cell associated barcodes",
                        type = int, default = 5000)
    parser.add_argument("--output", "-o", help = "File prefix for output files",
                        type = str, default = "OUT")
    parser.add_argument("--stats", "-s", help = "if set, true barcode statistics are run instead of barcode calling.")
    return parser.parse_args(args)
    
def set_logger(logger_instance):
    
    logger_instance.setLevel(logging.INFO)
    c_handler = logging.StreamHandler(stream=sys.stdout)
    c_handler.setLevel(logging.INFO)
    
    c_format = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    c_handler.setFormatter(c_format)
    
    logger_instance.addHandler(c_handler)

    logger_instance.info("Starting")

def main(args):
    args = parse_args(args)
    set_logger(logger)
    bc_len = 0
    if args.data_type == "10x":
        bc_len = 16
    elif args.data_type == "Double":
        bc_len = 20
    else:
        logger.error("Please specify the type of single cell data used. Options are 10x and Double.")
        exit(-3)
    # barcodes = pd.read_csv(args.bar_file, sep = "\t", header = None)
    # barcodes = barcodes.dropna()
    # barcodes = barcodes.iloc[:,0].tolist()
    true_barcodes = args.true_barcodes
    if true_barcodes:
        true_barcodes = pd.read_csv(true_barcodes, sep = "\t", header = None)
        #true_barcodes = pd.read_csv('10x_barcodes_5K.tsv', sep = "\t", header = None)
        #true_barcodes = pd.read_csv('barcodes.tsv', sep = "\t", header = None)
        true_barcodes = true_barcodes.iloc[:,0].tolist()
        if true_barcodes[0][-1] == '1':
            for i in range(len(true_barcodes)):
                true_barcodes[i] = true_barcodes[i][:-2]
        true_barcodes = set(true_barcodes)
    
    if args.barcode_list:
        list_file = open(args.barcode_list, "r")
        barcode_list = list_file.read()
        barcode_list = set(barcode_list.split("\n"))
        list_file.close()
    else:
        barcode_list = None
    
    out = args.output
    reads = pd.read_csv(args.reads, sep = "\t")
    ids = reads["#read_id"].tolist()
    observed = reads["barcode"]
    observed = observed.dropna()
    observed = observed.tolist()
    read_assignment = []
    barcodes = reads["barcode"]
    barcodes = barcodes.dropna()
    barcodes = barcodes[barcodes != "*"]
    barcodes = barcodes[barcodes != "barcode"]
    barcodes = barcodes.tolist()
    for i in range(len(ids)):
        if ids[i] != "#read_id":
            di = ids[i]
            o = observed[i]
            if o != "barcode":
                if len(o) == bc_len + 1:
                    o = o[:-1]
                read_assignment.append((di, o))
    logger.info("Imported barcodes from file")
    
    logger.info("Initializing Graph")
    graph = BarcodeGraph(args.threshold)
    graph.graph_construction(barcodes, bc_len)
    logger.info("Graph construction done")
    
    if not args.stats:
        graph.cluster(true_barcodes, barcode_list, args.n_cells)
        logger.info("Clustering done")
        
        
        graph.output_file(read_assignment, out, true_barcodes)
    
    disconnected = len(graph.counts.keys()) - len(graph.edges.keys())
    print(disconnected)
            
    if args.stats:
        logger.info("Statistics being calculated")
        #graph.graph_statistics(true_barcodes)
        graph.choose_true(true_barcodes, barcode_list, args.n_cells)
        #graph.visualize_graph()
        graph.true_barcode_stats(true_barcodes)
        #graph.large_component(true_barcodes)
        #graph.print_components(true_barcodes)
    
    
    if args.ground_truth is not None:
        truth = pd.read_csv(args.ground_truth, sep = "\t", header = None)
        #print(reads.iloc[:,0])
        ids = truth.iloc[1:,0].tolist()
        observed = truth.iloc[1:,1].tolist()
        read_assignment = []
        true_assignment = defaultdict(dict)
        seen = defaultdict(set)
        for i in range(len(ids)): 
            if ids[i] != "#read_id":
                true_bc = ids[i].split('_')[3]
                if true_bc == "PAR":
                    true_bc = ids[i].split('_')[5]
                observed_bc = observed[i]
                if observed_bc != "barcode" and observed_bc != "*":
                    read_assignment.append((ids[i], true_bc, observed_bc[:-1]))
                    observed_bc = observed_bc[:-1]
                    if true_bc in seen[observed_bc]:
                        true_assignment[observed_bc][true_bc] += 1
                    else:
                        true_assignment[observed_bc][true_bc] = 1
                        seen[observed_bc].add(true_bc)
        #for key in true_assignment.keys():
            #print(key, true_assignment[key])
        if true_barcodes:
            #graph.components_without_true(true_barcodes, true_assignment)
            #graph.compare_to_cluster(true_barcodes, true_assignment)
            #graph.isoquant_output(read_assignment, true_barcodes)
            graph.compare_results(true_assignment, true_barcodes)
    

if __name__ == "__main__":
    try:
        main(sys.argv[1:])
    except SystemExit:
        raise
    except KeyboardInterrupt:
        raise
    except:
        if logger.handlers:
            strout = StringIO()
            print_exc(file=strout)
            s = strout.getvalue()
            if s:
                logger.critical("Barcode Graph failed" +s)
            else:
                print_exc()
        else:
            sys.stderr.write("Barcode Graph failed")
            print_exc()
        sys.exit(-1)
