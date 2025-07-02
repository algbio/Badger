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

from barcode_extraction.molecule_structure import MoleculeStructure
from barcode_graph import BarcodeGraph
from extract_raw_barcodes import (
    process_in_parallel,
    process_single_thread,
    ListReadHandler,
    ListHandlerGenerator,
    BARCODE_CALLING_MODES,
    BarcodeCallingModes
)
import stats
from support import load_true_barcodes, load_extracted_barcodes

logger = logging.getLogger('Badger')


def parse_args(args):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    #parser.add_argument("--barcodes", "-b", help="tsv file containing the observed cell barcodes", ##this at some point I can just get from the barcode extraction output directly, but do that after thesis
    #                    type = str, dest = "bar_file", required = True)
    parser.add_argument("--threshold", "-t", help = "Maximal accepted difference between barcodes",
                        type = int, dest = "threshold", default = 1) 
    parser.add_argument("--input", "-i", help = "reads in FASTQ/FASTA (can be gzipped), BAM or TSV from barcode extraction",
                        type = str, dest = "input", required = True)
    parser.add_argument("--ground_truth", help = "File connecting each observed barcode to its read ID containing true barcode, only used for statistics",
                        type = str, default = None)
    parser.add_argument("--barcode_list", "-l", help = "List of all possible barcodes for the used method, helps identify correct barcodes",
                        type = str, dest = "barcode_list", default = None)
    parser.add_argument("--mode", type=str, help="mode to be used",
                        choices=[x.name for x in BarcodeCallingModes],
                        required=True)
    parser.add_argument("--molecule", type=str, help="file with molecule description format")
    parser.add_argument("--true_barcodes", help = "List of all true barcodes of the input data, for example obtained from short read data",
                        type = str, default = None)
    parser.add_argument("--n_cells", "-c", help = "expected number of cell associated barcodes",
                        type = int, default = 5000)
    parser.add_argument("--output", "-o", help = "File prefix for output files",
                        type = str, default = "OUT")
    parser.add_argument("--interval", "-i", help = "Percentage by which the number of cells is allowed to differ from estimated cell number, default 25%", default = 25, type = int)
    parser.add_argument("--stats", "-s", action='store_true', help = "if set, true barcode statistics are run instead of barcode calling.", default = False)
    parser.add_argument("--threads", "-tr", dest = "threads", default = 1, type = int)
    parser.add_argument("--high_sens", "-hs", action='store_true', help = "if set, Badger is run in high sensitivity mode. This increases recall but decreases precision", default = False)
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

    args.mode = BarcodeCallingModes[args.mode]
    if args.mode != BarcodeCallingModes.custom and args.molecule:
        logger.warning("You set %s mode, but also provided a molecule structure file %s. "
                       "Molecule structure file will have not effect, set mode to %s to use it." %
                       (args.mode, args.molecule, BarcodeCallingModes.custom))

    if args.mode == BarcodeCallingModes.custom:
        molecule_structure = MoleculeStructure(open(args.molecule))
        barcode_detector = BARCODE_CALLING_MODES[args.mode](molecule_structure)
    else:
        barcode_detector = BARCODE_CALLING_MODES[args.mode]()
    logger.info("Barcode caller created")

    if args.input.endswith("tsv"):
        # FIXME
        read_assignment = load_extracted_barcodes(args.input)
        logger.info("Imported barcodes from file")
    else:
        main_read_handler = ListReadHandler()
        if args.threads == 1:
            process_single_thread(args.input, barcode_detector, main_read_handler)
        else:
            handler_generator = ListHandlerGenerator()
            process_in_parallel(args, barcode_detector, handler_generator, main_read_handler)
        read_assignment = main_read_handler.read_storage
    # FIXME
    barcodes = list(filter(lambda x: x != "*", (ra[1] for ra in read_assignment)))


def detect_barcodes_simple(mode, args):
    true_barcodes = load_true_barcodes(args.true_barcodes)

    if args.barcode_list:
        list_file = open(args.barcode_list, "r")
        barcode_list = list_file.read()
        barcode_list = set(barcode_list.split("\n"))
        list_file.close()
    else:
        barcode_list = None
    
    out = args.output

    bc_len = 0
    if args.data_type.startswith("tenX"):
        bc_len = 16
    elif args.data_type == "Double":
        bc_len = 20
    else:
        logger.error("Please specify the type of single cell data used. Options are tenX_v2, tenX_v3 and Double.")
        exit(-3)

    logger.info("Initializing Graph")
    graph = BarcodeGraph(args.threshold)
    graph.graph_construction(barcodes, bc_len, args.threads)
    logger.info("Graph construction done")
    
    if not args.stats:
        graph.cluster(true_barcodes, barcode_list, args.n_cells, bc_len, args.interval)
        logger.info("Clustering done")
        
        
        graph.output_file(read_assignment, out, true_barcodes, bc_len, args.high_sens)
    
    disconnected = len(graph.counts.keys()) - len(graph.edges.keys())
    print(disconnected)
            
    if args.stats:
        logger.info("Statistics being calculated")
        tbcs = graph.get_cluster_centers(None, bc_len, barcode_list, args.n_cells, args.interval)
        stats.evaluate_centers(graph, tbcs, true_barcodes, barcode_list, bc_len)
        #stats.graph_statistics(graph, true_barcodes)
        #stats.choose_true(graph, true_barcodes, barcode_list, args.n_cells)
        #stats.visualize_graph(graph)
        stats.true_barcode_stats(graph, true_barcodes, bc_len)
        #stats.large_component(graph, true_barcodes)
        #stats.print_components(graph, true_barcodes)
    
    
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
