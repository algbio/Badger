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

from barcode_extraction.molecule_structure import MoleculeStructure, ElementType
from barcode_graph import BarcodeGraph
from extract_raw_barcodes import (
    process_in_parallel,
    process_single_thread,
    ListReadHandler,
    ListHandlerGenerator,
    BARCODE_CALLING_MODES,
    BarcodeCallingModes
)
from support import load_true_barcodes, load_extracted_barcodes

logger = logging.getLogger('Badger')


def parse_args(args):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", "-i", help = "reads in FASTQ/FASTA (can be gzipped), BAM "
                                                "or TSV from barcode extraction",
                        type = str, dest = "input", required = True)
    parser.add_argument("--data_type", "-d", type=str, help="protocol used",
                        choices=[x.name for x in BarcodeCallingModes],
                        required=True)
    parser.add_argument("--molecule", type=str, help="file with molecule description format"
                                                     "(effective only when data type is set to custom)")

    parser.add_argument("--barcode_list", "-l",
                        help="List of all possible barcodes for the used method, helps identify correct barcodes",
                        type=str, dest="barcode_list", default=None)
    parser.add_argument("--true_barcodes",
                        help="List of all true barcodes of the input data, for example obtained from short read data",
                        type=str, default=None)
    parser.add_argument("--n_cells", "-c", help="expected number of cell associated barcodes",
                        type=int, default=5000)
    parser.add_argument("--output", "-o", help = "File prefix for output files",
                        type = str)

    parser.add_argument("--threshold", help = "Maximal accepted difference between barcodes",
                        type = int, dest = "threshold", default = 1)
    parser.add_argument("--ground_truth", help = "File connecting each observed barcode to its read ID containing true barcode, only used for statistics",
                        type = str, default = None)

    parser.add_argument("--cell_count_variance",
                        help="Percentage by which the number of cells is allowed to differ from estimated cell number, default 25%%",
                        default=25, type=int)
    parser.add_argument("--stats", action='store_true', help = "if set, true barcode statistics are run instead of barcode calling.", default = False)
    parser.add_argument("--threads", "-t", dest = "threads", default = 1, type = int)
    parser.add_argument("--high_sens", action='store_true',
                        help = "if set, Badger is run in high sensitivity mode. This increases recall but decreases precision", default = False)
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

    args.data_type = BarcodeCallingModes[args.data_type]
    if args.data_type != BarcodeCallingModes.custom and args.molecule:
        logger.warning("You set %s mode, but also provided a molecule structure file %s. "
                       "Molecule structure file will have not effect, set mode to %s to use it." %
                       (args.data_type, args.molecule, BarcodeCallingModes.custom))

    molecule_structure = None
    if args.data_type == BarcodeCallingModes.custom:
        molecule_structure = MoleculeStructure(open(args.molecule))
        barcode_detector = BARCODE_CALLING_MODES[args.data_type](molecule_structure)
    else:
        barcode_detector = BARCODE_CALLING_MODES[args.data_type]()
    logger.info("Barcode caller created")

    if args.input.endswith("tsv"):
        read_assignments = load_extracted_barcodes(args.input)
        logger.info("Imported barcodes from file")
    else:
        main_read_handler = ListReadHandler()
        if args.threads == 1:
            process_single_thread(args.input, barcode_detector, main_read_handler)
        else:
            handler_generator = ListHandlerGenerator()
            process_in_parallel(args, barcode_detector, handler_generator, main_read_handler)
        read_assignments = main_read_handler.read_storage

    if args.data_type != BarcodeCallingModes.custom:
        detect_barcodes_simple(args, read_assignments, barcode_detector)
    elif molecule_structure is not None:
        for element in molecule_structure:
            if element.element_type == ElementType.VAR_FILE:
                detect_barcodes_simple(args, read_assignments, barcode_detector, element.element_name)
            elif element.element_type == ElementType.VAR_LIST:
                detect_barcodes_simple(args, read_assignments, barcode_detector, element.element_name)


def detect_barcodes_simple(args, read_assignments, barcode_detector, element_keyword="barcode"):
    true_barcodes = load_true_barcodes(args.true_barcodes)

    if args.barcode_list:
        list_file = open(args.barcode_list, "r")
        barcode_list = list_file.read()
        barcode_list = set(barcode_list.split("\n"))
        list_file.close()
    else:
        barcode_list = None

    extracted_barcodes = list(filter(lambda x: x != '*',
                                     [x.detected_results[element_keyword].seq
                                      if element_keyword in x.detected_results
                                      else '*' for x in read_assignments]))

    logger.info("Initializing Graph")
    bc_len = barcode_detector.barcode_length
    graph = BarcodeGraph(args.threshold)
    graph.graph_construction(extracted_barcodes, bc_len, args.threads)
    logger.info("Graph construction done")

    out = args.output
    if not args.stats:
        graph.cluster(true_barcodes, barcode_list, args.n_cells, bc_len, args.cell_count_variance)
        logger.info("Clustering done")

        graph.output_file(read_assignments, out, barcode_detector, args.high_sens, element_keyword)
    
    #disconnected = len(graph.counts.keys()) - len(graph.edges.keys())
    #print(disconnected)
            
    if args.stats:
        import stats
        logger.info("Statistics being calculated")
        tbcs = graph.get_cluster_centers(None, bc_len, barcode_list, args.n_cells, args.cell_count_variance)
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
