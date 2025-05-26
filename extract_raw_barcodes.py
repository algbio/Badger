#!/usr/bin/env python3
#
############################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################
import concurrent
import os
import random
import sys
import argparse
import gzip
from traceback import print_exc
import shutil
from concurrent.futures import ProcessPoolExecutor
from collections import defaultdict
from enum import Enum, unique

import pysam
from Bio import SeqIO
import logging

from barcode_extraction.barcode_callers import (
    TenXBarcodeExtractorV2,
    TenXBarcodeExtractorV3,
    ReadStats
)
from barcode_extraction.universal_extraction import (
    UniversalSingleMoleculeExtractor,
    MoleculeStructure
)

logger = logging.getLogger('BarcodeGraph')


@unique
class BarcodeCallingModes(Enum):
    custom = 0
    tenX_v2 = 12
    tenX_v3 = 13
    # tenX_5v2 = 15


READ_CHUNK_SIZE = 100000
BARCODE_CALLING_MODES = {BarcodeCallingModes.tenX_v2: TenXBarcodeExtractorV2,
                         BarcodeCallingModes.tenX_v3: TenXBarcodeExtractorV3,
                         BarcodeCallingModes.custom: UniversalSingleMoleculeExtractor}


class FileReadHandler:
    def __init__(self, outfile, formatter):
        self.output_table = outfile
        self.output_file = open(self.output_table, "w")
        self.formatter = formatter

    def add_header(self):
        self.output_file.write(self.formatter.header() + "\n")

    def add_read(self, barcode_result):
        self.output_file.write(self.formatter.format_result(barcode_result) + "\n")

    def dump_stats(self, read_stat):
        stat_out = open(self.output_table + ".stats", "w")
        stat_out.write(str(read_stat))
        stat_out.close()

    def __del__(self):
        self.output_file.close()


class ListReadHandler:
    def __init__(self):
        self.read_storage = []

    def add_header(self, header):
        pass

    def add_read(self, barcode_result):
        self.read_storage.append((barcode_result.read_id, barcode_result.barcode, barcode_result.UMI))

    def dump_stats(self, read_stat):
        pass


class BarcodeCaller:
    def __init__(self, barcode_detector, read_handler):
        self.barcode_detector = barcode_detector
        self.read_handler = read_handler
        self.read_handler.add_header()
        self.read_stat = ReadStats()

    def process(self, input_file):
        logger.info("Processing " + input_file)
        fname, outer_ext = os.path.splitext(os.path.basename(input_file))
        low_ext = outer_ext.lower()

        handle = input_file
        if low_ext in ['.gz', '.gzip']:
            handle = gzip.open(input_file, "rt")
            input_file = fname
            fname, outer_ext = os.path.splitext(os.path.basename(input_file))
            low_ext = outer_ext.lower()

        if low_ext in ['.fq', '.fastq']:
            self._process_fastx(SeqIO.parse(handle, "fastq"))
        elif low_ext in ['.fa', '.fasta']:
            self._process_fastx(SeqIO.parse(handle, "fasta"))
        elif low_ext in ['.bam', '.sam']:
            self._process_bam(pysam.AlignmentFile(input_file, "rb"))
        else:
            logger.error("Unknown file format " + input_file)
        logger.info("Finished " + input_file)

    def _process_fastx(self, input_read_iterator):
        counter = 0
        for r in input_read_iterator:
            if counter % 100 == 0:
                sys.stdout.write("Processed %d reads\r" % counter)
            counter += 1
            read_id = r.id
            seq = str(r.seq)
            self._process_read(read_id, seq)

    def _process_bam(self, input_read_iterator):
        counter = 0
        for r in input_read_iterator:
            if counter % 100 == 0:
                sys.stdout.write("Processed %d reads\r" % counter)
            counter += 1
            read_id = r.query_name
            seq = r.query_sequence
            self._process_read(read_id, seq)

    def _process_read(self, read_id, read_sequence):
        logger.debug("==== %s ====" % read_id)
        barcode_result = self.barcode_detector.find_patterns(read_id, read_sequence)
        self.read_handler.add_read(barcode_result)
        self.read_stat.add_read(barcode_result)

    def process_chunk(self, read_chunk):
        for read_id, seq in read_chunk:
            self._process_read(read_id, seq)


def fastx_file_chunk_reader(handler):
    current_chunk = []
    for r in handler:
        current_chunk.append((r.id, str(r.seq)))
        if len(current_chunk) >= READ_CHUNK_SIZE:
            yield current_chunk
            current_chunk = []
    yield current_chunk


def bam_file_chunk_reader(handler):
    current_chunk = []
    for r in handler:
        if r.is_secondary or r.is_supplementary:
            continue
        current_chunk.append((r.query_name, r.query_sequence))
        if len(current_chunk) >= READ_CHUNK_SIZE:
            yield current_chunk
            current_chunk = []
    yield current_chunk


def process_chunk(barcode_detector, read_chunk, output_file, num):
    output_file += "_" + str(num)
    read_handler = FileReadHandler(output_file, barcode_detector)
    barcode_caller = BarcodeCaller(barcode_detector, read_handler)
    barcode_caller.process_chunk(read_chunk)
    read_handler.dump_stats(barcode_caller.read_stat)
    return output_file


def process_single_thread(args):
    logger.info("Processing " + args.input)
    if args.mode == BarcodeCallingModes.custom:
        molecule_structure = MoleculeStructure(open(args.molecule))
        barcode_detector = BARCODE_CALLING_MODES[args.mode](molecule_structure)
    else:
        barcode_detector = BARCODE_CALLING_MODES[args.mode]()
    read_handler = FileReadHandler(args.output, barcode_detector)
    barcode_caller = BarcodeCaller(barcode_detector, read_handler)
    barcode_caller.process(args.input)

    read_handler.dump_stats(barcode_caller.read_stat)
    for l in str(barcode_caller.read_stat).split("\n"):
        if l:
            logger.info(l)
    logger.info("Finished barcode calling")


def process_in_parallel(args):
    input_file = args.input
    logger.info("Processing " + input_file)
    fname, outer_ext = os.path.splitext(os.path.basename(input_file))
    low_ext = outer_ext.lower()

    handle = input_file
    if low_ext in ['.gz', '.gzip']:
        handle = gzip.open(input_file, "rt")
        input_file = fname
        fname, outer_ext = os.path.splitext(os.path.basename(input_file))
        low_ext = outer_ext.lower()

    if low_ext in ['.fq', '.fastq']:
        read_chunk_gen = fastx_file_chunk_reader(SeqIO.parse(handle, "fastq"))
    elif low_ext in ['.fa', '.fasta']:
        read_chunk_gen = fastx_file_chunk_reader(SeqIO.parse(handle, "fasta"))
    elif low_ext in ['.bam', '.sam']:
        read_chunk_gen = bam_file_chunk_reader(pysam.AlignmentFile(input_file, "rb"))
    else:
        logger.error("Unknown file format " + input_file)
        exit(-1)

    tmp_dir = "barcode_calling_%x" % random.randint(0, 1 << 32)
    while os.path.exists(tmp_dir):
        tmp_dir = "barcode_calling_%x" % random.randint(0, 1 << 32)
    if args.tmp_dir:
        tmp_dir = os.path.join(args.tmp_dir, tmp_dir)
    os.makedirs(tmp_dir)

    tmp_barcode_file = os.path.join(tmp_dir, "bc")
    count = 0
    future_results = []
    output_files = []

    if args.mode == BarcodeCallingModes.custom:
        molecule_structure = MoleculeStructure(open(args.molecule))
        barcode_detector = BARCODE_CALLING_MODES[args.mode](molecule_structure)
    else:
        barcode_detector = BARCODE_CALLING_MODES[args.mode]()
    logger.info("Barcode caller created")

    with ProcessPoolExecutor(max_workers=args.threads) as proc:
        for chunk in read_chunk_gen:
            future_results.append(proc.submit(process_chunk, barcode_detector, chunk, tmp_barcode_file, count))
            count += 1
            if count >= args.threads:
                break

        reads_left = True
        while reads_left:
            completed_features, _ = concurrent.futures.wait(future_results, return_when=concurrent.futures.FIRST_COMPLETED)
            for c in completed_features:
                if c.exception() is not None:
                    raise c.exception()
                future_results.remove(c)
                output_files.append(c.result())
                if reads_left:
                    try:
                        chunk = next(read_chunk_gen)
                        future_results.append(proc.submit(process_chunk, barcode_detector, chunk, tmp_barcode_file, count))
                        count += 1
                    except StopIteration:
                        reads_left = False

        completed_features, _ = concurrent.futures.wait(future_results, return_when=concurrent.futures.ALL_COMPLETED)
        for c in completed_features:
            if c.exception() is not None:
                raise c.exception()
            output_files.append(c.result())

    outf = open(args.output, "w")
    stat_dict = defaultdict(int)
    for tmp_file in output_files:
        shutil.copyfileobj(open(tmp_file, "r"), outf)
        if not os.path.exists(tmp_file + ".stats"):
            logger.warning("Stats fiel %s was no found" % tmp_file + ".stats")
            continue
        for l in open(tmp_file + ".stats", "r"):
            v = l.strip().split("\t")
            if len(v) != 2:
                continue
            stat_dict[v[0]] += int(v[1])

    out_stats = open(args.output + ".stats", "w")
    for k, v in stat_dict.items():
        logger.info("%s %d" % (k, v))
        out_stats.write("%s %d\n" % (k, v))
    shutil.rmtree(tmp_dir)
    logger.info("Finished barcode calling")


def extract_barcodes_from_chunk(barcode_detector, read_chunk):
    read_handler = ListReadHandler()
    barcode_caller = BarcodeCaller(barcode_detector, read_handler)
    barcode_caller.process_chunk(read_chunk)
    return read_handler.read_storage


def extract_barcodes_single_thread(input_file, mode):
    logger.info("Extracting from " + input_file)
    barcode_detector = BARCODE_CALLING_MODES[mode]()
    read_handler = ListReadHandler()
    barcode_caller = BarcodeCaller(barcode_detector, read_handler)
    barcode_caller.process(input_file)
    logger.info("Finished barcode extraction")
    return read_handler.read_storage


def extract_barcodes_in_parallel(input_file, mode, threads):
    logger.info("Extracting from " + input_file)
    fname, outer_ext = os.path.splitext(os.path.basename(input_file))
    low_ext = outer_ext.lower()

    handle = input_file
    if low_ext in ['.gz', '.gzip']:
        handle = gzip.open(input_file, "rt")
        input_file = fname
        fname, outer_ext = os.path.splitext(os.path.basename(input_file))
        low_ext = outer_ext.lower()

    if low_ext in ['.fq', '.fastq']:
        read_chunk_gen = fastx_file_chunk_reader(SeqIO.parse(handle, "fastq"))
    elif low_ext in ['.fa', '.fasta']:
        read_chunk_gen = fastx_file_chunk_reader(SeqIO.parse(handle, "fasta"))
    elif low_ext in ['.bam', '.sam']:
        read_chunk_gen = bam_file_chunk_reader(pysam.AlignmentFile(input_file, "rb"))
    else:
        logger.error("Unknown file format " + input_file)
        exit(-1)

    count = 0
    future_results = []
    barcode_detector = BARCODE_CALLING_MODES[mode]()
    barcoded_reads = []
    logger.info("Barcode caller created")

    with ProcessPoolExecutor(max_workers=threads) as proc:
        for chunk in read_chunk_gen:
            future_results.append(proc.submit(extract_barcodes_from_chunk, barcode_detector, chunk))
            count += 1
            if count >= threads:
                break

        reads_left = True
        while reads_left:
            completed_features, _ = concurrent.futures.wait(future_results, return_when=concurrent.futures.FIRST_COMPLETED)
            for c in completed_features:
                if c.exception() is not None:
                    raise c.exception()
                future_results.remove(c)
                barcoded_reads += c.result()
                if reads_left:
                    try:
                        chunk = next(read_chunk_gen)
                        future_results.append(proc.submit(extract_barcodes_from_chunk, barcode_detector, chunk))
                        count += 1
                    except StopIteration:
                        reads_left = False

        completed_features, _ = concurrent.futures.wait(future_results, return_when=concurrent.futures.ALL_COMPLETED)
        for c in completed_features:
            if c.exception() is not None:
                raise c.exception()
            barcoded_reads += c.result()

    logger.info("Finished barcode extraction")
    return barcoded_reads


def load_barcodes(inf):
    barcode_list = []
    for l in open(inf):
        barcode_list.append(l.strip().split()[0])
    return barcode_list


def set_logger(logger_instance):
    logger_instance.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger_instance.addHandler(ch)


def parse_args(sys_argv):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", "-o", type=str, help="output prefix name", required=True)
    parser.add_argument("--mode", type=str, help="mode to be used",
                        choices=[x.name for x in BarcodeCallingModes],
                        required=True)
    parser.add_argument("--molecule", type=str, help="file with molecule description format")
    parser.add_argument("--input", "-i", type=str, help="input reads in [gzipped] FASTA, FASTQ, BAM, SAM",
                        required=True)
    parser.add_argument("--threads", "-t", type=int, help="threads to use (16)", default=16)
    parser.add_argument("--tmp_dir", type=str, help="folder for temporary files")

    args = parser.parse_args(sys_argv)
    return args


def main(sys_argv):
    args = parse_args(sys_argv)
    set_logger(logger)

    args.mode = BarcodeCallingModes[args.mode]
    if args.mode != BarcodeCallingModes.custom and args.molecule:
        logger.warning("You set %s mode, but also provided a molecule structure file %s. "
                       "Molecule structure file will have not effect, set mode to %s to use it." %
                       (args.mode, args.molecule, BarcodeCallingModes.custom))

    if args.threads == 1:
        process_single_thread(args)
    else:
        process_in_parallel(args)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main(sys.argv[1:])
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
