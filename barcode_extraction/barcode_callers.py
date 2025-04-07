###########################################################################
# Copyright (c) 2023 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
from collections import defaultdict
from enum import Enum, unique

from .kmer_indexer import KmerIndexer
from .common import find_polyt_start, reverese_complement, detect_exact_positions

logger = logging.getLogger('BarcodeGraph')


class BarcodeDetectionResult:
    NOSEQ = "*"

    def __init__(self, read_id, barcode=NOSEQ, UMI=NOSEQ, BC_score=-1, UMI_good=False, strand=".", additional_info=None):
        self.read_id = read_id
        self.barcode = barcode
        self.UMI = UMI
        self.BC_score = BC_score
        self.UMI_good = UMI_good
        self.strand = strand

    def is_valid(self):
        return self.barcode != BarcodeDetectionResult.NOSEQ

    def more_informative_than(self, that):
        raise NotImplemented()

    def get_additional_attributes(self):
        raise NotImplemented()

    def set_strand(self, strand):
        self.strand = strand

    def __str__(self):
        return "%s\t%s\t%s\t%d\t%s\t%s" % (self.read_id, self.barcode, self.UMI,
                                           self.BC_score, self.UMI_good, self.strand)

    def __getstate__(self):
        return (self.read_id,
                self.barcode,
                self.UMI,
                self.BC_score,
                self.UMI_good,
                self.strand)

    def __setstate__(self, state):
        self.read_id = state[0]
        self.barcode = state[1]
        self.UMI = state[2]
        self.BC_score = state[3]
        self.UMI_good = state[4]
        self.strand = state[5]

    @staticmethod
    def header():
        return "#read_id\tbarcode\tUMI\tBC_score\tvalid_UMI\tstrand"


class TenXBarcodeDetectionResult(BarcodeDetectionResult):
    def __init__(self, read_id, barcode=BarcodeDetectionResult.NOSEQ, UMI=BarcodeDetectionResult.NOSEQ,
                 BC_score=-1, UMI_good=False, strand=".",
                 polyT=-1, r1=-1, r1_score=0):
        BarcodeDetectionResult.__init__(self, read_id, barcode, UMI, BC_score, UMI_good, strand)
        self.r1 = r1
        self.polyT = polyT
        self.r1_score = r1_score

    def is_valid(self):
        return self.barcode != BarcodeDetectionResult.NOSEQ

    def more_informative_than(self, that):
        return self.r1_score > that.r1_score

    def get_additional_attributes(self):
        attr = []
        if self.polyT != -1:
            attr.append("PolyT detected")
        if self.r1 != -1:
            attr.append("R1 detected")
        return attr

    def set_strand(self, strand):
        self.strand = strand

    def __str__(self):
        return (BarcodeDetectionResult.__str__(self) +
                "\t%d\t%d" % (self.polyT, self.r1))

    def __getstate__(self):
        return (self.read_id,
                self.barcode,
                self.UMI,
                self.BC_score,
                self.UMI_good,
                self.strand,
                self.r1,
                self.polyT,
                self.r1_score)

    def __setstate__(self, state):
        self.read_id = state[0]
        self.barcode = state[1]
        self.UMI = state[2]
        self.BC_score = state[3]
        self.UMI_good = state[4]
        self.strand = state[5]
        self.r1 = state[6]
        self.polyT = state[7]
        self.r1_score = state[8]

    @staticmethod
    def header():
        return BarcodeDetectionResult.header() + "\tpolyT_start\tR1_end"


class ReadStats:
    def __init__(self):
        self.read_count = 0
        self.bc_count = 0
        self.umi_count = 0
        self.additional_attributes_counts = defaultdict(int)

    def add_read(self, barcode_detection_result):
        self.read_count += 1
        for a in barcode_detection_result.get_additional_attributes():
            self.additional_attributes_counts[a] += 1
        if barcode_detection_result.barcode != BarcodeDetectionResult.NOSEQ:
            self.bc_count += 1
        if barcode_detection_result.UMI_good:
            self.umi_count += 1

    def __str__(self):
        human_readable_str =  ("Total reads:\t%d\nBarcode detected:\t%d\nReliable UMI:\t%d\n" %
                      (self.read_count, self.bc_count, self.umi_count))
        for a in self.additional_attributes_counts:
            human_readable_str += "%s:\t%d\n" % (a, self.additional_attributes_counts[a])
        return human_readable_str


@unique
class TenXVersions(Enum):
     v2 = 2
     v3 = 3


class TenXBarcodeExtractor:
    TSO = "CCCATGTACTCTGCGTTGATACCACTGCTT"
    R1 = "CTACACGACGCTCTTCCGATCT" # 10x 3'
    BARCODE_LEN_10X = 16
    UMI_LENGTHS = {TenXVersions.v2: 10, TenXVersions.v3: 12}

    TERMINAL_MATCH_DELTA = 4
    STRICT_TERMINAL_MATCH_DELTA = 1

    def __init__(self, protocol_version=TenXVersions.v3):
        self.r1_indexer = KmerIndexer([TenXBarcodeExtractor.R1], kmer_size=6)
        self.UMI_LEN_10X = self.UMI_LENGTHS[protocol_version]

    def find_barcode_umi(self, read_id, sequence):
        read_result = self._find_barcode_umi_fwd(read_id, sequence)
        if read_result.polyT != -1:
            read_result.set_strand("+")

        rev_seq = reverese_complement(sequence)
        read_rev_result = self._find_barcode_umi_fwd(read_id, rev_seq)
        if read_rev_result.polyT != -1:
            read_rev_result.set_strand("-")

        if read_rev_result.is_valid() and read_result.is_valid():
            return read_result if read_result.more_informative_than(read_rev_result) else read_rev_result
        if read_rev_result.is_valid():
            return read_rev_result
        return read_result

    def _find_barcode_umi_fwd(self, read_id, sequence):
        logger.debug("== read id %s ==" % read_id)
        polyt_start = find_polyt_start(sequence)
        logger.debug("PolyT %d" % polyt_start)
        r1_start, r1_end, r1_score = None, None, 0
        if polyt_start != -1:
            # use relaxed parameters is polyA is found
            r1_occurrences = self.r1_indexer.get_occurrences(sequence[0:polyt_start + 1])
            r1_start, r1_end, r1_score = detect_exact_positions(sequence, 0, polyt_start + 1,
                                                                self.r1_indexer.k, self.R1,
                                                                r1_occurrences, min_score=9,
                                                                end_delta=self.TERMINAL_MATCH_DELTA)

        logger.debug("R1 %s" % str(r1_start))
        if r1_start is None:
            # if polyT was not found, or linker was not found to the left of polyT, look for linker in the entire read
            r1_occurrences = self.r1_indexer.get_occurrences(sequence)
            r1_start, r1_end, r1_score = detect_exact_positions(sequence, 0, len(sequence),
                                                                self.r1_indexer.k, self.R1,
                                                                r1_occurrences, min_score=17,
                                                                start_delta=self.STRICT_TERMINAL_MATCH_DELTA,
                                                                end_delta=self.STRICT_TERMINAL_MATCH_DELTA)

        if r1_start is None:
            return TenXBarcodeDetectionResult(read_id, polyT=polyt_start)
        logger.debug("LINKER: %d-%d" % (r1_start, r1_end))

        if polyt_start != -1 and polyt_start - r1_end < self.BARCODE_LEN_10X:
            return TenXBarcodeDetectionResult(read_id, polyT=polyt_start)

        if polyt_start == -1 or polyt_start - r1_end > self.BARCODE_LEN_10X + self.UMI_LEN_10X + 10:
            # if polyT was not detected earlier, use relaxed parameters once the linker is found
            presumable_polyt_start = r1_end + self.BARCODE_LEN_10X + self.UMI_LEN_10X
            search_start = presumable_polyt_start - 4
            search_end = min(len(sequence), presumable_polyt_start + 10)
            polyt_start = find_polyt_start(sequence[search_start:search_end], window_size=5, polya_fraction=1.0)
            if polyt_start != -1:
                polyt_start += search_start

        barcode_start = r1_end + 1
        barcode_end = r1_end + self.BARCODE_LEN_10X
        potential_barcode = sequence[barcode_start:barcode_end + 1]
        logger.debug("Barcode: %s" % (potential_barcode))
        potential_umi_start = barcode_end + 1
        potential_umi_end = polyt_start - 1
        if potential_umi_end - potential_umi_start <= 5:
            potential_umi_end = potential_umi_start + self.UMI_LEN_10X - 1
        potential_umi = sequence[potential_umi_start:potential_umi_end + 1]
        return TenXBarcodeDetectionResult(read_id, potential_barcode, potential_umi, BC_score=0, polyT=polyt_start, r1=r1_end, r1_score=r1_score)

    def find_barcode_umi_no_polya(self, read_id, sequence):
        logger.debug("===== FWD =====")
        read_result = self._find_barcode_umi_fwd(read_id, sequence)
        if read_result.polyT != -1:
            read_result.set_strand("+")
        if read_result.is_valid():
            return read_result

        rev_seq = reverese_complement(sequence)
        logger.debug("===== REV =====")
        read_rev_result = self._find_barcode_umi_fwd(read_id, rev_seq)
        if read_rev_result.polyT != -1:
            read_rev_result.set_strand("-")
        if read_rev_result.is_valid():
            return read_rev_result
        logger.debug("===== DONE =====")
        return read_result if read_result.more_informative_than(read_rev_result) else read_rev_result

    @staticmethod
    def result_type():
        return TenXBarcodeDetectionResult



class TenXBarcodeExtractorV2(TenXBarcodeExtractor):
    def __init__(self):
        TenXBarcodeExtractor.__init__(self, TenXVersions.v2)


class TenXBarcodeExtractorV3(TenXBarcodeExtractor):
    def __init__(self):
        TenXBarcodeExtractor.__init__(self, TenXVersions.v3)
