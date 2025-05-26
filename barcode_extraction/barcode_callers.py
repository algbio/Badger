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
from .extraction_result import ExtractionResult, DetectedElement
from .molecule_structure import MoleculeStructure, MoleculeElement, ElementType
logger = logging.getLogger('BarcodeGraph')


class ReadStats:
    def __init__(self):
        self.read_count = 0
        self.pattern_counts = defaultdict(int)

    def add_read(self, barcode_detection_result):
        self.read_count += 1
        for el in barcode_detection_result.detected_results:
            if barcode_detection_result.detected_results[el].start != -1:
                self.pattern_counts[el] += 1

    def __str__(self):
        human_readable_str =  ("Total reads:\t%d\n" % self.read_count)
        for a in self.pattern_counts:
            human_readable_str += "%s:\t%d\n" % (a, self.pattern_counts[a])
        return human_readable_str


@unique
class TenXMoleculeElements(Enum):
     r1 = 1
     barcode = 2
     umi = 3
     polyT = 4
     tso = 5


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
        self.molecule_structure = (
            MoleculeStructure.from_element_list([MoleculeElement(TenXMoleculeElements.r1.name, ElementType.CONST),
                                                 MoleculeElement(TenXMoleculeElements.barcode.name, ElementType.VAR_ANY, self.BARCODE_LEN_10X),
                                                 MoleculeElement(TenXMoleculeElements.umi.name, ElementType.VAR_ANY, self.UMI_LEN_10X),
                                                 MoleculeElement(TenXMoleculeElements.polyT.name, ElementType.PolyT)]))

    def header(self):
        return self.molecule_structure.header()

    def format_result(self, result: ExtractionResult):
        res_str = "%s\t%s" % (result.read_id, result.stand)
        for el in self.molecule_structure:
            if el == ElementType.PolyT:
                detected_element = result.detected_results[ElementType.PolyT.name]
                res_str += "\t%d\t%d" % (detected_element.start,detected_element.end)
            elif el.element_type == ElementType.CONST:
                detected_element = result.detected_results[el.element_name]
                res_str += "\t%d\t%d\t%d" % (detected_element.start,detected_element.end, detected_element.score)
            else:
                detected_element = result.detected_results[el.element_name]
                res_str += "\t%d\t%d\t%s" % (detected_element.start, detected_element.end, detected_element.seq)

    def find_patterns(self, read_id, sequence):
        detected_elements_fwd = self._find_patterns_fwd(read_id, sequence)
        rev_seq = reverese_complement(sequence)
        detected_elements_rev = self._find_patterns_fwd(read_id, rev_seq)

        fwd_has_polyt = ElementType.PolyT.name in detected_elements_fwd
        rev_has_polyt = ElementType.PolyT.name in detected_elements_rev
        if fwd_has_polyt and not rev_has_polyt:
            return ExtractionResult(read_id, '+', detected_elements_fwd)
        if not fwd_has_polyt and rev_has_polyt:
            return ExtractionResult(read_id, '-', detected_elements_rev)

        if len(detected_elements_fwd) >= len(detected_elements_rev):
            return ExtractionResult(read_id, '+', detected_elements_fwd)
        return ExtractionResult(read_id, '-', detected_elements_rev)

    def _find_patterns_fwd(self, read_id, sequence):
        logger.debug("== read id %s ==" % read_id)
        detected_elements = {}
        polyt_start = find_polyt_start(sequence)
        logger.debug("PolyT %d" % polyt_start)
        detected_elements[TenXMoleculeElements.polyT.name] = DetectedElement(polyt_start, score=0)
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
            return detected_elements
        logger.debug("LINKER: %d-%d" % (r1_start, r1_end))
        detected_elements[TenXMoleculeElements.r1.name] = DetectedElement(r1_start, r1_end, r1_score)

        if polyt_start != -1 and polyt_start - r1_end < self.BARCODE_LEN_10X:
            return detected_elements

        if polyt_start == -1 or polyt_start - r1_end > self.BARCODE_LEN_10X + self.UMI_LEN_10X + 10:
            # if polyT was not detected earlier, use relaxed parameters once the linker is found
            presumable_polyt_start = r1_end + self.BARCODE_LEN_10X + self.UMI_LEN_10X
            search_start = presumable_polyt_start - 4
            search_end = min(len(sequence), presumable_polyt_start + 10)
            polyt_start = find_polyt_start(sequence[search_start:search_end], window_size=5, polya_fraction=1.0)
            if polyt_start != -1:
                polyt_start += search_start
                detected_elements[TenXMoleculeElements.polyT.name] = DetectedElement(polyt_start, score=0)

        barcode_start = r1_end + 1
        barcode_end = r1_end + self.BARCODE_LEN_10X
        potential_barcode = sequence[barcode_start:barcode_end + 1]
        logger.debug("Barcode: %s" % (potential_barcode))
        detected_elements[TenXMoleculeElements.barcode.name] = DetectedElement(barcode_start, barcode_end, seq=potential_barcode)

        potential_umi_start = barcode_end + 1
        potential_umi_end = polyt_start - 1
        # FIXME: UMI length
        if potential_umi_end - potential_umi_start <= 5:
            potential_umi_end = potential_umi_start + self.UMI_LEN_10X - 1
        potential_umi = sequence[potential_umi_start:potential_umi_end + 1]
        detected_elements[TenXMoleculeElements.umi.name] = DetectedElement(potential_umi_start, potential_umi_end,
                                                                               seq=potential_umi)

        return detected_elements

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


class TenXBarcodeExtractorV2(TenXBarcodeExtractor):
    def __init__(self):
        TenXBarcodeExtractor.__init__(self, TenXVersions.v2)


class TenXBarcodeExtractorV3(TenXBarcodeExtractor):
    def __init__(self):
        TenXBarcodeExtractor.__init__(self, TenXVersions.v3)
