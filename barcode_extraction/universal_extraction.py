###########################################################################
# Copyright (c) 2025 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging

from .extraction_result import DetectedElement, ExtractionResult
from .kmer_indexer import KmerIndexer
from .common import find_polyt, reverese_complement, detect_exact_positions
from .molecule_structure import ElementType, MoleculeStructure

logger = logging.getLogger('Badger')


class UniversalSingleMoleculeExtractor:
    MIN_SCORE_COEFF = 0.75
    MIN_SCORE_COEFF_TERMMINAL = 0.5
    TERMINAL_MATCH_DELTA = 3
    MAX_LEN_DIFF = 0.2

    def __init__(self, molecule_structure: MoleculeStructure):
        self.molecule_structure = molecule_structure
        self.index_dict = {}
        self.has_polyt = False
        self.has_cdna = False
        self.elements_to_detect = set()
        self.elements_to_extract = set()

        for el in self.molecule_structure:
            if el.element_type == ElementType.CONST:
                self.index_dict[el.element_name] = KmerIndexer([el.element_value], max(6, int(el.element_length / 2) - 2))
            elif el.element_type == ElementType.PolyT:
                if self.has_polyt:
                    logger.error("Current version only supports a single polyT, extraction results may be suboptimal")
                self.has_polyt = True
            elif el.element_type == ElementType.cDNA:
                if self.has_cdna:
                    logger.error("Current version only supports a single cDNA, reads will not be split into molecules")
                self.has_cdna = True
            elif el.element_type in {ElementType.VAR_ANY, ElementType.VAR_FILE, ElementType.VAR_LIST}:
                self.elements_to_extract.add(el.element_name)
                if el.element_type in [ElementType.VAR_FILE, ElementType.VAR_LIST]:
                    self.elements_to_detect.add(el.element_name)

        if not self.has_cdna:
            logger.error("Molecule must include a cDNA")
            exit(-2)

    def header(self):
        return self.molecule_structure.header()

    def format_result(self, result: ExtractionResult):
        res_str = "%s\t%s" % (result.read_id, result.stand)
        for el in self.molecule_structure:
            if el.element_type == ElementType.cDNA: continue
            if el.element_type == ElementType.PolyT:
                if ElementType.PolyT.name not in result.detected_results:
                    res_str += "\t-1\t-1"
                    continue
                detected_element = result.detected_results[ElementType.PolyT.name]
                res_str += "\t%d\t%d" % (detected_element.start, detected_element.end)
            elif el.element_type == ElementType.CONST:
                if el.element_name not in result.detected_results:
                    res_str += "\t-1\t-1\t0"
                    continue
                detected_element = result.detected_results[el.element_name]
                res_str += "\t%d\t%d\t%d" % (detected_element.start, detected_element.end, detected_element.score)
            else:
                if el.element_name not in result.detected_results:
                    res_str += "\t-1\t-1\t*"
                    continue
                detected_element = result.detected_results[el.element_name]
                res_str += "\t%d\t%d\t%s" % (detected_element.start, detected_element.end, detected_element.seq)
        return res_str

    def find_patterns(self, read_id, sequence):
        detected_elements_fwd = self._find_patterns_fwd(read_id, sequence)
        rev_seq = reverese_complement(sequence)
        detected_elements_rev = self._find_patterns_fwd(read_id, rev_seq)

        if self.has_polyt:
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
        # TODO: make a list of DetectedElement rather that dict of str->DetectedElement for speed
        detected_elements = {}
        polyt_start, polyt_end = -1, -1
        if self.has_polyt:
            polyt_start, polyt_end = find_polyt(sequence)
            detected_elements[ElementType.PolyT.name] = DetectedElement(polyt_start, polyt_end, 0)
        logger.debug("PolyT %d" % polyt_start)

        self.detect_const_elements(sequence, polyt_start, polyt_end, detected_elements)
        self.extract_variable_elements5(detected_elements, sequence)
        self.extract_variable_elements3(detected_elements, sequence)
        logger.debug("== end read id %s ==" % read_id)
        return detected_elements

    def extract_variable_elements5(self, detected_elements, sequence):
        first_detected_const_element = None
        for i, el in enumerate(self.molecule_structure):
            if el.element_type == ElementType.cDNA:
                break
            if el.element_type in [ElementType.CONST, ElementType.PolyT] and el.element_name in detected_elements:
                first_detected_const_element = i
                break

        if first_detected_const_element is None: return

        # extracting elements preceding first detected const element
        current_pos = detected_elements[self.molecule_structure.ordered_elements[first_detected_const_element].element_name].start - 1
        for i in range(first_detected_const_element - 1, -1, -1):
            el = self.molecule_structure.ordered_elements[i]
            if not el.is_variable(): break
            potential_end = current_pos
            potential_start = potential_end - el.element_length + 1
            if potential_start < 0: break
            detected_elements[el.element_name] = DetectedElement(potential_start, potential_end, 0,
                                                                 seq=sequence[potential_start:potential_end+1])
            current_pos = potential_start - 1

        current_pos = detected_elements[self.molecule_structure.ordered_elements[first_detected_const_element].element_name].end + 1
        for i in range(first_detected_const_element + 1, len(self.molecule_structure.ordered_elements)):
            if current_pos >= len(sequence):
                break

            el = self.molecule_structure.ordered_elements[i]
            if el.element_type in [ElementType.cDNA, ElementType.PolyT]: break
            elif el.element_type == ElementType.CONST:
                if el.element_name in detected_elements:
                    current_pos = detected_elements[el.element_name].end + 1
                else:
                    current_pos += el.element_length
                continue

            potential_start = current_pos
            if i + 1 == len(self.molecule_structure.ordered_elements):
                potential_end = potential_start + el.element_length - 1
            else:
                next_el = self.molecule_structure.ordered_elements[i + 1]
                if next_el.element_name in detected_elements:
                    potential_end = detected_elements[next_el.element_name].start - 1
                else:
                    potential_end = potential_start + el.element_length - 1
            if potential_end >= len(sequence):
                potential_end = len(sequence) + 1

            potential_len = potential_end - potential_start + 1
            if abs(potential_len - el.element_length) <= self.MAX_LEN_DIFF * el.element_length:
                detected_elements[el.element_name] = DetectedElement(potential_start, potential_end, 0,
                                                                     seq=sequence[potential_start:potential_end+1])
            current_pos = potential_end + 1

    def extract_variable_elements3(self, detected_elements, sequence):
        last_detected_const_element = None
        for i in range(len(self.molecule_structure.ordered_elements) - 1, -1, -1):
            el = self.molecule_structure.ordered_elements[i]
            if el.element_type == ElementType.cDNA:
                break
            if el.element_type in [ElementType.CONST, ElementType.PolyT] and el.element_name in detected_elements:
                last_detected_const_element = i
                break

        if last_detected_const_element is None: return

        # extracting elements following last detected const element
        current_pos = detected_elements[self.molecule_structure.ordered_elements[last_detected_const_element].element_name].end + 1
        for i in range(last_detected_const_element + 1, len(self.molecule_structure.ordered_elements)):
            el = self.molecule_structure.ordered_elements[i]
            if not el.is_variable(): break
            potential_start = current_pos
            potential_end = potential_start + el.element_length - 1
            if potential_end >= len(sequence): break
            detected_elements[el.element_name] = DetectedElement(potential_start, potential_end, 0,
                                                                 seq=sequence[potential_start:potential_end+1])
            current_pos = potential_end + 1

        current_pos = detected_elements[self.molecule_structure.ordered_elements[last_detected_const_element].element_name].start - 1
        for i in range(last_detected_const_element - 1, -1, -1):
            if current_pos <= 0:
                break
            el = self.molecule_structure.ordered_elements[i]
            if el.element_type in [ElementType.cDNA, ElementType.PolyT]: break
            elif el.element_type == ElementType.CONST:
                if el.element_name in detected_elements:
                    current_pos = detected_elements[el.element_name].start - 1
                else:
                    current_pos -= el.element_length
                continue

            potential_end = current_pos
            if i == 0:
                potential_start = potential_end - el.element_length + 1
            else:
                prev_el = self.molecule_structure.ordered_elements[i - 1]
                if prev_el.element_name in detected_elements:
                    potential_start = detected_elements[prev_el].end + 1
                else:
                    potential_start = potential_end - el.element_length + 1
            if potential_start < 0:
                potential_start = 0

            potential_len = potential_end - potential_start + 1
            if abs(potential_len - el.element_length) <= self.MAX_LEN_DIFF * el.element_length:
                detected_elements[el.element_name] = DetectedElement(potential_start, potential_end, 0,
                                                                     seq=sequence[potential_start:potential_end+1])
            current_pos = potential_start - 1

    def detect_const_elements(self, sequence, polyt_start, polyt_end, detected_elements):
        # searching left of cDNA
        first_element_detected = False
        current_search_start = 0
        for i, el in enumerate(self.molecule_structure):
            if el.element_type in [ElementType.cDNA, ElementType.PolyT]:
                break

            if el.element_type != ElementType.CONST:
                continue

            search_start = current_search_start
            search_end = len(sequence) if polyt_start == -1 else polyt_start + 1
            element_occurrences = self.index_dict[el.element_name].get_occurrences_substr(sequence, search_start, search_end)
            if not first_element_detected:
                min_score = int(el.element_length * self.MIN_SCORE_COEFF_TERMMINAL)
                start_delta = -1
            else:
                min_score = int(el.element_length * self.MIN_SCORE_COEFF)
                start_delta = self.TERMINAL_MATCH_DELTA

            element_start, element_end, element_score = detect_exact_positions(sequence,
                                                                               search_start, search_end,
                                                                               self.index_dict[el.element_name].k,
                                                                               el.element_value,
                                                                               element_occurrences,
                                                                               min_score=min_score,
                                                                               start_delta=start_delta,
                                                                               end_delta=self.TERMINAL_MATCH_DELTA)

            if element_start is not None:
                current_search_start = element_end
                if not first_element_detected:
                    first_element_detected = True

                detected_elements[el.element_name] = DetectedElement(element_start, element_end, element_score)

        last_element_detected = False
        current_search_end = len(sequence)
        for i in range(len(self.molecule_structure.ordered_elements) - 1, -1, -1):
            el = self.molecule_structure.ordered_elements[i]
            if el.element_type in [ElementType.cDNA, ElementType.PolyT]:
                break

            if el.element_type != ElementType.CONST:
                continue

            search_start = 0 if polyt_start == -1 else polyt_end + 1
            search_end = current_search_end
            element_occurrences = self.index_dict[el.element_name].get_occurrences_substr(sequence, search_start, search_end)

            if not last_element_detected:
                min_score = int(el.element_length * self.MIN_SCORE_COEFF_TERMMINAL)
                end_delta = -1
            else:
                min_score = int(el.element_length * self.MIN_SCORE_COEFF)
                end_delta = self.TERMINAL_MATCH_DELTA

            element_start, element_end, element_score = detect_exact_positions(sequence,
                                                                               search_start, search_end,
                                                                               self.index_dict[el.element_name].k,
                                                                               el.element_value,
                                                                               element_occurrences,
                                                                               min_score=min_score,
                                                                               start_delta=self.TERMINAL_MATCH_DELTA,
                                                                               end_delta=end_delta)
            if element_start is not None:
                current_search_end = element_end
                if not last_element_detected:
                    last_element_detected = True

                detected_elements[el.element_name] = DetectedElement(element_start, element_end, element_score)



