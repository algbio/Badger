###########################################################################
# Copyright (c) 2025 University of Helsinki
# # All Rights Reserved
# See file LICENSE for details.
############################################################################

import logging
import math
from collections import defaultdict
from enum import Enum, unique

from .kmer_indexer import KmerIndexer
from .common import find_polyt, reverese_complement, detect_exact_positions

logger = logging.getLogger('BarcodeGraph')


@unique
class ElementType(Enum):
    PolyT = 1
    cDNA = 2
    CONST = 10
    VAR_ANY = 11
    VAR_LIST = 12
    VAR_FILE = 13


class MoleculeElement:
    def __init__(self, element_name: str, element_type: ElementType, element_value: str = ""):
        self.element_name = element_name
        self.element_type = element_type

        if self.element_type in [ElementType.PolyT, ElementType.cDNA]:
            self.element_value = None
            self.element_length = -1
        elif self.element_type == ElementType.CONST:
            self.element_value = element_value
            self.element_length = len(self.element_value)
        elif self.element_type == ElementType.VAR_FILE:
            self.element_value = element_value
            self.element_length = len(open(self.element_value, 'r').readline().strip())
        elif self.element_type == ElementType.VAR_LIST:
            self.element_value = element_value.split(',')
            self.element_length = len(element_value[0])
        elif self.element_type == ElementType.VAR_ANY:
            self.element_value = None
            self.element_length = int(element_value)
        else:
            assert False

    def is_variable(self):
        return self.element_type in {ElementType.VAR_ANY, ElementType.VAR_LIST, ElementType.VAR_FILE}


class MoleculeStructure:
    def __init__(self, str_iterator):
        self.ordered_elements: list[MoleculeElement] = []

        l = str_iterator.next()
        elements = list(map(lambda x: x.strip(), l.strip().split(':')))
        element_properties = {}
        for l in str_iterator:
            v = l.strip().split('\t')
            assert len(v) == 3
            element_properties[v[0]] = (v[1], v[2])

        for el in elements:
            if el not in element_properties:
                if el not in ElementType.__dict__:
                    logger.critical("Molecule element %s was not described in the format file" % el)
                    exit(-1)
                element_type = ElementType[el]
                self.ordered_elements.append(MoleculeElement(el, element_type))
            else:
                element_type, element_val = element_properties[el]
                if element_type not in ElementType.__dict__:
                    logger.critical("Molecule element type %s is not among the possible types" % element_type)
                    exit(-1)
                element_type = ElementType[element_type]
                self.ordered_elements.append(MoleculeElement(el, element_type, element_val))

    def __iter__(self):
        for e in self.ordered_elements:
            yield e


class DetectedElement:
    def __init__(self, start, end, score):
        self.start = start
        self.end = end
        self.score = score


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

    def find_patterns(self, read_id, sequence):
        read_result = self._find_patterns_fwd(read_id, sequence)
        if read_result.polyT != -1:
            read_result.set_strand("+")

        rev_seq = reverese_complement(sequence)
        read_rev_result = self._find_patterns_fwd(read_id, rev_seq)
        if read_rev_result.polyT != -1:
            read_rev_result.set_strand("-")

        if read_rev_result.is_valid() and read_result.is_valid():
            return read_result if read_result.more_informative_than(read_rev_result) else read_rev_result
        if read_rev_result.is_valid():
            return read_rev_result
        return read_result

    def _find_patterns_fwd(self, read_id, sequence):
        logger.debug("== read id %s ==" % read_id)
        detected_elements = {}
        polyt_start, polyt_end = -1, -1
        if self.has_polyt:
            polyt_start, polyt_end = find_polyt(sequence)
            detected_elements[ElementType.PolyT] = DetectedElement(polyt_start, polyt_end, 0)
        logger.debug("PolyT %d" % polyt_start)

        self.detect_const_elements(sequence, polyt_start, polyt_end, detected_elements)
        self.extract_variable_elements5(detected_elements)

    def extract_variable_elements5(self, detected_elements):
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
            detected_elements[el.element_name] = DetectedElement(potential_start, potential_end, 0)

        current_pos = detected_elements[self.molecule_structure.ordered_elements[first_detected_const_element].element_name].end + 1
        for i in range(first_detected_const_element + 1, len(self.molecule_structure.ordered_elements)):
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
            next_el = self.molecule_structure.ordered_elements[i + 1]
            if next_el.element_name in detected_elements:
                potential_end = detected_elements[next_el].start - 1
            else:
                potential_end = potential_start + el.element_length - 1
            potential_len = potential_end - potential_start + 1
            if abs(potential_len - el.element_length) <= self.MAX_LEN_DIFF * el.element_length:
                detected_elements[el.element_name] = DetectedElement(potential_start, potential_end, 0)

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
            if element_start:
                current_search_start = element_end
                if not first_element_detected:
                    first_element_detected = True

            detected_elements[el.element_name] = DetectedElement(element_start, element_end, element_score)

        # TODO: run the code above for reversed molecule structure
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
            if element_start:
                current_search_end = element_end
                if not last_element_detected:
                    last_element_detected = True

            detected_elements[el.element_name] = DetectedElement(element_start, element_end, element_score)



