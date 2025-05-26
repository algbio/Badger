from enum import unique, Enum

from barcode_extraction.universal_extraction import logger


@unique
class ElementType(Enum):
    PolyT = 1
    cDNA = 2
    CONST = 10
    VAR_ANY = 11
    VAR_LIST = 12
    VAR_FILE = 13


class MoleculeElement:
    def __init__(self, element_name: str, element_type: ElementType, element_value = ""):
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

        l = next(str_iterator)
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

    @classmethod
    def from_element_list(cls, element_list):
        ms = cls.__new__(cls)
        ms.ordered_elements = element_list

    def __iter__(self):
        for e in self.ordered_elements:
            yield e

    def header(self):
        header = "#read_id\tstrand"
        for el in self.ordered_elements:
            if el.element_type == ElementType.cDNA: continue
            if el.element_type == ElementType.PolyT:
                header += "\tpolyT_start\tpolyT_end"
            elif el.element_type == ElementType.CONST:
                header += "\t%s_start\t%s_end\t%s_score" % (el.element_name, el.element_name, el.element_name)
            else:
                header += "\t%s_start\t%s_end\t%s_seqeunce" % (el.element_name, el.element_name, el.element_name)
        return header
