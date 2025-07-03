import pandas as pd
from barcode_extraction.extraction_result import ExtractionResult, DetectedElement
from collections import defaultdict


def load_true_barcodes(inf):
    true_barcodes = inf
    if true_barcodes:
        true_barcodes = pd.read_csv(true_barcodes, sep = "\t", header = None)
        true_barcodes = true_barcodes.iloc[:,0].tolist()
        if true_barcodes[0][-1] == '1':
            for i in range(len(true_barcodes)):
                true_barcodes[i] = true_barcodes[i][:-2]
        true_barcodes = set(true_barcodes)
    return true_barcodes


def load_extracted_barcodes(in_tsv):
    reads = pd.read_csv(in_tsv, sep="\t", header=True)
    elements = defaultdict(list)
    for c in reads.columns:
        v = c.split('_')
        element_name = v[0]
        property_name = v[1]
        elements[element_name].append(property_name)

    read_assignments = []
    for r in reads:
        # if r["Barcode_sequence"] == '*':
        #    continue
        er = ExtractionResult(r["#read_id"], r["strand"], {})
        for e in elements.keys():
            de = DetectedElement(r[e + "_start"], r[e + "_end"])
            if e + "_sequence" in r:
                de.seq = r[e + "_sequence"]
            if e + "_score" in r:
                de.seq = r[e + "_score"]
            er.detected_results[e] = de
        read_assignments.append(er)
    return read_assignments
