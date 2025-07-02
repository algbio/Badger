import pandas as pd


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
    reads = pd.read_csv(in_tsv, sep="\t")
    ids = reads["#read_id"].tolist()
    observed = reads["barcode"]
    observed = observed.fillna('*')
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