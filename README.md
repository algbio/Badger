# Badger manual

1. [About Badger](#sec1) </br>
1.1. [Supported data types](#sec1.1)</br>
2. [Installation](#sec2)</br>
3. [Running Badger](#sec3)</br>
3.1. [Badger input](#sec3.1)</br>
3.2. [Command line options](#sec3.2)</br>
3.3. [Badger output](#sec3.3)</br>

<a name="sec1"></a>
# About Badger

Badger is a tool for long read barcode calling. For a given set of single cell
long reads, it extracts the barcodes, identifies cell-associated barcodes and
corrects extracted barcodes containing errors. The correction stage is based on
an edit distance graph.

Badger works in two stages. First, barcodes are extracted, which results in a file
containing readIDs, barcodes and further information about the reads. This is then
the input for the next step, in which the graph is constructed and the barcodes
corrected.

Currently supported protocols are 10x single cell and 10x visium.

<a name="sec1.1"></a>
## Supported data types

Badger support all kinds of long single cell RNA data:
* PacBio CCS
* ONT dRNA / ONT cDNA
* Assembled / corrected transcript sequences

Reads must be provided in FASTQ or FASTA format (can be gzipped) or
aligned reads as a BAM or SAM file.

<a name="sec2"></a>
# Installation

To obtain Badger you can download repository and install requirements.
Clone Badger repository and switch to the latest release:
```bash
git clone https://github.com/algbio/Badger.git
cd Badger
git checkout latest
```


<a name="sec3"></a>
# Running Badger
<a name="sec3.1"></a>
## Badger input
To run Badger, you should provide:
* Long single cell RNA reads (PacBio or Oxford Nanopore) in one of the following formats:
  * FASTA/FASTQ (can be gzipped);
  * Sorted and indexed BAM;
* Optionally a list of cell-associated barcodes
* Barcode whitelist for the used single cell sequencing protocol

<a name="sec3.2"></a>
## Badger command line options

### Barcode extraction step

#### Basic options

`--output` (or `-o`)
    Prefix for output files in respect to the current folder

`--help` (or `-h`)
    Prints help message.

#### Input options

`--barcodes` (or `-b`)
    Barcode whitelist for the used protocol

`--input` (or `-i`)
    Reads in FASTA, FASTQ, BAM or SAM format

`--mode`
    Extraction method to be used, currently only tenX

#### Algorithm parameters

`--threads` (or `-t`)
    Number of threads to use, default 16

### Barcode correction step

#### Basic options
`--output` (or `-o`)
    Prefix for output files in respect to the current folder

`--help` (or `-h`)
    Prints help message.

#### Input options

`--barcodes` (or `-b`)
    Barcodes extracted from the long reads in tsv format

`--reads` (or `-r`)
    This is the output file of the extraction step, used getting readIDs

`--barcode_list` (or `-l`)
    Barcode whitelist for the used protocol

`--true_barcodes`
    List of the cell-associated barcodes, optional

`--data_type` (or `-d`)
    Type of data to process, supported values are:  `10x` and `visium`

#### Algorithm parameters
<a name="params"></a>

`--threshold` (or `-t`)
    Maximal edit distance between barcodes to be connected in the graph, default 1

`--n_cells` (or `-c`)
    Expected number of cell-associated barcodes

### Examples
<a name="examples"></a>

* Extracting 10x single cell barcodes from reads

```bash
detect_barcodes.py --barcodes whitelist.txt --input scRNAseq_reads.fasta
  --mode tenX --output barcode_file
```

* Correcting extracted 10x single cell barcodes

```bash
barcodes.py --barcodes barcodes.tsv --reads barcode_file.tsv --barcode_list whitelist.txt
  -d 10x --output corrected_barcodes --n_cells 5000
```

<a name="sec3.3"></a>
## Badger output

### Output files

Both extraction and correction step have one output file. They will be shortly described here.

#### Barcode extraction output
OUTPUT_PREFIX.tsv - TSV file containing readID, barcodes and additional information

Columns are:
  * `#read_id` - readID
  * `barcode` - extracted barcode for the read, `*` if no barcode could be extracted
  * `BC_score` - 0 if a barcode could be extracted, -1 if not
  * `valid_UMI` - TRUE if a valid UMI was found, FALSE if not
  * `strand` - read direction in which a barcode was found as `+` or `-`, `.` if no barcode could be extracted
  * `polyT_start` - position in the read where the polyT sequence starts, -1 if no polyT sequence was found
  * `R1_end` - position in the read where the adapter ends, -1 if no adapter was found

#### Barcode correction output

OUTPUT_PREFIX_output_file.tsv - TSV file containing readID and assigned cell-associated barcode

Columns are:
  * `readID` - readID
  * `barcode` -  assigned barcode, `*` if no barcode could be assigned
