# CaTCH

These are the analysis scripts for the CaTCH experiments.
They should include everything needed to go from an unaligned .BAM file to barcode counts and an HTML report with various plots, some of which are interactive.

These instructions apply mainly to command-line use. However, the input formats apply equally to the Galaxy instance.


## Requirements


The 3rd-party packages below are needed and must be available/loaded/activated in the working environment before running any of the CaTCH scripts.
The software versions are provided for reference, most likely many other combinations of versions should work too.

- Python 3 (3.6.7):
    - argparse
    - pysam
    - Biopython
    - pandas
    - Levenshtein

- R (3.5.1):
    - tidyverse
    - data.table
    - plotly
    - ggrepel
    - RColorBrewer
    - patchwork

- The CaTCH python and R scripts provided here.



## Input


- One or multiple BAM files containing the CaTCH reads. Each BAM file can be multiplexed or a single sample. No alignment is needed.

- In the case of multiplexed BAMs, a tab-delimited file containing a table matching the multiplexing tags to sample names. The table should contain 2 named columns: "barcode" and "sample", for the sample tag and sample name respectively. Demultiplexing tags can be a mix of different lengths, but the shorter ones must not be substrings of the longer ones.
  ie.

```
barcode sample
TGAG    reference_1
TCGG    reference_2
ATCACG  notreatment_1
ATGGCG  notreatment_2
ATCTAG  treatmentA_1
CGAT    treatmentA_2
```

- A file containing the condition allocation for each sample to be included in the report and the assigned colour for the condition. It must consist of 3 named columns: "sample", "condition", "colour". Colour must be in an R-compatible string format, one colour per condition.
  ie.

```
sample  condition   colour
reference_1 reference   black
reference_2 reference   black
notreatment_1   notreatment blue
notreatment_2   notreatment blue
treatmentA_1    treatmentA  red
treatmentA_2    treatmentA  red
```




## Workflow

### Steps 1-4 : Pipeline


The repository contains a shell script (catch_workflow.sh) that automates execution of all the steps, in the presence of a SLURM HPC-cluster system.
If such is not available, follow the individual steps instead. Steps 1 and 2 do use quite a bit of memory, they are probably not laptop-friendly.

#### Parameters:

```
-b DIR          Directory with the BAM file(s) to be analysed.
-c FILE         Demultiplexing table, if applicable.
-o DIR          Output directory for the counts files.
-O DIR          Ouptut directory for the analysis report.
-v FILE         List of samples in the desired order, with corresponding condition and corresponding display colour.
-X DIR          Where to find all the scripts for this workflow (probably your local clone of the repository).
-m INT          Hamming distance at which to merge barcodes as likely sequencing errors, if applicable (0 ie. not applicable).
-n INT          Number of dark bases to allow in the pattern (0).
-A FLOAT        Barcode abundance threshold for the report (0.01).
-R INT_LIST     Comma-seperated list of rows in covars to be used as reference samples in the report, NOT counting the header line (1).
-r              Reverse complement the sample tags (False).
-i              Spike-in barcode was added (False). The spike sequence is hard-coded.
-s              Match the full pattern of semi-random barcodes instead of the short pattern (False).
```

### Step 1 : Count the barcodes

barcodingQuantifier.py


Two barcode designs are programmed into the script:

```
P7adapter_template_BARCODE_GenotypeTag_SampleTag_NNNNNN_P5adapter
                                                    <--- read direction
```

In this design the sample tag is included in read1 of the mate pair, and the samples can be demultiplexed by the counting script, if a sample tag table is provided (see input section above). The sample tags are expected to be found at a fixed interval from where the barcode is located.

```
P7adapter_SampleTag_template_BARCODE_GenotypeTag_P5adapter
                                             <--- read direction
```

In this design the sample tag is unlikely to be included in read1 containing the barcode, and therefore must be demultiplexed IN ADVANCE by other means, not provided here.


The genotype tag is not used in either case in the current implementation, except as spacer between the barcode and the sample tag in the first design.


Two versions of the barcode-matching pattern are programmed in, one more stringent than the other:

- stringent:

```
[TC]{1}[ATGC]{2}[AT]{1}[TG]{1}[TGC]{1}[AT]{1}[TGC]{1}[ATGC]{1}[TGC]{1}[ATGC]{1}[ATC]{1}[TGC]{1}[ATGC]{3}[CG]{1}[ATGC]{1}[CGA]{1}CGCCG[TC]{1}[AGTC]{2}[AT]{1}[TG]{1}[TCG]{1}[AT]{1}GCC[ATGC]{1}[ACT]{1}[ATGC]{4}[GC]{1}[AGTC]{1}[CGA]{1}CGCCG
```

- non-stringent:

```
CGCCG[ATGC]{7}GCC[ATGC]{9}CGCCG
```

The stringent pattern matches from the beginning of the barcode to the end. The non-stringent matches the less variable sequence in the middle of the barcode, and the start of the barcode is found with a fixed offset of 19nt before the pattern.

There are currently no parameters to specify a different barcode pattern or different offset of the pattern from the start of the barcode.

#### Parameters:

```
-f FILE           Demultiplexing tags table. (tab delimited: Barcode\tSample, with header line). Omit if data already demultiplexed.
-o DIR            Output directory for counts (./process/).
-r                Reverse complement the sample tags (default: false).
-s                Stringent barcode matching (default: false).
--n_dark INT      Number of dark bases to consider in the barcode pattern (0). This allows up to that number of undefined N bases within the barcode.
--spikein         The spike-in barcode was included (default: false).
```

### Step 2 : Merge barcodes within a certain edit distance (optional)

barcodingHammingMerge.py


The merge is based on the Hamming distance between barcodes. That means that only substitutions are considered, no insertions or deletions.
This step is very slow.
It is up to your judgement to apply this step or not, as well as to choose the edit distance.

#### Parameters:

```
-b FILE       A tab delimited file with barcode counts, as produced by the previous step.
-d INT        Hamming distance at which barcodes should be considered the same (Default: 1)
-o FILE       Output file.
```

### Step 3 : Collect all the counts into one table

mergeBCcounts.R


This is simply a matter of full-merging the tables from Step 1 (or from Step 2 instead, if you chose to apply it).
Please make one collective table for the barcode counts and one for the summaries. Both are needed for the next step.


#### Parameters (positional):

```
1. FILE               Output file.
2. FILE FILE ...      Input files.
```


### Step 4 : Compile the report

barcoding_results_run.R


This will compile a report using the barcoding_results_template.Rmd template.

#### Parameters:

```
-c FILE       The collective counts table from step 3.
-s FILE       The collective summaries table from step 3.
-d DIR        Directory in which to save all output files.
-T FILE       The path to barcoding_results_template.Rmd .
-v FILE       The file assigning samples to conditions and colours (see Input section above).
-p            Instead of an HTML report (more information and more interaction), output just the figures into a PDF file (for Illustrator or presentations).
-r INT_LIST   Comma-separated list (without spaces) of integers designating the samples to use as reference abundance. The numbers should correspond to the order in which  the samples appear in the collective summaries file, not counting the header line) (Default: 1).
-N INT        Count threshold for barcodes (Default: 50). The counting generates many low-count barcodes, as a result of sequencing errors. These inflate the number of barcodes, so this threshold is provided to cut out that noise.
-A FLOAT      Proportional abundance threshold to consider barcodes top hits (Default: 0.01). Range 0 - 1.
-B BC_LIST    Comma-separated list of additional barcode IDs that should be included in the report despite not being among the top shared ones. These are the IDs assigned to the barcodes by the report, so you'll need to run the report at least once before.
```
