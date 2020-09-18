# CaTCH

Catch is available as a command-line tool and as a VBC Galaxy tool.

## Input

- One or multiple BAM files containing the CaTCH reads. Each BAM file can be multiplexed or a single sample. No alignment is needed.

- In the case of multiplexed BAMs, a tab-delimited file containing a table matching the multiplexing tags to sample names.
The table should contain 2 named columns: "barcode" and "sample", for the sample tag and sample name respectively.
Demultiplexing tags can be a mix of different lengths, but the shorter ones must not be substrings of the longer ones.

```
Barcode Sample
TGAG    reference_1
TCGG    reference_2
ATCACG  notreatment_1
ATGGCG  notreatment_2
ATCTAG  treatmentA_1
CGAT    treatmentA_2
```

- A file containing the condition allocation for each sample to be included in the report and the assigned colour for the condition.
It must consist of 3 named columns: "sample", "condition", "colour". Colour must be in an R-compatible string format, one colour per condition.

```
Sample  Treatment   Colour
reference_1 reference   black
reference_2 reference   black
notreatment_1   notreatment blue
notreatment_2   notreatment blue
treatmentA_1    treatmentA  red
treatmentA_2    treatmentA  red
```

Alternatively, the above two tables can be combined into a single table for both uses. This is convenient when the input BAM contains all the samples
to be included in the report.

```
Barcode Sample  Treatment   Colour
TGAG    reference_1 reference   black
TCGG    reference_2 reference   black
ATCACG  notreatment_1   notreatment blue
ATGGCG  notreatment_2   notreatment blue
ATCTAG  treatmentA_1    treatmentA  red
CGAT    treatmentA_2    treatmentA  red
```

**In all cases, the order of columns is strictly required!**



## Workflow - VBC Galaxy

The VBCF Galaxy instance is not public. These instructions are for internal institute use only.

### Step 1 : Quantify barcodes

Each BAM file needs to be processed separately. Multiplexed BAMs additionally each require a sample tags table, as per the Input section above.
A counts table will be generated for each sample and a summary table for each BAM.

The original Obenauf-group barcode design is hard-coded into the script and is currently the only available option for Galaxy.
Both C.Umkehrer's and S.Cronin's variants can be used as single sample per BAM file. Umkehrer's variant can also be used with a
demultiplexing table (as per Input section above).

Some advanced options made available through Galaxy are discussed in the Command Line sections below.

### Step 2: Collect outputs

The sample-wise outputs need to be merged together. One collective table for the barcode counts, and one collective table for the summaries.
There are many ways to join/merge tables using the barcode sequence as row key, in R, Python, Excel, etc, or you can use the catch count merger tool on our Galaxy.

### Step 3: Visual report

This requires the two collective tables from step 2. Additionally it needs a conditions table, as per Input section above.
The output is and HTML report.

Some advanced options made available through Galaxy are discussed in the Command Line sections below.



## Workflow - Command line

This is applicable to anyone who obtains the source code.

### Requirements

The 3rd-party packages below are needed and must be available/loaded/activated in the working environment before running any of the CaTCH scripts.
The versions are provided for reference, most likely many other combinations of versions should work too.

- Python (3.6.7):
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


### Steps 1-4 : Pipeline

*catch_workflow.sh*

This is my Bash script to automate the steps of the pipeline. It requires a SLURM environment. Currenlty it only runs the original Obenauf group design.
I cannot promise it will work for you, but at the very least it can be used as template to write a new pipeline.

Alternatively, the steps can be executed individually, as described in the next sections.

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


### Step 1 : Quantify barcodes

*barcodingQuantifier.py*

Starting from a BAM file and an optional table with demultiplexing info, it uses regex matching to extract the barcode sequence.
It creates a barcode counts file (per sample), a demultiplication and extraction summary file (per BAM).

Additionally, for troubleshooting use, it outputs a BAM with reads where a barcode couldn't be identified (truncation or mutation in inflexible bases),
a BAM file for reads whose sample tags didn't match the demultiplex info, a text file with the tally of the unmatched sample tags, and a
BAM file with reads that matched the empty or spike-in patterns in the Obenauf-group design.

#### Parameters:

```
-f FILE         A single BAM file/
-b FILE         Demultiplexing tags table. (tab delimited: Barcode\tSample, with header line). Omit if data already demultiplexed.
-o DIR          Output directory for counts (./process/).
--revcomp       Reverse complement the sample tags (false).
--bc_len        Length of the barcode.
--gt_len        Number of bases between the start of the barcode and the end fo the sample tag. This is tested only for sample tags that precede the barcode
                and allows a spacer between barcode and sample tag, such as the genotype tag in the Obenauf-group design.
                Theoretically, negative values should work too, and allow the sample tag to be located after the barcode, but this is untested.
--spikein       Also search for the spike-in barcode (false). Its sequence is hard-coded.
--stringent     Use the full-length pattern (false). This is for the hard-coded Obenauf-group design only. When false, the shorter core pattern is used.
--n_dark INT    Number of dark bases to consider in the barcode pattern (0). This requiresthe pattern to explicitly specify Ns.
                There is a hard-coded version of the shortObenauf-group design for this.
--pattern       Barcode pattern as a regex string. This overrides the hard-coded patterns to enable other designs.
                The regex must begin at the start of the barcode, but doesn't need to reach the end of the barcode.
                bc_len, gt_len, n_dark are available to use with a custom pattern.
```

### Step 2 : Merge barcodes within a certain edit distance (optional)

*barcodingHammingMerge.py*

The idea is to reduce the noise from sequencing errors.
The merge is based on the Hamming distance between barcodes. That means that only substitutions are considered, no insertions or deletions.
It is up to your judgement to apply this step or not, as well as to choose the edit distance.

This step is very slow and resource hungry.

#### Parameters:

```
-b FILE       A tab delimited file with barcode counts, as produced by the previous step.
-d INT        Hamming distance at which barcodes should be considered the same (Default: 1)
-o FILE       Output file.
```

### Step 3 : Collect outputs

*mergeBCcounts.R*

Full-join/merge the tables from Step 1 (or from Step 2 instead, if you chose to apply it).
This can be done in many different ways in R, Python, Excels, etc.
A script is provided for convenience and consistency. You can provide it with as many input tables as you want.

Please make one collective table for the barcode counts and one for the summaries. Both are needed for the report.
That means running this script twice, once for all the barcode counts, and once for all the summaries.

#### Usage:

```
mergeBCcounts.R OUTPUT_FILE INPUT_FILE INPUT_FILE ...
```

### Step 4 : Visual report

*barcoding_results_run.R*

This will compile a report using the *barcoding_results_template.Rmd* template.

*IMPORTANT*
Please provide full paths to the files and directories. On Linux systems, you can wrap the relative path in `$(realpath FILE)` to avoid typing it in.

#### Parameters:

```
-c FILE         The collective counts table from step 3. Full path needed.
-s FILE         The collective summaries table from step 3. Full path needed.
-d DIR          Directory in which to save all output files. Full path needed.
-T FILE         The path to barcoding_results_template.Rmd.
-v FILE         The file assigning samples to conditions and colours (see Input section above). Full path needed.
-p              In addition to the HTML report, output static figures into a PDF file for Illustrator or presentations (Default: false).
-r INT_LIST     Comma-separated list (without spaces) of integers designating the samples to use as reference abundance.
                The numbers should correspond to the order in which  the samples appear in the collective summaries file, not counting the header line) (Default: 1).
-N INT          Count threshold for barcodes (Default: 50). The counting generates many low-count barcodes, as a result of sequencing errors.
                These inflate the number of barcodes, so this threshold is provided to cut out that noise.
-A FLOAT        Proportional abundance threshold to consider barcodes top hits (Default: 0.01). Range 0 - 1.
-B BC_LIST      Comma-separated list of additional barcode IDs that should be included in the report despite not being among the top shared ones.
                These are the IDs assigned to the barcodes by the report, so you'll need to run the report at least once before you know which extra IDs to provide.
```
