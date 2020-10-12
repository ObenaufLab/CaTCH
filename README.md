# CaTCH 0.7.2

Catch is available as a command-line tool for all, and as a VBC Galaxy tool for the members of the Vienna Biocenter.

## Input

- One or multiple BAM files containing the CaTCH reads. Each BAM file can be multiplexed or a single sample.

- In the case of multiplexed BAMs, each BAM file must be demultiplexed/counted *individually*. 
A tab-delimited file containing a table matching the multiplexing tags to sample names must be provided for each BAM file.
The table should contain 2 named columns: "Barcode" and "Sample", for the sample tag and sample name respectively.
Demultiplexing tags can be a mix of different lengths, but the short ones must *not* be substrings of the longer ones.

```
Barcode Sample
TGAG    reference_1
TCGG    reference_2
ATCACG  notreatment_1
ATGGCG  notreatment_2
ATCTAG  treatmentA_1
CGAT    treatmentA_2
```

- For colour-coding figures in the report, a file listing the samples to be included in the report (in the desired order), as well as the respective treatment/condition and the desired display colour. Unlike the demultiplexing table, the conditions table is *not optional*. 
It must consist of 3 named columns: "Sample", "Treatment", "Colour". Colour must be in an R-compatible string format, such as these: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf .

```
Sample  Treatment   Colour
reference_1 reference   black
reference_2 reference   black
notreatment_1   notreatment blue
notreatment_2   notreatment blue
treatmentA_1    treatmentA  red
treatmentA_2    treatmentA  red
```

For convenience, the above two tables can also be combined into a single table used both for demultiplexing and compiling the report. Both the quantifier and the report can recognize this format in place of theirs.

```
Barcode Sample  Treatment   Colour
TGAG    reference_1 reference   black
TCGG    reference_2 reference   black
ATCACG  notreatment_1   notreatment blue
ATGGCG  notreatment_2   notreatment blue
ATCTAG  treatmentA_1    treatmentA  red
CGAT    treatmentA_2    treatmentA  red
```


In all three cases, *the columns must appear in the specified order.*



## Workflow - VBC Galaxy

The VBC Galaxy instance is not public. These instructions are for internal institute use only. 
For general use, consult the next section: "Workflow - Command line".

### Step 1 : Quantify barcodes

* Each BAM file needs to be processed separately. 
* Multiplexed BAMs additionally *each* require a sample tags table, as per the Input section above. 
* A counts table will be generated for each sample and a summary table for each BAM.

The original Obenauf-group barcode design is hard-coded into the script and is currently the only 
available option for Galaxy. Both C.Umkehrer's and S.Cronin's variants can be used as single sample 
per BAM file. Umkehrer's variant can also be used as a multiplexed BAM with a demultiplexing table 
as per Input section above.

Some advanced options made available through Galaxy are discussed in the Command Line sections below.

### Step 2: Collect outputs

The individual sample outputs need to be combined: One collective table for the barcode counts, and one 
collective table for the summaries. Simply select the Step-1 Galaxy outputs to combine. Importing
external tables into Galaxy at this stage is *not well supported*. 

### Step 3: Visual report

This requires the two collective tables from step 2. Additionally it needs a conditions table, as per 
Input section above. The results can be obtained from the zipfile. The main output is an HTML report.
Tables of the correaltions, counts and proportions are also output. 
*The HTML does not display correctly within Galaxy*, so be sure to *save the zip* file and open the files on 
your computer instead.




## Workflow - Command line

This is applicable to anyone who obtains the source code.

### Requirements

The 3rd-party packages below are needed and must be available/loaded/activated in the working environment before running the CaTCH scripts.
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


### Step 1 : Quantify barcodes

**barcodingQuantifier.py**

Starting from a BAM file and an optional table with demultiplexing info, it uses regex matching 
to extract the barcode sequence. It creates:

* a barcode counts file per sample, and 
* an extraction summary file per BAM.

Additionally, for troubleshooting use, it outputs:

* a BAM with reads where a barcode couldn't be identified (truncation, or mutation in inflexible bases),
* a BAM file for reads whose sample tags didn't match the demultiplex info, 
* a text file with the tally of the sample tags of the reads that couldn't be demultiplexed, and 
* a BAM file with reads that matched the empty or spike-in patterns in the Obenauf-group design.


#### Parameters:

```
-f FILE         A single BAM file.
-b FILE         Demultiplexing tags table. (tab delimited: Barcode\tSample, with header line). Omit if data already demultiplexed.
-o DIR          Output directory for counts (./process/).
--revcomp       Reverse complement the sample tags (false).
--bc_len        Length of the barcode. Determines the size of the sequences extracted from the reads.
--gt_len        Number of bases between the start of the barcode and the end fo the sample tag. This is tested only for sample tags that precede the barcode
                and allows a spacer between barcode and sample tag, such as the genotype tag in the Obenauf-group design.
                Theoretically, negative values should work too and allow the sample tag to be located after the start of the barcode, but this is untested.
--spikein       Also search for the spike-in barcode (false). Its sequence is hard-coded.
--stringent     Use the full-length pattern (false). This is for the hard-coded Obenauf-group design only. When false, the shorter/core pattern is used.
--n_dark INT    Number of dark/unknown bases (N) to consider in the barcode pattern (0). 
                For this to have an effect, the provided --pattern must be designed to also match Ns.
                There is a hard-coded version of the short Obenauf-group design for this.
--pattern       Barcode pattern as a regex string. This overrides the hard-coded patterns to enable other designs.
                The regex must begin at the start of the barcode, but doesn't need to reach the end of the barcode.
                bc_len, gt_len, n_dark are available to use with a custom pattern.
```

### Step 2 : Merge barcodes within a certain edit distance (optional)

**barcodingHammingMerge.py**

The idea is to reduce the noise from sequencing errors. The merge is based on the Hamming distance 
between barcodes. That means that only substitutions are considered, no insertions or deletions. When
barcodes are merged, the sequence of the most abundant is the one that is kept.
It is up to your judgement to apply this step or not, as well as to choose the edit distance.

This step is very slow and resource hungry.

#### Parameters:

```
-b FILE       A tab delimited file with barcode counts, as produced by the previous step.
-d INT        Hamming distance at which barcodes should be considered the same (Default: 1)
-o FILE       Output file.
```

### Step 3 : Collect outputs

**mergeBCcounts.R**

Full-join/merge the tables from Step 1 (or from Step 2 instead, if applicable).
You can provide as many input tables as you want.

Please make one collective table for the barcode counts *and* one for the summaries. Both are needed 
for the report. That means running this script twice, once for all the barcode count files, and once for 
all the summary files.

#### Usage:

```
mergeBCcounts.R OUTPUT_FILE INPUT_FILE INPUT_FILE ...
```

### Step 4 : Visual report

**barcoding_results_run.R**

This will compile a report using the **barcoding_results_template.Rmd** template. It creates:

* an HTML report with interactive figures,
* optionally a PDF static copy of those figures,
* a table with the counts and percentages of barcodes and their assigned IDs, and
* a table with the pairwise sample correlations.

**IMPORTANT**
Please provide *full paths* to the files and directories.

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


### Steps 1-4 : Pipeline

**catch_workflow.sh**

This is my Bash script to automate the steps of the pipeline. It requires a SLURM environment. 
Currenlty it only runs the original Obenauf group design.

**I cannot promise it will work for you**, but at the very least it can be used as template to write 
a new pipeline. Alternatively, the workflow can be executed directly step by step, as described in the previous sections.

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


