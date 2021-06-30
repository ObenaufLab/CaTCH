# CaTCH 0.8.0.dev

CaTCH is a wetlab method for identifying cell clones from barcoded populations using CRISPRa-inducible reporters, developed in the Obenauf lab at the Institute for Molecular Pathology in Vienna.

*Isolating live cell clones from barcoded populations using CRISPRa-inducible reporters.*
**Nat Biotechnol. 2021 Feb;39(2):174-178.**

PMID: 32719478 
DOI: 10.1038/s41587-020-0614-0

This software package comprises the analysis steps for the sequencing data generated with the CaTCH method.
It is available as a command-line tool for everyone, and for the members of the Vienna Biocenter only it is also implemented in the VBC's Galaxy server.

This README file outlines the general workflow. For information on an individual step please consult the respective README file for that step.



## Input

### Sequencing data

The input data will be one or multiple BAM files containing the sequencing reads. 

The method has been tested with multiple BAMs containing one sample each, as well as with a single multiplexed BAM containing all the samples.
It is possible to use multiple multiplexed BAMs, but they may have to be quantified separately.

#### Two barcode designs have been tested:

Christian Umkehrer's original design allows the sample index to be contained in the first read. This design can be demultiplexed and quantified in a single step with the barcodingQuantifier.py script, so no additional step is required for pooled lanes.

*5'-P5adapter_NNNNNN_SampleINDEX_GenotypeTag_BARCODE_template_P7adapter-3'*

In an effort to move the barcode closer to the start of the read for higher quality, Shona Cronin's design moves the sample index to the end of the construct. It is unlikely that the first read will contain the sample index intact,
so pooled lanes with this design need to be demultiplexed in advance (outside of this package) using the second read of the pair.

*5'-P5adapter_GenotypeTag_BARCODE_template_SampleIndex_P7adapter-3'*

### Data description

You will need to compile a table with all the samples to be analyzed, containing the following columns and in that order:

1. `Tag` - The sample index. This should be a valid oligonucleotide (`[ATGC]`, no wildcards). A mix of Tag lengths is allowed, but the shorter ones should not be substrings of the longer ones, to prevent misidentifications. If mixing multiplexed and non-multiplexed BAMs, use a dummy value like `XXXX` for Samples that don't need demultiplexing.
2. `Sample` - A unique name for each sample.
3. `Group` - A coarse grouping of the samples. For example if you have multiple different treatments and/or mutliple controls, you can group them as treated/untreated.
4. `Treatment` - The specific treatment for that sample.
5. `Colour` - This is for display purposes in plots. You can use any colour name recognised by R. Most plain colour names (in lowercase) are recognised, though some are not as pleasant to look at. A suggested selection to choose from is: black, grey01 (dark) through grey99 (light), red, organge, darkgold, forestgreen, dodgerblue, steelblue, magenta, purple.

For columns 2, 3 and 4 the **values must be** plain ASCII alphanumeric strings, starting with a letter. Underscores are allowed. 
No accented or special letters, no spaces, no other symbols (including dashes). `^[A-Za-z][A-Za-z0-9_]*`

Column 1 is optional and only applies to quantifying multiplexed BAM files using Christian's design for the barcode construct.

```
Tag     Sample          Group       Treatment     Colour
TGAG    reference_1     untreated   reference     black
TCGG    reference_2     untreated   reference     black
ATCACG  notreatment_1   untreated   notreatment   blue
ATGGCG  notreatment_2   untreated   notreatment   blue
ATCTAG  treatmentA_1    treated     treatmentA    red
CGAT    treatmentA_2    treated     treatmentA    red
ATCCC   treatmentB_1    treated     treatmentB    orange
CAATT   treatmentB_2    treated     treatmentB    orange
```



## Workflow - VBC Galaxy

The VBC Galaxy instance is not public. These instructions are for internal institute use only. 
For general use, please consult the next section: "Workflow - Command line".

### Step 1 : Quantify barcodes

**Groups/Obenauf/catch barcode quantifier**

* Each BAM file needs to be quantified separately. 
* Each multiplexed BAM using Christian's design requires its own description table with just the relevant samples. Only the `Tag` and `Sample` columns are required, the other 3 columns are optional for this step.

Quantification will generate a read counts table for each sample. A summary table for each BAM will also be created.

### Step 2: Collect outputs

**Groups/Obenauf/catch barcode merger**

The individual sample outputs need to be combined: One collective table for all the barcode read counts, and one 
collective table for the summaries. Simply select the Step-1 Galaxy outputs to combine.

### Step 3: Visual report

**Groups/Obenaud/catch report**

* This requires the two collective tables (barcode counts and summaries) from Step 2. 
* Additionally it needs the description table with the samples to be included in the report. The `Tag` column is optional.




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
    - ggridges
    - RColorBrewer
    - patchwork


### Step 1 : Quantify barcodes

**barcodingQuantifier.py**

Starting from a BAM file and an optional table with demultiplexing info, it uses regex matching 
to extract the barcode sequence.

### Step 2 : Merge barcodes within a certain edit distance (optional)

**barcodingHammingMerge.py**

The idea is to reduce the noise from sequencing errors. The merge is based on the Hamming distance 
between barcodes. That means that only substitutions are considered, no insertions or deletions. When
barcodes are merged, the sequence of the most abundant is the one that is kept.
It is up to your judgement to apply this step or not, as well as to choose the edit distance.

This step is very slow and resource hungry.

### Step 3 : Collect outputs

**mergeBCcounts.R**

Full-join/merge the tables from Step 1 (or from Step 2 instead, if applicable).

Please make one collective table for the barcode counts and one collective table for the summaries. Both are needed 
for the report. so you will need to *run this step twice*, once for all the barcode count files, and once for 
all the summary files.

### Step 4 : Visual report

**barcoding_results_run.R**

This will compile a report using the **barcoding_results_template.Rmd** template. It creates:

* an HTML report with interactive figures,
* optionally a PDF static copy of those figures,
* a table with the counts and percentages of barcodes and their assigned IDs, and
* a table with the pairwise sample correlations.


### Steps 1-4 : Pipeline

**catch_workflow.sh**

The analysis can be executed step by step, as described in the previous sections. That is the most likely use case for you.

The Bash script included here is my rough automation of steps 1-4. The script depends on a SLURM environment, and assumes all the software dependencies have been made available in advance. This script is entirely optional. **I cannot promise it will work for you**, it probably won't. But at the very least it can serve as a template from which to adapt a new script that is suitable for your system.  
 



#### Parameters:

```
-b DIR          Directory with the BAM file(s) to be analysed.
-c FILE         Demultiplexing table, if applicable.
-o DIR          Output directory for the counts files.
-O DIR          Ouptut directory for the analysis report.
-v FILE         Description table.
-X DIR          Where to find all the scripts for this workflow (probably the path to your local clone of the repository).
-m INT          Hamming distance at which to merge barcodes as likely sequencing errors, if applicable (0 ie. not applicable).
-n INT          Number of dark bases to allow in the pattern (0).
-A FLOAT        Barcode abundance threshold for the report. Barcodes must make up at least this percentage of the sample's reads to be considered abundant (0.005 ie. 0.5% of reads).
-R INT_LIST     Comma-seperated list of rows in the description table to be used as reference samples in the report, NOT counting the header line (1).
-r              Reverse complement the sample tags (False).
-i              Spike-in barcode was added (False). The spike sequence is hard-coded in the quantifier script.
-s              Match the full pattern of the barcodes instead of the short pattern (False).
```


