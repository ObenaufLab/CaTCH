These are the analysis scripts for the CaTCH experiments.
They should include everything needed to go from an unaligned .BAM file to barcode counts and an HTML report with various plots, some of which are interactive.

================
# Requirements #
================

The exact software versions are provided for reference. Most likely many other combinations of versions should work too.

- Python 3 (3.6.7):
    - os
    - sys
    - re
    - argparse
    - collections
    - pysam
    - Biopython
    - pandas
    - Levenshtein

- R (3.5.1) :
    - tidyverse
    - data.table
    - plotly
    - ggrepel
    - RColorBrewer
    - pheatmap

- The python and R scripts provided here.


=========
# Input #
=========

- One or multiple BAM files containing the CaTCH reads. Each BAM file can be multiplexed or a single sample. No alignment is needed.

- In the case of multiplexed BAMs, a tab-delimited file containing a table matching the multiplexing tags to sample names. The table should contain 2 named columns: "barcode" and "sample", for the sample tag and sample name respectively.
  ie.

    barcode sample
    TGAG    reference_1
    TCGG    reference_2
    ATCACG  notreatment_1
    ATGGCG  notreatment_2
    ATCTAG  treatmentA_1
    CGAT    treatmentA_2

  Demultiplexing tags can be a mix of different lengths, but the shorter ones must not be substrings of the longer ones.

- A file containing the condition allocation for each sample to be included in the report and the assigned colour for the condition. It must consist of 3 named columns: "sample", "condition", "colour". Colour must be in an R-compatible string format, one colour per condition.
  ie.

    sample  condition   colour
    reference_1 reference   black
    reference_2 reference   black
    notreatment_1   notreatment blue
    notreatment_2   notreatment blue
    treatmentA_1    treatmentA  red
    treatmentA_2    treatmentA  red


============
# Workflow #
============


# Step 1 : Count the barcodes
=============================
barcodingQuantifier.py
----------------------

Two barcode designs are programmed into the script:

(1) P7adapter_template_BARCODE_GenotypeTag_SampleTag_NNNNNN_P5adapter  
                                                     <--- read direction

  In this design the sample tag is included in read1 of the mate pair, and the samples can be demultiplexed by the counting script, if a sample tag table is provided (see input section above). The sample tags are expected to be found at a fixed interval from where the barcode is located.

(2) P7adapter_SampleTag_template_BARCODE_GenotypeTag_P5adapter
                                              <--- read direction
                                                
  In this design the sample tag is unlikely to be included in read1 containing the barcode, and therefore must be demultiplexed _in_advance_ by other means, not provided here.


The genotype tag is not used in the current implementation, except as spacer between the barcode and the sample tag in the first design.

Two versions of the barcode-matching pattern are programmed in, one more stringent than the other:

(1) stringent: 
[TC]{1}[ATGC]{2}[AT]{1}[TG]{1}[TGC]{1}[AT]{1}[TGC]{1}[ATGC]{1}[TGC]{1}[ATGC]{1}[ATC]{1}[TGC]{1}[ATGC]{3}[CG]{1}[ATGC]{1}[CGA]{1}CGCCG[TC]{1}[AGTC]{2}[AT]{1}[TG]{1}[TCG]{1}[AT]{1}GCC[ATGC]{1}[ACT]{1}[ATGC]{4}[GC]{1}[AGTC]{1}[CGA]{1}CGCCG'

(2) nonstringent: 
'CGCCG[ATGC]{7}GCC[ATGC]{9}CGCCG'

The stringent pattern matches from the beginning of the barcode to the end. The non-stringent matches the less variable sequence in the middle of the barcode, and the start of the barcode is found with a fixed offset of 19nt before the pattern.

There are currently no parameters to specify a different barcode pattern or different offset of the pattern from the start of the barcode.

Parameters:
-----------
-f --file      A BAM file.
-b --barcodes  File with a demultiplexing table. Completely omit this if data is already demultiplexed.
-r --revcomp   Reverse complement the provided sample tags (default: false)
-i --spikein   Barcode GTGTGTGGAACGAGCACAGCGCCGAGAGACGGATATCACTAGTCGCCGCCATTTGCGCGCGCTCGCC was added as spiked-in (default: false). If true it will be matched explicitly.
-s --stringent Use the more stringent barcode pattern (default: false).
-o --outdir    Output directory for counts (default: ./process/)
--bc_len       Length of the barcode from start position of the match (default: 67).
--gt_len       Genotyping tag length (default: 20).


# Step 2 : Merge barcodes within a certain edit distance (optional)
===================================================================
barcodingHammingMerge.py
------------------------

The merge is based on the Hamming distance between barcodes. That means that only substitutions are considered, no insertions or deletions.
This step is very slow.
It is up to your judgement to apply this step or not, as well as to choose the edit distance.

Parameters:
-----------
-b --barcodes A tab delimited file with barcode counts, as produced by the previous step.
-d --hammDist Hamming distance at which barcodes should be considered the same (Default: 1)
-o --outfile  Output file.


# Step 3 : Collect all the counts into one table
================================================
mergeBCcounts.R
---------------

This is simply a matter of full-merging the tables from Step 1 (or from Step 2, if you choose to apply it) and can be achieved in many different ways. 
The way provided here with this script assumes that all the count files can be found in a common path and have a common pattern in the name.

For the barcode counts files, that pattern should be "_barcode-counts.txt"
For the step1 summaries, that pattern should be "_summary.txt"

Please make one collective table for the barcode counts and one for the summaries. Both are needed for the next step.


Parameters (positional):
------------------------
1. directory path in which to search for barcode counts files with "_barcode-counts.txt" in their name.
2. output file
3. file pattern


# Step 4 : Compile the report
=============================
barcoding_results_run.R
-----------------------

This will compile a report using the barcoding_results_template.Rmd template.

Parameters:
-----------
-c  The collective counts table from step 3.
-s  The collective summaries table from step 3.
-d  Directory in which to save all output files.
-T  The path to barcoding_results_template.Rmd .
-v  The file assigning samples to conditions and colours (see Input section above)
-p  Instead of an HTML report (more information and more interaction), output just the figures into a PDF file (for Illustrator or presentations).
-r  Comma-separated list (without spaces) of integers designating the samples to use as reference abundance. The numbers should correspond to the order in which  the samples appear in the collective summaries file) (Default: 1).
-N  Count threshold for barcodes (Default: 50). The counting generates many low-count barcodes, as a result of sequencing errors. These inflate the number of barcodes, so this threshold is provided to cut out that noise.
-A  Proportional abundance threshold to consider barcodes top hits (Default: 0.01). Range 0 - 1.

