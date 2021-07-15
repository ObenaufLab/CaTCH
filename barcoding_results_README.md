# Step 4 : Visual report

This creates:

* an HTML report with interactive figures,
* a PDF with static copies of those figures,
* a table with the counts and percentages of barcodes and their assigned IDs, and
* a table with the pairwise sample correlations.

In addition to the collective count and summary tables from the previous step, the report also needs the description table:

1. `Tag` - The sample index. This column is not used in this step, so it can be omitted or contain dummy values like XXXX.
2. `Sample` - The unique name of the sample.
3. `Group` - A coarse grouping of the samples. For example if you have multiple different treatments and/or mutliple controls, you can group them as treated/untreated.
4. `Treatment` - The specific treatment far that sample.
5. `Colour` - This is for display purposes in plots. You can use any colour name recognised by R. A suggested colourset to choose from could be: black, red, orange, dodgerblue, steelblue, forestgreen, darkgold, magenta, purple. For more precise control, hexadecimal RGB values in the format #RRGGBB should be possible as well.

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

* For columns 2, 3 and 4 the values must be plain ASCII alphanumeric strings, starting with a letter. Underscores `_` are allowed. No accented or special letters, no spaces, no other symbols (including dashes).
* Column 1 is optional.


**Important:** This step allows you to specify one or more samples as reference samples. If you provide this, barcodes that were not seen in the reference(s) will be filtered out from all the other samples as well. Use this only if you have samples collected at the beginning of the screen and you have sufficient sequencing depth to have confidence that no barcodes were missed.

The reference sample(s) will also be used to give short names to the barcodes, in order of decreasing abundance. If no references are specified, the barcodes will be named by order of abundance in the first sample. This does not affect the results per se, but can confuse comparison of results from different runs, if the barcode names are based on different samples. However, the results do also include the barcode sequences, which can be used to match the results across runs.


## Command line 

**barcoding_results_run.R**

This will compile an HTML report using an Rmd template, by default the provided **barcoding_results_template.Rmd**.

Please provide **full paths** to the files and directories.

#### Parameters:

```
-c FILE         The collective counts table from step 3.
-s FILE         The collective summaries table from step 3.
-d DIR          Directory in which to save all output files.
-T FILE         The path to barcoding_results_template.Rmd (or other suitable template that takes the same parameters).
-v FILE         The description table.
-p              In addition to the HTML report, output static figures into a PDF file for Illustrator or presentations (Default: false).
-r INT_LIST     Comma-separated list (without spaces) of integers designating the samples to use as reference abundance.
                The numbers should correspond to the order in which  the samples appear in the collective summaries file, not counting the header line) (Default: 1).
-N INT          Count threshold for barcodes (Default: 50). The counting generates many low-count barcodes, as a result of sequencing errors.
                These inflate the number of barcodes, so this threshold is provided to cut out that noise.
-A FLOAT        Proportional abundance threshold to consider barcodes top hits (Default: 0.01). Range 0 - 1.
-B BC_LIST      Comma-separated list of additional barcode IDs that should be included in the report despite not being among the top shared ones.
                These are the IDs assigned to the barcodes by the report, so you'll need to run the report at least once before you know which extra IDs to provide.
```

## VBC Galaxy (not public)

**GROUP/Obenauf/catch report**

### Cluster options

`Walltime`: Estimated completion time. This step takes a few minutes.

`Memory`: Estimated memory use. The default 4GB should be plenty. Increase this if you get slurm out-of-memory (oom) / memory-exceeded errors.

### Main options

`Collected Barcode Counts`: The collective table of barcode readcounts created in the previous step.

`Collected Barcode Summaries`: The collective tablee of summaries created in the previous step.

`Conditions and Colours`: The description table.

### Advanced options

`Count threshold for barcodes`: Minimum number of reads required for a barcode to be worth considering.

`Proportional abundance threshold`: Minimum percentage of a sample's reads required for a barcode to be considered enriched.

`Comma-separated list of integers designating the samples to use as reference abundance`: In order of appearance in the description table. Optional.

