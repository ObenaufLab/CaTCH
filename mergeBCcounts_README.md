### Step 3 : Collect outputs

Full-join/merge the tables from Step 1 (or from Step 2 instead, if applicable).
You can provide as many input tables as you want.


## Command line

**mergeBCcounts.R**

Please make one collective table for the barcode counts *and* one for the summaries. Both are needed 
for the report. That means running this script twice, once for all the barcode count files, and once for 
all the summary files.

### Usage:

```
mergeBCcounts.R OUTPUT_FILE INPUT_FILE INPUT_FILE ...
```

## VBC Galaxy (not public)

**GROUP/Obenauf/catch barcode count merger**

The Galaxy implementation of this step merges all the counts tables together and all the summary tables together.

### Cluster options

`Walltime`: Estimated completion time. This step takes a few seconds.

`Memory`:   Estimated memory use. The default 4GB should be plenty. Increase this if you get slurm out-of-memory (oom) / memory-exceeded errors. It will depend on the  number of samples and number of identified barcodes.


### Main options

`barcode counts`: Select the output from the quantification step. You can also upload samples from other quantification runs as individual table files or as a zipped bundle. **NOTE:** in order to be recognised, the tables must have filenames that end with `_barcode-counts.txt` and `_summary.txt`. If you use outputs from another CaTCH quantification step, this will be already done.


