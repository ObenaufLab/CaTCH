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

### Cluster options

```
Walltime: Estimated completion time. This step takes a few seconds.
Memory:   Estimated memory use. The default 4GB should be plenty. Increase this if you get slurm out-of-memory (oom) / memory-exceeded errors. It will depend on the  number of samples and number of identified barcodes.
```

### Main options

```
barcode counts: Either all the barcode counts files from Step 1, or all the summary files from Step 1.
```

You will need to run it separately for the barcode counts and for the summaries.
