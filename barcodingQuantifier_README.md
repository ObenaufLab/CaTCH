# Step 1 : Quantify barcodes

Starting from a BAM file and an optional table with demultiplexing info, it uses regex matching 
to extract the barcode sequence. It creates:

* a barcode counts file per sample, and 
* an extraction summary file per BAM.

Additionally, for troubleshooting use, it outputs:

* a BAM with reads where a barcode couldn't be identified (truncation, or mismatch in inflexible bases),
* a BAM file for reads whose sample tags didn't match the demultiplex info, 
* a text file with the tally of the sample tags found in reads that couldn't be matched to a sample, and 
* a BAM file with reads that matched the empty or spike-in patterns.

## Input

### Sequencing data

The input data is one BAM file containing sequencing reads. Multiple BAM files have to be quantified individually.

#### Two barcode designs have been tested:

Christian Umkehrer's original design allows the sample index to be contained in the first read. This design can be demultiplexed and quantified in a single step with the current script.

*5'-P5adapter_NNNNNN_SampleINDEX_GenotypeTag_BARCODE_template_P7adapter-3'*

In an effort to move the barcode closer to the start of the read for higher quality, Shona Cronin's design moves the sample index to the end of the construct. It is unlikely that the first read will contain the sample index intact,
so pooled lanes with this design need to be demultiplexed in advance using the second read of the pair.

*5'-P5adapter_GenotypeTag_BARCODE_template_SampleIndex_P7adapter-3'*

### Data description

A table with the following columns and in that order is helpful to have from the beginning covering all the samples to be analyzed:

1. `Tag` - The sample index. This should be a valid nucleotide [ATGC]. A mix of Tag lengths is allowed, but the shorter ones should not be substrings of the longer ones, to prevent misidentifications. If mixing multiplexed and non-multiplexed BAMs, use a dummy Tag like `XXXX` for Samples that are already in their own BAM.
2. `Sample` - A name for the sample.
3. `Group` - A coarse grouping of the samples. For example if you have multiple different treatments and/or mutliple controls, you can group them as treated/untreated.
4. `Treatment` - The specific treatment far that sample.
5. `Colour` - This is for display purposes in plots. You can use any colour name recognised by R.

For column 2 the values must be plain ASCII alphanumeric strings, starting with a letter. Underscores (`_`) are allowed. 
No accented or special letters, no spaces, no other symbols (including dashes).

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

* If the BAM file is not multiplexed, the description table is irrelevant and **should be ommitted** completely for quantifying that file.
* Each multiplexed BAM file will need its own demultiplexing table. The easiest way to do it is to create a duplicate of the description table and delete the rows that are not samples contained in the given BAM file.
* If none of the sample tags is reused for different samples in diffeerent lanes, you may be able to use the description table in its original form for each BAM file, without needing to create duplicates with subsets of the samples.


## Command line

**barcodingQuantifier.py**

### Parameters:

```
-f FILE         A single BAM file.
-b FILE         Optional demultiplexing table for use with Christian's design. This may be the description table or a subset of it.
-o DIR          Output directory for counts (./process/).
--revcomp       Reverse complement the sample tags.
--bc_len INT    Length of the barcode (not of the whole read) (68).
--gt_len INT    Number of nucleotides between the end fo the sample tag and the start of the barcode (20). It allows a spacer between barcode and sample tag.
                This is tested only for sample tags that precede the barcode. 
                Theoretically, negative values should work too and allow the sample tag to be located after the start of the barcode, as long as the end of the tag is at a fixed offset from the start of the barcode.
--spikein       Also search for the spike-in barcode. The Obenauf lab design for this is hard-coded.
--stringent     Use the full-length pattern instead of the core pattern. This is applicable to the hard-coded Obenauf-lab designs only.
--n_dark INT    Number of dark/unknown bases (N) to consider in the barcode pattern (0). 
                For this to have an effect, the provided --pattern must be designed to also match Ns.
                There is a hard-coded version of the short Obenauf-lab design for this. Too high Ns will make recognition of barodes unreliable.
--pattern       Barcode pattern as a regex string. This overrides the hard-coded patterns to enable other designs.
                The regex must begin at the start of the barcode, but doesn't need to reach the end of the barcode.
                bc_len, gt_len, n_dark are available to use with a custom pattern. 
                This features has not been extensively tested.
```



## VBC Galaxy (not public)

**GROUPS/Obenauf/catch barcode quantifier**

### Cluster options

`Walltime`: Estimated completion time. Quantification runs fairly fast and providing a shorter estimate will help the job get allocated sooner. Though that depends on the size of the BAM file, an hour should be more than adequate, often may be just a few minutes.
`Memory`:   Estimated memory use. Increase this if you get slurm out-of-memory (oom) / memory-exceeded errors. Again it will depend on the size of your BAM files.

### Main options

`read files`: A single BAM file.
`barcodes`:   Optional demultiplexing table for use with Christian's design. This may be the description table or a subset of it. Omit completely if the BAM contains a single sample.

### Advanced options

`rev-comp sample tags`:       Reverse complement the sample tags. Try this if you get unexpectedly low counts across multiple samples.
`Stringent barcode matching`: Use the longer more-explicit pattren to match barcodes.
`Number of dark bases`:       Number of undefined (N) bases to allow. This uses a special version of the short pattern and is quite low-tech, keep this number low to prevent unreliable recognition of the barcodes.
`Length of the barcode`:      Self-explanatory.
`Genotyping tag`:             Length of the sequence between the sample tag and the barcode. Intended as paceholder for genotyping tags, it is currently treated only as a spacer.

