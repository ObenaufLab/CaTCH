# Step 2 : Merge barcodes within a certain edit distance (optional)

The idea is to reduce the noise from sequencing errors. The merge is based on the Hamming distance between the already quantified barcodes. 
That means that only substitutions are considered, no insertions or deletions. When barcodes are merged, the sequence of the most abundant is the one that is kept.

It is up to your judgement whether to apply this step or not, as well as what edit distance to use.

This step is very slow and resource hungry.

In practice we don't use this step. For barcodes of the length we typically use, correct sequence barcodes are much more abundant than their versions with sequencing errors.
This also avoids accidentally merging valid barcodes. 
The problem of noise is handled by focusing on the most prominent sequences with the use of readcount and proportion filters instead.


## Command line

**barcodingHammingMerge.py**

### Parameters:

```
-b FILE       A tab delimited file with barcode counts, as produced by the quantifier in Step 1.
-d INT        Hamming distance at which barcodes should be considered the same sequence (Default: 1)
-o FILE       Output file.
```

## VBC Galaxy

This step is not included in Galaxy. 