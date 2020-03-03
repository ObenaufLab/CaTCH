#!/usr/bin/env python3

# Date located in: -
# from __future__ import print_function
import sys, os, re
from argparse import ArgumentParser, RawDescriptionHelpFormatter

from Levenshtein import hamming

# It would be cleaner to apply the distance limit directly while quantifying the barcodes.
# But that would explode the number of comparisons done, in exchange for having some fewer barcodes in the dictionary.
# Letting the quantifier count all barcodes at identity level and then filter the counted barcodes here in a separate set might be more efficient than measuring the distance for every read.

usage = "Hamming distance merge"
parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-b", "--barcodes", type=str, required=True, dest="barcodesFile", help="Tab delimited file: barcode_sequence \\t count")
parser.add_argument("-d", "--hammDist", type=int, default=1, help="Hamming distance at which barcodes should be considered the same (1)")
parser.add_argument("-o", "--outfile", type=str, required=True, help="Output file.")

args = parser.parse_args()

##############
# Read bc and merge
##############

# Doing it as barcodes are parsed, reduces the size of the dictionary.
bc = dict()
first = True
with open(args.barcodesFile, 'r') as f:
    for line in f:
        barcode, count = line.rstrip().split("\t")
        count = int(count)
        if first:
            bc[barcode] = count
            first = False
        else:
            found = False
            for b in list(bc):
                if len(b) == len(barcode) and hamming(b, barcode) <= args.hammDist:
                    if bc[b] >= count:
                        # If existing sequence is more abundant, update its count and discard the new sequence
                        # This is the most likely scenario, as the quantifier orders the barcodes by decreasing abundance.
                        bc[b] = bc[b] + count
                        found = True
                        break
                    else:
                        # Add the new sequence with combined count, and remove the old sequence
                        bc[barcode] = count + bc[b]
                        bc.pop(b)
                        found = True
                        break
            if not found:
                bc[barcode] = count

##############
# Output
##############

with open(args.outfile, 'w') as f:
    for barcode in bc:
        f.write(barcode + "\t" + str(bc[barcode]) + "\n")
