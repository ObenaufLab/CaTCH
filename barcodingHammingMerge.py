#!/usr/bin/env python3

# Date located in: -
# from __future__ import print_function
import sys, os, re
from collections import Counter
from argparse import ArgumentParser, RawDescriptionHelpFormatter

from Levenshtein import hamming

# It would be cleaner to apply the distance limit directly while quantifying the barcodes. 
# But for the benefit of keeping fewer barcodes in the dictionary, the distance must then be computed for every read.
# Instead, letting the quantifier tally barcodes at identity level and then measuring distances only for the distinct barcodes, reduces the number fo distances calculated.
# It is still a very slow process.

usage = "Hamming distance merge"
parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-b", "--barcodes", type=str, required=True, dest="barcodesFile", help="Tab delimited file: barcode_sequence \\t count")
parser.add_argument("-d", "--hammDist", type=int, default=1, help="Hamming distance at which barcodes should be considered the same (1)")
parser.add_argument("-o", "--outfile", type=str, required=True, help="Output file.")

args = parser.parse_args()

# Keep stats on BC distances.
hbcstat = Counter()  # Number of barcodes at given distance. Loses abundance information of the barcodes.
hrmin = Counter()    # Min read-count among the barcodes in a distance bin.
hrmax = Counter()    # Max read-count among the barcodes in a distance bin.

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
                # Hamming distance is only for equal lengths. Different lengths should not occur under current CaTCH design, but may be diagnostic counts or unforseen irregularities.
                # Explicitly exclude any diagnostic categories from having their distances compared.
                if len(b) == len(barcode) and b not in ["unknown", "SampleUnknown", "unmatched", "BCUnmatched", "empty", "EmptyVector", "spike", "SpikeIn"]:   
                    h = hamming(b, barcode)
                    sh = str(h)
                    hbcstat.update([sh])
                    if sh in hrmin:
                        hrmin[sh] = min([hrmin[sh], count, bc[b]])
                        hrmax[sh] = max([hrmax[sh], count, bc[b]])
                    else:
                        hrmin[sh] = count
                        hrmax[sh] = count
                    if len(b) == len(barcode) and h <= args.hammDist and not found:
                        if bc[b] >= count:
                            # If existing sequence is more abundant, update its count and discard the new sequence
                            # This is the most likely scenario, as the quantifier orders the barcodes by decreasing abundance.
                            bc[b] = bc[b] + count
                            found = True
                            # break
                        else:
                            # Add the new sequence with combined count, and remove the old sequence
                            bc[barcode] = count + bc[b]
                            bc.pop(b)
                            found = True
                            # break
            if not found:
                bc[barcode] = count

##############
# Output
##############

with open(args.outfile, 'w') as f:
    for barcode in bc:
        f.write(barcode + "\t" + str(bc[barcode]) + "\n")



def natural_sorted(l):
    """Sort list of numbers/strings in human-friendly order.

    Args:
        l(list): A list of strings.
    Returns:
        list
    """
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)

# print( args.outfile.replace('barcode-counts.txt', 'hammdist.txt') )
with open(args.outfile.replace('bardode-counts.txt', 'hammdist.txt'), 'w') as f:
    # f.write("# These are not based on all the pairwise distances, but rather on the distances encountered until a barcode hits a match.\n")
    f.write("# These are based on all the pairwise distances, each pair of barcodes measured once.\n")
    f.write("HammDist\tBCs\tminReads\tmaxReads\tSample\n")
    for sh in natural_sorted(hbcstat.keys()):
        # f.write(h[0] + "\t" + str(h[1]) + "\n")
        f.write(sh + "\t" + str(hbcstat[sh]) + "\t" + str(hrmin[sh]) + "\t" + str(hrmax[sh]) + "\t" + os.path.basename(args.outfile).replace('_bardode-counts.txt', '') + "\n")
