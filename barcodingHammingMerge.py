#!/usr/bin/env python

# Date located in: -
from __future__ import print_function
import sys, os, re

import subprocess

from argparse import ArgumentParser, RawDescriptionHelpFormatter

from os.path import basename

from Bio import SeqIO
from Bio.Seq import Seq
from string import upper
import Levenshtein

usage = "Hamming distance merge"
parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-b", "--barcodes", type=str, required=True, dest="barcodesFile", help="Tab delimited barcode counts file")

args = parser.parse_args()

##############
# Read bc
##############

bc = dict()

with(open(args.barcodesFile, 'r')) as f:
    for line in f:
        barcode, count = line.rstrip().split("\t")
        if int(count) > 1:
            bc[barcode] = int(count)

merge = True

print("Starting neighbor joining done.",file=sys.stderr)
restartCounter = 0

while merge:
    merge = False
    index = 0
    barcodes = bc.keys()
    barcodes.sort()
    while index < (len(barcodes) - 1):
        base = barcodes[index]
        test = barcodes[index + 1]

        if len(base) == len(test) :

            if Levenshtein.hamming(base, test) == 1:
                baseAb = bc[base]
                testAb = bc[test]

                if baseAb > testAb:
                    if testAb >= float(baseAb) / 8:
                        bc[base] += bc[test]
                        del bc[test]
                        merge = True
                else:
                    if baseAb >= float(testAb) / 8:
                        bc[test] += bc[base]
                        del bc[base]
                        merge = True

            if Levenshtein.hamming(base, test) == 2:
                baseAb = bc[base]
                testAb = bc[test]

                if baseAb > testAb:
                    if testAb >= float(baseAb) / 40:
                        bc[base] += bc[test]
                        del bc[test]
                        merge = True
                else:
                    if baseAb >= float(testAb) / 40:
                        bc[test] += bc[base]
                        del bc[base]
                        merge = True

        if merge:
            print("Merged. Restart " + str(restartCounter) + ".",file=sys.stderr)
            restartCounter += 1
            break

        index += 1

print("Lexicographical neighbor joining done.", file=sys.stderr)
print("Starting computationally intense n2 hamming search.",file=sys.stderr)

restartCounter = 0

merge = True

while merge:
    merge = False
    cur = 1
    prev = 0
    barcodes = bc.keys()
    barcodes.sort()
    while prev < (len(barcodes) - 1):
        while cur < len(barcodes):
            base = barcodes[prev]
            test = barcodes[cur]

            if len(base) == len(test) :

                if Levenshtein.hamming(base, test) == 1:
                    baseAb = bc[base]
                    testAb = bc[test]

                    if baseAb > testAb:
                        if testAb >= float(baseAb) / 8:
                            bc[base] += bc[test]
                            del bc[test]
                            merge = True
                    else:
                        if baseAb >= float(testAb) / 8:
                            bc[test] += bc[base]
                            del bc[base]
                            merge = True

                if Levenshtein.hamming(base, test) == 2:
                    baseAb = bc[base]
                    testAb = bc[test]

                    if baseAb > testAb:
                        if testAb >= float(baseAb) / 40:
                            bc[base] += bc[test]
                            del bc[test]
                            merge = True
                    else:
                        if baseAb >= float(testAb) / 40:
                            bc[test] += bc[base]
                            del bc[base]
                            merge = True

            if merge:
                print("Merged. Restart " + str(restartCounter) + ".",file=sys.stderr)
                restartCounter += 1
                break

            cur += 1

        if merge:
            break

        prev += 1

print("Computationally intense n2 hamming search done.",file=sys.stderr)

barcodes = bc.keys()
barcodes.sort()

for barcode in barcodes:
    print(barcode + "\t" + str(bc[barcode]))
