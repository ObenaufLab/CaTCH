#!/usr/bin/env python

# Date located in: -
from __future__ import print_function
import sys, os, re

import subprocess

from argparse import ArgumentParser, RawDescriptionHelpFormatter

from os.path import basename

from string import upper

usage = "Annotate annovar files with VAF and normal depth for maftools"
parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-i", "--input", type=str, required=True, dest="inputFile", help="Annovar input file")

args = parser.parse_args()

##############
# Read annovar file
##############

with(open(args.inputFile, 'r')) as f:

    header = next(f)
    header = header.rstrip()
    header = header.split("\t")
    numFields = len(header)
    #header += "\tNormalCoverage\tVAF"

    print('\t'.join(header),end="\t")
    print('VAF')

    for line in f:

        fields = line.rstrip().split("\t")
        gt = fields[23]
        VAF = gt.split(":")[6]

        listPointer = 0

        while listPointer < numFields:
            print(fields[listPointer], end = "\t")
            listPointer += 1

        print(VAF, end = "")

        while listPointer < len(fields):
            print("\t" + fields[listPointer], end = "")
            listPointer += 1

        print()
