#!/usr/bin/env python2.7

# Demultiplex samples and count the semi-random barcode occurrences in each.

# Date located in: -
from __future__ import print_function
import sys, os, re

import subprocess

from argparse import ArgumentParser, RawDescriptionHelpFormatter

import os.path

from Bio import SeqIO
from Bio.Seq import Seq
from string import upper
import Levenshtein

usage = "Barcoding"
parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-f", "--fastq", type=str, required=True, dest="fastqFile", help="Fastq file")
parser.add_argument("-b", "--barcodes", type=str, required=True, dest="barcodesFile", help="Demultiplexing 4nt-barcodes table (tab delimited: Barcode \\t Sample, with header line)")
parser.add_argument("-6", "--6bpbarcodes", type=str, required=False, dest="barcodes6bpFile", help="Demultiplexing 6nt-barcodes table (same format as above)")
parser.add_argument('-r', "--revcomp", action='store_true', dest="revcomp", help="Reverse complement barcodes (default: false)")
parser.add_argument('-i', "--spikein", action='store_true', dest="spikein", help="Barcode was spiked-in (default: false)")
parser.add_argument('-s', "--stringent", action='store_true', dest="stringent", help="Stringent barcode matching (default: false)")
parser.add_argument('-o', "--outdir", type=str, required=False, dest="outdir", default="./process", help="Output directory for counts (./process/)")

args = parser.parse_args()

##############
# Read bc
##############

bc = dict()
bcStats = dict()

with(open(args.barcodesFile, 'r')) as f:
    # skip header
    next(f)
    for line in f:
        barcode, sample = line.rstrip().split("\t")
        bc[barcode] = sample
        bcStats[barcode] = 0

bc["unmatched"] = "unmatched"
bcStats["unmatched"] = 0

##############
# Read 6bp barcodes
##############

if args.barcodes6bpFile:

    with(open(args.barcodes6bpFile, 'r')) as f:
    # skip header
        next(f)
        for line in f:
            barcode, sample = line.rstrip().split("\t")
            bc[barcode] = sample
            bcStats[barcode] = 0

if args.revcomp:

    bcRevCompStats = {}
    bcRevComp = {}

    for barcode in bc:
        if barcode != "unmatched":
            bcRevComp[str(Seq(barcode).reverse_complement())] = bc[barcode]
            bcRevCompStats[str(Seq(barcode).reverse_complement())] = 0
        else :
            bcRevComp["unmatched"] = "unmatched"
            bcRevCompStats["unmatched"] = 0

    bc = bcRevComp
    bcStats = bcRevCompStats

library = dict()

for multiplex in bc:
    library[multiplex] = dict()

if args.barcodes6bpFile:

    for multiplex in bc:
        library[multiplex] = dict()

fout = open(os.path.join(args.outdir, "unmatched.fastq"), "w")

for record in SeqIO.parse(args.fastqFile, "fastq"):

    # GBNSNNNVDNVNVWVMWNNRCGGCGBNSNNNNDNGGCWVMWNNRCGGCGBNSNNNVDNVNVWVMWNNR <- sequencing direction
    # RNNWMVWVNVNDVNNNSNB CGCCG RNNWMVW GCC NDNNNNSNB CGCCG

    if args.stringent:

        match = re.search('[TC]{1}[ATGC]{2}[AT]{1}[TG]{1}[TGC]{1}[AT]{1}[TGC]{1}[ATGC]{1}[TGC]{1}[ATGC]{1}[ATC]{1}[TGC]{1}[ATGC]{3}[CG]{1}[ATGC]{1}[CGA]{1}CGCCG[TC]{1}[AGTC]{2}[AT]{1}[TG]{1}[TCG]{1}[AT]{1}GCC[ATGC]{1}[ACT]{1}[ATGC]{4}[GC]{1}[AGTC]{1}[CGA]{1}CGCCG', str(record.seq))

        if match:
            barcode = str(record.seq[match.start(0):match.start(0) + 67])
            genotyping = str(record.seq[match.start(0) - 20:match.start(0)])
            multiplex = str(record.seq[match.start(0) - 24:match.start(0) - 20])
            multiplex6p = str(record.seq[match.start(0) - 26:match.start(0) - 20])

            #print(barcode)
            #print(genotyping)
            #print(record.seq)

            if multiplex in bc:
                bcStats[multiplex] += 1
                if not barcode in library[multiplex]:
                    library[multiplex][barcode] = 0
                library[multiplex][barcode] += 1
            elif multiplex6bp in bc:
                bcStats[multiplex6bp] += 1
                if not barcode in library[multiplex6bp]:
                    library[multiplex6bp][barcode] = 0
                library[multiplex6bp][barcode] += 1
            else:
                bcStats["unmatched"] += 1
                if not barcode in library["unmatched"]:
                    library["unmatched"][barcode] = 0
                library["unmatched"][barcode] += 1

        else :
            # Emtpy vector
            match = re.search('AGAGACGGATATCACTAGTCGTCTCCGTTCGCTCTAGACAGGGTACCCAGCATATGATAGGGTCCCCT', str(record.seq))

            if match:

                barcode = str(record.seq[match.start(0):match.start(0) + 67])
                genotyping = str(record.seq[match.start(0) - 20:match.start(0)])
                multiplex = str(record.seq[match.start(0) - 24:match.start(0) - 20])
                multiplex6bp = str(record.seq[match.start(0) - 26:match.start(0) - 20])

                if multiplex in bc:
                    bcStats[multiplex] += 1
                    if not barcode in library[multiplex]:
                        library[multiplex][barcode] = 0
                    library[multiplex][barcode] += 1
                elif multiplex6bp in bc:
                    bcStats[multiplex6bp] += 1
                    if not barcode in library[multiplex6bp]:
                        library[multiplex6bp][barcode] = 0
                    library[multiplex6bp][barcode] += 1
                else:
                    bcStats["unmatched"] += 1
                    if not barcode in library["unmatched"]:
                        library["unmatched"][barcode] = 0
                    library["unmatched"][barcode] += 1
            else :
                if args.spikein:
                # Spike-in barcode
                #  CCTAAAGCTTCTCCTGCCG GTGTGTGGAACGAGCACAGCgccgAGAGACGGATATCACTAGTCgccgCCATTTGCGCGCGCTCGCC
                # 4mer index 20nt genotyping GTGTGTGGAACGAGCACAGCgccgAGAGACGGATATCACTAGTCgccgCCATTTGCGCGCGCTCGCC

                    match = re.search('GTGTGTGGAACGAGCACAGCgccgAGAGACGGATATCACTAGTCgccgCCATTTGCGCGCGCTCGCC'.upper(), str(record.seq))

                    if match:

                        barcode = str(record.seq[match.start(0):match.start(0) + 67])
                        genotyping = str(record.seq[match.start(0) - 20:match.start(0)])
                        multiplex = str(record.seq[match.start(0) - 24:match.start(0) - 20])
                        multiplex6bp = str(record.seq[match.start(0) - 26:match.start(0) - 20])

                        if multiplex in bc:
                            bcStats[multiplex] += 1
                            if not barcode in library[multiplex]:
                                library[multiplex][barcode] = 0
                            library[multiplex][barcode] += 1
                        elif multiplex6bp in bc:
                            bcStats[multiplex6bp] += 1
                            if not barcode in library[multiplex6bp]:
                                library[multiplex6bp][barcode] = 0
                            library[multiplex6bp][barcode] += 1
                        else:
                            bcStats["unmatched"] += 1
                            if not barcode in library["unmatched"]:
                                library["unmatched"][barcode] = 0
                            library["unmatched"][barcode] += 1
                    else :
                        SeqIO.write(record, fout, "fastq")

                else:

                    SeqIO.write(record, fout, "fastq")


    else :

        match = re.search('CGCCG[ATGC]{7}GCC[ATGC]{9}CGCCG', str(record.seq))

        if match:
            barcode = str(record.seq[match.start(0) - 19:match.start(0) + 48])
            genotyping = str(record.seq[match.start(0) - 39:match.start(0) - 19])
            multiplex = str(record.seq[match.start(0) - 43:match.start(0) - 39])
            multiplex6bp = str(record.seq[match.start(0) - 45:match.start(0) - 39])

            if multiplex in bc:
                bcStats[multiplex] += 1
                if not barcode in library[multiplex]:
                    library[multiplex][barcode] = 0
                library[multiplex][barcode] += 1
            elif multiplex6bp in bc:
                bcStats[multiplex6bp] += 1
                if not barcode in library[multiplex6bp]:
                    library[multiplex6bp][barcode] = 0
                library[multiplex6bp][barcode] += 1
            else:
                bcStats["unmatched"] += 1
                if not barcode in library["unmatched"]:
                    library["unmatched"][barcode] = 0
                library["unmatched"][barcode] += 1

        else :
            # Emtpy vector
            match = re.search('AGAGACGGATATCACTAGTCGTCTCCGTTCGCTCTAGACAGGGTACCCAGCATATGATAGGGTCCCCT', str(record.seq))

            if match:
                barcode = str(record.seq[match.start(0):match.start(0) + 67])
                genotyping = str(record.seq[match.start(0) - 20:match.start(0)])
                multiplex = str(record.seq[match.start(0) - 24:match.start(0) - 20])
                multiplex6bp = str(record.seq[match.start(0) - 26:match.start(0) - 20])

                if multiplex in bc:
                    bcStats[multiplex] += 1
                    if not barcode in library[multiplex]:
                        library[multiplex][barcode] = 0
                    library[multiplex][barcode] += 1
                elif multiplex6bp in bc:
                    bcStats[multiplex6bp] += 1
                    if not barcode in library[multiplex6bp]:
                        library[multiplex6bp][barcode] = 0
                    library[multiplex6bp][barcode] += 1
                else:
                    bcStats["unmatched"] += 1
                    if not barcode in library["unmatched"]:
                        library["unmatched"][barcode] = 0
                    library["unmatched"][barcode] += 1
            else :

                if args.spikein:
                # Spike-in barcode
                #  CCTAAAGCTTCTCCTGCCG GTGTGTGGAACGAGCACAGCgccgAGAGACGGATATCACTAGTCgccgCCATTTGCGCGCGCTCGCC
                # 4mer index 20nt genotyping GTGTGTGGAACGAGCACAGCgccgAGAGACGGATATCACTAGTCgccgCCATTTGCGCGCGCTCGCC

                    match = re.search('GTGTGTGGAACGAGCACAGCgccgAGAGACGGATATCACTAGTCgccgCCATTTGCGCGCGCTCGCC'.upper(), str(record.seq))

                    if match:

                        barcode = str(record.seq[match.start(0):match.start(0) + 67])
                        genotyping = str(record.seq[match.start(0) - 20:match.start(0)])
                        multiplex = str(record.seq[match.start(0) - 24:match.start(0) - 20])
                        multiplex6bp = str(record.seq[match.start(0) - 26:match.start(0) - 20])

                        if multiplex in bc:
                            bcStats[multiplex] += 1
                            if not barcode in library[multiplex]:
                                library[multiplex][barcode] = 0
                            library[multiplex][barcode] += 1
                        elif multiplex6bp in bc:
                            bcStats[multiplex6bp] += 1
                            if not barcode in library[multiplex6bp]:
                                library[multiplex6bp][barcode] = 0
                            library[multiplex6bp][barcode] += 1
                        else:
                            bcStats["unmatched"] += 1
                            if not barcode in library["unmatched"]:
                                library["unmatched"][barcode] = 0
                            library["unmatched"][barcode] += 1
                    else :
                        SeqIO.write(record, fout, "fastq")

                else:

                    SeqIO.write(record, fout, "fastq")

fout.close()

multiplexes = library.keys()
multiplexes.sort()

for multiplex in multiplexes:

    fout = open(os.path.join(args.outdir, bc[multiplex] + "_barcode_counts.txt"), "w")

    barcodes = library[multiplex].keys()
    barcodes.sort()
    for barcode in barcodes:
        print(barcode + "\t" + str(library[multiplex][barcode]), file = fout)

    fout.close()

multiplexes = bc.keys()
multiplexes.sort()

for multiplex in multiplexes:
    print(bc[multiplex] + "\t" + str(bcStats[multiplex]))
