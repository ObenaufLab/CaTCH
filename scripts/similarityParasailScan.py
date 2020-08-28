#!/usr/bin/env python

from __future__ import print_function

import logging
import sys, pysam

from argparse import ArgumentParser, RawDescriptionHelpFormatter

try:
    import parasail
except ImportError as e:
    logging.error("Could not load parasail library. Please try reinstalling "
                  "parasail using pip.")
    sys.exit(1)

import operator


if parasail.can_use_sse2():
    parasail_sg_stat = parasail.sg_qx_stats_striped_32
    parasail_sg = parasail.sg_qx_striped_32
else:
    logging.warning("Warning: SSE not supported, falling back to standard alignment")
    parasail_sg_stat = parasail.sg_qx_stats
    parasail_sg = parasail.sg_qx


usage = "similarityParasailScan.py"
parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-b", "--barcodes", type=str, required=True, dest="barcodesFile", help="Tab delimited barcodes file: name sequence")
parser.add_argument("-r", "--reads", type=str, required=True, dest="readsFile", help="Bam file of reads")
parser.add_argument("-c", "--cutoff", type=float, required=False, default=0.8, dest="scoreCutoff", help="Minimum score cutoff")
parser.add_argument('-i', "--identity", action='store_true', dest="identity", help="Use identity instead of score (default: false)")

args = parser.parse_args()

barcodes = dict()
scoreMax = dict()

with open(args.barcodesFile, 'r') as f:
    for line in f:
        id, seq = line.rstrip().split('\t')
        barcodes[id] = seq.upper()

print("\t".join(barcodes.keys()))

samfile = pysam.AlignmentFile(args.readsFile, "rb", check_sq=False)

initScores = True

counter = 0

for read in samfile.fetch(until_eof=True):

    # Get the maximum achievable alignment score per barcode for reference
    # Randomly get first read, substitute substring with barcode
    # does not matter where for semi-global alignment
    # record score
    if initScores:
        initScores = False
        for barcode in barcodes:

            modRead = read.query_sequence
            modRead = ''.join([modRead[0:3], barcodes[barcode], modRead[3+len(barcodes[barcode]):]])

            alignment = parasail_sg_stat(s1=modRead,
                             s2=barcodes[barcode],
                             open=1,
                             extend=1,
                             matrix=parasail.pam100
            )
            scoreMax[barcode] = alignment.score

    scoreCollection = list()

    report = False

    for barcode in barcodes:
        alignment = parasail_sg_stat(s1=read.query_sequence,
                         s2=barcodes[barcode],
                         open=1,
                         extend=1,
                         matrix=parasail.pam100
        )
        score = alignment.score / scoreMax[barcode]

        identity = float(alignment.matches) / len(barcodes[barcode])

        if args.identity:
            if identity > args.scoreCutoff:
                report = True

            scoreCollection.append(str(identity))

        else:
            if score > args.scoreCutoff:
                report = True

            scoreCollection.append(str(score))

    if report:
        print("\t".join(scoreCollection))

    counter += 1

    if counter % 1000000 == 0:
        print(str(counter) + " reads processed.",file=sys.stderr)
