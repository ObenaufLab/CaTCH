#!/usr/bin/env python3

# Demultiplex samples and count the semi-random barcode occurrences in each.

# Date located in: -

# from __future__ import print_function

import sys, os, re
# import subprocess
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from collections import Counter

import pysam
from Bio.Seq import Seq

usage = "Barcoding"
parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-f", "--file", type=str, required=True, dest="bamFile", help="BAM file")
parser.add_argument("-b", "--barcodes", type=str, required=False, dest="barcodesFile", help="Demultiplexing tags table, regardless of tag length. (tab delimited: Barcode \\t Sample, with header line). Omit if data already demultiplexed.")
parser.add_argument('-r', "--revcomp", action='store_true', help="Reverse complement barcodes (default: false)")
parser.add_argument('-i', "--spikein", action='store_true', help="Barcode was spiked-in (default: false)")
parser.add_argument('-s', "--stringent", action='store_true', help="Stringent barcode matching (default: false)")
parser.add_argument('-o', "--outdir", type=str, required=False, default="./process", help="Output directory for counts (./process/)")
parser.add_argument('--bc_len', type=int, default=67, help="Length of the barcode from start position of the match (67).")
parser.add_argument('--gt_len', type=int, default=20, help="Genotyping tag length (20).")

args = parser.parse_args()

##############
# Read demultiplexing tags
##############

bc = dict()       # dictionary of sample tags
tagLens = None     # set of sample tag lengths
if args.barcodesFile:
    with(open(args.barcodesFile, 'r')) as f:
        # skip header
        next(f)
        for line in f:
            barcode, sample = line.rstrip().split("\t")
            if args.revcomp:
                barcode = str(Seq(barcode).reverse_complement())
            bc[barcode] = sample
    tagLens = {len(x) for x in bc}
else:
    bc["demuxed"] = os.path.basename(args.bamFile)

##############
# Barcode designs
##############

    # GBNSNNNVDNVNVWVMWNNRCGGCGBNSNNNNDNGGCWVMWNNRCGGCGBNSNNNVDNVNVWVMWNNR      <- sequencing direction
    # RNNWMVWVNVNDVNNNSNB CGCCG RNNWMVW GCC NDNNNNSNB CGCCG
# Christian's design:
    # P7adapter_template_BARCODE_GenotypeTag_SampleINDEX_NNNNNN_P5adapter       <--- read direction, unlikely to get past the template.
# Shona's design:
    # P7adapter_SampleIndex_template_BARCODE_GenotypeTag_P5adapter              <--- read direction, unlikely to contain the SampleIndex, demultiplexed externally by mate pair.

# Barcode pattern
bc_pattern_full = re.compile('[TC]{1}[ATGC]{2}[AT]{1}[TG]{1}[TGC]{1}[AT]{1}[TGC]{1}[ATGC]{1}[TGC]{1}[ATGC]{1}[ATC]{1}[TGC]{1}[ATGC]{3}[CG]{1}[ATGC]{1}[CGA]{1}CGCCG[TC]{1}[AGTC]{2}[AT]{1}[TG]{1}[TCG]{1}[AT]{1}GCC[ATGC]{1}[ACT]{1}[ATGC]{4}[GC]{1}[AGTC]{1}[CGA]{1}CGCCG')
bc_pattern_short = re.compile('CGCCG[ATGC]{7}GCC[ATGC]{9}CGCCG')
bc_pattern = None
offset = None
if args.stringent:
    bc_pattern = bc_pattern_full
    offset = 0
else:
    bc_pattern = bc_pattern_short
    offset = 19    # Compensate for the short pattern not being at the beginning of the strict pattern.
# Empty vector pattern
empty_pattern = re.compile('AGAGACGGATATCACTAGTCGTCTCCGTTCGCTCTAGACAGGGTACCCAGCATATGATAGGGTCCCCT')
# Spike-in barcode
    #            CCTAAAGCTTCTCCTGCCG GTGTGTGGAACGAGCACAGCgccgAGAGACGGATATCACTAGTCgccgCCATTTGCGCGCGCTCGCC
    # 4mer index     20nt genotyping GTGTGTGGAACGAGCACAGCgccgAGAGACGGATATCACTAGTCgccgCCATTTGCGCGCGCTCGCC
spike_pattern = re.compile('GTGTGTGGAACGAGCACAGCgccgAGAGACGGATATCACTAGTCgccgCCATTTGCGCGCGCTCGCC'.upper())

##############
# Initialize counters
##############

# Global stats
bcStats = Counter(unmatched=0)  # Counter of total matches per sample.

# Initialize a counter for each sample
bc["unmatched"] = "unmatched"       # unrecognized sample tags
library = dict()                    # dictionary of counters. One barcode counter per sample.
for multiplex in bc:                # including unmatched
    library[multiplex] = Counter()

##############
# Count the barcodes
##############

# Reduce repetitive code
# This function directly updates the counters in the global scope.
def recognize(mypattern, myread, samples=bc, bcl = args.bc_len, gtl = args.gt_len, sampleTagLens = tagLens, extra=offset):
    match = mypattern.search(myread)
    if match:
        start = match.start(0) - extra
        barcode = myread[start:(start + bcl)]          # The pattern is part of the barcode.
        genotyping = myread[(start - gtl):start]     # The genotyping tag is immediately before the barcode.

        found = False       # flag to report whether any of the different tag lengths had a match
        sampleTag=None
        if sampleTagLens:
            for k in sampleTagLens:                   # Flexible length.
                sampleTag = myread[(start - gtl - k):(start - gtl)]  # The sample tag is either immediately before the gentype tag, or not present.
                if sampleTag in bc:
                    found = True
                    bcStats.update([sampleTag])             # total hits for the sample
                    library[sampleTag].update([barcode])    # hits of the barode in the sample
                    break
            if not found:
                bcStats.update([bc["unmatched"]])
                library["unmatched"].update([barcode])
        else:
            bcStats.update([bc["demuxed"]])
            library["demuxed"].update([barcode])
        return True
    else:
        return False

with pysam.AlignmentFile(args.bamFile, 'rb', check_sq=False) as fin:      # Checking for reference sequences in the header has to be disabled for unaligned BAM.
    # Open a destination for reads that fail to match the demultiplexing tag
    # fout = open(os.path.join(args.outdir, "unmatched.sam"), "w")
    with pysam.AlignmentFile(os.path.join(args.outdir, "unmatched.bam"), 'wb', template=fin) as fout:
        # Scan the reads
        for record in fin:
            if recognize(bc_pattern, record.query_sequence):
                pass  # Counters are updated by the function directly, there is no additional logic to implement.
            elif recognize(empty_pattern, record.query_sequence) :
                pass
            elif args.spikein and recognize(spike_pattern, record.query_sequence):
                pass
            else :
                fout.write(record)

##############
# Output the barcode counts for each sample
##############

multiplexes = sorted(library.keys())

for multiplex in multiplexes:

    with open(os.path.join(args.outdir, bc[multiplex].replace('/', '.') + "_barcode_counts.txt"), "w") as fout:
        for v,c in library[multiplex].most_common():
            fout.write( v + "\t" + str(c) + "\n" )

##############
# Demultiplexing report
##############

multiplexes = sorted(bc.keys())

for multiplex in multiplexes:
    print(bc[multiplex] + "\t" + str(bcStats[multiplex]))
