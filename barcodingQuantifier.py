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

prefix = os.path.basename(args.bamFile).replace('.bam', '')

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
            multiplex, sample = line.rstrip().split("\t")
            if args.revcomp:
                multiplex = str(Seq(multiplex).reverse_complement())
            bc[multiplex] = sample
    tagLens = {len(x) for x in bc}
else:
    bc["demuxed"] = prefix

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
bcStats = Counter()  # Counter of total matches per sample.

# Initialize a counter for each sample
library = dict()                    # dictionary of counters. One barcode counter per sample.
for sample in bc.values():
    library[sample] = Counter()
bc["unknown"] = "sample-unknown"       # unrecognized sample tags
bc["unmatched"] = "bc-unmatched"       # barcode pattern didn't match

##############
# Count the barcodes
##############

# Reduce repetitive code
# This function directly updates the counters in the global scope.
def recognize(mypattern, record, fout, samples=bc, bcl = args.bc_len, gtl = args.gt_len, sampleTagLens = tagLens, extra=offset):
    myread = record.query_sequence
    match = mypattern.search(myread)
    if match:
        start = match.start(0) - extra
        barcode = myread[start:(start + bcl)]          # The pattern is part of the barcode.
        genotyping = myread[(start - gtl):start]     # The genotyping tag is immediately before the barcode.

        if sampleTagLens:
            found = False       # flag to report whether any of the tag lengths resulted in a match
            sampleTag=None
            for k in sampleTagLens:                   # Flexible length.
                sampleTag = myread[(start - gtl - k):(start - gtl)]  # The sample tag is either immediately before the gentype tag, or not present.
                if sampleTag in bc:
                    found = True
                    bcStats.update([bc[sampleTag]])             # total hits for the sample
                    library[bc[sampleTag]].update([barcode])    # hits of the barode in the sample
                    break
            if not found:
                bcStats.update([bc["unknown"]])
                ##library[bc["unknown"]].update([barcode])
                fout.write(record)
        else:
            bcStats.update([bc["demuxed"]])
            library[bc["demuxed"]].update([barcode])
        return True
    else:
        return False

if not os.path.exists(os.path.join(args.outdir, prefix)):
    os.mkdir(os.path.join(args.outdir, prefix))
with pysam.AlignmentFile(args.bamFile, 'rb', check_sq=False) as fin:      # Checking for reference sequences in the header has to be disabled for unaligned BAM.
    with pysam.AlignmentFile(os.path.join(args.outdir, prefix, prefix + "_unmatched.bam"), 'wb', template=fin) as unmatched_bc:
        with pysam.AlignmentFile(os.path.join(args.outdir, prefix, prefix + "_unknown.bam"), 'wb', template=fin) as unknown_sample:
            # Scan the reads
            for record in fin:
                if recognize(bc_pattern, record, unknown_sample):
                    pass  # Counters are updated by the function directly, there is no additional logic to implement.
                elif recognize(empty_pattern, record, unknown_sample) :
                    pass
                elif args.spikein and recognize(spike_pattern, record, unknown_sample):
                    pass
                else :
                    bcStats.update([bc["unmatched"]])
                    ##library[bc["unmatched"]].update([barcode])
                    unmatched_bc.write(record)

##############
# Output the barcode counts for each sample
##############

samples = sorted(library.keys())
for sample in samples:
    with open(os.path.join(args.outdir, prefix, sample.replace('/', '.')) + "_barcode-counts.txt", "w") as fout:
        for v,c in library[sample].most_common():
            fout.write( v + "\t" + str(c) + "\n" )

##############
# Demultiplexing report
##############

samples = sorted(bcStats.keys())
with open(os.path.join(args.outdir, prefix, prefix + "_summary.txt"), "w") as fout:
    for sample in samples:
        fout.write(sample + "\t" + str(bcStats[sample]) + "\n")
