#!/usr/bin/env python3

## version 0.8.0


# Demultiplex samples and count the semi-random barcode occurrences in each.

# Date located in: -

# from __future__ import print_function

import sys, os, re
# import subprocess
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from collections import Counter

import pysam
from Bio.Seq import Seq
import pandas as pd

usage = "Barcoding"
parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)

parser.add_argument("-f", "--file", type=str, required=True, dest="bamFile", help="A single BAM file.")
parser.add_argument("-b", "--barcodesFile", type=str, required=False, help="Optional tab-delimited demultiplexing table.")
parser.add_argument('-o', "--outdir", type=str, required=False, default="./process", help="Output directory (./process/).")

parser.add_argument('-r', "--revcomp", action='store_true', help="Reverse complement the sample tags (False).")
parser.add_argument('--bc_len', type=int, default=68, help="Length of the barcode (68).")
parser.add_argument('--gt_len', type=int, default=20, help="Number of bases between the end of the sample tag and the start of the barcode (20).")

parser.add_argument('-i', "--spikein", action='store_true', help="Also match the hard-coded spike-in barcode pattern.")
parser.add_argument('-s', "--stringent", action='store_true', help="Use the hard-coded long barcode pattern instead of the short one (False).")
parser.add_argument('--n_dark', type=int, default=0, help="Number of dark bases to allow in the barcode (0). For this to work, the barcode pattern must include N's for all the positions that are allowed to be dark.")
parser.add_argument("--pattern", type=str, required=False, help="Optional custom regex string to override the hard-coded design. It must include the start of the barcode.")

args = parser.parse_args()

prefix = os.path.basename(args.bamFile).replace('.bam', '')

##############
# Read demultiplexing tags
##############

bc = dict()       # dictionary of sample tags
tagLens = None     # set of sample tag lengths. When not none, demultiplexing will be applied.
if args.barcodesFile:
    with(open(args.barcodesFile, 'r')) as f:
        next(f) # skip header
        forbidden = re.compile('(^[0-9])|([^0-9a-zA-z_])')
        problems = False
        for line in f:
            try:
                multiplex, sample, treatgroup, treatment, colour = line.rstrip().split("\t")   # combined definitions table
            except ValueError:
                multiplex, sample = line.rstrip().split("\t")                      # demultiplex-only table
            # Sanitize
            match = forbidden.search(sample)
            if match:
                problems = True
                sys.stderr.write('Invalid sample name detected: ' + sample + ' .\n')
            if args.revcomp:
                multiplex = str(Seq(multiplex).reverse_complement())
            bc[multiplex] = sample
        if problems:
            sys.stderr.write('Sample names must: [1] start with a letter, [2] contain no spaces or symbols (except underscore).\n')
            sys.exit(0)
    tagLens = {len(x) for x in bc}
else:
    bc["demuxed"] = prefix

##############
# Barcode designs
##############

# Christian's design:
#     P7adapter_template_BARCODE_GenotypeTag_SampleINDEX_NNNNNN_P5adapter       <--- read direction. Sample index is contained in the read and can be used to demultiplex.
# Shona's design:
#     P7adapter_SampleIndex_template_BARCODE_GenotypeTag_P5adapter              <--- read direction, unlikely to contain the Sample index,
#                                                                                    reads must be demultiplexed in advance using the mate read.

#    20nt                |short end       short start|               19nt
# <--GBNSNNNVDNVNVWVMWNNRCGGCGBNSNNNNDNGGCWVMWNNRCGGCGBNSNNNVDNVNVWVMWNNR--<      <- sequencing direction <-
#
# >--RNNWMVWVNVNDVNNNSNBCGCCGRNNWMVWGCCNDNNNNSNBCGCCGRNNWMVWVNVNDVNNNSNBG-->      strict pattern
#    19nt               CGCCGRNNWMVWGCCNDNNNNSNBCGCCG                20nt         short pattern and offsets


# Umkehrer/Obenauf barcode pattern
bc_pattern_full = re.compile('[TC]{1}[ATGC]{2}[AT]{1}[TG]{1}[TGC]{1}[AT]{1}[TGC]{1}[ATGC]{1}[TGC]{1}[ATGC]{1}[ATC]{1}[TGC]{1}[ATGC]{3}[CG]{1}[ATGC]{1}[CGA]{1}CGCCG[TC]{1}[AGTC]{2}[AT]{1}[TG]{1}[TCG]{1}[AT]{1}GCC[ATGC]{1}[ACT]{1}[ATGC]{4}[GC]{1}[AGTC]{1}[CGA]{1}CGCCG')
bc_pattern_short = re.compile('CGCCG[ATGC]{7}GCC[ATGC]{9}CGCCG')
bc_pattern_dark = re.compile('[CN][GN][CN][CN][GN][ATGCN]{7}[GN][CN][CN][ATGCN]{9}[CN][GN][CN][CN][GN]')  # allow dark bases, based on the short pattern.
                                                                                                         # should be ok for high quality sequencing with very few Ns.
                                                                                                         # would be problematic if many Ns are present.
# Empty vector pattern
empty_pattern = re.compile('AGAGACGGATATCACTAGTCGTCTCCGTTCGCTCTAGACAGGGTACCCAGCATATGATAGGGTCCCCT')
# Spike-in barcode
#            CCTAAAGCTTCTCCTGCCG GTGTGTGGAACGAGCACAGCgccgAGAGACGGATATCACTAGTCgccgCCATTTGCGCGCGCTCGCC
# 4mer index     20nt genotyping GTGTGTGGAACGAGCACAGCgccgAGAGACGGATATCACTAGTCgccgCCATTTGCGCGCGCTCGCC
spike_pattern = re.compile('GTGTGTGGAACGAGCACAGCgccgAGAGACGGATATCACTAGTCgccgCCATTTGCGCGCGCTCGCC'.upper())

bc_pattern = None
offset = 0
if args.pattern:    # External pattern provided. Invalidates all parameters that apply only to the Umkehrer/Obenauf design.
    args.stringent = False
    args.spikein = False
    empty_pattern = None
    bc_pattern = re.compile(args.pattern)
elif args.stringent:
    bc_pattern = bc_pattern_full
elif args.n_dark > 0:
    bc_pattern = bc_pattern_dark
    offset = 19    # Compensate for the short pattern starting 20nt internally to the barcode.
else:
    bc_pattern = bc_pattern_short
    offset = 19    # Compensate for the short pattern starting 20nt internally to the barcode.




##############
# Initialize counters
##############

# Global stats
bcStats = Counter()  # Counter of total matches per sample.

# Initialize a counter for each sample
library = dict()                    # dictionary of counters. One barcode counter per sample.
for sample in bc.values():
    library[sample] = Counter()
bc["unknown"] = "SampleUnknown"       # unrecognized sample tags
bc["unmatched"] = "BCUnmatched"       # barcode pattern didn't match
bc["empty"] = "EmptyVector"           # Diagnostic 1
bc["spike"] = "SpikeIn"               # Diagnostic 2
unmatchedTags = Counter()

##############
# Count the barcodes
##############

# Reduce repetitive code
# This function directly updates the counters in the global scope.
# Reads that don't match the pattern return False (BCUnmatched),
# reads that match the pattern but do not match the demultiplication codes go to fout (SampleUnknown) and return True,
# reads that match the pattern and a demux code go to the library dict and return True.
# The diagnostic parameter is used to redirect reads that match special patterns (empty vector, spike-in), so they don't interfere with the library.
def recognize(mypattern, record, fout, bcl = args.bc_len, gtl = args.gt_len, sampleTagLens = tagLens, extra = 0, diagnostic = False, dark = 0):
    myread = record.query_sequence.upper()
    match = mypattern.search(myread)
    if match and not diagnostic:
        start = match.start(0) - extra               # The pattern is internal to the barcode.
        barcode = myread[start:(start + bcl)]
        if barcode.count('N') > dark:
            return True                     # Reject the match. Returning True means I don't need to try the diagnostic patterns.

        if sampleTagLens:                   # If there is a list of tag lengths, there is a demux table, and we gotta extract the sample tag sequence.
            found = False               # flag to report whether any of the tag lengths resulted in a match
            sampleTag = None
            for k in sampleTagLens:                   # Allows tags of mixed length
                sampleTag = myread[(start - gtl - k):(start - gtl)].upper()  # The sample tag ends at a fixed offset from the barcode.
                if sampleTag in bc:
                    found = True
                    if not diagnostic:
                        bcStats.update([ bc[sampleTag] ])             # total hits for the sample
                        library[ bc[sampleTag] ].update([ barcode ])    # hits of the barode in the sample
                    break
            if not found:
                fout.write(record)
                bcStats.update([ bc["unknown"] ])
                # Since the sample tag matches none of the known ones, what is it actually?
                unmatchedTags.update([ myread[(start - gtl - min(sampleTagLens)):(start - gtl)].upper() ])  # In a mixed tag length scenario, always use the shortest length,
                                                                                                            # to avoid going into random sequence and crowding the output.
        else:
            bcStats.update([ bc["demuxed"] ])           # Dummy sample to maintain the code logic and output format. Resolves to the file prefix.
            library[ bc["demuxed"] ].update([ barcode ])
        return True
    elif match:
        return True
    else:
        return False

# Parse BAM
if not os.path.exists(os.path.join(args.outdir, prefix)):
    os.mkdir(os.path.join(args.outdir, prefix))
with pysam.AlignmentFile(args.bamFile, 'rb', check_sq=False) as fin:      # Checking for reference sequences in the header has to be disabled for unaligned BAM.
    with pysam.AlignmentFile(os.path.join(args.outdir, prefix, prefix + "_unmatched.bam"), 'wb', template=fin) as unmatched_bc:
        with pysam.AlignmentFile(os.path.join(args.outdir, prefix, prefix + "_unknown.bam"), 'wb', template=fin) as unknown_sample:
            with pysam.AlignmentFile(os.path.join(args.outdir, prefix, prefix + "_diagnostic.bam"), 'wb', template=fin) as diagnostic_pattern:
                # Scan the reads
                for record in fin:
                    if recognize(bc_pattern, record, unknown_sample, extra=offset, dark=args.n_dark):
                        pass  # Counters and barcode library are updated by the function directly, there is no additional logic to implement.
                    elif args.spikein and recognize(spike_pattern, record, unknown_sample, diagnostic=True):
                        bcStats.update([ bc['spike'] ])
                        diagnostic_pattern.write(record)
                    elif (empty_pattern is not None) and recognize(empty_pattern, record, unknown_sample, diagnostic=True) :
                        bcStats.update([ bc['empty'] ])
                        diagnostic_pattern.write(record)
                    else :
                        # None of the patterns has matched
                        bcStats.update([bc["unmatched"]])
                        unmatched_bc.write(record)

##############
# Output the barcode counts for each sample
##############

samples = sorted(library.keys())
for sample in samples:
    with open(os.path.join(args.outdir, prefix, sample.replace('/', '.') + "_barcode-counts.txt"), "w") as fout:
        for v,c in library[sample].most_common():
            fout.write( v + "\t" + str(c) + "\n" )

###############
#  Output unmatched sample tag tallies, for troubleshooting
###############
with open(os.path.join(args.outdir, prefix, prefix + "_rogueTags.txt"), "w") as fout:
    for v,c in unmatchedTags.most_common():
        fout.write( v + "\t" + str(c) + "\n" )

##############
# Demultiplexing report
##############

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

samples = natural_sorted(bcStats.keys())
with open(os.path.join(args.outdir, prefix, prefix + "_summary.txt"), "w") as fout:
    for sample in samples:
        fout.write(sample + "\t" + str(bcStats[sample]) + "\n")
        sys.stdout.write(sample + "\t" + str(bcStats[sample]) + "\n")


assigned = pd.DataFrame(bcStats.most_common())
toprogue = pd.DataFrame(unmatchedTags.most_common()).head()
notdiagnostic = (assigned.iloc[:, 0] != 'BCUnmatched') & (assigned.iloc[:, 0] != 'SampleUnknown') & (assigned.iloc[:, 0] != 'EmptyVector') & (assigned.iloc[:, 0] != 'SpikeIn')
if len(unmatchedTags) > 0:
    if toprogue.iloc[0,1] > min(assigned[notdiagnostic].iloc[:,1]):
        sys.stdout.write("\nIt seems there are unassigned sample tags that have more reads than some of your samples. If this is unexpected, maybe the provided sample tags have typos in them or need to be reverse-complemented.\n")
        sys.stdout.write(str(toprogue))
        sys.stdout.write("\n...\n")
