#!/usr/bin/env sh
#
#SBATCH --get-user-env
#SBATCH -J CaTCH
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50000
#SBATCH --output=catch.out
#SBATCH --error=catch.err

function usage() {
    echo "Usage:"
    echo "      $0 -b BAM_FILE [[-c DEMUX_BC4nt | -6 DEMUX_BC6nt]] [-o OUTDIR] [-X PATH_TO_CATCH_SCRIPTS] [-K PATH/FILEUTILITIES.PY] [-r] [-i] [-s] [-m]"
    echo "Either or both of -c and -6 must be provided."
    exit 1
}

# Defaults
SCRIPTSPATH='/users/kimon.froussios/catch/scripts'
#FU='/users/kimon.froussios/utility_scripts'
FU=$SCRIPTSPATH
PYTHONPATH="$PYTHONPATH:$FU"
# Parse options.
while getopts 'b:c:6:o:u:X:K:rism' flag; do
  case "${flag}" in
    b) bam="${OPTARG}" ;;         # BAM file
    c) barcodes="${OPTARG}" ;;    # Demultiplexing table for 4nt barcodes
    6) barcodes2="${OPTARG}" ;;   # Demultiplexing table for 6nt barcodes
    o) outdir="${OPTARG}" ;;      # Output directory for the counts
    X) SCRIPTSPATH="${OPTARG}" ;; # Where to find all the scripts for this workflow
    K) FU="${OPTARG}" ;;          # Where to find my fileutilities.py and mylogs.py
    r) revcomp="${OPTARG}" ;;     # Reverse complement the barcodes
    i) spikedin="${OPTARG}" ;;    # Spike-in barcode was added (hard-coded barcode sequence)
    s) stringent="${OPTARG}" ;;   # Match full format of semi-random barcodes
    m) merge="${OPTARGS}" ;;      # Allow for sequencing errors in barcodes by merging very similar barcodes when one is common and the others are rare.
    *) usage ;;
  esac
done

if [[ -z "$outdir" ]]; then
  outdir='./process'
fi

if [[ ! -d $outdir ]]; then
  mkdir -p $outdir
fi

if [[ -z "$unknown" ]]; then
  unknown='./unknown.fastq'
fi

if [ $revcomp ]; then
  revcomp='-r'
else
  revcomp=''
fi

if [ $spikedin ]; then
  spikedin='-i'
else
  spikedin=''
fi

if [ $stringent ]; then
  stringent='-s'
else
  stringent=''
fi

set -e

# Modules
echo "${bam}: Loading modules"
module load samtools/1.9-foss-2017a
module load python/2.7.13-foss-2017a
module load biopython/1.70-foss-2017a-python-2.7.13
module load python-levenshtein/0.12.0-foss-2017a-python-2.7.13

prefix=$(basename $bam)
prefix=${prefix/.bam/}

# Extract to FASTQ
echo "${bam}: Extracting to FASTQ"
fastq=${bam/.bam/.fastq}
samtools fastq $bam -0 $fastq

# Match semi-random barcode format, demultiplex, and count the barcodes.
echo "${fastq}: Demultiplexing and counting barcodes"
if [[ -z "$barcodes2" ]]; then
  ${SCRIPTSPATH}/barcodingQuantifier.py -f $fastq -o $outdir -b $barcodes $revcomp $spikedin $stringent
else
  ${SCRIPTSPATH}/barcodingQuantifier.py -f $fastq -o $outdir -b $barcodes -6 $barcodes2 $revcomp $spikedin $stringent
fi
unknown="${outdir}/unmatched.fastq"
umcounts="${outdir}/unmatched_barcode_counts.txt"

# Save space.
echo "${prefix}: Compressing unmatched reads"
gzip $unknown
# echo "${fastq}: Re-compressing reads"
# gzip $fastq
rm $fastq

# Sort unmatched barcodes by count
echo "${prefix}: Sorting unmatched barcodes by count"
${SCRIPTSPATH}/sort_unmatched.R $umcounts
rm $umcounts      # Prevents it from being pulled into the merge in the next step. The sorted version has a different suffix and won't be pulled in.

# Merge sample-specific count tables into one collective table
echo "${prefix}: Merging sample-specific count tables into one"
${FU}/fileutilities.py T $outdir --dir _barcode_counts.txt | perl -e 'while(<>){~s/_barcode_counts$//;print}' | ${FU}/fileutilities.py P -i --merge outer no yes > ${outdir}/${prefix}_tmp.tsv
head -n 1 ${outdir}/${prefix}_tmp.tsv | perl -e 'while(<>){~s/_\|\d+//g; print}' > ${outdir}/${prefix}_barcode_counts.tsv   # Clean up the header
tail -n +2 ${outdir}/${prefix}_tmp.tsv >> ${outdir}/${prefix}_barcode_counts.tsv
rm ${outdir}/${prefix}_tmp.tsv
${FU}/fileutilities.py T $outdir --dir _barcode_counts.txt | ${FU}/fileutilities.py P --loop rm {abs}
# Rename barcode column
head -n 1 ${outdir}/${prefix}_barcode_counts.tsv | perl -e 'while(<>){~s/row_ID/barcode/;print}' > ${outdir}/tmp
tail -n +2 ${outdir}/${prefix}_barcode_counts.tsv >> ${outdir}/tmp
mv ${outdir}/tmp ${outdir}/${prefix}_barcode_counts.tsv

echo "${prefix}: Workflow finished!"
