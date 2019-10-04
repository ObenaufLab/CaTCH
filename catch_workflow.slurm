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
    echo "      $0 -b BAM_FILE [-c DEMUX_BC] [-o OUTDIR] [-X PATH_TO_CATCH_SCRIPTS] [-r] [-i] [-s]"
    exit 1
}

# Defaults
SCRIPTSPATH='/users/kimon.froussios/catch/scripts'
#FU='/users/kimon.froussios/utility_scripts'
PYTHONPATH="$PYTHONPATH:$SCRIPTSPATH"
# Parse options.
while getopts 'b:c:o:u:X:K:rism' flag; do
  case "${flag}" in
    b) bam="${OPTARG}" ;;         # BAM file
    c) barcodes="${OPTARG}" ;;    # Demultiplexing table (for all lengths af sample tags)
    o) outdir="${OPTARG}" ;;      # Output directory for the counts
    X) SCRIPTSPATH="${OPTARG}" ;; # Where to find all the scripts for this workflow
    r) revcomp="${OPTARG}" ;;     # Reverse complement the barcodes
    i) spikedin="${OPTARG}" ;;    # Spike-in barcode was added (hard-coded barcode sequence)
    s) stringent="${OPTARG}" ;;   # Match full format of semi-random barcodes
    m) merge="${OPTARGS}" ;;      # NOT IMPLEMENTED? Allow for sequencing errors in barcodes by merging very similar barcodes when one is common and the others are rare.
    *) usage ;;
  esac
done

if [[ ! -z "$barcodes" ]]; then
    barcodes="-b ${barcodes}"
fi

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
# echo "${bam}: Loading modules"
# module load samtools/1.9-foss-2017a
# module load python/2.7.13-foss-2017a
# module load biopython/1.70-foss-2017a-python-2.7.13
# module load pysam/0.14.1-foss-2017a-python-2.7.13
# module load python-levenshtein/0.12.0-foss-2017a-python-2.7.13

prefix=$(basename $bam)
prefix=${prefix/.bam/}

# Match semi-random barcode format, demultiplex, and count the barcodes.
echo "${prefix}: Demultiplexing and counting barcodes"
${SCRIPTSPATH}/barcodingQuantifier.py -f $bam -o $outdir $barcodes $revcomp $spikedin $stringent

##
## Merging very similar barcode sequences with vastly different abundances using barcodingHammingMerge.py would go somewhere here in the workflow.
## Currently such barcodes are kept separate and likely rejected when the count threshold is applied in the analysis afterwards.
##

# Merge sample-specific count tables into one collective table
echo "${prefix}: Merging sample-specific count tables into one"
${SCRIPTSPATH}/fileutilities.py T ${outdir}/${prefix} --dir _barcode-counts.txt | perl -e 'while(<>){~s/_barcode-counts$//;print}' | ${SCRIPTSPATH}/fileutilities.py P -i --appnd outer > ${outdir}/${prefix}_tmp.tsv
head -n 1 ${outdir}/${prefix}_tmp.tsv | perl -e 'while(<>){~s/_\|\d+//g; print}' > ${outdir}/${prefix}_barcode-counts.tsv   # Clean up the header
tail -n +2 ${outdir}/${prefix}_tmp.tsv >> ${outdir}/${prefix}_barcode-counts.tsv
${SCRIPTSPATH}/fileutilities.py T $outdir --dir _barcode-counts.txt | ${SCRIPTSPATH}/fileutilities.py P --loop rm {abs}
rm ${outdir}/${prefix}_tmp.tsv

# Rename barcode column
head -n 1 ${outdir}/${prefix}_barcode-counts.tsv | perl -e 'while(<>){~s/row_ID/barcode/;print}' > ${outdir}/tmp
tail -n +2 ${outdir}/${prefix}_barcode-counts.tsv >> ${outdir}/tmp
mv ${outdir}/tmp ${outdir}/${prefix}_barcode-counts.tsv

echo "${prefix}: Workflow finished!"
