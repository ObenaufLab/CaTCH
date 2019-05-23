#!/bin/bash
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
    echo "      $0 -b BAM_FILE [[-c DEMUX_BC4nt | -6 DEMUX_BC6nt]] [-o OUTDIR] [-u UNMATCHED_READS] [-r] [-i] [-s] [-m]"
    echo "Either or both of -c and -6 must be provided."
    exit 1
}

# Parse options.
while getopts 'b:c:6:o:u:rism' flag; do
  case "${flag}" in
    b) bam="${OPTARG}" ;;         # BAM file
    c) barcodes="${OPTARG}" ;;    # Demultiplexing table for 4nt barcodes
    6) barcodes2="${OPTARG}" ;;   # Demultiplexing table for 6nt barcodes
    o) outdir="${OPTARG}" ;;      # Output directory for the counts
    u) unknown="${OPTARG}" ;;     # Output FASTQ for reads that do not match the barcode pattern
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


# Modules
echo "${bam}: Loading modules"
module load samtools/1.9-foss-2017a
module load python/2.7.13-foss-2017a
module load biopython/1.70-foss-2017a-python-2.7.13
#source /users/kimon.froussios/miniconda3/bin/activate py27


# Extract to FASTQ
echo "${bam}: Converting to FASTQ"
fastq=${bam/.bam/.fastq}
samtools fastq $bam -0 $fastq
echo "${bam}: FASTQ done"


# Match semi-random barcode format, demultiplex, and count the barcodes.
echo "${fastq}: Demultiplexing and counting barcodes"
if [[ -z "$barcodes2" ]]; then
  python /users/kimon.froussios/catch/scripts/barcodingQuantifier.py -f $fastq -o $outdir -b $barcodes -u $unknown $revcomp $spikedin $stringent
else
  python /users/kimon.froussios/catch/scripts/barcodingQuantifier.py -f $fastq -o $outdir -b $barcodes -6 $barcodes2 -u $unknown $revcomp $spikedin $stringent
fi
echo "${fastq}: Counts done"

#counts=

# Collapse very similar barcodes
#if [[ ! -z "$merge" ]]; then
#  echo "${bam}: Absorbing barcodes with putative sequencing errors"
#
#
#  python /users/kimon.froussios/catch/scripts/barcodingHammingMerge.py -b $counts
#
#else
#  echo "${bam}: Not trying to account for sequencing errors in the barcodes"
#fi
