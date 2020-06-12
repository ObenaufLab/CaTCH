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
    echo "      $0 -b BAM_FILE [-c DEMUX_BC] [-o OUTDIR] [-X PATH_TO_CATCH_SCRIPTS] [-m HAMMING_DIST] [-r] [-i] [-s]"
    exit 1
}

# Defaults
SCRIPTSPATH='/users/kimon.froussios/catch'
#FU='/users/kimon.froussios/utility_scripts'
PYTHONPATH="$PYTHONPATH:$SCRIPTSPATH"
hamm=0
dxcnt=1
dist=1
post=1
dark=0
# Parse options.
while getopts 'b:c:o:X:m:n:ris123' flag; do
  case "${flag}" in
    b) bam="${OPTARG}" ;;         # BAM file
    c) barcodes="${OPTARG}" ;;    # Demultiplexing table (for all lengths of sample tags)
    o) outdir="${OPTARG}" ;;      # Output directory for the counts
    X) SCRIPTSPATH="${OPTARG}" ;; # Where to find all the scripts for this workflow
    m) hamm="${OPTARG}" ;;        # Hamming distance at which to merge barcodes as likely sequencing errors
    n) dark="${OPTARG}" ;;        # Number of dark bases to allow in the pattern (0)
    r) revcomp="${OPTARG}" ;;     # Reverse complement the barcodes
    i) spikedin="${OPTARG}" ;;    # Spike-in barcode was added (hard-coded barcode sequence)
    s) stringent="${OPTARG}" ;;   # Match full format of semi-random barcodes
    1) dxcnt=0 ;;                 # Skip demux and count
    2) dist=0 ;;                  # Skip hammind distance merge
    3) post=0 ;;                  # Skip table mergers
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

wait_for_jobs(){
  echo "waiting"
  sleep 60  # seconds, give time to the scheduler to put up the task
  sleeptime=120  # ask every 2 mins, for the first 10 mins
  n=1
  while true; do
    if [ $(squeue | grep kimon.fr | grep -c $1) -eq 0 ]; then
      break
    else
      echo "sleep another" $((sleeptime / 60)) "minutes..."
      sleep $sleeptime
    fi
    n=$((n + 1))
    if [ $n -eq 5 ]; then
      sleeptime=300  # if still running after 10 mins, ask every 5 mins
    fi
    if [ $n -eq 10 ]; then
      sleeptime=600  # if still running after 30 mins, ask every 10 mins
    fi
  done
}

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


if [ "$dxcnt" -eq 1 ]; then
  # Match semi-random barcode format, demultiplex, and count the barcodes.
  echo "${prefix}: Demultiplexing and counting barcodes"
  sbatch -J CaTCHc --mem=50G -o /dev/null -e /dev/null ${SCRIPTSPATH}/barcodingQuantifier.py -f $bam -o $outdir --n_dark $dark $barcodes $revcomp $spikedin $stringent
  wait_for_jobs CaTCHc
fi

if [ "$dist" -eq 1 ]; then
  # Merge barcodes with too similar sequences
  if [ "$hamm" -gt 0 ]; then
      echo "${prefix}: Merging similar barcodes at distance $hamm"
      ## preserve unmerged counts in a new name, so the merged can have the name that is used downstream
      ${SCRIPTSPATH}/fileutilities.py T ${outdir}/${prefix} --dir barcode-counts.txt | ${SCRIPTSPATH}/fileutilities.py P --loop mv {abs} {dir}/{bas}.raw.txt
      ${SCRIPTSPATH}/fileutilities.py T ${outdir}/${prefix} --dir barcode-counts.raw.txt | ${SCRIPTSPATH}/fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J CaTCHh ,--qos=medium ,--mem=20G ${SCRIPTSPATH}/barcodingHammingMerge.py ,-b {abs} ,-d $hamm ,-o {dir}/{cor}.txt
      # ${SCRIPTSPATH}/fileutilities.py T ${outdir}/${prefix} --dir barcode-counts.raw.txt | ${SCRIPTSPATH}/fileutilities.py P --loop ${SCRIPTSPATH}/barcodingHammingMerge.py ,-b {abs} ,-d $hamm ,-o {dir}/{cor}.txt
      wait_for_jobs CaTCHh
  fi
fi

if [ "$post" -eq 1 ]; then
  # Merge sample-specific count tables into one collective table
  echo "${prefix}: Merging sample-specific count tables into one"
  ${SCRIPTSPATH}/mergeBCcounts.R ${outdir}/${prefix} ${outdir}/${prefix}_barcode-counts.tsv _barcode-counts.txt

  # Also merge the summary reports.
  ${SCRIPTSPATH}/mergeBCcounts.R ${outdir}/${prefix} ${outdir}/${prefix}_summary.tsv _summary.txt
fi

# # Note to self how to merge the hammdist tables.  
# head -n 2 ${outdir}/${prefix}/*hammdist.txt | tail -n 1 > ${outdir}/${prefix}_hammdist.tsv
# tail -n +3 ${outdir}/${prefix}/*hammdist.txt | perl -e 'while($line=<STDIN>){print $line if $line !~/^[=\n]/}' >> ${outdir}/${prefix}_hammdist.tsv

echo "${prefix}: Workflow finished!"
