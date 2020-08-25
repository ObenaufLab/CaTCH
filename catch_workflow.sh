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
    echo "      $0 -b BAM_DIR -o PROCESS_DIR -O REPORT_DIR -v COVARS [-A ABUND_THRESH] [-R REF_ROWS] [-c DEMUX_BC] [-m HAMMING_DIST] [-X PATH_TO_CATCH_SCRIPTS] [-r] [-i] [-s]"
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
plot=1
abund=0.01
ref=1
revcomp=0
spikedin=0
stringent=0
# Parse options.
while getopts 'b:c:o:O:v:X:m:n:A:R:ris1234' flag; do
  case "${flag}" in
    b) bam="${OPTARG}" ;;         # BAM folder.
    c) barcodes="${OPTARG}" ;;    # Demultiplexing table (for all lengths of sample tags)
    o) outdir="${OPTARG}" ;;      # Output directory for the counts
    O) resdir="${OPTARG}" ;;      # Ouptut directory for the analysis report
    v) covars="${OPTARG}" ;;      # List of samples in desired order, corresponding condition, desired corresponding colour
    X) SCRIPTSPATH="${OPTARG}" ;; # Where to find all the scripts for this workflow
    m) hamm="${OPTARG}" ;;        # Hamming distance at which to merge barcodes as likely sequencing errors
    n) dark="${OPTARG}" ;;        # Number of dark bases to allow in the pattern (0)
    A) abund="${OPTARG}" ;;       # Barcode abundance threshold for the report (0.01).
    R) ref="${OPTARG}" ;;         # Comma seperated list of row numbers in covars to be used as reference samples in report, NOT counting the header line (1).
    r) revcomp=1 ;;               # Reverse complement the barcodes
    i) spikedin=1 ;;              # Spike-in barcode was added (hard-coded barcode sequence)
    s) stringent=1 ;;             # Match full format of semi-random barcodes
    1) dxcnt=0 ;;                 # Skip demux and count
    2) dist=0 ;;                  # Skip hammind distance merge
    3) post=0 ;;                  # Skip table mergers
    4) plot=0 ;;                  # Skip analysis report.
    *) usage ;;
  esac
done

if [[ -z "$outdir" ]]; then
  outdir='./process'
fi

if [[ ! -d $outdir ]]; then
  mkdir -p $outdir
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
    if [[ ! -z "$barcodes" ]]; then
        barcodes=",-b ${barcodes}"
    fi

    if [ "$revcomp" -eq 1 ]; then
      revcomp=',-r'
    else
      revcomp=''
    fi

    if [ "$spikedin" -eq 1 ]; then
      spikedin=',-i'
    else
      spikedin=''
    fi

    if [ "$stringent" -eq 1 ]; then
      stringent=',-s'
    else
      stringent=''
    fi

    # Match semi-random barcode format, demultiplex, and count the barcodes.
    echo "${prefix}: Demultiplexing and counting barcodes"
    ${SCRIPTSPATH}/fileutilities.py T $bam --dir 'bam$' | ${SCRIPTSPATH}/fileutilities.py P --loop sbatch ,-J CaTCHc ,--mem=50G ,-o /dev/null ,-e /dev/null ${SCRIPTSPATH}/barcodingQuantifier.py ,-f {abs} ,-o $outdir ,--n_dark $dark $barcodes $revcomp $spikedin $stringent
    wait_for_jobs CaTCHc
fi

if [ "$dist" -eq 1 ]; then
  # Merge barcodes with too similar sequences
  if [ "$hamm" -gt 0 ]; then
      echo "${prefix}: Merging similar barcodes at distance $hamm"
      ## preserve unmerged counts in a new name, so the merged can have the name that is used downstream
      ${SCRIPTSPATH}/fileutilities.py T ${outdir}/*/ --dir barcode-counts.txt | ${SCRIPTSPATH}/fileutilities.py P --loop mv {abs} {dir}/{bas}.raw.txt
      ${SCRIPTSPATH}/fileutilities.py T ${outdir}/*/ --dir barcode-counts.raw.txt | ${SCRIPTSPATH}/fileutilities.py P --loop sbatch ,-o /dev/null ,-e /dev/null ,-J CaTCHh ,--qos=medium ,--mem=20G ${SCRIPTSPATH}/barcodingHammingMerge.py ,-b {abs} ,-d $hamm ,-o {dir}/{cor}.txt
      wait_for_jobs CaTCHh
  fi
fi

if [ "$post" -eq 1 ]; then
  echo "${prefix}: Merging sample-specific count tables into one"
  # Merge sample-specific count tables into one collective table
  ${SCRIPTSPATH}/mergeBCcounts.R ${outdir}/${prefix}_barcode-counts.txt ${outdir}/*/*_barcode-counts.txt

  # Also merge the summary reports.
  ${SCRIPTSPATH}/mergeBCcounts.R ${outdir}/${prefix}_summaries.txt ${outdir}/*/*_summary.txt > /dev/null   # The stdout merge report is not sensible for the summaries file.
fi

# # Note to self how to merge the hammdist tables.
# head -n 2 ${outdir}/${prefix}/*hammdist.txt | tail -n 1 > ${outdir}/${prefix}_hammdist.tsv
# tail -n +3 ${outdir}/${prefix}/*hammdist.txt | perl -e 'while($line=<STDIN>){print $line if $line !~/^[=\n]/}' >> ${outdir}/${prefix}_hammdist.tsv


if [ "$plot" -eq 1 ]; then
  echo "${prefix}: Compiling barcoding report"
  srun --mem=9G ${SCRIPTSPATH}/barcoding_results_run.R -c ${outdir}/${prefix}_barcode-counts.txt -s ${outdir}/${prefix}_summaries.txt -d ${resdir}/${prefix} -v $covars -r $ref -A $abund
fi



echo "${prefix}: Workflow finished!"
