#!/usr/bin/env Rscript

library(getopt)

spec = matrix(c(
  'help'            , 'h', 0, "logical",   "Help",
  'countsFile'      , 'c', 1, "character", "Tab-separated table of counts with the collected counts for the samples.",
  'summariesFile'   , 's', 1, "character", "Tab-separated table with the collected summaries for the samples.",
  'hammingFile'     , 'e', 1, "character", "Tab-separated table with the collected hamming distances for the samples.",
  'resultsDir'      , 'd', 1, "character", "Directory in which to save all output.",
  'reportTemplate'  , 'T', 1, "character", "Template Rmd file for CaTCH report.",
  'covars'          , 'v', 1, "character", "File with DEseq2-style covariates table, listing samples in desired order, assigned condition (just one variable), and desired colour (in R-compatible string values, preferably avoid 'magenta' as it is used internally)",
  'topdf'           , 'p', 0, "logical",   "Output figures to pdf",
  'refsamps'        , 'r', 1, "character", "Comma-separated list of integers designating the samples to use as reference abundance (in order they appear in summariesFile)",
  'count_thresh'    , 'N', 1, "integer",   "Count threshold for barcodes to analyse",
  'abund_thresh'    , 'A', 1, "numeric",   "Proportional abundance threshold to consider barcodes top hits",
  'extrabc'         , 'B', 1, "character", "Comma-separated list of barcode IDs to include in addition to those chosen by the analysis."
), byrow=TRUE, ncol=5)

opt = getopt(spec)
# opt <- list(countsFile='/users/kimon.froussios/obenauf/shona/catch_R9739/process/000000000-G5Y7J_1_20200603B_20200604_barcode-counts.tsv', summariesFile='/users/kimon.froussios/obenauf/shona/catch_R9739/process/000000000-G5Y7J_1_20200603B_20200604_summary.tsv', samples='/users/kimon.froussios/obenauf/shona/catch_R9739/description/covars.txt', resultsDir='/users/kimon.froussios/obenauf/shona/catch_R9739/results/', refsamps='1', count_thresh=50, abund_thresh=0.01, extrabc='3149,3825,326,902,3321,2517,1511,1201,3902,263,468')

if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

if ( is.null(opt$reportFile) ) {
  opt$reportFile <- '~/catch/barcoding_results_template.Rmd'
}

if ( is.null(opt$topdf) ) {
  opt$topdf <- FALSE
}

if ( is.null(opt$count_thresh) ) {
  opt$count_thresh <- 50L
}
if ( is.null(opt$abund_thresh) ) {
  opt$abund_thresh <- 0.01
}

if ( is.null(opt$refsamps) ) {
  opt$resamps <- 1
} else {
  opt$refsamps <- as.integer(strsplit(opt$refsamps, ',', fixed=TRUE)[[1]])
}

if (! is.null(opt$extrabc) ) {
  opt$extrabc <- strsplit(opt$extrabc, ',', fixed=TRUE)[[1]]
  opt$extrabc <- opt$extrabc[order(opt$extrabc)]
}

rmarkdown::render(opt$reportFile,
                  output_file = sub('.txt|.tsv', '_report.html', basename(opt$countsFile)),
                  output_dir = opt$resultsDir,
                  params=list(counts = opt$countsFile,
                              summaries = opt$summariesFile,
                              hammdist = opt$hammingFile,
                              outpref = file.path(opt$resultsDir, sub('.txt|.tsv', '', basename(opt$countsFile))),
                              topdf = opt$topdf,
                              samples = opt$covars,
                              refsamps = opt$refsamps,
                              count_thresh = opt$count_thresh,
                              abund_thresh = opt$abund_thresh,
                              extra_bc = opt$extrabc)
)

