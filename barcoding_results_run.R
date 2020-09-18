#!/usr/bin/env Rscript

library(getopt)

spec = matrix(c(
  'help'            , 'h', 0, "logical",   "Help",
  'countsFile'      , 'c', 1, "character", "Tab-separated table of counts with the collected counts for the samples. Full path needed.",
  'summariesFile'   , 's', 1, "character", "Tab-separated table with the collected summaries for the samples. Full path needed.",
  'hammingFile'     , 'e', 1, "character", "Tab-separated table with the collected hamming distances for the samples. Full path needed.",
  'resultsDir'      , 'd', 1, "character", "Directory in which to save all output. Full path needed.",
  'reportTemplate'  , 'T', 1, "character", "Template Rmd file for CaTCH report.",
  'covars'          , 'v', 1, "character", "Table file listing sample names ('Sample') (in desired order), respective condition ('Treatment') and desired colour code ('Colour') (in R-compatible names or #RGB. Preferably avoid white as it is used internally). For convenience, an optional field 'Barcode' is allowed as the first column, to permit using one table for the whole pipeline. Full path needed.",
  'topdf'           , 'p', 0, "logical",   "Output figures to pdf",
  'refsamps'        , 'r', 1, "character", "Comma-separated list of integers designating the samples to use as reference abundance (in order they appear in summariesFile)",
  'count_thresh'    , 'N', 1, "integer",   "Count threshold for barcodes to analyse",
  'abund_thresh'    , 'A', 1, "numeric",   "Proportional abundance threshold to consider barcodes top hits",
  'extrabc'         , 'B', 1, "character", "Comma-separated list of barcode IDs to include in addition to those chosen by the analysis."
), byrow=TRUE, ncol=5)

opt = getopt(spec)

# opt <- list(countsFile='/Volumes/groups/obenauf/Kimon_Froussios/catch_test/process/data1_barcode-counts.txt', summariesFile='/Volumes/groups/obenauf/Kimon_Froussios/catch_test/process/data1_summary.txt', covars='/Volumes/groups/obenauf/Kimon_Froussios/catch_test/covars1.txt', resultsDir='/Volumes/groups/obenauf/Kimon_Froussios/catch_test/results', count_thresh=1, topdf=TRUE, reportTemplate='/Volumes/groups/obenauf/Kimon_Froussios/catch/barcoding_results_template.Rmd')

# opt <- list(countsFile='/groups/obenauf/Kimon_Froussios/catch_test/process/data1_barcode-counts.txt', summariesFile='/groups/obenauf/Kimon_Froussios/catch_test/process/data1_summary.txt', covars='/groups/obenauf/Kimon_Froussios/catch_test/covars1.txt', resultsDir='/groups/obenauf/Kimon_Froussios/catch_test/results', count_thresh=1, topdf=TRUE)

# opt <- list(countsFile='/Volumes/groups/obenauf/Kimon_Froussios/shona/catch_R10246/process/data_barcode-counts.txt', summariesFile='/Volumes/groups/obenauf/Kimon_Froussios/shona/catch_R10246/process/data_summaries.txt', covars='/Volumes/groups/obenauf/Kimon_Froussios/shona/catch_R10246/description/covars.txt', resultsDir='/Volumes/groups/obenauf/Kimon_Froussios/shona/catch_R10246/results')

if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

if (any( is.null(opt$countsFile), is.null(opt$summariesFile), is.null(opt$covars), is.null(opt$resultsDir) )) {
  stop('Missing input. Need counts file, summaries file, conditions file and output directory.')
}

if ( is.null(opt$reportTemplate) ) {
  opt$reportTemplate <- '~/catch/barcoding_results_template.Rmd'
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
  opt$refsamps <- 1
} else {
  opt$refsamps <- as.integer(strsplit(opt$refsamps, ',', fixed=TRUE)[[1]])
}

if (! is.null(opt$extrabc) ) {
  opt$extrabc <- strsplit(opt$extrabc, ',', fixed=TRUE)[[1]]
  opt$extrabc <- opt$extrabc[order(opt$extrabc)]
}

# Sanitize
DT <- read.delim(opt$covars, stringsAsFactors = FALSE)

if( (length(DT) == 3 && names(DT) != c('Sample', 'Treatment', 'Colour')) || 
    (length(DT) == 4 && names(DT) != c('Barcode', 'Sample', 'Treatment', 'Colour')) ) {
  stop('The covariates table must have fields named Barcode, Sample, Treatment, Colour. In that order. The Barcode column is optional, the others are required.')
}

bad <- grepl('(^[0-9])|([^0-9a-zA-z_])', DT$Sample, perl=TRUE)
if (any(bad)) {
  print(DT[bad, 'Sample'])
  stop('Invalid sample name(s) detected. Names must start with a letter, and contain no spaces and no symbols (except underscore).')
}
bad <- grepl('(^[0-9])|([^0-9a-zA-z_])', DT$Treatment, perl=TRUE)
if (any(bad)) {
  print(DT[bad, 'Treatment'])
  stop('Invalid treatment name(s) detected. Names must start with a letter, and contain no spaces and no symbols (except underscore).')
}

areColors <- function(x) {
  vapply(x, function(X) {
    tryCatch(is.matrix(col2rgb(X)), 
             error = function(e) FALSE)
  }, logical(1))
}
bad <- !areColors(DT$Colour)
if (any(bad)) {
  print(DT[bad, 'Colour'])
  stop("Invalid colour value(s). Choose among R's named colours.")
}

bad <- ! DT$Sample %in% names(read.delim(opt$countsFile))
if (any(bad)) {
  print(names(DT)[bad])
  stop("Some of the sample names in the covariates table were not found in the counts table.")
}

rmarkdown::render(opt$reportTemplate,
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

