#!/usr/bin/env Rscript

## version 0.8.2

library(getopt)

# Full paths are needed for input files.
spec = matrix(c(
  'help'            , 'h', 0, "logical",   "Help",
  'countsFile'      , 'c', 1, "character", "Tab-separated table of counts with the collected counts for the samples.",
  'summariesFile'   , 's', 1, "character", "Tab-separated table of counts with the collected summaries for the demultiplexing/counting.",
  'hammingFile'     , 'e', 1, "character", "Tab-separated table with the collected hamming distances for the samples.",
  'resultsDir'      , 'd', 1, "character", "Directory in which to save all output.",
  'reportTemplate'  , 'T', 1, "character", "Template Rmd file for CaTCH report.",
  'covars'          , 'v', 1, "character", "Table file listing sample names ('Sample') (in desired order), coarse sample grouping for identifying shared barcodes ('Group'), respective condition for each sample ('Treatment') and desired display colour for each sample ('Colour') (in R-compatible format, preferably avoid 'white'). For convenience, an optional field 'Tag' is allowed as the first column, to permit using one table for the whole pipeline. Full path needed.",
  'refsamps'        , 'r', 2, "character", "Optional, comma-separated list of rows designating the samples to use as reference abundance (in order they appear in covars, not counting the header row). Used for naming and filtering the barcodes. If omitted, the barcodes will be named by their abindance in thee first available sample, and no filtering based on reference count will take place.",
  'count_thresh'    , 'N', 1, "integer",   "Count threshold for barcodes to analyse",
  'abund_thresh'    , 'A', 1, "numeric",   "Proportional abundance threshold to consider barcodes top hits",
  'extrabc'         , 'B', 1, "character", "Comma-separated list of barcode IDs to include in addition to those chosen by the analysis. No spaces."
), byrow=TRUE, ncol=5)

opt = getopt(spec)

# opt <- list(countsFile='/Volumes/groups/obenauf/Kimon_Froussios/shona/catch_160598/process/data_barcode-counts.txt', summariesFile='/Volumes/groups/obenauf/Kimon_Froussios/shona/catch_160598/process/data_summaries.txt', covars='/Volumes/groups/obenauf/Kimon_Froussios/shona/catch_160598/description/samples.txt', resultsDir='/Volumes/groups/obenauf/Kimon_Froussios/shona/catch_160598/results', count_thresh=50, reportTemplate='~/catch/barcoding_results_template.Rmd')

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

dir.create(opt$resultsDir)


if ( is.null(opt$count_thresh) ) {
  opt$count_thresh <- 10L
}
if ( is.null(opt$abund_thresh) ) {
  opt$abund_thresh <- 0.005
}

if ( !is.null(opt$refsamps) ) {
  opt$refsamps <- as.integer(strsplit(opt$refsamps, ',', fixed=TRUE)[[1]])
}

if (! is.null(opt$extrabc) ) {
  opt$extrabc <- strsplit(opt$extrabc, ',', fixed=TRUE)[[1]]
  opt$extrabc <- opt$extrabc[order(opt$extrabc)]
}

# Sanitize
DT <- read.delim(opt$covars, stringsAsFactors = FALSE)

if( (length(DT) == 4 && names(DT) != c('Sample', 'Group', 'Treatment', 'Colour')) || 
    (length(DT) == 5 && names(DT) != c('Tag', 'Sample', 'Group', 'Treatment', 'Colour')) ) {
  stop('The covariates table must have fields named Tag, Sample, Group, Treatment, Colour, in that order. The Tag column is optional, the others are required.')
}

bad <- grepl('(^[0-9])|([^0-9a-zA-z_])', DT$Sample, perl=TRUE)
if (any(bad)) {
  stop(paste('Invalid sample name(s) detected (',
             paste(DT[bad, 'Sample'], collapse=", "),
             '). Names must start with a letter, and contain no spaces and no symbols (except underscore).'))
}

bad <- grepl('(^[0-9])|([^0-9a-zA-z_])', DT$Group, perl=TRUE)
if (any(bad)) {
  stop(paste('Invalid group name(s) detected (',
             paste(DT[bad, 'Group'], collapse=", "), 
             '). Names must start with a letter, and contain no spaces and no symbols (except underscore).'))
}

bad <- grepl('(^[0-9])|([^0-9a-zA-z_])', DT$Treatment, perl=TRUE)
if (any(bad)) {
  stop(paste('Invalid treatment name(s) detected (',
             paste(DT[bad, 'Treatment'], collapse=", "), 
             '). Names must start with a letter, and contain no spaces and no symbols (except underscore).'))
}

areColors <- function(x) {
  vapply(x, function(X) {
    tryCatch(is.matrix(col2rgb(X)), 
             error = function(e) FALSE)
  }, logical(1))
}

bad <- ! areColors(DT$Colour)
if (any(bad)) {
  stop( paste(paste(DT[bad, 'Colour'], collapse=", "), 
              "is/are invalid colour value(s).") )
}

bad <- ! DT$Sample %in% names(read.delim(opt$countsFile))
if (any(bad)) {
  stop(paste("Some sample names in the covariates table were not found in the counts table (", 
             paste(DT$Sample[bad], collapse=", "),") ."))
}

rmarkdown::render(opt$reportTemplate,
                  output_file = sub('.txt|.tsv', '_report.html', basename(opt$countsFile)),
                  output_dir = opt$resultsDir,
                  params=list(counts = opt$countsFile,
                              summaries = opt$summariesFile,
                              hammdist = opt$hammingFile,
                              outpref = file.path(opt$resultsDir, sub('.txt|.tsv', '', basename(opt$countsFile))),
                              samples = opt$covars,
                              refsamps = opt$refsamps,
                              count_thresh = opt$count_thresh,
                              abund_thresh = opt$abund_thresh,
                              extra_bc = opt$extrabc)
)

