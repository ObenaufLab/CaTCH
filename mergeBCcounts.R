#!/usr/bin/env Rscript

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
# args <- c('/Volumes/groups/obenauf/Kimon_Froussios/chris/CE83TANXX_review/process/CE83TANXX_7_20200221B_20200224', 'test.txt.')

files <- dir(args[1], pattern='barcode-counts.txt', full.names=TRUE)
out <- args[2]


# f <- files[1]
tables <- lapply(files, function(f) {
	DT <- fread(f)
	names(DT) <- c('barcode', sub('_barcode-counts.txt', '', basename(f)))
	DT
})

counts <- Reduce(function(df1, df2) merge(df1, df2, by='barcode', all=TRUE), tables)
counts[is.na(counts)] <- 0

# On-screen report of library sizes
print( colSums(counts[, names(counts) != 'barcode', with=FALSE]) )

# On-screen report of barcode counts
vapply(counts[, names(counts) != 'barcode', with=FALSE], 
			 function(x) { sum(x>0) }, 
			 integer(1) )

fwrite(counts, file=out, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
