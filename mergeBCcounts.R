#!/usr/bin/env Rscript

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
# args <- c('/Volumes/groups/obenauf/Kimon_Froussios/chris/CE83TANXX_review/process/CE83TANXX_7_20200221B_20200224', 'test.txt', '_barcode-counts.txt')

patt <-  args[3]
files <- dir(args[1], pattern=patt, full.names=TRUE)
out <- args[2]

# f <- files[1]
tables <- lapply(files, function(f) {
	DT <- fread(f)
	names(DT) <- c(ifelse(grepl('barcode', patt), 'barcode', 'sample'), sub(patt, '', basename(f)))
	DT
})

id <- names(tables[[1]])[1]
counts <- Reduce(function(df1, df2) merge(df1, df2, by=id, all=TRUE), tables)
counts[is.na(counts)] <- 0

# On-screen report of library sizes
print("Read Counts")
print( colSums(counts[, names(counts) != id, with=FALSE]) )

# On-screen report of barcode counts
print("Barcode Counts")
vapply(counts[, names(counts) != id, with=FALSE], 
			 function(x) { sum(x>0) }, 
			 integer(1) )

fwrite(counts, file=out, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
