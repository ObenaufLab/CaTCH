#!/usr/bin/env Rscript

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
# args <- c('/Volumes/groups/obenauf/Kimon_Froussios/shona/catch_R9739/testpro/barcode-counts.txt', '/Volumes/groups/obenauf/Kimon_Froussios/shona/catch_R9739/testpro/000000000-G5Y7J_1_20200603B_20200604/day0_barcode-counts.txt', '/Volumes/groups/obenauf/Kimon_Froussios/shona/catch_R9739/testpro/000000000-G5Y7J_1_20200603B_20200604/L11_treat_barcode-counts.txt')

out <- args[1]
files <- args[2:length(args)]

tables <- lapply(files, function(f) {
	# f <- files[1]
	DT <- fread(f)
	names(DT) <- c(ifelse(grepl('barcode', f), 'barcode', 'sample'), sub('_barcode-counts.txt|summary.txt', '', basename(f)))
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
print(vapply(counts[, names(counts) != id, with=FALSE], 
			 function(x) { sum(x>0) }, 
			 integer(1) ) )

fwrite(counts, file=out, quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
