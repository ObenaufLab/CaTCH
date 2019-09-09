#!/usr/bin/env Rscript

library(data.table)

args <- commandArgs(trailingOnly = TRUE)

d <- fread(args[1])

fwrite( d[order(V2, decreasing=TRUE),], file=gsub('.txt', '.tsv', args[1], fixed=TRUE),
        col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE )
