#!/usr/bin/env Rscript
library(reshape2)
library(Rtsne)
library(dbscan)
library(ggplot2)
library(getopt)

spec <- matrix(c(
  'help', 'h', 0, "logical", "Print help message",
  'inputfile', 'i', 1, "character", "pleae provide a table file with hmm alignment core ",
   
  'output_folder', 'o', 1, "character", "pleae provide path for output files "

), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

if ( !is.null(opt$help) ) {
        cat(getopt(spec, usage=TRUE))
        q(status=1)
}

if ( is.null(opt$inputfile) ) {
        cat("No input file specified!\n")
        cat(getopt(spec, usage=TRUE))
        q(status=1)
}

if ( is.null(opt$output_folder) ) {
        cat("No output path specified!\n")
        cat(getopt(spec, usage=TRUE))
        q(status=1)
}

data <- read.table(opt$inputfile,header = FALSE, sep = " ")

data_cast <- dcast(data, V1 ~ V2, value.var="V3")
data_cast[is.na(data_cast)] <- 1

type_index <- which(names(data_cast) == 'type')
total_cols <- ncol(data_cast)
if (type_index != total_cols) {
data_cast <- data_cast[, c(1:(type_index - 1), (type_index + 1):ncol(data_cast), type_index)]
}

file <- paste0(opt$output_folder,"/training_data.txt")
write.table(data_cast, file = file, sep = "\t", row.names = FALSE, col.names = TRUE,quote = FALSE)

