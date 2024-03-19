#!/usr/bin/env Rscript
library(reshape2)
library(Rtsne)
library(dbscan)
library(ggplot2)
library(getopt)

spec <- matrix(c(
  'help', 'h', 0, "logical", "Print help message",
  'inputfile', 'i', 1, "character", "pleae provide a table file with hmm alignment core ",
   'training_data', 't', 1, "character", "pleae provide a table file with hmm alignment core ",
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
if ( is.null(opt$training_data) ) {
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
train <- read.table(opt$training_data, header = TRUE, sep = "\t")

data_cast <- dcast(data, V1 ~ V2, value.var="V3")
data_cast[is.na(data_cast)] <- 1
result <- merge(data_cast, train, by = intersect(names(data_cast), names(train)), all = TRUE)
result[is.na(result)] <- 1
type_index <- which(names(result) == 'type')
total_cols <- ncol(result)
if (type_index != total_cols) {
result <- result[, c(1:(type_index - 1), (type_index + 1):ncol(result), type_index)]

}

file <- paste0(opt$output_folder,"/data.txt")
write.table(result, file = file, sep = "\t", row.names = FALSE, col.names = TRUE,quote = FALSE)

