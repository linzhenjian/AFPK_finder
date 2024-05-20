#!/usr/bin/env Rscript
library(Rtsne)
library(dbscan)
library(ggplot2)
library(getopt)
#library(RColorBrewer)
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

data <- read.delim(opt$inputfile,header=T)

col_count = ncol(data)
row_count <- nrow(data)
in_data <- as.matrix(data[1:row_count,2:(col_count-1)])
col_mean <- apply(in_data,2,mean)  
norm_data <- sweep(in_data,2,col_mean,"/")
number <- round(row_count / 100)/10
perplexities = ceiling(c(30,40,50,60,70,80,90,100) * number**0.3)

mytsne <- function(i, norm_data) {
  require(Rtsne)
  tsne <- Rtsne(norm_data, perplexity = i, check_duplicates = FALSE)
  return(tsne$Y)
}
file_name <- basename(opt$inputfile)
myhdbscan <- function(tsne){
	ds <- hdbscan(tsne,minPts=20)
	tsne_data <- data.frame(data[1:row_count,1],tsne,ds$cluster,data[1:row_count,col_count])
	colnames(tsne_data)=c("ID","X","Y","cluster","clade")
	datafile <- paste0(opt$output_folder,"/",file_name,"_tsne-db",i,".csv")
	write.table(tsne_data,file=datafile,sep="\t",row.names = FALSE, quote=FALSE)
	p <- ggplot(tsne_data) + 
	geom_point(aes(x=X,y=Y,colour=as.factor(clade),alpha=0.7))+ theme_bw() + theme(panel.grid=element_blank()) +
	scale_color_manual(values=c("red","blue","green", "purple","black","gray","pink","yellow","#8DEEEE", "#006400","#FFFF00","#191970"))
	return(p)
}

for(i in perplexities){
	tsne <- mytsne(i, norm_data)
	myhdbscan(tsne)
	pngfile <- paste0(opt$output_folder,"/",file_name,"_tsne",i,".png")
	ggsave(file=pngfile)
}
cat("AFLP analysis completed!!\n")
q()



