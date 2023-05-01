#!/bin/bash
#load R and install packages: install.packages(c("ggplot2","Rtsne","getopt","dbscan"))
#load HMMER3
#load prodigal

# export PATH="/path_to_/AFPK_finder/":$PATH
run_AFPK_finder.sh -p testing_data/ar-clade1.fa -o test_output 
