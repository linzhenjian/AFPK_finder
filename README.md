# AFLP_finder
fast identify Animal FAS Like PKS (AFLP)s
It utilizes the t-Distributed Stochastic Neighbor Embedding (t-SNE) technique for dimensionality reduction, which is based on the alignment scores to a series of ketosynthase related Hidden Markov Model (HMM). The HDBSCAN clustering method is then employed to identify the possible AFLP clusters for each protein sequence input.

## Requirements
* Bash
* R >= 3.2
* R packages: ggplot2, Rtsne, getopt, dbscan
*To install R packages just run R environment and execute command:
*install.packages(c("ggplot2","Rtsne","getopt","dbscan"))*

* HMMER 3.3.2


AFLP_finder was tested on Mac, Ubuntu 16.04 and Ubuntu 18.04.

## How to use
At first, add the folder with executables to your PATH variable.
then run "run_identifier.sh <KS-protein.fa> <output_path>"
