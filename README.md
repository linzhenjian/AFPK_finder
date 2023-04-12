# AFLP_finder
fast identify AFLPs
It's based on t-Distributed Stochastic Neighbor Embedding (t-SNE) (van der Maaten, 2014) technique for dimensionality reduction, that uses hmm align scores. HDBSCAN clustering method (Campello et.al., 2013) is used to identify the posible clade (AFLP) for each preotein sequence input.

## Requirements
* Bash
* Perl >= 5.10
* R >= 3.2
* R packages: ggplot2, Rtsne, getopt, dbscan

To install R packages just run R environment and execute command:
*install.packages(c("ggplot2","Rtsne","getopt","dbscan"))*


AFLP_finder was tested on Ubuntu 16.04 and Ubuntu 18.04.

## How to use
At first, add the folder with executables to your PATH variable.
then run "run_identifier.sh <KS-protein.fa> <output_path>"
