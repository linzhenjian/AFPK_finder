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


AFLP_finder was tested on MacOs 11, Ubuntu 16.04 and Ubuntu 18.04.

## How to use
At first, add the folder with executables to your PATH variable.
then run "run_identifier.sh <KS-protein.fa> <output_path>"

## results:
this plot shows the INPUT_KSs (purple) are adjacent to AFLPs (light blue),  so the INPUT_KSs are identified as AFLPs.
![tsne54](https://user-images.githubusercontent.com/52417798/231541965-086593ed-edc0-4f04-be00-dcdba9387a55.png)


output table example: it tells that INPUT_KS1 is clustered in cluster 4, cluster 4 contains 100% of the training AFLPs

![Capture](https://user-images.githubusercontent.com/52417798/231544396-0c8b1ea5-bd71-4c32-b5bc-76c3eb98ef97.PNG)
