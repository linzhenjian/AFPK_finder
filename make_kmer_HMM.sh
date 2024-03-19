#!/bin/bash

if [ $# -ne 3 ]; then
    echo "Usage: $0 <fasta_aln file> <kmer size> <output directry> # make kmer pHMM"
    usage
    exit 1
fi
file="$1"
kmer_size="$2"
output="$3"
BIN_PATH="`dirname \"$0\"`"
file_name=$(basename "$file" | awk -F . '{print $1}')
if [[ ! -f $file ]]; then
  echo "Error: inputfile does not exist or is not a regular file."
  exit 1
fi

if [[ ! $(file -b --mime-type "$file") == "text/plain" ]]; then
  echo "Error: inputfile is not a text file."
  exit 1
fi

if [[ ! -e $output ]]; then
  echo "Error: output does not exist or is not a regular file."
  exit 1
fi

if [[ !  $kmer_size ]]; then
  echo "Error: kmer_size does not exist or is not a regular file."
  exit 1
fi

if ! ( which seqkit> /dev/null ); then echo "You should install seqkit.";exit 1; fi
if ! ( which hmmbuild > /dev/null ); then echo "You should install hmmer3.";exit 1; fi


seqkit seq -w 0 $file | sed '$!N;s/\n/\t/' > $output/fasta_aln
awk -F '\t' '{print $1}' $output/fasta_aln > $output/myhead
awk -F '\t' '{print $2}' $output/fasta_aln > $output/myseq


#making Kmer p HMM
mkdir -p $output/hmm
seq=$(head -n 1 $output/myseq)
len=${#seq}
for ((i = 1; i <= len - kmer_size; i++)); do
	cut -c $i-$((i+kmer_size)) $output/myseq > $output/k_seq
	paste $output/myhead $output/k_seq | sed 's/\t/\n/g'  | sed "/^-\{$((kmer_size+1))\}$/d" | awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} {if ($2) print ">"$0}' > $output/hmm/$file_name\_kmer$i.align
	lines_file1=$(grep \^  $output/k_seq | wc -l)
	#echo $lines_file1
	lines_file1=$((lines_file1 / 2))
	#echo $lines_file1
	lines_file2=$(grep \^ $output/hmm/$file_name\_kmer$i.align | wc -l)
	if [ "$lines_file2" -gt "$lines_file1" ]; then
		hmmbuild $output/hmm/$file_name\_kmer$i.hmm $output/hmm/$file_name\_kmer$i.align
	fi
done
