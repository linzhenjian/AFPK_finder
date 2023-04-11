#!/bin/bash
ml R
if [ $# -ne 2 ]; then
    echo "Usage: $0 <fasta_file> <output directry>"
    usage
    exit 1
fi
file="$1"
output="$2"
rm -r $output
mkdir $output
if [[ ! -f $file ]]; then
  echo "Error: $file does not exist or is not a regular file."
  exit 1
fi

if [[ ! $(file -b --mime-type "$file") == "text/plain" ]]; then
  echo "Error: $file is not a text file."
  exit 1
fi

result=`cat $file | grep -v '^>' | grep -i -e [FPEJLZOIQ*X] |wc -l`

if (($result <= 0)); then
  echo "Error: $file is not a protein fasta file."
  exit 1	
fi

BIN_PATH="`dirname \"$0\"`"
#export PATH=$BIN_PATH

if ! ( which hmmsearch > /dev/null ); then echo "You should install hmmer3.";exit 1; fi

fasta_head=($(grep ^\> $file | awk '{print $1}' | sed 's/>//g' | awk '!seen[$1]++'))

for i in {1..31}; do
        $BIN_PATH/hmmsearch --noali --cpu 4 --tblout $output/$i.hmm_output $BIN_PATH/hmm/$i.hmm $file
done
echo -n "" > $output/hmm_matrix

for seq_name in ${fasta_head[*]}; do
	echo -n $seq_name" " >> $output/hmm_matrix
	for i in {1..31}; do
		score=$(awk '$1=="'"${seq_name}"'" {print $6}' $output/$i.hmm_output | sed -n 1p)
	if [ "$score" == "" ]; then 
		score=5
	fi
	echo -n $score" " >> $output/hmm_matrix
	done 
	echo  "INPUT_KS" >> $output/hmm_matrix
	
done
rm $output/*.hmm_output
awk '{print $1}' $file | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) " "; } $0 !~ ">" {c+=length($0);} END { print c; }' > $output/length
awk  'NR==FNR{a[$1]=$0;next} NR>FNR{print a[$1],$2}' $output/hmm_matrix $output/length | sed 's/ /\t/g' |  awk '{for(i=1;i<=NF;i++)if($i>300){print $0;next}}'>  $output/hmm_tab.txt
cat    $BIN_PATH/training_data.txt $output/hmm_tab.txt > $output/data.txt
$BIN_PATH/script.r -i $output/data.txt -o $output
