#!/usr/bin/env bash

if [[ -z "$1" || $1 == -h* ]]; then
    cat <<EOF
AFLP-finder Schmidt Lab, University of Utah
requirments: hmmer3 prodigal
Usage: $0  -p <protein fasta file>/-d <dna fasta file> [ -o <output folder> ] [ -t <CPU threads> ]
EOF
    exit
fi

# Parse command line options
while getopts "p:d:o:t:m:" opt; do
  case $opt in
    p)
      prot="$OPTARG"
      ;;
    d)
      nucl="$OPTARG"
      ;;
    o)
      output="$OPTARG"
      ;;
    t)
      THREADS="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

if [ -z $output ]; then
	mkdir ./AFLP_output
	output="./AFLP_output"
fi

if [ -d "$output" ]; then
    rm -rf "$output"  # delete existing folder
fi


if [ -z $THREADS ]; then THREADS=4; fi


mkdir $output

if ! ( [[ -f $prot ]]  || [[ -f $nucl ]] );then
        echo "Input files not found!"
        exit 1
fi

if ( [[ -f $prot ]] && [[ -z $nucl ]] ); then
	file=$prot
	result=`cat $prot | grep -v '^>' | grep -i -e [FPEJLZOIQ*X] |wc -l`
fi



if (($result <= 0)); then
  echo "Error: $file is not a protein fasta file."
  exit 1	
fi

if ( [[ -z $prot ]] && [[ -f $nucl ]] ); then
	prodigal -q -i $nucl  -a $nucl.prot  -p meta 	
	file=$nucl.prot
fi

BIN_PATH="`dirname \"$0\"`"


if ! ( which R > /dev/null ); then echo "You should install R.";exit 1; fi
if ! ( which prodigal > /dev/null ); then echo "You should install prodigal.";exit 1; fi
if ! ( which hmmsearch > /dev/null ); then echo "You should install hmmer3.";exit 1; fi


for i in {1..30}; do
        hmmsearch --noali --cpu $THREADS --tblout $output/$i.hmm_output $BIN_PATH/hmm/$i.hmm $file
done
fasta_head=($(cat $output/*.hmm_output | awk '$6>300 {print $1}' | sed '/#/d' | awk '!seen[$1]++'))
echo -n "" > $output/hmm_matrix

for seq_name in ${fasta_head[*]}; do
	echo -n $seq_name" " >> $output/hmm_matrix
	for i in {1..30}; do
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
awk  'NR==FNR{a[$1]=$0;next} NR>FNR{print a[$1],$2}' $output/hmm_matrix $output/length | sed 's/ /\t/g'  | awk '{ for (i=2; i<NF-1; i++) if ($i > 300) { print; next } }' >  $output/hmm_tab.txt
cat    $BIN_PATH/training_data.txt $output/hmm_tab.txt > $output/data.txt
$BIN_PATH/script.r -i $output/data.txt -o $output

for csv in $(ls $output/tsne-db*.csv); do
        name=$(echo $csv | awk -F '[/.]' '{print $(NF-1)}')
        echo "ID cluster clade percentage" | sed 's/ /\t/g' > $output/$name.id
        awk '$5=="INPUT_KS" {print $1,$4}' $csv > $output/temp_a
        awk '$5!="INPUT_KS" && $5!="clade" {print $4,$5}'  $csv | awk '{a[$0]++} END{for(i in a){print i,a[i] | "sort -n -r -k 2"}}' > $output/temp_b
        awk '$5!="INPUT_KS" && $5!="clade" {print $5}'  $csv | awk '{a[$0]++} END{for(i in a){print i,a[i] | "sort -n -r -k 2"}}' > $output/temp_c
        awk  'NR==FNR{a[$1]=$2;next} NR>FNR{print $0,a[$2]}'  $output/temp_c $output/temp_b | awk '{print $1,$2,($3/$4)*100"%"}'  | awk '{a[$1]=a[$1]" "$2" "$3} END {for (i in a) print i a[i]}'  > $output/temp_f
        awk  'NR==FNR{a[$1]=$0;next} NR>FNR{print $0,a[$2]}' $output/temp_f $output/temp_a | cut -d' ' -f1,3- | sed 's/ /\t/g' >> $output/$name.id
done
rm $output/temp_*

