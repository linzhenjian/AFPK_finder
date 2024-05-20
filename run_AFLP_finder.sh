#!/usr/bin/env bash

if [[ -z "$1" || $1 == -h* ]]; then
    cat <<EOF
AFLP-finder Schmidt Lab, University of Utah
requirments: hmmer3 prodigal
Usage: $0  -p <protein fasta file>/-d <dna fasta file> -h <HMM models in a folder> [ -o <output folder> ] [ -t <CPU threads> ] [ -m <training data> ] [ -l <training prot annotation table> ] [ -x <model> ]
example:  run_AFLP_finder.sh -p protein.fa -o ./output -m ./FAD_halogenase_finder_2-24-2024/traindata.fa   -l ./FAD_halogenase_finder_2-24-2024/train_anno.txt  -h ./FAD_halogenase_finder_2-24-2024/hmm
EOF
    exit
fi

BIN_PATH="`dirname \"$0\"`"

# Parse command line options
while getopts "p:d:o:t:m:l:h:x:" opt; do
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
    m)
      train_data="$OPTARG"
      ;;
    l)
      train_anno="$OPTARG"
      ;;
    h)
      hmms="$OPTARG"
      ;;
    x)
      model="$OPTARG"
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

# Create the output folder
mkdir $output

# Check if either a protein or DNA FASTA file exists

if ! ( [[ -f $prot ]]  || [[ -f $nucl ]] );then
        echo "Input files not found!"
        exit 1
fi

# Process protein FASTA file

if ( [[ -f $prot ]] && [[ -z $nucl ]] ); then
	file=$prot
	result=`cat $prot | grep -v '^>' | grep -i -e [FPEJLZOIQ*X] |wc -l`
fi


# Check if the file is a protein FASTA file

if (($result <= 0)); then
  echo "Error: $file is not a protein fasta file."
  exit 1	
fi

if ( [[ -z $prot ]] && [[ -f $nucl ]] ); then
	prodigal -q -i $nucl  -a $nucl.prot  -p meta 	
	file=$nucl.prot
fi


if ( [[ -z $hmms ]] ); then
        hmms="$BIN_PATH"/hmm
fi

if ( [[ -z $model ]] ); then
	model="$BIN_PATH"/training_data.txt
fi

if ! ( which R > /dev/null ); then echo "You should install R.";exit 1; fi
if ! ( which prodigal > /dev/null ); then echo "You should install prodigal.";exit 1; fi
if ! ( which hmmsearch > /dev/null ); then echo "You should install hmmer3.";exit 1; fi

# making the training data matraxi, if the training data provided, the model will be re-calculated

if ( [[ -f $train_data ]] && [[ -f $train_anno ]] ); then
	# Check if the training file is a protein FASTA file
	if ( [[ -f $train_data ]] ); then
	
		result=`cat $train_data | grep -v '^>' | grep -i -e [FPEJLZOIQ*X] |wc -l`
	fi

	if (($result <= 0)); then
  	echo "Error: $train_data is not a protein fasta file."
  	exit 1	
	fi
	# creat training matrix

	for hmm in "$hmms"/*.hmm; do
	 	file_name=$(basename "$hmm" | sed 's/[.]hmm//')	
       		hmmsearch --noali --cpu $THREADS --tblout $output/$file_name.hmm_output $hmm $train_data
	done
	echo "" > $output/train_temp
	for i in $output/*.hmm_output; do   awk  ' {print $1,"hmm_"$3,$9}' $i | sed '/#/d' | awk '!seen[$1]++' >> $output/train_temp; done
	#for i in $output/*.hmm_output; do name=$(basename $i | sed 's/[.]hmm_output//');  awk -v var="$name" ' {print $1,var,$9}' $i | sed '/#/d' | awk '!seen[$1]++' >> $output/train_temp; done
	
	cat $output/*.hmm_output | awk '{print $1}' | sed '/#/d' | awk '!seen[$1]++' > $output/list
	awk  'NR==FNR{a[$0]}NR>FNR{if ($1 in a) print $0}'  $output/list  $output/train_temp > $output/signifit_hit
	awk '{print $1,"type",$2}' $train_anno >> $output/signifit_hit
	
	$BIN_PATH/make_train.r -i $output/signifit_hit -o $output
	rm $output/*.hmm_output
	model=$output/training_data.txt
fi

#get the hmm scores for input KSs
for hmm in "$hmms"/*.hmm; do
	 file_name=$(basename "$hmm" | sed 's/[.]hmm//')	
        hmmsearch --noali --cpu $THREADS --tblout $output/$file_name.hmm_output $hmm $file
done

echo "" > $output/outputqq
for i in $output/*.hmm_output; do awk ' {print $1,"hmm_"$3,$9}' $i | sed '/#/d' | awk '!seen[$1]++' >> $output/outputqq; done
cat $output/*.hmm_output | awk '$9>100 {print $1}' | sed '/#/d' | awk '!seen[$1]++' > $output/list
awk  'NR==FNR{a[$0]}NR>FNR{if ($1 in a) print $0}'  $output/list  $output/outputqq > $output/signifit_hit
awk '{print $1,"type","INPUTseq"}' $output/list >> $output/signifit_hit
$BIN_PATH/combine_table.r -i $output/signifit_hit -t $model -o $output
rm $output/*.hmm_output

#normalize the hmm scores by row

awk 'BEGIN { FS="\t"; OFS="\t" } NR>1 { max=0; for (i=2; i<NF; i++) { if ($i>max) max=$i } for (i=2; i<NF; i++) $i=$i/max;   print}' $output/data.txt > $output/normal_data.txt  # normalized the data by the max in each row


$BIN_PATH/script.r -i $output/normal_data.txt -o $output
$BIN_PATH/script.r -i $output/data.txt -o $output
 
#summarize the clustering
for csv in $(ls $output/*tsne-db*.csv); do
        name=$(echo $csv | awk -F '/' '{print $NF}' | sed 's/[.]csv//')
        echo "ID cluster clade percentage" | sed 's/ /\t/g' > $output/$name.id
        awk '$5=="INPUTseq" {print $1,$4}' $csv > $output/temp_a
        awk '$5!="INPUTseq" && $5!="clade" {print $4,$5}'  $csv | awk '{a[$0]++} END{for(i in a){print i,a[i] | "sort -n -r -k 2"}}' > $output/temp_b
        awk '$5!="INPUTseq" && $5!="clade" {print $5}'  $csv | awk '{a[$0]++} END{for(i in a){print i,a[i] | "sort -n -r -k 2"}}' > $output/temp_c
        awk  'NR==FNR{a[$1]=$2;next} NR>FNR{print $0,a[$2]}'  $output/temp_c $output/temp_b | awk '{print $1,$2,($3/$4)*100"%"}'  | awk '{a[$1]=a[$1]" "$2" "$3} END {for (i in a) print i a[i]}'  > $output/temp_f
        awk  'NR==FNR{a[$1]=$0;next} NR>FNR{print $0,a[$2]}' $output/temp_f $output/temp_a | cut -d' ' -f1,3- | sed 's/ /\t/g' >> $output/$name.id
done
rm $output/temp_*

