#!/bin/bash

##############################################################################
# This script run the MCSC decontamination pipeline.
#
# To run this script, you need:
#    Perl > 5.18.2
#    DIAMOND blast (in your $PATH)
#    UNIREF90 and UNIREF100 databases
#
# Usage
# MCSC_decontamination.sh file.ini
#
#
##############################################################################

# parameters:
# $FASTA: sequences (REQUIRED; fasta file)
FASTA=""

# $LVL: clustering level (integer) 3 = 8 clusters, 4 = 16 clusters (default 32)
LVL=5

# "${OUT}": Directory to output the clusters (default PWD)
OUT="$PWD"

# $UNIREF90: path to the DIAMOND UNIREF90 database (REQUIRED) 
UNIREF90=""

# $UNIREF100: path to the UNIREF100 taxonomy list (REQUIRED)
UNIREF100=""

# $TAXDUMP: path to the NCBI taxonomy dump (REQUIRED)
TAXDUMP=""

# "${MCSC}": path to the MCSC_Decontamination folder (REQUIRED)
MCSC=""

# $TAXO_LVL: taxonomic level for the WR index (default: phylum)
TAXO_LVL="phylum"

# $WHITE_NAME: Name of the target taxon for the WR index (REQUIRED)
WHITE_NAME=""

# $T: number of threads (default 8)
T=8

# $AST: alignment score threshold for alignment correction. Default -1 value means not using it. To use it, use value over or equal to 0
AST=-1


### parameter validation

grep "=" $1 > "${OUT}"/param_file.txt

source "${OUT}"/param_file.txt

rm param_file.txt

## Check AST parameter
if [ "$AST" -ge 0 ];then
	echo "Using alignment correction for sequences with BIT score higher than $AST";
else
	echo "Not using alignment correction";
fi


FILE="$(basename "$FASTA")"
NAME="${FILE%.*}"


## skip if you only need to run cluster evaluation with a different taxon
if [ $# -eq 1 ]
then

    mkdir -p "${OUT}"

    #convert fastq to fasta if detected
    #check if input file is valid
    echo "Validating input file format..."
    if [ $(grep -F ".fastq" <<< "$FASTA") ]; then #is input fastq file? Check file extension fisrt.
        if [ "${FASTA##*.}" == "gz" ]; then #is it gzipped?
            LINE1="$(zcat "$FASTA" | head -n 1)" #assuming the first line is a sequence header
            CHAR1="${LINE1:0:1}"
             if [ "$CHAR1" == "@" ]; then
                echo "Converting fastq to fasta..."
                zcat "$FASTA" | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' \
                    > "${MCSC}"/data/file.fa
                 FASTA="${MCSC}"/data/file.fa
            else
                echo "Invalid fastq file."
                exit 1
            fi
        else
            LINE1="$(cat "$FASTA" | head -n 1)"
            CHAR1="${LINE1:0:1}"
            if [ "$CHAR1"0 == "@" ]; then
                echo "Converting fastq to fasta..."
                cat "$FASTA" | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' \
                    > "${MCSC}"/data/file.fa
                FASTA="${MCMC}"/data/file.fa
            else
                echo "Invalid fastq file."
                exit 1
            fi
        fi
    else #it's a fasta
        LINE1="$(head -n 1 $FASTA)"
        CHAR1="${LINE1:0:1}"
        if [ "$CHAR1" != ">" ]; then
            echo "Invalid fasta file."
            exit 1
        fi
    fi

    echo "Input file format seems to be OK."

    ## DIAMOND blast
    if [ ! -f "${OUT}"/"${NAME}".tsv ]
    then
        echo "Running DIAMOND blast..."
        diamond blastx \
            -d "$UNIREF90" \
            -q "$FASTA" \
            -o "${OUT}"/"${NAME}".tsv \
            -f 6 \
            -p "$T"
    fi

    ## Extract DIAMOND blast taxonomy
    echo "Extracting DIAMOND blast taxonomy..."
    perl "${MCSC}"/scripts/diamond_to_tagc.pl \
        "$UNIREF100" \
        "${OUT}"/"${NAME}".tsv
    
    ## Format extracted DIAMOND blast taxonomy
    echo "Formating the DIAMOND output..."
    perl "${MCSC}"/scripts/get_taxa_from_diamond.pl \
        "${OUT}"/"${NAME}".tsv.tagc \
        "${TAXDUMP}"/names.dmp \
        "${TAXDUMP}"/nodes.dmp | \
    sort -rnk3,3 | sort -uk1,1 | \
    sed "s/'\"//g" > "${OUT}"/taxo_uniq.txt

    ## MCSC clusters name        
    OUT_NAME=""${OUT}"/"${NAME}"_"

    ## MCSC algorithm
    echo "Running the MCSC algorithm..."
    "${MCSC}"/scripts/mcsc_fix "$FASTA" $((LVL+1)) "$OUT_NAME"

    sed -i 's/\s.*$//g' ${OUT_NAME}*.fasta


fi

## Create a temporary fasta to match the blast sequence name
sed 's/\s.*$//g' $FASTA > "${OUT}"/temp.fasta


## compute WR index and evaluate clusters 
echo "Computing the White-Ratio (WR) index and evaluating the clusters..."
perl "${MCSC}"/scripts/cluster_eval.pl \
    "$OUT" \
    "${OUT}"/taxo_uniq.txt \
    "$TAXO_LVL" \
    "$WHITE_NAME"    


## Check if white_name is in the taxonomy file
if [ -z $(grep -m 1 -o "${TAXO_LVL:0:2}_${WHITE_NAME}" $OUT/taxo_uniq.txt) ]; then
    echo "The "$TAXO_LVL" "$WHITE_NAME" is not present in taxo_uniq.txt file. Aborting."
    exit 1
fi



## print a plot of the cluster evaluation
echo -e "Printing MCSC output results..."
R CMD BATCH --no-save --no-restore \
    "--args $OUT/cluster_eval.tsv" \
    "${MCSC}"/scripts/clusters_evaluation.R \
    "${OUT}"/clusters_evaluation.Rout

## extract the file names from the R output
FILES=($(grep ".fasta" "${OUT}"/clusters_evaluation.Rout | sed "s/.*"$NAME"/"$NAME"/"))

## print the "good" cluster files on screen
echo "${FILES[@]}" | tr " " "\n"


if [ "$AST" -ge 0 ];then
	awk -F "\t" -v var="$WHITE_NAME" -v align_th="$AST" '{ if ($3 > align_th && $4 ~ var) print $1}' $OUT/taxo_uniq.txt > $OUT/keep.txt
	awk -F "\t" -v var="$WHITE_NAME" -v align_th="$AST" '{ if (($3 > align_th) && (!($4 ~ var))) print $1}' $OUT/taxo_uniq.txt > $OUT/reject.txt
	for f in "${FILES[@]}"; do
		grep "^>" $OUT/$f | sed 's/>//' | sed 's/;clust_.*//' | sed 's/;repSeq_.*//' >> $OUT/keep.txt
	done

	sort -uk1,1 $OUT/keep.txt > $OUT/temp.txt
	comm -23 $OUT/temp.txt $OUT/reject.txt > $OUT/keep_final.txt

	perl $MCSC/scripts/getFastaFromList.pl $OUT/temp.fasta list $OUT/keep_final.txt > "${OUT}"/"${NAME}"_decont.fasta
else
	## merge the "good" cluster files in a new fasta file labeled "decont"
	counter=0
	for f in "${FILES[@]}"; do
	let counter+=1
	if [ "$counter" -eq 1 ]; then
		cat $OUT/$f > "${OUT}"/"${NAME}"_decont.fasta
	else
		cat $OUT/$f >> "${OUT}"/"${NAME}"_decont.fasta
	fi
	done
fi

## move the cluster files in a sub directory
mkdir "${OUT}"/clusters
mv "${OUT}"/*cluster_*.fasta "${OUT}"/clusters/
rm $OUT/temp*


echo "Done"

