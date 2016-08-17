#!/bin/bash

##############################################################################
# This script run the MCSC decontamination pipeline.
#
# To run this script, you need:
#	Perl > 5.18.2
#	DIAMOND blast (in your $PATH)
#	UNIREF90 and UNIREF100 databases
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

# $OUT: Directory to output the clusters (default PWD)
OUT=$PWD

# $UNIREF90: path to the DIAMOND UNIREF90 database (REQUIRED) 
UNIREF90=""

# $UNIREF100: path to the UNIREF100 taxonomy list (REQUIRED)
UNIREF100=""

# $TAXDUMP: path to the NCBI taxonomy dump (REQUIRED)
TAXDUMP=""

# $MCSC: path to the MCSC_Decontamination folder (REQUIRED)
MCSC=""

# $TAXO_LVL: taxonomic level for the WR index (default: phylum)
TAXO_LVL="phylum"

# $WHITE_NAME: Name of the target taxon for the WR index (REQUIRED)
WHITE_NAME=""

# $T: number of threads (default 8)
T=8


### parameter validation



grep "=" $1 > $OUT/param_file.txt

source $OUT/param_file.txt

rm param_file.txt

## skip if you only need to run cluster evaluation with a different taxon
if [ $# -eq 1 ]
then

	mkdir -p "$OUT"	

	## get the fasta name
	FILE=$(basename "$FASTA")
	NAME="${FILE%.*}"

	## DIAMOND blast
	if ! [ -e $OUT/$NAME.daa ]
	then
		diamond blastx -d $UNIREF90 -q $FASTA -a $OUT/$NAME -t $OUT -p $T
	fi

	## Extract DIAMOND blast taxonomy
	perl $MCSC/scripts/daa_to_tagc.pl $UNIREF100 $OUT/${NAME}.daa
	
	## Format the DIAMOND output 
	perl $MCSC/scripts/get_taxa_from_diamond.pl $OUT/$NAME.daa.tagc \
	$TAXDUMP/names.dmp $TAXDUMP/nodes.dmp > $OUT/temp.txt
	sort -rk3,3 $OUT/temp.txt | sort -uk1,1 > $OUT/taxo_uniq.txt  
	sed "s/'\"//g" $OUT/taxo_uniq.txt > $OUT/temp.txt
	mv $OUT/temp.txt $OUT/taxo_uniq.txt

	## MCSC clusters name        
	OUT_NAME=${OUT}/${NAME}_ 

	## MCSC algorithm
	$MCSC/scripts/mcsc_fix $FASTA $((LVL+1)) $OUT_NAME

fi



FILE=$(basename "$FASTA")
NAME="${FILE%.*}"


## compute WR index and evaluate clusters 
perl $MCSC/scripts/cluster_eval.pl $OUT $OUT/taxo_uniq.txt $TAXO_LVL $WHITE_NAME	



## print a plot of the cluster evaluation
R CMD BATCH --no-save --no-restore "--args $OUT/cluster_eval.tsv" $MCSC/scripts/clusters_evaluation.R $OUT/clusters_evaluation.Rout

## extract the file names from the R output
FILES=($(grep ".fasta" "${OUT}"/clusters_evaluation.Rout | sed "s/.*"$NAME"/"$NAME"/"))

## print the "good" cluster files on scree
echo "${FILES[@]}" | tr " " "\n"



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

## move the cluster files in a sub directory
mkdir $OUT/clusters
mv $OUT/*cluster_*.fasta $OUT/clusters/

echo "Done"
