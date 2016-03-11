#!/bin/sh

##############################################################################
# This script run the MCSC decontamination pipeline.
#
# To run this script, you need Perl > 5.18.2
#
# If you do not have blast output file of your sequences on
# both databases you need blast 2.2.29+, a white-list and
# a black-list BLAST databases (See README file to get the right
# command to build those databases)
# 
# 
# if you do not have the blast output file of your sequences on both databases
# Usage: MCSC_decontamination.sh sequences.fasta MCSC_clusters_directory clustering_level(integer) new_blast_white_list_output.tsv new_blast_black_list_output.tsv blast_db_white blast_db_black
# 
# <OR>
# 
# if you already have the blast output on both databases
# MCSC_decontamination.sh sequences.fasta clustering_level(integer) MCSC_clusters_directory /full/Path/to/blast_white_list_output.tsv /full/path/to/blast_black_list_output.tsv
#
#
##############################################################################

# parameters:
# $FASTA: sequences (fasta file)
# $LVL: clustering level (integer) 3 = 8 clusters, 4 = 16 clusters ... 
# $OUT: Directory to output the clusters
# $PROGRAM: blastx or blastn
# $WHITE_TAB: Output file of the blast on the white list (.tsv file)
# $BLACK_TAB: Output file of the blast on the black list (.tsv file)
# $WHITE_DB: (optionnal) name of the white list blast DB
# $BLACK_DB: (optionnal) name of the black list blast DB

grep "=" $1 > param_file.txt

source ./param_file.txt

rm param_file.txt

## if you need to BLAST your sequences against databases
if ! [ -z "$WHITE_DB" ]
then
 	mkdir -p "$OUT"
	
	FILE=$(basename "$FASTA")
	NAME="${FILE%.*}"

	
	## Get the directory of the command.
	SOURCE="${BASH_SOURCE[0]}"
	while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  		DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  		SOURCE="$(readlink "$SOURCE")"
  		[[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
	done
	DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

	## call the MCSC algorithm to cluster the sequences
	OUT_NAME=${OUT}/${NAME}_

	$DIR/scripts/mcsc_fix $FASTA $LVL $OUT_NAME


	## blast on both datases
	$PROGRAM -db $WHITE_DB -query $FASTA -out $OUT/$WHITE_TAB -outfmt 6 -max_target_seqs 1 -num_threads 24 -evalue 1e-10
	$PROGRAM -db $BLACK_DB -query $FASTA -out $OUT/$BLACK_TAB -outfmt 6 -max_target_seqs 1 -num_threads 24 -evalue 1e-10
	
	eval "$a"
	eval "$b"

	## sort to get the best alignment
	sort -nk12,12 $OUT/$WHITE_TAB >$DIR/temp1
	sort -nk12,12 $OUT/$BLACK_TAB >$DIR/temp2
	mv $DIR/temp1 $OUT/$WHITE_TAB
	mv $DIR/temp2 $OUT/$BLACK_TAB
	
	## evaluate clusters
	perl $DIR/scripts/clusters_evaluation.pl $OUT $OUT/$WHITE_TAB $OUT/$BLACK_TAB

	## print a plot of the cluster evaluation
	R CMD BATCH --no-save --no-restore "--args $OUT/cluster_eval.tsv" $DIR/scripts/clusters_evaluation.R $OUT/clusters_evaluation.Rout

	
elif [ -z $WHITE_DB ]
then

	mkdir -p $OUT
        
        FILE=$(basename "$FASTA")
        NAME="${FILE%.*}"

        
        ## Get the directory of the command.
        SOURCE="${BASH_SOURCE[0]}"
        while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
                DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
                SOURCE="$(readlink "$SOURCE")"
                [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
        done
        DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

        ## call the MCSC algorithm to cluster the sequences

	OUT_NAME=${OUT}/${NAME}_
        $DIR/scripts/mcsc_fix $FASTA $LVL $OUT_NAME  


        ## sort to get the best alignment
        sort -nk12,12 $WHITE_TAB >$DIR/temp1
        sort -nk12,12 $BLACK_TAB >$DIR/temp2
        mv $DIR/temp1 $WHITE_TAB
        mv $DIR/temp2 $BLACK_TAB
        
        ## evaluate clusters
        perl $DIR/scripts/clusters_evaluation.pl $OUT $WHITE_TAB $BLACK_TAB

        ## print a plot of the cluster evaluation
        R CMD BATCH --no-save --no-restore "--args $OUT/cluster_eval.tsv" $DIR/scripts/clusters_evaluation.R $OUT/clusters_evaluation.Rout


else
	## If wrong number of arguments
	$STR=$'Wrong number of arguments!\n usage:\n clusters_eval.sh assembly.fasta DHCS_cluster_directory Blast_output_white Blast_output_black <blast_db_white> <blast_db_black>'
	echo $STR
fi
