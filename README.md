# MCSC Decontamination method


This is a pipeline for decontamination of biological sequences

you can use it as a pipeline with the MCSC_decontamination.sh script
or you can use only part of the pipeline (see below)

Author: Joel Lafond Lapalme
email for questions/support: joel.lafond.lapalme@gmail.com

For more information on the method, see the paper at
LINK

### Included ###
 
 ./MCSC_decontamination.sh             # script for pipeline
 ./example.ini                         # script with parameters and files required for the pipeline 
 ./scritps/mcsc_fix                    # mcsc algorithm
 ./scripts/clusters_evaluation.pl      # perl script to evaluate all clusters from the BLAST results
 ./scripts/clusters_evaluation.R       # R script that produces a jpeg of the White-ratio(WR) of every clusters


### Requirement ###

 In order to run the pipeline you the need the following
 - BLAST 2.2.26+ (or newer)
 - perl 5.18.2 (Pipeline was test with this version, older version could work)
 - R 3.1.1 (Pipeline was test with this version, older version could work)
 - White list BLAST database or a nucleotide/protein fasta file (for target species)
 - Black list BLAST database or a nucleotide/protein fasta file (for contaminants)
 - Fill the example.ini file


## hints and tips

 If you work with large databases the BLAST command will be long.
 You should consider make the blastx/blastn on a parrallel computer.

 Your white list should be a union of sequences from known organisms 
 that are close to your target species. Like C. elegans if you work
 on an unknown nematode.

 Your black list should be a union of sequences from known contaminants
 that are close the contaminants from your sample. If you have no idea
 which contaminants are in your sample, you should BLAST all (or sub-sample)
 your data on the NR databases of NCBI. By computing the species distribution
 of your best BLAST hit, you should have an idea of which type of organism
 are in your sample.






### Build databases

 If you do not already have your white or black list BLAST database you can use the following command:

 For nucleotide sequences
  $  makeblastdb -in white_list.fasta -dbtype 'nucl' -title white_list_db -out white_list_db -parse_seqids 

 For protein sequences
  $  makeblastdb -in white_list.fasta -dbtype 'prot' -title white_list_db -out white_list_db -parse_seqids 



### Example.ini

 The .ini file tell the pipeline your parameter and file needed to run the pipeline
 The original example.ini is available at https://github.com/Lafond-LapalmeJ/MCSC_Decontamination




### Pipeline 
 To run the pipeline simply call:

 $ sh MCSC_decontamination.sh file.ini

 The pipeline will output the following files:
 - n fasta files representing the n clusters
 - cluster_eval.tsv (tab with statistic on all clusters)
 - cluster_eval.jpeg (plot of the clustering)




### Step by step 


 Here are the commands to run the pipeline step by step
 1) clustering of the file to decontaminate
  $ ./dhcs_fix file.fasta Clustering_Level

 2) Blast on the white list and black list  databases
  $ blastx -db white_list -query file.fasta -out white_list_blast_output.tsv -outfmt 6 -max_target_seq 1 -evalue 1e-10
  $ blastx -db black_list -query file.fasta -out black_list_blast_output.tsv -outfmt 6 -max_target_seq 1 -evalue 1e-10
 
 3) Sort the BLast output
   $ sort -nk12,12 white_list_blast_output.tsv > white_list_blast_output_sorted.tsv
   $ sort -nk12,12 black_list_blast_output.tsv > black_list_blast_output_sorted.tsv

 4) evaluate clusters
   $ perl ./clusters_evaluation.pl /path/to/cluster/directory white_list_blast_output_sorted.tsv black_list_blast_output_sorted.tsv
   
 5)print a plot of the cluster evaluation 
  $ R CMD BATCH --no-save --no-restore "--args cluster_eval.tsv" clusters_evaluation.R
 
 

