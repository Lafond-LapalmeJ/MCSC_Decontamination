# MCSC Decontamination method


This is a pipeline for decontamination of biological sequences

Author: Joel Lafond Lapalme
email for questions/support: joel.lafond.lapalme@gmail.com


For more information on the method, see the paper at: https://academic.oup.com/bioinformatics/article-abstract/doi/10.1093/bioinformatics/btw793/2736361/A-new-method-for-decontamination-of-de-novo
### Included ###
 
 ./MCSC_decontamination.sh             # script for pipeline
 ./example.ini                         # parameters required for the pipeline 
 ./scritps/                            # all required scripts
 ./data/			       # required data for taxonomy
 ./test/			       # example

### Requirements ###

 In order to run the pipeline you the need the following
 - DIAMOND (v0.9.22; may work with older version) blast in your $PATH https://github.com/bbuchfink/diamond
 - perl 5.18.2 (Pipeline was test with this version, older version could work)
 - R 3.1.1 (Pipeline was test with this version, older version could work)
 - Uniref90 database build by DIAMOND makedb command (http://www.uniprot.org/downloads)
 - Uniref100 taxlist (https://github.com/GDKO/uniref_taxlist)


## Installation

1) Clone the MCSC_Decontamination repository
```
git clone https://github.com/Lafond-LapalmeJ/MCSC_Decontamination.git
```
2) Install DIAMOND
```
wget http://github.com/bbuchfink/diamond/releases/download/v0.9.22/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
#add this line to your .bashrc to add DIAMOND in your $PATH
export PATH=$PATH:PATH_TO_DIAMOND/diamondblast
```
3) Get the uniref100 taxlist
```
git clone https://github.com/GDKO/uniref_taxlist.git
cat uniref100.taxlist.gz.part-a* | gunzip > uniref100.taxlist
```
4) Get and build uniref90 database with DIAMOND
```
wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz | gunzip
diamond makedb --in uniref90.fasta --db uniref90
```
5) Download the ncbi taxonomy dmp
```
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz | tar -xvf
```




### Example.ini

 The .ini file tell the pipeline your parameters, path and file needed to run the pipeline
 The original example.ini is available at https://github.com/Lafond-LapalmeJ/MCSC_Decontamination




### Pipeline 
 To run the pipeline simply call:
```
 sh MCSC_decontamination.sh file.ini
```

### Output
 The pipeline will output the following files:
 - n fasta files representing the n clusters
 - cluster_eval.tsv (tab with statistic on all clusters)
 - cluster_eval.jpeg (plot of the clustering)
 - file_decont.fasta (decontaminated fasta file)



## Test data

If you follow the installation you should be able to decontaminate the test data.
First add the path to the uniref90 DIAMOND database and uniref100 in test.ini
Then, just run the command:

```
sh MCSC_decontamination.sh test/test.ini
```


### Time and performance
The longest part of the method is the DIAMOND blast of the fasta file.
If you already have a DIAMOND blast of this file, just insert the file.daa
in the output directory. The pipeline will skip the blastx and use the file.daa
as DIAMONd blast output.

Also if you misspelled the target taxon ($WHITE) you can run only the WR index calculation by changing your .ini file and call
```
sh MCSC_Decontamination file.ini --recalculate
```

### Alignment correction
The parameter AST (alignment score threshold) is disable by default (AST=-1).
With a AST value greater or equal than 0, Sequences with an alignment and a BIT score higher the AST value are classified by the taxonomy of that alignment instead of being classified with the MCSC cluster.
All the sequences without alignment to the database are classified according to the MCSC clustering.
Finally, sequences with a BIT score lower than the AST value are classified according to their MCSC cluster.

For example, if we set the AST to 100 and we want to decontaminate a nematode transcriptome. A sequence with an alignment to a non-nematode with a BIT score over 100 will be discarded whatever the clustering of this sequence and a sequence with an alignment to a nematode with a BIT score over 100 will be kept in the decontaminated transcriptome. Setting the AST to 0 means to classify all sequences that have an alignment by the taxonomy of this alignment and setting the AST to infinity (or any very large number) means that no sequence will be classified according to the alignment which comes back to the default decontamination method.

