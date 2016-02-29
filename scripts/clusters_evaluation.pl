#!/bin/perl
#########################################################################################
#
#	This perl script evaluate clusters from DHCS
#
#   usage: perl DHCS_clusters.pl DHCS_directory white_blast_output black_blast_output
#		
#
#########################################################################################
use strict;
use warnings;
use File::Basename;


sub log10 {
    my $n = shift;
    return log($n)/log(10);
}


;
# Get file
my $cluster_dir = shift;
my $white_blast_out = shift;
my $black_blast_out = shift;
# Get the directory, file name and extention
my ( $white_filename, $out_dir, $ext ) = fileparse( $white_blast_out, qr{ \. [^.]+ \z }msx );
my $outputName = "$out_dir" . "cluster_eval.tsv";
my %whitehash;
my %blackhash;
my %clusterhash;
# Open white list blast output file
open my $blast, $white_blast_out or die "Could not open $white_blast_out: $!";
# iterate line by line
while( my $line = <$blast>)  {   
    chomp $line;  
    
    #split line in array  
    my @lineArray = split(/\t|;|\s/, $line);
	my $str = shift @lineArray;
	@{ $whitehash{$str} } = @lineArray;
	
	#debug
	#print "$str\t$lineArray[4]\n;
    
}
close $blast;
# Open black list blast output file
open my $blast2, $black_blast_out or die "Could not open $black_blast_out: $!";
# iterate line by line
while( my $line = <$blast2>)  {   
    chomp $line;  
    
    #split line in array  
    my @lineArray = split(/\t|;|\s/, $line);
	
	my $str = shift @lineArray;
	@{ $blackhash{$str} } = @lineArray;
	
	#debug
	#print "$str\t$lineArray[4]\n;
    
}
close $blast2;
# Get list of cluster files
my @cluster_files;
opendir (my $DIR, $cluster_dir) or die "Cannot open $cluster_dir: $!";
@cluster_files = readdir($DIR);
@cluster_files = grep(/fasta/, @cluster_files);
closedir($DIR);
# Get cluster hash
foreach my $file (@cluster_files){
	my @contigs;
	my $str;
	open my $fasta, "$cluster_dir/$file" or die "Could not open $cluster_dir/$file: $!";
	while( my $line = <$fasta>)  {
		if( $line =~ /^>/ ){
			my @lineArray = split(/\t|;|\s/, $line);
			$str = $lineArray[0];
			$str = substr($str, 1);
			push @contigs, $str;
			
			#debug
			#print "$str\n";
			
		}
	}
	@{ $clusterhash{$file} } = @contigs;
	close $fasta;
}
# Create and open output file
open( OUTFILE, '>', "$outputName") or die "Can't create \"$outputName\": $!\n";

# Print header
print OUTFILE "Cluster\tNb_contigs\tWR\tNb_white_hit\tWhite_score_mean\tNb_black_hit\tBlack_score_mean\tNb_no_hit\n";


foreach my $cluster (keys %clusterhash) {
	my $white_hit = 0;
	my $white_score = 0;
	my $black_hit = 0;
	my $black_score = 0;
	my $wr = 0;
	my $no_hit = 0;
	my $beta = 0;
	my @cluster_contigs = @{ $clusterhash{$cluster} };
	foreach my $contig (@cluster_contigs){
		my $sw = 0;
		my $bw = 0;
		if (exists $whitehash{$contig}) {
			$sw = ${ whitehash{$contig} }[-1]; 
		}
		if (exists $blackhash{$contig}) { 
			$bw =  ${ blackhash{$contig} }[-1];
		}
		if ($sw > $bw) {
			$white_hit++;
			$white_score = $white_score + $sw; 
		} elsif($sw == 0 && $bw ==0){
			$no_hit++;
		} else{
			$black_hit++;
			$black_score = $black_score + $bw;
		}
	}
	my $total = $#cluster_contigs;

	$white_score = $white_score/($white_hit + 0.1);
	$black_score = $black_score/($black_hit + 0.1);
	
	#$wr = ($white_hit/$total) * $white_score;

	my $a = ( $white_score * ( $white_hit / ($total + 0.1) ));
	my $b = ( $black_score * ( $black_hit / ($total + 0.1) ));
	$wr = $a / ($a + $b + 0.1);

	#DEBUG
	#print "$cluster\t$total\t$white_hit\t$black_hit\t$both_hit\t$no_hit\n";
	
	print OUTFILE "$cluster\t$total\t$wr\t$white_hit\t$white_score\t$black_hit\t$black_score\t$no_hit\n";
}
close OUTFILE;

