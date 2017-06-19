#!/bin/perl

# Author : Joel Lafond-Lapalme

use strict;
use warnings;
use File::Basename;
use Getopt::Long;


############ FILEs validation ##################
die "Usage: perl cluster_eval.pl <cluster_directory> <blast_out> <taxo_level> <white_name>\n" if @ARGV != 4; 

################################################


# Get file
my $cluster_dir = shift;
my $blast_out = shift;
my $taxo_level = shift;
my $white_name = shift;


# Get the directory, file name and extention
my ( $white_filename, $out_dir, $ext ) = fileparse( $blast_out, qr{ \. [^.]+ \z }msx );
my $outputName = "$out_dir" . "cluster_eval.tsv";
my $whiteid = "$out_dir" . "white_id.txt";
my %whitehash;
my %clusterhash;




my $tx = "";
if($taxo_level =~ /kingdom/){$tx = "ki";}
elsif($taxo_level =~ /phylum/){$tx = "ph";}
elsif($taxo_level =~ /class/){$tx = "cl";}
elsif($taxo_level =~ /order/){$tx = "or";}
elsif($taxo_level =~ /family/){$tx = "fa";}
elsif($taxo_level =~ /genus/){$tx = "ge";}
elsif($taxo_level =~ /species/){$tx = "sp";}
else{die "Wrong taxonomic rank. Choices are kingdom, phylum, class, order, family, genus, species. You enter $taxo_level";}

print "The $taxo_level to keep is $white_name\n";
print "Top $taxo_level in $blast_out :\n";
my $cmd1 = "grep -o ";
my $cmd2 = "\"$tx" . "_[A-Za-z]*;\"";
my $cmd3 = " $blast_out | sort | uniq -c | sort -rnk1,1 | sed ";
my $cmd4 = "\"s/\\s*[0-9]*\\s*$tx" . "_//\"";
my $cmd5 = " | head -n 20 | grep -o \"[A-Za-z]*\"";
system("$cmd1$cmd2 $cmd3$cmd4$cmd5");

# Open blast output file
open my $blast, $blast_out or die "Could not open $blast_out: $!";
# iterate line by line
while( my $line = <$blast>)  {   
    chomp $line;  
    
    #split line in array  
    my @lineArray = split(/\t/, $line);
        my $str = shift @lineArray;
        $whitehash{$str} = join("\t", @lineArray);

        #debug
        #print "$str\t$lineArray[2]\n";
    
}
close $blast;


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
open( WHITEID, '>', "$whiteid") or die "Can't create \"$whiteid\": $!\n";

my @best_tax = `$cmd1$cmd2$cmd3$cmd4$cmd5`;
#print "$best_tax[1]";

chomp @best_tax;



my $taxa1 = shift @best_tax;
my $taxa2 = shift @best_tax;
my $taxa3 = shift @best_tax;
my $taxa4 = shift @best_tax;
my $taxa5 = shift @best_tax;
my $taxa6 = shift @best_tax;
my $taxa7 = shift @best_tax;
my $taxa8 = shift @best_tax;
my $taxa9 = shift @best_tax;
my $taxa10 = shift @best_tax;


# Print header
print OUTFILE "Cluster\tNb_contigs\tNb_hits\tWR\tNb_$taxa1\tNb_$taxa2\tNb_$taxa3\tNb_$taxa4\tNb_$taxa5\tNb_$taxa6\tNb_$taxa7\tNb_$taxa8\tNb_$taxa9\tNb_$taxa10\tother\tNb_no_hit\tBIT_score_sum_White\tBIT_score_sum_black\n";


foreach my $cluster (keys %clusterhash) {
        my $white_hit = 0;
        my $white_score = 0;
        my $black_hit = 0;
        my $black_score = 0;
		my $nb_taxa1 = 0;
		my $nb_taxa2 = 0;
		my $nb_taxa3 = 0;
		my $nb_taxa4 = 0;
		my $nb_taxa5 = 0;
		my $nb_taxa6 = 0;
		my $nb_taxa7 = 0;
		my $nb_taxa8 = 0;
		my $nb_taxa9 = 0;
		my $nb_taxa10 = 0;
		my $nb_other = 0;
        my $nb_hit = 0;
        my $no_hit = 0;
        my $wr = 0;
		
        my @cluster_contigs = @{ $clusterhash{$cluster} };
        foreach my $contig (@cluster_contigs){

                if (!exists $whitehash{$contig}) {
						$no_hit++;
                }
				else{
					$nb_hit++;
					if ($whitehash{$contig} =~ /$taxa1/) { 
							$nb_taxa1++;
					}
					elsif ($whitehash{$contig} =~ /$taxa2/) { 
							$nb_taxa2++;
					}
					elsif ($whitehash{$contig} =~ /$taxa3/) { 
							$nb_taxa3++;
					}
					elsif ($whitehash{$contig} =~ /$taxa4/) { 
							$nb_taxa4++;
					}
					elsif ($whitehash{$contig} =~ /$taxa5/) { 
							$nb_taxa5++;
					}
					elsif ($whitehash{$contig} =~ /$taxa6/) { 
					    		$nb_taxa6++;
					}
					elsif ($whitehash{$contig} =~ /$taxa7/) { 
					    		$nb_taxa7++;
					}
					elsif ($whitehash{$contig} =~ /$taxa8/) { 
					    		$nb_taxa8++;
					}
					elsif ($whitehash{$contig} =~ /$taxa9/) { 
					    		$nb_taxa9++;
					}
					elsif ($whitehash{$contig} =~ /$taxa10/) { 
					    		$nb_taxa10++;
					}
					else { 
							$nb_other++;
					}
                        my @lineArray = split(/\t/, $whitehash{$contig});
                        my $bit = $lineArray[1];	
						my $taxo = $lineArray[2];

					if($whitehash{$contig} =~ /$white_name/){
							$white_hit++;
							$white_score = ($bit + $white_score);
							print WHITEID "$contig\n";
					}
					else{
						$black_hit++;
						$black_score = ($bit + $black_score);
					}
				}

        }

	#debug
	#print "$white_score\n$white_hit\n$black_score\n$nb_hit\n$no_hit";
        if(($white_score+$black_score)>0){
		$wr = ($white_score/($white_score+$black_score));
	}
		
		my $nb_contig = $#cluster_contigs + 1;
		
		
        print OUTFILE "$cluster\t$nb_contig\t$nb_hit\t$wr\t$nb_taxa1\t$nb_taxa2\t$nb_taxa3\t$nb_taxa4\t$nb_taxa5\t$nb_taxa6\t$nb_taxa7\t$nb_taxa8\t$nb_taxa9\t$nb_taxa10\t$nb_other\t$no_hit\t$white_score\t$black_score\n";
}
close OUTFILE;
close WHITEID;
