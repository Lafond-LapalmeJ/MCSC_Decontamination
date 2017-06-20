#!/bin/perl

# Authors : Joel Lafond-Lapalme, Marc-Olivier Duceppe

use strict;
use warnings;
use File::Basename;
use Getopt::Long;
#use Data::Dumper;

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

# Get taxonomy level
my $tx = "";
if($taxo_level =~ m/kingdom/){$tx = "ki";}
elsif($taxo_level =~ m/phylum/){$tx = "ph";}
elsif($taxo_level =~ m/class/){$tx = "cl";}
elsif($taxo_level =~ m/order/){$tx = "or";}
elsif($taxo_level =~ m/family/){$tx = "fa";}
elsif($taxo_level =~ m/genus/){$tx = "ge";}
elsif($taxo_level =~ m/species/){$tx = "sp";}
else{die "Wrong taxonomic rank. Choices are: kingdom, phylum, class, order, family, genus, species. You enter $taxo_level";}

# Parse blast output file
my %whitehash;
open (my $blast_fh, '<', $blast_out) or die "Could not open $blast_out: $!\n";
# iterate line by line
while( my $line = <$blast_fh>)
{
    chomp($line); #remove carriage return at end of line
    $line =~ s/^\s+|\s+$//g;  # trim remove white space from both ends of a string
    next unless length($line);  # skip empty lines or lines with only whitespace characters

    #split line in array
    my($header, $bit, $taxid, $taxo) = split(/\t/, $line);
    (my $wanted_taxo) = $taxo =~ m/($tx\_[A-Za-z]*[^;])/;  # extract the wanted taxonomy level
    next unless $wanted_taxo;  #skip if wanted taxo is not present (e.g. not resolved to that level)
#    $wanted_taxo = "N/A" unless $wanted_taxo;
    $wanted_taxo =~ s/$tx\_//;  # remove headin taxonomy identifier
    $whitehash{$header}{'bitscore'} = $bit;
    $whitehash{$header}{'taxid'} = $taxid;
    $whitehash{$header}{'taxonomy'} = $wanted_taxo;
}
close ($blast_fh);

# Some stats
# Number of taxi
print "The $taxo_level to keep: $white_name\n\n";
my @best_tax;
my %taxa_stats;
foreach my $seq (sort keys %whitehash)
{
    my $taxon =  $whitehash{$seq}{'taxonomy'};
    $taxa_stats{$taxon}++;
    push(@best_tax, $taxon) unless grep(/$taxon/, @best_tax);
}
my $nb_taxa = keys (%taxa_stats);
print ("Total taxa in query: $nb_taxa\n\n");

# Print top taxa and number of times found
my @tax;
my @top_ten;
if ($nb_taxa <= 10)
{
    print ("Number of contigs for each $taxo_level found:\n");
    foreach my $taxon (sort keys %taxa_stats)
    {
        my $count = $taxa_stats{$taxon};
        print ("    $taxon $count\n");
    }
}
else
{
    print ("Top 10 $taxo_level in $white_filename$ext:\n");
#    foreach my $taxon (sort { @{ $taxa_stats{$b} } cmp @{ $taxa_stats{$a} } } keys %taxa_stats)
    foreach my $taxon (sort { $taxa_stats{$b} <=> $taxa_stats{$a} } keys %taxa_stats)
    {
        push(@tax, $taxon);  # ordered by number of hits
    }
    @top_ten = @tax[0 .. 9];  # only keep the top 10
    foreach my $top_taxon (@top_ten)
    {
        my $count = $taxa_stats{$top_taxon};
        print ("    $top_taxon $count\n");
    }
}

# Fetch fasta files in folder
my @cluster_files;
opendir (my $DIR, $cluster_dir) or die "Cannot open $cluster_dir: $!\n";
@cluster_files = grep(/fasta/, readdir($DIR));
closedir($DIR);

# Parse clustered fasta files
my %clusterhash;
foreach my $file (@cluster_files)
{
    open (my $fasta_fh, '<', "$cluster_dir/$file") or die "Could not open $cluster_dir/$file: $!\n";
    while( my $line = <$fasta_fh>)
    {
        chomp($line); #remove carriage return at end of line
        $line =~ s/^\s+|\s+$//g;  # trim remove white space from both ends of a string
        next unless length($line);  # skip empty lines or lines with only whitespace characters

        if( $line =~ m/^>/ )
        {
            my $header = substr((split(/\t|;|\s/, $line))[0], 1); # keep first field and remove ">"
            push ( @{ $clusterhash{$file} }, $header);  # Add to hash
        }
    }
    close ($fasta_fh);
}

# Create and open output file
open (my $outfile_fh, '>', $outputName) or die "Can't write to \"$outputName\": $!\n";
open (my$whiteID_fh, '>', $whiteid) or die "Can't write to \"$whiteid\": $!\n";

# Autodetect numer of best_tax and create right amount of variables.
# If more than 10, just make 10.
my $formated_taxa;
if ($nb_taxa > 10)
{
   $formated_taxa = join("\tNb_", @top_ten);
}
else
{
    $formated_taxa = join("\tNb_", @best_tax);
}

# Print header
print {$outfile_fh} ("Cluster\tNb_contigs\tNb_hits\tWR\tNb_$formated_taxa\tother\tNb_no_hit\tBIT_score_sum_White\tBIT_score_sum_black\n");

# Compute the WR of each fasta file (cluster)
foreach my $cluster (keys %clusterhash)
{
    if($cluster eq "2014-SEQ-0414_pilon_cluster_00000010.fasta")
    {
        print ("");
    }

    my $white_hit = 0;
    my $white_score = 0;
    my $black_hit = 0;
    my $black_score = 0;
    my $nb_hit = 0;
    my $no_hit = 0;
    my $nb_other = 0;
    my $wr = 0;
    my %taxo_counts;
    $taxo_counts{'no_hit'} = 0;

    foreach my $contig ( @{ $clusterhash{$cluster} } )
    {
        if( !exists $whitehash{$contig} )
        {
            $no_hit++;
        }
        else
        {
            $nb_hit++;
            my $t = $whitehash{$contig}{'taxonomy'};

            # if more thant 10, make "other" available
            if ($nb_taxa <= 10)
            {
                $taxo_counts{$t}++;
            }
            else
            {
                if (grep ($t, @top_ten))
                {
                    $taxo_counts{$t}++;
                }
                else
                {
                    $nb_other++;
#                    $taxo_counts{'other'}++;
                }
            }

            #If it's the taxonomic group we're looking for...
            if($whitehash{$contig}{'taxonomy'} eq $white_name)
            {
                $white_hit++;
                $white_score += $whitehash{$contig}{'bitscore'};
                print {$whiteID_fh} ("$contig\n");
            }
            else
            {
                $black_hit++;
                $black_score += $whitehash{$contig}{'bitscore'};
            }
        }
    }

    # Compute white ratio for cluster
    if(($white_score + $black_score) > 0)
    {
        $wr = ($white_score/($white_score+$black_score));
    }

    # Add number of counts for each taxon for all the contigs in the cluster
    my $nb_contig = keys $clusterhash{$cluster};
    my @counts_to_print;
    if ($nb_taxa <= 10)
    {
        foreach my $taxon (@best_tax)
        {
            if (exists $taxo_counts{$taxon})
            {
                push (@counts_to_print, $taxo_counts{$taxon});
            }
            else
            {
                push(@counts_to_print, 0);
            }
        }
        my $formated_counts = join("\t", @counts_to_print);
        print {$outfile_fh} ("$cluster\t$nb_contig\t$nb_hit\t$wr\t$formated_counts\t0\t$no_hit\t$white_score\t$black_score\n");
    }
    else  # ($nb_taxa > 10)
    {
        foreach my $taxon (@top_ten)
        {
            if (exists $taxo_counts{$taxon})
            {
                push (@counts_to_print, $taxo_counts{$taxon});
            }
            else
            {
                push(@counts_to_print, 0);
            }
        }
        my $formated_counts = join("\t", @counts_to_print);
        print {$outfile_fh} ("$cluster\t$nb_contig\t$nb_hit\t$wr\t$formated_counts\t$nb_other\t$no_hit\t$white_score\t$black_score\n");
    }
}

close ($outfile_fh);
close ($whiteID_fh);
