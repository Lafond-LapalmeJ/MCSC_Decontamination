#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;
use Data::Dumper;

#I\O files
my $diamondBlastFile = $ARGV[0];
my $namesFile = $ARGV[1];
my $nodesFile = $ARGV[2];

if (scalar(@ARGV) != 3)
{
        die "Usage: perl get_taxo.pl <blastFile> <names.dmp> <nodes.dmp>";
}



#2-Create hash with the names file where the keys are the taxids and the values are the names.
# 11 seconds to run on Intel� Core i7-3820 CPU @ 3.60GHz � 8, 62.1 GB ram
open(NAMES, "<", $namesFile) or die "Cannot open $namesFile: $!";
my %names;
while (my $line = <NAMES>)
{
        chomp($line);
        my @elem = split(/\|/, $line);
        #clean whitespaces before and after string of first 2 fields (taxid)
        for (my $i = 0; $i <= 1; $i++)
        #foreach (@elem)
        {
                #$_ =~ s/^\s+|\s+$//g;
                $elem[$i] =~ s/^\s+|\s+$//g;
        }

        if (exists $names{$elem[0]})
        {
                #Only use the value of the key the first time it's encountered
                next;
        }
        else
        {
                $names{$elem[0]} = $elem[1];
        }
}

close(NAMES);


#3-Create hash with the nodes file where the keys are the taxids and the values are the parents.
# 6 seconds to run on Intel� Core i7-3820 CPU @ 3.60GHz � 8, 62.1 GB ram
open(NODES, "<", $nodesFile) or die "Cannot open $nodesFile: $!";
my %nodes;
my %rank;
while (my $line = <NODES>)
{
        chomp($line);
        my @elem = split(/\|/, $line);
        #clean whitespaces before and after string of first 2 fields (taxid)
        for (my $i = 0; $i <= 2; $i++)
        {
                $elem[$i] =~ s/^\s+|\s+$//g;
        }
        $nodes{$elem[0]} = $elem[1];
        $rank{$elem[0]} = $elem[2];
}

close(NODES);


#4-Create hash with the blast file where the keys are the taxids and the values are kingdom.
open(BLAST, "<", $diamondBlastFile) or die "Cannot open $diamondBlastFile: $!";
my %taxids;
while (my $line = <BLAST>)
{
        chomp($line);
        my @elem = split(/\t/, $line);
        my $id = (split(";",$elem[1]))[0];
        if($id eq "N/A"){
                $id = 1;
        }
        
        $taxids{$id} = "";

}
close(BLAST);


#5-Get the taxonomy for each taxid present in blast output
foreach my $key (sort keys %taxids)
{
        my $taxid = $key;
        my @taxo;
        while ($taxid > 1)
        {
                my $tx_name = defined $names{$taxid} ? $names{$taxid} : '';
                my $tx_rank = defined $rank{$taxid} ? substr($rank{$taxid},0,2) : '';


                $taxids{$key} = "_" . $tx_rank . "_" . $tx_name . ";$taxids{$key}"; 
                if (defined $nodes{$taxid}) #if taxid present in nodes.dmp
                {
                        $taxid = $nodes{$taxid}; #grab the parent
                }
                else    #if taxid not in nodes.dmp, must have been merged
                {
                        #no parent return the root
                        $taxid = 1;
                }
        }
}




#6 print name gi and taxonomy
open(BLAST, "<", $diamondBlastFile) or die "Cannot open $diamondBlastFile: $!";
while (my $line = <BLAST>)
{
        chomp($line);
        my @elem = split(/\t/, $line);
  
        print "$line\t$taxids{$elem[1]}\n";


}
close(BLAST);

