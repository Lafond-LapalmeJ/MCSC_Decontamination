#!/usr/bin/env perl

# Author: Georgios Koutsovoulos

use strict;
use warnings;

my $tool_name = "diamond";
my $tool_path = `which $tool_name`;
die "No $tool_name command available\n" unless ( $tool_path );

my $uniref_file=$ARGV[0];
my $diamond_file=$ARGV[1];

die "Usage: daa_to_tagc.pl uniref100.taxlist [assembly_se_uniref.daa]\n" unless ( $uniref_file && $diamond_file );

open (IN,"$uniref_file");

my %uniref_taxid;

while (<IN>) {

chomp;
my @columns=split(/\t/,$_);
my ($db,$unirefid)=split(/_/,$columns[0]);
$uniref_taxid{$unirefid}=$columns[1];

}

close IN;
open my $cmd, "diamond view -a $diamond_file |";
open OUT, ">$diamond_file.tagc";

while (<$cmd>){

chomp;
my @columns=split(/\t/,$_);
my ($db,$unirefid)=split(/_/,$columns[1]);

if (exists $uniref_taxid{$unirefid}) {
	print OUT join("\t", $columns[0],$uniref_taxid{$unirefid}, $columns[11]) . "\n";
	}

}

close OUT;
