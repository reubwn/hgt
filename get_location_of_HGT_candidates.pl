#!/usr/bin/env perl

## author: reubwn May 2017

use strict;
use warnings;

my $usage = "
SYNOPSIS
  Takes a '*.HGT_results' file and a GFF and returns an '*.HGT_locations' file, specifying the location on each chromosome of HGT candidates.

OPTIONS:
  -i|--in      [FILE]   : taxified diamond/BLAST results file [required]
  -r|--results [FILE]   : *.HGT_results.txt file [required]
  -g|--gff     [FILE]   : GFF file [required]
  -s|--CHS     [FLOAT]  : Consesus Hits Support threshold [default>=0.9]
  -u|--hU      [INT]    : threshold for determining 'good' OUTGROUP (HGT) genes [default>=30]
  -U|--not     [INT]    : threshold for determining 'good' INGROUP genes [default<=0]
  -V|--heavy   [FLOAT]  : threshold for determining 'HGT heavy' scaffolds, with >= this proportion genes >= hU [default>=0.75]
  -h|--help             : prints this help message

OUTPUTS
  A '*.HGT_locations' file and a map file, and a list of scaffolds with 'too many' HGT candidates encoded on them to be believable ('*.HGT_heavy').
\n";

my ($infile,$resultsfile,$gfffile,$prefix,$help);
my $CHS = 0.9;
my $hU = 30;
my $not = 0;
my $heavy = 0.75;

GetOptions (
  'in|i=s'      => \$infile,
  'results|r=s' => \$resultsfile,
  'gff|g=s'     => \$gfffile,
  'CHS|s:f'     => \$CHS,
  'hU|p:i'      => \$hU,
  'not|U:i'     => \$not,
  'heavy|V:f'   => \$heavy,
  'help|h'      => \$help,
);

die $usage if $help;
die $usage unless ($infile && $resultsfile && $gfffile);

my $n = 1;
my (%query_names,%hgt_results,%saffolds,%gff);

## 1st parse HGT_results file to get query names:
open (my $RESULTS, $resultsfile) or die "[ERROR] Cannot open $resultsfile: $!\n";
while (<$RESULTS>) {
  chomp;
  next if /^\#/;
  my @F = split (/\s+/, $_);
  $query_names{$F[0]} = ();
  #$query_names{($F[0]=~s/.+\|//;)} = (); ##NB modification required!!! change as needed?
}
close $RESULTS;
print STDERR "[INFO] Number of queries: ".scalar(keys %query_names)."\n";

## parse GFF file:
open (my $GFF, $gfffile) or die "[ERROR] Cannot open $gfffile: $!\n";
while (<$GFF>) {
  chomp;
  my @F = split (/\s+/, $_);
  if ($F[8] =~ /ID\=(.+\|g\d+\.t\d+)[\.\;].+/)) {

  }
}

## parse HGT_results file:
open (my $RESULTS, $resultsfile) or die "[ERROR] Cannot open $resultsfile: $!\n";
LINE: while (<$RESULTS>) {

}


print STDERR "[INFO] Finished on ".`date`."\n";
