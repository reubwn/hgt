#!/usr/bin/env perl

## author: reubwn May 2017

use strict;
use warnings;

use Getopt::Long;

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

## parse HGT_results file:
open (my $RESULTS, $resultsfile) or die "[ERROR] Cannot open $resultsfile: $!\n";
while (<$RESULTS>) {
  chomp;
  next if /^\#/;
  my @F = split (/\s+/, $_);
  $hgt_results{qr/$F[0]/} = { ##HoH, key= query name as regex; val= {hu, ai, CHS, etc}
    'hU'  => $F[3],
    'AI'  => $F[6],
    'CHS' => $F[10],
    'TAX' => $F[11]
  };
}
close $RESULTS;
print STDERR "[INFO] Number of queries: ".scalar(keys %query_names)."\n";

## parse GFF file:
open (my $GFF, $gfffile) or die "[ERROR] Cannot open $gfffile: $!\n";
while (<$GFF>) {
  chomp;
  my @F = split (/\s+/, $_);
  #next unless $F[2] =~ /mrna/i; ##only look at mRNAs...NOPE doesnt work for some files...
  ## in the GFF line, want to find the appropriate result from the %hgt_results hash...
  if ($F[8] =~ (keys %hgt_results)) {
    print "$F[0]\n";
  }
}

## parse HGT_results file:
# open (my $RESULTS, $resultsfile) or die "[ERROR] Cannot open $resultsfile: $!\n";
# LINE: while (<$RESULTS>) {
#
# }


print STDERR "[INFO] Finished on ".`date`."\n";
