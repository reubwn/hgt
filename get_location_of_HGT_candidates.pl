#!/usr/bin/env perl

## author: reubwn May 2017

use strict;
use warnings;

use Getopt::Long;
use Data::Dumper qw(Dumper);

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

my ($infile,$gfffile,$prefix,$help);
my $CHS = 0.9;
my $hU = 30;
my $not = 0;
my $heavy = 0.75;

GetOptions (
  'in|i=s'      => \$infile,
  'gff|g=s'     => \$gfffile,
  'CHS|s:f'     => \$CHS,
  'hU|p:i'      => \$hU,
  'not|U:i'     => \$not,
  'heavy|V:f'   => \$heavy,
  'help|h'      => \$help,
);

die $usage if $help;
die $usage unless ($infile && $gfffile);

my $n = 1;
(my $insize = `wc -l $infile`) =~ s/\s.+//;
print STDERR "[INFO] Size of GFF: $insize\n";
my (%query_names,%hgt_results,%scaffolds,%gff);

## parse HGT_results file:
open (my $RESULTS, $infile) or die "[ERROR] Cannot open $infile: $!\n";
while (<$RESULTS>) {
  chomp;
  next if /^\#/;
  my @F = split (/\s+/, $_);
  my @gffline = split(/\s+/, `grep -m 1 -F $F[0] $gfffile`); ##grep 1st line from GFF containing query name
  print STDERR "[WARN] No scaffold name found for query $F[0]: something is wrong!\n" if (scalar(@gffline)==0);
  #my $scaffold = $gffline[0]; ##the scaffold will be the first element in @gffline, assuming a normal GFF

  $hgt_results{$F[0]} = { ##HoH, key= query name as regex; val= {hu, ai, CHS, etc}
    'hU'       => $F[3],
    'AI'       => $F[6],
    'CHS'      => $F[10],
    'taxonomy' => $F[11],
    'scaffold' => $gffline[0] ##the scaffold will be the first element in @gffline, assuming a normal GFF
  };
  print STDERR "\r[INFO] Working on query \#$n: $F[0] (".percentage($n,$insize)."\%)"; $|=1;
  $n++;
  last if percentage($n,$insize) == 10;
}
close $RESULTS;
print STDERR "[INFO] Number of queries: ".scalar(keys %hgt_results)."\n";

## parse GFF file:
# open (my $GFF, $gfffile) or die "[ERROR] Cannot open $gfffile: $!\n";
# GFF: while (<$GFF>) {
#   chomp;
#   next if /^\#/;
#   print "\r[INFO] Complete: ".(($n/$gffsize)*100);
#   my @F = split (/\s+/, $_);
#   #next unless $F[2] =~ /mrna/i; ##only look at mRNAs...NOPE doesnt work for some files...
#   ## in the GFF line, want to find the appropriate result from the %hgt_results hash...
#   QUERY: foreach my $query (keys %hgt_results) {
#     if (index($F[8], $query) >= 0) {
#       $hgt_results{$query}{'scaffold'} = $F[0];
#       next QUERY;
#     } else {
#       next GFF;
#     }
#   }
#   $n++;
# }

## parse HGT_results file:
# open (my $RESULTS, $resultsfile) or die "[ERROR] Cannot open $resultsfile: $!\n";
# LINE: while (<$RESULTS>) {
#
# }

print Dumper \%hgt_results;


print STDERR "[INFO] Finished on ".`date`."\n";

################################################################################

sub percentage {
    my $numerator = $_[0];
    my $denominator = $_[1];
    my $places = "\%.2f"; ## default is two decimal places
    if (exists $_[2]){$places = "\%.".$_[2]."f";};
    my $float = (($numerator / $denominator)*100);
    my $rounded = sprintf("$places",$float);
    return $rounded;
}

sub commify {
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text;
}
