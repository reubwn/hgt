#!/usr/bin/env perl

## author: reubwn May 2017

use strict;
use warnings;

use Getopt::Long;
use Sort::Naturally;
#use Tie::Hash::Regex;
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
  'hU|u:i'      => \$hU,
  'not|U:i'     => \$not,
  'heavy|V:f'   => \$heavy,
  'help|h'      => \$help,
);

die $usage if $help;
die $usage unless ($infile && $gfffile);

print STDERR "[INFO] Infile: $infile\n";
print STDERR "[INFO] GFF file: $gfffile\n";

my $n = 1;
(my $filesize = `wc -l $infile`) =~ s/\s.+\n//;
my (%query_names,%hgt_results,%scaffolds,%gff);

## parse HGT_results file:
open (my $RESULTS, $infile) or die "[ERROR] Cannot open $infile: $!\n";
while (<$RESULTS>) {
  chomp;
  next if /^\#/;
  my @F = split (/\s+/, $_);
  print STDERR "\r[INFO] Working on query \#$n: $F[0] (".percentage($n,$filesize)."\%)"; $|=1;

  my @gffline = split(/\s+/, `grep -m 1 -F $F[0] $gfffile`); ##grep 1st line from GFF containing query name
  die "[ERROR] No scaffold name found for query $F[0]\n" if (scalar(@gffline)==0);

  ## build HGT_results hash:
  $hgt_results{$F[0]} = { ##key= query; val= HoH
    'scaffold' => $gffline[0], ##the scaffold will be the first element in @gffline, assuming a normal GFF
    'hU'       => $F[3],
    'AI'       => $F[6],
    'CHS'      => $F[10],
    'taxonomy' => $F[11]
  };

  ## build scaffolds hash:
  push ( @{ $scaffolds{$gffline[0]} }, $F[0] ); ##key= scaffold; val= \@array of genes on that scaffold
  #push ( @{ $scaffolds{$gffline[0]}{'coords'} }, $F[0] ); ##key= scaffold; val= \@array of genes on that scaffold

  $n++;
  last if percentage($n,$filesize) == 1;
}
close $RESULTS;
print STDERR "\n";
print STDERR "[INFO] Number of queries: ".scalar(keys %hgt_results)."\n";
print STDERR "[INFO] Sorting scaffolds...\n";

# tie %hgt_results, 'Tie::Hash::Regex';
open (my $GFF, $gfffile) or die "[ERROR] Cannot open file $gfffile: $!\n";
while my $line (<$GFF>) {
  SEARCH: while (my ($k,$v) = each %hgt_results) {
    if ($line =~ /$k/) {
      print join ("\t", $hgt_results{$k}{'scaffold'}, $k, $hgt_results{$k}{'hU'}, $hgt_results{$k}{'AI'}, "\n");
      last SEARCH;
    }
  }
}
close $GFF;

## iterate across all genes per each scaffold:
# foreach my $scaff (nsort keys %scaffolds) {
#   my @genes = @{ $scaffolds{$scaff} };
#   foreach my $gene (@genes) {
#     print join ("\t", $scaff, $gene, $hgt_results{$gene}{'hU'}, $hgt_results{$gene}{'AI'}, $hgt_results{$gene}{'CHS'});
#   }
# }
#
# print Dumper \%scaffolds;

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

__END__
