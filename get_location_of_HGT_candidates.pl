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
  -i|--in     [FILE]  : *.HGT_results.txt file [required]
  -g|--gff    [FILE]  : GFF file [required]
  -u|--outgrp [INT]   : threshold hU score for determining 'good' OUTGROUP (HGT) genes [default>=30]
  -U|--ingrp  [INT]   : threshold hU score for determining 'good' INGROUP genes [default<=0]
  -y|--heavy  [FLOAT] : threshold for determining 'HGT heavy' scaffolds, with >= this proportion genes >= hU [default>=0.75]
  -h|--help           : prints this help message

OUTPUTS
  A '*.HGT_locations' file and a map file, and a list of scaffolds with 'too many' HGT candidates encoded on them to be believable ('*.HGT_heavy').
\n";

my ($infile,$namesfile,$gfffile,$prefix,$help);
my $outgrp = 30;
my $ingrp = 0;
my $heavy = 0.75;

GetOptions (
  'in|i=s'      => \$infile,
  'n|names=s' => \$namesfile,
  'gff|g=s'     => \$gfffile,
  'outgrp|u:i'      => \$outgrp,
  'ingrp|U:i'     => \$ingrp,
  'heavy|y:f'   => \$heavy,
  'help|h'      => \$help,
);

die $usage if $help;
die $usage unless ($infile && $gfffile);

print STDERR "[INFO] Infile: $infile\n";
print STDERR "[INFO] GFF file: $gfffile\n";
print STDERR "[INFO] Protein names file: $namesfile\n";

(my $infilesize = `wc -l $infile`) =~ s/\s.+\n//;
(my $namesfilesize = `wc -l $namesfile`) =~ s/\s.+\n//;
(my $bedfile = $infile) =~ s/HGT_results.+/HGT_locations.bed.txt/;
(my $locationsfile = $infile) =~ s/HGT_results.+/HGT_locations.txt/;
(my $summaryfile = $infile) =~ s/HGT_results.+/HGT_locations.summary.txt/;
(my $heavyfile = $infile) =~ s/HGT_results.+/HGT_locations.heavy.txt/;
my (%bed,%query_names,%hgt_results,%scaffolds,%gff,%seen);
my $n = 1;

## grep protein names from GFF and convert to BED:
open (my $NA, $namesfile) or die "[ERROR] Cannot open $namesfile: $!\n";
open (my $BED, ">$bedfile") or die "[ERROR] Cannot open $bedfile: $!\n";
while (my $gene = <$NA>) {
  chomp $gene;
  print STDERR "\r[INFO] Working on query \#$n: $gene (".percentage($n,$namesfilesize)."\%)"; $|=1;

  my ($start,$end,$chrom,$strand) = (1e+9,0,"NULL","NULL"); ##this will work so long as no start coord is ever >=1Gb!
  open (my $G, "grep -F $gene $gfffile |") or die "$!\n";
  while (<$G>) {
    chomp;
    my @F = split (/\s+/, $_);
    $start = $F[3] if $F[3] < $start;
    $end = $F[4] if $F[4] > $end;
    $chrom = $F[0];
    $strand = $F[6];
    $bed{$gene} = { ##key= gene; val= HoH
                    'chrom'  => $chrom,
                    'start'  => $start,
                    'end'    => $end,
                    'strand' => $strand
                   };
    # my @bedline = ("\t", $F[0],$start,$end,$gene,".",$F[6],"\n");
    # $bed{$F[0]} = \@bedline;
  }
  close $G;

  ## build scaffolds hash:
  push ( @{ $scaffolds{$chrom} }, $gene ); ##key= scaffold; val= \@array of genes on that scaffold

  $n++;
  last if $n == 50;
}
close $NA;
close $BED;
$n = 0;
#print Dumper \%bed;

## parse HGT_results file:
open (my $RESULTS, $infile) or die "[ERROR] Cannot open $infile: $!\n";
while (<$RESULTS>) {
  chomp;
  next if /^\#/;
  my @F = split (/\s+/);
  print STDERR "\r[INFO] Working on query \#$n: $F[0] (".percentage($n,$infilesize)."\%)"; $|=1;

#  my @gffline = split(/\s+/, `grep -m 1 -F $F[0] $gfffile`); ##grep 1st line from GFF containing query name
#  die "[ERROR] No scaffold name found for query $F[0]\n" if (scalar(@gffline)==0);

  ## build HGT_results hash:
  $hgt_results{$F[0]} = { ##key= query; val= HoH
  #  'scaffold' => $gffline[0], ##the scaffold will be the first element in @gffline, assuming a normal GFF
    'hU'       => $F[3],
    'AI'       => $F[6],
    'bbsumcat' => $F[9],
    'CHS'      => $F[10],
    'taxonomy' => $F[11]
  };

  ## build scaffolds hash:
  #push ( @{ $scaffolds{$gffline[0]} }, $F[0] ); ##key= scaffold; val= \@array of genes on that scaffold

  $n++;
  #last if percentage($n,$infilesize) == 0.1;
}
close $RESULTS;
print STDERR "\n";
print STDERR "[INFO] Number of queries: ".scalar(keys %hgt_results)."\n";
print STDERR "[INFO] Mapping results to GFF...\n";
$n=1;

## iterate through GFF:
open (my $LOC, ">$locationsfile") or die "[ERROR] Cannot open file $locationsfile: $!\n";
open (my $GFF, $gfffile) or die "[ERROR] Cannot open file $gfffile: $!\n";

foreach my $chrom (nsort keys %scaffolds) {
  foreach my $gene ( @{$scaffolds{$chrom}} ) { ##get genes on scaffold
    print "$chrom\t$gene\t$hgt_results{$gene}{hU}\n";
  }

}












# OUTER: while (<$GFF>) {
#   chomp;
#   my @F = split (/\s+/, $_);
#



  # INNER: foreach my $gene (nsort keys %hgt_results) {
  #   print STDERR "\r[INFO] Working on query \#$n: $gene"; $|=1;
  #   next OUTER if exists($seen{$gene});
  #
  #   ## iterate through all genes:
  #   if ( index($_, $gene)>=0 ) { ##search for substring match
  #     #unless (exists($seen{$gene})) {
  #       if ($hgt_results{$gene}{'hU'} >= $outgrp) {
  #         print $LOC join ("\t",$hgt_results{$gene}{'scaffold'},$gene,"OUTGROUP",$hgt_results{$gene}{'hU'},$hgt_results{$gene}{'AI'},$hgt_results{$gene}{'CHS'},$hgt_results{$gene}{'taxonomy'},"\n");
  #       } elsif ($hgt_results{$gene}{'hU'} <= $ingrp) {
  #         print $LOC join ("\t", $hgt_results{$gene}{'scaffold'},$gene,"INGROUP",$hgt_results{$gene}{'hU'},"\n");
  #       } else {
  #         print $LOC join ("\t", $hgt_results{$gene}{'scaffold'},$gene,"INTERMEDIATE",$hgt_results{$gene}{'hU'},"\n");
  #       }
  #
  #       $n++;
  #     #}
  #     $seen{$gene} = (); ##prevents printing again on another GFF line
  #     last INNER; ##quit the foreach loop
  #   }
  # }
# }
close $GFF;
close $LOC;

## iterate across all genes per each scaffold:
# foreach my $scaff (nsort keys %scaffolds) {
#   my @genes = @{ $scaffolds{$scaff} };
#   foreach my $gene (@genes) {
#     print join ("\t", $scaff, $gene, $hgt_results{$gene}{'hU'}, $hgt_results{$gene}{'AI'}, $hgt_results{$gene}{'CHS'});
#   }
# }
#
# print Dumper \%scaffolds;

print STDERR "\n";
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
