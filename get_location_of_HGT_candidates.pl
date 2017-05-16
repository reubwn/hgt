#!/usr/bin/env perl

## author: reubwn May 2017

use strict;
use warnings;

use Getopt::Long;
use Sort::Naturally;
use Data::Dumper qw(Dumper);

my $usage = "
SYNOPSIS
  Takes a '*.HGT_results' file and a GFF and returns an '*.HGT_locations' file, specifying the location on each chromosome of HGT candidates.

OPTIONS:
  -i|--in     [FILE]  : *.HGT_results.txt file [required]
  -g|--gff    [FILE]  : GFF file [required]
  -u|--outgrp [INT]   : threshold hU score for determining 'good' OUTGROUP (HGT) genes [default>=30]
  -U|--ingrp  [INT]   : threshold hU score for determining 'good' INGROUP genes [default<=0]
  -c|--CHS    [INT]   : threshold CHS score for determining 'good' OUTGROUP (HGT) genes [default>=90\%]
  -y|--heavy  [INT]   : threshold for determining 'HGT heavy' scaffolds, with >= this proportion genes >= hU [default>=75\%]
  -b|--bed            : also write bed file for 'good' HGT genes (eg, for intersection with RNASeq bamfile)
  -s|--subset [INT]   : test on a subset of proteins [for debug]
  -h|--help           : prints this help message

OUTPUTS
  A '*.HGT_locations' file and a map file, and a list of scaffolds with 'too many' HGT candidates encoded on them to be believable ('*.HGT_heavy').
\n";

my ($infile,$namesfile,$gfffile,$bed,$subset,$help);
my $outgrp_threshold = 30;
my $ingrp_threshold = 0;
my $CHS_threshold = 90;
my $heavy = 75;

GetOptions (
  'i|in=s'     => \$infile,
  'n|names=s'  => \$namesfile,
  'g|gff=s'    => \$gfffile,
  'u|outgrp:i' => \$outgrp_threshold,
  'U|ingrp:i'  => \$ingrp_threshold,
  'c|CHS:i'    => \$CHS_threshold,
  'y|heavy:i'  => \$heavy,
  'b|bed'      => \$bed,
  's|subset:i' => \$subset,
  'h|help'     => \$help,
);

die $usage if $help;
die $usage unless ($infile && $gfffile);

print STDERR "[INFO] Infile: $infile\n";
print STDERR "[INFO] GFF file: $gfffile\n";
print STDERR "[INFO] Proteins names file: $namesfile\n";
print STDERR "[INFO] hU threshold to determine strong evidence for OUTGROUP: >= $outgrp_threshold\n";
print STDERR "[INFO] hU threshold to determine strong evidence for INGROUP: <= $ingrp_threshold\n";
print STDERR "[INFO] CHS threshold to determine strong evidence for OUTGROUP: >= $CHS_threshold\%\n";
print STDERR "[INFO] Proportion of genes >= hU threshold to determine 'HGT heavy' scaffolds: $heavy\n";
print STDERR "[INFO] Write bedfile: TRUE\n" if $bed;

(my $locationsfile = $infile) =~ s/HGT_results.+/HGT_locations.txt/;
(my $summaryfile = $infile) =~ s/HGT_results.+/HGT_locations.summary.txt/;
(my $heavyfile = $infile) =~ s/HGT_results.+/HGT_locations.heavy.txt/;
(my $bedfile = $infile) =~ s/HGT_results.+/HGT_locations.OUTGROUP.bed/ if $bed;
my ($namesfilesize,$isfasta,%bed,%query_names,%hgt_results,%scaffolds,%gff,%seen);
my $n=1;

## autodetect if names are coming from fasta:
if ($namesfile =~ m/(fa|faa|fasta)$/) {
  print STDERR "[INFO] Proteins names file is fasta...\n";
  $namesfilesize = `grep -c ">" $namesfile`;
  $namesfilesize =~ s/\s.+\n//;
  $isfasta = 1;
} else {
  $namesfilesize = `wc -l $namesfile`;
  $namesfilesize =~ s/\s.+\n//;
}

## grep protein names from GFF and get coords of CDS:
open (my $NAMES, $namesfile) or die "[ERROR] Cannot open $namesfile: $!\n";
print STDERR "[INFO] Getting genomic coordinates of proteins from GFF file...\n";
while (my $gene = <$NAMES>) {
  chomp $gene;
  print STDERR "\r[INFO] Working on query \#$n: $gene (".percentage($n,$namesfilesize)."\%)"; $|=1;

  my ($start,$end,$chrom,$strand) = (1e+9,0,"NULL","NULL"); ##this will work so long as no start coord is ever >=1Gb!

  ## get coords of all items grepped by $gene
  open (my $G, "grep -F $gene $gfffile |") or die "$!\n";
  while (<$G>) {
    chomp;
    my @F = split (/\s+/, $_);
    $start = $F[3] if $F[3] < $start; ##then get ONLY the 1st
    $end = $F[4] if $F[4] > $end; ##... and last coords across all items
    $chrom = $F[0];
    $strand = $F[6];
    $bed{$gene} = { ##key= gene; val= HoH
                    'chrom'  => $chrom,
                    'start'  => $start, ##this should cover the 'gene region'
                    'end'    => $end, ##... encoded by the protein name
                    'strand' => $strand
                   };
  }
  close $G;

  ## build scaffolds hash:
  push ( @{ $scaffolds{$chrom} }, $gene ); ##key= scaffold; val= \@array of genes on that scaffold

  $n++;
  if ($subset) {last if $n == $subset};
}
close $NAMES;
print STDERR "\n";

## parse HGT_results file:
print STDERR "[INFO] Parsing HGT_results file...";
open (my $RESULTS, $infile) or die "[ERROR] Cannot open $infile: $!\n";
while (<$RESULTS>) {
  chomp;
  next if /^\#/;
  my @F = split (/\s+/);

  ## evaluate HGT evidence.
  ## assign a 'score' based on hU, bbsumcat and CHS,
  ## 0: hU <= 0; bbsumcat = INGROUP; CHS >= 90 [good INGROUP gene]
  ## 1: 30 < hU > 0 [intermediate score]
  ## 2: hU >= 30; bbsumcat = OUTGROUP; CHS >= 90 [good HGT candidate]
  my $HGT_evidence;
  if ($F[3] <= $ingrp_threshold && $F[9] eq "INGROUP" && $F[10] >= $CHS_threshold) {
    $HGT_evidence = 0; ##good INGROUP candidate
  } elsif ($F[3] >= $outgrp_threshold && $F[9] eq "OUTGROUP" && $F[10] >= $CHS_threshold) {
    $HGT_evidence = 2; ##good HGT candidate
  } else {
    $HGT_evidence = 1; ##intermediate
  }

  ## build HGT_results hash:
  $hgt_results{$F[0]} = { ##key= query; val= HoH
    'hU'       => $F[3],
    'AI'       => $F[6],
    'bbsumcat' => $F[9],
    'CHS'      => $F[10],
    'taxonomy' => $F[11],
    'evidence' => $HGT_evidence
  };
}
close $RESULTS;
print STDERR " found ".scalar(keys %hgt_results)." queries\n";
print STDERR "[INFO] Evaluating results...\n";
# $n=0;

## iterate through pseudo-GFF:
open (my $LOC, ">$locationsfile") or die "[ERROR] Cannot open file $locationsfile: $!\n";
open (my $SUM, ">$summaryfile") or die "[ERROR] Cannot open file $summaryfile: $!\n";
open (my $HEV, ">$heavyfile") or die "[ERROR] Cannot open file $heavyfile: $!\n";
print $SUM join ("\t", "SCAFFOLD","NUMGENES","UNASSIGNED","GOOD_INGRP","INTERMEDIATE","GOOD_OUTGRP","PROPORTION_OUTGRP","IS_LINKED","\n");
print $HEV join ("\t", "SCAFFOLD","NUMGENES","UNASSIGNED","GOOD_INGRP","INTERMEDIATE","GOOD_OUTGRP","PROPORTION_OUTGRP","IS_LINKED","\n");
my ($good_outgrp_total,$good_ingrp_total,$intermediate_total,$na_total,$is_linked_total,$is_heavy) = (0,0);

## iterate across scaffolds:
foreach my $chrom (nsort keys %scaffolds) {
  # print STDERR "\r[INFO] Working on scaffold \#$n: $chrom (".percentage($n,scalar(keys %scaffolds))."\%)"; $|=1;
  print $LOC "## Scaffold \#$n: $chrom\n";

  ## sort by start coord within the %bed hash:
  my ($good_outgrp,$good_ingrp,$intermediate,$na,$is_linked) = (0,0,0,0,0);
  foreach my $gene ( sort {$bed{$a}{start}<=>$bed{$b}{start}} @{$scaffolds{$chrom}} ) {
    if ( exists($hgt_results{$gene}{hU}) ) {
      if ($hgt_results{$gene}{evidence} == 2) {
        print $LOC join ("\t", $chrom,$bed{$gene}{start},$bed{$gene}{end},$gene,".",$bed{$gene}{strand},$hgt_results{$gene}{hU},$hgt_results{$gene}{evidence},$hgt_results{$gene}{taxonomy},"\n");
        $good_outgrp++;
      } else {
        print $LOC join ("\t", $chrom,$bed{$gene}{start},$bed{$gene}{end},$gene,".",$bed{$gene}{strand},$hgt_results{$gene}{hU},$hgt_results{$gene}{evidence},"\n");
        $good_ingrp++ if $hgt_results{$gene}{evidence} == 0;
        $intermediate++ if $hgt_results{$gene}{evidence} == 1;
      }
    } else {
      ## NA if either no hit or if hit to 'skipped' taxon (usually self-phylum):
      print $LOC join ("\t", $chrom,$bed{$gene}{start},$bed{$gene}{end},$gene,".",$bed{$gene}{strand},"NA","\n");
      $na++;
    }
  }

  ## sum for totals:
  $good_ingrp_total += $good_ingrp;
  $good_outgrp_total += $good_outgrp;
  $intermediate_total += $intermediate;
  $na_total += $na;

  ## evaluate if HGT candidate gene is encoded on a scaffold which also encodes a 'good_ingrp' gene:
  $is_linked = 1 if ($good_outgrp > 0 && $good_ingrp > 0); ##must have at least one strong evidence for both on the same scaffold
  print $SUM join ("\t", $chrom,scalar(@{$scaffolds{$chrom}}),$na,$good_ingrp,$intermediate,$good_outgrp,(percentage($good_outgrp,scalar(@{$scaffolds{$chrom}}))),$is_linked,"\n");
  $is_linked_total += $is_linked;

  ## evaluate proportion of HGT candidates per scaffold; print to 'heavy' if > threshold:
  if ( (percentage($good_outgrp,scalar(@{$scaffolds{$chrom}}))) > $heavy ) {
    print $HEV join ("\t", $chrom,scalar(@{$scaffolds{$chrom}}),$na,$good_ingrp,$intermediate,$good_outgrp,(percentage($good_outgrp,scalar(@{$scaffolds{$chrom}}))),$is_linked,"\n");
    $is_heavy++;
  }
$n++;
}
close $LOC;
close $SUM;
close $HEV;
print STDERR "\n";
print STDERR "[INFO] Number of good INGROUP genes: $good_ingrp_total\n";
print STDERR "[INFO] Number of good OUTGROUP genes (HGT candidate): $good_outgrp_total\n";
print STDERR "[INFO] Number of genes with intermediate score: $intermediate_total\n";
print STDERR "[INFO] Number of genes with no assignment (no-hitters or skipped-hitters): $na_total\n";
print STDERR "[INFO] Number of scaffolds with HGT proportion >= $heavy: $is_heavy\n";
print STDERR "[INFO] Number of 'good' HGT candidates encoded on same scaffold as 'good' INGROUP gene: $is_linked_total\n";
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
