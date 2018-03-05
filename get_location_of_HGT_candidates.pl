#!/usr/bin/env perl

## author: reubwn May 2017

use strict;
use warnings;

use Getopt::Long;
use Sort::Naturally;
use Data::Dumper qw(Dumper);

my $usage = "
SYNOPSIS
  Takes a '*.HGT_results' file and a GFF and returns '*.HGT_locations' files,
  specifying the location on each chromosome of HGT candidates.

OPTIONS:
  -i|--in     [FILE] : HGT_results.txt file [required]
  -g|--gff    [FILE] : GFF file [required] (accepts *.gz)
  -n|--names  [FILE] : names of proteins in GFF file, can be fasta of proteins used
  -r|--regex  [STR]  : optional regex to apply to seq headers if -n is a fasta file
  -u|--outgrp [INT]  : threshold hU score for determining 'good' OUTGROUP (HGT) genes [default>=30]
  -U|--ingrp  [INT]  : threshold hU score for determining 'good' INGROUP genes [default<=0]
  -c|--CHS    [INT]  : threshold CHS score for determining 'good' OUTGROUP (HGT) genes [default>=90\%]
  -y|--heavy  [INT]  : threshold for determining 'HGT heavy' scaffolds [default>=95\%]
  -b|--bed           : also write bed file for 'good' HGT genes
  -h|--help          : prints this help message

OUTPUTS
  (1) HGT_locations: reports gene-by-gene HGT results in pseudo-BED format. Gene
      positions on scaffolds inherited from gene coordinates in input GFF.
  (2) HGT_locations.summary: reports per-scaffold summary of number of HGT candidates
  (3) HGT_locations.heavy: reports scaffolds with a high proportion (dictated by --heavy)
      of HGT candidates
  (4) HGT_locations.bed: BED format file of HGT candidates. Useful for intersection
      with RNASeq mapping data.
\n";

my ($infile,$gfffile,$namesfile,$regexstr,$bed,$help);
my $outgrp_threshold = 30;
my $ingrp_threshold = 0;
my $CHS_threshold = 90;
my $heavy = 95;

GetOptions (
  'i|in=s'     => \$infile,
  'g|gff=s'    => \$gfffile,
  'n|names=s'  => \$namesfile,
  'r|regex:s'  => \$regexstr,
  'u|outgrp:i' => \$outgrp_threshold,
  'U|ingrp:i'  => \$ingrp_threshold,
  'c|CHS:i'    => \$CHS_threshold,
  'y|heavy:i'  => \$heavy,
  'b|bed'      => \$bed,
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
print STDERR "[INFO] Proportion of genes >= hU threshold to determine 'HGT heavy' scaffolds: $heavy\%\n";
print STDERR "[INFO] Write bedfile: TRUE\n" if ($bed);

(my $locationsfile = $infile) =~ s/HGT_results.+/HGT_locations.txt/;
(my $summaryfile = $infile) =~ s/HGT_results.+/HGT_locations.scaffold_summary.txt/;
(my $oversummaryfile = $infile) =~ s/HGT_results.+/HGT_locations.overall_summary.txt/;
(my $heavyfile = $infile) =~ s/HGT_results.+/HGT_locations.heavy.txt/;
(my $warningsfile = $infile) =~ s/HGT_results.+/HGT_locations.warnings.txt/;
(my $bedfile = $infile) =~ s/HGT_results.+/HGT_locations.bed/ if ($bed);
my ($namesfilesize,%orphans,%locations,%hgt_results,%scaffolds);
my $n=1;

## grep protein names from GFF and get coords of CDS:
if ($namesfile =~ m/(fa|faa|fasta)$/) { ##autodetect if names are coming from fasta
  ## get filesize:
  $namesfilesize = `grep -c ">" $namesfile`;
  chomp $namesfilesize;

  open (my $FAA, $namesfile) or die "[ERROR] Cannot open $namesfile: $!\n";
  print STDERR "[INFO] Getting genomic coordinates of proteins from GFF file...\n";
  print STDERR "[INFO] Proteins names file is from fasta (".commify($namesfilesize)." sequences)\n";

  my $regexvar = qr/$regexstr/ if ($regexstr);
  print STDERR "[INFO] Applying regex 's/$regexstr//' to fasta headers...\n" if ($regexstr);

  while (my $gene = <$FAA>) {
    if ($gene =~ m/^>/) {
      chomp $gene;
      $gene =~ s/^>//;
      $gene =~ s/$regexvar//ig if ($regexstr); ##apply regex if specified
      print STDERR "\r[INFO] Working on query \#$n: $gene (".percentage($n,$namesfilesize)."\%)"; $|=1;
      my ($start,$end,$introns) = (1e+9,0,-1); ##this will work so long as no start coord is ever >=1Gb!
      my ($chrom,$strand) = ("NULL","NULL");
      ## get coords of all items grepped by $gene
      my $G;
      if ($gfffile =~ m/gz$/) {
        open (my $G, "zcat $gfffile | grep -F CDS | grep -F \Q$gene\E |") or die "$!\n"; ##will return GFF lines matching "CDS" && $gene
      } else {
        open (my $G, "grep -F CDS $gfffile | grep -F \Q$gene\E |") or die "$!\n"; ##will return GFF lines matching "CDS" && $gene
      }
      while (<$G>) {
        chomp;
        my @F = split (/\s+/, $_);
        $start = $F[3] if $F[3] < $start; ##then get ONLY the 1st
        $end = $F[4] if $F[4] > $end; ##... and last coords across all CDS
        $introns++; ##the number of iterations of through <$G> corresponds to the num exons; therefore introns is -1 this
        $chrom = $F[0];
        $strand = $F[6];
        $locations{$gene} = { ##key= gene; val= HoH
                        'chrom'   => $chrom,
                        'start'   => $start, ##this should cover the 'gene region'
                        'end'     => $end, ##... encoded by the protein name
                        'strand'  => $strand,
                        'introns' => $introns
                       };
      }
      close $G;
      unless (exists($locations{$gene}{chrom})) { ##all genes should have a chrom...
        print STDERR "\n[WARN] Nothing found for gene $gene in GFF $gfffile!\n";
        $orphans{$gene} = ();
      }

      ## build scaffolds hash:
      push ( @{ $scaffolds{$chrom} }, $gene ); ##key= scaffold; val= \@array of genes on that scaffold
      $n++;
    } else {
      next;
    }
  }
  close $FAA;
} else {
  ## get filesize:
  $namesfilesize = `wc -l $namesfile`;
  $namesfilesize =~ s/\s.+\n//;

  open (my $NAMES, $namesfile) or die "[ERROR] Cannot open $namesfile: $!\n";
  print STDERR "[INFO] Getting genomic coordinates of proteins from GFF file...\n";

  ## iterate through genes in namesfile:
  GENE: while (my $gene = <$NAMES>) {
    chomp $gene;
    print STDERR "\r[INFO] Working on query \#$n: $gene (".percentage($n,$namesfilesize)."\%)"; $|=1;
    my ($start,$end,$introns) = (1e+9,0,-1); ##this will work so long as no start coord is ever >=1Gb!
    my ($chrom,$strand) = ("NULL","NULL");
    ## get coords of all items grepped by $gene
    open (my $G, "grep -F CDS $gfffile | grep -F \Q$gene\E |") or die "$!\n"; ##will return GFF lines matching "CDS" && $gene
    while (<$G>) {
      chomp;
      my @F = split (/\s+/, $_);
      $start = $F[3] if $F[3] < $start; ##then get ONLY the 1st
      $end = $F[4] if $F[4] > $end; ##... and last coords across all CDS
      $introns++; ##the number of iterations of through <$G> corresponds to the num exons; therefore introns is -1 this
      $chrom = $F[0];
      $strand = $F[6];
      $locations{$gene} = { ##key= gene; val= HoH
                      'chrom'   => $chrom,
                      'start'   => $start, ##this should cover the 'gene region'
                      'end'     => $end, ##... encoded by the protein name
                      'strand'  => $strand,
                      'introns' => $introns
                     };
    }
    close $G;

    ## evaluate GFF grep results - gene name may not be present if there is a mismatch in protein faa and GFF used:
    unless (exists($locations{$gene}{chrom})) { ##all genes should have a chrom...
      print STDERR "\n[WARN] Nothing found for gene $gene in GFF $gfffile!\n";
      $orphans{$gene} = (); ##count
      next GENE; ##move on
    }

    ## build scaffolds hash:
    push ( @{ $scaffolds{$chrom} }, $gene ); ##key= scaffold; val= \@array of genes on that scaffold
    $n++;
  }
  close $NAMES;
}
print STDERR "\n";

## print orphans to warnings file:
if (scalar(keys %orphans) > 0) {
  print STDERR "[WARN] Orphans found! ".scalar(keys %orphans)." gene names not in GFF\n";
  open (my $WARN, ">$warningsfile") or die "[ERROR] Cannot open file $warningsfile: $!\n";
  print $WARN "Following gene names not found in $gfffile (".scalar(keys %orphans)." genes):\n";
  foreach (nsort keys %orphans) {
    print $WARN "$_\n";
  }
  close $WARN;
}

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
print STDERR " found ".commify(scalar(keys %hgt_results))." queries\n";
print STDERR "[INFO] Evaluating results...\n";
 $n=1;

## iterate through pseudo-GFF:
open (my $LOC, ">$locationsfile") or die "[ERROR] Cannot open file $locationsfile: $!\n";
open (my $SUM, ">$summaryfile") or die "[ERROR] Cannot open file $summaryfile: $!\n";
open (my $HEV, ">$heavyfile") or die "[ERROR] Cannot open file $heavyfile: $!\n";
open (my $BED, ">$bedfile") or die "[ERROR] Cannot open file $bedfile: $!\n" if ($bed);
#print $LOC "## HGT_locations\n##".`date`."\n";
print $LOC join ("\t", "#SCAFFOLD","START","END","GENE","SCORE","STRAND","INTRONS","hU","EVIDENCE","TAXONOMY","\n");
print $SUM join ("\t", "#SCAFFOLD","NUMGENES","UNASSIGNED","GOOD_INGRP","INTERMEDIATE","GOOD_OUTGRP","PROPORTION_OUTGRP","IS_LINKED","\n");
print $HEV join ("\t", "#SCAFFOLD","NUMGENES","UNASSIGNED","GOOD_INGRP","INTERMEDIATE","GOOD_OUTGRP","PROPORTION_OUTGRP","IS_LINKED","\n");
my ($good_outgrp_total,$good_ingrp_total,$intermediate_total,$na_total,$intronized,$linked_total,$is_heavy,$num_genes_on_heavy_HGT,$num_genes_on_heavy_total) = (0,0,0,0,0,0,0,0,0);

## iterate across scaffolds:
foreach my $chrom (nsort keys %scaffolds) {
  #print STDERR "\r[INFO] Working on scaffold \#$n: $chrom (".percentage($n,scalar(keys %scaffolds))."\%)"; $|=1;
  #print $LOC "## Scaffold \#$n: $chrom\n";

  ## sort by start coord within the %locations hash:
  my ($good_outgrp,$good_ingrp,$intermediate,$na,$is_linked) = (0,0,0,0,0);
  foreach my $gene ( sort {$locations{$a}{start}<=>$locations{$b}{start}} @{$scaffolds{$chrom}} ) {
    ## evaluate if protein name was not found in the GFF (if protein file and GFF dont match up exactly):
    # unless (exists($locations{$gene}{chrom})) {
    #   open (my $WARN, ">$warningsfile") or die "[ERROR] Cannot open file $warningsfile: $!\n";
    #   print $WARN join ("\t", $gene,"No chrom found in GFF $gfffile","\n");
    #   next GENE; ##skip to next gene
    # }
    ## then evaluate if gene has associated hU:
    if ( exists($hgt_results{$gene}{hU}) ) {
      if ($hgt_results{$gene}{evidence} == 2) {
        print $LOC join ("\t", $chrom,$locations{$gene}{start},$locations{$gene}{end},$gene,".",$locations{$gene}{strand},$locations{$gene}{introns},$hgt_results{$gene}{hU},$hgt_results{$gene}{evidence},$hgt_results{$gene}{taxonomy},"\n");
        print $BED join ("\t", $chrom,$locations{$gene}{start},$locations{$gene}{end},$gene,".",$locations{$gene}{strand},"\n") if ($bed);
        $good_outgrp++;
        $intronized++ if $locations{$gene}{introns} > 0;
      } else {
        print $LOC join ("\t", $chrom,$locations{$gene}{start},$locations{$gene}{end},$gene,".",$locations{$gene}{strand},$locations{$gene}{introns},$hgt_results{$gene}{hU},$hgt_results{$gene}{evidence},"\n");
        $good_ingrp++ if $hgt_results{$gene}{evidence} == 0;
        $intermediate++ if $hgt_results{$gene}{evidence} == 1;
      }
    } else {
      ## NA if either no hit or if hit to 'skipped' taxon (usually self-phylum):
      print $LOC join ("\t", $chrom,$locations{$gene}{start},$locations{$gene}{end},$gene,".",$locations{$gene}{strand},$locations{$gene}{introns},"NA","\n");
      $na++;
    }
  }

  ## sum for totals:
  $good_ingrp_total += $good_ingrp;
  $good_outgrp_total += $good_outgrp;
  $intermediate_total += $intermediate;
  $na_total += $na;

  ## evaluate if HGT candidate gene is encoded on a scaffold which also encodes a 'good_ingrp' gene:
  if ( ($good_outgrp>0) && ($good_ingrp>0) ) { ##must have at least one strong evidence for both on the same scaffold
    print $SUM join ("\t", $chrom,scalar(@{$scaffolds{$chrom}}),$na,$good_ingrp,$intermediate,$good_outgrp,(percentage($good_outgrp,scalar(@{$scaffolds{$chrom}}))),"1","\n");
    $linked_total += $good_outgrp; ##sum total linked HGT candidates
    $is_linked = 1;
  } else {
    print $SUM join ("\t", $chrom,scalar(@{$scaffolds{$chrom}}),$na,$good_ingrp,$intermediate,$good_outgrp,(percentage($good_outgrp,scalar(@{$scaffolds{$chrom}}))),"0","\n");
  }
  # $is_linked = 1 if ($good_outgrp > 0 && $good_ingrp > 0);
  # print $SUM join ("\t", $chrom,scalar(@{$scaffolds{$chrom}}),$na,$good_ingrp,$intermediate,$good_outgrp,(percentage($good_outgrp,scalar(@{$scaffolds{$chrom}}))),$is_linked,"\n");
  # $is_linked_total += $is_linked;

  ## evaluate proportion of HGT candidates per scaffold; print to 'heavy' if > threshold:
  if ( (percentage($good_outgrp,scalar(@{$scaffolds{$chrom}}))) > $heavy ) {
    print $HEV join ("\t", $chrom,scalar(@{$scaffolds{$chrom}}),$na,$good_ingrp,$intermediate,$good_outgrp,(percentage($good_outgrp,scalar(@{$scaffolds{$chrom}}))),$is_linked,"\n");
    $num_genes_on_heavy_total += scalar(@{$scaffolds{$chrom}});
    $num_genes_on_heavy_HGT += scalar( grep { $hgt_results{$_}{evidence}==2 } @{$scaffolds{$chrom}} );
    $is_heavy++;
  }
$n++;
}
close $LOC;
close $SUM;
close $HEV;
close $BED if ($bed);

print STDERR "\n";
print STDERR "[RESULT] Number of good INGROUP genes: ".commify($good_ingrp_total)."\n";
print STDERR "[RESULT] Number of good OUTGROUP genes (HGT candidates): ".commify($good_outgrp_total)."\n";
print STDERR "[RESULT] Number of HGT candidates with (at least one) intron: ".commify($intronized)."\n";
print STDERR "[RESULT] Number of HGT candidates linked to good INGROUP gene: ".commify($linked_total)."\n";
print STDERR "[RESULT] Number of genes with intermediate score: ".commify($intermediate_total)."\n";
print STDERR "[RESULT] Number of genes with no assignment (no-hitters or hit-to-skippers): ".commify($na_total)."\n";
print STDERR "[RESULT] Number of 'HGT heavy' scaffolds: ".commify($is_heavy)."\n";
print STDERR "[RESULT] Number of genes on 'HGT heavy' scaffolds: ".commify($num_genes_on_heavy_total)." (total); ".commify($num_genes_on_heavy_HGT)." (HGT candidates)\n";
print STDERR "\n[INFO] Finished on ".`date`."\n";

open (my $OVER, ">$oversummaryfile") or die "[ERROR] Cannot open file $oversummaryfile: $!\n";
print $OVER "[INFO] Infile: $infile\n";
print $OVER "[INFO] GFF file: $gfffile\n";
print $OVER "[INFO] Proteins names file: $namesfile\n";
print $OVER "[RESULT] Number of good INGROUP genes: ".commify($good_ingrp_total)."\n";
print $OVER "[RESULT] Number of good OUTGROUP genes (HGT candidates): ".commify($good_outgrp_total)."\n";
print $OVER "[RESULT] Number of HGT candidates with (at least one) intron: ".commify($intronized)."\n";
print $OVER "[RESULT] Number of HGT candidates linked to good INGROUP gene: ".commify($linked_total)."\n";
print $OVER "[RESULT] Number of genes with intermediate score: ".commify($intermediate_total)."\n";
print $OVER "[RESULT] Number of genes with no assignment (no-hitters or hit-to-skippers): ".commify($na_total)."\n";
print $OVER "[RESULT] Number of 'HGT heavy' scaffolds: ".commify($is_heavy)."\n";
print $OVER "[RESULT] Number of genes on 'HGT heavy' scaffolds: ".commify($num_genes_on_heavy_total)." (total); ".commify($num_genes_on_heavy_HGT)." (HGT candidates)\n";
print $OVER "[INFO] Finished on ".`date`."\n";
close $OVER;

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
