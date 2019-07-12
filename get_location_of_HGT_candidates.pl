#!/usr/bin/env perl

## author: reubwn May 2017

use strict;
use warnings;

use Getopt::Long;
use Term::ANSIColor;
use Sort::Naturally;
use Data::Dumper qw(Dumper);
use List::MoreUtils qw(indexes);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

my $usage = "
SYNOPSIS
  Takes a '*.HGT_results' file and a GFF and returns '*.HGT_locations' files,
  specifying the location on each chromosome of HGT candidates.

OPTIONS
  -i|--in       [FILE] : *.HGT_results.txt file [required]
  -g|--gff      [FILE] : annotation file in GFF3 format (supports gzipped) [required]
  -p|--prot     [FILE] : proteins file in fasta format (see below) (supports gzipped) [required]
  -r|--regexp   [STR]  : optional regex to apply to seq headers if -n is a fasta file
  -u|--outgrp   [INT]  : hU threshold to determine strong evidence for 'outgroup' [default>=30]
  -U|--ingrp    [INT]  : hU threshold to determine strong evidence for 'ingroup' [default<=0]
  -c|--CHS      [INT]  : CHS threshold to determine strong evidence for 'outgroup' [default>=90\%]
  -y|--heavy    [INT]  : Proportion of genes >= hU threshold to find contaminant scaffolds [default>=95\%]
  -b|--bed             : also write bed file for HGTc
  -h|--help            : this help message

NOTES
  Fasta headers in <PROTEINS> must match with names in <GFF>
  If they don't, use option '-r' to apply regexp to fasta headers

OUTPUTS
  (1) *.HGT_locations: reports gene-by-gene HGT results in pseudo-BED format; HGTc coords inherited from GFF
  (2) *.HGT_locations.summary: reports per-scaffold summary of number of HGTc
  (3) *.HGT_locations.heavy: reports scaffolds with a suspicously high \%HGTc (possible contaminant)
  (4) *.HGT_locations.bed: BED format file of HGT candidates; useful for intersection with RNASeq mapping data
\n";

my ($in_file,$gff_file,$proteins_file,$regexp_option,$bed,$help,$debug);
my $outgrp_threshold = 30;
my $ingrp_threshold = 0;
my $CHS_threshold = 90;
my $heavy = 95;
my $mapping = 0;

GetOptions (
  'i|in=s'     => \$in_file,
  'g|gff=s'    => \$gff_file,
  'p|prot=s'   => \$proteins_file,
  'r|regexp:s' => \$regexp_option,
  'u|outgrp:i' => \$outgrp_threshold,
  'U|ingrp:i'  => \$ingrp_threshold,
  'c|CHS:i'    => \$CHS_threshold,
  'y|heavy:i'  => \$heavy,
  'b|bed'      => \$bed,
  'h|help'     => \$help,
  'd|debug'    => \$debug
);

die $usage if $help;
die $usage unless ($in_file && $gff_file);

print STDERR "[INFO] Infile: '$in_file'\n";
print STDERR "[INFO] GFF file: '$gff_file'\n";
print STDERR "[INFO] Protein ID file: '$proteins_file'\n";
print STDERR "[INFO] hU threshold to determine strong evidence for OUTGROUP: >= ".colored($outgrp_threshold, 'yellow')."\n";
print STDERR "[INFO] hU threshold to determine strong evidence for INGROUP: <= ".colored($ingrp_threshold, 'yellow')."\n";
print STDERR "[INFO] CHS threshold to determine strong evidence for OUTGROUP: >= ".colored("$CHS_threshold\%", 'yellow')."\n";
print STDERR "[INFO] Proportion of genes >= hU threshold to find contaminant scaffolds: ".colored("$heavy\%", 'yellow')."\n";
print STDERR "[INFO] Applying regex '$regexp_option' to names in '$proteins_file'\n" if ( $regexp_option );
print STDERR "[INFO] Write bedfile: TRUE\n" if ( $bed );

## decompress gzipped input files if necessary
my ($gff_is_gz, $names_is_gz) = (0,0);
if ($gff_file =~ m/gz$/) {
  print STDERR "[INFO] Gunzipping '$gff_file' ";
  $gff_file = decompress ($gff_file); ## $in_file inherits new filename with '.gz' extension removed
  print STDERR "to '$gff_file'\n";
  $gff_is_gz = 1;
}
if ($proteins_file =~ m/gz$/) {
  print STDERR "[INFO] Gunzipping '$proteins_file' ";
  $proteins_file = decompress ($proteins_file); ## $in_file inherits new filename with '.gz' extension removed
  print STDERR "to '$proteins_file'\n";
  $names_is_gz = 1;
}

## create output filenames
(my $locations_file = $in_file) =~ s/HGT_results.+/HGT_locations.txt/;
(my $summary_file = $in_file) =~ s/HGT_results.+/HGT_locations.scaffold_summary.txt/;
(my $genes_file = $in_file) =~ s/HGT_results.+/HGT_locations.HGT_linked.genelist.txt/;
(my $oversummary_file = $in_file) =~ s/HGT_results.+/HGT_locations.overall_summary.txt/;
(my $heavy_file = $in_file) =~ s/HGT_results.+/HGT_locations.heavy.txt/;
(my $warnings_file = $in_file) =~ s/HGT_results.+/HGT_locations.warnings.txt/;
(my $bed_file = $in_file) =~ s/HGT_results.+/HGT_locations.bed/ if ($bed);
my ($proteins_total,%orphans,%locations,%hgt_results,%scaffolds,%names_map,%protein_hash,%protein_hash_map);
my $regexp = qr/$regexp_option/ if ($regexp_option); ## this should store the command line argument as compiled regexp
my $n=1;

## parse CDS entries from GFF into memory
my ($GFF_fh, @GFF_array);
print STDERR "[INFO] Parsing '$gff_file' GFF file...\n";
open ($GFF_fh, $gff_file) or die "$!\n";
while (my $line = <$GFF_fh>) {
  chomp ($line);
  my @F = split (m/\s+/, $line);
  push (@GFF_array, $line) if $F[2] eq "CDS"; ## push CDS entries ONLY
}
close $GFF_fh;

## grep protein names from GFF and get coords of CDS:
print STDERR "[INFO] Parsing '$proteins_file' fasta file...\n";
chomp ($proteins_total = `grep -c ">" $proteins_file`); ## get filesize
print STDERR "[INFO] Number of sequences in '$proteins_file': ".commify($proteins_total)."\n";
print STDERR "[INFO] Applying regex 's/$regexp_option//' to fasta headers...\n" if ($regexp_option);

my $PROT_fh;
open ($PROT_fh, $proteins_file) or die "[ERROR] Cannot open $proteins_file: $!\n";
while (my $line = <$PROT_fh>) {
  if ($line =~ m/^>/) {
    chomp (my $gene = $line);
    $gene =~ s/^>//; ## trim ">"
    $gene =~ s/\s.*//; ## also trim anything after 1st whitespace
    $gene =~ /$regexp/ if ($regexp_option); ## apply regex if specified

    print STDERR "\r[INFO] Working on query \#$n: $gene (".percentage($n,$proteins_total)."\%)"; $|=1;

    ## initial coords of all items grepped by $gene
    my ($start,$end,$introns) = (1e+12,0,-1); ## this will work so long as no start coord is ever >=1Tb!
    my ($chrom,$strand) = ("NULL","NULL");

    ## grep the relevant CDS entries based on value of $gene anywhere in the GFF line
    foreach my $line ( grep { m/\Q$gene\E\;*/ } @GFF_array ) { ## assumes $gene is bounded by a ';', or nothing if EOL
      chomp ($line);
      print STDOUT "Grepped for '$gene': $line\n" if ( $debug );

      my @F = split (/\s+/, $line);
      $start = $F[3] if $F[3] < $start; ## then get ONLY the 1st
      $end = $F[4] if $F[4] > $end; ##... and last coords across all CDS
      $introns++; ## number of iterations corresponds to the num exons; therefore introns is -1 this
      $chrom = $F[0];
      $strand = $F[6];
      $locations{$gene} = { ##key= gene; val= HoH
                      'chrom'   => $chrom,
                      'start'   => $start, ## this should cover the 'gene region'
                      'end'     => $end, ##... encoded by the protein name
                      'strand'  => $strand,
                      'introns' => $introns
                     };
    }

    print STDOUT Dumper (%locations) if ( $debug );

    ## dynamically shrink @GFF_array so search should get faster as parsing progresses?
    ## first get indices...
    my @to_splice = indexes { m/\Q$gene\E\;*/ } @GFF_array; ## assumes $gene is bounded by a ';', or nothing if EOL
    ## and splice them out of @GFF_array
    print STDOUT "Deleting ".scalar(@to_splice)." indices (@to_splice)\n" if ( $debug );
    print STDOUT "Length of \@GFF_array: ".scalar(@GFF_array)."\n" if ( $debug );
    splice (@GFF_array, $to_splice[0], scalar(@to_splice)); ## this should be OK, as CDS will always be consecutive in GFF

    ## old way of grepping directly from GFF using system grep
      # my $G;
      # if ($gff_file =~ m/gz$/) {
      #   open ($G, "zcat $gff_file | LC_ALL=C grep -F CDS | LC_ALL=C grep -F \Q$gene\E |") or die "$!\n"; ##will return GFF lines matching "CDS" && $gene
      #   print STDERR "[WARN] File '$gff_file' is gzipped, but decompressing first will be much faster\n";
      # } else {
      #   open ($G, "LC_ALL=C grep -F CDS $gff_file | LC_ALL=C grep -F \Q$gene\E |") or die "$!\n"; ##will return GFF lines matching "CDS" && $gene
      # }
      # while (<$G>) {
      #   chomp;
      #   my @F = split (/\s+/, $_);
      #   $start = $F[3] if $F[3] < $start; ##then get ONLY the 1st
      #   $end = $F[4] if $F[4] > $end; ##... and last coords across all CDS
      #   $introns++; ##the number of iterations of through <$G> corresponds to the num exons; therefore introns is -1 this
      #   $chrom = $F[0];
      #   $strand = $F[6];
      #   $locations{$gene} = { ##key= gene; val= HoH
      #                   'chrom'   => $chrom,
      #                   'start'   => $start, ##this should cover the 'gene region'
      #                   'end'     => $end, ##... encoded by the protein name
      #                   'strand'  => $strand,
      #                   'introns' => $introns
      #                  };
      # }
      # close $G;

    unless (exists($locations{$gene}{chrom})) { ##all genes should have a chrom...
      print STDERR "\n[WARN] Nothing found for gene $gene in GFF $gff_file!\n";
      $orphans{$gene} = ();
    }

    ## build scaffolds hash:
    push ( @{ $scaffolds{$chrom} }, $gene ); ##key= scaffold; val= \@array of genes on that scaffold
    $n++;
  } else {
    next;
  }
}
close $PROT_fh;
print STDERR "\n";

## print orphans to warnings file:
if (scalar(keys %orphans) > 0) {
  print STDERR "[WARN] Orphans found! ".scalar(keys %orphans)." gene names not in GFF\n";
  open (my $WARN, ">$warnings_file") or die "[ERROR] Cannot open file $warnings_file: $!\n";
  print $WARN "Following gene names not found in $gff_file (".scalar(keys %orphans)." genes):\n";
  foreach (nsort keys %orphans) {
    print $WARN "$_\n";
  }
  close $WARN;
}

## parse HGT_results file:
print STDERR "[INFO] Parsing '$in_file' HGT results file...";
my ($RESULTS_fh, $HGTc_total);
open ($RESULTS_fh, $in_file) or die "[ERROR] Cannot open '$in_file': $!\n";
while (<$RESULTS_fh>) {
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
    $HGTc_total++;
  } else {
    $HGT_evidence = 1; ##intermediate
  }
  #
  # ## apply regex to protein ID in $F[0] if specifed:
  # my $protein_id; ## we want this to map to GFF IDs
  # if ($regexp_option) {
  #   ($protein_id = $F[0]) =~ $regexp; ## apply regex
  # } elsif ($mapping == 1) {
  #   $protein_id = $names_map{$F[0]}; ## inherit GFF ID
  # } else {
  #   $protein_id = $F[0]; ## do nowt
  # }

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
close $RESULTS_fh;
print STDERR " found ".commify(scalar(keys %hgt_results))." entries\n";
print STDERR "[INFO] Of these, ".colored(commify($HGTc_total), 'green')." (".colored(percentage($HGTc_total,$proteins_total), "\%", 'green')." total proteins) passed initial thresholds for HGTc\n";
print STDERR "[INFO] Evaluating results...\n";

## iterate through pseudo-GFF:
open (my $LOC, ">$locations_file") or die "[ERROR] Cannot open file $locations_file: $!\n";
open (my $SUM, ">$summary_file") or die "[ERROR] Cannot open file $summary_file: $!\n";
open (my $GENE, ">$genes_file") or die "[ERROR] Cannot open file $genes_file: $!\n";
open (my $HEV, ">$heavy_file") or die "[ERROR] Cannot open file $heavy_file: $!\n";
open (my $BED, ">$bed_file") or die "[ERROR] Cannot open file $bed_file: $!\n" if ($bed);
#print $LOC "## HGT_locations\n##".`date`."\n";
print $LOC join ("\t", "#SCAFFOLD","START","END","GENE","SCORE","STRAND","INTRONS","hU","EVIDENCE","TAXONOMY","\n");
print $SUM join ("\t", "#SCAFFOLD","NUMGENES","UNASSIGNED","GOOD_INGRP","INTERMEDIATE","GOOD_OUTGRP","PROPORTION_OUTGRP","IS_LINKED","\n");
print $HEV join ("\t", "#SCAFFOLD","NUMGENES","UNASSIGNED","GOOD_INGRP","INTERMEDIATE","GOOD_OUTGRP","PROPORTION_OUTGRP","IS_LINKED","\n");
my ($good_outgrp_total,$good_ingrp_total,$intermediate_total,$na_total,$intronized,$linked_total,$is_heavy,$num_HGTc_on_heavy,$num_genes_on_heavy_total) = (0,0,0,0,0,0,0,0,0);

## iterate across scaffolds:
foreach my $chrom (nsort keys %scaffolds) {

  my ($good_outgrp,$good_ingrp,$intermediate,$na,$is_linked) = (0,0,0,0,0);
  my (%good_outgrp_hash);
  ## sort by start coord within the %locations hash:
  foreach my $gene ( sort {$locations{$a}{start}<=>$locations{$b}{start}} @{$scaffolds{$chrom}} ) {

    ## evaluate if gene has associated hU:
    if ( exists($hgt_results{$gene}{hU}) ) {
      if ($hgt_results{$gene}{evidence} == 2) {
        $good_outgrp_hash{$gene} = (); ##add all HGTc to this hash
        $good_outgrp++;
        $intronized++ if $locations{$gene}{introns} > 0;
        print $LOC join ("\t", $chrom,$locations{$gene}{start},$locations{$gene}{end},$gene,".",$locations{$gene}{strand},$locations{$gene}{introns},$hgt_results{$gene}{hU},$hgt_results{$gene}{evidence},$hgt_results{$gene}{taxonomy},"\n");
        print $BED join ("\t", $chrom,$locations{$gene}{start},$locations{$gene}{end},$gene,".",$locations{$gene}{strand},"\n") if ($bed);
      } else {
        $good_ingrp++ if $hgt_results{$gene}{evidence} == 0;
        $intermediate++ if $hgt_results{$gene}{evidence} == 1;
        print $LOC join ("\t", $chrom,$locations{$gene}{start},$locations{$gene}{end},$gene,".",$locations{$gene}{strand},$locations{$gene}{introns},$hgt_results{$gene}{hU},$hgt_results{$gene}{evidence},"\n");
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
    $linked_total += $good_outgrp; ##sum total linked HGT candidates
    $is_linked = 1;
    print $GENE join ("\n", nsort @protein_hash_map{ keys %good_outgrp_hash }) . "\n"; ##print as a list all linked HGTc genes
    print $SUM join ("\t", $chrom,scalar(@{$scaffolds{$chrom}}),$na,$good_ingrp,$intermediate,$good_outgrp,(percentage($good_outgrp,scalar(@{$scaffolds{$chrom}}))),"1","\n");
  } else {
    print $SUM join ("\t", $chrom,scalar(@{$scaffolds{$chrom}}),$na,$good_ingrp,$intermediate,$good_outgrp,(percentage($good_outgrp,scalar(@{$scaffolds{$chrom}}))),"0","\n");
  }

  ## evaluate proportion of HGT candidates per scaffold; print to 'heavy' if > threshold:
  if ( (percentage($good_outgrp,scalar(@{$scaffolds{$chrom}}))) > $heavy ) {
    print $HEV join ("\t", $chrom,scalar(@{$scaffolds{$chrom}}),$na,$good_ingrp,$intermediate,$good_outgrp,(percentage($good_outgrp,scalar(@{$scaffolds{$chrom}}))),$is_linked,"\n");
    $num_genes_on_heavy_total += scalar(@{$scaffolds{$chrom}});
    $num_HGTc_on_heavy += scalar( grep { $hgt_results{$_}{evidence}==2 } @{$scaffolds{$chrom}} );
    $is_heavy++;
  }
$n++;
}
close $LOC;
close $SUM;
close $HEV;
close $BED if ($bed);

## gzip input files if that's how they came
`gzip $gff_file` if ($gff_is_gz == 1);
`gzip $proteins_file` if ($names_is_gz == 1);

print STDERR "\n";
print STDERR "[RESULT] Failed scaffolds: ".colored(commify($is_heavy), 'green')."\n";
print STDERR "[RESULT] Genes on failed scaffolds: ".colored(commify($num_genes_on_heavy_total), 'green')." (total); ".colored(commify($num_HGTc_on_heavy), 'green')." (HGT candidates)\n";
print STDERR "[RESULT] HGTc with introns: ".colored(commify($intronized), 'green bold')."\n";
print STDERR "[RESULT] HGTc linked to ingroup gene: ".colored(commify($linked_total), 'green bold underscore')."\n";
print STDERR "\n[INFO] Finished on ".`date`."\n";

open (my $OVER, ">$oversummary_file") or die "[ERROR] Cannot open file $oversummary_file: $!\n";
print $OVER "[INFO] Infile: $in_file\n";
print $OVER "[INFO] GFF file: $gff_file\n";
print $OVER "[INFO] Protein ID file: $proteins_file\n";
print $OVER "[INFO] Protein ID file: '$proteins_file'\n";
print $OVER "[INFO] hU threshold to determine strong evidence for OUTGROUP: >= ".colored($outgrp_threshold, 'yellow')."\n";
print $OVER "[INFO] hU threshold to determine strong evidence for INGROUP: <= ".colored($ingrp_threshold, 'yellow')."\n";
print $OVER "[INFO] CHS threshold to determine strong evidence for OUTGROUP: >= ".colored("$CHS_threshold\%", 'yellow')."\n";
print $OVER "[INFO] Proportion of genes >= hU threshold to find contaminant scaffolds: ".colored("$heavy\%", 'yellow')."\n";
print $OVER "[INFO] Applying regex '$regexp_option' to names in '$proteins_file'\n" if ( $regexp_option );
print $OVER "[INFO] Write bedfile: TRUE\n" if ( $bed );
print $OVER "[RESULT] Failed scaffolds: ".colored(commify($is_heavy), 'green')."\n";
print $OVER "[RESULT] Genes on failed scaffolds: ".colored(commify($num_genes_on_heavy_total), 'green')." (total); ".colored(commify($num_HGTc_on_heavy), 'green')." (HGT candidates)\n";
print $OVER "[RESULT] HGTc with introns: ".colored(commify($intronized), 'green bold')."\n";
print $OVER "[RESULT] HGTc linked to ingroup gene: ".colored(commify($linked_total), 'green bold underscore')."\n";
print $OVER "[INFO] Finished on ".`date`."\n";
close $OVER;


################################################################################ SUBS

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

sub decompress {
  die "[ERROR] Cannot find gunzip!\n" if ( system ("gunzip --help &>/dev/null") != 0 );
  my $in_file = $_[0];
  my $out_file = $in_file =~ s/\.gz//r;
  `gunzip $in_file`;
  return $out_file;
}

__END__
