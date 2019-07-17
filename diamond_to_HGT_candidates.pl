#!/usr/bin/env perl

## author: reubwn Nov 2016
## added hU functionality May 2017

use strict;
use warnings;

use Getopt::Long;
use Term::ANSIColor;
use Sort::Naturally;
use Data::Dumper qw(Dumper);
use List::Util qw(reduce sum min max);

my $usage = "
SYNOPSIS
  Goal is to take a taxified Diamond or BLAST file, and for each hit recurse up the tax tree until
  that hit can be categorised into **INGROUP** versus **OUTGROUP** (e.g., Metazoa vs non-Metazoa).

  1. Calculate HGT Index (hU): Calculate hU based on best bitscores to INGROUP vs OUTGROUP (default >= 30)
  2. Get Query Category: For each query, calculate the bitscoresum for ingroup vs outgroup across **all hits**;
     the category with the highest bitscoresum is the \"winner\"
  3. Get Consensus Hit Support: Assess support for the winning query taxid from secondary hits; winning taxid
     is well-supported if the rest of the hits agree with the INGROUP/OUTGROUP categorisation above --support_threshold
     (default = 90%)
  4. Print: General results printed to *.HGT_results, candidate HGT genes printed to *.HGT_candidates (default
     is to print candidates passing both hU and CHS thresholds)

DEFINITIONS
  HGT Index (hU) = (Best-hit bitscore for non-Metazoa) - (Best-hit bitscore for Metazoa)
  Alien Index (AI) = log10((Best-hit Evalue for Metazoa) + 1e-200) - log10((Best-hit Evalue for non-Metazoa) + 1e-200)
  Consesus Hit Support (CHS) = Proportion of all hits that agree with hU classification

INPUT
  Taxified diamond/blast text file, default is to have TaxID in the 13th column.
  OR a list of files to analyse, one per line. This is quicker for multiple files as taxonomy databases will only need build once.
  If list is provided, TaxID to skip (-k) can be specified as a second column in the list file and the path to the input proteins (-f)
  in the third column.

OUTPUTS
  A \"\*.HGT_results.txt\" file with the support values for each query; a \"\*.HGT_candidates.txt\" file with queries
  showing support over the specified thresholds.

NOTES
  Please supply input fasta file to get accurate proportion HGTc, as queries with no hits to UniRef90 are
  excluded from the diamond hits file

OPTIONS [* required]
  -i|--in              [FILE]*   : taxified diamond/blast results file (accepts gzipped)
  -l|--list            [FILE]    : list of diamond/blast files to analyse
  -p|--path            [PATH]*   : path to dir of taxonomy files
  -o|--nodes           [FILE]    : path to nodes.dmp
  -a|--names           [FILE]    : path to names.dmp
  -m|--merged          [FILE]    : path to merged.dmp
  -n|--nodesDB         [FILE]    : nodesDB.txt file from blobtools
  -f|--faa             [FILE]*   : fasta file of protein sequences under investigation (accepts gzipped)
  -t|--taxid_ingroup   [INT]     : NCBI taxid to define 'ingroup' [default=33208 (Metazoa)]
  -k|--taxid_skip      [INT]     : NCBI taxid to skip; hits to this taxid will not be considered
  -s|--CHS_threshold   [FLOAT]   : Consesus Hits Support threshold [default>=90\%]
  -u|--hU_threshold    [INT]     : hU threshold [default>=30]
  -e|--evalue_column   [INT]     : define evalue column [default=11]
  -b|--bitscore_column [INT]     : define bitscore column [default=12]
  -c|--taxid_column    [INT]     : define taxid column [default=13]
  -x|--prefix          [FILE]    : filename prefix for outfile [default=INFILE]
  -v|--verbose                   : say more things
  -h|--help                      : this help message
\n";

my ($in_file,$list,$nodes_file,$path,$names_file,$merged_file,$nodesDB_file,$proteins_file,$prefix,$outfile,$hgtcandidates_file,$warnings_file,$verbose,$help,$debug);
my $taxid_threshold = 33208; ##metazoa
my $taxid_skip_cmd = 0; ##default is 0, not a valid NCBI taxid and should not affect the tree recursion; NB Rotifera = 10190
my $support_threshold = 90;
my $hU_threshold = 30;
my $scoring = "sum";
my $evalue_column = 11;
my $bitscore_column = 12;
my $taxid_column = 13;

GetOptions (
  'i|infiles:s'           => \$in_file,
  'l|list:s'              => \$list,
  'p|path:s'              => \$path,
  'o|nodes:s'             => \$nodes_file,
  'a|names:s'             => \$names_file,
  'm|merged:s'            => \$merged_file,
  'n|nodesDB:s'           => \$nodesDB_file,
  'f|faa:s'               => \$proteins_file,
  't|taxid_threshold:i'   => \$taxid_threshold,
  'k|taxid_skip:i'        => \$taxid_skip_cmd,
  's|support_threshold:f' => \$support_threshold,
  'u|hU_threshold:i'      => \$hU_threshold,
  'e|evalue_column:i'     => \$evalue_column,
  'c|taxid_column:i'      => \$taxid_column,
  'b|bitscore_column:i'   => \$bitscore_column,
  'x|prefix:s'            => \$prefix,
  'v|verbose'             => \$verbose,
  'h|help'                => \$help,
  'd|debug'               => \$debug
);

die $usage if $help;
die "$usage\n[ERROR] Missing -i or -l argument!\n\n" unless ( $in_file || $list );
die "$usage\n[ERROR] Missing -p or -n or -oam argument!\n\n" unless ( $path || $nodesDB_file || ($nodes_file && $names_file && $merged_file) );
die "$usage\n[ERROR] Missing -f argument!\n\n" unless ( $proteins_file );

############################################## PARSE NODES

## parse nodes and names:
my (%nodes_hash, %names_hash, %rank_hash);
if ($path) {
  print STDERR "[INFO] Building taxonomy databases from tax files in '$path'...\n";
  open(my $NODES, "$path/nodes.dmp") or die $!;
  while (<$NODES>) {
    chomp;
    next if /\#/;
    my @F = map { s/^\s+|\s+$//gr } split (/\|/, $_); ## split nodes.dmp file on \s+|\s+ regex
    $nodes_hash{$F[0]} = $F[1]; ## key= child taxid; value= parent taxid
    $rank_hash{$F[0]} = $F[2]; ## key= taxid; value= rank
  }
  close $NODES;
  open (my $NAMES, "$path/names.dmp") or die $!;
  while (<$NAMES>) {
    chomp;
    next if /\#/;
    my @F = map { s/^\s+|\s+$//gr } split (/\|/, $_);
    $names_hash{$F[0]} = $F[1] if ($F[3] eq "scientific name"); ## key= taxid; value= species name
  }
  close $NAMES;
  if (-e "$path/merged.dmp") {
    open (my $MERGED, "$path/merged.dmp") or die $!;
    while (<$MERGED>) {
      chomp;
      next if /\#/;
      my @F = map { s/^\s+|\s+$//gr } split (/\|/, $_);
      $nodes_hash{$F[0]} = $F[1]; ## key= old taxid; value= new taxid
      ## this will behave as if old taxid is a child of the new one, which is OK I guess
    }
  }
} elsif ($nodes_file && $names_file) {
  print STDERR "[INFO] Building taxonomy databases from '$nodes_file' and '$names_file'...";
  open(my $NODES, $nodes_file) or die $!;
  while (<$NODES>) {
    chomp;
    next if /\#/;
    my @F = map { s/^\s+|\s+$//gr } split (/\|/, $_); ## split nodes.dmp file on \s+|\s+ regex
    $nodes_hash{$F[0]} = $F[1]; ## key= child taxid; value= parent taxid
    $rank_hash{$F[0]} = $F[2]; ## key= taxid; value= rank
  }
  close $NODES;
  open (my $NAMES, $names_file) or die $!;
  while (<$NAMES>) {
    chomp;
    next if /\#/;
    my @F = map { s/^\s+|\s+$//gr } split (/\|/, $_);
    $names_hash{$F[0]} = $F[1] if ($F[3] eq "scientific name"); ## key= taxid; value= species name
  }
  close $NAMES;
  if ($merged_file) {
    open (my $MERGED, $merged_file) or die $!;
    while (<$MERGED>) {
      chomp;
      next if /\#/;
      my @F = map { s/^\s+|\s+$//gr } split (/\|/, $_);
      $nodes_hash{$F[0]} = $F[1]; ## key= old taxid; value= new taxid
      ## this will behave as if old taxid is a child of the new one, which is OK I guess
    }
  }
} elsif ($nodesDB_file) {
  print STDERR "[INFO] Building taxonomy databases from '$nodesDB_file'...";
  open(my $NODES, $nodesDB_file) or die $!;
  while (<$NODES>) {
    chomp;
    next if /\#/;
    my @F = split (/\t/, $_);
    $nodes_hash{$F[0]} = $F[3]; ## key= child taxid; value= parent taxid
    $names_hash{$F[0]} = $F[2]; ## key= taxid; value= species name
    $rank_hash{$F[0]} = $F[1]; ## key= taxid; value= rank
  }
  close $NODES;
}
## print some info to STDERR:
print STDERR "[INFO] Nodes parsed: ".commify(scalar(keys %nodes_hash))."\n";
print STDERR "[INFO] Threshold taxid set to '$taxid_threshold' ($names_hash{$taxid_threshold})\n";
print STDERR "[INFO] INGROUP set to '$names_hash{$taxid_threshold}'; OUTGROUP is therefore 'non-$names_hash{$taxid_threshold}'\n";

############################################ INFILES

my @in_files;
my %taxid_skip_hash;
my %prots_file_hash;
if ($list) {
  open (my $LIST, $list) or die $!;
  while (<$LIST>) {
    chomp;
    my @F = split (/\s+/, $_);
    if ( -f $F[0] ) {
      push (@in_files, $F[0]);
    } else {
      die "[ERROR] File $F[0] does not exist!\n";
    }
    ## fetch the TaxID to skip from file:
    if ($F[1]) {
      $taxid_skip_hash{$F[0]} = $F[1];
    }
    ## count number of prots in each original input file
    if ( -f $F[2] ) {
      my $num_prots;
      if ($F[2] =~ m/gz$/) { ## can grep from gzipped
        chomp ($num_prots = `zgrep -c ">"`);
      } else {
        chomp ($num_prots = `grep -c ">" $F[2]`);
      }
      $prots_file_hash{$F[0]} = {'file' => $F[2], 'num' => $num_prots};
    }
  }
  close $LIST;
  print STDERR "[INFO] Number of files in '$list': ".@in_files."\n";
} else {
  @in_files = split (/\s+/, $in_file); ## split $in_file string
  print STDERR "[INFO] Number of files: ".@in_files."\n";
  ## get number of prots in each original input file
  my @prots_files = split (/\s+/, $proteins_file);
  if (scalar(@prots_files) == scalar(@in_files)) {
    for my $i (0 .. $#prots_files) {
      my $num_prots;
      if ( -f $prots_files[$i] ) {
        if ($prots_files[$i] =~ m/gz$/) { ## can grep from gzipped
          chomp ($num_prots = `zgrep -c ">" $prots_files[$i]`);
        } else {
          chomp ($num_prots = `grep -c ">" $prots_files[$i]`);
        }
        $prots_file_hash{$in_files[$i]} = {'file' => $prots_files[$i], 'num' => $num_prots};
      } else {
        die "[ERROR] File $prots_files[$i] does not exist!\n";
      }
    }
  } else {
    die "[ERROR] Number of hits files (".scalar(@in_files).") not equal to number of proteins files ".scalar(@prots_files)."!\n";
  }
}

############################################ OUTFILES

foreach my $in (@in_files) { ## iterate over multiple files if required
  print STDERR "[INFO] Analysing DIAMOND hits file '$in'\n";
  print STDERR "[INFO] Based on proteins in '$prots_file_hash{$in}{'file'}' (".commify($prots_file_hash{$in}{'num'})." total input queries)\n";

  ## define TaxID to skip:
  my $taxid_skip;
  if ($taxid_skip_cmd) { ## from command line input
    $taxid_skip = $taxid_skip_cmd;
    print STDERR "[INFO] Skipping any hits to taxid '$taxid_skip' ($names_hash{$taxid_skip})\n";
  } elsif ($taxid_skip_hash{$in}) {
    $taxid_skip = $taxid_skip_hash{$in};
    print STDERR "[INFO] Skipping any hits to taxid '$taxid_skip' ($names_hash{$taxid_skip})\n";
  } else {
    print STDERR "[WARN] TaxID to skip is not set!\n";
  }

  ## define outfiles:
  if ($prefix) {
    $outfile = "$prefix.HGT_results.$names_hash{$taxid_threshold}.txt";
    $hgtcandidates_file = "$prefix.HGT_candidates.$names_hash{$taxid_threshold}.hU$hU_threshold.CHS$support_threshold.txt";
    $warnings_file = "$prefix.HGT_warnings.txt";
  } else {
    $outfile = "$in.HGT_results.$names_hash{$taxid_threshold}.txt";
    $hgtcandidates_file = "$in.HGT_candidates.$names_hash{$taxid_threshold}.hU$hU_threshold.CHS$support_threshold.txt";
    $warnings_file = "$in.HGT_warnings.txt";
  }

  ## open outfiles:
  open (my $OUT, ">",$outfile) or die $!;
  open (my $HGT, ">",$hgtcandidates_file) or die $!;
  open (my $WARN, ">",$warnings_file) or die $!;
  my @header = (
    "#QUERY",
    "INGROUP",
    "OUTGROUP",
    "hU",
    "BIT_OUT",
    "BIT_IN",
    "AI",
    "EVAL_OUT",
    "EVAL_IN",
    "BBSUMCAT",
    "CHS",
    "TAXONOMY"
  );
  print $OUT join ("\t", @header, "\n");
  print $HGT join ("\t", @header, "\n");

  ############################################## PARSE DIAMOND

  ## parse Diamond file:
  my (%bitscores_per_query_hash, %evalues_per_query_hash);
  my ($total_entries,$skipped_entries_because_bad_taxid,$skipped_entries_because_skipped_taxid,$skipped_entries_because_unassigned) = (0,0,0,0);

  my $DIAMOND;
  if ($in =~ m/gz$/) {
    open ($DIAMOND, "zcat $in |") or die $!; ## if gzipped
    print STDERR "[INFO] File is gzipped\n";
  } else {
    open ($DIAMOND, $in) or die $!;
  }

  while (<$DIAMOND>) {
    chomp;
    next if /^\#/;
    $total_entries++;
    my @F = split (m/\s+/, $_);
    if (scalar(@F) <= 1) {
      die "[ERROR] File '$in' did not split properly: is it tab or space delimited?\n";
    } elsif ( (scalar(@F)<$taxid_column) or (scalar(@F)<$bitscore_column) or (scalar(@F)<$evalue_column) ) {
      die "[ERROR] File '$in' line $. did not split properly: number of columns is too small (".@F.")\n";
    }
    if ($F[($taxid_column-1)] !~ m/\d+/) {
      print $WARN join ("\t", $F[0], $., $F[($taxid_column-1)], "invalid/unrecognised/absent taxid", "\n");
      $skipped_entries_because_bad_taxid++;
      next;
    } elsif (check_taxid_has_parent($F[($taxid_column-1)]) == 1) {
      print $WARN join ("\t", $F[0], $., $F[($taxid_column-1)], "invalid/unrecognised parent taxid", "\n");
      $skipped_entries_because_bad_taxid++;
      next;
    } elsif ( tax_walk($F[($taxid_column-1)], $taxid_skip) eq "ingroup" ) { ## do not include any hits to within taxid $taxid_skip
      print $WARN join ("\t", $F[0], $., $F[($taxid_column-1)], "taxid within skipped ($taxid_skip)", "\n") if $verbose;
      $skipped_entries_because_skipped_taxid++;
      next;
    } elsif ( tax_walk($F[($taxid_column-1)]) eq "unassigned" ) {
      print $WARN join ("\t", $F[0], $., $F[($taxid_column-1)], "taxid unassigned/unclassified", "\n");
      $skipped_entries_because_unassigned++;
      next;
    } else {
      ## push all bitscores and evalues for every taxid into an array within a hash within a hash:
      push @{ $bitscores_per_query_hash{$F[0]}{$F[($taxid_column-1)]} }, $F[($bitscore_column-1)]; ## key= query; value= hash{ key= taxid; value= array[ bitscores ]}
      push @{ $evalues_per_query_hash{$F[0]}{$F[($taxid_column-1)]} }, $F[$evalue_column-1]; ## key= query; value= hash{ key= taxid; value= array [ evalues ]}
    }
  }
  close $DIAMOND;
  print STDERR "[INFO] Total number of hits parsed: ".commify($total_entries)."\n";
  print STDERR "[WARN] There were ".commify($skipped_entries_because_bad_taxid)." (".percentage($skipped_entries_because_bad_taxid,$total_entries)."\%) invalid taxid entries\n" if $skipped_entries_because_bad_taxid > 0;
  print STDERR "[WARN] There were ".commify($skipped_entries_because_skipped_taxid)." (".percentage($skipped_entries_because_skipped_taxid,$total_entries)."\%) skipped taxid entries\n" if $skipped_entries_because_skipped_taxid > 0;
  print STDERR "[WARN] There were ".commify($skipped_entries_because_unassigned)." (".percentage($skipped_entries_because_unassigned,$total_entries)."\%) unassigned/unclassified taxid entries\n" if $skipped_entries_because_unassigned > 0;

  ############################################ DEBUG

  print Dumper \%bitscores_per_query_hash if $debug;
  print Dumper \%evalues_per_query_hash if $debug;

  ############################################ MAIN

  ## get winning bitscore and taxid; calculate congruence among all taxids for all hits per query:
  my ($processed,$ingroup,$ingroup_supported,$outgroup,$outgroup_supported,$hU_supported,$AI_supported,$unassigned) = (0,0,0,0,0,0,0,0);
  my %hgt_candidates;
  print STDERR "[INFO] Calculating hU, AI and consensus hit support...\n";
  print STDERR "\n" if $verbose;
  foreach my $query (nsort keys %bitscores_per_query_hash) {
    my %bitscore_hash = %{ $bitscores_per_query_hash{$query} }; ## key= taxid; value= \@array of all bitscores for that taxid
    my %evalue_hash = %{ $evalues_per_query_hash{$query} }; ## key= taxid; value= \@array of all evalues for that taxid

    ## calculate alien index (AI):
    my ($ingroup_best_evalue, $outgroup_best_evalue) = (1,1);
    foreach my $taxid (keys %evalue_hash) {
      my $min_evalue = min( @{ $evalue_hash{$taxid} } );
      if (tax_walk($taxid) eq "ingroup") {
        $ingroup_best_evalue = $min_evalue if ($min_evalue < $ingroup_best_evalue); ## only accept it if it's LOWER (better) than current Evalue
      } elsif (tax_walk($taxid) eq "outgroup") {
        $outgroup_best_evalue = $min_evalue if ($min_evalue < $outgroup_best_evalue);
      }
    }
    ## AI = log(AI_in+constant) - log(AI_out+constant):
    my $AI = ( log10($ingroup_best_evalue + 1e-200) - log10($outgroup_best_evalue + 1e-200) );
    $AI_supported++ if $AI >= $hU_threshold;

    ## calculate HGT index (hU):
    my ($ingroup_best_bitscore, $outgroup_best_bitscore) = (0,0);
    foreach my $taxid (keys %bitscore_hash) {
      my $max_bitscore = max( @{ $bitscore_hash{$taxid} } );
      if (tax_walk($taxid) eq "ingroup") {
        $ingroup_best_bitscore = $max_bitscore if ($max_bitscore > $ingroup_best_bitscore); ## only accept it if it's HIGHER (better) than current bitscore
      } elsif (tax_walk($taxid) eq "outgroup") {
        $outgroup_best_bitscore = $max_bitscore if ($max_bitscore > $outgroup_best_bitscore);
      }
    }
    ## hU = B_out - B_in:
    my $hU = ($outgroup_best_bitscore - $ingroup_best_bitscore);
    $hU_supported++ if $hU >= $hU_threshold;

    ## calculate bitscoresums per taxid; get taxid of highest bitscoresum; get support for winning taxid from other hits:
    my (%bitscoresum_hash, %count_categories, %support_categories);
    my ($ingroup_bitscoresum, $outgroup_bitscoresum) = (0,0);

    foreach my $taxid (keys %bitscore_hash) {
      if ($scoring eq "sum") {
        if (tax_walk($taxid) eq "ingroup") {
          $ingroup_bitscoresum += sum( @{ $bitscore_hash{$taxid} } );
          #print STDERR join "\t", "\t", $query, $taxid, tax_walk($taxid), sum( @{ $bitscore_hash{$taxid} } ),"\n"; ## uncomment to see info for each each hit
        } elsif (tax_walk($taxid) eq "outgroup") {
          $outgroup_bitscoresum += sum( @{ $bitscore_hash{$taxid} } );
          #print STDERR join "\t", "\t", $query, $taxid, tax_walk($taxid), sum( @{ $bitscore_hash{$taxid} } ),"\n"; ## uncomment to see info for each each hit
        }
        my $bitscoresum = sum( @{ $bitscore_hash{$taxid} } );
        $bitscoresum_hash{$taxid} = $bitscoresum; ## key= taxid; value= bitscoresum
      } elsif ($scoring eq "individual") {
        my $bitscoresum = sum( @{ $bitscore_hash{$taxid} } );
        $bitscoresum_hash{$taxid} = $bitscoresum; ## key= taxid; value= bitscoresum
      }
      $count_categories{tax_walk($taxid)}++; ## count categories; if each hit's taxid falls within/outwith the $taxid_threshold
    }
    print "$query:\n" if $debug; ## debug
    print Dumper \%bitscoresum_hash if $debug; ## debug

    foreach my $cat (keys %count_categories) {
      $support_categories{$cat} = percentage($count_categories{$cat}, scalar(keys %bitscore_hash)); ## calculate proportion of support for the category of the winner
    }

    ## get taxid with highest bitscore:
    my ($taxid_with_highest_bitscore,$taxid_with_highest_bitscore_category,$taxid_with_highest_bitscore_category_support);
    if ($scoring eq "sum") {
      $taxid_with_highest_bitscore_category = $ingroup_bitscoresum > $outgroup_bitscoresum ? "ingroup" : "outgroup"; ## define query category based on bitscoresums of ingroup vs outgroup
      $taxid_with_highest_bitscore_category_support = $support_categories{$taxid_with_highest_bitscore_category}; ## % support from other hits
      $taxid_with_highest_bitscore = List::Util::reduce { $bitscoresum_hash{$b} > $bitscoresum_hash{$a} ? $b : $a } keys %bitscoresum_hash; ## winning taxid
      print STDERR "[INFO] [$query] Bitscoresum for INGROUP ($names_hash{$taxid_threshold}): $ingroup_bitscoresum\n" if $verbose;
      print STDERR "[INFO] [$query] Bitscoresum for OUTGROUP (non-$names_hash{$taxid_threshold}): $outgroup_bitscoresum\n" if $verbose;
      ## PRINT TO OUT:
      print $OUT join (
        "\t",
        $query, ##QUERY
        $names_hash{$taxid_threshold}, ##INGROUP
        "non-$names_hash{$taxid_threshold}", ##OUTGROUP
        $hU, ##hU
        $outgroup_best_bitscore, ##BIT_OUT
        $ingroup_best_bitscore, ##BIT_IN
        $AI, ##AI
        $outgroup_best_evalue, ##EVAL_OUT
        $ingroup_best_evalue, ##EVAL_IN
        ($ingroup_bitscoresum > $outgroup_bitscoresum ? "INGROUP" : "OUTGROUP"), ##BBSUMCAT
        $taxid_with_highest_bitscore_category_support, ##CHS
        tax_walk_to_get_rank_to_phylum($taxid_with_highest_bitscore), ##TAXONOMY
         "\n"
      );

    } elsif ($scoring eq "individual") {
      die "[ERROR] Sorry, individual scoring not supported atm, please change to 'sum'\n";
    }

    print STDERR "[INFO] [$query] Decision of bestsum bitscore: '$taxid_with_highest_bitscore_category' (support = $taxid_with_highest_bitscore_category_support)\n" if $verbose;
    print STDERR "[INFO] [$query] Best evalue for INGROUP ($names_hash{$taxid_threshold}): $ingroup_best_evalue\n" if $verbose;
    print STDERR "[INFO] [$query] Best evalue for OUTGROUP (non-$names_hash{$taxid_threshold}): $outgroup_best_evalue\n" if $verbose;
    print STDERR "[INFO] [$query] Alien Index = $AI\n[----]\n" if $verbose;

    ## count genes in various categories:
    if ( $taxid_with_highest_bitscore_category eq "unassigned" ) {
      $unassigned++;
    } elsif ( $taxid_with_highest_bitscore_category eq "ingroup" ) {
      $ingroup++;
      $ingroup_supported++ if ( $taxid_with_highest_bitscore_category_support >= $support_threshold );
    } elsif ( $taxid_with_highest_bitscore_category eq "outgroup" ) {
      $outgroup++;
      $outgroup_supported++ if ( $taxid_with_highest_bitscore_category_support >= $support_threshold );
    }

    ## PRINT TO HGT_CANDIDATES:
    ## print all queries with hU>=threshold AND CHS>=threshold to HGT_candidates file:
    if ( ($hU >= $hU_threshold) && ($taxid_with_highest_bitscore_category eq "outgroup") && ($taxid_with_highest_bitscore_category_support >= $support_threshold) ) {
      print $HGT join (
        "\t",
        $query, ##QUERY
        $names_hash{$taxid_threshold}, ##INGROUP
        "non-$names_hash{$taxid_threshold}", ##OUTGROUP
        $hU, ##hU
        $outgroup_best_bitscore, ##BIT_OUT
        $ingroup_best_bitscore, ##BIT_IN
        $AI, ##AI
        $outgroup_best_evalue, ##EVAL_OUT
        $ingroup_best_evalue, ##EVAL_IN
        ($ingroup_bitscoresum > $outgroup_bitscoresum ? "INGROUP" : "OUTGROUP"), ##BBSUMCAT
        $taxid_with_highest_bitscore_category_support, ##CHS
        tax_walk_to_get_rank_to_phylum($taxid_with_highest_bitscore), ##TAXONOMY
        "\n"
      );
      $hgt_candidates{$query} = ();
    }

    ## progress:
    $processed++;
    if ($processed % 1000 == 0){
      print STDERR "\r[INFO] Processed ".commify($processed)." queries...";
      $| = 1;
    }
  }
  close $OUT;
  close $HGT;
  close $WARN;

  print STDERR "\r[INFO] Processed ".commify($processed)." queries\n";
  print STDERR "[INFO] All results are printed to '$outfile'\n";
  print STDERR "[INFO] HGT candidates are printed to '$hgtcandidates_file'\n";
  #print STDERR "[INFO] Number of queries in unassigned/unclassified category: ".commify($unassigned)."\n" if $unassigned > 0;
  #print STDERR "[INFO] Number of queries in OUTGROUP category ('non-$names_hash{$taxid_threshold}'): ".commify($outgroup)."\n";
  #print STDERR "[INFO] Number of queries in OUTGROUP category ('non-$names_hash{$taxid_threshold}') with CHS >= $support_threshold\%: ".commify($outgroup_supported)."\n";
  print STDERR "[INFO] Number of queries with HGT Index >= $hU_threshold: ".colored(commify($hU_supported), 'green bold')."\n";
  print STDERR "[INFO] Number of queries with HGT Index >= $hU_threshold and CHS >= $support_threshold\% to non-$names_hash{$taxid_threshold}: ".colored(commify(scalar(keys(%hgt_candidates))), 'green bold underscore')." (".percentage(scalar(keys(%hgt_candidates)),$processed)."\% of ".commify($processed)." processed";
  if ($prots_file_hash{$in}) {
    print STDERR " or ".colored(percentage(scalar(keys(%hgt_candidates)),$prots_file_hash{$in}{'num'})."\%", 'green bold underscore')." of ".commify($prots_file_hash{$in}{'num'})." total input queries)\n";
  } else {
    print STDERR ")\n";
  }
  #print STDERR "[INFO] Number of queries with Alien Index (AI) >= $hU_threshold: ".commify($AI_supported)."\n";
  #print STDERR "[INFO] NUMBER OF HGT CANDIDATES: ".commify(scalar(keys(%hgt_candidates)))."\n";
  print STDERR "[INFO] Finished on ".`date`."\n";

}##in_files loop

############################################ SUBS

sub check_taxid_has_parent {
  my $taxid = $_[0];
  my $result = 0;
  unless ($nodes_hash{$taxid}) {
    $result = 1;
  }
  return $result; ## 0 = taxid exists; 1 = taxid does not exist
}

sub tax_walk {
    my $taxid = $_[0];
    my $walk_to;
    if (exists $_[1]) {
      $walk_to = $_[1];
    } else {
      $walk_to = $taxid_threshold; ## default is metazoa
    }

    ## first parent:
    my $parent = $nodes_hash{$taxid};
    my $result;

    ## return "unassigned" if hit has no valid taxid
    if ($parent !~ m/\d+/) {
      $result = "unassigned";
      return $result;
    }

    ## recurse the tree:
    while (1) {
      if ($parent == $walk_to) {
        $result = "ingroup"; ## is ingroup
        last;
      } elsif ($parent == 1) {
        $result = "outgroup"; ## root; i.e., the whole tree has been recursed without finding $threshold, therefore $taxid must reside in another part of the tax tree
        last;
      } elsif ($parent == 32644) {
        $result = "unassigned"; ## taxid for "unidentified"
        last;
      } elsif ($parent == 12908) {
        $result = "unassigned"; ## taxid for "unclassified sequences"
        last;
      } else { ## walk up the tree!
        $parent = $nodes_hash{$parent};
      }
    }
    return $result;
}

sub tax_walk_to_get_rank_to_phylum {
  my $taxid = $_[0];
  my $parent = $nodes_hash{$taxid};
  my $parent_rank = $rank_hash{$parent};
  my ($phylum,$kingdom,$superkingdom) = ("undef","undef","undef");

  while (1) {
    if ($parent_rank eq "phylum") {
      $phylum = $names_hash{$parent};
      #print "Found phylum: $phylum\n";
      $parent = $nodes_hash{$parent};
      $parent_rank = $rank_hash{$parent};
      next;
    } elsif ($parent_rank eq "kingdom") {
      $kingdom = $names_hash{$parent};
      #print "Found phylum: $kingdom\n";
      $parent = $nodes_hash{$parent};
      $parent_rank = $rank_hash{$parent};
      next;
    } elsif ($parent_rank eq "superkingdom") {
      $superkingdom = $names_hash{$parent};
      #print "Found phylum: $superkingdom\n";
      last;
    } elsif ($parent == 1) {
      last;
    } else {
      $parent = $nodes_hash{$parent};
      $parent_rank = $rank_hash{$parent};
    }
  }
  my $result = join (";",$superkingdom,$kingdom,$phylum);
  $result =~ s/\s+/\_/g; ## replace spaces with underscores
  return $result;
}

sub tax_walk_to_get_rank_to_species {
  my $taxid = $_[0];
  my $parent = $nodes_hash{$taxid};
  my $parent_rank = $rank_hash{$parent};
  my ($species,$genus,$family,$order,$class,$phylum,$kingdom,$superkingdom) = ("undef","undef","undef","undef","undef","undef","undef","undef");

  while (1) {
    if ($parent_rank eq "species") {
      $species = $names_hash{$parent};
      $parent = $nodes_hash{$parent};
      $parent_rank = $rank_hash{$parent};
      next;
    } elsif ($parent_rank eq "genus") {
      $genus = $names_hash{$parent};
      $parent = $nodes_hash{$parent};
      $parent_rank = $rank_hash{$parent};
      next;
    } elsif ($parent_rank eq "family") {
      $family = $names_hash{$parent};
      $parent = $nodes_hash{$parent};
      $parent_rank = $rank_hash{$parent};
      next;
    } elsif ($parent_rank eq "order") {
      $order = $names_hash{$parent};
      $parent = $nodes_hash{$parent};
      $parent_rank = $rank_hash{$parent};
      next;
    } elsif ($parent_rank eq "class") {
      $class = $names_hash{$parent};
      $parent = $nodes_hash{$parent};
      $parent_rank = $rank_hash{$parent};
      next;
    } elsif ($parent_rank eq "phylum") {
      $phylum = $names_hash{$parent};
      $parent = $nodes_hash{$parent};
      $parent_rank = $rank_hash{$parent};
      next;
    } elsif ($parent_rank eq "kingdom") {
      $kingdom = $names_hash{$parent};
      $parent = $nodes_hash{$parent};
      $parent_rank = $rank_hash{$parent};
      next;
    } elsif ($parent_rank eq "superkingdom") {
      $superkingdom = $names_hash{$parent};
      last;
    } elsif ($parent == 1) {
      last;
    } else {
      $parent = $nodes_hash{$parent};
      $parent_rank = $rank_hash{$parent};
    }
  }
  my $result = join (";",$superkingdom,$kingdom,$phylum,$class,$order,$family,$genus,$species);
  $result =~ s/\s+/\_/g; ## replace spaces with underscores
  return $result;
}

sub log10 {
    my $n = shift;
    return log($n)/log(10);
}

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
