#!/usr/bin/env perl

## author: reubwn July 2019

use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use Sort::Naturally;
use Data::Dumper qw(Dumper);

my $usage = "
SYNOPSIS
  Returns info on the taxon breakdown of the sequences in UniRef database.

OPTIONS [*]
  -i|--in   [FILE]* : UniRef fasta file [required]
  -p|--path [STR]*  : path to dir/ containing taxonomy files [required]
  -r|--rank [STR]   : taxonomic rank to report (default: 'phylum')
  -c|--count        : count number of specific rank values given by -r
  -h|--help         : this message
\n";

my ($in_file,$path,$count,$help,$debug);
my $rank_limit = "phylum";

GetOptions (
  'i|in=s'   => \$in_file,
  'p|path=s' => \$path,
  'r|rank:s' => \$rank_limit,
  'c|count'  => \$count,
  'h|help'   => \$help,
  'd|debug'  => \$debug
);

die $usage if $help;
die $usage unless ($in_file && $path);

############################################## PARSE NODES
## parse nodes and names:
my (%nodes_hash, %names_hash, %rank_hash, %rank_hash_reversed);
print STDERR "[INFO] Building taxonomy databases from tax files in '$path'...\n";
if ( -f "$path/nodes.dmp") {
  print STDERR "[INFO]   - $path/nodes.dmp\n";
} else {
  die "[ERROR] Cannot find $path/nodes.dmp\n";
}
open (my $NODES, "$path/nodes.dmp") or die $!;
while (<$NODES>) {
  chomp;
  next if /\#/;
  my @F = map { s/^\s+|\s+$//gr } split (/\|/, $_); ## split nodes.dmp file on \s+|\s+ regex
  $nodes_hash{$F[0]} = $F[1]; ## key= child taxid; value= parent taxid
  $rank_hash{$F[0]} = $F[2]; ## key= taxid; value= rank name
  $rank_hash_reversed{$F[2]} = $F[0]; ## key= rank name; value= taxid
}
close $NODES;
if ( -f "$path/names.dmp" ) {
  print STDERR "[INFO]   - $path/names.dmp\n";
} else {
  die "[ERROR] Cannot find $path/names.dmp\n";
}
open (my $NAMES, "$path/names.dmp") or die $!;
while (<$NAMES>) {
  chomp;
  next if /\#/;
  my @F = map { s/^\s+|\s+$//gr } split (/\|/, $_);
  $names_hash{$F[0]} = $F[1] if ($F[3] eq "scientific name"); ## key= taxid; value= rank name
}
close $NAMES;
if ( -f "$path/merged.dmp" ) {
  print STDERR "[INFO]   - $path/merged.dmp\n";
  open (my $MERGED, "$path/merged.dmp") or die $!;
  while (<$MERGED>) {
    chomp;
    next if /\#/;
    my @F = map { s/^\s+|\s+$//gr } split (/\|/, $_);
    $nodes_hash{$F[0]} = $F[1]; ## key= old taxid; value= new taxid
    ## this will behave as if old taxid is a child of the new one, which is OK I guess
  }
}
print STDERR "[INFO] Nodes parsed: ".scalar(keys %nodes_hash)."\n";

############################################## PARSE UNIREF
print STDERR "[INFO] Parsing NCBI:txids from '$in_file'...\n";
my ($IN, $proteins_total, $rank_count);
if ($in_file =~ m/gz$/) {
  open ($IN, "zcat $in_file |") or die "[ERROR] Cannot open UniRef file: $!\n"; ## grep the fasta headers directly
  # chomp ($proteins_total = `zgrep -c ">" $in_file`); ## get total number of proteins
} else {
  open ($IN, $in_file) or die "[ERROR] Cannot open UniRef file: $!\n";
  # chomp ($proteins_total = `grep -c ">" $in_file`); ## get total number of proteins
}

my %uniref_hash;
while (my $line = <$IN>) {
  if ($line =~ m/^\>/) {
    ## get the NCBI:txid from fasta header
    if ($line =~ m/TaxID=(\d+)\s/) {
      if ( $count ) {
        $rank_count += tax_walk_to_count_rank ($1);
      } else {
        my $rank = tax_walk_to_rank ($1);
        $uniref_hash{$rank}++;
      }
    }
    $proteins_total++;
    ## progress
    if ($proteins_total % 1000 == 0){
      print STDERR "\r[INFO] Processed ".commify($proteins_total)." queries...";
      $| = 1;
    }
  } else {
    next;
  }
}
close $IN;
print STDERR "\n[INFO] Number of sequences in '$in_file': ".commify($proteins_total)."\n";

if ( $count ) {
  print STDERR "[INFO] Rank was '$rank_limit'\n";
  print STDERR "[INFO] Found $rank_count proteins under '$rank_limit' (".percentage($rank_count,$proteins_total)."\%)\n";
} else {

  ############################################## OUTFILE
  my $out_file = basename ($in_file) . $rank_limit . ".txt";
  open (my $OUT, ">$out_file") or die "[ERROR] Cannot open outfile '$out_file': $!\n";

}

############################################ SUBS

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

## this sub returns the <RANK> of given taxid, up to $rank_limit
## built into a hash, can be used to compute distribution of ranks within database
sub tax_walk_to_rank {
  my $taxid = $_[0];
  my $parent = $nodes_hash{$taxid};
  my $parent_rank = $rank_hash{$parent};
  my $result = "undef";

  if ($rank_hash{$taxid} =~ m/^$rank_limit$/i) { ## no need to recurse tree
    $result = $rank_hash{$taxid};
    last;
  } else {
    while (1) { ## else recurse tree
      if ($parent_rank =~ m/^$rank_limit$/i) {
        $result = $rank_hash{$parent};
        last;
      } elsif ($parent == 1) {
        $result = "root";
        last;
      } else { ## climb the tree
        $parent = $nodes_hash{$parent};
        $parent_rank = $rank_hash{$parent};
      }
    }
  }

  $result =~ s/\s+/\_/g; ## replace spaces in some names with underscores
  return $result;
}

## this sub returns 1 if taxid is a child of $rank_limit, else 0
## so can be used to count number of <RANK> in database
sub tax_walk_to_count_rank {
  my $taxid = $_[0];
  my $parent = $nodes_hash{$taxid};
  my $parent_rank = $rank_hash{$parent};
  print STDERR "TaxID=$taxid; Rank=$rank_hash{$taxid}; Parent=$nodes_hash{$taxid}; ParentRank=$rank_hash{$parent}\n" if ( $debug );

  my $result = 0;
  # my $result = "$taxid_rank[$taxid]:$names_hash{$taxid}";

  if ($rank_hash{$taxid} =~ m/^$rank_limit$/i) { ## no need to recurse tree
    $result = 1;
    last;
  } else {
    while (1) { ## else recurse tree
      if ($rank_hash{$parent} =~ m/^$rank_limit$/i) {
        $result = 1;
        last;
      } elsif ($parent == 1) {
        last;
      } else { ## climb the tree
        $parent = $nodes_hash{$parent};
        $parent_rank = $rank_hash{$parent};
      }
    }
  }
  return $result;
}
