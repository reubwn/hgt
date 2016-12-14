#!/usr/bin/env perl

## author: reubwn Nov 2016

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;
use Sort::Naturally;

my $usage = "
SYNOPSIS:

OUTPUTS:

OPTIONS:
  -i|--in                [FILE]   : *HGT_candidates.txt file [required]
  -f|--fasta             [FILE]   : protein fasta file, e.g. UniRef90.fasta [required]
  -p|--path              [STRING] : path to dir/ containing tax files
  -g|--groups            [FILE]   : groups file, e.g. OrthologousGroups.txt from OrthoFinder
  -x|--prefix            [FILE]   : filename prefix for outfile [default = INFILE]
  -v|--verbose                    : say more things [default: be quiet]
  -h|--help                       : prints this help message

EXAMPLES:

\n";

my ($in,$fasta,$path,$groups,$prefix,$verbose,$help);

GetOptions (
  'in|i=s'              => \$in,
  'fasta|f=s'            => \$fasta,
  'path|p=s'         => \$path,
  'groups|g:s'           => \$groups,
  'prefix|x:s'          => \$prefix,
  'verbose|v'           => \$verbose,
  'help|h'              => \$help,
);

die $usage if $help;
die $usage unless ($in && $fasta && $path);

############################################## PARSE NODES

## parse nodes and names:
my (%nodes_hash, %names_hash, %rank_hash);
if ($path) {
  print STDERR "[INFO] Building taxonomy databases from tax files in '$path'...\n";
  open(my $NODES, "$path/nodes.dmp") or die "[ERROR] nodes.dmp not found in $path: $!\n";
  while (<$NODES>) {
    chomp;
    next if /\#/;
    my @F = map { s/^\s+|\s+$//gr } split (/\|/, $_); ## split nodes.dmp file on \s+|\s+ regex
    $nodes_hash{$F[0]} = $F[1]; ## key= child taxid; value= parent taxid
    $rank_hash{$F[0]} = $F[2]; ## key= taxid; value= rank
  }
  close $NODES;
  open (my $NAMES, "$path/names.dmp") or die "[ERROR] names.dmp not found in $path: $!\n";
  while (<$NAMES>) {
    chomp;
    next if /\#/;
    my @F = map { s/^\s+|\s+$//gr } split (/\|/, $_);
    $names_hash{$F[0]} = $F[1] if ($F[3] eq "scientific name"); ## key= taxid; value= species name
  }
  close $NAMES;
  if (-e "$path/merged.dmp") {
    open (my $MERGED, "$path/merged.dmp") or die "[ERROR] merged.dmp not found in $path: $!\n";
    while (<$MERGED>) {
      chomp;
      next if /\#/;
      my @F = map { s/^\s+|\s+$//gr } split (/\|/, $_);
      $nodes_hash{$F[0]} = $F[1]; ## key= old taxid; value= new taxid
      ## this will behave as if old taxid is a child of the new one, which is OK I guess
    }
  }
} else { die "[ERROR] Path '$path' not found\n" }

############################################## PARSE FASTA

print STDERR "[INFO] Building sequence database from '$fasta'...\n";
my %seq_hash;
my $seqio = Bio::SeqIO -> new( -file => $fasta, -format => 'fasta' );
while (my $seq_obj = $seqio -> next_seq() ) {
  $seq_hash{$seq_obj->display_id()} = $seq_obj->seq(); ## key= seqname; val= seqstring
}
print STDERR "[INFO] Read ".commify((scalar(keys %seq_hash)))." sequences\n";

############################################## PARSE INFILE

open (my $IN, $in) or die $!;
while (<$IN>) {
  chomp;
  next if ($_ =~ m/^\#/);
  my @F = split (/\s+/, $_);

}
close $IN;

############################################# SUBS

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
