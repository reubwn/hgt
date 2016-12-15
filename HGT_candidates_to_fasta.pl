#!/usr/bin/env perl

## author: reubwn Nov 2016

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;
use Sort::Naturally;
#use File::Grep qw( fgrep fmap fdo );

my $usage = "
SYNOPSIS:

1. Run 'makeblastdb -in uniref90.fasta -dbtype prot -parse_seqids' to allow for rapid sequence retrieval from large fasta DB like UniRef90

OUTPUTS:

OPTIONS:
  -i|--in              [FILE]   : HGT_candidates.txt file [required]
  -d|--diamond         [FILE]   : diamond/BLAST results file [required]
  -f|--fasta           [FILE]   : diamond/BLAST database fasta file, e.g. UniRef90.fasta [required]
  -p|--path            [STRING] : path to dir/ containing tax files
  -g|--groups          [FILE]   : groups file, e.g. OrthologousGroups.txt from OrthoFinder
  -t|--taxid_threshold [INT]    : NCBI taxid to recurse up to; i.e., threshold taxid to define 'ingroup' [default = 33208 (Metazoa)]
  -x|--prefix          [FILE]   : filename prefix for outfile [default = INFILE]
  -v|--verbose                  : say more things [default: be quiet]
  -h|--help                     : prints this help message

EXAMPLES:

\n";

my ($in,$diamond,$fasta,$path,$groups,$prefix,$verbose,$help);
my $taxid_threshold = 33208;

GetOptions (
  'in|i=s'              => \$in,
  'fasta|f=s'            => \$fasta,
  'diamond|d=s'      => \$diamond,
  'path|p=s'         => \$path,
  'groups|g:s'           => \$groups,
  'prefix|x:s'          => \$prefix,
  'verbose|v'           => \$verbose,
  'help|h'              => \$help,
);

die $usage if $help;
die $usage unless ($in && $fasta && $path && $diamond);

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

# print STDERR "[INFO] Building sequence database from '$fasta'...\n";
# my %seq_hash;
# my $processed = 0;
# my $seqio = Bio::SeqIO -> new( -file => $fasta, -format => 'fasta' );
# while (my $seq_obj = $seqio -> next_seq() ) {
#   $seq_hash{$seq_obj->display_id()} = $seq_obj->seq(); ## key= seqname; val= seqstring
#   ## progress
#   $processed++;
#   if ($processed % 1000 == 0){
#     print STDERR "\r[INFO] Read ".commify($processed)." sequences...";
#     $| = 1;
#   }
# }
# print STDERR "[INFO] Read ".commify((scalar(keys %seq_hash)))." sequences\n";

############################################## GET HGT CANDIDATE HITS

my (%hits_name_map, %hits_hash);

## get HGT candidates
chomp ( my @keys = `cut -f1 $in` );
my %hgt_candidates = map { $_ => 1 } @keys;
print "[INFO] Number of HGT candidates: ".commify((scalar(keys %hgt_candidates)))."\n";
open (my $DIAMOND, $diamond) or die $!;
while (<$DIAMOND>) {
  my @F = split (/\s+/, $_);
  if ( $hgt_candidates{$F[0]} ) {

    ## modify hit name based on taxid
    my ($hit_category,$new_hit_name); ## ingroup or outgroup
    if ( tax_walk($F[12]) eq "ingroup" ) {
      $new_hit_name = join ("_", $F[1], "IN");
    } elsif ( tax_walk($F[12]) eq "outgroup" ) {
      $new_hit_name = join ("_", $F[1], "OUT");
    } else {
      next;
    }
    $hits_name_map{$F[1]} = $new_hit_name; ## key= UniRef90 name; val= suffixed with IN|OUT
    push @{ $hits_hash{$F[0]} }, $F[1]; ## key= query name; val= [array of UniRef90 hit ids]
  }
}

############################################## PARSE INFILE

open (my $IN, $in) or die $!;
while (<$IN>) {
  chomp;
  next if ($_ =~ m/^\#/);
  my @F = split (/\s+/, $_);
  my @uniref_hits = @{ $hits_hash{$F[0]} };
  my $blastdbcmd_query_string = join (",", @{ $hits_hash{$F[0]} });

  my $file_name = $F[0];
  $file_name =~ s/\|/\_/;
  #open (my $FA, "<", $file_name) or die $!;
  if ( system ("blastdbcmd -db $fasta -dbtype prot -outfmt \%f -entry $blastdbcmd_query_string > $file_name") !=0 ) {
    die "[ERROR] blastdbcmd command did not work\n";
  }
  # foreach my $hit ( @uniref_hits ) {
  #   print $FA "\>$hits_name_map{$hit}\n$seq_hash{$hit}\n";
  # }
  #close $FA;

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
