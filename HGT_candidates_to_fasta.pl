#!/usr/bin/env perl

## author: reubwn Nov 2016

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;
use Sort::Naturally;
use Data::Dumper qw(Dumper);
#use File::Grep qw( fgrep fmap fdo );

my $usage = "
SYNOPSIS:

1. Run 'makeblastdb -in uniref90.fasta -dbtype prot -parse_seqids' to allow for rapid sequence retrieval from large fasta DB like UniRef90

OUTPUTS:

OPTIONS:
  -i|--in              [FILE]   : taxified diamond/BLAST results file [required]
  -c|--candidates      [FILE]   : *.HGT_candidates.txt file [required]
  -u|--uniref90        [FILE]   : diamond/BLAST database fasta file, e.g. UniRef90.fasta [required]
  -f|--fasta           [FILE]   : fasta file of query proteins [required]
  -p|--path            [STRING] : path to dir/ containing tax files
  -g|--groups          [FILE]   : groups file, e.g. OrthologousGroups.txt from OrthoFinder
  -t|--taxid_threshold [INT]    : NCBI taxid to recurse up to; i.e., threshold taxid to define 'ingroup' [default = 33208 (Metazoa)]
  -x|--prefix          [FILE]   : filename prefix for outfile [default = INFILE]
  -v|--verbose                  : say more things [default: be quiet]
  -h|--help                     : prints this help message

EXAMPLES:

\n";

my ($in,$candidates,$uniref90,$fasta,$path,$groups,$prefix,$verbose,$help);
my $taxid_threshold = 33208;

GetOptions (
  'in|i=s'        => \$in,
  'cadidates|c=S' => \$candidates,
  'uniref90|u=s'  => \$uniref90,
  'fasta|f=s'     => \$fasta,
  'path|p=s'      => \$path,
  'groups|g:s'    => \$groups,
  'prefix|x:s'    => \$prefix,
  'verbose|v'     => \$verbose,
  'help|h'        => \$help,
);

die $usage if $help;
die $usage unless ($in && $uniref90 && $fasta && $path && $cadidates);

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
  print STDERR "[INFO] Done $path/nodes.dmp\n";
  open (my $NAMES, "$path/names.dmp") or die "[ERROR] names.dmp not found in $path: $!\n";
  while (<$NAMES>) {
    chomp;
    next if /\#/;
    my @F = map { s/^\s+|\s+$//gr } split (/\|/, $_);
    $names_hash{$F[0]} = $F[1] if ($F[3] eq "scientific name"); ## key= taxid; value= species name
  }
  close $NAMES;
  print STDERR "[INFO] Done $path/names.dmp\n";
  if (-e "$path/merged.dmp") {
    open (my $MERGED, "$path/merged.dmp") or die "[ERROR] merged.dmp not found in $path: $!\n";
    while (<$MERGED>) {
      chomp;
      next if /\#/;
      my @F = map { s/^\s+|\s+$//gr } split (/\|/, $_);
      $nodes_hash{$F[0]} = $F[1]; ## key= old taxid; value= new taxid
      ## this will behave as if old taxid is a child of the new one, which is OK I guess
    }
    close $MERGED;
    print STDERR "[INFO] Done $path/merged.dmp\n";
  }
} else { die "[ERROR] Path '$path' not found\n" }

############################################## PARSE FASTA

print STDERR "[INFO] Building sequence database from '$fasta'...\n";
my %seq_hash;
my $processed = 0;
my $seqio = Bio::SeqIO -> new( -file => $fasta, -format => 'fasta' );
while (my $seq_obj = $seqio -> next_seq() ) {
  $seq_hash{$seq_obj->display_id()} = $seq_obj->seq(); ## key= seqname; val= seqstring

}
print STDERR "[INFO] Read ".commify((scalar(keys %seq_hash)))." sequences\n";

############################################## GET HGT CANDIDATE HITS

my (%hits_name_map, %hits_hash);

## get HGT candidates
chomp ( my @keys = `cut -f1 $candidates` );
my %hgt_candidates = map { $_ => 1 } @keys;
print "[INFO] Number of HGT candidates: ".commify((scalar(keys %hgt_candidates)))."\n";
open (my $DIAMOND, $in) or die "Cannot open file '$in': $!\n";
while (<$DIAMOND>) {
  my @F = split (/\s+/, $_);
  if ( $hgt_candidates{$F[0]} ) {

    ## modify hit name based on taxid
    my ($hit_category,$new_hit_name); ## ingroup or outgroup
    if ($F[12] !~ /\d+/) {
      ## invalid taxid (not a number)
      next;
    } elsif (check_taxid_has_parent($F[12]) == 1) {
      ## taxid doesn't have a parent
      next;
    } elsif ( tax_walk($F[12]) eq "unassigned" ) {
      ## taxid is unassigned
      next;
    } elsif ( tax_walk($F[12]) eq "ingroup" ) {
      my $phylum = tax_walk_to_get_rank_to_phylum($F[12]);
      $new_hit_name = join ("_", $F[1], "IN", $phylum);
      $hits_name_map{$F[1]} = $new_hit_name; ## key= UniRef90 name; val= suffixed with IN|OUT
      push @{ $hits_hash{$F[0]} }, $F[1]; ## key= query name; val= [array of UniRef90 hit ids]
      next;
    } elsif ( tax_walk($F[12]) eq "outgroup" ) {
      my $phylum = tax_walk_to_get_rank_to_phylum($F[12]);
      $new_hit_name = join ("_", $F[1], "OUT", $phylum);
      $hits_name_map{$F[1]} = $new_hit_name; ## key= UniRef90 name; val= suffixed with IN|OUT
      push @{ $hits_hash{$F[0]} }, $F[1]; ## key= query name; val= [array of UniRef90 hit ids]
      next;
    }
  }
}
close $DIAMOND;

print Dumper \%hits_hash if $verbose;

############################################## PARSE CANDIDATES

print STDERR "[INFO] Getting sequences from database '$uniref90'\n";
print STDERR "[INFO] Processing HGT_candidates file...\n";
open (my $CANDIDATES, $candidates) or die "Cannot open file '$candidates': $!\n";
while (<$CANDIDATES>) {
  chomp;
  next if ($_ =~ m/^\#/);
  my @F = split (/\s+/, $_);

  ## make query string for blastdbcmd
  my $blastdbcmd_query_string = join (",", @{ $hits_hash{$F[0]} });

  ## make filename based on query name:
  (my $file_name = $F[0]) =~ s/\|/\_/;
  open (my $FA, ">", $file_name) or die $!;
  open (my $CMD, "blastdbcmd -db $uniref90 -dbtype prot -outfmt \%f -entry $blastdbcmd_query_string |") or die $!;
  print $FA "\>$F[0]\n$seq_hash{$F[0]}\n";
  while (<$CMD>) {
    if ($_ =~ /^>/) {
      $_ =~ s/\>//;
      $_ =~ s/\s+//g;
      #print STDERR "$_-->$hits_name_map{$_}\n";
      print $FA "\>$hits_name_map{$_}\n";
    } else {
      print $FA $_;
    }
  }
  close $CMD;
  close $FA;

  ## progress
  $processed++;
  if ($processed % 10 == 0){
     print STDERR "\r[INFO] Processed ".commify($processed)." queries...";
     $| = 1;
  }
}
close $CANDIDATES;
print STDERR "\n[INFO] Finished\n\n";

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
