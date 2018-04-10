#!/usr/bin/env perl

## author: reubwn Nov 2016

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;
use Sort::Naturally;
use Data::Dumper qw(Dumper);
use File::Path;

my $usage = "
SYNOPSIS
  This program takes the results from the 'diamond_to_HGT_candidates.pl' script and returns a fasta file for each HGT candidate gene, containing
    1. the focal HGT candidate protein and
    2. the 'taxified' hits from UniRef90 (or other db) to that protein

  Taxified means that each protein is annotated as either ingroup or outgroup, based on the partition specified (default = 'Metazoan' vs 'non-Metazoan'), and taxonomic classification up to phylum level
  The resultant fasta files can then be aligned and used for phylogenetic reconstruction.

PREREQUISITES
  The script utilises the rapid lookup of large fasta files by 'blastdbcmd'. This requires the '-parse_seqids' flag to be specified during the creation of the blastdb: 'makeblastdb -in uniref90.fasta -dbtype prot -parse_seqids'

OUTPUTS
  A bunch of fasta files, one per HGT candidate.
  Default behaviour is to fetch up to 15 ingroup and 15 outgroup sequences, if available, to avoid uneccessarily large and uninformative alignments (set by --limit)

OPTIONS:
  -i|--in              [FILE]   : taxified diamond/BLAST results file (accepts gzipped) (required)
  -c|--candidates      [FILE]   : HGT_candidates.txt file [required]
  -u|--uniref90        [FILE]   : diamond/BLAST database fasta file, e.g. UniRef90.fasta (required)
  -f|--fasta           [FILE]   : fasta file of query proteins (required)
  -p|--path            [STRING] : path to dir/ containing tax files (required)
  -t|--taxid_threshold [INT]    : NCBI taxid to recurse up to; i.e., threshold taxid to define 'ingroup' [33208]
  -k|--taxid_skip      [INT]    : NCBI taxid to skip; hits to this taxid will not be considered in any calculations of support
  -d|--dir             [DIR]    : dir name to put fasta files into ['outdir']
  -l|--limit                    : maximum number of ingroup / outgroup sequences to fetch (if available) [15]
  -a|--alienness                : format sequence headers for alienness phylogeny testing instead [no]
  -v|--verbose                  : say more things
  -h|--help                     : print help
\n";

my ($in,$candidates,$uniref90,$fasta,$path,$alienness,$groups,$mafft,$raxml,$prefix,$verbose,$help);
my $taxid_threshold = 33208;
my $taxid_skip = 0; ## default is 0, which is not a valid NCBI taxid and should not affect the tree recursion
my $outdir = "outdir";
my $limit = 15;

GetOptions (
  'i|in=s'              => \$in,
  'c|cadidates=s'       => \$candidates,
  'u|uniref90=s'        => \$uniref90,
  'f|fasta=s'           => \$fasta,
  'p|path=s'            => \$path,
  't|taxid_threshold:i' => \$taxid_threshold,
  'k|taxid_skip:i'      => \$taxid_skip,
  'l|limit:i'           => \$limit,
  'd|dir:s'             => \$outdir,
  'a|alienness'         => \$alienness,
  'g|groups:s'          => \$groups,
  'm|mafft'             => \$mafft,
  'x|raxml'             => \$raxml,
  'v|verbose'           => \$verbose,
  'h|help'              => \$help,
);

die $usage if $help;
die $usage unless ($in && $candidates && $uniref90 && $fasta && $path);
die "Options $groups, $mafft and $raxml are not implemented yet!\n\n" if ($groups || $mafft || $raxml);

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
  print STDERR "[INFO] + ".$path."nodes.dmp\n";
  open (my $NAMES, "$path/names.dmp") or die "[ERROR] names.dmp not found in $path: $!\n";
  while (<$NAMES>) {
    chomp;
    next if /\#/;
    my @F = map { s/^\s+|\s+$//gr } split (/\|/, $_);
    $names_hash{$F[0]} = $F[1] if ($F[3] eq "scientific name"); ## key= taxid; value= species name
  }
  close $NAMES;
  print STDERR "[INFO] + ".$path."names.dmp\n";
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
    print STDERR "[INFO] + ".$path."merged.dmp\n";
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
my %hits_limit = ( "in_limit" => 0, "out_limit" => 0 );

## get HGT candidates
chomp ( my @keys = `cut -f1 $candidates` );
my %hgt_candidates = map { $_ => 1 } @keys;
print "[INFO] Number of HGT candidates: ".commify((scalar(keys %hgt_candidates)))."\n";
if ($taxid_skip) {
  print "[INFO] Skipping any hits to taxid '$taxid_skip' ($names_hash{$taxid_skip})\n";
} else {
  print "[WARN] Taxid to skip (-k) is not set! Suggest setting -k to the taxid of the phylum your organism comes from.\n";
}
print "[INFO] Printing alienness format\n" if ($alienness);

my $DIAMOND;
if ($in =~ m/gz$/) {
  open ($DIAMOND, "zcat $in |") or die "Cannot open file '$in': $!\n"; ## open from Gzip
} else {
  open ($DIAMOND, $in) or die "Cannot open file '$in': $!\n";
}
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
    } elsif ( tax_walk($F[12], $taxid_skip) eq "ingroup" ) {
      ## skip entries to close relatives as this may "confuse" the assignment of query+outgroup monophyly tests downstream
      next;
    } elsif ( tax_walk($F[12]) eq "ingroup" ) {
      ## add up to 15 INGROUP sequences (avoids uneccessarily large and uninformative alignments)
      $hits_limit{$F[0]}{"in_limit"}++; ## increment in_limit up to max set by $limit
      if ( $hits_limit{$F[0]}{"in_limit"} > $limit ) {
        next;
      } else {
        my $phylum = tax_walk_to_get_rank_to_phylum($F[12]);
        if ($alienness) {
          $new_hit_name = "$F[1]\@TOI-$names_hash{$F[12]}";
        } else {
          $new_hit_name = join ("|", $F[1], "IN", $phylum);
        }
        $hits_name_map{$F[1]} = $new_hit_name; ## key= UniRef90 name; val= suffixed with IN|OUT
        #push @{ $hits_hash{$F[0]} }, $F[1]; ## key= query name; val= [array of UniRef90 hit ids]
        $hits_hash{$F[0]}{$F[1]} = 1; ## key= query name; value= { hash of UniRef90 hit ids }; use hash to eliminate repeated hit names
        next;
      }
    } elsif ( tax_walk($F[12]) eq "outgroup" ) {
      ## add up to 15 OUTGROUP sequences
      $hits_limit{$F[0]}{"out_limit"}++; ## increment out_limit up to max set by $limit
      if ( $hits_limit{$F[0]}{"out_limit"} > $limit ) {
        next;
      } else {
        if ($alienness) {
          my $phylum = tax_walk_to_get_rank_to_phylum($F[12], 1);
          $new_hit_name = "$F[1]\@$phylum";
        } else {
          my $phylum = tax_walk_to_get_rank_to_phylum($F[12]);
          $new_hit_name = join ("|", $F[1], "OUT", $phylum);
        }
        $hits_name_map{$F[1]} = $new_hit_name; ## key= UniRef90 name; val= suffixed with IN|OUT
        #push @{ $hits_hash{$F[0]} }, $F[1]; ## key= query name; val= [array of UniRef90 hit ids]
        $hits_hash{$F[0]}{$F[1]} = 1; ## key= query name; value= { hash of UniRef90 hit ids }
        next;
      }
    }
  }
}
close $DIAMOND;
print Dumper \%hits_hash if $verbose;

############################################## PARSE CANDIDATES

## make outdir to save results:
if (-d $outdir) {
  rmtree([ "$outdir" ]);
  mkdir $outdir;
} else {
  mkdir $outdir;
}

print STDERR "[INFO] Getting sequences from database '$uniref90'\n";
print STDERR "[INFO] Writing results to dir: $outdir\n";
print STDERR "[INFO] Processing HGT_candidates file...\n";
open (my $CANDIDATES, $candidates) or die "Cannot open file '$candidates': $!\n";
while (<$CANDIDATES>) {
  chomp;
  next if ($_ =~ m/^\#/);
  my @F = split (/\s+/, $_);

  ## make query string for blastdbcmd
  my $blastdbcmd_query_string = join (",", keys %{ $hits_hash{$F[0]} });

  ## make filename based on query name:
  (my $file_name = $F[0]) =~ s/\|/\_/;
  open (my $FA, ">", "$file_name.fasta") or die $!;
  open (my $CMD, "blastdbcmd -db $uniref90 -dbtype prot -outfmt \%f -entry $blastdbcmd_query_string |") or die $!;
  ## print hits from UniRef90
  while (<$CMD>) {
    if ($_ =~ /^>/) {
      $_ =~ s/\>//;
      $_ =~ s/\s.*//g; ## Duncan Berger suggests edit to: "$_ =~ s/\s.*//g;", I guess in case there are other characters after the whitespace? Check.
      print $FA "\>$hits_name_map{$_}\n";
    } else {
      print $FA $_;
    }
  }
  ## print query sequence LAST!
  if ($alienness) {
    print $FA "\>$F[0]\@StudiedOrganism\n$seq_hash{$F[0]}\n";
  } else {
    print $FA "\>$F[0]\n$seq_hash{$F[0]}\n";
  }
  close $CMD;
  close $FA;

  ## mv fasta to outdir:
  system("mv $file_name.fasta $outdir");

  ## progress
  $processed++;
  if ($processed % 10 == 0){
     print STDERR "\r[INFO] Processed ".commify($processed)." queries...";
     $| = 1;
  }
}
close $CANDIDATES;
print STDERR "\n[INFO] Finished on ".`date`."\n";

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
  my $superkingdom_only = $_[1];
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
  my $result;
  if ($superkingdom_only) {
    $result = $superkingdom;
  } else {
    $result = join ("|",$superkingdom,$kingdom,$phylum);
  }
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
