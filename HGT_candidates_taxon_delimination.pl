#!/usr/bin/env perl

## author: reubwn July 2017

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;
use Sort::Naturally;
use Data::Dumper qw(Dumper);

my $usage = "
SYNOPSIS

OUTPUTS

OPTIONS:
  -i|--in              [FILE]   : taxified diamond/BLAST results file [required]
  -p|--path            [STRING] : path to dir/ containing tax files [required]
  -l|--list            [FILE]   : list of HGT candidate gene names, one per line [required]
  -t|--taxid_threshold [INT]    : NCBI taxid to recurse up to; i.e., threshold taxid to define 'ingroup' [default = 33208 (Metazoa)]
  -k|--taxid_skip      [INT]    : NCBI taxid to skip; hits to this taxid will not be considered in any calculations of support
  -v|--verbose                  : say more things [default: be quiet]
  -h|--help                     : prints this help message
\n";

my ($infile,$candidates,$listfile,$path,$nodesfile,$namesfile,$mergedfile,$nodesDBfile,$prefix,$noheader,$help);
my $taxid_threshold = 33208;
my $taxid_skip = 0; ## default is 0, which is not a valid NCBI taxid and should not affect the tree recursion

GetOptions (
  'i|in=s'              => \$infile,
  'l|list=s'            => \$listfile,
  'p|path:s'            => \$path,
  'o|nodes:s'           => \$nodesfile,
  'a|names:s'           => \$namesfile,
  'm|merged:s'          => \$mergedfile,
  'n|nodesDB:s'         => \$nodesDBfile,
  't|taxid_threshold:i' => \$taxid_threshold,
  'k|taxid_skip:i'      => \$taxid_skip,
  'e|nohead'            => \$noheader,
  'h|help'              => \$help,
);

die $usage if $help;
die $usage unless ($infile && $listfile && $path);

## outfile
my $outfile = "$infile.HGT_candidates.taxon_delimitation";
open (my $OUT, ">$outfile") or die "[ERROR] Cannot open outfile '$outfile': $!\n\n";
print $OUT join (
  "\t",
  "QUERY",
  "UNIREF90",
  "PERCID",
  "BITSCORE",
  "TAXID",
  "SUPERKING",
  "KINGDOM",
  "PHYLUM",
  "CLASS",
  "ORDER",
  "FAMILY",
  "GENUS",
  "SPECIES",
  "\n"
) unless ($noheader);

############################################## PARSE NODES
## parse nodes and names:
my (%nodes_hash, %names_hash, %rank_hash);
if ($path) {
  print STDERR "[INFO] Building taxonomy databases from tax files in '$path'...";
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
} elsif ($nodesfile && $namesfile) {
  print STDERR "[INFO] Building taxonomy databases from '$nodesfile' and '$namesfile'...";
  open(my $NODES, $nodesfile) or die $!;
  while (<$NODES>) {
    chomp;
    next if /\#/;
    my @F = map { s/^\s+|\s+$//gr } split (/\|/, $_); ## split nodes.dmp file on \s+|\s+ regex
    $nodes_hash{$F[0]} = $F[1]; ## key= child taxid; value= parent taxid
    $rank_hash{$F[0]} = $F[2]; ## key= taxid; value= rank
  }
  close $NODES;
  open (my $NAMES, $namesfile) or die $!;
  while (<$NAMES>) {
    chomp;
    next if /\#/;
    my @F = map { s/^\s+|\s+$//gr } split (/\|/, $_);
    $names_hash{$F[0]} = $F[1] if ($F[3] eq "scientific name"); ## key= taxid; value= species name
  }
  close $NAMES;
  if ($mergedfile) {
    open (my $MERGED, $mergedfile) or die $!;
    while (<$MERGED>) {
      chomp;
      next if /\#/;
      my @F = map { s/^\s+|\s+$//gr } split (/\|/, $_);
      $nodes_hash{$F[0]} = $F[1]; ## key= old taxid; value= new taxid
      ## this will behave as if old taxid is a child of the new one, which is OK I guess
    }
  }
} elsif ($nodesDBfile) {
  print STDERR "[INFO] Building taxonomy databases from '$nodesDBfile'...";
  open(my $NODES, $nodesDBfile) or die $!;
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
print STDERR " done\n";
print STDERR "[INFO] Nodes parsed: ".scalar(keys %nodes_hash)."\n";
print STDERR "[INFO] Threshold taxid set to '$taxid_threshold' ($names_hash{$taxid_threshold})\n";
print STDERR "[INFO] INGROUP set to '$names_hash{$taxid_threshold}'; OUTGROUP is therefore 'non-$names_hash{$taxid_threshold}'\n";
if ($taxid_skip) {
  print STDERR "[INFO] Skipping any hits to taxid '$taxid_skip' ($names_hash{$taxid_skip})\n";
} else {
  print STDERR "[WARN] Taxid to skip (-k) is not set! Suggest setting -k to the taxid of the phylum your organism comes from.\n";
}

############################################# PARSE list

print STDERR "[INFO] Getting HGT candidate names from '$listfile'...\n";
my (%candidate_list);
open (my $LIST, $listfile) or die "[ERROR] Cannot open '$listfile': $!\n\n";
while (<$LIST>) {
  chomp;
  $candidate_list{$_} = ();
}
print STDERR "[INFO] Read ".scalar(keys %candidate_list)." HGT candidate names from '$listfile'\n";

############################################ PARSE DIAMOND FILE

print STDERR "[INFO] Parsing hit info from '$infile'...\n";
my (%seen,$processed);
open (my $DIAMOND, $infile) or die "[ERROR] Cannot open '$infile': $!\n\n";
LINE: while (my $line = <$DIAMOND>) {

  ## progress:
  $processed++;
  if ($processed % 1000 == 0){
    print STDERR "\r[INFO] Processed ".commify($processed)." queries...";$| = 1;
  }

  ## code:
  chomp $line;
  my @F = split (/\s+/, $line);
  if (exists($seen{$F[0]})) {
    next LINE;
  } else {
    if (exists($candidate_list{$F[0]})) {
      if ($F[12] !~ m/\d+/) {
        next LINE;
      } elsif (check_taxid_has_parent($F[12]) == 1) {
        next LINE;
      } elsif ( tax_walk($F[12], $taxid_skip) eq "ingroup" ) { ## do not include any hits to within taxid $taxid_skip
        next LINE;
      } elsif ( tax_walk($F[12]) eq "unassigned" ) {
        next LINE;
      } else {
        print $OUT join (
          "\t",
          $F[0],
          $F[1],
          $F[2],
          $F[11],
          $F[12],
          tax_walk_to_get_rank_to_species($F[12]),
          "\n"
        );
        $seen{$F[0]} = ();
        next LINE;
      }
    } else {
      next LINE;
    }
  }
}
print STDERR "[INFO] Finished on ".`date`."\n";

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
  $result =~ s/\s+/\_/g; ## replace spaces in some names with underscores
  $result =~ s/\;/\t/g; ## then replace ; with tabs
  return $result;
}
