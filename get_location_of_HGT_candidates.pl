#!/usr/bin/env perl

## author: reubwn May 2017

use strict;
use warnings;

my $usage = "
SYNOPSIS
  Takes a '*.HGT_results' file and a GFF and returns an '*.HGT_locations' file, specifying the location on each chromosome of HGT candidates.

OUTPUTS
  A '*.HGT_locations' file and a map file, and a list of scaffolds with 'too many' HGT candidates encoded on them to be believable.

OPTIONS:
  -i|--in            [FILE]   : taxified diamond/BLAST results file [required]
  -r|--results       [FILE]   : *.HGT_results.txt file [required]
  -g|--gff           [FILE]   : GFF file [required]
  -s|--CHS_threshold [FLOAT]  : Consesus Hits Support threshold [default>=90\%]
  -u|--hU_threshold  [INT]    : hU threshold [default>=30]
  -v|--verbose                : say more things [default=quiet]
  -h|--help                   : prints this help message
\n";

my ($in,$candidates,$uniref90,$fasta,$path,$groups,$mafft,$raxml,$prefix,$verbose,$help);
my $taxid_threshold = 33208;
my $taxid_skip = 0; ## default is 0, which is not a valid NCBI taxid and should not affect the tree recursion
my $limit = 15;

GetOptions (
  'in|i=s'              => \$in,
  'cadidates|c=s'       => \$candidates,
  'uniref90|u=s'        => \$uniref90,
  'fasta|f=s'           => \$fasta,
  'path|p=s'            => \$path,
  'taxid_threshold|t:i' => \$taxid_threshold,
  'taxid_skip|k:i'      => \$taxid_skip,
  'limit|l:i'           => \$limit,
  'groups|g:s'          => \$groups,
  'mafft|m'             => \$mafft,
  'raxml|x'             => \$raxml,
  'verbose|v'           => \$verbose,
  'help|h'              => \$help,
);

die $usage if $help;
die $usage unless ($in && $candidates && $uniref90 && $fasta && $path);
die "Options $groups, $mafft and $raxml are not implemented yet!\n\n" if ($groups || $mafft || $raxml);
