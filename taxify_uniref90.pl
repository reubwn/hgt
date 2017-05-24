#!/usr/bin/env perl

## author: reubwn

use strict;
use warnings;

use Getopt::Long;
use Sort::Naturally;
use Data::Dumper qw(Dumper);
use List::Util qw(reduce sum min max);

my ($infile,$path);
my $taxid_skip = 42241;
my $taxid_threshold = 33208; ##metazoa

GetOptions (
  'i|in=s' => \$infile,
  'path|p:s'  => \$path,
);

my (%uniref);
my ($seqid,$taxid,$ncbitaxid) = ("NULL","NULL","NULL");
my $seqid_searchfor = '<entry id="';
my $taxid_searchfor = '<property type="common taxon ID" value="';
my $UniParc_searchfor = '<dbReference type="UniParc ID" id="';
# my $UniRef_searchfor =
my $NCBI_searchfor = '<property type="NCBI taxonomy" value="';
my $endquery_searchfor = '</entry>';

open (my $IN, $infile) or die $!;
while (my $line = <$IN>) {
  chomp $line;
  if ($line =~ m/^\Q$seqid_searchfor\E(\w+)\"/) {
    $seqid = $1;
  } elsif ($line =~ m/^\Q$taxid_searchfor\E(\d+)\"/) {
    $taxid = $1;
  } elsif ($line =~ m/\Q$NCBI_searchfor\E(\d+)\"/) {
    $ncbitaxid = $1;
  } elsif ($line =~ m/\Q$endquery_searchfor\E/) {
    print "[WARN] $seqid $taxid $ncbitaxid\n" if ($seqid eq "NULL" || $taxid eq "NULL" || $ncbitaxid eq "NULL");
    print "$seqid\t$taxid\t$ncbitaxid\n";
    ($seqid,$taxid,$ncbitaxid) = ("NULL","NULL","NULL");
  }
}
