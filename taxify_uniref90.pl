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

my %dbxref;
my ($seqid,$taxname,$taxid,$ncbitaxid) = ("NULL","NULL","NULL");
my $entry_begin = '<entry id="';
my $taxname_searchfor = '<property type="common taxon" value="';
my $taxid_searchfor = '<property type="common taxon ID" value="';
my $UniParc_searchfor = '<dbReference type="UniParc ID" id="';
my $UniRef_searchfor = '<property type="UniRef';
my $NCBI_searchfor = '<property type="NCBI taxonomy" value="';
my $entry_end = '</entry>';

open (my $IN, $infile) or die $!;
while (my $line = <$IN>) {
  chomp $line;
  if ($line =~ m/^\Q$entry_begin\E(\w+)\"/) {
    $seqid = $1;
  } elsif ($line =~ m/^\Q$taxname_searchfor\E(.+)\"/) {
    $taxname = $1;
  } elsif ($line =~ m/^\Q$taxid_searchfor\E(\d+)\"/) {
    $taxid = $1;
  } elsif ($line =~ m/^\Q$UniParc_searchfor\E(\w+)\"/) {
    $dbxref{$1} = ();
  } elsif ($line =~ m/^\Q$UniRef_searchfor\E.+(UniRef\w+)\"/) {
    $dbxref{$1} = ();
  } elsif ($line =~ m/\Q$NCBI_searchfor\E(\d+)\"/) {
    $ncbitaxid = $1;
  } elsif ($line =~ m/\Q$entry_end\E/) {
    if ($seqid eq "NULL" || $taxid eq "NULL" || $ncbitaxid eq "NULL") {
      print "[WARN] $seqid $taxid $ncbitaxid\n";
    } else {
      print join ("\t",
        $seqid,
        $taxid,
        $ncbitaxid,
        $taxname,
        (nsort keys %dbxref),
        "\n"
      );
      ($seqid,$taxid,$ncbitaxid,$taxname) = ("NULL","NULL","NULL","NULL");
      %dbxref = ();
    }
  }
}
