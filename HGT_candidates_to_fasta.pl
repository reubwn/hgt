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
  -g|--groups            [FILE]   : groups file, e.g. OrthologousGroups.txt from OrthoFinder
  -x|--prefix            [FILE]   : filename prefix for outfile [default = INFILE]
  -v|--verbose                    : say more things [default: be quiet]
  -h|--help                       : prints this help message

EXAMPLES:

\n";

my ($in,$fasta,$groups,$prefix,$verbose,$help);

GetOptions (
  'in|i=s'              => \$in,
  'fasta|f=s'            => \$fasta,
  'groups|g:s'           => \$groups,
  'prefix|x:s'          => \$prefix,
  'verbose|v'           => \$verbose,
  'help|h'              => \$help,
);

die $usage if $help;
die $usage unless ($in && $fasta);

############################################## PARSE FASTA

print STDERR "[INFO] Parsing sequences from '$fasta'...\n";
my %seq_hash;
my $seqio = Bio::SeqIO -> new( -file => $fasta, -format => 'fasta' );
while (my $seq_obj = $seqio -> next_seq() ) {
  $seq_hash{$seq_obj->display_id()} = $seq_obj->seq(); ## key= seqname; val= seqstring
}
print STDERR "[INFO] Read ".commify((scalar(keys %seq_hash)))." sequences\n";

############################################## PARSE INFILE

open (my $IN, $in) or die $!;
while (<$IN>) {
  
}

############################################# SUBS

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
