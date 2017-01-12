# hgt
scripts for the analysis of HGT in genome sequence data

## diamond_to_HGT_candidates.pl

### Synopsis

This script calculates two measures of "support" for identifying putative HGT candidate genes:

1. **Alien Index:** AI is a measure of how well a given sequence matches to one set of taxa (eg. Metazoa) relative to another, mutually exclusive set of taxa (eg. non-Metazoa). It uses only the best-hit to each category to calculate this: ```log((Best E-value for INGROUP) + 1e-200) - log((Best E-value for OUTGROUP) + 1e-200)```. See [Gladyshev et al. (2008) Science 320:1210 Suppl Mat](http://science.sciencemag.org/content/suppl/2008/05/29/320.5880.1210.DC1/Gladyshev.SOM.pdf) for more details.
2. **Hit Support Index:** The second measure uses information from all hits, as opposed to just the top 2 as above. Each query is designated *ingroup* (eg. Metazoa) versus *outgroup* (eg. non-Metazoa) based on the sum of bitscores across all hits; the category with the highest bitscore wins. Each individual hit is also designated *ingroup* or *outgroup*, based on its taxonomy. To be classified as a putative HGT candidate, a given protein must (a) be designated as *outgroup* and (b) have support for this designation from a given proportion of the rest of the hits (default is set to >=90%). See [Koutsovoulos et al. (2016) PNAS 113:5053](http://www.pnas.org/content/113/18/5053.abstract) for more details.

### Prerequisites

There's a bit of setup required before the script can be run:

1. **Download NCBI Taxonomy:** These scripts rely on the NCBI taxonomy databases to assign taxonomic info to proteins. In particular the files "nodes.dmp" and "names.dmp" are required, while "merged.dmp" is recommended. Download them like:

   ```
   >> wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
   >> tar xzvf taxdump.tar.gz
   >> export TAXPATH=/path/to/taxdump/
   ```

   It's also possible to use the "nodesDB.txt" file from the excellent [Blobtools](https://blobtools.readme.io/docs/what-is-blobtools) software.

2. **Taxify your BLAST/Diamond file:** Diamond is great for speed, but adding the taxid information to each hit requires an additional step. See [this Gist](https://gist.github.com/sujaikumar/9ad04e62449a2d7025b17144de67038b) by Sujai Kumar on how to set this up for the UniRef90 database.

   A typical Diamond script might then look like:

   ```
   ## run diamond
   >> diamond blastp \
   --sensitive \
   --index-chunks 1 \
   -e 1e-5 \
   -p 100 \
   -q $QUERY \
   -d $DB \
   -a $PREFIX

   ## taxify
   >> diamond view -a ${PREFIX}.daa \
   | perl -lne '
   BEGIN{open UT, "</path/to/uniref90.taxlist" or die $!; while (<UT>) { $ut{$1}=$2 if /^(\S+)\t(\S+)$/ } }
   {print "$_\t$ut{$1}" if /^\S+\t(\S+)/ and exists $ut{$1}}' \
   > ${PREFIX}.daa.taxid
   ```

   The taxified file has the usual 12 columns but with an additional 13th column containing the taxid of the hit protein.

### Options

Type ```diamond_to_HGT_candidates.pl -h``` to see help and options.

```
OPTIONS:
-i|--in                [FILE]   : tab formatted Diamond output file [required]
-p|--path              [STRING] : path to dir/ containing tax files [one of -p || (-o && -a) || -n is required]
-o|--nodes             [FILE]   : path to nodes.dmp
-a|--names             [FILE]   : path to names.dmp
-m|--merged            [FILE]   : path to merged.dmp
-n|--nodesDB           [FILE]   : nodesDB.txt file from blobtools
-g|--gff               [FILE]   : path to augustus-formatted GFF file [TODO]
-t|--taxid_threshold   [INT]    : NCBI taxid to recurse up to; i.e., threshold taxid to define 'ingroup' [default = 33208 (Metazoa)]
-k|--taxid_skip        [INT]    : NCBI taxid to skip; hits to this taxid will not be considered in any calculations of support
-r|--scoring           [STRING] : scoring strategy for calculating bestsum bitscore: 'sum' or 'individual' [default = 'sum']
-s|--support_threshold [FLOAT]  : Secondary Hits Support threshold for considering HGT candidates [default = 90\%]
-l|--alien_threshold   [INT]    : Alien Index threshold for considering HGT candidates (default >= 45)
-e|--evalue_column     [INT]    : define evalue column for --in (first column = 1) [default: 11]
-b|--bitscore_column   [INT]    : define bitscore column for --in (first column = 1) [default: 12]
-c|--taxid_column      [INT]    : define taxid column for --in (first column = 1) [default: 13]
-d|--delimiter         [STRING] : define delimiter to split --in (specify 'diamond' for Diamond files (\"\\s+\") or 'blast' for BLAST files (\"\\t\")) [default: diamond]
-x|--prefix            [FILE]   : filename prefix for outfile [default = INFILE]
-H|--header                     : don't print header [default: do print it]
-v|--verbose                    : say more things [default: be quiet]
-h|--help                       : prints this help message
```
