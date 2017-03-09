# hgt
scripts for the analysis of HGT in genome sequence data

## diamond_to_HGT_candidates.pl

### Synopsis

This script calculates two measures of "support" for identifying putative HGT candidate genes:

1. **Alien Index:** AI is a measure of how well a given sequence matches to one set of taxa (eg. Metazoa) relative to another, mutually exclusive set of taxa (eg. non-Metazoa). It uses only the best-hit to each category to calculate this: ```log((Best E-value for INGROUP) + 1e-200) - log((Best E-value for OUTGROUP) + 1e-200)```. See [Gladyshev et al. (2008) Science 320:1210 Suppl Mat](http://science.sciencemag.org/content/suppl/2008/05/29/320.5880.1210.DC1/Gladyshev.SOM.pdf) for more details.
2. **Hit Support Index:** The second measure uses information from all hits, as opposed to just the top 2 as above. Each query is designated *ingroup* (eg. Metazoa) versus *outgroup* (eg. non-Metazoa) based on the sum of bitscores across all hits; the category with the highest bitscore wins. Each individual hit is also designated *ingroup* or *outgroup*, based on its taxonomy. To be classified as a putative HGT candidate, a given protein must (a) be designated as *outgroup* and (b) have support for this designation from a given proportion of the rest of the hits (default is set to >=90%). This approach may be less prone to errors associated with contamination issues in genome data, as it does not rely on evidence of providence from only the top hit. See [Koutsovoulos et al. (2016) PNAS 113:5053](http://www.pnas.org/content/113/18/5053.abstract) for more details.

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
-s|--support_threshold [FLOAT]  : Secondary Hits Support threshold for considering HGT candidates [default = 90%]
-l|--alien_threshold   [INT]    : Alien Index threshold for considering HGT candidates (default >= 45)
-e|--evalue_column     [INT]    : define evalue column for --in (first column = 1) [default: 11]
-b|--bitscore_column   [INT]    : define bitscore column for --in (first column = 1) [default: 12]
-c|--taxid_column      [INT]    : define taxid column for --in (first column = 1) [default: 13]
-d|--delimiter         [STRING] : define delimiter to split --in (specify 'diamond' for Diamond files ("\s+") or 'blast' for BLAST files ("\t")) [default: diamond]
-x|--prefix            [FILE]   : filename prefix for outfile [default = INFILE]
-H|--header                     : don't print header [default: do print it]
-v|--verbose                    : say more things [default: be quiet]
-h|--help                       : prints this help message
```

### Examples

Assuming you have a taxified Diamond file with the hit taxids in the 13th column, and have downloaded the NCBI taxdump files to a location specified by $TAXPATH, running the script with default thresholds is as simple as:

```
>> diamond_to_HGT_candidates.pl -i diamond_results.daa.taxid -p $TAXPATH
```

To exclude hits to organisms that are closely related to your focal organism, use -k option (recommended); eg. if you're looking at HGT in a nematode genome:

```
>> diamond_to_HGT_candidates.pl -i diamond_results.daa.taxid -p $TAXPATH -k 6231
```

### Outputs

The program outputs 3 files, suffixed with the tags:

1. **HGT_results:** AI and HSI scores for all query proteins.
2. **HGT_candidates:** All queries which show evidence from *either* AI *or* HSI above the specifed thresholds (defaults are >45 for AI, >90% for HSI). Can be further parsed to get the intersection of these two lists.
3. **HGT_warnings:** Any non-fatal warnings detected during the run. Worth checking, but most warnings can probably safely be ignored.

## HGT_candidates_to_fasta.pl

### Synopsis

This script takes the results from the "HGT_candidates" output file, generated using ```diamond_to_HGT_candidates.pl```, and outputs a fasta file containing the sequence data for each HGT candidate + its hits from UniRef90, each of which has been annotated with *ingroup* or *outgroup* and some taxonomy information up to Phylum level. These fasta files can then be aligned, and phylogenies constructed from these alignments.

### Prerequisites

1. The real utility of the script is that it uses ```blastdbcmd``` from the BLAST suite for ultra-quick sequence retrieval. This requires you generate a BLAST database of the UniRef90 fasta sequences **with the -parse_seqids option specifed**, otherwise it can't retrieve the sequence data itself using the fasta header:

   ```
   >> makeblastdb -in uniref90.fasta -dbtype prot -parse_seqids
   ```

   Also requires that ```blastdbcmd``` is in your ```$PATH```.


### Options

Type ```HGT_candidates_to_fasta.pl -h``` to see the options (note some are not implemented yet...):

```
OPTIONS:
  -i|--in              [FILE]   : taxified diamond/BLAST results file [required]
  -c|--candidates      [FILE]   : *.HGT_candidates.txt file [required]
  -u|--uniref90        [FILE]   : diamond/BLAST database fasta file, e.g. UniRef90.fasta [required]
  -f|--fasta           [FILE]   : fasta file of query proteins [required]
  -p|--path            [STRING] : path to dir/ containing tax files [required]
  -t|--taxid_threshold [INT]    : NCBI taxid to recurse up to; i.e., threshold taxid to define 'ingroup' [default = 33208 (Metazoa)]
  -g|--groups          [FILE]   : groups file, e.g. OrthologousGroups.txt from OrthoFinder [TODO]
  -m|--mafft                    : run MAFFT alignment on each fasta file (default = no) [TODO]
  -x|--raxml                    : run RAxML phylogenetic reconstruction on each MAFFT alignment file [default = no] [TODO]
  -v|--verbose                  : say more things [default: be quiet]
  -h|--help                     : prints this help message
```

### Outputs

One fasta file per HGT candidate gene; each file is named after the focal species query gene name.

### Downstream analyses

1. Align fasta files using MAFFT (run from within dir of HGT candidate fasta files):

   ```
   >> mkdir ../mafft_alns
   >> for f in *.fasta; do echo $f; mafft --auto --quiet --thread 8 $f > ../mafft_alns/${f}.mafft; done
   ```

2. Construct phylogenies using RAxML (run from dir of alignments):

   ```
   >> cd ../mafft_alns/
   >> mkdir processed
   >> mkdir ../raxml_trees
   >> for file in *mafft; do
        echo $file;
        raxmlHPC-PTHREADS-AVX -f a -p 12345 -x 12345 -# 100 -m PROTGAMMAAUTO -T 16 -s ${file} -n ${file};
        raxmlHPC-PTHREADS-AVX -f b -t RAxML_bestTree.${file} -z RAxML_bootstrap.${file} -m PROTGAMMAAUTO -n RAxML_bestTree.${file};
        mv $file processed_files/;
        mv RAxML* ../raxml_trees/;
      done
   ```

   If it fails or you run out of cluster time you can just restart from within the mafft_alns dir, as the processed alignments should have been moved to processed/

## analyse_trees.R

Rscript to assess topological support for HGT candidates. For each tree, tests for monophyly with non-metazoan Eukayotes, Fungi, Plants, Bacteria and Archaea. Requires the {ape} R package.

### Inputs

Path to dir of (newick) trees.

### Outputs

Counts for monophyletic support for query with numerous groups; PDF of trees colour-coded by ingroup / outgroup categorisation.
