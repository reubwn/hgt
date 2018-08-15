# hgt
Scripts for the analysis of HGT in genome sequence data.

## Note

**Disclaimer.** These methods are predictive! Finding a candidate gene with high hU or AI or attached to a 'real' metazoan gene or with introns etc. does not necessarily mean that that gene has been acquired horizontally. Further tests would be required to provide further evidence for the evolutionary origin of each gene.

### Cite

If you use these scripts and find them useful, please cite [Nowell, R. W. et al. Comparative genomics of bdelloid rotifers: Insights from desiccating and nondesiccating species. PLoS Biol. 16, e2004830 (2018)](http://dx.doi.org/10.1371/journal.pbio.2004830).

The software also has a DOI:

[![DOI](https://zenodo.org/badge/76456664.svg)](https://zenodo.org/badge/latestdoi/76456664)

---

## diamond_to_HGT_candidates.pl

### Synopsis

This script analyses the output of Diamond/BLAST files and calculates 3 measures of support for identifying putative HGT candidate genes:

1. **HGT Index:** hU is a measure of how well a given sequence matches to one set of taxa (eg. Metazoa) relative to another, mutually exclusive set of taxa (eg. non-Metazoa). It uses best-hit bitscores and is defined: `(best-hit bitscore for OUTGROUP) - (best-hit bitsore for INGROUP)`. See [Boschetti et al. 2012](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003035) for more details.
2. **Alien Index:** AI is another measure based on E-values rather than bitscores: `log10((best-hit E-value for INGROUP) + 1e-200) - log10((best-hit E-value for OUTGROUP) + 1e-200)`. See [Gladyshev et al. 2008](http://science.sciencemag.org/content/suppl/2008/05/29/320.5880.1210.DC1/Gladyshev.SOM.pdf) for more details.
2. **Consensus Hit Support:** CHS uses information from all hits, as opposed to just the top 2 as for hU and AI. Each query is designated *ingroup* (eg. Metazoa) versus *outgroup* (eg. non-Metazoa) based on the sum of bitscores across all hits; the category with the highest bitscore wins. Each individual hit is also designated *ingroup* or *outgroup*, based on its taxonomy. To be classified as a putative HGT candidate, a given protein must (a) be designated as *outgroup* and (b) have support for this designation from a given proportion of the rest of the hits (default is set to >=90%). This approach may be less prone to errors associated with contamination issues in genome data, as it does not rely on evidence of providence from only the top hit. See [Koutsovoulos et al. 2016](http://www.pnas.org/content/113/18/5053.abstract) for more details.

### Prerequisites

1. Perl libraries: `Getopt::Long`, `Term::ANSIColor`, `Sort::Naturally`, `Data::Dumper`, and `List::Util`.

2. **Download NCBI Taxonomy:** These scripts rely on the NCBI taxonomy databases to assign taxonomic info to proteins. In particular the files "nodes.dmp" and "names.dmp" are required, while "merged.dmp" is recommended.
   ```
   >> wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
   >> tar xzvf taxdump.tar.gz
   >> export TAXPATH=/path/to/taxdump/
   ```

   It's also possible to use the "nodesDB.txt" file from [Blobtools](https://blobtools.readme.io/docs/what-is-blobtools) software.

3. [See update below] **Taxify your BLAST/Diamond file:** Diamond is great for speed, but adding the taxid information to each hit requires an additional step. See [this Gist](https://gist.github.com/sujaikumar/9ad04e62449a2d7025b17144de67038b) by Sujai Kumar on how to set this up for the UniRef90 database.

   A typical Diamond script might then look like:

   ```
   ## run diamond
   >> diamond blastp --sensitive --index-chunks 1 -k 500 -e 1e-5 -p 100 -q $QUERY -d $DB -a ${QUERY}.vs.uniref90.k500.1e5

   ## taxify output
   >> diamond view -a ${PREFIX}.daa | perl -lne 'BEGIN{open UT, "</path/to//uniref90.taxlist" or die $!; while (<UT>) { $ut{$1}=$2 if /^(\S+)\t(\S+)$/ } } {print "$_\t$ut{$1}" if /^\S+\t(\S+)/ and exists $ut{$1}}' > ${QUERY}.vs.uniref90.k500.1e5.daa.taxid;
   ```

   The taxified file has the usual 12 columns but with an additional 13th column containing the taxid of the hit protein.

---
### UPDATE Feb 2018
The latest UniRef90 database now contains the taxids in the fasta headers (thanks [@evolgenomology](https://github.com/evolgenomology)!), which makes generating the required taxlists a bit easier (no need to download and parse the XML file). Unfortunately you'll still need to add the taxid post-hoc as there is no option in Diamond to add this to the output file.

1. Download the UniRef90 fasta database (~15 Gb):
   ```
   wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
   ```
2. Generate taxlist somehow (can be as simple as):
   ```
   >> perl -lane 'if(/^>(\w+)\s.+TaxID\=(\d+)/){print "$1 $2"}' <(zcat uniref90.fasta.gz) | gzip > uniref90.fasta.taxlist.gz
   ```
3. Run Diamond as above:
   ```
   >> diamond blastp --sensitive --index-chunks 1 -k 500 -e 1e-5 -p 100 -q $QUERY -d /path/to/uniref90.fasta.dmnd -a ${QUERY}.vs.uniref90.k500.1e5
   ```
4. Taxify output:
   ```
   >> cat <(zcat /path/to/uniref90.fasta.taxlist.gz) <(diamond view -a ${QUERY}.vs.uniref90.k500.1e5.daa) \
   | perl -lane '
    if(@F==2){
      $tax{$F[0]}=$F[1];
    }else{
      if(exists($tax{$F[1]})){
        print join("\t",@F,$tax{$F[1]});
      } else {
        print join("\t",@F,"NA");
      }
    }
   ' | gzip > ${QUERY}.vs.uniref90.k500.1e5.daa.taxid.gz
   ```
   This will insert the taxid into a 13th column as above.
---

### Options

Type `-h` to see help and options.

```
  -i|--in              [FILE]   : taxified diamond/blast results file (accepts gzipped)
  -l|--list            [FILE]   : list of diamond/blast files to analyse [-i or -l required]
  -p|--path            [PATH]   : path to dir/ containing tax files
  -o|--nodes           [FILE]   : path to nodes.dmp
  -a|--names           [FILE]   : path to names.dmp
  -m|--merged          [FILE]   : path to merged.dmp
  -n|--nodesDB         [FILE]   : nodesDB.txt file from blobtools
  -t|--taxid_ingroup   [INT]    : NCBI taxid to define 'ingroup' [default=33208 (Metazoa)]
  -k|--taxid_skip      [INT]    : NCBI taxid to skip; hits to this taxid will not be considered
  -s|--CHS_threshold   [FLOAT]  : Consesus Hits Support threshold [default>=90\%]
  -u|--hU_threshold    [INT]    : hU threshold (default>=30)
  -e|--evalue_column   [INT]    : define evalue column [default=11]
  -b|--bitscore_column [INT]    : define bitscore column [default=12]
  -c|--taxid_column    [INT]    : define taxid column [default=13]
  -d|--delimiter       [STRING] : infile delimiter (diamond (\"\\s+\") or blast (\"\\t\")) [default=diamond]
  -x|--prefix          [FILE]   : filename prefix for outfile [default=INFILE]
  -v|--verbose                  : say more things
  -h|--help                     : this help message
```

### Usage

1. Assuming you have a taxified Diamond file with the hit taxids in the 13th column, and have downloaded the NCBI taxdump files to a location specified by $TAXPATH, running the script with default thresholds is as simple as:
  ```
  >> diamond_to_HGT_candidates.pl -i diamond_results.daa.taxid -p $TAXPATH
  ```
2. To exclude hits to organisms that are closely related to your focal organism, use -k option (recommended); eg. if you're looking at HGT in a nematode genome:
  ```
  >> diamond_to_HGT_candidates.pl -i diamond_results.daa.taxid -p $TAXPATH -k 6231
  ```
### Input

1. `-i` flag can read single file or string of filenames whitespace delimited: `-i "file1 file2 file3"`
2. `-l` flag can read a list of filenames in a textfile, one per line. In this case the `-k` flag can be specified on the command line (applied to all files) or as a second column in the list of filenames: `file1 6231` etc.

### Outputs

The program outputs 3 files, suffixed with the tags:

1. **HGT_results:** hU, AI and CHS scores for all query proteins.
2. **HGT_candidates:** All queries which show evidence from *both* hU (or AI) *and* CHS above the specifed thresholds (defaults are >=30 for hU, >=90% for CHS).
3. **HGT_warnings:** Any non-fatal warnings detected during the run. Worth checking, but most warnings can probably safely be ignored.

## HGT_candidates_to_fasta.pl

### Synopsis

This script takes the results from the "HGT_candidates" output file and outputs a fasta file containing the sequence data for each HGT candidate + its hits from UniRef90, each of which has been annotated with *ingroup* or *outgroup* and some taxonomy information up to Phylum level. These fasta files can then be aligned for tree building.

### Prerequisites

1. The script uses `blastdbcmd` from the BLAST suite for speedy sequence retrieval. This requires you generate a BLAST database of the UniRef90 fasta sequences **with the -parse_seqids option specifed**, otherwise it can't retrieve the sequence data itself using the fasta header:
   ```
   >> makeblastdb -in uniref90.fasta -dbtype prot -parse_seqids
   ```

   Also requires that `blastdbcmd` is in your `$PATH`.

### Options

Type `-h` to see the options:

```
-i|--in              [FILE]   : taxified diamond/BLAST results file [required]
-c|--candidates      [FILE]   : HGT_candidates.txt file [required]
-u|--uniref90        [FILE]   : diamond/BLAST database fasta file, e.g. UniRef90.fasta [required]
-f|--fasta           [FILE]   : fasta file of query proteins [required]
-p|--path            [STRING] : path to dir/ containing tax files [required]
-t|--taxid_threshold [INT]    : NCBI taxid to recurse up to; i.e., threshold taxid to define 'ingroup' [default = 33208 (Metazoa)]
-k|--taxid_skip      [INT]    : NCBI taxid to skip; hits to this taxid will not be considered in any calculations of support
-l|--limit                    : maximum number of ingroup / outgroup sequences to fetch (if available) [default = 15]
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

2. Construct phylogenies using iqtree:

   ```
   mkdir processed_files
   COUNT=1;

   ## iqtree commands
   for file in *mafft;
      do echo $COUNT: $file;
      iqtree-omp -s $file -st AA -nt 16 -quiet -bb 1000 -m TESTNEW -msub nuclear
      mv $file processed_files/;
      COUNT=$[COUNT+1];
   done

   mkdir iqtree treefiles
   mv *treefile treefiles/
   mv *bionj *gz *contree *iqtree *log *mldist *model *nex iqtree/
   ```

## get_locations_of_HGT_candidates.pl

### Synopsis

Takes the .HGT_results file from above, a GFF, and a list of protein names (can be parsed from fasta headers), and looks at the genomic location of all genes in light of their HGT candidacy. Physical "linkage" (here defined simply as genes encoded on the same scaffold) of HGT candidates (hU >= 30 and CHS >= 90% by default) with 'good' metazoan genes (hU <= 0 by default) is reported, as is the presence of introns in HGT candidates (but see note below). Scaffolds encoding a high proportion (default >= 95%) of HGT candidates are also flagged; these are likely derived from contaminant DNA and probably warrant further investigation and/or exclusion.

#### Notes

1. The presence of introns may not be a good measure of HGT support, see [Koutsovoulos et al](http://www.pnas.org/content/113/18/5053) and this [Gist](https://gist.github.com/GDKO/bc507bc9b620e6006a44). Number of introns is inferred by counting the number of CDS with a given protein ID from the input GFF - this works in most cases but may break if this correspondence does not fit your GFF.
2. A gene can have no information (annotated `NA`) if (a) it has no hit and thus no entry in the HGT_results file ("no-hitters"), or (b) it hits only to taxa that fall under the skipped category specified by `--taxid_skip` during the `diamond_to_HGT_candidates.pl` analysis ("hit-to-skippers"). Such genes are not considered in the current analysis - eg., if a HGT candidate is linked _only_ to hit-to-skippers it will still be counted as unlinked. This seems fair, since it is necessary to discount any association between each protein and its taxonomic annotation, and most genuine metazoan genes should have good hits to homologs across the Metazoa.
3. The sequence names specified in `--names` must correspond to the protein names in HGT_results and the GFF. An optional regex can be applied to fasta headers, this will be fed into a string substitution like: `s/$REGEX//ig`.

### Options

Type `-h` to see the options:
```
-i|--in     [FILE] : HGT_results.txt file [required]
-g|--gff    [FILE] : GFF file [required]
-n|--names  [FILE] : names of proteins in GFF file, can be fasta of proteins used
-r|--regex  [STR]  : optional regex to apply to seq headers if -n is a fasta file
-u|--outgrp [INT]  : threshold hU score for determining 'good' OUTGROUP (HGT) genes [default>=30]
-U|--ingrp  [INT]  : threshold hU score for determining 'good' INGROUP genes [default<=0]
-c|--CHS    [INT]  : threshold CHS score for determining 'good' OUTGROUP (HGT) genes [default>=90\%]
-y|--heavy  [INT]  : threshold for determining 'HGT heavy' scaffolds [default>=95\%]
-b|--bed           : also write bed file for 'good' HGT genes
-h|--help          : prints this help message
```

### Outputs

1. HGT_locations: reports gene-by-gene HGT results in pseudo-BED format. Gene positions on scaffolds inherited from gene coordinates in input GFF.
2. HGT_locations.summary: reports per-scaffold summary of number of HGT candidates
3. HGT_locations.heavy: reports scaffolds with a high proportion (dictated by `--heavy`) of HGT candidates
4. HGT_locations.bed (optional): BED format file of HGT candidates. Useful eg. for intersection with RNASeq mapping data.

## analyse_trees.R

Rscript to assess topological support for HGT candidates. For each tree, tests for monophyly with non-metazoan Eukayotes, Fungi, Plants, Bacteria and Archaea. Requires the {ape} R package.

Run with `-h` for options:
```
Usage: /home/reuben/software/tools/hgt/analyse_trees.R [options]

Options:
	-t CHARACTER, --trees=CHARACTER
		File containing multiple trees, one per line, newick format [default=NULL]
	-p CHARACTER, --path=CHARACTER
		Path to tree files [default=./]
	-a CHARACTER, --pattern=CHARACTER
		Pattern match to glob proper treefiles in --path [default=*]
	-q CHARACTER, --query=CHARACTER
		Value used to identify query sequence [default=NULL]
	-d CHARACTER, --delim=CHARACTER
		Character to delimit IN and OUT in sequence names [default=_]
	-o CHARACTER, --out=CHARACTER
		Tab delim outfile [default=results.tab]
	-f CHARACTER, --pdf=CHARACTER
		PDF trees outfile [default=results.pdf]
	-h, --help
		Show this help message and exit
```

### Inputs

Path to dir of (newick) trees.

### Outputs

Program calculates:
1. Counts for monophyletic support for query with numerous groups;
2. Assessment of topological evidence for HGT;
3. PDF of trees colour-coded by ingroup / outgroup categorisation;

Typical output might look like:
```
Mean number of taxa per tree:
  Total taxa: 12.87234
  Ingroup taxa: 0.893617
  Outgroup taxa: 10.97872
```
Showing the number of taxa designated as 'ingroup' (eg Metazoa) versus 'outgroup' (eg non-Metazoa) per analysed tree.

```
Number of trees with:
  Ingroup taxa only: 0
  Outgroup taxa only: 36
  Both: 11
```
Showing number of trees with _only_ ingroup or _only_ outgroup taxa present, or with _both_ ingroup and outgroup taxa present.

```
Number of 'Both' trees where Query can be monophyletic with:
  Ingroup taxa only: 2
  Outgroup taxa only: 2
  Both: 5
Average bootstrap support for monophyly of Query with:
  Ingroup: 82.71429
  Outgroup: 79.28571
```
Tests for monophyly of Query with in/outgroup explicitly ask the question: is it possible to draw a root somewhere on the tree that results in a monophyletic cluster containing the Query + (all of the) in/outgroup taxa? Note that Monophyly of Query+ingroup and monophyly of Query+outgroup are not mutually exclusive (counted as 'Both').

```
Assessment of HGT support catagories:
  1. Monophyly with outgroup taxa is possible (incl. clusters with only outgroup sequences): 43
  2. Monophyly with outgroup taxa is possible (excl. clusters with only outgroup sequences): 7
  3. Monophyly with outgroup taxa is possible, with ingroup impossible, at least 3 each in/outgroup taxa: 0
  4. Monophyly with outgroup taxa is possible (>=70% bootstrap support), with ingroup impossible, at least 3 each in/outgroup taxa: 0
```
1. It is possible root the tree such that Query is monophyletic with outgroup taxa, includes those trees with _only_ outgroup taxa;
2. It is possible root the tree such that Query is monophyletic with outgroup taxa, **excludes** those trees with _only_ outgroup taxa (ie, requires at least 1 ingroup + 1 outgroup + Query in tree);
3. More stringent; requires at least 3 members from ingroup and outgroup present in tree, and excludes trees where monophyly with ingroup is also possible;
4. More stringent again; as above but requires high bootstrap value.

```
Number of trees where Query+GROUP are monophyletic:
  Metazoa: 7 (0)
  Fungi: 13 (0)
  plants: 5 (0)
  other eukaryotes: 6 (0)
  Bacteria: 16 (0)
  Archaea: 0 (0)
```
Note that these counts require monophyly of Query + _all_ members of `GROUP`; if `GROUP` is paraphyletic, won't be counted.
