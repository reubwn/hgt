# hgt
scripts for the analysis of HGT in genome sequence data

## diamond_to_HGT_candidates.pl

Goal is to take a taxified Diamond or BLAST file, and for each hit recurse up the tax tree until that hit can be categorised into **INGROUP** versus **OUTGROUP** (e.g., Metazoan vs non-Metazoan etc.).

1. **Get Query Category:** For each query, calculate the bitscoresum for ingroup vs outgroup across **all hits**; the category with the highest bitscoresum is the \"winner\"
2. **Get Support:** Assess support for the winning query taxid from secondary hits; winning taxid is well-supported if the rest of the hits agree with the INGROUP/OUTGROUP categorisation above --support_threshold (default = 90%)
3. **Calculate AI:** Also calculate Alien Index based on best e-values to INGROUP vs OUTGROUP
4. **Print:** General results printed to HGT_results, candidate HGT genes printed to HGT_candidates (candidates printed if evidence of HGT from _either_ bitscoresum _or_ AI)

**Alien Index** is ```log((Best E-value for Metazoa) + 1e-200) - log((Best E-value for NonMetazoa) + 1e-200)``` (see Gladyshev et al., http://science.sciencemag.org/content/suppl/2008/05/29/320.5880.1210.DC1/Gladyshev.SOM.pdf)

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
