#!/usr/bin/Rscript

################################################################################
## hgt_tree #
#####################
## Script for testing HGT candidate phylogenies

################################################################################
## get options #
################

## load option parsing lib
suppressPackageStartupMessages(library(optparse))
option_list = list(
  make_option(c('-t','--trees'), type='character', default=NULL, help='File containing multiple trees, one per line, newick format [default=%default]', metavar='character'),
  make_option(c('-p','--path'), type='character', default="./", help='Path to tree files [default=%default]', metavar='character'),
  make_option(c('-a','--pattern'), type='character', default="*", help='Pattern match to glob proper treefiles in --path [default=%default]', metavar='character'),
  make_option(c('-q','--query'), type='character', default=NULL, help='Value used to identify query sequence [default=%default]', metavar='character'),
  make_option(c('-o','--out'), type='character', default="results.tab", help='Tab delim outfile [default=%default]', metavar='character'),
  make_option(c('-d','--pdf'), type='character', default="results.pdf", help='PDF trees outfile [default=%default]', metavar='character')
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(naturalsort))

################################################################################
## get trees #
##############

#if (exists(opt$trees)) {
#  tree.files <- read.trees(opt$trees)
#} else {
tree.files<-dir(path = opt$path, pattern = opt$pattern)
tree.files<-naturalsort(tree.files)
#}
#tree.files <- dir(path = "~/Google_Drive/Bdelloid/results/hgt/test/", pattern = "treefile$")
#opt$query<-"ARIC"
results.colnames<-c("filename","ntax.tot","ntax.in","ntax.out","only.in","only.out","both","in.mono","out.mono","in.q.mono","in.q.mono.sup","out.q.mono","out.q.mono.sup","euk.mono","euk.q.mono","euk.metazoa.q.mono","euk.fungi.q.mono","euk.plants.q.mono","euk.other.q.mono","bacteria.mono","bacteria.q.mono","archaea.mono","archaea.q.mono")
results<-data.frame(matrix(0,ncol=length(results.colnames),nrow=length(tree.files)),stringsAsFactors = F)
colnames(results)<-results.colnames

## open PDF 
pdf(file=opt$pdf,width=15,height=15)

for (i in (1:length(tree.files))) {
  cat("Treefile:", i, "\r")
  
  ## read tree
  tr <- read.tree(tree.files[i])
  results$filename[i]<-tree.files[i]
  results$ntax.tot[i]<-length(tr$tip.label)
  
  ## get ingroup, outgroup, query taxa
  taxa.in <- grep("_IN_", tr$tip.label)
  taxa.out <- grep("_OUT_", tr$tip.label)
  taxa.query <- grep(opt$query, tr$tip.label)
  
  ## how many innies and outies in tree?
  results$ntax.in[i]<-length(taxa.in)
  results$ntax.out[i]<-length(taxa.out)
  ## tip colours based on IN or OUT
  tipcols<-array(length(tr$tip.label))
  tipcols[taxa.in]<-1 ## black
  tipcols[taxa.out]<-"blue2" ## blue
  tipcols[taxa.query]<-"red2" ## red
  ## plot unrooted tree
  plot(tr,type="unrooted",lab4ut="axial",cex=0.75,show.node.label=T,font=1,tip.color=tipcols,main=tree.files[i])
  mtext(paste("ntax.total=",length(tr$tip.label),"; ntax.in=",length(taxa.in),"; ntax.out=",length(taxa.out),sep=""),col="grey30")
  
  ## populate with NA's to begin with
  results$in.q.mono.sup[i]<-NA
  results$out.q.mono.sup[i]<-NA
  
  if ( (length(taxa.in)>0)&&(length(taxa.out)==0) ){
    ## only innies in tree (should be rare/never)
    results$only.in[i]<-1
    results$in.mono[i]<-1
    results$in.q.mono[i]<-1
    
  } else if ((length(taxa.in)==0)&&(length(taxa.out)>0)) {
    ## only outies in tree (common)
    results$only.out[i]<-1
    results$out.mono[i]<-1
    results$out.q.mono[i]<-1
    
  } else if ((length(taxa.in)>0)&&(length(taxa.out)>0)) {
    ## both some innies and some outies in tree
    results$both[i]<-1
    
    ## are innies and/or outies monophyletic?
    tr.trim<-drop.tip(tr,taxa.query)
    results$in.mono[i]<-is.monophyletic(tr.trim,grep("_IN_",tr.trim$tip.label),reroot=T)
    results$out.mono[i]<-is.monophyletic(tr.trim,grep("_OUT_",tr.trim$tip.label),reroot=T)
    
    ## is query monophyletic with innies?
    results$in.q.mono[i]<-is.monophyletic(tr,c(taxa.query,taxa.in),reroot=T)
    ## reroot on any outie; get support value at in.query.mono node if TRUE
    if (is.monophyletic(tr,c(taxa.query,taxa.in),reroot=T)==T) {
      sr<-root(tr,taxa.out[1])
      ## sometimes there is no node label
      if (nchar(sr$node.label[getMRCA(sr,c(taxa.query,grep("_IN_",sr$tip.label)))-length(sr$tip.label)])>0){
        results$in.q.mono.sup[i]<-sr$node.label[getMRCA(sr,c(taxa.query,grep("_IN_",sr$tip.label)))-length(sr$tip.label)]
      }# else { results$in.q.mono.sup[i]<-NA }
      ## plot the monophyly
      is.monophyletic(tr,c(taxa.query,taxa.in),reroot=T,plot=T,edge.width=2,font=1,cex=0.75,show.node.label=F,main=tree.files[i])
      mtext("in.query.mono=TRUE",col="grey30")
      par(mfrow=c(1,1))
    }# else { results$in.q.mono.sup[i]<-NA }
    
    ## is query monophyletic with outies?
    results$out.q.mono[i]<-is.monophyletic(tr,c(taxa.query,taxa.out),reroot=T)
    ## reroot on any innie; get support value at out.query.mono node if TRUE
    if (is.monophyletic(tr,c(taxa.query,taxa.out),reroot=T)==T) {
      sr<-root(tr,taxa.in[1])
      ## sometimes there is no node label; so need to replace with 'NA'
      if (nchar(sr$node.label[getMRCA(sr,c(taxa.query,grep("_OUT_",sr$tip.label)))-length(sr$tip.label)])>0){
        results$out.q.mono.sup[i]<-sr$node.label[getMRCA(sr,c(taxa.query,grep("_OUT_",sr$tip.label)))-length(sr$tip.label)]
      }# else { results$out.q.mono.sup[i]<-NA }
      ## plot the monophyly
      is.monophyletic(tr,c(taxa.query,taxa.out),reroot=T,plot=T,edge.width=2,font=1,cex=0.75,show.node.label=F,main=tree.files[i])
      mtext("out.query.mono=TRUE",col="grey30")
      par(mfrow=c(1,1))
    }# else { results$out.q.mono.sup[i]<-NA }
  }
  
  ##############################################################################
  ##############################################################################
  ## Group-specific monophyly
  ## these all need to be evaluated individually as tree can contain a mixture of taxa, all of which need to be evaluated
  ## so if-else loops aren't appropriate
  
  ## EUKARYOTES
  if (length(grep("Eukaryota",tr$tip.label))>0){
    ## is euks monophyletic in the first place?
    tr.trim<-drop.tip(tr,taxa.query)
    results$euk.mono[i]<-is.monophyletic(tr.trim,c(grep("Eukaryota",tr.trim$tip.label)),reroot=T)
    ## is query monophyletic with euks?
    results$euk.q.mono[i]<-is.monophyletic(tr,c(taxa.query,grep("Eukaryota",tr$tip.label)),reroot=T)
    
    ## now look at kingdoms within Eukaryotes: Metazoa, Fungi, Plants, and 'other'
    if (length(grep("Eukaryota_Metazoa",tr$tip.label))>0){ results$euk.metazoa.q.mono[i]<-is.monophyletic(tr,c(taxa.query,grep("Eukaryota_Metazoa",tr$tip.label)),reroot=T) } 
    if (length(grep("Eukaryota_Fungi",tr$tip.label))>0){ results$euk.fungi.q.mono[i]<-is.monophyletic(tr,c(taxa.query,grep("Eukaryota_Fungi",tr$tip.label)),reroot=T) } 
    if (length(grep("Eukaryota_Viridiplantae",tr$tip.label))>0){ results$euk.plants.q.mono[i]<-is.monophyletic(tr,c(taxa.query,grep("Eukaryota_Viridiplantae",tr$tip.label)),reroot=T) } 
    if (length(grep("Eukaryota_undef",tr$tip.label))>0){ results$euk.other.q.mono[i]<-is.monophyletic(tr,c(taxa.query,grep("Eukaryota_undef",tr$tip.label)),reroot=T) }
  }
  
  ## BACTERIA
  if (length(grep("Bacteria",tr$tip.label))>0){
    ## are Bacteria monophyletic?
    tr.trim<-drop.tip(tr,taxa.query)
    results$bacteria.mono[i]<-is.monophyletic(tr.trim,c(grep("Bacteria",tr.trim$tip.label)),reroot=T)
    ## is query monophyletic with Bacteria?
    results$bacteria.q.mono[i]<-is.monophyletic(tr,c(taxa.query,grep("Bacteria",tr$tip.label)),reroot=T)
  } 
  
  ## ARCHAEA
  if (length(grep("Archaea",tr$tip.label))>0){
    ## are Archaea monophyletic?
    tr.trim<-drop.tip(tr,taxa.query)
    results$archaea.mono[i]<-is.monophyletic(tr.trim,c(grep("Archaea",tr.trim$tip.label)),reroot=T)
    ## is query monophyletic with Archaea?
    results$archaea.q.mono[i]<-is.monophyletic(tr,c(taxa.query,grep("Archaea",tr$tip.label)),reroot=T)
  }
}
dev.off()

## column totals / means where appropriate
cat("Number of trees analysed:",length(tree.files),"\n")
cat("Mean number taxa per tree:",mean(results$ntax.tot),"\n")
cat("Mean number ingroup taxa per tree:",mean(results$ntax.in),"\n")
cat("Mean number outgroup taxa per tree:",mean(results$ntax.out),"\n")
cat("Number of trees with only ingroup taxa present:",sum(results$only.in),"\n")
cat("Number of trees with only outgroup taxa present:",sum(results$only.out),"\n")
cat("Number of trees with both ingroup and outgroup taxa present:",sum(results$both),"\n")
cat("Number of trees where QUERY is monophyletic with INGROUP:",sum(results$in.q.mono),"\n")
cat("Number of trees where QUERY is monophyletic with OUTGROUP:",sum(results$out.q.mono),"\n")
cat("Average bootstrap support for monophyly of Q with INGROUP:",mean(as.numeric(results$in.q.mono.sup),na.rm=T),"\n")
cat("Average bootstrap support for monophyly of Q with OUTGROUP:",mean(as.numeric(results$out.q.mono.sup),na.rm=T),"\n")
cat("Number of trees where QUERY is monophyletic with METAZOA:",sum(results$euk.metazoa.q.mono),"\n")
cat("Number of trees where QUERY is monophyletic with FUNGI:",sum(results$euk.fungi.q.mono),"\n")
cat("Number of trees where QUERY is monophyletic with PLANTS:",sum(results$euk.plants.q.mono),"\n")
cat("Number of trees where QUERY is monophyletic with OTHER EUKARYOTES:",sum(results$euk.other.q.mono),"\n")
cat("Number of trees where QUERY is monophyletic with BACTERIA:",sum(results$bacteria.q.mono),"\n")
cat("Number of trees where QUERY is monophyletic with ARCHAEA:",sum(results$archaea.q.mono),"\n\n")
write.table(results, file=opt$out, quote=F, sep="\t", row.names=F)