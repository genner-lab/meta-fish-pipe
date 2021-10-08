#!/usr/bin/env Rscript

# install a bioc package with renv
# renv::install(packages="bioc::dada2")

# load libs
suppressMessages(library("here"))
suppressMessages(library("parallel"))
suppressMessages(library("tidyverse"))
suppressMessages(library("magrittr"))
suppressMessages(library("openssl"))
suppressMessages(library("ape"))
suppressMessages(library("optparse"))
suppressMessages(library("dada2"))

# function to reverse complement DNA
source("https://raw.githubusercontent.com/legalLab/protocols-scripts/master/scripts/revcomp_dna.R")


# function for making fasta files from tables
source("https://raw.githubusercontent.com/legalLab/protocols-scripts/master/scripts/tab2fas.R")


# function to subset the reference dataframe by marker and filter by sequence length
subset_by_marker <- function(prefix,df,thresh){
    df1 <- df %>% filter(!is.na(!!as.name(paste0("nucleotidesFrag.",prefix))))
    df1 %<>% filter(!!as.name(paste0("lengthFrag.",prefix)) >= (median(!!as.name(paste0("lengthFrag.",prefix)))*thresh))
    return(df1)
}


# collapses haplotypes (from a dataframe format to a dataframe format)
# need to specify columns that contain sequence lengths, and nucleotides
# hap_collapse_df(df=mydataframe,lengthcol="lengthFrag",nuccol="nucleotidesFrag")
# add a number of each haplotype
hap_collapse_df <- function(df,lengthcol,nuccol,cores){
    odf <- df[order(df[[lengthcol]],decreasing=TRUE),]
    reps <- mcmapply(FUN=function(x) which(str_detect(string=odf[[nuccol]], pattern=x) == TRUE)[1], odf[[nuccol]], SIMPLIFY=TRUE, USE.NAMES=FALSE, mc.cores=cores)
    ind <- unique(reps)
    dat <- odf[ind,]
    dat[["nHaps"]] <- as.numeric(table(reps))
    return(dat)
}


# function to get retrieve species names of sequences with an identical haplotype as your query 
# works on a dataframe
get_sames <- function(df,ids,nucs,sppVec,query){
    per.ind <- df[[sppVec]][str_detect(df[[nucs]], query)]
    return(per.ind)
}


# mode function
mode_avg <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}


# removes Ns from a DNAbin list
rm_ns <-function(bin){
    bin.char <- as.character(bin)
    bin.spread <- sapply(bin.char, paste, collapse="")
    bin.rep <- sapply(bin.spread, str_replace_all, "[^actg]", "")
    bin.split <- strsplit(bin.rep, "")
    bin.bin <- as.DNAbin(bin.split)
    return(bin.bin)
}


# function to create short file paths
cpath <- function(sense,step,r){
    path <- paste0(proj.path,"/processed-reads/",sense,"/",step,"-",r)
    return(path)
}


# function to query for contamination in an otu table
find_matches <- function(query,referencedb){
    res <- pull(referencedb,label)[which(grepl(query,pull(referencedb,nucleotides)))]
    res <- res[1]
    if(length(res)==0){
    res <- NA} else {res}
    return(res)
}
