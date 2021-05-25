#!/usr/bin/env Rscript

# load libs and funs
source(here::here("scripts/funs-libs.R"))

# get args
option_list <- list( 
    make_option(c("-p","--primer"), type="character")
    )

# set args
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))


# set prefix
if(opt$primer=="tele02") {
    prefix <- "12s.taberlet.noprimers"
} else if(opt$primer=="elas02" | opt$primer=="mifish-u" | opt$primer=="mifish-u-mod") {
    prefix <- "12s.miya.noprimers"
} else {
    stop("Primers must be 'tele02', 'elas02', 'mifish-u', or 'mifish-u-mod'")
}


# import reference library
custom.refs <- suppressMessages(suppressWarnings(read_csv(here("assets/custom-reference-library.csv"),guess_max=99999,col_types=cols())))

# subset the marker from the custom reference library reflib - using 50% length cutoff
reflib.sub <- subset_by_marker(prefix=prefix,df=custom.refs,thresh=0.5)

# collapse haps
reflib.red <- reflib.sub %>% group_by(sciNameValid) %>% group_modify(~ hap_collapse_df(df=.x,lengthcol=paste0("lengthFrag.",prefix),nuccol=paste0("nucleotidesFrag.",prefix),cores=1)) %>% ungroup()

# load up the queries
fish.queries <- read.FASTA(file=here("temp/taxonomic-assignment/asvs-fish.fasta"))

# set names as queries
names(fish.queries) <- paste0(names(fish.queries),".queries")

# concatenate with the reference and align
combined.fas <- c(fish.queries,tab2fas(reflib.red, seqcol=paste0("nucleotidesFrag.",prefix),namecol="dbid"))

# remove Ns
combined.fas <- rm_ns(bin=combined.fas)

# write out
write.FASTA(combined.fas,file=here("temp/taxonomic-assignment/epa/epa.input.fasta"))

# generate taxonomy file
reflib.red %>% mutate(sciNameValid=str_replace_all(sciNameValid," ", "_")) %>% 
    mutate(taxonomy=paste(subphylum,class,order,family,genus,sciNameValid,sep=";")) %>% 
    select(dbid,taxonomy) %>% 
    write_tsv(file=here("temp/taxonomic-assignment/epa/references.taxonomy.tsv"),col_names=FALSE)

# report
writeLines("\n...\nSequences prepared for EPA\n")
