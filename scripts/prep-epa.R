#!/usr/bin/env Rscript

# report
writeLines("\nPreparing sequences for EPA ...\n")

# load libs and funs
source(here::here("scripts/funs-libs.R"))

# import reduced reference library csv
reflib.red <- suppressMessages(suppressWarnings(read_csv(here("temp/taxonomic-assignment/custom-reference-library-reduced.csv"),guess_max=99999,col_types=cols())))

# import reduced reference library fasta
reflib.red.fas <- read.FASTA(file=here("temp/taxonomic-assignment/custom-reference-library.fasta"))

# load up the queries
fish.queries <- read.FASTA(file=here("temp/taxonomic-assignment/asvs-fish.fasta"))

# set names as queries
names(fish.queries) <- paste0(names(fish.queries),".queries")

# concatenate with the reference and align
combined.fas <- c(fish.queries,reflib.red.fas)

# remove Ns
combined.fas <- rm_ns(bin=combined.fas)

# write out
write.FASTA(combined.fas,file=here("temp/taxonomic-assignment/epa/epa.input.fasta"))

# generate taxonomy file
reflib.red %>% mutate(sciNameValid=str_replace_all(sciNameValid," ", "_")) %>% 
    mutate(taxonomy=paste(phylum,class,order,family,genus,sciNameValid,sep=";")) %>% 
    select(dbid,taxonomy) %>% 
    write_tsv(file=here("temp/taxonomic-assignment/epa/references.taxonomy.tsv"),col_names=FALSE)

# report
writeLines("\n...\nSequences prepared for EPA\n")
