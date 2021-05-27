#!/usr/bin/env Rscript

# report
writeLines("\nProcessing BLAST results ...\n")

# load libs and funs
source(here::here("scripts/funs-libs.R"))

# import reference library
custom.refs <- suppressMessages(suppressWarnings(read_csv(here("assets/custom-reference-library.csv"),guess_max=99999,col_types=cols())))

# load blast result
local.db.blast <- suppressMessages(suppressWarnings(read_tsv(here("temp/taxonomic-assignment/fish-blast-result.tsv"))))

# choose "best" hit based on bitscore
# also add scinames
local.db.blast.sorted <- local.db.blast %>%
    group_by(asv) %>%
    arrange(desc(blastBitscore),.by_group=TRUE) %>%
    filter(blastBitscore==max(blastBitscore)) %>%
    mutate(blastSpeciesID=pull(custom.refs,sciNameValid)[match(blastDbid,pull(custom.refs,dbid))]) %>%
    arrange(blastSpeciesID,.by_group=TRUE) %>%
    mutate(blastSpeciesID=paste(unique(blastSpeciesID),collapse="; ")) %>%
    slice(1) %>% 
    ungroup()

# write out for later 
local.db.blast.sorted %>% write_csv(here("temp/taxonomic-assignment/fish-blast-result-sorted.csv"))

writeLines("\n...\nBLAST results processed\n")
