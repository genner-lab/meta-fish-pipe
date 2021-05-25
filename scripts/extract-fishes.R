#!/usr/bin/env Rscript

# load libs and funs
source(here::here("scripts/funs-libs.R"))

# read in asv fasta
asvs.fasta <- read.FASTA(here("temp/taxonomic-assignment/asvs-clean-cat-relabel-derep.fasta"))

# read in taxonomy assignment
tax.ass.df <- suppressMessages(suppressWarnings(read_tsv(file=here("temp/taxonomic-assignment/sintax-output.tsv"),col_names=c("asv","idsProbs","strand","ids"),guess_max=999999)))

# filter the output
tax.ass.df %<>% 
    mutate(sintaxBS=str_replace_all(map(str_split(pull(tax.ass.df,idsProbs),"\\("), last),"\\)","")) %>%
    mutate(idsProbsTemp=str_replace_all(idsProbs,"[a-z]:",""),idsProbsTemp=str_replace_all(idsProbsTemp,"\\([0-9].[0-9][0-9]\\)","")) %>%
    separate(idsProbsTemp,into=c("kingdom","phylum","class","order","family","genus","sintaxSpeciesID"),sep=",") %>%
    mutate(sintaxSpeciesID=str_replace_all(sintaxSpeciesID,"_", " ")) %>%
    mutate(isFish=if_else(class=="Sarcopterygii" | class=="Cephalaspidomorphi" | class=="Elasmobranchii" | class=="Myxini" | class=="Holocephali" | class=="Actinopterygii" | class=="Chondrichthyes", TRUE, FALSE)) %>%
    select(asv,isFish,kingdom,phylum,class,order,family,genus,sintaxSpeciesID,sintaxBS)

# write out formatted sintax results
tax.ass.df %>% write_csv(file="temp/taxonomic-assignment/sintax-table.csv")

# extract the ASVs from the fasta
asvs.fish <- tax.ass.df %>% filter(isFish==TRUE) %>% distinct(asv) %>% pull()

# subset fasta
asvs.fasta.fish <- asvs.fasta[which(names(asvs.fasta) %in% asvs.fish)]

#write out
write.FASTA(asvs.fasta.fish,file=here("temp/taxonomic-assignment/asvs-fish.fasta"))

writeLines("\n...\nFish sequences extracted\n")
