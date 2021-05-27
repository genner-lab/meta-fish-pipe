#!/usr/bin/env Rscript

# report
writeLines("\nProcessing EPA results ...\n")

# load libs and funs
source(here::here("scripts/funs-libs.R"))

# load up epa results
epa.results <- suppressMessages(suppressWarnings(read_tsv(file=here("temp/taxonomic-assignment/epa/per_query.tsv"))))
epa.results %<>% rename(asv=name)
# LWR: likelihood weight that was assigned to this exact taxonomic path
# fract: LWR divided by the global total likelihood weight
# aLWR: accumulated likelihood weights that were assigned either to this taxonomic path or any taxonomic path below this
# afract: aLWR divided by the global total likelihood weight
# taxopath: the taxonomic path

# best the best id
epa.results.best <- epa.results %>% 
    group_by(asv) %>%
    arrange(desc(LWR),.by_group=TRUE) %>% 
    slice(1) %>% 
    mutate(epaID=sapply(str_split(taxopath,";"),last), epaID=str_replace_all(epaID,"_"," ")) %>% 
    ungroup()

# get species IDs
epa.results.species <- epa.results %>% 
    group_by(asv) %>% 
    filter(grepl("_",taxopath)) %>% 
    arrange(taxopath,.by_group=TRUE) %>% 
    mutate(epaBestSppID=sapply(str_split(taxopath,";"),last), epaBestSppID=str_replace_all(epaBestSppID,"_"," ")) %>%  
    mutate(epaAllSpp=paste(unique(epaBestSppID),collapse="; ")) %>% 
    arrange(desc(LWR),.by_group=TRUE) %>% 
    slice(1) %>% 
    ungroup() %>%
    select(asv,epaBestSppID,epaAllSpp)

# join and tidy
epa.results.filtered <- left_join(epa.results.best,epa.results.species) %>% 
    select(asv,LWR,epaID,epaBestSppID,epaAllSpp) %>% 
    rename(epaIdLWR=LWR)

# write out 
epa.results.filtered %>% write_csv(file=here("temp/taxonomic-assignment/epa/epa-results-filtered.csv"))

# report
writeLines("\n...\nEPA results processed\n")
