#!/usr/bin/env Rscript

# report
writeLines("\nAssembling results ...\n")

# load libs and funs
source(here::here("scripts/funs-libs.R"))

# get args
option_list <- list( 
    make_option(c("-c","--contam"), type="character")
    )
# set args
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))


############## LOAD DATA AND JOIN ##############
############## LOAD DATA AND JOIN ##############

# reload global blast, local blast, epa results
fish.otus <- suppressMessages(suppressWarnings(read_csv(file=here("temp/taxonomic-assignment/sintax-table.csv"))))
local.db.blast.sorted <- suppressMessages(suppressWarnings(read_csv(file=here("temp/taxonomic-assignment/fish-blast-result-sorted.csv"))))
epa.results.filtered <- suppressMessages(suppressWarnings(read_csv(file=here("temp/taxonomic-assignment/epa/epa-results-filtered.csv"))))

# read in contamination exclusions
contam.table <- suppressMessages(suppressWarnings(read_csv(here(opt$contam))))

# add dnas, get length and mode
dnas <- read.FASTA(file=here("temp/taxonomic-assignment/asvs-clean-cat-relabel-derep.fasta"))
dnas.table <- tibble(asv=names(dnas),nucleotides=mapply(paste,collapse="",as.character(dnas),USE.NAMES=FALSE)) %>% mutate(asvLength=str_length(nucleotides))

# combine
taxonomy.results <- purrr::reduce(list(fish.otus,local.db.blast.sorted,epa.results.filtered,dnas.table),left_join,by="asv") %>% rename(asvHash=asv)


############## ASSIGN TAXONOMY ##############
############## ASSIGN TAXONOMY ##############

# set a 90% of modal length threshold
minLen <- ceiling(mode_avg(pull(taxonomy.results,asvLength))*0.9)

# two step process - either (a) or (b) must be satisfied
# (a) epaAssign = highest likelihood EPA id and blast id the same, EPA likelihood > 0.9, and pident > 0.9
# (b) blastAssign = the best species-level EPA id (not necesarily best overall id) and blast id the same, and pident >0.97 and lengthLocal >0.9 of modal ASV length 
taxonomy.results %<>% 
    mutate(epaAssign=if_else(epaID==blastSpeciesID & epaIdLWR>=0.9 & blastPident>=90,TRUE,FALSE)) %>% 
    mutate(blastAssign=if_else(str_detect(blastSpeciesID,epaBestSppID) & blastLength>=minLen & blastPident>=97,TRUE,FALSE)) %>% 
    mutate(assigned=if_else(epaAssign==TRUE | blastAssign==TRUE,TRUE,FALSE)) %>% 
    mutate(assignedName=if_else(assigned==TRUE,if_else(epaAssign==TRUE,epaID,epaBestSppID),"NA")) %>%
    mutate(assigned=if_else(isFish==FALSE,FALSE,assigned))

# to assign other verts with sintax results
taxonomy.results %<>%
    dplyr::mutate(assigned=dplyr::if_else(isFish==FALSE & assigned==FALSE & sintaxBS>=0.95 & nreads>100,TRUE,assigned)) %>%
    dplyr::mutate(assignedName=dplyr::if_else(isFish==FALSE & assigned==TRUE & sintaxBS>=0.95 & nreads>100,sintaxSpeciesID,assignedName)) %>%
    dplyr::mutate(isFish=dplyr::if_else(isFish==FALSE & assigned==TRUE & sintaxBS>=0.95 & nreads>100,TRUE,isFish))
 

############## ADD CONTAMINANTS ##############
############## ADD CONTAMINANTS ##############

# make a label
contam.table %<>% mutate(label=paste(otuCode,bestId,sep="|"))

# add contam status
taxonomy.results %<>% mutate(contaminationID=purrr::map_chr(nucleotides,find_matches,referencedb=contam.table), isContaminated=if_else(is.na(contaminationID),FALSE,TRUE)) 


############## COMBINE ASVs AND EVENTS ##############
############## COMBINE ASVs AND EVENTS ##############

# load master asv table
asvs.all <- suppressMessages(suppressWarnings(read_csv(file=here("temp/taxonomic-assignment/asvs-all.csv"))))

# make id
asvs.all %<>% mutate(asvCode=paste(primer,lib,asv,sep="|")) %>% 
    select(-asv) %>% 
    rename(asvHash=md5,primerSet=primer,library=lib)

# get dirs for otu tables
libs.dirs <- list.dirs(path=here("temp/processing"),full.names=TRUE,recursive=FALSE)

# load up all the otu tables
tables.list <- suppressMessages(suppressWarnings(lapply(here(libs.dirs,"results/asv-table.tsv"),read_tsv)))

# pivot long and merge otu tables
asv.tables.df <- purrr::map_dfr(tables.list,pivot_longer,cols=-asv,names_to="sampleHash",values_to="nreads") %>% 
    filter(nreads > 0) %>% 
    rename(asvCode=asv)

# load up table of events
events.list <- suppressMessages(suppressWarnings(lapply(here(libs.dirs,"results/events-hashes.csv"),read_csv)))

# merge otu tables
events.df <- bind_rows(events.list) %>% rename(sampleHash=hashLabel)

# join the codes table and events
asv.tables.events <- asv.tables.df %>% left_join(events.df,by="sampleHash")

# join with the hash table 
asvs.by.sample <- asv.tables.events %>% left_join(asvs.all,by=c("asvCode","primerSet","library")) 

# write out the non-homologous reads
#asvs.by.sample %>% filter(is.na(asvHash)) %>% select(-asvHash) %>% write_csv(here("results/asvs-non-homologous.csv"))

# add reads to the taxonomy results
reads.by.asv <- asvs.by.sample %>% 
    filter(!is.na(asvHash)) %>% 
    group_by(asvHash) %>% 
    summarise(nreads=sum(nreads),.groups="drop")

# join
taxonomy.results %<>% left_join(reads.by.asv,by="asvHash")

# add taxonomy to the asvs.by.sample
taxonomy.results.reduced <- taxonomy.results %>% select(asvHash,isFish,assigned,assignedName,isContaminated)

# remove non-assigned and collapse asvs
fishes.by.sample <- asvs.by.sample %>% left_join(taxonomy.results.reduced,by="asvHash") %>% 
    filter(isFish==TRUE & assigned==TRUE & isContaminated==FALSE) %>%
    group_by(primerSet,library,eventID,sampleHash,replicateFilter,assignedName) %>% 
    summarise(nreads=sum(nreads),.groups="drop") %>%
    arrange(primerSet,library,eventID,sampleHash,desc(nreads),assignedName) %>%
    rename(nReads=nreads)


############## FILTER AND WRITE OUT ##############
############## FILTER AND WRITE OUT ##############

# write out fishes after removing blanks
fishes.by.sample %>% 
    filter(!grepl("Blank",replicateFilter)) %>% 
    select(-replicateFilter) %>%
    write_csv(file=here("results/fishes-by-sample.csv"))

# write out taxonomic identifications
taxonomy.results %>% select(asvHash,nreads,isFish,assigned,assignedName,
    sintaxSpeciesID,sintaxBS,kingdom,phylum,class,order,family,genus,
    blastAssign,blastSpeciesID,blastDbid,blastEvalue,blastLength,blastPident,blastNident,blastScore,blastBitscore,
    epaAssign,epaIdLWR,epaID,epaBestSppID,epaAllSpp,
    isContaminated,contaminationID,asvLength,nucleotides) %>%
    arrange(desc(isFish),desc(assigned),desc(nreads)) %>% 
    write_csv(file=here("results/taxonomic-assignments.csv"))


############## NEGATIVE CONTROLS ##############
############## NEGATIVE CONTROLS ##############

# isolate blanks with assigned fish reads
blanks.with.fishes <- fishes.by.sample %>% filter(grepl("Blank",replicateFilter))

# remove samples with reads from all events
events.df.blanks <- events.df %>% 
    filter(grepl("Blank",replicateFilter)) %>% 
    mutate(assignedName=NA,nReads=0) %>% 
    filter(!sampleHash %in% pull(blanks.with.fishes,sampleHash))

# combine and sort
blanks.all <- bind_rows(blanks.with.fishes,events.df.blanks) %>% 
    select(primerSet,library,replicateFilter,eventID,sampleHash,assignedName,nReads) %>%
    arrange(primerSet,library,replicateFilter,eventID,sampleHash,desc(nReads)) %>%
    rename(blankType=replicateFilter)

# write out
blanks.all %>% write_csv(file=here("results/controls-summary.csv"))


############## SUMMARISE STATS ##############
############## SUMMARISE STATS ##############

# load stats files
stats.list <- suppressMessages(suppressWarnings(lapply(here(libs.dirs,"logs/stats.csv"),read_csv)))

# merge
stats.merged <- bind_rows(stats.list)

# add assigned
assigned.sum <- fishes.by.sample %>% 
    mutate(library=paste(primerSet,library,sep="-")) %>% 
    group_by(library) %>% 
    summarise(nreads=sum(nReads),.groups="drop") %>% 
    mutate(stat="assigned")

# join and sort and write out
bind_rows(stats.merged,assigned.sum) %>% 
    arrange(library,desc(nreads)) %>% 
    write_csv(file=here("results/stats-summary.csv"))


############## REPORT ##############
############## REPORT ##############

# report
writeLines("\n...\nAssigned fish reads by sample written to: 'results/fishes-by-sample.csv'\n")
writeLines("\n...\nTaxonomic assignments by ASV written to: 'results/taxonomic-assignments.csv'\n")
writeLines("\n...\nSummary of negative controls written to: 'results/controls-summary.csv'\n")
writeLines("\n...\nSummary of library stats written to: 'results/stats-summary.csv'\n")
