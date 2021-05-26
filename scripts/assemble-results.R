#!/usr/bin/env Rscript

# load libs and funs
source(here::here("scripts/funs-libs.R"))


# reload global blast, local blast, epa results
fish.otus <- suppressMessages(suppressWarnings(read_csv(file=here("temp/taxonomic-assignment/sintax-table.csv"))))
local.db.blast.sorted <- suppressMessages(suppressWarnings(read_csv(file=here("temp/taxonomic-assignment/fish-blast-result-sorted.csv"))))
epa.results.filtered <- suppressMessages(suppressWarnings(read_csv(file=here("temp/taxonomic-assignment/epa/epa-results-filtered.csv"))))

# read in contamination exclusions
contam.table <- suppressMessages(suppressWarnings(read_csv(here("assets/contaminants-exclude.csv"))))

# add dnas, get length and mode
dnas <- read.FASTA(file=here("temp/taxonomic-assignment/asvs-clean-cat-relabel-derep.fasta"))
dnas.table <- tibble(asv=names(dnas),nucleotides=mapply(paste,collapse="",as.character(dnas),USE.NAMES=FALSE)) %>% mutate(asvLength=str_length(nucleotides))

# combine
taxonomy.results <- purrr::reduce(list(fish.otus,local.db.blast.sorted,epa.results.filtered,dnas.table),left_join, by="asv") %>% rename(asvHash=asv)

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

# filter contams and make a label
contam.table %<>% filter(primer=="tele02") %>% mutate(label=paste(otuCode,bestId,sep="|"))

# add contam status
taxonomy.results %<>% mutate(contaminationID=purrr::map_chr(nucleotides,find_matches,referencedb=contam.table), isContaminated=if_else(is.na(contaminationID),FALSE,TRUE)) 

# load master asv table
asvs.all <- suppressMessages(suppressWarnings(read_csv(file=here("temp/taxonomic-assignment/asvs-all.csv"))))

# make id
asvs.all %<>% mutate(asvCode=paste(primer,lib,asv,sep="-")) %>% select(-asv) %>% rename(asvHash=md5,primerSet=primer,library=lib)

# get dirs for otu tables
libs.dirs <- list.dirs(path=here("temp/processing"),full.names=TRUE,recursive=FALSE)

# load up all the otu tables
tables.list <- suppressMessages(suppressWarnings(lapply(here(libs.dirs,"results/asv-table.tsv"),read_tsv)))

# pivot long and merge otu tables
asv.tables.df <- purrr::map_dfr(tables.list,pivot_longer,cols=-asv,names_to="sampleHash",values_to="nreads") %>% filter(nreads > 0) %>% rename(asvCode=asv)

# load up table of events
events.list <- suppressMessages(suppressWarnings(lapply(here(libs.dirs,"results/events-hashes.csv"),read_csv)))

# merge otu tables
events.df <- bind_rows(events.list) %>% rename(sampleHash=hashLabel)

# join the codes table and events
asv.tables.events <- asv.tables.df %>% left_join(events.df,by="sampleHash")

# join with the hash table 
asvs.by.sample <- asv.tables.events %>% left_join(asvs.all,by=c("asvCode","primerSet","library")) 

# write out the non-homologous reads
#asvs.by.sample %>% filter(is.na(asvHash)) %>% select(-asvHash) %>% write_csv(here("results/asvs-non-hologous"))

# add reads to the taxonomy results
reads.by.asv <- asvs.by.sample %>% filter(!is.na(asvHash)) %>% group_by(asvHash) %>% summarise(nreads=sum(nreads))
taxonomy.results %<>% left_join(reads.by.asv,by="asvHash")

# add taxonomy to the asvs.by.sample
taxonomy.results.reduced <- taxonomy.results %>% select(asvHash,isFish,assigned,assignedName,isContaminated)

# remove non-assigned and collapse asvs
fishes.by.sample <- asvs.by.sample %>% left_join(taxonomy.results.reduced,by="asvHash") %>% 
    filter(isFish==TRUE & assigned==TRUE & isContaminated==FALSE) %>%
    group_by(primerSet,library,eventID,sampleHash,assignedName) %>% 
    summarise(nreads=sum(nreads),.groups="drop") %>%
    arrange(primerSet,library,eventID,sampleHash,desc(nreads),assignedName) %>%
    rename(nReads=nreads)

# write out
fishes.by.sample %>% write_csv(file=here("results/fishes-by-sample.csv"))

# report
writeLines("\n...\nOutput written to 'results/fishes-by-sample.csv'\n")

## NEED TO do
# controls
# taxonomic assignment
# add non-homol reads to assignments table?
# suppress messages



# write
# taxonomy.results %>% arrange(desc(isFish),assigned,kingdom,phylum,class,order,family,genus) %>% write_csv(file=here("temp/taxonomic-assignment/all.csv"))


# check
#taxonomy.results %>% select(asv,sintaxSpeciesID,assigned,assignedName,contaminationID) %>% arrange(desc(assigned),asv) %>% print(n=Inf)

# take a look
#taxonomy.results %>% 
#    select(-kingdom,-phylum,-class,-order,-family,-genus,-propTotal,-dnas,-blastDbid,-blastEvalue,-blastScore,-blastBitscore,-lib) %>%
#    arrange(desc(assigned),desc(epaAssign),desc(blastAssign)) %>% 
#    write_csv(file=paste0(proj.path,"/results/sanity-check.csv"))

# write out
#taxonomy.results %>% write_csv(file=paste0(proj.path,"/results/taxonomic-assignment.csv"))
