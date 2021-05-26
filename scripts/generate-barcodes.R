#!/usr/bin/env Rscript

# load libs funs
source(here::here("scripts/funs-libs.R"))

# get args
option_list <- list( 
    make_option(c("-p","--primer"), type="character"),
    make_option(c("-l","--lib"), type="character"),
    make_option(c("-f","--lenfwd"), type="numeric"),
    make_option(c("-r","--lenrev"), type="numeric")
    )
# set args
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))

# for testing
#opt <- NULL
#opt$primer <- "tele02"
#opt$lib <- "lib3"
#opt$lenfwd <- 18
#opt$lenrev <- 20

# make paths
proj.path <- here("temp/processing",paste0(opt$primer,"-",opt$lib))

# confirm proj path 
writeLines(paste0("\nOutput directory set to:\n",proj.path))

# load up the data
plates <- suppressMessages(suppressWarnings(read_csv(file=here("assets/sequencing-master.csv"))))

# check samples are all present in events-master - missing ones are blanks/ntc
events.master <- suppressMessages(suppressWarnings(read_csv(file=here("assets/events-master.csv"))))
writeLines("\nSample IDs not present in events table (should all be blanks):")
plates %>% filter(!eventID %in% pull(events.master,eventID)) %>% distinct(eventID) %>% pull(eventID)

# filter by marker
plates %<>% filter(primerSet==opt$primer & library==opt$lib)

# create the barcode tags
plates %<>% mutate(barcodesFwd=str_replace_all(oligoFwd,"N",""), 
    barcodesFwd=str_trunc(barcodesFwd, width=10, side="right", ellipsis=""),
    barcodesRev=str_replace_all(oligoRev,"N",""),
    barcodesRev=str_trunc(barcodesRev, width=10, side="right", ellipsis=""),
    primerFwd=str_trunc(oligoFwd, width=opt$lenfwd, side="left", ellipsis=""),
    primerRev=str_trunc(oligoRev, width=opt$lenrev, side="left", ellipsis=""),
    labelFwd=str_trunc(str_replace_all(oligoFwd,"N",""), width=unique(str_length(str_replace_all(oligoFwd,"N",""))-opt$lenfwd), side="right", ellipsis=""),
    labelRev=str_trunc(str_replace_all(oligoRev,"N",""), width=unique(str_length(str_replace_all(oligoRev,"N",""))-opt$lenrev), side="right", ellipsis="")
    )

writeLines("\nChecking all oligos and labels match (should all be TRUE):")

# check oligos match (should be true)
table(pull(plates,labelFwd) == pull(plates,labelRev))
# cat the labels
plates %<>% mutate(senseLabel=paste(eventID,primerSet,library,replicateFilter,replicatePCR,labelFwd,sep="."), antisenseLabel=paste(eventID,primerSet,library,replicateFilter,replicatePCR,labelRev,sep="."))
# check labels match (should be true)
table(pull(plates,senseLabel) == pull(plates,antisenseLabel))

# create the hashes
plates %<>% mutate(senseLabelMD5=str_trunc(md5(senseLabel),width=12,side="right",ellipsis=""), antisenseLabelMD5=str_trunc(md5(antisenseLabel),width=12,side="right",ellipsis="")) 
# check md5s match (should be true)
table(pull(plates,senseLabelMD5) == pull(plates,antisenseLabelMD5))
# check for collisions (should be true)
length(pull(plates,senseLabelMD5)) == length(unique(pull(plates,senseLabelMD5)))

# rename as hashLabel
plates %<>% mutate(hashLabel=if_else(senseLabelMD5 == antisenseLabelMD5, as.character(senseLabelMD5), "NA"))

# reverse comp the barcodes
plates %<>% mutate(barcodesFwdRevComp=revcomp_dna(dnacol=barcodesFwd),barcodesRevRevComp=revcomp_dna(dnacol=barcodesRev))

# make fasta 
bcF <- tab2fas(df=plates, seqcol="barcodesFwd", namecol="hashLabel")
bcR <- tab2fas(df=plates, seqcol="barcodesRev", namecol="hashLabel")

# write out in folder for running analyses
write.FASTA(bcF, file=here(proj.path,"barcodes-sense.fas"))
write.FASTA(bcR, file=here(proj.path,"barcodes-antisense.fas"))

# write out table
plates %>% select(primerSet,library,eventID,hashLabel) %>% write_csv(file=here(proj.path,"results/events-hashes.csv"))

writeLines(paste0("\nDone. Sample barcodes written to ",proj.path))
