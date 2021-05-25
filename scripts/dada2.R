#!/usr/bin/env Rscript

############## SET PARAMS ##############
############## SET PARAMS ##############

# load libs and funs
source(here::here("scripts/funs-libs.R"))

# get args
option_list <- list( 
    make_option(c("-p","--primer"), type="character"),
    make_option(c("-l","--lib"), type="character")
    )

# set args
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))

# make paths
proj.path <- here("temp/processing",paste0(opt$primer,"-",opt$lib))
dir.sense <- paste0(proj.path,"/processed-reads/sense")
dir.antisense <- paste0(proj.path,"/processed-reads/antisense")
# confirm proj path 
writeLines(paste0("\n...\nOutput directory set to:\n",proj.path,"\n"))

# trucLens
# for Miya MiFish - truncLen 105 gives at least 29 bp overlap for the longest amplicons (e.g. Raja clavata @ 181 bp), and 40 bp for the regular 170 bp
trucVal <- c(105,105)


############## QUALITY TRIM TRUNCATE ##############
############## QUALITY TRIM TRUNCATE ##############

# report
writeLines("\n...\nQuality trimming and truncating\n")
Sys.sleep(3)

# quality trim Ns and truncate
filterAndTrim(fwd=cpath("sense","trimmed","R1"), filt=cpath("sense","filtered","R1"), rev=cpath("sense","trimmed","R2"), filt.rev=cpath("sense","filtered","R2"), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, truncLen=trucVal, multithread=TRUE, verbose=TRUE, matchIDs=TRUE)
filterAndTrim(fwd=cpath("antisense","trimmed","R1"), filt=cpath("antisense","filtered","R1"), rev=cpath("antisense","trimmed","R2"), filt.rev=cpath("antisense","filtered","R2"), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, truncLen=trucVal, multithread=TRUE, verbose=TRUE, matchIDs=TRUE)


############## LEARN ERRORS ##############
############## LEARN ERRORS ##############

# report
writeLines("\n...\nLearning errors\n")
Sys.sleep(3)

# learn errors
set.seed(42)
sense.filt.R1.errs <- learnErrors(cpath("sense","filtered","R1"), multithread=TRUE, randomize=TRUE, nbases=1e+08, verbose=TRUE)
set.seed(42)
sense.filt.R2.errs <- learnErrors(cpath("sense","filtered","R2"), multithread=TRUE, randomize=TRUE, nbases=1e+08, verbose=TRUE)
set.seed(42)
antisense.filt.R1.errs <- learnErrors(cpath("antisense","filtered","R1"), multithread=TRUE, randomize=TRUE, nbases=1e+08, verbose=TRUE)
set.seed(42)
antisense.filt.R2.errs <- learnErrors(cpath("antisense","filtered","R2"), multithread=TRUE, randomize=TRUE, nbases=1e+08, verbose=TRUE)


############## MAKE PRIORS ##############
############## MAKE PRIORS ##############

# make some fish priors
if(opt$primer=="tele02") {
    prefix <- "12s.taberlet.noprimers"
} else if(opt$primer=="elas02" | opt$primer=="mifish-u" | opt$primer=="mifish-u-mod") {
    prefix <- "12s.miya.noprimers"
} else {
    stop("Primers must be 'tele02', 'elas02', 'mifish-u', or 'mifish-u-mod'")
}

# load reflib - first need to have run 'join-references.R' to create
reflib <- suppressMessages(suppressWarnings(read_csv(here("assets/custom-reference-library.csv"),guess_max=99999,col_types=cols())))

# subset the marker from the reflib and rev comp and trim
reflib.sub <- subset_by_marker(prefix=prefix,df=reflib,thresh=0.625) %>% 
    select(sciNameValid,!!as.name(paste0("nucleotidesFrag.",prefix))) %>% 
    rename(seqs=!!as.name(paste0("nucleotidesFrag.",prefix))) %>%
    filter(!str_detect(seqs,"[^actg]")) %>%
    mutate(fwd=toupper(str_trunc(seqs,width=trucVal[1],side="right",ellipsis=""))) %>%
    mutate(revcomp=toupper(str_trunc(revcomp_dna(seqs),width=trucVal[1],side="right",ellipsis=""))) 

# pull out the priors
fish.priors <- unique(c(pull(reflib.sub,fwd),pull(reflib.sub,revcomp)))


############## DENOISE ##############
############## DENOISE ##############

# report
writeLines("\n...\ndada2 denoising\n")
Sys.sleep(3)

# run dada denoising - takes time with pool=TRUE
sense.filt.R1.dada <- dada(cpath("sense","filtered","R1"), err=sense.filt.R1.errs, multithread=TRUE, pool=TRUE, priors=fish.priors)
sense.filt.R2.dada <- dada(cpath("sense","filtered","R2"), err=sense.filt.R2.errs, multithread=TRUE, pool=TRUE, priors=fish.priors)
antisense.filt.R1.dada <- dada(cpath("antisense","filtered","R1"), err=antisense.filt.R1.errs, multithread=TRUE, pool=TRUE, priors=fish.priors)
antisense.filt.R2.dada <- dada(cpath("antisense","filtered","R2"), err=antisense.filt.R2.errs, multithread=TRUE, pool=TRUE, priors=fish.priors)


############## DEREPLICATE ##############
############## DEREPLICATE ##############

# report
writeLines("\n...\nDereplication\n")
Sys.sleep(3)

# derep
sense.filt.R1.derep <- derepFastq(cpath("sense","filtered","R1"))
sense.filt.R2.derep <- derepFastq(cpath("sense","filtered","R2"))
antisense.filt.R1.derep <- derepFastq(cpath("antisense","filtered","R1"))
antisense.filt.R2.derep <- derepFastq(cpath("antisense","filtered","R2"))


############## MERGE ##############
############## MERGE ##############

# report
writeLines("\n...\nRead merging\n")
Sys.sleep(3)

# merge the R1 and R2
sense.merged <- mergePairs(dadaF=sense.filt.R1.dada, derepF=sense.filt.R1.derep, dadaR=sense.filt.R2.dada, derepR=sense.filt.R2.derep, verbose=TRUE, maxMismatch=0)
antisense.merged <- mergePairs(dadaF=antisense.filt.R1.dada, derepF=antisense.filt.R1.derep, dadaR=antisense.filt.R2.dada,  derepR=antisense.filt.R2.derep, verbose=TRUE, maxMismatch=0)

# make an OTU table
sense.seqtab <- makeSequenceTable(sense.merged)
antisense.seqtab <- makeSequenceTable(antisense.merged)

# reverse comp the antisense
colnames(antisense.seqtab) <- dada2::rc(colnames(antisense.seqtab))

# fix the names before merging
rownames(sense.seqtab) <- str_split_fixed(rownames(sense.seqtab),"\\.",4)[,1]
rownames(antisense.seqtab) <- str_split_fixed(rownames(antisense.seqtab),"\\.",4)[,1]

# merge the tables
merged.seqtab <- mergeSequenceTables(table1=sense.seqtab, table2=antisense.seqtab, repeats="sum")


############## REMOVE CHIMAERAS ##############
############## REMOVE CHIMAERAS ##############

# report
writeLines("\n...\nDetecting chimaeras\n")
Sys.sleep(3)

# remove chimaeras
merged.seqtab.nochim <- removeBimeraDenovo(merged.seqtab, method="consensus", multithread=TRUE, verbose=TRUE)


############## SAVE FILES ##############
############## SAVE FILES ##############

# report
writeLines("\n...\nSaving raw ASV files\n")
Sys.sleep(3)

# make df and fasta for IDs
otus.df <- tibble(names=paste0(opt$primer,"-",opt$lib,"-","asv",str_pad(seq_along(colnames(merged.seqtab.nochim)),width=4,side="left",pad="0")), dnas=colnames(merged.seqtab.nochim)) %>% mutate(len=str_length(dnas))

# write out
write.FASTA(tab2fas(df=otus.df, seqcol="dnas", namecol="names"), file=paste0(proj.path,"/results/asvs.fna"))

# save the OTU table as df
colnames(merged.seqtab.nochim) <- paste0(opt$primer,"-",opt$lib,"-","asv",str_pad(seq_along(colnames(merged.seqtab.nochim)),width=4,side="left",pad="0"))
write_tsv(as_tibble(t(merged.seqtab.nochim), rownames="asv"), file=paste0(proj.path,"/results/asv-table.tsv"))

# report
writeLines(paste0("\n...\nASVs written to:\n",paste0(proj.path,"/results/asvs.fna"),"\n"))
writeLines(paste0("\n...\nASV table written to:\n",paste0(proj.path,"/results/asv-table.tsv"),"\n"))


############## CLEAN ASVs WITH HMMs ##############
############## CLEAN ASVs WITH HMMs ##############

# report
writeLines("\n...\nCleaning ASVs\n")
Sys.sleep(3)

# get correct hmm
if(opt$primer=="tele02") {
    hmm <- here("assets/12s.taberlet.noprimers.hmm")
} else if(opt$primer=="elas02" | opt$primer=="mifish-u" | opt$primer=="mifish-u-mod") {
    hmm <- here("assets/12s.miya.noprimers.hmm")
} else {
    stop("Primers must be 'tele02', 'elas02', 'mifish-u', or 'mifish-u-mod'")
}

# set vars
infile <- paste0(proj.path,"/results/asvs.fna")
outfile <- paste0(proj.path,"/results/hmm.out")
cleanfile <- paste0(proj.path,"/results/asvs-clean.fasta")

# run hmm clean
string.hmmer <- paste0("hmmsearch -E 0.01 --incE 0.01 ",hmm," ",infile," | grep '>>' | sed -e 's/>> //g' -e 's/[[:space:]]//g' -e 's/$/$/' | sort | uniq"," > ",outfile)
system(command=string.hmmer,ignore.stdout=FALSE)

# run the grep clean 
string.grep <- paste0("grep -A 1 -f ",outfile," ",infile," | grep -v '^-' > ",cleanfile)
system(command=string.grep,ignore.stdout=FALSE)


# read asv table back in written table to get format correct
otu.tab <- as_tibble(read.table(file=paste0(proj.path,"/results/asv-table.tsv"), sep="\t", header=TRUE, as.is=TRUE, row.names=1, check.names=FALSE), rownames="asv")
dnas.curated <- read.FASTA(file=paste0(proj.path,"/results/asvs.fna"))

# load up data cleaned by the HMM
dnas.curated.clean <- read.FASTA(file=paste0(proj.path,"/results/asvs-clean.fasta"))
dnas.curated.dirty <- dnas.curated[which(!(names(dnas.curated) %in% names(dnas.curated.clean)))]
write.FASTA(dnas.curated.dirty, file=paste0(proj.path,"/results/asvs-dirty.fasta"))

# subset the clean sequences from the ASV table
otu.tab.clean <- otu.tab %>% filter(asv %in% names(dnas.curated.clean))
# write out
write_tsv(otu.tab.clean,file=paste0(proj.path,"/results/asv-table-clean.tsv"))

# get numbers of seqs lost during hmm search
total <- otu.tab %>% summarise_if(is.numeric, sum, na.rm=TRUE) %>% rowSums()
lost <- otu.tab %>% filter(!asv %in% names(dnas.curated.clean)) %>% summarise_if(is.numeric, sum, na.rm=TRUE) %>% rowSums()
retained <- otu.tab.clean %>% summarise_if(is.numeric, sum, na.rm=TRUE) %>% rowSums()

# report
writeLines(paste0("\n...\nCleaned ASVs written to:\n",paste0(proj.path,"/results/asvs-clean.fasta"),"\n"))
writeLines(paste0("\n...\nCleaned ASV table written to:\n",paste0(proj.path,"/results/asv-table-clean.tsv"),"\n"))

# write stats
stats <- paste0("merge,",sum(merged.seqtab),"\n","chim,",sum(merged.seqtab.nochim),"\n","homol,",retained)
writeLines(stats,paste0(proj.path,"/logs/dada-stats.csv"))
