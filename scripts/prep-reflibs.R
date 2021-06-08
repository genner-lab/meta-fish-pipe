#!/usr/bin/env Rscript

# report
writeLines("\nPreparing reference library sequences ...\n")

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

# to test
# prefix <- "12s.taberlet.noprimers" 

# import reference library
custom.refs <- suppressMessages(suppressWarnings(read_csv(here("temp/taxonomic-assignment/custom-reference-library.csv"),guess_max=99999,col_types=cols())))

# subset the marker from the custom reference library reflib - using 50% length cutoff
custom.refs.sub <- subset_by_marker(prefix=prefix,df=custom.refs,thresh=0.5)

# format custom
custom.refs.sub.clean <- custom.refs.sub %>% 
    mutate(source="custom",phylum="Chordata") %>%
    select(source,dbid,kingdom,phylum,class,order,family,genus,sciNameValid,paste0("lengthFrag.",prefix),paste0("nucleotidesFrag.",prefix)) %>%
    rename(length=paste0("lengthFrag.",prefix),nucleotides=paste0("nucleotidesFrag.",prefix))


# import refseq
refseq.refs <- suppressMessages(suppressWarnings(read_csv(here("temp/taxonomic-assignment/refseq-annotated.csv"),guess_max=99999,col_types=cols())))

# format refseq
refseq.refs.clean <- refseq.refs %>% 
    mutate(source="refseq") %>% 
    select(source,accession,kingdom,phylum,class,order,family,genus,scientificName,length,nucleotides) %>% 
    rename(dbid=accession,sciNameValid=scientificName)


# join custom and refseq
refseq.custom.combined <- bind_rows(custom.refs.sub.clean,refseq.refs.clean)

# dereplicate the combined refseq.custom and add labels
combined.derep <- refseq.custom.combined %>% 
    group_by(sciNameValid) %>% 
    group_modify(~ hap_collapse_df(df=.x,lengthcol="length",nuccol="nucleotides",cores=1)) %>% 
    ungroup() %>% 
    arrange(kingdom,phylum,class,order,family,genus,sciNameValid,source) %>% 
    mutate(label=paste0(dbid,";tax=k:",kingdom,",p:",phylum,",c:",class,",o:",order,",f:",family,",g:",genus,",s:",sciNameValid)) %>% 
    mutate(label=str_replace_all(label," ","_"))

# dereplicate the custom only
custom.refs.sub.clean.derep <- custom.refs.sub.clean %>% 
    group_by(sciNameValid) %>% 
    group_modify(~ hap_collapse_df(df=.x,lengthcol="length",nuccol="nucleotides",cores=1)) %>% 
    ungroup() %>% 
    arrange(kingdom,phylum,class,order,family,genus,sciNameValid)


# write out the combined reference library fasta for SINTAX
tab2fas(combined.derep, seqcol="nucleotides",namecol="label") %>% write.FASTA(file=here("temp/taxonomic-assignment/refseq-annotated.fasta"))

# write out the custom reference library fasta for BLAST and EPA
tab2fas(custom.refs.sub.clean.derep, seqcol="nucleotides",namecol="dbid") %>% write.FASTA(file=here("temp/taxonomic-assignment/custom-reference-library.fasta"))

# write out the custom reference library csv for BLAST and EPA
custom.refs.sub.clean.derep %>% write_csv(file=here("temp/taxonomic-assignment/custom-reference-library-reduced.csv"))

# report
writeLines("\n...\nReference library sequences prepared for assignment\n")
