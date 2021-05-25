#!/usr/bin/env Rscript

# load libs and funs
source(here::here("scripts/funs-libs.R"))

# load up the alignment
mafft.alignment <- read.FASTA(file=here("temp/taxonomic-assignment/epa/epa.aligned.fasta"))

# grep queries
mafft.alignment.queries <- mafft.alignment[grep("\\.queries",names(mafft.alignment))]
mafft.alignment.refs <- mafft.alignment[-grep("\\.queries",names(mafft.alignment))]

# remove suffix
names(mafft.alignment.queries) <- gsub("\\.queries","",names(mafft.alignment.queries))

# write DNAs

# write out
write.FASTA(mafft.alignment.queries,file=here("temp/taxonomic-assignment/epa/epa.queries.fasta"))
write.FASTA(mafft.alignment.refs,file=here("temp/taxonomic-assignment/epa/epa.references.fasta"))

# report
writeLines("\n...\nDNA split and formatted\n")
