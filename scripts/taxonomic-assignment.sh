#!/usr/bin/env sh

# set params #
while getopts t:p: option
do
case "${option}"
in
t) THREADS=${OPTARG};;
p) PRIMER=${OPTARG};;
esac
done

############## PREP REFERENCE LIBRARIES ##############
############## PREP REFERENCE LIBRARIES ##############

# format and subset reference libraries
scripts/prep-reflibs.R -p "$PRIMER"


############## MERGE AND DEREP ASVs ##############
############## MERGE AND DEREP ASVs ##############

# cat all the asvs
find temp/processing -name 'asvs-clean.fasta' -exec cat {} \; > temp/taxonomic-assignment/asvs-clean-cat.fasta

# add md5sums
vsearch --fastx_subsample temp/taxonomic-assignment/asvs-clean-cat.fasta --sample_pct 100 --relabel_md5 --relabel_keep --fasta_width 0 --fastaout temp/taxonomic-assignment/asvs-clean-cat-relabel.fasta

# make a table of asvs
echo "md5,primer,lib,asv" > temp/taxonomic-assignment/asvs-all.csv
grep ">" temp/taxonomic-assignment/asvs-clean-cat-relabel.fasta | sed -e 's/>//g' -e 's/ /,/g' -e 's/|/,/g' >> temp/taxonomic-assignment/asvs-all.csv

# derep the ASVs
vsearch --derep_fulllength temp/taxonomic-assignment/asvs-clean-cat-relabel.fasta --minuniquesize 1 --fasta_width 0 --output temp/taxonomic-assignment/asvs-clean-cat-relabel-derep.fasta


############## RUN SINTAX ON REFSEQ ##############
############## RUN SINTAX ON REFSEQ ##############

# run sintax
vsearch --threads "$THREADS" --sintax temp/taxonomic-assignment/asvs-clean-cat-relabel-derep.fasta --db temp/taxonomic-assignment/refseq-annotated.fasta --sintax_cutoff 0.7 --tabbedout temp/taxonomic-assignment/sintax-output.tsv


############## EXTRACT FISH ASVs ##############
############## EXTRACT FISH ASVs ##############

scripts/extract-fishes.R


############## RUN BLAST ON FISH ASVs ##############
############## RUN BLAST ON FISH ASVs ##############

# make blast db (only need to do this step once)
makeblastdb -in temp/taxonomic-assignment/custom-reference-library.fasta -dbtype nucl -blastdb_version 5

# blast the blast db
# get better hits with smaller word size
blastn -task blastn -num_threads "$THREADS" -evalue 1000 -word_size 7 -max_target_seqs 500 -db temp/taxonomic-assignment/custom-reference-library.fasta -outfmt "6 qseqid sseqid evalue length pident nident score bitscore" -out temp/taxonomic-assignment/fish-blast.out -query temp/taxonomic-assignment/asvs-fish.fasta

# join the header
echo -e 'asv\tblastDbid\tblastEvalue\tblastLength\tblastPident\tblastNident\tblastScore\tblastBitscore' | sed -e "s/-e //g" > temp/taxonomic-assignment/headers
cat temp/taxonomic-assignment/headers temp/taxonomic-assignment/fish-blast.out > temp/taxonomic-assignment/fish-blast-result.tsv
rm temp/taxonomic-assignment/fish-blast.out
rm temp/taxonomic-assignment/headers


############## PROCESS BLAST RESULTS ##############
############## PROCESS BLAST RESULTS ##############

# process blast
scripts/process-blast.R


############## ALIGN AND RUN EPA ##############
############## ALIGN AND RUN EPA ##############

# prep for EPA
scripts/prep-epa.R

# align with mafft (accurate)
mafft --thread "$THREADS" --retree 2 --maxiterate 2 temp/taxonomic-assignment/epa/epa.input.fasta > temp/taxonomic-assignment/epa/epa.aligned.fasta

# split up references and queries for epa
scripts/split.R

# make RAxML trees for all markers
# create a binary and set model
raxml-ng --parse --msa temp/taxonomic-assignment/epa/epa.references.fasta --model TN93+G --redo

# run a tree search
raxml-ng --search --msa temp/taxonomic-assignment/epa/epa.references.fasta.raxml.rba --tree pars{1} --seed 42 --threads 1 --redo

# optimise params
raxml-ng --evaluate --msa temp/taxonomic-assignment/epa/epa.references.fasta.raxml.rba --tree temp/taxonomic-assignment/epa/epa.references.fasta.raxml.rba.raxml.bestTree --prefix temp/taxonomic-assignment/epa/opt --redo

# run epa-ng
epa-ng --ref-msa temp/taxonomic-assignment/epa/epa.references.fasta --tree temp/taxonomic-assignment/epa/epa.references.fasta.raxml.rba.raxml.bestTree --query temp/taxonomic-assignment/epa/epa.queries.fasta --outdir temp/taxonomic-assignment/epa --model temp/taxonomic-assignment/epa/opt.raxml.bestModel --redo --preserve-rooting off

# run gappa
gappa examine assign --per-query-results --jplace-path temp/taxonomic-assignment/epa/epa_result.jplace --taxon-file temp/taxonomic-assignment/epa/references.taxonomy.tsv --out-dir temp/taxonomic-assignment/epa --allow-file-overwriting

# rename best tree to keep
mv temp/taxonomic-assignment/epa/epa.references.fasta.raxml.rba.raxml.bestTree temp/taxonomic-assignment/epa/references.tree.nwk

# clean up
rm temp/taxonomic-assignment/epa/*.log
rm temp/taxonomic-assignment/epa/*raxml*

# process results
scripts/process-epa.R

# report
sleep 3
printf "...\nTaxonomy assignments completed\n"
