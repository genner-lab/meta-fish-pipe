#!/usr/bin/env sh




############## MERGE AND DEREP ASVs ##############
############## MERGE AND DEREP ASVs ##############

# cat all the asvs
find temp/processing -name 'asvs-clean.fasta' -exec cat {} \; > temp/reference-library/asvs-clean-cat.fasta

# add md5sums
vsearch --fastx_subsample temp/reference-library/asvs-clean-cat.fasta --sample_pct 100 --relabel_md5 --relabel_keep --fasta_width 0 --fastaout temp/reference-library/asvs-clean-cat-relabel.fasta

# make a table of asvs
echo "md5,primer,lib,asv" > temp/reference-library/asvs-all.csv
grep ">" temp/reference-library/asvs-clean-cat-relabel.fasta | sed -e 's/>//g' -e 's/ /,/g' -e 's/-/,/g' >> temp/reference-library/asvs-all.csv

# derep the ASVs
vsearch --derep_fulllength temp/reference-library/asvs-clean-cat-relabel.fasta --minuniquesize 1 --fasta_width 0 --output temp/reference-library/asvs-clean-cat-relabel-derep.fasta


############## RUN SINTAX ON REFSEQ ##############
############## RUN SINTAX ON REFSEQ ##############

# run sintax
vsearch --threads 8 --sintax temp/reference-library/asvs-clean-cat-relabel-derep.fasta --db temp/reference-library/refseq-annotated.fasta --sintax_cutoff 0.7 --tabbedout temp/reference-library/sintax-output.tsv




############## MAKE BLAST DB ##############
############## MAKE BLAST DB ##############

# make blast db (only need to do this step once)
makeblastdb -in ../temp/reference-library/custom-refs.fasta -dbtype nucl -blastdb_version 5

# blast the blast db
# get better hits with smaller word size
blastn -task blastn -num_threads 8 -evalue 1000 -word_size 7 -max_target_seqs 500 -db ../reference-library/custom-refs.fasta -outfmt "6 qseqid sseqid evalue length pident nident score bitscore" -out "$DIR"/results/fish-blast.out -query "$DIR"/results/fishqueries.fasta

# join the header
echo -e "asv\tblastDbid\tblastEvalue\tblastLength\tblastPident\tblastNident\tblastScore\tblastBitscore" > "$DIR"/results/headers
cat "$DIR"/results/headers "$DIR"/results/fish-blast.out > "$DIR"/results/fish-blast-result.tsv
rm "$DIR"/results/fish-blast.out
rm "$DIR"/results/headers


