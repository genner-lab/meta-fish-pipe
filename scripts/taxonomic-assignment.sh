#!/usr/bin/env sh


# cat all the asvs
find temp/processing -name 'asvs.fna' -exec cat {} \; > temp/reference-library/asvs-cat.fna

# add md5sums
vsearch --fastx_subsample temp/reference-library/asvs-cat.fna --sample_pct 100 --relabel_md5 --relabel_keep --fasta_width 0 --fastaout temp/reference-library/asvs-relabel.fasta

# make a table of asvs
echo "md5,primer,lib,asv" > temp/reference-library/asvs.csv
grep ">" temp/reference-library/asvs-relabel.fasta | sed -e 's/>//g' -e 's/ /,/g' -e 's/-/,/g' >> temp/reference-library/asvs.csv

# derep the ASVs
vsearch --derep_fulllength temp/reference-library/asvs-relabel.fasta --minuniquesize 1 --fasta_width 0 --output temp/reference-library/asvs-derep.fasta


