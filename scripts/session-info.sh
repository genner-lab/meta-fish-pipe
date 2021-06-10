#!/usr/bin/env sh

# set params #
while getopts r:c: option
do
case "${option}"
in
r) REFSEQ=${OPTARG};;
c) CUSTOM=${OPTARG};;
esac
done

# make temp and results dirs
mkdir -p results temp/taxonomic-assignment/epa

# copy across refseq db
cp "$REFSEQ" temp/taxonomic-assignment/refseq-annotated.csv

# copy across custom db
cp "$CUSTOM" temp/taxonomic-assignment/custom-reference-library.csv

printf "...\nReference libraries copied\n"

# report R libs
Rscript scripts/get-r-libs.R

# add the installed software
{ printf "\n\n-----------------------\nLinux software versions\n-----------------------\n"; } >> temp/sessionInfo.txt
{ printf "\ncutadapt --version\n" & cutadapt --version; } >> temp/sessionInfo.txt
{ printf "\nvsearch --version\n" & vsearch --version 2>&1 | grep "linux"; } >> temp/sessionInfo.txt
{ printf "\nhmmsearch -h\n" & hmmsearch -h | grep "HMMER"; } >> temp/sessionInfo.txt
{ printf "\nblastn -version\n" & blastn -version | grep -v "Package"; } >> temp/sessionInfo.txt
{ printf "\nmafft --version\n" & mafft --version 2>&1; } >> temp/sessionInfo.txt
{ printf "\nraxml-ng --version\n" & raxml-ng --version | grep "RAxML-NG"; } >> temp/sessionInfo.txt
{ printf "\nepa-ng --version\n" & epa-ng --version; } >> temp/sessionInfo.txt
{ printf "\ngappa --version\n" & gappa --version; } >> temp/sessionInfo.txt
{ printf "\nseqkit version\n" & seqkit version | grep '^seqkit'; } >> temp/sessionInfo.txt

# copy to results
cp temp/sessionInfo.txt results/sessionInfo.txt

# report
printf "...\nSession info written to 'results/sessionInfo.txt'\n"
