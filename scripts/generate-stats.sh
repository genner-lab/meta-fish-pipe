#!/usr/bin/env sh

# set params #
while getopts p:l:t: option
do
case "${option}"
in
p) PRIMER=${OPTARG};;
l) LIB=${OPTARG};;
t) THREADS=${OPTARG};;
esac
done

# create marker lib combo
PROJ="$PRIMER""-""$LIB"
# make path
DIR="temp/processing/""$PROJ"

# get number of raw reads passing filter
PF="$(seqkit -j 8 stats -T -b "$DIR"/fastq/R1.fastq.gz | grep -v 'file' | cut -f4)"
echo "$PROJ",pf,"$PF" >> "$DIR"/logs/stats.csv

# get number reads with primer
ORIENT="$(cat "$DIR"/processed-reads/sense/R1.fastq.gz "$DIR"/processed-reads/antisense/R1.fastq.gz | seqkit -j 8 stats -T | grep -v 'file' | cut -f4)"
echo "$PROJ",primer,"$ORIENT" >> "$DIR"/logs/stats.csv

# get number reads with barcodes
BARCODE="$(cat "$DIR"/processed-reads/sense/dmplx/*.R1.fastq.gz "$DIR"/processed-reads/antisense/dmplx/*.R1.fastq.gz | seqkit -j 8 stats -T | grep -v 'file' | cut -f4)"
echo "$PROJ",barcode,"$BARCODE" >> "$DIR"/logs/stats.csv

# get number reads after trimming and minlen
TRIM="$(cat "$DIR"/processed-reads/sense/trimmed-R1/*.fastq.gz "$DIR"/processed-reads/antisense/trimmed-R1/*.fastq.gz | seqkit -j 8 stats -T | grep -v 'file' | cut -f4)"
echo "$PROJ",trim,"$TRIM" >> "$DIR"/logs/stats.csv

# get number reads after dada quality filtering
FILTER="$(cat "$DIR"/processed-reads/sense/filtered-R1/*.fastq.gz "$DIR"/processed-reads/antisense/filtered-R1/*.fastq.gz | seqkit -j 8 stats -T | grep -v 'file' | cut -f4)"
echo "$PROJ",filter,"$FILTER" >> "$DIR"/logs/stats.csv

sleep 3
printf "...\nStats generated\n"
