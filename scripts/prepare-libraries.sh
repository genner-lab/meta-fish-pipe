#!/usr/bin/env sh

# set params #
while getopts p:l: option
do
case "${option}"
in
p) PRIMER=${OPTARG};;
l) LIB=${OPTARG};;
esac
done


# create marker lib combo
PROJ="$PRIMER""-""$LIB"

# make path
DIR="temp/processing/""$PROJ"

# make required outputs
mkdir -p "$DIR"/fastq "$DIR"/logs "$DIR"/processed-reads/sense/dmplx "$DIR"/processed-reads/sense/trimmed-R1 "$DIR"/processed-reads/sense/trimmed-R2 "$DIR"/processed-reads/antisense/dmplx "$DIR"/processed-reads/antisense/trimmed-R1 "$DIR"/processed-reads/antisense/trimmed-R2 "$DIR"/results "$DIR"/trash

printf "...\nDirectories created for project '$PROJ'\n"
