#!/usr/bin/env sh

# set params #
while getopts p:l:f:r:t:m: option
do
case "${option}"
in
p) PRIMER=${OPTARG};;
l) LIB=${OPTARG};;
f) FWD=${OPTARG};;
r) REV=${OPTARG};;
t) THREADS=${OPTARG};;
m) MINLEN=${OPTARG};;
esac
done

# create marker lib combo
PROJ="$PRIMER""-""$LIB"
# make path
DIR="temp/processing/""$PROJ"
# set minlen
TRUCVAL="105"

# make grepable primers
FWDGR="$(echo "$FWD" | sed -e 's/W/[A\|T]/g' -e 's/R/[A\|G]/g' -e 's/Y/[C\|T]/g' -e 's/M/[A\|C]/g' -e 's/H/[A\|C\|T]/g' -e 's/N/[A\|C\|G\|T]/g' -e 's/$/\|\$/g')"
REVGR="$(echo "$REV" | sed -e 's/W/[A\|T]/g' -e 's/R/[A\|G]/g' -e 's/Y/[C\|T]/g' -e 's/M/[A\|C]/g' -e 's/H/[A\|C\|T]/g' -e 's/N/[A\|C\|G\|T]/g' -e 's/$/\|\$/g')"

# confirm
printf "...\nOutput directory set to '$PROJ'\n"
sleep 3

printf "...\nDetermining read orientation\n"
sleep 3



# sense
cutadapt -j "$THREADS" --no-indels --error-rate 0 --overlap "$MINLEN" --pair-adapters --action=none -g senseF="$FWD" -G senseR="$REV" --untrimmed-output "$DIR"/trash/untrimmed.R1.fastq.gz --untrimmed-paired-output "$DIR"/trash/untrimmed.R2.fastq.gz -o "$DIR"/processed-reads/sense/R1.fastq.gz -p "$DIR"/processed-reads/sense/R2.fastq.gz "$DIR"/fastq/R1.fastq.gz "$DIR"/fastq/R2.fastq.gz > "$DIR"/logs/cutadapt.sense.log

# antisense
cutadapt -j "$THREADS" --no-indels --error-rate 0 --overlap "$MINLEN" --pair-adapters --action=none -g antisenseF="$REV" -G antisenseR="$FWD" --untrimmed-output "$DIR"/trash/noprimer.R1.fastq.gz --untrimmed-paired-output "$DIR"/trash/noprimer.R2.fastq.gz -o "$DIR"/processed-reads/antisense/R1.fastq.gz -p "$DIR"/processed-reads/antisense/R2.fastq.gz "$DIR"/trash/untrimmed.R1.fastq.gz "$DIR"/trash/untrimmed.R2.fastq.gz > "$DIR"/logs/cutadapt.antisense.log



# check the primers
printf "...\nChecking forward primers\n"
sleep 3
gzip -cd "$DIR"/processed-reads/sense/R1.fastq.gz | sed -n '2~4p' | head -n 20 | GREP_COLOR="1;44" grep --color=always -E "$FWDGR" | GREP_COLOR="1;41" grep --color=always -E "$REVGR"
printf "...\nChecking reverse primers\n"
sleep 3
gzip -cd "$DIR"/processed-reads/sense/R2.fastq.gz | sed -n '2~4p' | head -n 20 | GREP_COLOR="1;44" grep --color=always -E "$FWDGR" | GREP_COLOR="1;41" grep --color=always -E "$REVGR"



# demuxing
printf "...\nDemultiplexing reads\n"
sleep 3

# combinatorial paired demultiplex on sense
cutadapt -j "$THREADS" --no-indels --error-rate 0 --overlap 10 --pair-adapters --action=none -g file:"$DIR"/barcodes-sense.fas -G file:"$DIR"/barcodes-antisense.fas --discard-untrimmed -o "$DIR"/processed-reads/sense/dmplx/{name}.R1.fastq.gz -p "$DIR"/processed-reads/sense/dmplx/{name}.R2.fastq.gz "$DIR"/processed-reads/sense/R1.fastq.gz "$DIR"/processed-reads/sense/R2.fastq.gz > "$DIR"/logs/cutadapt.dmplx.barcodes.sense.log

# combinatorial paired demultiplex on antisense
cutadapt -j "$THREADS" --no-indels --error-rate 0 --overlap 10 --pair-adapters --action=none -g file:"$DIR"/barcodes-antisense.fas -G file:"$DIR"/barcodes-sense.fas --discard-untrimmed -o "$DIR"/processed-reads/antisense/dmplx/{name}.R1.fastq.gz -p "$DIR"/processed-reads/antisense/dmplx/{name}.R2.fastq.gz "$DIR"/processed-reads/antisense/R1.fastq.gz "$DIR"/processed-reads/antisense/R2.fastq.gz > "$DIR"/logs/cutadapt.dmplx.barcodes.antisense.log



# trim
printf "...\nTrimming reads to minimum length '$TRUCVAL'\n"
sleep 3

# trim SENSE with cutadapt
sense="$(basename -a -s .R1.fastq.gz "$DIR"/processed-reads/sense/dmplx/*.R1.fastq.gz | uniq)"
# now run
for i in $sense; do
cutadapt -n 5 --error-rate 0.15 --minimum-length "$TRUCVAL" -g "$FWD" -G "$REV" -o "$DIR"/processed-reads/sense/trimmed-R1/"$i".R1.fastq.gz -p "$DIR"/processed-reads/sense/trimmed-R2/"$i".R2.fastq.gz --discard-untrimmed "$DIR"/processed-reads/sense/dmplx/"$i".R1.fastq.gz "$DIR"/processed-reads/sense/dmplx/"$i".R2.fastq.gz >> "$DIR"/logs/cutadapt.sense.trimming.log
done

# trim ANTISENSE with cutadapt
antisense="$(basename -a -s .R1.fastq.gz "$DIR"/processed-reads/antisense/dmplx/*.R1.fastq.gz | uniq)"
# now run
for i in $antisense; do
cutadapt -n 5 --error-rate 0.15 --minimum-length "$TRUCVAL" -g "$REV" -G "$FWD" -o "$DIR"/processed-reads/antisense/trimmed-R1/"$i".R1.fastq.gz -p "$DIR"/processed-reads/antisense/trimmed-R2/"$i".R2.fastq.gz --discard-untrimmed "$DIR"/processed-reads/antisense/dmplx/"$i".R1.fastq.gz "$DIR"/processed-reads/antisense/dmplx/"$i".R2.fastq.gz >> "$DIR"/logs/cutadapt.antisense.trimming.log
done

sleep 3
printf "...\nDone\n"
