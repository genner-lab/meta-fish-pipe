# meta-fish-pipe
Bioinformatics pipeline for fish metabarcoding

# get date and time
then=$(date)

### FOR LIB3 ###

# prep
scripts/prepare.sh -p tele02 -l lib3

# make symlinks
#ln -s ~/Projects/SeaDNA/temp-local-only/fastq/reduced/lib3/R1.fastq.gz ~/Projects/genner-lab/meta-fish-pipe/temp/processing/tele02-lib3/fastq/R1.fastq.gz
#ln -s ~/Projects/SeaDNA/temp-local-only/fastq/reduced/lib3/R2.fastq.gz ~/Projects/genner-lab/meta-fish-pipe/temp/processing/tele02-lib3/fastq/R2.fastq.gz
# uni
ln -s /media/1TB/rc16041/Projects-temp-local-only/seadna-temp-local-only/fastq/reduced/lib3/R1.fastq.gz ~/Projects/genner-lab/meta-fish-pipe/temp/processing/tele02-lib3/fastq/R1.fastq.gz
ln -s /media/1TB/rc16041/Projects-temp-local-only/seadna-temp-local-only/fastq/reduced/lib3/R2.fastq.gz ~/Projects/genner-lab/meta-fish-pipe/temp/processing/tele02-lib3/fastq/R2.fastq.gz

# generate barcodes
scripts/generate-barcodes.R -p tele02 -l lib3 -f 18 -r 20

# demux
scripts/demultiplex.sh -p tele02 -l lib3 -f AAACTCGTGCCAGCCACC -r GGGTATCTAATCCCAGTTTG -t 8 -m 18

# denoise with dada2
scripts/dada2.R -p tele02 -l lib3

# generate stats
scripts/generate-stats.sh -p tele02 -l lib3 -t 8


### FOR LIB4 ###

# prep
scripts/prepare.sh -p tele02 -l lib4

# make symlinks HOME
#ln -s ~/Projects/SeaDNA/temp-local-only/fastq/reduced/lib4/R1.fastq.gz ~/Projects/genner-lab/meta-fish-pipe/temp/processing/tele02-lib4/fastq/R1.fastq.gz
#ln -s ~/Projects/SeaDNA/temp-local-only/fastq/reduced/lib4/R2.fastq.gz ~/Projects/genner-lab/meta-fish-pipe/temp/processing/tele02-lib4/fastq/R2.fastq.gz
# uni
ln -s /media/1TB/rc16041/Projects-temp-local-only/seadna-temp-local-only/fastq/reduced/lib4/R1.fastq.gz ~/Projects/genner-lab/meta-fish-pipe/temp/processing/tele02-lib4/fastq/R1.fastq.gz
ln -s /media/1TB/rc16041/Projects-temp-local-only/seadna-temp-local-only/fastq/reduced/lib4/R2.fastq.gz ~/Projects/genner-lab/meta-fish-pipe/temp/processing/tele02-lib4/fastq/R2.fastq.gz

# generate barcodes
scripts/generate-barcodes.R -p tele02 -l lib4 -f 18 -r 20

# demux
scripts/demultiplex.sh -p tele02 -l lib4 -f AAACTCGTGCCAGCCACC -r GGGTATCTAATCCCAGTTTG -t 8 -m 18

# denoise with dada2
scripts/dada2.R -p tele02 -l lib4

# generate stats
scripts/generate-stats.sh -p tele02 -l lib4 -t 8



############## ALL LIBS - TAXONOMIC ASSIGNMENT ##############
############## ALL LIBS - TAXONOMIC ASSIGNMENT ##############

# run taxonomic assignment
scripts/taxonomic-assignment.sh -t 8 -p tele02 -r assets/refseq206-annotated-tele02.csv -c assets/custom-reference-library.csv

# assemble results
scripts/assemble-results.R


# print date and time again
now=$(date)
echo "$then"
echo "$now"

# clean up (optional) 
rm -r results temp



############## JOBS TO DO #####

- documentation
- update all the software and test


# grep a function
grep -r "hap_collapse_df" .


### refseq CHECK

# trim primers from the custom reference library
cutadapt -n 1 -e 0.3 -O 10 -g AAACTCGTGCCAGCCACC temp/reference-library/refseq-annotated.fasta | cutadapt --minimum-length 144 --maximum-length 216 -n 1 -e 0.3 -O 10 -a  CAAACTGGGATTAGATACCC -o temp/reference-library/refseq-annotated-trimmed.fasta -

# my files
mkdir -p reduced/lib3 reduced/lib4
gzip -cd SeaDNA_Tele02_Lib03v2_R1.fastq.gz | head -n 400000 | gzip > reduced/lib3/R1.fastq.gz
gzip -cd SeaDNA_Tele02_Lib03v2_R2.fastq.gz | head -n 400000 | gzip > reduced/lib3/R2.fastq.gz
gzip -cd SeaDNA_Teleo02_Lib-04_S2_L001_R1_001.fastq.gz | head -n 400000 | gzip > reduced/lib4/R1.fastq.gz
gzip -cd SeaDNA_Teleo02_Lib-04_S2_L001_R2_001.fastq.gz | head -n 400000 | gzip > reduced/lib4/R2.fastq.gz
