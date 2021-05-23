# meta-fish-pipe
Bioinformatics pipeline for fish metabarcoding



### FOR LIB3

# prep
scripts/prepare.sh -p tele02 -l lib3

# make symlinks
ln -s ~/Projects/SeaDNA/temp-local-only/fastq/reduced/lib3/R1.fastq.gz ~/Projects/genner-lab/meta-fish-pipe/temp/processing/tele02-lib3/fastq/R1.fastq.gz
ln -s ~/Projects/SeaDNA/temp-local-only/fastq/reduced/lib3/R2.fastq.gz ~/Projects/genner-lab/meta-fish-pipe/temp/processing/tele02-lib3/fastq/R2.fastq.gz

# generate barcodes
scripts/generate-barcodes.R -p tele02 -l lib3 -f 18 -r 20

# demux
scripts/demultiplex.sh -p tele02 -l lib3 -f AAACTCGTGCCAGCCACC -r GGGTATCTAATCCCAGTTTG -t 8 -m 18

# denoise with dada2
scripts/dada2.R -p tele02 -l lib3



### FOR LIB4

# prep
scripts/prepare.sh -p tele02 -l lib4

# make symlinks
ln -s ~/Projects/SeaDNA/temp-local-only/fastq/reduced/lib4/R1.fastq.gz ~/Projects/genner-lab/meta-fish-pipe/temp/processing/tele02-lib4/fastq/R1.fastq.gz
ln -s ~/Projects/SeaDNA/temp-local-only/fastq/reduced/lib4/R2.fastq.gz ~/Projects/genner-lab/meta-fish-pipe/temp/processing/tele02-lib4/fastq/R2.fastq.gz

# generate barcodes
scripts/generate-barcodes.R -p tele02 -l lib4 -f 18 -r 20

# demux
scripts/demultiplex.sh -p tele02 -l lib4 -f AAACTCGTGCCAGCCACC -r GGGTATCTAATCCCAGTTTG -t 8 -m 18

# denoise with dada2
scripts/dada2.R -p tele02 -l lib4





renv::install("Matrix@1.3-2")


# my files
gzip -cd SeaDNA_Tele02_Lib03v2_R1.fastq.gz | head -n 4000000 | gzip > reduced/lib3/R1.fastq.gz
gzip -cd SeaDNA_Tele02_Lib03v2_R2.fastq.gz | head -n 4000000 | gzip > reduced/lib3/R2.fastq.gz
gzip -cd SeaDNA_Teleo02_Lib-04_S2_L001_R1_001.fastq.gz | head -n 4000000 | gzip > reduced/lib4/R1.fastq.gz
gzip -cd SeaDNA_Teleo02_Lib-04_S2_L001_R2_001.fastq.gz | head -n 4000000 | gzip > reduced/lib4/R2.fastq.gz





FWD="AAACTCGTGCCAGCCACC"
REV="GGGTATCTAATCCCAGTTTG"