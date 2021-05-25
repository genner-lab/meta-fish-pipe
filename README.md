# meta-fish-pipe
Bioinformatics pipeline for fish metabarcoding



### FOR LIB3

# prep
scripts/prepare.sh -p tele02 -l lib3

# make symlinks
ln -s ~/Projects/SeaDNA/temp-local-only/fastq/reduced/lib3/R1.fastq.gz ~/Projects/genner-lab/meta-fish-pipe/temp/processing/tele02-lib3/fastq/R1.fastq.gz
ln -s ~/Projects/SeaDNA/temp-local-only/fastq/reduced/lib3/R2.fastq.gz ~/Projects/genner-lab/meta-fish-pipe/temp/processing/tele02-lib3/fastq/R2.fastq.gz
# uni
ln -s /media/1TB/rc16041/Projects-temp-local-only/seadna-temp-local-only/fastq/reduced/lib3/R1.fastq.gz ~/Projects/genner-lab/meta-fish-pipe/temp/processing/tele02-lib3/fastq/R1.fastq.gz
ln -s /media/1TB/rc16041/Projects-temp-local-only/seadna-temp-local-only/fastq/reduced/lib3/R2.fastq.gz ~/Projects/genner-lab/meta-fish-pipe/temp/processing/tele02-lib3/fastq/R2.fastq.gz

# generate barcodes
scripts/generate-barcodes.R -p tele02 -l lib3 -f 18 -r 20

# demux
scripts/demultiplex.sh -p tele02 -l lib3 -f AAACTCGTGCCAGCCACC -r GGGTATCTAATCCCAGTTTG -t 8 -m 18

# denoise with dada2
scripts/dada2.R -p tele02 -l lib3


### FOR LIB4

# prep
scripts/prepare.sh -p tele02 -l lib4

# make symlinks HOME
ln -s ~/Projects/SeaDNA/temp-local-only/fastq/reduced/lib4/R1.fastq.gz ~/Projects/genner-lab/meta-fish-pipe/temp/processing/tele02-lib4/fastq/R1.fastq.gz
ln -s ~/Projects/SeaDNA/temp-local-only/fastq/reduced/lib4/R2.fastq.gz ~/Projects/genner-lab/meta-fish-pipe/temp/processing/tele02-lib4/fastq/R2.fastq.gz
# uni
ln -s /media/1TB/rc16041/Projects-temp-local-only/seadna-temp-local-only/fastq/reduced/lib4/R1.fastq.gz ~/Projects/genner-lab/meta-fish-pipe/temp/processing/tele02-lib4/fastq/R1.fastq.gz
ln -s /media/1TB/rc16041/Projects-temp-local-only/seadna-temp-local-only/fastq/reduced/lib4/R2.fastq.gz ~/Projects/genner-lab/meta-fish-pipe/temp/processing/tele02-lib4/fastq/R2.fastq.gz

# generate barcodes
scripts/generate-barcodes.R -p tele02 -l lib4 -f 18 -r 20

# demux
scripts/demultiplex.sh -p tele02 -l lib4 -f AAACTCGTGCCAGCCACC -r GGGTATCTAATCCCAGTTTG -t 8 -m 18

# denoise with dada2
scripts/dada2.R -p tele02 -l lib4




############## ALL LIBS - TAXONOMIC ASSIGNMENT ##############
############## ALL LIBS - TAXONOMIC ASSIGNMENT ##############

# run taxonomic assignment
scripts/taxonomic-assignment.sh -t 8 -p tele02









#### SUBSET THE CUSTOM REFS
source(here::here("scripts/funs-libs.R"))
ref.lib <- read_csv(here("assets/custom-reference-library.csv"),guess_max=99999,col_types=cols())
ref.lib.filt <- ref.lib %>% filter(!is.na(nucleotidesFrag.12s.taberlet.noprimers))
write.FASTA(tab2fas(df=ref.lib.filt,seqcol="nucleotidesFrag.12s.taberlet.noprimers",namecol="dbid"),file=here("assets/custom-reference-library.fasta"))




### refseq CHECK

# trim primers from the custom reference library
cutadapt -n 1 -e 0.3 -O 10 -g AAACTCGTGCCAGCCACC temp/reference-library/refseq-annotated.fasta | cutadapt --minimum-length 144 --maximum-length 216 -n 1 -e 0.3 -O 10 -a  CAAACTGGGATTAGATACCC -o temp/reference-library/refseq-annotated-trimmed.fasta -

# derep
vsearch --derep_fulllength temp/reference-library/refseq-annotated-trimmed.fasta --minuniquesize 1 --fasta_width 0 --output temp/reference-library/refseq-annotated-trimmed-derep.fasta

vsearch --threads 8 --sintax temp/reference-library/asvs-clean-cat-relabel-derep.fasta --db temp/reference-library/refseq-annotated.fasta --sintax_cutoff 0.7 --tabbedout temp/reference-library/sintax-output.tsv


renv::install("Matrix@1.3-2")


# my files
mkdir -p reduced/lib3 reduced/lib4
gzip -cd SeaDNA_Tele02_Lib03v2_R1.fastq.gz | head -n 400000 | gzip > reduced/lib3/R1.fastq.gz
gzip -cd SeaDNA_Tele02_Lib03v2_R2.fastq.gz | head -n 400000 | gzip > reduced/lib3/R2.fastq.gz
gzip -cd SeaDNA_Teleo02_Lib-04_S2_L001_R1_001.fastq.gz | head -n 400000 | gzip > reduced/lib4/R1.fastq.gz
gzip -cd SeaDNA_Teleo02_Lib-04_S2_L001_R2_001.fastq.gz | head -n 400000 | gzip > reduced/lib4/R2.fastq.gz



FWD="AAACTCGTGCCAGCCACC"
REV="GGGTATCTAATCCCAGTTTG"



