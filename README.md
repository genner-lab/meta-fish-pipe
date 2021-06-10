[![DOI](https://zenodo.org/badge/xxx.svg)](https://zenodo.org/badge/latestdoi/xxx)

# meta-fish-pipe

A bioinformatics pipeline for fish metabarcoding. Currently supported 12S metabarcode markers are: 'tele02', 'mifish-u', 'elas02', and 'mifish-u-mod'.

### Setup (do once)

1. Install the following software on your system, and make them available on your $PATH: [cutadapt](https://github.com/marcelm/cutadapt) v3.4, [vsearch](https://github.com/torognes/vsearch) v2.17.0, [seqkit](https://github.com/shenwei356/seqkit) v0.16.1, [raxml-ng](https://github.com/amkozlov/raxml-ng) v1.0.2, [epa-ng](https://github.com/Pbdas/epa-ng) v0.3.8, [gappa](https://github.com/lczech/gappa) v0.7.1, [hmmer](http://hmmer.org/) v3.1b2, [ncbi-blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) v2.11.0, [mafft](http://mafft.cbrc.jp/alignment/software/) v7.271.

2. Clone the repository:

```
git clone https://github.com/genner-lab/meta-fish-pipe.git
```

3. Change directory:

```
cd meta-fish-pipe
```

4. Obtain correct R packages:

```
Rscript -e "renv::restore()"
```

### Add user data and setup directories (per project)

Users must provide 4 files for the project. They can be named however you like, but must be added to the `assets/` directory before running the pipeline. This files only need to be provided once at the beginning of the project.

1. Your custom 'curated fish' reference library in CSV format, which can be obtained from [https://github.com/genner-lab/meta-fish-lib](https://github.com/genner-lab/meta-fish-lib).

2. Your ReqSeq 'all taxa' reference library in CSV format, which can be obtained from [https://github.com/genner-lab/refseq-reflib](https://github.com/genner-lab/refseq-reflib).

3. A 'samples' file containing your sample/primer/tag information in CSV format. There must be 8 named columns, as shown in this example file: [samples-examples.csv](assets/examples/samples-examples.csv).

4. A file in CSV format containing DNA sequences that are known or suspected to be contaminants from other experiments in the lab. There must be 3 named columns, as shown in this example file: [contaminants-examples.csv](assets/examples/contaminants-examples.csv).

### Record project session and prepare project directories 

This script records the project session info (installed software and R package versions) and also sets up the project directories. Paths to the user submitted reference libraries need to be provided as follows. Script only needs to be run once at the beginning of the project.

- The '-r' flag is the path to the ReqSeq 'all taxa' reference library.

- The '-c' flag is the path to the custom 'curated fish' reference library.

```
scripts/session-info.sh -r assets/refseq206-annotated-tele02.csv -c assets/meta-fish-lib-v243.csv
```

### Prepare library directories and files

This script sets up your working area directories for each library.

- The '-p' flag is the primer set. Must be one of 'tele02', 'mifish-u', 'elas02', 'mifish-u-mod'. It must also match the 'primerSet' field in your provided samples file.

- The '-l' flag is the name of your library. It must also match the 'library' field in your provided samples file.

```
# prepare working area directories
scripts/prepare-libraries.sh -p tele02 -l lib1
```

- Now create a link to your R1 and R2 fastq read files from their current location on your computer into the location required for the pipeline (for each library). It is important that 'tele02' and 'lib1' paths match the flags used for `prepare.sh` and that the files are simply named 'R1.fastq.gz' and 'R2.fastq.gz'. 

```
# create symbolic link
ln -s /your/file/path/fileR1.fastq.gz temp/processing/tele02-lib1/fastq/R1.fastq.gz
ln -s /your/file/path/fileR2.fastq.gz temp/processing/tele02-lib1/fastq/R2.fastq.gz

```


### Generate barcodes (for each library)

This script creates the sample barcode tags for demultiplexing, using information in your sample sheet (for each library).

- The '-p' flag is the primer set. Must be one of 'tele02', 'mifish-u', 'elas02', 'mifish-u-mod'. It must also match the 'primerSet' field in your provided samples file.

- The '-l' flag is the name of your library. It must also match the 'library' field in your provided samples file.

- The '-f' flag is the length of your forward PCR primer in base pairs (bp).

- The '-r' flag is the length of your reverse PCR primer (bp).

- The '-m' flag is the path to your samples file in `assets/`.

```
# generate barcodes
scripts/generate-barcodes.R -p tele02 -l lib1 -f 18 -r 20 -m assets/sequencing-master-may2021.csv
```

### Demultiplex (for each library)

This script finds sequences with attached primers, demultiplexes the samples using the sample tags for each library, and then trims off the primers.

- The '-p' flag is the primer set. Must be one of 'tele02', 'mifish-u', 'elas02', 'mifish-u-mod'.

- The '-l' flag is the name of your library.

- The '-f' flag is the forward PCR primer sequence.

- The '-r' flag is the reverse PCR primer sequence.

- The '-t' flag is the number of processing threads.

- The '-m' flag is the length (bp) of the shorter primer (or length of both if same).

```
# demultiplex
scripts/demultiplex.sh -p tele02 -l lib1 -f AAACTCGTGCCAGCCACC -r GGGTATCTAATCCCAGTTTG -t 8 -m 18
```

### Denoise (for each library)

This script runs the quality trimming, dada2 denoising, dereplication, and chimaera removal. The resulting amplicon sequence variants (ASVs) are then cleaned with hidden Markov models to remove the non-homologous sequences.

- The '-p' flag is the primer set. Must be one of 'tele02', 'mifish-u', 'elas02', 'mifish-u-mod'.

- The '-l' flag is the name of your library.

```
# denoise with dada2
scripts/dada2.R -p tele02 -l lib1
```

### Generate stats (for each library)

This script generates reads numbers at each stage of the pipeline. 

- The '-p' flag is the primer set. Must be one of 'tele02', 'mifish-u', 'elas02', 'mifish-u-mod'.

- The '-l' flag is the name of your library.

- The '-t' flag is the number of processing threads.

```
# generate stats
scripts/generate-stats.sh -p tele02 -l lib1 -t 8
```

### Taxonomic assignment (once for all libraries)

Once all libraries have been processed, this script performs taxonomic assignment using: (a) sintax; (b) blastn; and (c) epa, using the ReqSeq 'all taxa' reference library and custom 'curated fish' reference library.

- The '-p' flag is the primer set. Must be one of 'tele02', 'mifish-u', 'elas02', 'mifish-u-mod'.

- The '-t' flag is the number of processing threads.

```
# run taxonomic assignment
scripts/taxonomic-assignment.sh -t 8 -p tele02
```

### Assemble results

This script assembles all assignments for all libraries, and generates results tables and writes them into `results/`.

- The '-c' flag is the path to your contaminants file in `assets/`.

```
# assemble results
scripts/assemble-results.R -c assets/contaminants-exclude-may2021.csv
```


### Clean up (optional)

This script deletes the `temp/` directory to save disk space. ONLY RUN THIS IF YOU ARE ABSOLUTELY SURE YOU DO NOT STILL NEED THE INTERMEDIATE FILES AND DATA IN THERE. 

```
# clean up
rm -r temp
```
