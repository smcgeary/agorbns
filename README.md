# AgoRBNS Analysis Pipeline

This repository contains analysis code for two publications:

1. McGeary SE, et al. "The biochemical basis of microRNA targeting efficacy." Science. 2019.
   [Link to paper](https://www.science.org/doi/10.1126/science.aav1741)

2. McGeary SE, et al. "MicroRNA 3'-compensatory pairing occurs through two binding modes, with affinity shaped by nucleotide identity and position." eLife. 2022.
   [Link to paper](https://elifesciences.org/articles/73188)

## Data Access

Raw sequencing data from 1. is available from GEO under accession number GSE140220.
Raw sequencing data from 2. is available from GEO under accession number GSE196458.



### Downloading the Data

#### Option 1: Using SRA toolkit
```bash
# Download individual samples
fastq-dump --split-files SRR10427205  # Replace with specific SRR number

# Or download all samples (24 total)
for i in {205..228}; do
    fastq-dump --split-files SRR104272$i
done