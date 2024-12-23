# AgoRBNS Analysis Pipeline

This repository contains analysis code for two publications:

1. McGeary SE, et al. "The biochemical basis of microRNA targeting efficacy." Science. 2019.
   [Link to paper](https://www.science.org/doi/10.1126/science.aav1741)

2. McGeary SE, et al. "MicroRNA 3'-compensatory pairing occurs through two binding modes, with affinity shaped by nucleotide identity and position." eLife. 2022.
   [Link to paper](https://elifesciences.org/articles/73188)

## Data Access

Raw sequencing data is available from GEO under the following accession numbers:
- GSE140220 (Science 2019 paper)
- GSE196458 (eLife 2022 paper)

## Requirements

- Conda (Miniconda or Anaconda) (tested using Miniconda 23.10)
- Git

## Installation and Setup

1. Clone this repository:
   ```bash
   git clone https://github.com/smcgeary/agorbns.git
   cd AGORBNS_working
   ```

2. Create the conda environment:
   ```bash
   conda env create -f environment.yml
   ```

3. Activate the conda environment:
   ```bash
   conda activate agorbns
   ```

4. Verify the environment is active:
   ```bash
   conda env list
   ```
   The active environment (agorbns) should be marked with *

5. When finished working, deactivate the environment:
   ```bash
   conda deactivate
   ```

Note: If the environment.yml file is updated, you can update your environment with:
```bash
   conda env update -f environment.yml
```

## Data acquisition

To download all data associated with both papers:
1. First assemble the metadata-table required for parsing the data:
```bash
   python SRA_download_scripts/make_SRA_and_sample_table.py 
```

2. Next download each dataset:
```bash
   bash SRA_download_scripts/download_SRA_files.sh 
```

## Repository Structure

Analyses proceed through three primary stages. First, the read data is extracted
using the functionality in the PreProcessReads directory. Then site types are
assigned to reads using the AssignSiteTypes directory. Then, biochemical
modeling that enables Kd estimation is performed using functionality in the
SolveForKds directory.


## Usage

The basic commands to perform from the command line are listed in order leading
up to each Figure sub-panel in the Figures.R script within the `general` directory.
This file is run by entering R from the terminal and loading it by calling:
```R
source("general/general.R")
```

## Citation

If you use this code in your research, please cite:

McGeary SE, et al. "The biochemical basis of microRNA targeting efficacy." Science. 2019.
McGeary SE, et al. "MicroRNA 3'-compensatory pairing occurs through two binding modes, with affinity shaped by nucleotide identity and position." eLife. 2022.