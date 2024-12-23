#!/bin/bash
#SBATCH -c 20                               # Request one core
#SBATCH -t 0-00:20                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem-per-cpu=5G                   # Memory per core in MiB
#SBATCH -o jobs/hostname_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e jobs/hostname_%j.err                 # File to which STDERR will be written, including job ID (%j)
                                           # You can change the filenames given with -o and -e to any filenames you'd like

################################################################################
#### NOTE: These lines are only necessary if submitting this script to a Slurm-
#### based job scheduler; specifically to access the sratoolkit commands.
#### If this script is being run directly in the terminal
#### `bash SRA_downnload_scripts/download_SRA_files.sh`, these can be commented
#### out as long as the agorbns environment has already been loaded.
module load miniconda3/23.1.0
source activate /home/sem689/.conda/envs/agorbns
################################################################################

python PreProcessReads/MakeReadFile.py $@