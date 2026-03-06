#!/bin/bash
#SBATCH -c 15                              # Request one core
#SBATCH -t 0-00:06                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem-per-cpu=300M                   # Memory per core in MiB
#SBATCH -o /n/groups/klein/mcgeary/gradschool/AgoRBNS/SRUN/out/%j.out  # File to which STDOUT will be written, with job ID (%j)
#SBATCH -e /n/groups/klein/mcgeary/gradschool/AgoRBNS/SRUN/err/%j.err  # File to which STDERR will be written, with job ID (%j)
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

python AssignSiteTypes/AssignMiRNAs.py $@