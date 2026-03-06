#!/bin/bash
#SBATCH -c 20                                             # Request twenty cores
#SBATCH -t 0-00:05                                   # Runtime in D-HH:MM format
#SBATCH -p short                                           # Partition to run in
#SBATCH --mem-per-cpu=1G                                # Memory per core in MiB
#SBATCH -o /n/groups/klein/mcgeary/gradschool/AgoRBNS/SRUN/out/%j.out  # File to which STDOUT will be written, with job ID (%j)
#SBATCH -e /n/groups/klein/mcgeary/gradschool/AgoRBNS/SRUN/err/%j.err  # File to which STDERR will be written, with job ID (%j)

################################################################################
#### NOTE: These lines are compatible with the O2 cluster at Harvard Medical
#### school, which uses the SLURM scheduler, and assume the script is being
#### called from the makefile with the `sbatch [script args]` command.
#### If this script is being run directly in the terminal using 
#### `bash [script args]`, these lines can be commented out, and the agorbns
#### conda environment should be activated.
module load miniconda3/23.1.0
source activate /home/sem689/.conda/envs/agorbns
################################################################################

python PreProcessReads/MakeMiRNAReadFile.py $@