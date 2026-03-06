#!/bin/bash
#SBATCH -c 20                                             # Request twenty cores
#SBATCH -t 0-00:20                                   # Runtime in D-HH:MM format
#SBATCH -p short                                           # Partition to run in
#SBATCH --mem-per-cpu=10G                               # Memory per core in MiB
#SBATCH -o /n/groups/klein/mcgeary/gradschool/AgoRBNS/SRUN/out/%j.out  # File to which STDOUT will be written, with job ID (%j)
#SBATCH -e /n/groups/klein/mcgeary/gradschool/AgoRBNS/SRUN/err/%j.err  # File to which STDERR will be written, with job ID (%j)

################################################################################
#### NOTE: These lines are necessary if submitting this script to a Slurm-
#### based job scheduler. In practice, submitting this while having already
#### activated the agorbns environment causes the script to fail.
module load miniconda3/23.1.0
source activate /home/sem689/.conda/envs/agorbns
################################################################################

which python
python --version
echo $PATH

python PreProcessReads/MakeReadFile.py $@