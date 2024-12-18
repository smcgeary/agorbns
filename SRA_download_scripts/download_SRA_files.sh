#!/bin/bash
#SBATCH -c 1                               # Request one core
#SBATCH -t 3-00:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                          # Partition to run in
#SBATCH --mem=1G                           # Memory total in MiB (for all cores)
#SBATCH -o jobs/hostname_%j_%t.out            # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e jobs/hostname_%j_%t.err            # File to which STDERR will be written, including job ID (%j)
# This script is for downloading the fastq files associated my graduate work.
# The script should go through the entire list of files and their associated
# metadata.

# Testing to make sure that the access of metadata is working properly.

file=SRA_download_scripts/SRA_and_metadata.txt

echo $file

lines=$(cut -f 1 $file)
# Step 1. Load the sratoolkit.

# Step 2. Iterate over the SRR lines and prefetch and then download the
# associated files.
for LINE in $lines
do
    echo "'$LINE'"
    time_0=$(date +%s)
    if [ ! -f "data/raw/sra/$LINE.sra" ]
    then
        echo "Prefetching data."
        prefetch $LINE -O data/raw/sra/$LINE
    else
        echo "Already prefetched."
    fi
    time_1=$(date +%s)
    elapsed=$(( time_1 - time_0 ))
    echo "Prefetch time: $elapsed s."
    if [ ! -f "data/raw/fastq/$LINE.fq.gz" ]
    then
        echo "Performing fasterq-dump."
        # The "-f" option enables the faster-q dump to overwrite any temporary
        # files that might have been made in a previous call to the script.
        fasterq-dump -f -O data/raw/fastq/ data/raw/sra/$LINE
        time_2=$(date +%s)
        elapsed=$(( time_2 - time_1 ))
        echo "Fasterq-dump time: $elapsed s."
        mv data/raw/fastq/$LINE.fastq data/raw/fastq/$LINE.fq
        gzip data/raw/fastq/$LINE.fq
        time_3=$(date +%s)
        elapsed=$(( time_3 - time_2 ))
        echo "Rename and gzip time: $elapsed s."
    else
        echo "Already downloaded and gzipped file."
    fi
    echo ""
done
# IFS="$OLDIFS"