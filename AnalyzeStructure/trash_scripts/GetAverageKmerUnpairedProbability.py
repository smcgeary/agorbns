################################################################################
#MakeSiteTypeReadFiles.py
################################################################################
import imp # Used to import general.py
import time # Used for time
import math
import csv
import itertools as it
import multiprocessing
import sys
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
imp.load_source("general",
                ("/lab/bartel1_ata/mcgeary/computation/"
                 "AgoRBNS/general/general.py"
                )
               )
from general import parse_arguments, print_time_elapsed, get_analysis_path, multiprocess_test, readline, readline_two, get_kmer_list
from general import *
# This script reads one read file and outputs multiple .txt files to be used in all
# downstream applications. Inputs are:

# mirna: The miRNA in the experiment to be processed.
# experiment: The type of experiment to be processed.
# condition: The sampel type within the experiment to be processed.

def get_kmer_probs(read_seqs, k, left_ext=0, right_ext=0, ind_start=0, ind_stop=0):
    """Takes a read sequence and identifies the position of all site types'.

    Args:
        read_sequence: The read sequence.
        sites: A list of sequences corresponding to each site type.

    Returns:
        A tuple containing the list of reads and a dictionary of read types.
    """

    # Sub-function to be used on each site:
    kmer_range = range(26-ind_start+left_ext,26+37+ind_stop+1-right_ext-k)
    # Defines the centered x coordinates for each tiling window
    # Defines the range of possible window positions relative to the
    kmers = get_kmer_list(k)
    kmer_counts_map = {kmer: [0, 0, 0] for kmer in kmers}
    for count,seq in enumerate(read_seqs):
        if count%100000 == 0:
            print(count)
        read, structure = (j.strip() for j in seq)
        for ind in kmer_range:
            kmer = read[ind:ind+k]
            logp = sum([math.log10(1-float(i)) for i in structure.split("\t")[(ind-left_ext):(ind+k+right_ext)]])
            p = 10**logp
            n = kmer_counts_map[kmer][2]
            kmer_counts_map[kmer][0] = (kmer_counts_map[kmer][0] * n + p)/float(n+1)
            kmer_counts_map[kmer][1] = (kmer_counts_map[kmer][1] * n + logp)/float(n+1)
            kmer_counts_map[kmer][2] += 1


    return kmer_counts_map


def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["miRNA", "experient", "condition", "k", "left", "right"]
    mirna, experiment, condition, k, left, right = parse_arguments(arguments)

    # Get the path to the read file and to that of where the site labels will
    # be written.
    full_reads_path = get_analysis_path(mirna, experiment, condition, "full_reads")
    structure_pairs_path = get_analysis_path(mirna,experiment,condition,"structures_bp_prob")
    kmer_unpaired_prob_path = get_analysis_path(mirna,experiment,condition,"average_kmer_unpaired_prob/%smer_%s-l_%s-r/" %(k, left, right))

    kmers = get_kmer_list(int(k))

    with open(full_reads_path,"rb") as file_in_reads:
        with open(structure_pairs_path, "rb") as file_in_structures:
            results = multiprocess_file([file_in_reads, file_in_structures],
                                        readline_two, 1000000, 
                                        get_kmer_probs,
                                        int(k),
                                        int(left),
                                        int(right))

    kmers = get_kmer_list(k)
    kmer_totals_map = {kmer: sum([i[kmer][2] for i in results]) for kmer in kmers}
    kmer_p_map = {kmer: sum([i[kmer][0]*i[kmer][2] for i in results])/kmer_totals_map[kmer] for kmer in kmers}
    kmer_logp_map = {kmer: sum([i[kmer][1]*i[kmer][2] for i in results])/kmer_totals_map[kmer] for kmer in kmers}

    total_kmers = sum(kmer_totals_map.values())

    with open(kmer_unpaired_prob_path,"wb") as file_out:
        file_out.write("\t".join(["kmer", "prob_unpaired", "log_prob_unpaired", "percentage"] + ["\n"]))
        file_out.write("".join(["\t".join([kmer,
                                str(kmer_p_map[kmer]),
                                str(kmer_logp_map[kmer]),
                                str(kmer_totals_map[kmer]/float(total_kmers)),
                                "\n"]) for kmer in kmers]))

    # Print the amount of time the script took to complete.
    print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()

