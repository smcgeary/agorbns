################################################################################
#MakeSiteTypeReadFiles.py
################################################################################
import imp # Used to import general.py
import time # Used for time
import numpy as np
import pandas as pd
import os
import subprocess
imp.load_source("general",
                ("/lab/bartel1_ata/mcgeary/computation/"
                 "AgoRBNS/general/general.py"
                )
               )
from general import parse_arguments, print_time_elapsed, get_analysis_path, multiprocess_file, multiprocess_test, readline, get_kmer_list
from sitetypes import get_seq_site_map
# This script reads one read file and outputs multiple .txt files to be used in all
# downstream applications. Inputs are:

# mirna: The miRNA in the experiment to be processed.
# experiment: The type of experiment to be processed.
# condition: The sampel type within the experiment to be processed.

# Functions from general this uses:


# Constants

# FUNCTIONS

def assign_kmer_to_read(read_seqs, kmer_weights, n_constant, read_length, k_length, mirna, experiment):
    """Takes a read sequence and identifies the position of all site types'.

    Args:
        read_sequence: The read sequence.
        sites: A list of sequences corresponding to each site type.

    Returns:
        A tuple containing the list of reads and a dictionary of read types.
    """

    # Sub-function to be used on each site:
    kmersXcounts = {barcode : {kmer : 0 for kmer in get_kmer_list(k_length)} for barcode in ["TGT", "ACA"]}
    # ticksite = 0
    # tick = 0
    for i, read in enumerate(read_seqs):
        read = read.strip()
        # print(read)
        # Define the barcode portion of the read.
        barcode = read[26 + read_length : 26 + read_length + 3]

        # Deals with the TCG instead of TGT ending for miR-1 equilibrium exp
        if barcode == "TCG":
            barcode = "TGT"
            if mirna != "miR-1" or experiment != "equilibrium":
                read = read[:26 + 37]+"TGTTCGTATGCCGTCTTCTGCTTG"
        if barcode == "TGT" and mirna == "miR-1" and experiment == "equilibrium":
            read = read[:26 + 37]+"TCGTATGCCGTCTTCTGCTTG"


        # Take the constant regiong + the number of constant sequences
        read_seq = read.strip()[(26 - n_constant) :
                                (26 + read_length + n_constant)]
        # Find all sites within the read
        kmers = [read_seq[i:i+k_length] for i in range(read_length +2*n_constant - k_length + 1)]
        weights = [kmer_weights[kmer] for kmer in kmers]
        total = sum(weights)
        # if "AACATTCC" in kmers:
        #     # print(kmersXcounts["AACATTCC"])
        #     # print(total)
        #     print(tick)
        #     print(ticksite)
        #     print(read_seq)
        #     ticksite +=1
        for i, kmer in enumerate(kmers):
            kmersXcounts[barcode][kmer] += float(weights[i])/total
        # tick += 1
    return kmersXcounts

def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["miRNA", "experiment", "condition", "n_constant", "kmer_length", "n_iter"]
    mirna, experiment, condition, n_constant, kmer_length, n_iter = parse_arguments(arguments)
    motifs_path = get_analysis_path(mirna, experiment, "motifs_PAPER", "sitekmers")
    print(motifs_path)
    # Load file with site types for the miRNA.
    if experiment in ["equil_mmseed_nb", "equil_seed_nb"]:
        read_length = 38
    else:
        read_length = 37
    n_iter = int(n_iter)
    if n_iter > 0:    

        new_file = get_analysis_path(mirna, experiment, condition, "full_reads", ext = "_iter%s" %(n_iter))

        if not os.path.exists(new_file):
            if n_iter == 1:
                old_file = get_analysis_path(mirna, experiment, condition, "full_reads")
            else:
                old_file = get_analysis_path(mirna, experiment, condition, "full_reads", ext = "_iter%s" %(n_iter-1))
            with open(motifs_path, "r") as motifs:
                motif = motifs.readlines()[n_iter].split("\t")[0]
            print(motif)
            with open(new_file,"w") as file_out:
                subprocess.call(["grep","-v",motif,old_file],stdout=file_out)
        extension = "_%s_%s_iter%s" %(n_constant, kmer_length, n_iter)
        print(new_file)

        reads_path = new_file
    else:
        extension = "_%s_%s" %(n_constant, kmer_length)
        reads_path = get_analysis_path(mirna, experiment, condition,
                                       "full_reads")
    kmer_counts_path = get_analysis_path(mirna, experiment, condition,
                                              "kmer_counts",
                                               ext=extension)
   
    kmer_list = get_kmer_list(int(kmer_length))

    kmer_weights = {kmer : 1 for kmer in kmer_list}
    with open(reads_path, "rb") as file_in:
        results = multiprocess_file(file_in,
                                    readline,
                                    int(1e6),
                                    assign_kmer_to_read,
                                    kmer_weights,
                                    int(n_constant),
                                    read_length,
                                    int(kmer_length),
                                    mirna,
                                    experiment)

    # Collect the output from each thread, which are the list
    # of reads and the dictionary of read types.


    # Flatten the list into one list (found on stackoverflow) and write
    # it to its output file.
    kmer_counts = {barcode : {kmer : sum([i[barcode][kmer] for i in results])
                    for kmer in kmer_list} for barcode in ["TGT", "ACA"]}
    # print(kmer_counts["AACATTCC"]*24)
    kmer_initial = kmer_counts
    # if condition not in ["I", "I_combined"]:
    #     for round in range(3):
    #         with open(reads_path, "rb") as file_in:
    #             results = multiprocess_test(file_in,
    #                                         readline,
    #                                         int(1e6),
    #                                         assign_kmer_to_read,
    #                                         1000000,
    #                                         kmer_counts,
    #                                         int(n_constant),
    #                                         read_length,
    #                                         int(kmer_length),
    #                                         mirna,
    #                                         experiment)
    #         kmer_counts = {kmer : sum([i[kmer] for i in results])
    #                            for kmer in kmer_list}

    #         kmer_final = {kmer : [kmer_initial[kmer], kmer_counts[kmer]] for kmer in kmer_list}


    #         kmers_total = pd.DataFrame.from_dict(kmer_final,orient="index").sort_index()
    #         kmers_total.columns = ["Initial","Counts"]
    #         print(kmers_total.sort_values(["Counts"],ascending=False)[:10])

    kmer_final = {kmer : [kmer_initial["TGT"][kmer], kmer_initial["ACA"][kmer]] for kmer in kmer_list}


    kmers_total = pd.DataFrame.from_dict(kmer_final,orient="index").sort_index()
    kmers_total.columns = ["TGT","ACA"]
    print(kmers_total[1:10])
    print(kmers_total.sort_values(["TGT"],ascending=False)[:10])

    kmers_total.to_csv(kmer_counts_path, sep = "\t")

    print(kmer_counts_path)

    # Print the amount of time the script took to complete.
    print_time_elapsed(time_start)

################################################################################

if __name__ == "__main__":
    main()

