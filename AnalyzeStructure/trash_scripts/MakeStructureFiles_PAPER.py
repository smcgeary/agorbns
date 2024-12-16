################################################################################
#MakeSiteTypeReadFiles.py
################################################################################
import imp # Used to import general.py
import time # Used for time
import itertools as it
import sys
import math
imp.load_source("general",
                ("/lab/bartel1_ata/mcgeary/computation/"
                 "AgoRBNS/general/general.py"
                )
               )
from general import parse_arguments, print_time_elapsed, get_analysis_path, multiprocess_file, multiprocess_test, readline
import RNA
import numpy as np
import pandas as pd
import numpy.matlib
np.set_printoptions(suppress=True,linewidth=np.nan,threshold=np.nan)
# This script reads one read file and outputs multiple .txt files to be used in all
# downstream applications. Inputs are:

# mirna: The miRNA in the experiment to be processed.
# experiment: The type of experiment to be processed.
# condition: The sampel type within the experiment to be processed.

# Functions from general this uses:


# Constants

# FUNCTIONS


def analyze_structure(read_seqs, n_constant, read_length, mirna, experiment):
    """Takes a read sequence and identifies the position of all site types'.

    Args:
        read_sequence: The read sequence.
        sites: A list of sequences corresponding to each site type.

    Returns:
        A tuple containing the list of reads and a dictionary of read types.
    """

    # Sub-function to be used on each site:

    # Initialize count dictionary with keys as site names and a value of
    # for each item.
    # count_site_map = {value: {"TGT" : 0, "ACA" : 0}
    #                   for value in site_seq_map.values()}
    # count_multisite_map = {value: {"TGT" : 0, "ACA" : 0}
    #                        for value in site_seq_map.values()}
    # count_site_map.update({"None": {"TGT" : 0, "ACA" : 0}})
    # count_multisite_map.update({"None": {"TGT" : 0, "ACA" : 0}})

    
    # Pre-allocate output with "None" for each line.
    if n_constant == "Full":
        pairs = list(it.combinations(range(len(read_seqs[0].strip())), 2))
        space = 0
        def CutRead(read_seq):
            return read.strip()
    else:
        print(read_length+2*int(n_constant))
        space = 26 - int(n_constant)
        pairs = list(it.combinations(range(read_length+2*int(n_constant)), 2))
        def CutRead(read_seq):
            read_cut = read.strip()[(26 - int(n_constant)) :
                                    (26 + read_length + int(n_constant))]
            return read_cut

    mfe_output = [0]*len(read_seqs)
    structure_output = ["NA"]*len(read_seqs)
    entropy_output = ["NA"] * len(read_seqs)
    time_start = time.time()
    for i, read in enumerate(read_seqs):
        if i % 100000 == 0:
            print(i)
            print_time_elapsed(time_start)
            sys.stdout.flush()
        # Define the barcode portion of the read.
        # Take the constant regiong + the number of constant sequences
        time_old = time.time()
        # Define the barcode portion of the read.
        # Take the constant regiong + the number of constant sequences
        read_seq = CutRead(read)
        mfe_output[i] = RNA.fold(read_seq)[1]
        pf = RNA.pf_fold(read_seq)[1]
        probs = [1]*len(read_seq)
        entropy = [0]*len(read_seq)
        for pair in pairs:
            [p1, p2] = [int(p) for p in pair]
            prob = float(RNA.get_pr(p1+1, p2+1))
            if prob > 0:
                probs[p1] -= prob
                probs[p2] -= prob
                entropy[p1] -= prob*math.log(prob)
                entropy[p2] -= prob*math.log(prob)
        for k, prob_free in enumerate(probs):
            entropy[k] -= prob_free*math.log(prob_free)
        structure_output[i] = "\t".join([str(j) for j in probs])
        entropy_output[i] = "\t".join([str(j) for j in entropy])


    return mfe_output, structure_output, entropy_output

        # Find all sites within the read
    #     site_maps = get_unique_sites(read_seq, site_seq_map)
    #     # if site_maps:
    #     #     print(site_maps)
    #     #     site_maps = [j for j in site_maps]
    #     #     print(site_maps)
    #     if site_maps:
    #         output_temp = ", ".join(["%s:%s-%s" % (site, index, end)
    #                                for site, index, end in site_maps])
    #         output[i] = output_temp
    #         for site_map in site_maps:
    #             count_site_map[site_map[0]][barcode] += 1
    #         if len(site_maps) > 1:
    #             key = ",".join(sorted([j[0] for j in site_maps]))
    #             if key in count_multisite_map:
    #                 count_multisite_map[key][barcode] += 1
    #             else:
    #                 count_multisite_map[key] = {"TGT" : 0, "ACA" : 0}
    #                 count_multisite_map[key][barcode] += 1
    #         else:
    #             count_multisite_map[site_maps[0][0]][barcode] += 1
    #     else:
    #         count_site_map["None"][barcode] += 1
    #         count_multisite_map["None"][barcode] += 1
    # return output, count_site_map, count_multisite_map

def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["miRNA", "experiment", "condition", "n_constant"]
    mirna, experiment, condition, n_constant = parse_arguments(arguments)

    # Load file with site types for the miRNA.

    # Initialize the dictionary with all sites, and two more dictionaries,
    # one mapping to the site name from the sequence, and one from site name
    # to the count.
    # Get the path to the read file and to that of where the site labels will
    # be written.
    if experiment in ["equil_mmseed_nb", "equil_seed_nb"]:
        read_length = 38
    else:
        read_length = 37    

    extension = "_%s" %(n_constant)

    reads_path = get_analysis_path(mirna, experiment, condition,
                                   "full_reads")
    bpprob_path = get_analysis_path(mirna, experiment, condition,
                                   "pairing_probabilities",ext=extension)
    posentropy_path = get_analysis_path(mirna, experiment, condition,
                                   "positional_entropy",ext=extension)

    mfe_path = get_analysis_path(mirna, experiment, condition,
                                         "mfe_files", ext=extension)

    with open(reads_path, "rb") as file_in:
        results = multiprocess_file(file_in,
                                    readline,
                                    int(10000),
                                    analyze_structure,
                                    n_constant,
                                    read_length,
                                    mirna,
                                    experiment)


    with open(bpprob_path, "wb") as bp_out:
        with open(mfe_path, "wb") as mfe_out:
            with open(posentropy_path, "wb") as entropy_out:
                for result in results:
                    for i,mfe in enumerate(result[0]):
                        mfe_out.write("%s\n" %(mfe))
                        bp_out.write("%s\n" %(result[1][i]))
                        entropy_out.write("%s\n" %(result[2][i]))


    print(bpprob_path)
    print(mfe_path)
    print(posentropy_path)
    print_time_elapsed(time_start)

################################################################################

if __name__ == "__main__":
    main()

