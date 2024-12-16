################################################################################
#MakeSiteTypeReadFiles.py
################################################################################
import imp # Used to import general.py
import time # Used for time
import numpy as np
import pandas as pd
import itertools as it
from more_itertools import unique_everseen
imp.load_source("general",
                ("/lab/bartel1_ata/mcgeary/computation/"
                 "AgoRBNS/general/general.py"
                )
               )
imp.load_source("RBNS_methods",
                ("/lab/bartel1_ata/mcgeary/computation/"
                 "AgoRBNS/general/RBNS_methods.py"
                )
               )
from general import *
from RBNS_methods import *
from sitetypes import get_seq_site_map

# FUNCTIONS
def count_read_kmers(read_seqs, _mirna, experiment, n_constant,
                     rand_length):
    time_start = time.time()
    mm_sites = []
    # Assign the strings for the 8mer and 6mer sites
    site_8mer = _mirna["8mer"]
    # Make a list of the nucleotides of the 8mer site.
    site_8mer_list = list(site_8mer)
    # Make a list of all possible single-nucleotide mismatches to the 8mer site.
    single_mm = [
        site_8mer[:i] + i_n + site_8mer[i+1:]
        for i in range(0, 8)
        for i_n in ["A", "C", "G", "T"] if i_n != site_8mer[i]
    ]
    # Make a list of all possible double mismatches to the 8mer site.
    double_mm = [
        site_8mer[:i] + i_n + site_8mer[i+1:j] + j_n + site_8mer[j+1:]
        for i in range(0, 7) for j in range(i + 1, 8)
        for i_n in ["A", "C", "G", "T"] if i_n != site_8mer[i]
        for j_n in ["A", "C", "G", "T"] if j_n != site_8mer[j]
    ]
    # Make the full list removing any redundant sites.
    full_list = [site_8mer] + single_mm + list(unique_everseen(double_mm)) + ["None"]
    site_dict = {i : 0 for i in full_list}
    none_kmers_dict = collections.defaultdict(int)
    print(none_kmers_dict)
    sys.stdout.flush()

    # Define the region where the programmed mismatch site should be within the
    # read.    
    for i_r, r in enumerate(read_seqs):
        programmed_kmer = r.strip()[25:25 + 8]
        if programmed_kmer in site_dict.keys():
            site_dict[programmed_kmer] += 1
        else:
            site_dict["None"] += 1
            none_kmers_dict[programmed_kmer] += 1
    sys.stdout.flush()
    site_dict_df = pd.DataFrame.from_dict(
        site_dict, orient="index"
    ).reindex(full_list)
    none_dict_df = pd.DataFrame.from_dict(
        none_kmers_dict, orient="index"
    )
    return [site_dict_df, none_dict_df]

def main():
    time_start = time.time()
    # Explanation of arguments:
    # n_constant : The number of constant sequence nucleotides appended to
    #   either side of the random/random+programmed portion of the library.
    # seedex_binary: 
    arguments = ["miRNA", "experiment", "condition", "n_constant", "-jobs",
                 "-temp_path", "-test_binary"]
    args = parse_arguments(arguments)
    (mirna, experiment, condition,
     n_constant, jobs, temp_path, test) = args
    # Determine length of random region.
    if (("let-7a" in mirna or "miR-1" in mirna) and
        mirna != "miR-155_let-7a" and
        experiment in ["equil_c_nb", "equil_s_nb", "equil_c_alt_nb",
                       "equil_sc_nb", "equil_sc_alt_nb", "equil_c2_nb",
                       "equil_c2_alt_nb"]):
        rand_length = 38
    else:
        rand_length = 37

    _mirna = Mirna(mirna)
    read_dir = "reads"
    if not jobs:
        jobs = 20
    args  = [_mirna, experiment, int(n_constant), int(rand_length)]
    ############################################################################
    ## PART 2: Set up the multiprocessing. #####################################
    reads_path = get_analysis_path(mirna, experiment, condition, read_dir)
    # Iterate through the read file via multiprocessing.
    threads = multiproc_file(reads_path, int(jobs), count_read_kmers, test,
                              *args)
    # Make a dataframe where each column is the results of one of the threads,
    # then sum across the rows to make the final, single column off output.
    output = pd.concat([thread[0] for thread in threads], axis=1).sum(axis=1)
    # Also make a second dataframe combining all of the kmers that were not one
    # or two mutations away from being a canonical 6mer site.
    none_df = pd.concat([thread[1] for thread in threads], axis=1).sum(axis=1)
    ## MAKE THE OUTPUT FILE ####################################################
    # Format the output path:    
    extension = "_%s" %n_constant
    if test:
        extension = "%s_test" %extension    
    if temp_path:
        kmers_path = "I_combined_temp/%s" %temp_path
    else:
        kmers_path = get_analysis_path(mirna, experiment, condition,
                                       "kmers_positional_programmed_only",
                                       ext=extension)
        none_path = get_analysis_path(mirna, experiment, condition,
                                       "kmers_positional_programmed_only_unexpected",
                                       ext=extension)
    # Print the path to allow for easy checking that the script worked, then
    # write the output file and the script is done.
    print(kmers_path)
    print(output)
    print(none_df)
    output.to_csv(kmers_path, index_label=False, sep="\t")
    none_df.to_csv(none_path, index_label=False, sep="\t")
    print("finished writing kmer file.")
    print_time_elapsed(time_start)
################################################################################

if __name__ == "__main__":
    main()

