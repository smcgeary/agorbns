################################################################################
#MakeSiteTypeReadFiles.py
################################################################################
import imp # Used to import general.py
import time # Used for time
import numpy as np
import pandas as pd
import itertools as it
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
def LCSuff(seq1, seq2):
    if len(seq1) == 0 or len(seq2) == 0:
        return("")
    elif seq1[-1] == seq2[-1]:
        return(LCSuff(seq1[:-1], seq2[:-1]) + seq1[-1])
    else:
        return("")

def LCSubString(seq1, seq2):
    pres1 = [seq1[:i] for i in range(1, len(seq1) + 1)]
    pres2 = [seq2[:i] for i in range(1, len(seq2) + 1)]
    pre_pairs = [(i, j) for i in pres1 for j in pres2]
    LCs = [LCSuff(i[0], i[1]) for i in pre_pairs]
    LCmax = max([len(i) for i in LCs])
    return [i for i in LCs if len(i) == LCmax]


def assign_flanks_to_read(read_seqs, _mirna, _sitelist, experiment,
                             n_constant, rand_length, buffer3p):
    time_start = time.time()
    for i_r, r in enumerate(read_seqs):
        # Get the read number:
        read = r.strip()
        # Make Read object.
        _read = Read(read, rand_length, _mirna, n_constant, experiment)
        _read.get_all_sites_in_read(_sitelist)
        add_sites_from_read(_sitelist, _read, buffer_=buffer3p)
        add_flanks_from_read(_sitelist, _read, buffer_=buffer3p)
    return _sitelist.top_counts_df(), _sitelist.flank_counts_df("TGT"), _sitelist.flank_counts_df("ACA")



def main():
    time_start = time.time()
    arguments = ["miRNA", "experiment", "condition", "n_constant", "sitelist",
                 "-buffer3p_binary", "-jobs", "-test_binary"]
    args = parse_arguments(arguments)
    (mirna, experiment, condition,
     n_constant,sitelist,buffer3p, jobs, test) = args
    _mirna = Mirna(mirna)

    if experiment in ["equil_mmseed_nb", "equil_seed_nb"]:
        read_length = 38
    else:
        read_length = 37    

    _sitelist = SiteList(_mirna, sitelist, int(read_length) + 2*int(n_constant))
    extension = "_%s_%s" %(n_constant, sitelist)

    if buffer3p:
        extension = "%s_buffer3p" %(extension)

    if condition == "I_combined":
        if "tp" in experiment:
            if mirna == "miR-7-24nt":
                input_list = [(mirna, experiment, condition)]
            else:
                input_list = INPUT_LIST_I_COMBINED_TP
        elif "miR-7" in mirna and experiment in ["equilibrium2_nb", "equilibrium3_nb"]:
            input_list = INPUT_LIST_I_COMBINED_NB
        else:
            input_list = INPUT_LIST_I_COMBINED
        if not jobs:
            jobs = 20
    elif condition == "0_combined":
        input_list = INPUT_LIST_0_COMBINED
        if not jobs:
            jobs = 20
    else:
        input_list = [(mirna, experiment, condition)]
        if not jobs:
            jobs = 20
    args  = [_mirna, _sitelist, experiment, int(n_constant),
                    read_length, buffer3p]
    threads = []
    for i_input in input_list:
        print(i_input)
        reads_path = get_analysis_path(i_input[0], i_input[1], i_input[2], "reads")
        if not jobs:
            jobs = 1
        threads += multiproc_file(reads_path, int(jobs), assign_flanks_to_read, test,
                                  *args)

    # for i_input in input_list:
    #     print(i_input)
    #     reads_path = get_analysis_path(i_input[0], i_input[1], i_input[2], "reads")
    #     with open(reads_path, "rb") as file_in:
    #         if test:
    #             total_reads = 10000
    #             process_size = total_reads/10
    #             threads += multiprocess_test(file_in,
    #                                         readline_one,
    #                                         process_size,
    #                                         assign_flanks_to_read,
    #                                         total_reads,
    #                                         _mirna,
    #                                         _sitelist,
    #                                         experiment,
    #                                         int(n_constant),
    #                                         read_length)
    #         else:
    #             threads += multiprocess_file(file_in,
    #                                         readline_one,
    #                                         JOB_SPLIT,
    #                                         assign_flanks_to_read,
    #                                         _mirna,
    #                                         _sitelist,
    #                                         experiment,
    #                                         int(n_constant),
    #                                         read_length)

    print("done reading files")
    flank_counts_path = get_analysis_path(mirna, experiment, condition,
                                              "flanks",
                                              ext=extension)

    site_counts = merge_data_frames([i[0] for i in threads])["TGT"]
    pulse_flank_counts = merge_data_frames([i[1] for i in threads])
    chase_flank_counts = merge_data_frames([i[2] for i in threads])

    print(pulse_flank_counts.sum(axis=0))
    print(chase_flank_counts.sum(axis=0))
    print(flank_counts_path)

    if test:
        print(pulse_flank_counts)
    if not test:
        pulse_flank_counts.to_csv(flank_counts_path, sep="\t", header=True)
    # print_time_elapsed(time_start)
################################################################################

if __name__ == "__main__":
    main()

