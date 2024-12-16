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
def count_read_kmers(read_seqs, kmer_len):
    kmer_dict = {kmer: 0 for kmer in get_kmer_list(kmer_len)}
    time_start = time.time()
    sys.stdout.flush()
    read_len = len(read_seqs[0].strip())
    for i_r, r in enumerate(read_seqs):
        read = r.strip()
        sys.stdout.flush()
        kmers = [read[i:i+kmer_len] for i in range(read_len - kmer_len + 1)]
        for kmer in kmers:
            kmer_dict[kmer] += 1
    return kmer_dict

def main():
    time_start = time.time()
    arguments = ["rbp", "condition", "-kmer_len", "-jobs", "-test_binary"]
    args = parse_arguments(arguments)
    (rbp, condition, kmer_len, jobs, test) = args

    
    if not kmer_len:
        kmer_len = 5
    if not jobs:
        jobs = 20

    kmers = get_kmer_list(kmer_len)
    args = [kmer_len]
    ############################################################################
    reads_path = get_analysis_path_burge(rbp, "equilibrium", condition, "reads")
    kmers_path = get_analysis_path_burge(rbp, "equilibrium", condition, "kmers",
                                         ext = "_%s-mers" %(kmer_len))

    threads = multiproc_file(reads_path, int(jobs), count_read_kmers, test,
                             *args)

    kmer_dict_final = {i: 0 for i in kmers}

    for kmer in kmer_dict_final.keys():
        kmer_dict_final[kmer] += sum([i[kmer] for i in threads])

    with open(kmers_path, "w") as file_out:
        file_out.write(
            "\n".join(
                ["%s\t%s" %(kmer, kmer_dict_final[kmer]) for kmer in kmers]
            )
        )

    print_time_elapsed(time_start)
################################################################################

if __name__ == "__main__":
    main()

