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

np.set_printoptions(linewidth=200)

# FUNCTIONS
def count_dinuc_freqs(read_seqs, _mirna, experiment, n_constant,
                     rand_length):
    time_start = time.time()
    sys.stdout.flush()
    read_len = rand_length + 2*n_constant - 1
    dinucs = np.zeros((16, read_len))
    for i_r, r in enumerate(read_seqs):
        read = r.strip()
        _read = Read(read, rand_length, _mirna, n_constant, experiment)
        indeces = [get_kmer_index(_read.seq[i:i+2]) for i in range(read_len)]
        for i, j in enumerate(indeces):
            dinucs[j, i] += 1
    return dinucs

def main():
    time_start = time.time()
    arguments = ["miRNA", "experiment", "condition", "n_constant", "-jobs", "-test_binary"]
    args = parse_arguments(arguments)
    mirna, experiment, condition, n_constant, jobs, test = args
    _mirna = Mirna(mirna)
    if experiment in ["equil_mmseed_nb", "equil_seed_nb"]:
        rand_length = 38
    else:
        rand_length = 37    
    mirna_base = "-".join(_mirna.name.split("-")[:2])
    if not jobs:
        jobs = 20
    args  = [_mirna, experiment, int(n_constant), int(rand_length)]
    ############################################################################
    reads_path = get_analysis_path(mirna, experiment, condition, "reads")
    threads = multiproc_file(reads_path, int(jobs), count_dinuc_freqs, test,
                             *args)
    print("done reading files")
    print("summing threads:")
    dinucs = sum(threads)
    print(dinucs)
    extension = "_%s" %(n_constant)

    out_path = get_analysis_path(mirna, experiment, condition, "dinucleotide_counts", ext=extension)
    print(out_path)
    np.savetxt(out_path, dinucs, fmt="%i", delimiter="\t")
    if not test:
        print("writing output:")
        _kmers.write(kmers_path)
        print("finished writing kmer file.")
        print(kmers_path)
    print_time_elapsed(time_start)
################################################################################

if __name__ == "__main__":
    main()

