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
def count_kmer_freqs(read_seqs, _mirna, experiment, n_constant,
                      rand_length, kmer_length):
    time_start = time.time()
    # kpl is the kmer position length.
    kpl = rand_length + 2*n_constant - kmer_length + 1
    kmers = [0]*(4**kmer_length*kpl)
    if experiment in ["kinetics", "kinetics_pilot"]:
        kmer_chase = [0]*(4**kmer_length*kpl)
    tot_reads = len(read_seqs)
    for i_r, r in enumerate(read_seqs):
        if i_r % (float(tot_reads)/10.0) == 0:
            print(float(i_r)/float(tot_reads))
            sys.stdout.flush()
        read = r.strip()
        _read = Read(read, rand_length, _mirna, n_constant, experiment)
        print(_read.seq)
        sys.stdout.flush()
        indeces = [get_kmer_index(_read.seq[i:i+kmer_length])
                   for i in range(kpl)]
        if (experiment in ["kinetics", "kinetics_pilot"] and
            _read.barcode == "ACA"):
            for i, j in enumerate(indeces):
                kmers_chase[j*kpl + i] += 1
        else:
            for i, j in enumerate(indeces):
                kmers[j*kpl + i] += 1
    print("done with loop")
    sys.stdout.flush()
    if experiment in ["kinetics", "kinetics_pilot"]:
        return kmers, kmers_chase
    else:
        print("about to return kmers")
        sys.stdout.flush()
        return kmers

def main():
    # Get the time for the starting position

    time_start = time.time()
    arguments = ["miRNA", "experiment", "condition", "n_constant",
                 "kmer_length", "-addcomp", "-jobs", "-test_binary"]
    args = parse_arguments(arguments)
    (mirna, experiment, condition,
     n_constant, kmer_length, addcomp, jobs, test) = args
    _mirna = Mirna(mirna)
    if experiment in ["equil_mmseed_nb", "equil_seed_nb"]:
        rand_length = 38
    else:
        rand_length = 37    
    mirna_base = "-".join(_mirna.name.split("-")[:2])
    if not jobs:
        jobs = 20
    args  = [_mirna, experiment, int(n_constant),int(rand_length),
             int(kmer_length)]
    # MULTIPROCESSING PART #####################################################
    reads_path = get_analysis_path(mirna, experiment, condition, "reads")
    # Determine the number of reads in the file in order to simulate more reads
    # that contain the competitor oligo.
    threads = multiproc_file(reads_path, int(jobs), count_kmer_freqs, test,
                             *args)



    # kpl = kmer position length
    kpl = int(rand_length) + 2*int(n_constant) - int(kmer_length) + 1

    if experiment in ["kinetics", "kinetics_pilot"]:
        threads_pulse = [np.array(i[0]).reshape(len(i[0])/kpl, kpl)
                         for i in threads]
        threads_chase = [np.array(i[1]).reshape(len(i[1])/kpl, kpl)
                         for i in threads]
        kmers = sum([i for i in threads_pulse])
        kmers_chase = sum([i for i in threads_chase])
    else:
        threads = [np.array(i).reshape(len(i)/kpl, kpl)
                   for i in threads]
        kmers = sum(threads)


    if addcomp:
        print(kmers)
        num_reads = kmers.sum(axis=0)[0]
        print(num_reads)
        comp_reads = int(float(num_reads)*float(addcomp)/100)

        additional_reads = simulate_competitor_reads(comp_reads, mirna,
                                                     "equilibrium", 11)
        thread = count_kmer_freqs(additional_reads, _mirna, experiment,
                                        int(n_constant),int(rand_length),
                                        int(kmer_length))
        thread_comp = np.array(thread).reshape(len(thread)/kpl, kpl)
        print(thread_comp)
        kmers = kmers + thread_comp
        print(kmers.sum(axis=0))




    if experiment in ["kinetics", "kinetics_pilot"]:
        extension = "_pulse_%s_k%s" %(n_constant, kmer_length)
    else:
        extension = "_%s_k%s" %(n_constant, kmer_length)
    if addcomp:
        extension = extension + "_addcomp%s" %(addcomp)
    if test:
        extension = extension + "_test"
    kmers_path = get_analysis_path(mirna, experiment, condition,
                                   "positional_kmers", ext=extension)
    print(kmers_path)
    np.savetxt(kmers_path, kmers, fmt="%i", delimiter="\t")

    if experiment in ["kinetics", "kinetics_pilot"]:
        extension_chase = "_chase_%s_k%s" %(n_constant, kmer_length)
        if test:
            extension_chase = extension_chase + "_test"
        kmers_path_chase = get_analysis_path(mirna, experiment, condition,
                                             "positional_kmers",
                                             ext=extension_chase)
        print(kmers_path_chase)
        np.savetxt(kmers_path_chase, kmers_chase, fmt="%i", delimiter="\t")
    print_time_elapsed(time_start)
################################################################################

if __name__ == "__main__":
    main()

