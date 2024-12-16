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
def count_read_kmers(read_seqs, _mirna, _kmerlist, experiment, n_constant,
                     rand_length, min_ex_str, buffer3p):
    time_start = time.time()
    sys.stdout.flush()
    for i_r, r in enumerate(read_seqs):
        read = r.strip()
        _read = Read(read, rand_length, _mirna, n_constant, experiment)
        exclude = sum([_read.seq.find(ex) for ex in min_ex_str])
        # print(exclude)
        # print(-1*len(min_ex_str))
        # print(_read)
        # print("."*n_constant + "_"*rand_length + "."*n_constant)

        if exclude == -1*len(min_ex_str):
            _kmerlist += _read
    return _kmerlist

def main():
    time_start = time.time()
    arguments = ["miRNA", "experiment", "condition", "n_constant", "kmer_len",
                 "-n_ex", "-uniq_binary", "-buffer3p_binary", "-final_binary",
                 "-weird", "-jobs", "-test_binary"]
    args = parse_arguments(arguments)
    (mirna, experiment, condition,
     n_constant, kmer_len, n_ex, uniq, buffer3p, final, weird, jobs, test) = args
    _mirna = Mirna(mirna)
    if experiment in ["equil_mmseed_nb", "equil_seed_nb"]:
        rand_length = 38
    else:
        rand_length = 37    
    _kmerlist = KmerList(int(kmer_len), int(n_constant), int(rand_length),
                         buffer_=buffer3p)

    # The new code by which the new sites are excluded from the analysis:#######
    mirna_base = "-".join(_mirna.name.split("-")[:2])
    if final:
        site_string = "paperfinal"
    elif weird:
        site_string = "papercutoffweirdsite%smer" %(weird)
    else:
        site_string = "papercutoff"


    _sitelist = SiteList(_mirna, site_string,
                         int(rand_length) + 2*int(n_constant))
    exclude_list = [i.name for i in _sitelist.sites if i]
    print(exclude_list)
    if n_ex:
        if n_ex[0] in ["A", "C", "T", "G"]:
            n_ex = n_ex.split(",")
            n_ex_strings = True
        else:
            n_ex = int(n_ex)
            n_ex_strings = False
    else:
        n_ex = len(exclude_list)
        n_ex_strings = True

    if uniq:
        read_dir = "reads_unique"
    else:
        read_dir = "reads"

    print(n_ex)
    if not n_ex_strings: 
        exclude_strings = [i.seq for i in _sitelist.sites][:n_ex]
    else:
        exclude_strings = n_ex
    print("exclude strings:")
    print(exclude_strings)
    min_ex_str = MinimalExcludeStringList(exclude_strings)


    if not jobs:
        jobs = 20

    if condition == "I_combined":
        if "tp" in experiment:
            input_list = INPUT_LIST_I_COMBINED_TP
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


    args  = [_mirna, _kmerlist, experiment, int(n_constant), int(rand_length),
             min_ex_str, buffer3p]
    ############################################################################

    threads = []


    for i_input in input_list:
        reads_path = get_analysis_path(i_input[0], i_input[1], i_input[2], read_dir)
        print(reads_path)
        if not jobs:
            jobs = 1
        threads += multiproc_file(reads_path, int(jobs), count_read_kmers, test,
                                  *args)


    
    kmers = get_kmer_list(int(kmer_len))[:10]
    for thread in threads:
        ns = thread.kmers["TGT"][:10]
        means = thread.pos_mean["TGT"][:10]
        sds = thread.pos_sd["TGT"][:10]
        print("thread:")
        for kmer, n, m, sd in zip(kmers, ns, means, sds):
            print(("%s\t%s\t%s\t%s" %(kmer, n, m, sd)))
    # print("done reading files")
    # print("summing threads:")
    _kmers = sum([thread for thread in threads])
    # print("Done summing threads.")
    # print(_kmers.pos_mean["TGT"][:10])
    pos_max = int(rand_length) + 2*int(n_constant) - int(kmer_len)
    # print("expected mean:")
    # print(pos_max/float(2))
    # print("expected sd:")
    expected_sd = (((pos_max + 1)**2 - 1)/12)**0.5
    # print(expected_sd)
    # for kmer, n, m, sd in zip(kmers, _kmers.kmers["TGT"][:10], _kmers.pos_mean["TGT"][:10],
    #                       _kmers.pos_sd["TGT"][:10]):
    #     print("%s\t%s\t%s\t%s" %(kmer, n, m, sd))
    extension = "_%s_k%s" %(n_constant, kmer_len)
    if uniq:
        extension = "%s_uniq" %(extension)
    if buffer3p:
        extension = "%s_buffer3p" %(extension)
    print((type(n_ex)))
    if type(n_ex) == list:
        extension = "%s_ex%s" %(extension, ",".join(exclude_strings))
    elif n_ex > 0:
        # extension_old = "%s_ex:%s" %(extension,
        #                              ",".join(exclude_list[:int(n_ex)]))
        extension = "%s_ex%s" %(extension, n_ex)
    if final:
        extension = "%s_final" %(extension)
    elif weird:
        extension = "%s_weirdsite%smer" %(extension, weird)
    kmers_path = get_analysis_path(mirna, experiment, condition,
                                   "kmers_cutoff_final", ext=extension)
    kmers_path_test = get_analysis_path(mirna, experiment, condition,
                                        "kmers_cutoff_final",
                                        ext=extension + "_test")
    print((sum(_kmers.kmers["TGT"])))
    print((sum(_kmers.kmers["ACA"])))
    if test:
        print(kmers_path_test)
        _kmers.write(kmers_path_test)
    else:
        print(kmers_path)
        print("writing kmer file:")
        _kmers.write(kmers_path)
        print("finished writing kmer file.")

    print_time_elapsed(time_start)
################################################################################

if __name__ == "__main__":
    main()

