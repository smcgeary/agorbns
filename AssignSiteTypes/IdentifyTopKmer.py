################################################################################
#GenerateSiteTypeKds.py
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

START_TIME = time.time()

def assign_kmers(read_seqs, _mirna, _sitelist, n_constant, rand_length,
                 kmer_length, n_sites, experiment):
    time_start = time.time()
    seq_length = rand_length + 2*n_constant
    _kmerlist = KmerList(kmer_length)
    for i, r in enumerate(read_seqs):
        if i % 100000 == 0 & i > 0:
            print(i)
            print_time_elapsed(time_start)
            sys.stdout.flush()
        # Get the read number:
        read = r.strip()

        # Make Read object.
        _read = Read(read, rand_length, _mirna, n_constant, experiment)
        _read.get_all_sites_in_read(_sitelist)
        all_sites = _read.all_site_names(overlap=False)
        topsites = _sitelist.list_site_names()[:n_sites]
        seq = _read.seq
        if len(list(set(topsites) & set(all_sites))) == 0:
            _kmerlist += _read
    return _kmerlist

def get_site_kds(mirna, experiment="equilibrium", n_constant=5,
                 sitelist="paper"):
    path = get_analysis_path(mirna, experiment,
                                 "%s_%s_PAPER" %(int(n_constant), sitelist),
                                 "kds_PAPER")
    return(pd.read_csv(path, sep="\t"))

def main():
    time_start = time.time()
    arguments = ["miRNA", "experiment", "condition", "n_constant", "sitelist",
                 "kmer_length", "n_sites", "-test_binary"]
    args = parse_arguments(arguments)
    mirna, experiment, condition, n_constant, sitelist, kmer_length, n_sites, test = args
    # Define mirna object:
    _mirna = Mirna(mirna)
    _sitelist = SiteList(_mirna, sitelist)

    kds = get_site_kds(mirna, experiment, n_constant, sitelist)
    _sitelist.reorder(kds)

    print(_sitelist.list_site_names()[int(n_sites)])

    if experiment in ["equil_mmseed_nb", "equil_seed_nb"]:
        rand_length = 38
    else:
        rand_length = 37    
    mirna_trim = "-".join(mirna.split("-")[:2])
    reads_path = get_analysis_path(mirna, experiment, condition, "reads")
    with open(reads_path, "rb") as file_in:
        if test:
            total_reads = 100000
            process_size = total_reads/10
            threads = multiprocess_test(file_in,
                                        readline_one,
                                        process_size,
                                        assign_kmers,
                                        total_reads,
                                        _mirna,
                                        _sitelist,
                                        int(n_constant),
                                        int(rand_length),
                                        int(kmer_length),
                                        int(n_sites),
                                        experiment)
        else:
            threads = multiprocess_file(file_in,
                                        readline_one,
                                        int(1e6),
                                        assign_kmers,
                                        _mirna,
                                        _sitelist,
                                        int(n_constant),
                                        int(rand_length),
                                        int(kmer_length),
                                        int(n_sites),
                                        experiment)
    extension = "_%s_k%s_%ssites" %(n_constant, kmer_length, n_sites)
    ############################################################################
    print("done reading file")
    kmer_counts_path = get_analysis_path(mirna, experiment, condition,
                                         "kmers", ext=extension)
    _kmers = threads[0]
    if len(threads) > 1:
        for _kmer in threads[1:]:
            _kmers += _kmer
    kmer_df = _kmers.df()
    print(kmer_df.loc["ATAACAAAA"])
    print(kmer_df.loc["ATGACAAAA"])
    if not test:
        print("Writing file:")
        _kmers.df().to_csv(kmer_counts_path, sep="\t", header=False)
    print("Done script.")
    print_time_elapsed(time_start)

################################################################################

if __name__ == "__main__":
    _trie = main()
