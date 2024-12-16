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



def assign_kmers(read_seqs, _mirna, _sitelist, n_constant, rand_length,
                 kmer_length, n_sites, experiment):
    time_start = time.time()
    seq_length = rand_length + 2*n_constant
    _kmerlist = KmerList(kmer_length)
    for i, r in enumerate(read_seqs):
        if i % 100 == 0:
            print(i)
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

def main():
    time_start = time.time()
    arguments = ["miRNA", "experiment", "condition", "n_constant", "sitelist",
                 "kmer_length", "n_sites", "-test_binary"]
    args = parse_arguments(arguments)
    mirna, experiment, condition, n_constant, sitelist, kmer_length, n_sites, test = args
    # Define mirna object:
    _mirna = Mirna(mirna)
    _sitelist = SiteList(_mirna, sitelist)
    if experiment in ["equil_mmseed_nb", "equil_seed_nb"]:
        rand_length = 38
    else:
        rand_length = 37    
    mirna_trim = "-".join(mirna.split("-")[:2])
    reads_path = get_analysis_path(mirna, experiment, condition, "reads")
    with open(reads_path, "rb") as file_in:
        if test:
            threads = multiprocess_test(file_in,
                                        readline_one,
                                        int(100),
                                        assign_kmers,
                                        100,
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

    print(_kmers.df().iloc[:20, :])

    _kmers.df().to_csv(kmer_counts_path, sep="\t", header=False)
################################################################################

if __name__ == "__main__":
    _trie = main()
