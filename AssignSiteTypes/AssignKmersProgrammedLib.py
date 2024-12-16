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
                     rand_length, seedex):
    time_start = time.time()
    mm_sites = []
    # Assign the strings for the 8mer and 6mer sites
    site_8mer = _mirna["8mer"]
    site_6mer = _mirna["6mer"]
    site_6merA1 = _mirna["6mer-A1"]
    site_6merm8 = _mirna["6mer-m8"]
    print(site_8mer)
    # Make a list of the nucleotides of the 8mer site.
    site_8mer_list = list(site_8mer)
    # Make a list of all of the strings representing the 18 possible internal
    # 8mer mismatch sites.
    mismatches = [[site_8mer[:i] + j + site_8mer[i+1:]
                   for j in ["A", "C", "G", "T"] if j != site_8mer[i]]
                  for i in range(1, 7)]
    # Simplify the list
    mismatches = [j for i in mismatches for j in i]
    kmer_len = _kmerlist.len
    # print(kmer_len)
    # sys.stdout.flush()
    fraction_each_map = [0, 0]
    # shift_r = 0
    # shift_not_used = 0
    # Define the region where the programmed mismatch site should be within the
    # read.    
    lib_stop = 25

    for i_r, r in enumerate(read_seqs):
        read = r.strip()
        _read = Read(read, rand_length, _mirna, n_constant, experiment)
        # print(_read.seq)
        # sys.stdout.flush()
        mismatch_pos = [find_all(_read.seq, mm) for mm in mismatches]
        mm = [mm for i, mm in enumerate(mismatches) if len(mismatch_pos[i]) != 0]
        # print(mm)
        pos = [pos[0] for pos in mismatch_pos if len(pos) != 0 and pos[0]]
        if not (seedex and (
            site_6mer in _read.seq or site_6merA1 in _read.seq or
            site_6merm8 in _read.seq)): 
            if lib_stop + n_constant in pos:
                # offset = 0
                kmers = [_read.seq[i:i+kmer_len]
                         for i in range(len(_read.seq) - kmer_len + 1)]
                for i, kmer in enumerate(kmers):
                    _kmerlist[_read.barcode][kmer][i] += 1
                fraction_each_map[0] += 1
            # elif 24 + n_constant in pos:
            #     offset = 1
            #     add = True
            #     shift_l += 1
            else:
                print(pos)
                sys.stdout.flush()
                fraction_each_map[1] += 1
                # add = False
                # shift_not_used += 1
            # if add:
    # NOT ESSENTIAL #
    # Diagnostic of how many reads are 1.) at the right position,
    # 2., at the wrong position, and 3.) not counted.
    # percentages = [float(i) for i in fraction_each_map]
    percentages = [float(i)/np.sum(fraction_each_map) for i in fraction_each_map]
    print(percentages)
    print(fraction_each_map)
    sys.stdout.flush()
    return _kmerlist

def main():
    time_start = time.time()
    # Explanation of arguments:
    # n_constant : The number of constant sequence nucleotides appended to
    #   either side of the random/random+programmed portion of the library.
    # seedex_binary: 
    arguments = ["miRNA", "experiment", "condition", "n_constant", "kmer_len",
                 "-seedex_binary", "-jobs", "-temp_path", "-test_binary"]
    args = parse_arguments(arguments)
    (mirna, experiment, condition,
     n_constant, kmer_len, seedex, jobs, temp_path, test) = args

    if (("let-7a" in mirna or "miR-1" in mirna) and
        mirna != "miR-155_let-7a" and
        experiment in ["equil_c_nb", "equil_s_nb", "equil_c_alt_nb",
                       "equil_sc_nb", "equil_sc_alt_nb", "equil_c2_nb",
                       "equil_c2_alt_nb"]):
        rand_length = 38
    else:
        rand_length = 37



    _mirna = Mirna(mirna)

    # ###################################################################
    _kmerlist = KmerPositionalList(int(kmer_len),
                                   2*int(n_constant) + rand_length - int(kmer_len) + 1 + 3)
    read_dir = "reads"
    if not jobs:
        jobs = 20
    args  = [_mirna, _kmerlist, experiment, int(n_constant),
             int(rand_length), seedex]
    ########################################################################

    reads_path = get_analysis_path(mirna, experiment, condition, read_dir)
    threads = multiproc_file(reads_path, int(jobs), count_read_kmers, test,
                              *args)

    # Get the list of all k-mers of length `kmer-len`.
    kmers = get_kmer_list(int(kmer_len))
    # Make the dictionary of the counts from all threads.
    _kmers = sum([thread for thread in threads])
    # Convert the dictionary into a pandas dataframe.
    # Transpose so that the indeces are the kmers and the columns are the
    # positions within the random library.
    output = pd.DataFrame.from_dict(_kmers["TGT"]).T
    # NOT REQUIRED #
    print(output.sum(axis=0))


    # Format the output path:    
    extension = "_%s_k%s" %(n_constant, kmer_len)
    if seedex:
        extension = "%s_seedex" %extension
    if test:
        extension = "%s_test" %extension    
    if temp_path:
        kmers_path = "I_combined_temp/%s" %temp_path
    else:
        kmers_path = get_analysis_path(mirna, experiment, condition,
                                       "kmers_positional_programmed", ext=extension)
    # Print the path to allow for easy checking that the script worked, then
    # write the output file and the script is done.
    print(kmers_path)
    output.to_csv(kmers_path, index_label=False, sep="\t")
    print("finished writing kmer file.")
    print_time_elapsed(time_start)
################################################################################

if __name__ == "__main__":
    main()

