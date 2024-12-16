################################################################################
#MakeSiteTypeReadFiles.py
################################################################################
import imp # Used to import general.py
import time # Used for time
import os
import subprocess
import numpy as np
import pandas as pd
import math
import sys
import time
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

p5_seq = "GGGCAGAGTTCTACAGTCCGACGATC"
p3_seq = "TGTTCGTATGCCGTCTTCTGCTTG"

p5_seq_RC = get_rc(p5_seq)
p3_seq_RC = get_rc(p3_seq)

OUTPUT_LEN = 2


def assign_site_type_to_read(read_seqs, _sitelist, n_constant, read_length,
                             _mirna, experiment, seq_site_map, site_seq_map):
    """Takes a read sequence and identifies the position of all site types'.

    Args:
        read_sequence: The read sequence.
        sites: A list of sequences corresponding to each site type.

    Returns:
        A tuple containing the list of reads and a dictionary of read types.
    """

    # Sub-function to be used on each site:
    # def get_unique_sites(read_seq, site_seq_map):
    #     site_maps = []
    #     # Iterate over keys in dictionary (site seqs)
    #     for key, value in site_seq_map.items():
    #         # Assign left-hand value of site location, if it exists.
    #         i_l_list = []
    #         start = 0
    #         i_l = read_seq[start:].find(key) + start
    #         while i_l != -1:
    #             i_l_list = i_l_list + [i_l + start]
    #             start += i_l + 1
    #             i_l = read_seq[start:].find(key)
    #         i_r_list = [i + len(key) - 1 for i in i_l_list]
    #         # Check if index is -1
    #         site_maps_temp = [(value, i_l, i_r) for i_l, i_r in zip(i_l_list, i_r_list)]
    #         for site_map in site_maps_temp:
    #             # Check if site_maps is empty (suggesting first map)
    #             if site_maps == []:
    #                 site_maps.append(site_map)
    #             else:
    #                 add = True
    #                 i_l, i_r = site_map[1:]
    #             # Else only add if the sites in the list do not overlap.
    #                 site_maps_copy = site_maps[:]
    #                 for site_map_j in site_maps_copy:
    #                     j_l, j_r = site_map_j[1: ]
    #                     # Checks if new map is an extension of old map
    #                     if i_l <= j_l and i_r >= j_r:                            
    #                         # Removes map
    #                         site_maps.remove(site_map_j)
    #                     elif ((i_l >= j_l and i_r < j_r)
    #                           or (i_l > j_l and i_r <= j_r)):
    #                         add = False
    #                 if add:
    #                     site_maps.append(site_map)
    #     return site_maps

    # Initialize count dictionary with keys as site names and a value of
    # for each item.
    count_site_map = {value: {"TGT" : 0, "ACA" : 0}
                      for value in _sitelist.sites}
    count_multisite_map = {value: {"TGT" : 0, "ACA" : 0}
                           for value in _sitelist.sites}
    count_site_map.update({"None": {"TGT" : 0, "ACA" : 0}})

    count_multisite_map.update({"None": {"TGT" : 0, "ACA" : 0}})
    sites_range = range(26 - n_constant, 26 + 37 + n_constant + 1)
    time_start = time.time()

    # if conex:
    #     count_site_map.update({"ConstantUGU": {"TGT" : 0, "ACA" : 0}})
    #     count_multisite_map.update({"ConstantUGU": {"TGT" : 0, "ACA" : 0}})
    #     p5_overlap = {key : 0 for key in range(20)}
    #     p3_overlap = {key : 0 for key in range(20)}

    # sitekmers_count_map = {value : {kmer: {"TGT" : 0, "ACA" : 0}
    #                                 for kmer in get_kmer_list(kmer_length)}
    #                        for value in site_seq_map.values() + ["None"]}
    # Pre-allocate output with "None" for each line.
    output = ["None"]*len(read_seqs)
    _break = False
    for i_r, r in enumerate(read_seqs):
        # Get the read number:
        read = r.strip()
        # NEW Get the _read object, work with this.
        _read = Read(read)
        print(_read.seq)
        _read.get_all_sites_in_read(_sitelist)
        print("sites:")
        print(_read.print_sites())
        # print(_read.topsite.seq)
        # Find all sites within the read
        top_site = _read.site_name()
        top_pos = _read.site_pos()
        print(top_site)
        print(top_pos)
    #     barcode = read[26 + read_length : 26 + read_length + 3]

    #     # Deals with the TCG instead of TGT ending for miR-1 equilibrium exp
    #     if barcode == "TCG":
    #         barcode = "TGT"
    #         if _mirna.name != "miR-1" or experiment != "equilibrium":
    #             read = read[:26 + 37]+"TGTTCGTATGCCGTCTTCTGCTTG"
    #     if barcode == "TGT" and _mirna.name == "miR-1" and experiment == "equilibrium":
    #         read = read[:26 + 37]+"TCGTATGCCGTCTTCTGCTTG"


    # #     # Take the constant regiong + the number of constant sequences
    #     read_seq = read.strip()[(26 - n_constant) :
    #                             (26 + read_length + n_constant)]
    #     # Find all sites within the read
    #     site_maps = get_unique_sites(read_seq, site_seq_map)
    #     # seq_kmers = [read_seq[j:j+kmer_length] for j in range(read_length + 2*n_constant - kmer_length)]
    #     # if len(p5_kmers[0]) < 6 and len(p3_kmers[0]) < 6 and conex:
    #     # If there are mapped sites:
        # if site_maps:
        #     print(" "*(26 - n_constant) + read_seq)
        #     seq_map = " "*len(read_seq)
        #     for site_map in site_maps:
        #         start, stop = site_map[1:]
        #         for pos in range(start, stop + 1):
        #             if seq_map[pos] == " ":
        #                 seq_map = seq_map[:pos] + "_" + seq_map[pos+1:]
        #             elif seq_map[pos] == "_":
        #                 seq_map = seq_map[:pos] + "=" + seq_map[pos+1:]

        #     print(" "*(26 - n_constant) + seq_map)
        #     print(site_maps)
        #     if "=" in seq_map:
        #         _break = True
        # if _break:
        #     return
    #         # Label the output file:
    #         output_temp = ", ".join(["%s:%s-%s" % (site, index, end)
    #                                for site, index, end in site_maps])
    #         # print(output_temp)
    #         output[i] = output_temp
    #         # for kmer in seq_kmers:
    #         #     topsite = site_maps[0][0]
    #         #     if kmer not in topsite and topsite not in kmer:
    #         #         sitekmers_count_map[topsite][kmer][barcode] += 1
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
    #         # for i, kmer in enumerate(seq_kmers):
    #         #     sitekmers_count_map["None"][kmer][barcode] += 1
    #         count_site_map["None"][barcode] += 1
    #         count_multisite_map["None"][barcode] += 1
    #     # else:
    #     #     output[i] = "ConstantUGU"
    # return count_site_map, count_multisite_map

def main():
    time_start = time.time()
    # Define all the relevant arguments.
    arguments = ["miRNA", "experiment", "condition", "n_constant", "sitelist",
                 "-test_binary"]
    args = parse_arguments(arguments)
    mirna, experiment, condition, n_constant, sitelist, test, = args
    _mirna = Mirna(mirna)
    _sitelist = SiteList(_mirna, sitelist)
    if experiment in ["equil_mmseed_nb", "equil_seed_nb"]:
        read_length = 38
    else:
        read_length = 37    
    extension = "_%s_%s" %(n_constant, sitelist)

    mirna_trim = "-".join(mirna.split("-")[:2])
    sites_file_name = ("/lab/bartel1_ata/mcgeary/computation/"
                       "AgoRBNS/AssignSiteTypes/sites.%s_%s.txt" % (mirna_trim, sitelist))
    with open(sites_file_name) as file_in:
        sites = file_in.read().split("\n")

    # Initialize the dictionary with all sites, and two more dictionaries,
    # one mapping to the site name from the sequence, and one from site name
    # to the count.

    seq_site_map = get_seq_site_map(mirna_trim, sitelist)
    site_seq_map = {value: key for key, value in seq_site_map.items()}


    # if conex:
    #     extension = extension + "_conex-%s" %(conex)
    reads_path = get_analysis_path(mirna, experiment, condition, "full_reads")
    with open(reads_path, "rb") as file_in:
        if test:
            results = multiprocess_test(file_in,
                                        readline_one,
                                        int(1e6),
                                        assign_site_type_to_read,
                                        1000,
                                        _sitelist,
                                        int(n_constant),
                                        read_length,
                                        _mirna,
                                        experiment,
                                        seq_site_map,
                                        site_seq_map)
        else:
            results = multiprocess_iter(file_in,
                                        readline_one,
                                        int(1e6),
                                        assign_site_type_to_read,
                                        _sitelist,
                                        int(n_constant),
                                        read_length,
                                        _mirna,
                                        experiment)
    ############################################################################
    #OUTPUT
    ############################################################################
    print("done reading file")
    site_counts_path = get_analysis_path(mirna, experiment, condition,
                                         "site_counts", ext=extension)
    multisite_counts_path = get_analysis_path(mirna, experiment, condition,
                                              "multisite_counts",
                                              ext=extension)

    # Collect the output from each thread, which are the list
    # of reads and the dictionary of read types.
    results_list = [[r2 for r in results for r2 in r[i_r]] for i_r in range(OUTPUT_LEN)]
    print(results_list)
    # count_threads = [i[0] for i in results] for j in seq(1)
    # multicount_threads = [i[1] for i in results]
    # if conex:
    #     p5_overlap_threads = [i[4] for i in results]
    #     p3_overlap_threads = [i[5] for i in results]
    #     print(p5_overlap_threads)
    #     print(p3_overlap_threads)

    #     p5_overlap_total = [sum([thread[i] for thread in p5_overlap_threads]) for i in range(20)]
    #     p3_overlap_total = [sum([thread[i] for thread in p3_overlap_threads]) for i in range(20)]
    #     for i in range(20):
    #         print("%s\t%s" %(i, p5_overlap_total[i]))
    #     for i in range(20):
    #         print("%s\t%s" %(i, p3_overlap_total[i]))
    # # Flatten the list into one list (found on stackoverflow) and write
    # # it to its output file.
    # site_assignments = [i for thread in site_threads for i in sublist]
    # with open(sites_path,"wb") as file_out:
    #     file_out.write("".join(["%s\n" % (i) for i in site_assignments]))

    # # # Construct the dictionary from the list of futures, and write it to
    # # # its output file.
    # keys = set([j for i in dict_threads for j in i.keys()])
    # bcs = list(set([k for i in dict_threads for j in i.values() for k in j.keys()]))[::-1]
    # count_values = [tuple(sum([thread[key][bc] for thread in dict_threads if key in thread])
    #                       for bc in bcs)
    #                 for key in keys]
    # if experiment != "kinetics":
    #     count_values = [sum(i) for i in count_values]
    # counts_sites_map = {key : value for (key, value) in zip(keys, count_values)}
    # # No longer have input our output files open. Now write summmary file
    # # to the same directory , ending in "summary.txt".
    # with open(site_counts_path,"wb") as file_out:
    #     if experiment == "kinetics":
    #         for key in sorted(keys):
    #             value = counts_sites_map[key]
    #             file_out.write("%s:\t%s\t%s\n" % (key, value[0], value[1]))
    #     else:
    #         for key in sorted(keys):
    #             value = counts_sites_map[key]
    #             file_out.write("%s:\t%s\n" % (key, value))

    # # Construct the dictionary from the list of futures, and write it to
    # # its output file.
    # keys_multi = set([j for i in multi_threads for j in i.keys()])
    # multi_values = [tuple(sum([thread[key][bc] for thread in multi_threads if key in thread])
    #             for bc in bcs)
    #           for key in keys_multi]
    # if experiment != "kinetics":
    #     multi_values = [sum(i) for i in multi_values]
    # multicounts_sites_map = {key : value for (key, value) in zip(keys_multi, multi_values)}
    # # No longer have input our output files open. Now write summmary file
    # # to the same directory , ending in "summary.txt".
    # with open(multisite_counts_path,"wb") as file_out:
    #     if experiment == "kinetics":
    #         for key in sorted(keys_multi):
    #             value = multicounts_sites_map[key]
    #             file_out.write("%s:\t%s\t%s\n" % (key, value[0], value[1]))
    #     else:
    #         for key in sorted(keys_multi):
    #             value = multicounts_sites_map[key]
    #             print("%s\t%s" %(key, value))
    #             file_out.write("%s:\t%s\n" % (key, value))
    # kmers = sitekmers_threads[0]["None"].keys()
    # sitekmers_pulse_dict = {key : {kmer : sum(
    #     [i[key][kmer]["TGT"] for i in sitekmers_threads]) 
    #                                for kmer in kmers}
    #                         for key in keys if key != "ConstantUGU"}
    # sitekmers_chase_dict = {key : {kmer : sum(
    #     [i[key][kmer]["ACA"] for i in sitekmers_threads]) 
    #                                for kmer in kmers}
    #                         for key in keys if key != "ConstantUGU"}
    # kmers_pulse_df = pd.DataFrame.from_dict(sitekmers_pulse_dict)
    # kmers_chase_df = pd.DataFrame.from_dict(sitekmers_chase_dict)

    # extension_k = extension + "_k%s" %(kmer_length)
    # if experiment == "kinetics":
    #     kmer_pulse_counts_path = get_analysis_path(mirna, experiment, condition,
    #                                            "sitekmer_counts_pulse",
    #                                            ext=extension_k)
    #     kmer_chase_counts_path = get_analysis_path(mirna, experiment, condition,
    #                                            "sitekmer_counts_chase",
    #                                            ext=extension_k)
    #     kmers_pulse_df.to_csv(kmer_pulse_counts_path, sep="\t")
    #     kmers_chase_df.to_csv(kmer_chase_counts_path, sep="\t")
    #     print(kmer_pulse_counts_path)
    #     print(kmer_chase_counts_path)
    # else:
    #     kmer_counts_path = get_analysis_path(mirna, experiment, condition,
    #                                          "sitekmer_counts",
    #                                          ext=extension_k)
    #     kmers_df = kmers_pulse_df + kmers_chase_df
    #     kmers_df.to_csv(kmer_counts_path, sep="\t")
    #     print(kmer_counts_path)
        # Print the amount of time the script took to complete.
    print_time_elapsed(time_start)
    print(multisite_counts_path)
    # print(file_out)
################################################################################

if __name__ == "__main__":
    main()

