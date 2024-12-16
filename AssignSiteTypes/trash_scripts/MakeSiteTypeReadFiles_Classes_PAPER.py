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
    # Initialize count dictionary with keys as site names and a value of
    # for each item.
    count_site_map = {value: {"TGT" : 0, "ACA" : 0}
                      for value in _sitelist.list_site_names()}
    count_multisite_map = {value: {"TGT" : 0, "ACA" : 0}
                           for value in _sitelist.list_site_names()}
    count_site_map.update({"None": {"TGT" : 0, "ACA" : 0}})

    count_multisite_map.update({"None": {"TGT" : 0, "ACA" : 0}})
    sites_range = range(26 - n_constant, 26 + 37 + n_constant + 1)
    time_start = time.time()
    # output = ["None"]*len(read_seqs)
    overlaps = 0
    for i_r, r in enumerate(read_seqs):
        # Get the read number:

        read_trim, read = [i.strip() for i in r]
        print(" "*26 + read_trim)
        print(read)
        # NEW Get the _read object, work with this.
        _read = Read(read, read_length, _mirna, experiment)
        _read.get_all_sites_in_read(_sitelist)
        # Find all sites within the read
        top_site = _read.site_name()
        top_pos = _read.site_pos()
        if top_site != "None":
            # if sum(_read.site_overlaps()):
            #     _read.print_sites()
            #     # print(_read.all_site_names())
            #     overlaps += 1
            #     print("top site")
            #     print(top_site)
            #     print("top pos")
            #     print(top_pos)
            count_site_map[top_site][_read.barcode] += 1
            multisite = ",".join(_read.all_site_names())
            # print(multisite)
            if multisite not in count_multisite_map.keys():
                count_multisite_map[multisite] = {"TGT": 0, "ACA": 0}
            count_multisite_map[multisite][_read.barcode] += 1
        else:
            count_site_map["None"][_read.barcode] += 1
            count_multisite_map["None"][_read.barcode] += 1
    return count_site_map, count_multisite_map

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
    reads_path = get_analysis_path(mirna, experiment, condition, "reads")
    full_reads_path = get_analysis_path(mirna, experiment, condition, "full_reads")

    with open(reads_path, "rb") as file_in:
        with open(full_reads_path, "rb") as file_in_2:
            if test:
                threads = multiprocess_test([file_in, file_in_2],
                                            readline_two,
                                            int(10000),
                                            assign_site_type_to_read,
                                            10000,
                                            _sitelist,
                                            int(n_constant),
                                            read_length,
                                            _mirna,
                                            experiment,
                                            seq_site_map,
                                            site_seq_map)
            else:
                threads = multiprocess_file([file_in, file_in_2],
                                            readline_two,
                                            int(1000),
                                            assign_site_type_to_read,
                                            _sitelist,
                                            int(n_constant),
                                            read_length,
                                            _mirna,
                                            experiment,
                                            seq_site_map,
                                            site_seq_map)
    ############################################################################
    #OUTPUT
    ############################################################################
    print("done reading file")
    site_counts_path = get_analysis_path(mirna, experiment, condition,
                                         "site_counts_classes", ext=extension)
    multisite_counts_path = get_analysis_path(mirna, experiment, condition,
                                              "multisite_counts_classes",
                                              ext=extension)

    dict_threads = [i[0] for i in threads]
    multi_threads = [i[1] for i in threads]


    keys = set([j for i in dict_threads for j in i.keys()])
    bcs = list(set([k for i in dict_threads for j in i.values() for k in j.keys()]))[::-1]
    count_values = [tuple(sum([thread[key][bc] for thread in dict_threads if key in thread])
                          for bc in bcs)
                    for key in keys]
    if experiment != "kinetics":
        count_values = [sum(i) for i in count_values]
    counts_sites_map = {key : value for (key, value) in zip(keys, count_values)}
    # No longer have input our output files open. Now write summmary file
    # to the same directory , ending in "summary.txt".
    with open(site_counts_path,"wb") as file_out:
        if experiment == "kinetics":
            for key in sorted(keys):
                value = counts_sites_map[key]
                # file_out.write("%s:\t%s\t%s\n" % (key, value[0], value[1]))
        else:
            for key in sorted(keys):
                value = counts_sites_map[key]
                # file_out.write("%s:\t%s\n" % (key, value))

    # Construct the dictionary from the list of futures, and write it to
    # its output file.
    keys_multi = set([j for i in multi_threads for j in i.keys()])
    multi_values = [tuple(sum([thread[key][bc] for thread in multi_threads if key in thread])
                for bc in bcs)
              for key in keys_multi]
    if experiment != "kinetics":
        multi_values = [sum(i) for i in multi_values]
    counts_multisites_map = {key : value for (key, value) in zip(keys_multi, multi_values)}

    counts = pd.DataFrame.from_dict(counts_sites_map, orient="index")
    counts = counts.reindex(_sitelist.list_site_names()+["None"])
    multicounts = pd.DataFrame.from_dict(counts_multisites_map, orient="index")
    print(counts)
    print(multicounts)
    # print(np.sum(counts))
    print(np.sum(multicounts))
    if not test:
        counts.to_csv(site_counts_path, sep="\t", header=False)
        multicounts.to_csv(multisite_counts_path, sep="\t", header=False)
    print_time_elapsed(time_start)
    # print(multisite_counts_path)
    # print(file_out)
################################################################################

if __name__ == "__main__":
    main()

