################################################################################
#MakeSiteTypeReadFiles.py
################################################################################
import imp # Used to import general.py
import time # Used for time
import math
import csv
import itertools as it
import multiprocessing
import sys
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
imp.load_source("general",
                ("/lab/bartel1_ata/mcgeary/computation/"
                 "AgoRBNS/general/general.py"
                )
               )
from general import parse_arguments, print_time_elapsed, get_analysis_path, multiprocess_test, readline, readline_three
from general import *
# This script reads one read file and outputs multiple .txt files to be used in all
# downstream applications. Inputs are:

# mirna: The miRNA in the experiment to be processed.
# experiment: The type of experiment to be processed.
# condition: The sampel type within the experiment to be processed.

def pull_structures(read_seqs, order_site_map, ind_start=5, ind_stop=5):
    """Takes a read sequence and identifies the position of all site types'.

    Args:
        read_sequence: The read sequence.
        sites: A list of sequences corresponding to each site type.

    Returns:
        A tuple containing the list of reads and a dictionary of read types.
    """

    # Sub-function to be used on each site:
    sites_range = range(26-ind_start,26+37+ind_stop+1)
    # Defines the centered x coordinates for each tiling window
    with open('general/site_5p_pos.txt',"r+") as site_pos_file:
        site_pos_map = dict()
        for line in site_pos_file:
            key, value = line.strip().split("\t")
            site_pos_map[key] = int(value)
    # Defines the range of possible window positions relative to the
    read_site_maps = {site: [] for site in order_site_map.keys() if site != "None"}
    struc_site_maps = {site: [] for site in order_site_map.keys() if site != "None"}
    for count,seq in enumerate(read_seqs):
        if count%100000 == 0:
            print(count)
        read, sites, structure = (j.strip() for j in seq)
        if sites != "None":
            coords = [i.split(":")[1] for i in sites.split(", ")]
            sites = [i.split(":")[0] for i in sites.split(", ")]
            coord_site_map = {site: coord for site, coord in zip(sites,coords)}
            ranks = [order_site_map[site] for site in sites]
            [(coord, site)] = [i for i in zip(coords, sites)
                               if order_site_map[i[1]] == min(ranks)]
            start = int(coord.split("-")[0])+26 - 5 -ind_start
            stop = int(coord.split("-")[1]) + 1 + 26 - 5 - ind_start
            l_1, l_2 = read[start-2], read[start-1]
            r_1, r_2 = read[stop], read[stop+1]

            if start in sites_range and stop in sites_range:
                out_row = [str(start), l_1, l_2, r_1, r_2, structure]
                struc_site_maps[site].append("\t".join(out_row))
                read_site_maps[site].append(read)

    return struc_site_maps, read_site_maps


def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["miRNA", "experient", "condition"]
    mirna, experiment, condition, = parse_arguments(arguments)

    # Load file with site types for the miRNA.
    sites_file_name = ("/lab/bartel1_ata/mcgeary/computation/"
                       "AgoRBNS/AssignSiteTypes/sites.%s.txt" % (mirna))
    add = True
    with open(sites_file_name) as file_in:
        sites = file_in.read().split("\n")
        # This swaps the four "Centered-..." for one single "Centered" entry in the list.
    sites_new = []
    for i in sites:
        if "Centered" not in i:
            sites_new.append(i)
        elif add:
            sites_new.append("Centered")
            add = False
    sites = sites_new   
    order_site_map = {site: i for i, site in enumerate(sites + ["None"])}
    ind_start = 0
    ind_stop = 0
    # Get the path to the read file and to that of where the site labels will
    # be written.
    full_reads_path = get_analysis_path(mirna, experiment, condition, "full_reads")
    full_sites_path = get_analysis_path(mirna,experiment,condition,"full_sites",ext = "_%s-%s" %(5, 5))
    structure_pairs_path = get_analysis_path(mirna,experiment,condition,"structures_bp_prob")

    with open(full_reads_path,"rb") as file_in_reads:
        with open(full_sites_path,"rb") as file_in_sites:
            with open(structure_pairs_path, "rb") as file_in_structures:
                results = multiprocess_file([file_in_reads, file_in_sites, file_in_structures],
                                            readline_three, 1000000, 
                                            pull_structures,
                                            order_site_map,
                                            int(ind_start),
                                            int(ind_stop))
    for site in sites:
        site_struct_path = get_analysis_path(mirna,experiment,condition,"structures_bp_prob/%s/" %(site),ext="_%s-%s" %(ind_start,ind_stop))
        site_reads_path = get_analysis_path(mirna,experiment,condition,"reads_by_site/%s/" %(site),ext="_%s-%s" %(ind_start,ind_stop))

        # Collect the output from each thread, which are the list
        # of reads and the dictionary of read types.
        output = [j for i in results for j in zip(i[0][site], i[1][site])]
        # print(reads)
        output.sort(key = lambda i: int(i[0].split("\t")[0]))
        structs = [i[0] for i in output]
        reads = [i[1] for i in output]

        struct_length = len(structs[0].split("\t"))-5
        with open(site_struct_path,"wb") as file_out:
            file_out.write("\t".join(["site pos", "left 1", "left 2", "right 1", "right 2"] + ["p%s" %(i) for i in range(1,struct_length+1)]+["\n"]))
            file_out.write("".join(["%s\n" % (i) for i in structs]))

        with open(site_reads_path,"wb") as file_out:
            file_out.write("Reads:\n")
            file_out.write("".join(["%s\n" % (i) for i in reads]))


    # Print the amount of time the script took to complete.
    print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()

