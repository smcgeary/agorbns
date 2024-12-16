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

def get_average_structure(read_seqs,order_site_map,win_size,ind_start=0, ind_stop=0,constant=False):
    """Takes a read sequence and identifies the position of all site types'.

    Args:
        read_sequence: The read sequence.
        sites: A list of sequences corresponding to each site type.

    Returns:
        A tuple containing the list of reads and a dictionary of read types.
    """

    # Sub-function to be used on each site:
    if constant == True:
        seq_range = range(len(read_seqs[0][0].strip())-win_size+1)
    else:
        seq_range = range(26,26+37-win_size+1)
    sites_range = range(26-ind_start,26+37+ind_stop+1)
    print(seq_range)
    print(sites_range)
    # Defines the centered x coordinates for each tiling window
    win_pos = [sum(range(i,i+win_size))/float(win_size) for i in range(26+37-win_size+1)]
    print(win_pos)
    with open('general/site_info.txt',"r+") as site_pos_file:
        site_pos_file.readline()
        site_pos_map = dict()
        for line in site_pos_file:
            line_strip = line.strip().split("\t")
            key = line_strip[0]
            value = line_strip[3]
            if key in order_site_map.keys() and key != "None":
                site_pos_map[key] = int(value)
    # Defines the range of possible window positions relative to the mirna
    win_pos_rel = sorted(set([i - j for i, j in it.product(win_pos,range(70))]))
    prob_averages = {site_key: {key: 0 for key in win_pos_rel} for site_key in order_site_map.keys()}
    prob_totals = {site_key: {key: 0 for key in win_pos_rel} for site_key in order_site_map.keys()}
    for count,seq in enumerate(read_seqs):
        if count%100000 == 0:
            print(count)
        read, sites, structure = (j.strip() for j in seq)
        print("hi")
        if sites != "None":
            print(read)
            print(sites)
            print(structure)
            coords = [i.split(":")[1] for i in sites.split(", ")]
            sites = [i.split(":")[0] for i in sites.split(", ")]
            coord_site_map = {site: coord for site, coord in zip(sites,coords)}
            ranks = [order_site_map[site] for site in sites]
            [(coord, site)] = [i for i in zip(coords, sites)
                               if order_site_map[i[1]] == min(ranks)]
            start = int(coord.split("-")[0])+26 - 5 - ind_start
            stop = int(coord.split("-")[1]) + 1 + 26 - 5 - ind_start
            if start in sites_range and stop in sites_range:
                prob_pos = [math.log10(1-float(j)) for j in structure.split("\t")]
                prob_win_pos = [sum(prob_pos[i:i+win_size])/win_size for i in range(len(read)-win_size+1)]
                for i in range(len(read)):
                    if i in seq_range:
                        rel_pos = win_pos[i]-(start+site_pos_map[site]-1)
                        prob_averages[site][rel_pos] += prob_win_pos[i]
                        prob_totals[site][rel_pos] += 1
    return prob_averages, prob_totals


def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["miRNA","experient","condition", "start", "stop", "sitelist"]
    mirna, experiment, condition, start, stop, sitelist = parse_arguments(arguments)

    # Load file with site types for the miRNA.
    sites_file_name = ("/lab/bartel1_ata/mcgeary/computation/"
                       "AgoRBNS/AssignSiteTypes/sites.%s_%s.txt" % (mirna, sitelist))
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
    print(sites)
    extension = "_%s-%s_%s" %(start, stop, sitelist)
    full_reads_path = get_analysis_path(mirna, experiment, condition, "full_reads")
    full_sites_path = get_analysis_path(mirna,experiment,condition,"full_sites",ext = extension)
    structure_pairs_path = get_analysis_path(mirna,experiment,condition,"structures_bp_prob")
    for win_size in range(1,21):
        site_averages_path = get_analysis_path(mirna,experiment,condition,"structures_bp_prob",ext="_%s-%s_geo_average_window_%s" %(ind_start,ind_stop,win_size))
        site_totals_path = get_analysis_path(mirna,experiment,condition,"structures_bp_prob",ext="_%s-%s_geo_totals_window_%s" %(ind_start,ind_stop,win_size))

        with open(full_reads_path,"rb") as file_in_reads:
            with open(full_sites_path,"rb") as file_in_sites:
                with open(structure_pairs_path, "rb") as file_in_structures:
                    results = multiprocess_test([file_in_reads, file_in_sites, file_in_structures],
                                                readline_three, 1000000, 
                                                get_average_structure, 1,
                                                order_site_map,
                                                int(win_size),
                                                int(ind_start),
                                                int(ind_stop),
                                                False)
        # Collect the output from each thread, which are the list
        # of reads and the dictionary of read types.
        # prob_threads = [i[0] for i in results]
        # totals_threads = [i[1] for i in results]
        # # print(bp_prob_threads[0])
        # # Flatten the list into one list (found on stackoverflow) and write
        # # it to its output file.
        # position_keys = ([j.keys() for i in prob_threads for j in i.values()])
        # # Flatten the list into one list (found on stackoverflow) get the unique
        # # number values within.
        # position_keys = sorted(list(set([i for j in position_keys for i in j])))
        # out_probs = pd.DataFrame.from_dict(
        #     {s_key : {p_key: sum([i[s_key][p_key] for i in prob_threads])
        #               for p_key in position_keys} for s_key in sites})
        # out_totals = pd.DataFrame.from_dict(
        #     {s_key : {p_key: sum([i[s_key][p_key] for i in totals_threads])
        #               for p_key in position_keys} for s_key in sites})
        # out_probs = out_probs/out_totals
        # out_probs.to_csv(site_averages_path,sep="\t",columns=out_probs.columns)
        # out_totals.to_csv(site_totals_path,sep="\t",columns=out_totals.columns)


        # Print the amount of time the script took to complete.
        print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()

