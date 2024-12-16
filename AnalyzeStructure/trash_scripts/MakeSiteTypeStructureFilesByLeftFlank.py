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

def get_average_structure(read_seqs,order_site_map,site_,win,ind_start=0, ind_stop=0):
    """Takes a read sequence and identifies the position of all site types'.

    Args:
        read_sequence: The read sequence.
        sites: A list of sequences corresponding to each site type.

    Returns:
        A tuple containing the list of reads and a dictionary of read types.
    """

    # Sub-function to be used on each site:
    # Defines the centered x coordinates for each tiling window
    with open('general/site_5p_pos.txt',"r+") as site_pos_file:
        site_pos_map = dict()
        for line in site_pos_file:
            key, value = line.strip().split("\t")
            site_pos_map[key] = int(value)
    with open('general/site_length.txt',"r+") as site_pos_file:
        site_len_map = dict()
        for line in site_pos_file:
            key, value = line.strip().split("\t")
            site_len_map[key] = int(value)

    kmers = get_kmer_list(4)
    # Defines the range of possible window positions relative to the
    if win%2==0:
        flank_averageprob_map = {kmer: [0]*(site_len_map[site_]+4+1) for kmer in kmers}
        flank_averagelogprob_map = {kmer: [0]*(site_len_map[site_]+4+1) for kmer in kmers}
    else:
        flank_averageprob_map = {kmer: [0]*(site_len_map[site_]+4) for kmer in kmers}
        flank_averagelogprob_map = {kmer: [0]*(site_len_map[site_]+4) for kmer in kmers}
    w_l = -(win-1)/2
    w_r = (win-1)/2+1
    print(w_l)
    print("w_r")
    print(w_r)
    flank_totals_map = {kmer: 0 for kmer in kmers}
    for count,seq in enumerate(read_seqs):
        if count%100000 == 0:
            print(count)
        read, sites, structure = (j.strip() for j in seq)
        if site_ in sites:
            coords = [i.split(":")[1] for i in sites.split(", ")]
            sites = [i.split(":")[0] for i in sites.split(", ")]
            coord_site_map = {site: coord for site, coord in zip(sites,coords)}
            ranks = [order_site_map[site] for site in sites]
            [(coord, site)] = [i for i in zip(coords, sites)
                               if order_site_map[i[1]] == min(ranks)]
            if site_ == site:
                start = int(coord.split("-")[0])+26 - 5 - ind_start
                stop = int(coord.split("-")[1]) + 1 + 26 - 5 - ind_start
                if start > (26 + 4) and stop < (26 + 37):
                    # print(read)
                    # print(" "*start+read[start:stop])
                    # for i in range(len(read)):
                    #     print(" "*i+read[i])
                    #     print(" "*i+str(structure.split("\t")[i]))
                    # print(" "*(start-4)+read[start-4:stop])
                    # print("\n".join(structure.split("\t")[start-4:stop]))
                    # print("--")
                    # print("\n".join(structure.split("\t")[start-4-(win-1)%2:stop]))
                    # print(structure.split("\t")[start])
                    # print([[math.log10(1 - float(j)) for j in structure.split("\t")[i+w_l:i+w_r]] for i in range(start-4-(win-1)%2,stop)])

                    # print([sum([math.log10(1 - float(j)) for j in structure.split("\t")[i+w_l:i+w_r]]) for i in range(start-4-(win-1)%2,stop)])
                    p = [10**sum([math.log10(1- float(j)) for j in structure.split("\t")[i+w_l:i+w_r]]) for i in range(start-4-(win-1)%2,stop)]
                    # print("p")
                    # print("p")
                    # print("\n".join([str(i) for i in prob]))
                    # print("\n".join([i for i in prob]))
                    logp = [math.log10(i) for i in p]
                    # print(logp)
                    # print("p")
                    # print(prob)
                    # print(site)
                    flank = read[start-4:start]
                    n = flank_totals_map[flank]
                    flank_averageprob_map[flank] = [(i[0]*n + i[1])/(n+1) for i in zip(flank_averageprob_map[flank],p)]
                    flank_averagelogprob_map[flank] = [(i[0]*n + i[1])/(n+1) for i in zip(flank_averagelogprob_map[flank],logp)]
                    flank_totals_map[flank] +=1
    return flank_averageprob_map, flank_averagelogprob_map, flank_totals_map


def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["miRNA","experient","condition", "site", "win"]
    mirna, experiment, condition, site, win = parse_arguments(arguments)
    win = int(win)
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
    site_flank_p_position_path = get_analysis_path(mirna,experiment,condition,"4ntflank_5p_unpaired_probs/p",ext= "_%s_%s" %(site, win))
    site_flank_logp_position_path = get_analysis_path(mirna,experiment,condition,"4ntflank_5p_unpaired_probs/logp",ext= "_%s_%s" %(site, win))
    site_flank_totals_path = get_analysis_path(mirna,experiment,condition,"4ntflank_5p_unpaired_probs/totals",ext= "_%s_%s" %(site, win))

    with open(full_reads_path,"rb") as file_in_reads:
        with open(full_sites_path,"rb") as file_in_sites:
            with open(structure_pairs_path, "rb") as file_in_structures:
                results = multiprocess_file([file_in_reads, file_in_sites, file_in_structures],
                                            readline_three, 1000000, 
                                            get_average_structure,
                                            order_site_map,
                                            site,
                                            win,
                                            int(ind_start),
                                            int(ind_stop))
        # Collect the output from each thread, which are the list
        # of reads and the dictionary of read types.
    # print(results[0])
    with open('general/site_length.txt',"r+") as site_pos_file:
        site_len_map = dict()
        for line in site_pos_file:
            key, value = line.strip().split("\t")
            site_len_map[key] = int(value)
    if win%2 == 0:
        range_site = range(site_len_map[site]+4+1)
        range_labels = range(0,site_len_map[site]+1)[::-1]
    else:
        range_site = range(site_len_map[site]+4)
        range_labels = range(1,site_len_map[site]+1)[::-1]        
    kmers = get_kmer_list(4)
    kmer_totals_map = {kmer: sum([i[2][kmer] for i in results]) for kmer in kmers}
    kmer_p_map = {kmer: [sum([i[0][kmer][j]*i[2][kmer]/kmer_totals_map[kmer] for i in results]) for j in range_site] for kmer in kmers if kmer_totals_map[kmer] > 0}
    kmer_logp_map = {kmer: [sum([i[1][kmer][j]*i[2][kmer]/kmer_totals_map[kmer] for i in results]) for j in range_site] for kmer in kmers if kmer_totals_map[kmer] > 0}
    # kmer_logp_map = {kmer: sum([i[kmer][1]*i[kmer][2] for i in results])/kmer_totals_map[kmer] for kmer in kmers}

    # total_kmers = sum(kmer_totals_map.values())
    print(range_labels)
    with open(site_flank_p_position_path,"wb") as file_out:
        file_out.write("\t".join(["kmer"] + ["flank_+%s" %(str(i+(win-1)%2*0.5).split(".0")[0]) for i in range(4,0,-1)]+ ["position_%s" %(str(i+(win-1)%2*0.5).split(".0")[0]) for i in range_labels]) + "\n")
        file_out.write("".join(["\t".join([kmer] + [str(i) for i in kmer_p_map[kmer]])+"\n" for kmer in kmers if kmer_totals_map[kmer] > 0]))

    with open(site_flank_logp_position_path,"wb") as file_out:
        file_out.write("\t".join(["kmer"] + ["flank_+%s" %(str(i+(win-1)%2*0.5).split(".0")[0]) for i in range(4,0,-1)]+ ["position_%s" %(str(i+(win-1)%2*0.5).split(".0")[0]) for i in range_labels]) + "\n")
        file_out.write("".join(["\t".join([kmer] + [str(i) for i in kmer_logp_map[kmer]])+"\n" for kmer in kmers if kmer_totals_map[kmer] > 0]))

    with open(site_flank_totals_path,"wb") as file_out:
        file_out.write("\t".join(["kmer"] + ["total"]) + "\n")
        file_out.write("".join(["\t".join([kmer] + [str(kmer_totals_map[kmer])])+"\n" for kmer in kmers if kmer_totals_map[kmer] > 0]))





        # Print the amount of time the script took to complete.
        print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()

