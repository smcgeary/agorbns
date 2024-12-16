################################################################################
#MakeSiteTypeReadFiles.py
################################################################################
import imp # Used to import general.py
import time # Used for time
import math
import csv
import ast
import itertools as it
import sys
import numpy as np
import pandas as pd
imp.load_source("general",
                ("/lab/bartel1_ata/mcgeary/computation/"
                 "AgoRBNS/general/general.py"
                )
               )
from general import parse_arguments, print_time_elapsed, get_analysis_path, multiprocess_file, multiprocess_test, readline, readline_three
from general import *
# This script reads one read file and outputs multiple .txt files to be used in all
# downstream applications. Inputs are:

# mirna: The miRNA in the experiment to be processed.
# experiment: The type of experiment to be processed.
# condition: The sampel type within the experiment to be processed.


def pull_structures(read_seqs, order_site_map, ind_start, ind_stop,
                    win_start, win_length, constant, mir_center):
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
    with open('general/site_info.txt',"r+") as site_pos_file:
        site_pos_map = dict()
        for line in site_pos_file:
            line_strip = line.strip().split("\t")
            key = line_strip[0]
            if key in order_site_map.keys() and key != "None":
                if mir_center:
                    # print("mirna center!")
                    value = int(line_strip[2]) + int(line_strip[3]) - 1
                else:
                    value = int(line_strip[2])
                site_pos_map[key] = value
    # print(site_pos_map)
    # Defines the range of possible window positions relative to the
    strucAndFlank_map = {site: [] for site in order_site_map.keys() if site != "None"}
    strucGeoMeanAndFlank_map = {site: [] for site in order_site_map.keys() if site != "None"}


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
            # print(read)
            # print(" "*26 + read[26 : 26 + 37])
            # print(" "*(26-ind_start) + "_"*(37 +ind_start + ind_stop))

            # print(site)
            start = int(coord.split("-")[0])+26 - ind_start
            stop = int(coord.split("-")[1]) + 1 + 26 - ind_start
            # print(" "*start + read[start:stop])
            if start in sites_range and stop in sites_range:
                l_1, l_2 = read[start-2], read[start-1]
                r_1, r_2 = read[stop], read[stop+1]
                flank = "".join([l_1, l_2, r_1, r_2])
                mir_win = range(start + site_pos_map[site] - win_start - win_length + 1, start + site_pos_map[site] + 1 - win_start)
                if "b" in site and mir_center:
                    win_l = start + site_pos_map[site] - win_start - win_length 
                else:
                    win_l = start + site_pos_map[site] - win_start - win_length + 1
                win_r = start + site_pos_map[site] + 1 - win_start
                if constant == False:
                    win_l = max(win_l,26)
                    win_r = min(win_r,26+37)
                # for i in range(win_l,win_r):
                #     print(" "*i+str(1 - float(structure.split("\t")[i])))
                # print(" "*win_l + read[win_l:win_r])
                if win_l < win_r:
                    strucMean_win = 10**np.mean([math.log10(1 - float(i)) for i in structure.split("\t")[win_l:win_r]])
                    struc_win = 10**np.sum([math.log10(1 - float(i)) for i in structure.split("\t")[win_l:win_r]])
                    # print(struc_win)
                    # print(win_r - win_l)
                    # print(win_length)
                    if win_r - win_l == win_length:
                        # print("ok")
                        strucAndFlank_map[site].append((flank,struc_win))
                    strucGeoMeanAndFlank_map[site].append((flank,strucMean_win))

    return([strucAndFlank_map, strucGeoMeanAndFlank_map])


def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["miRNA", "experient", "condition", "ind_start", "ind_stop", "sitelist", "win_start", "win_length", "constant", "mir_center"]
    mirna, experiment, condition, ind_start, ind_stop, sitelist, win_start, win_length, constant, mir_center = parse_arguments(arguments)

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
    # Get the path to the read file and to that of where the site labels will
    # be written.
    full_reads_path = get_analysis_path(mirna, experiment, condition, "full_reads")
    full_sites_path = get_analysis_path(mirna,experiment,condition,"full_sites",ext = "_%s-%s_%s" %(ind_start, ind_stop, sitelist))
    structure_pairs_path = get_analysis_path(mirna,experiment,condition,"structures_bp_prob")

    # print(int(win_stop))
    # print(order_site_map)
    constant = ast.literal_eval(constant)
    mir_center = ast.literal_eval(mir_center)
    with open(full_reads_path,"rb") as file_in_reads:
        with open(full_sites_path,"rb") as file_in_sites:
            with open(structure_pairs_path, "rb") as file_in_structures:
                results = multiprocess_file([file_in_reads, file_in_sites, file_in_structures],
                                            readline_three,  1000000,
                                            pull_structures,
                                            order_site_map,
                                            int(ind_start),
                                            int(ind_stop),
                                            int(win_start),
                                            int(win_length),
                                            constant,
                                            mir_center)
    kmers = ["".join(kmer)for kmer in list(it.product(["A","C","G","T"],repeat=4))]
    output_structures = {site: [j for i in results for j in i[0][site]]for site in order_site_map.keys() if site != "None"}
    output_structures_geo_mean = {site: [j for i in results for j in i[1][site]] for site in order_site_map.keys() if site != "None"}



    for site in sites:
        structure_flanks_path = get_analysis_path(mirna,experiment,condition,"probability_unpaired_by_site/%s/%s" %(sitelist, site), ext = "_%s-%s_%s_%s-%s_constant-%s_mir_center-%s" %(ind_start, ind_stop, sitelist, win_start, win_length, constant, mir_center))
        structure_flanks_geo_mean_path = get_analysis_path(mirna,experiment,condition,"probability_unpaired_geometric_mean_by_site/%s/%s" %(sitelist, site), ext = "_%s-%s_%s_%s-%s_constant-%s_mir_center-%s" %(ind_start, ind_stop, sitelist, win_start, win_length, constant, mir_center))

        out = pd.DataFrame(output_structures[site])
        if site == "8mer":
            print(out)
            print(structure_flanks_path)
        n = out.size
        out.to_csv(structure_flanks_path,sep="\t",index = False, header = False)
        if site == "8mer":
            print(out)
            print(structure_flanks_geo_mean_path)
        out = pd.DataFrame(output_structures_geo_mean[site])
        n = out.size
        out.to_csv(structure_flanks_geo_mean_path,sep="\t",index = False, header = False)



    # parameters_path = get_analysis_path(mirna,experiment,condition,"site_flank_unpaired_data",ext=file_ext)
    # summary_out = pd.DataFrame.from_dict(summary)
    # summary_out.to_csv(parameters_path,sep="\t")

    # Print the amount of time the script took to complete.
    print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()

