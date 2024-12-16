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

AU_weights_left = [1.0 / i for i in range(30,0,-1)]
AU_weights_right = [1.0 / i for i in range(1,31)]

def pull_structures(read_seqs, order_site_map, ind_start = 5, ind_stop = 5,
                    win_start = 1, win_stop = 15, constant = True):
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
                value = int(line_strip[2]) + int(line_strip[3]) - 1
                site_pos_map[key] = value
    # Defines the range of possible window positions relative to the
    struc_site_flanks_map = {site: {"".join(kmer): [] for kmer in list(it.product(["A","C","G","T"],repeat=4))} for site in order_site_map.keys() if site != "None"}
    AU_win_map = {site: {"".join(kmer): [] for kmer in list(it.product(["A","C","G","T"],repeat=4))} for site in order_site_map.keys() if site != "None"}
    AU_read_map = {site: {"".join(kmer): [] for kmer in list(it.product(["A","C","G","T"],repeat=4))} for site in order_site_map.keys() if site != "None"}
    AU_context_map = {site: {"".join(kmer): [] for kmer in list(it.product(["A","C","G","T"],repeat=4))} for site in order_site_map.keys() if site != "None"}

    tally_8mer = 0
    count = 0
    for count, seq in enumerate(read_seqs):
        seq = read_seqs[count]
        # if count%100000 == 0:
        #     print(count)
        read, sites, structure = (j.strip() for j in seq)
        read = "GGGCAGAGTTCTACAGTCCGACGATC" + read + "TATGCCGTCTTCTGCTTG"
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
            l_1, l_2 = read[start-2], read[start-1]
            r_1, r_2 = read[stop], read[stop+1]
            flank = "".join([l_1, l_2, r_1, r_2])

            if start in sites_range and stop in sites_range and site == "8mer" and flank == "AAAA":
                mir_win = range(start + site_pos_map[site] - win_stop, start + site_pos_map[site] + 1 - win_start)
                AU_score_L_l = max(0, start - 30)
                AU_score_L_r = start
                AU_score_R_l = stop 
                AU_score_R_r = min(len(read),stop + 30)
                AU_flanks_left = read[AU_score_L_l:AU_score_L_r]
                AU_flanks_right = read[AU_score_R_l:AU_score_R_r]
                # if site == "8mer":
                #     print(read)
                #     print(" "*(start)+read[start:stop])
                #     print(" "*(AU_score_L_l) + "."*(AU_score_L_r - AU_score_L_l) + " "*(stop - start) + "."*(AU_score_R_r - AU_score_R_l))
                #     print(" "*(AU_score_L_l) + read[AU_score_L_l:AU_score_L_r] + " "*(stop - start) + read[AU_score_R_l:AU_score_R_r])

                #     print(len(read[AU_score_L_l:AU_score_L_r]))
                #     print(len(read[AU_score_R_l:AU_score_R_r]))
                #     print(" "*(start+site_pos_map[site]-1)+read[start+site_pos_map[site]-1])
                #     print(AU_weights_left)
                #     print(flank)
                if "b" in site:
                    win_l = start + site_pos_map[site] - win_stop - 1
                else:
                    win_l = start + site_pos_map[site] - win_stop
                win_r = start + site_pos_map[site] + 1 - win_start
                if constant == False:
                    win_l = max(win_l,26)
                    win_r = min(win_r,26+37)
                # if site == "8mer":
                #     print(" "*win_l + "_"*(win_r-win_l))
                    # print(" "*win_l + read[win_l:win_r])
                AUleftoffset = 30 - len(AU_flanks_left)
                AUrightoffset = len(AU_flanks_right)
                AU_score_left = sum(AU_weights_left[:AUleftoffset]) + sum([score for i, score in enumerate(AU_weights_left[AUleftoffset : ]) if AU_flanks_left[i] in ["A","T"]])
                AU_score_right = sum([score for i, score in enumerate(AU_weights_right[ : AUrightoffset]) if AU_flanks_right[i] in ["A","T"]]) + sum(AU_weights_right[AUrightoffset:])
                # if site == "8mer":
                #     print("AU_score_left")
                #     print(AU_score_left)
                #     print("AU_score_right")
                #     print(AU_score_right)
                AU_score = AU_score_left + AU_score_right
                AU_win = float((read[win_l:win_r].count("A") + read[win_l:win_r].count("T")))/len(read[win_l : win_r])
                AU_read = float((read[26 : 26 + 37].count("A") + read[26 : 26 + 37].count("T")))/len(read[26 : 26 + 37])
                struc_win = 10**np.mean([math.log10(1 - float(i)) for i in structure.split("\t")[win_l:win_r]])

                if site == "8mer":
                    # print(AU_score / (sum(AU_weights_left) + sum(AU_weights_right)))
                    # print(AU_win)
                    # print(AU_read)
                    print("%s\t%s\t%s" %(count,struc_win, AU_win))
                    print(read)
                    print(" "*start + read[start:])

                    # for i in range(win_l,win_r):
                    #     print(" "*i+str(1 - float(structure.split("\t")[i])))
                    # print(" "*win_l + read[win_l:win_r])


                struc_site_flanks_map[site][flank].append(struc_win)
                AU_win_map[site][flank].append(AU_win)
                AU_read_map[site][flank].append(AU_read)
                AU_context_map[site][flank].append(AU_score)


            # if start in sites_range and stop in sites_range:
            #     out_row = [str(start), l_1, l_2, r_1, r_2, structure]
            #     struc_site_maps[site].append("\t".join(out_row))
            #     read_site_maps[site].append(read)

    return [struc_site_flanks_map, AU_win_map, AU_read_map, AU_context_map]


def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["miRNA", "experient", "condition", "ind_start", "ind_stop", "win_start", "win_stop"]
    mirna, experiment, condition, ind_start, ind_stop, win_start, win_stop = parse_arguments(arguments)

    # Load file with site types for the miRNA.
    sites_file_name = ("/lab/bartel1_ata/mcgeary/computation/"
                       "AgoRBNS/AssignSiteTypes/sites.%s_current.txt" % (mirna))
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
    full_reads_path = get_analysis_path(mirna, experiment, condition, "reads")
    full_sites_path = get_analysis_path(mirna,experiment,condition,"full_sites",ext = "_%s-%s" %(ind_start, ind_stop))
    structure_pairs_path = get_analysis_path(mirna,experiment,condition,"structures_bp_prob")
    AU_win_path = get_analysis_path(mirna,experiment,condition,"AU_win")
    AU_read_path = get_analysis_path(mirna,experiment,condition,"AU_read")
    AU_score_path = get_analysis_path(mirna,experiment,condition,"AU_cs")

    constant = True
    print(int(win_stop))
    # print(order_site_map)
    with open(full_reads_path,"rb") as file_in_reads:
        with open(full_sites_path,"rb") as file_in_sites:
            with open(structure_pairs_path, "rb") as file_in_structures:
                results = multiprocess_file([file_in_reads, file_in_sites, file_in_structures],
                                            readline_three,  10000,
                                            pull_structures,
                                            order_site_map,
                                            int(ind_start),
                                            int(ind_stop),
                                            int(win_start),
                                            int(win_stop),
                                            constant)
    # kmers = ["".join(kmer)for kmer in list(it.product(["A","C","G","T"],repeat=4))]
    # output_structures = {site: {kmer: [j for i in results for j in i[0][site][kmer]] for kmer in kmers} for site in order_site_map.keys() if site != "None"}
    # output_AU_win = {site: {kmer: [j for i in results for j in i[1][site][kmer]] for kmer in kmers} for site in order_site_map.keys() if site != "None"}
    # output_AU_read = {site: {kmer: [j for i in results for j in i[2][site][kmer]] for kmer in kmers} for site in order_site_map.keys() if site != "None"}
    # output_AU_cs = {site: {kmer: [j for i in results for j in i[3][site][kmer]] for kmer in kmers} for site in order_site_map.keys() if site != "None"}


    # # print(output_structures["8mer"])
    # # print(output_AU_win["8mer"])
    # # print(output_AU_read["8mer"])
    # print(output_AU_cs["8mer"])

    # summary = {site: {kmer: "0" for kmer in kmers} for site in order_site_map.keys() if site != "None"}
    # file_ext = "_%s-%s_%s-%s" %(ind_start,ind_stop,win_start,win_stop)
    # if constant == False:
    #     file_ext = file_ext + "_noconstant"
    # for site in ["8mer", "7mer-m8", "7mer-A1", "6mer"]:
    #     print(site)
    #     for flank in kmers:

    #         site_flank_path = get_analysis_path(mirna,experiment,condition,"site_flank_unpaired_data/%s/%s" %(site,flank),ext=file_ext)
    #         out = pd.DataFrame(output_structures[site][flank])
    #         n = out.size
    #         if n > 0:
    #             summary[site][flank] = ",".join([str(i) for i in [n, float(out.mean()), float(out.var())]])
    #         out.to_csv(site_flank_path,sep="\t",index = False, header = False)

    #         site_flank_AU_win_path = get_analysis_path(mirna,experiment,condition,"AU_win/%s/%s" %(site,flank),ext=file_ext)
    #         out = pd.DataFrame(output_AU_win[site][flank])
    #         out.to_csv(site_flank_AU_win_path,sep="\t",index = False, header = False)

    #         site_flank_AU_read_path = get_analysis_path(mirna,experiment,condition,"AU_read/%s/%s" %(site,flank),ext=file_ext)
    #         out = pd.DataFrame(output_AU_read[site][flank])
    #         out.to_csv(site_flank_AU_read_path,sep="\t",index = False, header = False)

    #         site_flank_AU_context_path = get_analysis_path(mirna,experiment,condition,"AU_cs/%s/%s" %(site,flank),ext=file_ext)
            
    #         out = pd.DataFrame(output_AU_cs[site][flank])
    #         out.to_csv(site_flank_AU_context_path,sep="\t",index = False, header = False)



    # parameters_path = get_analysis_path(mirna,experiment,condition,"site_flank_unpaired_data",ext=file_ext)
    # summary_out = pd.DataFrame.from_dict(summary)
    # summary_out.to_csv(parameters_path,sep="\t")

    # # # Print the amount of time the script took to complete.
    # print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()

