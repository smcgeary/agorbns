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
import itertools as it
import RNA
imp.load_source("general",
                ("/lab/bartel1_ata/mcgeary/computation/"
                 "AgoRBNS/general/general.py"
                )
               )

imp.load_source("sitetypes",
                ("/lab/bartel1_ata/mcgeary/computation/"
                 "AgoRBNS/AssignSiteTypes/sitetypes.py"
                )
               )
from general import parse_arguments, print_time_elapsed, get_kmer_list, get_analysis_path, multiprocess_file, multiprocess_test, readline_two, geo_mean, get_site_flanks
from sitetypes import get_seq_site_map
import random, string

def randomword(length):
   return ''.join(random.choice(string.lowercase) for i in range(length))

# This script reads one read file and outputs multiple .txt files to be used in all
# downstream applications. Inputs are:

# mirna: The miRNA in the experiment to be processed.
# experiment: The type of experiment to be processed.
# condition: The sampel type within the experiment to be processed.

# Functions from general this uses:


# Constants

# FUNCTIONS

AU_weights_left = [1.0 / i for i in range(30,0,-1)]
AU_weights_right = [1.0 / i for i in range(1,31)]
# AU_weights_right[0] = 0.5

def RNAplfold(read, win, file_num, logprob=False):
    # This should give the probability of the region in the read complementary
    # to position "mir_start" to "mir_stop" in the miRNA sequence.
    read_name = randomword(20)
    temp_fa_filename = "%s%s.fa" %(read_name, file_num)
    with open(temp_fa_filename, "wb") as f:
        f.write(">%s%s\n%s" %(read_name, file_num, read))
    # call RNAplfold:
    l_param = len(read)
    w_param = l_param
    u_param = win
    if logprob:
        log_param = "-O "
        temp_plfold_filename = '%s%s_openen' %(read_name, file_num)
    else:
        log_param = ""
        temp_plfold_filename = '%s%s_lunp' %(read_name, file_num)

    call = 'RNAplfold -L %s -W %s -u %s %s< %s%s.fa' %(l_param, w_param,
                                                       u_param, log_param,
                                                       read_name, file_num)
    subprocess.call([call], shell=True, stdout = subprocess.PIPE)
    rnaplfold_data = pd.read_csv(temp_plfold_filename, sep = '\t',
                                 header = 1).set_index(' #i$')
    rnaplfold_data = rnaplfold_data.iloc[:, :u_param]
    os.remove(temp_fa_filename)
    os.remove(temp_plfold_filename)
    os.remove('%s%s_dp.ps' %(read_name, file_num))
    return rnaplfold_data


def calculate_local_au(utr, site_start, site_end, site_type):
    """
    Calculate the local AU score

    Parameters
    ----------
    utr: string, utr sequence

    site_type: string, site type

    site_start: int, start of site

    site_end: int, end of site

    Output
    ------
    float: local AU score
    """
    # find A, U and weights upstream of site
    up_site_adder = int(site_type not in ['7mer-m8', '8mer'])
    
    upstream = utr[max(0, site_start - 30): site_start]
    upstream_str = upstream

    l= max(0, site_start - 30)
    upstream = [int(x in ['A', 'T']) for x in upstream]
    inv_upweights = [(x + 1 + up_site_adder)
                 for x in range(len(upstream))][::-1]

    upweights = [1.0 / (x + 1 + up_site_adder)
                 for x in range(len(upstream))][::-1]

    # find A,U and weights downstream of site
    down_site_adder = int(site_type in ['7mer-A1', '8mer'])
    downstream = utr[site_end:min(len(utr), site_end + 30)]
    downstream_str = downstream
    downstream = downstream + "A"*(30 - len(downstream))
    downstream = [int(x in ['A', 'T']) for x in downstream]

    inv_downweights = [(x + 1 + down_site_adder)
                   for x in range(len(downstream))]
    downweights = [1.0 / (x + 1 + down_site_adder)
                   for x in range(len(downstream))]
    # print("_"*l + upstream_str + "."*(site_end - site_start) + downstream_str)
    # print(" "*l + "".join([str(i) for i in upstream]) + "."*(site_end - site_start) + "".join([str(i) for i in downstream]))
    # for num, weight in enumerate(inv_upweights):
    #     print(" "*(l + num) + str(weight) + " "*(30 + site_end - site_start - len(str(weight))) + str(inv_downweights[num]))

    weighted = np.dot(upstream, upweights) + np.dot(downstream, downweights)
    total = float(sum(upweights) + sum(downweights))

    return weighted / total




def get_local_au_score(read, start, stop, sitetype):
    if sitetype not in ["8mer", "7mer-m8", "7mer-A1", "6mer"]:
        return False
    else:
        upstream_site = int(sitetype not in ["7mer-m8", "8mer"])
        upstream_seq = read[max(0, start - 30):start]

        upstream = [int(x in ["A", "T"]) for x in upstream_seq]
        upweights = [1.0 / (x + 1 + upstream_site)
                     for x in range(len(upstream))][::-1]
        upweights_total = [1.0 / (x + 1 + upstream_site)
                     for x in range(30)][::-1]
        downstream_site = int(sitetype in ["7mer-A1", "8mer"])
        downstream_seq = read[stop:min(len(read), stop + 30)]
        downstream = [int(x in ["A", "T"]) for x in downstream_seq]
        downstream = downstream + [1 for i in range(30 - len(downstream))]
        downweights = [1.0 / (x + 1 + downstream_site)
                       for x in range(len(downstream))]
        downweights_total = [1.0 / (x + 1 + downstream_site)
                       for x in range(30)]

        weighted = np.dot(upstream, upweights) + np.dot(downstream, downweights)
        total = float(sum(upweights_total) + sum(downweights_total))

        return(weighted / total)
    return

def get_read_structural_data(read_seqs, order_site_map, site_pos_map, n_constant,
                             read_length, mirna, experiment, win_start = -1, win_stop = 15, win = 5):
    """Takes a read sequence and identifies the position of all site types'.

    Args:
        read_sequence: The read sequence.
        sites: A list of sequences corresponding to each site type.

    Returns:
        A tuple containing the list of reads and a dictionary of read types.
    """
    # These are the slope and intercept for converting the logplfold to plfold.
    # (I had to figure this out manually using lm in R!)
    b = -0.00001798580677720675
    m = -1.62250426156107407927
    # Initialize count dictionary with keys as site names and a value of
    # for each item.
    sites_range = range(26 - n_constant, 26 + 37 + n_constant + 1)
    site_flanks_map = {site: {"".join(kmer): [] for kmer in get_kmer_list(4)} 
                       for site in ["8mer"]}
    # Pre-allocate output with "None" for each line.
    pl_pos_5p = ["f%s" %(i) for i in range(1, win_stop - 8 + 1)]
    pl_pos_3p = ["t%s" %(i) for i in range(1, -1*win_start + 1)]
    print(pl_pos_3p)
    pl_pos_seed = ["s%s" %(i) for i in range(1, 8 + 1)[::-1]]
    pl_pos_keys = pl_pos_5p + pl_pos_seed + pl_pos_3p
    pl_win_keys = ["w%s" %(i) for i in range(1, win+1)]
    pos_win_plfold_map = {i: {j : [] for j in pl_win_keys} for i in pl_pos_keys}
    pos_win_logplfold_map = {i: {j : [] for j in pl_win_keys} for i in pl_pos_keys}
    flanks = []
    time_start = time.time()
    tick = 0
    for num_i, i in enumerate(read_seqs):
        tick +=1
        if tick % 10000 == 0:
            print(tick)
            print_time_elapsed(time_start)
            sys.stdout.flush()
        read, sites = (j.strip() for j in i)
        # Define the barcode portion of the read.
        barcode = read[26 + read_length : 26 + read_length + 3]
        # Deals with the TCG instead of TGT ending for miR-1 equilibrium exp
        if barcode == "TCG":
            barcode = "TGT"
            if mirna != "miR-1" or experiment != "equilibrium":
                read = read[:26 + 37] + "TGTTCGTATGCCGTCTTCTGCTTG"
        if barcode == "TGT" and mirna == "miR-1" and experiment == "equilibrium":
            read = read[:26 + 37] + "TCGTATGCCGTCTTCTGCTTG"
        # Take the constant region + the number of constant sequences
        # Find all sites within the read
        if sites != "None":
            # CONVERT ALL COORDINATES TO 0 is 1, start of the read!
            coords = [i.split(":")[1] for i in sites.split(", ")]
            sites = [i.split(":")[0] for i in sites.split(", ")]
            coord_site_map = {site: coord for site, coord in zip(sites,coords)}
            ranks = [order_site_map[site] for site in sites]
            [(coord, site)] = [i for i in zip(coords, sites)
                               if order_site_map[i[1]] == min(ranks)]
            # Get the start and stop positions of the site, in pythonic
            # coordinates.
            start = int(coord.split("-")[0]) + 26 - n_constant
            stop = int(coord.split("-")[1]) + 1 + 26 - n_constant
            if start in sites_range and stop in sites_range and site in ["8mer"]:
                mirp1_span = site_pos_map[site] - 1
                # print("mirp1_span")
                # print(mirp1_span)
                mirp1 = start + mirp1_span
                win_r = mirp1 - win_start + 1
                win_l = mirp1 - win_stop + 1
                flank = read[start - 2:start] + "." + read[stop:stop + 2]
                flanks.append(flank)
                if "b" in site:
                    win_l -= 1
                # Perform plfold on the entire read:
                plfold = RNAplfold(read, win, num_i)
                logplfold = RNAplfold(read, win, num_i, logprob=True)
                test = plfold.iloc[:3, :3]
                test_log = logplfold.iloc[:3, :3]
                # print(plfold.iloc[:3, :3])
                # print(logplfold.iloc[:3, :3])
                # seed_inds = range(start, stop)[::-1]
                full_inds = range(win_l, win_r)
                # print(seed_inds)
                # print(full_inds)
                # print(len(full_inds))
                # print(pos_win_plfold_map.keys())
                # print(len(pos_win_plfold_map.keys()))

                struc = RNA.fold(read)[0]
                pos_zips = zip(pl_pos_keys, full_inds)
                # print(pos_zips)
                # print(plfold.shape)
                # print(" "*win_l + "."*(win_r - win_l))
                # for i, pos in enumerate(full_inds):
                #     if pl_pos_keys[i] == "seed8":
                #         print(" "*start + read[start:stop])
                #     if pl_pos_keys[i] == "threep1":
                #         print(" "*start + read[start:stop])
                #     print(" "*pos + pl_pos_keys[i])
                # print(" "*win_l + "."*(win_r - win_l))
                for (i_c, win_key) in enumerate(pl_win_keys):
                    first_position = pos_zips[0]
                    shift = i_c/2
                    i_r_list = [i + i_c/2 for i in full_inds]
                    left_nas = [float('nan')]*max(0, -1*min(i_r_list))
                    pl_left_bound = max(0, min(i_r_list))
                    pl_right_bound = min(plfold.shape[0], max(i_r_list) + 1)
                    # print(full_inds)
                    # print(pl_left_bound)
                    # print(pl_right_bound)
                    # print(plfold.shape)
                    pl_row = list(plfold.iloc[pl_left_bound:pl_right_bound,i_c])
                    logpl_row = list(logplfold.iloc[pl_left_bound:pl_right_bound,i_c])
                    right_nas = [float('nan')]*max(0, max(i_r_list)-plfold.shape[0]+1)
                    output_row = left_nas + pl_row + right_nas
                    logoutput_row = left_nas + logpl_row + right_nas
                    # print(len(left_nas))
                    # print(len(right_nas))
                    # if flank == "CC.AG":
                    #     print(win_key)
                    #     print(len(output_row))
                    # print(len(pl_pos_keys))
                    for (pos_key, out, logout) in zip(pl_pos_keys, output_row, logoutput_row):
                        # if win_key == "w1":
                        #     print("%s\t%s" %(pos_key, out))
                        # print(i_r)
                        # print(plfold.shape)
                        # print(" "*(ind + shift) + "*" + " "*(len(read)-ind - shift -1) + pos)
                        pos_win_plfold_map[pos_key][win_key].append(out)
                        pos_win_logplfold_map[pos_key][win_key].append(logout)
    return [pos_win_plfold_map, pos_win_logplfold_map, flanks]

def main():
    time_start = time.time()
    # Define all the relevant arguments.
    arguments = ["miRNA", "experiment", "condition", "n_constant", "sitelist", "win_start", "win_stop", "win_max"]
    mirna, experiment, condition, n_constant, sitelist, win_start, win_stop, win_max = parse_arguments(arguments)
    # Load file with site types for the miRNA.
    sites_file_name = ("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/"
                       "AssignSiteTypes/sites.%s_%s.txt" % (mirna.split("-alt")[0], sitelist))
    # Initialize the dictionary with all sites, and two more dictionaries,
    # one mapping to the site name from the sequence, and one from site name
    # to the count.
    with open(sites_file_name) as file_in:
        sites = file_in.read().split("\n")
    order_site_map = {site: i for i, site in enumerate(sites)}
    with open('general/site_info.txt',"r+") as site_pos_file:
        site_pos_map = dict()
        for line in site_pos_file:
            line_strip = line.strip().split("\t")
            key = line_strip[0]
            if key in order_site_map.keys() and key != "None":
                value = int(line_strip[2]) + int(line_strip[3]) - 1
                site_pos_map[key] = value
    # Get the path to the read file and to that of where the site labels will
    # be written.
    if experiment in ["equil_mmseed_nb", "equil_seed_nb"]:
        read_length = 38
    else:
        read_length = 37    
    extension_sites = "_%s_%s" %(n_constant, sitelist)
    extension_struc = "_%s" %("Full")
    reads_path = get_analysis_path(mirna, experiment, condition,
                                   "full_reads")
    sites_path = get_analysis_path(mirna, experiment, condition,
                                   "sites", ext=extension_sites)
    with open(reads_path, "rb") as file_in_reads:
        with open(sites_path, "rb") as file_in_sites:
            cwd = os.getcwd()
            folder_name = randomword(20)
            os.makedirs("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/AnalyzeStructure/RNAplfold_temp/%s" %(folder_name))
            os.chdir("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/AnalyzeStructure/RNAplfold_temp/%s" %(folder_name))

            results = multiprocess_test([file_in_reads, file_in_sites],
                                        readline_two,
                                        100000,
                                        get_read_structural_data,
                                        20000,
                                        order_site_map,
                                        site_pos_map,
                                        int(n_constant),
                                        read_length,
                                        mirna,
                                        experiment,
                                        int(win_start),
                                        int(win_stop),
                                        int(win_max))
            os.chdir(cwd)
            os.rmdir("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/AnalyzeStructure/RNAplfold_temp/%s" %(folder_name))
    # pl_pos_5p = ["f%s" %(i) for i in range(1, win_stop - 8 + 1)]
    # pl_pos_3p = ["t%s" %(i) for i in range(1, -1*win_start + 1)]
    # print(pl_pos_3p)
    # pl_pos_seed = ["s%s" %(i) for i in range(1, 8 + 1)[::-1]]
    # pl_pos_keys = pl_pos_5p + pl_pos_seed + pl_pos_3p
    # pl_win_keys = ["w%s" %(i) for i in range(1, win+1)]



    pl_pos_5p = ["f%s" %(i) for i in range(1, int(win_stop) - 8 + 1)]
    pl_pos_3p = ["t%s" %(i) for i in range(1, -1*int(win_start) + 1)]
    pl_pos_seed = ["s%s" %(i) for i in range(1, 8 + 1)[::-1]]
    pl_pos = pl_pos_5p + pl_pos_seed + pl_pos_3p
    pl_win = ["w%s" %(i) for i in range(1, int(win_max) + 1)]

    keys_new = list(it.product(pl_pos, pl_win))
    key_names_new = ["|".join(i) for i in keys_new]
    flanks_list = [k for j in results for k in j[2]]
    indeces = ["flanks"] + key_names_new
    file_ext = "_%s_%s-%s-%s" %(n_constant, win_start, win_stop, win_max)
    for i_out, outname in enumerate(["plfold", "logplfold"]):
        out_dict = {"|".join(i) : [k for j in results for k in j[i_out][i[0]][i[1]]] for i in keys_new} 
        out_df = pd.DataFrame.from_dict(out_dict)[key_names_new]
        out_df["flanks"] = flanks_list
        out_df = out_df[indeces]
        print(out_df["f1|w22"])
        print(out_df.iloc[:10, :10])
        out_path = get_analysis_path(mirna, experiment, condition, "%s_PAPER/%s" %(outname, "8mer"),ext=file_ext)
        print(out_path)
        out_df.to_csv(out_path,sep="\t", header=True, index=False, na_rep="NaN")
 
    # plfold_dict 
    # logplfold_dict = {"|".join(i) : [k for j in results for k in j[1][i[0]][i[1]]] for i in keys_new}

    # logplfold_df = pd.DataFrame.from_dict(logplfold_dict)[key_names_new]
    # logplfold_df["flanks"] = flanks_list
    # logplfold_df = logplfold_df[indeces]
    # print(plfold_df)
    # print(logplfold_df)
    # logplfold_path = get_analysis_path(mirna,experiment,condition,"logplfold_PAPER/%s" %("8mer"),ext=file_ext)
    # print(plfold_path)
    # print(logplfold_path)
    # logplfold_df.to_csv(logplfold_path,sep="\t", header=True, index=False)



    # flanks = get_kmer_list(4)
    # output = {pos: {win: [j for i in results for j in i[pos][win]] for pos in results[0].values()[0]} for pos in ["8mer"] if site != "None"}

  
    # summary = {site: {flank: "0" for flank in flanks} for site in order_site_map.keys() if site != "None"}
    # for site in ["8mer"]:
        # with open(site_flank_path, "wb") as file_flank:
        #     file_flank.write("\t".join(["Flank", "plfold", "accessibility", "AU_cs", "AU_win", "AU_read", "AU_win_wo", "AU_read_wo"])+"\n")
        # with open(site_flank_path, "ab") as file_flank:
        #     output_rows = "\n".join([
        #         "\n".join([
        #             "\t".join([flank]+[str(j) for j in i]) for i
        #              in output[site][flank]
        #          ]) for flank in flanks if len(output[site][flank]) > 0])
        #     file_flank.write(output_rows)
        # print(site_flank_path)


    #         # out = pd.DataFrame(output_structures[site][flank],columns=["sa", "H"])
    #         # n = out.size
    #         # if n > 0:
    #         #     summary[site][flank] = ",".join([str(i) for i in [n, float(out.mean()), float(out.var())]])
    #         # out.to_csv(site_flank_path,sep="\t",index = False, header = False)

    #         # site_flank_AU_win_path = get_analysis_path(mirna,experiment,condition,"AU_content_PAPER/%s/%s" %(site,flank),ext=file_ext)
    #         # out = pd.DataFrame(output_AU_win[site][flank])
    #         # out.to_csv(site_flank_AU_win_path,sep="\t",index = False, header = False)

    #         # site_flank_AU_read_path = get_analysis_path(mirna,experiment,condition,"AU_contentinread_PAPER/%s/%s" %(site,flank),ext=file_ext)
    #         # out = pd.DataFrame(output_AU_read[site][flank])
    #         # out.to_csv(site_flank_AU_read_path,sep="\t",index = False, header = False)

    #         # site_flank_AU_context_path = get_analysis_path(mirna,experiment,condition,"AU_contextscore_PAPER/%s/%s" %(site,flank),ext=file_ext)
    #         # out = pd.DataFrame(output_AU_cs[site][flank])
    #         # out.to_csv(site_flank_AU_read_path,sep="\t",index = False, header = False)



    # # Print the amount of time the script took to complete.
    # print_time_elapsed(time_start)

################################################################################

if __name__ == "__main__":
    main()

