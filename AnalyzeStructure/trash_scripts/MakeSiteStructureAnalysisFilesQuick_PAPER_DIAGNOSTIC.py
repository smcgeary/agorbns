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
from general import parse_arguments, print_time_elapsed, get_analysis_path, multiprocess_file, multiprocess_test, readline_two
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




def RNAplfold(read,mir_start, mir_stop,file_num):
    # This should give the probability of the region in the read complementary
    # to position "mir_start" to "mir_stop" in the miRNA sequence.


    read_name = randomword(20)
    temp_fa_filename = "%s%s.fa" %(read_name, file_num)

    with open(temp_fa_filename, "wb") as f:
        f.write(">%s%s\n%s" %(read_name, file_num, read))

    # call RNAplfold:
    l_param = len(read)
    w_param = l_param
    u_param = mir_stop - mir_start + 1

    call = 'RNAplfold -L %s -W %s -u %s < %s%s.fa' %(l_param, w_param, u_param,
                                                     read_name, file_num)
    subprocess.call([call], shell=True, stdout = subprocess.PIPE)

    temp_lunp_filename = '%s%s_lunp' %(read_name, file_num)
    rnaplfold_data = pd.read_csv(temp_lunp_filename, sep = '\t', header = 1).set_index(' #i$')
    os.remove(temp_fa_filename)
    os.remove(temp_lunp_filename)
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
    print(site_type)
    print(utr)
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
    print("_"*l + upstream_str + "."*(site_end - site_start) + downstream_str)
    print(" "*l + "".join([str(i) for i in upstream]) + "."*(site_end - site_start) + "".join([str(i) for i in downstream]))
    for num, weight in enumerate(inv_upweights):
        print(" "*(l + num) + str(weight) + " "*(30 + site_end - site_start - len(str(weight))) + str(inv_downweights[num]))

    weighted = np.dot(upstream, upweights) + np.dot(downstream, downweights)
    total = float(sum(upweights) + sum(downweights))

    return weighted / total




def get_local_au_score(read, start, stop, sitetype):
    if sitetype not in ["8mer", "7mer-m8", "7mer-A1", "6mer"]:
        return False
    else:
        upstream_site = int(sitetype not in ["7mer-m8", "8mer"])
        upstream_seq = read[max(0, start - 30): start]

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

        # print(" "*(max(0, start - 30)) + upstream_seq + " "*(stop - start) + downstream_seq)
        # print(" "*(max(0, start - 30)) + "".join([str(i) for i in upstream]) + " "*(stop - start) + "".join([str(i) for i in downstream]))
        # for i, weight in enumerate(upweights):
        #     if upstream[i] == 1:
        #         print(" "*(max(0, start - 30) + i) + str(weight))
        #     else:
        #         print("")
        # print(" "*(max(0, start - 30)) + upstream_seq + " "*(stop - start) + downstream_seq)
        # for i, weight in enumerate(downweights):
        #     if downstream[i] == 1:
        #         print(" "*(stop + i) + str(weight))
        #     else:
        #         print("")
        # print("kathy_AU_left:")
        # print(np.dot(upstream, upweights))
        # print("kathy_AU_right:")
        # print(np.dot(downstream, downweights))
        weighted = np.dot(upstream, upweights) + np.dot(downstream, downweights)
        # if sitetype in ["8mer", "7mer-m8", "7mer-A1", "6mer"]:
        #     print("downstream score:")
        #     print(np.dot(downstream, downweights))
        #     print("upstream score:")
        #     print(np.dot(upstream, upweights))
        #     print("upstream total:")
        #     print(sum(upweights_total))
        #     print("downstream total:")
        #     print(sum(downweights_total))
        # print("up_totals")
        # print(sum(upweights_total))
        # print("down_totals:")
        # print(sum(downweights_total))
        total = float(sum(upweights_total) + sum(downweights_total))

        return(weighted / total)
    return

def get_read_structural_data(read_seqs, order_site_map, site_pos_map, n_constant,
                             read_length, mirna, experiment, win_start = -1,
                             win_stop = 15, constant = True):
    """Takes a read sequence and identifies the position of all site types'.

    Args:
        read_sequence: The read sequence.
        sites: A list of sequences corresponding to each site type.

    Returns:
        A tuple containing the list of reads and a dictionary of read types.
    """

    # Initialize count dictionary with keys as site names and a value of
    # for each item.
    sites_range = range(26-n_constant,26+37+n_constant+1)

    site_flanks_map = {site: {"".join(kmer): [] for kmer 
        in list(it.product(["A","C","G","T"],repeat=4))} 
        for site in ["8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-A1", "6mer-m8"]}

    # Pre-allocate output with "None" for each line.
    time_start = time.time()
    tick = 0
    tick = 0
    for num_i, i in enumerate(read_seqs):

        tick +=1
        # if tick % 10000 == 0:
        #     print(tick)
        #     print_time_elapsed(time_start)
        #     sys.stdout.flush()

        read, sites = (j.strip() for j in i)
        # Define the barcode portion of the read.
        barcode = read[26 + read_length : 26 + read_length + 3]
        if read == "GGGCAGAGTTCTACAGTCCGACGATCACTTGTGCTATGCTTGGGAACAACGTGGAAATGGCAGTCGTATGCCGTCTTCTGCTTG":
            print("FOUND IT %s" %(num_i))
            print(sites)
            print(read)
        # Deals with the TCG instead of TGT ending for miR-1 equilibrium exp
        if barcode == "TCG":
            barcode = "TGT"
            if mirna != "miR-1" or experiment != "equilibrium":
                read = read[:26 + 37]+"TGTTCGTATGCCGTCTTCTGCTTG"
        if barcode == "TGT" and mirna == "miR-1" and experiment == "equilibrium":
            read = read[:26 + 37]+"TCGTATGCCGTCTTCTGCTTG"


        # Take the constant regiong + the number of constant sequences

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
            l_1, l_2 = read[start-2], read[start-1]
            r_1, r_2 = read[stop], read[stop+1]
            flank = "".join([l_1, l_2, r_1, r_2])

            if start in sites_range and stop in sites_range and site == "8mer" and flank == "AAAA":
                mirp1_span = site_pos_map[site] - 1
                mirp1 = start + mirp1_span

                win_r = mirp1 - win_start + 1 + 1 
                win_l = mirp1 - win_stop + 1
                if "b" in site:
                    win_l -= 1

                # Perform plfold on the entire read:
                plfold = RNAplfold(read,win_start, win_stop, num_i)
                df_cols = plfold.shape[1]

                # Get average secondary structure
                p_acc_pl = plfold.iloc[win_r-1, df_cols-2]**(1/(float(win_r-win_l)))
                p_acc_geo = 10**np.mean(
                    plfold.iloc[win_l:win_r, 0].apply(math.log10))

                # Get the flanking nucleotides:
                l_1, l_2 = read[start-2], read[start-1]
                r_1, r_2 = read[stop], read[stop+1]
                flank = "".join([l_1, l_2, r_1, r_2])

                AU_cs = get_local_au_score(read, start, stop, site)


                AU = [int(x in ["A", "T"]) for x in read]
                AU_win = AU[win_l:win_r]
                AU_read = AU[26 : 26 + 37]
                AU_win_wo_site = [AU[i_au] for i_au in range(win_l, win_r)
                                  if i_au not in range(start, stop)]
                AU_read_wo_site = [AU[i_au] for i_au in range(26,26+37)
                                  if i_au not in range(start, stop)]
                AU_f_win, AU_f_read, AU_f_win_wo_site, AU_f_read_wo_site = [
                    float(sum(i))/len(i) for i in [
                        AU_win, AU_read, AU_win_wo_site, AU_read_wo_site]]

                if site in ["8mer"]:
                    print("%s\t%s\t%s" %(num_i,p_acc_geo, AU_f_win))
                    print(read)
                    print(" "*start + read[start:])
                # site_flanks_map[site][flank].append((pl_accessibility, AU_context_score, AU_win, AU_read))
        num_i +=1
    return site_flanks_map





def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["miRNA", "experiment", "condition", "n_constant", "sitelist"]
    mirna, experiment, condition, n_constant, sitelist = parse_arguments(arguments)

    # Load file with site types for the miRNA.
    sites_file_name = ("/lab/bartel1_ata/mcgeary/computation/"
                       "AgoRBNS/AssignSiteTypes/sites.%s_%s.txt" % (mirna.split("-alt")[0], sitelist))


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
                                   "sites",ext=extension_sites)

    constant = True
    win_start = 1
    win_stop = 15
    with open(reads_path, "rb") as file_in_reads:
        with open(sites_path, "rb") as file_in_sites:
            cwd = os.getcwd()
            folder_name = randomword(20)
            os.makedirs("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/AnalyzeStructure/RNAplfold_temp/%s" %(folder_name))
            os.chdir("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/AnalyzeStructure/RNAplfold_temp/%s" %(folder_name))

            results = multiprocess_file([file_in_reads, file_in_sites],
                                        readline_two,10000,
                                        get_read_structural_data,
                                        order_site_map,
                                        site_pos_map,
                                        int(n_constant),
                                        read_length,
                                        mirna,
                                        experiment,
                                        win_start,
                                        win_stop,
                                        constant)
            os.chdir(cwd)
            os.rmdir("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/AnalyzeStructure/RNAplfold_temp/%s" %(folder_name))
    # flanks = ["".join(flank)for flank in list(it.product(["A","C","G","T"],repeat=4))]
    # output = {site: {flank: [j for i in results for j in i[site][flank]] for flank in flanks} for site in ["8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-A1", "6mer-m8"] if site != "None"}

  
    # # summary = {site: {flank: "0" for flank in flanks} for site in order_site_map.keys() if site != "None"}
    # file_ext = "_%s_%s-%s" %(n_constant,win_start,win_stop)
    # if constant == False:
    #     file_ext = file_ext + "_noconstant"
    # for site in ["8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-A1", "6mer-m8"]:
    #     site_flank_path = get_analysis_path(mirna,experiment,condition,"structural_analysis_PAPER_absolute/%s" %(site),ext=file_ext)
    #     with open(site_flank_path, "wb") as file_flank:
    #         file_flank.write("\t".join(["Flank", "plfold","AU_cs","AU_win", "AU_read"])+"\n")
    #     with open(site_flank_path, "ab") as file_flank:
    #         output_rows = "\n".join([
    #             "\n".join([
    #                 "\t".join([flank]+[str(j) for j in i]) for i
    #                  in output[site][flank]
    #              ]) for flank in flanks if len(output[site][flank]) > 0])
    #         file_flank.write(output_rows)
    #     print(site_flank_path)


            # out = pd.DataFrame(output_structures[site][flank],columns=["sa", "H"])
            # n = out.size
            # if n > 0:
            #     summary[site][flank] = ",".join([str(i) for i in [n, float(out.mean()), float(out.var())]])
            # out.to_csv(site_flank_path,sep="\t",index = False, header = False)

            # site_flank_AU_win_path = get_analysis_path(mirna,experiment,condition,"AU_content_PAPER/%s/%s" %(site,flank),ext=file_ext)
            # out = pd.DataFrame(output_AU_win[site][flank])
            # out.to_csv(site_flank_AU_win_path,sep="\t",index = False, header = False)

            # site_flank_AU_read_path = get_analysis_path(mirna,experiment,condition,"AU_contentinread_PAPER/%s/%s" %(site,flank),ext=file_ext)
            # out = pd.DataFrame(output_AU_read[site][flank])
            # out.to_csv(site_flank_AU_read_path,sep="\t",index = False, header = False)

            # site_flank_AU_context_path = get_analysis_path(mirna,experiment,condition,"AU_contextscore_PAPER/%s/%s" %(site,flank),ext=file_ext)
            # out = pd.DataFrame(output_AU_cs[site][flank])
            # out.to_csv(site_flank_AU_read_path,sep="\t",index = False, header = False)



    # Print the amount of time the script took to complete.
    print_time_elapsed(time_start)

################################################################################

if __name__ == "__main__":
    main()

