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



def RNAplfold(read, file_num):
    # This should give the probability of the region in the read complementary
    # to position "mir_start" to "mir_stop" in the miRNA sequence.


    read_name = randomword(20)
    temp_fa_filename = "%s%s.fa" %(read_name, file_num)

    with open(temp_fa_filename, "wb") as f:
        f.write(">%s%s\n%s" %(read_name, file_num, read))

    # call RNAplfold:
    l_param = len(read)
    w_param = l_param

    call = 'RNAplfold -L %s -W %s -u 20 < %s%s.fa' %(l_param, w_param, 
                                                     read_name, file_num)
    subprocess.call([call], shell=True, stdout = subprocess.PIPE)

    temp_lunp_filename = '%s%s_lunp' %(read_name, file_num)
    rnaplfold_data = pd.read_csv(temp_lunp_filename, sep = '\t', header = 1).set_index(' #i$')
    os.remove(temp_fa_filename)
    os.remove(temp_lunp_filename)
    os.remove('%s%s_dp.ps' %(read_name, file_num))


    return rnaplfold_data



def get_plfold_data(read_seqs, order_site_map, site_pos_map, n_constant,
                             read_length, mirna, experiment):
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
    labels_left = ["l%s" %(i) for i in range(15,0,-1)]
    labels_mid = ["seed%s" %(i) for i in range(7,1,-1)]
    labels_right = ["r%s" %(i) for i in range(1,16)]
    labels_all = labels_left + labels_mid + labels_right
    heights = ["len%s" %(i) for i in range(1, 21)]

    site_flanks_map = {"".join(kmer): {
        j : pd.DataFrame(0, index=heights, columns=labels_all) for j in ["N", "lin", "log"]} for kmer 
        in list(it.product(["A","C","G","T"],repeat=4))} 
    print(site_flanks_map["AAAA"])
    # Pre-allocate output with "None" for each line.
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
            if start in sites_range and stop in sites_range and site in ["8mer"]:
                print(read)
                mirp1_span = site_pos_map[site] - 1
                mirp1 = start + mirp1_span
                print(" "*mirp1+"*")
                print(site_flanks_map.columns())

                if "b" in site:
                    win_l -= 1

                # Perform plfold on the entire read:
                plfold = RNAplfold(read, num_i)
                print(plfold)

                df_cols = plfold.shape[1]

                # Get average secondary structure
                p_acc_pl = plfold.iloc[win_r-1, df_cols-2]**(1/(float(win_r-win_l)))
                p_acc_geo = 10**np.mean(
                    plfold.iloc[win_l:win_r, 0].apply(math.log10))

                # Get the flanking nucleotides:
                l_1, l_2 = read[start-2], read[start-1]
                r_1, r_2 = read[stop], read[stop+1]
                flank = "".join([l_1, l_2, r_1, r_2])
                print(site_flanks_map[flank])

    


                site_flanks_map[site][flank].append((p_acc_pl, p_acc_geo, AU_cs, AU_win, AU_read, AU_win_wo, AU_read_wo))
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

    with open(reads_path, "rb") as file_in_reads:
        with open(sites_path, "rb") as file_in_sites:
            cwd = os.getcwd()
            folder_name = randomword(20)
            os.makedirs("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/AnalyzeStructure/RNAplfold_temp/%s" %(folder_name))
            os.chdir("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/AnalyzeStructure/RNAplfold_temp/%s" %(folder_name))

            results = multiprocess_test([file_in_reads, file_in_sites],
                                        readline_two,
                                        100000,
                                        get_plfold_data,
                                        10000,
                                        order_site_map,
                                        site_pos_map,
                                        int(n_constant),
                                        read_length,
                                        mirna,
                                        experiment)
            os.chdir(cwd)
            os.rmdir("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/AnalyzeStructure/RNAplfold_temp/%s" %(folder_name))
    flanks = ["".join(flank)for flank in list(it.product(["A","C","G","T"],repeat=4))]
    output = {site: {flank: [j for i in results for j in i[site][flank]] for flank in flanks} for site in ["8mer"] if site != "None"}

  
    # summary = {site: {flank: "0" for flank in flanks} for site in order_site_map.keys() if site != "None"}
    # file_ext = "_%s_%s-%s" %(n_constant,win_start,win_stop)
    # for site in ["8mer"]:
    #     site_flank_path = get_analysis_path(mirna,experiment,condition,"plfoldwindowanalysis_PAPER/%s" %(site),ext=file_ext)
    #     with open(site_flank_path, "wb") as file_flank:
    #         file_flank.write("\t".join(["Flank", "plfold", "accessibility", "AU_cs", "AU_win", "AU_read", "AU_win_wo", "AU_read_wo"])+"\n")
    #     with open(site_flank_path, "ab") as file_flank:
    #         output_rows = "\n".join([
    #             "\n".join([
    #                 "\t".join([flank]+[str(j) for j in i]) for i
    #                  in output[site][flank]
    #              ]) for flank in flanks if len(output[site][flank]) > 0])
    #         file_flank.write(output_rows)
    #     print(site_flank_path)



    # Print the amount of time the script took to complete.
    print_time_elapsed(time_start)

################################################################################

if __name__ == "__main__":
    main()

