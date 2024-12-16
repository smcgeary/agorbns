################################################################################
#GenerateSiteTypeKds.py
################################################################################
import imp # Used to import general.py
import time # Used for time
import numpy as np
import pandas as pd
import math
from scipy import optimize
from scipy import special
import csv
import itertools as it

# from special import gammaln
# from special import xlogy
pd.set_option('display.width', 1000)
imp.load_source("general",
                ("/lab/bartel1_ata/mcgeary/computation/"
                 "AgoRBNS/general/general.py"
                )
               )
from general import parse_arguments, print_time_elapsed, get_analysis_path, ensure_directory

# This script reads one read file and outputs multiple .txt files to be used in all
# downstream applications. Inputs are:

# mirna: The miRNA in the experiment to be processed.
# experiment: The type of experiment to be processed.
# condition: The sampel type within the experiment to be processed.

# Functions from general this uses:


# Constants

# FUNCTIONS

def get_flanking_table(mirna,site,sitelist,start,stop):
    # Load the list of sites
    # Determine the order of the sites, to pick the the one at the top of list.
    # Initialize the the counts dictionary.
 
    # Initialize the columns / data headings, and the matrix.
    conditions = ["I", "I_combined", "40", "12.6", "4", "1.26", "0.4", "0"]
    flanks = ["".join(kmer[:2] + (".", ) + kmer[2:])
              for kmer in it.product(["A", "C", "G", "T"], repeat=4)]
    flanks_single = ["".join(i) for i in it.product(["A", "C", "G", "T"], repeat=2)]
    matrix = pd.DataFrame(0, index=flanks, columns=conditions)
    lmatrix = pd.DataFrame(0, index=flanks_single, columns=conditions)
    rmatrix = pd.DataFrame(0, index=flanks_single, columns=conditions)
    for cond in conditions:

        counts_flank_sites_map = {flank: 0 for flank in flanks}
        counts_lflank_sites_map = {flank: 0 for flank in flanks_single}
        counts_rflank_sites_map = {flank: 0 for flank in flanks_single}

    	site_cond_name = ("/lab/solexa_bartel/mcgeary/AgoRBNS/%s"
    					  "/equilibrium/full_flanking/%s_%s-%s_%s.txt" % (mirna, cond, start, stop, sitelist))

        with open(site_cond_name, "rb") as file_in:
            temp = enumerate(file_in.readline().strip().split("\t"))
            for i, t in temp:
                print(i)
                print(t)
                print(site)
                if site == t:
                    ind = int(i)
            # [ind] = [int(i) for i, s in temp
            #             if s == site]
            print(ind)
            for line in file_in:
                all = line.strip().split("\t")
                flank = all[0].split(":")[0]
                lflank, rflank = flank.split(".")
                val = int(all[ind + 1])
                counts_flank_sites_map[flank] += val
                counts_lflank_sites_map[lflank] += val
                counts_rflank_sites_map[rflank] += val
        matrix[cond] = [counts_flank_sites_map[i] for i in flanks]
        lmatrix[cond] = [counts_lflank_sites_map[i] for i in flanks_single]
        rmatrix[cond] = [counts_rflank_sites_map[i] for i in flanks_single]      
        print(matrix)  
    return([matrix, lmatrix, rmatrix])

def main():
    time_start = time.time()

    # Define all the relevant arguments.

    arguments = ["miRNA", "sitelist","start", "stop"]
    [mirna, sitelist,start, stop] = parse_arguments(arguments)
    if sitelist:
        print("hi")
        sites_file_name = ("/lab/bartel1_ata/mcgeary/computation/"
                           "AgoRBNS/AssignSiteTypes/sites.%s_%s.txt" % (mirna.split("-alt")[0], sitelist))
    else:
        sites_file_name = ("/lab/bartel1_ata/mcgeary/computation/"
                   "AgoRBNS/AssignSiteTypes/sites.%s.txt" % (mirna.split("-alt")[0]))


    print(sitelist)
    with open(sites_file_name) as file_in:
        sites = file_in.read().split("\n")
    print(sites)
    add = True
    sites_new = []
    for i in sites:
        if "Centered" not in i:
            sites_new.append(i)
        elif add:
            sites_new.append("Centered")
            add = False
    sites = sites_new
    for site in sites:
        print(site)
        print(site)

        sfXc, slfXc, srfXc = get_flanking_table(mirna,site,sitelist,start,stop)
        sfXc.columns=["I", "I_combined", "A40","A12.65","A4","A1.265","A0.4","A0"]
        slfXc.columns = sfXc.columns
        srfXc.columns = sfXc.columns

        sfXc = sfXc[(sfXc.T != 0).any()]
        slfXc = slfXc[(slfXc.T != 0).any()]
        srfXc = srfXc[(srfXc.T != 0).any()]

        site_count_tables_path = get_analysis_path(
            mirna,"equilibrium", "%s_flanking" % (site), "full_site_count_tables", ext = "_%s-%s_%s" %(start, stop, sitelist))
        site_left_count_tables_path = get_analysis_path(
            mirna,"equilibrium", "%s_leftflanking" % (site), "full_site_count_tables", ext = "_%s-%s_%s" %(start, stop, sitelist))
        site_right_count_tables_path = get_analysis_path(
            mirna,"equilibrium", "%s_rightflanking" % (site), "full_site_count_tables", ext = "_%s-%s_%s" %(start, stop, sitelist))       
        print(site_count_tables_path)
        with open(site_count_tables_path,"wb") as file_out:
            sfXc.to_csv(site_count_tables_path,sep="\t")
        with open(site_left_count_tables_path,"wb") as file_out:
            slfXc.to_csv(site_left_count_tables_path,sep="\t")
        with open(site_right_count_tables_path,"wb") as file_out:
            srfXc.to_csv(site_right_count_tables_path,sep="\t")

    print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()
