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
pd.set_option('display.width', 1000)
imp.load_source("general",
                ("/lab/bartel1_ata/mcgeary/computation/"
                 "AgoRBNS/general/general.py"
                )
               )
imp.load_source("sitetools",
                ("/lab/bartel1_ata/mcgeary/computation/"
                 "AgoRBNS/AssignSiteTypes/sitetypes.py"
                )
               )
from general import parse_arguments
from general import print_time_elapsed
from general import get_analysis_path
from general import ensure_directory
from sitetools import get_site_seq
from sitetools import get_seq_site_map

# This script reads one read file and outputs multiple .txt files to be used in all
# downstream applications. Inputs are:

# mirna: The miRNA in the experiment to be processed.
# experiment: The type of experiment to be processed.
# condition: The sampel type within the experiment to be processed.

# Functions from general this uses:


# Constants

# FUNCTIONS

def get_site_table(mirna, exp, start, stop, sitelist):
    # Load the list of sites
    sites_file_name = ("/lab/bartel1_ata/mcgeary/computation/"
                       "AgoRBNS/AssignSiteTypes/sites.%s_%s.txt" % (mirna, sitelist))
    # This ensures that both the final data file and the place where the figure will go are made.
    ensure_directory("/lab/solexa_bartel/mcgeary/AgoRBNS/%s/%s/kds" % (mirna, exp)) 
    ensure_directory("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/figures/kds/%s/%s" % (mirna, exp)) 

    with open(sites_file_name) as file_in:
        sites = file_in.read().split("\n")+["None"]
    print(sites)
    site_seq_map = get_seq_site_map(mirna, sitelist)
    seqs = [site_seq_map[site] for site in sites if site != "None"]+["None"]
    # This swaps the four "Centered-..." for one single "Centered" entry in the list.
    add = True
    sites_new = []
    for i in sites:
        if "Centered" not in i:
            sites_new.append(i)
        elif add:
            sites_new.append("Centered")
            add = False
    sites = sites_new


    # Determine the order of the sites, to pick the the one at the top of list.
    order_site_map = {site: i for i, site in enumerate(sites)}
    # Initialize the the counts dictionary.
 
    # Initialize the columns / data headings, and the matrix.
    if exp == "equil_seed_nb":
        conditions = ["I", "12.6", "12.6_2", "4", "4_2", "1.26", "0.4", "0"]
    else:
        conditions = ["I", "40", "12.6", "4", "1.26", "0.4", "0"]
    matrix = pd.DataFrame(0,index=sites,columns=["seq"] + conditions)
    matrix["seq"] = seqs
    for cond in conditions:

        counts_sites_map = {site: 0 for site in sites}

    	site_cond_name = ("/lab/solexa_bartel/mcgeary/AgoRBNS/%s"
    					  "/%s/full_multisite_counts/%s_%s-%s_%s.txt" % (mirna, exp, cond, start, stop, sitelist))
        with open(site_cond_name) as file_in:
    		multisites = [i.split(":\t")[0].split(",") for i in file_in.readlines()]
        with open(site_cond_name) as file_in:
            counts = [int(i.split(":\t")[1]) for i in file_in.readlines()]
        for j, multisite in enumerate(multisites):
            ranks = [order_site_map[k] for k in multisite]
            [site] = [k for k in multisite if order_site_map[k] == min(ranks)]
            counts_sites_map[site] += counts[j]
        matrix[cond] = [counts_sites_map[i] for i in sites]
    return(matrix)

def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["miRNA", "experiment", "sitelist", "start", "stop"]
    mirna, exp, sitelist, start, stop = parse_arguments(arguments)
    siteXcount = get_site_table(mirna, exp, start, stop, sitelist)
    print(siteXcount)
    if exp == "equil_seed_nb":
        siteXcount.columns = ["Seq", "I", "A12.6", "A12.6_2", "A4", "A4_2", "A1.26", "A0.4", "A0"]
    else:
        siteXcount.columns=["Seq","I","A40","A12.65","A4","A1.265","A0.4","A0"]
    site_count_tables_path = get_analysis_path(mirna, exp, "all_sites_%s-%s_%s" %(start, stop, sitelist),"full_site_count_tables")
    print(site_count_tables_path)
    with open(site_count_tables_path,"wb") as file_out:
        siteXcount.to_csv(site_count_tables_path,sep="\t")
    
    print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()
