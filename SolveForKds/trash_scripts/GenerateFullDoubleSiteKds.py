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

def get_site_table(mirna, start, stop):
    # Load the list of sites
    sites_file_name = ("/lab/bartel1_ata/mcgeary/computation/"
                       "AgoRBNS/AssignSiteTypes/sites.%s.txt" % (mirna))
    # This ensures that both the final data file and the place where the figure will go are made.
    ensure_directory("/lab/solexa_bartel/mcgeary/AgoRBNS/%s/equilibrium/kds" % (mirna)) 
    ensure_directory("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/figures/kds/%s" % (mirna)) 

    with open(sites_file_name) as file_in:
        sites = file_in.read().split("\n")+["None"]
    print(sites)
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
    sites_bipartite = []

    # Determine the order of the sites, to pick the the one at the top of list.
    order_site_map = {site: i for i, site in enumerate(sites)}
    print(order_site_map)
    # Initialize the the counts dictionary.
 
    # Initialize the columns / data headings, and the matrix.
    conditions = ["I", "40", "12.6", "4", "1.26", "0.4", "0"]
    for cond in conditions:
        site_cond_bipart_name = ("/lab/solexa_bartel/mcgeary/AgoRBNS/%s"
                                 "/equilibrium/full_double_bipartite_counts/%s_%s-%s.txt" % (mirna, cond, start, stop))

        with open(site_cond_bipart_name) as file:
            keys = [line.split(":\t")[0] for line in file]
            for key in keys:
                if key not in sites:
                    sites.append(key)
    print(sites)
    print(sites_bipartite)

    matrix = pd.DataFrame(0,index=sites,columns=conditions)
    for cond in conditions:

        counts_sites_map = {site: 0 for site in sites}

    	site_cond_name = ("/lab/solexa_bartel/mcgeary/AgoRBNS/%s"
    					  "/equilibrium/full_double_multisite_counts/%s_%s-%s.txt" % (mirna, cond, start, stop))
        site_cond_bipart_name = ("/lab/solexa_bartel/mcgeary/AgoRBNS/%s"
                                 "/equilibrium/full_double_bipartite_counts/%s_%s-%s.txt" % (mirna, cond, start, stop))

        with open(site_cond_name) as file_in:
    		multisites = [i.split(":\t")[0].split(",") for i in file_in.readlines()]
        with open(site_cond_name) as file_in:
            counts = [int(i.split(":\t")[1]) for i in file_in.readlines()]
        for j, multisite in enumerate(multisites):
            ranks = [order_site_map[k] for k in multisite]
            # print(multisite)
            # print(ranks)
            [site] = [k for k in multisite if order_site_map[k] == min(ranks)]
            counts_sites_map[site] += counts[j]
        with open(site_cond_bipart_name) as file_in:
            biparts = [i.split(":\t")[0] for i in file_in.readlines()]
        with open(site_cond_bipart_name) as file_in:
            counts = [sum([int(j) for j in i.split("\t")[1:]]) for i in file_in.readlines()]
        for j, bipart in enumerate(biparts):
            counts_sites_map[bipart] += counts[j]


        matrix[cond] = [counts_sites_map[i] for i in sites]
    

    return(matrix)



def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["miRNA","start","stop"]
    [mirna, start, stop] = parse_arguments(arguments)
    siteXcount = get_site_table(mirna, start, stop)
    print(siteXcount)
    siteXcount.columns=["I","A40","A12.65","A4","A1.265","A0.4","A0"]
    site_count_tables_path = get_analysis_path(mirna,"equilibrium","all_sites_%s-%s" %(start,stop),"full_double_site_count_tables")
    print(site_count_tables_path)
    with open(site_count_tables_path,"wb") as file_out:
        siteXcount.to_csv(site_count_tables_path,sep="\t")
    
    print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()
