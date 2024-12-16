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
from general import get_kmer_list_constant_insert
from sitetools import get_site_seq
from sitetools import get_seq_site_map


# FUNCTIONS

def get_site_table(mirna, exp, start, stop, mir_start, mir_stop, kmer_list):
    # Load the list of sites
    # This swaps the four "Centered-..." for one single "Centered" entry in the list.


    # Initialize the columns / data headings, and the matrix.
    if exp == "equil_mmseed_nb":
        conditions = ["I", "12.6", "12.6_2", "4", "4_2", "1.26", "0.4", "0"]
    else:
        conditions = ["I", "40", "12.6", "4", "1.26", "0.4", "0"]
    matrix = pd.DataFrame(0, index = kmer_list, columns = conditions)
    for cond in conditions:

        counts_sites_map = {kmer: 0 for kmer in kmer_list}

    	site_cond_name = ("/lab/solexa_bartel/mcgeary/AgoRBNS/%s"
    					  "/%s/full_multisite_counts/%s_%s-%s_12mers_%s-%s.txt" % (mirna, exp, cond, start, stop, mir_start, mir_stop))
        with open(site_cond_name) as file_in:
    		multisites = [i.split(":\t")[0].split(",") for i in file_in.readlines()]
        with open(site_cond_name) as file_in:
            counts = [int(i.split(":\t")[1]) for i in file_in.readlines()]
        for j, multisite in enumerate(multisites):
            sites = [k for k in multisite]
            for site in sites:
                counts_sites_map[site] += counts[j]/float(len(multisite))
        matrix[cond] = [counts_sites_map[kmer] for kmer in kmer_list]
    return(matrix)

def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["miRNA", "experiment", "start", "stop", "miR_start", "miR_stop"]
    mirna, exp, start, stop, mir_start, mir_stop = parse_arguments(arguments)


    kmer_list = get_kmer_list_constant_insert(mirna, int(mir_start),
                                              int(mir_stop))
    kmer_list.append("None")

    siteXcount = get_site_table(mirna, exp, start, stop, int(mir_start),
                                int(mir_stop), kmer_list)
    print(siteXcount)
    if exp == "equil_mmseed_nb":
        siteXcount.columns = ["I", "A12.6", "A12.6_2", "A4", "A4_2", "A1.26", "A0.4", "A0"]
    else:
        siteXcount.columns=["I", "A40", "A12.65", "A4", "A1.265",
                            "A0.4", "A0"]

    site_count_tables_path = get_analysis_path(
        mirna, exp, "all_sites_%s-%s_12mers_%s-%s_update1" %(start, stop, mir_start,
                                                     mir_stop),
        "full_site_count_tables"
    )
    print(siteXcount.sum(axis=0))
    print(site_count_tables_path)

    with open(site_count_tables_path, "wb") as file_out:
        siteXcount.to_csv(site_count_tables_path, sep = "\t")
    
    print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()
