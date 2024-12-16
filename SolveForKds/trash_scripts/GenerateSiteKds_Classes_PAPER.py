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
imp.load_source("RBNS_methods",
                ("/lab/bartel1_ata/mcgeary/computation/"
                 "AgoRBNS/general/RBNS_methods.py"
                )
               )
imp.load_source("sitetools",
                ("/lab/bartel1_ata/mcgeary/computation/"
                 "AgoRBNS/AssignSiteTypes/sitetypes.py"
                )
               )
from general import *
from sitetools import get_site_seq
from sitetools import get_seq_site_map
from RBNS_methods import *

# This script reads one read file and outputs multiple .txt files to be used in all
# downstream applications. Inputs are:

# FUNCTIONS

def get_site_table(_sitelist, experiment, n_constant):
    # Load the list of sites
    sites = _sitelist.list_site_names() + ["None"]
    mirna = _sitelist.mirna
    # Determine the order of the sites, to pick the the one at the top of list.
    if experiment == "equil_seed_nb":
        conditions = ["I", "12.6", "12.6_2", "4", "4_2", "1.26", "0.4", "0"]
    else:
        conditions = ["I", "I_combined", "40", "12.6", "4", "1.26", "0.4", "0"]
    sXc = pd.DataFrame(0, index=sites, columns=conditions)
    for cond in conditions:
    	cond_path = "%s%s/%s/site_counts_classes/%s_%s_%s.txt" %(
            SOLEXA_DIR, mirna, experiment, cond, n_constant, _sitelist.name
        )
        s_c = pd.read_csv(cond_path, header=None, index_col=0,sep="\t")
        sXc[cond] = s_c
        print(sXc)
    return(sXc)

def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["miRNA", "experiment", "n_constant", "sitelist"]
    mirna, experiment, n_constant, sitelist = parse_arguments(arguments)

    _mirna = Mirna(mirna)
    _sitelist = SiteList(_mirna, sitelist)
    print(_sitelist)
    sXc = get_site_table(_sitelist, experiment, n_constant)
    print(sXc)
    if experiment == "equil_seed_nb":
        sXc.columns = ["I", "A12.6", "A12.6_2", "A4", "A4_2", "A1.26", "A0.4", "A0"]
    else:
        sXc.columns=["I", "I_combined", "A40","A12.65","A4","A1.265","A0.4","A0"]
    sXc_path = get_analysis_path(mirna, experiment,
                                 "all_sites_%s_%s" %(n_constant, sitelist),
                                 "site_count_tables_classes")
    print(sXc)
    print(sXc_path)
    sXc.to_csv(sXc_path, sep="\t")
    
    ensure_directory("%s%s/%s/kds" % (SOLEXA_DIR, mirna, experiment)) 
    ensure_directory("%sfigures/kds/%s/%s" % (HOME_DIR, mirna, experiment)) 


    print_time_elapsed(time_start)

################################################################################

if __name__ == "__main__":
    main()
