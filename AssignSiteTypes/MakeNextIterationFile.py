################################################################################
#GenerateSiteTypeKds.py
################################################################################
import imp # Used to import general.py
import time # Used for time
import numpy as np
import pandas as pd
import math
import subprocess
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
from general import get_kmer_list
from sitetools import get_site_seq
from sitetools import get_seq_site_map

# This script reads one read file and outputs multiple .txt files to be used in all
# downstream applications. Inputs are:

# FUNCTIONS


def MakeNextIterationFile(mirna, experiment, condition, motif, n_iter):
    print(motif)
    file_in = ("/lab/solexa_bartel/mcgeary/AgoRBNS/%s/%s/full_reads/"
               "%s_iter%s.txt") %(mirna, experiment, condition, n_iter-1)
    print(file_in)

    file_out = ("/lab/solexa_bartel/mcgeary/AgoRBNS/%s/%s/full_reads/"
               "%s_iter%s.txt") %(mirna, experiment, condition, n_iter)
    print(file_out)

    call = "bsub 'grep -v %s %s > %s'" %(motif, file_in, file_out)
    print(call)

    subprocess.call([call], shell=True, stdout = subprocess.PIPE)

def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["miRNA", "experiment", "n_constant", "n_iter"]
    mirna, exp, n_constant, n_iter = parse_arguments(arguments)
    
    n_iter = int(n_iter)

    motifs_path = get_analysis_path(mirna, exp, "motifs_PAPER", "sitekmers")
    
    with open(motifs_path,"r") as file_in:
        motif = file_in.readlines()[n_iter].split("\t")[0]
    for condition in ["I", "I_combined", "0", "0.4", "1.26", "4", "12.6", "40"]:
        MakeNextIterationFile(mirna, exp, condition, motif, n_iter)

################################################################################

if __name__ == "__main__":
    main()
