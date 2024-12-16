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
from general import get_kmer_list
from sitetools import get_site_seq
from sitetools import get_seq_site_map

# This script reads one read file and outputs multiple .txt files to be used in all
# downstream applications. Inputs are:

# FUNCTIONS

def get_kmer_enrichments(mirna, exp, n_constant, kmer_length, n_iter):
    # Load the list of sites

 
    # Initialize the columns / data headings, and the matrix.

    conditions = ["I_combined", "4"]
    matrix1 = pd.DataFrame(0,index=get_kmer_list(kmer_length), columns=conditions)
    matrix2 = pd.DataFrame(0,index=get_kmer_list(kmer_length), columns=conditions)
    print(n_iter)
    for cond in conditions:
        print(cond)
        if n_iter == 0:
            site_cond_name = ("/lab/solexa_bartel/mcgeary/AgoRBNS/%s"
                              "/%s/kmer_counts/%s_%s_%s.txt" % (mirna, exp, cond, n_constant, kmer_length))
        else:
            site_cond_name = ("/lab/solexa_bartel/mcgeary/AgoRBNS/%s"
                              "/%s/kmer_counts/%s_%s_%s_iter%s.txt" % (mirna, exp, cond, n_constant, kmer_length, n_iter))
        with open(site_cond_name) as file_in:
            file_in.readline()
            for num, line in enumerate(file_in.readlines()):
                linesplit = line.strip().split("\t")
                kmer = linesplit[0]
                count1 = linesplit[1]
                count2 = linesplit[2]
                matrix1.at[kmer,cond] = float(count1)
                matrix2.at[kmer,cond] = float(count2)
    if kmer_length == 8:
        print(matrix1.loc["AACATTCC"])
        print(matrix1.sum(axis=0))
    matrix1 = matrix1.div(matrix1.sum(axis=0))
    matrix2 = matrix2.div(matrix2.sum(axis=0))
    mat_I = matrix1["I_combined"]
    if kmer_length == 8:
        print(matrix1.loc["AACATTCC"])

    mat_Ago1 = matrix1[["4"]].sum(axis=1)
    mat_Ago2 = matrix2[["4"]].sum(axis=1)
    R1 = pd.DataFrame(mat_Ago1 / mat_I)
    R2 = pd.DataFrame(mat_Ago2 / mat_I)
    R = pd.concat([R1, R2], axis = 1)
    R.columns = ["R","skaR"]
    return(R)


def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["miRNA", "experiment", "n_constant", "n_iter"]
    mirna, experiment, n_constant, n_iter = parse_arguments(arguments)
    n_iter = int(n_iter)
    for kmer_length in range(5,12):
        print(kmer_length)
        if n_iter > 0:
            R_value_path = get_analysis_path(mirna, experiment, str(kmer_length),"kmer_R_values",ext="_%s_iter%s" %(n_constant, n_iter))

        else:
            n_iter = False
            R_value_path = get_analysis_path(mirna, experiment, str(kmer_length),"kmer_R_values",ext="_%s" %(n_constant))
        print(R_value_path)

        print(n_iter)



        kmer_length = int(kmer_length)
        R_values = get_kmer_enrichments(mirna, experiment, n_constant, kmer_length, n_iter)
        R_values_sorted = R_values.sort_values(by="R",ascending=False)
        skaR_values_sorted = R_values.sort_values(by="skaR",ascending=False)

        print(R_values_sorted[0:10])
        print(skaR_values_sorted[0:10])

        print(R_value_path)
        R_values_sorted.to_csv(R_value_path,sep="\t")
        print_time_elapsed(time_start)

################################################################################

if __name__ == "__main__":
    main()
