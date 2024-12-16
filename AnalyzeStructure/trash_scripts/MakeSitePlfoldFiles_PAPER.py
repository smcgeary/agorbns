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
import time
import itertools as it
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
from general import *
from RBNS_methods import *

#23456789!123456789@123456789#123456789$123456789%123456789^123456789&123456789*
# SCRIPT SPECIFIC VARIABLES
PL_WIN = 30
LABELS_5P, LABELS_3P = [["f%sp%s" %(i, j + 1) for j in range(20)]
                        for i in [5, 3]]
LABELS_SEED = ["mi%s" %(i + 1) for i in range(8)]
FULL_LABELS = (["Read"] + ["Flank"] + LABELS_5P[::-1] +
               LABELS_SEED[::-1] + LABELS_3P)
print(FULL_LABELS)
OUT_HEADER = "\t".join(FULL_LABELS) + "\n"
PRINT_BUFFER = 25

# FUNCTIONS
def get_read_structural_data(read_seqs, job_ind, _sitelist, site_mir5pnt_map,
                             n_constant, read_length, _mirna, experiment,
                             condition, site):
    # These are the slope and intercept for converting the logplfold to plfold.
    # (I had to figure this out manually using lm in R!)
    print("job index:\t%s" %(job_ind))
    b = -0.00001798580677720675
    m = -1.62250426156107407927
    # Initialize count dictionary with keys as site names and a value of
    # for each item.
    sites_range = range(26 - n_constant, 26 + 37 + n_constant + 1)
    # These are the position labels, no longer keys as I'm doing it.
    time_start = time.time()
    ##############PL
    pl_inds = [-20, 28]
    if "b" in site:
        pl_inds[1] += 1
   
    ##############
    # SET UP THE FILES BEING WRITTEN TO:
    exts = ["_%s_win%s_thread%s_%s" %(site, i_f, job_ind, randomword(8))
            for i_f in range(PL_WIN)]
    paths = [get_analysis_path(_mirna.name, experiment, condition,
                               "plfold_2018_temp", ext=i)
             for i in exts]
    outs = [open(i, "w") for i in paths]
    # This writes the header, closes it, and then reopens the files:
    for out in outs:
        # if job_ind == 0:
        #     out.write("\t".join(OUT_HEADER) + "\n")
        out.write("")
        out.close()
    outs = [open(i, "a") for i in paths]
    # Iterate through the reads:
    for r in read_seqs:
        # Get the read number:
        [i_r, read] = (r[0], r[1].strip())
        # NEW Get the _read object, work with this.
        _read = Read(read)
        _read.get_all_sites_in_read(_sitelist)
        # Find all sites within the read
        top_site = _read.site_name()
        top_pos = _read.site_pos()
        if top_site == site and (p in sites_range for p in top_pos):
            site_l, site_r = top_pos
            # Get the left-most paired position of the miRNA at the site:
            i_mirnt1 = site_l + site_mir5pnt_map[site] - 1
            flank = _read.site_flank()
            plfold = RNAplfold(read, PL_WIN, i_r)
            for i_win in range(PL_WIN):
                win_r, win_l = (i_mirnt1 - i + 1 + i_win/2 for i in pl_inds)
                pl_l = max(0, win_l)
                pl_r = min(win_r, len(read))
                Na_l = [float('nan')]*max(0, -win_l)
                Na_r = [float('nan')]*max(0, win_r - len(read))                 ###
                # Perform plfold on the entire read:
                pl_output = Na_l + list(plfold.iloc[pl_l:pl_r, i_win]) + Na_r
                output = [i_r] + [flank] + pl_output
                outs[i_win].write("\t".join([str(i) for i in output]) + "\n")
                struc = RNA.fold(read)
                # print(" "*PRINT_BUFFER + read + "\t" + site)
                # print(" "*(PRINT_BUFFER + site_l) + read[site_l:site_r])
                # for num in range(0,5):
                #     print(" "*(PRINT_BUFFER+win_l-i_win/2*2)+"_"*(i_win/2)+"".join([(str(k) + " "*(10 - len(str(k))))[num] for k in output]))
                # print(" "*PRINT_BUFFER + struc[0])
                # for num in range(1,4):
                #     print(" "*(PRINT_BUFFER+win_l-i_win/2)+"".join([(k + " "*(10-len(k)))[num] for k in file_header[2:]]))
    return paths

def main():
    time_start = time.time()
    # Define all the relevant arguments.
    arguments = ["miRNA", "experiment", "condition", "n_constant", "sitelist",
                 "site"]
    args = parse_arguments(arguments)   
    mirna, experiment, condition, n_constant, sitelist, site = args
    # Load the _mirna and _sitelist classes:
    _mirna = Mirna(mirna)
    _sitelist = SiteList(_mirna, sitelist)
    # Make a stupid dictionary:
    with open('general/site_info.txt',"r+") as site_pos_file:   ##
        site_mir5pnt_map = dict()                                                   ##
        for line in site_pos_file:                                              ##
            line_strip = line.strip().split("\t")                               ##
            key = line_strip[0]                                                 ##
            if key in _sitelist.list_site_names() and key != "None":            ##
                value = int(line_strip[2]) + int(line_strip[3]) - 1             ##
                site_mir5pnt_map[key] = value                                       ##
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
    # sites_path = get_analysis_path(mirna, experiment, condition,
    #                                "sites", ext=extension_sites)
    with open(reads_path, "rb") as file_in_reads:
        # with open(sites_path, "rb") as file_in_sites:
        results = multiprocess_iter(file_in_reads,
                                         readline_one,
                                         int(5e4),
                                         get_read_structural_data,
                                         _sitelist,
                                         site_mir5pnt_map,
                                         int(n_constant),
                                         read_length,
                                         _mirna,
                                         experiment,
                                         condition,
                                         site)
    # Iterate over the range of windows (each is a different out file):
    base_extension = "_%s_%s_%s" %(n_constant, sitelist, site)
    for i_f in range(PL_WIN):
        extension = "%s_win%s" %(base_extension, i_f)
        # Construct path and open file:
        path = get_analysis_path(_mirna.name, experiment, condition,
                "plfold_2018_PAPER", ext=extension)
        with open(path, "w") as out_file:
            # Write the header:
            out_file.write(OUT_HEADER)
            # Iterate over the results:
            for result in results:
                # Write the entire tempfile to output:
                with open(result[i_f], "rb") as in_file:
                    out_file.write(in_file.read())
                # Remove the old file.
                os.remove(result[i_f])

    print_time_elapsed(time_start)

################################################################################

if __name__ == "__main__":
    main()

