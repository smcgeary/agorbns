################################################################################
#MakeSiteTypeReadFiles.py
################################################################################
import imp # Used to import general.py
import time # Used for time
import numpy as np
import pandas as pd
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
from sitetypes import get_seq_site_map


def main():
    time_start = time.time()
    arguments = ["path"]
    args = parse_arguments(arguments)
    [path] = args
    with open(path, "r+") as file_in:
    	header = file_in.readline()
    	total = 0
    	line = file_in.readline()
    	while line:
	    	counts = int(line.strip().split("\t")[1]) + 1
	    	total += counts
	    	line = file_in.readline()
	print(total)


################################################################################

if __name__ == "__main__":
    main()

