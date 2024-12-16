################################################################################
#MakeSiteTypeReadFiles.py
################################################################################
import imp # Used to import general.py
import time # Used for time
import numpy as np
import pandas as pd

imp.load_source("general",
                ("/lab/bartel1_ata/mcgeary/computation/"
                 "AgoRBNS/general/general.py"
                )
               )
from general import get_rc, seq_mirna_map, parse_arguments, print_time_elapsed, get_analysis_path, multiprocess_file, readline

# This script reads one read file and outputs multiple .txt files to be used in all
# downstream applications. Inputs are:

# mirna: The miRNA in the experiment to be processed.
# experiment: The type of experiment to be processed.
# condition: The sampel type within the experiment to be processed.

# Functions from general this uses:


# Constants

# FUNCTIONS

def get_site_positions(site_seqs, site_order, startpos, stoppos):
    """Takes a read sequence and identifies the position of all site types'.

    Args:
        read_sequence: The read sequence.
        sites: A list of sequences corresponding to each site type.

    Returns:
        A tuple containing the list of reads and a dictionary of read types.
    """

    # Initialize count dictionary with keys as site names and a value of
    # for each item.
    count_site5p_map = {key: np.array([0]*(37+startpos+stoppos)) for key in site_order.keys()}
    count_site3p_map = {key: np.array([0]*(37+startpos+stoppos)) for key in site_order.keys()}


    for i, site_seq in enumerate(site_seqs):
        # Remove trailing \n
        site_seq = site_seq.strip()
        if site_seq != "None":
            site_seq = site_seq.strip().split(", ")
            if len(site_seq) > 1:
                ranks = [site_order[i.split(":")[0]] for i in site_seq]
                site_seq = [site_seq[i] for i in range(len(site_seq)) if ranks[i] == min(ranks)]
            site, ends = site_seq[0].split(":")
            start, stop = [int(i) for i in ends.split("-")]
            count_site5p_map[site][start] += 1
            count_site3p_map[site][stop] += 1
    return count_site5p_map, count_site3p_map


def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["miRNA", "experient", "condition", "-start", "-stop"]
    mirna, experiment, condition, start, stop = parse_arguments(arguments)
    if start == None:
        start = 0
    else:
        start = int(start)
    if stop == None:
        stop = 0
    else:
        stop = int(stop)
    # Load file with site types for the miRNA.
    sites_file_name = ("/lab/bartel1_ata/mcgeary/computation/"
                       "AgoRBNS/AssignSiteTypes/sites.%s.txt" % (mirna))
    with open(sites_file_name) as file_in:
        sites = file_in.read().split("\n")
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


    # Get the path to the read file and to that of where the site labels will
    # be written.
    if start != 0 or stop !=0:
        sites_path = get_analysis_path(mirna,experiment,condition,"full_sites",ext="_%s-%s" %(start, stop))
        sites_position5p_counts_path = get_analysis_path(mirna,experiment,condition,"site_position5p_counts",ext="_%s-%s" %(start, stop))
        sites_position3p_counts_path = get_analysis_path(mirna,experiment,condition,"site_position3p_counts",ext="_%s-%s" %(start, stop))
    else:
        sites_path = get_analysis_path(mirna,experiment,condition,"sites")
        sites_position5p_counts_path = get_analysis_path(mirna,experiment,condition,"site_position5p_counts")
        sites_position3p_counts_path = get_analysis_path(mirna,experiment,condition,"site_position3p_counts")
    print(sites_path)
    print(sites_position5p_counts_path)
    print(sites_position3p_counts_path)
    with open(sites_path,"rb") as file_in:
        results = multiprocess_file(
            file_in, readline, 1000000, get_site_positions, order_site_map, start, stop)
    # Collect the output from each thread, which are the list
    # of reads and the dictionary of read types.
    count_pos5p_threads = [i[0] for i in results]
    count_pos3p_threads = [i[1] for i in results]
    counts_pos5p = pd.DataFrame(
        [np.array([thread[key] for thread in count_pos5p_threads ]).sum(axis=0)
         for key in sites],index=sites, columns = range(1,37+start+stop+1))
    counts_pos3p = pd.DataFrame(
        [np.array([thread[key] for thread in count_pos3p_threads ]).sum(axis=0)
         for key in sites],index=sites, columns = range(1,37+start+stop+1))

    counts_pos5p.to_csv(sites_position5p_counts_path, sep = "\t")
    counts_pos3p.to_csv(sites_position3p_counts_path, sep = "\t")

    print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()

