################################################################################
#MakeSiteTypeReadFiles.py
################################################################################
import imp # Used to import general.py
import time # Used for time
import math
imp.load_source("general",
                ("/lab/bartel1_ata/mcgeary/computation/"
                 "AgoRBNS/general/general.py"
                )
               )
from general import parse_arguments, print_time_elapsed, get_analysis_path, multiprocess_file, multiprocess_test, readline
from sitetypes import get_seq_site_map
# This script reads one read file and outputs multiple .txt files to be used in all
# downstream applications. Inputs are:

# mirna: The miRNA in the experiment to be processed.
# experiment: The type of experiment to be processed.
# condition: The sampel type within the experiment to be processed.

# Functions from general this uses:


# Constants

# FUNCTIONS

# TODO: comment this.
def assign_site_type_to_read(read_seqs, site_seq_map, site, start, stop):
    """Takes a read sequence and identifies the position of all site types'.

    Args:
        read_sequence: The read sequence.
        sites: A list of sequences corresponding to each site type.

    Returns:
        A tuple containing the list of reads and a dictionary of read types.
    """

    # Sub-function to be used on each site:
    def get_unique_sites(read_seq, site_seq_map):
        site_maps = []
        # Iterate over keys in dictionary (site seqs)
        for key, value in site_seq_map.items():
            # Assign left-hand value of site location, if it exists.
            i_l = read_seq.find(key)

            # Check if index is -1
            if i_l != -1:
                i_r = i_l + len(key) - 1
                site = (value, i_l, i_r)
                # if "CCATTCCA" in read_seq:
                #     print(site_seq_map)
                # Check if site_maps is empty (suggesting first map)
                if site_maps == []:

                    site_maps.append(site)

                # Else only add if the sites in the list do not overlap.
                else:
                    add = True
                    site_maps_temp = [i for i in site_maps]
                    for site_map in site_maps:
                        j_l, j_r = site_map[1: ]

                        # Checks if new map is an extension of old map
                        if i_l <= j_l and i_r >= j_r:
                            
                            # Removes map
                            site_maps_temp.remove(site_map)

                        elif ((i_l >= j_l and i_r < j_r)
                              or (i_l > j_l and i_r <= j_r)):
                            add = False
                    if add:
                        site_maps_temp.append((value, i_l, i_r))
                    site_maps = site_maps_temp

        return site_maps



    output = []
    for i, read_seq in enumerate(read_seqs):
        read_seq = read_seq.strip()[(26-5):(26+37+5)]
        site_maps = get_unique_sites(read_seq, site_seq_map)
        if site_maps:
            site_maps = [j for j in site_maps if (j[1] - 5 >= -1*start) and (j[2] - 5 < 37 + stop)]
        if site_maps:
            if site in [j[0] for j in site_maps]:
                pos = [j[1] for j in site_maps if site in j[0]][0]
                output.append((read_seq,pos))

    return output

def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["miRNA", "experient", "condition", "site","-start", "-stop","-min","-max"]
    mirna, experiment, condition, site_interest, start, stop, min_space, max_space = parse_arguments(arguments)
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

    with open('general/site_info.txt',"r+") as site_info_file:
        site_region_map = dict()
        site_info_file.readline()
        for line in site_info_file:
            row = line.strip().split("\t")
            key = row[0]
            value = row[5]
            site_region_map[key] = value

    with open(sites_file_name) as file_in:
        sites = file_in.read().split("\n")

    # Initialize the dictionary with all sites, and two more dictionaries,
    # one mapping to the site name from the sequence, and one from site name
    # to the count.
    seq_site_map = get_seq_site_map(mirna)
    site_seq_map = {value: key for key, value in seq_site_map.items()}
    for key in site_seq_map.keys():
        if "Centered" in site_seq_map[key]:
            site_seq_map[key] = "Centered"
    # Get the path to the read file and to that of where the site labels will
    # be written.
    reads_path = get_analysis_path(mirna,experiment,condition,"full_reads")
    sites_path = get_analysis_path(mirna,experiment,condition,"full_double_sites",ext="_%s-%s" %(start, stop))
    output_list_path = get_analysis_path(mirna, experiment, condition, "aligned_reads_by_site",ext = "_%s_%s-%s" %(site_interest, start, stop))

    with open(reads_path,"rb") as file_in:
        results = multiprocess_file(
            file_in, readline, 1000000, assign_site_type_to_read, site_seq_map, site_interest, start, stop)
    # Collect the output from each thread, which are the list
    # of reads and the dictionary of read types.
    list = [j for i in results for j in i]
    list = [i for i in list if int(i[1]) > int(min_space) & int(i[1]) < int(max_space)]
    maximum = max([int(i[1]) for i in list])
    minimum = min([int(i[1]) for i in list])
    out = []
    for i in list:
        out_i = i[0][i[1]-minimum:-(maximum - i[1])]
        if len(out_i) > 0:
            out.append(out_i)
    # Flatten the list into one list (found on stackoverflow) and write
    # it to its output file.
    # site_assignments = [i for sublist in site_threads for i in sublist]
    with open(output_list_path,"wb") as file_out:
        file_out.write("".join(["%s\n" % (i) for i in out]))

    # # Construct the dictionary from the list of futures, and write it to
    # # its output file.
    # keys = set([j for i in dict_threads for j in i.keys()])
    # values = [sum([thread[key] for thread in dict_threads if key in thread])
    #           for key in keys]
    # counts_sites_map = {key : value for (key,value) in zip(keys,values)}

    # # No longer have input our output files open. Now write summmary file
    # # to the same directory , ending in "summary.txt".
    # with open(site_counts_path,"wb") as file_out:
    #     for key in sorted(keys):
    #         value = counts_sites_map[key]
    #         file_out.write("%s:\t%s\n" % (key, value))

    # # Construct the dictionary from the list of futures, and write it to
    # # its output file.
    # keys = set([j for i in multi_threads for j in i.keys()])
    # values = [sum([thread[key] for thread in multi_threads if key in thread])
    #           for key in keys]
    # counts_sites_map = {key : value for (key,value) in zip(keys,values)}

    # # No longer have input our output files open. Now write summmary file
    # # to the same directory , ending in "summary.txt".
    # with open(multisite_counts_path,"wb") as file_out:
    #     for key in sorted(keys):
    #         value = counts_sites_map[key]
    #         file_out.write("%s:\t%s\n" % (key, value))

    # # Construct the dictionary from the list of futures, and write it to
    # # its output file.
    # keys = set([j for i in bipartite_threads for j in i.keys()])
    # values = [[len(j) for thread in bipartite_threads if key in thread.keys() for j in thread[key]]
    #           for key in keys]
    # max_dist = max([i for j in values for i in j])

    # counts_sites_map = {key : "\t".join([str(value.count(i)) for i in range(max_dist)]) for (key,value) in zip(keys,values)}


    # # No longer have input our output files open. Now write summmary file
    # # to the same directory , ending in "summary.txt".
    # with open(bipartite_counts_path,"wb") as file_out:
    #     for key in sorted(keys):
    #         value = counts_sites_map[key]
    #         file_out.write("%s:\t%s\n" % (key, value))


    # Print the amount of time the script took to complete.
    print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()

