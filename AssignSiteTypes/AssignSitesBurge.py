################################################################################
#MakeSiteTypeReadFiles.py
################################################################################
import imp # Used to import general.py
import time # Used for time
imp.load_source("general",
                ("/lab/bartel1_ata/mcgeary/computation/"
                 "AgoRBNS/general/general.py"
                )
               )
from general import parse_arguments, print_time_elapsed, get_analysis_path, multiprocess_file, multiprocess_test, readline
from sitetypes import get_seq_site_map
# This script reads one read file and outputs multiple .txt files to be used in all
# downstream applications. Inputs are:

# prot: The miRNA in the experiment to be processed.
# experiment: The type of experiment to be processed.
# condition: The sampel type within the experiment to be processed.

# Functions from general this uses:


# Constants

# FUNCTIONS

def assign_site_type_to_read(read_seqs, sites, start, stop):
    """Takes a read sequence and identifies the position of all site types'.

    Args:
        read_sequence: The read sequence.
        sites: A list of sequences corresponding to each site type.

    Returns:
        A tuple containing the list of reads and a dictionary of read types.
    """

    # Sub-function to be used on each site:
    def get_unique_sites(read_seq, sites):
        site_maps = []
        # Iterate over keys in dictionary (site seqs)
        for key in sites:
            # Assign left-hand value of site location, if it exists.
            i_l = read_seq.find(key)
            # Check if index is -1
            if i_l != -1:
                i_r = i_l + len(key) - 1
                site = (key, i_l, i_r)
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
                        site_maps_temp.append((key, i_l, i_r))
                    site_maps = site_maps_temp
        return site_maps

    # Initialize count dictionary with keys as site names and a value of
    # for each item.
    count_site_map = {value: 0 for value in sites}
    count_multisite_map = {value: 0 for value in sites}
    count_site_map.update({"None": 0})
    count_multisite_map.update({"None": 0})
    # Pre-allocate output with "None" for each line.
    output = ["None"]*len(read_seqs)

    for i, read in enumerate(read_seqs):
        # Define the barcode portion of the read.
        # print(read.strip())
        read_seq = read.strip()[(23 - start) : (23 + 40 + stop)]
        # print("_"*(23-start)+read_seq)
        site_maps = get_unique_sites(read_seq, sites)
        if site_maps:
            # print(site_maps)
            site_maps = [j for j in site_maps]
        if site_maps:
            output_temp = ", ".join(["%s:%s-%s" % (site, index, end)
                                   for site, index, end in site_maps])
            output[i] = output_temp
            for site_map in site_maps:
                # print("."*(23-start+int(site_map[1]))+site_map[0])
                count_site_map[site_map[0]] += 1
            if len(site_maps) > 1:
                key = ",".join(sorted([j[0] for j in site_maps]))
                if key in count_multisite_map:
                    count_multisite_map[key] += 1
                else:
                    count_multisite_map[key] = 0
                    count_multisite_map[key] += 1
            else:
                count_multisite_map[site_maps[0][0]] += 1
        else:
            count_site_map["None"] += 1
            count_multisite_map["None"] += 1
    return output, count_site_map, count_multisite_map

def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["protein", "condition", "sitelist","start", "stop"]
    prot, condition, sitelist, start, stop = parse_arguments(arguments)

    # Load file with site types for the prot.
    sites_file_name = ("/lab/bartel1_ata/mcgeary/computation/"
                       "AgoRBNS/AssignSiteTypes/sites.%s_%s.txt" % (prot.split("-alt")[0], sitelist))

    print(sites_file_name)
    with open(sites_file_name) as file_in:
        sites = file_in.read().split("\n")

    # Initialize the dictionary with all sites, and two more dictionaries,
    # one mapping to the site name from the sequence, and one from site name
    # to the count.
    print(sites)
    # Get the path to the read file and to that of where the site labels will
    # be written.
    extension = "_%s-%s_%s" %(start, stop, sitelist)

    reads_path = get_analysis_path(prot,"equilibrium",condition,"full_reads")
    sites_path = get_analysis_path(prot,"equilibrium",condition,"full_sites",ext=extension)
    site_counts_path = get_analysis_path(prot,"equilibrium",condition,"full_site_counts",ext=extension)
    multisite_counts_path = get_analysis_path(prot,"equilibrium",condition,"full_multisite_counts",ext=extension)

    with open(reads_path,"rb") as file_in:
        results = multiprocess_file(file_in,
                                    readline,
                                    int(1e6),
                                    assign_site_type_to_read,
                                    sites,
                                    int(start),
                                    int(stop))

    # Collect the output from each thread, which are the list
    # of reads and the dictionary of read types.
    site_threads = [i[0] for i in results]
    dict_threads = [i[1] for i in results]
    multi_threads = [i[2] for i in results]


    # Flatten the list into one list (found on stackoverflow) and write
    # it to its output file.
    site_assignments = [i for sublist in site_threads for i in sublist]
    with open(sites_path,"wb") as file_out:
        file_out.write("".join(["%s\n" % (i) for i in site_assignments]))

    # # Construct the dictionary from the list of futures, and write it to
    # # its output file.
    keys = set([j for i in dict_threads for j in i.keys()])
    values = [sum([thread[key] for thread in dict_threads if key in thread])
              for key in keys]
    counts_sites_map = {key : value for (key,value) in zip(keys,values)}
    # No longer have input our output files open. Now write summmary file
    # to the same directory , ending in "summary.txt".
    with open(site_counts_path,"wb") as file_out:
        for key in sorted(keys):
            value = counts_sites_map[key]
            file_out.write("%s:\t%s\n" % (key, value))

    # Construct the dictionary from the list of futures, and write it to
    # its output file.
    keys = set([j for i in multi_threads for j in i.keys()])
    values = [sum([thread[key] for thread in multi_threads if key in thread])
              for key in keys]

    counts_sites_map = {key : value for (key,value) in zip(keys,values)}

    # No longer have input our output files open. Now write summmary file
    # to the same directory , ending in "summary.txt".
    with open(multisite_counts_path,"wb") as file_out:
        for key in sorted(keys):
            value = counts_sites_map[key]
            file_out.write("%s:\t%s\n" % (key, value))


    # Print the amount of time the script took to complete.
    print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()

