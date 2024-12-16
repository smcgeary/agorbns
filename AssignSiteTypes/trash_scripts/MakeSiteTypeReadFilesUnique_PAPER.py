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

# mirna: The miRNA in the experiment to be processed.
# experiment: The type of experiment to be processed.
# condition: The sampel type within the experiment to be processed.

# Functions from general this uses:


# Constants

# FUNCTIONS

def assign_site_type_to_read(read_seqs, site_seq_map, start, stop, read_length):
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

    # Initialize count dictionary with keys as site names and a value of
    # for each item.
    count_site_map = {value: {"TGT" : 0, "ACA" : 0} for value in site_seq_map.values()}
    count_multisite_map = {value: {"TGT" : 0, "ACA" : 0} for value in site_seq_map.values()}
    count_site_map.update({"None": {"TGT" : 0, "ACA" : 0}})
    count_multisite_map.update({"None": {"TGT" : 0, "ACA" : 0}})
    # Pre-allocate output with "None" for each line.
    output = ["None"]*len(read_seqs)
    output_single = ["None"]*len(read_seqs)

    for i, read in enumerate(read_seqs):
        # Define the barcode portion of the read.
        [count, read] = read.strip().split(" ")
        barcode = read[26 + read_length : 26 + read_length + 3]
        if barcode == "TCG":
            barcode = "TGT"
        read_seq = read.strip()[(26 - start) : (26 + read_length + stop)]
        site_maps = get_unique_sites(read_seq, site_seq_map)
        if site_maps:
            site_maps = [j for j in site_maps]
        if site_maps:
            output_temp = ", ".join(["%s:%s-%s" % (site, index, end)
                                   for site, index, end in site_maps])
            output[i] = output_temp
            for site_map in site_maps:
                count_site_map[site_map[0]][barcode] += 1
            if len(site_maps) > 1:
                key = ",".join(sorted([j[0] for j in site_maps]))
                if key in count_multisite_map:
                    count_multisite_map[key][barcode] += 1
                else:
                    count_multisite_map[key] = {"TGT" : 0, "ACA" : 0}
                    count_multisite_map[key][barcode] += 1
            else:
                count_multisite_map[site_maps[0][0]][barcode] += 1
        else:
            count_site_map["None"][barcode] += 1
            count_multisite_map["None"][barcode] += 1
    return output, count_site_map, count_multisite_map

def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["miRNA", "experiment", "condition", "sitelist","start", "stop"]
    mirna, experiment, condition, sitelist, start, stop = parse_arguments(arguments)

    # Load file with site types for the miRNA.
    sites_file_name = ("/lab/bartel1_ata/mcgeary/computation/"
                       "AgoRBNS/AssignSiteTypes/sites.%s_%s.txt" % (mirna.split("-alt")[0], sitelist))

    print(sites_file_name)
    with open(sites_file_name) as file_in:
        sites = file_in.read().split("\n")

    # Initialize the dictionary with all sites, and two more dictionaries,
    # one mapping to the site name from the sequence, and one from site name
    # to the count.
    seq_site_map = get_seq_site_map(mirna, sitelist)
    site_seq_map = {value: key for key, value in seq_site_map.items()}
    for key in site_seq_map.keys():
        if "Centered" in site_seq_map[key]:
            site_seq_map[key] = "Centered"
    # Get the path to the read file and to that of where the site labels will
    # be written.
    if experiment in ["equil_mmseed_nb", "equil_seed_nb"]:
        read_length = 38
    else:
        read_length = 37    

    extension = "_%s-%s_%s" %(start, stop, sitelist)

    reads_path = get_analysis_path(mirna,experiment,condition,"full_reads_unique")
    sites_path = get_analysis_path(mirna,experiment,condition,"full_sites_unique",ext=extension)
    site_counts_path = get_analysis_path(mirna,experiment,condition,"full_site_counts_unique",ext=extension)
    multisite_counts_path = get_analysis_path(mirna,experiment,condition,"full_multisite_counts_unique",ext=extension)

    with open(reads_path,"rb") as file_in:
        results = multiprocess_file(file_in,
                                    readline,
                                    int(1e6),
                                    assign_site_type_to_read,
                                    site_seq_map,
                                    int(start),
                                    int(stop),
                                    read_length)

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
    bcs = list(set([k for i in dict_threads for j in i.values() for k in j.keys()]))[::-1]
    if experiment == "kinetics":
        values = [tuple(sum([thread[key][bc] for thread in dict_threads if key in thread])
                    for bc in bcs)
                  for key in keys]
    else:
        values = [sum((sum([thread[key][bc] for thread in dict_threads if key in thread])
                    for bc in bcs))
                  for key in keys]
    counts_sites_map = {key : value for (key,value) in zip(keys,values)}
    # No longer have input our output files open. Now write summmary file
    # to the same directory , ending in "summary.txt".
    with open(site_counts_path,"wb") as file_out:
        if experiment == "kinetics":
            for key in sorted(keys):
                value = counts_sites_map[key]
                file_out.write("%s:\t%s\t%s\n" % (key, value[0], value[1]))
        else:
            for key in sorted(keys):
                value = counts_sites_map[key]
                file_out.write("%s:\t%s\n" % (key, value))

    # Construct the dictionary from the list of futures, and write it to
    # its output file.
    keys = set([j for i in multi_threads for j in i.keys()])
    if experiment == "kinetics":
        values = [tuple(sum([thread[key][bc] for thread in multi_threads if key in thread])
                    for bc in bcs)
                  for key in keys]
    else:
        values = [sum((sum([thread[key][bc] for thread in multi_threads if key in thread])
                    for bc in bcs))
                  for key in keys]

    counts_sites_map = {key : value for (key,value) in zip(keys,values)}

    # No longer have input our output files open. Now write summmary file
    # to the same directory , ending in "summary.txt".
    with open(multisite_counts_path,"wb") as file_out:
        if experiment == "kinetics":
            for key in sorted(keys):
                value = counts_sites_map[key]
                file_out.write("%s:\t%s\t%s\n" % (key, value[0], value[1]))
        else:
            for key in sorted(keys):
                value = counts_sites_map[key]
                file_out.write("%s:\t%s\n" % (key, value))


    # Print the amount of time the script took to complete.
    print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()

