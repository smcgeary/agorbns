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

# TODO: comment this.
def assign_site_type_to_read(read_seqs, site_seq_map, site_region_map, start, stop):
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



    # Initialize count dictionary with keys as site names and a value of
    # for each item.
    count_site_map = {value: 0 for value in site_seq_map.values()}
    count_multisite_map = {value: 0 for value in site_seq_map.values()}
    count_bipartite_map = dict()
    count_site_map.update({"None": 0})
    count_multisite_map.update({"None": 0})
    # Pre-allocate output with "None" for each line.
    output = ["None"]*len(read_seqs)
    # memory = False
    for i, read_seq in enumerate(read_seqs):
        read_seq = read_seq.strip()[(26-5):(26+37+5)]
        site_maps = get_unique_sites(read_seq, site_seq_map)
        if site_maps:
            site_maps = [j for j in site_maps if (j[1] - 5 >= -1*start) and (j[2] - 5 < 37 + stop)]
        if site_maps:
            output_temp = ", ".join(["%s:%s-%s" % (site, index - 5 + start, end - 5 + start)
                                   for site, index, end in site_maps if (index - 5 >= -1*start) and (end - 5 < 37 + stop)])
            output[i] = output_temp
            dist = False
            regions = [site_region_map[site_map[0]] for site_map in site_maps]
            if "Seed" in regions and "3prime" in regions:
                seed_inds = [i for i, x in enumerate(regions) if x == "Seed"]
                threep_inds = [i for i, x in enumerate(regions) if x == "3prime"]
                seed_starts = [site_maps[i][1] for i in seed_inds]
                threep_ends = [site_maps[i][2] for i in threep_inds]
                ind_seed = False
                ind_threep = False
                for i, seed_start in enumerate(seed_starts):
                    for j, threep_end in enumerate(threep_ends):
                        if seed_start - threep_end > 0:
                            if dist:
                                if seed_start - threep_end < dist:
                                    dist = seed_start - threep_end
                                    ind_seed = i
                                    ind_threep = j
                            else:
                                dist = seed_start - threep_end
                                ind_seed = i
                                ind_threep = j
                if dist:
                    threep_final = site_maps[threep_inds[ind_threep]]
                    seed_final = site_maps[seed_inds[ind_seed]]
                    dist_seq = read_seq[threep_final[2]+1:seed_final[1]]
                    bipartite_key = ("%s|%s" %(seed_final[0], threep_final[0]))
                    if bipartite_key in count_bipartite_map.keys():
                        count_bipartite_map[bipartite_key].append(dist_seq)
                    else:
                        count_bipartite_map[bipartite_key] = [dist_seq]
                    site_maps.remove(seed_final)
                    site_maps.remove(threep_final)
            for site_map in site_maps:
                count_site_map[site_map[0]] += 1
            if dist == False:
                key = ",".join(sorted([j[0] for j in site_maps]))
                if key in count_multisite_map:
                    count_multisite_map[key] += 1
                else:
                    count_multisite_map[key] = 1
            elif site_maps:
                count_multisite_map[site_maps[0][0]] += 1
        else:
            count_site_map["None"] += 1
            count_multisite_map["None"] += 1

    return output, count_site_map, count_multisite_map, count_bipartite_map

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
    site_counts_path = get_analysis_path(mirna,experiment,condition,"full_double_site_counts",ext="_%s-%s" %(start, stop))
    multisite_counts_path = get_analysis_path(mirna,experiment,condition,"full_double_multisite_counts",ext="_%s-%s" %(start, stop))
    bipartite_counts_path = get_analysis_path(mirna,experiment,condition,"full_double_bipartite_counts",ext="_%s-%s" %(start, stop))

    with open(reads_path,"rb") as file_in:
        results = multiprocess_file(
            file_in, readline, 1000000, assign_site_type_to_read, site_seq_map, site_region_map, start, stop)
    # Collect the output from each thread, which are the list
    # of reads and the dictionary of read types.
    site_threads = [i[0] for i in results]
    dict_threads = [i[1] for i in results]
    multi_threads = [i[2] for i in results]
    bipartite_threads = [i[3] for i in results]

    # Flatten the list into one list (found on stackoverflow) and write
    # it to its output file.
    site_assignments = [i for sublist in site_threads for i in sublist]
    with open(sites_path,"wb") as file_out:
        file_out.write("".join(["%s\n" % (i) for i in site_assignments]))

    # Construct the dictionary from the list of futures, and write it to
    # its output file.
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

    # Construct the dictionary from the list of futures, and write it to
    # its output file.
    keys = set([j for i in bipartite_threads for j in i.keys()])
    values = [[len(j) for thread in bipartite_threads if key in thread.keys() for j in thread[key]]
              for key in keys]
    max_dist = max([i for j in values for i in j])

    counts_sites_map = {key : "\t".join([str(value.count(i)) for i in range(max_dist)]) for (key,value) in zip(keys,values)}


    # No longer have input our output files open. Now write summmary file
    # to the same directory , ending in "summary.txt".
    with open(bipartite_counts_path,"wb") as file_out:
        for key in sorted(keys):
            value = counts_sites_map[key]
            file_out.write("%s:\t%s\n" % (key, value))


    # Print the amount of time the script took to complete.
    print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()

