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
def assign_site_type_to_read(read_seqs, site_seq_map, start, stop, mirna, experiment):
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
    # memory = False
    for i, read in enumerate(read_seqs):
        read = read.strip()
        barcode = read[26 + 37 : 26 + 37 + 3]
        if barcode == "TCG":
            barcode = "TGT"
            if mirna != "miR-1" or experiment != "equilibrium":
                read_seq = read_seq[:26 + 37]+"TGTTCGTATGCCGTCTTCTGCTTG"
        if barcode == "TGT" and mirna == "miR-1" and experiment == "equilibrium":
            read_seq = read_seq[:26 + 37]+"TCGTATGCCGTCTTCTGCTTG"

        read_seq = read.strip()[(26 - start) : (26 + 37 + stop)]

        site_maps = get_unique_sites(read_seq, site_seq_map)
        if site_maps:
            site_maps = [j for j in site_maps if (j[1] - 5 >= -1*start) and (j[2] - 5 < 37 + stop)]
        if site_maps:
            output_temp = ", ".join(["%s:%s-%s" % (site, index, end)
                                   for site, index, end in site_maps])
            output[i] = output_temp
            print(site_maps)
            print(read.strip())
            print(read.strip()[:26] + "_"*37 + read.strip()[(26 + 37):])
            print(" "*(26 - start) + read_seq)
            for site_map in site_maps:
                print(site_map[0])
                print("."*(26 - start + site_map[1]) + read_seq[site_map[1]:site_map[2]+1])

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
    arguments = ["miRNA", "experient", "sitelist","start", "stop", "-nb_binary", "-chase_barcode"]
    print(parse_arguments(arguments))
    mirna, experiment, sitelist, start, stop, nb, chase_barcode = parse_arguments(arguments)
    if chase_barcode == None:
        chase_barcoe = False
    else:
        chase_barcode = True

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
    if nb:
        experiment = experiment + "_nb"
    extension = "_%s-%s_%s" %(start, stop, sitelist)
    sites_path = get_analysis_path(mirna, experiment, "I_combined", "full_sites", ext=extension)
    site_counts_path = get_analysis_path(mirna, experiment, "I_combined", "full_site_counts", ext=extension)
    multisite_counts_path = get_analysis_path(mirna, experiment, "I_combined", "full_multisite_counts", ext=extension)
    reads_output_path = get_analysis_path(mirna, experiment, "I_combined", "full_reads", ext=extension)
    with open(reads_output_path, "w") as file_out_write:
        file_out_write.write("")

    pairs = [("let-7a", "equilibrium"),
             ("miR-1", "equilibrium"),
             ("miR-124", "equilibrium"),
             ("let-7a", "kinetics"),
             ("miR-1", "kinetics"),
             ("miR-124", "kinetics"),
             ("lsy-6", "kinetics"),
             ("miR-1","kinetics_pilot")]
    
    site_assignments = []
    if experiment == "kinetics":
        counts_sites_map = {key : [0, 0] for key in sites+["None"]}
    else:
        counts_sites_map = {key : 0 for key in sites+["None"]}
    counts_multisites_map = {}


    for pair in pairs:
        mirna_input, experiment_input = pair
        if experiment_input == "kinetics_pilot":
            condition_temp = "I_TGT"
        else:
            condition_temp = "I"
        reads_path = get_analysis_path(mirna_input,experiment_input,condition_temp,"full_reads")

        print(pair)
        with open(reads_path,"rb") as file_in:
            results = multiprocess_test(file_in,
                                        readline,
                                        int(1e6),
                                        assign_site_type_to_read,
                                        site_seq_map,
                                        100,
                                        int(start),
                                        int(stop),
                                        mirna,
                                        experiment)
        # with open(reads_path,"rb") as file_in:
        #     with open(reads_output_path, "a") as file_out_write:
        #         file_out_write.write(file_in.read())



    #     # Collect the output from each thread, which are the list
    #     # of reads and the dictionary of read types.
    #     site_threads = [i[0] for i in results]
    #     dict_threads = [i[1] for i in results]
    #     multi_threads = [i[2] for i in results]

    # # Flatten the list into one list (found on stackoverflow) and write
    # # it to its output file.
    #     site_assignments += [i for sublist in site_threads for i in sublist]
    #     keys = set([j for i in dict_threads for j in i.keys()])
    #     bcs = list(set([k for i in dict_threads for j in i.values() for k in j.keys()]))[::-1]
    #     print(bcs)
    #     if experiment == "kinetics":
    #         values = [tuple(sum([thread[key][bc] for thread in dict_threads if key in thread])
    #                     for bc in bcs)
    #                   for key in keys]
    #         for (key,value) in zip(keys, values):
    #             counts_sites_map[key][0] += value[0]
    #             counts_sites_map[key][1] += value[1]                

    #     else:
    #         check = [tuple(sum([thread[key][bc] for thread in dict_threads if key in thread])
    #                     for bc in bcs)
    #                  for key in keys]
    #         values = [tuple(sum([thread[key][bc] for thread in dict_threads if key in thread])
    #                     for bc in bcs)[0]
    #                   for key in keys]
    #         for (key,value) in zip(keys, values):
    #             counts_sites_map[key] += value

    #     keys = set([j for i in multi_threads for j in i.keys()])
    #     if experiment == "kinetics":
    #         values = [tuple(sum([thread[key][bc] for thread in multi_threads if key in thread])
    #                     for bc in bcs)
    #                   for key in keys]
    #         for (key,value) in zip(keys, values):
    #             if key in counts_multisites_map.keys():
    #                 counts_multisites_map[key][0] += value[0]
    #                 counts_multisites_map[key][1] += value[1]                
    #             else:
    #                 counts_multisites_map[key] = list(value)

    #     else:
    #         values = [tuple(sum([thread[key][bc] for thread in multi_threads if key in thread])
    #                     for bc in bcs)[0]
    #                   for key in keys]

    #         for (key,value) in zip(keys, values):
    #             if key in counts_multisites_map:
    #                 counts_multisites_map[key] += value
    #             else:
    #                 counts_multisites_map[key] = value


    # with open(sites_path,"wb") as file_out:
    #     file_out.write("".join(["%s\n" % (i) for i in site_assignments]))

    # # # Construct the dictionary from the list of futures, and write it to
    # # # its output file.
    # # No longer have input our output files open. Now write summmary file
    # # to the same directory , ending in "summary.txt".
    # with open(site_counts_path,"wb") as file_out:
    #     if experiment == "kinetics":
    #         for key in sorted(counts_sites_map.keys()):
    #             value = counts_sites_map[key]
    #             file_out.write("%s:\t%s\t%s\n" % (key, value[0], value[1]))
    #     else:
    #         for key in sorted(counts_sites_map.keys()):
    #             value = counts_sites_map[key]
    #             file_out.write("%s:\t%s\n" % (key, value))

    # # Construct the dictionary from the list of futures, and write it to
    # # its output file.

    # # No longer have input our output files open. Now write summmary file
    # # to the same directory , ending in "summary.txt".
    # with open(multisite_counts_path,"wb") as file_out:
    #     if experiment == "kinetics":
    #         for key in sorted(counts_multisites_map.keys()):
    #             value = counts_multisites_map[key]
    #             file_out.write("%s:\t%s\t%s\n" % (key, value[0], value[1]))
    #     else:
    #         for key in sorted(counts_multisites_map.keys()):
    #             value = counts_multisites_map[key]
    #             file_out.write("%s:\t%s\n" % (key, value))


    # # Print the amount of time the script took to complete.
    # print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()

