################################################################################
#MakeSiteTypeReadFiles.py
################################################################################
import imp # Used to import general.py
import sys
import time # Used for time
import math
imp.load_source("general",
                ("/lab/bartel1_ata/mcgeary/computation/"
                 "AgoRBNS/general/general.py"
                )
               )
from general import get_kmer_list
from general import parse_arguments
from general import print_time_elapsed
from general import get_analysis_path
from general import multiprocess_file
from general import multiprocess_test
from general import readline
from general import get_rc
from general import get_kmer_list_constant_insert
from general import seq_mirna_map
from sitetypes import get_site_seq
from sitetypes import get_seq_site_map

# FUNCTIONS

def assign_site_type_to_read(read_seqs, kmer_list, start, stop, kmer_length, mir_start, mir_stop):
    """Takes a read sequence and identifies the position of all site types'.

    Args:
        read_sequence: The read sequence.
        sites: A list of sequences corresponding to each site type.

    Returns:
        A tuple containing the list of reads and a dictionary of read types.
    """
    constant_len = mir_stop - mir_start + 1
    offset_index = int(math.ceil((kmer_length - 8)/2.0))
    start_ind = 8 - mir_stop + offset_index
    stop_ind = 8 - mir_start + offset_index + 1
    constant_seq = kmer_list[0][start_ind : stop_ind]
    print(constant_seq)
    # Sub-function to be used on each site:
    def get_unique_sites(read_seq, kmer_list):
        site_maps = []
        # Iterate over keys in dictionary (site seqs)
        # print(read_seq)
        for index in range(0, len(read_seq) - kmer_length + 1):
            kmer = read_seq[index : index + kmer_length]
            seed = kmer[start_ind : stop_ind]
            # Check if index is -1
            if seed==constant_seq:
                site_maps.append((kmer, index))
        return site_maps

    # Initialize count dictionary with keys as site names and a value of
    # for each item.
    count_site_map = {kmer: {"TGT" : 0, "ACA" : 0} for kmer in kmer_list}
    count_multisite_map = {kmer: {"TGT" : 0, "ACA" : 0} for kmer in kmer_list}
    count_site_map.update({"None": {"TGT" : 0, "ACA" : 0}})
    count_multisite_map.update({"None": {"TGT" : 0, "ACA" : 0}})
    # Pre-allocate output with "None" for each line.
    output = ["None"]*len(read_seqs)
    time_old = time.clock()
    for i, read in enumerate(read_seqs):
        if i % 100000 == 0:
            print(i)
            time_new = time.clock()
            print(time_new - time_old)
            time_old = time_new

            sys.stdout.flush()
        # Define the barcode portion of the read.
        barcode = read[26 + 37 : 26 + 37 + 3]
        if barcode == "TCG":
            barcode = "TGT"
        read_seq = read.strip()[(26 - start) : (26 + 37 + stop)]
        site_maps = get_unique_sites(read_seq, kmer_list)
        if site_maps:
            site_maps = [j for j in site_maps]
            # print(read.strip())
            # for site in site_maps:
            #     print(" "*(26 - start + site[1]) + site[0])

        if site_maps:
            output_temp = ", ".join(["%s:%s" % (site, index)
                                   for site, index in site_maps])
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
    arguments = ["miRNA", "experient", "condition", "start", "stop",
                 "kmer","miR_start", "miR_stop", "-nb_binary"]
    (mirna, experiment, condition, start, stop, kmer_length, mir_start, 
     mir_stop,nb) = parse_arguments(arguments)

    kmer_list = get_kmer_list_constant_insert(mirna, int(kmer_length), int(mir_start),
                                              int(mir_stop))
    # Initialize the dictionary with all sites, and two more dictionaries,
    # one mapping to the site name from the sequence, and one from site name
    # to the count.
    # Get the path to the read file and to that of where the site labels will
    # be written.
    if nb:
        experiment = experiment + "_nb"
    extension = "_%s-%s_%smers_%s-%s" %(start, stop, kmer_length, mir_start, mir_stop)

    reads_path = get_analysis_path(mirna,experiment,condition,"full_reads")
    sites_path = get_analysis_path(mirna,experiment,condition,"full_sites",ext=extension)
    site_counts_path = get_analysis_path(mirna,experiment,condition,"full_site_counts",ext=extension)
    multisite_counts_path = get_analysis_path(mirna,experiment,condition,"full_multisite_counts",ext=extension)
    print(sites_path)
    print(site_counts_path)
    print(multisite_counts_path)
    with open(reads_path,"rb") as file_in:
        results = multiprocess_file(file_in,
                                    readline,
                                    int(1e6),
                                    assign_site_type_to_read,
                                    kmer_list,
                                    int(start),
                                    int(stop),
                                    int(kmer_length),
                                    int(mir_start),
                                    int(mir_stop))

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


###############################################################################

if __name__ == "__main__":
    main()

