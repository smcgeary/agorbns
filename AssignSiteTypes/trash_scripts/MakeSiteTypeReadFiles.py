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
from general import parse_arguments, print_time_elapsed, get_analysis_path, multiprocess_file, readline
from sitetypes import get_seq_site_map, assign_site_type_to_read
# This script reads one read file and outputs multiple .txt files to be used in all
# downstream applications. Inputs are:

# mirna: The miRNA in the experiment to be processed.
# experiment: The type of experiment to be processed.
# condition: The sampel type within the experiment to be processed.



def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["miRNA","experient","condition", "-nb_binary"]
    mirna, experiment, condition, nb = parse_arguments(arguments)

    # Load file with site types for the miRNA.
    sites_file_name = ("/lab/bartel1_ata/mcgeary/computation/"
                       "AgoRBNS/AssignSiteTypes/sites.%s.txt" % (mirna))
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
    if nb:
        experiment = experiment + "_nb"
    reads_path = get_analysis_path(mirna,experiment,condition,"reads")
    sites_path = get_analysis_path(mirna,experiment,condition,"sites")
    site_counts_path = get_analysis_path(mirna,experiment,condition,"site_counts")
    multisite_counts_path = get_analysis_path(mirna,experiment,condition,"multisite_counts")

    with open(reads_path,"rb") as file_in:
        results = multiprocess_file(
            file_in, readline, int(1e6), assign_site_type_to_read, site_seq_map)

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


    # Print the amount of time the script took to complete.
    print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()

