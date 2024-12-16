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
import itertools as it
from general import get_rc, seq_mirna_map, parse_arguments, print_time_elapsed, get_analysis_path, multiprocess_file, readline, readline_two

# This script reads one read file and outputs multiple .txt files to be used in all
# downstream applications. Inputs are:

# mirna: The miRNA in the experiment to be processed.
# experiment: The type of experiment to be processed.
# condition: The sampel type within the experiment to be processed.

# Functions from general this uses:


# Constants

# FUNCTIONS



def get_flanking_nucleotide(read_seqs, order_site_map, startpos, stoppos):
    """Takes a read sequence and identifies the position of all site types'.

    Args:
        read_sequence: The read sequence.
        sites: A list of sequences corresponding to each site type.

    Returns:
        A tuple containing the list of reads and a dictionary of read types.
    """

    # Sub-function to be used on each site:
    count_flank_site_map = {
        site: {"".join(kmer[:2]+(".", )+kmer[2:]): 0
               for kmer in it.product(["A","C","G","T"],repeat=4)}
        for site in order_site_map.keys()}
    for i in read_seqs:
        seq, sites = (j.strip() for j in i)
        seq = seq[(26-startpos):(26+37+stoppos)]
        # print(seq)
        if sites != "None":
            coords = [i.split(":")[1] for i in sites.split(", ")]
            sites = [i.split(":")[0] for i in sites.split(", ")]
            coord_site_map = {site: coord for site, coord in zip(sites,coords)}
            ranks = [order_site_map[site] for site in sites]
            [(coord, site)] = [i for i in zip(coords, sites)
                               if order_site_map[i[1]] == min(ranks)]
            start = int(coord.split("-")[0])
            stop = int(coord.split("-")[1]) + 1
            # print(site)
            # print("_"*start+seq[start:stop])
            left_start = start - 2
            right_stop = stop + 2
            if left_start >= 0 and right_stop <= len(seq):
                # print("_"*left_start+seq[left_start:right_stop])
                flank = seq[left_start:start]+"."+ seq[stop:right_stop]
                count_flank_site_map[site][flank] += 1

    return count_flank_site_map

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
    # This swaps the four "Centered-..." for one single "Centered" entry in the list.
    add = True
    sites_new = []
    for i in sites:
        if "Centered" not in i:
            sites_new.append(i)
        elif add:
            sites_new.append("Centered")
            add = False
    sites = sites_new  
    # print(sites) 
    order_site_map = {site: i for i, site in enumerate(sites)}
    print(order_site_map)
    # Initialize the dictionary with all sites, and two more dictionaries,
    # one mapping to the site name from the sequence, and one from site name
    # to the count.

    # Get the path to the read file and to that of where the site labels will
    # be written.
    if start != 0 or stop !=0:
        reads_path = get_analysis_path(mirna, experiment, condition, "full_reads")
        sites_path = get_analysis_path(mirna, experiment, condition, "full_sites",ext="_%s-%s" %(start, stop))
        flanks_path = get_analysis_path(mirna, experiment, condition, "full_flanking",ext="_%s-%s" %(start, stop))
    else:
        reads_path = get_analysis_path(mirna, experiment, condition, "full_reads")
        sites_path = get_analysis_path(mirna, experiment, condition, "sites")
        flanks_path = get_analysis_path(mirna, experiment, condition, "flanking")

    with open(reads_path,"rb") as file_in_reads:
        with open(sites_path,"rb") as file_in_sites:
            results = multiprocess_file([file_in_reads, file_in_sites],
                                        readline_two, 1000000, 
                                        get_flanking_nucleotide,
                                        order_site_map,
                                        start, stop)

    # # Collect the output from each thread, which are the list
    # # of reads and the dictionary of read types.
    # site_threads = [i[0] for i in results]

    dict_threads = [i for i in results]
    # Construct the dictionary from the list of futures, and write it to
    # its output file.
    keys = set([j for i in dict_threads for j in i.keys()])

    kmers = [i for i in dict_threads[0].values()[0].keys()]

    values = [{kmer: sum([thread[key][kmer] for thread in dict_threads]) for kmer in kmers}
              for key in keys]

    count_flank_site_map = {key : value for (key,value) in zip(keys,values)}

    # No longer have input our output files open. Now write summmary file
    # to the same directory , ending in "summary.txt".
    with open(flanks_path,"wb") as file_out:
        file_out.write("\t%s\n" % ("\t".join(sorted(keys))))
        for kmer in sorted(kmers):
            file_out.write("%s:\t%s\n" % (
                kmer, "\t".join([str(count_flank_site_map[key][kmer]) for key in sorted(keys)])))


    # # Print the amount of time the script took to complete.
    print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()

