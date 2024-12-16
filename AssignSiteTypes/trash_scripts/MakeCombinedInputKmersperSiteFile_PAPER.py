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
from general import get_rc, seq_mirna_map, parse_arguments, print_time_elapsed, get_analysis_path, multiprocess_file, multiprocess_test, readline, readline_two

# This script reads one read file and outputs multiple .txt files to be used in all
# downstream applications. Inputs are:

# mirna: The miRNA in the experiment to be processed.
# experiment: The type of experiment to be processed.
# condition: The sampel type within the experiment to be processed.

# Functions from general this uses:


# Constants

# FUNCTIONS



def get_kmers(read_seqs, order_site_map, k, start, stop, mirna, experiment):
    """Takes a read sequence and identifies the position of all site types'.

    Args:
        read_sequence: The read sequence.
        sites: A list of sequences corresponding to each site type.

    Returns:
        A tuple containing the list of reads and a dictionary of read types.
    """

    # Sub-function to be used on each site:
    count_kmer_site_map = {
        site: {"".join(kmer): {"TGT" : 0, "ACA" : 0}
               for kmer in it.product(["A","C","G","T"],repeat=k)}
        for site in order_site_map.keys()}
    for i in read_seqs:
        seq, sites = (j.strip() for j in i)
        read = seq.strip()
        barcode = read[26 + 37 : 26 + 37 + 3]
        if barcode == "TCG":
            barcode = "TGT"
            if mirna != "miR-1" or experiment != "equilibrium":
                read = read[:26 + 37]+"TGTTCGTATGCCGTCTTCTGCTTG"
        if barcode == "TGT" and mirna == "miR-1" and experiment == "equilibrium":
            read = read[:26 + 37]+"TCGTATGCCGTCTTCTGCTTG"

        read_seq = read[(26 - start - 2) : (26 + 37 + stop + 2)]
        kmers = [read_seq[i:(i+k)] for i in range(2,len(read_seq)-k-2+1)]

        if sites != "None":
            coords = [i.split(":")[1] for i in sites.split(", ")]
            sites = [i.split(":")[0] for i in sites.split(", ")]
            coord_site_map = {site: coord for site, coord in zip(sites,coords)}
            ranks = [order_site_map[site] for site in sites]
            [(coord, site)] = [i for i in zip(coords, sites)
                               if order_site_map[i[1]] == min(ranks)]
            # print(read)
            # print(read[:26] + "_"*37 + read[(26 + 37):])
            # print(" "*(26 - start) + read_seq[2:-2])
            # print(" "*(26 - start - 2 + index) + read_seq[index : end])
            # print(" "*(26 - start - 2 + left_start) + read_seq[left_start : index] + "_"*(end - index) + read_seq[end : right_stop])
        else:
            site = "None"
        for i, kmer in enumerate(kmers):
            # print(" "*(26 - start + i) + kmer)
            count_kmer_site_map[site][kmer][barcode] += 1


            # print(site)
            # print(flank)
            # count_flank_site_map[site][flank][barcode] += 1

    return count_kmer_site_map

def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["miRNA","experient", "condition", "sitelist", "start", "stop", "k", "-nb_binary"]
    mirna, experiment, condition, sitelist, start, stop, k, nb = parse_arguments(arguments)

    # Load file with site types for the miRNA.
    sites_file_name = ("/lab/bartel1_ata/mcgeary/computation/"
                       "AgoRBNS/AssignSiteTypes/sites.%s_%s.txt" % (mirna.split("-alt")[0], sitelist))

    print(sites_file_name)

    with open(sites_file_name) as file_in:
        sites = file_in.read().split("\n")
        # This swaps the four "Centered-..." for one single "Centered" entry in the list.
    sites_new = []
    for i in sites:
        if "Centered" not in i:
            sites_new.append(i)
        elif add:
            sites_new.append("Centered")
            add = False
    sites = sites_new + ["None"]  
    order_site_map = {site: i for i, site in enumerate(sites)}
    # Initialize the dictionary with all sites, and two more dictionaries,
    # one mapping to the site name from the sequence, and one from site name
    # to the count.

    # Get the path to the read file and to that of where the site labels will
    # be written.
    extension = "_%s-%s_%s" %(start, stop, sitelist)

    # Get the path to the read file and to that of where the site labels will
    # be written.
    reads_path = get_analysis_path(mirna, experiment, condition, "full_reads")
    sites_path = get_analysis_path(mirna, experiment, condition, "full_sites", ext = extension)
    kmers_path = get_analysis_path(mirna, experiment, condition, "full_sitekmers", ext = extension + "_k%s" %(k))

    with open(reads_path,"rb") as file_in_reads:
        with open(sites_path,"rb") as file_in_sites:
            results = multiprocess_file([file_in_reads, file_in_sites],
                                        readline_two,
                                        int(1e6), 
                                        get_kmers,
                                        order_site_map,
                                        int(k),
                                        int(start),
                                        int(stop),
                                        mirna,
                                        experiment)

    # # Collect the output from each thread, which are the list
    # # of reads and the dictionary of read types.

    dict_threads = [i for i in results]

    # Construct the dictionary from the list of futures, and write it to
    # its output file.
    keys = set([j for i in dict_threads for j in i.keys()])

    kmers = [i for i in dict_threads[0].values()[0].keys()]

    values = [{kmer: (tuple([sum([thread[key][kmer]["TGT"] for thread in dict_threads]),
                       sum([thread[key][kmer]["ACA"] for thread in dict_threads])])) for kmer in kmers}
              for key in keys]
    count_kmer_site_map = {key : value for (key,value) in zip(keys,values)}
    # No longer have input our output files open. Now write summmary file
    # to the same directory , ending in "summary.txt".
    print(kmers_path)
    with open(kmers_path,"wb") as file_out:
        file_out.write("\t%s\n" % ("\t".join(sorted(keys))))
        for kmer in sorted(kmers):
            file_out.write("%s:\t%s\n" % (
                kmer, "\t".join([str(count_kmer_site_map[key][kmer]) for key in sorted(keys)])))


    # # Print the amount of time the script took to complete.
    print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()

