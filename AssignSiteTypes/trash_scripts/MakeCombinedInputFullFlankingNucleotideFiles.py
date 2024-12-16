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



def get_flanking_nucleotide(read_seqs, order_site_map, start, stop, mirna, experiment):
    """Takes a read sequence and identifies the position of all site types'.

    Args:
        read_sequence: The read sequence.
        sites: A list of sequences corresponding to each site type.

    Returns:
        A tuple containing the list of reads and a dictionary of read types.
    """

    # Sub-function to be used on each site:
    count_flank_site_map = {
        site: {"".join(kmer[:2]+(".", )+kmer[2:]) : {"TGT" : 0, "ACA" : 0}
               for kmer in it.product(["A","C","G","T"],repeat=4)}
        for site in order_site_map.keys()}

    for i in read_seqs:
        print(i)
        print('was a read')
        seq, sites = (j.strip() for j in i)
        read_seq = seq[(26-5):(26+37+5)]
        barcode = read_seq[5+37:5+37+3]
        barcode_original = barcode
        if barcode == "TCG":
            barcode = "TGT"
            if mirna != "miR-1" or experiment != "equilibrium":
                read_seq = read_seq[:5+37]+"TGTTC"
        if barcode == "TGT" and mirna == "miR-1" and experiment == "equilibrium":
             read_seq = read_seq[:5+37]+"TCGTA"

        if sites != "None":
            coords = [i.split(":")[1] for i in sites.split(", ")]
            sites = [i.split(":")[0] for i in sites.split(", ")]
            coord_site_map = {site: coord for site, coord in zip(sites,coords)}
            ranks = [order_site_map[site] for site in sites]
            [(coord, site)] = [i for i in zip(coords, sites)
                               if order_site_map[i[1]] == min(ranks)]

            start = int(coord.split("-")[0])
            stop = int(coord.split("-")[1]) + 1
            left_start = start - 2
            right_stop = stop + 2
            if left_start >= 0 and right_stop <= 37 + start + stop:
                flank = read_seq[left_start:start]+"."+ read_seq[stop:right_stop]
                if site == "7mer-m8bT3":
                    print(site)
                    print(flank)
                    print(read_seq)
                    print(seq)
                    print(" "*(26-5) + read_seq)
                    print(" "*(26-5+start) + read_seq[start:stop])
                    print(" "*(26-5 + left_start) + read_seq[left_start:start]+"_"*(stop - start) + read_seq[stop:right_stop])
                    print(barcode_original)
                    print(barcode)
                count_flank_site_map[site][flank][barcode] += 1

    return count_flank_site_map

def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["miRNA","experient", "-sitelist","-start", "-stop", "-nb_binary"]
    mirna, experiment, sitelist, start, stop, nb = parse_arguments(arguments)
    if start == None:
        start = 0
    else:
        start = int(start)
    if stop == None:
        stop = 0
    else:
        stop = int(stop)

    # Load file with site types for the miRNA.
    if sitelist:
        print("hi")
        sites_file_name = ("/lab/bartel1_ata/mcgeary/computation/"
                           "AgoRBNS/AssignSiteTypes/sites.%s_%s.txt" % (mirna.split("-alt")[0], sitelist))
    else:
        sites_file_name = ("/lab/bartel1_ata/mcgeary/computation/"
                   "AgoRBNS/AssignSiteTypes/sites.%s.txt" % (mirna.split("-alt")[0]))


    print(sitelist)
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
    sites = sites_new   
    order_site_map = {site: i for i, site in enumerate(sites)}
    # Initialize the dictionary with all sites, and two more dictionaries,
    # one mapping to the site name from the sequence, and one from site name
    # to the count.

    if sitelist:
        print("hi")
        sites_file_name = ("/lab/bartel1_ata/mcgeary/computation/"
                           "AgoRBNS/AssignSiteTypes/sites.%s_%s.txt" % (mirna.split("-alt")[0], sitelist))
    else:
        sites_file_name = ("/lab/bartel1_ata/mcgeary/computation/"
                   "AgoRBNS/AssignSiteTypes/sites.%s.txt" % (mirna.split("-alt")[0]))


    extension = "_%s-%s" %(start, stop)
    if sitelist:
        extension = extension + "_" + sitelist


    # Get the path to the read file and to that of where the site labels will
    # be written.
    reads_path = get_analysis_path(mirna, experiment, "I_combined", "full_reads", ext = extension)
    sites_path = get_analysis_path(mirna, experiment, "I_combined", "full_sites", ext = extension)
    flanks_path = get_analysis_path(mirna, experiment, "I_combined", "full_flanking", ext = extension)

    with open(reads_path,"rb") as file_in_reads:
        with open(sites_path,"rb") as file_in_sites:
            results = multiprocess_test([file_in_reads, file_in_sites],
                                        readline_two, 1000000, 
                                        get_flanking_nucleotide,
                                        10,
                                        order_site_map,
                                        start,
                                        stop,
                                        mirna,
                                        experiment)

    # # Collect the output from each thread, which are the list
    # # of reads and the dictionary of read types.
    # site_threads = [i[0] for i in results]

    dict_threads = [i for i in results]

    # Construct the dictionary from the list of futures, and write it to
    # its output file.
    sites = list(set([j for i in dict_threads for j in i.keys()]))


    flanks = sorted(list(set([k for i in dict_threads for j in i.values() for k in j.keys()])))

    bcs = list(set([l for i in dict_threads for j in i.values() for k in j.values() for l in k.keys()]))
    print(bcs)
    print(sites)
    print(flanks)
    values = [{sites: sum([thread[site][flank]["TGT"] for thread in dict_threads]) for flank in flanks}
              for site in sites]

    count_flank_site_map = {site : value for (site,value) in zip(sites,values)}

    # # No longer have input our output files open. Now write summmary file
    # # to the same directory , ending in "summary.txt".
    # print(flanks_path)
    # with open(flanks_path,"wb") as file_out:
    #     file_out.write("\t%s\n" % ("\t".join(sorted(keys))))
    #     for kmer in sorted(kmers):
    #         file_out.write("%s:\t%s\n" % (
    #             kmer, "\t".join([str(count_flank_site_map[key][kmer]) for key in sorted(keys)])))


    # # # Print the amount of time the script took to complete.
    # print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()

