################################################################################
#LibraryDesign.py
################################################################################
import imp # Used to import general.py
imp.load_source("general",
                ("/lab/bartel1_ata/mcgeary/computation/"
                 "AgoRBNS/general/general.py"
                )
               )
imp.load_source("RBNS_methods",
                ("/lab/bartel1_ata/mcgeary/computation/"
                 "AgoRBNS/general/RBNS_methods.py"
                )
               )
from general import *
from RBNS_methods import *
from sitetypes import get_seq_site_map
from operator import add

pd.set_option('max_columns', 600)
pd.set_option("max_colwidth", 600)
# FUNCTIONS

MIRNAS = ["miR-1", "let-7a", "miR-155", "miR-124", "lsy-6", "miR-7"]


#FINAL DESIGN CONSTRAINTS
kLEN = 180


def assign_variant(read_seqs, mirna_site_map, seq100nt_variant_map):
    time_start = time.time()
    output_dic = make_variant_count_dictionary(mirna_site_map)
    unmapped = 0
    for i_r, r in enumerate(read_seqs):
        # Get the read number:
        variant = r.strip()
        if variant in seq100nt_variant_map:
            mirna, site, i, sequence = seq100nt_variant_map[variant]
            output_dic[mirna][site][int(i)] += 1
        else:
            unmapped += 1
    return([output_dic, unmapped])




def main():
    time_start = time.time()
    # Parse arguments:
    path = "ReporterScreen/Lim_UTR_sequences.txt"

    with open(path) as file_in:
        line = file_in.readline()
        line = file_in.readline()
        while line:
            (utr, mirna, mir1_rep, mir124_rep, sequence) = line.strip().split("\t")
            line = file_in.readline()
            print("______________________")
            for mirna in ["miR-1", "miR-124"]:
                _mirna = Mirna(mirna)
                _sitelist = SiteListRep(_mirna, "paperfinal")
                _utr = Utr(sequence)
                _utr.get_all_sites_in_utr(_sitelist)
                add_sites_from_utr(_sitelist, _utr)
                seqs = [_site.seq for _site in _utr.sites]
                names = [_site.name for _site in _utr.sites]
                ls = [_site.l for _site in _utr.sites]
                rs = [_site.r for _site in _utr.sites]
                overlaps = _utr.site_overlaps()
                delete_overlap = False
                for i, overlap in enumerate(overlaps):
                    if overlap > 0:
                        delete_overlap = True
                        overlap_pair = names[i:i+2]
                        if overlap_pair[0] == "8mer":
                            i_del = i + 1
                        else:
                            i_del = i
                if delete_overlap:
                    full_len = len(seqs)
                    seqs = [seqs[i] for i in range(full_len) if i != i_del]
                    names = [names[i] for i in range(full_len) if i != i_del]
                    ls = [ls[i] for i in range(full_len) if i != i_del]
                    rs = [rs[i] for i in range(full_len) if i != i_del]
                distances = [l - r for (l, r) in zip(ls[1:], rs[:-1])]
                out_string = str(ls[0]) + "_"
                for i, name in enumerate(names[:-1]):
                    out_string = out_string + name + "_" + str(distances[i]) + "_"
                if rs[-1]:
                    out_string = out_string + names[-1] + "_" + str(len(sequence) - rs[-1])
                print(out_string)


    # mirna, experiment, condition, rep, jobs, test = args
    # # Allocate job argument:
    # if not jobs:
    #     jobs = 20
    # # Update output extensions depending on test argument:
    # ext = ",%s" %(rep)
    # if test:
    #     ext_out = "%s_test" %(ext)
    # else:
    #     ext_out = ext
    # reads_path = get_analysis_path(mirna, experiment, condition, "reads",
    #                                ext=ext)

    # variant_path = get_analysis_path(mirna, experiment, condition, "counts",
    #                                  ext=ext_out)


    # print(reads_path)

    # # Get the variant dictionary:
    # mirna_site_map = get_site_sequence_map_dictionary()
    # seq100nt_variant_map = get_variant_dictionary(mirna_site_map)
    # # Assign the arguments for the function.
    # args = [mirna_site_map, seq100nt_variant_map]
    # # Use multiprocessing on the function assign variants:
    # threads = multiproc_file(reads_path, int(jobs), assign_variant, test,
    #                               *args)

    # print("_"*20 + "Finished with multiprocessing." + "_"*20)
    # print("")
    # output = make_variant_count_dictionary(mirna_site_map)

    # # Preallocate the mapped variable:
    # mapped = 0
    # for mirna in MIRNAS:
    #     for site in output[mirna]:
    #         for i in range(184):
    #             # Sum all of the instances of one variant across all the
    #             # microprocessed threads:
    #             mapped_i = sum([thread[0][mirna][site][i]
    #                             for thread in threads])
    #             # Add this sum to the apporpirate part of the totalled
    #             # dictionary.
    #             output[mirna][site][i] = mapped_i
    #             # Add to the tally of mapped variants.
    #             mapped += mapped_i
    # # Add the second multiprocesses variable together to get the unmapped
    # # variable.
    # unmapped = sum([thread[1] for thread in threads])
    # print("mapped: %s" %(mapped))
    # print("unmapped: %s" %(unmapped))

    # print(variant_path)
    # with open(variant_path, "w+") as file_out:
    #     for mirna in MIRNAS:
    #         site_order = get_mirna_site_order(mirna)
    #         for site in site_order:
    #             text_out = "\t".join(
    #                 [mirna, site] + [str(i) for i in output[mirna][site]]
    #             )
    #             print(text_out)
    #             file_out.write(text_out + "\n")




    print_time_elapsed(time_start)

    return


################################################################################

if __name__ == "__main__":
    main()

