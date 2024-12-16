################################################################################
#MakeSiteTypeReadFiles.py
################################################################################
import imp # Used to import general.py
import time # Used for time
import numpy as np
import pandas as pd
import itertools as it
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

# FUNCTIONS

LIB_5p =     "GGGCAGAGTTCTACAGTCCGACGATC"
LIB_3p =          "TCGTATGCCGTCTTCTGCTTG" # Does not include TGT or ACA


def count_read_kmers(read_seqs, kmer_len, n_constant, rand_length):
    kmer_dict = {"TGT" : {kmer : 0 for kmer in get_kmer_list(kmer_len)},
                 "ACA" : {kmer : 0 for kmer in get_kmer_list(kmer_len)}}
    time_start = time.time()
    sys.stdout.flush()
    read_len = rand_length + 2*n_constant + 3
    for i_r, r in enumerate(read_seqs):
        read = r.strip()
        barcode = read[rand_length:]
        full_read = LIB_5p[-n_constant:] + read + LIB_3p[:n_constant]
        kmers = [full_read[i:i + kmer_len] for i in range(read_len - kmer_len + 1)]
        sys.stdout.flush()
        for kmer in kmers:
            kmer_dict[barcode][kmer] += 1
    return kmer_dict


def count_read_kmers_mask(read_seqs, kmer_len, n_constant, rand_length, mkmers):
    kmer_dict = {"TGT" : {kmer : 0 for kmer in get_kmer_list(kmer_len)},
                 "ACA" : {kmer : 0 for kmer in get_kmer_list(kmer_len)}}
    time_start = time.time()
    sys.stdout.flush()
    read_len = rand_length + 2*n_constant + 3
    for i_r, r in enumerate(read_seqs):
        read = r.strip()
        barcode = read[rand_length:]
        full_read = LIB_5p[-n_constant:] + read + LIB_3p[:n_constant]
        kmers = [full_read[i:i+kmer_len] for i in range(read_len - kmer_len + 1)]
        kmer_intersect = set(kmers).intersection(set(mkmers))
        instances = 0
        if kmer_intersect:
            read_list = list(full_read)
            for i, kmer in enumerate(kmers):
                if kmer in kmer_intersect:
                    # instances += 1
                    read_list[i:i+kmer_len] = "x"*(kmer_len)

            # if instances > 1:
            #     print('more than 1')
            full_read = "".join(read_list)
            kmers = [full_read[i:i+kmer_len] for i in range(read_len - kmer_len + 1)]
            kmers = [kmer for kmer in kmers if "x" not in kmer]
        for kmer in kmers:
            kmer_dict[barcode][kmer] += 1
    return kmer_dict






def main():
    time_start = time.time()
    arguments = ["miRNA", "condition_1", "condition_2", "n_constant",
                 "kmer_len", "-sitelist", "-z", "-jobs", "-test_binary"]
    args = parse_arguments(arguments)
    (mirna, condition_1, condition_2, n_constant, kmer_len, sitelist, z, jobs, test) = args
    _mirna = Mirna(mirna)

    if not kmer_len:
        kmer_len = 5
    if not jobs:
        jobs = 30
    if not sitelist:
        sitelist = "paperfinal"
    kmers = get_kmer_list(kmer_len)
    args = (kmer_len,)

    rand_length = 37

    if not jobs:
        jobs = 20
    args  = [int(kmer_len), int(n_constant), int(rand_length)]

    if z:
        z = int(z)
    else:
        z = 3    

    out_ext = "_%s_%s_k%s_z%s" %(n_constant, sitelist, kmer_len, z)

    if test:
        out_ext = "%s_test" %(ext)


    ############################################################################


    reads_path_input = get_analysis_path(mirna, "kinetics", "I", "reads")
    reads_path_1 = get_analysis_path(mirna, "kinetics", condition_1,
                                     "reads")
    if condition_2 == "pulse_v_chase":
        cond_out = "%s_pulse_v_chase" %(condition_1)
    else:
        reads_path_2 = get_analysis_path(mirna, "kinetics", condition_2,
                                         "reads")
        cond_out = "%s_v_%s_pulse" %(condition_1, condition_2)

    fullkmers_path = get_analysis_path(mirna, "kinetics", cond_out,
                                            "full_kmers", ext=out_ext)


    topkmers_path = get_analysis_path(mirna, "kinetics", cond_out,
                                            "top_kmers", ext=out_ext)


    print(reads_path_input)
    print(reads_path_1)
    if condition_2 != "pulse_v_chase":
        print(reads_path_2)
    print(topkmers_path)
   
    # Temporary while troubleshooting:

    # ________________________First round________________________________________
    threads_i = multiproc_file(reads_path_input,
                                   int(jobs), count_read_kmers, test, *args)
    threads_cond_1 = multiproc_file(reads_path_1, int(jobs), count_read_kmers,
                                    test, *args)
    
    if condition_2 != "pulse_v_chase":
        threads_cond_2 = multiproc_file(reads_path_2, int(jobs), count_read_kmers,
                                        test, *args)

    print(len(threads_i))
    print(threads_i[0]["TGT"].items()[:10])


    kmer_list = get_kmer_list(kmer_len)

    # Make the kmer dictionaries for the two halves of the data:
    kmer_dict_all = {
        "TGT_I" : {kmer : sum([i["TGT"][kmer] for i in threads_i])
                   for kmer in kmer_list},
        "ACA_I" : {kmer : sum([i["ACA"][kmer] for i in threads_i])
                   for kmer in kmer_list},
        "TGT_1" : {kmer : sum([i["TGT"][kmer] for i in threads_cond_1])
                   for kmer in kmer_list},
        "ACA_1" : {kmer : sum([i["ACA"][kmer] for i in threads_cond_1])
                   for kmer in kmer_list}
    }

    counts_df = pd.DataFrame(kmer_dict_all)
    counts_df = counts_df[["TGT_I", "ACA_I", "TGT_1", "ACA_1"]]

    print(counts_df.iloc[:10, :])


    if condition_2 != "pulse_v_chase":
        kmer_dict_all = {
            "TGT_2" : {kmer : sum([i["TGT"][kmer] for i in threads_cond_2])
                       for kmer in kmer_list},
            "ACA_2" : {kmer : sum([i["ACA"][kmer] for i in threads_cond_2])
                       for kmer in kmer_list}
        }
        counts_df_added = pd.DataFrame(kmer_dict_all)
        counts_df_added = counts_df_added[["TGT_2", "ACA_2"]]
        print(counts_df_added.iloc[:10, :])

        counts_df = pd.concat([counts_df, counts_df_added], axis=1,
                                 sort=False)





    freq_df = counts_df/counts_df.sum()
    R_df = pd.concat([freq_df["TGT_1"]/freq_df["TGT_I"],
                      freq_df["ACA_1"]/freq_df["ACA_I"]], axis=1, sort=False)
    R_df.columns = ["TGT_1", "ACA_1"]
    print(R_df.iloc[:10, :])

    if condition_2 != "pulse_v_chase":
        R_df = pd.concat(
            [R_df, freq_df["TGT_2"]/freq_df["TGT_I"],
             freq_df["ACA_2"]/freq_df["ACA_I"]],
            axis=1, sort=False
        )
        R_df.columns = ["TGT_1", "ACA_1", "TGT_2", "ACA_2"]
        R_df_final = R_df["TGT_1"]/R_df["TGT_2"]
    else:
        R_df_final = R_df["TGT_1"]/R_df["ACA_1"]

    R_df_final.columns = ["0"]
    R_df_sorted = R_df_final.sort_values(axis=0, ascending=False)

    R_mean = R_df_final.mean()
    R_std  = R_df_final.std()

    z_df = (R_df_sorted - R_mean)/R_std

    print(R_df_sorted.iloc[:20])

    z_cutoff = R_mean + z*R_std
    print(R_mean)
    print(R_std)
    print(z_cutoff)

    top_R = R_df_sorted.iloc[0]
    top_kmer = R_df_sorted.index[0]

    top_z = z_df.iloc[0]

    print(top_R)
    print(top_kmer)

    R_df_final.to_csv(fullkmers_path, sep="\t")

    topkmers = []
    
    # ________________________Subsequent round__________________________________
    while top_R >= z_cutoff:
        topkmers.append([top_kmer, top_R, top_z])
        print(topkmers)

        mkmers = [topkmer[0] for topkmer in topkmers]
        args  = [int(kmer_len), int(n_constant), int(rand_length), mkmers]

        threads_i = multiproc_file(reads_path_input, int(jobs),
                                   count_read_kmers_mask, test, *args)
        threads_cond_1 = multiproc_file(reads_path_1, int(jobs),
                                        count_read_kmers_mask, test, *args)
        
        if condition_2 != "pulse_v_chase":
            threads_cond_2 = multiproc_file(reads_path_2, int(jobs),
                                            count_read_kmers_mask, test, *args)

        # Make the kmer dictionaries for the two halves of the data:
        kmer_dict_all = {
            "TGT_I" : {kmer : sum([i["TGT"][kmer] for i in threads_i])
                       for kmer in kmer_list if kmer not in mkmers},
            "ACA_I" : {kmer : sum([i["ACA"][kmer] for i in threads_i])
                       for kmer in kmer_list if kmer not in mkmers},
            "TGT_1" : {kmer : sum([i["TGT"][kmer] for i in threads_cond_1])
                       for kmer in kmer_list if kmer not in mkmers},
            "ACA_1" : {kmer : sum([i["ACA"][kmer] for i in threads_cond_1])
                       for kmer in kmer_list if kmer not in mkmers}
        }

        counts_df = pd.DataFrame(kmer_dict_all)
        counts_df = counts_df[["TGT_I", "ACA_I", "TGT_1", "ACA_1"]]

        if condition_2 != "pulse_v_chase":
            kmer_dict_all = {
                "TGT_2" : {kmer : sum([i["TGT"][kmer] for i in threads_cond_2])
                           for kmer in kmer_list if kmer not in mkmers},
                "ACA_2" : {kmer : sum([i["ACA"][kmer] for i in threads_cond_2])
                           for kmer in kmer_list if kmer not in mkmers}
            }
            counts_df_added = pd.DataFrame(kmer_dict_all)
            counts_df_added = counts_df_added[["TGT_2", "ACA_2"]]
            print(counts_df_added.iloc[:10, :])

            counts_df = pd.concat([counts_df, counts_df_added], axis=1,
                                     sort=False)

        freq_df = counts_df/counts_df.sum()
        R_df = pd.concat([freq_df["TGT_1"]/freq_df["TGT_I"],
                          freq_df["ACA_1"]/freq_df["ACA_I"]], axis=1,
                          sort=False)

        R_df.columns = ["TGT_1", "ACA_1"]
        print(R_df.iloc[:10, :])

        if condition_2 != "pulse_v_chase":
            R_df = pd.concat(
                [R_df, freq_df["TGT_2"]/freq_df["TGT_I"],
                 freq_df["ACA_2"]/freq_df["ACA_I"]],
                axis=1, sort=False
            )
            R_df.columns = ["TGT_1", "ACA_1", "TGT_2", "ACA_2"]
            R_df_final = R_df["TGT_1"]/R_df["TGT_2"]
        else:
            R_df_final = R_df["TGT_1"]/R_df["ACA_1"]

        R_df_final.columns = ["0"]
        R_df_sorted = R_df_final.sort_values(axis=0, ascending=False)

        R_mean = R_df_final.mean()
        R_std  = R_df_final.std()

        z_df = (R_df_sorted - R_mean)/R_std


        print(R_df_sorted.iloc[:20])
        top_kmer = R_df_sorted.index[0]

        top_R = R_df_sorted.iloc[0]
        top_z = z_df.iloc[0]
        print(top_kmer)
        print(top_R)

        print("TOP KMERS")
        print("\n".join(["\t".join([str(i) for i in topkmer])
                        for topkmer in topkmers]))

    print("TOP KMERS FINAL")
    print("\n".join(["\t".join([str(i) for i in topkmer])
                    for topkmer in topkmers]))


    with open(topkmers_path, "w+") as file_out:
        file_out.write("\n".join(["\t".join([str(i) for i in topkmer])
                                 for topkmer in topkmers]))
        file_out.write("\n")

    print(topkmers_path)

    print("finished writing kmer file.")

    print_time_elapsed(time_start)

################################################################################

if __name__ == "__main__":
    main()

