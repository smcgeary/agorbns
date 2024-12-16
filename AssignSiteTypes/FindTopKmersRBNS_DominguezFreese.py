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
def count_read_kmers(read_seqs, kmer_len):
    kmer_dict = {kmer: 0 for kmer in get_kmer_list(kmer_len)}
    time_start = time.time()
    sys.stdout.flush()
    read_len = len(read_seqs[0].strip())
    for i_r, r in enumerate(read_seqs):
        read = r.strip()
        sys.stdout.flush()
        kmers = [read[i:i+kmer_len] for i in range(read_len - kmer_len + 1)]
        for kmer in kmers:
            kmer_dict[kmer] += 1
    return kmer_dict

def count_read_kmers_mask(read_seqs, kmer_len, mkmers):
    kmer_dict = {kmer: 0 for kmer in get_kmer_list(kmer_len)}
    time_start = time.time()
    sys.stdout.flush()
    read_len = len(read_seqs[0].strip())
    for i_r, r in enumerate(read_seqs):
        read = r.strip()
        kmers = [read[i:i+kmer_len] for i in range(read_len - kmer_len + 1)]
        kmer_intersect = set(kmers).intersection(set(mkmers))
        instances = 0
        if kmer_intersect:
            read_list = list(read)
            for i, kmer in enumerate(kmers):
                if kmer in kmer_intersect:
                    # instances += 1
                    read_list[i:i+kmer_len] = "x"*(kmer_len)

            # if instances > 1:
            #     print('more than 1')
            read = "".join(read_list)
            kmers = [read[i:i+kmer_len] for i in range(read_len - kmer_len + 1)]
            kmers = [kmer for kmer in kmers if "x" not in kmer]
        for kmer in kmers:
            kmer_dict[kmer] += 1
    return kmer_dict



def main():
    time_start = time.time()
    arguments = ["rbp", "condition", "-kmer_len", "-jobs", "-test_binary"]
    args = parse_arguments(arguments)
    (rbp, condition, kmer_len, jobs, test) = args    
    if not kmer_len:
        kmer_len = 5
    if not jobs:
        jobs = 30
    kmers = get_kmer_list(kmer_len)
    args = (kmer_len,)
    ############################################################################
    if test:
        out_ext = "_test"
    else:
        out_ext = ""
    reads_path_input = get_analysis_path_burge(rbp, "equilibrium", "I", "reads")
    reads_path = get_analysis_path_burge(rbp, "equilibrium", condition, "reads")
    topkmers_path = get_analysis_path_burge(rbp, "equilibrium", condition,
                                            "top_kmers", ext=out_ext)

    # Temporary while troubleshooting:

    kwargs = {"halve_reads" : True}
    # ________________________First round________________________________________
    threads_i = multiproc_file(reads_path_input,
                                   int(jobs), count_read_kmers, test, *args, **kwargs)
    threads_rbp = multiproc_file(reads_path, int(jobs), count_read_kmers,
                                   test, *args, **kwargs)

    split_i = len(threads_i)/2
    split_rbp = len(threads_rbp)/2

    # Make the kmer dictionaries for the two halves of the data:
    kmer_dict_all = [
        {kmer : [sum([i[kmer] for i in threads_i]),
                 sum([i[kmer] for i in threads_rbp])]
         for kmer in kmers},
        {kmer : [sum([i[kmer] for i in threads_i[:split_i]]),
                 sum([i[kmer] for i in threads_rbp[:split_rbp]])]
         for kmer in kmers},
        {kmer : [sum([i[kmer] for i in threads_i[split_i:]]),
                 sum([i[kmer] for i in threads_rbp[split_rbp:]])]
         for kmer in kmers}
    ]

    counts_df_list = [pd.DataFrame.from_dict(kmer_dict, orient="index",
                                             columns=["I", "RBP"]).astype(np.float64)
                      for kmer_dict in kmer_dict_all]

    for counts_df in counts_df_list:
        counts_df.sort_index(inplace=True)


    total_kmers_list = [counts_df.sum().astype(np.float64) for counts_df in counts_df_list]

    print(total_kmers_list)
    freq_df_list = [counts_df/total_kmers
                    for counts_df, total_kmers in zip(counts_df_list, total_kmers_list)]



    R_df_list = [freq_df.iloc[:, 1]/freq_df.iloc[:, 0]
                 for freq_df in freq_df_list]

    Rcutoff_df_list = [R_df.mean() + 3*R_df.std()
                 for R_df in R_df_list]

    best_kmer_list = [R_df.idxmax() for R_df in R_df_list]

    kmer_weight_list = [R_df.loc[best_kmer] - 1
                        for R_df, best_kmer in zip(R_df_list, best_kmer_list)]
    

    print(kmer_weight_list)
    best_kmer_tup_list = []

    round = 2
    #____________________________Next rounds____________________________________
    while (kmer_weight_list[1] >= Rcutoff_df_list[1] - 1 and
           kmer_weight_list[2] >= Rcutoff_df_list[2] - 1):
        if best_kmer_list[1] != best_kmer_list[2]:
            print("someting is wrong")
            return
        best_kmer_tup_list.append((best_kmer_list[1],
                                   np.mean(kmer_weight_list[1:])))

        print(best_kmer_tup_list)
        print("ROUND %s____________________________________" %(round))
        masked_kmers = [i[0] for i in best_kmer_tup_list]
        args = [kmer_len, masked_kmers]

        threads_i = multiproc_file(reads_path_input, int(jobs),
                                   count_read_kmers_mask, test, *args, **kwargs)
        threads_rbp = multiproc_file(reads_path, int(jobs),
                                     count_read_kmers_mask, test, *args,
                                     **kwargs)

        # Make the kmer dictionaries for the two halves of the data:
        kmer_dict_all = [
            {kmer : [sum([i[kmer] for i in threads_i]),
                     sum([i[kmer] for i in threads_rbp])]
             for kmer in kmers if kmer not in masked_kmers},
            {kmer : [sum([i[kmer] for i in threads_i[:split_i]]),
                     sum([i[kmer] for i in threads_rbp[:split_rbp]])]
             for kmer in kmers if kmer not in masked_kmers},
            {kmer : [sum([i[kmer] for i in threads_i[split_i:]]),
                     sum([i[kmer] for i in threads_rbp[split_rbp:]])]
             for kmer in kmers if kmer not in masked_kmers}
        ]

        counts_df_list = [pd.DataFrame.from_dict(kmer_dict, orient="index",
                                             columns=["I", "RBP"]).astype(np.float64)
                          for kmer_dict in kmer_dict_all]

        for counts_df in counts_df_list:
            counts_df.sort_index(inplace=True)

        print([counts_df.sum() for counts_df in counts_df_list])
        print(total_kmers_list)
        freq_df_list_1 = [counts_df/total_kmers
                        for counts_df, total_kmers in zip(counts_df_list, total_kmers_list)]

        freq_df_list_2 = [counts_df/counts_df.sum()
                        for counts_df in counts_df_list]


        print(freq_df_list_1[1].iloc[:10, :])
        print(freq_df_list_2[1].iloc[:10, :])
        R_df_list_1 = [freq_df.iloc[:, 1]/freq_df.iloc[:, 0]
                     for freq_df in freq_df_list_1]

        R_df_list_2 = [freq_df.iloc[:, 1]/freq_df.iloc[:, 0]
                     for freq_df in freq_df_list_2]


        print(R_df_list_1[1].iloc[:10])
        print(R_df_list_2[1].iloc[:10])


        best_kmer_list = [R_df.idxmax() for R_df in R_df_list_1]

        kmer_weight_list_1 = [R_df.loc[best_kmer] - 1
                            for R_df, best_kmer in zip(R_df_list_1, best_kmer_list)]

        kmer_weight_list_2 = [R_df.loc[best_kmer] - 1
                            for R_df, best_kmer in zip(R_df_list_2, best_kmer_list)]
        
        print(kmer_weight_list_1)
        print(kmer_weight_list_2)
        kmer_weight_list = kmer_weight_list_2

        print("\n".join(["\t".join([str(i) for i in best_kmer_tup])
                         for best_kmer_tup in best_kmer_tup_list]))



        round += 1

    print(best_kmer_tup_list)

    topkmers_df = pd.DataFrame(best_kmer_tup_list)
    print(topkmers_df)
    topkmers_df.index = topkmers_df[0]
    print(topkmers_df)
    del topkmers_df.index.name
    print(topkmers_df)
    topkmers_df = topkmers_df.iloc[:, 1:]
    print(topkmers_df)
    topkmers_df.columns = ["R-1"]
    print(topkmers_df)

    print(topkmers_path)
    topkmers_df.to_csv(topkmers_path, sep="\t")


    print_time_elapsed(time_start)
################################################################################

if __name__ == "__main__":
    main()

