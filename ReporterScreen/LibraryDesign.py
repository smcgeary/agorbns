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

pd.set_option('max_columns', 600)
pd.set_option("max_colwidth", 600)
# FUNCTIONS

MIRNAS = ["miR-1", "let-7a", "miR-155", "miR-124", "lsy-6", "miR-7"]

def count_utr_sequence(utrs, split=False):
    """Calculates the dinucleotides frequencies of a DNA sequence.

    Args:
        utrs (list): A list of strings, containing only the characters
            A, C, G, and U.
    Returns:
        dict: A dictionary with keys that are the 16 DNA dinucleotides and the 
            keys are the corresponding frequencies within the entire list of
            sequences. The dictionary sums to 1.0.
    """
    dinucs = get_kmer_list(2)
    counts_dinuc_map = {dinuc : 0 for dinuc in dinucs}
    for utr in utrs:
        for i in range(len(utr) - 1):
            counts_dinuc_map[utr[i:i + 2]] += 1
    totals = sum(counts_dinuc_map.values())
    freq_dinuc_map = {dinuc : float(counts_dinuc_map[dinuc])/totals
                       for dinuc in dinucs}
    return freq_dinuc_map

def count_each_utr_sequence(utrs):
    """Calculates the dinucleotides frequencies of a DNA sequence.

    Args:
        utrs (list): A list of strings, containing only the characters
            A, C, G, and U.
    Returns:
        dict: A dictionary with keys that are the 16 DNA dinucleotides and the 
            keys are the corresponding frequencies within the entire list of
            sequences. The dictionary sums to 1.0.
    """
    dinucs = get_kmer_list(2)
    utrs_DinucFreq_map = ["None"] * len(utrs)
    tally = 0
    for utr_i, utr in enumerate(utrs):
        # print("utr_i")
        # print(utr_i)
        # print('utr')
        # print(utr)
        counts_dinuc_map = {dinuc : 0 for dinuc in dinucs}
        for i in range(len(utr) - 1):
            counts_dinuc_map[utr[i:i + 2]] += 1
        totals = sum(counts_dinuc_map.values())
        utrs_DinucFreq_map[utr_i] = {
            dinuc : float(counts_dinuc_map[dinuc])/totals for dinuc in dinucs
        }
        tally += 1
    utrs_DinucFreq_df = pd.DataFrame.from_dict(utrs_DinucFreq_map)
    return utrs_DinucFreq_df


def generate_utr_sequences(num_seq, len_seq, frac_dinuc_dict):
    """Generates a DNA sequences using the provided dinucleotide frequencies.

    Args:
        utrs (list): A list of strings, containing only the characters
            A, C, G, and U.
    Returns:
        dict: A dictionary with keys that are the 16 DNA dinucleotides and the 
            keys are the corresponding frequencies within the entire list of
            sequences. The dictionary sums to 1.0.
    """

    out = []
    a = frac_dinuc_dict.keys()
    p = frac_dinuc_dict.values()
    len_seq_start = len_seq

    for i in range(num_seq):
        # print(i)
        start = True
        # print(i)
        len_seq = len_seq_start
        utr = ""
        # print("out")
        # print(out)
        while len_seq != 0:
            # print("utr")
            # print(utr)
            if start:
                dinuc = np.random.choice(a, p=p)
                start = False
                utr = dinuc
            else: 
                all_use = [i for i in zip(a, p) if i[0][0] == utr[-1]]
                a_use = [i[0] for i in all_use]
                p_use = [i[1] for i in all_use]
                p_sum = np.sum(p_use)
                p_use = [i/p_sum for i in p_use]
                dinuc = np.random.choice(a_use, p=p_use)
                utr += dinuc[1]
            len_seq -= 1
        out += [utr]
    if num_seq == 1:
        return out[0]
    else:
        return out

def generate_library_sequence(site, n_seq, len_seq, freq_dinuc_map, max_attempts=50):
    """Generates library sequences for the reporter library.

    Args:
        site (str): A string giving the DNA sequence corresponding to the miRNA
            target site.
        n_seq (int): The number of sequences requested
    """
    df_cols = ["Variant", "Other sites"]
    sequences_df = pd.DataFrame(data=None, index=range(n_seq), columns=df_cols)
    non_AA = get_kmer_list(2)[1:]

    x_nucs_map = {"b" : ["C", "G", "T"],
                  "d" : ["A", "G", "T"],
                  "h" : ["A", "C", "T"],
                  "v" : ["A", "C", "G"]}
    time_start = time.time()
    for row_i in range(n_seq):
        # time_start = time.time()
        # print("variants made: %s" %(row_i))
        # Establish stopping condition, and reset attempts
        hold = True
        attempts = 0
        trialseqs_df = pd.DataFrame(data=None, index=range(max_attempts),
                                    columns=df_cols)
        while hold:
            sites_found = []
            site_list = list(site)
            for i, nuc in enumerate(site_list):
                if nuc == "B":
                    dinuc_list = list(random.choice(non_AA))
                    site_list[i] = dinuc_list[0]
                    site_list[i + 1] = dinuc_list[1]
                elif nuc in x_nucs_map.keys():
                    site_list[i] = random.choice(x_nucs_map[nuc])
            site_new = "".join(site_list)
            # print(site)
            # print(site_new)

            len_non_site = len_seq - len(site_new)
            len_l = len_non_site / 2
            len_r = len_non_site - len_l
            seq_l = generate_utr_sequences(1, len_l, freq_dinuc_map)
            seq_r = generate_utr_sequences(1, len_r, freq_dinuc_map)

            total_seq = seq_l + site_new + seq_r
            _reporter_variant = ReporterVariant(total_seq)

            # Check variant for other miRNA sequences
            for mirna in MIRNAS:
                _mirna = Mirna(mirna)
                _sitelist = SiteList(_mirna, "paperfinal", len_seq)
                _reporter_variant.get_all_sites_in_reporter_variant(_sitelist, ties=True)
                site_names =["%s|%s" %(_mirna.name, i.name) for i in _reporter_variant.sites]
                if site_names != ["%s|None" %(_mirna.name)]:
                    lr = ["%s.%s" %(i.l, i.r) for i in _reporter_variant.sites]
                    site_names_lr = ["|".join(i) for i in zip(site_names, lr)]
                    # print(site_names_lr)
                    sites_found += site_names_lr
            if len(sites_found) == 1:
                hold = False
                sites_found_str = "None"
                sequence_row = pd.Series([total_seq, sites_found_str],
                                            index=df_cols)
                # print("sequence_row")
                # print(sequence_row)
                # print("sequences_df, pre")
                # print(sequences_df)
                sequences_df.loc[row_i, ] = sequence_row 
                # print("sequences_df, post")
                # print(sequences_df)
            elif attempts == max_attempts:
                hold = False

                other_sites = list(trialseqs_df["Other sites"])
                # First filter by having the fewest number of sites:
                sites = [i.split(",") for i in other_sites]
                len_sites = [len(i) for i in sites]
                if min(len_sites) != 1:
                    print("reset")
                    attempts = 0
                    break
                min_sites = [i for i in range(len(len_sites)) if len_sites[i] == min(len_sites)]
                # print("min_sites:")
                # print(min_sites)
                trimmed_sites = [other_sites[i] for i in min_sites]
                trialseqs_df = trialseqs_df.iloc[min_sites, ]

                # Next filter by having a different miRNA than this one:
                mirnas_other = [[j.split("|")[0] for j in i.split(",")] for i in trimmed_sites]
                mir_index = [i for i, mir in enumerate(mirnas_other) if mir != mirna]
                trialseqs_df = trialseqs_df.iloc[mir_index, ]


                other_sites_2 = [trimmed_sites[i] for i, mir in enumerate(mirnas_other) if mir != mirna]
                mirnas_other2 = [j.split("|")[0] for i in other_sites_2 for j in i.split(",")]

                sites_other = [j.split("|")[1] for i in other_sites_2 for j in i.split(",")]
                # print(mirnas_other2)
                site_lengths = [len(Mirna(i[0])[i[1]]) for i in zip(mirnas_other2, sites_other)]
                # print(site_lengths)
                len_index = [i for i, len_i in enumerate(site_lengths) if len_i == min(site_lengths)][0]
                # print(len_index)
                sequence_row = trialseqs_df.iloc[len_index,]
                # print(sequence_row)
                sequences_df.loc[row_i,] = sequence_row
            else:
                # print(total_seq)
                sites_found_str = ",".join(sites_found[1:])
                # print(sites_found_str)
                row_df = pd.Series([total_seq, sites_found_str],
                                   index=df_cols)
                # print("row_df:")
                # print(row_df)
                # print("trialseqs_df")
                # print(trialseqs_df)
                # print(attempts)
                trialseqs_df.loc[attempts,] = row_df
                attempts += 1
                # print("trialseqs_df post:")
                # print(trialseqs_df)
        # print(sequences_df)
    time_stop = time.time()
    time_del = time_stop - time_start
    time_min = int(time_del / 60)
    time_sec = time_min
    print(time_del)

    return(sequences_df)


def main():
    plt.ion()
    plt.cla()
    plt.close()

    # 1: Read in the UTRs from Kathy's table, tally up all the dinucleotide
    # frequencies, to generate a probability table to generate random 3'UTR
    # sequence.

    # Make the list of UTRS:
    utr_path = ("/lab/bartel4_ata/kathyl/RNA_Seq/analysis/" + 
                "data/no_baseline_analysis/merged.txt")

    with open(utr_path) as file:
        # Drop the first row because this is the header file.
        print(file.readline())
        print(file.readline())

    with open(utr_path) as file:
        # Drop the first row because this is the header file.
        mrna_utrs = [i.split("\t")[18] for i in file.readlines()[1:]]
    with open(utr_path) as file:
        # Drop the first row because this is the header file.
        mrna_names = [i.split("\t")[0] for i in file.readlines()[1:]]

    with open(utr_path) as file:
        # Drop the first row because this is the header file.
        mrna_exps = [math.exp(np.mean([float(j) for j
                                        in i.split("\t")[2:17]]))
                     for i in file.readlines()[1:]]

    with open(utr_path) as file:
        # Drop the first row because this is the header file.
        mrna_utrl = [float(i.split("\t")[20]) for i in file.readlines()[1:]]

    mrna_exp_df = pd.Series(mrna_exps, index=mrna_names)
    mrna_utrl_df = pd.Series(mrna_utrl, index=mrna_names)
    # for i, ex in enumerate(mrna_exps):
    #     print(mrna_names[i])
    #     print(ex)
    #     print(float(ex))


    # 2: Test that the two functions, count_utr_sequences(), and
    # generate_utr_sequences(), work, by using generate_utr_sequences() to
    # simulate a large number of utr sequences, then calculate the dinucleotide
    # freuqncies within the simulated utrs, and compare them to the frequencies
    # originally calculated for the mRNA utr sequences.
    freq_MrnaDinuc_map = count_utr_sequence(mrna_utrs)

    # freq_MrnaDinucEach_df = count_each_utr_sequence(mrna_utrs)
    # freq_MrnaDinucEach_df.index = mrna_names

    # freq_MrnaDinucEachExp_df = freq_MrnaDinucEach_df.multiply(mrna_exp_df, axis=0)
    # freq_MrnaDinucEachExp_df = freq_MrnaDinucEachExp_df / sum(mrna_exp_df)

    # freq_MrnaDinucEachUtrLen_df = freq_MrnaDinucEach_df.multiply(mrna_utrl_df, axis=0)
    # freq_MrnaDinucEachUtrLen_df = freq_MrnaDinucEachUtrLen_df / sum(mrna_utrl_df)

    # mrna_utrexp_df = mrna_utrl_df*mrna_exp_df
    # freq_MrnaDinucEachExpUtrLen_df = freq_MrnaDinucEach_df.multiply(mrna_utrexp_df, axis=0)
    # freq_MrnaDinucEachExpUtrLen_df = freq_MrnaDinucEachExpUtrLen_df / sum(mrna_utrexp_df)


    # sim_utrs = generate_utr_sequences(10, 1000, freq_MrnaDinuc_map)
    # freq_SimDinuc_map = count_utr_sequence(sim_utrs)
    # frac_sim = pd.DataFrame.from_dict(freq_SimDinuc_map, orient="index",
    #                                   columns=["sim"])

    # dinuc_freq_df = pd.concat(
    #     [pd.DataFrame.from_dict(i[0], orient="index", columns=[i[1]]) for i in
    #      zip([freq_MrnaDinuc_map, freq_SimDinuc_map], ["mrna_l", "sim_l"])] + 
    #     [freq_MrnaDinucEach_df.mean(axis=0),
    #      freq_MrnaDinucEachExp_df.sum(axis=0),
    #      freq_MrnaDinucEachUtrLen_df.sum(axis=0),
    #      freq_MrnaDinucEachExpUtrLen_df.sum(axis=0)],
    #     axis=1, sort=True)
    # dinuc_freq_df.columns.values[2:] = ["mrna_", "mrna_e", "mrna_l2", "mrna_el"]
    # print(dinuc_freq_df)

    # fig, axarr = plt.subplots(2, 2)
    # lty1, = axarr[0, 0].plot([0, 0.15], [0, 0.15], label="x=y")
    # pts1, = axarr[0, 0].plot(dinuc_freq_df["mrna_"], dinuc_freq_df["mrna_e"], "bo",
    #                  label="By mRNA expression")
    # lty2, = axarr[0, 1].plot([0, 0.15], [0, 0.15], label="x=y")
    # pts2, = axarr[0, 1].plot(dinuc_freq_df["mrna_"], dinuc_freq_df["mrna_l"], "ro",
    #                  label="By UTR length")

    # lty3, = axarr[1, 0].plot([0, 0.15], [0, 0.15], label="x=y")
    # pts3, = axarr[1, 0].plot(dinuc_freq_df["mrna_"], dinuc_freq_df["mrna_el"], "go",
    #                  label="By both")

    # axarr[1, 1].axis("off")

    # axarr[1, 0].set_xlabel("Equal")
    # axarr[1, 0].set_ylabel("Weighted")
    # axarr[1, 1].legend(handles=[pts1, pts2, pts3, lty1], )


    mirnas = ["miR-1", "let-7a", "miR-155", "miR-124", "lsy-6", "miR-7"]
    total_sites = 0

    x_nuc_map = {"A" : "b",
                 "C" : "d",
                 "G" : "h",
                 "T" : "v"}

    site_seqs_all = []
    seq_site_mirna_map = {mirna : None for mirna in mirnas}
    # A dictionary to tabulate the weird sites for each miRNA.
    unknownsites_mirna_map = {mirna : [] for mirna in mirnas}

    # Fill dataframe with miRNA, name, and sequence of each site.
    for mirna in mirnas:
        _mirna = Mirna(mirna)
        _sitelist = SiteList(_mirna, "paperfinal", 37)
        site_seqs = [(_mirna.name, i.name, i.seq) for i in _sitelist.sites]
        seq_site_map = {i.name : i.seq for i in _sitelist.sites}
        unknownsites = [i.name for i in _sitelist.sites if "mer" not in i.name]
        unknownsites_dict = {i : i for i in unknownsites}
        for i in unknownsites:
            for j in unknownsites:
                if i in j and len(i) != len(j):
                    outside = j.split(i)
        unknownsites_mirna_map[mirna] = unknownsites
        seq_site_mirna_map[mirna] = seq_site_map
        site_seqs_all += site_seqs
        for _site in _sitelist.sites:
            total_sites += 1
    sites_df = pd.DataFrame(site_seqs_all, columns=["miRNA", "site", "seq"])

    for i_row in range(sites_df.shape[0]):
        mirna, site, core = list(sites_df.iloc[i_row, ])
        s8mer = seq_site_mirna_map[mirna]["8mer"]
        # Part of script where the non-pairing characters are added to the 
        # sequence. Necessary for 1.) Seed sites and 2.) 3' sites, and 3.) the
        # AA-[site] sites of miR-124.
        # 1.) Seed sites:
        if "7mer-m8" in site:
            xA1 = x_nuc_map[s8mer[-1]]
            sites_df.loc[i_row, "seq"] =       core + xA1
        if "7mer-A1" in site:
            xm8 = x_nuc_map[s8mer[0]]
            sites_df.loc[i_row, "seq"] = xm8 + core
        if site in ["6mer", "6mer-bG7"]:
            xA1 = x_nuc_map[s8mer[-1]]
            xm8 = x_nuc_map[s8mer[0]]
            sites_df.loc[i_row, "seq"] = xm8 + core + xA1
        if site == "6mer-A1":
            xm7 = x_nuc_map[s8mer[1]]
            sites_df.loc[i_row, "seq"] = xm7 + core
        if site == "6mer-m8":
            xm2 = x_nuc_map[s8mer[-2]]
            sites_df.loc[i_row, "seq"] =       core + xm2
        if site == "5mer-A1":
            xm6 = x_nuc_map[s8mer[2]]
            sites_df.loc[i_row, "seq"] = xm6 + core
        if site == "5mer-m2.6":
            xA1 = x_nuc_map[s8mer[-1]]
            xm7 = x_nuc_map[s8mer[1]]
            sites_df.loc[i_row, "seq"] = xm7 + core + xA1
        if site == "5mer-m3.7":
            xm2 = x_nuc_map[s8mer[-2]]
            xm8 = x_nuc_map[s8mer[0]]
            sites_df.loc[i_row, "seq"] = xm8 + core + xm2
        if site == "5mer-m8":
            xm3 = x_nuc_map[s8mer[-3]]
            sites_df.loc[i_row, "seq"] =       core + xm3
        # 2.) 3'-only sites of miR-155, miR-124, and lsy-6:
        if "mer-" in site:
            k, lr = site.split("mer-")
            if k in ["9", "10", "11"]:
                # Split schematic: m$(l).$(r)w$(p)
                #                 0|1
                #                   $(l).$(r)w$(p)
                #                           0|1
                #                   $(l).$(r)
                #                      0|1
                l, r = [int(i) for i
                        in lr.split("m")[1].split("w")[0].split(".")]
                _mirna = Mirna(mirna)
                if r < len(_mirna):
                    x_5p = x_nuc_map[get_rc(_mirna.seq[r])]
                else:
                    x_5p = ""
                x_3p = x_nuc_map[get_rc(_mirna.seq[l - 2])]
                site_full = x_5p + core + x_3p
                sites_df.loc[i_row, "seq"] = site_full
        # 3.) AA-[site] sites of miR-124:
        if mirna == "miR-124":
            _mirna = Mirna(mirna)
            # Check x positions at 3p end of AA sites:
            if site[:3] == "AA-":
                site_strip = site[3:]
                if "7mer-m8" in site_strip:
                    x_3p = x_nuc_map[get_rc(_mirna.seq[0])]
                elif "6mer-m8" in site_strip:
                    x_3p = x_nuc_map[get_rc(_mirna.seq[1])]
                elif "5mer-m8" in site_strip:
                    x_3p = x_nuc_map[get_rc(_mirna.seq[2])]
                else:
                    x_3p = ""
                sites_df.loc[i_row, "seq"] = core + x_3p
            # Check that non AA-[site] sites of miR-124 that have a paired
            # AA-[site] are not allowed to have AA in the flanking dinucleotide
            # context.
            elif "AA-" + site in seq_site_mirna_map["miR-124"].keys():
                sites_df.loc[i_row, "seq"] = "Bb" + sites_df.loc[i_row, "seq"]



    # Stuff to print to make sure that each row is good.
    # print(sites_df.iloc[:50,])
    # print("\n\n")
    # print(sites_df.iloc[51:100,])
    # print("\n\n")
    # print(sites_df.iloc[101:150,])
    # print("\n\n")
    # print(sites_df.iloc[151:,])
    sites_df.to_csv("180930_reporter_sites_check.txt", sep="\t")

    each = int(3e4 / total_sites)
    # print(each)
    # for i_row in range(sites_df.shape[0]):
    # each = 2
    sequence_length = 150
    print(each)
    out_df = pd.DataFrame(None, columns=["miRNA", "Site", "Variant", "Other sites"])
    for i_row in range(sites_df.shape[0]):
        print(i_row)
        # print("out_df pre:")
        # print(out_df)
        mirna, name, seq = list(sites_df.iloc[i_row, :])
        print(mirna)
        print(name)
        _mirna = Mirna(mirna)
        SeqsOthersites_df = generate_library_sequence(seq, each, sequence_length,
                                                  freq_MrnaDinuc_map)
        MirnaSite_df = pd.DataFrame.from_dict({"miRNA" : [mirna]*each,
                                              "Site"  : [name]*each})
        out_df_i = pd.concat([MirnaSite_df, SeqsOthersites_df], axis=1)
        ensure_directory("variants/%s" %(mirna))
        out_df_i.to_csv("variants/%s/%s_library_variants_test.txt" %(mirna, name), sep="\t")
        # print("out_df_i:")
        # print(out_df_i)
        out_df = out_df.append(out_df_i, ignore_index=True, sort=False)
        # print("out_df:")
        # print(out_df)
        # print("\n".join(mir_site_seqs))
        out_df.to_csv("variants/total_test.txt", sep="\t")
    # print(total_sites)

################################################################################

if __name__ == "__main__":
    main()

