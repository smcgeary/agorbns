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




def check_library_sequence(sequences):
    """Generates library sequences for the reporter library.

    Args:
        site (str): A string giving the DNA sequence corresponding to the miRNA
            target site.
        n_seq (int): The number of sequences requested
    """
    df_cols = ["Variant", "Other sites", "Utr sites"]
    sequences_df = pd.DataFrame(data=None, index=range(len(sequences)), columns=df_cols)
    non_AA = get_kmer_list(2)[1:]
    splice_donor = r"GGGT[GA]AGT"
    bstxii_site = r"CCA[ACGT]{6,6}TGG"
    polya_site = "AATAAA"
    pumilio_site = r"TGTA[ACGT]ATAA"
    # splice_acceptor = "CAGG"

                #87654321#              # Distance from 5p end to 5p end of site.
    # const_5p = "CAAGTAAAGCGGCCGCACTCC"
    # const_5p = "TACTGCCTCCACGCTGATGG"      # Designed by Charlie and Dave
    # const_5p = "CTGCCTCCACGCTGATGGCG"      # Newest version, drop 2 nt on 5' end and add 3' CG.
                                     #123456789!1234567#         # Length = 17 nt         
    const_5p   =                     "TCTACAGTCCGACGATC"         # October 11th, has the sequencing primer
    rep_5p_utr =     "CCACGCTGATGGTTCAGAGTTCTACAGTCCGACGATC"
                     #123456789!123456789@123456789#1234567#
    new_rep_5p_utr = "CCACGCTG"     +      "CTACAGTCCGACGATC"
    new_rep_5p_utr = "CCACGCTG""GGTTCAGAGTTCTACAGTCCGACGATC"
    # const_3p = "GGCCAATGCCCTGGCTCACAAATAC"
    # const_3p = "CCAATGCCCTGGCTCACAAA"      # Designed by Charlie and Dave
    # const_3p = "GCCCAATGCCCTGGCTCACA"      # Newest version, add 5' GC and drop 2 nt on 3' end.
    const_3p    = "TACCAATGCCCTGGCTC"         # October 11th, 17nt with GC changed to TA.
    rep_3p_utr0     = "TACCAATGCCCTGGCTC"    +    "ACAAATACCACTGAGATCTTTTTCCCTCTGCCAAAAATTATGGGGACATCATGAAGCCCCTTGAGCATCTGACTTCTGGCT"
    new_rep_3p_utr0 = "TACCAATGCCCTGGCTCTGCCCTGGCTCACAAATACCACTGAGATCTTTTTCCCTCTGCCAAAAATTATGGGGACATCATGAAGCCCCTTGAGCATCTGACTTCTGGCT"
    # new_rep_3p_charlie = "TACCAATGCCCTGGCTCTGCCCTGGCTCACAAATACCACTGAGATCTTTTTCCCTCTGCCAAAAATTATGGGGACATCATGAAGCCCCTTGAGCATCTGACTTCTGGCTAATAAAGGAAATTTATTTTCATTGCAATAGTGTGTTGGAATTTTTTGTGTCTCTCA"


    rep_pA  = "AATAAAGGAAATTTATTTTCATTGCA" + "A"*20
    
    rep_3p_utr = rep_3p_utr0 + rep_pA
    new_rep_3p_utr = new_rep_3p_utr0 + rep_pA
                 #123456789!1234567#         # Length = 17 nt


    x_nucs_map = {"b" : ["C", "G", "T"],
                  "d" : ["A", "G", "T"],
                  "h" : ["A", "C", "T"],
                  "v" : ["A", "C", "G"]}
    time_start = time.time()
    time_start_i = time_start
    # mirna_site = "%s|%s" %(mirna, name)
    # min_acheived_other_sites = 2

    for row_i in range(len(sequences)):
        if row_i % 1000 == 0:
            print(row_i)
        total_seq = sequences[row_i][0]
        mir_site = sequences[row_i][1]
        # print(total_seq)
        # print(total_seq[:len(const_5p)])
        # print(const_5p)
        variant_seq = total_seq[len(const_5p):len(total_seq) - len(const_3p)]

        # print(" "*len(const_5p) + variant_seq)
        # print(new_rep_5p_0)
        # print(new_rep_3p_1)
        utr_seq_old = rep_5p_utr + variant_seq + rep_3p_utr
        utr_seq_new = new_rep_5p_utr + variant_seq + new_rep_3p_utr
        # print(utr_seq_old)
        # print(utr_seq_new)
        sites_found_utr_old = []
        sites_found_utr_new = []

            # utr_seq   = rep_5p_utr + seq_l + site_new + seq_r + rep_3p_utr
            # _reporter_variant = ReporterVariant(total_seq)
        _utr_variant_old = ReporterVariant(utr_seq_old)
        _utr_variant_new = ReporterVariant(utr_seq_new)
            # # Check variant for other miRNA sequences
        for mirna_i in MIRNAS:
            _mirna = Mirna(mirna_i)
            _sitelist = SiteList(_mirna, "paperfinal", 180)
            # _reporter_variant.get_all_sites_in_reporter_variant(_sitelist, ties=True)
            _utr_variant_old.get_all_sites_in_reporter_variant(_sitelist, ties=True)
            _utr_variant_new.get_all_sites_in_reporter_variant(_sitelist, ties=True)
            # site_names =["%s|%s" %(_mirna.name, i.name) for i in _reporter_variant.sites]
            site_names_utr_old =["%s|%s" %(_mirna.name, i.name) for i in _utr_variant_old.sites]
            site_names_utr_new =["%s|%s" %(_mirna.name, i.name) for i in _utr_variant_new.sites]
            # if site_names != ["%s|None" %(_mirna.name)]:
            #     lr = ["%s.%s" %(i.l, i.r) for i in _reporter_variant.sites]
            #     site_names_lr = ["|".join(i) for i in zip(site_names, lr)]
            #     sites_found += site_names_lr
            # if 
            # print(site_names_utr)
            # if site_names_utr != ["%s|None" %(_mirna.name)]:
            #     lr = ["%s.%s" %(i.l, i.r) for i in _utr_variant.sites]
            #     site_names_lr = ["|".join(i) for i in zip(site_names_utr, lr)]
            #     sites_found_utr += site_names_lr
            sites_found_utr_old += site_names_utr_old
            sites_found_utr_new += site_names_utr_new
        # print(sites_found_utr_old)
        # print(sites_found_utr_new)
        AAT_site = re.search("AATAAAG", utr_seq_new).span()
        # print(AAT_site)
        # print(utr_seq_new)
        # print(utr_seq_new[:AAT_site[1] - 2])
        utr_seq_polyA_l = utr_seq_new[:AAT_site[1] - 2]
        utr_seq_polyA_r = utr_seq_new[AAT_site[0] + 1:]
        # print(" "*(AAT_site[0] + 1) + utr_seq_new[AAT_site[0] + 1:])

        if sites_found_utr_old != sites_found_utr_new:
            print("it's different!")
            return
        # if (
        #     polya_site in utr_seq_polyA_l or
        #     polya_site in utr_seq_polyA_r):
        #     print("polyA site!")
        #     print(utr_seq_new)
        #     site_span = re.search(polya_site, utr_seq_new).span()

        #     print(" "*(site_span[0]) + polya_site)
        #     print(mir_site)


        #     return
        if re.search(splice_donor, utr_seq_new):
            print("splice donor site!")
            print(utr_seq_new)
            splice_span = re.search(splice_donor, utr_seq_new).span()

            print(" "*(splice_span[0] + 1) + splice_donor)
            return
        if re.search(pumilio_site, utr_seq_new):
            print("pumilio site!")
            return
            # span = re.search(bstxii_site, total_seq).span()
            # l_bound, r_bound = list(map(add, span, [1, -1]))
            # l_seq_check = total_seq[:r_bound]
            # r_seq_check = total_seq[l_bound:]
            # site_span_r = len(new_rep_5p_utr) + len(seq_l) + len(site_new) - 2
            # site_span_l = len(new_rep_5p_utr) + len(seq_l) + 1
            # alt_polyA_split_l = total_seq[:site_span_r]
            # alt_polyA_split_r = total_seq[site_span_l:]

            # if(polya_site in total_seq[:AAT_site[1] - 1])

            # if(polya_site in total_seq or (
            #             name == "AATAAAG" and
            #             polya_site not in alt_polyA_split_l and
            #             polya_site not in alt_polyA_split_r)) and
            #         not re.search(bstxii_site, l_seq_check) and
            #         not re.search(bstxii_site, r_seq_check)
            #                         ):




    #         # Successful stopping condition:
    #         # print(sites_found[0])
    #         # print(sites_found[0].split("|")[1])
    #         # print(name)
    #         if ( not proofread or (
    #                 len(sites_found) == 1 and
    #                 sites_found[0].split("|")[1] == name and    
    #                 (polya_site not in total_seq or (
    #                     name == "AATAAAG" and
    #                     polya_site not in alt_polyA_split_l and
    #                     polya_site not in alt_polyA_split_r)) and
    #                 not re.search(splice_donor, total_seq) and
    #                 not re.search(pumilio_site, total_seq) and
    #                 not re.search(bstxii_site, l_seq_check) and
    #                 not re.search(bstxii_site, r_seq_check)
    #                               )
    #            ):
    #             hold = False
    #             min_acheived_other_sites = 1
    #             if len(sites_found) == 1:
    #                 sites_found_str = sites_found[0]
    #             else:
    #                 sites_found_str = ",".join(sites_found)

    #             sites_found_utr_str = ",".join(sites_found_utr)
    #             sequence_row = pd.Series([total_seq, sites_found_str,
    #                                       sites_found_utr_str],
    #                                         index=df_cols)
    #             sequences_df.loc[row_i, ] = sequence_row 
    #         # Stopping condition at which point the best sequence is picked.
    #         elif attempts == max_attempts:
    #             print("at max_attempts")
    #             other_sites = list(trialseqs_df["Other sites"])
    #             # First filter by having the fewest number of sites:
    #             sites = [i.split(",") for i in other_sites]
    #             len_sites = [len(i) for i in sites]
    #             # print("len_sites")
    #             # print(len_sites)
    #             # print('sites')
    #             # print("\n".join([" ".join(i) for i in sites]))
    #             if (
    #                     min(len_sites) != min_acheived_other_sites or
    #                     (polya_site in total_seq and name != "AATAAAG") or
    #                     (name == "AATAAAG" and
    #                     (polya_site in alt_polyA_split_l or
    #                      polya_site in alt_polyA_split_r)) and

    #                     re.search(splice_donor, total_seq) or
    #                     re.search(pumilio_site, total_seq) or
    #                     re.search(bstxii_site, l_seq_check) or
    #                     re.search(bstxii_site, r_seq_check)
    #                ):
    #                 attempts = 0
    #                 # print("reject")
    #                 # print(min(len_sites))
    #                 # print(polya_site in total_seq)
    #                 # print(re.search(splice_donor, total_seq))
    #                 # print(re.search(pumilio_site, total_seq))
    #                 # print(re.search(bstxii_site, l_seq_check))
    #                 # print(re.search(bstxii_site, r_seq_check))
    #             else:
    #                 #
    #                 # print("proceed")
    #                 # print(trialseqs_df)
    #                 min_sites = [i for i in range(len(len_sites)) if len_sites[i] == 2]
    #                 # print("min sites")
    #                 # print(min_sites)
    #                 # trimmed_sites = [other_sites[i] for i in min_sites]
    #                 trialseqs_df = trialseqs_df.iloc[min_sites, ]
    #                 other_sites = list(trialseqs_df["Other sites"])

    #                 # print("trialseqs_df.shape")
    #                 # print(trialseqs_df.shape)

    #                 # Next filter by having a different miRNA than this one:
    #                 # mirnas_all = [[j.split("|")[0] for j in i.split(",")] for i in other_sites]
    #                 # sites_other = [[j.split("|")[1] for j in i.split(",")] for i in other_sites]
    #                 mirna_sites = [["|".join(j.split("|")[:2]) for j in i.split(",")] for i in other_sites]
    #                 mirna_sites_other = [[j for j in i if j != mirna_site] for i in mirna_sites]
    #                 # print("mirna_sites")
    #                 # print(mirna_sites)
    #                 # print('mirna_site')
    #                 # print(mirna_site)
    #                 # print("mirna_sites_other")
    #                 # print(mirna_sites_other)
    #                 mirna_other_index_1 = [i for i, mir_site_i in enumerate(mirna_sites_other) if mir_site_i != []]
    #                 mirna_other_index = [i for i in mirna_other_index_1 if mirna_sites_other[i][0].split("|")[0] != mirna]
    #                 # print("mir_index")
    #                 # print(mirna_other_index)
    #                 trialseqs_df = trialseqs_df.iloc[mirna_other_index, ]
    #                 # print("mirna_sites_other")
    #                 # print(mirna_sites_other)
    #                 mirna_sites_other_new = [mirna_sites_other[i][0] for i in mirna_other_index]
    #                 # print("mirna_sites_other_new")
    #                 # print(mirna_sites_other_new)
    #                 # other_sites = list(trialseqs_df["Other sites"])

    #                 # sites_new = [i.split(",") for i in other_sites]

    #                 # print(trialseqs_df)
    #                 # for i in mirna_sites_other_new:
    #                 #     print(i)
    #                 #     print(i.split("|"))
    #                 site_lengths = [len(Mirna(i.split("|")[0])[i.split("|")[1]]) for i in mirna_sites_other_new]
    #                 # print(site_lengths)
    #                 len_index = [i for i, len_i in enumerate(site_lengths) if len_i == min(site_lengths)][0]
    #                 # print(len_index)
    #                 # print(trialseqs_df.shape)
    #                 sequence_row = trialseqs_df.iloc[len_index,]
    #                 # print(sequence_row)
    #                 sequences_df.loc[row_i,] = sequence_row
    #                 hold = False

    #         elif ( not proofread or (
    #                 len(sites_found) == 2 and
    #                 (polya_site not in total_seq or (
    #                     name == "AATAAAG" and
    #                     polya_site not in alt_polyA_split_l and
    #                     polya_site not in alt_polyA_split_r)) and
    #                 not re.search(splice_donor, total_seq) and
    #                 not re.search(pumilio_site, total_seq) and
    #                 not re.search(bstxii_site, l_seq_check) and
    #                 not re.search(bstxii_site, r_seq_check)
    #                                 )
    #              ):

    #             # print(total_seq)
    #             sites_found_str = ",".join(sites_found)
    #             # print(sites_found_str)
    #             sites_found_utr_str = ",".join(sites_found_utr)

    #             row_df = pd.Series([total_seq, sites_found_str, sites_found_utr_str],
    #                                index=df_cols)
    #             # print("row_df:")
    #             # print(row_df)
    #             # print("trialseqs_df")
    #             # print(trialseqs_df)
    #             # print(attempts)
    #             trialseqs_df.loc[attempts,] = row_df
    #             attempts += 1
    #             # print("trialseqs_df post:")
    #             # print(trialseqs_df)
    #     # print(sequences_df)
    # time_stop = time.time()
    # time_del = time_stop - time_start
    # time_min = int(time_del / 60)
    # time_sec = time_min
    # print(time_del)

    return(sequences_df)


def main():
    time_start = time.time()
    arguments = ["miRNA"]
    args = parse_arguments(arguments)
    [mirna] = args

    suffix = "oct12_final"

    _mirna = Mirna(mirna)

    print(_mirna)

    # out_df_i = pd.concat([MirnaSite_df, SeqsOthersites_df], axis=1)
    # out_df_noproofread_i = pd.concat([MirnaSite_df, SeqsOthersites_noproofread_df], axis=1)

    out_dir = "variants/%s/%s" %(suffix, _mirna.name)
    ensure_directory(out_dir)
    
    in_path = "variants/%s/final_order.fa" %(suffix)

    print(in_path)

    with open(in_path) as file:
        variants = list(file)

    # print(variants[:10])


    sequences = [variants[2*i + 1].strip() for i in range(29992)]
    sites = [variants[2*i].strip()[2:] for i in range(29992)]


    # print(sequences)
    # print(sites)

    seq_site_tuples = [[sequences[i], sites[i]] for i in range(len(sequences))]


    check_library_sequence(seq_site_tuples)

    mirnas = ["miR-1", "let-7a", "miR-155", "miR-124", "lsy-6", "miR-7"]
    total_sites = 0

    x_nuc_map = {"A" : "b",
                 "C" : "d",
                 "G" : "h",
                 "T" : "v"}

    site_seqs_all = []
    seq_site_mirna_map = {mirna_i : None for mirna_i in mirnas}
    # A dictionary to tabulate the weird sites for each miRNA.
    unknownsites_mirna_map = {mirna_i : [] for mirna_i in mirnas}
    # Fill dataframe with miRNA, name, and sequence of each site.
    for mirna_i in mirnas:
        _mirna_i = Mirna(mirna_i)
        _sitelist = SiteList(_mirna_i, "paperfinal", 37)
        site_seqs = [(_mirna_i.name, i.name, i.seq) for i in _sitelist.sites]
        seq_site_map = {i.name : i.seq for i in _sitelist.sites}
        unknownsites = [i.name for i in _sitelist.sites if "mer" not in i.name]
        unknownsites_dict = {i : i for i in unknownsites}
        for i in unknownsites:
            for j in unknownsites:
                if i in j and len(i) != len(j):
                    outside = j.split(i)
        unknownsites_mirna_map[mirna_i] = unknownsites
        seq_site_mirna_map[mirna_i] = seq_site_map
        site_seqs_all += site_seqs
        for _site in _sitelist.sites:
            total_sites += 1
    sites_df = pd.DataFrame(site_seqs_all, columns=["miRNA", "site", "seq"])
    for i_row in range(sites_df.shape[0]):
        mirna_i, site, core = list(sites_df.iloc[i_row, ])
        s8mer = seq_site_mirna_map[mirna_i]["8mer"]
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



    sites_df.to_csv("180930_reporter_sites_check.txt", sep="\t")
    if len_vars:
        len_vars = int(len_vars)
        suffix = "len%snt_%s" %(len_vars, suffix)
    else:
        len_vars = 180
    each = int(3e4 / total_sites)
    # out_df = pd.DataFrame(None, columns=["miRNA", "Site", "Variant", "Other sites"])
    for i_row in range(sites_df.shape[0]):
        mirna_i, name, seq = list(sites_df.iloc[i_row, :])
        if mirna == mirna_i and (not single_site or name == single_site):
            print(name)
            # _mirna = Mirna(mirna)
            SeqsOthersites_df = generate_library_sequence(
                mirna, name, seq, each, freq_MrnaDinuc_map, len_vars=len_vars
            )

            SeqsOthersites_noproofread_df = generate_library_sequence(
                mirna, name, seq, each, freq_MrnaDinuc_map, len_vars=len_vars,
                proofread=False
            )


            MirnaSite_df = pd.DataFrame.from_dict({"miRNA" : [_mirna.name]*each,
                                                  "Site"  : [name]*each})
            out_df_i = pd.concat([MirnaSite_df, SeqsOthersites_df], axis=1)
            out_df_noproofread_i = pd.concat([MirnaSite_df, SeqsOthersites_noproofread_df], axis=1)

            out_dir = "variants/%s/%s" %(suffix, _mirna.name)
            ensure_directory(out_dir)
            
            out_path = "variants/%s/%s/%s_library_variants_%s.txt" %(suffix,
                                                              _mirna.name,
                                                              name, suffix)
            out_path_noproofread = "variants/%s/%s/%s_library_variants_%s_noproofread.txt" %(suffix,
                                                              _mirna.name,
                                                              name, suffix)

            print(out_path)
            # out_df_i.to_csv(out_path, sep="\t")
            # out_df_noproofread_i.to_csv(out_path_noproofread, sep="\t")

################################################################################

if __name__ == "__main__":
    main()

