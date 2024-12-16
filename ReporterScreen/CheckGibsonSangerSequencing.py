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


def check_library_sequence(i_seq, sequence, seq_name_map, name_seq_map, seq_type="gibson"):
    """Generates library sequences for the reporter library.

    Args:
        site (str): A string giving the DNA sequence corresponding to the miRNA
            target site.
        n_seq (int): The number of sequences requested
    """
    # df_cols = ["Variant", "Other sites", "Utr sites"]
    # sequences_df = pd.DataFrame(data=None, index=range(len(sequences)), columns=df_cols)
    # non_AA = get_kmer_list(2)[1:]
    # splice_donor = r"GGGT[GA]AGT"
    # bstxii_site = r"CCA[ACGT]{6,6}TGG"
    # polya_site = "AATAAA"
    # pumilio_site = r"TGTA[ACGT]ATAA"
    # splice_acceptor = "CAGG"

                #87654321#              # Distance from 5p end to 5p end of site.
    # const_5p = "CAAGTAAAGCGGCCGCACTCC"
    # const_5p = "TACTGCCTCCACGCTGATGG"      # Designed by Charlie and Dave
    # const_5p = "CTGCCTCCACGCTGATGGCG"      # Newest version, drop 2 nt on 5' end and add 3' CG.
                                     #123456789!1234567#         # Length = 17 nt         
    const_5p   =                         "TCTACAGTCCGACGATC"         # October 11th, has the sequencing primer
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
    if seq_type == "gibson":
        if i_seq == "6":
            sequence = sequence[:56] + "-"*36 + sequence[56:]
        elif i_seq == "8":
            sequence = sequence[:157] + "-"*16 + sequence[157:]

    sequence_old = sequence
    sanger_sequence = sequence[10:250]
    _variant = ReporterVariant(sanger_sequence)

    # print(i_seq + "." + "_"*170)
    # print("")
    num_sites = 0

    if re.search(const_5p, sanger_sequence):
        [span_5p_l, span_5p_r] = re.search(const_5p, sanger_sequence).span()
        pre_seq = " "*span_5p_l + const_5p
    else:
        span_5p_r = 0
        pre_seq = ""

    if re.search(const_3p, sanger_sequence):
        [span_3p_l, span_3p_r] = re.search(const_3p, sanger_sequence).span()
        seq_3p_sanger = sanger_sequence[span_3p_l:span_3p_l + 20]
    else:
        span_3p_l = len(sanger_sequence)
        span_3p_r = len(sanger_sequence)
        seq_3p_sanger = ""


    site_sequences_all = []
    ls_all = []
    rs_all = []

    for mirna in MIRNAS:
        _mirna = Mirna(mirna)
        _sitelist = SiteList(_mirna, "paperfinal", len(sanger_sequence))
        _variant.get_all_sites_in_reporter_variant(_sitelist, ties=True)
        if _variant.sites[0].name != "None":
            sites_all = list(_variant.sites)
            names = [i.name for i in sites_all]
            ls = [i.l for i in sites_all]
            rs = [i.r for i in sites_all]
            seqs = [i.seq for i in sites_all]

            site_sequences_all += seqs
            ls_all += ls
            rs_all += rs
            # for i in range(len(ls)):
            #     print(pre_seq +
            #           "n"*(ls[i] - len(pre_seq)) +
            #           seqs[i] +
            #           "n"*(span_3p_l - rs[i]) + 
            #           seq_3p_sanger + 
            #           " "*(250 - span_3p_r) +
            #           mirna + " " + names[i])
            num_sites += len(sites_all)
    num_sites = len(set(site_sequences_all))

    if num_sites == 0:
        print("\n"*8)
        return
    elif i_seq == "2" and seq_type == "gibson":
        print("Sanger sequence:")
        best_variant = seq_name_map["miR-155_7mer-A1_25"]
        print(sanger_sequence[:150] + "...")
        tups = zip(list(sanger_sequence), list(" "*span_5p_l + best_variant))
        tup_check = [i == j.upper() for i, j in tups]
        matches = "".join([[" ", "|"][0 + i ] for i in tup_check])
        print(matches[:150] + "...")
        k_len = 54
        variant_sequence = " "*span_5p_l + best_variant
        print(variant_sequence[:150] + "...")
        print("best sequence:")
        print(" ")

        buffer_str = "\t"*2 + "..."
        print(buffer_str + sanger_sequence[120:])
        print(buffer_str + matches[120:])
        print(buffer_str + variant_sequence[120:])
     

        print("TACCAATGCCCTGGCTC" in sanger_sequence)
        print("TACCAATGCCCTGGCTC" in sequence_old)
        print(sanger_sequence[span_3p_l:span_3p_r])

    else:
        query_sequence = sanger_sequence[ls_all[0] - 10:rs_all[0] + 10]

        styled_sequence = (
            sanger_sequence[:span_5p_r] +
            sanger_sequence.lower()[span_5p_r:ls_all[0]] + 
            sanger_sequence[ls_all[0]:rs_all[0]] + 
            sanger_sequence.lower()[rs_all[0]:span_3p_l] + 
            sanger_sequence[span_3p_l:])


        for key in name_seq_map:
            if re.search(query_sequence, key):
                key_old = key
                if i_seq == "9":
                    key = key[:27] + "-" + key[27:]
                span_l, span_r = re.search(query_sequence, key).span()
                # print("."*(ls_all[0]))
                # print("*"*(span_l + 10))
                left_space = " "*(ls_all[0] - span_l - 10)
                # left_space = " "*0
                key_left = key[:(span_5p_r - span_5p_l)]
                rand_left = key[(span_5p_r - span_5p_l):(span_l + 10)]
                site_seq = key[(span_l + 10) : (span_r - 10)]
                var_3p_span = re.search("TACCAATGCCCTGGCTC", key).span()[0]
                rand_right = key[(span_r - 10): var_3p_span]
                if seq_type == "gibson":
                    key_right = new_rep_3p_utr0
                else:
                    key_right = rep_3p_utr0
                # print(left_space + key + " "*(250 - len(left_space) - len(key)) + name_seq_map[key])
                aligned_variant = (
                    left_space + key_left + rand_left.lower() + site_seq +
                    rand_right.lower() + key_right
                )
                
                aligned_variant = aligned_variant[:len(styled_sequence)] + "\t" + name_seq_map[key_old]


                tups = zip(list(_variant.seq), list(aligned_variant))
                tup_check = [i == j.upper() for i, j in tups]
                matches = "".join([[" ", "|"][0 + i ] for i in tup_check])
                print(styled_sequence[:120] + "...")
                print(matches[:120])
                print(aligned_variant[:120] + "...")
                print("")
                buffer_str = "\t"*2 + "..."
                print(buffer_str + styled_sequence[120:])
                print(buffer_str + matches[120:])
                print(buffer_str + aligned_variant[120:])
                print("")
                # print(left_space + " "*span_l + query_sequence)
        return



def main():
    time_start = time.time()

    path = "ReporterScreen/GibsonSangerSequencing_181126"

    print(path)

    files = os.listdir(path)
    print(files)

    new_rep_5p_utr = "CCACGCTG""GGTTCAGAGTTCTACAGTCCGACGATC"


    const_3p    = "TACCAATGCCCTGGCTC"         # October 11th, 17nt with GC changed to TA.


    print("R" + "\t" + "_"*(40 - len(new_rep_5p_utr)) + new_rep_5p_utr + "."*30 + const_3p)

    num_i = 1

    reporter_path = "ReporterScreen/final_order_twist.fa"

    with open(reporter_path) as fa_file_in:
        file_fa = list(fa_file_in)
        names = [file_fa[i*2].strip()[2:] for i in range(len(file_fa)/2)]
        seqs = [file_fa[i*2 + 1].strip() for i in range(len(file_fa)/2)]

    print(names[:10])
    print(seqs[:10])
    seq_name_map = {names[i] : seqs[i] for i in range(len(names))}
    name_seq_map = {seqs[i] : names[i] for i in range(len(names))}

    files = [file for file in files if ".seq" in file]
    file_order = [int(file.split("TRL-")[1].split("_")[0]) for file in files]

    files_sorted = [file for _, file in
                    sorted(zip(file_order,files), key=lambda pair: pair[0])]


    for file in files_sorted:
        if ".seq" in file:
            full_path = "%s/%s" %(path, file)
            num_i = file.split("TRL-")[1].split("_")[0]
            with open(full_path) as file_in:
                sequence_full = "".join([i.strip() for i in list(file_in)])
                # print(sequence_full)
                check_library_sequence(str(num_i), sequence_full, seq_name_map,
                                       name_seq_map)
                # num_i += 1


    path = "ReporterScreen/LigationSangerSequencing_181126"

    print(path)

    files = os.listdir(path)
    print(files)

    new_rep_5p_utr = "CCACGCTG""GGTTCAGAGTTCTACAGTCCGACGATC"


    const_3p    = "TACCAATGCCCTGGCTC"         # October 11th, 17nt with GC changed to TA.


    print("R" + "\t" + "_"*(40 - len(new_rep_5p_utr)) + new_rep_5p_utr + "."*30 + const_3p)

    num_i = 1

    reporter_path = "ReporterScreen/final_order_twist.fa"

    with open(reporter_path) as fa_file_in:
        file_fa = list(fa_file_in)
        names = [file_fa[i*2].strip()[2:] for i in range(len(file_fa)/2)]
        seqs = [file_fa[i*2 + 1].strip() for i in range(len(file_fa)/2)]

    print(names[:10])
    print(seqs[:10])
    seq_name_map = {names[i] : seqs[i] for i in range(len(names))}
    name_seq_map = {seqs[i] : names[i] for i in range(len(names))}

    files = [file for file in files if ".seq" in file]
    file_order = [int(file.split("clone")[1].split("_")[0]) for file in files]

    files_sorted = [file for _, file in
                    sorted(zip(file_order,files), key=lambda pair: pair[0])]


    for file in files_sorted:
        if ".seq" in file:
            full_path = "%s/%s" %(path, file)
            num_i = file.split("clone")[1].split("_")[0]
            with open(full_path) as file_in:
                sequence_full = "".join([i.strip() for i in list(file_in)])
                # print(sequence_full)
                check_library_sequence(str(num_i), sequence_full, seq_name_map,
                                       name_seq_map, seq_type="ligation")






    return


################################################################################

if __name__ == "__main__":
    main()

