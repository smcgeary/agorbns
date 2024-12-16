################################################################################
#GenerateSiteTypeKds.py
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

################################################################################
# FUNCTIONS
def main():
    time_start = time.time()
    # Define all the relevant arguments.
    arguments = ["miRNA", "condition", "kmer_len", "-n_ex", "-alt_mirna", "-alt_exp",
                 "-alt_n_constant",
                 "-n", "-check", "-weird", "-uniq_binary", "-buffer3p_binary",
                 "-sub_binary", "-I0_binary", "-final_binary",
                 "-combined_binary", "-cutoffdil_binary"]
    (mirna, condition, kmer_len, n_ex,
     alt_mirna,  alt_exp, alt_n_constant,
     n, check, weird, uniq, buffer, sub,
     I0, final, combined, cutoffdil     ) = parse_arguments(arguments)
    print(sub)
    # Initialize objects:
    # if mirna == "miR-7" and experiment == "equilibrium2_nb":
    #     _mirnas = [Mirna(mirna)  for mirna in
    #                ["miR-7-23nt", "miR-7-24nt", "miR-7-25nt"]]
    # else:
    _mirnas = [Mirna(mirna)]
    if (alt_n_constant):
        n_constant = alt_n_constant
    else:
        n_constant = "5"
    if n:
        n = int(n)
    else:
        n = 10
    if final:
        site_string = "paperfinal"
    elif weird:
        site_string = "papercutoffweirdsite%smer" %(weird)
    elif cutoffdil:
        site_string = "papercutoffdil"
    else:
        site_string = "papercutoff"
    _sitelist = SiteList(_mirnas[0], site_string, 1)
    exclude_list = [i.name for i in _sitelist.sites]
    if n_ex:
        n_ex = int(n_ex)
    else:
        n_ex = len([i for i in exclude_list if i != ""])
# if experiment in ["equil_mmseed_nb", "equil_seed_nb"]:
#     rand_length = 38
# else:
    rand_length = 37    

    print("condition")
    print(condition)
    print("n_constant")
    print(n_constant)

    if mirna == "miR-7-23nt":
        experiment = "equilibrium2_nb"
        if condition == "4":
            condition = "12.6"
    if mirna == "miR-7-24nt":
        experiment = "equilibrium_tp"
    elif mirna == "miR-124":
        experiment = "equilibrium_2_tp"
    elif mirna == "miR-1":
        experiment = "equilibrium_tp"
    else:
        experiment = "equilibrium"
    if alt_exp:
        experiment = alt_exp
    # if mirna == "miR-1":
    #     buffer = True
    # else:
    #     buffer = False

    print(buffer)
    _kmerlist_I = KmerList(int(kmer_len), int(n_constant), int(rand_length),
                           fast=True)
    if combined:
        _kmerlist_I.read(_mirnas[0], experiment, "I_combined", n_constant, n_ex,
                         uniq=uniq, buffer=buffer, final=final,
                         cutoffdil=cutoffdil)
    elif I0:
        _kmerlist_I.read(_mirnas[0], experiment, "I", n_constant, 0,
                         uniq=uniq, buffer=buffer, final=final,
                         cutoffdil=cutoffil)
    else:
        _kmerlist_I.read(_mirnas[0], experiment, "I", n_constant, n_ex,
                         uniq=uniq, buffer=buffer, final=final,
                         cutoffdil=cutoffdil)
    print("loop")

    for _mirna in _mirnas:
        _kmerlist = KmerList(int(kmer_len), int(n_constant), int(rand_length), fast=True)
        _kmerlist.read(_mirna, experiment, condition, n_constant, n_ex,
                       uniq=uniq, buffer=buffer, final=final,
                       cutoffdil=cutoffdil)
        _kmerlist_R_I = _kmerlist/_kmerlist_I
        kmer_tups = zip(get_kmer_list(int(kmer_len)),_kmerlist.kmers, _kmerlist_I.kmers)

        if alt_mirna:
            _mirna = Mirna(alt_mirna)
        print("%s:" %(_mirna.name))

        _kmerlist_R_I.head_with_mirna_sites(_mirna, n=n, sub=sub)
        _kmerlist_R_I.top_and_adjacent()
        _kmerlist_R_I.top_substitute()

        print("right match")
        print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("TTCACCGCTGT")])
        print("right mismatch")
        print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("TTCACCGCTGA")])
        print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("TTCACCGCTGC")])
        print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("TTCACCGCTGG")])
        print("left match")
        print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("ATTCACCGCTG")])
        print("left mismatch")
        print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("CTTCACCGCTG")])
        print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("GTTCACCGCTG")])
        print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("TTTCACCGCTG")])

        # print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("TATCGCCTTTA")])
        # print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("CATCGCCTTTA")])
        # print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("GATCGCCTTTA")])

        # print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("ATCGCCTTTAA")])
        # print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("ATCGCCTTTAT")])
        # print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("ATCGCCTTTAC")])
        # print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("ATCGCCTTTAG")])


        # print("CTTCCGCTGC")
        # print("TTCCGCTGCA")
        # print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("GAATTCACCG")])
        # print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("GTATTCACCG")])
        # print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("GCATTCACCG")])
        # print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("GGATTCACCG")])


        # print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("AATTCACCGC")])
        # print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("TATTCACCGC")])
        # print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("CATTCACCGC")])
        # print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("GATTCACCGC")])

        # print("TTCCGCTGCC")
        # print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("TTCCGCTGCC")])
        # print("TTCCGCTGCG")
        # print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("TTCCGCTGCG")])
        # print("TTCCGCTGCT")
        # print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("TTCCGCTGCT")])
        # print("ACTTCCGCTG")
        # print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("ACTTCCGCTG")])
        # print("CCTTCCGCTG")
        # print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("CCTTCCGCTG")])
        # print("GCTTCCGCTG")
        # print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("GCTTCCGCTG")])
        # print("TCTTCCGCTG")
        # print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("TCTTCCGCTG")])

                
        if check:
            print(check)
            check_kmers = [check[i:i+_kmerlist_I.len] for i in range(len(check) - _kmerlist_I.len + 1)]
            for i, check_i in enumerate(check_kmers):
                print(" "*i + check_i)
                print(_kmerlist_R_I.kmers["TGT"][get_kmer_index(check_i)])


################################################################################

if __name__ == "__main__":
    main()
