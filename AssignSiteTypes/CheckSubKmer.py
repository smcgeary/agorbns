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
    arguments = ["miRNA", "condition", "kmer",
                 "-n_ex", "-alt_mirna", "-buffer3p_binary"]
    (mirna, condition,
     kmer, n_ex, alt_mirna, buffer) = parse_arguments(arguments)
    # Initialize objects:
    _mirnas = [Mirna(mirna)]
    print(_mirnas)
    n_constant = "5"
    _sitelist = SiteList(_mirnas[0], "papercutoff", 1)
    exclude_list = [i.name for i in _sitelist.sites]
    if n_ex:
        n_ex = int(n_ex)
    else:
        n_ex = len([i for i in exclude_list if i != ""])
    rand_length = 37    
    if mirna == "miR-7-23nt":
        experiment = "equilibrium2_nb"
        if condition == "4":
            condition = "12.6"
    else:
        experiment = "equilibrium"

    kmer_range = range(6, len(kmer) + 1)
    print(kmer)
    for kmer_len in kmer_range:
        sub_kmers = [kmer[i:i + kmer_len] for i in range(len(kmer) - kmer_len + 1)]
        try: 
            _kmerlist_I = KmerList(int(kmer_len), int(n_constant), int(rand_length), fast=True)
            _kmerlist_I.read(_mirnas[0], experiment, "I", n_constant, n_ex, buffer=buffer, final=True)
            _kmerlist = KmerList(int(kmer_len), int(n_constant), int(rand_length), fast=True)
            _kmerlist.read(_mirnas[0], experiment, condition, n_constant, n_ex, buffer=buffer, final=True)
            _kmerlist_R_I = _kmerlist/_kmerlist_I
            tops = sorted(_kmerlist_R_I.kmers["TGT"])[::-1][:2]
            print(tops)
            for i, sub_kmer in enumerate(sub_kmers):
                if i > 0 and i < len(sub_kmers) - 1:
                    l_ = _kmerlist_R_I.kmers["TGT"][get_kmer_index(sub_kmers[i - 1])]
                    k_ = _kmerlist_R_I.kmers["TGT"][get_kmer_index(sub_kmers[i])]
                    r_ = _kmerlist_R_I.kmers["TGT"][get_kmer_index(sub_kmers[i + 1])]
                    lin_R = k_ / ((float(l_) + float(r_))/2)
                    print("."*i + sub_kmer + "."*(len(kmer) - len(sub_kmer) - i) + " %2.2f" %(_kmerlist_R_I.kmers["TGT"][get_kmer_index(sub_kmer)]) + "\t%2.2f" %(lin_R))
                else:
                    print("."*i + sub_kmer + "."*(len(kmer) - len(sub_kmer) - i) + " %2.2f" %(_kmerlist_R_I.kmers["TGT"][get_kmer_index(sub_kmer)]) + "\tNA")

        except Exception:
            pass



    return

    # if I0:
    #     _kmerlist_I.read(_mirnas[0], experiment, "I", n_constant, 0, buffer=buffer, final=True)
    # else:
    # print("loop")

    # for _mirna in _mirnas:
    #     _kmerlist = KmerList(int(kmer_len), int(n_constant), int(rand_length), fast=True)
    #     _kmerlist.read(_mirna, experiment, condition, n_constant, n_ex, buffer=buffer, final=True)
    #     _kmerlist_R_I = _kmerlist/_kmerlist_I
    #     kmer_tups = zip(get_kmer_list(kmer_len),_kmerlist.kmers, _kmerlist_I.kmers)

    #     if alt_mirna:
    #         _mirna = Mirna(alt_mirna)
    #     print("%s:" %(_mirna.name))

    #     _kmerlist_R_I.head_with_mirna_sites(_mirna, n=n, sub=sub)
    #     _kmerlist_R_I.top_and_adjacent()
    #     _kmerlist_R_I.top_substitute()
    #     print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("GCTTCCGCTA")])
    #     print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("GCTTCCGCTC")])
    #     print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("GCTTCCGCTG")])
    #     print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("GCTTCCGCTT")])

        # print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("AACATTCG")])
        # print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("CACATTCG")])
        # print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("GACATTCG")])
        # print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("TACATTCG")])

        # print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("AACATTCT")])
        # print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("CACATTCT")])
        # print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("GACATTCT")])
        # print(_kmerlist_R_I.kmers["TGT"][get_kmer_index("TACATTCT")])


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

                
        # if check:
        #     print(check)
        #     check_kmers = [check[i:i+_kmerlist_I.len] for i in range(len(check) - _kmerlist_I.len + 1)]
        #     for i, check_i in enumerate(check_kmers):
        #         print(" "*i + check_i)
        #         print(_kmerlist_R_I.kmers["TGT"][get_kmer_index(check_i)])


################################################################################

if __name__ == "__main__":
    main()
