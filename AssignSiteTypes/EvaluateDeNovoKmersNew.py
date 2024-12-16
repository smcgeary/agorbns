################################################################################
#EvaluateDeNovoKmersNew.py
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
    arguments = ["miRNA", "experiment", "condition", "kmer_len", "-n_ex", 
                 "-n", "-kmer_check", "-uniq_binary", "-buffer3p_binary",
                 "-sub_binary", "-I0_binary", "-final_binary",
                 "-combined_binary", "-cutoffdil_binary"]
    (mirna, experiment, condition, kmer_len, n_ex,
     n, kmer_check, uniq, buffer, sub,
     I0, final, combined, cutoffdil     ) = parse_arguments(arguments)
    # Make the miRNA object.
    _mirna = Mirna(mirna)
    n_constant = "5"
    if n:
        n = int(n)
    else:
        n = 10
    if "tp" in experiment:
        site_string = "papercutoff2"
    else:
        site_string = "papercutoff"
    rand_length = 37    

    _sitelist = SiteList(_mirna, site_string, 1)
    exclude_list = [i.name for i in _sitelist.sites]
    if n_ex:
        n_ex = int(n_ex)
        single_n_ex = True
    else:
        n_ex = len([i for i in exclude_list if i != ""])
        single_n_ex = False

    if mirna == "miR-1" and experiment == "equilibrium":
        buffer = True
    else:
        buffer = False

    start = True
    print(n_ex)
    if single_n_ex:
        # This just prints out a single n_ex.
        _kmerlist_I = KmerList(int(kmer_len), int(n_constant), int(rand_length),
                               fast=True)
        _kmerlist_I.read(_mirna, experiment, "I", n_constant, n_ex, uniq=uniq,
                         buffer=buffer, final=final, cutoffdil=cutoffdil,
                         resub=True)
        _kmerlist = KmerList(int(kmer_len), int(n_constant), int(rand_length),
                             fast=True)
        _kmerlist.read(_mirna, experiment, condition, n_constant, n_ex,
                       uniq=uniq, buffer=buffer, final=final,
                       cutoffdil=cutoffdil, resub=True)
        _kmerlist_R_I = _kmerlist/_kmerlist_I

        if kmer_check:
            if len(kmer_check) >= int(kmer_len):
                sub_kmers = [kmer_check[i : i + int(kmer_len)]
                             for i in range(len(kmer_check) - int(kmer_len) + 1)]
                print_string = ""
                for i, sub_kmer in enumerate(sub_kmers):
                    added_string = sub_kmer + ": "
                    print_string += added_string
                    added_thing = "%6.3f" %_kmerlist_R_I.kmers["TGT"][get_kmer_index(sub_kmer)]
                    print_string += added_thing
                    print_string += ", "
                print_string = print_string[:-2]
                print(print_string)
            elif len(kmer_check) == int(kmer_len) - 1:
                sub_kmers_l = [nuc + kmer_check for nuc in ["A", "C", "G", "T"]]
                sub_kmers_r = [kmer_check + nuc for nuc in ["A", "C", "G", "T"]]
                sub_kmers = sub_kmers_l + sub_kmers_r
                print(sub_kmers)
                R_values = [_kmerlist_R_I.kmers["TGT"][get_kmer_index(sub_kmer)]
                            for sub_kmer in sub_kmers]
                print(R_values)
                ind_max = [i for i in range(len(sub_kmers))
                           if R_values[i] == max(R_values)][0]
                R_values_no_max = [R_value for R_value in R_values
                                   if R_value != max(R_values)]

                ind_2nd_max = [i for i in range(len(sub_kmers))
                               if R_values[i] == max(R_values_no_max)][0]
                print(ind_max)
                kmer_sort = [sub_kmers[ind_max], R_values[ind_max], 0, 0, 0, 0]
                kmer_sort_2 = [sub_kmers[ind_2nd_max], R_values[ind_2nd_max], 0, 0, 0, 0]

                _mirna.get_best_name(kmer_sort, suboptimal=False) 
                _mirna.get_best_name(kmer_sort_2, suboptimal=False) 
        else:
            _kmerlist_R_I.head_with_mirna_sites(_mirna, n=n, sub=sub)
            _kmerlist_R_I.top_and_adjacent()
            _kmerlist_R_I.top_substitute()
    else:
        for n_ex_i in range(int(n_ex) + 1):
            _kmerlist_I = KmerList(int(kmer_len), int(n_constant), int(rand_length),
                                   fast=True)
            _kmerlist_I.read(_mirna, experiment, "I", n_constant, n_ex_i,
                                 uniq=uniq, buffer=buffer, final=final,
                                 cutoffdil=cutoffdil, resub=True)

            _kmerlist = KmerList(int(kmer_len), int(n_constant),
                                 int(rand_length), fast=True)
            _kmerlist.read(_mirna, experiment, condition, n_constant, n_ex_i,
                           uniq=uniq, buffer=buffer, final=final,
                           cutoffdil=cutoffdil, resub=True)
            _kmerlist_R_I = _kmerlist/_kmerlist_I
            _kmerlist_R_I.output_for_table(_mirna, n_ex_i+1, n=n, sub=sub,
                                           start=start)
            start = False

################################################################################

if __name__ == "__main__":
    main()
