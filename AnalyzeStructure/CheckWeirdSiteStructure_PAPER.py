################################################################################
#MakeSiteTypeReadFiles.py
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
import RNA

def randomword(length):
   return ''.join(random.choice(string.lowercase) for i in range(length))

# FUNCTIONS

AU_weights_left = [1.0 / i for i in range(30,0,-1)]
AU_weights_right = [1.0 / i for i in range(1,31)]
# AU_weights_right[0] = 0.5


WC = {
    "A" : "T",
    "T" : "A", 
    "C" : "G",
    "G" : "C"
}

wob_WC = {
    "A" : ["T"],
    "T" : ["A", "G"], 
    "C" : ["G"],
    "G" : ["C", "T"]
}


def check_GCTTCCG(read_seqs, _mirna, _sitelist, experiment,
                  n_constant, rand_length, buffer3p,
                  num_i):
    """
    """

    out_lengths = []
    out_wob_lengths = []
    out_percentages = []
    out_left_stem = []
    out_right_stem = []

    for i_r, r in enumerate(read_seqs):
        print("_____________________________")
        # Get the read number:
        read = r.strip()
        # Make Read object.
        _read = Read(read, rand_length, _mirna, n_constant, experiment)
        # Find all sites within the read
        add_sites_from_read(_sitelist, _read, buffer_=buffer3p)
        inds = [i for i, site in enumerate(_read.sites) if site.name == "GCTTCCGC"]
        for ind in inds:
            l = _read.sites[ind].l
            r = _read.sites[ind].r
            wc_extend = True
            extend = 1
            l_new = l
            r_new = r
            while wc_extend and l_new >= 0 and r_new <= len(_read.seq):
                l_new = l - extend
                r_new = r + extend - 1
                l_nuc = _read.fullseq[l_new + 21]
                r_nuc = _read.fullseq[r_new + 21]
                # print(" "*(l_new + 21) + l_nuc + "-"*(r_new - l_new - 1) + r_nuc + "(%s)" %(WC[l_nuc]))
                if WC[l_nuc] != r_nuc:
                    wc_extend = False
                else:
                    extend += 1
            extend_final = extend - 1
            print("extend_final")
            print(extend_final)

            wob_wc_extend = True
            wob_extend = 1
            l_new = l
            r_new = r
            while wob_wc_extend and l_new >= 0 and r_new <= len(_read.seq):
                l_new = l - wob_extend
                r_new = r + wob_extend - 1
                l_nuc = _read.fullseq[l_new + 21]
                r_nuc = _read.fullseq[r_new + 21]
                # print(" "*(l_new + 21) + l_nuc + "-"*(r_new - l_new - 1) + r_nuc + "(%s)" %(WC[l_nuc]))
                if r_nuc not in wob_WC[l_nuc]:
                    wob_wc_extend = False
                else:
                    wob_extend += 1
            wob_extend_final = wob_extend - 1
            print("wob_extend_final")
            print(wob_extend_final)
            if wob_extend_final < extend_final:
                print("this doesn't make sense.")
                return

            str_pf_fold, pf_dG = RNA.pf_fold(_read.fullseq)
            site_str_pf_fold = str_pf_fold[(l + 21):(r + 21)]
            left_pf_fold = str_pf_fold[(l + 18):(l + 23)]
            right_pf_fold = str_pf_fold[(r + 19):(r + 19 + 5)]

            string_TTCC_open = sum([i == "." for i in site_str_pf_fold[2:6]])/float(4)
            string_TTCC_stem_l = sum([i in ["(", "{"] for i in left_pf_fold])/float(5)
            string_TTCC_stem_r = sum([i in [")", "}"] for i in right_pf_fold])/float(5)
            sys.stdout.flush()

            out_lengths.append(extend_final)
            out_wob_lengths.append(wob_extend_final)
            out_percentages.append(string_TTCC_open)
            out_left_stem.append(string_TTCC_stem_l)
            out_right_stem.append(string_TTCC_stem_r)


            # print("%s %s\t%s\t%s\t%s\t%s" %(num_i, i_r, extend_final,
            #                                 string_TTCC_open,
            #                                 string_TTCC_stem_l,
            #                                 string_TTCC_stem_r))
            sys.stdout.flush()

    print("done job %s" %(num_i))
    sys.stdout.flush()
    return out_lengths, out_wob_lengths, out_percentages, out_left_stem, out_right_stem



def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["n_constant", "-buffer3p_binary", "-jobs", "-test_binary"]
    (
        n_constant, buffer3p, jobs, test
    ) = parse_arguments(arguments)

    print(buffer3p)
    mirna = "miR-1"
    experiment = "equilibrium"
    n_constant = 5
    sitelist = "paperfinal"
    _mirna = Mirna(mirna)
    if experiment in ["equil_mmseed_nb", "equil_seed_nb"]:
        read_length = 38
    else:
        read_length = 37
    if jobs:
        jobs = int(jobs)
    else:
        jobs = 2

    _mirna = Mirna(mirna)

    conditions = ["I", "0.4", "1.26", "4", "12.6", "40"]
    out_averages = [None, None, None, None, None, None]
    out_stds = [None, None, None, None, None, None]
    _sitelist = SiteList(_mirna, "paperfinal",
                     int(read_length) + 2*int(n_constant))

    cwd = os.getcwd()
    folder_name = randomword(20)
    plfold_path = "%s/AnalyzeStructure/RNAplfold_temp/%s" %(HOME_DIR,
                                                                folder_name)
    os.makedirs(plfold_path)
    os.chdir(plfold_path)

    for i_c, condition in enumerate(conditions):
        print(condition)
        i_input = (mirna, experiment, condition)
        args  = [_mirna, _sitelist, experiment, int(n_constant),
                 read_length, buffer3p]

    # Move to the temporary folder for making the pl fold files:

        reads_path = get_analysis_path(i_input[0], i_input[1], i_input[2],
                                       "GCTTCCGC_reads")
        threads = multiproc_file(reads_path, int(jobs), check_GCTTCCG,
                                    test, *args)

        averages = [np.mean([k for i in threads for k in i[j]]) for j in range(5)]
        stds = [np.std([k for i in threads for k in i[j]]) for j in range(5)]
        ns = [len([k for i in threads for k in i[j]]) for j in range(5)]
        print(stds)
        print(ns)
        ses = [i[0]/float(math.sqrt(i[1])) for i in zip(stds, ns)]

        out_averages[i_c] = averages
        out_stds[i_c] = ses

    os.chdir(cwd)
    os.rmdir(plfold_path)

    print(out_averages)
    print(out_stds)

    out_averages_df = pd.DataFrame(out_averages, index=conditions,
                                   columns=["hairpin_nt", "hairpin_wob_nt",
                                            "% TTCC open", "len left GC",
                                            "len right GC"])
    out_stds_df = pd.DataFrame(out_stds, index=conditions,
                               columns=["hairpin_nt", "hairpin_wob_nt",
                                        "% TTCC open", "len left GC",
                                        "len right GC"])

    print(out_averages_df)
    print(out_stds_df)

    print_time_elapsed(time_start)

################################################################################

if __name__ == "__main__":
    main()

