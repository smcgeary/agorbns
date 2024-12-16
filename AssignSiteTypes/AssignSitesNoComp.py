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

# FUNCTIONS
def LCSuff(seq1, seq2):
    if len(seq1) == 0 or len(seq2) == 0:
        return("")
    elif seq1[-1] == seq2[-1]:
        return(LCSuff(seq1[:-1], seq2[:-1]) + seq1[-1])
    else:
        return("")

def LCSubString(seq1, seq2):
    pres1 = [seq1[:i] for i in range(1, len(seq1) + 1)]
    pres2 = [seq2[:i] for i in range(1, len(seq2) + 1)]
    pre_pairs = [(i, j) for i in pres1 for j in pres2]
    LCs = [LCSuff(i[0], i[1]) for i in pre_pairs]
    LCmax = max([len(i) for i in LCs])
    return [i for i in LCs if len(i) == LCmax]

def assign_site(read_seqs, _mirna, _sitelist, experiment, n_constant,
                rand_length, min_ex_str, buffer3p, competitor_seq, c_k, m_k):
    time_start = time.time()
    counter = [0, 0, 0]
    sys.stdout.flush()
    comp_rc = get_rc(competitor_seq)
    mir_rc = get_rc(_mirna.seq)
    comp_kmers = [comp_rc[i:i + c_k]
                  for i in range(len(competitor_seq) - c_k + 1)]

    # print(_mirna.seq)
    # for i in range(8 - m_k + 1):
    #     print(" "*i + _mirna.seq[i:i + m_k])
    mirna_seqs = [get_rc(_mirna.seq[i:i + m_k])
                  for i in range(8 - m_k + 1)]
    # print(mirna_seqs)

    sys.stdout.flush()
    for i, kmer in enumerate(comp_kmers):
        print(" "* i + kmer)
        sys.stdout.flush()
        print(sum([m_k in kmer for m_k in mirna_seqs]))

    bools = [sum([m_k in kmer for m_k in mirna_seqs]) for kmer in comp_kmers]
    # print(bools)
    comp_exclude_kmers = [c_kmer for i, c_kmer in enumerate(comp_kmers)
                          if not bools[i]]
    
    for i, kmer in enumerate(comp_kmers):
        if kmer in comp_exclude_kmers:
            print(" "* i + kmer)
        else:
            print(" "* i + kmer.lower())
    sys.stdout.flush()


    k_6mers = [_mirna[site] for site in ["6mer", "6mer-m8", "6mer-A1"]]

    for i_r, r in enumerate(read_seqs):
        # Get the read number:
        read = r.strip()
        # Make Read object.
        _read = Read(read, rand_length, _mirna, n_constant, experiment)
        _read.get_all_sites_in_read(_sitelist, buffer_=buffer3p)

        sys.stdout.flush()
        # Find all sites within the read


        if buffer3p:
            seq_use = _read.seq[:n_constant + 36]
        else:
            seq_use = _read.seq
        exclude = sum([_read.seq.find(ex) for ex in min_ex_str])
        exclude_2 = sum([_read.seq.find(comp_k) for comp_k in comp_exclude_kmers])
        mirs = sum([seq_use.find(mir_k) for mir_k in k_6mers])

        # if exclude == -1*len(min_ex_str):
        #     add_sites_from_read(_sitelist, _read, buffer_=buffer3p)
        add_ = False
        if exclude == -1*len(min_ex_str) and (
            exclude_2 == -1*len(comp_exclude_kmers) or mirs != -1*len(k_6mers)
        ):
            add_ = True
            all_sites = _read.sites
            add_sites_from_read(_sitelist, _read, buffer_=buffer3p)
        if ((_read.topsite.name == "CACACACA" and add_) or 
            (_read.topsite.name == "GCTTCCGC" and add_)):
            print("                     ")
            print([i.name for i in _read.sites])
            print(exclude)
            print(exclude_2)
            print(mirs)
            print(_read)
            print([_read.seq.find(mir_k) for mir_k in k_6mers])
            print([i.l for i in _read.sites])

            sys.stdout.flush()

        #         print([i.name for i in _read.sites])
        #         starts = [_read.seq.find(mir_k) for mir_k in k_6mers]
        #         print(_read.seq)
        #         for i in starts:
        #             print(" "*i + _read.seq[i:i + 6])
        #         sys.stdout.flush()

        #         return
        #     sys.stdout.flush()

    if _sitelist.__class__.__name__ in ["TwelvemerSiteList", "SixteenmerSiteList"]:
        output = (_sitelist.top_counts_df(),
                  _sitelist.multi_counts_df())
    else:
        output = (_sitelist.top_counts_df(),
            _sitelist.multi_counts_df(),
            _sitelist.single_pos_counts_df(),
            _sitelist.top_pos_counts_df(),
            _sitelist.single_pos_counts_df(barcode="ACA"),
            _sitelist.top_pos_counts_df(barcode="ACA"))
    time_done = time.time()
    print("Finished job future in %f3.2 seconds." %(time_done - time_start))
    sys.stdout.flush()

    return (output)

def main():
    time_start = time.time()
    arguments = ["miRNA", "experiment", "condition", "n_constant", "sitelist",
                 "-mir_start", "-split16", "-n_lag", "-uniq_binary",
                 "-buffer3p_binary", "-jobs", "-test_binary", "-c_k", "-m_k"]
    args = parse_arguments(arguments)
    (mirna, experiment, condition, n_constant, sitelist, mir_start, split16,
     n_lag, uniq, buffer3p, jobs, test, c_k, m_k) = args

    competitor_seq = seq_competitor_map[mirna]

    _mirna = Mirna(mirna)

    if experiment in ["equil_mmseed_nb", "equil_seed_nb"]:
        read_length = 38
    else:
        read_length = 37

    if c_k:
        c_k = int(c_k)
    else:
        c_k = 6



    if m_k:
        m_k = int(m_k)
    else:
        m_k = 6

    m_k = min(c_k, m_k)





    if sitelist == "12mers":
        _sitelist = TwelvemerSiteList(_mirna, mir_start)
        extension = "_%s_%s_%s-%s" %(n_constant, sitelist, mir_start,
                                  int(mir_start) + 3)
    elif sitelist == "16mers":
        _sitelist = SixteenmerSiteList(_mirna, mir_start, split16)
        extension = "_%s_%s_%s-%s_%s" %(n_constant, sitelist, mir_start,
                                  int(mir_start) + 3, split16)        
    else:
        _sitelist = SiteList(_mirna, sitelist,
                             int(read_length) + 2*int(n_constant))
        extension = "_%s_%s_comp_c%sm%s" %(n_constant, sitelist, c_k, m_k)
    if test:
        extension = "%s_test" %(extension)
    # exclude_list = [i.name for i in _sitelist.sites if i]
    if n_lag:
        print(n_lag)
        n_lag = int(n_lag)
        with open("SolveForOffRates/LaggingKmers.txt") as file_in:
            kmers = [i.split("\t")[0] for i in list(file_in)[1:]]
        print(kmers)
        exclude_strings = kmers[:n_lag]
        print(exclude_strings)
        min_ex_str = MinimalExcludeStringList(exclude_strings)
        print(min_ex_str)
        extension = "%s_ex%s" %(extension, n_lag)
    else:
        min_ex_str = []


    if uniq:
        read_dir = "reads_unique"
        extension = "%s_uniq" %(extension)
    else:
        read_dir = "reads"


    if buffer3p:
        extension = "%s_buffer3p" %(extension)


    if condition == "I_combined":
        input_list = INPUT_LIST_I_COMBINED
        if not jobs:
            jobs = 20
    elif condition == "0_combined":
        input_list = INPUT_LIST_0_COMBINED
        if not jobs:
            jobs = 20
    else:
        input_list = [(mirna, experiment, condition)]
        if not jobs:
            jobs = 20



    args  = [_mirna, _sitelist, experiment, int(n_constant), read_length,
             min_ex_str, buffer3p, competitor_seq, c_k, m_k]
    print(args)
    ############################################################################
    threads = []

    # For each file in the input list, perform the multiprocessing:
    for i_input in input_list:
        reads_path = get_analysis_path(i_input[0], i_input[1], i_input[2], read_dir)
        print(reads_path)
        if not jobs:
            jobs = 1
        threads += multiproc_file(reads_path, int(jobs), assign_site, test,
                                  *args)


    # The summed data from each of the multiprocessing threads:
    equil_exps = ["equil_pilot", "equilibrium", "equil_flowthrough",
                  "equilibrium_nb", "equilibrium2_nb", "equilibrium3_nb",
                  "equilibrium_mmseed_nb"]
    kinetics_exps = ["kinetics"]
    if experiment in equil_exps:
        counts = merge_data_frames([i[0] for i in threads])["TGT"]
        multicounts = merge_data_frames([i[1] for i in threads])["TGT"]
    elif experiment in kinetics_exps:
        counts_pulse = merge_data_frames([i[0] for i in threads])["TGT"]
        counts_chase = merge_data_frames([i[0] for i in threads])["ACA"]
        counts_pulse.index = ["%s_p" %(i) for i in counts_pulse.index]
        counts_chase.index = ["%s_c" %(i) for i in counts_chase.index]
        counts = pd.concat([counts_pulse, counts_chase])
        multicounts_pulse = merge_data_frames([i[1] for i in threads])["TGT"]
        multicounts_chase = merge_data_frames([i[1] for i in threads])["ACA"]
        multicounts_pulse.index = ["%s_p" %(i) for i in multicounts_pulse.index]
        multicounts_chase.index = ["%s_c" %(i) for i in multicounts_chase.index]
        multicounts = pd.concat([multicounts_pulse, multicounts_chase])
    else:
        raise ValueError("Undefined experiment") 
    if sitelist not in ["12mers", "16mers"]:
        if experiment in equil_exps:
            single_pos_counts = merge_data_frames([i[2] for i in threads])
            top_pos_counts = merge_data_frames([i[3] for i in threads])
        elif experiment == "kinetics":
            single_pos_counts_pulse = merge_data_frames([i[2] for i in threads])
            single_pos_counts_chase = merge_data_frames([i[4] for i in threads])
            single_pos_counts_pulse.index = ["%s_p" %(i) for i in single_pos_counts_pulse.index]
            single_pos_counts_chase.index = ["%s_c" %(i) for i in single_pos_counts_chase.index]
            single_pos_counts = pd.concat([single_pos_counts_pulse, single_pos_counts_chase])
            top_pos_counts_pulse = merge_data_frames([i[3] for i in threads])
            top_pos_counts_chase = merge_data_frames([i[5] for i in threads])
            top_pos_counts_pulse.index = ["%s_p" %(i) for i in top_pos_counts_pulse.index]
            top_pos_counts_chase.index = ["%s_c" %(i) for i in top_pos_counts_chase.index]
            top_pos_counts = pd.concat([top_pos_counts_pulse, top_pos_counts_chase])

    # The names of the output files:
    site_counts_path = get_analysis_path(mirna, experiment, condition,
                                       "site_counts", ext=extension)
    multisite_counts_path = get_analysis_path(mirna, experiment, condition,
                                            "multisite_counts",
                                            ext=extension)
    single_pos_counts_path = get_analysis_path(mirna, experiment, condition,
                                            "single_pos_counts",
                                            ext=extension)
    top_pos_counts_path = get_analysis_path(mirna, experiment, condition,
                                            "top_pos_counts",
                                            ext=extension)
    outputpaths = [get_analysis_path(mirna, experiment, condition, i,
                                     ext=extension)
                   for i in ["site_counts", "multisite_counts",
                             "single_pos_counts", "top_pos_counts"]]

    for outputpath in outputpaths:
        print(outputpath)

    # Names for the positional data tables:
    if sitelist not in ["12mers", "16mers"]:
        positional_names = (["5C%s" %(i) for i in range(int(n_constant), 0, -1)] +
                            ["N%s" %(i) for i in range(37)] + 
                            ["3C%s" %(i) for i in range(1, int(n_constant) + 1)])
    print("done reading files")

    # if test:
    #     print(counts)
    #     print(" ")
    #     print(counts.sum(0))
    #     print(" ")
    #     print(single_pos_counts.iloc[:, 30:36])
    #     print(" ")
    #     print(top_pos_counts.iloc[:, 30:36])

    #     print(single_pos_counts_path)
    #     single_pos_counts_path = single_pos_counts_path.split(".")[0] + "_test.txt"
    #     top_pos_counts_path = top_pos_counts_path.split(".")[0] + "_test.txt"
    #     print(single_pos_counts_path)
    #     single_pos_counts.to_csv(single_pos_counts_path, sep="\t", header=positional_names)
    #     top_pos_counts.to_csv(top_pos_counts_path, sep="\t", header=positional_names)


    # if not test:
    print(counts)
    # print(counts.sum(0))
    # print(multisite_counts_path)
    multicounts = multicounts.astype(int)
    # print(multicounts)
    # print(multicounts.sum(0))
    print(site_counts_path)
    counts.to_csv(site_counts_path, sep="\t", header=False)
    multicounts.to_csv(multisite_counts_path, sep="\t", header=False)
    if sitelist not in ["12mers", "16mers"]:
        single_pos_counts.to_csv(single_pos_counts_path, sep="\t", header=positional_names)
        top_pos_counts.to_csv(top_pos_counts_path, sep="\t", header=positional_names)

    print_time_elapsed(time_start)
################################################################################

if __name__ == "__main__":
    main()

