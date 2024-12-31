################################################################################
# AssignSites.py
################################################################################
import imp # Used to import general.py
imp.load_source("general", "general/general.py")
imp.load_source("RBNS_methods", "general/RBNS_methods.py")
from general import *
from RBNS_methods import *

# FUNCTIONS
def assign_site(read_seqs, _mirna, _sitelist, experiment,
                             n_constant, rand_length, min_ex_str, buffer3p):
    time_start = time.time()
    for i_r, r in enumerate(read_seqs):
        # Get the read number:
        read = r.strip()
        # Make Read object.
        _read = Read(read, rand_length, _mirna, n_constant, experiment)
        # Find all sites within the read
        exclude = sum([_read.seq.find(ex) for ex in min_ex_str])
        if exclude == -1*len(min_ex_str):
            add_sites_from_read(_sitelist, _read, buffer_=buffer3p)
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
    print(("Finished job future in %f3.2 seconds." %(time_done - time_start)))
    sys.stdout.flush()

    return (output)

def main():
    time_start = time.time()
    arguments = ["miRNA", "experiment", "condition", "n_constant", "sitelist",
                 "-mir_start", "-split16", "-n_lag", "-alt_mirna", 
                 "-uniq_binary", "-buffer3p_binary",
                 "-jobs", "-test_binary"]
    ## Consider getting rid of n_lag argument.
    args = parse_arguments(arguments)
    (mirna, experiment, condition, n_constant, sitelist, mir_start, split16,
     n_lag, alt_mirna, uniq, buffer3p, jobs, test) = args


    if alt_mirna:
        _mirna = Mirna(alt_mirna)
    else:
        _mirna = Mirna(mirna)

    if (("let-7a" in mirna or "miR-1" in mirna) and
        mirna != "miR-155_let-7a" and
        experiment in ["equil_c_nb", "equil_s_nb", "equil_c_alt_nb",
                       "equil_sc_nb", "equil_sc_alt_nb", "equil_c2_nb",
                       "equil_c2_alt_nb"]):
        read_length = 38
    else:
        read_length = 37

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
        extension = "_%s_%s" %(n_constant, sitelist)
    if alt_mirna:
        extension= "%s_alt-%s" %(extension, alt_mirna)
    if test:
        extension = "%s_test" %(extension)
    print([i.name for i in _sitelist.sites])
    # if n_lag:
    #     # print(n_lag)
    #     n_lag = int(n_lag)
    #     with open("SolveForOffRates/LaggingKmers.txt") as file_in:
    #         kmers = [i.split("\t")[0] for i in list(file_in)[1:]]
    #     print(kmers)
    #     exclude_strings = kmers[:n_lag]
    #     print(exclude_strings)
    #     min_ex_str = MinimalExcludeStringList(exclude_strings)
    #     print(min_ex_str)
    #     extension = "%s_ex%s" %(extension, n_lag)
    # else:
    min_ex_str = [] # This can all go away

    if uniq:
        read_dir = "reads_unique"
        extension = "%s_uniq" %(extension)
    else:
        read_dir = "reads"

    if buffer3p:
        extension = "%s_buffer3p" %(extension)

    if condition == "I_combined":
        if "tp" in experiment:
            if mirna == "miR-7-24nt":
                input_list = [(mirna, experiment, condition)]
            else:
                input_list = INPUT_LIST_I_COMBINED_TP
        elif "miR-7" in mirna and experiment in ["equilibrium2_nb", "equilibrium3_nb"]:
            input_list = INPUT_LIST_I_COMBINED_NB
        else:
            input_list = INPUT_LIST_I_COMBINED
        if not jobs:
            jobs = 20
    elif condition == "0_combined":                 # Presumably this can go away
        input_list = INPUT_LIST_0_COMBINED          # Presumably this can go away
        if not jobs:
            jobs = 20
    else:
        input_list = [(mirna, experiment, condition)]
        if not jobs:
            jobs = 20
    args  = [_mirna, _sitelist, experiment, int(n_constant), read_length,
             min_ex_str, buffer3p]
    print(args)
    ############################################################################
    threads = []

    # For each file in the input list, perform the multiprocessing:
    for i_input in input_list:
        reads_path = get_analysis_path(i_input[0], i_input[1], i_input[2], read_dir)
        if not jobs:
            jobs = 1
        threads += multiproc_file(reads_path, int(jobs), assign_site, test,
                                  *args)
    # The summed data from each of the multiprocessing threads:
    # equil_exps = ["equil_pilot", "equilibrium", "equil_flowthrough",
    #               "equilibrium_nb", "equilibrium2_nb", "equilibrium3_nb",
    #               "equil_c_nb", "equil_s_nb", "equil_c2_nb", "equil_c2_alt_nb",
    #               "equil_c_alt_nb", "equil_sc_alt_nb", "equil_sc_nb",
    #               "equilibrium_tp", "equilibrium_met_tp", "equilibrium_2_tp"]
    # kinetics_exps = ["kinetics"]
    # if experiment in equil_exps:
    counts = merge_data_frames([i[0] for i in threads])["TGT"]
    multicounts = merge_data_frames([i[1] for i in threads])["TGT"]
    # elif experiment in kinetics_exps:
    #     counts_pulse = merge_data_frames([i[0] for i in threads])["TGT"]
    #     counts_chase = merge_data_frames([i[0] for i in threads])["ACA"]
    #     counts_pulse.index = ["%s_p" %(i) for i in counts_pulse.index]
    #     counts_chase.index = ["%s_c" %(i) for i in counts_chase.index]
    #     counts = pd.concat([counts_pulse, counts_chase])
    #     multicounts_pulse = merge_data_frames([i[1] for i in threads])["TGT"]
    #     multicounts_chase = merge_data_frames([i[1] for i in threads])["ACA"]
    #     multicounts_pulse.index = ["%s_p" %(i) for i in multicounts_pulse.index]
    #     multicounts_chase.index = ["%s_c" %(i) for i in multicounts_chase.index]
    #     multicounts = pd.concat([multicounts_pulse, multicounts_chase])
    # else:
    #     raise ValueError("Undefined experiment") 
    if sitelist not in ["12mers", "16mers"]:
        # if experiment in equil_exps:
        single_pos_counts = merge_data_frames([i[2] for i in threads])
        top_pos_counts = merge_data_frames([i[3] for i in threads])
        # elif experiment == "kinetics":
        #     single_pos_counts_pulse = merge_data_frames([i[2] for i in threads])
        #     single_pos_counts_chase = merge_data_frames([i[4] for i in threads])
        #     single_pos_counts_pulse.index = ["%s_p" %(i) for i in single_pos_counts_pulse.index]
        #     single_pos_counts_chase.index = ["%s_c" %(i) for i in single_pos_counts_chase.index]
        #     single_pos_counts = pd.concat([single_pos_counts_pulse, single_pos_counts_chase])
        #     top_pos_counts_pulse = merge_data_frames([i[3] for i in threads])
        #     top_pos_counts_chase = merge_data_frames([i[5] for i in threads])
        #     top_pos_counts_pulse.index = ["%s_p" %(i) for i in top_pos_counts_pulse.index]
        #     top_pos_counts_chase.index = ["%s_c" %(i) for i in top_pos_counts_chase.index]
        #     top_pos_counts = pd.concat([top_pos_counts_pulse, top_pos_counts_chase])

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
    multicounts = multicounts.astype(int)
    counts.to_csv(site_counts_path, sep="\t", header=False)
    multicounts.to_csv(multisite_counts_path, sep="\t", header=False)
    if sitelist not in ["12mers", "16mers"]:
        single_pos_counts.to_csv(single_pos_counts_path, sep="\t", header=positional_names)
        top_pos_counts.to_csv(top_pos_counts_path, sep="\t", header=positional_names)

    print_time_elapsed(time_start)
################################################################################

if __name__ == "__main__":
    main()

