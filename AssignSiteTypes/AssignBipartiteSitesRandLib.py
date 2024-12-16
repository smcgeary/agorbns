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
from collections import defaultdict
np.core.arrayprint._line_width = 200

## GLOBAL CONSTANTS ############################################################

# FUNCTIONS
def assign_site(read_seqs, _mirna, experiment, n_con, n_lib, buffer3p):
    time_start = time.time()
    ## SET UP THE DICTIONARIES REQUIRED FOR CONVERTING THE THREEPRIME ##########
    ## SEQUENCES TO THEIR NAMES, FOR PUTTING THEM INTO THE OUTPUT ##############
    ## DICTIONARIES. ###########################################################
    # Make the list of all possible programmed sites in each read.
    seq_8mer = _mirna["8mer"]
    _sitelistseed = SiteList(_mirna, "programmedbase", n_lib + 2*n_con)
    _sitelistthrp = SiteList(_mirna, "threeprime", n_lib + 2*n_con)
    # print(_sitelistthrp)
    # print([site.seq for site in _sitelistthrp.sites])
    # print([site.name for site in _sitelistthrp.sites])
    # dict_sites = defaultdict(int)
    # for site_seq in [site.seq for site in _sitelistthrp.sites]:
    #     dict_sites[site_seq] += 1
    # print(dict_sites)
    # return
    # 1. Make the seed site 3′ nucleotide dictionaries:
    seed_end_positions = [8, 8, 7, 7, 6, 8] + [8]*18
    endPos_seedName_map = {name : pos for (name, pos)
                           in zip([site.name for site in _sitelistseed.sites],
                                  seed_end_positions)}

    counts = defaultdict(int)
    for i_r, r in enumerate(read_seqs):
        # check_read = False
        # print_line = False
        read = r.strip()
        _read = Read(read, n_lib, _mirna, n_con, experiment)
        add_sites_from_read(_sitelistseed, _read, buffer_=buffer3p)
        seed_sites = _read.sites
        seed_names = [site.name for site in seed_sites]
        add_sites_from_read(_sitelistthrp, _read, buffer_=buffer3p)
        thrp_sites = _read.sites
        thrp_names = [site.name for site in thrp_sites]
        thrp_l = [site.l for site in thrp_sites]
        thrp_r = [site.r for site in thrp_sites]
        seed_l = [site.l for site in seed_sites]
        seed_r = [site.r for site in seed_sites]
        match = re.search("ACAACTACCTCA", _read.seq)
        bipartite_sites = []
        if (seed_names != ["None"] and thrp_names != ["None"]):
            # First exclude any threep sites that overlap any of the seed_sites.
            thrp_sites_final = []
            for thrp_site in thrp_sites:
                # Define the ends of the thrp_site.
                add_seq = True
                t_l, t_r  = (thrp_site.l, thrp_site.r)
                for seed_site in seed_sites:
                    # Define the ends of the seed_site.
                    s_l, s_r = (seed_site.l, seed_site.r)
                    # Check if the thrp site overlaps with the seed site.
                    s_l_int = s_l >= t_l and s_l < t_r
                    s_r_int = s_r > t_l and s_r <= t_r
                    t_l_int = t_l >= s_l and t_l < s_r
                    t_r_int = t_r > s_l and t_r <= s_r
                    if s_l_int or s_r_int or t_l_int or t_r_int:
                        # This is overlap in which the three-prime site is 5-p of the
                        # seed site, such that the three-prime site could potentially
                        # be clipped.
                        if s_l_int and t_r_int:
                            len_thrp = int(thrp_site.name.split("mer")[0])
                            thrp_limits = [int(i) for i in thrp_site.name.split("mer-m")[1].split(".")]
                            overlap = t_r - s_l
                            len_thrp_new = len_thrp - overlap
                            thrp_limits_new = list(thrp_limits)
                            thrp_limits_new[0] = thrp_limits_new[0] + overlap
                            if len_thrp_new >= 4:
                                thrp_name = "%smer-m%s.%s" %(len_thrp_new, thrp_limits_new[0], thrp_limits_new[1])
                                thrp_site = Read.ReadSite(_sitelistthrp, thrp_name, t_l, s_l)
                            else:
                                add_seq = False
                        else:   
                            add_seq = False
                if add_seq:
                    thrp_sites_final += [thrp_site]
            thrp_sites = thrp_sites_final
            if len(thrp_sites) != 0:
                # makes sure that only the maximum length three-prime site is used.
                thrp_sites_max_length = max([int(thrp_site.name.split("mer")[0]) for thrp_site in thrp_sites])
                thrp_sites = [thrp_site for thrp_site in thrp_sites if int(thrp_site.name.split("mer")[0]) == thrp_sites_max_length]
                for seed_site in seed_sites:
                    s_l, s_r = (seed_site.l, seed_site.r)
                    t_r = [thrp_site.r for thrp_site in thrp_sites]
                    thrp_sites_temp = [thrp_site for thrp_site in thrp_sites if thrp_site.r <= s_l]
                    seed_pos_end = endPos_seedName_map[seed_site.name]
                    for thrp_site in thrp_sites_temp:
                        # Calculate the actual distance between the seed site and
                        # thrp site.
                        t_l, t_r = (thrp_site.l, thrp_site.r)
                        pos = seed_site.l - thrp_site.r
                        # Calculate the adjustment for where the seed site ends
                        # (i.e., The 7mer-A1 ends at position 7, not 8), and add
                        # 1 to give the position at which the 3′ site begins.
                        pos += endPos_seedName_map[seed_site.name] + 1
                        full_site_name = "%s|%s|%s" %(thrp_site.name, pos, seed_site.name)
                        bipartite_sites += [full_site_name]
                # Condition in which at least one bipartite site has bene found:
                if len(bipartite_sites) != 0:
                    # print(" ".join(bipartite_sites))
                    # sys.stdout.flush()
                    # Add all of the bipartite sites.
                    for bipartite_site in bipartite_sites:
                        counts[bipartite_site] += 1./len(bipartite_sites)
                # Otherwise, add each of the seed and thrp sites to the count
                # dictionary, dividing the read by the total number of sites.
                else:
                    if print_read:
                            print("added separately, split counts in %s" %(len(seed_sites) + len(thrp_sites)))
                            sys.stdout.flush()
                    for seed_site in seed_sites:
                        counts[seed_site.name] += 1./(len(seed_sites) + len(thrp_sites))
                    for thrp_site in thrp_sites:
                        counts[thrp_site.name] += 1./(len(seed_sites) + len(thrp_sites))
            # The condition in which none of the thrp sites survive the overlap
            # with seed site criteria, in which case there will only be seed
            # sites left. So, add the seed sites to the dictionary.
            else:
                for seed_site in seed_sites:
                    counts[seed_site.name] += 1./len(seed_sites)
        # Conditions in which wasn't potential seed and thrp sites to begin
        # with, such that the bipartite site analysis portion of the loop was
        # never entered. So, that means there are only either seed or thrp sites
        # within the read, so check for each "None" condition and add whichever
        # one is correct.
        elif seed_names != ["None"]:
            for seed_site in seed_sites:
                counts[seed_site.name] += 1./len(seed_sites)
        elif thrp_names != ["None"]:
            for thrp_site in thrp_sites:
                counts[thrp_site.name] += 1./len(thrp_sites)
        else:
            counts["None"] += 1
    print("i_r: %s" %i_r)
    print("summed dict: %s" %sum(counts.values()))
    return(counts)

def main():
    time_start = time.time()
    arguments = [
        "miRNA", "experiment", "condition", "n_constant", "-jobs",
        "-buffer3p_binary", "-test_binary"
    ]
    args = parse_arguments(arguments)
    (
        mirna, experiment, condition, n_constant, jobs, buffer3p, test
    ) = args
    # Determine the length of the "random" portion.
    # `Random` is in quotes because in the programmed libraries positions 26-31
    # are not actually random.
    # if ("let-7a" in mirna or "miR-1" == mirna) and mirna != "miR-155_let-7a":
    #     rand_length = 38
    # else:
    rand_length = 37
    _mirna = Mirna(mirna)
    # Assign three different extensions for the three different count files
    # being made.
    extension = "_%s_randthrp" %n_constant
    # extension_collapsed = "%s_collapsed" %extension
    extension_comp = "%s_comp" %extension
    extension_suppcomp = "%s_suppcomp" %extension

    if buffer3p:
        extension = "%s_buffer3p" %extension
        # extension_collapsed = "%s_buffer3p" %extension_collapsed
        extension_comp = "%s_buffer3p" %extension_comp
        extension_suppcomp = "%s_buffer3p" %extension_suppcomp

    if test:
        extension = "%s_test" %extension
        # extension_collapsed = "%s_test" %extension_collapsed
        extension_comp = "%s_test" %extension_comp
        extension_suppcomp = "%s_test" %extension_suppcomp

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
    elif condition == "0_combined":
        input_list = INPUT_LIST_0_COMBINED
        if not jobs:
            jobs = 20
    else:
        input_list = [(mirna, experiment, condition)]
        if not jobs:
            jobs = 20

    args  = [_mirna, experiment, int(n_constant), rand_length, buffer3p]
    print(args)
    ############################################################################
    threads = []
    # For each file in the input list, perform the multiprocessing:
    for i_input in input_list:
        reads_path = get_analysis_path(
            i_input[0], i_input[1], i_input[2], "reads"
        )
        if not jobs:
            jobs = 1
        threads += multiproc_file(reads_path, int(jobs), assign_site, test, *args)
    # Initialize an output dictionary for which all of the value across the
    # threads will be combined.
    output_merged = defaultdict(int)
    for thread in threads:
        for (key, value) in thread.items():
            output_merged[key] += value
    # MAKE THE OUTPUT DATAFRAME, WITH NOTHING COLLAPSED
    counts = pd.DataFrame.from_dict(output_merged, orient="index")
    # Make the dataframe in which the supplementary and compensatory sites are
    # averaged together, along with the shifted seed sites.
    output_comp = defaultdict(int)
    output_suppcomp = defaultdict(int)
    # Iterate over each row in the output dataframe, checking if it is a
    # bipartite site (due to the presence of "|X|" where X is an integer), and
    # if it is, split up the site and change the seed site to either Comp for
    # compensatory site, Supp for supplementary site, or Offset for an offset
    # 6mer site (being the 6mer-m8 or the 6mer-A1).
    for name in counts.index.values:
        if "|" in name:
            thrp_site, num, base_site = name.split("|")
            if base_site in ["8mer", "7mer-m8", "7mer-A1", "6mer"]:
                new_base = "Supp"
                new_comp_base = base_site
            elif base_site in ["6mer-m8", "6mer-A1"]:
                new_base = "Offset"
                new_comp_base = base_site
            else:
                new_base = "Comp"
                new_comp_base = new_base
            new_name = "|".join([thrp_site, num, new_base])
            new_comp_name = "|".join([thrp_site, num, new_comp_base])
        else:
            new_name = name
            new_comp_name = name
        output_suppcomp[new_name] += counts.loc[name]
        output_comp[new_comp_name] += counts.loc[name]
    # Make a pandas dataframe from the suppcomp dictionary.
    counts_suppcomp = pd.DataFrame.from_dict(output_suppcomp, orient="index")
    counts_comp = pd.DataFrame.from_dict(output_comp, orient="index")
    # DEFINE THE PATHS TO OUTPUT THE FILES.
    site_counts_path = get_analysis_path(
        mirna, experiment, condition, "site_counts", ext=extension
    )
    site_counts_suppcomp_path = get_analysis_path(
        mirna, experiment, condition, "site_counts",
        ext=extension_suppcomp
    )
    site_counts_comp_path = get_analysis_path(
        mirna, experiment, condition, "site_counts",
        ext=extension_comp
    )    
    print(site_counts_path)
    # print(site_counts_collapsed_path)
    print(site_counts_suppcomp_path)
    print(site_counts_comp_path)
    print(counts.sum())
    print(counts_comp.sum())
    print(counts_suppcomp.sum())
    # Output the dataframes to the three paths.
    counts.to_csv(site_counts_path, sep="\t", header=False)
    # counts_collapsed.to_csv(site_counts_collapsed_path, sep="\t", header=False)
    counts_comp.to_csv(site_counts_comp_path, sep="\t", header=False)
    counts_suppcomp.to_csv(site_counts_suppcomp_path, sep="\t", header=False)
    print_time_elapsed(time_start)
################################################################################

if __name__ == "__main__":
    main()

