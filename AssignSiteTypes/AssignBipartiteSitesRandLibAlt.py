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
    names_can = ["8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-m8", "6mer-A1"]
    _sitelistseed = SiteList(_mirna, "programmedbase", n_lib + 2*n_con)
    _sitelistthrp = SiteList(_mirna, "threeprime", n_lib + 2*n_con)
    print([site.name for site in _sitelistseed.sites])
    print([site.name for site in _sitelistthrp.sites])
    # 1. Make the seed site 3′ nucleotide dictionaries:
    seed_end_positions = [8, 8, 7, 7, 6, 8] + [8]*18
    endPos_seedName_map = {name : pos for (name, pos)
                           in zip([site.name for site in _sitelistseed.sites],
                                  seed_end_positions)}

    # Initialize the site count dictionary that will allow the kds to be fit.
    counts = defaultdict(int)

    for i_r, r in enumerate(read_seqs):
        # check_read = False
        # print_line = False
        read = r.strip()
        _read = Read(read, n_lib, _mirna, n_con, experiment)
        add_sites_from_read(_sitelistseed, _read, buffer_=buffer3p)
        seed_sites = _read.sites
        # seed_names = [site.name for site in seed_sites]
        add_sites_from_read(_sitelistthrp, _read, buffer_=buffer3p)
        thrp_sites = _read.sites
        # thrp_names = [site.name for site in thrp_sites]
        # print("_____________________________")
        # print(thrp_sites)
        # print(seed_sites)
        # print(thrp_sites[0])
        # print(seed_sites[0])
        # print(thrp_sites[0].name == "None")
        # print(seed_sites[0].name == "None")
        sys.stdout.flush()
        bipartite_sites = []
        if (seed_sites[0].name != "None" and thrp_sites[0].name != "None"):
            # First exclude any threep sites that overlap any of the seed_sites.
            thrp_sites_final = []
            for tp in thrp_sites:
                # Define the ends of the thrp_site.
                add_thrp = True
                for s in seed_sites:
                    # Check if the thrp site overlaps with the seed site.
                    if (s.r >= tp.l and s.r <= tp.r) or (s.l >= tp.l and s.l <= tp.r):
                        add_thrp = False
                # 6.iii.h If it passes both tests, add to threeprime list.
                if add_thrp:
                    thrp_sites_final += [tp]
            thrp_sites = thrp_sites_final
            if len(thrp_sites) != 0:
                thrp_sites_pre = thrp_sites
                for s in seed_sites:
                    # s_l, s_r = (seed_site.l, seed_site.r)
                    thrp_sites = [ts for ts in thrp_sites_pre if ts.r <= s.l]
                    seed_pos_end = endPos_seedName_map[s.name]
                    for ts in thrp_sites:
                        # Calculate the actual distance between the seed site and
                        # thrp site.
                        # t_l, t_r = (thrp_site.l, thrp_site.r)
                        pos = s.l - ts.r + endPos_seedName_map[s.name] + 1
                        # Calculate the adjustment for where the seed site ends
                        # (i.e., The 7mer-A1 ends at position 7, not 8), and add
                        # 1 to give the position at which the 3′ site begins.
                        # pos += endPos_seedName_map[s.name] + 1
                        # full_site_name = "%s|%s|%s" %(ts.name, pos, s.name)
                        bipartite_sites += ["%s|%s|%s" %(ts.name, pos, s.name)]
                # Condition in which at least one bipartite site has bene found:
                if len(bipartite_sites) != 0:
                    v = 1./len(bipartite_sites)
                    # Add all of the bipartite sites.
                    for s in bipartite_sites:
                        counts[s] += v
                # Otherwise, add each of the seed and thrp sites to the count
                # dictionary, dividing the read by the total number of sites.
                else:
                    v = 1./(len(seed_sites) + len(thrp_sites))
                    for s in seed_sites + thrp_sites:
                        counts[s.name] += v
            # The condition in which none of the thrp sites survive the overlap
            # with seed site criteria, in which case there will only be seed
            # sites left. So, add the seed sites to the dictionary.
            else:
                for s in seed_sites:
                    counts[s.name] += 1./len(seed_sites)
        # Conditions in which wasn't potential seed and thrp sites to begin
        # with, such that the bipartite site analysis portion of the loop was
        # never entered. So, that means there are only either seed or thrp sites
        # within the read, so check for each "None" condition and add whichever
        # one is correct.
        elif seed_sites[0].name != "None":
            v = 1./len(seed_sites)
            for s in seed_sites:
                counts[s.name] += v
        elif thrp_sites[0].name != "None":
            v = 1./len(thrp_sites)
            for s in thrp_sites:
                counts[s.name] += v
        else:
            counts["None"] += 1
            # print(bipartite_sites)
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
    extension = "_%s_new" %n_constant
    extension_collapsed = "%s_collapsed_new" %extension
    extension_suppcomp = "%s_suppcomp_new" %extension

    if buffer3p:
        extension = "%s_buffer3p" %extension
        extension_collapsed = "%s_buffer3p" %extension_collapsed
        extension_suppcomp = "%s_buffer3p" %extension_suppcomp

    if test:
        extension = "%s_test" %extension
        extension_collapsed = "%s_test" %extension_collapsed
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
            elif base_site in ["6mer-m8", "6mer-A1"]:
                new_base = "Offset"
            else:
                new_base = "Comp"
            new_name = "|".join([thrp_site, num, new_base])
        else:
            new_name = name
        output_suppcomp[new_name] += counts.loc[name]
    # Make a pandas dataframe from the suppcomp dictionary.
    counts_suppcomp = pd.DataFrame.from_dict(output_suppcomp, orient="index")
    # DEFINE THE PATHS TO OUTPUT THE FILES.
    site_counts_path = get_analysis_path(
        mirna, experiment, condition, "bipartite_site_counts", ext=extension
    )
    site_counts_suppcomp_path = get_analysis_path(
        mirna, experiment, condition, "bipartite_site_counts",
        ext=extension_suppcomp
    )
    print(site_counts_path)
    # print(site_counts_collapsed_path)
    print(site_counts_suppcomp_path)
    print(counts.sum())
    print(counts_suppcomp.sum())
    # Output the dataframes to the three paths.
    counts.to_csv(site_counts_path, sep="\t", header=False)
    # counts_collapsed.to_csv(site_counts_collapsed_path, sep="\t", header=False)
    counts_suppcomp.to_csv(site_counts_suppcomp_path, sep="\t", header=False)
    print_time_elapsed(time_start)
################################################################################

if __name__ == "__main__":
    main()

