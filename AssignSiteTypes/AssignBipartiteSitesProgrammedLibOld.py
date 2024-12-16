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
def assign_site(read_seqs, _mirna, experiment, n_con, n_lib):
    time_start = time.time()
    # Define the left and right hand sides of the programmed site.
    prog_l = 25 + n_con
    prog_r = prog_l + 8
    # if end3prand:
    thrp_r = prog_l
    ## SET UP THE DICTIONARIES REQUIRED FOR CONVERTING THE THREEPRIME ##########
    ## SEQUENCES TO THEIR NAMES, FOR PUTTING THEM INTO THE OUTPUT ##############
    ## DICTIONARIES. ###########################################################
    # Make the list of all possible programmed sites in each read.
    seq_8mer = _mirna["8mer"]
    seq_8mer_list = list(seq_8mer)
    names_can = ["8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-m8", "6mer-A1"]
    _sitelist = SiteList(_mirna, "programmed", n_lib + 2*n_con)

    # Make the list of each of the sequences and names that are all of the 6
    # canonical site and the 16 mismatch possibilities.
    # Get the sequences of all canonical and single mismatch sites.
    seqs_7m8 = [seq_8mer[:7] + i for i in DNTS if i != seq_8mer[7]]
    seqs_7A1 = [i + seq_8mer[1:] for i in DNTS if i != seq_8mer[0]]
    seqs_6   = [i + seq_8mer[1:7] + j
                for i in DNTS for j in DNTS
                if i != seq_8mer[0] if j != seq_8mer[7]]
    seqs_6m8 = [seq_8mer[:6] + i + j
                for i in DNTS for j in DNTS
                if i != seq_8mer[6] if j != seq_8mer[7]]
    seqs_6A1 = [i + j + seq_8mer[2:]
                for i in DNTS for j in DNTS
                if i != seq_8mer[0] if j != seq_8mer[1]]
    seqs_8merMm = [seq_8mer[:i] + j + seq_8mer[i+1:]
                   for i in range(1, 7)
                   for j in ["A", "C", "G", "T"] if j != seq_8mer[i]]
    # Define the names and sequence of the canonical sites.
    names_can_rep = (
        ["8mer"] + 3*["7mer-m8"] + 3*["7mer-A1"] +
        9*["6mer"] + 9*["6mer-m8"] + 9*["6mer-A1"]
    )
    seqs_can = [seq_8mer] + seqs_7m8 + seqs_7A1 + seqs_6 + seqs_6m8 + seqs_6A1
    # Get the names of the 8mer 1-nt mismatch sites.
    names_8merMm = ["8mer-mm%s%s" %(j, 8 - i)
                    for i in range(1, 7)
                    for j in ["A", "C", "G", "T"] if j != seq_8mer[i]]
    # Assin the names and sequences of all supplemental and compensatory sites
    # in the library (which may or may not all be included in the library, but
    # are expected to be there at least at a low frequency due to 1-nt
    # mismatches 
    names_supp_comp = names_can_rep + names_8merMm
    seqs_supp_comp = seqs_can + seqs_8merMm
    # Make the dictionary relating each string to the name of that seed mismatch
    # type. 
    programmed_seq_map = {seq_str : seq_name for (seq_str, seq_name)
                          in zip(seqs_supp_comp, names_supp_comp)}
    # Generate a list of the three prime site names, site sequences, and a
    # dictionary with them as values and keys, respectively.
    # if (outto11mers):
    kmer_max = 12
    # else:
    #     kmer_max = 10
    # if downto2mers:
    #     kmer_min = 2
    # elif downto3mers:
    #     kmer_min = 3
    # elif downto4mers:
    kmer_min = 4
    # else:
    #     kmer_min = 5
    names_thrp = ["%smer-m%s.%s" %(k, p, p + k - 1)
                  for k in range(kmer_min, kmer_max)
                  for p in range(9, len(_mirna) - k + 2)]
    seqs_thrp = [_mirna[site_name] for site_name in names_thrp]
    thrpName_seq_map = {thrp_seq : thrp_name for (thrp_seq, thrp_name)
                        in zip(seqs_thrp, names_thrp)}
    # Make string of mirna 3-prime region to compare for each read.
    mirna_3p_seq_rc = get_rc(_mirna.seq[8:])
    len_mir3p = len(mirna_3p_seq_rc)

    ## MAKE THE OUTPUT DICTIONARIES ############################################
    # Make the output dictionary for instances with only the programmed site.
    progOnly_map = {site : 0.0 for site in names_supp_comp}
    # Make the output dictionary for when seed sites, but not three prime sites,
    # are in the random region.
    singleSeed_map = {
        prog_site : {              # Lenth of read          #   # len site #  # correction#
            site : np.array([[0.0]*(n_lib + 2*n_con - int(site[0]) + 1)])
            for site in names_can + names_8merMm
        }
        for prog_site in names_supp_comp
    }

    # Make the output dictionary for when three prime sites, but not seed sites,
    # are found in the random region.
    singleThrP_map = {
        prog_site : {
            "%smer" %n_k : np.array(
                        # Lenth of read #
                [[0.0]*(n_lib + 2*n_con - n_k + 1)]*(len_mir3p - n_k + 1)
            )
            for n_k in range(kmer_min, kmer_max)
        }
        for prog_site in names_supp_comp
    }
    # make the output dictionary for when both seed sites and three prime sites
    # are found in the random region. Note that this output dictionary does not
    # hold onto positional information with respect to where the two sites are
    # found within the library.
    thrPAndSeed_map = {
        prog_site : {
            name_seedAndMm : {
                name_thrp : 0.0 for name_thrp in names_thrp
            }
            for name_seedAndMm in names_can + names_8merMm
        }
        for prog_site in names_supp_comp
    }
    # Counter for the five possible ways a read is counted:
    # 1. Programmed site only, 2. programmed and random-position seed sites, 
    # 3. programmed and random-position three prime sites, 4. programmed and
    # both seed and three prime sites in the random-position region, and 5.
    # A read that doesn't have a mmseed in the programmed position (suggesting
    # that the seed is shifted with respect to the intended position)
    fraction_each_map = [0, 0, 0, 0, 0, 0]

    average_per_read_map = {kmer : (0, 0) for kmer in range(kmer_min, kmer_max)}
    print(average_per_read_map)
    sys.stdout.flush()

    # Define the region where the programmed mismatch site should be within the
    # read.
    ## ITERATE OVER THE READS TO GET THE VALUES FOR EACH OF THE FOUR OUTPUT ####
    ## DICTIONARIES ############################################################
    counter_test = 0
    # variable defining regions to be excluded for non-programmed, seed sites.
    ex_pos = [i + n_lib for i in [25, 26, 27]]
    count_8mer = 0
    for i_r, r in enumerate(read_seqs):
        check_read = False
        print_line = False
        read = r.strip()
        _read = Read(read, n_lib, _mirna, n_con, experiment)
        # Grab the sequence of the programmed site.
        prog_kmer = _read.seq[prog_l:prog_r]
        if prog_kmer in seqs_supp_comp:
            # 1.A Get the name of the programmed site. #######################
            name_prog = programmed_seq_map[prog_kmer]
            if name_prog == "6mer-A1":
                prog_l_s = prog_l + 2
            elif name_prog in ["7mer-A1", "6mer"]:
                prog_l_s = prog_l + 1
            else:
                prog_l_s = prog_l
            prog_tuple = (name_prog, prog_l_s)
            # 1.B Get all seed + mm seed sites in the random region. ###########
            _read.get_all_sites_in_read(_sitelist)
            seed_tuples = [(site.name, site.l) for site in _read.sites]
            # Reassign the 
            seed_tuples = [i for i in seed_tuples if i != prog_tuple]
            seed_names = [i[0] for i in seed_tuples]
            seed_starts = [i[1] for i in seed_tuples]
            seed_stops = [seed_start + int(seed_name[0]) for (seed_start, seed_name) in zip(seed_starts, seed_names)]
            seed_stops = [i[1] + int(i[0][0]) for i in seed_tuples]
            # 1.C Look for all instances of 3p sites. ##########################
            [len_match, (p_mir_rc_temp, p_lib_temp)] = LCSubString(mirna_3p_seq_rc,
                                                      _read.seq[:thrp_r])
            # This portion of the script was added in order to check for short
            # three-prime sites that overlap with the seed or seed-mismatch
            # sites. They work by making sure that both of the ends of
            # each of the mismatch sites are at least 2 nucleotides away from the ends of the 
            p_mir_rc = []
            p_lib = []
            # Loop to get rid of three prime sites that are hovering over the
            # programmed region.
            for i, p_lib_temp_i in enumerate(p_lib_temp):
                if abs(p_lib_temp_i - prog_l) > 2 and abs(p_lib_temp_i + len_match - prog_r) > 2:
                    add_sequence = True
                    # if check_read:
                    #     print(p_mir_rc_temp[i])
                    #     print(p_lib_temp_i)
                    #     print(seed_tuples)
                    #     print(_read.seq)
                    #     print(' '*p_lib_temp_i + "_"*len_match)
                    #     sys.stdout.flush()
                    p_r_temp_i = p_lib_temp_i + len_match
                    for (seed_start, seed_stop) in zip(seed_starts, seed_stops):
                        # if check_read:
                        #     print(" "*seed_start + "*"*(seed_stop - seed_start))
                        if p_lib_temp_i >= seed_start and p_r_temp_i <= seed_stop:
                            print("this is true")
                            sys.stdout.flush()
                            add_sequence = False
                    if add_sequence:
                        p_mir_rc.append(p_mir_rc_temp[i])
                        p_lib.append(p_lib_temp[i])
            # Use all three prime sites 4 nt or longer.
            if (len_match >= kmer_min and len(p_mir_rc) != 0 and
                len_match < kmer_max):
                thrp_tuples = list(zip(p_mir_rc, p_lib))
                # if check_read:
                #     print(thrp_tuples)
                #     print(len_match)
                #     sys.stdout.flush()
            else:
                thrp_tuples = []

            if name_prog == "8mer-mmA7" and "8mer" in [i[0] for i in seed_tuples]:
                print("_________________________________")
                print(seed_tuples)
                print(thrp_tuples)
                if len(thrp_tuples) == 0:
                    print(i_r)
                sys.stdout.flush()
            # 2. Add the counts of this read to the appropriate output #########
            # dictionary, depending on whether there are seed+mm sites, threep #
            # sites, both, or neither. #########################################
            # 2.A Check for programmed site only. ##############################
            if len(seed_tuples) + len(thrp_tuples) == 0:
                progOnly_map[name_prog] += 1
                fraction_each_map[0] += 1
                # if check_read:
                #     print("in 1; programmed only, and len_match is %s." %(len_match))
                #     sys.stdout.flush()
            # 2.B Check for programmed and seed site only. #####################
            elif len(thrp_tuples) == 0:
                # if check_read:
                #     print("in 2; programmed and seed sites only.")
                #     print(seed_tuples)
                #     print("___________________________________________________")
                #     sys.stdout.flush()
                count = 1./len(seed_tuples)
                name_sites = [i[0] for i in seed_tuples]
                if "8mer" in name_sites and check_read:
                    print("%s used in count" %i_r)
                    sys.stdout.flush()
                for (name_site, pos_site) in seed_tuples:
                    singleSeed_map[name_prog][name_site][0, pos_site] += count
                fraction_each_map[2] += 1
            # 2.C Check for programmed and three prime site only. ##############
            elif len(seed_tuples) == 0:
                # if check_read:
                #     print("in 3; programmed and threep only")
                #     print(len_match)
                #     print(thrp_tuples)
                #     print(1/len(thrp_tuples))
                #     print("___________________________________________________")
                #     sys.stdout.flush()
                count = 1./len(thrp_tuples)
                (len_ave, len_num) = average_per_read_map[len_match]
                len_ave_new = (len_ave*len_num + len(thrp_tuples))/(len_num + 1.)
                average_per_read_map[len_match] = (len_ave_new, len_num + 1.)
                for (p_mir_rc, p_lib) in thrp_tuples:
                    n_thrp = "%smer" %(len_match)
                    singleThrP_map[name_prog][n_thrp][p_mir_rc, p_lib] += count
                fraction_each_map[1] += 1
            # 2.D Check for programmed and seed, and three prime site only. ####
            else:
                # if check_read:
                #     print("in 4; programmed, threep, and seed sites.")
                #     print(seed_tuples)
                #     print(thrp_tuples)
                #     print("length of match is %s" %len_match)
                #     print("___________________________________________________")
                #     sys.stdout.flush()
                count = 1./(len(seed_tuples) * len(thrp_tuples))
                for (p_mir_rc, p_lib) in thrp_tuples:
                    p_l = len_mir3p + 8 + 1 - len_match - p_mir_rc
                    p_r = p_l + len_match - 1
                    name_thrP = "%smer-m%s.%s" %(len_match, p_l, p_r)
                    for (name_site, pos_site) in seed_tuples:
                        thrPAndSeed_map[name_prog][name_site][name_thrP] += count
                fraction_each_map[3] += 1
        else:
            # if check_read:
            #     print("in 5")
            #     sys.stdout.flush()
            fraction_each_map[4] += 1

    print("got to the end of the loop")
    print(fraction_each_map)
    fraction_each_map = [i/float(i_r + 1) for i in fraction_each_map]
    fract_names = ["Programmed_Only", "Seed in random region",
                   "ThreeP in random region",
                   "Seed and ThreeP in random region",
                   "Incorrect programmed region",
                   "Larger than 9 nt of 3p end complementarity"]
    fract_tuples = list(zip(fract_names, fraction_each_map))

    for fract_tup_i in fract_tuples:
        print("%s: %1.0f%%" %(fract_tup_i[0], fract_tup_i[1]*100))
    sys.stdout.flush()
    print("8mer count: %s" %(count_8mer))
    print(average_per_read_map)
    print("map at end:")
    print(singleSeed_map["8mer-mmA7"]["8mer"])

    sys.stdout.flush()
    return (progOnly_map, singleSeed_map, singleThrP_map, thrPAndSeed_map)

def main():
    time_start = time.time()
    arguments = [
        "miRNA", "experiment", "condition", "n_constant", "-jobs", "-test_binary"
    ]
    args = parse_arguments(arguments)
    (
        mirna, experiment, condition, n_constant, jobs, test
    ) = args
    # Determine the length of the "random" portion.
    # `Random` is in quotes because in the programmed libraries positions 26-31
    # are not actually random.
    if ("let-7a" in mirna or "miR-1" == mirna) and mirna != "miR-155_let-7a":
        rand_length = 38
    else:
        rand_length = 37
    _mirna = Mirna(mirna)
    # Assign three different extensions for the three different count files
    # being made.
    extension = "_%s_old" %n_constant
    extension_collapsed = "%s_collapsed_old" %extension
    extension_suppcomp = "%s_suppcomp_old" %extension
    if test:
        extension = "%s_test" %extension
        extension_collapsed = "%s_test" %extension_collapsed
        extension_suppcomp = "%s_test" %extension_suppcomp

    input_list = [(mirna, experiment, condition)]
    if not jobs:
        jobs = 20
    args  = [_mirna, experiment, int(n_constant), rand_length]
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
        threads += multiproc_file(reads_path, int(jobs), assign_site, test,
                                  *args)

    ## INITIALIZE THE OUTPUT DICTIONARY KEYS AND OBJECTS TO COMBINE THE ########
    ## RESULTS FROM EACH MULTIPROCESSING THREAD ################################

    progOnly_map, singleSeed_map, singleThrP_map, thrPAndSeed_map = threads[0]
    # Use the first thread to get all of the dictionary keys (this is probably
    # bad coding).
    # Get the mismatch keys:
    prog_keys = list(progOnly_map.keys())
    # Get the seed and mismatch keys:
    seed_keys = list(singleSeed_map[prog_keys[0]].keys())
    thrL_keys = list(singleThrP_map[prog_keys[0]].keys())
    thrP_keys = list(thrPAndSeed_map[prog_keys[0]][seed_keys[0]].keys())

    ## ADD UP ALL THE COUNTS FROM THE THREADS FOR ALL FOUR DICTIONARIES ########
    for thread in threads[1:]:
        # Get the programmed key for the first dictionaries.
        for prog_key in prog_keys:
            val = thread[0][prog_key]
            progOnly_map[prog_key] += val
            # Get the seed key for the second and fourth dictionaries.
            for seed_key in seed_keys:
                val = thread[1][prog_key][seed_key]
                singleSeed_map[prog_key][seed_key] += val
                # Get the threeP site key for the fourth dictionary.
                for thrP_key in thrP_keys:
                    val = thread[3][prog_key][seed_key][thrP_key]
                    thrPAndSeed_map[prog_key][seed_key][thrP_key] += val
            # Get the threeP length key for the third dictionary.
            for thrL_key in thrL_keys:
                val = thread[2][prog_key][thrL_key]
                singleThrP_map[prog_key][thrL_key] += val


    ## MAKE THE FINAL SINGLE COUNT WHERE THE VARIOUS KEYS AND DIMENSIONS OF ####
    ## FOUR OUTPUT DICTIONARIES ARE TRANSFORMED INTO SITE NAMES ################
    # 1. Initialize the final dataframe with the counts of reads with only
    # programmed mismatch sites.
    counts = pd.DataFrame.from_dict(progOnly_map, orient="index")
    counts_collapsed = counts.copy()
    counts_suppcomp = counts.copy()
    base_sites = progOnly_map.keys()
    # Make a dictionary that maps each of the 24 possible programmed site to
    # either "Supp" or "Comp".
    progSite_suppCompSite_map = {base_site : collapsed_site for
                                 (base_site, collapsed_site) in
                                 zip(base_sites, ["Supp"]*6 + ["Comp"]*18)}

    counts.index = ["NA|NA|%s" %i for i in counts.index]
    progAndSeed_map = {}
    progAndSeed_collapsed_map = {}
    progAndSeed_suppcomp_map = defaultdict(int)


    # 2. Add the counts from each of the 16 seed-in-random-region dictionaries.
    for prog_key in prog_keys:
        base_key = progSite_suppCompSite_map[prog_key]
        for seed_key in seed_keys:
            collapsed_val = 0
            key_collapsed = "%s|%s" %(seed_key, prog_key)
            df_i = singleSeed_map[prog_key][seed_key]
            len_seed = int(seed_key.split("mer")[0])
            mir_offset = 25 + int(n_constant) - len_seed + 8 + 1
            for j in range(df_i.shape[1]):
                site_pos = -j + mir_offset
                if site_pos < 1:
                    site_pos -= 1
                key = "%s|%s|%s" %(seed_key, site_pos, prog_key)
                key_suppcomp = "%s|%s|%s" %(seed_key, site_pos, base_key)
                val = df_i[0, j]
                progAndSeed_map[key] = val
                progAndSeed_suppcomp_map[key_suppcomp] += val
                collapsed_val += val
            progAndSeed_collapsed_map[key_collapsed] = collapsed_val
    counts = counts.append(
        pd.DataFrame.from_dict(progAndSeed_map, orient="index")
    )
    counts_collapsed = counts_collapsed.append(
        pd.DataFrame.from_dict(progAndSeed_collapsed_map, orient="index")
    )
    counts_suppcomp = counts_suppcomp.append(
        pd.DataFrame.from_dict(progAndSeed_suppcomp_map, orient="index")
    )

    print("added df 2")
    # 3. Add the counts from each of the 16X5 threep-in-random-region
    # dictionaries.
    progAndThrP_map = {}
    progAndThrP_collapsed_map = {}
    progAndThrP_suppcomp_map = defaultdict(int)

    for prog_key in prog_keys:
        base_key = progSite_suppCompSite_map[prog_key]
        for thrL_key in thrL_keys:
            len_thrP = int(thrL_key.split("mer")[0])
            array_thrP_map = singleThrP_map[prog_key][thrL_key]
            # Get the dimensions of the 2D array.
            n_mir, n_lib = array_thrP_map.shape
            # Flatten the array into one list.
            array_list = [j for i in array_thrP_map.tolist() for j in i]
            # Calculate the names that correspond to each of the count values
            # within the list.
            # The list was flattened maintaining the ordering of the rows,
            # therefore the inner list comprehension iterates over the columns,
            # i.e., the positions within the library.

            mir_len = n_mir + 8
            off_l = mir_len
            off_r = mir_len + len_thrP - 1
            off_mirp = 25 + int(n_constant) - len_thrP + 1 + 8
            name_base = ["%s-m%s.%s" %(thrL_key, off_l - i, off_r - i)
                         for i in range(n_mir)]
            names = ["%s|%s|%s" %(base, off_mirp - j, prog_key)
                     for base in name_base for j in range(n_lib)]
            names_collapsed = ["%s|%s" %(base, prog_key)
                               for base in name_base for j in range(n_lib)]
            names_suppcomp = ["%s|%s|%s" %(base, off_mirp - j, base_key)
                               for base in name_base for j in range(n_lib)]
            for name in names_collapsed:
                progAndThrP_collapsed_map[name] = 0.0
            # Iterate over zipped tuples of the names and counts to add them
            zip_iter = zip(names, names_collapsed, names_suppcomp,
                               array_list)
            for name, name_collapsed, name_suppcomp, count in zip_iter:
                progAndThrP_map[name] = count
                progAndThrP_collapsed_map[name_collapsed] += count
                progAndThrP_suppcomp_map[name_suppcomp] += count
    print("collapsed 3 length:")
    counts = counts.append(
        pd.DataFrame.from_dict(progAndThrP_map, orient="index")
    )
    counts_collapsed = counts_collapsed.append(
        pd.DataFrame.from_dict(progAndThrP_collapsed_map, orient="index")
    )
    counts_suppcomp = counts_suppcomp.append(
        pd.DataFrame.from_dict(progAndThrP_suppcomp_map, orient="index")
    )

    print("added df 3")
    progAndSeedAndThrP_map = {}
    progAndSeedAndThrP_suppcomp_map = defaultdict(int)

    for prog_key in prog_keys:
        base_key = progSite_suppCompSite_map[prog_key]
        for seed_key in seed_keys:
            for thrP_key in thrP_keys:
                name = "%s&%s|NA|%s" %(seed_key, thrP_key, prog_key)
                suppcomp_name = "%s&%s|NA|%s" %(seed_key, thrP_key, base_key)
                count = thrPAndSeed_map[prog_key][seed_key][thrP_key]
                progAndSeedAndThrP_map[name] = count
                progAndSeedAndThrP_suppcomp_map[suppcomp_name] += count
    print("collapsed 4 length:")
    print(len(progAndSeedAndThrP_map))
    counts = counts.append(
        pd.DataFrame.from_dict(progAndSeedAndThrP_map, orient="index")
    )
    counts_collapsed = counts_collapsed.append(
        pd.DataFrame.from_dict(progAndSeedAndThrP_map, orient="index")
    )
    counts_suppcomp = counts_suppcomp.append(
        pd.DataFrame.from_dict(progAndSeedAndThrP_suppcomp_map, orient="index")
    )
    print("added df 4")
    # The names of the output file:
    site_counts_path = get_analysis_path(
        mirna, experiment, condition, "programmed_site_counts", ext=extension
    )
    site_counts_collapsed_path = get_analysis_path(
        mirna, experiment, condition, "programmed_site_counts",
        ext=extension_collapsed
    )
    site_counts_suppcomp_path = get_analysis_path(
        mirna, experiment, condition, "programmed_site_counts",
        ext=extension_suppcomp
    )
    print(site_counts_path)
    print(site_counts_collapsed_path)
    print(site_counts_suppcomp_path)
    # Output the dataframes to the three paths.
    counts.to_csv(site_counts_path, sep="\t", header=False)
    counts_collapsed.to_csv(site_counts_collapsed_path, sep="\t", header=False)
    counts_suppcomp.to_csv(site_counts_suppcomp_path, sep="\t", header=False)
    print_time_elapsed(time_start)
################################################################################

if __name__ == "__main__":
    main()

