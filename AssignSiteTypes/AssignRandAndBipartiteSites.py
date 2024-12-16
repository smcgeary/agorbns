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
np.core.arrayprint._line_width = 200

# FUNCTIONS
def assign_site(read_seqs, _mirna, experiment, n_constant, rand_length):
    time_start = time.time()
    ## SET UP THE DICTIONARIES REQUIRED FOR CONVERTING THE THREEPRIME ##########
    ## SEQUENCES TO THEIR NAMES, FOR PUTTING THEM INTO THE OUTPUT ##############
    ## DICTIONARIES. ###########################################################
    # Make the list of mismatch sites to identify in each read.
    seq_8mer = _mirna["8mer"]
    seq_8mer_list = list(seq_8mer)
    names_can = ["8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-m8", "6mer-A1"]
    _sitelist = SiteList(_mirna, "programmed", rand_length + 2*n_constant)

    # Make the list of each of the sequences and names that are all 16 mismatch
    # possibilities.
    seqs_8merMm = [[seq_8mer[:i] + j + seq_8mer[i+1:]
                    for j in ["A", "C", "G", "T"] if j != seq_8mer[i]]
                   for i in range(1, 7)]
    seqs_8merMm = [j for i in seqs_8merMm for j in i]
    names_8merMm = [["8mer-mm%s%s" %(j, 8 - i) for j in ["A", "C", "G", "T"]
                     if j != seq_8mer[i]]
                    for i in range(1, 7)]
    names_8merMm = [j for i in names_8merMm for j in i]
    # Make the dictionary relating each string to the name of that seed mismatch
    # type. 
    mmname_seq_map = {seq_str : seq_name for (seq_str, seq_name)

                      in zip(seqs_8merMm, names_8merMm)}
    # Generate a list of the three prime site names.
    names_thrP = [["%smer-m%s.%s" %(k, p, p + k - 1)
                   for p in range(9, len(_mirna) - k + 2)]
                  for k in range(5, 10)]
    names_thrP = [j for i in names_thrP for j in i]
    thrp_seqs = [_mirna[site_name] for site_name in names_thrP]
    thrpName_seq_map = {thrp_seq : thrp_name for (thrp_seq, thrp_name)
                        in zip(thrp_seqs, names_thrP)}
    # Make string of mirna 3-prime regoin to compare for each read.
    mirna_3p_seq_rc = get_rc(_mirna.seq[8:])
    len_mir3p = len(mirna_3p_seq_rc)

    ## MAKE THE OUTPUT DICTIONARIES ############################################
    # Make the output dictionary for instances with only the programmed site.
    progOnly_map = {name_mm : 0.0 for name_mm in names_8merMm}
    # Make the output dictionary for when seed sites, but not three prime sites,
    # are in the random region.
    singleSeed_map = {
        name_mm : {
            site : np.array([[0.0]*24])
            for site in names_can + names_8merMm
        }
        for name_mm in names_8merMm
    }
    # Make the output dictionary for when three prime sites, but not seed sites,
    # are found in the random region.
    singleThrP_map = {
        name_mm : {
            "%smer" %j : np.array([[0.0]*(25 - j + 1)]*(len_mir3p - j + 1))
            for j in [5, 6, 7, 8, 9]
        }
        for name_mm in names_8merMm
    }
    # make the output dictionary for when both seed sites and three prime sites
    # are found in the random region. Note that this output dictionary does not
    # hold onto positional information with respect to where the two sites are
    # found within the library.
    thrPAndSeed_map = {
        name_mm : {
            name_seedAndMm : {
                name_thrP : 0.0 for name_thrP in names_thrP
            }
            for name_seedAndMm in names_can + names_8merMm
        }
        for name_mm in names_8merMm
    }
    # Counter for the five possible ways a read is counted:
    # 1. Programmed site only, 2. programmed and random-position seed sites, 
    # 3. programmed and random-position three prime sites, 4. programmed and
    # both seed and three prime sites in the random-position region, and 5.
    # A read that doesn't have a mmseed in the programmed position (suggesting
    # that the seed is shifted with respect to the intended position)
    fraction_each_map = [0, 0, 0, 0, 0]
    # Define the region where the programmed mismatch site should be within the
    # read.
    lib_stop = 25
    ## ITERATE OVER THE READS TO GET THE VALUES FOR EACH OF THE FOUR OUTPUT ####
    ## DICTIONARIES ############################################################
    counter_test = 0
    for i_r, r in enumerate(read_seqs):
        print_line = False
        read = r.strip()
        _read = Read(read, rand_length, _mirna, n_constant, experiment)
        if _read.seq[lib_stop:lib_stop+8] in seqs_8merMm:
            add = True
        else:
            add = False
            fraction_each_map[4] += 1
        if add:
            # 1.A Get the name of the mismatch position. #######################
            name_mm = mmname_seq_map[_read.seq[lib_stop:lib_stop+8]]
            # if name_mm == "8mer-mmA7":
            #     print(_read.seq + "TG")
            #     counter_test += 1
            # 1.B Get all seed + mm seed sites in the random region. ###########
            _read.get_all_sites_in_read(_sitelist)
            seed_site_tuples = []
            for site in _read.sites:
                if site.name != "None":
                    if site.l <= 24:
                        seed_site_tuples += [(site.name, site.l)]
            # 1.C Look for all instances of 3p sites. ##########################
            match_3p = LCSubString(mirna_3p_seq_rc, _read.seq[:lib_stop])
            [len_match, (p_mir, p_lib)] = match_3p
            num_sites_mult = 0
            p_mir_lib_tuples = [] 
            if 5 <= len_match <= 9:
                print_line = True
                p_mir_lib_tuples = list(zip(p_mir, p_lib))
            # 2. Add the counts of this read to the appropriate output #########
            # dictionary, depending on whether there are seed+mm sites, threep #
            # sites, both, or neither. #########################################
            # 2.A Check for programmed site only. ##############################
            if len(seed_site_tuples) + len(p_mir_lib_tuples) == 0:
                # if name_mm == "8mer-mmA7":
                #     print("progOnly adding")
                #     sys.stdout.flush()
                progOnly_map[name_mm] += 1
                fraction_each_map[0] += 1
            # 2.B Check for programmed and seed site only. #####################
            elif len(p_mir_lib_tuples) == 0:
                # if name_mm == "8mer-mmA7":
                #     print("seed only")
                #     print(seed_site_tuples)

                # print("singleSeed adding")
                # sys.stdout.flush()
                count = 1./len(seed_site_tuples)
                # count = 1.
                for (name_site, pos_site) in seed_site_tuples:
                    # if name_site == "8mer-mmA7":
                    #     print("single " + name_site)
                    #     sys.stdout.flush()
                    j = pos_site
                    singleSeed_map[name_mm][name_site][0, j] += count
                fraction_each_map[2] += 1
            # 2.C Check for programmed and three prime site only. ##############
            elif len(seed_site_tuples) == 0:
                # print("singleThreeP adding")
                # sys.stdout.flush()
                # if name_mm == "8mer-mmA7":
                #     print("thrp only")
                #     print(p_mir_lib_tuples)

                count = 1./len(p_mir_lib_tuples)
                for (p_mir_i, p_lib_i) in p_mir_lib_tuples:
                    i, j = p_mir_i, p_lib_i
                    key_thrP = "%smer" %(len_match)
                    singleThrP_map[name_mm][key_thrP][i, j] += count
                fraction_each_map[1] += 1
            # 2.D Check for programmed and seed, and three prime site only. ####
            else:
                # print("ThreePAndSeed adding")
                # sys.stdout.flush()
                # print(name_mm)
                # if name_mm == "8mer-mmA7" and "8mer" in [i[0] for i in seed_site_tuples]:
                #     print("__________________________________")
                    # sys.stdout.flush()
                count = 1./(len(seed_site_tuples) * len(p_mir_lib_tuples))
                # count = 1.
                for (p_mir_i, p_lib_i) in p_mir_lib_tuples:
                    p_l = len_mir3p + 8 + 1 - len_match - p_mir_i
                    p_r = p_l + len_match - 1
                    name_thrP = "%smer-m%s.%s" %(len_match, p_l, p_r)
                    # print(name_thrP)
                    # sys.stdout.flush()
                    for (name_site, pos_site) in seed_site_tuples:
                        # if name_site == "8mer-mmA7":
                        #     print(seed_site_tuples)
                        #     print(p_mir_lib_tuples)
                        #     print(name_site)
                        #     print(name_thrP)
                        #     print(count)
                        #     sys.stdout.flush()
                        thrPAndSeed_map[name_mm][name_site][name_thrP] += count
                fraction_each_map[3] += 1
    print("got to the end of the loop")
    print(fraction_each_map)
    fraction_each_map = [i/float(i_r + 1) for i in fraction_each_map]
    fract_names = ["Programmed_Only", "Seed in random region",
                   "ThreeP in random region",
                   "Seed and ThreeP in random region",
                   "Incorrect programmed region"]
    fract_tuples = list(zip(fract_names, fraction_each_map))

    for fract_tup_i in fract_tuples:
        print("%s: %1.0f%%" %(fract_tup_i[0], fract_tup_i[1]*100))
    
    return(progOnly_map, singleSeed_map, singleThrP_map, thrPAndSeed_map)

def main():
    time_start = time.time()
    arguments = ["miRNA", "experiment", "condition", "n_constant", "-jobs",
                 "-test_binary"]
    args = parse_arguments(arguments)
    (mirna, experiment, condition, n_constant, jobs, test) = args

    # Determine the length of the "random" portion.
    # `Random` is in quotes because in the programmed libraries positions 26-31
    # are not actually random.
    if "let-7a" in mirna and experiment in ["equilibrium_mmseed_nb",
                                            "equilibrium_seed_nb"]:
        read_length = 38
    else:
        read_length = 37
    _mirna = Mirna(mirna)
    extension = "_%s" %n_constant
    if test:
        extension = "%s_test" %(extension)

    input_list = [(mirna, experiment, condition)]
    if not jobs:
        jobs = 20
    args  = [_mirna, experiment, int(n_constant), read_length]
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
    mm_keys = list(progOnly_map.keys())
    # Get the seed and mismatch keys:
    seed_keys = list(singleSeed_map[mm_keys[0]].keys())
    thrL_keys = list(singleThrP_map[mm_keys[0]].keys())
    thrP_keys = list(thrPAndSeed_map[mm_keys[0]][seed_keys[0]].keys())

    ## ADD UP ALL THE COUNTS FROM THE THREADS FOR ALL FOUR DICTIONARIES ########
    for thread in threads[1:]:
        for mm_key in mm_keys:
            progOnly_map[mm_key] += thread[0][mm_key]
            for seed_key in seed_keys:
                singleSeed_map[mm_key][seed_key] += thread[1][mm_key][seed_key]
                for thrP_key in thrP_keys:
                    thrPAndSeed_map[mm_key][seed_key][thrP_key] += thread[3][mm_key][seed_key][thrP_key]
            for thrL_key in thrL_keys:
                singleThrP_map[mm_key][thrL_key] += thread[2][mm_key][thrL_key]

    print(progOnly_map["8mer-mmA7"])
    print(np.sum(np.concatenate(list(singleSeed_map["8mer-mmA7"].values()))))
    # print([np.sum(i) for i in list(singleThrP_map["8mer-mmA7"].values())])
    print(np.sum([np.sum(i) for i in list(singleThrP_map["8mer-mmA7"].values())]))
    print(np.sum([np.sum(j) for i in list(thrPAndSeed_map["8mer-mmA7"].values()) for j in list(i.values())]))

    print("double 8mers:")
    print(np.sum(list(thrPAndSeed_map["8mer-mmA7"]["8mer"].values())))

    ## MAKE THE FINAL SINGLE COUNT WHERE THE VARIOUS KEYS AND DIMENSIONS OF ####
    ## FOUR OUTPUT DICTIONARIES ARE TRANSFORMED INTO SITE NAMES ################
    # 1. Initialize the final dataframe with the counts of reads with only
    # programmed mismatch sites.
    counts = pd.DataFrame.from_dict(progOnly_map, orient="index")
    counts.index = ["NA|NA|%s" %i for i in counts.index]
    print("added df 1")
    # print(counts)
    progAndSeed_map = {}
    # 2. Add the counts from each of the 16 seed-in-random-region dictionaries.
    for mm_key in mm_keys:
        for seed_key in seed_keys:
            len_seed = int(seed_key.split("mer")[0])
            for j in range(24):
                key = "%s|%s|%s" %(seed_key, 25 - j - len_seed + 8 + 1, mm_key)
                val = singleSeed_map[mm_key][seed_key][0,j]
                progAndSeed_map[key] = val
    counts = counts.append(
        pd.DataFrame.from_dict(progAndSeed_map, orient="index")
    )
    print("added df 2")
    # print(counts)
    # 3. Add the counts from each of the 16X5 threep-in-random-region dictionaries.
    progAndThrP_map = {}
    for mm_key in mm_keys:
        for thrL_key in thrL_keys:
            len_thrP = int(thrL_key.split("mer")[0])
            array_thrP_map = singleThrP_map[mm_key][thrL_key]
            # Get the dimensions of the 2D array.
            n_mir, n_lib = array_thrP_map.shape
            # Flatten the array into one list.
            array_list = [j for i in array_thrP_map.tolist() for j in i]
            # Calculate the names that correspond to each of the count values
            # within the list.
            # The list was flattened maintaining the ordering of the rows,
            # therefore the inner list comprehension iterates over the columns,
            # i.e., the positions within the library.
            rownames = ["%s-m%s.%s" %(thrL_key, n_mir + 8 - i, n_mir + 8 - i + len_thrP - 1)
                        for i in range(n_mir)]
            names = [["%s-m%s.%s|%s|%s" %(thrL_key, # 5mer"
                                          n_mir + 8 - i,    # 9 
                                          n_mir + 8 - i + len_thrP - 1, # 9 + 5 -1 = 13
                                          25 - j - len_thrP + 1 + 8, # 25 - 0 - (5 + 1) + 8 = 29
                                          mm_key) for j in range(n_lib)]
                     for i in range(n_mir)]
            names = [j for i in names for j in i]
            if thrL_key == "7mer":
                print(singleThrP_map[mm_key][thrL_key])
                print(rownames)
                # print(zip(names, array_list))
                # print(names)
                # print("length of array list:")
                # print(len(array_list))
                # print("length of names:")
                # print(len(names))
            # Iterate over zipped tuples of the names and counts to add them
            for name, count in zip(names, array_list):
                progAndThrP_map[name] = count
    counts = counts.append(
        pd.DataFrame.from_dict(progAndThrP_map, orient="index")
    )

    print("added df 3")
    # print(counts)
    # 4. Add the counts from each of the seed_and_threep-in-random-region
    # dictionaries.
    progAndSeedAndThrP_map = {}
    for mm_key in mm_keys:
        for seed_key in seed_keys:
            for thrP_key in thrP_keys:
                name = "%s&%s|NA|%s" %(seed_key, thrP_key, mm_key)
                count = thrPAndSeed_map[mm_key][seed_key][thrP_key]
                progAndSeedAndThrP_map[name] = count
    # print(pd.DataFrame.from_dict(progAndSeedAndThrP_map, orient="index"))
    counts = counts.append(
        pd.DataFrame.from_dict(progAndSeedAndThrP_map, orient="index")
    )
    print("added df 4")
    # print(counts)
    # count_names = list(counts.index)
    # for string in count_names[:10]:
    #     name_1 = string.split("|")[2]
    #     print(name_1)
    # inds_keep = [i for i, string in enumerate(count_names) if string.split("|")[2] == "8mer-mmA7"]
    # print(inds_keep[:10])
    # count_sum = counts.iloc[inds_keep, 0]
    # tally = np.sum(count_sum)
    # print(tally)
    # print(counts.loc["NA|NA|8mer-mmA7"])
    print(counts.loc[["7mer-m11.17|11|8mer-mmA7", "7mer-m11.17|12|8mer-mmA7",
                      "7mer-m11.17|13|8mer-mmA7", "7mer-m11.17|14|8mer-mmA7",
                      "7mer-m11.17|15|8mer-mmA7", "7mer-m11.17|16|8mer-mmA7",
                      "7mer-m11.17|17|8mer-mmA7", "7mer-m11.17|18|8mer-mmA7",
                      "7mer-m11.17|19|8mer-mmA7", "7mer-m11.17|20|8mer-mmA7",
                      "7mer-m11.17|21|8mer-mmA7", "7mer-m11.17|22|8mer-mmA7",
                      "7mer-m11.17|23|8mer-mmA7", "7mer-m11.17|24|8mer-mmA7"]])

    print(counts.loc[[ "8mer|9|8mer-mmA7", "8mer|10|8mer-mmA7",
                      "8mer|11|8mer-mmA7", "8mer|12|8mer-mmA7",
                      "8mer|13|8mer-mmA7", "8mer|14|8mer-mmA7",
                      "8mer|15|8mer-mmA7", "8mer|16|8mer-mmA7",
                      "8mer|17|8mer-mmA7", "8mer|18|8mer-mmA7",
                      "8mer|19|8mer-mmA7", "8mer|20|8mer-mmA7",
                      "8mer|21|8mer-mmA7", "8mer|22|8mer-mmA7",
                      "8mer|23|8mer-mmA7", "8mer|24|8mer-mmA7",
                      "8mer|25|8mer-mmA7", "8mer|26|8mer-mmA7"]])

    # inds_keep_double = [i for i, string in enumerate(count_names) if "&" in string]
    # counts_double = counts.iloc[inds_keep_double, 0]
    # count_names_double = list(counts_double.index)

    # print("Counts sum:")
    # inds_keep_double_8mer = [i for i, string in enumerate(count_names_double) if string.split("&")[0] == "8mer" and string.split("|")[2] == "8mer-mmA7"]
    # counts_double_8mer = counts_double.iloc[inds_keep_double_8mer]

    # print(counts_double_8mer)
    # tally2 = np.sum(counts_double_8mer)
    # print(tally2)


    # return
    # The names of the output files:
    site_counts_path = get_analysis_path(mirna, experiment, condition,
                                       "programmed_site_counts", ext=extension)
    # multisite_counts_path = get_analysis_path(mirna, experiment, condition,
    #                                         "multisite_counts",
    #                                         ext=extension)
    # single_pos_counts_path = get_analysis_path(mirna, experiment, condition,
    #                                         "single_pos_counts",
    #                                         ext=extension)
    # top_pos_counts_path = get_analysis_path(mirna, experiment, condition,
    #                                         "top_pos_counts",
    #                                         ext=extension)
    # outputpaths = [get_analysis_path(mirna, experiment, condition, i,
    #                                  ext=extension)
    #                for i in ["site_counts", "multisite_counts",
    #                          "single_pos_counts", "top_pos_counts"]]

    # for outputpath in outputpaths:
    #     print(outputpath)
    print(site_counts_path)
    counts.to_csv(site_counts_path, sep="\t", header=False)
    # # multicounts.to_csv(multisite_counts_path, sep="\t", header=False)
    # if sitelist not in ["12mers", "16mers"]:
    #     single_pos_counts.to_csv(single_pos_counts_path, sep="\t", header=positional_names)
    #     top_pos_counts.to_csv(top_pos_counts_path, sep="\t", header=positional_names)
    print_time_elapsed(time_start)
################################################################################

if __name__ == "__main__":
    main()

