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


# Added part to extract the Kd values associated with the reporter library
# designs.




## GLOBAL CONSTANTS ############################################################

# FUNCTIONS
def assign_site(read_seqs, _mirna, experiment, n_con, n_lib):
    time_start = time.time()

    # 1. DEFINE FUNCTION CONSTANTS.
    # Left and right hand sides of the programmed site.
    l_libprog = 25 + n_con
    r_libprog = l_libprog + 8
    # The minimum and maximum threeprime kmer lengths queried (the max is
    # written as 12 and not 11 due to pythonic indexing.
    kmer_min = 4
    kmer_max = 12
    # Make string of mirna 3-prime region to compare for each read.
    mirna_3p_seq_rc = get_rc(_mirna.seq[8:])
    len_mir3p = len(mirna_3p_seq_rc)
    len_mir = len(_mirna)
    ## 2. SET UP THE DICTIONARIES REQUIRED FOR QUERYING THE PROGRAMMED REGION ##
    ## AND THE THREEPRIME COMPLEMENTARITY WITHIN EACH READ. ####################
    # Make the list of all possible programmed sites in each read.
    seq_8mer = _mirna["8mer"]
    _sitelist = SiteList(_mirna, "programmed", n_lib + 2*n_con)
    # 2.i.a Define the sequences of all canonical and single mismatch sites.
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
    seqs_can = [seq_8mer] + seqs_7m8 + seqs_7A1 + seqs_6 + seqs_6m8 + seqs_6A1
    seqs_comp = [seq_8mer[:i] + j + seq_8mer[i+1:] for i in range(1, 7)
                 for j in DNTS if j != seq_8mer[i]]
    # 2.i.b Define the names the canonical sites and single mismatch sites.
    names_can = (
        ["8mer"] + 3*["7mer-m8"] + 3*["7mer-A1"] +
        9*["6mer"] + 9*["6mer-m8"] + 9*["6mer-A1"]
    )
    names_comp = ["8mer-mm%s%s" %(j, 8 - i)
                    for i in range(1, 7)
                    for j in DNTS if j != seq_8mer[i]]
    names_suppcomp = names_can + names_comp
    seqs_suppcomp = seqs_can + seqs_comp
    # 2.i.c Make the programmed site dictionary.
    programmedName_seq_map = {seq_suppcomp : name_suppcomp
                              for (seq_suppcomp, name_suppcomp)
                              in zip(seqs_suppcomp, names_suppcomp)}
    # 2.ii.a Define the threeprime site names.
    names_thrp = ["%smer-m%s.%s" %(k, p, p + k - 1)
                  for k in range(kmer_min, kmer_max)
                  for p in range(9, len(_mirna) - k + 2)]
    # 2.ii.b Make the threprime site dictionary.
    thrpName_seq_map = {_mirna[name_thrp] : name_thrp
                        for name_thrp in names_thrp}

    ## 3. MAKE OUTPUT DICTIONARIES #############################################
    # 3.i Make the output dictionary for instances with only the programmed
    # site.
    progOnly_map = {site : 0.0 for site in names_suppcomp}
    # 3.ii Make the output dictionary for when seed sites, but not three prime
    # sites, are in the random region.
    singleSeed_map = {
        prog_site : {
                                    #          Num. of columns       #
                                    #    Num. of library positions   #
            site : np.array([[0.0]*(n_lib + 2*n_con - int(site[0]) + 1)])
            for site in names_suppcomp
        }
        for prog_site in names_suppcomp
    }
    # 3.iii Make the output dictionary for when three prime sites, but not seed
    # sites, are found in the random region.
    ## DIMENSION OF EACH ARRAY: number
    singleThrP_map = {
        prog_site : {
            "%smer" %n_k : np.array(
                     #       Num. of columns       |       Num. of rows      #
                     #  Num. of library positions  | Num. of miRNA positions #
                [[0.0]*(n_lib + 2*n_con - n_k + 1)]*(len_mir3p - n_k + 1)
            )
            for n_k in range(kmer_min, kmer_max)
        }
        for prog_site in names_suppcomp
    }
    # 3.iv Make the output dictionary for when both seed sites and threeprime
    # sites are found in the random region. Note that this output dictionary
    # does not hold onto positional information with respect to where the two
    # sites are found within the library.
    thrPAndSeed_map = {
        prog_site : {
            name_seedAndMm : {
                name_thrp : 0.0 for name_thrp in names_thrp
            }
            for name_seedAndMm in names_suppcomp
        }
        for prog_site in names_suppcomp
    }
    ## 4. CREATE DIAGNOSTIC DICTIONARIES (NOT ESSENTIAL, BUT USEFUL). ##########
    # 4.i Make a counter for the five possible ways a read is counted:
    # 1. Programmed site only,
    # 2. programmed and random-position seed sites, 
    # 3. programmed and random-position three prime sites,
    # 4. programmed and both seed and three prime sites in the random-position
    #    region
    # 5. No programmed site in the programmed region (suggesting that the seed
    # is shifted with respect to the intended position)
    counts_outputType_map = {
        i : 0 for i in ["Programmed", "Seed", "Threeprime", "SeedAndThreePrime",
                        "None"]
     }
    # 4.ii Make a dictionary to count the average number of occurences of
    # threeprime sites per read.
    perReadFreq_len_map = {kmer : (0, 0) for kmer in range(kmer_min, kmer_max)}
    print_counter = 0.0
    # NEW FOR COUNTING THE LOOP SEQUENCES:
    # Make a dictionary that records all instances of loop sequences, enabling
    # later sequence based analysis of the enrichment of various sequences.
    loop_sequences = {
        mir_pos : {
            pair_len : {
                offset : [] for offset in range(6)
            }
            for pair_len in range(4, 10)
        }
        for mir_pos in range(11, 14)
    }
    #### New 
    ReporterSeq_map = {
        prog_site : defaultdict(int) for prog_site in names_comp
    }


    ## 5. ITERATE OVER THE READS ###############################################
    for i_r, r in enumerate(read_seqs):
        sys.stdout.flush()
        read = r.strip()
        _read = Read(read, n_lib, _mirna, n_con, experiment)
        ## 6. CHECK FOR PROGRAMMED, SEED, AND THREEPRIME SITES #################
        # 6.i.a Determine the site at the programmed site.
        name_prog = programmedName_seq_map.get(_read.seq[l_libprog:r_libprog])
        if name_prog:
            # 6.i.b Get the name and left-hand position of the programmed site.
            if name_prog == "6mer-A1":
                l_prog = l_libprog + 2
            elif name_prog in ["7mer-A1", "6mer"]:
                l_prog = l_libprog + 1
            else:
                l_prog = l_libprog
            namel_prog = (name_prog, l_prog)
            # 6.ii.a Get all seed in the read region.
            _read.get_all_sites_in_read(_sitelist)
            namelr_seeds_pre = [(site.name, site.l, site.r)
                                for site in _read.sites]
            # 6.ii.b Remove the seed site that is at the programmed region. 
            namelr_seeds = [i for i in namelr_seeds_pre if i[:2] != namel_prog]
            # 6.iii.a Look for 3p sites.
            [len_thrp, l_mirAndThrPs_pre] = LCSubString(mirna_3p_seq_rc,
                                                      _read.seq[:l_libprog])
            # 6.iii.b Make lengths into list of tuples rather than a tuple of
            # lists.
            l_mirAndThrPs_pre = [(i, j) for (i, j) in zip(l_mirAndThrPs_pre[0], l_mirAndThrPs_pre[1])]
            # 6.iii.c Pre-allocate the final threeprime site object.
            l_mirAndThrPs = []
            # 6.iii.d If the length of the threeprime site is within the range.
            if len_thrp >= kmer_min and len_thrp < kmer_max:
                # 6.iii.e Iterate each of the potential three prime sites.
                for i, l_mirAndThrP in enumerate(l_mirAndThrPs_pre):
                    l_thrp = l_mirAndThrP[1]
                    r_thrp = l_thrp + len_thrp
                    # 6.iii.f Check that both of its ends are more than two
                    # nucleotides away from the programmed region.
                    if (abs(l_thrp - l_libprog) <= 2 or abs(r_thrp - r_libprog) <= 2):
                        print("THIS DOES MATTER")
                        sys.stdout.flush()
                        return()
                    if abs(l_thrp - l_libprog) > 2 and abs(r_thrp - r_libprog) > 2:
                        add_thrp = True
                        for (name_seed, l_seed, r_seed) in namelr_seeds:
                            # 6.iii.g Makes sure it isn't contained within any seed
                            # site.
                            if l_thrp >= l_seed and r_thrp <= r_seed:
                                add_thrp = False
                        # 6.iii.h If it passes both tests, add to threeprime list.
                        if add_thrp:
                            l_mirAndThrPs.append((l_mirAndThrP))

            ## 7. ADD ANY PROGRAMMED, SEED, OR THREEPRIME SITES TO THEIR #######
            ## RESPECTIVE OUTPUT DICTIONARIES. #################################
            # 7.i A Check for programmed site only.
            if ("TATACAACCAATTTTCAACCTCA" in read or
                "TATACAACCGATTATCCACCTCA" in read or
                "TATACAACCAAAAATCTACCACA" in read or
                "TATACAACCGATAAGCTCCCTCA" in read or
                "TATACAACCAATATACTACCTGA" in read):
                print("READ HAPPENS")
                print(len(namelr_seeds))
                print(len(l_mirAndThrPs))
                print(namelr_seeds)
                print(l_mirAndThrPs)
            if len(namelr_seeds) + len(l_mirAndThrPs) == 0:
                progOnly_map[name_prog] += 1
                counts_outputType_map["Programmed"] += 1
            # 7.ii Check for programmed and seed site only.
            elif len(l_mirAndThrPs) == 0:
                count = 1./len(namelr_seeds)
                for (name_seed, l_seed, r_seed) in namelr_seeds:
                    singleSeed_map[name_prog][name_seed][0, l_seed] += count
                counts_outputType_map["Seed"] += 1
            # 7.iii Check for programmed and threeprime site only.
            elif len(namelr_seeds) == 0:
                count = 1./len(l_mirAndThrPs)

                # 7.iii.b Update the preReadFreq object used for diagnostic
                # purposes.
                (len_ave, len_num) = perReadFreq_len_map[len_thrp]
                # Calculate the new weighted average.
                len_ave_new = (len_ave*len_num + len(l_mirAndThrPs))/(len_num + 1.)
                perReadFreq_len_map[len_thrp] = (len_ave_new, len_num + 1.)

                n_thrp = "%smer" %len_thrp
                for (l_mirThrP, l_thrP) in l_mirAndThrPs:
                    # Define the left and right indeces with respect to miRNA position 1
                    r_mir = len_mir - l_mirThrP
                    l_mir = r_mir - len_thrp + 1
                    # Define the left and right indeces with respect to library
                    # programmed site 1.
                    r_lib = 25 + n_con + 8 - l_thrP
                    l_lib = r_lib - len_thrp + 1
                    name_thrP = "%s-m%s.%s" %(n_thrp, l_mir, r_mir)
                    offset = l_lib - l_mir
                    name_check = "%s|%s|%s" %(name_thrP, l_lib, name_prog)
                    # print("###################################")
                    # print("name_check: %s" %name_check)
                    # print(_read.seq)
                    # print(' '*l_thrP + "_"*len_thrp + " "*(25 + n_con - l_thrP - len_thrp) + "_"*8)
                    # print(mirna_3p_seq_rc)
                    # print(" "*l_mirThrP + "_"*len_thrp)
                    # sys.stdout.flush()

                    # This boolean is in order to exclude sites that are separated
                    # due to their inclusion in the let-7a reporter assay experiments.
                    add_to_normal_dict = True
                    # ADDED CODE WHEN RECORDING LOOP SEQUENCES.
                    # Added code: Boolean for analyses looking at loop length
                    # for 3p-sites prior to design of the high-throughput
                    # reporter assay.
                    # Going to record loop sequences for pairing of lengths
                    # 4–9 bp in length, offsets of 0–5 nt, and with pairing
                    # beginning at miRNA position 11.
                    if (11 <= l_mir <= 13 and
                        4 <= len_thrp <= 9 and
                        0 <= offset <= 5):
                        loop_sequence = _read.seq[l_thrP + len_thrp:25 + n_con]
                        loop_sequences[l_mir][len_thrp][offset].append(loop_sequence)
                    # Added code when assigning Kds to the sites used directly
                    # in the reporter library experiment.
                    if ((l_mir == 11 and len_thrp == 9 or
                         l_mir == 13 and len_thrp in [4, 9]) and
                        offset in [1, 4] and name_prog in names_comp):
                        if i_r == 192791:
                            print(name_prog)
                            print(name_prog in names_comp)
                            print(((l_mir == 11 and len_thrp == 9) or
                        (l_mir == 13 and len_thrp in [4, 9]) and
                        offset in [1, 4] and
                        name_prog in names_comp))
                            sys.stdout.flush()
                        # Extract the loop sequence, which begins at the left-
                        # hand position of the library index + the length of
                        # pairing. It extends to the beginning of the programmed
                        # site, which is at a defined position (25 + the number
                        # of constant nucleotides appended to the read).
                        loop_sequence = _read.seq[l_thrP + len_thrp:25 + n_con]
                        # Trim to just get the region equivalent to the region
                        # that was inserted in the let-7a reporter assay:
                        #         NNNAAAAN|
                        #            xxxx |
                        #            1234 |    <- offset = 4
                        #           -5  -1|
                        added_loop = loop_sequence[-(offset + 1):-1]
                        if added_loop in ["A", "T", "AAAA", "TTTT"]:
                            name_reporter = "%s|%s|%s" %(name_thrP, l_lib, added_loop)
                            # print(name_reporter)
                            if name_prog == "6mer-A1":
                                print(i_r)
                                print(name_reporter)
                                print(read)
                                sys.stdout.flush()
                            ReporterSeq_map[name_prog][name_reporter] += 1
                            add_to_normal_dict = False
                        elif len(added_loop) == 4 and len_thrp == 9:
                            if l_mir == 13:
                                print(loop_sequence)
                                sys.stdout.flush()
                            all_AU_bool = len(
                                set(added_loop).difference(set("AT"))
                            ) == 0
                            if all_AU_bool:
                                print("_______________")
                                print(loop_sequence)
                                print(added_loop)
                                sys.stdout.flush()
                                name_reporter = "%s|%s|%s" %(name_thrP, l_lib, "WWWW")
                                if name_prog == "6mer-A1":
                                    print(i_r)
                                    print(name_reporter)
                                    print(read)
                                    print(names_comp)
                                    print(name_prog)
                                    sys.stdout.flush()
                                ReporterSeq_map[name_prog][name_reporter] += 1
                                add_to_normal_dict = False
                        sys.stdout.flush()
                    if "10mer-m13.22|13|8mer-mm" in name_check:
                        print(name_check + "\t" + str(i_r))
                        print_counter += count
                        print(print_counter)
                        sys.stdout.flush()
                    if add_to_normal_dict:
                        singleThrP_map[name_prog][n_thrp][l_mirThrP, l_thrP] += count
                counts_outputType_map["Threeprime"] += 1
            # 2.D Check for programmed and seed, and three prime site only. ####
            else:
                count = 1./(len(namelr_seeds) * len(l_mirAndThrPs))
                for (l_mirRCThrP, l_thrp) in l_mirAndThrPs:
                    # Converting the position along the reverse-complemented
                    # Mirna threeprime sequence to the position along the miRNA
                    # Example: pos = 1, len = 6, let-7a-21nt
                    # 2         1
                    #10987654321098765431           13 + 8
                    # 1                              - 1
                    # 123456                         - 6
                    #      |                         + 1
                    #     15
                          # miRNA len #   # rc pos #   # site 
                    p_l = len_mir3p + 8 - l_mirRCThrP - len_thrp + 1 
                    p_r = p_l + len_thrp - 1
                    name_thrP = "%smer-m%s.%s" %(len_thrp, p_l, p_r)
                    for (name_seed, l_seed, r_seed) in namelr_seeds:
                        thrPAndSeed_map[name_prog][name_seed][name_thrP] += count
                counts_outputType_map["SeedAndThreePrime"] += 1
        else:
            counts_outputType_map["None"] += 1
    print("got to the end of the loop")
    print("Output-type tallies:")
    for (key, value) in counts_outputType_map.items():
        print("%s: %1.0f%%" %(key, value/float(i_r + 1)*100))
    print(ReporterSeq_map)
    for key in ReporterSeq_map.keys():
        print(ReporterSeq_map[key])
    sys.stdout.flush()
    return (progOnly_map, singleSeed_map, singleThrP_map, thrPAndSeed_map,
            print_counter, ReporterSeq_map, loop_sequences)

def main():
    # 1. PARSE ARGUMENTS AND INITIALIZE VARIABLES ##############################
    time_start = time.time()
    arguments = [
        "miRNA", "experiment", "condition", "n_constant", "-jobs",
        "-test_binary"
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
    # Make the nput files.
    input_list = [(mirna, experiment, condition)]
    if not jobs:
        jobs = 20
    # Define the arguments to be given to the iterated function.
    args  = [_mirna, experiment, int(n_constant), rand_length]
    ############################################################################
    ## 2. PERFORM MMULTIPROCESSING OF THE READ FILE OR FILES (DEPENDING ON
    ## WHETHER THE COMBINED INPUT IS BEING USED).
    # 2.i Initialize the threads.
    threads = []
    # 2.ii For each file in the input list, perform the multiprocessing.
    for (mirna, experiment, condition) in input_list:
        reads_path = get_analysis_path(mirna, experiment, condition, "reads")
        print(reads_path)
        threads += multiproc_file(reads_path, int(jobs), assign_site, test,
                                  *args)
    ## 3. COLLECT THE THE THREADS AND MAKE THE OUTPUT FILES.
    # 3.i Initialized each summed map using the first thread.
    progOnly_map, singleSeed_map, singleThrP_map, thrPAndSeed_map, counter_all, loop_sequences, reporter_seq_map = threads[0]
    # 3.ii Get the keys required to iterate over the rest of the threads:
    prog_keys = list(progOnly_map.keys())
    seed_keys = list(singleSeed_map[prog_keys[0]].keys())
    thrL_keys = list(singleThrP_map[prog_keys[0]].keys())
    thrP_keys = list(thrPAndSeed_map[prog_keys[0]][seed_keys[0]].keys())
    # 3.iii Iterate over the rest of the threads to add their dictionaries to
    # the initialized one.
    for t in threads[1:]:
        counter_all += t[-3]
    print(counter_all)
    for t in threads[1:]:
        # Get the programmed key for the first dictionaries.
        for pk in prog_keys:
            progOnly_map[pk] += t[0][pk]
            # Get the seed key for the second and fourth dictionaries.
            for sk in seed_keys:
                singleSeed_map[pk][sk] += t[1][pk][sk]
                # Get the threeP site key for the fourth dictionary.
                for tpk in thrP_keys:
                    thrPAndSeed_map[pk][sk][tpk] += t[3][pk][sk][tpk]
            # Get the threeP length key for the third dictionary.
            for tkl in thrL_keys:
                singleThrP_map[pk][tkl] += t[2][pk][tkl]
    ## 7. TRANSFORM THESE FOUR DICTIONARIES INTO THE OUTPUT DATAFRAMES. ########
    # 7.i. Initialize the final dataframe for both the full site counts and the
    # collapsed suppcomp counts.
    counts = pd.DataFrame.from_dict(progOnly_map, orient="index")
    # counts_collapsed = counts.copy()
    counts_suppcomp = counts.copy()
    # 7.ii Make a dictionary that maps each of the 24 possible programmed sites to
    # either "Supp" or "Comp".
    progSite_suppCompSite_map = {key : val for (key, val) in
                                 zip(prog_keys, ["Supp"]*6 + ["Comp"]*18)}

    # counts.index = ["NA|NA|%s" %i for i in counts.index]
    progAndSeed_map = {}
    progAndSeed_suppcomp_map = defaultdict(int)
    # 2. Add the counts from each of the 16 seed-in-random-region dictionaries.
    for pk in prog_keys:
        sck = progSite_suppCompSite_map[pk]
        for sk in seed_keys:
            df_i = singleSeed_map[pk][sk]
            len_s = int(sk.split("mer")[0])
            # Gives the position of index 0 in the dataframe with respect to the
            # miRNA position 1.
            # Example, n_constant = 1, 8mer seed site, at position 0
            # 1                                           1     n_constant
            #  1234567890123456789012345               + 25     N
            # 12345678                                 -  8     len seed
            #                           87654321       +  8
            #        7654321098765432109               +  1
            #        2      2         1                = 27                        
            p0 = 25 + int(n_constant) - len_s + 8 + 1
            for p in range(df_i.shape[1]):
                v = df_i[0, p]
                mir_p = p0 - p
                if mir_p < 1:
                    mir_p -= 1
                prefix = "%s|%s|" %(sk, mir_p)
                progAndSeed_map[prefix + pk] = v
                progAndSeed_suppcomp_map[prefix + sck] += v
    counts = counts.append(
        pd.DataFrame.from_dict(progAndSeed_map, orient="index")
    )
    counts_suppcomp = counts_suppcomp.append(
        pd.DataFrame.from_dict(progAndSeed_suppcomp_map, orient="index")
    )

    print("added df 2")
    # 3. Add the counts from each of the 16X5 threep-in-random-region
    # dictionaries.
    progAndThrP_map = {}
    # progAndThrP_collapsed_map = {}
    progAndThrP_suppcomp_map = defaultdict(int)

    for pk in prog_keys:
        sck = progSite_suppCompSite_map[pk]
        for tlk in thrL_keys:
            len_thrP = int(tlk.split("mer")[0])
            df_i = singleThrP_map[pk][tlk]
            # Get the dimensions of the 2D array.
            n_mir, n_lib = df_i.shape
            # Flatten the array into one list.
            vs = [j for i in df_i.tolist() for j in i]
            # Gives the position of index 0 in the dataframe with respect to the
            # miRNA position 1.
            # Example, n_constant = 1, 8mer seed site, at position 0
            #     2         1
            #    109876543210987654321  l               r
            # 0  1234567890   X        12 (8 + 4 - 0)  21 (8 + 4 + 10 -1 - 0)                             1     n_constant
            # 1   1234567890  X        11 (8 + 4 - 1)  20 (8 + 4 + 10 -1 - 1)
            # 2    1234567890 X        10 (8 + 4 - 2)  19 (8 + 4 + 10 -1 - 2)
            # 3     1234567890X         9 (8 + 4 - 3)  18 (8 + 4 + 10 -1 - 3)
            # The list was flattened maintaining the ordering of the rows,
            # therefore the inner list comprehension iterates over the columns,
            # i.e., the positions within the library.
            off_l = n_mir + 8
            off_r = off_l + len_thrP - 1
            off_mirp = 25 + int(n_constant) - len_thrP + 8 + 1
            name_thrPs = ["%s-m%s.%s" %(tlk, off_l - i, off_r - i)
                          for i in range(n_mir)]
            prefix = ["%s|%s|" %(n, off_mirp - j)
                      for n in name_thrPs
                      for j in range(n_lib)]
            for pre, z in zip(prefix, vs):
                progAndThrP_map[pre + pk] = z
                # progAndThrP_collapsed_map[name_collapsed] += count
                progAndThrP_suppcomp_map[pre + sck] += z
    print("collapsed 3 length:")
    counts = counts.append(
        pd.DataFrame.from_dict(progAndThrP_map, orient="index")
    )
    counts_suppcomp = counts_suppcomp.append(
        pd.DataFrame.from_dict(progAndThrP_suppcomp_map, orient="index")
    )

    print("added df 3")
    progAndSeedAndThrP_map = {}
    progAndSeedAndThrP_suppcomp_map = defaultdict(int)

    for pk in prog_keys:
        sck = progSite_suppCompSite_map[pk]
        for sk in seed_keys:
            for tpk in thrP_keys:
                v = thrPAndSeed_map[pk][sk][tpk]
                prefix = "%s&%s|NA|" %(sk, tpk)
                progAndSeedAndThrP_map[prefix + pk] = v
                progAndSeedAndThrP_suppcomp_map[prefix + sck] += v
    print("collapsed 4 length:")
    print(len(progAndSeedAndThrP_map))
    counts = counts.append(
        pd.DataFrame.from_dict(progAndSeedAndThrP_map, orient="index")
    )
    counts_suppcomp = counts_suppcomp.append(
        pd.DataFrame.from_dict(progAndSeedAndThrP_suppcomp_map, orient="index")
    )
    print("added df 4")
    ##
    # Assign for the two different count files being made.
    extension = "_%s_progthrp" %n_constant
    extension_suppcomp = "%s_suppcomp" %extension
    extension_save = extension
    if test:
        extension = "%s_test" %extension
        extension_suppcomp = "%s_test" %extension_suppcomp

    # The names of the output file:
    site_counts_path = get_analysis_path(
        mirna, experiment, condition, "site_counts", ext=extension + "_new2"
    )
    site_counts_suppcomp_path = get_analysis_path(
        mirna, experiment, condition, "site_counts",
        ext=extension_suppcomp + "_new2"
    )
    print(site_counts_path)
    print(site_counts_suppcomp_path)

    # Write the output files.
    # Looking at loop lengths from 4 to 9
    # positions 11 to 13
    # offsets of 0 to +5.
    for l_mir in range(11, 14):
        for len_thrp in range(4, 10):
            for offset in range(0, 6):
                extension = extension_save + "_pos%s_len%s_off%s" %(l_mir, len_thrp, offset)
                if test:
                    extension = "%s_test" %extension
                loop_sequences_path = get_analysis_path(
                    mirna, experiment, condition, "loop_sequences",
                    ext=extension + "_new2"
                )
                print(loop_sequences_path)
                loops = []
                for t in threads:
                    loops += t[-1][l_mir][len_thrp][offset]
                with open(loop_sequences_path, "w+") as file:
                    file.write("\n".join(loops))


    # NEW PORTION WHERE THE SCRIPT MAKES THE REPORTER SEQUENCES.
    counts_reporterSite_map = defaultdict(int)
    countsSuppComp_reporterSite_map = defaultdict(int)

    for t in threads:
        thrpSiteDict_progSite_map = t[-2]
        prog_sites = thrpSiteDict_progSite_map.keys()
        print(prog_sites)
        for prog_site in prog_sites:
            counts_thrpSite_map = thrpSiteDict_progSite_map[prog_site]
            for thrp_site in counts_thrpSite_map.keys():
                print(thrp_site)
                counts_site = counts_thrpSite_map[thrp_site]
                print(counts_site)
                counts_reporterSite_map[thrp_site + "|" + prog_site] += counts_site
                countsSuppComp_reporterSite_map[thrp_site + "|Comp"] += counts_site

    print(counts_reporterSite_map)
    print(countsSuppComp_reporterSite_map)

    print(counts.iloc[-30:, :])
    print(counts_suppcomp.iloc[-30:, :])


    counts = counts.append(
        pd.DataFrame.from_dict(counts_reporterSite_map, orient="index")
    )
    counts_suppcomp = counts_suppcomp.append(
        pd.DataFrame.from_dict(countsSuppComp_reporterSite_map, orient="index")
    )

    print(counts.iloc[-30:, :])
    print(counts_suppcomp.iloc[-30:, :])
    print(site_counts_path)
    print(site_counts_suppcomp_path)

    # Output the dataframes to the two paths.
    counts.to_csv(site_counts_path, sep="\t", header=False)
    counts_suppcomp.to_csv(site_counts_suppcomp_path, sep="\t", header=False)
    print_time_elapsed(time_start)
################################################################################

if __name__ == "__main__":
    main()

