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
def assign_site(read_seqs, _mirna, experiment, n_con, n_lib, start_mm, stop_mm, new):
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

    ## SET UP THE DICTIONARIES REQUIRED FOR CONVERTING THE THREEPRIME ##########
    ## SEQUENCES TO THEIR NAMES, FOR PUTTING THEM INTO THE OUTPUT ##############
    ## DICTIONARIES. ###########################################################
    # Make the list of all possible programmed sites in each read.
    seq_8mer = _mirna["8mer"]
    seq_8mer_list = list(seq_8mer)
    names_can = ["8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-m8", "6mer-A1"]
    _sitelist = SiteList(_mirna, "programmed", n_lib + 2*n_con)
    ###### BEGIN THE PORTION SPECIFIC TO THE MISMATCH CODE (NOT FOUND IN THE ###
    ###### BIPARTITESITESPROGRAMMEDLIBRARY.PY CODE) ############################
    # Get the name of the base site, as well as its sequence.
    len_mm = stop_mm - start_mm + 1
    name_3p_site_mms = "%smer-m%s.%s" %(len_mm, start_mm, stop_mm)
    mm_seq = _mirna[name_3p_site_mms]
    # Assign the nucleotides of a target that is one nucleotide longer on either
    # side than the site for which the mismatches, bulges, and deletions are
    # being counted. These will be used to make sure that the site that is found
    # is not actually longer on either side than intended (this is required
    # because I am not looking for the *longest* comlementary stretch for each
    # mismatch type, I am rather looking for each defined site, without respect
    # to the end of the site).
    if start_mm != 9:
        tar_right = get_rc(_mirna.seq[start_mm - 2])
    else:
        tar_right = ""
    if stop_mm != len(_mirna):
        tar_left = get_rc(_mirna.seq[stop_mm])
    else:
        tar_left = ""
    # Build the mismatch dictionary to be used in addition to the seed, ########
    # threeprime, and seed+threeprime dictionaries for counting purposes. The ##
    # if conditional statement for each of the three site types is to make sure
    # that the sequence of the mismatch, bulge, or deletion site is not itself
    # complementary to another region of the miRNA. This will inherently not ###
    # apply to any sites that are deletions on either side of the site, since
    # it will just be one nucleotide shorter as a result.
    site_mm_list = [i for i in get_all_mismatches(mm_seq, internal_only=False)
                    if i[0] not in get_rc(_mirna.seq)]
    site_b_list = [i for i in get_all_bulges(mm_seq)
                   if i[0] not in get_rc(_mirna.seq)]
    site_d_list = [i for i in get_all_deletions(mm_seq)
                   if i[0] not in get_rc(_mirna.seq)]
    # Convert the tuple entries in each of the mismatch, bulge, and deletion ###
    # lists into a name, which requires different string formating for each type
    # of site.
    site_mm_names = ["%smm%s%s" %(name_3p_site_mms, name[3], stop_mm - name[1])
                     for name in site_mm_list]
    site_b_names = ["%s%s%s" %(name_3p_site_mms, name[2], stop_mm - name[1] + 1)
                    for name in site_b_list]
    site_d_names = ["%sd%s" %(name_3p_site_mms, stop_mm - name[1])
                    for name in site_d_list]
    # 2.A Deal with the reduncancy of the bulge sites, caused by the
    # bulging of nucleotides adjacent to complementary nucleotides.
    site_b_seqs = [i[0] for i in site_b_list]
    # Make a dictionary where the keys are the sequences, and the values are all
    # of the names of sites that have that sequence.
    seq_site_map = defaultdict(list)
    for (name_b, seq_b) in zip(site_b_names, site_b_seqs):
        seq_site_map[seq_b] += [name_b]
    # Make a dictionary where the keys are now names corrected by using the
    # values of the prior dictionary.
    siteb_seq_map = defaultdict(list)
    for item in seq_site_map.items():
        sys.stdout.flush()
        (key, value) = item
        if len(value) > 1:
            # Here is where the name e.g. "9mer-m11.19bT13.15" is constructed
            # from the list of [ "9mer-m11.19b13", "9mer-m11.19b14", 
            # "9mer-m11.19b15"]. 
            nums = [int(i.split("b")[1][1:]) for i in value]
            min_str = str(min(nums))
            max_str = str(max(nums))
            base = value[0].split("b")[0]
            nuc = value[0].split("b")[1][0]
            name = "%sb%s(%s.%s)" %(base, nuc, min_str, max_str)
        else:
            name = value[0]
            sys.stdout.flush()
        siteb_seq_map[name] = key
    # 2.B Deal with the reduncancy of the deletion sites, caused by the
    # removal of nucleotides within a streth of identical nucleotides.
    site_d_seqs = [i[0] for i in site_d_list]
    # Make a dictionary where the keys are the sequences, and the values are all
    # of the names of sites that have that sequence.
    seq_site_map = defaultdict(list)
    for (name_d, seq_d) in zip(site_d_names, site_d_seqs):
        seq_site_map[seq_d] += [name_d]
    # Make a dictionary where the keys are now names corrected by using the
    # values of the prior dictionary.
    sited_seq_map = defaultdict(list)
    for item in seq_site_map.items():
        (key, value) = item
        if len(value) > 1:
            # Here is where the name e.g. "9mer-m11.19bT13.15" is constructed
            # from the list of [ "9mer-m11.19b13", "9mer-m11.19b14", 
            # "9mer-m11.19b15"]. 
            nums = [int(i.split("d")[1]) for i in value]
            min_str = str(min(nums))
            max_str = str(max(nums))
            base = value[0].split("d")[0]
            name = "%sd(%s.%s)" %(base, min_str, max_str)
        else:
            name = value[0]
        sited_seq_map[name] = key
    # 3. Make the final list of mismatch and bulge site, where each element is a
    # tuple, where the first entry is the sequence, and the second is the name.
    # Use these names for the final application of the site to the dictionary.
    # Keep these lists split up to allow for first querying the bulge sites
    # before next querying the mismatch sites, since the bulge sites are one
    # nucleotide longer, which means that if a bulge site is found it should
    # take precedence over a mismatch site.
    site_mm_list = list(zip([i[0] for i in site_mm_list], site_mm_names))
    site_b_list = list(i[::-1] for i in siteb_seq_map.items())
    site_d_list = list(i[::-1] for i in sited_seq_map.items())
    # Make the concatenated list of names to use in the output dictionary.
    site_mmbd_names = [i[1] for i in site_mm_list + site_b_list + site_d_list]

    ###### HERE IT RETURNS TO LOOKING LIKE THE BIPARTITE SITE SCRIPT ###########
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
    seqs_can = [seq_8mer] + seqs_7m8 + seqs_7A1 + seqs_6 + seqs_6m8 + seqs_6A1
    seqs_comp = [seq_8mer[:i] + j + seq_8mer[i+1:] for i in range(1, 7)
                 for j in DNTS if j != seq_8mer[i]]

    # seqs_8merMm = [seq_8mer[:i] + j + seq_8mer[i+1:]
    #                for i in range(1, 7)
    #                for j in ["A", "C", "G", "T"] if j != seq_8mer[i]]
    # 2.i.b Define the names the canonical sites and single mismatch sites.
    names_can = (
        ["8mer"] + 3*["7mer-m8"] + 3*["7mer-A1"] +
        9*["6mer"] + 9*["6mer-m8"] + 9*["6mer-A1"]
    )
    # Get the names of the 8mer 1-nt mismatch sites.
    names_comp = ["8mer-mm%s%s" %(j, 8 - i)
                    for i in range(1, 7)
                    for j in DNTS if j != seq_8mer[i]]
    # Assin the names and sequences of all supplemental and compensatory sites
    # in the library (which may or may not all be included in the library, but
    # are expected to be there at least at a low frequency due to 1-nt
    # mismatches 
    names_suppcomp = names_can + names_comp
    seqs_suppcomp = seqs_can + seqs_comp
    # 2.i.c Make the programmed site dictionary.
    programmedName_seq_map = {seq_suppcomp : name_suppcomp
                              for (seq_suppcomp, name_suppcomp)
                              in zip(seqs_suppcomp, names_suppcomp)}
    # Generate a list of the three prime site names, site sequences, and a
    # dictionary with them as values and keys, respectively.
    names_thrp = ["%smer-m%s.%s" %(k, p, p + k - 1)
                  for k in range(kmer_min, kmer_max)
                  for p in range(9, len(_mirna) - k + 2)]
    # 2.ii.b Make the threprime site dictionary.
    thrpName_seq_map = {_mirna[name_thrp] : name_thrp
                        for name_thrp in names_thrp}

    ## MAKE THE OUTPUT DICTIONARIES ############################################
    # Make the output dictionary for instances with only the programmed site.
    progOnly_map = {site : 0.0 for site in names_suppcomp}
    # 3.ii Make the output dictionary for when seed sites, but not three prime
    # sites, are in the random region.
    singleSeed_map = {
        prog_site : {              # Lenth of read # # len site #  # correction#
            site : np.array([[0.0]*(n_lib + 2*n_con - int(site[0]) + 1)])
            for site in names_suppcomp
        }
        for prog_site in names_suppcomp
    }
    ## THIS IS DIFFERENT THAN THE ASSIGNBIPARTITESITESPROGRAMMEDLIB SCRIPT #####
    # 3.ii Make the output dictionary for when seed sites, but not three prime
    # sites, are in the random region.
    singleThrPMm_map = {
        prog_site : {              # Lenth of read # # len site #  # correction#
            site : np.array([[0.0]*(n_lib + 2*n_con - int(site[0]) + 1)])
            for site in site_mmbd_names
        }
        for prog_site in names_suppcomp
    }
    ####### BACK TO NORMAL FOR BIPARTITESITESPROGRAMMEDLIB SCRIPT ##############
    # 3.iii Make the output dictionary for when three prime sites, but not seed
    # sites, are found in the random region.
    singleThrP_map = {
        prog_site : {
            "%smer" %n_k : np.array(
                        # Lenth of read #
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
            name_seed : {
                name_thrp : 0.0 for name_thrp in names_thrp
            }
            for name_seed in names_suppcomp
        }
        for prog_site in names_suppcomp
    }
    # 3.v Make the output dictionary for when both seed sites and threeprime
    # sites are found in the random region. Note that this output dictionary
    # does not hold onto positional information with respect to where the two
    # sites are found within the library.
    thrPMmAndSeed_map = {
        prog_site : {
            name_seed : {
                name_mmbd : 0.0 for name_mmbd in site_mmbd_names
            }
            for name_seed in names_suppcomp
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
        i : 0 for i in ["Programmed", "Seed", "ThreePrime", "ThreePrimeMmbd",
                        "SeedAndThreePrime", "SeedAndThreePrimeMmbd", "None"]
    }
    # 4.ii Make a dictionary to count the average number of occurences of
    # threeprime sites per read.
    perReadFreq_len_map = {kmer : (0, 0) for kmer in range(kmer_min, kmer_max)}
    ## 5. ITERATE OVER THE READS ###############################################
    for i_r, r in enumerate(read_seqs):
        read = r.strip()
        _read = Read(read, n_lib, _mirna, n_con, experiment)
        # Grab the sequence of the programmed site.
        # prog_kmer = _read.seq[l_libprog:r_libprog]
        # `allow_thrp` is a boolean that I might remove, because I cant' come up
        # with a good reason for using it. Right now it is used to not allow any
        # three prime sites to be counted if a mismatch site has been found that
        # extends beyond the limits of the intended site in the script.
        allow_thrp = True
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
            ######### HERE THE SCRIPT BREAKS WITH THE BIPARTITE RANDOM SCRIPT
            ####################################################################
            # 1.C Look for all instances of 3p sites in read:
            mm_tuples = [i for i in site_mm_list if i[0] in _read.seq[:l_libprog]]
            b_tuples = [i for i in site_b_list if i[0] in _read.seq[:l_libprog]]
            d_tuples = [i for i in site_d_list if i[0] in _read.seq[:l_libprog]]
            # if i_r in [13659, 189094, 228260, 270515, 370380, 616595, 672184, 754111, 793677, 814271, 900607, 979929]:
            #     print(mm_tuples)
            #     print(b_tuples)
            #     print(d_tuples)
            #     sys.stdout.flush()
            # The conditions before the "else" statement check whether only
            # one category of mismatch site is present, in which case checking
            # for either deletions within mismatch sites, or mismatch and
            # deletion sites with the bulges is not necessary.
            l_mm, l_b, l_d = [len(i) for i in [mm_tuples, b_tuples, d_tuples]]
            if l_mm + l_b + l_d == 0:
                mmbd_tuples = []
            elif l_b + l_d == 0:
                mmbd_tuples = mm_tuples
            elif l_d + l_mm == 0:
                mmbd_tuples = b_tuples
            elif l_b + l_mm == 0:
                mmbd_tuples = d_tuples
            else:
                # If this statement is reached, that that means at least two
                # categories of mmbd type sites are present, so the categories
                # have to be checked for smaller sites within larger sites.
                # Remove any deletion sites found within the mismatch sites.
                d_seqs = [i[0] for i in d_tuples]
                mm_seqs = [i[0] for i in mm_tuples]
                d_tuples_pass = [i for i in d_tuples if np.sum([i[0] in j for j in mm_seqs]) == 0]
                # Remove any mismatch or deletion sites that are found within
                # the bulge sites.
                b_seqs = [i[0] for i in b_tuples]
                mmd_tuples_pass = [i for i in d_tuples_pass + mm_tuples if np.sum([i[0] in j for j in b_seqs]) == 0]
                mmbd_tuples = b_tuples + mmd_tuples_pass
            # if i_r in [13659, 189094, 228260, 270515, 370380, 616595, 672184, 754111, 793677, 814271, 900607, 979929]:
            #     print("____________________________")
            #     print(mmbd_tuples)
            #     sys.stdout.flush()
            if len(mmbd_tuples) != 0:
                mmbd_tuples_new = []
                mmbd_names = [i[1] for i in mmbd_tuples]
                mmbd_tuples = [(i[1], _read.seq[:l_libprog].find(i[0])) for i
                               in mmbd_tuples]

                mmbd_starts = [i[1] for i in mmbd_tuples]
                mmbd_stops = [i[0] + len_mm + int("b" in i[1]) - int("d" in i[1])
                             for i in zip(mmbd_starts, mmbd_names)]
                bs_global = False
                for mmbd_name, mmbd_start, mmbd_stop in zip(mmbd_names, mmbd_starts, mmbd_stops):
                    ## 200829 I'm not sure why I was looking for the "edge" type
                    ## beforehand. It seems that even if the mismatch is on the
                    ## edge, I would want to still make sure that the mismatch
                    ## site should not have complementarity to the miRNA at the
                    ## positions yet further out. Also, it isn't obvious to me
                    ## that only the mismatch sites should be queried– shouldn't
                    ## the bulge and deletion sites *also* require that there is
                    ## no pairing on either side of the mismatch site? I will
                    ## check this and name everything "new"
                    # Check if the mm/b/d site found is a misatch site.
                    ## 201009 I am checking this formally by making a "new" and a non
                    ## new versin of this analysis, where the `new` analysis will
                    ## allow pairing to a yet further external nucleotide if it is
                    ## a mismatch site and it's on the end, where as the normal
                    ## version of the analysis will not allow any sites to have
                    ## complementarity to either next-nucleotide out, no matter
                    ## whether it is a bulged-site, a deletion-site, an edge-mismatch
                    ## site, or an internal-mismatch site.
                    mms = re.findall(r'(?:mm)(?:A|C|G|T)(?:\(\d+\.\d+\)|\d+)', mmbd_name)

                    if mms:
                        mm_pos = int(mms[0][3:])
                        bool_edge_l = (mm_pos == start_mm)
                        bool_edge_r = (mm_pos == stop_mm)
                    else:
                        bool_edge_l = False
                        bool_edge_r = False
                    if new:
                        bool_internal = (
                            (mmbd_start == 0  or bool_edge_r or
                             _read.seq[mmbd_start - 1] != tar_left) and
                            (mmbd_stop == l_libprog or bool_edge_l or
                             _read.seq[mmbd_stop] != tar_right)
                        )
                    else:
                        bool_internal = (
                            (mmbd_start == 0  or 
                             _read.seq[mmbd_start - 1] != tar_left) and
                            (mmbd_stop == l_libprog or
                             _read.seq[mmbd_stop] != tar_right)
                        )

                    # if mms:
                    #     print("_______________________________")
                    #     print(mmbd_name)
                    #     print(i_r)
                    #     print(_read.seq[:l_libprog])
                    #     print(" "*mmbd_start + _mirna[mmbd_name])
                    #     print(" "*mmbd_start + _mirna[mmbd_name.split("mm")[0]])
                    #     print(" "*(mmbd_start - 1) + tar_left)
                    #     print(" "*(mmbd_stop) + tar_right)
                    #     print(_mirna["21mer-m1.21"])
                    #     print(bool_internal)
                    #     print(bool_internal_alt)
                    #     if i_r == 99657:
                    #         print(_mirna["10mer-m11.20bA20"])
                    #     sys.stdout.flush()
                    if bool_internal:
                        mmbd_tuples_new.append((mmbd_name, mmbd_start))
                    else:
                        allow_thrp = False
                mmbd_tuples = mmbd_tuples_new
            ##### THIS RETURNS TO ORIGINAL BIPARTITESITE VERSION OF THE SCRIPT #
            ####################################################################
             # 6.iii.a Look for 3p sites.
            [len_thrp, l_mirAndThrPs_pre] = LCSubString(mirna_3p_seq_rc,
                                                      _read.seq[:l_libprog])
            # 6.iii.b Make lengths in to list of tuples rather than a tuple of
            # lists.
            l_mirAndThrPs_pre = [(i, j) for (i, j) in zip(l_mirAndThrPs_pre[0], l_mirAndThrPs_pre[1])]
            # 6.iii.c Pre-allocate the final threeprime site object.
            l_mirAndThrPs = []
            # 6.iii.d If the length of the threeprime site is within the range.
            # Loop to get rid of three prime sites that are hovering over the
            # programmed region.

            for i, l_mirAndThrP in enumerate(l_mirAndThrPs_pre):
                l_thrp = l_mirAndThrP[1]
                r_thrp = l_thrp + len_thrp
                # 6.iii.f Check that both of its ends are more than two
                # nucleotides away from the programmed region.
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

            # If the three prime pairing is longer than the mismatch/bulge don't
            # count the mismatch/bulge..
            # If the threep match exists and is longer than the min length:
            if len_thrp >= kmer_min and len(l_mirAndThrPs) != 0 and allow_thrp:
                # If the threep match is longer than the mismatch length, remove
                # the mismatch sites, because there is a longer fully
                # complementary site within.
                if len_thrp >= len_mm:
                    mmbd_tuples = []
                    # If the match is too long, don't count the site. (This is
                    # internal to the outer loop because in either case the
                    # mismatch site shouldn't be kept, because the read is
                    # enriched not because of the mismatch site).
                    if len_thrp >= kmer_max:
                        l_mirAndThrPs = []
                # If here, that means that the three-prime site is either
                # shorter than the mismatch site, or that it doesn't exist. If
                # the mismatch site exists, it should take precedence:
                elif len(mmbd_tuples) != 0:
                    l_mirAndThrPs = []
            # Here the three-prime site is either too short or it doesn't exist, so
            # assign it as `[]` in either case.
            else:
                l_mirAndThrPs = []
            # # Count the threeprime pairing if it is within the bounds hasn't
            # # been removed (via the p_mir_rc condition).
            # if (len_thrp >= kmer_min and len(p_mir_rc) != 0 and
            #     len_thrp < kmer_max):
            #     # If either there is no threep site or if the length of the
            #     # match to the threep is equal to or longer than the mismatch
            #     # site, don't consider that mismatch site.
            #     if len_thrp >= len_mm or len(mmbd_tuples) == 0:
            #         l_mirAndThrPs = list(zip(p_mir_rc, p_lib))
            #     else:
            #         l_mirAndThrPs = []
            # else:
            #     l_mirAndThrPs = []
            # if len(mmbd_tuples) != 0:
            #     if len_thrp >= 4 and allow_thrp:
            #         numThrPSites_count_map[str(len(l_mirAndThrPs))] += 1
            #         lenThrPSites_count_map[str(len_thrp)] += 1
            #     else:
            #         numThrPSites_count_map["0"] += 1

            # 2. Add the counts of this read to the appropriate output #########
            # dictionary, depending on whether there are seed+mm sites, threep #
            # sites, both, or neither. #########################################
            # 2.A Check for programmed site only. ##############################
            if len(l_mirAndThrPs) > 0 and len(mmbd_tuples) > 0:
                print(l_mirAndThrPs)
                print(mmbd_tuples)
                sys.stdout.flush()
                return()
            if len(namelr_seeds) + len(l_mirAndThrPs) + len(mmbd_tuples) == 0:
                progOnly_map[name_prog] += 1
                counts_outputType_map["Programmed"] += 1
            # 2.B Check for programmed and seed site only. #####################
            elif len(l_mirAndThrPs) + len(mmbd_tuples) == 0:
                count = 1./len(namelr_seeds)
                name_sites = [i[0] for i in namelr_seeds]
                for (name_seed, l_seed, r_seed) in namelr_seeds:
                    # if print_line:
                    #     print(" ".join([name_prog, name_seed, str(l_seed)]))
                    #     sys.stdout.flush()
                    singleSeed_map[name_prog][name_seed][0, l_seed] += count
                counts_outputType_map["Seed"] += 1
            # 2.C Check for programmed and three prime site only. ##############
            elif len(namelr_seeds) + len(mmbd_tuples) == 0:
                count = 1./len(l_mirAndThrPs)
                # 7.iii.b Update the preReadFreq object used for diagnostic
                # purposes.
                (len_ave, len_num) = perReadFreq_len_map[len_thrp]
                # Calculate the new weighted average.
                len_ave_new = (len_ave*len_num + len(l_mirAndThrPs))/(len_num + 1.)
                perReadFreq_len_map[len_thrp] = (len_ave_new, len_num + 1.)
                for (l_mirThrP, l_thrP) in l_mirAndThrPs:
                    n_thrp = "%smer" %len_thrp
                    singleThrP_map[name_prog][n_thrp][l_mirThrP, l_thrP] += count
                counts_outputType_map["ThreePrime"] += 1
            # 2.D Check for programmed and seed, and three prime site only. ####
            elif len(namelr_seeds) == 0 and len(l_mirAndThrPs) == 0:
                count = 1./len(mmbd_tuples)
                for (mmbd_name, mmbd_pos) in mmbd_tuples:
                    singleThrPMm_map[name_prog][mmbd_name][0, mmbd_pos] += count
                counts_outputType_map["ThreePrimeMmbd"] += 1
            elif len(l_mirAndThrPs) == 0:
                count = 1./(len(namelr_seeds) * len(mmbd_tuples))
                for (mmbd_name, mmbd_pos) in mmbd_tuples:
                    for (name_seed, l_seed, r_seed) in namelr_seeds:
                        thrPMmAndSeed_map[name_prog][name_seed][mmbd_name] += count
                counts_outputType_map["SeedAndThreePrimeMmbd"] += 1
            else:
                count = 1./(len(namelr_seeds) * len(l_mirAndThrPs))
                for (l_mirRCThrP, l_thrp) in l_mirAndThrPs:
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
    sys.stdout.flush()

    return (progOnly_map, singleSeed_map, singleThrP_map, singleThrPMm_map,
            thrPAndSeed_map, thrPMmAndSeed_map)

def main():
    time_start = time.time()
    arguments = [
        "miRNA", "experiment", "condition", "n_constant", "start_mm", "stop_mm",
        "-jobs", "-new_binary", "-test_binary"
    ]
    args = parse_arguments(arguments)
    (
        mirna, experiment, condition, n_constant, start_mm, stop_mm, jobs, new, test
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
    extension = "_%s_progthrp_m%s.%smmsd" %(n_constant, start_mm, stop_mm)
    extension_suppcomp = "_%s_progthrp_suppcomp_m%s.%smmsd" %(n_constant, start_mm, stop_mm)
    if test:
        extension = "%s_test" %extension
        extension_suppcomp = "%s_test" %extension_suppcomp
    if new:
        extension = "%s_new" %extension
        extension_suppcomp = "%s_new" %extension_suppcomp
    input_list = [(mirna, experiment, condition)]
    if not jobs:
        jobs = 20
    args  = [_mirna, experiment, int(n_constant), rand_length, int(start_mm),
             int(stop_mm), new]
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

    (progOnly_map, singleSeed_map, singleThrP_map,
     singleThrPMm_map, thrPAndSeed_map, thrPMmAndSeed_map) = threads[0]
    # Use the first thread to get all of the dictionary keys (this is probably
    # bad coding).
    # Get the mismatch keys:
    prog_keys = list(progOnly_map.keys())
    # Get the seed and mismatch keys:
    seed_keys = list(singleSeed_map[prog_keys[0]].keys())
    thrL_keys = list(singleThrP_map[prog_keys[0]].keys())
    thrPMm_keys = list(singleThrPMm_map[prog_keys[0]].keys())
    thrP_keys = list(thrPAndSeed_map[prog_keys[0]][seed_keys[0]].keys())

    thrPMm_keys_alt = list(thrPMmAndSeed_map[prog_keys[0]][seed_keys[0]].keys())
    ## ADD UP ALL THE COUNTS FROM THE THREADS FOR ALL FOUR DICTIONARIES ########
    for thread in threads[1:]:
        # Get the programmed key for the first dictionaries.
        for prog_key in prog_keys:
            val = thread[0][prog_key]
            progOnly_map[prog_key] += val
            # Get the seed key for the second and fifth dictionaries.
            for seed_key in seed_keys:
                val = thread[1][prog_key][seed_key]
                singleSeed_map[prog_key][seed_key] += val
                # Get the threeP site key for the fifth dictionary.
                for thrP_key in thrP_keys:
                    val = thread[4][prog_key][seed_key][thrP_key]
                    thrPAndSeed_map[prog_key][seed_key][thrP_key] += val
                for thrPMm_key in thrPMm_keys:
                    val = thread[5][prog_key][seed_key][thrPMm_key]
                    thrPMmAndSeed_map[prog_key][seed_key][thrPMm_key] += val
            # Get the threeP length key for the third dictionary.
            for thrL_key in thrL_keys:
                val = thread[2][prog_key][thrL_key]
                singleThrP_map[prog_key][thrL_key] += val
            # Get the threePMM name key for the fourth dictionary.
            for thrPMm_key in thrPMm_keys:
                val = thread[3][prog_key][thrPMm_key]
                singleThrPMm_map[prog_key][thrPMm_key] += val


    ## MAKE THE FINAL SINGLE COUNT WHERE THE VARIOUS KEYS AND DIMENSIONS OF ####
    ## FOUR OUTPUT DICTIONARIES ARE TRANSFORMED INTO SITE NAMES ################
    # 1. Initialize the final dataframe with the counts of reads with only
    # programmed mismatch sites.
    counts = pd.DataFrame.from_dict(progOnly_map, orient="index")
    counts_suppcomp = counts.copy()
    base_sites = progOnly_map.keys()
    # Make a dictionary that maps each of the 24 possible programmed site to
    # either "Supp" or "Comp".
    progSite_suppCompSite_map = {base_site : collapsed_site for
                                 (base_site, collapsed_site) in
                                 zip(base_sites, ["Supp"]*6 + ["Comp"]*18)}

    # counts.index = ["NA|NA|%s" %i for i in counts.index]
    progAndSeed_map = {}
    progAndSeed_suppcomp_map = defaultdict(int)


    # 2. Add the counts from each of the 16 seed-in-random-region dictionaries.
    for prog_key in prog_keys:
        base_key = progSite_suppCompSite_map[prog_key]
        for seed_key in seed_keys:
            # collapsed_val = 0
            # key_collapsed = "%s|%s" %(seed_key, prog_key)
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
                # collapsed_val += val
            # progAndSeed_collapsed_map[key_collapsed] = collapsed_val
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
            # names_collapsed = ["%s|%s" %(base, prog_key)
            #                    for base in name_base for j in range(n_lib)]
            names_suppcomp = ["%s|%s|%s" %(base, off_mirp - j, base_key)
                               for base in name_base for j in range(n_lib)]
            # Iterate over zipped tuples of the names and counts to add them
            for name, name_suppcomp, count in zip(names, names_suppcomp,
                               array_list):
                progAndThrP_map[name] = count
                progAndThrP_suppcomp_map[name_suppcomp] += count
    print("collapsed 3 length:")
    counts = counts.append(
        pd.DataFrame.from_dict(progAndThrP_map, orient="index")
    )
    counts_suppcomp = counts_suppcomp.append(
        pd.DataFrame.from_dict(progAndThrP_suppcomp_map, orient="index")
    )
    print("added df 3")


    progAndThrPMm_map = {}
    progAndThrPMm_suppcomp_map = defaultdict(int)
    for prog_key in prog_keys:
        base_key = progSite_suppCompSite_map[prog_key]
        for thrPMm_key in thrPMm_keys:
            df_i = singleThrPMm_map[prog_key][thrPMm_key]
            len_site = int(thrPMm_key.split("mer")[0])
            mir_offset = 25 + int(n_constant) - len_site + 8 + 1 - int("b" in thrPMm_key) + int("d" in thrPMm_key)
            for j in range(df_i.shape[1]):
                site_pos = -j + mir_offset
                if site_pos < 1:
                    site_pos -= 1
                key = "%s|%s|%s" %(thrPMm_key, site_pos, prog_key)
                key_suppcomp = "%s|%s|%s" %(thrPMm_key, site_pos, base_key)
                val = df_i[0, j]
                progAndThrPMm_map[key] = val
                progAndThrPMm_suppcomp_map[key_suppcomp] += val
    counts = counts.append(
        pd.DataFrame.from_dict(progAndThrPMm_map, orient="index")
    )
    counts_suppcomp = counts_suppcomp.append(
        pd.DataFrame.from_dict(progAndThrPMm_suppcomp_map, orient="index")
    )
    print("added df 4")

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
    counts = counts.append(
        pd.DataFrame.from_dict(progAndSeedAndThrP_map, orient="index")
    )
    counts_suppcomp = counts_suppcomp.append(
        pd.DataFrame.from_dict(progAndSeedAndThrP_suppcomp_map, orient="index")
    )
    print("added df 5")


    progAndSeedAndThrPMm_map = {}
    progAndSeedAndThrPMm_suppcomp_map = defaultdict(int)
    for prog_key in prog_keys:
        base_key = progSite_suppCompSite_map[prog_key]
        for seed_key in seed_keys:
            for thrPMm_key in thrPMm_keys:
                name = "%s&%s|NA|%s" %(seed_key, thrPMm_key, prog_key)
                suppcomp_name = "%s&%s|NA|%s" %(seed_key, thrPMm_key, base_key)
                count = thrPMmAndSeed_map[prog_key][seed_key][thrPMm_key]
                progAndSeedAndThrPMm_map[name] = count
                progAndSeedAndThrPMm_suppcomp_map[suppcomp_name] += count
    print("collapsed 4 length:")
    counts = counts.append(
        pd.DataFrame.from_dict(progAndSeedAndThrPMm_map, orient="index")
    )
    counts_suppcomp = counts_suppcomp.append(
        pd.DataFrame.from_dict(progAndSeedAndThrPMm_suppcomp_map, orient="index")
    )
    print("added df 6")



    ############### Make the output files ######################################
    # The names of the output file:
    site_counts_path = get_analysis_path(
        mirna, experiment, condition, "site_counts", ext=extension
    )
    site_counts_suppcomp_path = get_analysis_path(
        mirna, experiment, condition, "site_counts",
        ext=extension_suppcomp
    )
    print(site_counts_path)
    print(site_counts_suppcomp_path)
    # Output the dataframes to the three paths.
    counts.to_csv(site_counts_path, sep="\t", header=False)
    counts_suppcomp.to_csv(site_counts_suppcomp_path, sep="\t", header=False)
    print_time_elapsed(time_start)
################################################################################

if __name__ == "__main__":
    main()

