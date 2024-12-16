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
def assign_site(read_seqs, _mirna, experiment, n_con, n_lib, start_mm, stop_mm,
                buffer3p):
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


    counts = defaultdict(int)
    for i_r, r in enumerate(read_seqs):
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

        ##################### INSERTION OF MISMATCH CODE #######################
        mm_tuples = [i for i in site_mm_list if i[0] in _read.seq]
        b_tuples = [i for i in site_b_list if i[0] in _read.seq]
        d_tuples = [i for i in site_d_list if i[0] in _read.seq]
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
        follow_read = False
        if len(mmbd_tuples) != 0:
            print_later = False
            mmbd_names = [i[1] for i in mmbd_tuples]
            mmbd_tuples = [(i[1], _read.seq.find(i[0])) for i
                           in mmbd_tuples]

            mmbd_starts = [i[1] for i in mmbd_tuples]
            mmbd_stops = [i[0] + len_mm + int("b" in i[1]) - int("d" in i[1])
                         for i in zip(mmbd_starts, mmbd_names)]
            mmbd_tuples_new = []
            for mmbd_name, mmbd_start, mmbd_stop in zip(mmbd_names, mmbd_starts, mmbd_stops):
                # First exclusion criterion: If either of the two nucleotides
                # on either side of the mismatch/bulge/deltion site are
                # complementary to the extended position of the miRNA, exclude
                # the site (because it is actually a longer type of site).
                bool_non_wc_edge = (
                    (mmbd_start == 0  or 
                     _read.seq[mmbd_start - 1] != tar_left) and
                    (mmbd_stop == len(_read.seq) or
                     _read.seq[mmbd_stop] != tar_right)
                )
                # Since it passes this criteria, now check whether or not it
                # overlaps with any seed sites.
                if bool_non_wc_edge:
                    bool_overlap_seed = False
                    if seed_names != ["None"]:
                        # print("Also seed sites:")
                        for seed_site in seed_sites:
                        # Define the ends of the seed_site.
                            s_l, s_r = (seed_site.l, seed_site.r)
                            # Check if the thrp site overlaps with the seed site.
                            s_l_int = s_l >= mmbd_start and s_l < mmbd_stop
                            s_r_int = s_r > mmbd_start and s_r <= mmbd_stop
                            t_l_int = mmbd_start >= s_l and mmbd_start < s_r
                            t_r_int = mmbd_stop > s_l and mmbd_stop <= s_r
                            bool_overlap_seed_i = bool(sum([s_l_int, s_r_int,
                                                        t_l_int, t_r_int]))
                            if bool_overlap_seed_i:
                                bool_overlap_seed = True
                    if not bool_overlap_seed:
                        mmbd_tuples_new.append((mmbd_name, mmbd_start))
            mmbd_tuples = mmbd_tuples_new

        bipartite_sites = []
        if seed_names != ["None"] and (thrp_names != ["None"] or mmbd_tuples): #______POTENTIAL BIPARTITE SITES_______|
            # First exclude any threep sites that overlap any of the seed_sites.                                     #|  
            thrp_sites_final = []                                                                                    #| 
            if thrp_names != ["None"]: #___________________CHECK 3p SITES FOR SEED OVERLAP________________________|  #|
                for thrp_site in thrp_sites:                                                                     #|  #|
                    # Define the ends of the thrp_site.                                                          #|  #|
                    add_seq = True                                                                               #|  #|
                    t_l, t_r  = (thrp_site.l, thrp_site.r)                                                       #|  #|
                    for seed_site in seed_sites:                                                                 #|  #|
                        # Define the ends of the seed_site.                                                      #|  #|
                        s_l, s_r = (seed_site.l, seed_site.r)                                                    #|  #|
                        # Check if the thrp site overlaps with the seed site.                                    #|  #|
                        s_l_int = s_l >= t_l and s_l < t_r                                                       #|  #|
                        s_r_int = s_r > t_l and s_r <= t_r                                                       #|  #|
                        t_l_int = t_l >= s_l and t_l < s_r                                                       #|  #|
                        t_r_int = t_r > s_l and t_r <= s_r                                                       #|  #|
                        # If any amount of overlap of the seed and thrp site:                                    #|  #|
                        if s_l_int or s_r_int or t_l_int or t_r_int:                                             #|  #|
                            # If the threeprime site overlaps the seed site from the                             #|  #|
                            # 5′ end:                                                                            #|  #|                             
                            if s_l_int and t_r_int:                                                              #|  #|
                                # Define the overlap:                                                            #|  #|
                                len_overlap = t_r - s_l                                                          #|  #|
                                # Get the length and range of pairing                                            #|  #|
                                len_thrp = int(thrp_site.name.split("mer")[0])                                   #|  #|
                                thrp_pairing = [                                                                 #|  #|
                                    int(i) for i                                                                 #|  #|
                                    in thrp_site.name.split("mer-m")[1].split(".")                               #|  #|
                                ]                                                                                #|  #|
                                # Update the length and range of pairing according                               #|  #|                          
                                # to the extent of overlap with the seed site.                                   #|  #|                      
                                len_thrp_new = len_thrp - len_overlap                                            #|  #|              
                                thrp_pairing[0] = thrp_pairing[0] + len_overlap                                  #|  #|                      
                                # If the remaining site is still at least 4 nt in                                #|  #|                          
                                # length, redefine the site as the shorter site.                                 #|  #|
                                if len_thrp_new >= 4:                                                            #|  #|
                                    thrp_name = "%smer-m%s.%s" %(                                                #|  #|
                                        len_thrp_new, thrp_pairing[0],                                           #|  #|
                                        thrp_pairing[1]                                                          #|  #|
                                    )                                                                            #|  #|
                                    thrp_site = Read.ReadSite(_sitelistthrp,                                     #|  #|
                                                              thrp_name, t_l, s_l)                               #|  #|
                                # Otherwise discard the site.                                                    #|  #|
                                else:                                                                            #|  #|
                                    add_seq = False                                                              #|  #|
                            # This else occurs if there is overlap but with the                                  #|  #|
                            # threeprime site overlapping the seed site from the                                 #|  #|
                            # 5'-end.                                                                            #|  #|
                            else:                                                                                #|  #|
                                add_seq = False                                                                  #|  #|
                    if add_seq:                                                                                  #|  #|
                        thrp_sites_final += [thrp_site] #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx|  #|
            thrp_sites = thrp_sites_final                                                                            #|
            check_read = False                                                                                       #|
            if len(thrp_sites) != 0: #______________THERE ARE 3P SITES (MIGHT BE MMBD SITES)______________________|  #|
                ######### EXTRA PART FOR MISMATCH SCRIPT #######################                                 #|  #|
                ##### Goal here is to check whether, if there are both 3p sites                                  #|  #|
                ##### and mismatch sites, to check if the mismatch site is of                                    #|  #|
                ##### greater length than the three-prime sites. If so, the                                      #|  #|
                ##### "threeprime site" is going to be the mismatch site.                                        #|  #|
                for seed_site in seed_sites:                                                                     #|  #|
                    s_l, s_r = (seed_site.l, seed_site.r)                                                        #|  #|
                    t_r = [thrp_site.r for thrp_site in thrp_sites]                                              #|  #|
                    # Just check the thrp_sites that are to the left of the                                      #|  #|
                    # seed site in question.                                                                     #|  #|
                    # IMPORTANT TO USE TEMP, BECAUSE A SEED SITE ALL THE WAY AT                                  #|  #|
                    # THE 5' END OF THE READ SHOULDN"T STOP A THRP SITE FROM                                     #|  #|
                    # FORMING A BIPARTITE SITE WITH ANOTHER SEED FURTHER 3'                                      #|  #|
                    # WITHIN THE READ.                                                                           #|  #|
                    thrp_sites_temp = [thrp_site for thrp_site in thrp_sites if thrp_site.r <= s_l]              #|  #|
                    # Get the max length thrp_site at this length.                                               #|  #|
                    thrp_sites_len = [int(thrp_site.name.split("mer")[0])                                        #|  #|
                                      for thrp_site in thrp_sites_temp]                                          #|  #|
                    thrp_sites_temp = [                                                                          #|  #|
                        thrp_site for i, thrp_site in enumerate(thrp_sites_temp)                                 #|  #|
                        if thrp_sites_len[i] == max(thrp_sites_len)                                              #|  #|
                    ]                                                                                            #|  #|
                    if mmbd_tuples and (len(thrp_sites_len) == 0 or len_mm > max(thrp_sites_len)): #_____|       #|  #|       There are either site  mmbd sites longer than 3p sites or there are no 3p sites
                        mmbd_r = [                                                                      #|       #|  #|
                            mmbd_l + len_mm + ("b" in mmbd_name)                                        #|       #|  #|
                            - ("d" in mmbd_name)                                                        #|       #|  #|
                            for (mmbd_name, mmbd_l) in mmbd_tuples                                      #|       #|  #|
                        ]                                                                               #|       #|  #|
                        mmbd_temp = [                                                                   #|       #|  #|       Check the mmbd sites for seed overlap
                            mmbd for (i_m, mmbd) in enumerate(mmbd_tuples)                              #|       #|  #|
                            if mmbd_r[i_m] <= s_l                                                       #|       #|  #|
                        ]                                                                               #|       #|  #|
                    else:                                                                               #|       #|  #|
                        mmbd_temp = [] #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx|       #|  #|
                    if mmbd_temp: #______________________________________________________________________|       #|  #|       There are mmbd sites that fully 5prime of a seed site  
                        for (mmbd_name, mmbd_l) in mmbd_temp:                                           #|       #|  #|             So add them
                            # Calculate the actual distance between the seed site and                   #|       #|  #|
                            # thrp site.                                                                #|       #|  #|
                            mmbd_r = mmbd_l + len_mm + ("b" in mmbd_name) - ("d" in mmbd_name)          #|       #|  #|
                            pos = seed_site.l - mmbd_r                                                  #|       #|  #|
                            pos += endPos_seedName_map[seed_site.name] + 1                              #|       #|  #|
                            full_site_name = "%s|%s|%s" %(mmbd_name, pos, seed_site.name)               #|       #|  #|
                            bipartite_sites += [full_site_name]                                         #|       #|  #|
                    else:                                                                               #|       #|  #|
                        for thrp_site in thrp_sites_temp:                                               #|       #|  #|       There aren't mbd sites that are fully 5 prime of a seed site
                            # Calculate the actual distance between the seed site and                   #|       #|  #|             so add threep contiguous site if they exist
                            # thrp site.                                                                #|       #|  #|
                            t_l, t_r = (thrp_site.l, thrp_site.r)                                       #|       #|  #|
                            pos = seed_site.l - thrp_site.r                                             #|       #|  #|
                            # Calculate the adjustment for where the seed site ends                     #|       #|  #|
                            # (i.e., The 7mer-A1 ends at position 7, not 8), and add                    #|       #|  #|
                            # 1 to give the position at which the 3′ site begins.                       #|       #|  #|
                            pos += endPos_seedName_map[seed_site.name] + 1                              #|       #|  #|
                            full_site_name = "%s|%s|%s" %(thrp_site.name, pos, seed_site.name)          #|       #|  #|
                            bipartite_sites += [full_site_name] #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx|       #|  #|
                if len(bipartite_sites) != 0: #_____________________MADE SOME BIPARTITE SITES_________________|  #|  #|
                    # Add all of the bipartite sites.                                                        #|  #|  #|
                    for bipartite_site in bipartite_sites:                                                   #|  #|  #|
                        counts[bipartite_site] += 1./len(bipartite_sites)                                    #|  #|  #|
                # Otherwise, add each of the seed and thrp sites to the count                                #|  #|  #|
                # dictionary, dividing the read by the total number of sites.                                #|  #|  #|
                else: #.........................................................NO BIPARTITE SITES............|  #|  #|
                    max_thrp_sites_len = max(                                                                #|  #|  #|
                        [int(thrp_site.name.split("mer")[0])                                                 #|  #|  #|
                         for thrp_site in thrp_sites]                                                        #|  #|  #|
                    )                                                                                        #|  #|  #|
                    if mmbd_tuples and len_mm > max_thrp_sites_len: #________ADD SEED AND MMBD SOLO_______|  #|  #|  #|
                        for seed_site in seed_sites:                                                     #|  #|  #|  #|
                            counts[seed_site.name] += 1./(len(seed_sites) + len(mmbd_tuples))            #|  #|  #|  #|
                        for mmbd_tuple in mmbd_tuples:                                                   #|  #|  #|  #|
                            counts[mmbd_tuple[0]] += 1./(len(seed_sites) + len(mmbd_tuples))             #|  #|  #|  #|
                    else: #.................................................ADD SEED AND THREEPRIME SOLO..|  #|  #|  #|
                        for seed_site in seed_sites:                                                     #|  #|  #|  #|
                            counts[seed_site.name] += 1./(len(seed_sites) + len(thrp_sites))             #|  #|  #|  #|
                        for thrp_site in thrp_sites:                                                     #|  #|  #|  #|
                            counts[thrp_site.name] += 1./(len(seed_sites) + len(thrp_sites)) #xxxxxxxxxxxx|xxx|  #|  #|
            # The condition in which none of the thrp sites survive the overlap                                  #|  #|
            # with seed site criteria, in which case there will only be seed                                     #|  #|
            # sites left. So, add the seed sites to the dictionary.                                              #|  #|
            elif mmbd_tuples: #.............................................SEED AND MMBD, NO THREEPRIME..........|  #|
                mmbd_r = [                                                                                       #|  #|
                    mmbd_l + len_mm + ("b" in mmbd_name) - ("d" in mmbd_name)                                    #|  #|
                    for (mmbd_name, mmbd_l) in mmbd_tuples                                                       #|  #|
                ]                                                                                                #|  #|
                for seed_site in seed_sites:                                                                     #|  #|
                    s_l, s_r = (seed_site.l, seed_site.r)                                                        #|  #|
                    # Just check the thrp_sites that are to the left of the                                      #|  #|
                    # seed site in question.                                                                     #|  #|
                    mmbd_temp = [                                                                                #|  #|
                        mmbd for (i_m, mmbd) in enumerate(mmbd_tuples)                                           #|  #|
                        if mmbd_r[i_m] <= s_l                                                                    #|  #|
                    ]                                                                                            #|  #|
                    for (mmbd_name, mmbd_l) in mmbd_temp:                                                        #|  #|
                        # Calculate the actual distance between the seed site and                                #|  #|
                        # thrp site.                                                                             #|  #|
                        mmbd_r = mmbd_l + len_mm + ("b" in mmbd_name) - ("d" in mmbd_name)                       #|  #|
                        pos = seed_site.l - mmbd_r                                                               #|  #|
                        pos += endPos_seedName_map[seed_site.name] + 1                                           #|  #|
                        full_site_name = "%s|%s|%s" %(mmbd_name, pos, seed_site.name)                            #|  #|
                        bipartite_sites += [full_site_name]                                                      #|  #|
                if len(bipartite_sites) != 0: #_________________________SOME BIPARTITE (MMBD) SITES__________|   #|  #|
                    # Add all of the bipartite sites.                                                       #|   #|  #|
                    for bipartite_site in bipartite_sites:                                                  #|   #|  #|
                        counts[bipartite_site] += 1./len(bipartite_sites)                                   #|   #|  #|
                # Otherwise, add each of the seed and thrp sites to the count                               #|   #|  #|
                # dictionary, dividing the read by the total number of sites.                               #|   #|  #|
                else: #................................................NO BIPARTITE (MMBD) SITES.............|   #|  #|
                    count = 1./(len(seed_sites) + len(mmbd_tuples))                                         #|   #|  #|
                    for seed_site in seed_sites:                                                            #|   #|  #|
                        counts[seed_site.name] += count                                                     #|   #|  #|
                    for mmbd_tuple in mmbd_tuples:                                                          #|   #|  #|
                        counts[mmbd_tuple[0]] += count #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx|   #|  #|
            else: #..........................................................ONLY SEED SITES......................|  #|
                for seed_site in seed_sites:                                                                     #|  #|
                    counts[seed_site.name] += 1./len(seed_sites) #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx|  #|
        # Conditions in which wasn't potential seed and thrp sites to begin                                          #|
        # with, such that the bipartite site analysis portion of the loop was                                        #|
        # never entered. So, that means there are only either seed or thrp sites                                     #|
        # within the read, so check for each "None" condition and add whichever                                      #|
        # one is correct.                                                                                            #|
        elif seed_names != ["None"]: #.......................................ONLY SEED SITES..........................|
            # if follow_read:                                                                                        #|
            #     print("Only seed and maybe mmbd")                                                                  #|
            #     sys.stdout.flush()                                                                                 #|
            print("____________________")                                                                            #|
            print(i_r)                                                                                               #|
            print(seed_sites)                                                                                        #|
            print(mmbd_tuples)                                                                                       #|
            print(thrp_sites)                                                                                        #|
            print(thrp_sites[0].name)                                                                                #|
            sys.stdout.flush()                                                                                       #|
            if len(mmbd_tuples) != 0 or thrp_sites[0].name != "None":                                                #|
                print("ISSUE")                                                                                       #|
                sys.stdout.flush()                                                                                   #|
                return()                                                                                             #|
            for seed_site in seed_sites:                                                                             #|
                counts[seed_site.name] += 1./len(seed_sites)                                                         #|
        elif thrp_names != ["None"]: #...................................ONLY THREE PRIME (AND MAYBE MMBD)............|
            max_thrp_sites_len = max(                                                                                #|
                        [int(thrp_site.name.split("mer")[0])                                                         #|
                         for thrp_site in thrp_sites]                                                                #|
            )                                                                                                        #|
            if mmbd_tuples and len_mm > max_thrp_sites_len: #________________USE MMBD_____________________________|  #|
                for mmbd_tuple in mmbd_tuples:                                                                   #|  #|
                    counts[mmbd_tuple[0]] += 1./len(mmbd_tuples)                                                 #|  #|
            else: #........................................................USE THRP...............................|  #|
                for thrp_site in thrp_sites:                                                                     #|  #|
                    counts[thrp_site.name] += 1./len(thrp_sites) #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx|  #|
        elif mmbd_tuples: #........................................................ONLY MMBD..........................|                                                                                           #|
            for mmbd_tuple in mmbd_tuples:                                                                           #|
                counts[mmbd_tuple[0]] += 1./len(mmbd_tuples)                                                         #|
            sys.stdout.flush()                                                                                       #|
        else: #.................................................................NO SITE...............................|
            if follow_read:                                                                                          #|
                print("No site counted")                                                                             #|
                sys.stdout.flush()                                                                                   #|
            counts["None"] += 1 #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx|
        # Conditional check:#|
        if seed_names != ["None"] and (thrp_names == ["None"] or len(thrp_sites) == 0) and mmbd_tuples:
            print("line: %s" %i_r)
            # print("seed_names:")
            # print(seed_names)
            # print("thrp_names:")
            # print(thrp_names)
            # print("mmbd_tuples:")
            # print(mmbd_tuples)
            print("bipartite_sites:")
            print(bipartite_sites)
            if bipartite_sites:
                print("HAS BIPARTITE SITES")
            else:
                print("NO BIPARTITE SITES")
            print("__________________________________________________________")
            sys.stdout.flush()
    print("i_r: %s" %i_r)
    print("summed dict: %s" %sum(counts.values()))
    return(counts)

def main():
    time_start = time.time()
    arguments = [
        "miRNA", "experiment", "condition", "n_constant", "start_mm", "stop_mm",
        "-jobs", "-buffer3p_binary", "-test_binary"
    ]
    args = parse_arguments(arguments)
    (
        mirna, experiment, condition, n_constant, start_mm, stop_mm, jobs,
        buffer3p, test
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
    extension = "_%s_randthrp_m%s.%smmsd" %(n_constant, start_mm, stop_mm)
    extension_comp = "_%s_randthrp_comp_m%s.%smmsd" %(n_constant, start_mm, stop_mm)
    extension_suppcomp = "_%s_randthrp_suppcomp_m%s.%smmsd" %(n_constant, start_mm, stop_mm)

    if buffer3p:
        extension = "%s_buffer3p" %extension
        extension_comp = "%s_buffer3p" %extension_comp
        extension_suppcomp = "%s_buffer3p" %extension_suppcomp

    extension = "%s_new" %extension
    extension_comp = "%s_new" %extension_comp
    extension_suppcomp = "%s_new" %extension_suppcomp



    if test:
        extension = "%s_test" %extension
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

    args  = [_mirna, experiment, int(n_constant), rand_length, int(start_mm),
             int(stop_mm), buffer3p]
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

