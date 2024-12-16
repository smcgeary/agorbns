import operator
import itertools

def hamming_distance( str1, str2 ):
    """
    - Returns the Hamming distance of str1 and str2
    """
    assert len(str1) == len(str2)
    ne = operator.ne
    return sum(itertools.imap(ne, str1, str2))





def is_homopolymer( kmer ):
    """
    - Returns True or False, depending on if kmer is a homopolymer
    """
    nts_S = set( [x for x in kmer] )
    if (len( nts_S ) == 1):
        return True
    else:
        return False



def is_homopolymer_w_1_intervening( kmer ):
    """
    - Returns a dictionary depending on if kmer is a homopolymer with
        exactly 1 intervening nt
    """

    nts_S = set( [x for x in kmer] )
    if (len( nts_S ) == 2):
        num_by_nt_D = {}
        for nt in kmer:
            try:
                num_by_nt_D[nt] += 1
            except KeyError:
                num_by_nt_D[nt] = 1

        nt_numtimes_T_L = [(nt, num_by_nt_D[nt]) for nt in num_by_nt_D]
        nt_numtimes_T_L.sort( key = lambda x: -1*x[1] )
        if (nt_numtimes_T_L[0][1] == len( kmer ) - 1) and\
                (nt_numtimes_T_L[1][1] == 1):
            major_nt = nt_numtimes_T_L[0][0]
            minor_nt = nt_numtimes_T_L[1][0]
            pos_of_minor_nt = kmer.find( minor_nt )
            return_D = {"is_homopolymer_w_1_intervening": True,
                    "major_nt": major_nt,
                    "minor_nt": minor_nt,
                    "pos_of_minor_nt": pos_of_minor_nt}
        else:
            return_D = {"is_homopolymer_w_1_intervening": False}
    else:
        return_D = {"is_homopolymer_w_1_intervening": False}
    return return_D



def get_best_match_of_kmer_to_foundingkmer(
        kmer_to_align,
        founding_kmer,
        possible_comps_by_foundingk_alignk_D ):
    """
    - Will try to align the kmer_to_align to the founding_kmer (trying all
        possible sliding combinations), with the number of mismatches
        allowed specified by possible_comps_by_foundingk_alignk_D

    - Returns:
        return_D = {"best_match_offset": best_match_offset,
            "best_match_side": best_match_side,
            "best_match": best_match}
    """
    best_match_offset = None
    best_match_side = None
    best_match = None
    best_match_position_in_allowedmatchesL = 100
    ####    allowed_matches_L is like:
    ####        ["side1_mismatch0", "side0_mismatch1", "side2_mismatch0"],
    ####        where side is the number of unaligned (overhang) positions,
    ####        and mismatch is the # of mismatched among the aligned positions
    ####    - This list is ordered from best -> worst, so if multiple offsets
    ####        are in the list, we'll use the one that has the lowest
    ####        best_match_position_in_allowedmatchesL
    allowed_matches_L = possible_comps_by_foundingk_alignk_D\
            [len(founding_kmer)][len(kmer_to_align)]
    sides_allowed_L = [int(x.split("side")[-1][0]) for x in allowed_matches_L]

    for offset in range( -3, len( founding_kmer ) ):

        #### Get the # of nt hanging of the "side" of the founding kmer
        if (offset < 0):
            side = abs( offset )
            pos_to_align = len( kmer_to_align ) - side
            kmer_to_align_for_mismatches = kmer_to_align[(-1*pos_to_align):]
            founding_kmer_for_mismatches = founding_kmer[:pos_to_align]

        else:
            side = max( offset + len( kmer_to_align ) - len( founding_kmer ), 0 )

            if (side == 0):
                kmer_to_align_for_mismatches = kmer_to_align
            else:
                kmer_to_align_for_mismatches = kmer_to_align[:-1*side]

            founding_kmer_for_mismatches = founding_kmer[
                    offset:offset+len( kmer_to_align_for_mismatches )]

        if side not in sides_allowed_L:
            continue

        #### Get the # of mismatches between the kmer_to_align_for_mismatches
        ####    and founding_kmer_for_mismatches
        mismatches = hamming_distance(
                kmer_to_align_for_mismatches,
                founding_kmer_for_mismatches )
        #### See if this side/mismatch combination is allowed
        this_match = "side{0}_mismatch{1}".format( side, mismatches )

        if (this_match in allowed_matches_L):
            pos_in_allowedmatchesL = allowed_matches_L.index( this_match )
            #### If this offset is better than (i.e., has a lover index in
            ####    allowed_matches_L) the previous best one, record it
            if (pos_in_allowedmatchesL < best_match_position_in_allowedmatchesL):
                best_match_position_in_allowedmatchesL = pos_in_allowedmatchesL
                best_match_offset = offset
                best_match_side = side
                best_match = this_match

    return_D = {"best_match_offset": best_match_offset,
            "best_match_side": best_match_side,
            "best_match": best_match}

    return return_D











def add_kmer_to_existing_motif_class_or_start_new(
        kmer,
        excess_R,
        kmer_offset_R_T_Ls_by_classnum_D ):
    """
    - For a kmer & its excess_R, will try to:
        1. Add it to an existing motif_class (in kmer_offset_R_T_Ls_by_classnum_D)
        2. Start a new motif_class if that's not possible

    - Returns:        return_D =
        {"kmer_offset_R_T_Ls_by_classnum_D": kmer_offset_R_T_Ls_by_classnum_D,
         "motif_class_added_to": motif_class_added_to,
         "offset_used": offset_used}

    See possible_comps_by_foundingk_alignk_D below for the different
        alignment/mismatch combinations that allow a kmer to be added to a logo
        based on its relationship to the founding kmer of an existing
        motif_class
    """
    possible_comps_by_foundingk_alignk_D = {
            4: {4: ["side1_mismatch0", "side0_mismatch1", "side2_mismatch0"]},

            5: {5: ["side1_mismatch0", "side0_mismatch1", "side2_mismatch0",
                    "side1_mismatch1"],
                4: ["side1_mismatch0", "side0_mismatch1", ]},

            6: {6: ["side1_mismatch0", "side0_mismatch1", "side2_mismatch0",
                    "side1_mismatch1", "side0_mismatch2", "side2_mismatch1"],
                5: ["side1_mismatch0", "side0_mismatch1", "side1_mismatch1",
                    "side2_mismatch0"],
                4: ["side1_mismatch0", "side0_mismatch1"]},

            7: {7: ["side1_mismatch0", "side2_mismatch0", "side0_mismatch1",
                "side1_mismatch1", "side3_mismatch0", "side0_mismatch2",
                "side2_mismatch1", "side1_mismatch2"],
                6: ["side1_mismatch0", "side0_mismatch1", "side2_mismatch0",
                    "side1_mismatch1", "side2_mismatch1", "side3_mismatch0",
                    "side0_mismatch2"],
                5: ["side1_mismatch0", "side0_mismatch1", "side1_mismatch1",
                    "side2_mismatch0"],
                4: ["side1_mismatch0", "side0_mismatch1"]},
            8: {8: ["side0_mismatch1",
                    "side0_mismatch2",
                    "side1_mismatch0",
                    "side1_mismatch1",
                    "side1_mismatch2",
                    "side2_mismatch0",
                    "side2_mismatch1",
                    "side3_mismatch0",
                    "side3_mismatch1",
                    "side4_mismatch0"],
                7: ["side1_mismatch0",
                    "side2_mismatch0",
                    "side0_mismatch1",
                    "side1_mismatch1",
                    "side3_mismatch0",
                    "side0_mismatch2",
                    "side2_mismatch1",
                    "side1_mismatch2"],
                6: ["side1_mismatch0", "side0_mismatch1", "side2_mismatch0",
                    "side1_mismatch1", "side2_mismatch1", "side3_mismatch0",
                    "side0_mismatch2"],
                5: ["side1_mismatch0", "side0_mismatch1", "side1_mismatch1",
                    "side2_mismatch0"],
                4: ["side1_mismatch0", "side0_mismatch1"]} }

    score_D = {"match": 1,
        "mismatch": -1,
        "side": -0.7}

    name_score_T_Ls_by_alignk_D = {
        8: [],
        7: [],
        6: [],
        5: [],
        4: []}

    #### Go through and score all of the possibilities in
    ####    possible_comps_by_foundingk_alignk_D
    for founding_k, founding_k_D in possible_comps_by_foundingk_alignk_D.iteritems():
        for align_k, L in founding_k_D.iteritems():
            for alignment in L:
                side = int( alignment.split("side")[-1][0] )
                mismatch = int( alignment[-1] )
                match = align_k - side - mismatch

                score = match*score_D["match"] +\
                    mismatch * score_D["mismatch"] +\
                    side * score_D["side"]
                name = "{0}_{1}".format( founding_k, alignment )
                name_score_T_Ls_by_alignk_D[align_k].append(
                        (name, score) )

    #### Go through all of the align_ks, sort by descending score
    pref_order_Ls_by_alignk_D = {}
    for align_k, name_score_T_L in name_score_T_Ls_by_alignk_D.iteritems():
        name_score_T_L.sort( key = lambda x: -1* x[1] )
        pref_order_Ls_by_alignk_D[align_k] = [x[0] for x in name_score_T_L]



    existing_motif_class_nums_L = range( len( kmer_offset_R_T_Ls_by_classnum_D ) )
    #### the length of the kmer to be aligned
    align_k = len( kmer )

    #### Go through each of the existing motif classes &
    matches_to_existing_L = []
    pos_of_best_match_in_orderedL = 100
    best_existing_motif_class = None
    best_offset = None
    for existing_motif_class in existing_motif_class_nums_L:
        founding_kmer = kmer_offset_R_T_Ls_by_classnum_D[existing_motif_class][0][0]
        print("founding kmer: %s" %(founding_kmer))

        #### If the founding_kmer is a homopolymer, deal with it differently
        if ( is_homopolymer( founding_kmer ) == True ):

            #### Determine if this is a homopolymer off-by-1; if so, include
            ####    it with the minor_nt at the center position
            homo_off_by_1_returned_D = is_homopolymer_w_1_intervening(
                    kmer )

            print(homo_off_by_1_returned_D)
            print(homo_off_by_1_returned_D["is_homopolymer_w_1_intervening"] == False)

            if (homo_off_by_1_returned_D["is_homopolymer_w_1_intervening"] == False):
                continue
            print("past first conditional")

            #### Make sure that the homopolymer nt is the same between the
            ####    founding_kmer and kmer
            if (founding_kmer[0] != homo_off_by_1_returned_D["major_nt"]):
                continue
            print("here in loop")
            #### If it's made it this far, kmer is a homopolymer off-by-1; if
            ####    NO add'l kmers have been added to this motif class, add
            ####    this kmer to have offset 0
            if (len( kmer_offset_R_T_Ls_by_classnum_D[existing_motif_class] ) == 1):

                pos_of_best_match_in_orderedL = -1
                best_existing_motif_class = existing_motif_class
                best_offset = 0

            #### Otherwise, if an off-by-1 homopolymer has already been added
            ####    to this motif class, match up the non-homopolymer nt
            ####    positions
            else:
                non_homo_already_added_kmer = kmer_offset_R_T_Ls_by_classnum_D[existing_motif_class][1][0]
                num_nt_before_non_homo_of_already_added = is_homopolymer_w_1_intervening(
                        non_homo_already_added_kmer )["pos_of_minor_nt"]
                num_nt_before_non_homo_this_kmer = is_homopolymer_w_1_intervening(
                        kmer )["pos_of_minor_nt"]
                offset_to_align_nonhomo = num_nt_before_non_homo_of_already_added -\
                        num_nt_before_non_homo_this_kmer

                pos_of_best_match_in_orderedL = -1
                best_existing_motif_class = existing_motif_class
                best_offset = offset_to_align_nonhomo

        #### If the founding_kmer is NOT a homopolymer, do the normal alignment
        else:
            #### best_match_returned_D has keys:
            ####    best_match_offset (is None if there's no match)
            ####    best_match_side (the number of unaligned positions; None if
            ####        theres no match)
            ####    best_match: like "side0_mismatch1"
            best_match_returned_D = get_best_match_of_kmer_to_foundingkmer(
                    kmer,
                    founding_kmer,
                    possible_comps_by_foundingk_alignk_D )
            print(best_match_returned_D)
            if type( best_match_returned_D["best_match_offset"] ) is int:
                foundingk_match_name = "{0}_{1}".format(
                        len(founding_kmer), best_match_returned_D["best_match"] )
                try:
                    pos_of_this_match_in_orderedL = pref_order_Ls_by_alignk_D[
                            align_k].index( foundingk_match_name )
                    if (pos_of_this_match_in_orderedL < pos_of_best_match_in_orderedL):
                        pos_of_best_match_in_orderedL = pos_of_this_match_in_orderedL
                        best_existing_motif_class = existing_motif_class
                        best_offset = best_match_returned_D["best_match_offset"]
                except ValueError:
                    pass

    print("best_existing_motif_class: %s" %(best_existing_motif_class))
    #### If the kmer aligned to any existing motif class
    if (type( best_existing_motif_class ) is int):
        kmer_offset_R_T_Ls_by_classnum_D[best_existing_motif_class].append(
                ( kmer, best_offset, excess_R ) )
        offset_used = best_offset
        motif_class_added_to = best_existing_motif_class
    #### Otherwise, start a new class
    else:
        motif_class_added_to = len( kmer_offset_R_T_Ls_by_classnum_D )
        kmer_offset_R_T_Ls_by_classnum_D[ motif_class_added_to ] =\
                [(kmer, 0, excess_R)]
        offset_used = 0

    return_D = {"kmer_offset_R_T_Ls_by_classnum_D":
            kmer_offset_R_T_Ls_by_classnum_D,
        "motif_class_added_to": motif_class_added_to,
            "offset_used": offset_used}
    return return_D


test_data = [("UAUAU", 4.05),
             ("UUUUU", 3.12),
             ("UUUAU", 2.82),
             ("AUUAU", 2.38),
             ("UAAAU", 2.22),
             ("UUAAU", 2.05),
             ("UAUGU", 1.98),
             ("UAUUU", 1.86)]



output_dict = {}

for kmer, R in test_data[:3]:
    print("_________" + kmer)
    print(R)


    output = add_kmer_to_existing_motif_class_or_start_new(kmer, R, output_dict)
    output_dict = output["kmer_offset_R_T_Ls_by_classnum_D"]
    print(output_dict)
