def calc_3p_score(pos_l, pos_r, offset):
    ## Note, this function does not use python indexing, such that pairng to
    ## positions 13 to 16 would be put in with the argument (13, 16, 0).

    # First get the bone of the score for overlap with the 3p region.
    score_supp = max(0, min(pos_r, 16) - max(pos_l, 13) + 1)/2.
    # Next calculate the portion of the score that is agnostic to whether it
    # pairs to nucleotides 13–16 (i.e., just the complementarity).
    score_rest = (pos_r - pos_l + 1)/2.
    # Calculate the offset penalty, but not allowing it to be less than 1,
    # since the offset penalty should never *help* the pairing.
    offset_penalty = max((abs(offset) - 2)/2., 0)
    # Add up the score, but do not allow the score to be less than zero.
    score_tot = max(score_rest + score_supp - offset_penalty, 0)
    # print("%s\t%s\t%s\t%s" %(pos_l, pos_r, offset, score_tot))
    return(score_tot)



def calculate_threep_score(df_row, feature_map, mirna_seq, upstream_limit):
    """
    Calculate the three-prime pairing score
    Parameters
    NOTE: This is a modified version of Kathy's function, parameterized in order
    to directly compare with her output tables.
    ----------
    df_row: a row from a dataframe taken generated from one of the 12mer files.
    feature_map: a nexted dictionary containing the ORF and UTR sequences,
        as well as their lengths (even though this is redundant by having the
        sequences)
    mirna_seq: string: the sequence of the miRNA.
    Output
    ------
    float: 3' pairing score
    """
    # First get the transcript name.
    transcript_use = df_row["transcript"]
    # Assign the sequence of the ORF-UTR.
    orfutr_seq = feature_map["orf_utr3"][transcript_use]
    # Get the position of the 12mer, by taking the ORF length, the location of
    # 12mer relative to the start of the 3'UTR (this is how it is given in 
    # Kathy's file), and converting this in to where wthin the ORF-UTR string
    # this location corresponds to.
    len_orf = int(feature_map["orf_length"][transcript_use])
    str_loc = int(df_row["utr3_loc"])
    site_start = len_orf - 3 + str_loc

    # utr = orf_utr_sequence
    # upstream_limit = 15
    if site_start <= 0:
        return 0

    # get the 3' region of the mirna and the corresponding utr seq
    mirna_seq_3p = mirna_seq[8:]  # miRNA sequence from position 9 onward
    trailing = orfutr_seq[max(0, site_start - upstream_limit): site_start + 2]  # site sequence up to edges of possible 8mer site
    utr_5p = get_rc(trailing, rna=True)
    # print(mirna_seq_3p)
    # print(utr_5p)

    # initiate array for dynamic programming search
    scores = np.empty((len(utr_5p) + 1, len(mirna_seq_3p) + 1))
    scores.fill(np.nan)
    possible_scores = [0]

    # fill in array
    for i, nt1 in enumerate(utr_5p):
        for j, nt2 in enumerate(mirna_seq_3p):
            if nt1 == nt2:
                new_score = 0.5 + 0.5 * ((j > 3) & (j < 8))
                if not np.isnan(scores[i, j]):
                    new_score += scores[i, j]
                    scores[i + 1, j + 1] = new_score
                    possible_scores.append(new_score)
                else:
                    offset_penalty = max(0, (abs(i - j) - 2) * 0.5)
                    scores[i + 1, j + 1] = new_score - offset_penalty
            else:
                scores[i + 1, j + 1] = float('NaN')
    # print(scores)
    # print(np.nanmax(possible_scores))
    # print(np.nanargmax(possible_scores))
    return np.nanmax(possible_scores)



TP = dict()
TP["min_pairing"] = 2
TP["supp_l"] = 13
TP["supp_r"] = 16
TP["offset_opt"] = 0
TP["offset_tol"] = 2
TP["w_pairing"] = 0.5
TP["w_supp"] = 0.5
TP["w_offset"] = 0.5



def calculate_threep_score_alt(df_row, feature_map, mirna_seq, upstream_limit,
                               two_bmodes=False, print_check=False, 
                               target_scan=False, **kwargs):
    """
    Calculate the three-prime pairing score
    Parameters
    NOTE: This is a modified version of Kathy's function, parameterized in order
    to directly compare with her output tables.
    ----------
    df_row: a row from a dataframe taken generated from one of the 12mer files.
    feature_map: a nexted dictionary containing the ORF and UTR sequences,
        as well as their lengths (even though this is redundant by having the
        sequences)
    mirna_seq: string: the sequence of the miRNA.
    Output
    ------
    float: 3' pairing score
    """
    for key in kwargs.keys():
        TP[key] = kwargs[key]
    # First get the transcript name.
    transcript_use = df_row["transcript"]
    # Assign the sequence of the ORF-UTR.
    orfutr_seq = feature_map["orf_utr3"][transcript_use]
    # Get the position of the 12mer, by taking the ORF length, the location of
    # 12mer relative to the start of the 3'UTR (this is how it is given in 
    # Kathy's file), and converting this in to where wthin the ORF-UTR string
    # this location corresponds to.
    len_orf = int(feature_map["orf_length"][transcript_use])
    str_loc = int(df_row["utr3_loc"])
    site_start = len_orf - 3 + str_loc
    # utr = orf_utr_sequence
    # upstream_limit = 15
    if site_start <= 0:
        if print_check:
            print("here")
            print("site_start: %s" %site_start)
            print("str_loc: %s" %str_loc)
        return 0

    # get the 3' region of the mirna and the corresponding utr seq
    mirna_seq_3p = mirna_seq[8:]  # miRNA sequence from position 9 onward
    trailing = orfutr_seq[max(0, site_start - upstream_limit): site_start + 2]  # site sequence up to edges of possible 8mer site
    utr_5p = get_rc(trailing, rna=True)

    if mirna_seq_3p[2] == "G" and two_bmodes:
        bonus_score = 1 + (mirna_seq_3p[3] == "G")
    else:
        two_bmodes = False
    # initiate array for dynamic programming search
    # scores = np.empty((len(utr_5p) + 1, len(mirna_seq_3p) + 1))
    # scores.fill(np.nan)

    lens = np.empty((len(utr_5p) + 1, len(mirna_seq_3p) + 1))
    lens.fill(0)
    possible_scores = [0]

    ### These are the parameters that define the current Threep score.
    # fill in array
    terminate_bool = False
    for i, nt1 in enumerate(utr_5p):
        for j, nt2 in enumerate(mirna_seq_3p):
            if nt1 == nt2:
                lens[i + 1, j + 1] = lens[i, j] + 1
            if lens[i + 1, j + 1] >= TP["min_pairing"]:
                len_pairing = lens[i + 1, j + 1]
                # Convert the indeces into the coordinates in the miRNA.
                pos_r = j + 1 + 8
                pos_l = pos_r - len_pairing + 1
                # Calculate the Threep score for these limits.
                score_rest = TP["w_pairing"]*len_pairing
                score_supp = TP["w_supp"]*max(0, min(pos_r, TP["supp_r"]) - max(pos_l, TP["supp_l"]) + 1)
                if pos_l == 11 and pos_r >= 14 and two_bmodes:
                    print("found one")
                    print("offset: %s" %(i - j))
                    score_supp += bonus_score
                    offset_penalty = max(0, TP["w_offset"]*(abs(i - j - 4) - TP["offset_tol"]))
                else:
                    offset_penalty = max(0, TP["w_offset"]*(abs(i - j - TP["offset_opt"]) - TP["offset_tol"]))
                threep_score = score_supp + score_rest - offset_penalty
                # Add the score to the list.
                possible_scores.append(threep_score)
                # scores[i + 1, j + 1] = float('NaN')

    # Miscellaneous diagnostic things to be printed to make sure the score works
    # correctly.

    if print_check:
        # scores[0, 1:] = list(range(1, len(mirna_seq_3p) + 1))
        # scores[1:, 0] = list(range(1, len(utr_5p) + 1))
        lens[0, 1:] = list(range(1, len(mirna_seq_3p) + 1))
        lens[1:, 0] = list(range(1, len(utr_5p) + 1))
        print(utr_5p)
        print(mirna_seq_3p)
        # print(scores)
        print(lens)
        print(possible_scores)
        # print(possible_scores_2)
        print("__________")
        print(np.nanmax(possible_scores))
        print("site_start: %s" %(site_start))

        # print(np.nanmax(possible_scores_2))
    return np.nanmax(possible_scores)

