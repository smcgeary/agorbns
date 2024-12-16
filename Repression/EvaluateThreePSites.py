################################################################################
#rDevelopThreePScore.py
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
imp.load_source("repression_functions",
                ("/lab/bartel1_ata/mcgeary/computation/"
                 "AgoRBNS/Repression/repression_functions.py"
                )
               )

from general import *
from RBNS_methods import *
from repression_functions import *
from sitetypes import get_seq_site_map


np.set_printoptions(edgeitems=10, linewidth=100000)

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


# K_supp_min, K_supp_max = [13, 16] # The range defining the 3p supplemental region
# K_offset_opt = 0 # The optimal offset.
# K_offset_tol = 2 # The width of tolerance for deviation from the optimal offset.
# K_pair_weight = 0.5 # The weight given to any contiguous pairing.
# K_supp_weight = 0.5 # The weight given to any pairing in the 3p supplemental region.
# K_offset_weight = 0.5 # The negative weight given for offsets that deviate more than `offset_tol` away from `offset_opt`.
# K_min_pairing = 2 # The minimum length of contiguous pairing required for something to be counted as having threeprime pairing.





def calculate_threep_score_alt(df_row, feature_map, mirna_seq, upstream_limit,
                               two_bmodes=False, print_check=False, **kwargs):
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
        return [0, 0]

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
    possible_lens = [0]
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
                possible_lens.append(int(len_pairing))

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
    ind_use = np.nanargmax(possible_scores)
    len_used = possible_lens[ind_use]
    return([possible_scores[ind_use], possible_lens[ind_use]])




def main():
    time_start = time.time()
    # Parse command line arguments:
    arguments = ["cell_line", "mirna", "kds", "upstream_limit",
                 "-min_pairing", "-supp_l", "-supp_r", "-offset_opt",
                 "-offset_tol", "-w_pairing", "-w_supp", "-w_offset", 
                 "-use_kd_model_binary", "-two_bmodes_binary"]
    args = parse_arguments(arguments)
    (cell_line, mirna, kds, upstream_limit,
     min_pairing, supp_l, supp_r, offset_opt,
     offset_tol, w_pairing, w_supp, w_offset,
     use_kds, two_bmodes                     ) = args
    # Get the file that has the UTR sequences to be fed into the threeprime
    # threeprime score function.
    feature_map = make_feature_dictionary(cell_line)
    # This line gets rid of the "\" character that is included when using the
    # make file to process the jobs. (THIS IS A BAND-AID BECAUSE I CAN'T FIGURE
    # OUT HOW TO GET TH EMAKE OUTPUT TO TREAT A "\" CHARACTER LIKE AN ESCAPE
    # CHARACTER AS IF I SUBMITTED THE JOB DIRECTLY FROM THE TERMINAL).
    mirna = mirna.replace("\\", "")
    # Gets the input datafile to be modified in terms of its threeprime score.
    input_df = make_input_data_table(mirna, kds, threep_canon_only=False,
                                    test=False)

    mirna_name_kl = get_kathy_mirna_name(mirna)
    path = "/lab/solexa_bartel/mcgeary/transfections/%s" %cell_line
    path += "/ThrP_counts/%s_kds/" %kds
    kwargs= dict()
    if min_pairing:
        kwargs["min_pairing"] = int(min_pairing)
    if supp_l:
        kwargs["supp_l"] = int(supp_l)
    if supp_r:
        kwargs["supp_r"] = int(supp_r)
    if offset_opt:
        kwargs["offset_opt"] = int(offset_opt)
    if offset_tol:
        kwargs["offset_tol"] = int(offset_tol)
    if w_pairing:
        kwargs["w_pairing"] = float(w_pairing)
    if w_supp:
        kwargs["w_supp"] = float(w_supp)
    if w_offset:
        kwargs["w_offset"] = float(w_offset)

    path_special = ",".join(["%s_%s" %(key, kwargs[key]) for key in kwargs.keys()])
    if two_bmodes:
        if path_special == "":
            path_special = "two_bmodes"
        else:
            path_special += ",two_bmodes"
    if path_special != "":
        path += "%s/" %path_special
    ensure_directory(path)
    path += "%s_len.txt" %mirna_name_kl
    # Get the sequence of the miRNA.
    mirseq = Mirna(mirna).seq
    # Pre-allocate the output array.
    output_df = input_df.copy()
    # Get the 8mer mismatch site names:
    _mirna = Mirna(mirna)
    seq_8mer = _mirna["8mer"]
    names_8mer_mm = ["8mer-mm%s%s" %(j, 8 - i)
                    for i in range(1, 7)
                    for j in DNTS if j != seq_8mer[i]]
    seqs_8mer_mm = [_mirna[name] for name in names_8mer_mm]
    print(names_8mer_mm)
    print(seqs_8mer_mm)
    name_seq_dict = {seq : name for (seq, name) in zip(seqs_8mer_mm, names_8mer_mm)}
    print(name_seq_dict)

    # This matrix is not strictly necessary; it's there in order to compare that
    # the threepriime pairing score hasn't been changed in any unintended ways
    # whille modifying the code to allow the various features of the score to be
    # varied.
    output_check = np.zeros([input_df.shape[0], 2])
    # The output of this file, is a dictionary reporting on the lengths of each
    # of the three-p pairing configurations used for each three-p score.
    count_ThrPLen_map = {"seed" : collections.defaultdict(int),
                         "8mer-mm" : collections.defaultdict(int),
                         "no site" : collections.defaultdict(int)}
    # Iterate over the rows.
    for i_row in range(input_df.shape[0]):
    # for i_row in range(1000):
        df_row = input_df.iloc[i_row, ]
        threep_kl = df_row["Threep"]
        kmer = df_row["12mer"]
        kmer_8 = kmer[2:10]
        stype = df_row["stype"]
        # Conditional logic to assign the threep site length to either "seed",
        # "8mer-mm", or "no site".
        if kmer_8 in seqs_8mer_mm:
            stype = "8mer-mm"
        elif stype != "no site":
            stype = "seed"
        threep_sm = calculate_threep_score_alt(
            df_row, feature_map, mirseq, int(upstream_limit), two_bmodes,
            **kwargs
        )
        output_check[i_row, ] = [threep_kl, threep_sm[0]]
        count_ThrPLen_map[stype][str(threep_sm[1])] += 1
        if threep_kl != threep_sm[0]:
            print(output_check[i_row, ])
        output_df.loc[i_row, "Threep"] = threep_sm[0]
    # Convert the dictionary into a pandas dataframe.
    out_df = pd.DataFrame.from_dict(count_ThrPLen_map, orient="columns")
    out_df.fillna(0, inplace=True)
    out_df = out_df.astype(int)
    out_df.index = out_df.index.astype(int)

    # Sort the datagrame by the indeces (the lengths of 3p pairing)
    out_df.sort_index(axis=0, inplace=True, ascending=True)
    print(path)
    # Write the file.
    out_df.to_csv(path, sep="\t", index_label=False)
    # Print time taken exit script.
    time_done = time.time()
    print(("Finished in %3.2f seconds." %(time_done - time_start)))
    return()

if __name__ == "__main__":
    main()
