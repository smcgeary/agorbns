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


TP = dict()
TP["min_pairing"] = 2
TP["supp_l"] = 13
TP["supp_r"] = 16
TP["offset_opt"] = 0
TP["offset_tol"] = 2
TP["w_pairing"] = 0.5
TP["w_supp"] = 0.5
TP["w_offset"] = 0.5



def calculate_kd_threep_score(df_row, feature_map, mirna_seq,
                               pars_use, print_check=False, **kwargs):
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
    upstream_limit = 20
    if site_start <= 0:
        return([0, 0, 0, ""])

    # get the 3' region of the mirna and the corresponding utr seq
    mirna_seq_3p = mirna_seq[8:]  # miRNA sequence from position 9 onward
    trailing = orfutr_seq[max(0, site_start - upstream_limit): site_start + 2]  # site sequence up to edges of possible 8mer site
    utr_5p = get_rc(trailing, rna=True)

    # if mirna_seq_3p[2] == "G" and two_bmodes:
    #     bonus_score = 1 + (mirna_seq_3p[3] == "G")
    # else:
    #     two_bmodes = False
    # initiate array for dynamic programming search
    scores = np.empty((len(utr_5p) + 1, len(mirna_seq_3p) + 1))
    scores.fill(np.nan)
    # Initialize array for dynamic programming search with Kd scores.
    kd_scores = np.empty((len(utr_5p) + 1, len(mirna_seq_3p) + 1))
    kd_scores.fill(np.nan)
    # Initialize array to track the lengths of the ThreePrime scores.
    lens = np.empty((len(utr_5p) + 1, len(mirna_seq_3p) + 1))
    lens.fill(0)
    possible_scores = [0]
    possible_kd_scores = [0]
    possible_lens = [0]
    possible_pairings = ['']
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
                offset = i - j
                # if pos_l == 11 and pos_r >= 14 and two_bmodes:
                #     print("found one")
                #     print("offset: %s" %offset)
                #     score_supp += bonus_score
                #     offset_penalty = max(0, TP["w_offset"]*(abs(offset - 4) - TP["offset_tol"]))
                # else:
                offset_penalty = max(0, TP["w_offset"]*(abs(offset - TP["offset_opt"]) - TP["offset_tol"]))
                # Calculate the original ThreeP score.
                threep_score = score_supp + score_rest - offset_penalty
                # Define the pairing and offset labels to calcultae the Kd-based
                # ThreeP score.
                pairing_label = "%i|%i" %(pos_l, pos_r)
                offset_label = str(int(offset))
                # Look up the pairing and offset terms.
                if pairing_label in pars_use[0].index.values:
                    pairing_param = pars_use[0].loc[pairing_label]
                else:
                    print(pairing_label)
                    print(len_pairing - 8 + 1)
                    labels_try = ["%i|%i" %(pos_l + i, pos_l + i + 8 - 1) for i in range(int(len_pairing) - 8 + 1)]
                    print(labels_try)
                    pairings_try = [pars_use[0].loc[pairing_i] for pairing_i in labels_try]
                    print(pairings_try)
                    pairing_param = max(pairings_try)

                if offset_label in pars_use[1].index.values:
                    offset_param = pars_use[1].loc[offset_label]
                else:
                    offset_param = 0
                # Calculate the Kd-based ThreeP score.
                kd_score = pairing_param*offset_param
                # Mark the specific pairing that is beign used.
                pairing_str = "%s_%s" %(pairing_label, offset_label)
                # Add the score to the list.
                possible_scores.append(threep_score)
                possible_kd_scores.append(kd_score)
                possible_lens.append(int(len_pairing))
                possible_pairings.append(pairing_str)
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
        print(lens)
        print(possible_scores)
        print(possible_kd_scores)
        print(possible_pairings)
    ind_threep = np.nanargmax(possible_scores)
    len_threep = possible_lens[ind_threep]
    ind_kd = np.nanargmax(possible_kd_scores)
    pairing_kd = possible_pairings[ind_kd]
    pairing_threep = possible_pairings[ind_threep]
    return([possible_scores[ind_threep], len_threep,
            possible_kd_scores[ind_kd], pairing_kd, pairing_threep])




def main():
    time_start = time.time()
    # Parse command line arguments:
    arguments = ["cell_line", "mirna", "kds", "mirna_pars", "experiment_pars",
                 "-decomp_binary", "-fit_max_len", "-keep_data_binary"]
    args = parse_arguments(arguments)
    (cell_line, mirna, kds, mirna_pars, experiment_pars, decomp, 
     fit_max_len, keep_data) = args
    # Get the file that has the UTR sequences to be fed into the threeprime
    # threeprime score function.
    feature_map = make_feature_dictionary(cell_line)
    ## LOAD THE KD PARAMETERS FROM THE MODELING.
    if decomp:
        decomp_str = "_decomp_max%s" %fit_max_len
        if keep_data:
            keepdata_str = "_keepdata"
        else:
            keepdata_str = ""
    else:
        decomp_str = ""
        keepdata_str = ""
    kd_pars_path = (
        "/lab/solexa_bartel/mcgeary/AgoRBNS"
        "/%s/%s/threep_models/pairing_and_offset%s%s.txt" %(
            mirna_pars, experiment_pars, decomp_str, keepdata_str
        )
    )
    pars = pd.DataFrame.from_csv(kd_pars_path, sep="\t")["MLE"]
    # Spilt up the parameters into pairing and offset coefficients.
    pars_pairs = pars.filter(regex="\\|", axis=0)
    pars_offset = pars.filter(regex="^((?!\\|).)*$", axis=0)
    pars_all = [pars_pairs, pars_offset]
    # This line gets rid of the "\" character that is included when using the
    # make file to process the jobs. (THIS IS A BAND-AID BECAUSE I CAN'T FIGURE
    # OUT HOW TO GET THE MAKEFILE OUTPUT TO TREAT A "\" CHARACTER LIKE AN ESCAPE
    # CHARACTER AS IF I SUBMITTED THE JOB DIRECTLY FROM THE TERMINAL).
    mirna = mirna.replace("\\", "")
    if mirna == "let-7a":
        mirna_use = "let-7a-21nt"
    else:
        mirna_use = mirna
    # Gets the input datafile to be modified in terms of its threeprime score.
    input_df = make_input_data_table(mirna, kds, threep_canon_only=False,
                                     test=False)

    mirna_name_kl = get_kathy_mirna_name(mirna)
    out_feat_path = "/lab/solexa_bartel/mcgeary/transfections/%s" %cell_line
    out_feat_path += "/kd_feature_files/%s_kds/" %kds
    if experiment_pars == "equilibrium" or experiment_pars == "equilibrium2_nb":
        out_feat_path += "random/"
    else:
        out_feat_path += "programmed/"
    ensure_directory(out_feat_path)
    out_feat_path += "%s.txt" %mirna_name_kl
    # Get the sequence of the miRNA.
    mirseq = Mirna(mirna_use).seq
    # Get the 8mer mismatch site names:
    _mirna = Mirna(mirna)
    # Get the list of names of the 8mer-mm sites.
    seq_8mer = _mirna["8mer"]
    names_8mer_mm = ["8mer-mm%s%s" %(j, 8 - i)
                    for i in range(1, 7)
                    for j in DNTS if j != seq_8mer[i]]
    seqs_8mer_mm = [_mirna[name] for name in names_8mer_mm]
    name_seq_dict = {seq : name for (seq, name)
                     in zip(seqs_8mer_mm, names_8mer_mm)}
    # Pre-allocate the output matrix.
    output_df = input_df.copy()
    output_df["Threep_pairing"] = 'none'
    output_df["Kd_score"] = 0
    output_df["Kd_pairing"] = 'none'
    output_df["Threep_len"] = 0
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
    # Output dictionary for counting the lengths of each of the chosen sites
    # based on the modeling coefficients.
    count_ThrPKdLen_map = {"seed" : collections.defaultdict(int),
                         "8mer-mm" : collections.defaultdict(int),
                         "no site" : collections.defaultdict(int)}
    # Iterate over the rows.
    # for i_row in range(input_df.shape[0]):
    for i_row in range(20):
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
        threep_sm = calculate_kd_threep_score(
            df_row, feature_map, mirseq, pars_all
        )
        output_check[i_row, ] = [threep_kl, threep_sm[2]]
        count_ThrPLen_map[stype][str(threep_sm[1])] += 1
        count_ThrPKdLen_map[stype][str(threep_sm[3])] += 1
        # if threep_kl != threep_sm[0]:
        #     print(output_check[i_row, ])
        output_df.loc[i_row, "Threep"] = threep_sm[0]
        output_df.loc[i_row, "Threep_pairing"] = threep_sm[1]
        output_df.loc[i_row, "Kd_score"] = threep_sm[2]
        output_df.loc[i_row, "Kd_pairing"] = threep_sm[3]
        output_df.loc[i_row, "Threep_len"] = threep_sm[4]
    print(output_df.iloc[:20, :])
    return()
    # Convert the dictionary into a pandas dataframe.
    out_df = pd.DataFrame.from_dict(count_ThrPKdLen_map, orient="columns")
    out_df.fillna(0, inplace=True)
    out_df = out_df.astype(int)
    out_df.index = out_df.index.astype(int)

    count_path = "/lab/solexa_bartel/mcgeary/transfections/%s" %cell_line
    count_path += "/ThrP_counts/%s_kds/" %kds
    count_path += "%skd_pairings.txt" %mirna_name_kl


    # Sort the datagrame by the indeces (the lengths of 3p pairing)
    # out_df.sort_index(axis=0, inplace=True, ascending=True)
    print(path)
    # Write the file.
    output_df.to_csv(path, sep="\t", index=False, index_label=False)
    # Print time taken exit script.
    # Write the count file.
    out_df.to_csv(count_path, sep="\t", index_label=False)

    time_done = time.time()
    print(("Finished in %3.2f seconds." %(time_done - time_start)))
    return()

if __name__ == "__main__":
    main()
