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




def main():
    time_start = time.time()
    # Parse command line arguments:
    arguments = ["cell_line", "mirna", "kds", "-upstream_limit",
                 "-min_pairing", "-supp_l", "-supp_r", "-offset_opt",
                 "-offset_tol", "-w_pairing", "-w_supp", "-w_offset", 
                 "-use_kd_model_binary", "-two_bmodes_binary",
                 "-oligoG_binary"]
    args = parse_arguments(arguments)
    (cell_line, mirna, kds, upstream_limit,
     min_pairing, supp_l, supp_r, offset_opt,
     offset_tol, w_pairing, w_supp, w_offset,
     use_kds, two_bmodes, oligoG             ) = args
    # Get the file that has the UTR sequences to be fed into the threeprime
    # threeprime score function.
    feature_map = make_feature_dictionary(cell_line)
    if not upstream_limit:
        if kds == "targetscan":
            upstream_limit = 13
        else:
            upstream_limit = 15
    # This line gets rid of the "\" character that is included when using the
    # make file to process the jobs. (THIS IS A BAND-AID BECAUSE I CAN'T FIGURE
    # OUT HOW TO GET TH EMAKE OUTPUT TO TREAT A "\" CHARACTER LIKE AN ESCAPE
    # CHARACTER AS IF I SUBMITTED THE JOB DIRECTLY FROM THE TERMINAL).
    mirna = mirna.replace("\\", "")
    # Gets the input datafile to be modified in terms of its threeprime score.
    input_df = make_input_data_table(mirna, kds, threep_canon_only=False,
                                    test=False)


    ts_off_dict = {"6mer" : -1, "7mer-1a": -1, "8mer-1a": 0, "7mer-m8": 0}
  
    total_diff = 0
    correct = 0
    wrong = 0
    output_df = input_df.copy()
    # Define the differences that are constant when fitting the biochemical
    # model and target scan.
    if kds == "targetscan":
        trans_col_str = "Gene ID"
        trans_feat_str = "utr3" # can take out
        threep_col_str = "Threep score"
    else:
        trans_col_str = "transcript"
        trans_feat_str = "orf_utr3"
        threep_col_str = "Threep"
        mirna_seq = Mirna(get_sean_mirna_name(mirna)).seq


    mirna_name_kl = get_kathy_mirna_name(mirna)
    path = "/lab/solexa_bartel/mcgeary/transfections/%s" %cell_line
    if kds == "targetscan":
        path += "/target_scan/new_threep_scores/"
    else:
        path += "/new_feature_files/%s_kds/" %kds
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
    print(path_special)
    # TWO BINDING MODES
    if two_bmodes:
        if path_special == "":
            path_special = "two_bmodes"
        else:
            path_special += ",two_bmodes"

    # NONLINEAR BENFIT TO G
    if oligoG:
        if path_special == "":
            path_special = "oligoG"
        else:
            path_special += ",oligoG"
    if path_special != "":
        path += "%s/" %path_special
    print(path)
    ensure_directory(path)
    if kds == "targetscan":
        path += "features.txt"
    else:
        path += "%s.txt" %mirna_name_kl
    print(path)
    # Get the sequence of the miRNA.
 
    # Pre-allocate the output array.
    output_df = input_df.copy()
    output_df["Threep_pairing"] = np.zeros([input_df.shape[0]])
    # This matrix is not strictly necessary; it's there in order to compare that
    # the threepriime pairing score hasn't been changed in any unintended ways
    # whille modifying the code to allow the various features of the score to be
    # varied.
    # output_check = np.zeros([input_df.shape[0], 2])
    # for i_row in range(input_df.shape[0]):
    # # for i_row in range(10):
    #     df_row = input_df.iloc[i_row, ]
    #     # print(df_row)
    #     threep_kl = df_row["Threep"]
    #     threep_sm = calculate_threep_score_alt(
    #         df_row, feature_map, mirseq, int(upstream_limit), two_bmodes,
    #         **kwargs
    #     )
    #     output_check[i_row, ] = [threep_kl, threep_sm]
    #     # print(output_check[i_row, ])
    #     # if threep_kl != threep_sm:
    #     #     print(output_check[i_row, ])
    #     output_df.loc[i_row, "Threep"] = threep_sm




    # Iterate over the rows:
    # for row_i in range(100):
    for row_i in range(input_df.shape[0]):
        # Get the row.
        row_use = input_df.iloc[row_i, ]
        # If using target scan, get the specific miRNA and site type in order to
        # get the ranges of the positions that are used.
        if kds == "targetscan":
            site = row_use["Site type"]
            mirna = get_sean_mirna_name(row_use["miRNA family"])
            mirna_seq = Mirna(mirna).seq
            start_str = "Site start"
            offset_start = ts_off_dict[site] - 2
            offset_limit = ts_off_dict[site]
        else:
            start_str = "loc"
            offset_start = -3
            offset_limit = 0
        # Define the offset positions and the upstream limit.
        start_pos = int(row_use[start_str]) + offset_start
        upstream_use = upstream_limit + offset_limit
        # Get Kathy's score (only really necessary when proofreading.)
        threep_kl = row_use[threep_col_str]
        # Get the transcript sequence.
        transcript = feature_map[trans_feat_str][row_use[trans_col_str]]
        threep_sm_tup = calculate_threep_score_new(
            transcript, mirna_seq, start_pos, upstream_use, two_bmodes, oligoG,
            **kwargs
        )
        threep_sm = threep_sm_tup[0]
        if threep_kl == threep_sm:
            correct += 1
        if threep_kl != threep_sm:
            wrong += 1
            # print("%s\t%s" %(threep_kl, threep_sm))
        diff = (threep_kl - threep_sm)**2
        output_df.loc[row_i, threep_col_str] = threep_sm
        output_df.loc[row_i, "Threep_pairing"] = threep_sm_tup[1]
        total_diff += diff
    print(input_df.head())
    print(output_df.head())
    print(total_diff)
    print("%s\t%s" %(correct, wrong))
    print(path)
    output_df.to_csv(path, sep="\t", index=False, index_label=False)




    time_done = time.time()
    print(("Finished in %3.2f seconds." %(time_done - time_start)))
    return()

if __name__ == "__main__":
    main()
