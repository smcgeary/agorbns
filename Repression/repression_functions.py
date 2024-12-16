################################################################################
#repression_test.py
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
from scipy.optimize import minimize
from scipy import stats


pd.set_option("display.precision", 6)



base_exp_dir = "/lab/solexa_bartel/mcgeary/transfections/"
star_index = ("/nfs/genomes/human_gp_feb_09_no_random/STAR_2.6.0/"
              "GRCh37.75_overhang_100/")
annot_file = "%smetadata/annotations.gtf" %base_exp_dir


kds_base_dir = ("/lab/solexa_bartel/klin/miRNA_models_data_old/"
                "model_inputs/biochem/")
measured_kds_dir = "measured_kds/"
predicted_kds_dir = ("predicted_kds/feat1_allgenes_lr003_nodropout_batch50"
                     "_rebalancekds2k_noaugment_repweight095_mask_w3_netpred/"
                     "hela/")



TP = dict()
TP["min_pairing"] = 2
TP["supp_l"] = 13
TP["supp_r"] = 16
TP["offset_opt"] = 0
TP["offset_tol"] = 2
TP["w_pairing"] = 0.5
TP["w_supp"] = 0.5
TP["w_offset"] = 0.5



def mean_center_df(df):
    return(df.sub(df.mean(axis=1), axis=0))


def sigmoid(x):
    # x = x.astype(float)
    return(1 / (1 + np.exp(-x.astype(float))))


def get_kathy_mirna_name(mirna):
    mirna = mirna.replace("*", "_pass").replace("-", "").replace("R", "r")
    mirna = mirna.replace("let7a", "let7")
    return(mirna)

def get_sean_mirna_name(mirna):
    mirna = mirna.replace("_pass", "*").replace("mir", "miR-")
    mirna = mirna.replace("lsy6", "lsy-6").replace("let7", "let-7a")
    return(mirna)

def make_feature_dictionary(cell_line):
    transcript_file_dir = ("/lab/solexa_bartel/mcgeary/transfections"
                           "/metadata/%s_transcripts_2.txt" %(cell_line))
    # Make a dictionary that has the information on each transcript in order to
    # ensure I am performing the analysis the same way that Kathy was doing it.
    with open(transcript_file_dir) as transcript_file:
        # Make the outer layer of keys.
        outer_keys = transcript_file.readline().strip().split("\t")[1:]
        features_mRNA_map = {i : collections.defaultdict(lambda: "Key does not exist") for i in outer_keys}
        line = transcript_file.readline()
        while line:
            line_list = line.strip().split("\t")
            transcript_id = line_list[0]
            features = line_list[1:]
            for i, item in enumerate(features):
                features_mRNA_map[outer_keys[i]][transcript_id] = item
            line = transcript_file.readline()
    return(features_mRNA_map)

def make_input_data_table(mirna, kds, threep_canon_only=True, test=False):
    # if 'kds' == targetscan, bypass the miRNA designation and output the single
    # file that has the targetscan features.
    if kds == "targetscan":
        path = ("/lab/bartel4_ata/kathyl/RNA_Seq/data/TS7_output_files/"
                "refseq_utr90_fixed_ribosome_shadow/features.txt")
        print(path)
        input_df = pd.read_csv(path, sep="\t")
        return(input_df)
    else:
        mirna_file = get_kathy_mirna_name(mirna) + ".txt"
        if kds == "measured":
            path = kds_base_dir + measured_kds_dir + mirna_file
        elif kds == "predicted":
            path = kds_base_dir + predicted_kds_dir + mirna_file
        else:
            raise ValueError("Incorrect 'kds' assignment.")
        print(path)
        input_df = pd.read_csv(path, sep="\t")
        # Give all noncanonical sites an structural score equal to the mean over
        # all of the canonical sites.
        logSA_mean = np.nanmean(input_df["logSA_diff"])
        input_df["logSA_diff"].fillna(logSA_mean, inplace=True)
        # Apply the Kd cutoff to the data:
        input_df = input_df[input_df["log_kd"] < 0].reset_index(drop=True)
            # Assign a threep score of zero to all noncanonical sites.
        if threep_canon_only:
            input_df.loc[(input_df["stype"] == "no site"), "Threep"] = 0
        if test:
            input_df = input_df.iloc[:50, ]
        return(input_df)

def make_full_input_data_table(mirnas, kds, threep_canon_only=True, test=False):
    # Iterate over the list of miRNAs to get each table.
    input_df = pd.concat(
        [make_input_data_table(mirna, kds, threep_canon_only=threep_canon_only,
                               test=test) for mirna in mirnas],
        axis=0
    )
    # Reformat the strings in the miRNA column.
    input_df["mir"] = input_df["mir"].apply(lambda x: x.replace("*", "_pass"))
    input_df["mir_gp"] = input_df["mir"].apply(lambda x: x.replace("_pass", ""))
    return(input_df)

def make_transfection_data_table(cell_line, mirnas):
    # Path to transfection data.
    input_data_path = ("/lab/solexa_bartel/mcgeary/transfections/%s/"
                       "count_tables/logtpm_batchnormalized.txt" %cell_line)
    #Input data.
    data_df = pd.read_csv(input_data_path, sep="\t", index_col=0)
    # Subset the data (will make this an option later on, now just trying
    # Figure 5).
    if mirnas == "five":
        data_df = data_df[["lsy-6", "miR-1", "miR-124", "miR-155", "miR-7"]]
    elif mirnas == "sixteen":
        data_df = data_df[["lsy-6", "miR-1", "miR-124", "miR-137",
                           "miR-139", "miR-143", "miR-144", "miR-153",
                           "miR-155", "miR-182", "miR-199a", "miR-204",
                           "miR-205", "miR-216b", "miR-223", "miR-7", ]]
    # Mean center these predictions:
    return(mean_center_df(data_df))

def get_predictions(pars, in_df, mirnas_use, biochemplus=False):
    # Assign the parameters:
    # Assign the free miRNA concentrations.
    ag_s = pars[0] + np.append([0], pars[1:len(mirnas_use)])
    # Assign the repression and orf-site penalty parameters.
    b, c = pars[len(mirnas_use):len(mirnas_use) + 2]
    # Assign the feature weights with either the input parameters for the
    # biochemplus model, or  0 for the biochem model.
    if biochemplus:
        cSA, cThrP, cPCT = pars[len(mirnas_use) + 2:]
    else:
        cSA, cThrP, cPCT = [0, 0, 0]
    out_df = in_df.copy()
    ag_mirna_map = pd.DataFrame.from_dict(
        {mirna : ag_i for mirna, ag_i in zip(mirnas_use, ag_s)},
        orient="index"
    )
    # Calculate the Kd offset from the features (will be 0 for model="biochem")
    out_df["offset"] = cThrP*in_df["Threep"] + cPCT*in_df["PCT"] 
    out_df["offset_all"] = cSA*in_df["logSA_diff"] + c*in_df["in_ORF"]
    # Subtract offset from log_kd
    # out_df["log_kd_site"] = in_df["log_kd"] - (cSA*in_df["logSA_diff"] + cThrP*in_df["Threep"] + cPCT*in_df["PCT"])
    # out_df["log_kd_bg"] = -1*cSA*in_df["logSA_diff"]
    out_df["a_g"] = ag_mirna_map.loc[out_df["mir"]].reset_index(drop=True)
    # Get the specifically bound sites:
    out_df["N_g"]    = 1/(1 + np.exp(in_df["log_kd"] - out_df["a_g"] - out_df["offset"] - out_df["offset_all"]))
    # Get the background binding:
    out_df["N_g,bg"] = 1/(1 + np.exp(-1*out_df["offset_all"] - out_df["a_g"]))
    # Sum the N_g and N_g,bg values across each transcript-and-miRNA combination.
    # print(out_df[["mir", "Threep", "PCT", "in_ORF", "logSA_diff", "a_g", "log_kd", "offset", "offset_all", "N_g", "N_g,bg"]])
    N_full = out_df.groupby(["transcript", "mir_gp"]).agg({'N_g': np.sum,
                                                        'N_g,bg': np.sum})
    N_full.reset_index(level=[0,1], inplace=True)
    # Convert these occupancies into the final predicted repression values:
    N_full["pred"] = np.log((1 + np.exp(b)*N_full["N_g,bg"])/(1 + np.exp(b)*N_full["N_g"]))
    # Convert data table back into a matrix with dimensions tranxcript X miRNA,
    # and replace the NaN values (which are transcripts for which there were no
    # 12 mers against a particular miRNA) with 0, both the background binding
    # and specific binding would be the same in this case.
    pred_df = N_full.pivot(index="transcript",
                           columns="mir_gp",
                           values="pred").fillna(0)
    # Change column names to be consistent with data_df.
    pred_df.columns = [get_sean_mirna_name(mirna) for mirna 
                       in pred_df.columns.tolist()]
    return(pred_df)

def get_predictions_final(pars, in_df, mirnas_use, biochemplus=False):
    # Assign the parameters:
    # Assign the free miRNA concentrations.
    ag_s = pars[0] + np.append([0], pars[1:len(mirnas_use)])
    # Assign the repression and orf-site penalty parameters.
    b, c = pars[len(mirnas_use):len(mirnas_use) + 2]
    # Assign the feature weights with either the input parameters for the
    # biochemplus model, or  0 for the biochem model.
    if biochemplus:
        cSA, cThrP, cPCT = pars[len(mirnas_use) + 2:]
    else:
        cSA, cThrP, cPCT = [0., 0., 0.]
    mir_dict = {mirna : ag_i for mirna, ag_i in zip(mirnas_use, ag_s)}

    ag_mirna_map = pd.DataFrame.from_dict(mir_dict, orient="index")

    # Calculate the Kd offset from the features (will be 0 for model="biochem")
    in_df.loc[:, "offset"] = cThrP*in_df["Threep"] + cPCT*in_df["PCT"]
    # offset = cThrP*in_df_np[:, col_name_map["Threep"]] + cPCT*in_df_np[:, col_name_map["PCT"]]
    in_df.loc[:, "offset_all"] = cSA*in_df["logSA_diff"] + c*in_df["in_ORF"]
    # offset_all = cSA*in_df_np[:, col_name_map["logSA_diff"]] + c*in_df_np[:, col_name_map["in_ORF"]]

    # Subtract offset from log_kd
    in_df["a_g"] = ag_mirna_map.loc[in_df["mir"]].reset_index(drop=True)
    # a_g = np.array([mir_dict[i] for i in in_df_np[:, col_name_map["mir"]]])
    # Get the specifically bound sites:
    in_df.loc[:, "N_g"]    = sigmoid(in_df["log_ka"] + in_df["a_g"] + in_df["offset"] + in_df["offset_all"])
    # N_g = sigmoid(in_df_np[:, col_name_map["log_ka"]] + a_g + offset + offset_all)
    in_df.loc[:, "N_g,bg"] = sigmoid(in_df["offset_all"] + in_df["a_g"])
    # N_gbg = sigmoid(offset_all + a_g)
    # Sum the N_g and N_g,bg values across each transcript-and-miRNA combination.
    # print(out_df[["mir", "Threep", "PCT", "in_ORF", "logSA_diff", "a_g", "log_kd", "offset", "offset_all", "N_g", "N_g,bg"]])
    N_full = in_df.groupby(["mir_gp", "transcript"]).agg({'N_g': np.sum,
                                                        'N_g,bg': np.sum})
    # N_full_alt = np.array([[np.sum(N_g[l_i:r_i]), np.sum(N_gbg[l_i:r_i])] for (l_i, r_i) in limits])


    N_full.reset_index(level=[0,1], inplace=True)
    # Convert these occupancies into the final predicted repression values:
    # pred_df = np.log1p(np.exp(b)*N_full["N_g,bg"]) - np.log1p(np.exp(b)*N_full["N_g"])


    # pred_df = np.log1p(np.exp(b)*N_full_alt[:, 1]) - np.log1p(np.exp(b)*N_full_alt[:, 0])

    # # print(pred_df)
    # # print(pred_df_alt)
    # return(pred_df)

    N_full["pred"] = np.subtract(np.log1p(np.exp(b)*N_full["N_g,bg"]), np.log1p(np.exp(b)*N_full["N_g"]))
    # Convert data table back into a matrix with dimensions tranxcript X miRNA,
    # and replace the NaN values (which are transcripts for which there were no
    # 12 mers against a particular miRNA) with 0, both the background binding
    # and specific binding would be the same in this case.
    pred_df = N_full.pivot(index="transcript",
                           columns="mir_gp",
                           values="pred").fillna(0)
    # # Change column names to be consistent with data_df.
    pred_df.columns = [get_sean_mirna_name(mirna) for mirna 
                       in pred_df.columns.tolist()]
    return(pred_df)




def calculate_threep_score_new(target_seq, mirna_seq, site_start, upstream_limit,
                               two_bmodes=False, oligoG=False, print_check=False, 
                               **kwargs):
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
    # for key in kwargs.keys():
    #     TP[key] = kwargs[key]
    # # First get the transcript name.
    # transcript_use = df_row["transcript"]
    # # Assign the sequence of the ORF-UTR.
    # orfutr_seq = feature_map["orf_utr3"][transcript_use]
    # Get the position of the 12mer, by taking the ORF length, the location of
    # 12mer relative to the start of the 3'UTR (this is how it is given in 
    # Kathy's file), and converting this in to where wthin the ORF-UTR string
    # this location corresponds to.
    # len_orf = int(feature_map["orf_length"][transcript_use])
    # str_loc = int(df_row["utr3_loc"])
    # site_start = len_orf - 3 + str_loc
    # utr = orf_utr_sequence
    # upstream_limit = 15
    for key in kwargs.keys():
        TP[key] = kwargs[key]
    if site_start <= 0:
        return([0, ""])

    # get the 3' region of the mirna and the corresponding utr seq
    mirna_seq_3p = mirna_seq[8:]  # miRNA sequence from position 9 onward
    # trailing = target_seq[max(0, site_start - upstream_limit - 2 ): site_start]  # site sequence up to edges of possible 8mer site
    trailing = target_seq[max(0, site_start - upstream_limit): site_start + 2]  # site sequence up to edges of possible 8mer site
    utr_5p = get_rc(trailing, rna=True)

    # print(trailing + " " + str(len(trailing)))
    # print(" "*len(trailing) + mirna_seq[7::-1])
    # print(trailing + target_seq[site_start+2:site_start + 10])
    # print(target_seq[max(0, site_start - upstream_limit):max(0, site_start - upstream_limit) + 50])
    min_reach = 11
    if (mirna_seq_3p[1] == "G" or mirna_seq_3p[2] == "G") and two_bmodes:
        if mirna_seq_3p[1] == "G":
            min_reach = 10
        bonus = 1*((mirna_seq_3p[1] == "G") + (mirna_seq_3p[2] == "G") + (mirna_seq_3p[3] == "G"))
    else:
        two_bmodes = False
    # initiate array for dynamic programming search

    lens = np.empty((len(utr_5p) + 1, len(mirna_seq_3p) + 1))
    lens.fill(0)
    possible_scores = [0]
    possible_pairings = ['']

    ### These are the parameters that define the current Threep score.
    # fill in array
    terminate_bool = False
    print_pairing = False
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
                if pos_l == min_reach and len_pairing >= 5 and two_bmodes:
                    print("found one")
                    print("offset: %s" %(i - j))
                    score_supp += bonus_score
                    offset_penalty = max(0, TP["w_offset"]*(abs(offset - 4) - TP["offset_tol"]))
                else:
                    # offset_penalty = max(0, TP["w_offset"]*(abs(i - j - TP["offset_opt"]) - TP["offset_tol"]))
                    offset_penalty = max(0, TP["w_offset"]*(abs(offset - TP["offset_opt"]) - TP["offset_tol"]))
                # Calculate the threep score.
                threep_score = score_supp + score_rest - offset_penalty
                # Define the pairing and offset labels to calcultae the Kd-based
                # ThreeP score.
                pairing_label = "%i|%i" %(pos_l, pos_r)
                mir_seq = mirna_seq_3p[int(j - len_pairing + 1):j]
                if oligoG:
                    G_bonus = 0
                    G_it = 0
                    for nuc in mir_seq:
                        if nuc == "G":
                            G_it += 1
                            G_bonus += G_it
                        else:
                            G_it = 0
                    if G_bonus > 1:
                        print("%s\t%s\t%s\t%s" %(mir_seq, G_bonus, threep_score, threep_score + G_bonus))
                    threep_score += G_bonus
                offset_label = str(int(offset))
                # Mark the specific pairing that is beign used.
                pairing_str = "%s_%s" %(pairing_label, offset_label)
                # Add the score and the pairing to the respective list.
                possible_scores.append(threep_score)
                possible_pairings.append(pairing_str)

    # Miscellaneous diagnostic things to be printed to make sure the score works
    # correctly.
    if print_check:
        # scores[0, 1:] = list(range(1, len(mirna_seq_3p) + 1))
        # scores[1:, 0] = list(range(1, len(utr_5p) + 1))
        lens[0, 1:] = list(range(1, len(mirna_seq_3p) + 1))
        lens[1:, 0] = list(range(1, len(utr_5p) + 1))
        print(utr_5p)
        print(mirna_seq_3p)
        print(lens)
        print(possible_scores)
        print("__________")
        print(np.nanmax(possible_scores))
        # print(np.nanmax(possible_scores_2))
        print("site_start: %s" %(site_start))
    ind_threep = np.nanargmax(possible_scores)

    return([possible_scores[ind_threep], possible_pairings[ind_threep]])





def grad_model(pars, input_df, data_df, mirnas, biochemplus):
    return()


import itertools

def calc_rsquared(xs, ys):
    r = stats.linregress(xs, ys)[2]
#     return r**2
    return np.sign(r) * (r**2)

def get_r2s(filename, mirs, model_name, batch=None, transcript_info=None, batch_size=None):
    all_r2s = []
    all_df = []
    for mir in mirs:
        temp = pd.read_csv(filename.replace('MIR', mir), sep='\t')
        temp = temp[temp['mir'] == mir]
        if batch is not None:
            temp['transcript'] = [x.replace('\'','').replace('b','') for x in temp['transcript']]
            temp['batch'] = transcript_info.loc[temp['transcript'].values]['batch'].values
            temp = temp[temp['batch'] == batch]
            assert len(temp) == batch_size
        all_df.append(temp)
        all_r2s.append([mir, model_name, calc_rsquared(temp['pred_normed'], temp['label_normed'])])
    all_df = pd.concat(all_df)
    if len(mirs) > 1:
        all_r2s.append([f'{len(mirs)} miRNAs', model_name, calc_rsquared(all_df['pred_normed'], all_df['label_normed'])])

    return pd.DataFrame(all_r2s, columns=['mir', 'model', 'r2'])




def main():
    # # This is to test the gradient function to see if it is faster.
    # model = "biochem"
    # cell_line = "HeLa"
    # mirnas = "five"
    # kds = "predicted"
    # test = False
    # pars_ag = [-4.7] + 31*[0]
    # pars_model = [-1.898, 0.54]
    # pars = pars_ag + pars_model
    # biochemplus = False
    # mirnas_use = ["lsy-6", "miR-1", "miR-124"]

    # mirnas_use = ["lsy-6", "lsy-6_pass", "miR-1", "miR-1_pass",
    #               "miR-124", "miR-124_pass", "miR-137", "miR-137_pass",
    #               "miR-139", "miR-139_pass", "miR-143", "miR-143_pass",
    #               "miR-144", "miR-144_pass", "miR-153", "miR-153_pass",
    #               "miR-155", "miR-155_pass", "miR-182", "miR-182_pass",
    #               "miR-199a", "miR-199a_pass", "miR-204", "miR-204_pass",
    #               "miR-205", "miR-205_pass", "miR-216b", "miR-216b_pass",
    #               "miR-223", "miR-223_pass", "miR-7", "miR-7_pass"]


    # input_df = make_full_input_data_table(mirnas_use, kds, test=test)
    # # Make the data table with the transfection values.
    # trans_mc_df = make_transfection_data_table(cell_line, mirnas)
    # # Only use the rows in the input data dataframe for which there are
    # # transfection data.
    # input_df = input_df[input_df["transcript"].isin(trans_mc_df.index)]
    # input_df.reset_index(drop=True, inplace=True)

    # input_df_new = input_df[["transcript", "mir", "Threep", "logSA_diff", "in_ORF", "PCT", "mir_gp"]]

    # input_df_new.loc[:, "log_ka"] = -1*input_df["log_kd"]
    # # Get the Kathy-formatted miRNA labels so as to use her input tables.
    # mirnas_use_kl = [get_kathy_mirna_name(mirna) for mirna in mirnas_use]

    # ag_s = pars_ag
    # ag_mirna_map = pd.DataFrame.from_dict(
    #     {mirna : ag_i for mirna, ag_i in zip(mirnas_use_kl, ag_s)},
    #     orient="index"
    # )
    # check = np.array(ag_mirna_map.loc[input_df_new["mir"]])
    # print(check)

    # mirnas_col = np.array(input_df_new["mir"])
    # transcripts_col = np.array(input_df_new["transcript"])

    # starts = [0]
    # m_prev = mirnas_col[0]
    # t_prev = transcripts_col[0]
    # for i, (m, t) in enumerate(zip(mirnas_col, transcripts_col)):
    #     if m != m_prev or t != t_prev:
    #         starts.append(i)
    #     m_prev = m
    #     t_prev = t
    # stops = starts[1:] + [input_df_new.shape[0]]

    # limits = [i for i in zip(starts, stops)]

    # print(limits[:10])


    # input_df_np = input_df_new.to_numpy()
    # col_name_map = dict(zip(input_df_new.columns, list(range(len(input_df_new.columns)))))
    # time_1 = time.time()
    # out = get_predictions(pars, input_df, mirnas_use_kl, biochemplus=biochemplus)
    # time_2 = time.time()
    # out_new = get_predictions_final(pars, input_df_new, mirnas_use_kl, biochemplus=biochemplus)
    # time_3 = time.time()
    # out_3 = get_predictions_final(pars, input_df_new, mirnas_use_kl, biochemplus=biochemplus)
    # time_4 = time.time()
    # out_4 = get_predictions(pars, input_df, mirnas_use_kl, biochemplus=biochemplus)
    # time_5 = time.time()




    # print(("Old: %3.3f seconds." %(time_2 - time_1)))
    # print(("New: %3.3f seconds." %(time_3 - time_2)))
    # print(("Old: %3.3f seconds." %(time_5 - time_4)))
    # print(("New: %3.3f seconds." %(time_4 - time_3)))

    # print(out)
    # print(out_new)

    mirs5 = ['mir1','mir124','mir155','mir7','lsy6']

    nsites = pd.read_csv('/lab/bartel4_ata/kathyl/RNA_Seq/outputs/biochem/resubmission/nsites.txt', sep='\t', index_col=0)
    print(nsites.head())

    nsites_long = nsites.reset_index().melt(id_vars=['transcript'])
    nsites_long.columns = ['transcript','mir','utr_nsites_canon']
    nsites_long = nsites_long.set_index(['transcript','mir'])
    print(nsites_long.head())

    TS7_original_predictions = pd.read_csv('/lab/bartel4_ata/kathyl/RNA_Seq/outputs/ts7/no_xval_results/original_predictions.txt', sep='\t')
    TS7_multisite_predictions = pd.read_csv('/lab/bartel4_ata/kathyl/RNA_Seq/outputs/ts7/no_xval_results/multisite_predictions.txt', sep='\t')
    TS7_original_predictions = TS7_original_predictions.rename(columns={
        'pred_normed':'pred_normed_unbounded',
        'pred_bounded_normed': 'pred_normed'
    }).set_index(['transcript','mir'])
    TS7_original_predictions = pd.concat([TS7_original_predictions, nsites_long], axis=1, join='inner').sort_index()
    TS7_multisite_predictions = TS7_multisite_predictions.rename(columns={
        'pred_normed':'pred_normed_unbounded',
        'pred_bounded_normed': 'pred_normed'
    }).set_index(['transcript','mir'])
    TS7_multisite_predictions = pd.concat([TS7_multisite_predictions, nsites_long], axis=1, join='inner').sort_index()
    print(TS7_original_predictions.iloc[:20, :])

    temp0 = TS7_original_predictions.query('mir in @mirs5')
    temp = TS7_multisite_predictions.query('mir in @mirs5')
    # print(temp0.head())
    # print(temp.head())
    # print(temp0.shape)
    print(calc_rsquared(temp0['pred_normed'], temp0['label_normed']))
    print(calc_rsquared(temp['pred_normed'], temp['label_normed']))

if __name__ == "__main__":
    main()
