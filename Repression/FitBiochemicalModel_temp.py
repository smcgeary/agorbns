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

# import statsmodels.formula.api as smf
# import regex


from scipy.optimize import minimize

# import bisect
# import json

# import numpy.ma as ma

# import utils





pd.set_option("display.precision", 6)

# I first used the above address, but the bottom one gives near-perfect (but
# not perfect) agreement with Kathy's files, in the folder:
# /lab/solexa_bartel/klin/miRNA_models_data/transfections/hela/compiled/ .

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





# def make_feature_dictionary(cell_line):
#     transcript_file_dir = ("/lab/solexa_bartel/mcgeary/transfections"
#                            "/metadata/%s_transcripts_2.txt" %(cell_line))
#     # Make a dictionary that has the information on each transcript in order to
#     # ensure I am performing the analysis the same way that Kathy was doing it.
#     with open(transcript_file_dir) as transcript_file:
#         # Make the outer layer of keys.
#         outer_keys = transcript_file.readline().strip().split("\t")[1:]
#         features_mRNA_map = {i : collections.defaultdict(lambda: "Key does not exist") for i in outer_keys}
#         line = transcript_file.readline()
#         while line:
#             line_list = line.strip().split("\t")
#             transcript_id = line_list[0]
#             features = line_list[1:]
#             for i, item in enumerate(features):
#                 features_mRNA_map[outer_keys[i]][transcript_id] = item
#             line = transcript_file.readline()
#     return(features_mRNA_map)

def mean_center_df(df):
    return(df.sub(df.mean(axis=1), axis=0))




def get_kathy_mirna_name(mirna):
    atypical_name_dict = {"let-7a" : "let7" , "lsy-6" : "lsy6"}
    if mirna in atypical_name_dict.keys():
        return(atypical_name_dict[mirna])
    else:
        return "ir".join(mirna.split("iR-"))

def get_sean_mirna_name(mirna):
    atypical_name_dict = {"let7" : "let-7a" , "lsy6" : "lsy-6"}
    if mirna in atypical_name_dict.keys():
        return(atypical_name_dict[mirna])
    else:
        return "iR-".join(mirna.split("ir"))



def make_input_data_table(mirna, kds, threep_canon_only=True, test=False):
    mirna_file = get_kathy_mirna_name(mirna) + ".txt"
    if kds == "measured":
        path = kds_base_dir + measured_kds_dir + mirna_file
    elif kds == "predicted":
        path = kds_base_dir + predicted_kds_dir + mirna_file
    input_df = pd.read_csv(path, sep="\t")
    # Apply the Kd cutoff to the data:
    input_df = input_df[input_df["log_kd"] < 0].reset_index(drop=True)
    # Give all noncanonical sites an structural score equal to the mean over
    # all of the canonical sites.
    logSA_mean = np.nanmean(input_df["logSA_diff"])
    input_df["logSA_diff"].fillna(logSA_mean, inplace=True)
        # Assign a threep score of zero to all noncanonical sites.
    if threep_canon_only:
        input_df.loc[(input_df["stype"] == "no site"), "Threep"] = 0
    if test:
        input_df = input_df.iloc[:50, ]
    return(input_df)

def make_transfection_data_table(cell_line):
    # Path to transfection data.
    input_data_path = ("/lab/solexa_bartel/mcgeary/transfections/%s/"
                       "count_tables/logtpm_batchnormalized.txt" %cell_line)
    #Input data.
    data_df = pd.read_csv(input_data_path, sep="\t", index_col=0)
    # Subset the data (will make this an option later on, now just trying
    # Figure 5).
    data_df = data_df[["lsy-6", "miR-1", "miR-124", "miR-155", "miR-7"]]
    # Mean center these predictions:
    # data_mean_centered_df = data_df.sub(data_df.mean(axis=1), axis=0)
    return(mean_center_df(data_df))




def get_predictions(pars, in_df, mirnas_use, biochemplus=False):
    # Assign the parameters:
    print(pars)
    print(biochemplus)
    print(mirnas_use)
    # Assign the free miRNA concentrations.
    # log_ag_s = pars[0] + np.append([0], pars[1:len(mirnas_use)])
    ag_s = np.exp(pars[0] + np.append([0], pars[1:len(mirnas_use)]))
    print(ag_s)
    # ag_s = np.exp(pars[:len(mirnas_use)])
    # Assign the repression and orf-site penalty parameters.
    b, c = np.exp(pars[len(mirnas_use):len(mirnas_use) + 2])
    print(b)
    print(c)
    # Assign the feature weights with either the input parameters for the
    # biochemplus model, or  0 for the biochem model.
    if biochemplus:
        cSA, cThrP, cPCT = pars[len(mirnas_use) + 2:]
    else:
        cSA, cThrP, cPCT = [0, 0, 0]
    print(cSA)
    print(cThrP)
    print(cPCT)
    out_df = in_df.copy()
    ag_mirna_map = pd.DataFrame.from_dict(
        {mirna : ag_i for mirna, ag_i in zip(mirnas_use, ag_s)},
        orient="index"
    )
    # Calculate the Kd offset from the features (will be 0 for model="biochem")
    out_df["offset"] = cSA*in_df["logSA_diff"] + cThrP*in_df["Threep"] + cPCT*in_df["PCT"]
    print(out_df)
    # Subtract offset from log_kd
    out_df["log_kd"] = in_df["log_kd"] - out_df["offset"]
    out_df["log_kd_bg"] = -1*out_df["offset"]
    print(out_df)
    out_df["a_g"] = ag_mirna_map.loc[out_df["mir"]].reset_index(drop=True)
    print(out_df)
    # Get the specifically bound sites:
    out_df["N_g"] = out_df["a_g"]/(out_df["a_g"] + np.exp(out_df["log_kd"])*c**(-1*in_df["in_ORF"]))
    # Get the background binding:
    out_df["N_g,bg"] = out_df["a_g"]/(out_df["a_g"] + np.exp(out_df["log_kd_bg"])*c**(-1*in_df["in_ORF"]))

    # Sum the N_g and N_g,bg values across each transcript-and-miRNA combination.
    N_full = out_df.groupby(["transcript", "mir"]).agg({'N_g': np.sum,
                                                        'N_g,bg': np.sum})
    print(N_full)
    N_full.reset_index(level=[0,1], inplace=True)
    print(N_full)
    return()
    # Convert these occupancies into the final predicted repression values:
    N_full["pred"] = np.log((1 + b*N_full["N_g,bg"])/(1 + b*N_full["N_g"]))
    # Convert data table back into a matrix with dimensions tranxcript X miRNA,
    # and replace the NaN values (which are transcripts for which there were no
    # 12 mers against a particular miRNA) with 0, both the background binding
    # and specific binding would be the same in this case.
    pred_df = N_full.pivot(index="transcript",
                           columns="mir",
                           values="pred").fillna(0)
    # Change column names to be consistent with data_df.
    pred_df.columns = [get_sean_mirna_name(mirna) for mirna 
                       in pred_df.columns.tolist()]
    return(pred_df)


def objective_function(pars, input_df, data_df, mirnas, biochemplus):
    # Sums the squared-difference of the prediction dataframe and the
    # transfection dataframe.
    print(biochemplus)
    pred_df = mean_center_df(get_predictions(pars, input_df, mirnas,
                                             biochemplus=biochemplus))
    N = pred_df.shape[0]*pred_df.shape[1]
    output = pred_df.sub(data_df).pow(2).sum().sum()/N
    print(output)
    return(output)



def main():
    time_start = time.time()
    # Parse command line arguments:
    arguments = ["model", "cell_line", "mirnas", "kds", "-test_binary"]
    args = parse_arguments(arguments)
    (model, cell_line, mirnas, kds, test) = args

    # Make the features dictionary that has the various sequences and lengths
    # associated with the annotations of each mRNA.
    # features_mRNA_map = make_feature_dictionary(cell_line)

    # Scripting to generate the name of the file that contains the input data
    # for the biochemical model.
    if model == "biochemplus":
        biochemplus = True
    else:
        biochemplus = False
    if mirnas == "five":
        mirnas_use = ["lsy-6", "miR-1", "miR-124", "miR-155", "miR-7"]
    input_df = pd.concat([make_input_data_table(mirna, kds, test=test)
                          for mirna in mirnas_use],
                         axis=0)
    test_bool = False
    print(input_df)
    # Make the data table with the transfection values.
    data_mean_centered_df = make_transfection_data_table(cell_line)

    input_df = input_df[input_df["transcript"].isin(data_mean_centered_df.index)]

    input_df.index = [i for i in range(input_df.shape[0])]


    ag_mirna_map = pd.DataFrame.from_dict(
        {"lsy6"   : -5.224992752075195,
         "mir1"   : -4.526782989501953,
         "mir124" : -5.244050979614258,
         "mir155" : -4.7737507820129395,
         "mir7"   : -5.3175435066223145},
         orient="index"
     )
    mirnas_use_kl = [get_kathy_mirna_name(mirna) for mirna in mirnas_use]

    pars = [-5.224992752075195, -4.526782989501953, -5.244050979614258,
            -4.7737507820129395, -5.3175435066223145, 0.8655766248703003,
            -1.848806619644165]

    if biochemplus:
        pars = [-5, -5, -5, -5, -5, 0, 0, 1, 2, 3]
        # pars = [-5.55907344,  0.62406557, -0.57499982,  0.30788662, -0.26999334,
        #         0.04727272, -1.93243033,  0.18722666,  0.2685592 ,  1.69639351]
        pars = [-5.562954425811768, -4.947486877441406, -6.119287490844727,
                -5.257762908935547, -5.826164722442627, 0.2, -2, 0.18821396,
                0.2381905, 1.4934731]
        # pars = -1*np.ones(10)
    else:
        # pars = [-5, -5, -5, -5, -5, 0, 0]
        # pars = -1*np.ones(7)
        pars = [-5.224992752075195, -4.526782989501953, -5.244050979614258,
                -4.7737507820129395, -5.3175435066223145, 0.8655766248703003,
                -1.848806619644165]
    print(pars)
    # pars[1:5] = [i - pars[0] for i in pars[1:5]]
    print(pars)
    pars = [-4.947486877441406, 0.5735143423080444, -1.7090857028961182, 0.11533742398023605, 0.09985063970088959, 0.9219825267791748]
    input_df = input_df.loc[input_df["mir"] == "mir1", ]
    mirnas_use_kl = ["mir1"]
    print(input_df)
    print(pars)
    preds = get_predictions(pars, input_df, mirnas_use_kl,
                            biochemplus=biochemplus)
    return()
    loss = objective_function(pars, input_df, data_mean_centered_df,
                              mirnas_use_kl, biochemplus=biochemplus)
    # print(loss)
    # # print(preds)
    # # print(ag_mirna_map)
    # # input_df = input_df.iloc[:15, ]
    # # print(input_df)
    # # print(input_df["mir"])
    # # print(np.exp(ag_mirna_map.loc[input_df["mir"]].reset_index(drop=True)))
    # # input_df["ag"] = np.exp(ag_mirna_map.loc[input_df["mir"]].reset_index(drop=True))
    # # print(input_df)
    # return()
    # ag_df = pd.DataFrame.from_dict(ag_mirna_map, orient="index")
    # ag_df_list = input_df["mir"].tolist()
    # ag_input = np.exp(ag_df.loc[ag_df_list])
    # ag_input.index = input_df.index

    # input_df["ag"] = ag_input


    # pars = [0.8655766248703003, -1.848806619644165]

    # pars = rep()
    print("About to try the minimization.")
    res = minimize(objective_function, pars, method="BFGS",
                   args=(input_df, data_mean_centered_df, mirnas_use_kl,
                         biochemplus))
    print("Done the minimization.")
    print(res)

    # Get the predictions based on the model parameters
    out_df = get_predictions(res.x, input_df, mirnas_use_kl, biochemplus=biochemplus)
    print(out_df.iloc[:5,])
    ##### OUTPUT ##########################
    path = "%s%s/count_tables/pred_%s_%s_%s_new.txt" %(base_exp_dir, cell_line,
                                                       model, mirnas, kds)
    print(path)
    out_df.to_csv(path, sep="\t", index_label=False)


    time_done = time.time()
    print(("Finished job future in %3.2f seconds." %(time_done - time_start)))
    return()





if __name__ == "__main__":
    main()
