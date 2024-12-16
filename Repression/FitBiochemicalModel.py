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
imp.load_source("repression_functions",
                ("/lab/bartel1_ata/mcgeary/computation/"
                 "AgoRBNS/Repression/repression_functions.py"
                )
               )


from general import *
from RBNS_methods import *
from repression_functions import *
from sitetypes import get_seq_site_map


from scipy.optimize import minimize

pd.set_option("display.precision", 6)

# I first used the above address, but the bottom one gives near-perfect (but
# not perfect) agreement with Kathy's files, in the folder:
# /lab/solexa_bartel/klin/miRNA_models_data/transfections/hela/compiled/ .




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




def objective_function(pars, input_df, data_df, mirnas, biochemplus):
    # Sums the squared-difference of the prediction dataframe and the
    # transfection dataframe.
    pred_df = mean_center_df(get_predictions(pars, input_df, mirnas,
                                             biochemplus=biochemplus))
    N = pred_df.shape[0]*pred_df.shape[1]
    output = pred_df.sub(data_df).pow(2).sum().sum()/N
    global tick
    if tick%100 == 0:
        print(pars)
        print(output)
    tick += 1
    return(output)







def main():
    time_start = time.time()
    # Parse command line arguments:
    arguments = ["model", "cell_line", "mirnas", "kds", "-test_binary",
                 "-checkpars_binary"]
    args = parse_arguments(arguments)
    (model, cell_line, mirnas, kds, test, checkpars) = args

    # Make the features dictionary that has the various sequences and lengths
    # associated with the annotations of each mRNA.
    # features_mRNA_map = make_feature_dictionary(cell_line)

    # Scripting to generate the name of the file that contains the input data
    # for the biochemical model.
    if model == "biochemplus":
        biochemplus = True
    else:
        biochemplus = False
    if mirnas == "six":
        mirnas_use = ["let-7a", "lsy-6", "miR-1", "miR-124", "miR-155", "miR-7"]
    if mirnas == "five":
        mirnas_use = ["lsy-6", "miR-1", "miR-124", "miR-155", "miR-7"]
    elif cell_line == "HeLa" and mirnas == "sixteen" and kds == "predicted":
        mirnas_use = ["lsy-6", "lsy-6_pass", "miR-1", "miR-1_pass",
                      "miR-124", "miR-124_pass", "miR-137", "miR-137_pass",
                      "miR-139", "miR-139_pass", "miR-143", "miR-143_pass",
                      "miR-144", "miR-144_pass", "miR-153", "miR-153_pass",
                      "miR-155", "miR-155_pass", "miR-182", "miR-182_pass",
                      "miR-199a", "miR-199a_pass", "miR-204", "miR-204_pass",
                      "miR-205", "miR-205_pass", "miR-216b", "miR-216b_pass",
                      "miR-223", "miR-223_pass", "miR-7", "miR-7_pass"]
    input_df = make_full_input_data_table(mirnas_use, kds, test=test)
    print(input_df)
    # Make the data table with the transfection values.
    trans_mc_df = make_transfection_data_table(cell_line, mirnas)
    # Only use the rows in the input data dataframe for which there are
    # transfection data.
    input_df = input_df[input_df["transcript"].isin(trans_mc_df.index)]
    input_df.reset_index(drop=True, inplace=True)
    # Get the Kathy-formatted miRNA labels so as to use her input tables.
    mirnas_use_kl = [get_kathy_mirna_name(mirna) for mirna in mirnas_use]
    if cell_line == "HeLa":
        if mirnas == "six":
            pars_ag = [-5] + 5*[0]
        if mirnas == "five":
            pars_ag = [-5] + 4*[0]
        elif mirnas == "sixteen":
            pars_ag = [-3*(i%2) for i in range(32)]
            pars_ag[0] = -5
            # The miRNA guide sequences that are found in the file associated
            # with Kathy's optimization of the model with the neural-net
            # predicted Kd values.
            pars_ag = [-5.270323753356934, -8.370148658752441,
                       -5.458196640014648, -6.255456924438477,
                       -6.0075178146362305, -8.565059661865234,
                       -5.497497081756592, -5.923926830291748,
                       -5.677675724029541, -6.931764125823975,
                       -5.691253185272217, -6.096938133239746,
                       -5.80425500869751,  -6.091868877410889,
                       -6.030906677246094, -6.6293768882751465,
                       -5.477324962615967, -7.008443832397461,
                       -5.775873184204102, -6.593321800231934,
                       -5.8094611167907715, -6.7227020263671875,
                       -5.835589408874512, -6.089090824127197,
                       -5.693760395050049, -7.3324480056762695,
                       -5.310624122619629, -6.137135028839111,
                       -5.7454915046691895, -6.138209819793701,
                       -5.9393744468688965, -8.204909324645996]
            pars_ag[1:32] = [i - pars_ag[0] for i  in pars_ag[1:32]]
        if biochemplus:
            pars_model = [-2, 0.5, 0, 0, 0]
            pars_model = [0.5735143423080444, -1.7090857028961182, 0.11533742398023605, 0.09985063970088959, 0.9219825267791748]
        else:
            pars_model = [-2, 0.5]

        #     # pars = [-5.55907344,  0.62406557, -0.57499982,  0.30788662, -0.26999334,
        #     #         0.04727272, -1.93243033,  0.18722666,  0.2685592 ,  1.69639351]
        #     pars = [-5.562954425811768, -4.947486877441406, -6.119287490844727,
        #             -5.257762908935547, -5.826164722442627, -2, 0.2, 0.18821396,
        #             0.2381905, 1.4934731]
        #     # pars = -1*np.ones(10)
        # else:
        #     # pars = [-5, -5, -5, -5, -5, 0, 0]
        #     # pars = -1*np.ones(7)
        #     pars = [-5.224992752075195, -4.526782989501953, -5.244050979614258,
        #             -4.7737507820129395, -5.3175435066223145, 0.8655766248703003,
        #             -1.848806619644165]

    pars = pars_ag + pars_model


    ####### Perform the optimization. ##########################################

    global tick
    tick = 0

    print("About to try the minimization.")
    res = minimize(objective_function, pars, method="BFGS",
                   args=(input_df, trans_mc_df, mirnas_use_kl,
                         biochemplus))
    print("Done the minimization.")
    print(res)

    # Get the predictions based on the model parameters
    out_df = get_predictions(res.x, input_df, mirnas_use_kl, biochemplus=biochemplus)

    # Assign the optimized parameters, and convert the miRNA guide strands such
    # that they are all their own values (rather than being the difference
    # between themselves and the first one).
    pars_out = res.x
    pars_out[1:len(mirnas_use_kl)] = [i + pars_out[0] for i
                                      in pars_out[1:len(mirnas_use_kl)]]
    # Assign the names of each of the parameters.
    names_pars = ["ln_%s_ag" %i for i in mirnas_use] + ["ln_b", "ln_c_orf"]
    if biochemplus:
        names_pars += ["ln_c_SA", "ln_c_ThreeP", "ln_c_PCT"]
    print(names_pars)
    out_pars_df = pd.DataFrame.from_dict({i : j for (i, j)
                                          in zip(names_pars, pars_out)},
                                          orient="index")
    ##### OUTPUT ###############################################################
    # Assign the paths for both the predictions and the parameter files.
    path = "%s%s/count_tables/pred_%s_%s_%s.txt" %(
        base_exp_dir, cell_line, model, mirnas, kds
    )
    path_pars = "%s%s/model_parameters/pred_%s_%s_%s.txt" %(
        base_exp_dir, cell_line, model, mirnas, kds
    )

    # Print both the predictions and the parameter files to their respective
    # paths.
    out_df.to_csv(path, sep="\t", index_label=False)
    out_pars_df.to_csv(path_pars, sep="\t", index_label=False, header=False)
    # Done with the script.
    time_done = time.time()
    print(("Finished job future in %3.2f seconds." %(time_done - time_start)))
    return()

if __name__ == "__main__":
    main()
