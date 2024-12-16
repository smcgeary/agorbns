################################################################################
# sean_fit_ts_params.py                                                        #
################################################################################
# This is where I try to re-implement the fitting of the target scan model using
# Kathy's python notebooks.
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

# import time

# import numpy as np
# import pandas as pd

import statsmodels.api as sm
import statsmodels.formula.api as smf

from scipy.optimize import minimize
# from collections import defaultdict

tick = 0

sites = ["6mer", "7mer-1a", "7mer-m8", "8mer-1a"]


mirs16 = sorted(['mir1','mir124','mir155','mir7','lsy6','mir153','mir139','mir144','mir223','mir137',
            'mir205','mir143','mir182','mir199a','mir204','mir216b'])

upper_bound_dict = {
    '8mer-1a': -0.03,
    '7mer-m8': -0.02,
    '7mer-1a': -0.01,
    '6mer': 0.0
}

feature_list = ['TA','SPS','SA', #1–3
                'siRNA_1A', 'siRNA_1C', 'siRNA_1G', #4
                'siRNA_8A', 'siRNA_8C', 'siRNA_8G', #5
                'site_8A', 'site_8G', 'site_8C', #6
                'Local_AU_score', 'Threep_score', 'Min_dist_score', #7–9
                'UTR_length_score','ORF_length', 'ORF_8mers', 'Off6m_score', #10–13
                'PCT' # 14
               ]

scaled_feature_list = ["Threep_score", "SPS", "TA", "UTR_length_score", #1–4
                       "ORF_length", "Min_dist_score", "Local_AU_score", #5–7
                       "SA", "PCT"] #8–9


def mean_center(df):
    return(df.sub(df.mean(axis=1), axis=0))

def get_tpm_data(linear=True, mean_centered=False):
    TPM_FILE = '/lab/bartel4_ata/kathyl/RNA_Seq/outputs/biochem/final_fixed/merged.txt'
    # Read in the tpm file.
    transXmir_tpm = pd.read_csv(TPM_FILE, sep='\t', index_col=0)
    # Names the index names as transcript.
    transXmir_tpm.index.name = 'transcript'
    transXmir_tpm = transXmir_tpm.loc[:, transXmir_tpm.columns.isin(mirs16)]
    transXmir_tpm.sort_index(axis=1, inplace=True)
    # Potentially mean-center the matrix.
    if mean_centered:
        transXmir_tpm = mean_center(transXmir_tpm)
    # Potentially linearize the matrix.
    if linear:
        vec_tpm = transXmir_tpm.reset_index().melt(id_vars=['transcript'])
        vec_tpm = vec_tpm.rename(columns={'variable': 'mir', 'value': 'label'})
        vec_tpm = vec_tpm.set_index(keys=['transcript','mir']).sort_index()
        out = vec_tpm
    else:
        out = transXmir_tpm
    # Return the matrix.
    return(out)


def expand_feats_stypes(features, stypes, expand_vars, single_vars):
    expanded_features = []
    for stype in stypes:
        temp = features[expand_vars]
        stype_filter = features[[stype]].values
        temp.columns = [x + '_' + stype for x in temp.columns]
        temp *= stype_filter
        expanded_features.append(temp)

    expanded_features.append(features[single_vars])
    expanded_features = pd.concat(expanded_features, axis=1, join='inner')

    # get rid of columns of all zeros, for example 6mer PCT
    for col in expanded_features.columns:
        if np.std(expanded_features[col].values) < 0.00001:
            expanded_features = expanded_features.drop(columns=[col])

    return expanded_features

def rescale_features(feat_df, scale_df):
    # Convert the indeces to being numeric (otherwise the ind_use command
    # doesn't work).
    feat_df.reset_index(level=[0,1], inplace=True)
    for site in sites:
        # Select only the rows associated with that site type.
        inds_use = feat_df[feat_df[f'Stype_{site}']==1].index
        # Select only those columns associated with features corresponding to
        # that site type.
        cols_use = ["%s_%s" %(feat_i, site) for feat_i in scaled_feature_list]
        # Make sure the column is in the feature matrix.
        cols_use = [i for i in cols_use if i in feat_df.columns.values]
        # Iterate over the columns.
        for col_use_i in cols_use:
            # Grab the values.
            vals = feat_df.loc[inds_use, col_use_i]
            # Determine the 5th and 95th percentiles.
            coef_min = np.percentile(vals, 5)
            coef_max = np.percentile(vals, 95)
            # Identify the row to use.
            scale_row_use = col_use_i.replace("_" + site, "")
            # Update the scale matrix with percentiles.
            scale_df.loc[scale_row_use, site + " min"] = coef_min
            scale_df.loc[scale_row_use, site + " max"] = coef_max
            # re-scale the values and use them to replace the values in the
            # original feature matrix.
            vals = (vals - coef_min)/(coef_max - coef_min)
            feat_df.loc[inds_use, col_use_i] = vals
    # Convert the indeces back into the multi-index type.
    feat_df = feat_df.set_index(keys=['transcript','mir']).sort_index()
    return(feat_df, scale_df)



def make_design_and_scale_dfs(transcripts_use=[], aggregate=True, upper_bound=False,
                       rescale_feat=False):
    feat_df = pd.read_csv(
        ('/lab/bartel4_ata/kathyl/RNA_Seq/data/TS7_output_files/'
         'refseq_utr90_fixed_ribosome_shadow/features.txt'), 
        sep='\t'
    )
    # Rename "miRNA_family" column as "mir".
    feat_df.rename(columns={'miRNA family': 'mir'}, inplace=True)
    # Swap " " for "_" in column names.
    feat_df.columns = [x.replace(' ','_') for x in feat_df.columns]
    # Rename "Gene_ID" column as "transcript".
    feat_df.rename(columns={'Gene_ID': 'transcript'}, inplace=True)
    # Selects on those rows that have miRNAs within the sixteen HeLa miRNAs.
    # NOTE, this will have to be changed in order to allow testing on the test
    # set.
    feat_df = feat_df[feat_df['mir'].isin(mirs16)]
    # Makes sure transcript is in the TPM data.
    if len(transcripts_use) != 0:
        feat_df = feat_df[feat_df['transcript'].isin(transcripts_use)]
    # Sort the rows by transcript and then miRNA.
    feat_df = feat_df.set_index(keys=['transcript','mir']).sort_index()
    # Define the site types.
    global ts7_stypes
    ts7_stypes = list(feat_df['Site_type'].unique())
    # Rename the site-type columns to "Site_type_6mer" from "6mer".
    for stype in ts7_stypes:
        feat_df[stype] = (feat_df['Site_type'] == stype).astype(float)
    single_vars = ts7_stypes.copy()
    if upper_bound:
        feat_df["upper_bound"] = [upper_bound_dict[x] for x in feat_df['Site_type']]
        single_vars += ["upper_bound"]
    # Expand the matrix by constructing one column for each feature-and-site-
    # type pair.
    feat_df = expand_feats_stypes(feat_df, ts7_stypes, feature_list,
                                  single_vars)
    # Rename the site coefficient columns (e.g., rename "6mer" to "Stype_6mer",
    # etc.).
    feat_df.rename(columns={x:'Stype_{}'.format(x) for x in ts7_stypes}, inplace=True)
    # Define the scaling matrix, which gives the min and max values for
    # rescaling each of the 9 features which can be rescaled.
    scale_df = pd.DataFrame(0, index=feature_list,
        columns=[j for i in ["min", "max"] for j in ["%s %s" %(site, i) for site in sites]])
        # columns=[j for i in sites for j in [i + " min", i + " max"]])
    # Allocate 1 to all of the maximum values.
    scale_df.loc[:, scale_df.columns.str.endswith('max')] = 1
    # Optionally rescale the features.
    if rescale_feat:
        feat_df, scale_df = rescale_features(feat_df, scale_df)
    # Replace dashes with no spaces.
    feat_df.columns = [x.replace('-','') for x in feat_df.columns]
    feat_df = feat_df.astype(np.float64)
    if aggregate:
        feat_df = feat_df.groupby(['transcript','mir']).agg(np.sum)
    return(feat_df, scale_df)


def fit_original(params = [], rescale_feat=False):
    tpms = get_tpm_data(linear=True, mean_centered=False)
    # Get the list of transcripts to use.
    transcripts_use = list(set(tpms.index.get_level_values(0).values))
    # Get design matrix
    design_df, scale_df = make_design_and_scale_dfs(transcripts_use=transcripts_use,
                                   upper_bound=False, aggregate=True,
                                   rescale_feat=rescale_feat)
    design_df = design_df.reindex(tpms.index).fillna(0.0)
    # Get features
    all_features = design_df.columns
    # Remove the labels (do I need to do this? [I do, to get transcript index.])
    tpms.reset_index(inplace=True)
    design_df.reset_index(inplace=True)
    # Define the list of features to use in the training.
    train_feats = [i for i in all_features if np.std(design_df[i].values) != 0]
    print("Training the model.")
    formula = 'label ~ {} - 1'.format(' + '.join(train_feats))
    # Add the data ('label') into the design matrix.
    design_df["label"] = tpms["label"]
    # Define the model.
    if len(params) == 0:
        mod = smf.mixedlm(formula=formula, data=design_df,
                          groups=design_df['transcript'])
        # Fit the model
        mdf = mod.fit()
        # Assign the parameters to an output dataframe.
        params = pd.DataFrame(mdf.params).drop(index=['Group Var'])
    return(params, scale_df)

def fit_with_bounds(params=[], rescale_feat=False):
    tpms_mc = get_tpm_data(linear=False, mean_centered=True)
    # Get the list of transcripts to use.
    transcripts_use = list(tpms_mc.index.values)
    # Get design matrix
    design_df, scale_df = make_design_and_scale_dfs(
        transcripts_use=transcripts_use, upper_bound=True, aggregate=False
    )
    # Split the design matrix into an upper bounds vector and the rest of the
    # design matrix.
    upper_bounds = design_df["upper_bound"]
    design_df.drop(columns=['upper_bound'], inplace=True)
    # Get features
    all_features = design_df.columns
    train_feats = [i for i in all_features if np.std(design_df[i].values) != 0]
    # Assign the final design matrix.
    design_df = design_df[train_feats]
    # Initialize the parameters.
    if len(params) == 0:
        params = pd.DataFrame(-1, index=train_feats)
    def objective_function(pars):
        # Get the predicted effect of each site, using the upper bounds.
        preds = pd.concat([design_df.dot(pars), upper_bounds], axis=1).min(axis=1)
        # Sum over the sites within each transcript.
        preds = preds.groupby(['transcript','mir']).agg(np.sum).reset_index(level=[0,1])
        # Reset the index and convert the matrix into a mirna-by-transcript
        # matrix.
        preds = preds.pivot(index="transcript", columns="mir", values=0).fillna(0)
        preds = preds.reindex(tpms_mc.index).fillna(0.0)
        preds_mc = mean_center(preds)
        obj = preds_mc.sub(tpms_mc).pow(2).sum().sum()
        global tick
        if tick%1000 == 0:
            print(pars)
            print(obj)
        tick += 1

        return(obj)
    params = params.iloc[0:design_df.shape[1], :]
    res = minimize(objective_function, params, method="BFGS")
    # res = minimize(objective_function, params, options={"maxiter" : 1})
    params = res.x
    params = pd.DataFrame(params)
    params.index = train_feats
    return(params, scale_df)


def main():
    # Define the miRNAs.
    time_start = time.time()
    # Parse command line arguments:
    arguments = ["-bound_binary", "-rescale_binary", "-alt_binary"]
    args = parse_arguments(arguments)
    (bound, rescale, alt) = args
    if bound:
        bound_str = "_bounded"
        alt_str = ""
    elif alt:
        bound_str = ""
        alt_str ="_alt"
    else:
        bound_str = ""
        alt_str = ""
    if rescale:
        rescale_str = "_rescale"
    else:
        rescale_str = ""

    # The base fit, used to initialize the parameters.
    params, scale_df = fit_original(rescale_feat=rescale)
    # If the upper bound flag is chosen, this refits with the new objective
    # function that utilizes the upper bounds.
    if bound:
        params, scale_df = fit_with_bounds(params, rescale_feat=rescale)
    output_df = pd.DataFrame(
        0, 
        index=["Intercept"] + feature_list,
        columns= ["%s %s" %(j, i) for i in ["coeff", "min", "max"]
                  for j in ["6mer", "7mer-m8", "7mer-1a", "8mer-1a"]]
    )
    # Name the indeces of the output matrix.
    output_df.index.name = "Feature"
    # Identify the indeces and columns to use to add the scale_df values into
    # the output_df.
    scale_ind = scale_df.index.values
    scale_col = scale_df.columns.values
    # Add the scale values.
    output_df.loc[scale_ind, scale_col] = scale_df.loc[scale_ind, scale_col]
    # Iterate over the sites.
    for site in sites:
        # First, assign the intercept coefficients:
        output_col = "%s coeff" %site
        params_ind = "Stype_%s" %site.replace("-", "")
        output_df.loc["Intercept", output_col] = params.loc[params_ind, 0]
        for feat in feature_list:
            params_ind = "%s_%s" %(feat, site.replace("-", ""))
            if params_ind in params.index.values:
                output_df.loc[feat, output_col] = params.loc[params_ind, 0]
    path_1 = "/lab/solexa_bartel/mcgeary/transfections/HeLa/target_scan/params/"
    path_2 = "params%s%s%s.txt" %(bound_str, alt_str, rescale_str)
    path = path_1 + path_2
    print(path)
    output_df.to_csv(path, sep='\t')

    # # print(params.loc["Stype_6mer"])
    # # print(params.loc[["Stype_6mer", "Stype_7mer1a"]])
    # # print(["Stype_%s" %i.replace("-", "") for i in sites])
    # # vals = params.loc[["Stype_%s" %i.replace("-", "") for i in sites]]
    # # print(vals)
    # # output_df.loc["Intercept", ["%s coeff" %i for i in sites]] = vals
    # # print(output_df)
    # return()


    # # del pars_in.index.name



    # time_done = time.time()
    # print(("Finished first optimization in %3.2f seconds." %(time_done - time_start)))


    # params_new = fit_with_bounds(params)
    # print(params_new)
    # params_new = pd.DataFrame(params_new)

    # time_done = time.time()
    # print(("Finished second optimization in %3.2f seconds." %(time_done - time_start)))


    # params_for_output = []
    # row = ['Intercept']
    # print(row)
    # for stype in ts7_stypes:
    #     stype = stype.replace('-','')
    #     row.append(params_new.loc[f'Stype_{stype}'][0])
    #     print(row)
    # row += [0] * 8
    # print(row)
    # params_for_output.append(row)
    # print(params_for_output)
    # for feature in feature_list:
    #     row = [feature]
    #     for stype in ts7_stypes:
    #         stype = stype.replace('-','')
    #         try:
    #             row.append(params_new.loc[f'{feature}_{stype}'][0])
    #         except:
    #             print(feature, stype)
    #             row.append(0)
    #     row += [0] * 4 + [1] * 4
    #     params_for_output.append(row)
    # cols = ['Feature'] + [f'{stype} coeff' for stype in ts7_stypes] + [f'{stype} min' for stype in ts7_stypes] + [f'{stype} max' for stype in ts7_stypes]
    # params_for_output = pd.DataFrame(params_for_output, columns=cols).set_index('Feature')
    # params_for_output.to_csv('Repression/bounded_parameters.txt', sep='\t')
    # params_for_output.head()


    # # params_for_output.to_csv("Repression/new_parameters_mean_centered.txt", sep="\t")



    # return()
    # # # Makes the prediction file that is as large as the data file.
    # # canon_features_all_ix = canon_features_expanded_agg.reindex(tpm_long.index).fillna(0.0)
    # # canon_features_all_ix = pd.concat([tpm_long, canon_features_all_ix], axis=1)
    # # canon_features_all_ix['label'] = canon_features_all_ix['label'].astype(np.float64)
    # # print("canon_features_all_ix.head()")
    # # print(canon_features_all_ix.head())



    # # return()


    # # Part that seems to do the actual fitting.    

    # design_matrix = canon_features_all_ix[train_feats]
    # tpm_dim = tpm_long_df_norm.shape
    # pars_init = pd.DataFrame([1]*73)
    # pars_init.index = design_matrix.columns
    # tick = 0
    # # def loss_fun(pars):
    # def objective_function(pars):
    #     preds = design_matrix.dot(pars).reset_index(level=[0,1])
    #     pred_df = preds.pivot(index="transcript", columns="mir").fillna(0)
    #     pred_df_norm = pred_df.sub(pred_df.mean(axis=1), axis=0)
    #     output = pred_df_norm.sub(tpm_long_df_norm).pow(2).sum().sum()
    #     global tick
    #     if tick%100 == 0:
    #         print(pars)
    #         print(output)
    #     tick += 1

    #     return(output)

    # def bounded_objective_function(pars):
    #     preds_1 = pd.concat([design_matrix_bound.dot(pars), upper_bounds],
    #                         axis=1)
    #     return(preds_1)

    # print(design_matrix_bound.head(15))
    # print(upper_bounds.head(15))

    # preds = bounded_objective_function(pars_init)    
    # print(preds.head(15))





    # return()




    # # Makes the prediction file that is as large as the data file.
    # canon_features_all_ix = canon_features_expanded_agg.reindex(tpm_long.index).fillna(0.0)



    # print(objective_function(pars_init))

    # # Alternative way of fitting the data, in which the data are mean-centered.
    # # res = minimize(objective_function, pars_init, method="BFGS")
    # # pars = res.x
    # # pars_out = {i: j for (i, j) in zip(design_matrix.columns.values, pars)}
    # # pars_out = pd.DataFrame.from_dict(pars_out, orient="index")
    # # print(pars_out)
    # # pars_out.to_csv("Repression/brand_new_pars.txt", sep="\t", index_label=False, header=False)



    # formula2 = 'label ~ {}'.format(' + '.join(train_feats))

    # mod = smf.mixedlm(formula=formula, data=train, groups=train['transcript'])
    # print(mod)
    # mdf = mod.fit()
    # print("Finished training the model.")

    # params = pd.DataFrame(mdf.params)
    # print(params)

    # mod2 = smf.mixedlm(formula=formula2, data=train, groups=train['transcript'])
    # print(mod2)
    # mdf2 = mod2.fit()
    # params2 = pd.DataFrame(mdf2.params)

    # print(params)
    # print(params2)
    # return()


    # params_alt = pars_in
    # print(params_alt)
    # params = params_alt
    # params_for_output = []
    # row = ['Intercept']
    # print(row)
    # for stype in ts7_stypes:
    #     stype = stype.replace('-','')
    #     row.append(params.loc[f'Stype_{stype}'][0])
    #     print(row)
    # row += [0] * 8
    # print(row)
    # params_for_output.append(row)
    # print(params_for_output)
    # for feature in feature_list:
    #     row = [feature]
    #     for stype in ts7_stypes:
    #         stype = stype.replace('-','')
    #         try:
    #             row.append(params.loc[f'{feature}_{stype}'][0])
    #         except:
    #             print(feature, stype)
    #             row.append(0)
    #     row += [0] * 4 + [1] * 4
    #     params_for_output.append(row)
    # cols = ['Feature'] + [f'{stype} coeff' for stype in ts7_stypes] + [f'{stype} min' for stype in ts7_stypes] + [f'{stype} max' for stype in ts7_stypes]
    # params_for_output = pd.DataFrame(params_for_output, columns=cols).set_index('Feature')
    # params_for_output.to_csv('Repression/new_parameters.txt', sep='\t')
    # params_for_output.head()


    # params_for_output = []
    # row = ['Intercept']
    # print(row)
    # for stype in ts7_stypes:
    #     stype = stype.replace('-','')
    #     row.append(params.loc[f'Stype_{stype}'][0])
    #     print(row)
    # row += [0] * 8
    # print(row)
    # params_for_output.append(row)
    # print(params_for_output)
    # for feature in feature_list:
    #     row = [feature]
    #     for stype in ts7_stypes:
    #         stype = stype.replace('-','')
    #         try:
    #             row.append(params.loc[f'{feature}_{stype}'][0])
    #         except:
    #             print(feature, stype)
    #             row.append(0)
    #     row += [0] * 4 + [1] * 4
    #     params_for_output.append(row)
    # cols = ['Feature'] + [f'{stype} coeff' for stype in ts7_stypes] + [f'{stype} min' for stype in ts7_stypes] + [f'{stype} max' for stype in ts7_stypes]
    # params_for_output = pd.DataFrame(params_for_output, columns=cols).set_index('Feature')
    # params_for_output.to_csv('Repression/new_parameters.txt', sep='\t')
    # params_for_output.head()


    # params_for_output.to_csv("Repression/new_parameters_mean_centered.txt", sep="\t")



    # time_done = time.time()
    # print(("Finished in %3.2f seconds." %(time_done - time_start)))
    # return()


if __name__ == "__main__":
    main()
