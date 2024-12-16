################################################################################
# sean_fit_ts_params.py                                                        #
################################################################################
# This is where I try to re-implement the fitting of the target scan model using
# Kathy's python notebooks.
import time

import numpy as np
import pandas as pd

import statsmodels.api as sm
import statsmodels.formula.api as smf

from scipy.optimize import minimize

import matplotlib.pyplot as plt

tick = 0

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

def make_design_matrix(transcripts_use=[], aggregate=True, upper_bound=False):
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
    # print(feat_df.iloc[0:30, 0:8])
    # print(feat_df.iloc[0:30, 8:16])
    # print(upper_bound_vec.iloc[0:30, :])
    # Here the original matrix function is replaced with a new matrix in which
    # each of the columns is one of the possible features of targetscan. There
    # are 14 features in addition to which site type it is, which means that
    # there is an expectation that there would be 15*4 = 60 columns. There are
    # actually more because the features of siRNA nucleotide 1, 8, and site
    # nucleotide 8 are each 3 columns each. Lastly, there is a upper bound
    # column, which gives the maximum amount of repression that the row can
    # have.
    feat_df = expand_feats_stypes(feat_df, ts7_stypes, feature_list,
                                  single_vars)
    # Rename the site coefficient columns (e.g., rename "6mer" to "Stype_6mer",
    # etc.).
    feat_df.rename(columns={x:'Stype_{}'.format(x) for x in ts7_stypes}, inplace=True)
    # Replace dashes with no spaces.
    feat_df.columns = [x.replace('-','') for x in feat_df.columns]
    print(feat_df.head(n=20))
    return()
    feat_df = feat_df.astype(np.float64)
    if aggregate:
        feat_df = feat_df.groupby(['transcript','mir']).agg(np.sum)
    return(feat_df)
    # print(canon_features_stype_expanded.head())
    # Rename the columns "6mer", "7mer-1a", etc., to "Stype_6mer", "Stype_7mer-1a", etc.
    # canon_features_stype_expanded = canon_features_stype_expanded.rename(columns={x:'Stype_{}'.format(x) for x in ts7_stypes})
    # canon_features_stype_expanded.columns = [x.replace('-','') for x in canon_features_stype_expanded.columns]
    # Add the batch numbers to each row.
    # canon_features_stype_expanded['batch'] = all_tpms.reindex(canon_features_stype_expanded.reset_index()['transcript'].values)['batch'].values
    # Make a new object, that drops the batches and the upper bounds, and
    # aggregates multiple sites per mRNA.

    ## TO DO, RE-SCALE THE FEATURES SO THAT THEY ARE BOUNDED FROM ZERO TO 1.
    # canon_features_norm = rescale_features(canon_features_stype_expanded)
    # return()



    # all_features = canon_features_expanded_agg.columns
    # canon_features_expanded_agg.columns = all_features
    # print("Post the feature column reassignment")
    # print(canon_features_expanded_agg.head())

    # canon_features_expanded_agg = canon_features_expanded_agg.astype(np.float64)
    # print("Post the np.float64 assignment")
    # print(canon_features_expanded_agg.head())

    # SEAN SPECIFIC: Make modified design matrix and upper bound matrices, to 
    # perform optimization in a way that still allows for the upper bounds.
    # design_matrix_bound = canon_features_stype_expanded.drop(columns=['batch', 'upper_bound'])
    # upper_bounds = canon_features_stype_expanded["upper_bound"]
    # print(design_matrix_bound.head(n=15))
    # print(upper_bounds.head(n=15))


def fit_original():
    tpms = get_tpm_data(linear=True, mean_centered=False)
    # Get the list of transcripts to use.
    transcripts_use = list(set(tpms.index.get_level_values(0).values))
    # Get design matrix
    design_df = make_design_matrix(transcripts_use=transcripts_use,
                                       upper_bound=False, aggregate=True)
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
    mod = smf.mixedlm(formula=formula, data=design_df,
                      groups=design_df['transcript'])
    # Fit the model
    mdf = mod.fit()
    # Assign the parameters to an output dataframe.
    params = pd.DataFrame(mdf.params)
    return(params)

def fit_with_bounds(params):
    tpms_mc = get_tpm_data(linear=False, mean_centered=True)
    print(tpms_mc.head())
    # Get the list of transcripts to use.
    transcripts_use = list(tpms_mc.index.values)
    # Get design matrix
    design_df = make_design_matrix(transcripts_use=transcripts_use,
                                   upper_bound=True, aggregate=False)
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

    # pars_init = pd.DataFrame([-1]*len(train_feats))
    # pars_init.index = train_feats

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
        # fig = plt.figure()
        # ax = fig.add_subplot(1, 1, 1)
        # ax.spines['left'].set_position('center')
        # ax.spines['bottom'].set_position('zero')
        # ax.spines['right'].set_color('none')
        # ax.spines['top'].set_color('none')
        # ax.xaxis.set_ticks_position('bottom')
        # ax.yaxis.set_ticks_position('left')

        # # plot the function
        # x = preds_mc
        # y = tpms_mc
        # plt.plot(x,y, 'o', color="black")

        # # show the plot
        # plt.show(block=False)
        obj = preds_mc.sub(tpms_mc).pow(2).sum().sum()
        global tick
        if tick%100 == 0:
            print(pars)
            print(obj)
        tick += 1

        return(obj)
    params = params.iloc[0:design_df.shape[1], :]
    res = minimize(objective_function, params, method="BFGS", options={"maxiter" : 1})
    # res = minimize(objective_function, params, options={"maxiter" : 1})
    params = res.x
    params = pd.DataFrame(params)
    params.index = train_feats
    return(params)
    # design_df = design_df.reindex(tpms.index).fillna(0.0)
    # print(design_df.head(n=30))
    # # Get features
    # # Remove the labels (do I need to do this? [I do, to get transcript index.])
    # tpms.reset_index(inplace=True)
    # design_df.reset_index(inplace=True)
    # # Define the list of features to use in the training.
    # train_feats = [i for i in all_features if np.std(design_df[i].values) != 0]
    # print("Training the model.")
    # formula = 'label ~ {} - 1'.format(' + '.join(train_feats))
    # # Add the data ('label') into the design matrix.
    # design_df["label"] = tpms["label"]
    # # Define the model.
    # mod = smf.mixedlm(formula=formula, data=design_df,
    #                   groups=design_df['transcript'])
    # # Fit the model
    # mdf = mod.fit()
    # # Assign the parameters to an output dataframe.
    # params = pd.DataFrame(mdf.params)
    # return(params)


def rescale_features(feature_df):
    feature_df.reset_index(level=[0,1], inplace=True)

    print(feature_df.head())
    sites = ["6mer"]
    for site in sites:
        inds_use = feature_df[feature_df[f'Stype_{site}']==1].index
        print(inds_use)
        cols_use = [col_i for col_i in feature_df.columns.values
                    if "_%s" %site in col_i and col_i != "Stype_%s" %site]
        print(cols_use)
        cols_use
        for col_use_i in cols_use:
            vals = feature_df.loc[inds_use, col_use_i]
            print(vals)
            break
    return()

def main():
    # Define the miRNAs.
    time_start = time.time()
    ## TODO. Make this mirs16 list an option based on the cell type; i.e., 


    pars_in = pd.read_csv("Repression/brand_new_pars.txt", sep="\t", index_col=0, header=None, names=[0])
    del pars_in.index.name


    params = fit_original()
    print(params)
    return
    time_done = time.time()
    print(("Finished first optimization in %3.2f seconds." %(time_done - time_start)))


    params_new = fit_with_bounds(params)
    print(params_new)
    params_new = pd.DataFrame(params_new)

    time_done = time.time()
    print(("Finished second optimization in %3.2f seconds." %(time_done - time_start)))


    params_for_output = []
    row = ['Intercept']
    print(row)
    for stype in ts7_stypes:
        stype = stype.replace('-','')
        row.append(params_new.loc[f'Stype_{stype}'][0])
        print(row)
    row += [0] * 8
    print(row)
    params_for_output.append(row)
    print(params_for_output)
    for feature in feature_list:
        row = [feature]
        for stype in ts7_stypes:
            stype = stype.replace('-','')
            try:
                row.append(params_new.loc[f'{feature}_{stype}'][0])
            except:
                print(feature, stype)
                row.append(0)
        row += [0] * 4 + [1] * 4
        params_for_output.append(row)
    cols = ['Feature'] + [f'{stype} coeff' for stype in ts7_stypes] + [f'{stype} min' for stype in ts7_stypes] + [f'{stype} max' for stype in ts7_stypes]
    params_for_output = pd.DataFrame(params_for_output, columns=cols).set_index('Feature')
    params_for_output.to_csv('Repression/bounded_parameters.txt', sep='\t')
    params_for_output.head()


    # params_for_output.to_csv("Repression/new_parameters_mean_centered.txt", sep="\t")



    return()
    # # Makes the prediction file that is as large as the data file.
    # canon_features_all_ix = canon_features_expanded_agg.reindex(tpm_long.index).fillna(0.0)
    # canon_features_all_ix = pd.concat([tpm_long, canon_features_all_ix], axis=1)
    # canon_features_all_ix['label'] = canon_features_all_ix['label'].astype(np.float64)
    # print("canon_features_all_ix.head()")
    # print(canon_features_all_ix.head())



    # return()


    # Part that seems to do the actual fitting.    

    design_matrix = canon_features_all_ix[train_feats]
    tpm_dim = tpm_long_df_norm.shape
    pars_init = pd.DataFrame([1]*73)
    pars_init.index = design_matrix.columns
    tick = 0
    # def loss_fun(pars):
    def objective_function(pars):
        preds = design_matrix.dot(pars).reset_index(level=[0,1])
        pred_df = preds.pivot(index="transcript", columns="mir").fillna(0)
        pred_df_norm = pred_df.sub(pred_df.mean(axis=1), axis=0)
        output = pred_df_norm.sub(tpm_long_df_norm).pow(2).sum().sum()
        global tick
        if tick%100 == 0:
            print(pars)
            print(output)
        tick += 1

        return(output)

    def bounded_objective_function(pars):
        preds_1 = pd.concat([design_matrix_bound.dot(pars), upper_bounds],
                            axis=1)
        return(preds_1)

    print(design_matrix_bound.head(15))
    print(upper_bounds.head(15))

    preds = bounded_objective_function(pars_init)    
    print(preds.head(15))





    return()




    # Makes the prediction file that is as large as the data file.
    canon_features_all_ix = canon_features_expanded_agg.reindex(tpm_long.index).fillna(0.0)



    print(objective_function(pars_init))

    # Alternative way of fitting the data, in which the data are mean-centered.
    # res = minimize(objective_function, pars_init, method="BFGS")
    # pars = res.x
    # pars_out = {i: j for (i, j) in zip(design_matrix.columns.values, pars)}
    # pars_out = pd.DataFrame.from_dict(pars_out, orient="index")
    # print(pars_out)
    # pars_out.to_csv("Repression/brand_new_pars.txt", sep="\t", index_label=False, header=False)



    formula2 = 'label ~ {}'.format(' + '.join(train_feats))

    mod = smf.mixedlm(formula=formula, data=train, groups=train['transcript'])
    print(mod)
    mdf = mod.fit()
    print("Finished training the model.")

    params = pd.DataFrame(mdf.params)
    print(params)

    mod2 = smf.mixedlm(formula=formula2, data=train, groups=train['transcript'])
    print(mod2)
    mdf2 = mod2.fit()
    params2 = pd.DataFrame(mdf2.params)

    print(params)
    print(params2)
    return()


    params_alt = pars_in
    print(params_alt)
    params = params_alt
    params_for_output = []
    row = ['Intercept']
    print(row)
    for stype in ts7_stypes:
        stype = stype.replace('-','')
        row.append(params.loc[f'Stype_{stype}'][0])
        print(row)
    row += [0] * 8
    print(row)
    params_for_output.append(row)
    print(params_for_output)
    for feature in feature_list:
        row = [feature]
        for stype in ts7_stypes:
            stype = stype.replace('-','')
            try:
                row.append(params.loc[f'{feature}_{stype}'][0])
            except:
                print(feature, stype)
                row.append(0)
        row += [0] * 4 + [1] * 4
        params_for_output.append(row)
    cols = ['Feature'] + [f'{stype} coeff' for stype in ts7_stypes] + [f'{stype} min' for stype in ts7_stypes] + [f'{stype} max' for stype in ts7_stypes]
    params_for_output = pd.DataFrame(params_for_output, columns=cols).set_index('Feature')
    params_for_output.to_csv('Repression/new_parameters.txt', sep='\t')
    params_for_output.head()


    params_for_output = []
    row = ['Intercept']
    print(row)
    for stype in ts7_stypes:
        stype = stype.replace('-','')
        row.append(params.loc[f'Stype_{stype}'][0])
        print(row)
    row += [0] * 8
    print(row)
    params_for_output.append(row)
    print(params_for_output)
    for feature in feature_list:
        row = [feature]
        for stype in ts7_stypes:
            stype = stype.replace('-','')
            try:
                row.append(params.loc[f'{feature}_{stype}'][0])
            except:
                print(feature, stype)
                row.append(0)
        row += [0] * 4 + [1] * 4
        params_for_output.append(row)
    cols = ['Feature'] + [f'{stype} coeff' for stype in ts7_stypes] + [f'{stype} min' for stype in ts7_stypes] + [f'{stype} max' for stype in ts7_stypes]
    params_for_output = pd.DataFrame(params_for_output, columns=cols).set_index('Feature')
    params_for_output.to_csv('Repression/new_parameters.txt', sep='\t')
    params_for_output.head()


    params_for_output.to_csv("Repression/new_parameters_mean_centered.txt", sep="\t")



    time_done = time.time()
    print(("Finished in %3.2f seconds." %(time_done - time_start)))
    return()


if __name__ == "__main__":
    main()
