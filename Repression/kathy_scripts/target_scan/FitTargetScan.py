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

import statsmodels.api as sm
import statsmodels.formula.api as smf

from scipy.optimize import minimize
from scipy.stats import linregress

# CONSTANTS OF THE SCRIPT
tick = 0

sites = ["6mer", "7mer-1a", "7mer-m8", "8mer-1a"]

mirs17 = sorted([  'mir1',  'mir124', 'mir155',    'mir7',
                   'lsy6',  'mir153', 'mir139',  'mir144',
                 'mir223',  'mir137', 'mir205',  'mir143',
                 'mir182', 'mir199a', 'mir204', 'mir216b',
                                                   'let7'])

mirs16 = sorted([  'mir1',  'mir124', 'mir155',    'mir7',
                   'lsy6',  'mir153', 'mir139',  'mir144',
                 'mir223',  'mir137', 'mir205',  'mir143',
                 'mir182', 'mir199a', 'mir204', 'mir216b'])

mirs6 = sorted(['mir1', 'mir124', 'mir155', 'mir7',
                                    'lsy6', 'let7'])

mirs5 = sorted(['mir1', 'mir124', 'mir155', 'mir7',
                                            'lsy6'])

list_mirs_map = {
    "five"      :  mirs5,
    "six"       :  mirs6,
    "sixteen"   : mirs16,
    "seventeen" : mirs17
}

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

def get_tpm_data(linear=True, mirs=mirs16, mean_centered=False):
    TPM_FILE = '/lab/bartel4_ata/kathyl/RNA_Seq/outputs/biochem/final_fixed/merged.txt'
    # Read in the tpm file.
    transXmir_tpm = pd.read_csv(TPM_FILE, sep='\t', index_col=0)
    # Names the index names as transcript.
    transXmir_tpm.index.name = 'transcript'
    transXmir_tpm = transXmir_tpm.loc[:, transXmir_tpm.columns.isin(mirs)]
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

def rescale_features(feat_df, scale_df, center_feat=False):
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
            print(col_use_i)
            # Grab the values.
            vals = feat_df.loc[inds_use, col_use_i]
            # Determine the 5th and 95th percentiles.
            coef_0 = np.percentile(vals, 0)
            coef_min = np.percentile(vals, 5)
            coef_max = np.percentile(vals, 95)
            coef_100 = np.percentile(vals, 100)
            print([coef_0, coef_min, coef_max, coef_100])
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



def make_design_and_scale_dfs(transcripts_use=[], mirs=mirs16, aggregate=True,
                              upper_bound=False, rescale_feat=False,
                              center_feat=False, two_bmodes=False, oligoG=False,
                              **kwargs):
    if len(kwargs) != 0:
        print("kwargs")
        print(kwargs)
    else:
        print("no kwargs")
    path_special = ",".join(["%s_%s" %(key, kwargs[key]) for key in kwargs.keys()])
    print(path_special)
    if two_bmodes:
        if path_special == "":
            path_special = "two_bmodes"
        else:
            path_special = ",two_bmodes"

    if oligoG:
        if path_special == "":
            path_special = "oligoG"
        else:
            path_special = ",oligoG"



    if path_special != "":
        path = "/lab/solexa_bartel/mcgeary/transfections/HeLa/target_scan/new_threep_scores/%s/features.txt" %path_special
    else:
        path = ('/lab/bartel4_ata/kathyl/RNA_Seq/data/TS7_output_files/'
                'refseq_utr90_fixed_ribosome_shadow/features.txt')
    print(path)
    feat_df = pd.read_csv(path,  sep='\t')


    # Rename "miRNA_family" column as "mir".
    feat_df.rename(columns={'miRNA family': 'mir'}, inplace=True)
    # Swap " " for "_" in column names.
    feat_df.columns = [x.replace(' ','_') for x in feat_df.columns]
    # Rename "Gene_ID" column as "transcript".
    feat_df.rename(columns={'Gene_ID': 'transcript'}, inplace=True)
    # Selects on those rows that have miRNAs within the sixteen HeLa miRNAs.
    # NOTE, this will have to be changed in order to allow testing on the test
    # set.
    feat_df = feat_df[feat_df['mir'].isin(mirs)]
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
    if rescale_feat or center_feat:
        feat_df, scale_df = rescale_features(feat_df, scale_df,
                                             center_feat=center_feat)
    # Replace dashes with no spaces.
    feat_df.columns = [x.replace('-','') for x in feat_df.columns]
    feat_df = feat_df.astype(np.float64)
    if aggregate:
        feat_df = feat_df.groupby(['transcript','mir']).agg(np.sum)
    return(feat_df, scale_df)

def calculate_bound_predictions(pars, design_df, tpms_mc):
    # Separate the upper bounds from the design matrix.
    upper_bounds = design_df["upper_bound"]
    design_df.drop(columns="upper_bound", inplace=True)
    # Calculate the predictions using the design matrix
    pars = pars.loc[design_df.columns.values, :]
    # preds = pd.concat([design_df.dot(pars), upper_bounds], axis=1).min(axis=1)
    preds = design_df.dot(pars)
    # Sum over the sites within each transcript.
    preds = preds.groupby(['transcript','mir']).agg(np.sum).reset_index(level=[0,1])
    # Reset the index and convert the matrix into a mirna-by-transcript
    # matrix.
    preds = preds.pivot(index="transcript", columns="mir", values=0).fillna(0)
    preds = preds.reindex(tpms_mc.index).fillna(0.0)
    return(preds)

def fit_original(params = [], mirs=mirs16, rescale_feat=False, two_bmodes=False, **kwargs):
    print(two_bmodes)
    print(kwargs)
    tpms = get_tpm_data(mirs=mirs, linear=True, mean_centered=False)
    # Get the list of transcripts to use.
    transcripts_use = list(set(tpms.index.get_level_values(0).values))
    # Get design matrix
    design_df, scale_df = make_design_and_scale_dfs(mirs=mirs, transcripts_use=transcripts_use,
                                   upper_bound=False, aggregate=True,
                                   rescale_feat=rescale_feat,
                                   two_bmodes=two_bmodes, **kwargs)
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

def fit_alt(params=[], mirs=mirs16, rescale_feat=False, center_feat=False, two_bmodes=False, **kwargs):
    tpms_mc = get_tpm_data(linear=False, mean_centered=True)
    # Get the list of transcripts to use.
    transcripts_use = list(tpms_mc.index.values)
    # Get design matrix
    design_df, scale_df = make_design_and_scale_dfs(mirs=mirs16,
        transcripts_use=transcripts_use, upper_bound=False, aggregate=True,
        two_bmodes=two_bmodes, center_feat=center_feat, **kwargs
    )
    # Get features
    all_features = design_df.columns
    train_feats = [i for i in all_features if np.std(design_df[i].values) != 0]
    print(train_feats)
    # Assign the final design matrix.
    design_df = design_df[train_feats]
    print(design_df.head())
    # Sum over the sites within each transcript.
    # Initialize the parameters.
    if len(params) == 0:
        params = pd.DataFrame([0]*len(train_feats), index=train_feats)
    def objective_function(pars):
        # Get the predicted effect of each site, using the upper bounds.
        preds = design_df.dot(pars).reset_index(level=[0,1])
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


def fit_with_bounds(params=[], mirs=mirs16, rescale_feat=False, **kwargs):
    tpms_mc = get_tpm_data(linear=False, mean_centered=True)
    # Get the list of transcripts to use.
    transcripts_use = list(tpms_mc.index.values)
    # Get design matrix
    design_df, scale_df = make_design_and_scale_dfs(
        mirs=mirs, transcripts_use=transcripts_use, upper_bound=True,
        aggregate=False
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
        params = pd.DataFrame([-1]*len(train_feats), index=train_feats)
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

    # # Calculate the final predictions:
    # preds = pd.concat([design_df.dot(params), upper_bounds], axis=1).min(axis=1)
    # # Sum over the sites within each transcript.
    # preds = preds.groupby(['transcript','mir']).agg(np.sum).reset_index(level=[0,1])
    # # Reset the index and convert the matrix into a mirna-by-transcript
    # # matrix.
    # preds = preds.pivot(index="transcript", columns="mir", values=0).fillna(0)
    # preds = preds.reindex(tpms_mc.index).fillna(0.0)
    # preds_mc = mean_center(preds)


    return(params, scale_df)


def main():
    # Define the miRNAs.
    time_start = time.time()
    # Parse command line arguments:
    arguments = ["-mirnas", "-bounded_binary", "-rescaled_binary", "-centered_binary","-alt_binary", "-inpars_binary",
                 "-upstream_limit", "-min_pairing", "-supp_l", "-supp_r", "-offset_opt",
                 "-offset_tol", "-w_pairing", "-w_supp", "-w_offset", 
                 "-use_kd_model_binary", "-two_bmodes_binary", "-oligoG_binary"]
    args = parse_arguments(arguments)
    (mirnas, bounded, rescaled, centered, alt, inpars, upstream_limit,
     min_pairing, supp_l, supp_r, offset_opt,
     offset_tol, w_pairing, w_supp, w_offset,
     use_kds, two_bmodes, oligoG) = args
    if bounded:
        bounded_str = "_bounded"
        alt_str = ""
    elif alt:
        bounded_str = ""
        alt_str ="_alt"
    else:
        bounded_str = ""
        alt_str = ""

    if centered:
        centered_str = "_centered"
    else:
        centered_str = ""

    if rescaled:
        rescaled_str = "_rescaled"
    else:
        rescaled_str = ""

    if not mirnas:
        mirnas = "sixteen"

    mirs = list_mirs_map[mirnas]
    mirs_str = "_%s" %mirnas

    cell_line="HeLa"

    print(mirs)



    path = "/lab/solexa_bartel/mcgeary/transfections/%s" %cell_line
    path += "/target_scan/new_threep_scores/"
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
    if two_bmodes:
        if path_special == "":
            path_special = "two_bmodes"
        else:
            path_special += ",two_bmodes"

    if oligoG:
        if path_special == "":
            path_special = "oligoG"
        else:
            path_special += ",oligoG"

    if path_special != "":
        path_special = "_" + path_special


    print(path_special)

    base_out_dir = "/lab/solexa_bartel/mcgeary/transfections/HeLa/target_scan/"

    path_condition = "%s%s%s%s%s" %(mirs_str, bounded_str, alt_str, rescaled_str, centered_str)

    # Assign the path of the eventual output dataframe.
    # params_dir = params/"
    # params_file = "params%s%s%s%s%s%s.txt" %(path_special, mirs_str, bound_str, alt_str, rescaled_str, centered_str)
    # params_path = params_dir + params_file
    params_path = base_out_dir + "params/params" + path_special + path_condition + ".txt"
    path_preds = base_out_dir + "predictions/predictions" + path_special + path_condition + ".txt"
    path_preds_mc = base_out_dir + "predictions/predictions_mc" + path_special + path_condition + ".txt"
    path_preds_sep = base_out_dir + "predictions_separated/predictions_separated" + path_special + path_condition + ".txt"
    path_tpms = base_out_dir + "tpms/tpms" + path_special + path_condition + ".txt"
    path_tpms_mc = base_out_dir + "tpms/tpms_mc" + path_special + path_condition + ".txt"
    path_r2 = base_out_dir + "r2/r2" + path_special + path_condition + ".txt"
 
    print(params_path)
 
    if inpars:
        params_in = pd.read_csv(params_path, sep="\t", index_col=0)
        params = pd.DataFrame(
            [0]*4*params_in.shape[0],
            index=["%s_%s" %(i, j.replace("-", "")) for i in feature_list for j in sites] + ["Stype_%s" %i.replace("-", "") for i in sites])

        for site in sites:
            # First, assign the intercept coefficients:
            input_col = "%s coeff" %site
            params_ind = "Stype_%s" %site.replace("-", "")
            # output_df.loc["Intercept", output_col] = params.loc[params_ind, 0]
            params.loc[params_ind, 0] = params_in.loc["Intercept", input_col] 
            for feat in feature_list:
                params_ind = "%s_%s" %(feat, site.replace("-", ""))
                # if params_ind in params.index.values:
                params.loc[params_ind, 0] = params_in.loc[feat, input_col]

    else:
        # The base fit, used to initialize the parameters.
        if alt:
            params, scale_df = fit_alt(mirs=mirs, rescale_feat=rescaled, center_feat=centered,
                                       two_bmodes=two_bmodes, oligoG=oligoG, **kwargs)
        else:
            params, scale_df = fit_original(mirs=mirs, rescale_feat=rescaled, center_feat=centered,
                                            two_bmodes=two_bmodes, oligoG=oligoG, **kwargs)
        # If the upper bounded flag is chosen, this refits with the new objective
        # function that utilizes the upper bounds.
        if bounded:
            params, scale_df = fit_with_bounds(params, mirs=mirs, rescale_feat=rescaled, center_feat=centered,
                                               two_bmodes=two_bmodes, oligoG=oligoG, **kwargs)
    ## CALCULATE THE BOUNDED PREDICTIONS AND PRINT THE R^2 FOR THE PARAMETERS ##
    # Get the mean-centered transfection tpm data.
    tpms = get_tpm_data(linear=False, mirs=mirs, mean_centered=False)
    tpms_16 = get_tpm_data(linear=False, mirs=mirs16, mean_centered=False)
    tpms_17 = get_tpm_data(linear=False, mirs=mirs17, mean_centered=False)
    tpms_mc = mean_center(tpms)
    tpms_16_mc = mean_center(tpms_16)
    tpms_17_mc = mean_center(tpms_17)

    design_df, scale_df = make_design_and_scale_dfs(
        mirs=mirs, transcripts_use=tpms_mc.index.values, aggregate=False,
        upper_bound=True, rescale_feat=rescaled, center_feat=centered,
        two_bmodes=two_bmodes, oligoG=oligoG, **kwargs
    )
    design_16_df, scale_16_df = make_design_and_scale_dfs(
        mirs=mirs16, transcripts_use=tpms_mc.index.values, aggregate=False,
        upper_bound=True, rescale_feat=rescaled, center_feat=centered,
        two_bmodes=two_bmodes, oligoG=oligoG, **kwargs
    )
    design_17_df, scale_17_df = make_design_and_scale_dfs(
        mirs=mirs17, transcripts_use=tpms_mc.index.values, aggregate=False,
        upper_bound=True, rescale_feat=rescaled, center_feat=centered,
        two_bmodes=two_bmodes, oligoG=oligoG, **kwargs
    )
    print(design_17_df.head())
    # Generate the design matrix and the scaling factors for the parameters.
    sites_short = [i.replace("-", "") for i in sites]
    dot_separate = []
    for feat_i in feature_list + ["Stype"]:
        print(feat_i)
        labels = ["%s_%s" %(feat_i, site_i) for site_i in sites_short]
        labels = [i for i in labels if i in design_df.columns.values]
        if len(labels) != 0:
            design_feat = design_df[labels]
            params_feat = params.loc[labels]
            feat_dot = design_feat.dot(params_feat)
            feat_dot.columns = [feat_i]
            dot_separate.append(feat_dot)
    dot_full = pd.concat(dot_separate, axis=1)
    print(dot_full.head())


    # dir_preds_sep = "/lab/solexa_bartel/mcgeary/transfections/HeLa/target_scan/predictions_separated/"
    # file_preds_sep = "predictions_separated%s%s%s%s%s.txt" %(path_special, bound_str, alt_str, rescaled_str, centered_str)
    # path_preds_sep = dir_preds_sep + file_preds_sep

    dot_full.to_csv(path_preds_sep, sep='\t')



    preds = calculate_bound_predictions(params, design_df, tpms_mc)
 
    preds_16 = calculate_bound_predictions(params, design_16_df, tpms_16_mc)
    preds_17 = calculate_bound_predictions(params, design_17_df, tpms_17_mc)

    print(preds.head(n=2))
    preds_mc = mean_center(preds)
    preds_16_mc = mean_center(preds_16)
    preds_17_mc = mean_center(preds_17)
    # dir_preds = "/lab/solexa_bartel/mcgeary/transfections/HeLa/target_scan/predictions/"
    # file_preds = "predictions%s%s%s%s%s.txt" %(path_special, bound_str, alt_str, rescaled_str, centered_str)
    # file_preds_mc = "predictions_mc%s%s%s%s%s.txt" %(path_special, bound_str, alt_str, rescaled_str, centered_str)
    # path_preds = dir_preds + file_preds
    # path_preds_mc = dir_preds + file_preds_mc
    # print(path_preds)
    # print(path_preds_mc)

    # dir_tpms = "/lab/solexa_bartel/mcgeary/transfections/HeLa/target_scan/tpms/"
    # file_tpms = "tpms%s%s%s%s%s.txt" %(path_special, bound_str, alt_str, rescaled_str, centered_str)
    # file_tpms_mc = "tpms_mc%s%s%s%s%s.txt" %(path_special, bound_str, alt_str, rescaled_str, centered_str)
    # path_tpms = dir_tpms + file_tpms
    # path_tpms_mc = dir_tpms + file_tpms_mc
    # print(path_tpms)
    # print(path_tpms_mc)
    preds.to_csv(path_preds, sep='\t')
    preds_mc.to_csv(path_preds_mc, sep='\t')
    tpms.to_csv(path_tpms, sep='\t')
    tpms_mc.to_csv(path_tpms_mc, sep='\t')

    print(preds_17_mc.columns.values)
    print(tpms_17_mc.columns.values)

    # Calculate the r^2 with the sixteen miRNAs.
    r2_list = []
    for mir_list in [mirs5, mirs6, mirs16, mirs17]:
        if "let7" in mir_list:
            preds_mc_use = preds_17_mc
            tpms_mc_use = tpms_17_mc
        else:
            preds_mc_use = preds_16_mc
            tpms_mc_use = tpms_16_mc
        print("_________________")
        print(preds_mc_use.columns.values)
        print(tpms_mc_use.columns.values)
        x_val = pd.melt(preds_mc_use, value_vars=mir_list)["value"]
        y_val = pd.melt(tpms_mc_use, value_vars=mir_list)["value"]
        r2 = (linregress(x_val, y_val)[2])**2
        print(mir_list)
        print(r2)
        r2_list.append(r2)

    pd.DataFrame([r2_list], index=["r2"], columns=["five", "six", "sixteen", "seventeen"]).to_csv(path_r2, sep="\t")
    print(path_r2)


    ## CONSTRUCT THE OUTPUT MATRIX FOR THE PARAMETERS ##########################
    if not inpars:
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

        # Output the parameters.
        output_df.to_csv(params_path, sep='\t')




if __name__ == "__main__":
    main()
