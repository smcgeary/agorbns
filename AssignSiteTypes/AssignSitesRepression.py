################################################################################
#MakeSiteTypeReadFiles.py
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

from scipy import stats


SITES_3p = ["8mer-m9.16",
            "8mer-m10.17", "9mer-m9.17",
            "8mer-m11.18", "9mer-m10.18", "10mer-m9.18",
            "8mer-m12.19", "9mer-m11.19", "10mer-m10.19", "11mer-m9.19",
            "8mer-m13.20", "9mer-m12.20", "10mer-m11.20", "11mer-m10.20",
            "8mer-m14.21", "9mer-m13.21", "10mer-m12.21", "11mer-m11.21",
            "8mer-m15.22", "9mer-m14.22", "10mer-m13.22", "11mer-m12.22",
            "8mer-m16.23", "9mer-m15.23", "10mer-m14.23", "11mer-m13.23"]

SITES_SEED = ["8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-A1", "6mer-m8",
              "8mer-mmT5", "5mer-m3.7", "8mer-mmA3", "5mer-m8"]


def update_mean(mean_1, n_1, mean_2, n_2):
    weight_1 = mean_1*n_1
    weight_2 = mean_2*n_2
    mean_overall = (weight_1 + weight_2)/float(n_1 + n_2)
    return(mean_overall)

def GetAU_freq(sequence):
    return(np.mean([nuc in ["A", "T"] for nuc in sequence]))




def analyze_repression_file(file_in, _mirna, _sitelist, bg_method, exclude_ribosome, new):
    # Convert Kathy-style miRNA naming to my style naming.
    if new:
        nosite_list = ["no_6mer", "no_guide_canon", "no_guide_canon_pass_6mer"]
    else:
        nosite_list = ["nosite1", "nosite2", "nosite3"]
    mirna_str = _mirna.name
    mirna_str = re.sub("miR-", "mir", mirna_str)
    mirna_str = re.sub("lsy-6", "lsy6", mirna_str)
    mirna_str = re.sub("let-7a", "let7", mirna_str)
    mirna_str = re.sub("-.*nt$", "", mirna_str)
    # Take in entire file, check for the total number of lines that have a
    # defined repression for "nosite3".
    lines = file_in.read()
    utrs = csv.DictReader(lines.splitlines(), delimiter="\t")
    totalrows = 0
    for utr in utrs:
        if utr[nosite_list[bg_method]] != "":
            totalrows += 1
    # Construct all the variables and objects that will be updated: the
    # repression data frame, that has gene_name as first column,
    # the site indeces as the second column, and the the fold-change value as
    # the third column, the positional matrix, for where each miRNA site is
    # within each 3'UTR, and the 3p autonomous and 3p supplementary UTRs.
    # Also define i_r, which is the index of each row when iterating over the
    # 3' UTRs, check is the variable that defines which rows get printed.
    # utr_AU_list is a list of all the individual AU contents values per each 
    # 3' UTR, for which the mean and sem are calculated at the end of the
    # script.
    sites = [site.name for site in _sitelist.sites]
    rep_df = [[None] + [0 for i in sites] + [0] for j in range(totalrows)]

    # rep_site_df = [[None] + [0 for i in sites] + [0] for j in range(totalrows)]

    flank_df = [[[] for i in range(256)] for j in range(totalrows)]
    flank_num_df = [[0 for i in range(256)] for j in range(totalrows)]

    flank_site_df = {site.name : {flank : 0 for flank in get_kmer_list(4)}
                     for site in _sitelist.sites}
    position_matrix = [[None] + [[] for i in sites] for j in range(totalrows)]
    rep_3p_auto_utrs = []
    rep_3p_supp_utrs = []
    i_r = 0
    check = 0
    AUfreq_utr_list = []
    AUfreq_totalutr_mean = 0
    AUfreq_totalutr_length = 0
    # Define the list of sites, and then make the dictionary that has the
    # list of AU frequency within the flanking dinucleotides surrounding each
    # site instance. 
    AUfreq_site_map = {site : [] for site in sites}
    # TODO: how does this work?
    site_cols = {site : i + 1 for i, site in enumerate(sites)}
    flank_cols = {kmer : i for i, kmer in enumerate(get_kmer_list(4))}


    for flank in get_kmer_list(4):
        print(flank_cols[flank])
    # Iterate again over the utrs, checking the sites actually.

    utrs = csv.DictReader(lines.splitlines(), delimiter="\t")
    # Iterate over the utrs:
    for utr in utrs:
        check += 1
        if check % 100 == 0:
            print(check)
        if utr[nosite_list[bg_method]] != "":
            # Make UTR object.
            _utr = Utr(utr["sequence"])
            # Get the AU frequency of each flanking dinucleotide
            len_utr = len(utr["sequence"])
            AUfreq_utr = GetAU_freq(utr["sequence"])
            AUfreq_utr_list += [AUfreq_utr]
            AUfreq_totalutr_mean = update_mean(AUfreq_totalutr_mean,
                                               AUfreq_totalutr_length,
                                               AUfreq_utr, len_utr)
            AUfreq_totalutr_length += len_utr
            # Add all the sites to the UTR.
            _utr.get_all_sites_in_utr(_sitelist)
            add_sites_from_utr(_sitelist, _utr)
            site_names = [_site.name for _site in _utr.sites]
            # Assign the repression
            rep = float(utr[mirna_str]) - float(utr[nosite_list[bg_method]])
            # Assign the UTR name.
            rep_df[i_r][0] = utr[""]
            if site_names != ["None"]:
                # AU frequency in the flanking nucleotides determination.
                ends = [(_site.l, _site.r) for _site in _utr.sites]
                flanks = [_utr.seq[end[0]-2:end[0]] + _utr.seq[end[1]:end[1]+2]
                          for end in ends]
                AUfreq_sites = [GetAU_freq(flank) for flank in flanks]
                for site_name, flank, end, AUfreq_site in zip(site_names, flanks, ends, AUfreq_sites):
                    if len(flank) == 2:
                        if end[1] == len(_utr.seq):
                            flank += "AA"
                        elif end[0] == 0:
                            if utr[""] in ["NM_000700", "NM_001004",
                                           "NM_001137550", "NM_002160"]:
                                flank = "AA" + flank
                            elif utr[""] in ["NM_004515"]:
                                flank = "GA"
                    elif len(flank) == 3:
                        if end[1] == len(_utr.seq) - 1:
                            flank += "A"
                    if (end[0] >= 15 or not exclude_ribosome) and len(flank) == 4:
                        rep_df[i_r][site_cols[site_name]] += 1
                        position_matrix[i_r][site_cols[site_name]] += end
                        flank_df[i_r][flank_cols[flank]] += [site_name]
                        flank_num_df[i_r][flank_cols[flank]] += 1
                        AUfreq_site_map[site_name] += [AUfreq_site]
                        flank_site_df[site_name][flank] += 1

                if set(site_names).intersection(set(SITES_3p)):
                    ends_3p = [[site] + list(ends[i])
                               for i, site in enumerate(site_names)
                               if site in SITES_3p]
                    ends_non_3p = [[site] + list(ends[i])
                                   for i, site in enumerate(site_names)
                                   if site not in SITES_3p]
                    paired_3p_non3p_sites = []
                    for i in ends_3p:
                        site_3p, l_3p, r_3p = i
                        for j in ends_non_3p:
                            site_seed, l_seed, r_seed = j
                            dist = l_seed - r_3p - 1
                            paired_3p_non3p_sites.append([utr[""], site_seed,
                                                          site_3p, l_seed,
                                                          r_seed, l_3p, r_3p,
                                                          dist, rep])
                    if (paired_3p_non3p_sites):
                        min_dist = float("inf")
                        min_tup = paired_3p_non3p_sites[0]
                        for tup in paired_3p_non3p_sites:
                            dist = tup[-2]
                            if ((min_dist > 0 and dist < min_dist) or
                                (min_dist < 0 and dist > 0)):
                                min_tup = tup
                                min_dist = dist
                        if min_dist >= 20 or min_dist < 0:
                            rep_3p_auto_utrs.append(min_tup)
                        else:
                            rep_3p_supp_utrs.append(min_tup)
                    else:
                        rep_3p_auto_utrs.append([utr[""], None, None,
                                                None, None,
                                                None, None,
                                                None, rep])
            rep_df[i_r][len(sites) + 1] = rep
            i_r += 1
    
    rep_df = pd.DataFrame(rep_df).set_index([0])
    rep_df.index.name = None
    rep_df.columns = sites + ["fc"]

    flank_df = pd.DataFrame(flank_df)
    flank_df.index = rep_df.index
    flank_df.index.name = None
    flank_df.columns = [kmer[:2] + "." + kmer[2:] for kmer in get_kmer_list(4)]

    flank_num_df = pd.DataFrame(flank_num_df)
    flank_num_df.index = rep_df.index
    flank_num_df.index.name = None
    flank_num_df.columns = flank_df.columns


    rep_3p_auto_df = pd.DataFrame(rep_3p_auto_utrs)

    if rep_3p_auto_df.shape[0] > 0:
        rep_3p_auto_df.set_index([0], inplace=True)
        rep_3p_auto_df.index.name = None
        rep_3p_auto_df.columns = ["Seed", "ThreeP", "s_l", "s_r", "thr_l",
                                  "thr_r", "dist", "fc"]

    rep_3p_supp_df = pd.DataFrame(rep_3p_supp_utrs)
    if rep_3p_supp_df.shape[0] > 0:
        rep_3p_supp_df.set_index([0], inplace=True)
        rep_3p_supp_df.index.name = None
        rep_3p_supp_df.columns = ["Seed", "ThreeP", "s_l", "s_r", "thr_l",
                                  "thr_r", "dist", "fc"]
    AUfreq_eachutr_mean = np.mean(AUfreq_utr_list)
    AUfreq_eachutr_std = np.std(AUfreq_utr_list)

    return(rep_df, flank_df, flank_num_df, rep_3p_auto_df, rep_3p_supp_df, i_r, AUfreq_site_map,
           AUfreq_eachutr_mean, AUfreq_eachutr_std, AUfreq_totalutr_mean, flank_site_df)

def main():
    # Start timer.
    time_start = time.time()
    # Parse arguments.
    arguments = ["miRNA", "sitelist", "bg_method", "-exrib_binary",
                 "-new_binary", "-test_binary"]
    mirna, sitelist, bg_method, exrib, new, test = parse_arguments(arguments)
    # Make miRNA and sitelist objects.
    _mirna = Mirna(mirna)
    # Get Kds for Kd conditional
    # if mirna == "miR-7-23nt":
    #     kd_exp = "equilibrium2_nb"
    # else:
    #     kd_exp = "equilibrium"
    # if mirna in ["miR-1", "miR-7-23nt"]:
    #     combined_str = "_nocombInput"
    #     if mirna == "miR-1":
    #         buffer3p_str = "_buffer3p"
    #     else:
    #         buffer3p_str = ""
    # else:
    #     combined_str = ""
    # ext = "%s%s_PAPER" %(buffer3p_str, combined_str)
    # kds_path = get_analysis_path(mirna, kd_exp, "5_%s" %(sitelist), "kds_PAPER", ext=ext)
    # kds = pd.read_csv(kds_path, sep="\t")



    _sitelist = SiteListRep(_mirna, sitelist)

    args  = [_mirna, _sitelist, new]
    ############################################################################
    # For each file in the input list, perform the multiprocessing:
    if new:
        rep_file = "/lab/bartel4_ata/kathyl/RNA_Seq/outputs/biochem/final/merged.txt"
    else:
        rep_file = "/lab/bartel4_ata/kathyl/RNA_Seq/analysis/data/no_baseline_analysis/merged.txt"


    print("about to open file")
    with open(rep_file) as file_in:
        (
            rep_df, flank_df, flank_num_df, rep_3p_auto_df, rep_3p_supp_df, i_r, AUfreq_site_map,
            AUfreq_eachutr_mean, AUfreq_eachutr_std, AUfreq_totalutr_mean, flank_site_df
        ) = analyze_repression_file(file_in, _mirna, _sitelist,
                                    int(bg_method) - 1, exrib, new)

    sites = [site.name for site in _sitelist.sites]

    # Make the repression matrix:
    rep_df_each_site = pd.DataFrame.from_dict(
        {site : [np.mean(rep_df.loc[rep_df[site] >= 1, "fc"]/math.log(2)),
                 np.std(rep_df.loc[rep_df[site] >= 1, "fc"]/math.log(2)),
                 rep_df.loc[rep_df[site] >= 1].shape[0]]
         for site in sites}, orient="index").reindex(sites, copy=False)
    # make the AU frequency matrix:
    AUfreq_site_df = pd.DataFrame.from_dict(
        {site : [np.mean(AUfreq_site_map[site]),
                 np.std(AUfreq_site_map[site]),
                 len(AUfreq_site_map[site])]
         for site in sites}, orient="index").reindex(sites, copy=False)
    
    AUfreq_site_df.columns = ["mean", "std", "n"]
    AUfreq_utrs_df = pd.DataFrame([[AUfreq_eachutr_mean, AUfreq_eachutr_std,
                                    rep_df.shape[0]],
                                   [AUfreq_totalutr_mean, None, None]])
    AUfreq_utrs_df.columns = ["mean", "std", "n"]
    AUfreq_utrs_df.index = ["each utr", "total utr"]
    # Merge the two dataframes into one dataframe.
    AUfreq_df = AUfreq_site_df.append(AUfreq_utrs_df)

    rep_df_no_site = pd.DataFrame.from_dict({"None" : [np.mean(rep_df.loc[rep_df.iloc[:,:-1].sum(axis=1) == 0, "fc"]/math.log(2)),
                                                      np.std(rep_df.loc[rep_df.iloc[:,:-1].sum(axis=1) == 0, "fc"]/math.log(2)),
                                                      rep_df.loc[rep_df.iloc[:,:-1].sum(axis=1) == 0].shape[0]]}, orient="index")
    
    site_flank_df = pd.DataFrame.from_dict(flank_site_df).transpose()
    rep_sites_all = pd.concat([rep_df.iloc[:, :-1], flank_num_df, rep_df.iloc[:, -1]], axis=1, join='inner')


    print("no site average:")
    print(rep_df_no_site)

    rep_df_each_site.columns = ["mean", "sd", "n"]
    fc_nosite_vec = rep_df.loc[rep_df.iloc[:,:-1].sum(axis=1) == 0, "fc"]
    extension = sitelist + "_" + bg_method

    if exrib:
        extension += "_exrib"

    if new:
        extension = "%s_new" %(extension)
    else:
        extension = "%s_old" %(extension)
    rep_df_path = get_analysis_path(mirna, "repression_hela_cs", extension, "lin_model_df")
    rep_flank_df_path = get_analysis_path(mirna, "repression_hela_cs", extension, "lin_flank_model_df")
    AUfreq_df_path = get_analysis_path(mirna, "repression_hela_cs", extension, "AUfreq_df")
    site_flank_count_path = get_analysis_path(mirna, "repression_hela_cs", extension, "site_flank_count_df")
    if not test:
        rep_df.to_csv(rep_df_path, sep="\t")
        rep_sites_all.to_csv(rep_flank_df_path, sep="\t")
        site_flank_df.to_csv(site_flank_count_path, sep="\t")
        print(rep_flank_df_path)
        print(site_flank_count_path)
        print("wrote flank dataframe")
    else:
        print(rep_df.shape)
        print(rep_df.iloc[0:6, ])
        print(rep_df.sum(axis=0))
    AUfreq_df.to_csv(AUfreq_df_path, sep="\t")
    print(AUfreq_df_path)
    print(rep_df_path)
    print_time_elapsed(time_start)
################################################################################

if __name__ == "__main__":
    main()

