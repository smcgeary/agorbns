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

# FUNCTIONS
# def LCSuff(seq1, seq2):
#     if len(seq1) == 0 or len(seq2) == 0:
#         return("")
#     elif seq1[-1] == seq2[-1]:
#         return(LCSuff(seq1[:-1], seq2[:-1]) + seq1[-1])
#     else:
#         return("")

# def LCSubString(seq1, seq2):
#     pres1 = [seq1[:i] for i in range(1, len(seq1) + 1)]
#     pres2 = [seq2[:i] for i in range(1, len(seq2) + 1)]
#     pre_pairs = [(i, j) for i in pres1 for j in pres2]
#     LCs = [LCSuff(i[0], i[1]) for i in pre_pairs]
#     LCmax = max([len(i) for i in LCs])
#     return [i for i in LCs if len(i) == LCmax]

def assign_site(file_in, _mirna, _sitelist):
    print(SITES_3p)
    threeP_only_fc_dict = {i : [] for i in SITES_3p}
    mirna_str = _mirna.name
    mirna_str = re.sub("miR-", "mir", mirna_str)
    mirna_str = re.sub("lsy-6", "lsy6", mirna_str)
    mirna_str = re.sub("let-7a", "let7", mirna_str)
    mirna_str = re.sub("-.*nt$", "", mirna_str)
    print(mirna_str)
    lines = file_in.read()
    lines_new = []

    utrs = csv.DictReader(lines.splitlines(), delimiter="\t")
    totalrows = 0
    sites = [site.name for site in _sitelist.sites]
    flanks_dict = {site : [] for site in sites}
    site_cols = {site : i + 1 for i, site in enumerate(sites)}
    for utr in utrs:
        if utr["nosite3"] != "":
            totalrows += 1
    print(totalrows)
    utrs = csv.DictReader(lines.splitlines(), delimiter="\t")
    # output_list = [None]*totalrows
    rep_df = [[None] + [0 for i in sites] + [0] for j in range(totalrows)]
    print(rep_df[0])
    position_matrix = [[None] + [[] for i in sites] for j in range(totalrows)]
    rep_3p_auto_utrs = []
    rep_3p_supp_utrs = []
    i_r = 0
    check = 0
    utr_AU_list = []
    for utr in utrs:
        check += 1
        if check % 100 == 0:
            print(check)
        if utr["nosite3"] != "":
        # Make UTR object.
            _utr = Utr(utr["sequence"])
            utr_AU = [nuc in ["A", "T"] for nuc in utr["sequence"]]
            utr_AU = sum(utr_AU)/float(len(utr_AU))
            utr_AU_list += [utr_AU]
            _utr.get_all_sites_in_utr(_sitelist)
            add_sites_from_utr(_sitelist, _utr)
            site_names = [_site.name for _site in _utr.sites]
            rep = float(utr[mirna_str]) - float(utr["nosite3"])
            rep_df[i_r][0] = utr[""]
            if site_names != ["None"]:
                ls = [_site.l for _site in _utr.sites]
                rs = [_site.r for _site in _utr.sites]
                ends = zip(ls, rs)
                flanks = [_utr.seq[end[0]-2:end[0]] + _utr.seq[end[1]:end[1]+2]
                          for end in ends]
                AU_cs = [sum([i in ["A", "T"] for i in flank])/4.0
                         for flank in flanks]
                for site_name in site_names:
                    # full_sites[site_name].append(float(rep))
                    rep_df[i_r][site_cols[site_name]] += 1
                    position_matrix[i_r][site_cols[site_name]] += ends
                for site_name, flank in zip(site_names, flanks):
                    flanks_dict[site_name] += AU_cs
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

    rep_3p_auto_df = pd.DataFrame(rep_3p_auto_utrs)
    if rep_3p_auto_df.shape[0] > 0:
        rep_3p_auto_df.set_index([0], inplace=True)
        rep_3p_auto_df.index.name = None
        rep_3p_auto_df.columns = ["Seed", "ThreeP", "s_l", "s_r", "thr_l", "thr_r", "dist", "fc"]

    rep_3p_supp_df = pd.DataFrame(rep_3p_supp_utrs)
    if rep_3p_supp_df.shape[0] > 0:
        rep_3p_supp_df.set_index([0], inplace=True)
        rep_3p_supp_df.index.name = None
        rep_3p_supp_df.columns = ["Seed", "ThreeP", "s_l", "s_r", "thr_l", "thr_r", "dist", "fc"]
    print("here")
    AU_utr_mean = np.mean(utr_AU_list)
    AU_utr_sem = stats.sem(utr_AU_list)
    return(rep_df, rep_3p_auto_df, rep_3p_supp_df, i_r, flanks_dict, AU_utr_mean, AU_utr_sem)

def main():
    time_start = time.time()
    arguments = ["miRNA", "sitelist", "-test_binary", "-new_binary"]
    mirna, sitelist, test, new = parse_arguments(arguments)
    _mirna = Mirna(mirna)

    _sitelist = SiteListRep(_mirna, sitelist)
    sites = [site.name for site in _sitelist.sites]
    args  = [_mirna, _sitelist]
    ############################################################################
    # For each file in the input list, perform the multiprocessing:
    if new:
        rep_file = "/lab/bartel4_ata/kathyl/RNA_Seq/outputs/biochem/final/merged.txt"
    else:
        rep_file = "/lab/bartel4_ata/kathyl/RNA_Seq/analysis/data/no_baseline_analysis/merged.txt"
    

    print("about to open file")
    with open(rep_file) as file_in:
        rep_df, rep_3p_auto_df, rep_3p_supp_df, i_r, flanks_dict, AU_utr_mean, AU_utr_sem = assign_site(file_in, _mirna,
                                                            _sitelist)


    for site in _sitelist.sites:
        site_str = site.name
        print(site_str)
        print(np.mean(flanks_dict[site_str]))
        print(len(flanks_dict[site_str]))
    print(AU_utr_mean)
    print(AU_utr_sem)



    rep_df_each_site = pd.DataFrame.from_dict({site : [np.mean(rep_df.loc[rep_df[site] >= 1, "fc"]/math.log(2)),
                                                       np.std(rep_df.loc[rep_df[site] >= 1, "fc"]/math.log(2)),
                                                       rep_df.loc[rep_df[site] >= 1].shape[0]]
                        for site in sites}, orient="index").reindex(sites, copy=False)

    flank_df_each_site = pd.DataFrame.from_dict({site : [np.mean(flanks_dict[site]),
                                                         np.std(flanks_dict[site]),
                                                         len(flanks_dict[site])]
                        for site in sites}, orient="index").reindex(sites, copy=False)
    print(flank_df_each_site)

    rep_df_no_site = pd.DataFrame.from_dict({"None" : [np.mean(rep_df.loc[rep_df.iloc[:,:-1].sum(axis=1) == 0, "fc"]/math.log(2)),
                                                      np.std(rep_df.loc[rep_df.iloc[:,:-1].sum(axis=1) == 0, "fc"]/math.log(2)),
                                                      rep_df.loc[rep_df.iloc[:,:-1].sum(axis=1) == 0].shape[0]]}, orient="index")
    rep_df_each_site.columns = ["mean", "sd", "n"]
    fc_nosite_vec = rep_df.loc[rep_df.iloc[:,:-1].sum(axis=1) == 0, "fc"]
    extension = sitelist
    if new:
        extension = "%s_new" %(extension)
    else:
        extension = "%s_old" %(extension)
    rep_df_path = get_analysis_path(mirna, "repression_hela_cs", extension, "lin_model_df")
    if not test:
        rep_df.to_csv(rep_df_path, sep="\t")
    else:
        print(rep_df.shape)
        print(rep_df.iloc[0:6, ])
        print(rep_df.sum(axis=0))
    print(rep_df_path)
    print_time_elapsed(time_start)
################################################################################

if __name__ == "__main__":
    main()

