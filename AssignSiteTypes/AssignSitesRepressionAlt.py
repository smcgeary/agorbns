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

# FUNCTIONS
def LCSuff(seq1, seq2):
    if len(seq1) == 0 or len(seq2) == 0:
        return("")
    elif seq1[-1] == seq2[-1]:
        return(LCSuff(seq1[:-1], seq2[:-1]) + seq1[-1])
    else:
        return("")

def LCSubString(seq1, seq2):
    pres1 = [seq1[:i] for i in range(1, len(seq1) + 1)]
    pres2 = [seq2[:i] for i in range(1, len(seq2) + 1)]
    pre_pairs = [(i, j) for i in pres1 for j in pres2]
    LCs = [LCSuff(i[0], i[1]) for i in pre_pairs]
    LCmax = max([len(i) for i in LCs])
    return [i for i in LCs if len(i) == LCmax]

def assign_site(file_in, _mirna, _sitelist):
    mirna_str = _mirna.name
    mirna_str = re.sub("miR-", "mir", mirna_str)
    mirna_str = re.sub("lsy-6", "lsy6", mirna_str)
    mirna_str = re.sub("let-7a", "let7", mirna_str)
    mirna_str = re.sub("-.*nt$", "", mirna_str)
    print(mirna_str)
    time_start = time.time()
    lines = file_in.read()
    lines_new = []
    utrs = csv.DictReader(lines.splitlines(), delimiter="\t")
    totalrows = 0
    sites = [site.name for site in _sitelist.sites]
    site_cols = {site : i for i, site in enumerate(sites)}
    site_pairs = [",".join(i) for i in it.product(sites, repeat=2)]
    single_sites = {site : [] for site in sites + ["None"]}
    full_sites = {site : [] for site in sites + ["None"]}
    double_sites = {pair : [] for pair in site_pairs + ["None"]}
    for utr in utrs:
        if utr["nosite3"] != "":
            totalrows += 1
    utrs = csv.DictReader(lines.splitlines(), delimiter="\t")
    output_list = [None]*totalrows
    print([i for i in sites])
    output_matrix = [[0 for i in sites] + [0] for j in range(totalrows)]
    print(output_matrix[0])
    print(output_matrix[1])
    i_r = 0
    sites_8mer = []
    for utr in utrs:
        if utr["nosite3"] != "":
        # Make UTR object.
            _utr = Utr(utr["sequence"])
            _utr.get_all_sites_in_utr(_sitelist)
            add_sites_from_utr(_sitelist, _utr)
            site_names = [_site.name for _site in _utr.sites]
            rep = float(utr[mirna_str]) - float(utr["nosite3"])
            if "8mer" in site_names:
                sites_8mer += [i_r]
            if site_names != ["None"]:
                ls = [_site.l for _site in _utr.sites]
                rs = [_site.r for _site in _utr.sites]
                dists = [l - r for l, r in zip(ls[1:], rs[:-1])] + [rs[-1]]
                dist_strings = ["(%s)" %(i) for i in dists]
                site_dist_pars = ["%s,%s" %(site, dist_strings[i]) for i, site in enumerate(site_names)]
                output_site = "(%s)," %(ls[0]) + ",".join(site_dist_pars)
                if len(site_names) == 1:
                    single_sites[site_names[0]].append(float(rep))
                if len(site_names) == 2:
                    site_pair_key = ",".join(site_names)
                    double_sites[site_pair_key].append(float(rep))
                for site_name in site_names:
                    full_sites[site_name].append(float(rep))
                    output_matrix[i_r][site_cols[site_name]] += 1
            else:
                output_site = "None,(" + str(len(_utr)) + ")" 
                single_sites["None"].append(float(rep))
                double_sites["None"].append(float(rep))
                full_sites["None"].append(float(rep))
            output_list[i_r] = (output_site, rep)
            output_matrix[i_r][len(sites)] = rep
            i_r += 1
    return output_list, single_sites, double_sites, full_sites, output_matrix

def main():
    time_start = time.time()
    arguments = ["miRNA", "sitelist", "-test_binary"]
    args = parse_arguments(arguments)
    mirna, sitelist, test = args
    _mirna = Mirna(mirna)

    _sitelist = SiteListRep(_mirna, sitelist)

    print(_sitelist)
    args  = [_mirna, _sitelist]
    ############################################################################
    # For each file in the input list, perform the multiprocessing:
    rep_file = "/lab/bartel4_ata/kathyl/RNA_Seq/analysis/data/no_baseline_analysis/merged.txt"
    with open(rep_file) as file_in:
        out = assign_site(file_in, _mirna, _sitelist)

    print(out)

    all_sites = out[0]
    dist = [0]*30
    for row in all_sites:
        sites = row[0].split(",")
        if sites[0] == "None":
            dist[0] +=1
        else:
            sites = sites[1:-1]
            num_sites = (len(sites) + 1)/2
            if num_sites >= len(dist):
                dist_new = [0]*(num_sites + 1)
                for i, d in enumerate(dist):
                    dist_new[i] = d
                dist = dist_new
            dist[num_sites] +=1
    print(dist)
    sites = [site.name for site in _sitelist.sites]
    total_list = out[0]
    single_sites, double_sites, full_sites = out[1:4]
    output_matrix = out[4]

    print(output_matrix[0])
    print(output_matrix[1])
    output_matrix_df = pd.DataFrame(output_matrix)
    output_matrix_df.columns = sites + ["fc"]

    rep_file = "/lab/bartel4_ata/kathyl/RNA_Seq/analysis/data/no_baseline_analysis/merged.txt"
    with open(rep_file, "r") as file_in:
        lines = file_in.read()
    utrs = csv.DictReader(lines.splitlines(), delimiter="\t")
    mRNA_names = []
    for utr in utrs:
        if utr["nosite3"] != "":
            mRNA_names.append(utr[""])

    output_matrix_df.index = mRNA_names


    single_sites_summary = {site: [np.mean(single_sites[site]),
                                np.std(single_sites[site]),
                                len(single_sites[site])]
                         for site in single_sites.keys()}
    double_sites_summary = {site: [np.mean(double_sites[site]),
                                np.std(double_sites[site]),
                                len(double_sites[site])]
                         for site in double_sites.keys()}

    double_sites_isolated = {site: []
                         for site in sites + ["None"]}
    # print(double_sites_isolated)
    for key in double_sites.keys():
        if key != "None":
            key1, key2 = key.split(",")
            # print([key1, key2])
            # print(double_sites[key])
            # print(double_sites_isolated[key1])
            # print(double_sites_isolated[key2])
            double_sites_isolated[key1] += double_sites[key]
            double_sites_isolated[key2] += double_sites[key]
        else:
            double_sites_isolated["None"] += double_sites[key]

    full_sites_summary = {site: [np.mean(full_sites[site]),
                                np.std(full_sites[site]),
                                len(full_sites[site])]
                         for site in full_sites.keys()}

    double_sites_isolated_summary = {site: [np.mean(double_sites_isolated[site]),
                                np.std(double_sites_isolated[site]),
                                len(double_sites_isolated[site])]
                         for site in double_sites_isolated.keys()}


    single_sites_df = pd.DataFrame.from_dict(single_sites_summary, orient="index").reindex(sites + ["None"], copy=False)
    double_sites_df = pd.DataFrame.from_dict(double_sites_summary, orient="index")
    double_sites_isolated_df = pd.DataFrame.from_dict(double_sites_isolated_summary, orient="index").reindex(sites + ["None"], copy=False)

    full_sites_df = pd.DataFrame.from_dict(full_sites_summary, orient="index").reindex(sites + ["None"], copy=False)
    single_sites_df.columns = ["mean", "sd", "n"]
    double_sites_df.columns = ["mean", "sd", "n"]
    double_sites_isolated_df.columns = ["mean", "sd", "n"]
    full_sites_df.columns = ["mean", "sd", "n"]
    
    extension = "_" + sitelist
    outputpaths = [get_analysis_path(mirna, "repression_hela_cs", "", i,
                                     ext=extension)
                   for i in ["singlesite_rep", "doublesite_rep",
                             "fullsite_rep", "doublesite_isolated_rep"]]
    print(outputpaths)
    print(output_matrix_df.iloc[1:10,1:10])
    print(output_matrix_df.sum(axis=0))
    output_matrix_path = get_analysis_path(mirna, "repression_hela_cs", sitelist + "_alt" , "lin_model_df")
    print(output_matrix_path)
    output_matrix_df.to_csv(output_matrix_path, sep="\t")

    single_sites_df.to_csv(outputpaths[0], sep="\t", header=["mean", "sd", "n"])
    double_sites_df.to_csv(outputpaths[1], sep="\t", header=["mean", "sd", "n"])
    full_sites_df.to_csv(outputpaths[2], sep="\t", header=["mean", "sd", "n"])
    double_sites_isolated_df.to_csv(outputpaths[3], sep="\t", header=["mean", "sd", "n"])

    print_time_elapsed(time_start)
################################################################################

if __name__ == "__main__":
    main()

