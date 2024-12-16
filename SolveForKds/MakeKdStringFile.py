################################################################################
#GenerateSiteTypeKds.py
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

################################################################################
# FUNCTIONS
def make_data_table(_sitelist, experiment, extension):
    sites = _sitelist["names"] + ["None"]
    mirna = _sitelist.mirna
    conditions = CONDITIONS_PY[experiment][mirna]
    sXc = pd.DataFrame(0, index=sites, columns=conditions)
    for condition in conditions:
        path = get_analysis_path(mirna, experiment, condition, "site_counts",
                                 ext=extension)
        sXc[condition] = pd.read_csv(path, header=None, index_col=0, sep="\t")
    return(sXc)

def main():
    time_start = time.time()
    # Define all the relevant arguments.
    arguments = ["miRNA", "experiment", "n_constant", "sitelist", "-mir_start"]
    mirna, experiment, n_constant, sitelist, mir_start = parse_arguments(arguments)
    # Initialize objects:
    _mirna = Mirna(mirna)
    if sitelist == "12mers":
        _sitelist = TwelvemerSiteList(_mirna, mir_start)
        condition = "_%s_%s_%s-%s" %(n_constant, sitelist, mir_start,
                                  int(mir_start) + 3)
    else:
        _sitelist = SiteList(_mirna, sitelist)
    # Make site table:
    kds_full = []
    seqs_full = []
    for site, seq in zip(_sitelist["names"], _sitelist["seqs"]):
        print(site)
        print(seq)
        condition = "%s_%s_%s_PAPER" %(n_constant, sitelist, site)
        path = get_analysis_path(mirna, experiment, condition,
                                 "kds_PAPER")
        pars = pd.read_csv(path, sep="\t")
        flanks_l = [re.sub(r'^([^.]*).([^.]*)$', r'\g<1>', i) for
                    i in pars.index.values]
        flanks_r = [re.sub(r'^([^.]*).([^.]*)$', r'\g<2>', i) for 
                    i in pars.index.values]
        sequence_full = ["".join([flank_l, seq, flank_r]) for
                         flank_l, flank_r in zip(flanks_l, flanks_r)]

        seqs_full += sequence_full
        kds = list(logit(pars["Mean"], 10))
        kds_full += kds
        if site == "8mer-bT(4.6)":
            print(pars.iloc[159,:])
            print(kds[159])
    out_df = pd.DataFrame(kds_full)
    out_df.index = seqs_full
    path = get_analysis_path(mirna, experiment, "site_and_flank_kds",
                             "kds_PAPER",
                             ext="_%s_%s" %(n_constant, _sitelist.name))
    print(path)
    out_df.transpose().to_csv(path, sep="\t", index=False)
    print_time_elapsed(time_start)

################################################################################

if __name__ == "__main__":
    main()
