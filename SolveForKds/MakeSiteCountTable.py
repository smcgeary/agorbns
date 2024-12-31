################################################################################
#MakeSiteCountTable.py
################################################################################
import imp # Used to import general.py
imp.load_source("general", "general/general.py")
imp.load_source("RBNS_methods", "general/RBNS_methods.py")
from general import *
from RBNS_methods import *

################################################################################
# FUNCTIONS
def make_data_table(_sitelist, experiment, extension,
                    buffer3p=False, comp=False, include_zeros=False):
    count_dir = "site_counts"
    mirna = _sitelist.mirna
    if _sitelist.name == "programmedbase":
        if experiment == "equilibrium" or experiment == "equilibrium2_nb":
            sites_path = get_analysis_path(mirna, experiment, "I_combined",
                                     count_dir, ext=extension)
        else:
            sites_path = get_analysis_path(mirna, experiment, "I",
                                     count_dir, ext=extension)
        counts = pd.read_csv(sites_path, header=None, index_col=0, 
                             filter_na=False, sep="\t")
        sites = list(counts.index)
    else:
        sites = _sitelist["names"] + ["None"]
    if experiment == "kinetics":
        sites = ["%s_p" %(i) for i in sites] + ["%s_c" %(i) for i in sites]
    conditions = CONDITIONS_PY[experiment][mirna]
    if comp:
        conditions = [i for i in conditions if i != "I_combined"]

    if include_zeros:
        sXc = pd.DataFrame(0., index=[], columns=[])
    else:
        sXc = pd.DataFrame(0., index=sites, columns=conditions)
    for condition in conditions:
        if (condition in ["I", "I_combined", "0", "0_combined"]
                and mirna in ["miR-7-23nt", "miR-7-24nt", "miR-7-25nt"]
                and experiment == "equilibrium2_nb"):
            mirna_temp = "miR-7-23nt"
        else:
            mirna_temp = mirna

        path = get_analysis_path(mirna_temp, experiment, condition,
                                 count_dir, ext=extension)
        try:
            print(path)
            counts = pd.read_csv(path, header=None, na_filter=False,
                                 index_col=0, sep="\t")
            print(counts)
            if include_zeros:
                sXc = pd.concat([sXc, counts], axis=1)
            else:
                sXc[condition] = counts.iloc[:, 0]

        except ValueError:
            print("error")
    sXc.fillna(value=0, inplace=True)
    return(sXc)

def main():
    time_start = time.time()
    # Define all the relevant arguments.
    arguments = [
        "miRNA", "experiment", "n_constant", "sitelist", "-mir_start",
        "-start_mm", "-stop_mm", "-split16", "-uniq_binary", "-buffer3p_binary",
        "-comp_binary", "-c_k", "-m_k", "-new_binary", "-new2_binary",
        "-zeros_binary"
    ]
    (
        mirna, experiment, n_constant, sitelist, mir_start, start_mm, stop_mm,
        split16, uniq, buffer3p, comp, c_k, m_k, new, new2, zeros
    ) = parse_arguments(arguments)
    ### MOD The `new2` argument is there for the modification for the let-7a
    # reporter library experiments, in which the kds associated with the AAAA
    # UUUU sites need to be corrected.
    # Initialize objects:
    _mirna = Mirna(mirna)
    output_dir = "site_count_tables"

    if sitelist == "12mers":
        _sitelist = TwelvemerSiteList(_mirna, mir_start)
        extension = "_%s_%s_%s-%s" %(n_constant, sitelist, mir_start,
                                  int(mir_start) + 3)
    elif sitelist == "16mers":
        _sitelist = SixteenmerSiteList(_mirna, mir_start, split16)
        extension = "_%s_%s_%s-%s_%s" %(n_constant, sitelist, mir_start,
                                  int(mir_start) + 3, split16)
    elif "randthrp" in sitelist:
        _sitelist = SiteList(_mirna, "programmedbase", 1)
        extension = "_%s_%s" %(n_constant, sitelist)
    elif "progthrp" in sitelist:
        _sitelist = SiteList(_mirna, "programmedbase", 1)
        extension = "_%s_%s" %(n_constant, sitelist)
    else:
        _sitelist = SiteList(_mirna, sitelist, 1)
        extension = "_%s_%s" %(n_constant, sitelist)

    if start_mm and stop_mm and ("randthrp" in sitelist or "progthrp" in sitelist):
        extension = "%s_m%s.%smmsd" %(extension, start_mm, stop_mm)
    if uniq:
        extension = "%s_uniq" %(extension)

    if comp:
        extension = "%s_comp" %(extension)
        if c_k:
            c_k = int(c_k)
        else:
            c_k = 6
        if m_k:
            m_k = int(m_k)
        else:
            m_k = 6

        m_k = min(m_k, c_k)
        extension = "%s_c%sm%s" %(extension, c_k, m_k)
    if buffer3p:
        extension = "%s_buffer3p" %(extension)
    if new:
        extension = "%s_new" %(extension)
    elif new2:
        extension = "%s_new2" %(extension)
    # Make site table:
    sXc = make_data_table(_sitelist, experiment, extension, buffer3p=buffer3p,
                          comp=comp, include_zeros=zeros)
    colnames = CONDITIONS_R[experiment][mirna]
    if comp:
        colnames = [i for i in colnames if i != "I_combined"]
    sXc.columns = colnames
    if _sitelist.__class__.__name__ == "SiteList":
        print(sXc)
        print(sXc.sum(axis=0))
    else:
        print(sXc.iloc[:10, :])
    if zeros:
        extension = "%s_zeros" %(extension)

    path = get_analysis_path(mirna, experiment, "all_sites%s" %(extension),
                             output_dir)
    # print(sXc.iloc[:, 2:].sum(axis=1))
    # sXc = sXc[sXc.iloc[:, 2:].sum(axis=1) != 0]
    print(path)
    sXc.to_csv(path, sep="\t") 
    ensure_directory("data/processed/%s/%s/kds" % (mirna, experiment)) 
    ensure_directory("figures/kds/%s/%s" % (mirna, experiment)) 
    print_time_elapsed(time_start)

################################################################################

if __name__ == "__main__":
    main()
