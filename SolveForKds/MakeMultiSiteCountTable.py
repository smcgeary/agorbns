################################################################################
#MakeMultiSiteCountTable.py
################################################################################
import imp # Used to import general.py
imp.load_source("general", "general/general.py")
imp.load_source("RBNS_methods", "general/RBNS_methods.py")
from general import *
from RBNS_methods import *

################################################################################
# FUNCTIONS
def make_data_table(_sitelist, experiment, extension):
    sites = _sitelist["names"] + ["None"]
    mirna = _sitelist.mirna
    conditions = CONDITIONS_PY[experiment][mirna]
    sXc = pd.DataFrame(0, index=sites, columns=conditions)
    sXc_list = [None]*len(conditions)
    for i, condition in enumerate(conditions):
        path = get_analysis_path(mirna, experiment, condition, "multisite_counts",
                                 ext=extension)
        print(path)
        print(pd.read_csv(path, header=None, index_col=0, sep="\t"))
        sXc_list[i] = pd.read_csv(path, header=None, index_col=0,
                                  na_filter=False, sep="\t")
    sXc = pd.concat(sXc_list, axis=1).fillna(0)
    print(sXc)
    return(sXc)

def main():
    time_start = time.time()
    # Define all the relevant arguments.
    arguments = ["miRNA", "experiment", "n_constant", "sitelist", "-uniq_binary",
                 "-buffer3p_binary"]
    mirna, experiment, n_constant, sitelist, uniq, buffer3p = parse_arguments(arguments)
    # Initialize objects:
    _mirna = Mirna(mirna)
    _sitelist = SiteList(_mirna, sitelist, 1)
    extension = "_%s_%s" %(n_constant, sitelist)
    if uniq:
        extension = "%s_uniq" %(extension)
    if buffer3p:
        extension = "%s_buffer3p" %(extension)
    # Make site table:
    sXc = make_data_table(_sitelist, experiment, extension)
    sXc.columns = CONDITIONS_R[experiment][mirna]
    path = get_analysis_path(mirna, experiment, "all_sites%s" %(extension),
                             "multisite_count_tables")
    sXc.to_csv(path, sep="\t")
    print(sXc)
    print(path)
    ensure_directory("data/processed/%s/%s/kds" % (mirna, experiment)) 
    ensure_directory("figures/kds/%s/%s" % (mirna, experiment)) 
    print_time_elapsed(time_start)

################################################################################

if __name__ == "__main__":
    main()
