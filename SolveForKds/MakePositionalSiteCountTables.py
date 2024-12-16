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
        if (condition in ["I", "I_combined", "0", "0_combined"]
                and mirna in ["miR-7-23nt", "miR-7-24nt", "miR-7-25nt"]):
            mirna_temp = "miR-7-23nt"
        else:
            mirna_temp = mirna
        path = get_analysis_path(mirna_temp, experiment, condition,
                                 "site_counts", ext=extension)
        print(path)
        try:
            counts = pd.read_csv(path, header=None, index_col=0, sep="\t")
            print(counts)
            sXc[condition] = counts
        except ValueError:
            print("error")
    return(sXc)

def main():
    time_start = time.time()
    # Define all the relevant arguments.
    arguments = ["miRNA", "experiment", "n_constant", "sitelist", "-mir_start",
                 "-split16", "-buffer3p_binary"]
    mirna, experiment, n_constant, sitelist, mir_start, split16, buffer3p = (
        parse_arguments(arguments)
    )
    # Initialize objects:
    
    _mirna = Mirna(mirna)

    if sitelist == "12mers":
        _sitelist = TwelvemerSiteList(_mirna, mir_start)
        extension = "_%s_%s_%s-%s" %(n_constant, sitelist, mir_start,
                                  int(mir_start) + 3)
    elif sitelist == "16mers":
        _sitelist = SixteenmerSiteList(_mirna, mir_start, split16)
        extension = "_%s_%s_%s-%s_%s" %(n_constant, sitelist, mir_start,
                                  int(mir_start) + 3, split16)
    else:
        _sitelist = SiteList(_mirna, sitelist)
        extension = "_%s_%s" %(n_constant, sitelist)
    if buffer3p:
        extension = "%s_buffer3p" %(extension)
    print(_sitelist.sites)
    # Make site table:
    sXc = make_data_table(_sitelist, experiment, extension, buffer3p)
    sXc.columns = CONDITIONS_R[experiment][mirna]
    if _sitelist.__class__.__name__ == "SiteList":
        print(sXc)
    else:
        print(sXc.iloc[:10, :])
    path = get_analysis_path(mirna, experiment, "all_sites%s" %(extension),
                             "site_count_tables")
    print(path)
    sXc.to_csv(path, sep="\t") 
    ensure_directory("%s%s/%s/kds" % (SOLEXA_DIR, mirna, experiment)) 
    ensure_directory("%sfigures/kds/%s/%s" % (HOME_DIR, mirna, experiment)) 
    print_time_elapsed(time_start)

################################################################################

if __name__ == "__main__":
    main()
