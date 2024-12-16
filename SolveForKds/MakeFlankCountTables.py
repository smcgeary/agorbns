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
def make_data_table(mirna, experiment, extension, site):
    flanks = SiteList.FLANKS
    conditions = CONDITIONS_PY[experiment][mirna]
    fXc = pd.DataFrame(0, index=flanks, columns=conditions)
    for condition in conditions:
        path = get_analysis_path(mirna, experiment, condition, "flanks",
                                 ext=extension)
        fXc[condition] = pd.read_csv(path, index_col=0, sep="\t")[site]
    return(fXc)

def main():
    time_start = time.time()
    # Define all the relevant arguments.
    arguments = ["miRNA", "experiment", "n_constant", "sitelist"]
    mirna, experiment, n_constant, sitelist = parse_arguments(arguments)
    extension = "_%s_%s" %(n_constant, sitelist, 1)
    # Initialize objects:
    _mirna = Mirna(mirna)
    _sitelist = SiteList(_mirna, sitelist)
    for site in _sitelist["names"]:
        print(site)
        # Get site-flanking count table:
        sfXc = make_data_table(mirna, experiment, extension, site)
        sfXc.columns = CONDITIONS_R[experiment][mirna]
        sfXc = sfXc[(sfXc.T != 0).any()]
        path = get_analysis_path(mirna, experiment, "%s_flanking" %(site),
                                 "site_count_tables", ext=extension)
        print(path)
        # sfXc.to_csv(path, sep="\t")
    print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()
