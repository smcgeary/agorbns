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

CONDITIONS1_NB = ["I", "12.6", "12.6_2", "4", "1.26", "0.4", "0"]
CONDITIONS2_NB = ["I", "A12.65", "A12.65_2", "A4", "A1.265", "A0.4", "A0"]
CONDITIONS1 = ["I", "I_combined", "40", "12.6", "4", "1.26", "0.4", "0"]
CONDITIONS2 = ["I", "I_combined", "A40", "A12.65", "A4", "A1.265", "A0.4", "A0"]

################################################################################
# FUNCTIONS
def make_data_table(_sitelist, site, experiment, extension):
    flanks = SiteList.FLANKS
    mirna = _sitelist.mirna
    conditions = CONDITIONS_PY[experiment][mirna]
    matrix = pd.DataFrame(0, index=flanks, columns=conditions)
    for cond in conditions:

        counts_flank_sites_map = {flank: 0 for flank in flanks}
        flank_path = get_analysis_path(mirna, experiment, cond, "flanks",
                                       extension)

        alt = pd.DataFrame.from_csv(flank_path, sep="\t")
        if cond == "I_combined" and site == "8mer":
            print(site)
            print(cond)
            print(flank_path)
            print(alt)
            print(alt.sum(axis=0))

        matrix[cond] = alt[site]
    return(matrix)

def main():
    time_start = time.time()
    # Define all the relevant arguments.
    arguments = ["miRNA", "experiment", "n_constant", "sitelist",
    "-buffer3p_binary"]
    (mirna, experiment,
     n_constant, sitelist, buffer3p) = parse_arguments(arguments)
    # Initialize objects:
    _mirna = Mirna(mirna)
    _sitelist = SiteList(_mirna, sitelist, 1)
    if experiment == "equil_seed_nb":
        conditions1 = CONDITIONS1_NB
        conditions2 = CONDITIONS2_NB
    else:
        conditions1 = CONDITIONS1
        conditions2 = CONDITIONS2
    extension = "_%s_%s" %(n_constant, sitelist)
    if buffer3p:
        extension = "%s_buffer3p" %(extension)
    for site in _sitelist["names"]:
        # print(site)
        # Get site-flanking count table:
        sfXc = make_data_table(_sitelist, site, experiment, extension)
        # print(sfXc)
        sfXc.columns = CONDITIONS_R[experiment][mirna]
        sfXc = sfXc[(sfXc.T != 0).any()]
        path = get_analysis_path(mirna, experiment, "%s_flanking" %(site),
                                 "site_count_tables", ext=extension)
        sfXc.to_csv(path, sep="\t")
    print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()
