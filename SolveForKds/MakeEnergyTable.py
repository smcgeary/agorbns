################################################################################
#MakeEnergyTable.py
################################################################################
import imp # Used to import general.py
imp.load_source("general", "general/general.py")
imp.load_source("RBNS_methods", "general/RBNS_methods.py")
from general import *
from RBNS_methods import *

################################################################################
# FUNCTIONS
def make_data_table(mirna):

    return()

def main():
    time_start = time.time()
    # Define all the relevant arguments.
    # arguments = [
    #     "miRNA", "experiment", "n_constant", "sitelist", "-mir_start",
    #     "-start_mm", "-stop_mm", "-split16", "-uniq_binary", "-buffer3p_binary",
    #     "-comp_binary", "-c_k", "-m_k", "-new_binary", "-new2_binary",
    #     "-zeros_binary"
    # ]
    # (
    #     mirna, experiment, n_constant, sitelist, mir_start, start_mm, stop_mm,
    #     split16, uniq, buffer3p, comp, c_k, m_k, new, new2, zeros
    # ) = parse_arguments(arguments)

    mirnas = ["miR-1", "let-7a", "miR-155", "miR-124", "lsy-6", "miR-7"]
    site_types = ["8mer", "7mer-m8", "7mer-A1", "6mer"]

    dGdf = pd.DataFrame(index=site_types, columns=mirnas, dtype=float)
    for mirna in mirnas:
        print("_____________________")
        _mirna = Mirna(mirna)
        seed = _mirna.seq[:8]

        # Get the base sites:
        seq_8mer = _mirna["8mer"]
        seq_7merm8 = _mirna["7mer-m8"]
        seq_7merA1 = _mirna["7mer-A1"]
        seq_6mer = _mirna["6mer"]
        # For each site type, calculate all possible sites, then remove larger
        # site-types within.
        all_seq_7merm8 = ["%s%s" %(seq_7merm8, i) for i in DNTS]
        all_seq_7merm8 = list(set(all_seq_7merm8) - set([seq_8mer]))
        all_seq_7merA1 = ["%s%s" %(i, seq_7merA1) for i in DNTS]
        all_seq_7merA1 = list(set(all_seq_7merA1) - set([seq_8mer]))
        all_seq_6mer = ["%s%s%s" %(i[0], seq_6mer, i[1]) for i in it.product(DNTS, DNTS)]
        all_seq_6mer = list(
            set(all_seq_6mer) - set([seq_8mer]) - set(all_seq_7merm8)
            - set(all_seq_7merA1)
        )
        # Calculate the duplex energy for the seed with each site, and average
        # them.
        dG_8mer = np.array([get_mfe(seed, i) for i in [seq_8mer]]).mean()
        dG_7merm8 = np.array([get_mfe(seed, i) for i in all_seq_7merm8]).mean()
        dG_7merA1 = np.array([get_mfe(seed, i) for i in all_seq_7merA1]).mean()
        dG_6mer = np.array([get_mfe(seed, i) for i in all_seq_6mer]).mean()
        # Allocate each to its position in the table.
        dGdf.loc["8mer", mirna]= dG_8mer
        dGdf.loc["7mer-m8", mirna] = dG_7merm8
        dGdf.loc["7mer-A1", mirna] = dG_7merA1
        dGdf.loc["6mer", mirna] = dG_6mer
    # Write the table.
    path = "general/miRNA_canonical_site_energies.txt"
    print(dGdf)
    dGdf.to_csv(path, sep="\t") 
    print_time_elapsed(time_start)

################################################################################

if __name__ == "__main__":
    main()
