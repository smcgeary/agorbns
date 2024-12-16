################################################################################
#LibraryDesign.py
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
from operator import add

pd.set_option('max_columns', 600)
pd.set_option("max_colwidth", 600)
# FUNCTIONS

MIRNAS = ["miR-1", "let-7a", "miR-155", "miR-124", "lsy-6", "miR-7"]


#FINAL DESIGN CONSTRAINTS
kLEN = 180



def main():
    # Fill dataframe with miRNA, name, and sequence of each site.
    suffix = "oct12_final"
    out_path = "variants/%s/final_order.fa" %(suffix)
    print(out_path)
    with open(out_path, "w+") as file_out:
        for mirna in MIRNAS:
            print(mirna)
            _mirna = Mirna(mirna)
            _sitelist = SiteList(_mirna, "paperfinal", 37)
            for _site in _sitelist.sites:
                name = _site.name
                input_path = "variants/%s/%s/%s_library_variants_%s.txt" %(
                    suffix, _mirna.name, name, suffix)
                with open(input_path) as file:
                    seqs = [i.split("\t")[3] for i in list(file)][1:]
                num = 1
                for seq in seqs:
                    fa_header = "> %s_%s_%s\n" %(mirna, name, num)
                    file_out.write(fa_header)
                    file_out.write(seq + "\n")
                    num += 1



################################################################################

if __name__ == "__main__":
    main()

