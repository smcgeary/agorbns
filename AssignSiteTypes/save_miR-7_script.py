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


def main():
	path = "/lab/solexa_bartel/mcgeary/AgoRBNS/miR-7-23nt/equilibrium2_nb/kmers_cutoff_final/*"
	ls_list = ls(path)
	for path_in in ls_list:
		# This gets the "_ex8" part of the file out, to check the conditionals.
		file_name = path_in.split(".")[0]
		path_out = "%s_pre_8mer-bA8.txt" %(file_name)
		print(path_in)
		print(path_out)
		with open(path_in , "r+") as file_in:
			with open(path_out, "w+") as file_out:
				file_out.write(file_in.read())
		


################################################################################

if __name__ == "__main__":
    main()
