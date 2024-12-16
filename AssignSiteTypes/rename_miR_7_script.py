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
		if "pre_8mer-bA8" in path_in:
			file_split = path_in.split("/")
			pre_ex = file_split[:-1]
			file = file_split[-1]
			ex = file.split("_pre_8mer-bA8")[0].split("_")[-1].split(".")[0][2:]
			if ex:
				ex = int(ex)
			else:
				ex = 0
			if ex > 3:
				print("will be renamed")
				new_ex = str(ex + 1)
				file_base = file.split("_ex")[0]
				path_out = "/".join(pre_ex) + "/" + file_base + "_ex" + new_ex + ".txt"
				print(path_in)
				print(path_out)
				with open(path_in , "r+") as file_in:
					with open(path_out, "w+") as file_out:
						file_out.write(file_in.read())

		# # This gets the "_ex8" part of the file out, to check the conditionals.
		# file_name = path_in.split(".")[0]
		# path_out = "%s_pre_8mer-bA8.txt" %(file_name)
		# print(path_in)
		# print(path_out)
		


################################################################################

if __name__ == "__main__":
    main()
