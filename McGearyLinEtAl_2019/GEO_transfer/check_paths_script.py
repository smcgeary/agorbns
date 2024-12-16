################################################################################
#CheckPathScripts.py
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


geo_path = "2017_Paper/GEO_transfer/GEO_table.txt"

file = open(geo_path)

# Iterate over the lines 1-95
header = file.readline()
print("Loop #1, through small RNA-seq, AGO-RBNS, and massively parallel reporter assay.")
for i in range(96):
	line = file.readline()
	# print(line)
	line_split = line.split("\t")
	print(line_split)
	raw_path = line_split[4]
	print(raw_path)
	if not os.path.exists(raw_path):
		print(line)
	# get the processed reads path
	reads_path = get_analysis_path(i_input[0], i_input[1], i_input[2], "reads")


print("Loop #2, through HeLa and HEK283FT transfection data.")
for line in file.readlines():
	line_split = line.split("\t")
	raw_path = line_split[4]
	print(raw_path)
	if not os.path.exists(raw_path):
		print(line)


