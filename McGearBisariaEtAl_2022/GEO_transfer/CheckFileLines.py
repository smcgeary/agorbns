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

geo_table_path = "ThreePrimeTargetingPaper/GEO_transfer/GEO_table_new2.txt"

geo_out_path = "/lab/solexa_bartel/mcgeary/GEO_McGearyBisaria2021/raw_files/"

geo_table_file = open(geo_table_path, encoding="ISO-8859-1")


header = geo_table_file.readline()
i = 0
line = geo_table_file.readline()
raw_file_names = []
while line:
# while i < 24:
    line_split = line.strip().split("\t")
    path_out_1 = line_split[-1]
    if path_out_1[0] != '"':
        path_out_1 = '"' + path_out_1 + '"'
    path_out = geo_out_path + path_out_1
    shell_call = "zcat %s | head" %(path_out)

    p1 = subprocess.Popen([shell_call], shell=True, stdout = subprocess.PIPE)
    # Parse the output
    out_pre = p1.communicate()[0].decode('utf-8').split("\n")
    output = len(out_pre[1])
    sequencer = out_pre[0]
    # print(sequencer)
    print(path_out_1 + "\t" + str(output))
    i += 1
    line = geo_table_file.readline()


