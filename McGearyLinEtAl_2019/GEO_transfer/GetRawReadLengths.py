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

raw_data_path = "/lab/solexa_bartel/mcgeary/AgoRBNS/GEO_McGearyLin2019/raw_files/"


file = open(geo_path)

# Iterate over the lines 1-95
header = file.readline()
print("Loop #1, through small RNA-seq, AGO-RBNS, and massively parallel reporter assay.")

len_check = 2
# len_check = 2

def move_tar_file(in_path, outpath):
    call = 'tar -xOzf %s | gzip -f > %s' %(in_path, outpath)
    subprocess.call([call], shell=True)

def move_gzip_file(in_path, outpath):
    call = 'cp %s %s' %(in_path, outpath)    
    subprocess.call([call], shell=True)


i = 0
line = file.readline()
while line:
    ## GET 100 RAW READS FROM THE RAW FILE PATH.
    # Get the path that will be used for the raw data, from the `.txt` file.
    line_split = line.split("\t")
    path_out = raw_data_path + line_split[3]
    # print(path_out)
    call = 'zcat %s | head -2 | tail -1 | wc -m' %(path_out)
    check = subprocess.call([call], shell=True)
    i += 1
    line = file.readline()



