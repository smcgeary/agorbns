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

geo_table_path = "2017_Paper/GEO_transfer/GEO_table.txt"

geo_out_path = "/lab/solexa_bartel/mcgeary/GEO_McGearyLin2019/raw_files/"


geo_table_file = open(geo_path)

# Two alternate functions for moving the file to the correct directory, one that
# untars and then gunzips the file at a new place, and the other that copies the
# file to a new directory.
def move_tar_file(in_path, outpath):
    call = 'tar -xOzf %s | gzip -f > %s' %(in_path, outpath)
    subprocess.call([call], shell=True)

def move_gzip_file(in_path, outpath):
    call = 'cp %s %s' %(in_path, outpath)    
    subprocess.call([call], shell=True)


header = geo_table_file.readline()
i = 0
line = geo_table_file.readline()
while line:
    print("%s " %(i + 1) + "_"*80)
    line_split = line.split("\t")
    path_in = line_split[4]
    path_out = raw_data_path + line_split[3]
    print(path_in)
    print(path_out)
    if not os.path.exists(path_in):
        print(line)
        break
    if path_in[-7:] == ".tar.gz":
        tar_func = True
    else:
        tar_func = False
    # Write the file as a `.txt.gz` file.
    time_start = time.time()
    if tar_func:
        "Untar-ing:"
        move_tar_file(path_in, path_out)
    else:
        "Moving:"
        move_gzip_file(path_in, path_out)
    print_time_elapsed(time_start)
    i += 1
    line = geo_table_file.readline()



