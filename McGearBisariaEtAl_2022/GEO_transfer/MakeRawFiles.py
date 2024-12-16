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


def move_tar_file(in_path, outpath):
    call = 'tar -xOzf %s | gzip -f > %s' %(in_path, outpath)
    subprocess.call([call], shell=True)

def move_gzip_file(in_path, outpath):
    call = 'cp %s %s' %(in_path, outpath)    
    subprocess.call([call], shell=True)


header = geo_table_file.readline()
i = 0
line = geo_table_file.readline()
raw_file_names = []
while line:
# while i < 24:
    # print("%s " %(i + 1) + "_"*80)
    # print(line.strip())
    line_split = line.strip().split("\t")
    path_in = line_split[0]
    # print(line_split[-2:])
    path_out_1 = line_split[-1]
    # print(path_out_1)
    if path_out_1[0] != '"':
        path_out_1 = '"' + path_out_1 + '"'
    path_out = geo_out_path + path_out_1
    if i == 35:
        print(path_in)
        print(path_out)
    if not os.path.exists(path_in):
        print("doesn't exist")
        print(i)
        print(path_in)
        print(path_out)
    if path_in not in raw_file_names:
        raw_file_names.append(path_in)
    else:
        print("redundant name!")
        break
    if path_in[-7:] == ".tar.gz":
        tar_func = True
    else:
        tar_func = False
    # Write the file as a `.txt.gz` file.
    time_start = time.time()
    if i == 35:
        if tar_func:
            "Untar-ing:"
            move_tar_file(path_in, path_out)
        else:
            "Moving:"
            move_gzip_file(path_in, path_out)
        print_time_elapsed(time_start)
    i += 1
    line = geo_table_file.readline()


