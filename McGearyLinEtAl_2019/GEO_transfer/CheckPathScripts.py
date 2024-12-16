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

len_check = 71
# len_check = 2

len_skip = 67
for i in range(len_skip):
    line = file.readline()

for i in range(len_skip, len_check):
    line = file.readline()
    print(line)
    ## GET 100 RAW READS FROM THE RAW FILE PATH.
    # Get the path that will be used for the raw data, from the `.txt` file.
    print("%s " %(i + 1) + "_"*80)
    line_split = line.split("\t")
    path = line_split[4]
    if not os.path.exists(path):
        print(line)
        break
    # Open the file object and extract the first 400 lines of the file
    if path[-7:] == ".tar.gz":
        tar_in = tarfile.open(path, "r:gz")
        member = tar_in.getmembers()[0]
        file_in = tar_in.extractfile(member)
    elif path[-3:] == ".gz":
        file_in = gzip.open(path, "rt")
    else:
        file_in = open(path, "r")
    if i > 39:
        range_reads = 400000
    else:
        range_reads = 2000
    raw_reads = []
    # raw_reads = file_in.readlines()
    while len(raw_reads) < 100:
        lines = [file_in.readline() for j in range(4)]
        try:
            lines = [line.decode("utf-8") for line in lines]
        except (UnicodeDecodeError, AttributeError):
            pass
        seq_line = lines[1]
        if "N" not in seq_line:
            raw_reads.append(seq_line[:37])
    print("raw_reads")
    for read in raw_reads:
        print(read)
    # Reiterate through the raw reads to pick the 2nd line and covert it from a
    # byte object to a character object.
    # raw_reads = [line.strip() 
    #              for (i, line) in enumerate(raw_reads) if i%4 == 1]
    # raw_reads = [line[:37] for line in raw_reads if "N" not in line]
    # raw_reads = raw_reads[:100]
    if path[-7:] == ".tar.gz":
        tar_in.close()
    else:
        file_in.close()
    ## 2. GET 10 PROCESSED READS.
    # Get the path to the processed reads in the AGO-RBNS folder system.
    experiment, condition, rep, mirna = line_split[11:15]
    if len(rep) != 0:
        condition = "%s,%s" %(condition, rep)
    mirna = mirna.strip().split("_")[0]
    reads_path = get_analysis_path(mirna, experiment, condition, "reads")
    # Get the first ten processed reads.
    file_in = open(reads_path, "r")
    reads_check = [file_in.readline().strip()[:37] for i in range(10)]
    print("processed reads:")
    for read in reads_check:
        print(read)
    # Check that all ten processed reads are found in the first 100 reads
    # of the raw reads.
    num_found_in_raw_reads = np.sum([read in raw_reads for read in reads_check])
    if num_found_in_raw_reads == 10:
        print("SAMPLE %s: OK!" %(i + 1))
    else:
        print("SAMPLE %s: FAILED!" %(i + 1))
        print(num_found_in_raw_reads)
        print(path)
        print("First 10 raw reads:")
        print("\n".join(raw_reads[:10]))
        print(reads_path)
        print("Processed reads:")
        print("\n".join(reads_check))


# print("Loop #2, through HeLa and HEK283FT transfection data.")
# for line in file.readlines():
#   line_split = line.split("\t")
#   raw_path = line_split[4]
#   print(raw_path)
#   if not os.path.exists(raw_path):
#       print(line)


