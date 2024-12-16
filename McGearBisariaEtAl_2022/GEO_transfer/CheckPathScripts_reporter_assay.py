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
from subprocess import Popen, PIPE

geo_path = "ThreePrimeTargetingPaper/GEO_transfer/GEO_table_reporter_assay.txt"
geo_out_path = "/lab/solexa_bartel/mcgeary/GEO_McGearyBisaria2021_3/raw_files/"
    


file = open(geo_path, encoding="ISO-8859-1")

# Iterate over the lines 1-95
header = file.readline()
header_names = header.strip().split("\t")
print(header_names)

header_names_dict = {j : i for (i, j) in enumerate(header_names)}

print(header_names_dict)


print("Loop #1, through small RNA-seq, AGO-RBNS, and massively parallel reporter assay.")

len_check = 71
len_check = 2
len_check = 8

for i in range(0, len_check):
    line = file.readline()
    line_split = line.strip().split("\t")
    ## GET 100 RAW READS FROM THE RAW FILE PATH.
    # Get the path that will be used for the raw data, from the `.txt` file.
    print("%s " %(i + 1) + "_"*80)
    line_split = line.split("\t")
    path_original = line_split[0]
    file_path = line_split[3].strip()
    out_path = geo_out_path + file_path
    path = out_path
    if not os.path.exists(path):
        print("doesn't exist")
        print("|" + path + "|")
        break
    ##############################1 
    # Open the file object and extract the first 400 lines of the file
    if path[-7:] == ".tar.gz":
        tar_in = tarfile.open(path, "r:gz")
        member = tar_in.getmembers()[0]
        file_in = tar_in.extractfile(member)
    elif path[-3:] == ".gz":
        file_in = gzip.open(path, "rt")
    else:
        file_in = open(path, "r")
    range_reads = 3000
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
            raw_reads.append(seq_line[:120])
    # # Reiterate through the raw reads to pick the 2nd line and covert it from a
    # # byte object to a character object.
    # raw_reads = [line.strip() 
    #              for (i, line) in enumerate(raw_reads) if i%4 == 1]
    # raw_reads = [line[:120] for line in raw_reads if "N" not in line]
    # raw_reads = raw_reads[:1000]
    if path[-7:] == ".tar.gz":
        tar_in.close()
    else:
        file_in.close()
    #################################### 2
    # Make sure the file on solexa_public has the same reads.
    path = path_original
    if not os.path.exists(path):
        print("doesn't exist")
        print("|" + path + "|")
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
    range_reads = 2000
    raw_solexa_reads = []
    while len(raw_solexa_reads) < 100:
        lines = [file_in.readline() for j in range(4)]
        try:
            lines = [line.decode("utf-8") for line in lines]
        except (UnicodeDecodeError, AttributeError):
            pass
        seq_line = lines[1]
        if "N" not in seq_line:
            raw_solexa_reads.append(seq_line[:120])
    
    if path[-7:] == ".tar.gz":
        tar_in.close()
    else:
        file_in.close()

    geo_test = raw_reads == raw_solexa_reads
    print("Are the GEO reads the same as the solexa_reads? -%s" %geo_test)
    ## 2. GET 10 PROCESSED READS.
    # Get the path to the processed reads in the AGO-RBNS folder system.
    experiment = "twist_reporter_assay_3p_2_tp"
    mirna, rep = out_path.split("/")[-1].split(".")[:-2][0].split("_")[1:]
    rep = rep.split("rep")[1]
    if mirna == "mock":
        mirna = "miR-1"
        condition = "no_duplex_parallel"
    else:
        condition = "duplex_parallel"
    reads_path = get_analysis_path(mirna, experiment, condition, "reads", ext="," + rep)
    test_paths = reads_path == line_split[1].replace('"', '')
    print("Is the local path the same as the manually determined path? -%s" %test_paths)
    if not test_paths:
        print(reads_path)
        print(line_split[1])
        break
    # Get the first ten processed reads.
    file_in = open(reads_path, "r")
    reads_check = [file_in.readline().strip()[:120] for i in range(10)]
    # Check that all ten processed reads are found in the first 100 reads
    # of the raw reads.
    num_found_in_raw_reads = np.sum([read in raw_reads for read in reads_check])
    if num_found_in_raw_reads == 10:
        print("SAMPLE %s: OK! - First 10 processed reads are within %s raw reads" %(i + 1, len(raw_reads)))
    else:
        print("SAMPLE %s: FAILED!" %(i + 1))
        print(num_found_in_raw_reads)
        print(path)
        print("First 10 raw reads:")
        print("\n".join(raw_reads[:10]))
        print(reads_path)
        print("Processed reads:")
        print("\n".join(reads_check))
        print(out_path)
        grep_ints = []
        for read_check in reads_check:
            # command = "zgrep -n %s %s | head -1" %(read_check, out_path)
            p1 = Popen(["zgrep -n -m 1 %s %s" %(read_check, out_path)], shell=True, stdout=PIPE)
            p2 = Popen(['head -1'], shell=True, stdin=p1.stdout, stdout=PIPE)
            final_output = p2.communicate()[0].decode().split(":")[0]
            print(final_output)
            grep_ints += [int(final_output)]
        print(grep_ints)
        for i, int_i in enumerate(grep_ints[:-1]):
            if int_i >= grep_ints[i+ 1]:
                print("out of order")
                break
        print("indeces were in order.")


# print("Loop #2, through HeLa and HEK283FT transfection data.")
# for line in file.readlines():
#   line_split = line.split("\t")
#   raw_path = line_split[4]
#   print(raw_path)
#   if not os.path.exists(raw_path):
#       print(line)


