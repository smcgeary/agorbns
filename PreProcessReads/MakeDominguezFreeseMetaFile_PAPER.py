################################################################################
#MakeReadFile.py
################################################################################
import csv # Used to open Librarys_SL.csv
import gzip # Used to open the fa.gz file
import imp # Used to import general.py
import time
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor
from itertools import repeat
imp.load_source("general",
                ("/lab/bartel1_ata/mcgeary/computation/"
                 "AgoRBNS/general/general.py"
                )
               )
from general import get_analysis_path, multiproc_file, print_time_elapsed
from general import *

# Functions from general this uses:
def check_read(reads):
    """Takes a list of reads and outputs a list of reads for which all reads
        pass the criteria specificed by the arguments 'spikes', 'barcode',
        'tag' and 'k_phi_x'.

    Args:
        file: The file object being read.

    Returns:
        A tuple containing the list of reads and a dictionary of read types.
    """
    output_reads = []
    indeces = [False]*len(reads)
     # Assign each of the four lines of the read.
    for header, full_seq, header2, score in izip(*[iter(reads)]*4):
        if "N" not in full_seq:
            output_reads.append(full_seq.strip())
    return output_reads

def main():
    time_start = time.time()
    path = "/lab/solexa_bartel/mcgeary/RBNS_DominguezFreese/raw_data/"
    # This is the official metadata file from the ENCODE website.
    csv_file = path + "metadata.tsv"
    # This is the extra file I downloaded to deal with the lack of annotations
    # for the input and no-protein samples.
    report_file = path + "Experiment Report 2018_12_6.tsv"

    # First use the report file to assign file names to each of the input and
    # no-protein controls for each of the RBPs.

    # Initialize dictionaries:
    input_rbp_map = {}
    mock_rbp_map = {}

    with open(report_file, 'U') as master:
        next(master, None)
        library_list = csv.DictReader(master, delimiter="\t")
        for row in library_list:
            # Get the textual description of each file from the file, and the 
            # the string of file names.
            desc = row["Description"].split("RBNS) ")[1]
            files = row["Files"].split("/")

            # Conditional making sure the file list isn't empty:
            if len(files) > 1:
                # Manually identified this is the right index.
                files = files[2]
                # String parsing of the description to get the names of the
                # RBPs and as well whether the row is an input sample or a mock
                # sample.
                desc_split = desc.split(" ")
                desc_cond = desc_split[0]
                if desc_cond == "input":
                    rbp = desc_split[-1]

                    # Conditional checking for 'RBP1/RBP2' formatting, which is
                    # actually two different RBPs.
                    if "/" in rbp:
                        rbps = rbp.split("/")
                    else:
                        rbps = [rbp]
                    for rbp in rbps:
                        input_rbp_map[rbp] = files

                    # Conditional looking for an "and" in the description,
                    # because this also is evidence that it is two different
                    # RBPs.
                    desc_split = desc_split[:-1]
                    if desc_split[-1] in ["or", "and"]:
                        rbp2 = desc_split[-2]
                        if "/" in rbp2:
                            rbps = rbp2.split("/")
                        else:
                            rbps = [rbp2]
                        for rbp in rbps:
                            input_rbp_map[rbp] = files


                if desc_cond == "mock":
                    rbp = desc_split[-1]
                    # Conditional checking for 'RBP1/RBP2' formatting, which is
                    # actually two different RBPs.
                    if "/" in rbp:
                        rbps = rbp.split("/")
                    else:
                        rbps = [rbp]
                    for rbp in rbps:
                        mock_rbp_map[rbp] = files
                    # Conditional looking for an "and" in the description,
                    # because this also is evidence that it is two different
                    # RBPs.
                    desc_split = desc_split[:-1]
                    if desc_split[-1] == "and":
                        rbp2 = desc_split[-2]
                        if "/" in rbp2:
                            rbps = rbp2.split("/")
                        else:
                            rbps = [rbp2]
                        for rbp in rbps:
                            mock_rbp_map[rbp] = files

    samples_rbp_map = {}
    rbp_exp_map = {}
    exp_rbp_map = {}

    # Load the database of experiments:
    with open(csv_file, 'U') as master:
        library_list = csv.DictReader(master, delimiter="\t")
        # Determine the experiment to be analyzed
        exp_file = None
        i = 0
        for row in library_list:
            if row["File format"] == "fastq":
                rbp, exp = [row[key] for key in
                            ["Experiment target", "Experiment accession"]]

                if rbp not in rbp_exp_map and " " not in rbp:
                    rbp = rbp.split("-")[0]
                    exp_rbp_map[rbp] = exp
                    rbp_exp_map[exp] = rbp


        library_list = csv.DictReader(master, delimiter="\t")
        # Determine the experiment to be analyzed
        exp_file = None
        i = 0

    samples_rbp_map = {i : [] for i in exp_rbp_map.keys()}

    with open(csv_file, 'U') as master:
        library_list = csv.DictReader(master, delimiter="\t")
        for row in library_list:
            if row["File format"] == "fastq":
                conc, exp, file = [row[key] for key in
                                   ["RBNS protein concentration",
                                    "Experiment accession",
                                    "File accession"]]
                if exp in rbp_exp_map:
                    rbp = rbp_exp_map[exp].split("-")[0]
                    samples_rbp_map[rbp] += [[conc, file]]

    # Define all the relevant arguments
    arguments = ["rbp","condition", "-test_binary"]
    rbp, condition, test = parse_arguments(arguments)


    if condition == "I":
        file = input_rbp_map[rbp]
        conc = 0
    else:
        sample_list = samples_rbp_map[rbp]
        sample_list = [[i[0].split(" nM")[0], i[1]] for i in sample_list]
        sample_list_sorted = sorted(sample_list,
                          key=lambda sample: int(sample[0]))
        file = sample_list_sorted[int(condition) - 1][1]
        conc = int(sample_list_sorted[int(condition) - 1][0])

    file_path = path + file + ".fastq.gz"
    reads_path = get_analysis_path_burge(rbp, "equilibrium", condition, "reads")

    # Take in the file using the multiprocess_file function in general.py.
    results = multiproc_file(file_path, 10, check_read, test)
    # Flatten the list into one list (found on stackoverflow) and write
    # it to its output file.
    read_seqs = [i for sublist in results for i in sublist]
    print(len(read_seqs))
    with open(reads_path,"wb") as file_out:
        file_out.write("".join(["%s\n" %(i) for i in read_seqs]))

    print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()

