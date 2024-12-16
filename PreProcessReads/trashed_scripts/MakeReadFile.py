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
from general import get_analysis_path, multiprocess_file, print_time_elapsed
from general import *

# This script reads one library and outputs a .txt file to be used in all
# downstream applications. Inputs are:

# mirna: The miRNA in the experiment to be processed.
# experiment: The type of experiment to be processed.
# condition: The sampel type within the experiment to be processed.

# Functions from general this uses:


# Constants

k_phi_x = "".join(open(("/lab/bartel1_ata/mcgeary/computation/"
                        "AgoRBNS/general/phi_x.txt")).read().split("\r"))


# FUNCTIONS

def get_fasta_path(exp,barcode,lane):
    """Retrieves the complete path to the zipped fastq file 
        in solexa_bartel/mcgeary.

        Args:
            exp: The folder in which the data is, representing the experiment.
            barcode: The multiplexing barcode associated with the sample.
            lane: The lane the sample sequenced on within the sequencing run.

        Returns:
            The string representing the full path to the zipped fastq file.
    """
    path = "/lab/solexa_bartel/mcgeary/%s/fastq/" % (exp)
    file = "%s-s_%s_1_sequence.fa.gz" % (barcode, lane)
    full_path = path + file
    return full_path

def get_read(file):
    """Retrieves the next four lines of the file.

    Args:
        file: The file object being read.

    Returns:
        A list of four reads, corresponding to the header, the sequence,
        the header, and the quality scores corresponding to each base.
        As a consequence, the read file is updated to be four lines ahead.
    """
    out = []
    for i in range(4):
        line = file.readline().strip()
        if line:
            out.append(line)
        else:
            out = False
            break
    return out


def check_read(reads, spikes, barcode, tag, k_phi_x = k_phi_x):
    """Takes a list of reads and outputs a list of reads for which all reads
        pass the criteria specificed by the arguments 'spikes', 'barcode',
        'tag' and 'k_phi_x'.

    Args:
        file: The file object being read.

    Returns:
        A tuple containing the list of reads and a dictionary of read types.
    """
    k_phi_x_rc = get_rc(k_phi_x)
    counts_readtypes_map = dict()
    for i in ["0","N", "B", "M", "X", "T", "S","P"]:
        counts_readtypes_map[i] = 0
    indeces = [False]*len(reads)

     # Assign each of the four lines of the read.
    for i, read in enumerate(reads):
        [header, full_seq, header2, score] = read
        seq = full_seq[:37]
        seq_tag = full_seq[37:40]

        # Assign the multiplexing barcode sequence
        read_barcode = header.split("#")[1].split("/")[0]
        # List of conditionals to except the read:
        # 1. That the first line doesn't end with 1;0
        if header[-1] == 0:
            counts_readtypes_map["0"] += 1
        # 2. That "N" is not in the read sequence.
        elif "N" in seq:
            counts_readtypes_map["N"] += 1
        # 3. That "B" is not in the quality score.
        elif "B" in score:
            counts_readtypes_map["B"] += 1
        # 4. That the read has the correct multiplexing bardcode.
        elif read_barcode != barcode:
            counts_readtypes_map["M"] += 1
        # 5. That the sequence is not within a 1 nt mismatch to
        #   the phi-x genome, or its reverse complement.
        elif (seq in k_phi_x or seq in k_phi_x_rc):
            counts_readtypes_map["X"] += 1
        # 6. That the seqence is not a spike sequence.
        elif True in [test_match(seq, j, mis = 1) for j in spikes]:
            print([test_match(seq, j, mis = 1) for j in spikes])
            counts_readtypes_map["S"] += 1
        elif seq_tag != tag:
            counts_readtypes_map["T"] += 1
        else:
            indeces[i] = True
            # print(seq)
            counts_readtypes_map["P"] += 1

    output_reads = [read[1][:40] for i, read in enumerate(reads) if indeces[i]]

    return output_reads, counts_readtypes_map


def main():
    time_start = time.time()

    # Define all the relevant arguments
    arguments = ["miRNA","experient","condition", "-nb_binary"]
    mirna, experiment, condition, nb = parse_arguments(arguments)
    if nb:
        csv_file = "Libraries_NB.csv"
    else:
        csv_file = "Libraries_SL.csv"
    # Load the database of experiments:
    with open(csv_file, 'U') as master:
        library_list = csv.DictReader(master)
        # Determine the experiment to be analyzed
        exp_file = None
        for row in library_list:
            if (mirna in row["miRNA"].split(",") and
                    row["Exp_type"] in experiment and
                    row["Sample type"] == condition and
                    row["Exp"] != "CANCELLED"):
                exp_file = row

    # Retrieve the path to the file, from the data in the appropriate
    # row of the database.
    exp, barcode, lane = [exp_file[i] for i in ["Exp", "Barcode", "Lane"]]
    file_path = get_fasta_path(exp, barcode, lane)
    if nb:
        experiment = experiment + "_nb"
    reads_path = get_analysis_path(mirna,experiment,condition,"reads")
    summary_path = get_analysis_path(mirna,experiment,condition,"summary")

    # Retrieve the spike used, and the 3' tag used to ensure that the
    # reads pass quality control.
    spikes, tag = [exp_file[i] for i in ["Spike", "3prime"]]
    spike_seqs = [seq_spike_map[i] for i in spikes.split(",")]

    # Take in the file using the multiprocess_file function in general.py.
    with gzip.open(file_path, "rb") as file_in:
        results = multiprocess_file(file_in, get_read, int(1e6), check_read,
                                    spike_seqs, barcode, tag)

    # Collect the output from each thread, which are the list
    # of reads and the dictionary of read types.
    read_threads = [i[0] for i in results]
    dict_threads = [i[1] for i in results]

    # Flatten the list into one list (found on stackoverflow) and write
    # it to its output file.
    read_seqs = [i for sublist in read_threads for i in sublist]
    with open(reads_path,"wb") as file_out:
        file_out.write("".join(["%s\n" % (i) for i in read_seqs]))

    # Construct the dictionary from the list of futures, and write it to
    # its output file.
    keys = dict_threads[0].keys()
    values = [sum([thread[key] for thread in dict_threads]) for key in keys]
    counts_readtypes_map = {key : value for (key,value) in zip(keys,values)}

    # No longer have input our output files open. Now write summmary file
    # to the same directory , ending in "summary.txt".
    with open(summary_path,"wb") as file_out:
        keys = ["0","N", "B", "M", "X", "T", "S","P"]
        for key in keys:
            value = counts_readtypes_map[key]
            file_out.write("%s:\t%s\n" % (key, value))

    # Print the amount of time the script took to complete.
    print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()

