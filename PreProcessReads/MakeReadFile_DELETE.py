################################################################################
#MakeReadFile.py
################################################################################
import csv # Used to open Librarys_SL.csv
import gzip # Used to open the fa.gz file
import imp # Used to import general.py
import time
import tarfile
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor
from itertools import repeat
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
import RBNS_methods
from RBNS_methods import *

# This script reads one library and outputs a .txt file to be used in all
# downstream applications. Inputs are:

# mirna: The miRNA in the experiment to be processed.
# experiment: The type of experiment to be processed.
# condition: The sampel type within the experiment to be processed.

# Functions from general this uses:


# Constants

k_phi_x = "".join(open(("/lab/bartel1_ata/mcgeary/computation/"
                        "AgoRBNS/general/phi_x.txt")).read().split("\r"))
k_phi_x_rc = get_rc(k_phi_x)


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

def get_fasta_path_new(exp,barcode,lane):
    """Retrieves the complete path to the zipped fastq file 
        in solexa_bartel/mcgeary.

        Args:
            exp: The folder in which the data is, representing the experiment.
            barcode: The multiplexing barcode associated with the sample.
            lane: The lane the sample sequenced on within the sequencing run.

        Returns:
            The string representing the full path to the zipped fastq file.
    """
    path = "/lab/solexa_public/Bartel/%s/QualityScore/" % (exp)
    file = "%s-s_%s_1_sequence.txt.tar.gz" % (barcode, lane)
    full_path = path + file
    return full_path


def get_fasta_path_burge(exp,srr_file):
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
    file = "%s.fastq.gz" % (srr_file)
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
    out = [file.readline().strip() for i in range(4)]
    return out


def check_read(reads, spikes, barcode, tag, experiment):
    """Takes a list of reads and outputs a list of reads for which all reads
        pass the criteria specificed by the arguments 'spikes', 'barcode',
        'tag' and 'k_phi_x'.

    Args:
        file: The file object being read.

    Returns:
        A tuple containing the list of reads and a dictionary of read types.
    """
    counts_readtypes_map = dict()
    counts_readtypes_map = {i:0 for i in ["0", "N", "B", "M", "X", "T", "S","P"]}
    indeces = [False]*len(reads)
    read_len = 40 - max([len(i) for i in tag.split(",")])
     # Assign each of the four lines of the read.
    for i, read in enumerate(reads):
        print(read.header1)
        print(read.N)
        print(read.bin_score)
        print(read.B)
        print(read.qc)
        print(read.multi)
        print(read.seq)
        print(tag)
        print(read.check_library_barcode(read_len, tag))
        print("read will pass:")
        print(read.check_multiplex_barcode(barcode))
        # print(k_phi_x)
        print(read.seq)
        print(LCSubString("ABAABABAAB", "BAABABBAABA"))
        check = LCSubString(read.seq, spikes[0])
        print(read.seq)
        print(spikes[0])
        check_ln = check[0]
        pos_1 = check[1][0]
        pos_2 = check[1][1]
        print([read.seq[i:i+check_ln] for i in pos_1])
        print([spikes[0][i:i+check_ln] for i in pos_2])
        # Assign the multiplexing barcode sequence
        # List of conditionals to except the read:
        # 1. That the first line doesn't end with 1;0

        print(read.check_if_passes(barcode, read_len, tag))
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
            counts_readtypes_map["S"] += 1
        elif seq_tag not in tag.split(","):
            counts_readtypes_map["T"] += 1
        else:
            indeces[i] = True
            counts_readtypes_map["P"] += 1
    output_reads = [read[1][:40] for i, read in enumerate(reads) if indeces[i]]

    return output_reads, counts_readtypes_map


def main():
    time_start = time.time()

    # Define all the relevant arguments
    arguments = ["miRNA","experiment","condition", "-nb_binary", "-cb_binary","-rep"]
    mirna, experiment, condition, nb, cb, rep = parse_arguments(arguments)
    if nb:
        csv_file = "Libraries_NB.csv"
    elif cb:
        csv_file = "Libraries_CB.csv"
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
                    row["Exp"] != "CANCELLED" and
                    (rep == None or row["Rep"] == rep)):
                print("check!")
                exp_file = row
    # Retrieve the path to the file, from the data in the appropriate
    # row of the database.
    exp, barcode, lane = [exp_file[i] for i in ["Exp", "Barcode", "Lane"]]
    file_path = get_fasta_path(exp, barcode, lane)
    if nb:
        experiment = experiment + "_nb"
    if rep:
        reads_path = get_analysis_path(mirna,experiment,condition,"full_reads", ext = "," + rep)
        summary_path = get_analysis_path(mirna,experiment,condition,"full_summary", ext = "," + rep)

    else:
        reads_path = get_analysis_path(mirna,experiment,condition,"full_reads")
        summary_path = get_analysis_path(mirna,experiment,condition,"full_summary")

    # Retrieve the spike used, and the 3' tag used to ensure that the
    # reads pass quality control.
    spikes, tag = [exp_file[i] for i in ["Spike", "3prime"]]
    spike_seqs = [seq_spike_map[i] for i in spikes.split(",")]
    print(file_path)
    # Take in the file using the multiprocess_file function in general.py.


    # # with tarfile.open(file_path, "r") as tar_file:
    #     # print(file_in.readline())
    #     tar_info = tar_file.getmembers()[0]
    #     file_in = tar_file.extractfile(tar_info)
    with gzip.open(file_path, "rb") as file_in:


        results = multiprocess_test(file_in,
                                    extract_raw_read,
                                    int(1e6),
                                    check_read,
                                    10,
                                    spike_seqs,
                                    barcode,
                                    tag,
                                    experiment)

    # # Collect the output from each thread, which are the list
    # # of reads and the dictionary of read types.
    # read_threads = [i[0] for i in results]
    # dict_threads = [i[1] for i in results]

    # # Flatten the list into one list (found on stackoverflow) and write
    # # it to its output file.
    # read_seqs = [i for sublist in read_threads for i in sublist]
    # print(len(read_seqs))

    # p5_seq = "GGGCAGAGTTCTACAGTCCGACGATC"
    # if mirna == "miR-1" and experiment == "equilibrium":
    #     p3_seq = "TATGCCGTCTTCTGCTTG"
    # elif tag == "TG":
    #     p3_seq = "TTCGTATGCCGTCTTCTGCTTG"
    # else:
    #     p3_seq = "TCGTATGCCGTCTTCTGCTTG"
    # with open(reads_path,"wb") as file_out:
    #     file_out.write("".join(["%s%s%s\n" % (p5_seq, i, p3_seq) for i in read_seqs]))

    # # Construct the dictionary from the list of futures, and write it to
    # # its output file.
    # keys = dict_threads[0].keys()
    # values = [sum([thread[key] for thread in dict_threads]) for key in keys]
    # counts_readtypes_map = {key : value for (key,value) in zip(keys,values)}

    # # No longer have input our output files open. Now write summmary file
    # # to the same directory , ending in "summary.txt".
    # with open(summary_path,"wb") as file_out:
    #     keys = ["0","N", "B", "M", "X", "T", "S","P"]
    #     for key in keys:
    #         value = counts_readtypes_map[key]
    #         file_out.write("%s:\t%s\n" % (key, value))

    # # Print the amount of time the script took to complete.
    # print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()

