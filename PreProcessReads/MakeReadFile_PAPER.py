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

def get_fasta_path(exp, barcode, lane):
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


# print(int(JOB_SPLIT/10))
def check_read(reads, spikes, barcode, tag, experiment, k_phi_x = k_phi_x):
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
    output_reads = []
    output_spikes = []
    tag_list = tag.split(",")
    read_len = 40 - len(tag_list[0])
     # Assign each of the four lines of the read.
    for i, read in enumerate(reads):
        if i%int(JOB_SPLIT/10) == 0:
            print(i)
            sys.stdout.flush()
        [header, full_seq, header2, score] = read
        seq = full_seq[:read_len]
        print(seq)
        seq_tag = full_seq[read_len:40]
        # Assign the multiplexing barcode sequence
        read_barcode = header.split("#")[1].split("/")[0]
        # List of conditionals to except the read:
        # 1. That the first line doesn't end with 1;0
        # q0  = header[-1] == 0
        # qN  = "N" in seq
        # qB  = "B" in score
        # qM = read_barcode != barcode
        # qX = (seq in k_phi_x or seq in k_phi_x_rc)
        # qS = sum([test_read_match(seq, j, mis = 1, indel = 1) for j in spikes])
        # qT = seq_tag not in tag_list
        # print([int(i) for i in [q0, qN, qB, qM, qX, qS, qT]])
        # print(seq_tag)
        # print(tag.split(","))
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
        elif True in [test_read_match(seq, j, mis = 1, indel = 1) for j in spikes]:
            output_spikes.append(full_seq[:40])
            # print([test_read_match(seq, j, mis = 1, indel = 1) for j in spikes])
            counts_readtypes_map["S"] += 1
        elif seq_tag not in tag_list:
            counts_readtypes_map["T"] += 1
        else:
            output_reads.append(full_seq[:40])
            counts_readtypes_map["P"] += 1
    summary_df = pd.DataFrame.from_dict(
        counts_readtypes_map, orient="index").reindex(
        ["0", "N", "B", "M", "X", "S", "T", "P"], copy=False)
    return output_reads, output_spikes, summary_df


def main():
    time_start = time.time()

    # Define all the relevant arguments
    arguments = ["miRNA","experiment","condition", "-nb_binary", "-rep",
                 "-test_binary"]
    mirna, experiment, condition, nb, rep, test = parse_arguments(arguments)
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
                    row["Exp"] != "CANCELLED" and
                    (rep == None or row["Rep"] == rep)):
                exp_file = row
    # Retrieve the path to the file, from the data in the appropriate
    # row of the database.
    exp, barcode, lane = [exp_file[i] for i in ["Exp", "Barcode", "Lane"]]
    file_path = get_fasta_path(exp, barcode, lane)
    print(file_path)
    if nb:
        experiment = experiment + "_nb"
    if rep:
        ext = "," + rep
    else:
        ext = ""
    reads_path = get_analysis_path(mirna,experiment,condition,"reads", ext = ext)
    reads_path_new = get_analysis_path(mirna,experiment,condition,"reads", ext = ext)

    spikes_path = get_analysis_path(mirna,experiment,condition,"spike_reads", ext = ext)
    summary_path = get_analysis_path(mirna,experiment,condition,"summary", ext = ext)

    # Retrieve the spike used, and the 3' tag used to ensure that the
    # reads pass quality control.
    spikes, tag = [exp_file[i] for i in ["Spike", "3prime"]]
    spike_seqs = [seq_spike_map[i] for i in spikes.split(",")]
    # Take in the file using the multiprocess_file function in general.py.
    with tarfile.open(file_path, "r:gz") as tar_in:
        for member in tar_in.getmembers():
            print(member)
            file_in = tar_in.extractfile(member)
            if file_in:
                if test:
                    results = multiprocess_test(file_in,
                                                get_read,
                                                100,
                                                check_read,
                                                1000,
                                                spike_seqs,
                                                barcode,
                                                tag,
                                                experiment)
                else:
                    results = multiprocess_file(file_in,
                                                get_read,
                                                int(1e6),
                                                check_read,
                                                spike_seqs,
                                                barcode,
                                                tag,
                                                experiment)
    # Collect the output from each thread, which are the list
    # of reads and the dictionary of read types.
    read_seqs = [j for i in results for j in i[0]]
    spike_seqs = [j for i in results for j in i[1]]
    pds = [i[2] for i in results]
    pd_summary = merge_data_frames(pds)
    # Flatten the list into one list (found on stackoverflow) and write
    # it to its output file.
    with open(reads_path,"wb") as file_out:
        file_out.write("".join([i + "\n" for i in read_seqs]))

    with open(spikes_path,"wb") as file_out:
        file_out.write("".join([i + "\n" for i in spike_seqs]))

    print(pd_summary)
    pd_summary.to_csv(summary_path, sep="\t", header=False)

    print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()

