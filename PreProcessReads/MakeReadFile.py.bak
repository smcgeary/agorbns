################################################################################
#MakeReadFile.py
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

# This script reads one library and outputs a .txt file to be used in all
# downstream applications. Inputs are:

# mirna: The miRNA in the experiment to be processed.
# experiment: The type of experiment to be processed.
# condition: The sampel type within the experiment to be processed.

# Functions from general this uses:


# Constants

k_phi_x = "".join(open(("/lab/bartel1_ata/mcgeary/computation/"
                        "AgoRBNS/general/phi_x.txt")).read().split("\r"))





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
    print("in check_read")
    sys.stdout.flush()
    print(k_phi_x)
    k_phi_x_rc = get_rc(k_phi_x)
    counts_readtypes_map = dict()
    for i in ["0","N", "B", "M", "X", "T", "S","P"]:
        counts_readtypes_map[i] = 0
    output_reads = []
    output_spikes = []
    tag_list = tag.split(",")
    print(tag_list)
    sys.stdout.flush()
    print(experiment)
    if experiment == "twist_reporter_assay":
        read_len = 100
    else:
        read_len = 40 - len(tag_list[0])
     # Assign each of the four lines of the read.
    i = 0
    total_reads = len(reads)/4
    print("Total reads: %s" %total_reads)
    sys.stdout.flush()
    for header, full_seq, header2, score in izip(*[iter(reads)]*4):
        if i%(total_reads/10) == 0:
            print(i/float(total_reads))
            sys.stdout.flush()
        seq = full_seq[:read_len]
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
        elif True in [test_read_match(seq, j, mis = 1, indel = 1)
                      for j in spikes]:
            output_spikes.append(full_seq[:40])
            # print([test_read_match(seq, j, mis = 1, indel = 1) for j in spikes])
            counts_readtypes_map["S"] += 1
        elif seq_tag not in tag_list and experiment != "twist_reporter_assay":
            counts_readtypes_map["T"] += 1
        else:
            output_reads.append(full_seq[:40])
            counts_readtypes_map["P"] += 1
        i += 1
    summary_df = pd.DataFrame.from_dict(
        counts_readtypes_map, orient="index").reindex(
        ["0", "N", "B", "M", "X", "S", "T", "P"], copy=False)
    print("about to return output")
    return output_reads, output_spikes, summary_df


def main():
    time_start = time.time()

    # Define all the relevant arguments
    arguments = ["miRNA","experiment","condition", "-rep",
                 "-jobs", "-test_binary"]
    mirna, experiment, condition, rep, jobs, test = parse_arguments(arguments)

    nb_tp_flag = experiment.split("_")[-1]
    if nb_tp_flag == "nb":
        nb = True
        tp = False
    elif nb_tp_flag == "tp":
        tp = True
        nb = False
    else:
        nb = False
        tp = False


    # Retrieve the path to the file, from the data in the appropriate
    # row of the database.
    file_path, spikes, tag, barcode = get_exp_info(mirna, experiment, condition,
                                                   rep=rep, nb=nb, tp=tp)
    spike_seqs = [seq_spike_map[i] for i in spikes.split(",") if i != "NA"]
    if rep:
        ext = "," + rep
    else:
        ext = ""
    if test:
        ext = "%s_test" %(ext)
    if jobs:
        jobs = int(jobs)
    else:
        jobs = 20
    reads_path, spikes_path, summary_path = [get_analysis_path(mirna,
                                                               experiment,
                                                               condition,
                                                               i, ext = ext)
                                             for i in ["reads", "spike_reads",
                                                       "summary"]]
    # Retrieve the spike used, and the 3' tag used to ensure that the
    # reads pass quality control.
    # Take in the file using the multiprocess_file function in general.py.
    args  = [spike_seqs, barcode, tag, experiment]
    if not jobs:
        jobs = 19

    print(file_path)
    print(reads_path)
    print(spikes_path)
    print(summary_path)
    threads = multiproc_file(file_path, jobs, check_read, test, *args)

    # Collect the output from each thread, which are the list
    # of reads and the dictionary of read types.
    read_seqs = [j for i in threads for j in i[0]]
    spike_seqs = [j for i in threads for j in i[1]]
    pds = [i[2] for i in threads]
    pd_summary = merge_data_frames(pds)
    # Flatten the list into one list (found on stackoverflow) and write
    # it to its output file.
    print(reads_path)
    with open(reads_path,"wb") as file_out:
        file_out.write("".join([i + "\n" for i in read_seqs]))
    with open(spikes_path,"wb") as file_out:
        file_out.write("".join([i + "\n" for i in spike_seqs]))
    print(pd_summary)
    pd_summary.to_csv(summary_path, sep="\t", header=False)

    print_time_elapsed(time_start)
    print("done script")


################################################################################

if __name__ == "__main__":
    main()

