################################################################################
#MakeReadFile.py
################################################################################
import imp # Used to import general.py
imp.load_source("general", "general/general.py")
# imp.load_source("RBNS_methods", "general/RBNS_methods.py")
from general import *
# from RBNS_methods import *

# This script reads one library and outputs a .txt file to be used in all
# downstream applications. Inputs are:

# mirna: The miRNA in the experiment to be processed.
# experiment: The type of experiment to be processed.
# condition: The sampel type within the experiment to be processed.

# Functions from general this uses:


# Constants
k_phi_x = "".join(open(("general/phi_x.txt")).read().split("\n"))

def check_read(reads, spikes, barcode, tag, experiment, raw_original,
               filt_indices, byte_convert, global_index, k_phi_x = k_phi_x):
    ################################################################################
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
    for i in ["0", "N", "B", "M", "X", "T", "S", "P"]:
        counts_readtypes_map[i] = 0
    output_reads = []
    output_spikes = []
    output_inds = []
    tag_list = tag.split(",")
    if "twist_reporter_assay" in experiment:
        read_len = 100
    else:
        read_len = 40 - len(tag_list[0])
     # Assign each of the four lines of the read.
    i = 0
    total_reads = len(reads)/4
    print(("Total reads: %s" %total_reads))
    sys.stdout.flush()
    if byte_convert:
        def get_four_reads(read_tup):
            return(i.decode("utf-8") for i in read_tup)
    else:
        def get_four_reads(read_tup):
            return(i for i in read_tup)
    for read_tup in zip(*[iter(reads)]*4):
        header, full_seq, header2, score = get_four_reads(read_tup)
        sys.stdout.flush()
        if i%(total_reads/10) == 0:
            print((i/float(total_reads)))
            sys.stdout.flush()
        seq = full_seq[:read_len]
        seq_tag = full_seq[read_len:40]
        if filt_indices is not None:
            # print("hi")
            # print(i)
            # print(global_index + i)
            # print(min(filt_indices))
            # print(sys.stdout.flush())
            if (global_index + i) in filt_indices:
                if "twist_reporter_assay_3p" in experiment:
                    output_reads.append(full_seq[:120])
                elif "twist_reporter_assay" in experiment:
                    output_reads.append(full_seq[:100])
                else:
                    output_reads.append(full_seq[:40])
                counts_readtypes_map["P"] += 1
            else:
                counts_readtypes_map["0"] += 1
        else:
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
            #
            # NOTE: The data as hosted by the SRA does not contain the original
            # information from the reads from the Whitehead Genome Technology CORE.
            # The original reads had the following structure:
            #
            # @WIGTC-HISEQ:4:1107:1514:1981#ACAGTG/1;0
            # NGCAGNNTTAGCATTCACAGACAACACCGCACACGAGGGC
            # +WIGTC-HISEQ:4:1107:1514:1981#ACAGTG/1;0
            # BOZ\_BBQY]^_]_^^_^_^^[^^_]^^^^[^[^^^^^^^
            #
            # where the "ACAGTG" in lines 1 and 3 refers to the sample index
            # barcode, and the "0" at the end of lines 1 and 3 refers to whether
            # the read passed the Illumina 1.5 pipeline QC filter, with 0 indicating
            # the read failed and 1 indicating it passed. Furthermore, the reads
            # used the Illumina 1.5-era Phred+64 quality scores, in which an
            # upper-case "B" indicates either a stretch of low-quality bases at the
            # 3′ end of the read, or low-quality base of unknown quality within the
            # middle. The SRA-downloaded reads contain neither the sample index
            # barcode nor the binary filtering information, and they use the
            # Phred+33 quality scores, in which B no longer indicates an anomalously
            # poor-quality read. For this reason, the associated lines of code
            # relating to these three filtering criteria are commented out in the
            # Preprocessing scripts, with "#*#"" on the right-hand-side of each
            # relevant line. The code is preserved in case original data files are
            # being used. Finally, in the directories with figures corresponding to
            # each of the two manuscripts (McGearyLinEtAl2019 and
            # McGearyBisariaEtAl2022), those figures with the suffix "_old"
            # correspond to the original figures, with data processed as originally
            # described, and those without the "_old" suffix correspond to figures
            # regenerated using the SRA-hosted files for which the filtering steps
            # were omitted.

            # Assign the multiplexing barcode sequence

            # List of conditionals to except the read:
            # 1. That the first line doesn't end with 1;0                       #**#
            read_barcode = header.split("#")[1].split("/")[0]                   #**#
            if raw_original and header[-1] == 0:                                #**#
                    counts_readtypes_map["0"] += 1                              #**#
            # 2. That "N" is not in the read sequence.
            elif "N" in seq:
                counts_readtypes_map["N"] += 1
            # 3. That "B" is not in the quality score.                          #**#
            elif raw_original and "B" in score:                                 #**#
                counts_readtypes_map["B"] += 1                                  #**#
            # 4. That the read has the correct multiplexing bardcode.           #**#
            elif raw_original and read_barcode != barcode:                      #**#
                counts_readtypes_map["M"] += 1                                  #**#
            # 5. That the sequence is not within a 1 nt mismatch to
            #   the phi-x genome, or its reverse complement.
            elif (seq in k_phi_x or seq in k_phi_x_rc):
                counts_readtypes_map["X"] += 1
            # 6. That the seqence is not a spike sequence.
            elif True in [test_read_match(seq, j, mis = 1, indel = 1)
                        for j in spikes]:
                output_spikes.append(full_seq[:40])
                counts_readtypes_map["S"] += 1
            elif seq_tag not in tag_list and "twist_reporter_assay" not in experiment:
                counts_readtypes_map["T"] += 1
            else:
                if "twist_reporter_assay_3p" in experiment:
                    output_reads.append(full_seq[:120])
                elif "twist_reporter_assay" in experiment:
                    output_reads.append(full_seq[:100])
                else:
                    output_reads.append(full_seq[:40])
                counts_readtypes_map["P"] += 1
                output_inds += [global_index + i]
        i += 1
    summary_df = pd.DataFrame.from_dict(
        counts_readtypes_map, orient="index").reindex(
        ["0", "N", "B", "M", "X", "S", "T", "P"], copy=False)
    print("about to return output")
    sys.stdout.flush()
    return output_reads, output_spikes, summary_df, output_inds


def main():
    time_start = time.time()

    # Define all the relevant arguments
    arguments = ["miRNA","experiment","condition", "-rep",
                 "-jobs", "-raw_original_binary", "-filt_by_original_binary",
                 "-test_binary"]
    mirna, experiment, condition, rep, jobs, raw_original, filt_by_original, test = parse_arguments(arguments)

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
    print("mirna: " + mirna)
    print("experiment: " + experiment)
    print("condition: " + condition)
    print("rep: %s" %rep)
    print("nb: " + str(nb))
    print("tp: " + str(tp))
    file_args = get_exp_info(
        mirna, experiment, condition, rep=rep, nb=nb, tp=tp,
        raw_original=raw_original, filt_by_original=filt_by_original,
    )
    if raw_original or filt_by_original:
        file_path, spikes, tag, barcode, int_path = file_args
    else:
        file_path, spikes, tag, barcode = file_args
    spike_seqs = [seq_spike_map[i] for i in spikes.split(",") if i != "NA"]

    if rep:
        ext = "," + rep
    else:
        ext = ""
    
    if raw_original:
        ext = "%s_raw_original" %ext

    if filt_by_original:
        ext = "%s_filtbyoriginal" %ext

    if test:
        ext = "%s_test" %ext

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
    print(file_path)
    print(reads_path)
    print(spikes_path)
    print(summary_path)
    if raw_original:
        if test:
            int_path = "%s_test.txt" %int_path
        else:
            int_path = "%s.txt" %int_path
        print(int_path)
    elif filt_by_original:
        int_path = "%s.txt" %int_path
        print(int_path)
    # Retrieve the spike used, and the 3' tag used to ensure that the
    # reads pass quality control.
    # Take in the file using the multiprocess_file function in general.py.
    if filt_by_original:
        with open(int_path, 'r') as f:
            filt_indices = set(int(line.strip()) for line in f if line.strip())
    else:
        filt_indices = None
    args = [spike_seqs, barcode, tag, experiment, raw_original, filt_indices]
    if not jobs:
        jobs = 19

    print(file_path + "\t" + reads_path)
    print(spikes_path)
    print(summary_path)
    
    threads = multiproc_file(file_path, jobs, check_read, test, *args)

    # Collect the output from each thread, which are the list
    # of reads and the dictionary of read types.
    read_seqs = [j for i in threads for j in i[0]]
    spike_seqs = [j for i in threads for j in i[1]]
    pds = [i[2] for i in threads]
    pd_summary = merge_data_frames(pds)
    if raw_original:
        ints_used = [j for i in threads for j in i[3]]
    # Flatten the list into one list and write it to its output file.
    print(reads_path)
    with open(reads_path,"w") as file_out:
        file_out.write("".join([i + "\n" for i in read_seqs]))
    with open(spikes_path,"w") as file_out:
        file_out.write("".join([i + "\n" for i in spike_seqs]))
    if raw_original:
        with open(int_path, "w") as file_out:
            print(file_out)
            file_out.write("".join([str(i) + "\n" for i in ints_used]))
    print(pd_summary)
    print(summary_path)
    pd_summary.to_csv(summary_path, sep="\t", header=False)

    print_time_elapsed(time_start)
    print("done script")


################################################################################

if __name__ == "__main__":
    main()

