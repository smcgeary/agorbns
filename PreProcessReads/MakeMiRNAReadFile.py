################################################################################
# MakeMiRNAReadFile.py
################################################################################
import imp  # Used to import general.py

imp.load_source("general", "general/general.py")
from general import *

# Constants
k_phi_x = "".join(open(("general/phi_x.txt")).read().split("\n"))


# Function
def check_read(
    reads,
    barcode,
    experiment,
    raw_original,
    filt_indices,
    byte_convert,
    global_index,
    k_phi_x=k_phi_x,
):
    """Takes a list of reads and outputs a list of reads which
       pass the criteria specified by the arguments.

    Args:
        reads: The list of reads to process
        barcode: The barcode to match
        experiment: The experiment identifier
        raw_original: Whether to apply additional filters and track indices
        filt_indices: Set of indices to filter by (or None)
        byte_convert: Whether to convert from bytes
        global_index: Starting index for this chunk
        k_phi_x: PhiX contamination sequences to filter out

    Returns:
        A tuple containing (list of reads, list of spikes, summary dataframe, list of indices)
    """
    k_phi_x_rc = get_rc(k_phi_x)
    counts_readtypes_map = dict()
    for i in ["0", "N", "B", "M", "X", "T", "S", "P"]:
        counts_readtypes_map[i] = 0
    output_reads = []
    output_spikes = []
    output_inds = []
    i = 0
    print("entering loop")
    if byte_convert:

        def get_four_reads(read_tup):
            return (i.decode("utf-8") for i in read_tup)

    else:

        def get_four_reads(read_tup):
            return (i for i in read_tup)

    total_reads = len(reads) / 4
    for read_tup in zip(*[iter(reads)] * 4):
        header, full_seq, header2, score = get_four_reads(read_tup)
        seq = full_seq.strip()
        sys.stdout.flush()
        if i % (total_reads / 10) == 0:
            print((i / float(total_reads)))
            sys.stdout.flush()

        # If filtering by original indices
        if filt_indices is not None:
            if (global_index + i) in filt_indices:
                output_reads.append(seq)
                counts_readtypes_map["P"] += 1
            else:
                counts_readtypes_map["0"] += 1
        else:
            # Apply normal filters
            if header[-1] == 0:
                counts_readtypes_map["0"] += 1
            elif "N" in seq:
                counts_readtypes_map["N"] += 1
            elif seq in k_phi_x or seq in k_phi_x_rc:
                counts_readtypes_map["X"] += 1
            elif "B" in score:
                counts_readtypes_map["B"] += 1
            else:
                output_reads.append(seq)
                counts_readtypes_map["P"] += 1
                # Track indices when in raw_original mode
                if raw_original:
                    output_inds.append(global_index + i)
        i += 1

    summary_df = pd.DataFrame.from_dict(counts_readtypes_map, orient="index").reindex(
        ["0", "N", "B", "M", "X", "S", "T", "P"], copy=False
    )
    return output_reads, output_spikes, summary_df, output_inds


def main():
    time_start = time.time()

    # Define all the relevant arguments
    arguments = [
        "miRNA",
        "experiment",
        "condition",
        "-rep",
        "-jobs",
        "-raw_original_binary",
        "-filt_by_original_binary",
        "-test_binary",
    ]

    mirna, experiment, condition, rep, jobs, raw_original, filt_by_original, test = (
        parse_arguments(arguments)
    )

    nb_check = experiment.split("_")[-1]
    if nb_check == "nb":
        nb = True
    else:
        nb = False

    # Retrieve the path to the file - UPDATED to match MakeReadFile.py
    file_args = get_exp_info(
        mirna,
        experiment,
        condition,
        rep=rep,
        nb=nb,
        raw_original=raw_original,
        filt_by_original=filt_by_original,
    )

    # Extract int_path if provided
    if raw_original or filt_by_original:
        file_path, spikes, tag, barcode, int_path = file_args
    else:
        file_path, spikes, tag, barcode = file_args

    if spikes:
        spike_seqs = [seq_spike_map[i] for i in spikes.split(",") if i != "NA"]
    else:
        spike_seqs = []

    # Build file extension
    if rep:
        ext = "," + rep
    else:
        ext = ""

    if raw_original:
        ext = "%s_raw_original" % ext

    if filt_by_original:
        ext = "%s_filtbyoriginal" % ext

    if test:
        ext = "%s_test" % ext

    reads_path, spikes_path, summary_path = [
        get_analysis_path(mirna, experiment, condition, i, ext=ext)
        for i in ["reads", "spike_reads", "summary"]
    ]

    # Update int_path based on test flag - MATCHES MakeReadFile.py exactly
    if raw_original:
        if test:
            int_path = "%s_test.txt" % int_path
        else:
            int_path = "%s.txt" % int_path
        print(f"Indices will be written to: {int_path}")
    elif filt_by_original:
        int_path = "%s.txt" % int_path
        print(f"Reading indices from: {int_path}")

    # Load indices if filtering by original
    if filt_by_original:
        with open(int_path, "r") as f:
            try:
                filt_indices = set(int(line.strip()) for line in f if line.strip())
                print(f"Loaded {len(filt_indices)} indices for filtering")
            except ValueError as e:
                print(f"Error reading indices file: {e}")
                sys.exit(1)
    else:
        filt_indices = None

    # Set up arguments for multiprocessing
    args = [barcode, experiment, raw_original, filt_indices]

    if not jobs:
        jobs = 19

    print(f"Processing file: {file_path}")
    print(f"Output paths: {reads_path}, {spikes_path}, {summary_path}")

    # Process file
    threads = multiproc_file(file_path, int(jobs), check_read, test, *args)

    # Collect the output from each thread
    read_seqs = [j for i in threads for j in i[0]]
    spike_seqs = [j for i in threads for j in i[1]]
    pds = [i[2] for i in threads]
    pd_summary = merge_data_frames(pds)
    print(len(read_seqs), "reads processed")
    # sys.exit()
    # Write indices to file if raw_original is True
    if raw_original:
        ints_used = [j for i in threads for j in i[3]]
        print(f"Writing {len(ints_used)} indices to {int_path}")
        with open(int_path, "w") as file_out:
            file_out.write("".join([str(i) + "\n" for i in ints_used]))

    # Write reads and spikes to files
    with open(reads_path, "w") as file_out:
        file_out.write("".join([i + "\n" for i in read_seqs]))

    with open(spikes_path, "w") as file_out:
        file_out.write("".join([i + "\n" for i in spike_seqs]))

    print(pd_summary)
    pd_summary.to_csv(summary_path, sep="\t", header=False)

    print_time_elapsed(time_start)
    print("done script")


if __name__ == "__main__":
    main()
