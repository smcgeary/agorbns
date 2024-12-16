################################################################################
#MakeSiteTypeReadFiles.py
################################################################################
import imp # Used to import general.py
import time # Used for time
import numpy as np
import pandas as pd
import itertools as it
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
from sitetypes import get_seq_site_map

np.set_printoptions(linewidth=200)

# FUNCTIONS
def count_kmer_freqs(read_seqs, _mirna, experiment, n_constant,
                      rand_length, kmer_length, comp_rc_kmers):
    time_start = time.time()
    # kpl is the kmer position length.
    kpl = rand_length + 2*n_constant - kmer_length + 1
    # Make output dictionary.
    kmer_dict = {
        "TGT" : {kmer : 0 for kmer in comp_rc_kmers},
        "ACA" : {kmer : 0 for kmer in comp_rc_kmers}
    }
    # Make the total reads dictionary (allows multiple 3' end barcodes.)
    count_reads = {"TGT" : 0, "ACA" : 0}
    tot_reads = len(read_seqs)
    # Loop through the list of reads.
    for i_r, r in enumerate(read_seqs):

        if i_r % (tot_reads/10) == 0:
            print(float(i_r)/float(tot_reads))
            sys.stdout.flush()
        read = r.strip()
        _read = Read(read, rand_length, _mirna, n_constant, experiment)
        # print(_read.seq)
        # sys.stdout.flush()
        kmers = [_read.seq[i:i+kmer_length] for i in range(kpl)]
        for kmer in kmers:
            if kmer in comp_rc_kmers:
                # if kmer == "TCTTCCTCCGCACC":
                #     str_ind = _read.seq.find(kmer)
                #     print(" "*str_ind + kmer)
                #     print(_read.seq)
                #     sys.stdout.flush()
                kmer_dict[_read.barcode][kmer] += 1
        # Increment the pulse or chase dicationary value up by 1.
        count_reads[_read.barcode] += 1
    # Convert dictionaries to lists for output.
    counts_pulse = [kmer_dict["TGT"][kmer] for kmer in comp_rc_kmers]
    counts_chase = [kmer_dict["ACA"][kmer] for kmer in comp_rc_kmers]
    num_reads = [count_reads["TGT"], count_reads["ACA"]]
    time_done = time.time()
    print("Finished job future in %f3.2 seconds." %(time_done - time_start))
    sys.stdout.flush()
    return (counts_pulse, counts_chase, num_reads)

def main():
    # Get the time for the starting position.
    time_start = time.time()
    # Define and parse the arguments.
    arguments = ["miRNA", "experiment", "condition", "n_constant",
                 "kmer_length", "seed_offset", "-addcomp", "-wob_binary",
                 "-alt_mirna_comp", "-jobs", "-test_binary"]
    (mirna, experiment, condition,
     n_constant, kmer_length,
     seed_offset, addcomp, wob,
     alt_mirna_comp, jobs, test) = parse_arguments(arguments)
    # Assign length of random sequence.
    if (
        (experiment in ["equilibrium_mmseed_nb",
                        "equilibrium_seed_nb",
                        "equilibrium_mmseed_2_nb"]) and
        (mirna in ["miR-1", "let-7a-21nt", "let-7a"])
    ):
        rand_length = 38
    else:
        rand_length = 37    
    # Assign Mirna object.
    _mirna = Mirna(mirna)
    if not jobs:
        jobs = 20
    # Assign kmers complementary to the competitor oligo.
    if alt_mirna_comp:
        mirna_comp = alt_mirna_comp
    elif "miR-7" in mirna:
        if "tp" in experiment:
            mirna_comp = "miR-7.tp"
        elif "nb" in experiment:
            mirna_comp = "miR-7.nb"
    else:
        mirna_comp = mirna
    comp_rc_seq = get_rc(seq_competitor_map[mirna_comp][13+int(seed_offset):])
    print(comp_rc_seq)

    comp_rc_kmers = [comp_rc_seq[i:i+int(kmer_length)]
                     for i in range(len(comp_rc_seq) - int(kmer_length) + 1)]
    if wob:
        comp_rc_kmers += get_wobble_kmers(comp_rc_kmers)

    # MULTIPROCESSING PART #####################################################
    # Define the path for the analysis, and the arguments to be given to the
    # function.

    reads_path = get_analysis_path(mirna, experiment, condition, "reads")
    args  = [_mirna, experiment, int(n_constant),int(rand_length),
             int(kmer_length), comp_rc_kmers]
    # Use the multiproc_file function from general with the path.
    threads = multiproc_file(reads_path, int(jobs), count_kmer_freqs, test,
                             *args)
    # kpl = kmer position length
    kpl = int(rand_length) + 2*int(n_constant) - int(kmer_length) + 1
    total_kmers = [sum([thread[2][i] for thread in threads])*kpl
                   for i in range(2)]
    print(total_kmers)
    # n_kmers = float(n_reads*kpl)

    counts_pulse = [sum([thread[0][i] for thread in threads])
                    for i in range(len(comp_rc_kmers))]
    counts_chase = [sum([thread[1][i] for thread in threads])
                    for i in range(len(comp_rc_kmers))]
    print("before adding comp")
    print(counts_pulse)
    print(counts_chase)
    print(total_kmers)
    if addcomp:
        n_reads = sum(total_kmers)
        print(n_reads)
        comp_reads_14 = int(float(n_reads)*float(addcomp)*0.3/100)
        comp_reads_13 = int(float(n_reads)*float(addcomp)*1/100)
        comp_reads_12 = int(float(n_reads)*float(addcomp)*4/100)
        comp_reads_11 = int(float(n_reads)*float(addcomp)*8/100)
        comp_reads_10 = int(float(n_reads)*float(addcomp)*4/100)
        comp_reads_9 = int(float(n_reads)*float(addcomp)*4/100)
        comp_reads_8 = int(float(n_reads)*float(addcomp)*4/100)
        comp_reads_7 = int(float(n_reads)*float(addcomp)*4/100)
        comp_reads_6 = int(float(n_reads)*float(addcomp)*40/100)



        print(comp_reads_14)
        print("proportion:")
        print(float(comp_reads_14)/n_reads)
        print("Simulating reads")
        additional_reads_14 = simulate_competitor_reads(comp_reads_14, mirna,
                                                     "equilibrium", 14)
        additional_reads_13 = simulate_competitor_reads(comp_reads_13, mirna,
                                                        "equilibrium", 13)
        additional_reads_12 = simulate_competitor_reads(comp_reads_12, mirna,
                                                        "equilibrium", 12)
        additional_reads_11 = simulate_competitor_reads(comp_reads_11, mirna,
                                                "equilibrium", 11)
        additional_reads_10 = simulate_competitor_reads(comp_reads_10, mirna,
                                                "equilibrium", 10)
        additional_reads_9 = simulate_competitor_reads(comp_reads_9, mirna,
                                                "equilibrium", 9)
        additional_reads_8 = simulate_competitor_reads(comp_reads_8, mirna,
                                                "equilibrium", 8)
        additional_reads_7 = simulate_competitor_reads(comp_reads_7, mirna,
                                                "equilibrium", 7)

        additional_reads_6 = simulate_competitor_reads(comp_reads_6, mirna,
                                                "equilibrium", 6)


        additional_reads = (additional_reads_14 + additional_reads_13 +
                            additional_reads_12 + additional_reads_11 +
                            additional_reads_10 + additional_reads_9 +
                            additional_reads_8 + additional_reads_8 +
                            additional_reads_7 + additional_reads_6)

        thread = count_kmer_freqs(additional_reads, _mirna, experiment,
                                  int(n_constant), int(rand_length),
                                  int(kmer_length), comp_rc_kmers)
        counts_pulse = [counts_pulse[i] + thread[0][i]
                        for i in range(len(comp_rc_kmers))]
        counts_chase = [counts_chase[i] + thread[1][i]
                        for i in range(len(comp_rc_kmers))]

        total_kmers = [total_kmers[i] + thread[2][i]*kpl for i in range(2)]
        print("after adding comp")
        print(counts_pulse)
        print(counts_chase)
        print(total_kmers)


    df_index = comp_rc_kmers + ["Total"]

    kmers = pd.DataFrame(counts_pulse + [total_kmers[0]], index=df_index)
    kmers.columns = [""]

    if experiment in ["kinetics", "kinetics_pilot"]:
        extension = "_pulse_%s_k%s_off%s" %(n_constant, kmer_length, seed_offset)
    else:
        extension = "_%s_k%s_off%s" %(n_constant, kmer_length, seed_offset)
    if addcomp:
        extension = extension + "_addcomp%s" %addcomp
    if wob:
        extension = extension + "_wob"
    if alt_mirna_comp:
        extension = extension + "_altcomp_%s" %alt_mirna_comp
    if test:
        extension = extension + "_test"
    kmers_path = get_analysis_path(mirna, experiment, condition,
                                   "competitor_kmers", ext=extension)

    print(kmers_path)

    kmers.to_csv(kmers_path, sep="\t", header=False)

    if experiment in ["kinetics", "kinetics_pilot"]:
        kmers_chase = pd.DataFrame(counts_chase + [total_kmers[1]], index=df_index)
        extension_chase = "_chase_%s_k%s" %(n_constant, kmer_length)
        if addcomp:
            extension = extension + "_addcomp%s" %(addcomp)
        if wob:
            extension = extension + "_wob"
        if test:
            extension_chase = extension_chase + "_test"
        kmers_path_chase = get_analysis_path(mirna, experiment, condition,
                                             "competitor_kmers",
                                             ext=extension_chase)
        print(kmers_path_chase)
        kmers_chase.to_csv(kmers_path_chase, sep="\t", header=False)

    print_time_elapsed(time_start)
################################################################################

if __name__ == "__main__":
    main()

