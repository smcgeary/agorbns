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

# FUNCTIONS
def count_read_kmers(read_seqs, _mirna, _kmerlist, experiment, n_constant,
                     rand_length, min_ex_str, buffer3p):
    time_start = time.time()
    sys.stdout.flush()
    for i_r, r in enumerate(read_seqs):
        read = r.strip()
        _read = Read(read, rand_length, _mirna, n_constant, experiment)
        exclude = sum([_read.seq.find(ex) for ex in min_ex_str])
        # time_3 = time.time()
        if exclude == -1*len(min_ex_str):
            # time_4 = time.time()
            # print("passed conditional: %f" %(time_4 - time_3))
            _kmerlist += _read
        # time_5 = time.time()
        # print("Passed conditional and added kmers: %f" %(time_5 - time_3))
        # sys.stdout.flush()

    return _kmerlist

def main():
    time_start = time.time()
    arguments = ["miRNA", "experiment", "condition", "n_constant", "kmer_len",
                 "-n_ex", "-uniq_binary", "-buffer3p_binary", "-final_binary",
                 "-weird", "-cutoffdil_binary", "-alt_miRNA", "-jobs",
                 "-temp_path", "-test_binary"]
    args = parse_arguments(arguments)
    (mirna, experiment, condition,
     n_constant, kmer_len, n_ex, uniq,
     buffer3p, final, weird, cutoffdil, 
     alt_miRNA, jobs, temp_path, test  ) = args

    if experiment in ["equil_mmseed_nb", "equil_seed_nb"]:
        rand_length = 38
    else:
        rand_length = 37

    if alt_miRNA:
        _mirna = Mirna(alt_miRNA)
    else:
        _mirna = Mirna(mirna)

    if final:
        site_string = "paperfinal"
    elif weird:
        site_string = "papercutoffweirdsite%smer" %(weird)
    elif cutoffdil:
        site_string = "papercutoffdil"
    else:
        site_string = "papercutoff"


    _sitelist = SiteList(_mirna, site_string,
                         int(rand_length) + 2*int(n_constant))
    exclude_list = [i.name for i in _sitelist.sites if i]
    print(exclude_list)
    if n_ex:
        if n_ex[0] in ["A", "C", "T", "G"]:
            n_ex = n_ex.split(",")
            n_ex_strings = True
        else:
            n_ex = int(n_ex)
            n_ex_strings = False
    else:
        n_ex = len(exclude_list)
        n_ex_strings = True

    ############################################################################
    if (condition == "I_combined" or 
        (condition == "I" and mirna == "miR-7-24nt"
         and experiment == "equilibrium_tp")):
        args_master = args
        # This updates the argument list.
        for i, argument in enumerate(arguments):
            if argument == "-alt_miRNA":
                args[i] = mirna
        if (condition == "I" and mirna == "miR-7-24nt"
            and experiment == "equilibrium_tp"):
            condition_list = ["I_%s" %(i + 1) for i in range(6)]
            input_list = [(mirna, experiment, i) for i in condition_list]
        elif "tp" in experiment:
            input_list = INPUT_LIST_I_COMBINED_TP
        else:
            input_list = INPUT_LIST_I_COMBINED
        # Create lists for the paths and job_ids of the bsub jobs used to
        # perform the combined kmer analysis.
        temp_paths = []
        job_ids = []
        print(input_list)
        # Iterate over the input tuples in `input_list`:   
        for input_tuple in input_list:
            args = args_master
            # Replace the `mirna`, `experiment`, and `condition` arguments:
            args[:3] = input_tuple
            # Iterate through the `arguments` list to find the index associated
            # with the `-temp_path` flag and assign a temporary path at that
            # location.
            for i, argument in enumerate(arguments):
                if argument == "-temp_path":
                    path = "%s_%s.txt" %(time.time(), randomword(20))
                    args[i] = path
                    temp_paths.append(path)
            # Make the string of arguments to be used in the bsub:
            arg_str = make_args_string(args, arguments)
            # Make the full call to the shell.
            shell_call = ("bsub -q 18 -n 20 -R span[hosts=1] python "
                          "AssignSiteTypes/AssignKmers.py %s" %arg_str)
            print(shell_call)
            # Open up the subprocess that submits this call to the shell.
            p1 = subprocess.Popen([shell_call], shell=True,
                                  stdout = subprocess.PIPE)
            # Get the output of the subprocess, convert it to a utf-8 string
            # (which is required since python 3).
            p1_output = str(p1.communicate()[0], "utf-8")
            # Cut the string to get the job id (between < and >).
            job_id = p1_output.split("\n")[0].split("<")[1].split(">")[0]
            # Add this job id to the list.
            job_ids.append(job_id)
        job_ids_current = [i for i in get_job_ids() if i != ""]
        print(job_ids_current)
        print(intersection(job_ids, job_ids_current))
        num_jobs_remaining = len(intersection(job_ids, job_ids_current))
        while num_jobs_remaining >= 1:
            print("Waiting 10 seconds")
            sys.stdout.flush()
            time.sleep(10)
            job_ids_current = [i for i in get_job_ids() if i != ""]
            jobs_current = intersection(job_ids, job_ids_current)
            print(jobs_current)
            sys.stdout.flush()
            num_jobs_remaining = len(jobs_current)

        threads = []
        for temp_path_i in temp_paths:
            _kmer_i = KmerList(int(kmer_len), int(n_constant), rand_length, fast=True)
            _kmer_i.read(mirna, experiment, condition, int(n_constant),
                         int(n_ex), temp_path=temp_path_i)
            threads.append(_kmer_i)

    ###################################################################
    else:
        _kmerlist = KmerList(int(kmer_len), int(n_constant), rand_length,
                             buffer_=buffer3p)

        # The new code by which the new sites are excluded from the analysis:###
        if uniq:
            read_dir = "reads_unique"
        else:
            read_dir = "reads"
        if not n_ex_strings: 
            exclude_strings = [i.seq for i in _sitelist.sites][:n_ex]
        else:
            exclude_strings = n_ex
        print("exclude strings:")
        print(exclude_strings)
        min_ex_str = MinimalExcludeStringList(exclude_strings)

        if not jobs:
            jobs = 20
        args  = [_mirna, _kmerlist, experiment, int(n_constant),
                 int(rand_length), min_ex_str, buffer3p]
        ########################################################################

        reads_path = get_analysis_path(mirna, experiment, condition, read_dir)
        print(reads_path)
        threads = multiproc_file(reads_path, int(jobs), count_read_kmers, test,
                                  *args)



    kmers = get_kmer_list(int(kmer_len))[:10]
    for thread in threads:
        ns = thread.kmers["TGT"][:10]
        means = thread.pos_mean["TGT"][:10]
        sds = thread.pos_sd["TGT"][:10]
        print("thread:")
        for kmer, n, m, sd in zip(kmers, ns, means, sds):
            print(("%s\t%s\t%s\t%s" %(kmer, n, m, sd)))
    # print("done reading files")
    # print("summing threads:")
    _kmers = sum([thread for thread in threads])
    # print("Done summing threads.")
    # print(_kmers.pos_mean["TGT"][:10])
    pos_max = int(rand_length) + 2*int(n_constant) - int(kmer_len)
    # print("expected mean:")
    # print(pos_max/float(2))
    # print("expected sd:")
    expected_sd = (((pos_max + 1)**2 - 1)/12)**0.5
    # print(expected_sd)
    # for kmer, n, m, sd in zip(kmers, _kmers.kmers["TGT"][:10], _kmers.pos_mean["TGT"][:10],
    #                       _kmers.pos_sd["TGT"][:10]):
    #     print("%s\t%s\t%s\t%s" %(kmer, n, m, sd))
    extension = "_%s_k%s" %(n_constant, kmer_len)
    if uniq:
        extension = "%s_uniq" %(extension)
    if buffer3p:
        extension = "%s_buffer3p" %(extension)
    print((type(n_ex)))
    if type(n_ex) == list:
        extension = "%s_ex%s" %(extension, ",".join(exclude_strings))
    elif n_ex > 0:
        # extension_old = "%s_ex:%s" %(extension,
        #                              ",".join(exclude_list[:int(n_ex)]))
        extension = "%s_ex%s" %(extension, n_ex)
    
    if final:
        extension = "%s_final" %(extension)
    elif weird:
        extension = "%s_weirdsite%smer" %(extension, weird)

    if cutoffdil:
        extension = "%s_dil" %extension

    if test:
        extension = "%s_test" %extension
    
    if temp_path:
        kmers_path = "I_combined_temp/%s" %temp_path
    else:
        kmers_path = get_analysis_path(mirna, experiment, condition,
                                       "kmers_cutoff_final", ext=extension)
    print((sum(_kmers.kmers["TGT"])))
    print((sum(_kmers.kmers["ACA"])))
    print(kmers_path)
    _kmers.write(kmers_path)
    print("finished writing kmer file.")

    print_time_elapsed(time_start)
################################################################################

if __name__ == "__main__":
    main()

