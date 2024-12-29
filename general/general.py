################################################################################
#general.py
################################################################################
import argparse
import csv
import gzip
import itertools as it
import multiprocessing
import os
import pandas as pd
import random
import re
import string
import sys
import tarfile
import time
from subprocess import Popen, PIPE
from concurrent.futures import ProcessPoolExecutor

# HOME_DIR = "/lab/bartel1_ata/mcgeary/computation/AgoRBNS/"
DATA_DIR = "/data/processed/"


# Map of the sequences associated with the quantitative spikes put into some of
# the libraries (that are actually used in the analyses).
seq_spike_map = {
    "s124_8mer" :         "TCTCATCTACCTCCCGGTTTTAATGAATAGTGCCTTAGAC",
    "s430a_8mer" :        "TCTCATCTACCTCCCGGTTTTAATGAATAAGCACTTAGAC",
    "s427_8mer" :         "TCTCATCTACCTCCCGGTTTTAATGAATAGCACTTTCGAC",
    "s14_8mer" :          "TCTCATCTACCTCCCGGTTTTAATGAATATCGCTCCCGAC",
    "s155_8mer" :         "TCTCATCTACCTCCCGGTTTTAATGAATAAGCATTAAGAC",
    "mir1_libconstruct1" :    "ATCTACCTCCCGGTTTTAATGAATAACATTCCAGACTCGT",
    "mir1_libconstruct2" :    "ATCTACCTCCCGGTTTTAATGAATACTACCTCAGACTCGT"
}


# Various conditions among the experiments performed:
conditions_py = ["I", "I_combined", "40", "12.6", "4", "1.26", "0.4", "0"]
# conditions_py = ["I", "I", "40", "12.6", "4", "1.26", "0.4", "0"]

# conditions_py_tp = ["I", "I_combined", "40", "12.6", "4", "1.26", "0.4",
#                     "0.126", "0"]

conditions_py_nocombI = ["I", "I", "40", "12.6", "4", "1.26", "0.4", "0"]

# conditions_kinetics_py = ["I", "0", "2", "5", "8", "15", "30", "90", "300",
#                           "900","2700", "7200", "Equil", "Equil-", "0-"]

conditions_R = ["I", "I_combined", "A40", "A12.65", "A4", "A1.265", "A0.4",
                "A0"]
# conditions_R = ["I", "I", "A40", "A12.65", "A4", "A1.265", "A0.4",
#                 "A0"]
conditions_R_tp = ["I", "I_combined", "A40", "A12.65", "A4", "A1.265", "A0.4",
                   "A0.126", "A0"]

conditions_R_nocombI = ["I", "I", "A40", "A12.65", "A4", "A1.265", "A0.4", "A0"]


# conditions_kinetics_R = ["I", "T0", "T2", "T5", "T8", "T15", "T30", "T90",
#                          "T300", "T900", "T2700", "T7200", "Equil", "Equil-",
#                          "T0-"]

# dups_let_7 = ["2", "T2", "5", "T5",  "8", "T8", "15", "T15",  "7200", "T7200"]
# dups_rest  = ["2", "T2", "5", "T5",  "8", "T8", "15", "T15", "30", "T30"]


# def dup_function(cond, mirna):
#     if mirna == "let-7a":
#         dups = dups_let_7
#     else:
#         dups = dups_rest
#     if cond in dups:
#         return ["%s,1" %(cond), "%s,2" %(cond)]
#     else:
#         return [cond]

CONDITIONS_PY = {
    # BEGIN EQUILIBRIUM PAPER LIBRARIES
    "equilibrium3_nb" : {
        "miR-7-24nt" : conditions_py
    },
    "equilibrium2_nb" : {
        "miR-7-23nt" : [i for i in conditions_py if i != "4"],
        "miR-7-24nt" : conditions_py,
        "miR-7-25nt" : conditions_py
    },
    "equilibrium_nb" : {
        "let-7a-21nt" : conditions_py_nocombI,    
        "miR-7-22nt"  : conditions_py,
        "miR-7-23nt"  : conditions_py,
        "miR-7-24nt"  : conditions_py,
        "miR-7-25nt"  : conditions_py
    },
    "equilibrium_tp" : {
        "miR-1"      : conditions_py,
        # "miR-122"    : conditions_py_tp,
        "miR-124"    : ["I", "I_combined", "40", "12.6", "4,1", "4,2", "1.26,1",
                        "1.26,2", "0.4", "0.126", "0.04", "0"],
        "miR-7-24nt" : ["I", "I_combined", "40", "12.6", "4", "1.26", "0.4",
                        "0"]
    },
    "equilibrium_2_tp" : {
        "miR-124" : ["I", "I_combined", "40", "12.6", "4", "1.26", "0.4",
                     "0.126", "0"]
    },
    # "equilibrium_met_tp" : {
    #     "miR-1"   : conditions_py_tp
    # },    
    "equilibrium" : {
        "miR-1"   : conditions_py,
        "let-7a"  : conditions_py,
        "miR-155" : conditions_py,
        "miR-124" : conditions_py,
        "lsy-6"   : conditions_py
    },
    # "equil_flowthrough" : {
    #     "let-7a"  : ["I", "40_2h", "40_4h", "40_2h_nc1", "40_2h_nc2",
    #                  "40_2h_ft", "0_2h", "0_4h", "0_2h_nc1", "0_2h_nc2",
    #                  "0_2h_ft"],
    #     "miR-124" : ["I", "40_2h_ft", "0_2h", "0_4h", "0_2h_nc1", "0_2h_nc2",
    #                  "0_2h_ft"]
    # },
    "equil_pilot" : {
        "miR-1"   : ["I", "L100A10"]
        # "miR-155" : ["I", "L100A10", "L10A10", "L100A0", "L10A0"]
    },
    # END EQUILIBRIUM PAPER LIBRARIES
    # BEGIN THREEPRIME PAPER LIBRARIES
    "equil_c_nb" : {
        "let-7a-21nt"    : conditions_py_nocombI,                               # 1
        "miR-1"          : conditions_py_nocombI,                               # 7
        "let-7a_plus1"   : conditions_py_nocombI,                               # 9
        "let-7a_minus1"  : conditions_py_nocombI,                               # 10
        "let-7a_miR-155" : conditions_py_nocombI,
        "miR-155_let-7a" : conditions_py_nocombI
    },
    "equil_s_nb" : {                                                            # 2
        "let-7a-21nt" : ["I", "I", "12.6", "12.6_2", "4", "4_2", "1.26", "0.4",
                         "0"]
    },
    "equil_c2_alt_nb" : {                                                       # 3
        "let-7a-21nt"  : ["I", "I", "40", "12.6", "12.6_2", "4", "1.26",
                          "0.4", "0"]
    },
    "equil_c_alt_nb" : {                                                        # 4
        "miR-1"  : conditions_py_nocombI
    },
    "equil_sc_alt_nb" : {                                                       # 5
        "miR-1"  : conditions_py_nocombI
    },
    "equil_c2_nb" : {                                                           # 6
        "let-7a-21nt"  : ["I", "I", "40", "12.6", "12.6_2", "4", "1.26",
                          "0.4", "0"]
    },
    "equil_sc_nb" : {
        "let-7a-21nt" : conditions_py_nocombI,                                  # 11
        "miR-1"       : conditions_py_nocombI,                                  # 8
        "miR-155"     : conditions_py_nocombI                                   # 12
    },
    # END THREEPRIME PAPER LIBRARIES
}

CONDITIONS_R = {
    # BEGIN EQUILIBRIUM PAPER LIBRARIES
    "equilibrium3_nb" : {
        "miR-7-24nt" : conditions_R
    },
    "equilibrium2_nb" : {
        "miR-7-23nt" : [i for i in conditions_R if i != "A4"],
        "miR-7-24nt" : conditions_R,
        "miR-7-25nt" : conditions_R
    },
    "equilibrium_nb" : {
        "let-7a-21nt" : conditions_R,
        "miR-7-22nt"  : conditions_R,
        "miR-7-23nt"  : conditions_R,
        "miR-7-24nt"  : conditions_R,
        "miR-7-25nt"  : conditions_R
    },
    "equilibrium_2_tp" : {
        "miR-124" : ["I_combined", "I", "A40", "A12.65", "A4", "A1.265", "A0.4",
                     "A0.1265", "A0"]
    },
    "equilibrium_tp" : {
        "miR-1"      : conditions_R_tp,
        "miR-122"    : conditions_R_tp,
        "miR-124"    : ["I", "I_combined", "A40", "A12.65", "A4,1", "A4,2",
                        "A1.265,1", "A1.265,2", "A0.4", "A0.1265", "A0.04",
                        "A0"],
        "miR-7-24nt" : ["I", "I_combined", "A40", "A12.65", "A4", "A1.265",
                        "A0.4", "A0"]
    },
    # # "equilibrium_met_tp" : {
    # #     "miR-1"   : conditions_R_tp
    # },    
    "equilibrium" : {
        "miR-1"   : conditions_R,
        "let-7a"  : conditions_R,
        "miR-155" : conditions_R,
        "miR-124" : conditions_R,
        "lsy-6"   : conditions_R
    },
    # "equil_flowthrough" : {
    #     "let-7a"  : ["I", "A40_2h", "A40_4h", "A40_2h_nc1", "A40_2h_nc2",
    #                  "A40_2h_ft", "A0_2h", "A0_4h", "A0_2h_nc1", "A0_2h_nc2",
    #                  "A0_2h_ft"],
    #     "miR-124" : ["I", "A40_2h_ft", "A0_2h", "A0_4h", "A0_2h_nc1",
    #                  "A0_2h_nc2", "A0_2h_ft"]
    # },
    "equil_pilot" : {
        "miR-1"   : ["I", "L100A10"],
        # "miR-155" : ["I", "L100A10", "L10A10", "L100A0", "L10A0"]
    },
    # END EQUILIBRIUM PAPER LIBRARIES
    # THESE LIBRARIES MATTER FOR THE THREEPRIME PAPER.
    "equil_c_nb" : {
        "let-7a-21nt"    : conditions_R,
        "miR-1"          : conditions_R,
        "let-7a_plus1"   : conditions_R,
        "let-7a_minus1"  : conditions_R,
        "let-7a_miR-155" : conditions_R,
        "miR-155_let-7a" : conditions_R
    },
    "equil_s_nb" : {
        "let-7a-21nt" : ["I", "I_combined", "A12.65", "A12.65_2", "A4", "A4_2",
                         "A1.265", "A0.4", "A0"]
    },
    "equil_c2_alt_nb" : {
        "let-7a-21nt"  : ["I", "I_combined", "A40", "A12.65", "A12.65_2", "A4",
                          "A1.265", "A0.4", "A0"]
    },
    "equil_c_alt_nb" : {
        "miR-1" : conditions_R
    },
    "equil_sc_alt_nb" : {
        "miR-1" : conditions_R
    },
    "equil_c2_nb" : {
        "let-7a-21nt"  : ["I", "I_combined", "A40", "A12.65", "A12.65_2", "A4",
                          "A1.265", "A0.4", "A0"]
    },
    "equil_sc_nb" : {
        "let-7a-21nt" : conditions_R,
        "miR-1"       : conditions_R,
        "miR-155"     : conditions_R
    },
    # END THREEPRIME PAPER LIBRARIES
}

# # INPUT_LIST_I_COMBINED = [("let-7a",  "equilibrium",     "I"),
# #                          ("miR-124", "equilibrium",     "I"),
# #                          ( "let-7a",    "kinetics",     "I"),
# #                          (  "miR-1",    "kinetics",     "I"),
# #                          ("miR-124",    "kinetics",     "I"),
# #                          (  "lsy-6",    "kinetics",     "I"),
# #                          (  "miR-1",   "kin_pilot", "I_TGT"),
# #                          (  "miR-1",   "kin_pilot", "I_ACA")]


INPUT_LIST_I_COMBINED = [( "let-7a", "equilibrium",     "I"),
                         ("miR-124", "equilibrium",     "I"),
                         (  "miR-1",   "kin_pilot", "I_TGT")]

INPUT_LIST_I_COMBINED_TP = [(  "miR-1",     "equilibrium_tp",   "I"),
                            (  "miR-1", "equilibrium_met_tp",   "I"),
                            ("miR-122",     "equilibrium_tp",   "I"),
                            ("miR-124",     "equilibrium_tp",   "I"),
                            ("miR-124",   "equilibrium_2_tp", "I,1"),
                            ("miR-124",   "equilibrium_2_tp", "I,2")]


INPUT_LIST_I_COMBINED_NB = [("miR-7-23nt", "equilibrium2_nb",   "I"),
                            ("miR-7-24nt", "equilibrium3_nb",   "I")]


# INPUT_LIST_0_COMBINED = [("let-7a",  "equilibrium",     "0"),
#                          ("miR-1",   "equilibrium",     "0"),
#                          ("miR-155", "equilibrium",     "0"),
#                          ("miR-124", "equilibrium",     "0"),
#                          ("let-7a",     "kinetics",     "0-"),
#                          ("miR-1",      "kinetics",     "0-"),
#                          ("miR-124",    "kinetics",     "0-"),
#                          ("lsy-6",      "kinetics",     "0-"),
#                          ("let-7a",     "kinetics",     "Equil-"),
#                          ("miR-1",      "kinetics",     "Equil-"),
#                          ("miR-124",    "kinetics",     "Equil-"),
#                          ("lsy-6",      "kinetics",     "Equil-"),
#                          ("miR-1",     "kin_pilot", "I_TGT"),
#                          ("miR-1",     "kin_pilot", "I_ACA")]





def print_time_elapsed(time_start):
    time_sec = time.time() - time_start
    time_h = time_sec // 3600
    time_m = (time_sec % 3600) // 60
    time_s = time_sec % 60
    print(("%.0f:%.0f:%.2f" % (time_h, time_m, time_s)))


def get_rc(seq, rna=False):
    """
    Get the reverse complement of a DNA/RNA sequence.

    Takes a DNA or RNA sequence and returns its reverse complement. Can handle both
    DNA and RNA inputs through the rna parameter.

    Args:
        seq (str): Input sequence containing A,C,G,T/U,N nucleotides
        rna (bool, optional): If True, uses U instead of T in complement. Defaults to False.
    
    Returns:
        str: Reverse complement of the input sequence. Will contain U if rna=True, T if rna=False.
    
    Examples:
        >>> get_rc("ACGT")
        'ACGT'
        >>> get_rc("ACGU", rna=True)
        'ACGU'
    """
    nucleotide_complement_map = {
    "C" : "G",
    "G" : "C",
    "T" : "A",
    "U" : "A",
    "N" : "N",  
    }
    # Assigns A to U or T depending on 'rna' flag.
    if rna:
        nucleotide_complement_map["A"] = "U"
    else:
        nucleotide_complement_map["A"] = "T"
    # Constructs a list of the characters substituted and reversed in order.
    rev_complement_list = [nucleotide_complement_map[i] for i in seq][::-1]
    rev_complement = "".join(rev_complement_list)
    return(rev_complement)



def parse_arguments(arguments):
    """
    Uses argparse to parse command line arguments and gives a list with
    each corresponding variable.

    Args:
        arguments: A list of strings representing the arguments needed in 
        the command line.

    Returns:
        A list of values corresponding to every argument, including
        flags arguments. Flags for which no argument was entered return a
        None object.
    """
    parser = argparse.ArgumentParser()
    for argument in arguments:
        if "_binary" in argument:
            parser.add_argument(argument.split("_binary")[0],
                                action='store_true')
        else:
            parser.add_argument(argument)
    args = vars(parser.parse_args())
    arguments_no_dash = []
    for i in arguments:
        if i[0] == "-":
            arguments_no_dash.append(i.split("_binary")[0][1:])
        else:
            arguments_no_dash.append(i)
    return [args[i] for i in arguments_no_dash]


def get_exp_info(mirna, experiment, condition, rep=None, nb=None, tp=None):
    """
    Retrieves the complete path to the zipped fastq file 
    in solexa_bartel/mcgeary.

    Args:
        exp: The folder in which the data is, representing the experiment.
        barcode: The multiplexing barcode associated with the sample.
        lane: The lane the sample sequenced on within the sequencing run.

    Returns:
        The string representing the full path to the zipped fastq file.
    """
    if nb:
        csv_file = "Libraries_NB.csv"
    elif tp:
        csv_file = "Libraries_TP.csv"
    else:
        csv_file = "Libraries_SL.csv"
    csv_file = "Libraries.csv"
    # Load the database of experiments:
    with open(csv_file, 'U') as master:
        library_list = csv.DictReader(master)
        # Determine the experiment to be analyzed
        exp_file = None
        for row in library_list:
            ## THIS CAN BE UN-COMMENTED AS A DIAGNOSTIC OF WHY THE PREPROCESSING
            ## ISN'T WORKING.
            # if row["SL"] == "15":
            #     print(row)
            #     print((mirna in row["miRNA"].split(",")))
            #     print(experiment)
            #     print((row["Exp_type"]))
            #     print((experiment == row["Exp_type"]))
            #     print((row["Sample type"]))
            #     print(condition)
            #     print((row["Sample type"] == condition))
            #     print((row["Rep"]))
            if (mirna in row["miRNA"].split(",") and
                    row["Exp_type"] == experiment and
                    row["Sample type"] == condition and
                    row["Exp"] != "CANCELLED" and
                    (rep == None or row["Rep"] == rep)):
                exp_file = row
    # Retrieve the path to the file, from the data in the appropriate
    # row of the database.
    exp, barcode, lane = [exp_file[i] for i in ["Exp", "Barcode", "Lane"]]
    SRR_number, spikes, tag = [exp_file[i] for i in ["SRR number", "Spike", "3prime"]]

    date = int(exp.split("_")[0])
    if date <= 180500:
        tar_string = ".tar"
    else:
        tar_string = ""

    path = "/lab/solexa_public/Bartel/%s/QualityScore/" % (exp)
    path = "data/raw/fastq/"
    print(path)
    file = "%s-s_%s_1_sequence.txt%s.gz" % (barcode, lane, tar_string)
    file = "%s.fq.gz" % (SRR_number)
    full_path = path + file
    return full_path, spikes, tag, barcode


def ensure_directory(directory):
    """Checks for a directory and creates it if it does not exist.

        Args:
            directory: The directory in question.
            string: The string being searched for the 'key' in. Note that this
                string must be longer than the 'key' argument.
            mismatch: the number of tolerated mismatches in the match.

        Returns:
            Does not return output, but a function call generates the directory
            if it doesn't exist.
    """
    if not os.path.exists(directory):
        os.makedirs(directory)

def get_analysis_path(mirna, experiment, condition, analysis_type, ext = "",
                      suffix = "txt"):
    """Generates the intended path and name reads extracted from the fastq
        file.

        Args:
            mirna: The miRNA in the experiment to be processed.
            experiment: The type of experiment to be processed.
            condition: The sample type within the experiment to be processed.
            analysis_type: The type of analysis file to be read/written.
            ext: Any extension information required.
            suffix: The file type.
        Returns:
            The string representing the full path to the output file.
    """
    path_outer = "data/processed/"
    path_inner = "%s/%s/%s" % (mirna, experiment, analysis_type)
    directory = path_outer + path_inner
    print(directory)
    # sys.stdout.flush()
    ensure_directory(directory)

    full_path = "%s/%s%s.%s" % (directory, condition, ext, suffix)
    return full_path


def multiproc_file(path, n_jobs, func, test, *args, **kwargs):
    if "halve_reads" in kwargs and kwargs["halve_reads"]:
      halve_reads = True
    else:
      halve_reads = False
    print("number of cores:")
    print((multiprocessing.cpu_count()))
    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        # Initialize the list of futures (the actual processes), and the list
        # of the iterable (called the iterator) that we be repeatedly cleared
        # when the list is greater than the intended process size.
        futures = []

        if ("_sequence.txt.tar.gz" in path or
            "_sequence.txt.gz" in path or
            ".fastq.gz" in path or
            ".fq.gz" in path):
            print("in here")
            mult = 4
            wc_proc = Popen("zcat %s | wc -l" %(path), shell=True, stdout=PIPE)
        else:
            mult = 1
            print("here")
            print("path exists using python command:")
            print(os.path.exists(path))
            wc_proc = Popen("wc -l %s" %(path), shell=True, stdout=PIPE)
            # print(wc_proc.communicate()[0])
        print("path:")
        print(path)
        if test:
            process_size = 300000
            n_jobs = 2
        else:
            file_lines = int(wc_proc.communicate()[0].split()[0])
            print(("file_lines: %s" %(file_lines)))
            if halve_reads:
                if (file_lines/mult) % 2 == 0:
                    line_mod_best = 0
                else:
                    line_mod_best = 1
                line_mod = ((file_lines/mult) % n_jobs)
                job_mod = n_jobs % 2
                while line_mod != line_mod_best or job_mod != 0:
                    n_jobs -= 1
                    line_mod = ((file_lines/mult) % n_jobs)
                    job_mod = n_jobs % 2
            process_size = file_lines//mult//n_jobs
        # Conditional for starting new job.
        i_jobs = 0
        remaining_jobs = n_jobs
        tick = 0
        print(n_jobs)
        if path[-7:] == ".tar.gz":
            print("This is a tar file")
            tar_in = tarfile.open(path, "r:gz")
            member = tar_in.getmembers()[0]
            file_in = tar_in.extractfile(member)
        elif path[-3:] == ".gz":
            file_in = gzip.open(path, "rt")
        else:
            file_in = open(path, "r")
        time_start = time.time()
        job_num = 0
        # Conditional for those scripts that require the index of the job
        # within the script itself, so pass an additional argument.
        if func.__name__ == "check_read":
            if "tar.gz" in path:
                byte_convert = True
            else:
                byte_convert = False
            args = [i for i in args] + [byte_convert]
        if func.__name__ in ["get_read_structural_data", "get_plfold_matrix",
                             "check_GCTTCCG"]:
            args = [i for i in args] + [i_jobs]
        # Boolean conditional used for preprocessing, because the `.txt.tar.gz`
        # files need to be converted from bytes to characters, while the
        # `.txt.gz` files are alreay read as characters.

        while n_jobs > 1:
            print(process_size*mult)
            print(process_size)
            print(mult)
            job_reads = [file_in.readline() for i in range(process_size*mult)]
            i_jobs += 1
            if func.__name__ in ["get_read_structural_data", "get_plfold_matrix",
                                 "check_GCTTCCG"]:
                args[-1] += 1
            sys.stdout.flush()
            futures.append(executor.submit(func, job_reads, *args))
            n_jobs -= 1
        if not test:
            line = file_in.readline()
            final_job = []
            while line:
                final_job.append(line)
                line = file_in.readline()
            if func.__name__ in ["get_read_structural_data", "get_plfold_matrix",
                                 "check_GCTTCCG"]:
                args[-1] += 1
            print(len(final_job))
            futures.append(executor.submit(func, final_job, *args))
        print("through multiprocess")
        time.sleep(10)
        output = [i.result() for i in futures]
        print("made output")
    if path[-7:] == ".tar.gz":
        tar_in.close()
        # print("closed tar")
    else:
        file_in.close()
    # Return a list of the results from each job, to be parsed uniquely
    # depending on the repeated_function.
    print(len(output))
    return output


def test_read_match(key, string, mis=0, indel=0):
    """
    Determines if the key matches the string tolerating a certain number of
    matches

    Args:
        key: The string to be searched for.
        string: The string being searched for the 'key' in. Note that this
            string must be longer than the 'key' argument.
        mis: the number of tolerated mismatches in the match.

    Returns:
        A list of values corresponding to every argument, including
        flags arguments. Flags for which no argument was entered return a
        'None' object.
    """
    # Iterate along possible windows of key within string
        # Set mismatch counter
    mm = [0]*(1+ 2 * (indel))
    # Iterate along the string within this window 
    for j in range(0,len(key)):
        if string[j] != key[j]:
            if mm[indel] > 0:

                if string[j] != key[j-1]:
                    mm[indel-1] += 1
                if string[j] != key[j-1]:
                    mm[indel+1] +=1
            mm[indel] += 1
        # If there are more mismatches than the limit,
        # break loop to pass to the next window.
        if True not in [m < mis for m in mm]:
            return(False)
    print(key)
    print(string)
    print(mis)
    print(indel)
    return True


def merge_data_frames(dfs):
    df_out = dfs[0]
    if len(dfs) > 1:
        for df_i in dfs[1:]:
            df_out = df_out.add(df_i, fill_value=0)
    return(df_out)


def get_kmer_list(k, rna=False):
    if not rna:
        NTS = DNTS
    if k > 0:
        return ["".join(i) for i in list(it.product(NTS,repeat=int(k)))]
    else:
        return []




def main():
    return

if __name__ == "__main__":
    main()







# import argparse
# import concurrent.futures
# import collections
# import csv
# import itertools as it
# import math
# import matplotlib.pyplot as plt
# import multiprocessing
# import numpy as np
# import os
# import pandas as pd
# import random
# import re
# import subprocess
# import sys
# import string
# import tarfile
# from importlib import reload
# from subprocess import Popen, PIPE

# pd.set_option('display.width', 1000)
# pd.set_option('display.precision', 3)

# from concurrent.futures import ProcessPoolExecutor






# seq_mirna_map = {
#     "miR-1"       : "UGGAAUGUAAAGAAGUAUGUAU",
#     "miR-155"     : "UUAAUGCUAAUCGUGAUAGGGGU",
#     "miR-124"     : "UAAGGCACGCGGUGAAUGCCAA",
#     "let-7a-21nt" : "UGAGGUAGUAGGUUGUAUAGU",
#     "let-7a"      : "UGAGGUAGUAGGUUGUAUAGUU",
#     "lsy-6"       : "UUUUGUAUGAGACGCAUUUCGA",
#     "miR-7"       : "UGGAAGACUAGUGAUUUUGUUGU",
#     "miR-7-22nt"  : "UGGAAGACUAGUGAUUUUGUUG",
#     "miR-7-23nt"  : "UGGAAGACUAGUGAUUUUGUUGU",
#     "miR-7-24nt"  : "UGGAAGACUAGUGAUUUUGUUGUU",
#     "miR-7-25nt"  : "UGGAAGACUAGUGAUUUUGUUGUUU",
#     "miR-671"     : "AGGAAGCCCUGGAGGGGCUGG",
#     "miR-21"      : "UAGCUUAUCAGACUGAUGUUGA"
# }

# seq_capture_map = {
# 	"miR-1"       : "UCUUCCUCCGCACCACACACAUUCCAACCUUACACAC",
#     "let-7a-21nt" : "UCUUCCUCCGCACCACACCUACCUCAACCUUACACAC",
#     "let-7a"      : "UCUUCCUGCGCACCAAGCCUACCUCAACUUUACACAC",
#     "miR-155"     : "AAAUAAAGACGACAACUCAGCAUUAAACCUUACACAC",
#     "miR-124"     : "UCUUCCUCCGCCACAGAAGUGCCUUAACCUUACACAC",
#     "lsy-6"       : "UCUUCCUCCGCCACAGAAAUACAAAAACCUUACACAC",
#     "miR-7.nb"    : "UCUUCCUCCGCACCACACGUCUUCCAACCUUACACAC",
#     "miR-7.tp  "  : "ACAUCGUCCGCACCACACGUCUUCCAACCUUACACAC"
# }

# seq_competitor_map = {
# 	"miR-1"       : "AAGGTTGGAATGTGTGTGGTGCGGAGGAAGA",
# 	"let-7a-21nt" : "AAGGTTGAGGTAGGTGTGGTGCGGAGGAAGA",
# 	"let-7a"      : "AAAGTTGAGGTAGGCTTGGTGCGCAGGAAGA",
# 	"miR-155"     : "AAGGTTTAATGCTGAGTTGTCGTCTTTATTT",
# 	"miR-124"     : "AAGGTTAAGGCACTTCTGTGGCGGAGGAAGA",
# 	"lsy-6"       : "AAGGTTTTTGTATTTCTGTGGCGGAGGAAGA",
# 	"miR-7.nb"    : "AAGGTTGGAAGACGTGTGGTGCGGAGGAAGA",
# 	"miR-7.tp"    : "AAGGTTGGAAGACGTGTGGTGCGGACGATGT",
# }

# seq_lib_map = {
#     "5p" : {
#         "miR-1_pilot"   : "GGGUUCAGAGUUCUACAGUCCGACGAUC",
#         "miR-155_pilot" : "GGGUUCAGAGUUCUACAGUCCGACGAUC",
#         "miR-1"         : "GGGCAGAGUUCUACAGUCCGACGAUC", 
#         "let-7a"        : "GGGCAGAGUUCUACAGUCCGACGAUC",
#         "miR-155"       : "GGGCAGAGUUCUACAGUCCGACGAUC",
#         "miR-124"       : "GGGCAGAGUUCUACAGUCCGACGAUC",
#         "lsy-6"         : "GGGCAGAGUUCUACAGUCCGACGAUC",
#         "miR-7"         : "GGGCAGAGUUCUACAGUCCGACGAUC", 
#     },
#     "3p" : {
#         "miR-1_pilot"   : "UAUGCCGUCUUCUGCUUG",
#         "miR-155_pilot" : "UAUGCCGUCUUCUGCUUG",
#         "miR-1"         : "UAUGCCGUCUUCUGCUUG", 
#         "let-7a"        : "UGUUCGUAUGCCGUCUUCUGCUUG",
#         "miR-155"       : "UGUUCGUAUGCCGUCUUCUGCUUG",
#         "miR-124"       : "UGUUCGUAUGCCGUCUUCUGCUUG",
#         "lsy-6"         : "UGUUCGUAUGCCGUCUUCUGCUUG",
#         "miR-7.nb1"     : "UGUUCGUAUGCCGUCUUCUGCUUG", 
#         "miR-7.nb2"     : "UGUUCGUAUGCCGCUGUGUGCUUG",
#         "miR-7.tp"      : "UGUUCGUAUGCCGCUGUGUGCUUG",
#     }
# }

# read3p_mirna_map = {
#     "miR-1"      : "TCG",
#     "let-7a"     : "TGT",
#     "miR-155"    : "TGT",
#     "miR-124"    : "TGT",
#     "lsy-6"      : "TGT",
#     "miR-7-23nt" : "TGT"
# }


# marker_18nt = "AGCGUGUAGGGAUCCAAA"
# weird   =   "GGAGCGTGTAGGATCCAAATCGTATGCC"
# marker_18nt_d = "AGCGTGTAGGGATCCAAA"
# marker_30nt = "GGCAUUAACGCGGCCGCUCUACAAUAGUGA"
#                # GGTCGTTTCCCGGCCCATGCACCATCGT
# marker_30nt_d = "GGCATTAACGCGGCCGCTCTACAATAGTGA"

# adapter_5p = "GUUCAGAGUUCUACAGUCCGACGAUC"
# adapter_5p_d = "GTTCAGAGTTCTACAGTCCGACGATC"
# adapter_3p_sRS = "TCGTATGCCGTCTTCTGCTTG"
# dme_miR_14_5p = "GGGAGCGAGACGGGGACUCACU"
# dme_miR_14_5p_d = "GGGAGCGAGACGGGGACTCACT"
# xtr_miR_427 = "GAAAGUGCUUUCUGUUUUGGGCG"
# xtr_miR_427_d = "GAAAGTGCTTTCTGTTTTGGGCG"

# strange_sequence = "AGCGTGTAGGGATCCAAA"

# temp="GGCATTAACGCGGCCG"


# # [marker_18nt_d, marker_30nt_d, adapter_5p_d, dme_miR_14_5p_d, xtr_miR_427_d] = [
# #     re.sub("U", "T", i) for i in [marker_18nt, marker_30nt, adapter_5p,
# #                                   dme_miR_14_5p, xtr_miR_427]
# # ]

# constant_seqs = [marker_18nt_d,
#                  marker_30nt_d,
#                  adapter_5p_d,
#                  adapter_3p_sRS,
#                  dme_miR_14_5p_d,
#                  xtr_miR_427_d]



# weird_rev =                                              "3p-GCTTCGACCCCGTCTCCTTCTCCAAGCACCACTAGTAGGT-p5"
# adapter_3p =                                                                                 "TCGTATGCCGTCTTCTGCTTG"
# pilot_library =             "GGGUUCAGAGUUCUACAGUCCGACGAUCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNUAUGCCGUCUUCUGCUUG"
# miR_1_library =               "GGGCAGAGUUCUACAGUCCGACGAUCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNUCGUAUGCCGUCUUCUGCUUG"
# pulselibrary =             "GGGCAGAGUUCUACAGUCCGACGAUCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNUGUUCGUAUGCCGUCUUCUGCUUG"
# rev_155 =                               "UGGGGAUAGUGCUAAUCGUAAUU"
# num_155 =                                              "87654321"

# chaselibrary =             "GGGCAGAGUUCUACAGUCCGACGAUCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNACAUCGUAUGCCGUCUUCUGCUUG"
# newlibrary =               "GGGCAGAGUUCUACAGUCCGACGAUCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNUGUUCGUAUGCCGCUGUGUGCUUG"

# weird =                   "5p-TGGATGATCACCACGAACCTCTTCCTCTGCCCCAGCTTCG-p3"
# rev5pconstat=                    "NCUAG CAGCCUGACAUCUUGAGACGGG"
# weird_rev = "3p-GCTTCGACCCCGTCTCCTTCTCCAAGCACCACTAGTAGGT-p5"

# NTS = ["A", "C", "G", "U"]
DNTS = ["A", "C", "G", "T"]

# WOBBLE_NTS = ["A", "C", "G", "U", "H", "V"]

# RNA_comp_map = {
#     "C" : "G",
#     "G" : "C",
#     "A" : "U",
#     "U" : "A",
#     "V" : "H",
#     "H" : "V",  
# }



# def get_job_ids():
#     return(Popen("bjobs | cut -d ' ' -f 1", shell=True,
#           stdout=PIPE, stderr=PIPE).communicate()[0].decode().split("\n"))[1:-1]

# def get_job_hosts():
#     return(Popen("bjobs | cut -c 46- | cut -d ' ' -f 1", shell=True,
#           stdout=PIPE, stderr=PIPE).communicate()[0].decode().split("\n"))[1:-1]


# def get_job_statuses():
#     return(Popen("bjobs | cut -c 17- | cut -d ' ' -f 1", shell=True,
#           stdout=PIPE, stderr=PIPE).communicate()[0].split("\n"))[1:-1]

# def count_lines(path):
#     return(int(Popen("wc -l %s" %(path), shell=True,
#           stdout=PIPE).communicate()[0].split()[0]))

# def ls(path):
#     return(Popen("ls %s" %(path), shell=True,
#           stdout=PIPE, stderr=PIPE).communicate()[0].strip().split())

# def logit(vector, maximum):
#     return(-np.log(maximum/vector - 1))

# def make_args_string(args, arguments):
#     # Returns the command line string to use after calling a particular sript
#     # given:
#     # arg = This is the list of argument values give by the user.
#     # arguments = This is the string of argument labels used by argparse.
#     # 1. First make a list of tuples associating each argument label with its
#     # argument value.
#     arg_tuples = [(i[1], i[0]) for i in zip(args, arguments) if i[0]]
#     # 2. Enumerate over the list of argument tuples making two changes:
#     for i, arg_i in enumerate(arg_tuples):
#         arg_name, arg = arg_i
#         # 2.1. If the argument name doesn't begin with "-", remove it from the
#         # list, because it is a positional argument, rather than flag, and
#         # therefor doesn't need a label in the command line calling of the
#         # script.
#         if arg_name[0] != "-":
#             arg_name = ""
#         # If the argument name does have a "-", it is a flag, and any "_binary"
#         # suffix should be removed because that suffix should only exist within
#         # the python script.
#         elif "_binary" in arg_name:
#             arg_name = "".join(arg_name.split("_binary"))
#             arg = ""
#         arg_tuples[i] = (arg_name, arg)
#     # 3. Remove the nested tuple structure of the argument list, removing
#     # any blank entries (should only be blank) due to the mising name.
#     arg_tuples = [str(j) for i in arg_tuples for j in i if j != ""]
#     arg_str = " ".join(arg_tuples)
#     return(arg_str)



# ## FUNCTIONS
# def logit(vector, maximum):
#     return(-np.log(maximum/vector - 1))


# def pd_dict(dict, sortindex=True, sortcolumn=True):
#     df = pd.DataFrame.from_dict(dict, orient="index")
#     if sortindex:
#         df.reindex(sorted(df.index), copy=False)
#     if sortcolumn:
#         df.sort_index(axis=1, inplace=True)
#     return(df)


def randomword(length):
    return ''.join(random.choice(string.ascii_lowercase) for i in range(length))

# def randomseq(length, rna=False):
#     if not rna:
#         NTS = DNTS
#     return ''.join(random.choice(NTS) for i in range(length))

# def geo_mean(x_vector):
#     return np.exp(np.mean(np.log(x_vector)))


# def get_site_flanks(read, site_pos_5p, site_pos_3p, num_flanks=2):
#     flanks_5p = read[site_pos_5p - num_flanks:site_pos_5p]
#     flanks_3p = read[site_pos_3p:site_pos_3p + num_flanks]
#     return "".join([flanks_5p, flanks_3p])


# def CorrectReadSequence(read, mirna, read_length, experiment):
#     barcode = read[26 + read_length : 26 + read_length + 3]
#     # Deals with the TCG instead of TGT ending for miR-1 equilibrium exp
#     if barcode == "TCG":
#         barcode = "TGT"
#         if mirna != "miR-1" or experiment != "equilibrium":
#             read = read[:26 + 37] + "TGTTCGTATGCCGTCTTCTGCTTG"
#     if barcode == "TGT" and mirna == "miR-1" and experiment == "equilibrium":
#         read = read[:26 + 37] + "TCGTATGCCGTCTTCTGCTTG"

# def LCSuffOld(seq1, seq2):
#     if len(seq1) == 0 or len(seq2) == 0:
#         return("")
#     elif seq1[-1] == seq2[-1]:
#         return(LCSuffOld(seq1[:-1], seq2[:-1]) + seq1[-1])
#     else:
#         return("")

# def LCSubStringOld(seq1, seq2):
#     pres1 = [seq1[:i] for i in range(1, len(seq1) + 1)]
#     pres2 = [seq2[:i] for i in range(1, len(seq2) + 1)]
#     pre_pairs = [(i, j) for i in pres1 for j in pres2]
#     LCs = [LCSuffOld(i[0], i[1]) for i in pre_pairs]
#     LCmax = max([len(i) for i in LCs])
#     return [i for i in LCs if len(i) == LCmax]

# def MinimalExcludeStringList(s_list):
#     if len(s_list) > 0:
#         s_list.sort(key = lambda s: len(s))
#         min_len = len(s_list[0])
#         max_len = len(s_list[-1])
#         removes = [0]*len(s_list)
#         for s_1 in s_list:
#             for i, s_2 in enumerate(s_list):
#                 if s_2.find(s_1) != -1 and s_1 != s_2:
#                     removes[i] +=1
#         s_list = [s for (i, s) in enumerate(s_list) if removes[i] == 0]
#     return(s_list)

# def LCSuff(seq1, seq2):
#     if len(seq1) > 0 and len(seq2) > 0 and seq1[-1] == seq2[-1]:
#         return(LCSuff(seq1[:-1], seq2[:-1]) + 1)
#     else:
#         return(0)

# def LCSubString(seq1, seq2, alt_len=None):
#     # Make all prefixes of both strings
#     if len(seq1) == 0 or len(seq2) == 0:
#         return([0, [[], []]])
#     else:
#         pre_1 = [seq1[:i] for i in range(1, len(seq1) + 1)]
#         pre_2 = [seq2[:i] for i in range(1, len(seq2) + 1)]
#         # All paired combinations between the prefixes:
#         pre_pairs = [(i, j) for i in pre_1 for j in pre_2]
#         # The Longest common suffixes of the pairs
#         LCs = [LCSuff(i[0], i[1]) for i in pre_pairs]
#         if len(LCs) == 0:
#             print(seq1)
#             print(seq2)
#         if alt_len:
#             LCmax = alt_len
#         else:
#             LCmax = max(LCs)
#         if LCmax == 0:
#             return([0, [[], []]])
#         else:
#             inds_max = [ind for (ind, lc) in enumerate(LCs) if lc == LCmax]
#             # Find the starting positions of the longest substrings in the first
#             # and second sequence, respectively:

#             inds_1 = [i//len(pre_2) - LCmax + 1 for i in inds_max]
#             inds_2 = [i%len(pre_2) - LCmax + 1 for i in inds_max]
#             return([LCmax, [inds_1, inds_2]])


# def intersection(lst1, lst2): 
#     return list(set(lst1) & set(lst2))


# temp_text = ["CATTCCCTAGGTCAAATCCTCCATTCCACATTCCCACTCG",
#              "ATGCAAATCCCTTTGACATTCCCTTCCGCATCCTAAATCG",
#              "CAACATTCCATGATCTGCGTGAAAATGTTCAATAGGTTCG",
#              "CAATTGACATCTACATACCATATTTTAACTATCCTATTCG",
#              "CTAAGACATGGACATACAACGATGATCAACATCCGTGTCG",
#              "ACAATTACCCACATTCCGAGTCATTCCTCCAACCAGCTCG",
#              "CATTCCATTCTCCCAAAACTCTCCAAACATAGGTGGGTCG",
#              "TGAGCCATTCCGACATCCGTCTGACATTCCGCTCTACTCG",
#              "TTATGTGATACTTAGTCTTGTTAACGTGAAATGCGACTCG",
#              "TGTTGAAGGTCCTTGTCGACGCTAAATTCTTTGGGTATCG",]

# # def find_all(seq, key, start=0):
# #     match = seq.find(key, start)
# #     if match >= 0:
# #         return([match] + find_all(seq, key, match + 1))
# #     else:
# #         return([])

def find_all(seq, key, start=0, iteration=0):
    match = re.search(key, seq[start:])
    if iteration > 3:
        return([])
    if match:
        # Update the position for the substring:
        pos = match.span()[0] + start
        pos_output = pos
        key_update = key
        # Updates position 1 with respect to "[^G]:"
        while re.match("\[\^.*\]", key_update):
            # Shift position 1 to the left if [^N] at start:
            pos_output += 1
            # Update key_update to remove "[^N]" from beginning of string:
            key_update = re.sub("^\[\^.*\]", "", key_update)
        return([pos_output] + find_all(seq, key, pos + 1))
    else:
        return([])


# def calculate_local_au(utr, site_start, site_end, site_type):
#     """
#     Calculate the local AU score

#     Parameters
#     ----------
#     utr: string, utr sequence

#     site_type: string, site type

#     site_start: int, start of site

#     site_end: int, end of site

#     Output
#     ------
#     float: local AU score
#     """
#     # find A, U and weights upstream of site
#     up_site_adder = int(site_type not in ['7mer-m8', '8mer'])
    
#     upstream = utr[max(0, site_start - 30): site_start]
#     upstream_str = upstream

#     l= max(0, site_start - 30)
#     # print(site_type)
#     # print(utr[0:100])
#     upstream = [int(x in ['A', 'T']) for x in upstream]
#     inv_upweights = [(x + 1 + up_site_adder)
#                  for x in range(len(upstream))][::-1]

#     upweights = [1.0 / (x + 1 + up_site_adder)
#                  for x in range(len(upstream))][::-1]

#     # find A,U and weights downstream of site
#     down_site_adder = int(site_type in ['7mer-A1', '8mer'])
#     downstream = utr[site_end:min(len(utr), site_end + 30)]
#     downstream_str = downstream
#     downstream = downstream + "A"*(30 - len(downstream))
#     downstream = [int(x in ['A', 'T']) for x in downstream]

#     inv_downweights = [(x + 1 + down_site_adder)
#                    for x in range(len(downstream))]
#     downweights = [1.0 / (x + 1 + down_site_adder)
#                    for x in range(len(downstream))]
#     # print("_"*l + upstream_str + "."*(site_end - site_start) + downstream_str)
#     # print(" "*l + "".join([str(i) for i in upstream]) + "."*(site_end - site_start) + "".join([str(i) for i in downstream]))
#     # for num, weight in enumerate(inv_upweights):
#     #     print(" "*(l + num) + str(weight) + " "*(30 + site_end - site_start - len(str(weight))) + str(inv_downweights[num]))

#     weighted = np.dot(upstream, upweights) + np.dot(downstream, downweights)
#     total = float(sum(upweights) + sum(downweights))

#     return weighted / total


# def CheckRandomLibrary(lib_seq):
#     lib_seq_rc = get_rc(lib_seq)
#     print(lib_seq_rc)
#     for key, seq in list(seq_mirna_map.items()):
#         print("________________")
#         print(key)
#         seq_rc = get_rc(seq, rna = True)
#         seq_rev = seq[::-1]
#         matches = LCSubString(seq_rc, lib_seq)
#         if len(matches) > 0:
#             print(matches)
#             for match in matches:
#                 lib_seq_temp = lib_seq
#                 end_pos = lib_seq_temp.find(match)
#                 while end_pos != -1:
#                     print((" "*10 + lib_seq))
#                     mir_pos = seq.find(get_rc(match, rna=True))
#                     print((" "*(10 + end_pos - len(str(mir_pos + len(match))))
#                           + str(mir_pos + len(match)) + "."*len(match) + str(mir_pos + 1)))
#                     # print(seq)
#                     # print(end_pos)
#                     # print(len(seq))
#                     # print(mir_pos)
#                     # print(len(match))
#                     mirna_adjust_pos = 10 + end_pos - (len(seq) - mir_pos - len(match))
#                     # print(mirna_adjust_pos)
#                     print((" "*mirna_adjust_pos + seq_rev))
#                     print("\n")
#                     lib_seq_temp = lib_seq_temp[end_pos+1:]
#                     end_pos = lib_seq_temp.find(match)
#                     mir_pos = seq.find(get_rc(match))


#     return



# def get_mfe_heteroduplex(rna_seq, dna_seq):
#     # This should give the delta G from the variation in RNA fold, for use with
#     # the competitor oligo.

#     # Make filename that is unambiguously unique:
#     temp_filename = "%s_%s.fa" %(time.time(), randomword(20))
#     with open(temp_filename, "wb") as f:
#         f.write(">RNA\n%s\n>DNA\n%s" %(rna_seq, dna_seq))
#     shell_call = "/nfs/apps/ViennaRNA/bin/RNAduplex < %s" %(temp_filename)
#     p1 = subprocess.Popen([shell_call], shell=True, stdout = subprocess.PIPE)
#     # Parse the output 
#     output = p1.communicate()[0].split("\n")[2].split(" ")[-1]
#     deltaG = output[1:-1]
#     os.remove(temp_filename)
#     return(float(deltaG))



# def get_competitor_oligo_mfe(mirna, len_k):
#     # Get competitor oligo sequence
#     comp_oligo = seq_competitor_map[mirna]
#     comp_seqs = [comp_oligo[i:i+len_k]
#                  for i in range(len(comp_oligo) - len_k + 1)]
#     rna_seqs = [get_rc(i, rna=True) for i in comp_seqs]

#     print(comp_seqs)
#     print(rna_seqs)
#     dGs = [get_mfe_heteroduplex(rna_seq, dna_seq)
#                 for (rna_seq, dna_seq) in zip(rna_seqs, comp_seqs)]
#     print(dGs)
#     out_path = "CompetitorOligoMFEs/%s_k%s.txt" %(mirna, len_k)
#     with open(out_path, "wb") as file:
#         file.write("\n".join(["%s\t%s" %(comp_seq, dG)
#                               for (comp_seq, dG) in zip(comp_seqs, dGs)]))
#     return


# def get_competitor_oligo_mfe_with_lib(mirna, len_k, lib="5p"):
#     # Get competitor oligo sequence
#     comp_oligo = seq_competitor_map[mirna]
#     comp_seqs = [comp_oligo[i:i+len_k]
#                  for i in range(len(comp_oligo) - len_k + 1)]

#     lib_seq = seq_lib_map[lib][mirna]
#     print(lib_seq)
#     if lib == "5p":
#         rna_seqs = [lib_seq + get_rc(i, rna=True) for i in comp_seqs]
#     else:
#         rna_seqs = [get_rc(i, rna=True) + lib_seq for i in comp_seqs]

#     print(comp_seqs)
#     print(rna_seqs)
#     dGs = [get_mfe_heteroduplex(rna_seq, comp_oligo)
#                 for rna_seq in rna_seqs]
#     print(dGs)
#     out_path = "CompetitorOligoMFEs/%s_k%s_lib%s.txt" %(mirna, len_k, lib)
#     with open(out_path, "wb") as file:
#         file.write("\n".join(["%s\t%s" %(comp_seq, dG)
#                               for (comp_seq, dG) in zip(comp_seqs, dGs)]))
#     return




def get_mfe(seq1, seq2, wobble=True):
    '''Returns the ∆G in kcal/mol for the two sequences. This function calls the
        viennarna library to use RNAduplex fold. Which itself calls the nearest-
        neighbor model. The default parameters assume that the characters
        correspond to RNA, irrespective of whether the T or U is used. The
        default setting of this function is that wobble pairing is considered
        (rather than considered a mismatch).

    Args:
        seq1: The first RNA sequence.
        seq2: The second RNA sequence.
        wobble: Boolean for whether or not wobbles are considered.

    Returns:
        The ∆G of pairing, in kcal/mol.


    '''
    import RNA
    if wobble:
        RNA.cvar.noGU = 0
        RNA.cvar.temperature = 37
    else:
        RNA.cvar.noGU = 1
        RNA.cvar.temperature = 37
    duplex = RNA.duplexfold(seq1,seq2)
    out = duplex.energy
    return(out)




def get_site_seq(mirna, start, stop, wobble=False,
                           mismatch=False, bulge=False, rna=False):
    """Returns the sequence (in DNA form) based on the complementarity of
        the site.

    Args:
        mirna: The miRNA.
        start: The first nucleotide contained in the site (1 = first position).
        stop: The last nucleotide contained in the site.
        wobble: If not False, a number specifying the identify of the wobble
            paired site.

    Returns:
        A string giving the sequence corresponding to the requisite site.
    """
    # Convert miRNA sequence to list (which can be updated)
    rc_list = [i for i in seq_mirna_map[mirna]]
    print(rc_list)
    # Make the sequence consistent with the A1 rule.
    rc_list[0] = "T"
    # Update for wobbles.
    if wobble:
        wobble_nucleotide_map = {"G" : "A", "U" : "C"}
        if type(wobble) == list:
            for wob in wobble:
                rc_list[wob - 1] = wobble_nucleotide_map[rc_list[wob - 1]]
        else:
            wobble_nucleotide_map = {"G" : "A", "U" : "C"}
            rc_list[wobble - 1] = wobble_nucleotide_map[rc_list[wobble - 1]]
    # Update for mismatches.
    if mismatch:
        rc_list[mismatch[0] - 1] = get_rc(mismatch[1])
    # Define the end before bulge introduction.
    rc_list = rc_list[ : stop]
    # Define bulge.
    if bulge:
        rc_list.insert(bulge[0] - 1, get_rc(bulge[1]))
    # Define starting position.
    rc_list = rc_list[start - 1 : ]
    print(rc_list)
    # Get and return the reverce complement of the merged string.
    site_seq = get_rc("".join(rc_list), rna=rna)
    return site_seq

# def get_miRNA_mfes(mirna, start, stop, subseq=True, alt_start=False, alt_stop=False):
#     print(mirna)
#     print(start)
#     print(stop)
#     t_seq = get_site_seq(mirna, start, stop, rna=True)
#     print(t_seq)
#     print("hi")
#     m_seq = seq_mirna_map[mirna]
#     mstart = 0
#     mstop = len(m_seq)
#     if subseq:
#         mstart = start
#         mstop = stop
#     if alt_start:
#         mstart = alt_start
#     if alt_stop:
#         mstop = alt_stop
#     m_seq = m_seq[(mstart - 1):mstop]
#     print(m_seq)
#     return(get_mfe(m_seq, t_seq))
    

# def GetAllSeedSites(mirna, start, stop, mm=[]):
#     mir_seq = seq_mirna_map[mirna]
#     mir_mm_r = list(range(stop, 8))
#     mir_mm_l = list(range(0, start - 1))
#     t_seq_base = get_rc(mir_seq[start - 1:stop], rna=True)
#     if len(mir_mm_r) > 0:
#         t_seq_l = [n for n in nts if n != get_rc(mir_seq[mir_mm_r[0]], rna=True)]
#         if len(mir_mm_r) > 1:
#             t_seq_ll = get_kmer_list(len(mir_mm_r)-1, rna=True)
#             t_seq_l = ["".join(i) for i in list(it.product(t_seq_ll, t_seq_l))]
#     else:
#         t_seq_l = [""]
#     if len(mir_mm_l) > 0:
#         t_seq_r = [n for n in nts if n != get_rc(mir_seq[mir_mm_l[::-1][0]], rna=True)]
#         if len(mir_mm_l) > 1:
#             t_seq_rr = get_kmer_list(len(mir_mm_l)-1, rna=True)
#             t_seq_r = ["".join(i) for i in list(it.product(t_seq_r, t_seq_rr))]
#     else:
#         t_seq_r = [""]
#     return(["".join(i) for i in list(it.product(t_seq_l, [t_seq_base], t_seq_r))])













# def get_analysis_path_burge(rbp, experiment, condition, analysis_type, ext = "",
#                       suffix = "txt"):
#     """Generates the intended path and name reads extracted from the fastq
#         file.

#         Args:
#             mirna: The miRNA in the experiment to be processed.
#             experiment: The type of experiment to be processed.
#             condition: The sample type within the experiment to be processed.
#             analysis_type: The type of analysis file to be read/written.
#         Returns:
#             The string representing the full path to the output file.
#     """
#     path_outer = "/lab/solexa_bartel/mcgeary/RBNS_DominguezFreese/"
#     path_inner = "%s/%s/%s" % (rbp, experiment, analysis_type)
#     directory = path_outer + path_inner
#     sys.stdout.flush()
#     ensure_directory(directory)

#     full_path = "%s/%s%s.%s" % (directory, condition, ext, suffix)
#     return full_path



# def get_func_arg_names(func):
#     return func.__code__.co_varnames[:func.__code__.co_argcount]





# # TODO fix this
# def readline_one(file_in):
#     return file_in.readline()

# def readline_two(file_in_pair):
#     out1 = readline_one(file_in_pair[0])
#     out2 = readline_one(file_in_pair[1])
#     if out1 and out2:
#         return [out1, out2]
#     else:
#         return False

# def readline_three(file_in_triplet):
#     out1 = readline_one(file_in_triplet[0])
#     out2 = readline_one(file_in_triplet[1])
#     out3 = readline_one(file_in_triplet[2])
#     if out1 and out2 and out3:
#         return [out1, out2, out3]
#     else:
#         return False

# def readline_four(file_in_quartet):
#     out1 = readline_one(file_in_quartet[0])
#     out2 = readline_one(file_in_quartet[1])
#     out3 = readline_one(file_in_quartet[2])
#     out4 = readline_one(file_in_quartet[3])
#     if out1 and out2 and out3 and out4:
#         return [out1, out2, out3, out4]
#     else:
#         return False

# def test_match(key, string, mis=0):
#     """Determines if the key matches the string tolerating a certain number of
#         matches

#         Args:
#             key: The string to be searched for.
#             string: The string being searched for the 'key' in. Note that this
#                 string must be longer than the 'key' argument.
#             mis: the number of tolerated mismatches in the match.

#         Returns:
#             Either returns the first index within the string where the key
#             matches it, or 'False' if it never matches.
#     """
#     # Iterate along possible windows of key within string
#     for i in range(0, len(string)-len(key)+1):
#         # Set mismatch counter
#         mm = 0
#         # Iterate along the string within this window 
#         for j in range(0,len(key)):
#             if string[i+j] != key[j]:
#                 mm += 1
#             # If there are more mismatches than the limit,
#             # break loop to pass to the next window.
#             if mm > mis:
#                 break
#             # Check if it made it to the end, and if so,
#             # return the match index.
#             if j == len(key)-1:
#                 return i
#     return False








# def longest_match(key, string, len_min, mis=0):
#     """Determines if the key matches the string tolerating a certain number of
#         matches

#         Args:
#             key: The string to be searched for.
#             string: The string being searched for the 'key' in. Note that this
#                 string must be longer than the 'key' argument.
#             len_min: The minimum length of complementarity, under which the
#                 function returns [False, False, False].
#             mis: the number of tolerated mismatches in the match.

#         Returns:
#             Returns a list, where the first item is the max length, and the
#             second and third items are the starting positions of complementarity
#             in the key and string, respectively. 
#     """
#     # Iterate along possible windows of key within string
#     len_best = False
#     key_l_best = False
#     str_l_best = False
#     for len_use in range(len_min, len(key) + 1):
#         sub_keys = [key[i:i + len_use] for i in range(len(key) - len_use + 1)]
#         str_l = [test_match(sub_key, string, mis=mis) for sub_key in sub_keys]
#         if sum([type(i) == bool for i in str_l]) != len(str_l):
#             len_best = len_use
#             key_l_best = [i for i in range(len(sub_keys)) if type(str_l[i]) != bool] 
#             str_l_best = [i for i in str_l if type(i) != bool]
#         else:
#             break
#     return([len_best, key_l_best, str_l_best])




# def get_kmer_list_constant_insert(mirna, kmer_length, mir_start, mir_stop):
#     # Define the length of the full kmer.

#     # Get the constant sequence that is identical to a region of the miRNA.
#     kmer_const = get_site_seq(mirna, mir_start, mir_stop)

#     # Get the list of random kmers.
#     kmer_length_rand = kmer_length - len(kmer_const)
#     kmer_rand_list = get_kmer_list(kmer_length_rand)

#     # position to insert the constant sequence within the random kmer:
#     insert_index = int(8 - mir_stop + math.ceil((kmer_length - 8)/2.0))
#     print(insert_index)
#     # Make the list of full kmers, which have inserted constant sequence.
#     kmer_list = [i[ : insert_index] + kmer_const +
#                  i[insert_index: ] for i in kmer_rand_list]

#     return(kmer_list)

# def make_trie(kmer_length, seq_length):
#     """Recursive make a full trie based on the given kmer length.
#     Each node has a child called "stop" that stores position counts.
#     """
#     if kmer_length == 0:
#         return {'stop': np.zeros(seq_length)}
#     else:
#         new_dict = {'stop': np.zeros(seq_length)}
#         for nt in ['A','C','G','T']:
#             new_dict[nt] = make_trie(kmer_length - 1, seq_length) 
#         return  new_dict

# def add_to_trie(seq, trie, pos):
#     """Loops through characters of a sequence and recursively update trie counts."""
#     trie['stop'][pos] += 1
#     if len(seq) > 0:
#         char = seq.pop(0)
#         add_to_trie(seq, trie[char], pos)

# def count_trie(trie, current_kmer):
#     """Traverses trie and returns a dictionary with all the substring and position counts"""
#     kmer_dict = {current_kmer: trie['stop']}
#     if len(trie) > 1:
#         for nt in ['A','C','G','T']:
#             kmer_dict.update(count_trie(trie[nt], current_kmer + nt))

#     return kmer_dict

# def count_kmers(seqs, kmer_length, seq_length):
#     """Given a list of sequences:
#     1. Make a full trie
#     2. Add each sequence to the trie
#     3. Aggregate counts into a dictionary and return dict
#     """
#     mytrie = make_trie(kmer_length, seq_length)

#     for seq in seqs:
#         for pos in range(len(seq)):
#             subseq = list(seq[pos:min(len(seq),pos+kmer_length)])
#             add_to_trie(subseq, mytrie, pos)

#     return count_trie(mytrie, '')


# def count_positional_kmers(seqs, len_kmer):
#     """Given a list of sequences:
#     1. Make a full trie
#     2. Add each sequence to the trie
#     3. Aggregate counts into a dictionary and return dict
#     """
#     len_seq = len(seqs[0])
#     num_pos = len_seq - len_kmer + 1
#     kmer_dict = {kmer : [0]*num_pos for kmer in get_kmer_list(2)}
#     for seq in seqs:
#         for pos in range(num_pos):
#             kmer_dict[seq[pos:pos+len_kmer]][pos] += 1
#     kmer_df = pd.DataFrame.from_dict(kmer_dict).transpose()
#     kmer_df.index = get_kmer_list(2)
#     return kmer_df



# def sequential(kmers, seqs, kmer_length, seq_length):
#     """Naive method that just loops through all the sequences and updates counts
#     """
#     master_kmer_dict = {kmer: np.zeros(seq_length) for kmer in kmers}
#     for seq in seqs:
#         for i in range(seq_length):
#             for j in range(i+1, min(seq_length, i+kmer_length) + 1):
#                 subseq = seq[i:j]
#                 if subseq in kmers:
#                     master_kmer_dict[subseq][i] += 1

#     return master_kmer_dict


# def GetThreePrimeBarcodeOverlap(read, adapt=adapter_3p):
#     window = len(adapt)
#     if window == 0:
#         return(None)
#     else:
#         end = read[-window:]
#         # print(end)
#         # print(" "*(len(read) - window) + adapt)
#         # print(" "*(len(read) - window) + end)
#         position = read.find(adapt)
#         # print(position)
#         if adapt != end:
#             return(GetThreePrimeBarcodeOverlap(read, adapt=adapt[:-1]))
#         else:
#             return(len(read) - window)


# def simulate_competitor_read(n, mirna, exp, len_comp):
#     """Generates reads with regions of complementarity to a particular miRNA's
#         competitor oligo, that are otherwise constructed using the dinucleotide
#         frequencies of the input-sequencing from the RBNS experiments.
#         matches

#         Args:
#             n: The number of reads to be generated.
#             mirna: The miRNA name.
#             exp: The particular RBNS experiment.
#             len_comp: The length of complementarity to the competitor oligo.

#         Returns:
#             A list of read sequences 40-nt in length, which can be used in
#             combination with the actual reads from the AGO-RBNS experiments.
#     """
#     time_start = time.time()
#     # Define constants:
#     len_rand = 37
#     dinucs = get_kmer_list(2)
#     read3p = read3p_mirna_map[mirna]
#     # Get the dinucleotide frequencies for that miRNA's input sequence:
#     dinuc_path = get_analysis_path(mirna, exp, "I", "positional_kmers",
#                                    ext="_0_k2")
#     df_dinucs = pd.read_csv(dinuc_path, sep="\t", header=None)
#     df_dinucs.index = dinucs
#     df_dinucs = df_dinucs/df_dinucs.sum()
#     # Get the competitor oligo sequence, and all possible kmers that are
#     # reverse complementary to the competitor oligo:
#     comp_rc = get_rc(seq_competitor_map[mirna])
#     comp_kmers = [comp_rc[i:i + len_comp]
#                   for i in range(len(comp_rc) - len_comp + 1)]
#     # Pre-allocate the output list:
#     reads = [None]*n
#     for i in range(n):
#         if i % (n/10) == 0:
#             print((float(i)/n))
#         # Draw random position and random kmer of complementarity to the
#         # competitor oligo.
#         pos_kmer = random.choice(list(range(len_rand - len_comp + 1)))
#         kmer = random.choice(comp_kmers)

#         # Build read starting with this random kmer at this random position.
#         read = kmer
#         # Extend the sequence 5-prime and 3-prime of the kmer to the full read
#         # length using the positional dinucleotide probabilities.
#         # Five prime:
#         for j in range(0, pos_kmer)[::-1]:
#             probs = df_dinucs.filter(regex="%s$" %(read[0]), axis=0)[j]
#             probs = probs/probs.sum()
#             dinuc = np.random.choice(probs.index.tolist(), p=probs)
#             read = dinuc[0] + read
#         # Three prime:
#         for j in range(pos_kmer + len_comp - 1, df_dinucs.shape[1]):
#             probs = df_dinucs.filter(regex="^%s" %(read[-1]), axis=0)[j]
#             probs = probs/probs.sum()
#             dinuc = np.random.choice(probs.index.tolist(), p=probs)
#             read = read + dinuc[1]
#         # Add 3-prime-most constant sequence in the read, and update the list of
#         # reads.
#         read = read + read3p
#         reads[i] = read 
#     print_time_elapsed(time_start)
#     return(reads)


# def simulate_competitor_reads(n, mirna, exp, len_comp):
#     """Generates reads with regions of complementarity to a particular miRNA's
#         competitor oligo, that are otherwise constructed using the dinucleotide
#         frequencies of the input-sequencing from the RBNS experiments.
#         matches

#         Args:
#             n: The number of reads to be generated.
#             mirna: The miRNA name.
#             exp: The particular RBNS experiment.
#             len_comp: The length of complementarity to the competitor oligo.

#         Returns:
#             A list of read sequences 40-nt in length, which can be used in
#             combination with the actual reads from the AGO-RBNS experiments.
#     """
#     time_start = time.time()
#     # Define constants and call external functions:
#     len_rand = 37
#     nucs = get_kmer_list(1)
#     dinucs = get_kmer_list(2)
#     read3p = read3p_mirna_map[mirna]
#     comp_rc = get_rc(seq_competitor_map[mirna])

#     # Get the dinucleotide frequencies for that miRNA's input sequence:
#     dinuc_path = get_analysis_path(mirna, exp, "I", "positional_kmers",
#                                    ext="_0_k2")
#     df_dinucs = pd.read_csv(dinuc_path, sep="\t", header=None)
#     df_dinucs.index = dinucs

#     # Create dictionaries for adding 5-prime and 3-prime nucleotides to the
#     # initial kmer.
#     prob_pos_5pnuc_map = {}
#     prob_pos_3pnuc_map = {}

#     # Create dictionaries for the nucleotide to choose based on the adjacent
#     # nucleotide sequence.
#     for nuc in nucs:
#         prob_pos_5pnuc_map[nuc] = {}
#         prob_pos_3pnuc_map[nuc] = {}
#         for pos in range(df_dinucs.shape[1]):
#             df_5p = df_dinucs.filter(regex="%s$" %(nuc), axis=0)[pos]
#             df_3p = df_dinucs.filter(regex="^%s" %(nuc), axis=0)[pos]
#             prob_pos_5pnuc_map[nuc][pos] = list(df_5p/df_5p.sum())
#             prob_pos_3pnuc_map[nuc][pos] = list(df_3p/df_3p.sum())
    
#     # Get the competitor oligo sequence, and all possible kmers that are
#     # reverse complementary to the competitor oligo:
#     comp_kmers = [comp_rc[i:i + len_comp]
#                   for i in range(len(comp_rc) - len_comp + 1)]
#     # Pre-allocate the output list:
#     reads = [None]*n

#     # Define a function to use for each read, that initializes the starting
#     # kmer. if the length of complementarity to the competitor is zero, the read
#     # is initialized by picking the 5-prime most dinucleotide of the read based
#     # on their dinucleotide frequencies. Other wise, the read is initialized
#     if len_comp == 0:
#         pos_kmer = 0
#         prob_start = list(df_dinucs[pos_kmer]/df_dinucs[pos_kmer].sum())
#         len_comp = 2
#         def kmer_start():
#             read = np.random.choice(dinucs, p=prob_start)
#             return((pos_kmer, read))
#     else:
#         def kmer_start():
#             pos_kmer = random.choice(list(range(len_rand - len_comp + 1)))
#             read = random.choice(comp_kmers)
#             return((pos_kmer, read))
#     for i in range(n):
#         if i != 0 and i % (n/10) == 0:
#             print((float(i)/n))
#             if i != 0:
#                 time_cur = time.time()
#                 d_time = time.time() - time_start
#                 time_sec = d_time/i*(n-i)
#                 time_h = time_sec // 3600
#                 time_m = (time_sec % 3600) // 60
#                 time_s = time_sec % 60
#                 print(("%.0f:%.0f:%.2f remaining" % (time_h, time_m, time_s)))


#         # Draw random position and random kmer of complementarity to the
#         # competitor oligo.
#         pos_kmer, read = kmer_start()
#         # else
#         # Build read starting with this random kmer at this random position.
#         # Extend the sequence 5-prime and 3-prime of the kmer to the full read
#         # length using the positional dinucleotide probabilities.
#         # Five prime:
#         for j in range(0, pos_kmer)[::-1]:
#             # probs = df_dinucs.filter(regex="%s$" %(read[0]), axis=0)[j]
#             # probs = probs/probs.sum()
#             nuc = np.random.choice(nucs, p=prob_pos_5pnuc_map[read[0]][j])
#             read = nuc + read
#         # Three prime:
#         for j in range(pos_kmer + len_comp - 1, df_dinucs.shape[1]):
#             # probs = df_dinucs.filter(regex="^%s" %(read[-1]), axis=0)[j]
#             # probs = probs/probs.sum()
#             nuc = np.random.choice(nucs, p=prob_pos_3pnuc_map[read[-1]][j])
#             read = read + nuc
#         # Add 3-prime-most constant sequence in the read, and update the list of
#         # reads.
#         read = read + read3p
#         reads[i] = read 
#     print_time_elapsed(time_start)
#     return(reads)


# def get_wobble_kmers(kmer_list):
#     # Make the output list of wobble kmers.
#     wob_kmer_list = []
#     # Define a dictionary that returns the appropriate wobble nucleotide.
#     wob_nuc_map = {"C" : "T", "A" : "G"}
#     # Iterate over the kmers
#     for kmer in kmer_list:
#         # Make a sublist of wobble kmers pertaining to this single kmer.
#         wob_kmer_sublist = [kmer]
#         # Iterate over the nucleotides of the kmer.
#         for i, nuc in enumerate(kmer):
#             # If a particular nucleotide position is a wobble position
#             if nuc in wob_nuc_map:
#                 # Make temporary list with all the current kmers.
#                 wob_kmer_sublist_temp = list(wob_kmer_sublist)
#                 # Add each of the pre-existing kmers in the list back in with
#                 # The additional wobble nucleotide at this position.
#                 for wob_kmer in wob_kmer_sublist:
#                     wob_kmer_sublist_temp.append(
#                         wob_kmer[:i] + wob_nuc_map[nuc] + wob_kmer[i+1:]
#                     )
#                 # Replace the old wobble sublist with the new wobble sublst.
#                 wob_kmer_sublist = list(wob_kmer_sublist_temp)
#         # Add each of these kmers to the new list, except for the first entry,
#         # Which is the original kmer.
#         wob_kmer_list += wob_kmer_sublist[1:]
#     return(wob_kmer_list)




def main():
    # print(get_miRNA_mfes("miR-1", 1, 8, ))
    return


################################################################################

if __name__ == "__main__":
    main()

