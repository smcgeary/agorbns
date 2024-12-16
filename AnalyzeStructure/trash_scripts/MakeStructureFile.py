################################################################################
#MakeSiteTypeReadFiles.py
################################################################################
import imp # Used to import general.py
import time # Used for time
import math
import itertools as it
import multiprocessing
import sys
imp.load_source("general",
                ("/lab/bartel1_ata/mcgeary/computation/"
                 "AgoRBNS/general/general.py"
                )
               )
from general import parse_arguments, print_time_elapsed, get_analysis_path, multiprocess_file, readline
# This script reads one read file and outputs multiple .txt files to be used in all
# downstream applications. Inputs are:

# mirna: The miRNA in the experiment to be processed.
# experiment: The type of experiment to be processed.
# condition: The sampel type within the experiment to be processed.

def get_structures(args_sub):
    time_start = time.time()
    temp_file, bp_prob_path, mfe_path, read_seqs = args_sub
    bp_prob_temp = bp_prob_path.split(".txt")[0]+"_"+str(temp_file)+".txt"
    print(bp_prob_temp)
    sys.stdout.flush()
    mfe_prob_temp = mfe_path.split(".txt")[0]+"_"+str(temp_file)+".txt"
    with open(bp_prob_temp,"w") as file_out:
        print(file_out)
        sys.stdout.flush()
        file_out.write("")
    with open(mfe_prob_temp,"w") as file_out:
        file_out.write("")
    import RNA

    """Takes a read sequence and calculates the base pairing probabilities
    at each position.

    Args:
        read_sequence: The read sequence.

    Returns:
        A tuple containing the list the basepairing probabilites.
    """

    # Sub-function to be used on each site:
    probs_out = range(len(read_seqs))
    mfe_out = range(len(read_seqs))
    # print(read_seqs)
    pairs = list(it.combinations(range(1,len(read_seqs[0].strip())+1), 2))
    with open(bp_prob_temp,"a") as bp_temp:
        with open(mfe_prob_temp,"a") as mfe_temp:
            for i, read_seq in enumerate(read_seqs):
                if i % 1000 == 0:
                    print("%sk" %(i/1000))
                    print_time_elapsed(time_start)
                    sys.stdout.flush()
                read_seq = read_seq.strip()

                mfe = RNA.fold(read_seq)[1]
                mfe_out[i] = mfe
                pf = RNA.pf_fold(read_seq)[1]
                mfe_temp.write("%s\n" %(str(mfe)))
                probs = [0]*len(read_seq)
                for pair in pairs:
                    [p1, p2] = pair
                    prob = RNA.get_pr(p1, p2)
                    probs[int(p1) - 1] += prob
                    probs[int(p2) - 1] += prob
                probs_out[i] = probs
                bp_temp.write("%s\n" %("\t".join([str(j) for j in probs])))
    return probs_out, mfe_out


def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["miRNA","experient","condition"]
    mirna, experiment, condition = parse_arguments(arguments)

    # Get the path to the read file and to that of where the site labels will
    # be written.
    reads_path = get_analysis_path(mirna,experiment,condition,"full_reads")
    bp_prob_path = get_analysis_path(mirna,experiment,condition,"structures_bp_prob")
    mfe_path = get_analysis_path(mirna,experiment,condition,"structures_mfe")

    with open(reads_path,"rb") as file_in:
        print("about to read lines")
        file_reads = file_in.readlines()
        # print(file_reads)
        print("about to split reads")
        file_sub_reads = [file_reads[x:x+1000000] for x in range(0, len(file_reads), 1000000)]
        args_sub = [(i,bp_prob_path,mfe_path,sub_reads) for i, sub_reads in enumerate(file_sub_reads)]
        print(len(file_sub_reads))
        print([len(i) for i in file_sub_reads])
        print("making pools")
        pool = multiprocessing.Pool(22)
        # print(file_sub_reads)
        print("starting multiprocess")
        results = pool.map(
            get_structures, args_sub)
    # print(results)
    # Collect the output from each thread, which are the list
    # of reads and the dictionary of read types.
    bp_prob_threads = [i[0] for i in results]
    mfe_threads = [i[1] for i in results]
    # print(bp_prob_threads[0])
    # Flatten the list into one list (found on stackoverflow) and write
    # it to its output file.
    # print(mfe_threads[0])
    bp_prob_list = [i for sublist in bp_prob_threads for i in sublist]
    with open(bp_prob_path,"wb") as file_out:
        file_out.write("".join(["\t".join([str(j) for j in i])+"\n" for i in bp_prob_list]))

    # Flatten the list into one list (found on stackoverflow) and write
    # it to its output file.
    mfe_list = [i for sublist in mfe_threads for i in sublist]
    with open(mfe_path,"wb") as file_out:
        file_out.write("".join(["%s\n" % (i) for i in mfe_list]))


    # Print the amount of time the script took to complete.
    print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()

