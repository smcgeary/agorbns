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

def main():

    arguments = ["i","start_stop","bp_prob_path","mfe_path","reads_path"]
    temp_file, start_stop, bp_prob_path, mfe_path, reads_path = parse_arguments(arguments)
    start, stop = start_stop
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

    with open(read_seqs,"rb") as file_in:
        read_seqs = file_in.readlines()[start:stop]
    probs_out = range(len(read_seqs))
    mfe_out = range(len(read_seqs))
    # print(read_seqs)
    pairs = list(it.combinations(range(len(read_seqs[0])), 2))
    print(pairs)
    with open(bp_prob_temp,"a") as bp_temp:
        with open(mfe_prob_temp,"a") as mfe_temp:
            for i, read_seq in enumerate(read_seqs):
                if i % 1000 == 0:
                    print("%sk" %(i/1000))
                    print_time_elapsed(time_start)
                    sys.stdout.flush()
                read_seq = read_seq.strip()
                fix = True
                while fix:
                    mfe = RNA.fold(read_seq)[1]
                    mfe_out[i] = mfe
                    pf = RNA.pf_fold(read_seq)[1]
                    mfe_temp.write("%s\n" %(str(mfe)))
                    probs = [0]*len(read_seq)
                    for pair in pairs:
                        [p1, p2] = pair
                        prob = RNA.get_pr(p1, p2)
                        prob = float(prob)
                        probs[int(p1) - 1] += prob
                        probs[int(p2) - 1] += prob
                    fix = False
                    attempt = 1
                    for i_p, p in enumerate(probs):
                        if p > 1 or p < 0 or math.isnan(p):
                            attempt +=1
                            print("failed")
                            print(p)
                            print(read_seq)
                            print(" "*i_p+read_seq[i_p])
                            print(p)
                            sys.stdout.flush()
                            fix = True
                            break
                probs_out[i] = probs
                bp_temp.write("%s\n" %("\t".join([str(j) for j in probs])))
    return probs_out, mfe_out




    # Print the amount of time the script took to complete.
    print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()

