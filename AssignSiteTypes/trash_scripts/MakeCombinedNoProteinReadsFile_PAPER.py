################################################################################
#MakeSiteTypeReadFiles.py
################################################################################
import imp # Used to import general.py
import time # Used for time
imp.load_source("general",
                ("/lab/bartel1_ata/mcgeary/computation/"
                 "AgoRBNS/general/general.py"
                )
               )
from general import parse_arguments, print_time_elapsed, get_analysis_path, multiprocess_file, multiprocess_test, readline
from operator import add
from sitetypes import get_seq_site_map
# This script reads one read file and outputs multiple .txt files to be used in all
# downstream applications. Inputs are:

# mirna: The miRNA in the experiment to be processed.
# experiment: The type of experiment to be processed.
# condition: The sampel type within the experiment to be processed.

# Functions from general this uses:


# Constants

# FUNCTIONS

# TODO: comment this.

def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["miRNA", "experient"]
    print(parse_arguments(arguments))
    mirna, experiment = parse_arguments(arguments)
    reads_output_path = get_analysis_path(mirna, experiment, "I_combined", "full_reads")
    with open(reads_output_path, "w") as file_out_write:
        file_out_write.write("")

    pairs = [("let-7a", "equilibrium"),
             ("miR-1", "equilibrium"),
             ("miR-124", "equilibrium"),
             ("let-7a", "kinetics"),
             ("miR-1", "kinetics"),
             ("miR-124", "kinetics"),
             ("lsy-6", "kinetics"),
             ("miR-1","kinetics_pilot")]

    # Initialize the site assignments.

    for pair in pairs:
        # Define the mirna and experiment variables for analyzing the right
        # input experiment.
        print(pair)
        mirna_input, experiment_input = pair
        if experiment_input == "kinetics_pilot":
            condition_temp = "I_TGT"
        else:
            condition_temp = "I"

        # Assign the path for the input reads to be read (NOT to be written)
        reads_path = get_analysis_path(mirna_input,experiment_input,condition_temp,"full_reads")


        with open(reads_path,"rb") as file_in:
            with open(reads_output_path, "a") as file_out_write:
                if mirna == "miR-1" and experiment == "equilibrium":
                    line = file_in.readline()
                    tick = 0
                    while line:
                        barcode = line[25+37+1:25+37+4]
                        if barcode in ["TGC", "TGT"]:
                            line = line[:25+37+1]+"TCGTATGCCGTCTTCTGCTTG\n"
                            file_out_write.write(line)
                        line = file_in.readline()
                        tick +=1
                        if tick % 1000000 == 0:
                            print(tick)
                elif experiment == "equilibrium":
                    line = file_in.readline()
                    while line:
                        barcode = line[25+37+1:25+37+4]
                        if barcode in ["TGC", "TGT"]:
                            line = line[:25+37+1]+"TGTTCGTATGCCGTCTTCTGCTTG\n"
                            file_out_write.write(line)
                        line = file_in.readline()

  

    # Print the amount of time the script took to complete.
    print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()

