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
    if experiment == "kinetics":
        pulse_reads_output_path = get_analysis_path(mirna, experiment, "I_TGT_combined", "full_reads")
        chase_reads_output_path = get_analysis_path(mirna, experiment, "I_ACA_combined", "full_reads")

        pulse_file_out_write = open(pulse_reads_output_path, "w")
        chase_file_out_write = open(chase_reads_output_path, "w")
        pulse_file_out_write.write("")
        chase_file_out_write.write("")
        pulse_file_out_write.close()
        chase_file_out_write.close()
    else:
        reads_output_path = get_analysis_path(mirna, experiment, "I_combined", "full_reads")

        file_out_write = open(reads_output_path, "w")
        file_out_write.write("")
        file_out_write.close()
    pairs = [("let-7a", "equilibrium"),
             ("miR-1", "equilibrium"),
             ("miR-124", "equilibrium"),
             ("let-7a", "kinetics"),
             ("miR-1", "kinetics"),
             ("miR-124", "kinetics"),
             ("lsy-6", "kinetics"),
             ("miR-1","kin_pilot","I_TGT"),
             ("miR-1","kin_pilot","I_ACA")]

    # Initialize the site assignments.

    for pair in pairs:
        # Define the mirna and experiment variables for analyzing the right
        # input experiment.
        print(pair)
        if len(pair) == 3:
            mirna_input, experiment_input, condition_temp = pair
        else:
            mirna_input, experiment_input = pair
            condition_temp = "I"

        # Assign the path for the input reads to be read (NOT to be written)
        reads_path = get_analysis_path(mirna_input,experiment_input,condition_temp,"full_reads")

        if experiment == "kinetics":
            with open(reads_path,"rb") as file_in:
                with open(pulse_reads_output_path, "a") as pulse_file_out_write:
                    with open(chase_reads_output_path, "a") as chase_file_out_write:
                        line = file_in.readline()
                        tick = 0
                        while line:
                            barcode = line[25+37+1:25+37+4]
                            if barcode == "TGC" and experiment_input == "equilibrium" and mirna_input == "miR-1":
                                line = line[:25+37+1]+"TGTTCGTATGCCGTCTTCTGCTTG\n"
                                pulse_file_out_write.write(line)
                            elif barcode == "TGT":
                                pulse_file_out_write.write(line)
                            elif barcode == "ACA":
                                chase_file_out_write.write(line)
                            if tick % 1000000 == 0:
                                print(tick)

                            tick += 1
                            line = file_in.readline()
        else:
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

