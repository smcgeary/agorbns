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
from general import get_rc, seq_mirna_map, parse_arguments, print_time_elapsed, get_analysis_path, multiprocess_file, readline

# This script reads one read file and outputs multiple .txt files to be used in all
# downstream applications. Inputs are:

# mirna: The miRNA in the experiment to be processed.
# experiment: The type of experiment to be processed.
# condition: The sampel type within the experiment to be processed.

# Functions from general this uses:


# Constants

# FUNCTIONS



def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["miRNA","experient","condition"]
    mirna, experiment, condition = parse_arguments(arguments)

 
    reads_path = get_analysis_path(mirna,experiment,condition,"reads")
    fasta_path = get_analysis_path(mirna,experiment,condition,"reads_fasta","fa")

    with open(reads_path,"rb") as file_in:
        with open(fasta_path,"wb") as file_out:
            for i, line in enumerate(file_in.readlines()):
                file_out.write(">%s\n%s" %(str(i), line))

    # Print the amount of time the script took to complete.
    print_time_elapsed(time_start)


################################################################################

if __name__ == "__main__":
    main()

