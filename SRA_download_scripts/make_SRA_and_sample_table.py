import gzip
import csv

# This script is for downloading the fastq files associated my graduate work.
# The script should go through the entire list of files and their associated
# metadata.

# Prior steps. I found the "SRP" files associated with each of the five GEO
# SubSeries associated with the SuperSeries GSE140220 (The biochemical basis of
# microRNA targeting efficacy): GSE140214, GSE140215, GSE140216, GSE140217,
# and GSE140218. These are SRP229580, SRP233292, SRP235216, SRP234771, and
# SRP235217. Searching each of these five accessions at
# 'https://www.ncbi.nlm.nih.gov/sra/' enables the download of a meta-data table
# by clicking "Send to:", "File", and "Summary." I named each file
# `$SRP_sra_result.csv`, where `$SRP` refers to the SRP accession. I also
# downloaded the accession list itself, using "Send to:", "File", and
# "Accession", naming each file "$SRP_SraAccList.txt". Finally, I downloaded
# each of the series matrix files from the GEO, using thie following command:
#
# `wget --recursive --no-parent -nd
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE140nnn/GSE140214/matrix/
# -P SRA_download_scripts'
#
# so that I would have the original file name associated with each file, which
# is important for identifying what the miRNA, experiment, and otherwise
# condition of each file is.

script_dir = "/n/groups/klein/mcgeary/gradschool/SRA_download_scripts"

# These two lists require manual curation currently in order for the script to
# work.
GSE_IDs = [
    "GSE140214",
    "GSE140215",
    "GSE140216",
    "GSE140217",
    "GSE140218",
    "GSE174715",
    "GSE196457",
]
SRP_IDs = [
    "SRP229580",
    "SRP233292",
    "SRP235216",
    "SRP234771",
    "SRP235217",
    "SRP320549",
    "SRP359105",
]
SRP_from_GSE_map = dict(zip(GSE_IDs, SRP_IDs))

## Functions required for the parsing the data.


def construct_directory(string):
    """Makes the directory associated with a SRA sample name.

    Makes the directory to be used in testing the GitHub scripts and
    re-downloading the data from McGeary et al., 2019 and 2022. This file does
    not deal with the ambiguities associated with input read files that are
    applicable to multiple miRNAs, such as those of let-7a.

    Args:
        string: The input string from which the directory information must be
        generated.

    Returns:
        A new string resembling the final directory structure of the final read
        file.

    """
    if "smallRNA-seq_purifiedAGO2" in string:
        mir = string.split("smallRNA-seq_purifiedAGO2-")[1]
        exp = "smallRNA_seq"
        cond = "AGO2_pur"
    elif ", " in string:
        exp, mir, cond = string.split(", ")
        exp = exp + "_let73p"
    elif "mpra" in string:
        exp, mir, cond = string.split("_")
    elif ("HeLa" in string) | ("HEK" in string):
        cell_line, mir, exp, rep, batch = string.split("_")
        exp = "%s_%s" % (exp, cell_line)
        cond = "%s_%s" % (rep, batch)
    else:
        str_split = string.split("_")
        mir, pur, cond = str_split[:3]
        pur = pur.replace("purification", "pur")
        cond = cond.replace("percent", "%")
        cond = cond.replace("inputlibrary", "input")
        cond = cond.replace("input2", "input_reseq")
        if len(str_split) == 4:
            cond += "_%s" % str_split[3]
        cond = "%s_%s" % (pur, cond)
        exp = "agorbns"
    dir_string = "%s/%s/reads/%s.txt" % (mir, exp, cond)

    return dir_string


def main():
    # Pre-allocate the list of strings to be used for the output file.
    output_list = []
    # GENERATE OUTPUT.__________________________________________________________
    # Iterate over the GSE IDs (each of the five accessions related to the five
    # data types in the 2019 paper).
    for GSE_use in GSE_IDs:
        # Define file name for the GEO matrix file used to get the GSM
        # accessions of each SRA file, and the name of sample (which enables
        # parsing what the sample actually corresponds to).
        SRP_use = SRP_from_GSE_map[GSE_use]
        GSE_file = "%s/%s_series_matrix.txt.gz" % (script_dir, GSE_use)
        SRP_summary_file = "%s/%s_sra_result.csv" % (script_dir, SRP_use)
        SRP_acc_file = "%s/%s_SraAccList.txt" % (script_dir, SRP_use)
        # Pre-allocate the list of names and GSM IDs, which have positional
        # correspondence.
        names = []
        GSM_IDs = []
        # Iterate through the file and find the two relevant rows, being the
        # sample name or the GSM ID row.
        with gzip.open(GSE_file, "rt") as f:
            for line in f:
                line_split = line.strip().split("\t")
                if line_split[0] == "!Sample_title":
                    names = [i.strip('"') for i in line_split[1:]]
                elif line_split[0] == '"ID_REF"':
                    GSM_IDs = [i.strip('"') for i in line_split[1:]]
        # Construct dictionary for retrival of the names using the GSM IDs.
        name_from_GSM_map = dict(zip(GSM_IDs, names))

        # Pre-allocate the GSM accession string. This step is necessary to map
        # the GSM accessions from the GEO metadata file to the SRR accessions in
        # the summary list. In principle, one might expect that the rows of the
        # two different files are the same, but I am not making that assumption.
        GSM_acc = []
        with open(SRP_summary_file) as csv_file:
            csv_scan = csv.DictReader(csv_file)
            for row in csv_scan:
                GSM_acc += [row["Experiment Title"].split(":")[0]]
        # Through manual inspection, identified that the GSM file is in opposite
        # order compared to accession list. I am therefore reversing the order
        # of this list.
        GSM_acc = GSM_acc[::-1]
        # Now opening the actual Accession List, which includes a list of SRP
        # accessions that can actually be used to download the files using the
        # prefetch and fastq-dump programs.
        with open(SRP_acc_file) as acc_file:
            SRR_acc = [i.strip() for i in acc_file.readlines() if len(i.strip()) != 0]
        # Construct a dictionary mapping the GSM accessions to the SRR
        # accessions.
        SRR_from_GSM_map = dict(zip(GSM_acc, SRR_acc))
        # Iterate over each GSM accession and get the corresponding SRR
        # accession and also determine the directory placement of the
        # corresponding file.
        for GSM_i in name_from_GSM_map.keys():
            SRR_use = SRR_from_GSM_map[GSM_i]
            dir_string = construct_directory(name_from_GSM_map[GSM_i])
            output_list += ["%s\t%s" % (SRR_use, dir_string)]

    # WRITE OUTPUT _____________________________________________________________
    # Assign the name of output file, and write the concatenated file
    output_dir = "%s/SRA_and_metadata.txt" % script_dir
    with open(output_dir, "w+") as file_out:
        file_out.write("\n".join(output_list) + "\n")


################################################################################
if __name__ == "__main__":
    main()
