################################################################################
#repression_test.py
################################################################################
import imp # Used to import general.py
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

# I first used the above address, but the bottom one gives near-perfect (but
# not perfect) agreement with Kathy's files, in the folder:
# /lab/solexa_bartel/klin/miRNA_models_data/transfections/hela/compiled/ .

base_exp_dir = "/lab/solexa_bartel/mcgeary/transfections/"
star_index = ("/nfs/genomes/human_gp_feb_09_no_random/STAR_2.6.0/"
              "GRCh37.75_overhang_100/")
annot_file = "%smetadata/annotations.gtf" %base_exp_dir
temp_dir   = "%sSTAR_alignment_temp" %base_exp_dir
raw_files_dir = "/lab/solexa_bartel/mcgeary/GEO_McGearyLin2019/raw_files/"



def make_data_table(cell_line):
    count_dir = "%s%s/STAR_htseq/" %(base_exp_dir, cell_line)

    # Get the orf lengths for the TPM calculation:
    orf_dir = "%s/metadata/%s_orf_lengths.txt" %(base_exp_dir,
                                                 cell_line)
    orf_lengths = pd.read_csv(orf_dir, header=None, index_col=0, sep="\t")
    shell_string = "ls %s*" %count_dir
    file_str = subprocess.run(
        shell_string,
        check=True,
        shell=True,
        stdout=PIPE).stdout.decode("utf-8")
    file_list = file_str.strip().split("\n")
    print(file_list[0])
    mirnas = [file.split("/")[-1].split("_")[1] for file in file_list]
    reps = [file.split("/")[-1].split("_")[3].split("rep")[1] for file in file_list]
    batches = [file.split("/")[-1].split("_")[4].split("batch")[1].split(".")[0]
               for file in file_list]
    # Make the dataframe with metadata reporting on the miRNA, rep, and batch
    # associated with each sample.
    metadata = pd.DataFrame(np.array([mirnas, reps, batches]),
                            index=["miRNA", "rep", "batch"],
                            columns=list(range(len(file_list))))
    print(metadata)
    # Make the data table by iterating over all of the columns.
    gXc = pd.concat(
        [pd.read_csv(file, header=None, index_col=0, sep="\t")
         for file in file_list],
        axis=1
    )
    # Add column labels to the count data:
    gXc.columns = list(range(len(file_list)))
    # First exclude rows based on not having a "NM" number.
    rows_exclude = list(filter(lambda x: re.search(r'__', x), gXc.index))
    gXc.drop(index=rows_exclude, inplace=True)
    # Next exclude mRNAs that do not have at least 10 reads in each sample.
    gXc_check = gXc >= 10
    print(orf_lengths.shape)
    gXc = gXc.loc[gXc_check.all(axis=1), ]
    print(gXc.shape)
    orf_lengths = orf_lengths.loc[gXc.index, ]
    print(orf_lengths)
    orf_lengths.columns.name = None
    orf_lengths.rename_axis(None, axis = 1)
    print(orf_lengths)
    print(gXc)
    print(gXc.iloc[:4, :4])
    print(orf_lengths.iloc[:4, :])
    gXc_length_norm = gXc.divide(orf_lengths.to_numpy(), axis=0)
    print(gXc_length_norm.iloc[:4, :4])
    gXc = np.log(gXc_length_norm.divide(gXc_length_norm.sum(axis=0), axis=1)*1e6)
    # Append the meta data containing the miRNA, rep, and batch information to
    # the data table.
    gXc = pd.concat([metadata, gXc])
    return(gXc)




def main():
    time_start = time.time()
    # Parse command line arguments:
    arguments = ["cell_line"]
    args = parse_arguments(arguments)
    print(args)
    (cell_line,) = args
    print(cell_line)
    # Assign the cell-line argument.
    gXc = make_data_table(cell_line)
    print(gXc)
    # Define the path and output the count table.
    path = "%s%s/count_tables/raw.txt" %(base_exp_dir, cell_line)
    print(path)
    gXc.to_csv(path, sep="\t") 
    return()





if __name__ == "__main__":
    main()
