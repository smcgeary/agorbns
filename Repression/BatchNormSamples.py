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

import statsmodels.formula.api as smf

# I first used the above address, but the bottom one gives near-perfect (but
# not perfect) agreement with Kathy's files, in the folder:
# /lab/solexa_bartel/klin/miRNA_models_data/transfections/hela/compiled/ .

base_exp_dir = "/lab/solexa_bartel/mcgeary/transfections/"
star_index = ("/nfs/genomes/human_gp_feb_09_no_random/STAR_2.6.0/"
              "GRCh37.75_overhang_100/")
annot_file = "%smetadata/annotations.gtf" %base_exp_dir
temp_dir   = "%sSTAR_alignment_temp" %base_exp_dir
raw_files_dir = "/lab/solexa_bartel/mcgeary/GEO_McGearyLin2019/raw_files/"



def batch_normalize(gXtpm):
    # Define the miRNA and batch variables to be used with each row.
    miRNA = gXtpm.iloc[0, ]
    batch = gXtpm.iloc[2, ]
    # Remove the metadata from the top of the dataframe.
    gXtpm_data = gXtpm.iloc[3:, ]
    # Define the output data frame, which will have the batch normalization
    # coefficients for the miRNA and the batches.
    gXtpm_bn = pd.DataFrame(
        np.zeros(shape=(gXtpm_data.shape[0],
        len(miRNA.unique()) + len(batch.unique()) - 1))
    )
    # Name the indeces and the columns of the pre-allocated output matrix.
    # Omitting the "1" batch, as this is the reference point for the linear
    # model fitting. The dimensions are given by the remaining rows of the data
    # matrix (each gene) and the total number of unique miRNA categories and the
    # total number of batches minus the reference batch of 1.
    gXtpm_bn.index = gXtpm_data.index
    gXtpm_bn.columns = (list(np.sort(miRNA.unique())) +
                        list(np.sort(batch.unique()))[1:])
    # Iterate over each row of the data matrix:
    num_rows = gXtpm_data.shape[0]
    print(num_rows/10)
    print(int(num_rows/10.))
    for i_row in range(num_rows):
        if i_row % (int(num_rows/10.)) == 0:
            print(int(10*i_row/float(num_rows)))
        # Make a specific dataframe that has the logtpm values, the miRNA, and
        # the batch.
        ltpm = pd.to_numeric(gXtpm_data.iloc[i_row, ])
        ltpm.name = "ltpm"
        df_fit = pd.concat([ltpm, miRNA, batch], axis=1)
        # Fit the model, extract the coefficients.
        params = smf.ols(formula='ltpm ~ miRNA + C(batch) - 1',
                         data=df_fit).fit().params
        # Assign the coefficients to the output row.
        gXtpm_bn.iloc[i_row, ] = np.array(params)
    # Output the batch-normalized, log-tpm matrix.
    return(gXtpm_bn)



def main():
    time_start = time.time()
    # Parse command line arguments:
    arguments = ["cell_line"]
    args = parse_arguments(arguments)
    (cell_line,) = args
    # Define path to input data and read in table.
    path = "%s%s/count_tables/raw.txt" %(base_exp_dir, cell_line)
    gXtpm = pd.read_csv(path, index_col=0, sep="\t")
    # Perform the batch normalization.
    gXtpm_bn = batch_normalize(gXtpm)
    # Get the output path, and write the table to thie path.
    path = "%s%s/count_tables/logtpm_batchnormalized.txt" %(base_exp_dir, cell_line)
    print(path)
    gXtpm_bn.to_csv(path, sep="\t") 
    # Print the total amount of time the script took.
    time_done = time.time()
    print(("Finished job future in %f3.2 seconds." %(time_done - time_start)))
    return()





if __name__ == "__main__":
    main()
