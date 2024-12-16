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
import regex

# I first used the above address, but the bottom one gives near-perfect (but
# not perfect) agreement with Kathy's files, in the folder:
# /lab/solexa_bartel/klin/miRNA_models_data/transfections/hela/compiled/ .

base_exp_dir = "/lab/solexa_bartel/mcgeary/transfections/"
star_index = ("/nfs/genomes/human_gp_feb_09_no_random/STAR_2.6.0/"
              "GRCh37.75_overhang_100/")
annot_file = "%smetadata/annotations.gtf" %base_exp_dir
temp_dir   = "%sSTAR_alignment_temp" %base_exp_dir
raw_files_dir = "/lab/solexa_bartel/mcgeary/GEO_McGearyLin2019/raw_files/"

measured_kds_dir = "/lab/solexa_bartel/klin/miRNA_models_data_old/model_inputs/biochem/measured_kds/"


CONVERT_MIRNA_NAMES_SPECIAL = {"let-7a" : "let7" ,
                               "lsy-6" : "lsy6"}



def make_feature_dictionary(cell_line):
    transcript_file_dir = ("/lab/solexa_bartel/mcgeary/transfections"
                           "/metadata/%s_transcripts_2.txt" %(cell_line))
    # Make a dictionary that has the information on each transcript in order to
    # ensure I am performing the analysis the same way that Kathy was doing it.
    with open(transcript_file_dir) as transcript_file:
        # Make the outer layer of keys.
        outer_keys = transcript_file.readline().strip().split("\t")[1:]
        features_mRNA_map = {i : collections.defaultdict(lambda: "Key does not exist") for i in outer_keys}
        line = transcript_file.readline()
        while line:
            line_list = line.strip().split("\t")
            transcript_id = line_list[0]
            features = line_list[1:]
            for i, item in enumerate(features):
                features_mRNA_map[outer_keys[i]][transcript_id] = item
            line = transcript_file.readline()
    return(features_mRNA_map)


def convert_mirna_name_kathy(mirna):
    if mirna in CONVERT_MIRNA_NAMES_SPECIAL.keys():
        return(CONVERT_MIRNA_NAMES_SPECIAL[mirna])
    else:
        return "ir".join(mirna.split("iR-"))


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
            print(int(10.*i_row/float(num_rows)))
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
    arguments = ["mirna", "cell_line"]
    args = parse_arguments(arguments)
    (mirna, cell_line) = args

    # Make the features dictionary that has the various sequences and lengths
    # associated with the annotations of each mRNA.
    features_mRNA_map = make_feature_dictionary(cell_line)

    # Scripting to generate the name of the file that contains the input data
    # for the biochemical model.
    mirna_file = convert_mirna_name_kathy(mirna) + ".txt"
    measured_kds_file = open(measured_kds_dir + mirna_file)

    # Grab the first line of this file because this has the names of each of the
    # features in the model.
    site_features = measured_kds_file.readline().strip().split("\t")

    # THIS PART WILL GO INTO A LOOP OR A COMPREHENSION EVENTUALLY. #############
    # Make sure that all of the orf and UTR sequence features that I am using
    # agree with those generated by Kathy. These should be, 1.) that the 12-nt
    # kmer is indeed in the mRNA at the location that it is, and 2.) that the
    # ORF length and 3'UTR length is also the same.
    line = measured_kds_file.readline()
    i_line = 0
    tally_good = 0
    tally_multi = 0
    tally_edge = 0
    tally_none = 0

    multi_dist = collections.defaultdict(lambda : 0)
    N_b = collections.defaultdict(lambda : 0)
    N_bg = collections.defaultdict(lambda : 0)
    # TODO Make sure this is coded in a better way.
    ag_mirna_map = {"lsy-6"   : -5.224992752075195,
                    "miR-1"   : -4.526782989501953,
                    "miR-124" : -5.244050979614258,
                    "miR-155" : -4.7737507820129395,
                    "miR-7"   : -5.3175435066223145}
    ag = math.exp(ag_mirna_map[mirna])
    b = math.exp(0.8655766248703003)
    c = math.exp(-1.848806619644165)
    transcript_pre = line.strip().split("\t")[0]
    while line:
        print_next = False
        line_dict = {i : j for (i, j) in zip(site_features, line.strip().split("\t"))}
        transcript = line_dict["transcript"]
        kmer12 = line_dict["12mer"]
        in_orf = line_dict["in_ORF"]
        print(transcript)
        print(line_dict["utr3_loc"])
        print(in_orf)
        kd = math.exp(float(line_dict["log_kd"]))
        if in_orf:
            N_b[transcript] += ag/(ag + c*kd)
            N_bg[transcript] += ag/(ag + c)
        else:
            N_b[transcript] += ag/(ag + kd)
            N_bg[transcript] += ag/(ag + 1)
        #     print(N_b)
        #     print(N_bg)
        # elif transcript == "NM_000031.5":
        #     break


        # orf_seq = features_mRNA_map["orf"][transcript]
        # utr3_seq = features_mRNA_map["utr3"][transcript]
        # orf_utr3_seq = features_mRNA_map["orf_utr3"][transcript]
        # orf_length = int(features_mRNA_map["orf_length"][transcript])
        # ind_orf_utr3 = orf_utr3_seq.find(kmer12)
        # utr3_loc_sm = ind_orf_utr3 - orf_length + 3
        # utr3_loc_kl = int(line_dict["utr3_loc"])
        # if utr3_loc_sm != utr3_loc_kl:
        #     alt_find = find_all(orf_utr3_seq, kmer12)
        #     if len(alt_find) > 1:
        #         check = [utr3_loc_kl == i_try - orf_length + 3 for i_try in alt_find]
        #         somewhere_within = np.multiply(check, 1).sum()
        #     else:
        #         somewhere_within = 0
        #     edge_kmer = "X" in kmer12
        #     # print(num_instances_orf_utr)
        #     if somewhere_within:
        #         # multi_dist[str(num_instances_orf_utr)] += 1
        #         tally_multi += 1
        #     elif edge_kmer:
        #         tally_edge += 1
        #     else:
        #         print("XXXXXXXXXXXXXXXXXXXXXXXX")
        #         print(kmer12)
        #         print(i_line)
        #         print(num_instances_orf_utr)
        #         print(edge_kmer)
        #         print(utr3_loc_sm)
        #         print(utr3_loc_kl)
        #         tally_none += 1
        #         break
        # else:
        #     if "X" in kmer12:
        #         print(utr3_loc_sm)
        #         print(utr3_loc_kl)
        #         print(orf_length)
        #         print(kmer12)
        #     tally_good += 1

        # # print(type(utr3_loc_kl))
        # # print(type(utr3_loc_sm))
        # # print(orf_seq[-3:])
        # # print(utr3_seq[:3])
        # # print(str(utr3_loc_sm) + "\t" + str(utr3_loc_kl))
        # transcript_pre = transcript
        # line_dict_pre = line_dict.copy()
        # utr3_loc_sm_pre = utr3_loc_sm

        line = measured_kds_file.readline()
        i_line += 1
    out_put = dict({trans : math.log((1 + b*N_bg[trans])/(1 + b*N_b[trans]))
                    for trans in N_b.keys()})
    output_df = pd.DataFrame.from_dict(out_put, orient="index")
    print(output_df)
    
    return()
    print(N_b)
    print(N_bg)
    print(math.log((1 + b_param*N_bg)/(1 + b_param*N_b)))

    print("Done script")
    print(tally_good)
    print(tally_multi)
    print(tally_edge)
    print(tally_none)
    print(multi_dist)

    total_kmers = tally_good + tally_multi + tally_edge + tally_none

    print(tally_good/float(total_kmers))
    print(tally_multi/float(total_kmers))
    print(tally_edge/float(total_kmers))
    print(tally_none/float(total_kmers))

    # # Define path to input data and read in table.
    # path = "%s%s/count_tables/raw.txt" %(base_exp_dir, cell_line)
    # gXtpm = pd.read_csv(path, index_col=0, sep="\t")
    # # Perform the batch normalization.
    # gXtpm_bn = batch_normalize(gXtpm)
    # # Get the output path, and write the table to thie path.
    # path = "%s%s/count_tables/logtpm_batchnormalized.txt" %(base_exp_dir, cell_line)
    # print(path)
    # gXtpm_bn.to_csv(path, sep="\t") 
    # # Print the total amount of time the script took.
    # time_done = time.time()
    # print(("Finished job future in %f3.2 seconds." %(time_done - time_start)))
    return()





if __name__ == "__main__":
    main()
