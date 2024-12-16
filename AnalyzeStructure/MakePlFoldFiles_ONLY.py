################################################################################
#MakeSiteTypeReadFiles.py
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

def randomword(length):
   return ''.join(random.choice(string.lowercase) for i in range(length))

# This script reads one read file and outputs multiple .txt files to be used in all
# downstream applications. Inputs are:

# mirna: The miRNA in the experiment to be processed.
# experiment: The type of experiment to be processed.
# condition: The sampel type within the experiment to be processed.

# Functions from general this uses:


# Constants
kLFLANK, kRFLANK, kWIN = 16, 16, 30

# FUNCTIONS
def RNAplfold(read, win, mir1p, file_num, logprob=False):
    # This should give the probability of the region in the read complementary
    # to position "mir_start" to "mir_stop" in the miRNA sequence.
    read_name = randomword(20)
    temp_fa_filename = "%s%s.fa" %(read_name, file_num)
    with open(temp_fa_filename, "wb") as f:
        f.write(">%s%s\n%s" %(read_name, file_num, read))
    # call RNAplfold:
    l_param = len(read)
    w_param = l_param
    u_param = win
    if logprob:
        log_param = "-O "
        temp_plfold_filename = '%s%s_openen' %(read_name, file_num)
    else:
        log_param = ""
        temp_plfold_filename = '%s%s_lunp' %(read_name, file_num)

    call = 'RNAplfold -L %s -W %s -u %s %s< %s%s.fa' %(l_param, w_param,
                                                       u_param, log_param,
                                                       read_name, file_num)
    subprocess.call([call], shell=True, stdout = subprocess.PIPE)
    rnaplfold_data = pd.read_csv(temp_plfold_filename, sep = '\t',
                                 header = 1).set_index(' #i$').iloc[:, :win].transpose()
    for row in range(rnaplfold_data.shape[0]):
        rnaplfold_data.iloc[row, ] = (rnaplfold_data.iloc[row, ].shift(-((row + 1)/2)))**(1/float(row + 1))
    os.remove(temp_fa_filename)
    os.remove(temp_plfold_filename)
    os.remove('%s%s_dp.ps' %(read_name, file_num))
    return rnaplfold_data


def get_plfold_matrix(read_seqs, _mirna, _sitelist, experiment, site,
                      n_constant, rand_length, buffer3p, num_i):
    """Takes a read sequence and identifies the position of all site types'.

    Args:
        read_sequence: The read sequence.
        sites: A list of sequences corresponding to each site type.

    Returns:
        A tuple containing the list of reads and a dictionary of read types.
    """


    flanks = get_kmer_list(4)

    df_columns = (["f5p%s" %(i) for i in range(kLFLANK, 0, -1)] +
                  ["sp%s" %(i) for i in range(8, 0, -1)] + 
                  ["f3p%s" %(i) for i in range(1, kRFLANK + 1)])

    df_index = [str(i + 1) for i in range(kWIN)]

    plmean_flanks_map = {
        flank: pd.DataFrame(0.0, index=df_index, columns=df_columns)
        for flank in flanks
    }

    n_flanks_map = {
        flank: pd.DataFrame(0.0, index=df_index, columns=df_columns)
        for flank in flanks
    }

    time_start = time.time()
    tick = 0
    stop_max = 0
    for i_r, r in enumerate(read_seqs):
        # Get the read number:
        read = r.strip()
        # Make Read object.
        _read = Read(read, rand_length, _mirna, n_constant, experiment)
        sys.stdout.flush()
        # Find all sites within the read
        add_sites_from_read(_sitelist, _read, buffer_=buffer3p)
        if _read.sites[0].name != "None":
            read_sites = [i.name for i in _read.sites]
            if site in read_sites:
                sys.stdout.flush()
                ind_site = [i for i, site_i in enumerate(read_sites)
                             if site_i == site][0]
                start, stop = _read.site_pos(full=True)
                site_i = _read.sites[ind_site]
                mirp1_span = _sitelist.m1pos[site_i.name] - 1
                mirp1 = start + mirp1_span
                # print(mirp1_span)
                # print(_read.fullseq)
                # print(" "*start + "_"*(stop - start))
                # print(" "*mirp1 + _read.fullseq[mirp1])
                sys.stdout.flush()
                # Perform plfold on the entire read:
                headers = n_flanks_map["AAAA"].columns
                plfold = RNAplfold(_read.fullseq, kWIN, mirp1, num_i)
                l_pl = mirp1 - 8 + 1 - kLFLANK
                r_pl = mirp1 + kRFLANK + 1
                plfold_subset = plfold.iloc[:, l_pl:r_pl]
                plfold_subset.index = df_index
                plfold_subset.columns = df_columns
                flank = "".join([_read.fullseq[i] for i in
                 [start - 2, start - 1, stop, stop + 1]])
                # print(" "*(start - 2) +flank[:2] + " "*(stop - start) + flank[2:])
                # if flank == "ATTA":
                #     print(" "*l_pl + "6543210987654321" + "O.,..,.O" + "1234567890123456")
                #     for row in range(plfold_subset.shape[0]):
                #         strings = "".join([str(i)[0] for i in plfold_subset.iloc[row, :]])
                #         print(" "*l_pl + strings)
                #     # print(plfold_subset)
                #     sys.stdout.flush()
                mean_old = plmean_flanks_map[flank] * n_flanks_map[flank]
                mean_i = np.log(plfold_subset)
                n_i = (mean_i * 0).replace(0, 1).fillna(0)
                n_new = n_flanks_map[flank] + n_i
                mean_new = (mean_old + mean_i.fillna(0)) / n_new.replace(0, 1)
                plmean_flanks_map[flank] = mean_new
                n_flanks_map[flank] = n_new
                # if flank == "ATTA":
                #     print(old_weight)
                #     print(new_mean*0)
                #     # non_zero = non_zero.replace(NaN, 0)
                #     print(non_zero)
                #     print(n_flanks_map_alt["ATTA"])
                #     sys.stdout.flush()
    return [plmean_flanks_map, n_flanks_map]




def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["miRNA", "experiment", "condition", "n_constant", "sitelist",
                 "site", "-buffer3p_binary", "-jobs", "-test_binary"]

    # Parse arguments and assign script parameters accordingly:
    (
        mirna, experiment, condition, n_constant, sitelist, site, buffer3p,
        jobs, test
    ) = parse_arguments(arguments)

    if experiment in ["equil_mmseed_nb", "equil_seed_nb"]:
        read_length = 38
    else:
        read_length = 37
    if condition == "I_combined":
        input_list = INPUT_LIST_I_COMBINED
        if not jobs:
            jobs = 20
    elif condition == "0_combined":
        input_list = INPUT_LIST_0_COMBINED
        if not jobs:
            jobs = 20
    else:
        input_list = [(mirna, experiment, condition)]
        if not jobs:
            jobs = 20

    # Call instance variables:
    _mirna = Mirna(mirna)
    _sitelist = SiteList(_mirna, sitelist,
                         int(read_length) + 2*int(n_constant))

    # Assign arguments to be passed to the function
    # "get_plfold_matrix.""
    args  = [_mirna, _sitelist, experiment, site, int(n_constant), read_length,
             buffer3p]

    # Move to the temporary folder for making the pl fold files:
    cwd = os.getcwd()
    folder_name = randomword(20)
    plfold_path = "%s/AnalyzeStructure/RNAplfold_temp/%s" %(HOME_DIR,
                                                            folder_name)
    os.makedirs(plfold_path)
    os.chdir(plfold_path)

    # For each file in the input list, perform the multiprocessing:
    threads = []
    for i_input in input_list:
        print(i_input)
        reads_path = get_analysis_path(i_input[0], i_input[1], i_input[2],
                                       "reads")
        if not jobs:
            jobs = 1
        threads += multiproc_file(reads_path, int(jobs), get_plfold_matrix,
                                  test, *args)
    # Move back to the home directory and delete the folder just made
    os.chdir(cwd)
    os.rmdir(plfold_path)


    flanks = get_kmer_list(4)


    df_columns = (["f5p%s" %(i) for i in range(kLFLANK, 0, -1)] +
                ["sp%s" %(i) for i in range(8, 0, -1)] + 
                ["f3p%s" %(i) for i in range(1, kRFLANK + 1)])

    df_index = [str(i + 1) for i in range(kWIN)]

    plmean_flanks_map_final = {
        flank: pd.DataFrame(0.0, index=df_index, columns=df_columns)
        for flank in flanks
    }
    n_flanks_map_final = {
        flank: pd.DataFrame(0, index=df_index, columns=df_columns)
        for flank in flanks
    }
    print(len(threads))
    for flank in flanks:
        # print(flank)
        for thread in threads:
            mean_old = plmean_flanks_map_final[flank]
            n_old = n_flanks_map_final[flank]

            mean_i = thread[0][flank]
            n_i = thread[1][flank]
            n_new = n_old + n_i

            mean_new = (mean_old*n_old + mean_i*n_i)/n_new.replace(0, 1)
            plmean_flanks_map_final[flank] = mean_new
            n_flanks_map_final[flank] = n_new



    df_final_index = ["%s_w%s" %(i[1], i[0])
                      for i in it.product(list(plmean_flanks_map_final["AAAA"].index),
                                   list(plmean_flanks_map_final["AAAA"].columns))]
    df_final_columns = flanks

    outputmatrix = pd.DataFrame(0.0, index=df_final_index, columns=df_final_columns)

    for flank in flanks:
        flank_matrix = pd.Series(plmean_flanks_map_final[flank].values.flatten())
        flank_matrix.index = df_final_index
        outputmatrix.loc[df_final_index, flank] = flank_matrix


    # Name the output file extension.
    extension = "_%s_%s_%s" %(n_constant, sitelist, site)
    if buffer3p:
        extension = "%s_buffer3p" %(extension)
    if test:
        extension = "%s_test" %(extension)

    output_path = get_analysis_path(
        mirna, experiment, condition,
            "plfold_2018_PAPER", ext=extension + "_ONLY"
    )
    print(output_path)

    print(outputmatrix)

    outputmatrix.to_csv(output_path, sep="\t", header=flanks)
    return

    print("get here")
    for site in ["8mer"]:
        site_flank_path = get_analysis_path(
            mirna, experiment, condition,
            "structural_analysis_PAPER_realfinal/%s" %(site), ext=extension
        )
        print(site_flank_path)
        # if not test:
        with open(site_flank_path, "wb") as file_flank:
            file_flank.write("\t".join(["Flank", "plfold", "accessibility", "AU_cs", "AU_win", "AU_read", "AU_win_wo", "AU_read_wo"])+"\n")
        with open(site_flank_path, "ab") as file_flank:
            output_rows = "\n".join([
                "\n".join([
                    "\t".join([flank]+[str(j) for j in i]) for i
                     in output[site][flank]
                 ]) for flank in flanks if len(output[site][flank]) > 0])
            file_flank.write(output_rows)

    # Print the amount of time the script took to complete.
    print_time_elapsed(time_start)

################################################################################

if __name__ == "__main__":
    main()

