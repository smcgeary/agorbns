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

# FUNCTIONS

AU_weights_left = [1.0 / i for i in range(30,0,-1)]
AU_weights_right = [1.0 / i for i in range(1,31)]
# AU_weights_right[0] = 0.5



# def RNAplfold(read, mir_start, mir_stop, file_num):
#     # This should give the probability of the region in the read complementary
#     # to position "mir_start" to "mir_stop" in the miRNA sequence.


#     read_name = randomword(20)
#     temp_fa_filename = "%s%s.fa" %(read_name, file_num)

#     with open(temp_fa_filename, "wb") as f:
#         f.write(">%s%s\n%s" %(read_name, file_num, read))

#     # call RNAplfold:
#     l_param = len(read)
#     w_param = l_param
#     u_param = mir_stop - mir_start + 1

#     call = 'RNAplfold -L %s -W %s -u %s < %s%s.fa' %(l_param, w_param, u_param,
#                                                      read_name, file_num)
#     subprocess.call([call], shell=True, stdout = subprocess.PIPE)

#     temp_lunp_filename = '%s%s_lunp' %(read_name, file_num)
#     rnaplfold_data = pd.read_csv(temp_lunp_filename, sep = '\t', header = 1).set_index(' #i$')
#     os.remove(temp_fa_filename)
#     os.remove(temp_lunp_filename)
#     os.remove('%s%s_dp.ps' %(read_name, file_num))


#     return rnaplfold_data


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
    # print(read)
    for row in range(rnaplfold_data.shape[0]):
        rnaplfold_data.iloc[row, ] = (rnaplfold_data.iloc[row, ].shift(-((row + 1)/2)))**(1/float(row + 1))
    # print(" "*(mir1p - 8 + 1) + read[mir1p - 8 + 1: mir1p + 1])
    # print("".join([(str(i) + "   ")[0] for i in rnaplfold_data.iloc[0, :]]))
    # print("".join([(str(i) + "   ")[1] for i in rnaplfold_data.iloc[0, :]]))
    # print("".join([(str(i) + "   ")[2] for i in rnaplfold_data.iloc[0, :]]))
    # print("".join([(str(i) + "   ")[3] for i in rnaplfold_data.iloc[0, :]]))
    # print("".join([(str(i) + "   ")[4] for i in rnaplfold_data.iloc[0, :]]))
    # print("".join([(str(i) + "   ")[5] for i in rnaplfold_data.iloc[0, :]]))

    # print(" "*(mir1p - 8 + 1) + read[mir1p - 8 + 1: mir1p + 1])
    # print("".join([(str(i) + "   ")[0] for i in rnaplfold_data.iloc[1, :]]))
    # print("".join([(str(i) + "   ")[1] for i in rnaplfold_data.iloc[1, :]]))
    # print("".join([(str(i) + "   ")[2] for i in rnaplfold_data.iloc[1, :]]))
    # print("".join([(str(i) + "   ")[3] for i in rnaplfold_data.iloc[1, :]]))
    # print("".join([(str(i) + "   ")[4] for i in rnaplfold_data.iloc[1, :]]))
    # print("".join([(str(i) + "   ")[5] for i in rnaplfold_data.iloc[1, :]]))

    # print(" "*(mir1p - 8 + 1) + read[mir1p - 8 + 1: mir1p + 1])
    # print("".join([(str(i) + "   ")[0] for i in rnaplfold_data.iloc[2, :]]))
    # print("".join([(str(i) + "   ")[1] for i in rnaplfold_data.iloc[2, :]]))
    # print("".join([(str(i) + "   ")[2] for i in rnaplfold_data.iloc[2, :]]))
    # print("".join([(str(i) + "   ")[3] for i in rnaplfold_data.iloc[2, :]]))
    # print("".join([(str(i) + "   ")[4] for i in rnaplfold_data.iloc[2, :]]))
    # print("".join([(str(i) + "   ")[5] for i in rnaplfold_data.iloc[2, :]]))
    os.remove(temp_fa_filename)
    os.remove(temp_plfold_filename)
    os.remove('%s%s_dp.ps' %(read_name, file_num))
    return rnaplfold_data


def get_plfold_matrix(read_seqs, _mirna, _sitelist,
                       experiment, n_constant, rand_length,
                       buffer3p, lflank, rflank, win, num_i):
    """Takes a read sequence and identifies the position of all site types'.

    Args:
        read_sequence: The read sequence.
        sites: A list of sequences corresponding to each site type.

    Returns:
        A tuple containing the list of reads and a dictionary of read types.
    """

    # Initialize count dictionary with keys as site names and a value of
    # for each item.
    # sites_range = range(26 - n_constant, 26 + 37 + n_constant + 1)

    # site_flanks_map = {site: {"".join(kmer): [] for kmer 
    #     in list(it.product(["A","C","G","T"],repeat=4))} 
    #     for site in ["8mer"]}

    # pl_pos_5p = ["f%s" %(i) for i in range(1, win_stop - 8 + 1)]
    # pl_pos_3p = ["t%s" %(i) for i in range(1, -1*win_start + 1)]
    # print(pl_pos_3p)
    # pl_pos_seed = ["s%s" %(i) for i in range(1, 8 + 1)[::-1]]
    # pl_pos_keys = pl_pos_5p + pl_pos_seed + pl_pos_3p
    # pl_win_keys = ["w%s" %(i) for i in range(1, win+1)]
    # pos_win_plfold_map = {i: {j : [] for j in pl_win_keys} for i in pl_pos_keys}
    # pos_win_logplfold_map = {i: {j : [] for j in pl_win_keys} for i in pl_pos_keys}
    # flanks = []
    # time_start = time.time()
    # tick = 0


    flanks = get_kmer_list(4)

    df_columns = (["f5p%s" %(i) for i in range(lflank, 0, -1)] +
                ["sp%s" %(i) for i in range(8, 0, -1)] + 
                ["f3p%s" %(i) for i in range(1, rflank + 1)])

    df_index = [str(i + 1) for i in range(win)]

    # test rows:
    # df_columns = (["f5p%s" %(i) for i in range(2, 0, -1)] +
    #             ["sp%s" %(i) for i in range(8, 0, -1)] + 
    #             ["f3p%s" %(i) for i in range(1, 3)])

    plmean_flanks_map = {
        flank: pd.DataFrame(0.0, index=df_index, columns=df_columns)
        for flank in flanks
    }
    n_flanks_map = {
        flank: pd.DataFrame(0, index=df_index, columns=df_columns)
        for flank in flanks
    }


    # print(plmean_flanks_map["AAAA"])
    # print(n_flanks_map["AAAA"])
    # Pre-allocate output with "None" for each line.
    time_start = time.time()
    tick = 0
    num_8mers = 0
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
            if "8mer" in read_sites:
                ind_8mer = [i for i, site in enumerate(read_sites)
                             if site == "8mer"][0]
                (l_f, r_f) = _read.site_pos(full=True)
                l_r_list = _read.all_sites_pos()
                # _read.print_sites(_sitelist)        
                site = _read.sites[ind_8mer]
                start, stop = (l_f, r_f)
                if stop > stop_max:
                    stop_max = stop
                mirp1_span = _sitelist.m1pos[site.name] - 1
                mirp1 = start + mirp1_span

                # win_r = mirp1 - win_start + 1 + 1 
                # win_l = mirp1 - win_stop + 1
                # if "b" in site.name:
                #     win_l -= 1

                # Perform plfold on the entire read:
                headers = n_flanks_map["AAAA"].columns
                # print(list(headers))
                # print(" "*(mirp1 - 8 + 1 - lflank) + "".join([i[0] for i in headers]))

                plfold = RNAplfold(_read.fullseq, win, mirp1, num_i)

                l_pl = mirp1 - 8 + 1 - lflank
                r_pl = mirp1 + rflank + 1
                # print(l_pl)
                # print(r_pl)
                plfold_subset = plfold.iloc[:, l_pl:r_pl]
                # print(" "*(mirp1 - 8 + 1 - lflank) + "".join([i[0] for i in headers]))

                # print(" "*(l_pl) + "".join([str(i)[0] for i in plfold_subset.iloc[0, :]]))
                # print(" "*(l_pl) + "".join([str(i)[1] for i in plfold_subset.iloc[0, :]]))
                # print(" "*(l_pl) + "".join([str(i)[2] for i in plfold_subset.iloc[0, :]]))
                # print(" "*(l_pl) + "".join([str(i)[3] for i in plfold_subset.iloc[0, :]]))
                # print(" "*(l_pl) + "".join([str(i)[4] for i in plfold_subset.iloc[0, :]]))

                # print(" "*(l_pl) + "".join([(str(i) + "  ")[0] for i in plfold_subset.iloc[1, :]]))
                # print(" "*(l_pl) + "".join([(str(i) + "  ")[1] for i in plfold_subset.iloc[1, :]]))
                # print(" "*(l_pl) + "".join([(str(i) + "  ")[2] for i in plfold_subset.iloc[1, :]]))
                # print(" "*(l_pl) + "".join([(str(i) + "  ")[3] for i in plfold_subset.iloc[1, :]]))
                # print(" "*(l_pl) + "".join([(str(i) + "  ")[4] for i in plfold_subset.iloc[1, :]]))


                plfold_subset.index = df_index
                plfold_subset.columns = df_columns
                # print(plfold_subset)
                # print(plfold_subset.notnull())
                flank = "".join([_read.fullseq[i] for i in
                 [start - 2, start - 1, stop, stop + 1]])
                # print("flank")
                # print(flank)
                # update old_matrix
                # print(plmean_flanks_map[flank])
                for row in range(plfold_subset.shape[0]):
                    # print(row)
                    # print(plfold_subset.iloc[row, :].notnull())
                    col_inds = plfold_subset.iloc[row, :].notnull()
                    mean_old = plmean_flanks_map[flank].iloc[row, :][col_inds]
                    n_old = n_flanks_map[flank].iloc[row, :][col_inds]
                    # print("original")
                    # print(mean_old)
                    # print("n")
                    # print(n_old)
                    mean_i = np.log(plfold_subset.iloc[row, :][col_inds])
                    # print("new")
                    # print(mean_i)
                    n_tots = n_old + 1
                    # print("n new")
                    # print(n_tots)

                    mean_new = (mean_old*n_old + mean_i)/(n_tots)
                    # print("mean_new")
                    # print(mean_new)

                    plmean_flanks_map[flank].iloc[row, :][col_inds] = mean_new
                    n_flanks_map[flank].iloc[row, :][col_inds] = n_tots

                    # print(plmean_flanks_map[flank])
                    # print(n_flanks_map[flank])


                    # row_num = (
                    #     np.log(plmean_flanks_map[flank].iloc[row, :][col_inds]*n_flanks_map[flank].iloc[row, :][col_inds]) +
                    #     np.log(plfold_subset.iloc[row, :][col_inds])
                    # )
                    # print("row_num")
                    # print(row_num)
                    # print(np.log(plmean_flanks_map[flank].iloc[row, :][col_inds]))
                    # row_new = np.exp(
                    #     (
                    #         np.log(plmean_flanks_map[flank].iloc[row, :][col_inds]*n_flanks_map[flank].iloc[row, :][col_inds]) +
                    #         np.log(plfold_subset.iloc[row, :][col_inds])
                    #     ) / (
                    #         (n_flanks_map[flank].iloc[row, :][col_inds] + 1)
                    #     )
                    # )
                    # print(row_new)
                # print(plmean_flanks_map[flank])


                # struc = RNA.fold(read)[0]
                # pos_zips = zip(pl_pos_keys, full_inds)
                # print(pos_zips)
                # # print(plfold.shape)
                # # print(" "*win_l + "."*(win_r - win_l))
                # # for i, pos in enumerate(full_inds):
                # #     if pl_pos_keys[i] == "seed8":
                # #         print(" "*start + read[start:stop])
                # #     if pl_pos_keys[i] == "threep1":
                # #         print(" "*start + read[start:stop])
                # #     print(" "*pos + pl_pos_keys[i])
                # # print(" "*win_l + "."*(win_r - win_l))
                # for (i_c, win_key) in enumerate(pl_win_keys):
                #     first_position = pos_zips[0]
                #     shift = i_c/2
                #     i_r_list = [i + i_c/2 for i in full_inds]
                #     left_nas = [float('nan')]*max(0, -1*min(i_r_list))
                #     pl_left_bound = max(0, min(i_r_list))
                #     pl_right_bound = min(plfold.shape[0], max(i_r_list) + 1)
                #     # print(full_inds)
                #     # print(pl_left_bound)
                #     # print(pl_right_bound)
                #     # print(plfold.shape)
                #     pl_row = list(plfold.iloc[pl_left_bound:pl_right_bound,i_c])
                #     logpl_row = list(logplfold.iloc[pl_left_bound:pl_right_bound,i_c])
                #     right_nas = [float('nan')]*max(0, max(i_r_list)-plfold.shape[0]+1)
                #     output_row = left_nas + pl_row + right_nas
                #     logoutput_row = left_nas + logpl_row + right_nas
                #     # print(len(left_nas))
                #     # print(len(right_nas))
                #     # if flank == "CC.AG":
                #     #     print(win_key)
                #     #     print(len(output_row))
                #     # print(len(pl_pos_keys))
                #     for (pos_key, out, logout) in zip(pl_pos_keys, output_row, logoutput_row):
                #         # if win_key == "w1":
                #         #     print("%s\t%s" %(pos_key, out))
                #         # print(i_r)
                #         # print(plfold.shape)
                #         # print(" "*(ind + shift) + "*" + " "*(len(read)-ind - shift -1) + pos)
                #         pos_win_plfold_map[pos_key][win_key].append(out)
                #         pos_win_logplfold_map[pos_key][win_key].append(logout)
    return [plmean_flanks_map, n_flanks_map]




def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["miRNA", "experiment", "condition",
                 "n_constant", "sitelist", "-buffer3p_binary",
                 "-win", "-lflank", "-rflank",
                 "-jobs", "-test_binary"]

    # Parse arguments and assign script parameters accordingly:
    (
        mirna, experiment, condition,
        n_constant, sitelist, buffer3p,
        win, lflank, rflank,
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

    if lflank:
        lflank = int(lflank)
    else:
        lflank = 16
    if rflank:
        rflank = int(rflank)
    else:
        rflank = 16
    if win:
        win = int(win)
    else:
        win = 30

    # Call instance variables:
    _mirna = Mirna(mirna)
    _sitelist = SiteList(_mirna, sitelist,
                         int(read_length) + 2*int(n_constant))

    # Assign arguments to be passed to the function
    # "get_plfold_matrix.""
    args  = [_mirna, _sitelist, experiment,
             int(n_constant), read_length, buffer3p,
             lflank, rflank, win]

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


    df_columns = (["f5p%s" %(i) for i in range(lflank, 0, -1)] +
                ["sp%s" %(i) for i in range(8, 0, -1)] + 
                ["f3p%s" %(i) for i in range(1, rflank + 1)])

    df_index = [str(i + 1) for i in range(win)]

    plmean_flanks_map_final = {
        flank: pd.DataFrame(0.0, index=df_index, columns=df_columns)
        for flank in flanks
    }
    n_flanks_map_final = {
        flank: pd.DataFrame(0, index=df_index, columns=df_columns)
        for flank in flanks
    }
    # print(n_flanks_map_final["AAAA"])
    # print(plmean_flanks_map_final["AAAA"])
    for flank in flanks:
        # print(flank)
        for row in range(threads[0][0][flank].shape[0]):
            for thread in threads:
                row_ijk = thread[1][flank].iloc[row, :]
                # print(row_ijk)
                col_inds = (row_ijk != 0)
                # if flank == "TTTC":
                #     print(col_inds)
                #     print(row_ijk)
            # print(row)
            # print(plfold_subset.iloc[row, :].notnull())
                mean_old = plmean_flanks_map_final[flank].iloc[row, :][col_inds]
                n_old = n_flanks_map_final[flank].iloc[row, :][col_inds]
                # print("original")
                # print(mean_old)
                # print("n")
                # print(n_old)
                mean_i = thread[0][flank].iloc[row, :][col_inds]
                # print("new")
                # print(mean_i)
                n_tots = n_old + 1
                # print("n new")
                # print(n_tots)

                mean_new = (mean_old*n_old + mean_i)/(n_tots)
                # print("mean_new")
                # print(mean_new)

                plmean_flanks_map_final[flank].iloc[row, :][col_inds] = mean_new
                n_flanks_map_final[flank].iloc[row, :][col_inds] = n_tots

                # print(plmean_flanks_map_final[flank])
                # print(n_flanks_map_final[flank])


    print(plmean_flanks_map_final["TTTC"])
    print(n_flanks_map_final["TTTC"])


    df_final_index = ["%s_w%s" %(i[1], i[0])
                      for i in it.product(list(plmean_flanks_map_final["AAAA"].index),
                                   list(plmean_flanks_map_final["AAAA"].columns))]
    df_final_columns = flanks



    outputmatrix = pd.DataFrame(0.0, index=df_final_index, columns=df_final_columns)

    for flank in flanks:
        for df_ind in df_index:
            for df_col in df_columns:
                df_final_ind = "%s_w%s" %(df_col, df_ind)
                outputmatrix.loc[df_final_ind, flank] = plmean_flanks_map_final[flank].loc[df_ind, df_col]

    print(outputmatrix)




    # Name the output file extension.
    extension = "_%s_%s_%s_l%s_r%s_w%s" %(n_constant, sitelist, "8mer", lflank, rflank, win)
    if buffer3p:
        extension = "%s_buffer3p" %(extension)
    if test:
        extension = "%s_test" %(extension)

    output_path = get_analysis_path(
        mirna, experiment, condition,
            "plfold_2018_PAPER", ext=extension
    )
    print(output_path)



    outputmatrix.to_csv(output_path, sep="\t", header=flanks)
    return

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

