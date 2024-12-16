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



def RNAplfold(read, mir_start, mir_stop, file_num):
    # This should give the probability of the region in the read complementary
    # to position "mir_start" to "mir_stop" in the miRNA sequence.


    read_name = randomword(20)
    temp_fa_filename = "%s%s.fa" %(read_name, file_num)

    with open(temp_fa_filename, "wb") as f:
        f.write(">%s%s\n%s" %(read_name, file_num, read))

    # call RNAplfold:
    l_param = len(read)
    w_param = l_param
    u_param = mir_stop - mir_start + 1

    call = 'RNAplfold -L %s -W %s -u %s < %s%s.fa' %(l_param, w_param, u_param,
                                                     read_name, file_num)
    subprocess.call([call], shell=True, stdout = subprocess.PIPE)

    temp_lunp_filename = '%s%s_lunp' %(read_name, file_num)
    rnaplfold_data = pd.read_csv(temp_lunp_filename, sep = '\t', header = 1).set_index(' #i$')
    os.remove(temp_fa_filename)
    os.remove(temp_lunp_filename)
    os.remove('%s%s_dp.ps' %(read_name, file_num))


    return rnaplfold_data


def calculate_local_au(utr, site_start, site_end, site_type):
    """
    Calculate the local AU score

    Parameters
    ----------
    utr: string, utr sequence

    site_type: string, site type

    site_start: int, start of site

    site_end: int, end of site

    Output
    ------
    float: local AU score
    """
    # find A, U and weights upstream of site
    up_site_adder = int(site_type not in ['7mer-m8', '8mer'])
    
    upstream = utr[max(0, site_start - 30): site_start]
    upstream_str = upstream

    l= max(0, site_start - 30)
    print(site_type)
    print(utr)
    upstream = [int(x in ['A', 'T']) for x in upstream]
    inv_upweights = [(x + 1 + up_site_adder)
                 for x in range(len(upstream))][::-1]

    upweights = [1.0 / (x + 1 + up_site_adder)
                 for x in range(len(upstream))][::-1]

    # find A,U and weights downstream of site
    down_site_adder = int(site_type in ['7mer-A1', '8mer'])
    downstream = utr[site_end:min(len(utr), site_end + 30)]
    downstream_str = downstream
    downstream = downstream + "A"*(30 - len(downstream))
    downstream = [int(x in ['A', 'T']) for x in downstream]

    inv_downweights = [(x + 1 + down_site_adder)
                   for x in range(len(downstream))]
    downweights = [1.0 / (x + 1 + down_site_adder)
                   for x in range(len(downstream))]
    print("_"*l + upstream_str + "."*(site_end - site_start) + downstream_str)
    print(" "*l + "".join([str(i) for i in upstream]) + "."*(site_end - site_start) + "".join([str(i) for i in downstream]))
    for num, weight in enumerate(inv_upweights):
        print(" "*(l + num) + str(weight) + " "*(30 + site_end - site_start - len(str(weight))) + str(inv_downweights[num]))

    weighted = np.dot(upstream, upweights) + np.dot(downstream, downweights)
    total = float(sum(upweights) + sum(downweights))

    return weighted / total




def get_local_au_score(read, start, stop, sitetype):
    if sitetype not in ["8mer", "7mer-m8", "7mer-A1", "6mer"]:
        return False
    else:
        upstream_site = int(sitetype not in ["7mer-m8", "8mer"])
        upstream_seq = read[max(0, start - 30): start]

        upstream = [int(x in ["A", "T"]) for x in upstream_seq]
        upweights = [1.0 / (x + 1 + upstream_site)
                     for x in range(len(upstream))][::-1]
        upweights_total = [1.0 / (x + 1 + upstream_site)
                     for x in range(30)][::-1]
        downstream_site = int(sitetype in ["7mer-A1", "8mer"])
        downstream_seq = read[stop:min(len(read), stop + 30)]
        downstream = [int(x in ["A", "T"]) for x in downstream_seq]
        downstream = downstream + [1 for i in range(30 - len(downstream))]
        downweights = [1.0 / (x + 1 + downstream_site)
                       for x in range(len(downstream))]
        downweights_total = [1.0 / (x + 1 + downstream_site)
                       for x in range(30)]

        weighted = np.dot(upstream, upweights) + np.dot(downstream, downweights)
        total = float(sum(upweights_total) + sum(downweights_total))

        return(weighted / total)
    return





# def assign_site(read_seqs, _mirna, _sitelist, experiment,
#                              n_constant, rand_length, buffer3p):
#     time_start = time.time()
#     counter = [0, 0, 0]
#     for i_r, r in enumerate(read_seqs):
#         # Get the read number:
#         read = r.strip()
#         # Make Read object.
#         _read = Read(read, rand_length, _mirna, n_constant, experiment)
#         sys.stdout.flush()
#         # Find all sites within the read
#         add_sites_from_read(_sitelist, _read, buffer_=buffer3p)
#     if _sitelist.__class__.__name__ in ["TwelvemerSiteList", "SixteenmerSiteList"]:
#         print("hi")
#         output = (_sitelist.top_counts_df(),
#                   _sitelist.multi_counts_df())
#     else:
#         output = (_sitelist.top_counts_df(),
#             _sitelist.multi_counts_df(),
#             _sitelist.single_pos_counts_df(),
#             _sitelist.top_pos_counts_df(),
#             _sitelist.single_pos_counts_df(barcode="ACA"),
#             _sitelist.top_pos_counts_df(barcode="ACA"))
#     return (output)




def get_read_structural_data(read_seqs, _mirna, _sitelist, experiment,
                             n_constant, rand_length, buffer3p, win_start,
                             win_stop, num_i):
    """Takes a read sequence and identifies the position of all site types'.

    Args:
        read_sequence: The read sequence.
        sites: A list of sequences corresponding to each site type.

    Returns:
        A tuple containing the list of reads and a dictionary of read types.
    """

    # Initialize count dictionary with keys as site names and a value of
    # for each item.
    sites_range = range(26 - n_constant, 26 + 37 + n_constant + 1)

    site_flanks_map = {site: {"".join(kmer): [] for kmer 
        in list(it.product(["A","C","G","T"],repeat=4))} 
        for site in ["8mer"]}

    # Pre-allocate output with "None" for each line.
    time_start = time.time()
    tick = 0
    num_8mers = 0
    max_stop = 0
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
                # print("------------------------")
                ind_8mer = [i for i, site in enumerate(read_sites)
                             if site == "8mer"][0]
                # print(ind_8mer)
                # print(_read.site_pos())
                # print(_read.topsite)
                (l, r)  = _read.site_pos()
                stop_rand = r - 5
                if stop_rand > max_stop:
                    max_stop = stop_rand
                (l_f, r_f) = _read.site_pos(full=True)
                l_r_list = _read.all_sites_pos()
                # print("l_r_list:")
                # print(l_r_list)
                # _read.print_sites(_sitelist)        
                # print(_read.fullseq)
                # for (l_f, r_f) in l_r_list:
                #     print(" "*l_f + _read.fullseq[l_f:r_f])

                site = _read.sites[ind_8mer]
                # CONVERT ALL COORDINATES TO 0 is 1, start of the read!
                # coords = [i.split(":")[1] for i in sites.split(", ")]
                # sites = [i.split(":")[0] for i in sites.split(", ")]
                # coord_site_map = {site: coord for site, coord in zip(sites,coords)}
                # ranks = [order_site_map[site] for site in sites]
                # [(coord, site)] = [i for i in zip(coords, sites)
                #                    if order_site_map[i[1]] == min(ranks)]

                # # Get the start and stop positions of the site, in pythonic
                # # coordinates.
                # start = int(coord.split("-")[0]) + 26 - n_constant
                # stop = int(coord.split("-")[1]) + 1 + 26 - n_constant
                # if start in sites_range and stop in sites_range and site in ["8mer"]:
                start, stop = (l_f, r_f)
                # print(start)
                # print(stop)
                mirp1_span = _sitelist.m1pos[site.name] - 1
                mirp1 = start + mirp1_span

                win_r = mirp1 - win_start + 1 + 1 
                win_l = mirp1 - win_stop + 1
                if "b" in site.name:
                    win_l -= 1
                plfold = RNAplfold(_read.fullseq, win_start, win_stop, num_i)
                df_cols = plfold.shape[1]
                # # print(plfold)
                # print(plfold.iloc[win_r - 1, ])
                # print("         1         2         3         4         5         6         7         8         9"[:len(_read.fullseq)])
                # print("123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890"[:len(_read.fullseq)])
                # print("".join([(str(i) + "   ")[0] for i in plfold.iloc[:, df_cols - 2]]))
                # print("".join([(str(i) + "   ")[1] for i in plfold.iloc[:, df_cols - 2]]))
                # print("".join([(str(i) + "   ")[2] for i in plfold.iloc[:, df_cols - 2]]))
                # print("".join([(str(i) + "   ")[3] for i in plfold.iloc[:, df_cols - 2]]))
                # print("".join([(str(i) + "   ")[4] for i in plfold.iloc[:, df_cols - 2]]))
                # Get average secondary structure
                p_acc_pl = plfold.iloc[win_r-1, df_cols-2]**(1/(float(win_r - win_l)))
                p_acc_geo = 10**np.mean(
                    plfold.iloc[win_l:win_r, 0].apply(math.log10)
                )
                l_1, l_2 = _read.fullseq[start-2], _read.fullseq[start-1]
                r_1, r_2 = _read.fullseq[stop], _read.fullseq[stop+1]
                flank = "".join([_read.fullseq[i] for i in
                                 [start - 2, start - 1, stop, stop + 1]])
                # print(flank)

                AU_cs = get_local_au_score(read, start, stop, site)

                AU = [int(x in ["A", "T"]) for x in _read.fullseq]
                # print(" "*l_f + _read.fullseq[l_f : r_f])
                # print(" "*win_l + _read.fullseq[win_l:win_r])
                # print(_read.fullseq)
                # print(" "*(start - 2)+ l_1 + l_2 + "_"*(stop - start) + r_1 + r_2)
                # print("".join([str(i) for i in AU]))
                AU_win = AU[win_l:win_r]
                AU_read = AU[26 : 26 + 37]
                AU_win_wo = [AU[i_au] for i_au in range(win_l, win_r)
                                  if i_au not in range(start, stop)]
                AU_read_wo = [AU[i_au] for i_au in range(26,26+37)
                                  if i_au not in range(start, stop)]
                # Normalize the AU strings such that they are between 0 and 1.
                AU_win, AU_read, AU_win_wo, AU_read_wo = [
                    float(sum(i))/len(i) for i in [
                        AU_win, AU_read, AU_win_wo, AU_read_wo]]
                # Add read info to flank content
                site_flanks_map[site.name][flank].append((p_acc_pl, p_acc_geo, AU_cs, AU_win, AU_read, AU_win_wo, AU_read_wo))
    return site_flanks_map, max_stop





def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["miRNA", "experiment", "condition", "n_constant", "sitelist",
                 "win_start", "win_stop", "-buffer3p_binary", "-jobs",
                 "-test_binary"]
    (
        mirna, experiment, condition,
        n_constant,sitelist, win_start,
        win_stop, buffer3p, jobs, test
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


    _mirna = Mirna(mirna)
    _sitelist = SiteList(_mirna, sitelist,
                         int(read_length) + 2*int(n_constant))

    # Assign arguments to be passed to the function
    # "get_read_structural_data.""
    args  = [_mirna, _sitelist, experiment,
             int(n_constant), read_length, buffer3p,
             int(win_start), int(win_stop)]

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
        threads += multiproc_file(reads_path, int(jobs),
                                  get_read_structural_data, test, *args)

    # Move back to the home directory and delete the folder just made
    os.chdir(cwd)
    os.rmdir(plfold_path)

    # Make a list of flanking dinucleotides
    flanks = get_kmer_list(4)
    output_matrix = {site: {flank: [j for i in threads
                             for j in i[0][site][flank]]
                             for flank in flanks}
                             for site in ["8mer"] if site != "None"}

    max_stops = [i[1] for i in threads]
    print(max_stops)

    # Name the output file extension.
    extension = "_%s_%s_%s-%s" %(n_constant, sitelist, win_start, win_stop)
    if buffer3p:
        extension = "%s_buffer3p" %(extension)
    if test:
        extension = "%s_test" %(extension)
        print(output_matrix)
    # Write the output files. Here it is only written for the 8mer site type.
    for site in ["8mer"]:
        output_path = get_analysis_path(
            mirna, experiment, condition,
            "structural_analysis_PAPER_realfinal/%s" %(site), ext=extension
        )
        with open(output_path, "wb") as file_flank:
            file_flank.write("\t".join(["Flank", "plfold", "accessibility",
                                        "AU_cs", "AU_win", "AU_read",
                                        "AU_win_wo", "AU_read_wo"]) + "\n")
        with open(output_path, "ab") as file_flank:
            output_str = "\n".join([
                "\n".join([
                    "\t".join([flank]+[str(j) for j in i]) for i
                     in output_matrix[site][flank]
                 ]) for flank in flanks if len(output_matrix[site][flank]) > 0])
            file_flank.write(output_str)
        print(output_path)

    # Print the amount of time the script took to complete.
    print_time_elapsed(time_start)

################################################################################

if __name__ == "__main__":
    main()

