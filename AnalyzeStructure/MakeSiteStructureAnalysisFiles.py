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
   return ''.join(random.choice(string.ascii_lowercase) for i in range(length))


# mirna: The miRNA in the experiment to be processed.
# experiment: The type of experiment to be processed.
# condition: The sampel type within the experiment to be processed.

# Functions from general this uses:


# Constants
AU_weights_left = [1.0 / i for i in range(30,0,-1)]
AU_weights_right = [1.0 / i for i in range(1,31)]

# FUNCTIONS

# AU_weights_right[0] = 0.5



def RNAplfold(read, mir_start, mir_stop, file_num):
    # This should give the probability of the region in the read complementary
    # to position "mir_start" to "mir_stop" in the miRNA sequence.


    read_name = randomword(20)
    temp_fa_filename = "%s%s.fa" %(read_name, file_num)

    with open(temp_fa_filename, "wb") as f:
        f.write((">%s%s\n%s" %(read_name, file_num, read)).encode())

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

def get_read_structural_data(read_seqs, _mirna, _sitelist, site, experiment,
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
    col_names = ["flank", "plfold", "accessibility", "AU_cs", "AU_win",
                 "AU_read", "AU_win_wo", "AU_read_wo"]

    sites_range = range(26 - n_constant, 26 + 37 + n_constant + 1)

    # site_flanks_map = {site: {"".join(kmer): [] for kmer 
    #     in list(it.product(["A","C","G","T"],repeat=4))} 
    #     for site in ["8mer"]}

    flank_struct_map = {
        flank: pd.DataFrame(0, index=[], columns=col_names)
        for flank in get_kmer_list(4)}

    # Pre-allocate output with "None" for each line.
    time_start = time.time()
    tick = 0
    num_8mers = 0
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
                ind_site = [i for i, site_i in enumerate(read_sites)
                             if site_i == site][0]
                start, stop = _read.site_pos(full=True)
                l_r_list = _read.all_sites_pos()
                # site_i = _read.sites[ind_site]
                mirp1_span = _sitelist.m1pos[site] - 1
                mirp1 = start + mirp1_span

                win_r = mirp1 - win_start + 1 + 1 
                win_l = mirp1 - win_stop + 1
                if "b" in site:
                    win_l -= 1
                plfold = RNAplfold(_read.fullseq, win_start, win_stop, num_i)
                df_cols = plfold.shape[1]
                # Get average secondary structure
                p_acc_pl = plfold.iloc[win_r-1, df_cols-2]**(1/(float(win_r - win_l)))
                p_acc_geo = 10**np.mean(
                    plfold.iloc[win_l:win_r, 0].apply(math.log10)
                )
                flank = "".join([_read.fullseq[i] for i in
                                 [start - 2, start - 1, stop, stop + 1]])
                AU_cs = get_local_au_score(read, start, stop, site)

                AU = [int(x in ["A", "T"]) for x in _read.fullseq]
                AU_win = AU[win_l:win_r]
                AU_read = AU[26:26 + 37]
                AU_win_wo = [AU[i_au] for i_au in range(win_l, win_r)
                                  if i_au not in range(start, stop)]
                AU_read_wo = [AU[i_au] for i_au in range(26, 26 + 37)
                                  if i_au not in range(start, stop)]
                # Normalize the AU strings such that they are between 0 and 1.
                AU_win, AU_read, AU_win_wo, AU_read_wo = [
                    float(sum(i))/len(i) for i in [
                        AU_win, AU_read, AU_win_wo, AU_read_wo]]
                # Add read info to flank content
                row = pd.DataFrame([flank, p_acc_pl, p_acc_geo, AU_cs, AU_win,
                                    AU_read, AU_win_wo, AU_read_wo],
                                   index=col_names).T
                # Add read info to flank content
                flank_struct_map[flank] = flank_struct_map[flank].append(row)

    return flank_struct_map





def main():
    time_start = time.time()

    # Define all the relevant arguments.
    arguments = ["miRNA", "experiment", "condition", "n_constant", "sitelist",
                 "site", "win_start", "win_stop", "-alt_mir_exp_cond",
                 "-buffer3p_binary", "-jobs", "-test_binary"]
    (
        mirna, experiment, condition, n_constant, sitelist, site, win_start,
        win_stop, alt_mir_exp_cond, buffer3p, 
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
        if alt_mir_exp_cond:
            str_split = alt_mir_exp_cond.split("_")
            print(str_split)
            if str_split[1] == "kin":
                str_split[1] = "_".join(str_split[1:3])
                str_split = str_split[:2] + str_split[3:]
            print(str_split)
            if str_split[-1] == "TGT":
                str_split[-2] = "_".join(str_split[-2:])
                str_split = str_split[:-1]
            print(str_split)

            mirna_use, experiment_use, condition_use = str_split
            input_list = [(mirna_use, experiment_use, condition_use)]
        else:
            input_list = [(mirna, experiment, condition)]
        if not jobs:
            jobs = 20

    _mirna = Mirna(mirna)
    _sitelist = SiteList(_mirna, sitelist,
                         int(read_length) + 2*int(n_constant))

    # Assign arguments to be passed to the function
    # "get_read_structural_data.""
    args  = [_mirna, _sitelist, site, experiment, int(n_constant), read_length, buffer3p,
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
    # output_matrix = {site: {flank: [j for i in threads
    #                          for j in i[0][site][flank]]
    #                          for flank in flanks}
    #                          for site in ["8mer"] if site != "None"}

    output = pd.concat(
        [pd.concat([thread[flank] for thread in threads])
         for flank in flanks])


    print(output)
    extension = "_%s_%s_%s-%s" %(n_constant, sitelist, win_start, win_stop)
    if alt_mir_exp_cond:
        extension = "%s_alt_%s" %(extension, alt_mir_exp_cond)
    if buffer3p:
        extension = "%s_buffer3p" %(extension)
    if test:
        extension = "%s_test" %(extension)
    # Write the output files. Here it is only written for the 8mer site type.
    out_path = get_analysis_path(
        mirna, experiment, condition,
        "structural_analysis/%s" %(site), ext=extension
    )
    print(out_path)
    output.to_csv(out_path, sep="\t", header=True, index=False)
    print_time_elapsed(time_start)

################################################################################

if __name__ == "__main__":
    main()

