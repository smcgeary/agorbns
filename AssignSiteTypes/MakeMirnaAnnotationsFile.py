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
from collections import defaultdict

# FUNCTIONS
spike_seqs = [marker_18nt_d,
              marker_30nt_d[:18],
              adapter_5p_d,
              adapter_3p_sRS,
              dme_miR_14_5p_d,
              xtr_miR_427_d,
              "TGGAATGTAAAGAAGTATGTAT",
              "TTAATGCTAATCGTGATAGGGGT"]
spike_names = ["18nt_marker", "30nt_marker", "5p_adapter",
                               "3p_adapter", "dme-miR-14-5p",
                               "xtr-miR-427", "miR-1", "miR-155"]






# def assign_mirna(read_seqs):
#     time_start = time.time()
#     counter = [0, 0, 0]
#     mirnas = defaultdict(list)
#     spike_dict = defaultdict(list)
#     for spike_seq in spike_seqs:
#         spike_dict[spike_seq] = {i : 0 for i in range(40)}

#     spike_dict["None"] = {i : 0 for i in range(40)}
#     unassigned = defaultdict(int)
#     for i_r, r in enumerate(read_seqs):
#         # Get the read number:
#         seq = r.strip()
#         # print(seq)

#         adapter_pos = GetThreePrimeBarcodeOverlap(seq)
#         # print(adapter_pos)
#         barcode =  seq[:14]
#         mirna = seq[14:adapter_pos]
#         mirnas[mirna].append(barcode)
#         assigned = False
#         for spike_seq in spike_seqs:
#             spike_pos = seq.find(spike_seq)
#             i_alt = test_match(spike_seq,seq,mis=4)
#             # print(i_alt)
#             # print(spike_pos)
#             if i_alt:
#                 spike_dict[spike_seq][i_alt] += 1
#                 assigned = True
#                 break

#         if not assigned:
#             spike_dict["None"][0] += 1
#             unassigned[seq[12:]] += 1


#         # print(seq[14:adapter_pos])
#         # # Make Read object.
#         # _read = Read(read, rand_length, _mirna, n_constant, experiment)
#         # _read.get_all_sites_in_read(_sitelist)
#         # site_names = [_site.name for _site in _read.sites]
#         sys.stdout.flush()
#         # # Find all sites within the read
#         # add_sites_from_read(_sitelist, _read)
#     return (mirnas, spike_dict, unassigned)

def main():
    file_path = "AssignSiteTypes/mature.fa"

    with open(file_path, "r+") as file_in:
        lines = list(file_in)
        lines = lines[:2000]
        for line in izip(*[iter(lines)]*2):
            name, seq = line
            mirna = name.split(" ")[0][1:]
            seq = re.sub("U", "T", seq.strip())
            species = mirna.split("-")[0]
            if species == "hsa":
                print("%s\t%s" %(mirna, seq))

    return
    # with 
    # time_start = time.time()
    # arguments = ["miRNA", "experiment", "condition", "-jobs", "-test_binary"]
    # args = parse_arguments(arguments)
    # mirna, experiment, condition, jobs, test = args



    # input_list = [(mirna, experiment, condition)]
    # if not jobs:
    #     jobs = 20

    # ############################################################################
    # threads = []
    # # For each file in the input list, perform the multiprocessing:
    # reads_path = get_analysis_path(mirna, experiment, condition, "reads")
    # if not jobs:
    #     jobs = 1
    # threads = multiproc_file(reads_path, int(jobs), assign_mirna, test)

    # threads_summed = defaultdict(list)
    # spikes_summed = defaultdict(list)
    # unassigned_all = defaultdict(int)

    # spike_dict_all = defaultdict(list)
    # for spike_seq in spike_seqs:
    #     spike_dict_all[spike_seq] = {i : 0 for i in range(40)}
    # spike_dict_all["None"] = {i : 0 for i in range(40)}

    # print(threads[0][1])
    # for thread in threads:
    #     for key in thread[0]:
    #         threads_summed[key] += thread[0][key]
    #     for key in thread[1]:
    #         print(key)
    #         print(thread[1][key])
    #         for pos in range(40):
    #             spike_dict_all[key][pos] += thread[1][key][pos]
    #     for key in thread[2]:
    #         unassigned_all[key] += thread[2][key]

    # seq_spike_dict = {i[0] : i[1] for i in zip(spike_names, spike_seqs)}
    # spike_seq_dict = {i[1] : i[0] for i in zip(spike_names, spike_seqs)}
    # spike_seq_dict["None"] = "None"


    # unassigned_all_df =pd.DataFrame.from_dict(unassigned_all, orient="index")
    # unassigned_all_df.sort_values([0], ascending=[False], inplace=True)

    # print(unassigned_all_df.iloc[:10,])
    # print(unassigned_all_df.sum())

    # # return

    # spike_dict_all_pd = pd.DataFrame.from_dict(spike_dict_all)
    # spike_dict_all_pd.columns = [spike_seq_dict[i] for i in list(spike_dict_all_pd.columns.values)]
    # print(spike_dict_all_pd)

    # table_path = get_analysis_path(mirna, experiment, condition,
    #                                  "AGO_pur_counts")
    # unassigned_path = get_analysis_path(mirna, experiment, condition,
    #                                  "AGO_pur_counts", ext="_unassigned_reads")


    # # print(unassigned_reads[:10])


    # print(table_path)
    # print(unassigned_path)
    # if not test:
    #     unassigned_all_df.to_csv(unassigned_path, sep="\t")
    #     spike_dict_all_pd.transpose().to_csv(table_path, sep="\t")

    # return

    # seq_spike_dict = {i[0] : i[1] for i in zip(spike_names, spike_seqs)}

    # spikes_dist_paths = [get_analysis_path(mirna, experiment, condition,
    #                                      "AGO_pur_counts",
    #                                      ext="_%s_positions" %(i))
    #                      for i in ["18nt_marker", "30nt_marker", "5p_adapter",
    #                                "3p_adapter", "dme-miR-14-5p",
    #                                "xtr-miR-427"]]
    # spike_paths_dict = {i[0] : i[1]
    #                     for i in zip(spike_names, spikes_dist_paths)}
    # for spike_name in spike_names:
    #     positions = spikes_summed[seq_spike_dict[spike_name]]
    #     print(positions)
    #     positions_text = ",".join([str(i) for i  in positions])
    #     print(positions_text)
    #     out_file = spike_paths_dict[spike_name]
    #     print(positions)
    #     print(out_file)
    #     if not test:
    #         with open(out_file, "w+") as file_out:
    #             file_out.write(positions_text)


    # print(spikes_dist_paths)

    # # print(threads_summed)
    # # for key, val in threads_summed.items():
    # #     print(key)
    # #     print(val[0])

    # seqs_df = pd.DataFrame([(key, len(val[0]), len(list(set(val[0]))))
    #                         for key, val in threads_summed.items()]).set_index([0])
    # seqs_df.columns = ["total", "collapsed"]
    # seqs_df.sort_values(["total", "collapsed"], ascending=[False, False], inplace=True)
    # seqs_df.index.name = None
    # print(seqs_df.iloc[:10,:])
    # seq_counts_path = get_analysis_path(mirna, experiment, condition, "AGO_pur_counts")
    # if not test:
    #     seqs_df.to_csv(seq_counts_path, sep="\t", header=False)
    # # return
    # # # The summed data from each of the multiprocessing threads:
    # # counts = merge_data_frames([i[0] for i in threads])["TGT"]
    # # multicounts = merge_data_frames([i[1] for i in threads])["TGT"]
    # # if sitelist not in ["12mers", "16mers"]:
    # #     single_pos_counts = merge_data_frames([i[2] for i in threads])
    # #     top_pos_counts = merge_data_frames([i[3] for i in threads])

    # # # The names of the output files:
    # # site_counts_path = get_analysis_path(mirna, experiment, condition,
    # #                                    "site_counts", ext=extension)
    # # multisite_counts_path = get_analysis_path(mirna, experiment, condition,
    # #                                         "multisite_counts",
    # #                                         ext=extension)
    # # single_pos_counts_path = get_analysis_path(mirna, experiment, condition,
    # #                                         "single_pos_counts",
    # #                                         ext=extension)
    # # top_pos_counts_path = get_analysis_path(mirna, experiment, condition,
    # #                                         "top_pos_counts",
    # #                                         ext=extension)
    # # outputpaths = [get_analysis_path(mirna, experiment, condition, i,
    # #                                  ext=extension)
    # #                for i in ["site_counts", "multisite_counts",
    # #                          "single_pos_counts", "top_pos_counts"]]

    # # print(outputpaths[0])

    # # print("done reading files")

    # # if test:
    # #     print(counts)
    # # if not test:
    # #     counts.to_csv(site_counts_path, sep="\t", header=False)
    # #     multicounts.to_csv(multisite_counts_path, sep="\t", header=False)
    # #     if sitelist not in ["12mers", "16mers"]:
    # #         single_pos_counts.to_csv(single_pos_counts_path, sep="\t", header=positional_names)
    # #         top_pos_counts.to_csv(top_pos_counts_path, sep="\t", header=positional_names)

    # print_time_elapsed(time_start)
################################################################################

if __name__ == "__main__":
    main()

