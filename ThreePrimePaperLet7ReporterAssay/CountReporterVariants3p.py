################################################################################
#LibraryDesign.py
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
from operator import add

pd.set_option('max_columns', 600)
pd.set_option("max_colwidth", 600)
# FUNCTIONS

MIRNAS = ["miR-1", "let-7a", "miR-155", "miR-124", "lsy-6", "miR-7"]


#FINAL DESIGN CONSTRAINTS
kLEN = 180

# def get_mirna_site_order(mirna_var):
#     # Fill dataframe with miRNA, name, and sequence of each site.
#     variants_path = "oligo_df_210726_library0.txt"    
#     # Pre-allocate the dictionary of miRNA sites:
#     output = []
#     with open(variants_path, "r") as file_in:
#         header = file_in.readline().strip()
#         while header:
#             mirna, site, index = header.split("> ")[1].split("_")
#             sequence = file_in.readline().strip()
#             # If site not yet in the subdictionary, add it:
#             if mirna == mirna_var and site not in output:
#                 output += [site]
#             header = file_in.readline().strip()
#     return output



# def get_site_sequence_map_dictionary():
#     # Fill dataframe with miRNA, name, and sequence of each site.
#     variants_path = "ThreePrimePaperLet7ReporterAssay/libraries_210726/oligo_df_210726_library0.txt"    
#     # Pre-allocate the dictionary of miRNA sites:
#     return()
#     with open(variants_path, "r") as file_in:
#         header = file_in.readline().strip()
#         print(header)
#         return()
#         while header:
#             mirna, site, index = header.split("> ")[1].split("_")
#             sequence = file_in.readline().strip()
#             # If site not yet in the subdictionary, add it:
#             if site not in mirna_site_map[mirna]:
#                 mirna_site_map[mirna][site] = [sequence]
#             # Otherwise, add the sequence.
#             else:
#                 mirna_site_map[mirna][site] += [sequence]
#             header = file_in.readline().strip()
#     return mirna_site_map
#     # Now this is a dictionary of dictionaries, where each ultimate value is
#     # a list of the 184 variants ascribed to that site.
#     # Make a new dictionary in which the key is the first 100 nt of the variant,
#     # and the value is a tuple with the mirna, site, the variant index, and the
#     # the full sequence.


def get_variant_dictionary(read_len=120):
    # Read in the variants file as a pandas DataFrame.
    var_path = "ThreePrimePaperLet7ReporterAssay/libraries_210726/oligo_df_210726_library0.txt"    
    var_df = pd.read_csv(var_path, sep="\t", index_col=0)
    # Add a new column to the dataframe .
    var_df["Read_seq"] = var_df["Sequence"].apply(
        lambda x: x[17:].split("TACCAATGCCCTGGCTC")[0][:read_len]
    )
    # Define the columns to use in the final dataframe, and make a 0-row
    # dataframe in which the non-redundant rows will be modified and stored.
    columns_final = ["Seed", "ThreePrime", "Bulge",
                       "Location", "Context", "Read_seq"]
    var_df_uniq = pd.DataFrame(columns=columns_final)
    # Definition of function to be used for the string collapsation of the rows
    # with multiple pieces of information.
    def collapse_function(x):
        if x.isna().all():
            return(np.nan)
        elif x.isna().any():
            return("|").join(sorted([str(i_x) for i_x in x]))
        else:
            return("|".join(sorted(list(set(x)))))
    # Iterate over the rows of the unique Read sequences within the dataframe
    # to collapse the variants.
    for i, read_i in enumerate(list(set(var_df["Read_seq"]))):
        var_df_i = var_df[var_df["Read_seq"] == read_i][columns_final]
        if var_df_i.shape[0] > 1:
            var_df_uniq_i = var_df_i.apply(collapse_function, axis=0)
            var_df_uniq_i.name = var_df_i.index[0]
            # print_round = True
        else:
            var_df_uniq_i = var_df_i
            # print_round = False
        # var_df_uniq_i.name = i
        var_df_uniq = var_df_uniq.append(var_df_uniq_i)
    # Sort the dataframe based on the (now non-consecutive) indeces, add a
    # column for the counts, and rename the row indeces according to the
    # expected sequencing reads.
    var_df_uniq.sort_index(inplace=True)
    # Add a row for the unmapped data.
    var_df_unmapped = pd.Series(
        [np.nan]*(len(columns_final) - 1) + ["Unmapped"], index=columns_final,
        name = "Unmapped"
    )
    var_df_uniq = var_df_uniq.append(var_df_unmapped)
    var_df_uniq["Counts"] = [0]*var_df_uniq.shape[0]
    var_df_uniq.index = var_df_uniq["Read_seq"]
    return(var_df_uniq)

# def make_variant_count_dictionary():
#     output = {mirna : {} for mirna in MIRNAS}
#     for mirna in MIRNAS:
#         for site in mirna_site_map[mirna]:
#             output[mirna][site] = [0]*184
#     return(output)

def assign_variant(read_seqs, var_df_uniq):
    time_start = time.time()
    # print(read_seqs[:10])
    # print(var_df_uniq)
    # sys.stdout.flush()
    unmapped = 0
    print(len(read_seqs))
    for i_r, r in enumerate(read_seqs):
        # Get the read number:
        r = r.strip()
        if r not in var_df_uniq.index:
            r = "Unmapped"
        var_df_uniq.loc[r, "Counts"] += 1
    print(var_df_uniq[var_df_uniq["Counts"] > 0])
    return(var_df_uniq)

# def assign_variant_mm(read_seqs, mirna_site_map, seq100nt_variant_map):
#     time_start = time.time()
#     output_dic = make_variant_count_dictionary(mirna_site_map)
#     unmapped = 0
#     tot_mapped = 0
#     mirnas = ["miR-1", "let-7a", "miR-155", "miR-124", "lsy-6", "miR-7-23nt"]
#     sitelist = "paperfinal"
#     seq25nt_left = {key[12:37] : item
#                     for key, item in seq100nt_variant_map.items()}
#     seq25nt_mid = {key[37:62] : item
#                     for key, item in seq100nt_variant_map.items()}
#     seq25nt_right = {key[62:87] : item
#                     for key, item in seq100nt_variant_map.items()}

#     for i_r, r in enumerate(read_seqs):
#         if i_r % 100000 == 0:
#             print(i_r)
#         # Get the read number:
#         r = r.strip()
#         mapped = False
#         sys.stdout.flush()
#         if r in seq100nt_variant_map:
#             mirna, site, i, sequence = seq100nt_variant_map[r]
#             output_dic[mirna][site][int(i)] += 1
#             mapped = True
#             tot_mapped += 1
#         else:
#             sys.stdout.flush()
#             for i in range(len(r) - 25 + 1):
#                 r_i = r[i:i + 25]
#                 map_ = False
#                 if r_i in seq25nt_left:
#                     map_ = seq25nt_left[r_i]
#                 elif r_i in seq25nt_mid:
#                     map_ = seq25nt_mid[r_i]
#                 elif r_i in seq25nt_right:
#                     map_ = seq25nt_right[r_i]
#                 sys.stdout.flush()
#                 if map_:
#                     mirna, site, i, sequence = map_
#                     _mirna = Mirna(mirna)
#                     _sitelist = SiteList(_mirna, "paperfinal", 100)
#                     site_seq = _mirna[site]
#                     _variant = ReporterVariant(r)
#                     _variant.get_all_sites_in_reporter_variant(_sitelist)
#                     ts = _variant.topsite
#                     sys.stdout.flush()
#                     site_bool = _variant.topsite
#                     other_topsites = []
#                     other_toppos = []
#                     for mirna_i in MIRNAS:
#                         if mirna_i != mirna:
#                             _mirna = Mirna(mirna_i)
#                             _sitelist = SiteList(_mirna, "paperfinal", 100)
#                             _variant = ReporterVariant(r)
#                             _variant.get_all_sites_in_reporter_variant(_sitelist)
#                             ts_other = _variant.topsite.name
#                             ts_pos = _variant.topsite.l
#                             other_topsites.append(ts_other)
#                             other_toppos.append(ts_pos)
#                     # print(other_topsites)
#                     # print(other_toppos)
#                     unique_other_topsites = list(set(other_topsites))
#                     # print(unique_other_topsites)
#                     # print(ts.name)
#                     # print(mirna)
#                     # print(unique_other_topsites == ["None"])
#                     sys.stdout.flush()
#                     if ts.name == site and unique_other_topsites == ["None"]:
#                         # var = sequence[17:117]
#                         # mismatches = [int(i == j) for i, j in zip(list(r), list(var))]
#                         output_dic[mirna][site][int(i)] += 1
#                         mapped = True
#                         tot_mapped += 1
#                         sys.stdout.flush()
#                         break
#         if not mapped:
#             unmapped += 1
#     return([output_dic, unmapped])





def main():
    time_start = time.time()
    # Parse arguments:
    arguments = ["miRNA","experiment","condition", "rep", "-mm_binary",
                 "-read_len", "-jobs", "-test_binary"]
    args = parse_arguments(arguments)
    print(args)
    mirna, experiment, condition, rep, mm, read_len, jobs, test = args
    # Allocate job argument:
    if not jobs:
        jobs = 20
    if not read_len:
        read_len = 120

    # Update output extensions depending on test argument:
    if rep:
        ext = ",%s" %(rep)
    else:
        ext = ""
    ext_out = ext
    ext_out = "%s_%snt" %(ext_out, read_len)
    if mm:
        ext_out = "%s_mm" %(ext_out)
    if test:
        ext_out = "%s_test" %(ext_out)
    reads_path = get_analysis_path(mirna, experiment, condition, "reads",
                                   ext=ext)

    variant_path = get_analysis_path(mirna, experiment, condition, "counts",
                                     ext=ext_out)

    print(reads_path)
    print(variant_path)
    # Get the variant dictionary:
    var_count_df = get_variant_dictionary(read_len=read_len)
    print(var_count_df[:10])

    args = [var_count_df]

    threads = multiproc_file(reads_path, int(jobs), assign_variant, test,
                                  *args)
    print("_"*20 + "Finished with multiprocessing." + "_"*20)
    var_df_all = threads[0]
    var_df_all["Counts"] = sum([thread["Counts"] for thread in threads])
    print(315)
    print(var_df_all[var_df_all["Counts"] > 0])
    
    print(variant_path)
    var_df_all.to_csv(variant_path, sep="\t", index=False)

    print_time_elapsed(time_start)
    return


################################################################################

if __name__ == "__main__":
    main()

