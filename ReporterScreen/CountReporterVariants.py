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

def get_mirna_site_order(mirna_var):
    # Fill dataframe with miRNA, name, and sequence of each site.
    variants_path = "variants/oct12_final/final_order.fa"    
    # Pre-allocate the dictionary of miRNA sites:
    output = []
    with open(variants_path, "r") as file_in:
        header = file_in.readline().strip()
        while header:
            mirna, site, index = header.split("> ")[1].split("_")
            sequence = file_in.readline().strip()
            # If site not yet in the subdictionary, add it:
            if mirna == mirna_var and site not in output:
                output += [site]
            header = file_in.readline().strip()
    return output



def get_site_sequence_map_dictionary():
    # Fill dataframe with miRNA, name, and sequence of each site.
    variants_path = "variants/oct12_final/final_order.fa"    
    # Pre-allocate the dictionary of miRNA sites:
    mirna_site_map = {mirna : {} for mirna in MIRNAS}
    with open(variants_path, "r") as file_in:
        header = file_in.readline().strip()
        while header:
            mirna, site, index = header.split("> ")[1].split("_")
            sequence = file_in.readline().strip()
            # If site not yet in the subdictionary, add it:
            if site not in mirna_site_map[mirna]:
                mirna_site_map[mirna][site] = [sequence]
            # Otherwise, add the sequence.
            else:
                mirna_site_map[mirna][site] += [sequence]
            header = file_in.readline().strip()
    return mirna_site_map
    # Now this is a dictionary of dictionaries, where each ultimate value is
    # a list of the 184 variants ascribed to that site.
    # Make a new dictionary in which the key is the first 100 nt of the variant,
    # and the value is a tuple with the mirna, site, the variant index, and the
    # the full sequence.
def get_variant_dictionary(mirna_site_map, trim=True):
    if trim:
        start = 17
        stop = 117
    else:
        start = 14
        stop = 120
    seq100nt_variant_map = {}
    for mirna, sites in mirna_site_map.items():
        for site, sequences in sites.items():
            # Conditional making sure each mirna-and-site combination has
            # 184 different variants.
            if len(sequences) != 184:
                print("not correct number of variants!")
                return
            # This i variable is there to keep track of iterating through each
            # sequencing list.
            i = 0
            for sequence in sequences:
                if "CAGAAACTGTAGCCCTTCGAACTTGAGGCTACTTGGGGATTGGTTCGAAGAAATACATCCAAACCGAAAGTGCCATTAGCTTGTGGATCTCAGCCTTTGG" in sequence:
                    print("CAGAAACTGTAGCCCTTCGAACTTGAGGCTACTTGGGGATTGGTTCGAAGAAATACATCCAAACCGAAAGTGCCATTAGCTTGTGGATCTCAGCCTTTGG")
                    print(sequence[start:stop])
                seq100nt_variant_map[sequence[start:stop]] = (mirna, site, i,
                                                        sequence)
                i += 1
    return seq100nt_variant_map

def make_variant_count_dictionary(mirna_site_map):
    output = {mirna : {} for mirna in MIRNAS}
    for mirna in MIRNAS:
        for site in mirna_site_map[mirna]:
            output[mirna][site] = [0]*184
    return(output)

def assign_variant(read_seqs, mirna_site_map, seq100nt_variant_map):
    time_start = time.time()
    output_dic = make_variant_count_dictionary(mirna_site_map)
    print(output_dic.values[0])
    sys.stdout.flush()
    return
    unmapped = 0
    for i_r, r in enumerate(read_seqs):
        # Get the read number:
        variant = r.strip()
        if variant in seq100nt_variant_map:
            mirna, site, i, sequence = seq100nt_variant_map[variant]
            output_dic[mirna][site][int(i)] += 1
        else:
            unmapped += 1
    return([output_dic, unmapped])

def assign_variant_mm(read_seqs, mirna_site_map, seq100nt_variant_map):
    time_start = time.time()
    output_dic = make_variant_count_dictionary(mirna_site_map)
    unmapped = 0
    tot_mapped = 0
    mirnas = ["miR-1", "let-7a", "miR-155", "miR-124", "lsy-6", "miR-7-23nt"]
    sitelist = "paperfinal"
    seq25nt_left = {key[12:37] : item
                    for key, item in seq100nt_variant_map.items()}
    seq25nt_mid = {key[37:62] : item
                    for key, item in seq100nt_variant_map.items()}
    seq25nt_right = {key[62:87] : item
                    for key, item in seq100nt_variant_map.items()}

    for i_r, r in enumerate(read_seqs):
        if i_r % 100000 == 0:
            print(i_r)
        # Get the read number:
        r = r.strip()
        mapped = False
        sys.stdout.flush()
        if r in seq100nt_variant_map:
            mirna, site, i, sequence = seq100nt_variant_map[r]
            output_dic[mirna][site][int(i)] += 1
            mapped = True
            tot_mapped += 1
        else:
            sys.stdout.flush()
            for i in range(len(r) - 25 + 1):
                r_i = r[i:i + 25]
                map_ = False
                if r_i in seq25nt_left:
                    map_ = seq25nt_left[r_i]
                elif r_i in seq25nt_mid:
                    map_ = seq25nt_mid[r_i]
                elif r_i in seq25nt_right:
                    map_ = seq25nt_right[r_i]
                sys.stdout.flush()
                if map_:
                    mirna, site, i, sequence = map_
                    _mirna = Mirna(mirna)
                    _sitelist = SiteList(_mirna, "paperfinal", 100)
                    site_seq = _mirna[site]
                    _variant = ReporterVariant(r)
                    _variant.get_all_sites_in_reporter_variant(_sitelist)
                    ts = _variant.topsite
                    sys.stdout.flush()
                    site_bool = _variant.topsite
                    other_topsites = []
                    other_toppos = []
                    for mirna_i in MIRNAS:
                        if mirna_i != mirna:
                            _mirna = Mirna(mirna_i)
                            _sitelist = SiteList(_mirna, "paperfinal", 100)
                            _variant = ReporterVariant(r)
                            _variant.get_all_sites_in_reporter_variant(_sitelist)
                            ts_other = _variant.topsite.name
                            ts_pos = _variant.topsite.l
                            other_topsites.append(ts_other)
                            other_toppos.append(ts_pos)
                    # print(other_topsites)
                    # print(other_toppos)
                    unique_other_topsites = list(set(other_topsites))
                    # print(unique_other_topsites)
                    # print(ts.name)
                    # print(mirna)
                    # print(unique_other_topsites == ["None"])
                    sys.stdout.flush()
                    if ts.name == site and unique_other_topsites == ["None"]:
                        # var = sequence[17:117]
                        # mismatches = [int(i == j) for i, j in zip(list(r), list(var))]
                        output_dic[mirna][site][int(i)] += 1
                        mapped = True
                        tot_mapped += 1
                        sys.stdout.flush()
                        break
        if not mapped:
            unmapped += 1
    return([output_dic, unmapped])





def main():
    time_start = time.time()
    # Parse arguments:
    arguments = ["miRNA","experiment","condition", "-rep", "-mm_binary", "-jobs",
                 "-test_binary"]
    args = parse_arguments(arguments)
    print(args)
    mirna, experiment, condition, rep, mm, jobs, test = args
    # Allocate job argument:
    if not jobs:
        jobs = 20
    # Update output extensions depending on test argument:
    if rep:
        ext = ",%s" %(rep)
    else:
        ext = ""
    ext_out = ext
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
    mirna_site_map = get_site_sequence_map_dictionary()
    if mm:
        trim = False
    else:
        trim = True
    seq100nt_variant_map = get_variant_dictionary(mirna_site_map, trim=True)

    # seq25nt_variant_map_1 = collections.defaultdict(str)
    # seq25nt_variant_map_2 = collections.defaultdict(str)
    # seq25nt_variant_map_3 = collections.defaultdict(str)
    # seq25nt_variant_map_4 = collections.defaultdict(str)
    # seq25nt_variant_map_5 = collections.defaultdict(str)

    # for key, value in seq100nt_variant_map.items():
    #     seq25nt_variant_map_1[key[:25]] = value
    #     seq25nt_variant_map_2[key[25:50]] = value
    #     seq25nt_variant_map_3[key[50:75]] = value
    #     seq25nt_variant_map_4[key[75:]] = value
    #     seq25nt_variant_map_5[key[10:35]] = value

    # Assign the arguments for the function.
    args = [mirna_site_map, seq100nt_variant_map]
    # Use mismatch function if "mm" argument is true.
    if mm:
        assign_variant_func = assign_variant_mm
        # args = args + [seq25nt_variant_map_1, seq25nt_variant_map_2,
        #                seq25nt_variant_map_3, seq25nt_variant_map_4,
        #                seq25nt_variant_map_5]
    else:
        assign_variant_func = assign_variant
    print(assign_variant_func)
    print(list(seq100nt_variant_map.keys())[0])
    return
    # Use multiprocessing on the function assign variants:
    threads = multiproc_file(reads_path, int(jobs), assign_variant_func, test,
                                  *args)
    return
    print("_"*20 + "Finished with multiprocessing." + "_"*20)
    print("")
    output = make_variant_count_dictionary(mirna_site_map)

    # Preallocate the mapped variable:
    mapped = 0
    for mirna in MIRNAS:
        for site in output[mirna]:
            for i in range(184):
                # Sum all of the instances of one variant across all the
                # microprocessed threads:
                mapped_i = sum([thread[0][mirna][site][i]
                                for thread in threads])
                # Add this sum to the apporpirate part of the totalled
                # dictionary.
                output[mirna][site][i] = mapped_i
                # Add to the tally of mapped variants.
                mapped += mapped_i
    # Add the second multiprocesses variable together to get the unmapped
    # variable.
    print("mapped: %s" %(mapped))
    unmapped = sum([thread[1] for thread in threads])
    print("unmapped: %s" %(unmapped))
    
    print(variant_path)
    with open(variant_path, "w+") as file_out:
        for mirna in MIRNAS:
            site_order = get_mirna_site_order(mirna)
            for site in site_order:
                text_out = "\t".join(
                    [mirna, site] + [str(i) for i in output[mirna][site]]
                )
                file_out.write(text_out + "\n")




    print_time_elapsed(time_start)

    return


################################################################################

if __name__ == "__main__":
    main()

