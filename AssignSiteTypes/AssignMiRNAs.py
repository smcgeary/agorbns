################################################################################
#AssignMiRNAs.py
################################################################################
import imp # Used to import general.py
imp.load_source("general", "general/general.py")
imp.load_source("RBNS_methods", "general/RBNS_methods.py")
from general import *
from RBNS_methods import *
# from sitetypes import get_seq_site_map
from collections import defaultdict

# marker_18nt = "AGCGUGUAGGGAUCCAAA"
# weird   =   "GGAGCGTGTAGGATCCAAATCGTATGCC"
marker_18nt_d = "AGCGTGTAGGGATCCAAA"
# marker_30nt = "GGCAUUAACGCGGCCGCUCUACAAUAGUGA"
#                # GGTCGTTTCCCGGCCCATGCACCATCGT
marker_30nt_d = "GGCATTAACGCGGCCGCTCTACAATAGTGA"

# adapter_5p = "GUUCAGAGUUCUACAGUCCGACGAUC"
adapter_5p_d = "GTTCAGAGTTCTACAGTCCGACGATC"
adapter_3p_sRS = "TCGTATGCCGTCTTCTGCTTG"
# dme_miR_14_5p = "GGGAGCGAGACGGGGACUCACU"
dme_miR_14_5p_d = "GGGAGCGAGACGGGGACTCACT"
# xtr_miR_427 = "GAAAGUGCUUUCUGUUUUGGGCG"
xtr_miR_427_d = "GAAAGTGCTTTCTGTTTTGGGCG"



# FUNCTIONS
marker_seqs = [marker_18nt_d, marker_30nt_d]
marker_names = ["18nt_marker", "30nt_marker"]

adapter_seqs = [adapter_5p_d, adapter_3p_sRS]
adapter_names = ["5p_adapter", "3p_adapter"]

spike_seqs = [dme_miR_14_5p_d, xtr_miR_427_d]
spike_names = ["dme-miR-14-5p", "xtr-miR-427"]


file_path = "AssignSiteTypes/mature.fa"
name_file_path = "AssignSiteTypes/mirna_name_conversion_table.txt"

lit_base_map = dict()


def GetThreePrimeBarcodeOverlap(read, adapt=adapter_3p_sRS):
    window = len(adapt)
    if window == 0:
        return(None)
    else:
        end = read[-window:]
        # print(end)
        # print(" "*(len(read) - window) + adapt)
        # print(" "*(len(read) - window) + end)
        position = read.find(adapt)
        # print(position)
        if adapt != end:
            return(GetThreePrimeBarcodeOverlap(read, adapt=adapt[:-1]))
        else:
            return(len(read) - window)



with open(name_file_path , "r+") as file_in:
    lines = list(file_in)
    lines = lines[1:]
    for line in lines:
        # Split the rows by the tabs
        mir_lit, mir_base, mir_ID, family, species = line.strip().split("\t")
        # Get the human species
        spec_mirbase = mir_base[:3]
        if spec_mirbase == "hsa":
            key = mir_base[4:]
            # Check if name is already in the dictionary:
            if key in lit_base_map:
                # This checks if it is only there once
                if type(lit_base_map[key]) == str:
                    lit_base_map[key] = [lit_base_map[key]] + [mir_lit[4:]]
                # This checks if it is a list:
                elif type(lit_base_map[key]) == list:
                    if sum(["miR-500" in i for i in lit_base_map[key]]) != 0:
                        print(lit_base_map[key])
                    lit_base_map[key].append(mir_lit[4:])
                else:
                    raise (ValueError("Not a string or list"))
            else:
                lit_base_map[mir_base[4:]] = mir_lit[4:]
print(lit_base_map["miR-500a-5p"])

print("here")
tick = 0
for key, value in lit_base_map.items():
    value_orig = value
    if type(value) != str:
        keeps = [i[-3:] != "-3p" and i[-3:] != "-5p" for i in value]
        value = [i[0] for i in zip(value, keeps) if i[1]]
        if len(value) == 0:
            value = value_orig[0]
        elif len(value) == 1:
            value = value[0]
        elif len(value) > 1:
            if value[0] == "miR-128":
                value = value[0]
            elif sum(["miR-500" in i for i in value]) != 0:
                value = value[0]
                print(value)
            else:
                len_values = [len(i) for i in value]
                [value] = [i for i in value if len(i) == max(len_values)]
        if type(value) != str:
            raise(ValueError("value not a string"))
        lit_base_map[key] = value
    tick += 1


seqs = []
mirnas = []
print("third loop")
with open(file_path, "r+") as file_in:
    lines = list(file_in)
    for line in zip(*[iter(lines)]*2):
        name, seq = line
        mirna = name.split(" ")[0][1:]
        seq = re.sub("U", "T", seq.strip())[:18]
        species = mirna.split("-")[0]
        if species == "hsa":
            mirna = mirna[4:]
            if mirna in lit_base_map:
                mirna = lit_base_map[mirna]
            mirnas += [mirna]
            seqs += [seq]

print("fourth loop")
mirna_seq_map = defaultdict(list)
for seq, mirna in zip(seqs, mirnas):
    mirna_seq_map[seq].append(mirna)

mirna_seq_map_new = dict()
mirna_seq_map = {key : "/".join(value) for key, value in mirna_seq_map.items()}


seq_mirna_map = {value : key for key, value in mirna_seq_map.items()}

for name, seq in zip(spike_names, spike_seqs):
    seq_mirna_map[name] = seq[:18]


def assign_mirna(read_seqs, seq_mirna_map, test):
    time_start = time.time()
    counter = [0, 0, 0]
    seq_dict = defaultdict(list)
    for i_r, r in enumerate(read_seqs):
        # Get the read number:
        seq = r.strip()
        adapter_pos = GetThreePrimeBarcodeOverlap(seq)
        barcode =  seq[:14]
        seq_segment = seq[14:adapter_pos]
        if seq in ["AGGCAGGCAGGCGGTGGAATGTAAAGAAGTATGTATTCGT",
                   "GGCTGCTGGTTGGATGGAATGTAAAGAAGTATGTATTCGT",
                   "CGACATCGATCAGTTGGAATGTAAAGAAGTATCGTATGCC"]:
            print(seq)
            print("-"*14 + "."*len(seq_segment) + "-"*(len(seq) - adapter_pos))
            print(" "*14 + "123456789!12345678")
            # print("adapter_pos:")
            # print(adapter_pos)
            # print("barcode:")
            # print(barcode)
            # print("seq_segment:")
            # print(seq_segment)
            # print(len(seq_segment) >= 18)
            # print(len(seq_dict[seq_segment[:18]]))
            sys.stdout.flush()
        if len(seq_segment) >= 18:
            seq_dict[seq_segment[:18]].append(barcode[:18])
        if seq in ["AGGCAGGCAGGCGGTGGAATGTAAAGAAGTATGTATTCGT",
                   "GGCTGCTGGTTGGATGGAATGTAAAGAAGTATGTATTCGT",
                   "CGACATCGATCAGTTGGAATGTAAAGAAGTATCGTATGCC"]:
            print(len(seq_dict[seq_segment[:18]]))
            sys.stdout.flush()
        # if adapter_5p_d[:18] == seq[14:14+18] and test:
        #     print(seq[:adapter_pos])
        #     print(" "*14 + adapter_5p_d)
        sys.stdout.flush()
    return seq_dict

def main():
    time_start = time.time()
    arguments = ["miRNA", "experiment", "condition", "-jobs", "-test_binary",
                 "-no_marker_binary", "-no_adapter_binary"]
    args = parse_arguments(arguments)
    mirna, experiment, condition, jobs, test, no_marker, no_adapter = args
    if not jobs:
        jobs = 20


    ext = ""
    if not no_adapter:
        for name, seq in zip(adapter_names, adapter_seqs):
            seq_mirna_map[name] = seq[:18]
    else:
        ext += "_noadapter"
    if not no_marker:
        for name, seq in zip(marker_names, marker_seqs):
            seq_mirna_map[name] = seq[:18]
    else:
        ext += "_nomarker"
    if test:
        ext += "_test"

    args = [seq_mirna_map, test]

    reads_path = get_analysis_path(mirna, experiment, condition, "reads")
    threads = multiproc_file(reads_path, int(jobs), assign_mirna, test, *args)

    mirna_count_dict = {name : sum([len(thread[seq_mirna_map[name]])
                                    for thread in threads])
                        for name in seq_mirna_map}
    mirna_count_unique_dict = {name : sum([len(list(set(thread[seq_mirna_map[name]])))
                                    for thread in threads])
                        for name in seq_mirna_map}

    uncounted = sum([sum([len(thread[name]) for name in thread
                          if name not in seq_mirna_map.values()])
                     for thread in threads])
    mirna_count_dict["Unmapped"] = uncounted
    mirna_count_unique_dict["Unmapped"] = uncounted

    print(mirna_count_dict["miR-1"])
    print("unmapped: %s" %(uncounted))
    uncounted_plus = uncounted
    print(no_marker)
    if not no_marker:
        print(mirna_count_dict["18nt_marker"])
        print(mirna_count_dict["30nt_marker"])
        uncounted_plus += mirna_count_dict["18nt_marker"]
        uncounted_plus += mirna_count_dict["30nt_marker"]
    if not no_adapter:
        print(mirna_count_dict["5p_adapter"])
        print(mirna_count_dict["3p_adapter"])
        uncounted_plus += mirna_count_dict["5p_adapter"]
        uncounted_plus += mirna_count_dict["3p_adapter"]
    print("unmapped plus: %s" %(uncounted_plus))

    mirna_count_df = pd.DataFrame.from_dict(mirna_count_dict, orient="index")
    mirna_count_unique_df = pd.DataFrame.from_dict(mirna_count_unique_dict, orient="index")

    table_path = get_analysis_path(mirna, experiment, condition,
                                     "AGO_pur_counts", ext=ext)
    unique_table_path = get_analysis_path(mirna, experiment, condition,
                                     "AGO_pur_counts", ext=ext + "_unique")

    print(table_path)
    print(unique_table_path)
    mirna_count_df.to_csv(table_path, sep="\t")
    mirna_count_unique_df.to_csv(unique_table_path, sep="\t")

    print_time_elapsed(time_start)
################################################################################

if __name__ == "__main__":
    main()

