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
np.core.arrayprint._line_width = 200

## GLOBAL CONSTANTS ############################################################

# FUNCTIONS
def get_all_8mer_mm_sites(_mirna):
  # First get the 8mer site and find the fist nucleotide.
  site_8mer_vec = list(_mirna["8mer"])[::-1]
  names = ["8mer-mm%s%s" %(j, i + 1)
           for i in range(1, 7)
           for j in ["A", "C", "G", "T"] if j != site_8mer_vec[i]]
  # Output `site_names` as a 1-D vector.
  return(names)






def assign_site(read_seqs, _mirna):
    time_start = time.time()
    # Define the names of the possible seed sites.
    names_8nt = get_all_8mer_mm_sites(_mirna)
    sites_8nt = names_8nt
    sites_7nt = ["7mer-m8", "7mer-A1"]
    sites_6nt = ["6mer-m8", "6mer", "6mer-A1"]
    mir_thrp_seq = get_rc(_mirna.seq[8:])
    counts = defaultdict(int)
    mm_counts = defaultdict(int)
    full_names_all = []
    mm_names_all = []
    # exit_flag = False
    for i_r, r in enumerate(read_seqs):
        # Hierarchically look for each of the relevant seed sites, in descending
        # order of the 8mer, an 8mer mismatch site, a 7-nt site (7mer-m8 or
        # 7mer-A1, or a 6mer, 6mer-A1 or 6mer-m8).
        site_found = [i for i in ["8mer"] if _mirna[i] in r]
        if not site_found:
            site_found = [i for i in sites_8nt if _mirna[i] in r]
        if not site_found:
            site_found = [i for i in sites_7nt if _mirna[i] in r]
        if not site_found:
            site_found = [i for i in sites_6nt if _mirna[i] in r]
        # If found at least one seed site, look for threeprime sites compatible
        # with each of the seed sites (since this could be different).
        if site_found:
            full_names = []
            mm_names = []
            # Get the left-most starting position of each seed site.
            l_seed_pos = [r.find(_mirna[i]) for i in site_found]
            for (site_found_i, l_seed_pos_i) in zip(site_found, l_seed_pos):
                # print("r: %s" %r)
                # print("r: %s" %r[:l_seed_pos_i])
                # Correct the starting distance for not beginning with seed
                # nucleotide 8.
                if site_found_i in ["7mer-A1", "6mer"]:
                    spacing_adj = 1
                elif site_found_i == "6mer-A1":
                    spacing_adj = 2
                else:
                    spacing_adj = 0
                # Get the longest substring of the miRNA 3-prime end.
                [len_thrp, l_mirAndThrPs_pre] = LCSubString(mir_thrp_seq, r[:l_seed_pos_i])
                # Grab out the starting position of the substring within the
                # miRNA 3-prime end and within the target sequence.
                l_mirAndThrPs_pre = [(i, j) for (i, j) in zip(l_mirAndThrPs_pre[0], l_mirAndThrPs_pre[1])]
                # This just serves as a check, so that len_thrp_alt == len_thrp.
                [len_thrp_alt, key_l, str_l] = longest_match(mir_thrp_seq, r[:l_seed_pos_i], 4, mis=0)
                # Also get the longest sequences allowing for single and double
                # mismatches (not used in the paper).
                # [len_thrp_mm, key_l_mm, str_l_mm] = longest_match(mir_thrp_seq, r[:l_seed_pos_i], 4, mis=1)
                # [len_thrp_mm2, key_l_mm2, str_l_mm2] = longest_match(mir_thrp_seq, r[:l_seed_pos_i], 4, mis=2)
                # If a contiguous site is greater than four nucleotides in
                # length, proceed to define the 3-prime pairing.
                if len_thrp >= 4:
                    # Iterate over the starting positions (paired tuple giving
                    # the miRNA and target starting position).
                    for (l_mir_pos, l_tar_pos) in l_mirAndThrPs_pre:
                        # Get ends.
                        #  1@987654321!987654321
                        #  NNNNNNNNNNNNNSSSSSSSS
                        #     xxxxxxx
                        #  |                       21
                        #  012|                        -  3
                        #     123456|                         - 7 + 1 = 12
                        #  |                       21
                        #  012|                        -  3           = 19
                        l_mirThrp_pos = len(_mirna) - l_mir_pos - len_thrp + 1
                        r_mirThrp_pos = len(_mirna) - l_mir_pos
                        # Construct name.
                        thrp_name = "%smer-m%s.%s" %(len_thrp, l_mirThrp_pos, r_mirThrp_pos)
                        # Get target position.
                        #  0123456789!123456789@1234
                        #  NNNNNNNNNNNNNNNNNSSSSSSSS
                        #       xxxxxxx            .
                        #  01234|                  .  5
                        #       123456|            .    + 7 = 12
                        #  012345678!1|            .          12 
                        #  0123456789!123456|      .  17 - 12 + 9 - spacing _adj
                        #             |321!987654321          = 14
                        r_tar_pos = l_tar_pos + len_thrp
                        thrp_mir_pos = l_seed_pos_i - r_tar_pos + 9 - spacing_adj
                        full_name = "%s|%s|%s" %(thrp_name, thrp_mir_pos, site_found_i)
                        full_names.append(full_name)
                        if len_thrp != len_thrp_alt:
                            print(len_thrp)
                            print(len_thrp_alt)
                            print("discrepancy between methods")
                            return()
                        # If double mismatch greater than 11, use it.
                        # if len_thrp_mm2 > 11:
                        #     print("this happens")
                        #     print(full_names)
                        #     print(i_r)
                        #     print(r)
                        #     return()
                        #     len_thrp_mm = len_thrp_mm2
                        #     key_l_mm = key_l_mm2
                        #     str_l_mm = str_l_mm2
                        # if len_thrp_mm > len_thrp and (len_thrp_mm2 > 11 or len_thrp > 5):
                        #     for i in range(len(str_l_mm)):
                        #         mir_l_mm = len(mir_thrp_seq) + 8 - key_l_mm[i] - len_thrp_mm + 1
                        #         mir_r_mm = len(mir_thrp_seq) + 8 - key_l_mm[i]
                        #         thrp_mm_name =  "%smer-m%s.%s" %(
                        #             len_thrp_mm, mir_l_mm, mir_r_mm
                        #         )
                        #         spacing = l_seed_pos_i - (str_l_mm[i] + len_thrp_mm) + 9 - spacing_adj
                        #         mir_seq = mir_thrp_seq[key_l_mm[i]:key_l_mm[i]+ len_thrp_mm][::-1]
                        #         tar_seq = r[str_l_mm[i]:str_l_mm[i] + len_thrp_mm][::-1]
                        #         suffix = ""
                        #         for i in range(len(mir_seq)):
                        #             mir_nuc = mir_seq[i]
                        #             tar_nuc = tar_seq[i]
                        #             if mir_nuc != tar_nuc:
                        #                 pos = i + mir_l_mm
                        #                 suffix += "mm%s%s" %(tar_nuc, pos)
                        #         thrp_mm_name += "%s|%s|%s" %(suffix, spacing, site_found_i)
                        #         mm_names.append(thrp_mm_name)
                        # else:
                        mm_names.append(full_name)
                else:
                    full_names.append(site_found_i)
                    mm_names.append(site_found_i)
                # print("XXXXXXX")
                # print(full_names)
                # print(mm_names)
                # print("YYYYYYYY")
                if len_thrp >= 12:
                    print("12 or more")
                    print(mm_names)
                    print(full_names)
                if len_thrp < 4:
                    print("less than 4")
                    print(mm_names)
                    print(full_names)
            # # This part makes no sense.
            # if len(full_names) > 1:
            #     print("more than 1")
            #     print(i_r)
            #     print(r)
            #     exit_flag = True
            #     for (l_seed_pos_i, site_name_i) in zip(l_seed_pos, site_found):
            #         seed_len = int(site_name_i.split("mer")[0])
            #         print(" "*l_seed_pos_i + "_"*seed_len + " "*(50 - l_seed_pos_i - seed_len) + site_name_i)

            #     print(" "*l_tar_pos + "_"*len_thrp + " "*(50 - l_tar_pos - len_thrp) + thrp_name)
            #     print(site_found)
            #     print(full_names)
            # if len(mm_names) == 4 or len(full_names) == 4:
            #     print(full_names)
            #     print(mm_names)
            #     return()
        else:
            full_names = []
            mm_names = []
        full_names_set = set(full_names)
        mm_names_set = set(mm_names)
        if len(full_names_set) != len(full_names):
            full_names = list(full_names_set)
        if len(mm_names_set) != len(mm_names):
            mm_names = list(mm_names_set)
        # if exit_flag:
        #     print(site_found)
        #     print(full_names)
        #     print(mm_names)
        #     print(i_r)
        #     print(r)
        #     return()
        # if site_found and len(site_found) > 1:
        #     print(site_found)
        #     print(l_seed_pos)
        #     print(full_names)
        #     print(mm_names)
        #     print(i_r)
        #     print(r)
        #     return()

        counts[str(len(full_names))] += 1
        mm_counts[str(len(mm_names))] += 1
        if full_names == []:
            full_names = "No seed"
            mm_names = "No seed"
        else:
            full_names = ";".join(full_names)
            mm_names = ";".join(mm_names)
        full_names_all.append(full_names)
        mm_names_all.append(mm_names)
    # print(counts)
    # print(mm_counts)
    # print(i_r)
    # print(sum(mm_counts.values()))
    # print(sum(mm_counts.values()))
    return([full_names_all, mm_names_all])

def main():
    time_start = time.time()
    arguments = [
        "miRNA", "-test_binary"
    ]
    args = parse_arguments(arguments)
    mirna, test = args
    if mirna == "let-7a": 
        _mirna = Mirna("let-7a-21nt")
    else:
        _mirna = Mirna(mirna)
    # FIX THE FACT THAT THE BECKER SITES GET SKIPPED OVER IF THEY ARE AT
    # POSITION [0] OF THE READ. This seems to be a hack that I did to deal with
    # the error that occurs when I give the LCSubString function a string of
    # length 0, which necessarily must happen if the seed site is at position 0.
    # The correct way to deal with this would be to edit the LCSubString
    # function to automatically return two empty lists as output alongside the
    # integer 0 if either string is of length 0. (i.e., it would return
    # [0, [ [], [] ] ]).

    # This may change the results of the Becker re-analysis but the likelihood
    # that it gives a different conclusion seems essentially 0. All it means is
    # that the Kd fold changes might change due to not being influenced by sites
    # that occur at the first position of the library.
    # Assign three different extensions for the three different count files
    # being made.
    # THIS HAS BEEN FIXED, made two "no seeds" become 8mers, added a 6mer-m8
    # site to a multi site, and made four "no seed"'s become 6mer-A1's.
    path = "Becker,Ober-Reynolds_et_al_data/%s_v2.txt" %mirna
    output_path = "Becker,Ober-Reynolds_et_al_data/%s_v2_with_names_no_mm.txt" %mirna

    file = open(path, "r")

    header = file.readline()
    sequences_out = []
    line = file.readline().strip()
    tick = 0
    kd_med_all = []
    kd_lCI_all = []
    kd_uCI_all = []
    kon_all = []
    R2_kon_all = []
    num_points_kon_all = []
    kcleave_all = []
    R2_kcleave_all = []
    dG_all = []
    while line:
        line_strip = line.split("\t")
        (sequence, kd_med, kd_lCI, kd_uCI, kon,
         R2_kon, num_points_kon, kcleave, R2_kcleave, dG) = line_strip
        kd_med_all.append(kd_med)
        kd_lCI_all.append(kd_lCI)
        kd_uCI_all.append(kd_uCI)
        kon_all.append(kon)
        R2_kon_all.append(R2_kon)
        num_points_kon_all.append(num_points_kon)
        kcleave_all.append(kcleave)
        R2_kcleave_all.append(R2_kcleave)
        dG_all.append(dG)

        sequences_out += [sequence]
        line = file.readline().strip()
        tick += 1
    thread = assign_site(sequences_out, _mirna)
    print(thread[0][0])
    output = list(zip(sequences_out, thread[0], thread[1], kd_med_all,
                      kd_lCI_all, kd_uCI_all, kon_all, R2_kon_all,
                      num_points_kon_all, kcleave_all, R2_kcleave_all, dG_all))
    output_df = pd.DataFrame(
        output,
        columns=["Sequence", "Site name", "Mismatch site name", "Kd_median",
                 "Kd_lCI", "Kd_uCI", "kon", "R2_kon", "N_kon", "kcleave",
                 "R2_kcleave", "dG"]
    )
    print(output_path)
    # Output the dataframes to the three paths.
    output_df.to_csv(output_path, sep="\t", index=False)
    # counts_collapsed.to_csv(site_counts_collapsed_path, sep="\t", header=False)
    print_time_elapsed(time_start)
################################################################################

if __name__ == "__main__":
    main()

