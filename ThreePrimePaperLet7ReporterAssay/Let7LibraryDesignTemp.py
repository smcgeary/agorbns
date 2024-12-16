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



#FINAL DESIGN CONSTRAINTS
kLEN = 180

def count_utr_sequence(utrs, split=False):
    """Calculates the dinucleotides frequencies of a DNA sequence.

    Args:
        utrs (list): A list of strings, containing only the characters
            A, C, G, and U.
    Returns:
        dict: A dictionary with keys that are the 16 DNA dinucleotides and the 
            keys are the corresponding frequencies within the entire list of
            sequences. The dictionary sums to 1.0.
    """
    dinucs = get_kmer_list(2)
    counts_dinuc_map = {dinuc : 0 for dinuc in dinucs}
    for utr in utrs:
        for i in range(len(utr) - 1):
            counts_dinuc_map[utr[i:i + 2]] += 1
    totals = sum(counts_dinuc_map.values())
    freq_dinuc_map = {dinuc : float(counts_dinuc_map[dinuc])/totals
                       for dinuc in dinucs}
    return freq_dinuc_map


def generate_utr_sequences(num_seq, len_seq, frac_dinuc_dict):
    """Generates a DNA sequences using the provided dinucleotide frequencies.

    Args:
        utrs (list): A list of strings, containing only the characters
            A, C, G, and U.
    Returns:
        dict: A dictionary with keys that are the 16 DNA dinucleotides and the 
            keys are the corresponding frequencies within the entire list of
            sequences. The dictionary sums to 1.0.
    """

    out = []
    a = list(frac_dinuc_dict.keys())
    p = list(frac_dinuc_dict.values())
    print(a)
    print(p)
    print(frac_dinuc_dict["GG"])
    len_seq_start = len_seq
    for i in range(num_seq):
        # print(i)
        start = True
        # print(i)
        len_seq = len_seq_start
        utr = ""
        # print("out")
        # print(out)
        while len_seq != 1:
            # print("utr")
            # print(utr)
            if start:
                dinuc = np.random.choice(a, p=p)
                start = False
                utr = dinuc
            else: 
                all_use = [i for i in zip(a, p) if i[0][0] == utr[-1]]
                a_use = [i[0] for i in all_use]
                p_use = [i[1] for i in all_use]
                p_sum = np.sum(p_use)
                p_use = [i/p_sum for i in p_use]
                dinuc = np.random.choice(a_use, p=p_use)
                utr += dinuc[1]
            len_seq -= 1
        out += [utr]
    if num_seq == 1:
        return out[0]
    else:
        return out


def add_sequence(seq, len_seq, frac_dinuc_dict):
    """Generates a DNA sequences using the provided dinucleotide frequencies.

    Args:
        seq (string): A DNA sequence containing only the characters
            A, C, G, and T.
        len_seq (int): The number of additional nucleotides to add to `seq`.
        frac_dinuc_dict (dict): A dictionary containing each of the dinucleotide
            frequencies from which the nucleotide sequence is to be sampled.
    Returns:
        seq (string): A DNA sequence with `len_seq` nucleotides added to its
            3-prime end. 
    """

    a = list(frac_dinuc_dict.keys())
    p = list(frac_dinuc_dict.values())
    while len_seq != 0:
        all_use = [i for i in zip(a, p) if i[0][0] == seq[-1]]
        a_use = [i[0] for i in all_use]
        p_use = [i[1] for i in all_use]
        dinuc = np.random.choice(a_use, p=p_use/np.sum(p_use))
        seq += dinuc[1]
        len_seq -= 1
    return seq

def check_sequence_for_sites(seq, _mirna, offset, position, seed_pos):
    """Checks a reporter sequence for the extent of its contiguous
        complementarity to either the seed or 3-prime end of a miRNA sequence.

    Args:
        seq (string): A DNA sequence containing only the characters A, C, G, and
            T.
        _mirna (Mirna): An object representing a miRNA.
        thrp (bool): A True/False boolean determining whether the 3-prime end of
            of the miRNA is checked, or the seed. The default is `False`, such
            that it looks for seed complementarity.
    Returns:
        None. The function is intended to print the greatest complementarity for
        visual assessment of the reporter sequence.
    """
    # 1. Make sure that there is no longer complementarity than 6-nt to let-7a,
    # either all together, or just the seed, or the three prime end.
    mir_seq = get_rc(_mirna.seq[1:]) + "A"
    mir_thrp_seq = get_rc(_mirna.seq[8:])
    mir_seed_seq = get_rc(_mirna.seq[1:8]) + "A"
    # Get the miR-1 8mer sequence.
    mir1_seed_seq = get_rc(Mirna("miR-1").seq[1:8]) + "A"
    mir_seqs = [mir_seq, mir_thrp_seq, mir_seed_seq, mir1_seed_seq]
    mir_lims = [6, 4, 6, 6]
    pos_pairing = []
    full_pairing = [" "]*len(seq)
    char = ["f", "t", "s", "1"]
    for i in range(4):
        # First define the maximum overlap of the sequence with the substring.
        [len_comp, l_mirAndThrPs] = LCSubString(mir_seqs[i], seq)
        if i == 3:
            print(len_comp)
        # If the length is larger than the maximum allowed limit for that type
        # of sequence, iterate from this length to the minimum allowed limit,
        # marking each of the sequences with pairing.
        if len_comp >= mir_lims[i]:
            for len_i in range(len_comp, mir_lims[i] - 1, -1):
                [len_comp, l_mirAndThrPs] = LCSubString(mir_seqs[i], seq, alt_len=len_i)
                for j in range(len(l_mirAndThrPs[0])):
                    mir_pos = l_mirAndThrPs[0][j]
                    seq_pos = l_mirAndThrPs[1][j]
                    pairing_range = [k for k in range(seq_pos, seq_pos + len_comp)]
                    for k in pairing_range:
                        if char[i] == "f":
                            full_pairing[k] = char[i].upper()
                        else:
                            full_pairing[k] = char[i]
                    pos_pairing += pairing_range


    pos_pairing = list(set(pos_pairing))

    full_pairing = "".join(full_pairing)
    # if "1" in full_pairing:
    print(seq)
    print(full_pairing)
    pos_pairing.sort()
    return(pos_pairing)


def mutate_seq(seq, allowed_pos, paired_pos, offset, position, seed_pos):
    """Mutates a DNA sequence as a function of the list of all of the paired
        positions found by the `check_sequence_for_sites` function, the list of
        allowed positions outputed by the 'add_sites_to_seq' function, as well
        as the 'offset' 'position' and 'seed_pos' arguments.

    Args:
        seq (string): A DNA sequence containing only the characters A, C, G, and
            T.
        allowed_pos (list): Lists all of the positions that pairing is allowed,
            as put into the sequence using the `add_sites_to_seq` function.
        paired_pos (list): Lists all of the positions that were determined to
            contribute in some way to pairing in the `check_sequence_for_sites`
            function.
        offset (str): A str of the offset sequence that was added by the
            `add_sites_to_seq` function. This is used to convert the indeces
            from their positions in the context of site-added sequence to that
            of the starting contextual sequence.
        position (str): Also used for the conversion of the indeces from the
            site-added sequence to that of the starting contextual sequence,
            since the left, right, or both sites might be shifted by the offset
            length.
        seed_pos (tup): A tuple or list giving the right-hand definition of each
            of the two site positions, used to consistently position the seed and
            3-prime sites within the reporter construct.
    Returns:
        seq_new (str): A new DNA sequence that has been mutated in the center of
        each of the contiguous regions of complementarity as determined by the
        `check_sequences_for_sites` function.
    """
    dif_seq = [i for i in paired_pos if i not in allowed_pos]
    if position in ["left", "dual"]:
        for i in range(len(dif_seq)):
            pos = dif_seq[i]
            if pos > seed_pos[0] - 9:
                dif_seq[i] = pos - len(offset)
    if position in ["right", "dual"]:
        for i in range(len(dif_seq)):
            pos = dif_seq[i]
            if pos > seed_pos[1] + len(offset)*(position == "dual") - 9:
                dif_seq[i] = pos - len(offset)
    seq_list = list(seq)
    contig_nucs = []
    for i in dif_seq:
        if len(contig_nucs) == 0:
            contig_nucs.append(i)
        elif i == contig_nucs[-1] + 1:
            contig_nucs.append(i)
        else:
            # Mutate the central nucleotide within this contiguous segment.
            pos_mut = math.ceil(sum(contig_nucs)/len(contig_nucs))
            nuc = seq_list[pos_mut]
            mut_nuc = np.random.choice([i for i in ["A", "C", "G", "T"] if i != nuc])
            seq_list[pos_mut] = mut_nuc
            # Re-assign the contiguous set of nucleotides.
            contig_nucs = [i]
    # Mutate the last segment.
    pos_mut = math.ceil(sum(contig_nucs)/len(contig_nucs))
    nuc = seq_list[pos_mut]
    mut_nuc = np.random.choice([i for i in ["A", "C", "G", "T"] if i != nuc])
    seq_list[pos_mut] = mut_nuc
    # Form the new sequence from the list.
    seq_new = "".join(seq_list)

    return(seq_new)

def add_sites_to_seq(seq, _mirna, seed, thrp, offset, position, seed_pos):
    """ Adds a potentially seed-and-3-prime site to a library DNA molecule, at
        either of the two positions defined by the lin-41 3 3-prime UTR when
        adusting for the two +1 offsets and the bulged nucleotide in the 5-prime
        seed site.

    Args:
        seq (str): A DNA sequence containing only the characters A, C, G, and
            T.
        _mirna (Mirna): An object representing a miRNA.
        seed (str): A name of a seed site type. Must be one of "8mer",
            "7mer-m8", "7mer-A1", or "6mer".
        thrp (str): A name of a thrp site, of the form "9mer-m11.19". Can
            also be "", which corresponds to no 3-prime site being added.
        offset (str): The sequence of offset nucleotides to be inserted between
            the two positions across from miRNA positions 9 and 10.
        position (str): One of "left", "right", or "dual", referring to which of
            the two site positions the site is put into.
        seed_pos (tup): A tuple or list giving the right-hand definition of each
            of the two site positions, used to consistently position the seed and
            3-prime sites within the reporter construct.
    TODO: add functionality for 
    """
    # 1a. Get the seed site from the miRNA object.
    seed_site = _mirna[seed]
    # 1b. Get the left- and right-hand positions of where to replace the seed
    #     site relative to the position of the A1.
    if seed in ["8mer", "7mer-m8", "8mer-w6"]:
        seed_l = -8
    elif seed in ["7mer-A1", "6mer"]:
        seed_l = -7
    if seed in ["8mer", "7mer-A1", "8mer-w6"]:
        seed_r = 0
    elif seed in ["7mer-m8", "6mer"]:
        seed_r = -1
    # 2a. Get the thrp site and the left- and right-hand positions of where to
    #     put the site relative to the position of the A1.
    if thrp != "":
        thrp_site = _mirna[thrp]
        thrp_l = -1*int(thrp.split("mer-m")[1].split(".")[1])
        thrp_r = -1*int(thrp.split("mer-m")[1].split(".")[0]) + 1
    else:
        thrp_site = ""
        thrp_l = 0
        thrp_r = 0
    paired_nucs = []
    # Add the seed, threep, and offset sequences into the reporter.
    if position in ["left", "dual"]:
        seed_0 = seed_pos[0]
        seq = seq[:seed_0 + seed_l] + seed_site + seq[seed_0 + seed_r:]
        paired_nucs += [i for i in range(seed_0 + seed_l + len(offset), seed_0 + seed_r + len(offset))]
        seq = seq[:seed_0 + thrp_l] + thrp_site + seq[seed_0 + thrp_r:]
        paired_nucs += [i for i in range(seed_0 + thrp_l, seed_0 + thrp_r)]
        seq = seq[:seed_0 - 9] + offset + seq[seed_0 - 9:]
    if position in ["right", "dual"]:
        seed_0 = seed_pos[1] + len(offset)
        seq = seq[:seed_0 + seed_l] + seed_site + seq[seed_0 + seed_r:]
        paired_nucs += [i for i in range(seed_0 + seed_l + len(offset), seed_0 + seed_r + len(offset))]
        seq = seq[:seed_0 + thrp_l] + thrp_site + seq[seed_0 + thrp_r:]
        paired_nucs += [i for i in range(seed_0 + thrp_l, seed_0 + thrp_r)]
        seq = seq[:seed_0 - 9] + offset + seq[seed_0 - 9:]
    paired_nucs.sort()
    buffer_seq = [" "]*len(seq)
    for i in paired_nucs:
        buffer_seq[i] = "x"
    buffer_seq = "".join(buffer_seq)
    print(seq)
    print(buffer_seq)

    # Return the reporter sequence.
    return(seq, paired_nucs)




def main():
    time_start = time.time()
    # arguments = ["len_vars", "-single_site"]
    # args = parse_arguments(arguments)
    # len_vars, single_site = args
    mirna = "let-7a"
    _mirna = Mirna(mirna)

    ## 1: GENERATE UTR FREQUENCY COUNTS ########################################
    # Read in the UTRs from Kathy's table, tally up all the dinucleotide
    # frequencies, to generate a probability table to generate random 3'UTR
    # sequence. ################################################################

    # Load UTR dataframe:
    utr_path = ("/lab/bartel4_ata/kathyl/RNA_Seq/analysis/" + 
                "data/no_baseline_analysis/merged.txt")
    utr_df = pd.read_csv(utr_path, sep="\t", header=0, index_col=0)
    # 2: Test that the two functions, count_utr_sequences(), and
    # generate_utr_sequences(), work, by using generate_utr_sequences() to
    # simulate a large number of utr sequences, then calculate the dinucleotide
    # freuqncies within the simulated utrs, and compare them to the frequencies
    # originally calculated for the mRNA utr sequences.
    freq_MrnaDinuc_map = count_utr_sequence(utr_df["sequence"])

    freq_MrnaDinuc_df = pd.DataFrame.from_dict(freq_MrnaDinuc_map,
                                               orient="index",
                                               columns=["mRNA freq"])

    # 3: Load the lin-41 3'-UTR sequence, to use for the site counting.
    with open("ThreePrimePaperLet7ReporterAssay/lin41_3pUTR.txt") as file:
        lin41_utr = file.read().strip().upper()

    # Using meta-data, get the sequences of the two seed sites and the shared
    # 3p site that is found in the lin41 3 prime utr.
    lin41_seed_sites = [_mirna[i] for i in ["8mer-bA5", "7mer-m8mmG6"]]
    lin41_3p_site = _mirna["9mer-m11.19"]

    # Define the left-hand and right-hand limits of the lin41 site, so as to
    # calculate the middle and thus define the limits that should be put into
    # the 146-nt length region of the reporter library.
    # The left-hand site is given by the first exaple of the 3-prime site.
    lin41_site_l = find_all(lin41_utr, lin41_3p_site)[0]
    # The right-hand site is given by the last position of the right-hand
    # seed site, itself given by its left-hand occurence plus its length.
    lin41_site_r = max([lin41_utr.find(i) + len(i) for i in lin41_seed_sites])

    print(lin41_utr[lin41_site_l:lin41_site_r])
    # Get the center and length of the lin41 site.
    lin41_site_center = (lin41_site_l + lin41_site_r)/2
    lin41_site_len = lin41_site_r - lin41_site_l
    # The added nucleotides are the length of the plasmid fragment minus the
    # the full length of the lin41 dual site.
    len_added_nucs = 149 - lin41_site_len
    l_flank = math.floor(len_added_nucs/2)
    r_flank = math.ceil(len_added_nucs/2)
    # Define the lin41 reporter sequence to be used.
    lin41_rep_seq = lin41_utr[lin41_site_l - l_flank:lin41_site_r + r_flank]
    # Calculate the length of the bridge between the two dual sites.
    seed_site_r = [lin41_rep_seq.find(i) + len(i) for i in lin41_seed_sites]
    thrp_site_l = find_all(lin41_rep_seq, lin41_3p_site)
    n_bridge = thrp_site_l[1] - seed_site_r[0]
    print(n_bridge) 
    # Define the limits of the let7a site type:
    # DIST          l_flank       11          3        9
    # FEATURE     left flank | 3p site 1 | loop 1 | seed 1 | 
    # DIST            28          11          3       8*        r_flank - 1*
    # FEATURE     intersite  | 3p site 2 | loop 2 | seed 2 |  right flank
    # * The reason the r_flank is one nucleotide shorter is that the seed2 site
    # is actually a 7mer-m8, not an 8mer, but I am including the non-A
    # nucleotide from the intervening sequence as part of the site.
    # Define the limits of the lin41 sequence:
    l41_l_thrp1 = l_flank
    l41_r_thrp1 = l41_l_thrp1 + 9
    l41_l_seed1 = l41_r_thrp1 + 3               # This has an offset of 1
    l41_r_seed1 = l41_l_seed1 + 9               # This has 1 bulged nucleotide
    l41_l_thrp2 = l41_r_seed1 + n_bridge
    l41_r_thrp2 = l41_l_thrp2 + 9
    l41_l_seed2 = l41_r_thrp2 + 3               # This has an offset of 1
    l41_r_seed2 = l41_l_seed2 + 8

    # Get the starting positions for the seed A1 target nucleotides in a
    # non lin-41 site (i.e., one with no offsets and no bulged target site for
    # its first site).
    rep_r_seed1 = l41_r_seed1 - 2
    rep_r_seed2 = l41_r_seed2 - 3

    l_seq     = lin41_rep_seq[            : l41_l_thrp1]
    thrp1_seq = lin41_rep_seq[l41_l_thrp1 : l41_r_thrp1]
    loop1_seq = lin41_rep_seq[l41_r_thrp1 : l41_l_seed1]
    seed1_seq = lin41_rep_seq[l41_l_seed1 : l41_r_seed1]
    inter_seq = lin41_rep_seq[l41_r_seed1 : l41_l_thrp2]
    thrp2_seq = lin41_rep_seq[l41_l_thrp2 : l41_r_thrp2]
    loop2_seq = lin41_rep_seq[l41_r_thrp2 : l41_l_seed2]
    seed2_seq = lin41_rep_seq[l41_l_seed2 : l41_r_seed2]
    r_seq     = lin41_rep_seq[l41_r_seed2 :        ]

    # Check that the addition of the parts recapitulates the original site.
    bi1_seq = thrp1_seq + loop1_seq + seed1_seq
    bi2_seq = thrp2_seq + loop2_seq + seed2_seq
    full_seq_check = l_seq + bi1_seq + inter_seq + bi2_seq + r_seq
    print(lin41_rep_seq == full_seq_check)

    # Make the lin-41 sequence context, by replacing the entire sequence of the
    # two dual sites with sequence sampled from human 3-prime UTR dinucleotide
    # frequencies.
    l_seq_context = add_sequence(l_seq, len(bi1_seq) - 2, freq_MrnaDinuc_map)

    inter_context = add_sequence(inter_seq, len(bi2_seq) - 1, freq_MrnaDinuc_map)

    lin41_bg_context = (l_seq_context + inter_context + r_seq)
    check_sequence_for_sites(lin41_bg_context, _mirna, "", "", [])

    check_sequence_for_sites(lin41_rep_seq, _mirna, "", "", [])





    seqs = [generate_utr_sequences(1, 146, freq_MrnaDinuc_map),
            lin41_bg_context]

    site_list = [[       "",            "",     "",  "left"],
                 [   "8mer",            "",     "",  "left"],
                 [   "8mer",            "",     "", "right"],
                 [   "8mer",            "",     "",  "dual"],
                 ["7mer-m8",            "",     "",  "left"],
                 ["7mer-m8",            "",     "", "right"],
                 ["7mer-m8",            "",     "",  "dual"],
                 ["7mer-A1",            "",     "",  "left"],
                 ["7mer-A1",            "",     "", "right"],
                 ["7mer-A1",            "",     "",  "dual"],
                 [   "6mer",            "",     "",  "left"],
                 [   "6mer",            "",     "", "right"],
                 [   "6mer",            "",     "",  "dual"],
                 ["8mer-w6",            "",     "",  "left"],
                 ["8mer-w6",            "",     "", "right"],
                 ["8mer-w6",            "",     "",  "dual"],
                 ["8mer-w6", "4mer-m13.16",     "",  "left"],
                 ["8mer-w6", "4mer-m13.16",     "", "right"],
                 ["8mer-w6", "4mer-m13.16",     "",  "dual"],
                 ["8mer-w6", "4mer-m13.16",    "T",  "left"],
                 ["8mer-w6", "4mer-m13.16",    "T", "right"],
                 ["8mer-w6", "4mer-m13.16",    "T",  "dual"],
                 ["8mer-w6", "4mer-m13.16", "TTTT",  "left"],
                 ["8mer-w6", "4mer-m13.16", "TTTT", "right"],
                 ["8mer-w6", "4mer-m13.16", "TTTT",  "dual"],
                 ["8mer-w6", "9mer-m11.19",     "",  "left"],
                 ["8mer-w6", "9mer-m11.19",     "", "right"],
                 ["8mer-w6", "9mer-m11.19",     "",  "dual"],
                 ["8mer-w6", "9mer-m11.19",    "T",  "left"],
                 ["8mer-w6", "9mer-m11.19",    "T", "right"],
                 ["8mer-w6", "9mer-m11.19",    "T",  "dual"],
                 ["8mer-w6", "9mer-m11.19", "TTTT",  "left"],
                 ["8mer-w6", "9mer-m11.19", "TTTT", "right"],
                 ["8mer-w6", "9mer-m11.19", "TTTT",  "dual"],
                 ["8mer-w6", "9mer-m13.21",     "",  "left"],
                 ["8mer-w6", "9mer-m13.21",     "", "right"],
                 ["8mer-w6", "9mer-m13.21",     "",  "dual"],
                 ["8mer-w6", "9mer-m13.21",    "T",  "left"],
                 ["8mer-w6", "9mer-m13.21",    "T", "right"],
                 ["8mer-w6", "9mer-m13.21",    "T",  "dual"],
                 ["8mer-w6", "9mer-m13.21", "TTTT",  "left"],
                 ["8mer-w6", "9mer-m13.21", "TTTT", "right"],
                 ["8mer-w6", "9mer-m13.21", "TTTT",  "dual"]]


    # site_list = [[   "8mer", "", "",  "left"],
    #              ["8mer-w6", "", "",  "left"]]



    for i in range(2):
        seq = generate_utr_sequences(1, 146, freq_MrnaDinuc_map)
        complete = False
        attempts = 0
        while complete == False:
            variants = []
            print('start of loop')
            complete = True
            for site_type_i in site_list:
                print(site_type_i)
                site, thrp, offset, position = site_type_i
                variant, allowed_pos = add_sites_to_seq(seq, _mirna, site, thrp, offset, position, [rep_r_seed1, rep_r_seed2])
                paired_pos = check_sequence_for_sites(variant, _mirna, offset, position, [rep_r_seed1, rep_r_seed2])
                union_list = list(set(paired_pos) | set(allowed_pos))
                union_list.sort()
                if union_list != allowed_pos:
                    complete = False
                    attempts +=1
                    if attempts == 10:
                        print("Try a new sequence.")
                        seq = generate_utr_sequences(1, 146, freq_MrnaDinuc_map)
                    else:
                        print("Mutate sequence attempt %s." %attempts)
                        seq = mutate_seq(seq, allowed_pos, paired_pos, offset, position, [rep_r_seed1, rep_r_seed2])
                    break
                variants.append(variant)
        print("made it to the end")

        print("\n".join(variants))
        print(len(variants))
        return()

        variant = add_sites_to_seq(seq, _mirna, "8mer", "", "", "right", [rep_r_seed1, rep_r_seed2])
        check_sequence_for_sites(variant, _mirna, "", "", "right", [rep_r_seed1, rep_r_seed2])

        variant = add_sites_to_seq(seq, _mirna, "8mer", "", "", "dual", [rep_r_seed1, rep_r_seed2])
        check_sequence_for_sites(variant, _mirna, "", "", "dual", [rep_r_seed1, rep_r_seed2])


        variant = add_sites_to_seq(seq, _mirna, "8mer", "4mer-m13.16", "", "left", [rep_r_seed1, rep_r_seed2])
        check_sequence_for_sites(variant, _mirna, "4mer-m13.16", "", "left", [rep_r_seed1, rep_r_seed2])

        variant = add_sites_to_seq(seq, _mirna, "8mer", "4mer-m13.16", "", "right", [rep_r_seed1, rep_r_seed2])
        check_sequence_for_sites(variant, _mirna, "4mer-m13.16", "", "left", [rep_r_seed1, rep_r_seed2])

        variant = add_sites_to_seq(seq, _mirna, "8mer", "4mer-m13.16", "", "dual", [rep_r_seed1, rep_r_seed2])
        check_sequence_for_sites(variant, _mirna, "4mer-m13.16", "", "left", [rep_r_seed1, rep_r_seed2])


        # add_sites_to_seq(seq, _mirna, "8mer", "4mer-m13.16", "", "left", [rep_r_seed1, rep_r_seed2])
        # add_sites_to_seq(seq, _mirna, "8mer", "4mer-m13.16", "CCCC", "left", [rep_r_seed1, rep_r_seed2])
    return()

    

    total_sites = 0

    x_nuc_map = {"A" : "b",
                 "C" : "d",
                 "G" : "h",
                 "T" : "v"}

    _mirna_i = Mirna(mirna_i)
    _sitelist = SiteList(_mirna_i, "paperfinal", 37)
    site_seqs = [(_mirna_i.name, i.name, i.seq) for i in _sitelist.sites]
    seq_site_map = {i.name : i.seq for i in _sitelist.sites}
    unknownsites = [i.name for i in _sitelist.sites if "mer" not in i.name]
    unknownsites_dict = {i : i for i in unknownsites}
    for i in unknownsites:
        for j in unknownsites:
            if i in j and len(i) != len(j):
                outside = j.split(i)
    unknownsites_mirna_map[mirna_i] = unknownsites
    seq_site_mirna_map[mirna_i] = seq_site_map
    site_seqs_all += site_seqs
    for _site in _sitelist.sites:
        total_sites += 1
    sites_df = pd.DataFrame(site_seqs_all, columns=["miRNA", "site", "seq"])
    for i_row in range(sites_df.shape[0]):
        mirna_i, site, core = list(sites_df.iloc[i_row, ])
        s8mer = seq_site_mirna_map[mirna_i]["8mer"]
        # Part of script where the non-pairing characters are added to the 
        # sequence. Necessary for 1.) Seed sites and 2.) 3' sites, and 3.) the
        # AA-[site] sites of miR-124.
        # 1.) Seed sites:
        if "7mer-m8" in site:
            xA1 = x_nuc_map[s8mer[-1]]
            sites_df.loc[i_row, "seq"] =       core + xA1
        if "7mer-A1" in site:
            xm8 = x_nuc_map[s8mer[0]]
            sites_df.loc[i_row, "seq"] = xm8 + core
        if site in ["6mer", "6mer-bG7"]:
            xA1 = x_nuc_map[s8mer[-1]]
            xm8 = x_nuc_map[s8mer[0]]
            sites_df.loc[i_row, "seq"] = xm8 + core + xA1
        if site == "6mer-A1":
            xm7 = x_nuc_map[s8mer[1]]
            sites_df.loc[i_row, "seq"] = xm7 + core
        if site == "6mer-m8":
            xm2 = x_nuc_map[s8mer[-2]]
            sites_df.loc[i_row, "seq"] =       core + xm2
        if site == "5mer-A1":
            xm6 = x_nuc_map[s8mer[2]]
            sites_df.loc[i_row, "seq"] = xm6 + core
        if site == "5mer-m2.6":
            xA1 = x_nuc_map[s8mer[-1]]
            xm7 = x_nuc_map[s8mer[1]]
            sites_df.loc[i_row, "seq"] = xm7 + core + xA1
        if site == "5mer-m3.7":
            xm2 = x_nuc_map[s8mer[-2]]
            xm8 = x_nuc_map[s8mer[0]]
            sites_df.loc[i_row, "seq"] = xm8 + core + xm2
        if site == "5mer-m8":
            xm3 = x_nuc_map[s8mer[-3]]
            sites_df.loc[i_row, "seq"] =       core + xm3
        # 2.) 3'-only sites of miR-155, miR-124, and lsy-6:
        if "mer-" in site:
            k, lr = site.split("mer-")
            if k in ["9", "10", "11"]:
                # Split schematic: m$(l).$(r)w$(p)
                #                 0|1
                #                   $(l).$(r)w$(p)
                #                           0|1
                #                   $(l).$(r)
                #                      0|1
                l, r = [int(i) for i
                        in lr.split("m")[1].split("w")[0].split(".")]
                _mirna = Mirna(mirna)
                if r < len(_mirna):
                    x_5p = x_nuc_map[get_rc(_mirna.seq[r])]
                else:
                    x_5p = ""
                x_3p = x_nuc_map[get_rc(_mirna.seq[l - 2])]
                site_full = x_5p + core + x_3p
                sites_df.loc[i_row, "seq"] = site_full
        # 3.) AA-[site] sites of miR-124:
        if mirna == "miR-124":
            _mirna = Mirna(mirna)
            # Check x positions at 3p end of AA sites:
            if site[:3] == "AA-":
                site_strip = site[3:]
                if "7mer-m8" in site_strip:
                    x_3p = x_nuc_map[get_rc(_mirna.seq[0])]
                elif "6mer-m8" in site_strip:
                    x_3p = x_nuc_map[get_rc(_mirna.seq[1])]
                elif "5mer-m8" in site_strip:
                    x_3p = x_nuc_map[get_rc(_mirna.seq[2])]
                else:
                    x_3p = ""
                sites_df.loc[i_row, "seq"] = core + x_3p
            # Check that non AA-[site] sites of miR-124 that have a paired
            # AA-[site] are not allowed to have AA in the flanking dinucleotide
            # context.
            elif "AA-" + site in seq_site_mirna_map["miR-124"].keys():
                sites_df.loc[i_row, "seq"] = "Bb" + sites_df.loc[i_row, "seq"]



    # sites_df.to_csv("180930_reporter_sites_check.txt", sep="\t")
    print(sites_df)
    return()
    if len_vars:
        len_vars = int(len_vars)
        suffix = "len%snt_%s" %(len_vars, suffix)
    else:
        len_vars = 180
    each = int(3e4 / total_sites)
    for i_row in range(sites_df.shape[0]):
        mirna_i, name, seq = list(sites_df.iloc[i_row, :])
        if mirna == mirna_i and (not single_site or name == single_site):
            print(name)
            # _mirna = Mirna(mirna)
            SeqsOthersites_df = generate_library_sequence(
                mirna, name, seq, each, freq_MrnaDinuc_map, len_vars=len_vars
            )

            # SeqsOthersites_noproofread_df = generate_library_sequence(
            #     mirna, name, seq, each, freq_MrnaDinuc_map, len_vars=len_vars,
            #     proofread=False
            # )


            MirnaSite_df = pd.DataFrame.from_dict({"miRNA" : [_mirna.name]*each,
                                                  "Site"  : [name]*each})
            out_df_i = pd.concat([MirnaSite_df, SeqsOthersites_df], axis=1)
            out_df_noproofread_i = pd.concat([MirnaSite_df, SeqsOthersites_noproofread_df], axis=1)

            # out_dir = "variants/%s/%s" %(suffix, _mirna.name)
            # ensure_directory(out_dir)
            
            # out_path = "variants/%s/%s/%s_library_variants_%s.txt" %(suffix,
            #                                                   _mirna.name,
            #                                                   name, suffix)
            # out_path_noproofread = "variants/%s/%s/%s_library_variants_%s_noproofread.txt" %(suffix,
            #                                                   _mirna.name,
            #                                                   name, suffix)

            # print(out_path)
            # out_df_i.to_csv(out_path, sep="\t")
            # out_df_noproofread_i.to_csv(out_path_noproofread, sep="\t")

################################################################################


def generate_library_sequence(name, site, n_seq, freq_dinuc_map,
                              proofread=True, max_attempts=10, len_vars=146):
    """Generates library sequences for the reporter library.
    Args:
        site (str): A string giving the DNA sequence corresponding to the miRNA
            target site.
        n_seq (int): The number of sequences requested
    """
    df_cols = ["Variant", "Other sites", "Utr sites"]
    sequences_df = pd.DataFrame(data=None, index=range(n_seq), columns=df_cols)

    splice_donor = r"GGGT[GA]AGT"
    bstxii_site = r"CCA[ACGT]{6,6}TGG"
    polya_site = "AATAAA"
    pumilio_site = r"TGTA[ACGT]ATAA"
                                             #123456789!1234567#         # Length = 17 nt         
    const_5p   =                             "TCTACAGTCCGACGATC"         # October 11th, has the sequencing primer
    rep_5p_utr = "GTACAAATAACCACGCTGGGTTCAGAGTTCTACAGTCCGACGATC"


    const_3p    = "TACCAATGCCCTGGCTC"         # October 11th, 17nt with GC changed to TA.
    # The following sequence is the actual library sequence, owing to the Gibson
    # scheme that was used for the actual library preparation. The UTR differs
    # according to having an extra "TGCCCTGGCTC" in tandem with the 5′-most
    # instance of it.
    #                    TGCCCTGGCTC
    #                    xxxxxxxxxxxTGCCCTGGCTC
    #                               xxxxxxxxxxx
    rep_3p_utr0 = "TACCAATGCCCTGGCTCTGCCCTGGCTCACAAATACCACTGAGATCTTTTTCCCTCTGCCAAAAATTATGGGGACATCATGAAGCCCCTTGAGCATCTGACTTCTGGCT"

    rep_pA  = "AATAAAGGAAATTTATTTTCATTGCA" + "A"*20
    
    rep_3p_utr = rep_3p_utr0 + rep_pA
                 #123456789!1234567#         # Length = 17 nt

    x_nucs_map = {"b" : ["C", "G", "T"],
                  "d" : ["A", "G", "T"],
                  "h" : ["A", "C", "T"],
                  "v" : ["A", "C", "G"]}
    time_start = time.time()
    time_start_i = time_start
    mirna_site = "%s|%s" %(mirna, name)
    min_acheived_other_sites = 2

    for row_i in range(n_seq):
        # print("___________________________________________________________________")
        # print(row_i)
        # time_start = time.time()
        if row_i % 10 == 0 or row_i == 1:
            print("Finished made: %s" %(row_i))
            print("min_site: %s" %(min_acheived_other_sites))
            time_new = time.time()
            print(time_new - time_start_i)
            time_start_i = time_new
        # Establish stopping condition, and reset attempts
        hold = True
        attempts = 0
        trialseqs_df = pd.DataFrame(data=None, index=range(max_attempts),
                                    columns=df_cols)
        while hold:
            sites_found = []
            sites_found_utr = []
            site_list = list(site)
            for i, nuc in enumerate(site_list):
                if nuc == "B":
                    dinuc_list = list(random.choice(non_AA))
                    site_list[i] = dinuc_list[0]
                    site_list[i + 1] = dinuc_list[1]
                elif nuc in x_nucs_map.keys():
                    site_list[i] = random.choice(x_nucs_map[nuc])
            site_new = "".join(site_list)
            len_non_site = len_vars - len(site_new + const_5p + const_3p)
            len_l = len_non_site / 2
            len_r = len_non_site - len_l
            seq_l = generate_utr_sequences(1, len_l, freq_dinuc_map)
            seq_r = generate_utr_sequences(1, len_r, freq_dinuc_map)
            total_seq = const_5p + seq_l + site_new + seq_r + const_3p
            utr_seq   = rep_5p_utr + seq_l + site_new + seq_r + rep_3p_utr
            _reporter_variant = ReporterVariant(total_seq)
            _utr_variant = ReporterVariant(utr_seq)
            # Check variant for other miRNA sequences
            for mirna_i in MIRNAS:
                _mirna = Mirna(mirna_i)
                _sitelist = SiteList(_mirna, "paperfinal", len_vars)
                _reporter_variant.get_all_sites_in_reporter_variant(_sitelist, ties=True)
                _utr_variant.get_all_sites_in_reporter_variant(_sitelist, ties=True)
                site_names =["%s|%s" %(_mirna.name, i.name) for i in _reporter_variant.sites]
                site_names_utr =["%s|%s" %(_mirna.name, i.name) for i in _utr_variant.sites]
                if site_names != ["%s|None" %(_mirna.name)]:
                    lr = ["%s.%s" %(i.l, i.r) for i in _reporter_variant.sites]
                    site_names_lr = ["|".join(i) for i in zip(site_names, lr)]
                    sites_found += site_names_lr
                if site_names_utr != ["%s|None" %(_mirna.name)]:
                    lr = ["%s.%s" %(i.l, i.r) for i in _utr_variant.sites]
                    site_names_lr = ["|".join(i) for i in zip(site_names_utr, lr)]
                    sites_found_utr += site_names_lr

            # Successful stopping condition:
            span = re.search(bstxii_site, total_seq).span()
            l_bound, r_bound = list(map(add, span, [1, -1]))
            l_seq_check = total_seq[:r_bound]
            r_seq_check = total_seq[l_bound:]
            site_span_r = len(const_5p) + len(seq_l) + len(site_new) - 2
            site_span_l = len(const_5p) + len(seq_l) + 1
            alt_polyA_split_l = total_seq[:site_span_r]
            alt_polyA_split_r = total_seq[site_span_l:]
            # print(sites_found[0])
            # print(sites_found[0].split("|")[1])
            # print(name)
            if ( not proofread or (
                    len(sites_found) == 1 and
                    sites_found[0].split("|")[1] == name and    
                    (polya_site not in total_seq or (
                        name == "AATAAAG" and
                        polya_site not in alt_polyA_split_l and
                        polya_site not in alt_polyA_split_r)) and
                    not re.search(splice_donor, total_seq) and
                    not re.search(pumilio_site, total_seq) and
                    not re.search(bstxii_site, l_seq_check) and
                    not re.search(bstxii_site, r_seq_check)
                                  )
               ):
                hold = False
                min_acheived_other_sites = 1
                if len(sites_found) == 1:
                    sites_found_str = sites_found[0]
                else:
                    sites_found_str = ",".join(sites_found)

                sites_found_utr_str = ",".join(sites_found_utr)
                sequence_row = pd.Series([total_seq, sites_found_str,
                                          sites_found_utr_str],
                                            index=df_cols)
                sequences_df.loc[row_i, ] = sequence_row 
            # Stopping condition at which point the best sequence is picked.
            elif attempts == max_attempts:
                print("at max_attempts")
                other_sites = list(trialseqs_df["Other sites"])
                # First filter by having the fewest number of sites:
                sites = [i.split(",") for i in other_sites]
                len_sites = [len(i) for i in sites]
                if (
                        min(len_sites) != min_acheived_other_sites or
                        (polya_site in total_seq and name != "AATAAAG") or
                        (name == "AATAAAG" and
                        (polya_site in alt_polyA_split_l or
                         polya_site in alt_polyA_split_r)) and

                        re.search(splice_donor, total_seq) or
                        re.search(pumilio_site, total_seq) or
                        re.search(bstxii_site, l_seq_check) or
                        re.search(bstxii_site, r_seq_check)
                   ):
                    attempts = 0
                else:
                    min_sites = [i for i in range(len(len_sites)) if len_sites[i] == 2]
                    trialseqs_df = trialseqs_df.iloc[min_sites, ]
                    other_sites = list(trialseqs_df["Other sites"])

                    mirna_sites = [["|".join(j.split("|")[:2]) for j in i.split(",")] for i in other_sites]
                    mirna_sites_other = [[j for j in i if j != mirna_site] for i in mirna_sites]
                    mirna_other_index_1 = [i for i, mir_site_i in enumerate(mirna_sites_other) if mir_site_i != []]
                    mirna_other_index = [i for i in mirna_other_index_1 if mirna_sites_other[i][0].split("|")[0] != mirna]
                    trialseqs_df = trialseqs_df.iloc[mirna_other_index, ]
                    mirna_sites_other_new = [mirna_sites_other[i][0] for i in mirna_other_index]
                    site_lengths = [len(Mirna(i.split("|")[0])[i.split("|")[1]]) for i in mirna_sites_other_new]
                    # print(site_lengths)
                    len_index = [i for i, len_i in enumerate(site_lengths) if len_i == min(site_lengths)][0]
                    # print(len_index)
                    # print(trialseqs_df.shape)
                    sequence_row = trialseqs_df.iloc[len_index,]
                    # print(sequence_row)
                    sequences_df.loc[row_i,] = sequence_row
                    hold = False

            elif ( not proofread or (
                    len(sites_found) == 2 and
                    (polya_site not in total_seq or (
                        name == "AATAAAG" and
                        polya_site not in alt_polyA_split_l and
                        polya_site not in alt_polyA_split_r)) and
                    not re.search(splice_donor, total_seq) and
                    not re.search(pumilio_site, total_seq) and
                    not re.search(bstxii_site, l_seq_check) and
                    not re.search(bstxii_site, r_seq_check)
                                    )
                 ):

                # print(total_seq)
                sites_found_str = ",".join(sites_found)
                # print(sites_found_str)
                sites_found_utr_str = ",".join(sites_found_utr)

                row_df = pd.Series([total_seq, sites_found_str, sites_found_utr_str],
                                   index=df_cols)
                # print("row_df:")
                # print(row_df)
                # print("trialseqs_df")
                # print(trialseqs_df)
                # print(attempts)
                trialseqs_df.loc[attempts,] = row_df
                attempts += 1
                # print("trialseqs_df post:")
                # print(trialseqs_df)
        # print(sequences_df)
    time_stop = time.time()
    time_del = time_stop - time_start
    time_min = int(time_del / 60)
    time_sec = time_min
    print(time_del)

    return(sequences_df)


if __name__ == "__main__":
    main()

