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



#FINAL DESIGN CONSTRAINTS
kLEN = 146

# FUNCTIONS ####################################################################
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

def generate_utr_sequence(len_seq, frac_dinuc_dict):
    """Generates a DNA sequences using the provided dinucleotide frequencies.

    Args:
        utrs (list): A list of strings, containing only the characters
            A, C, G, and U.
    Returns:
        dict: A dictionary with keys that are the 16 DNA dinucleotides and the 
            keys are the corresponding frequencies within the entire list of
            sequences. The dictionary sums to 1.0.
    """
    # Define the starting output list.
    a = list(frac_dinuc_dict.keys())
    p = list(frac_dinuc_dict.values())
    # Initialize the sequence with a dinucleotide sequence.
    utr = np.random.choice(a, p=p)
    # Iterate over the rest of the sequences.
    for j in range(len_seq - 2):
        # Take only those dinucleotides compatible with the last nucleotide of
        # the current sequence.
        all_use = [k for k in zip(a, p) if k[0][0] == utr[-1]]
        # Makes sure that only A and T will be inserted at positions 9–13 in
        # flanking contexts.
        if (46 <= j < 49 or 92 <= j < 96):
            all_use = [k for k in all_use if k[0][1] in ["A", "T"]]
        a_use = [k[0] for k in all_use]
        p_use = [k[1] for k in all_use]
        p_use = [i/sum(p_use) for i in p_use]
        dinuc = np.random.choice(a_use, p=p_use)
        utr += dinuc[1]
    return utr



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
    # Get the length of prior sequence for positional information.
    seq_len_orig = len(seq)
    for pos in range(len_seq):
        # Take only those dinucleotides compatible with the last nucleotide of
        # the current sequence.
        ap_use = [i for i in zip(a, p) if i[0][0] == seq[-1]]
        # Makes sure that only A and T will be inserted at positions 9–13 in
        # flanking contexts.
        if (44 < pos + seq_len_orig - 2 < 49 or 91 < pos + seq_len_orig - 2 < 96):
            ap_use = [i for i in ap_use if i[0][1] in ["A", "T"]]
        a_use = [i[0] for i in ap_use]
        p_use = [i[1] for i in ap_use]
        dinuc = np.random.choice(a_use, p=p_use/np.sum(p_use))
        seq += dinuc[1]
    return seq

def check_sequence_for_sites(seq, _mirna, site, offset, position, seed_pos,
                             len_const_5p, print_pairing=False):
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
    # 1. Define the relevant sequences, being the full let-7 sequence, and the 
    #    miR-1 seed.
    mir_seq = get_rc(_mirna.seq[1:]) + "A"
    mir1_seed_seq = Mirna("miR-1")["8mer"]

    mir_seqs = [mir_seq, mir1_seed_seq]
    char = ["f", "1"]

    pos_pairing = []
    full_pairing = [" "]*len(seq)
    for i in range(2):
        # First define the maximum overlap of the sequence with the substring.
        [len_comp, l_mirAndThrPs] = LCSubString(mir_seqs[i], seq)
        # If the length is larger than the maximum allowed limit for that type
        # of sequence, iterate from this length to the minimum allowed limit,
        # marking each of the sequences with pairing.
        if len_comp >= 6:
            for len_i in range(len_comp, 5, -1):
                [len_comp, l_mirAndThrPs] = LCSubString(mir_seqs[i], seq, alt_len=len_i)
                for j in range(len(l_mirAndThrPs[0])):
                    mir_pos = l_mirAndThrPs[0][j]
                    seq_pos = l_mirAndThrPs[1][j]
                    pairing_range = [k for k in range(seq_pos, seq_pos + len_comp)]
                    for k in pairing_range:
                        if full_pairing[k] in ["F", "f"] and char[i] == "1":
                            full_pairing[k] = "F"
                        else:
                            full_pairing[k] = char[i]
                    pos_pairing += pairing_range

    # 2. Define the let-7a 3-prime region, to check for 4-nt of complementarity
    #    anywhere within the 22 nt stargint at position 1 of the miRNA.
    # Define miRNA 3-prime region.
    mir_thrp_seq = get_rc(_mirna.seq[8:])
    # Calculate the shift for both the window of sequence checked, and the
    # right-hand reference point owing to having bulged nucleotides or the
    # bulged seed site in lin-41.
    shift_l = len(offset)*(position in ["left", "dual"]) + (site == "lin_41")
    shift_r = len(offset)*(position in ["right", "dual"])
    win_l = 22 + shift_l
    win_r = 22 + shift_r
    seed_l = len_const_5p + seed_pos[0] + shift_l + 1
    seed_r = len_const_5p + seed_pos[1] + shift_l + shift_r + 1
    win_lr = [win_l, win_r]
    seed_lr = [seed_l, seed_r]
    # Determine which sites to check for.
    if position == "dual":
        thrp_check = [0, 1]
    elif position == "left":
        thrp_check = [0]
    elif position == "right":
        thrp_check = [1]
    else:
        thrp_check = []
    # Iterate over sites.
    for i in thrp_check:
        seq_use = seq[seed_lr[i] - win_lr[i]:seed_lr[i]]
        # Get max complementarity.
        len_comp = LCSubString(mir_thrp_seq, seq_use)[0]
        # If longer than 4, check all lengths between this length and 4.
        if len_comp >= 4:
            for len_i in range(len_comp, 3, -1):
                [len_comp, (mir_ls, seq_ls)] = LCSubString(mir_thrp_seq, seq_use,
                                                        alt_len=len_i)
                for seq_l in seq_ls:
                    seq_pos = seq_l + seed_lr[i] - win_lr[i]
                    pairing_range = [k for k in range(seq_pos,
                                                      seq_pos + len_comp)]
                    pos_pairing += pairing_range     
                    for k in pairing_range:
                        if full_pairing[k] in ["f", "F", "T", "1"]:
                            full_pairing[k] = "T"
                        else:
                            full_pairing[k] = "t"

    print("".join(full_pairing))


    # Important sequences from prior library generation.
    splice_donor = r"GGGT[GA]AGT"
    bstxii_site = r"CCA[ACGT]{6,6}TGG"
    bsrgi_site = r"TGTACA"
    polya_site = "AATAAA"
    pumilio_site = r"TGTA[ACGT]ATAA"

    non_mir_seqs = [splice_donor, bstxii_site, bsrgi_site,
                    polya_site, pumilio_site]

    for non_mir_seq in non_mir_seqs:
        res = re.search(non_mir_seq, seq)
        if res:
            pairing_range = [k for k in range(res.span()[0], res.span()[1])]
            for k in pairing_range:
                if full_pairing[k] != " " and full_pairing[k] != "x":
                    full_pairing[k] = "X"
                else:
                    full_pairing[k] = "x"
            pos_pairing += pairing_range
    # Remove redundant entries and order `pos_pairing`.
    pos_pairing = list(set(pos_pairing))
    # Sort the pairing entries.
    pos_pairing.sort()
    # Make `full_pairing` diagnostic string.
    full_pairing = "".join(full_pairing)
    if print_pairing:
        print(seq)
        print(full_pairing)
    return(pos_pairing)


def mutate_seq(seq, site, allowed_pos, paired_pos, offset, position, seed_pos,
               new_method=False):
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
    # First get the indeces of all the positions with nucleotides that form
    # undesired complementarity.
    dif_seq = [i for i in paired_pos if i not in allowed_pos]
    # Correct for the offset of the left-hand site.
    if position in ["left", "dual"]:
        for i in range(len(dif_seq)):
            pos = dif_seq[i]
            if pos > seed_pos[0] - 9:
                dif_seq[i] = pos - len(offset)
    if site == "lin-41":
        for i in range(len(dif_seq)):
            pos = dif_seq[i]
            if pos > seed_pos[0] - 6:
                dif_seq[i] = pos - 1
    if position in ["right", "dual"]:
        for i in range(len(dif_seq)):
            pos = dif_seq[i]
            if pos > seed_pos[1] - 9:
                dif_seq[i] = pos - len(offset)
    seq_list = list(seq)
    contig_nucs = []
    dif_seq += [dif_seq[-1] + 2]
    for i in dif_seq:
        if len(contig_nucs) == 0:
            contig_nucs.append(i)
        elif i == contig_nucs[-1] + 1:
            contig_nucs.append(i)
        else:
            # Mutate the central nucleotide within this contiguous segment.
            pos_mut = math.ceil(sum(contig_nucs)/len(contig_nucs))
            nuc = seq_list[pos_mut]
            print(" "*pos_mut + nuc)
            if (47 <= pos_mut < 51 or 94 <= pos_mut < 98):
                mut_nuc = np.random.choice([i for i in ["A", "T"] if i != nuc])
            else:
                mut_nuc = np.random.choice([i for i in ["A", "C", "G", "T"] if i != nuc])
            seq_list[pos_mut] = mut_nuc
            # Re-assign the contiguous set of nucleotides.
            contig_nucs = [i]
    # Form the new sequence from the list.
    seq_new = "".join(seq_list)
    return(seq_new)

def add_sites_to_seq(seq, _mirna, seed, thrp, offset, position, seed_pos,
                     print_pairing=False):
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
        print_pairing (bool): A boolean indicating whether or not to print
            strings reporting on the positions at which the sequence was
            changed.
    """
    # 1a. Get the seed site from the miRNA object.
    if seed == "lin-41":
        seed_site = _mirna["8mer-bA5"]
        seed_site_2 = _mirna["7mer-m8w6"]
    else:
        seed_site = _mirna[seed]
    # 1b. Get the left- and right-hand positions of where to replace the seed
    #     site relative to the position of the A1.
    if seed == "lin-41": # Required for lin-41 due to different seed sites.
        seed_l = -8
        seed_r = 0
        seed_r_2 = -1
    elif seed == "No_site_0": # Doesn't matter for no-sites.
        seed_l = 0
        seed_r = 0
    else:
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
        thrp_site = "" # Doesn't matter for no-sites.
        thrp_l = 0
        thrp_r = 0
    paired_nucs = []
    # Add the seed, threep, and offset sequences into the reporter.
    if position in ["left", "dual"]:
        # Add the seed pairing to the sequence and record the positions.
        seed_0 = seed_pos[0]
        seq = seq[:seed_0 + seed_l] + seed_site + seq[seed_0 + seed_r:]
        paired_nucs += [
            i for i in range(seed_0 + seed_l + len(offset),
                             seed_0 + seed_r + len(offset) + (seed == "lin-41"))
        ]
        # Add the 3-prime pairing to the sequence and record the positions.
        seq = seq[:seed_0 + thrp_l] + thrp_site + seq[seed_0 + thrp_r:]
        # Record the paired positions, adding the offset owing to the right-hand
        # offset length. Recording with the offset reflects the offset yet that
        # has not yet been added, and happens in 9 lines.
        paired_nucs += [i for i in range(seed_0 + thrp_l, seed_0 + thrp_r)]
        # Add the offset nucleotides into the sequence.
        seq = seq[:seed_0 - 9] + offset + seq[seed_0 - 9:]
        # Record the offset sequence.
        paired_nucs += [i for i in range(seed_0 - 9, seed_0 - 9 + len(offset))]

    if seed == "lin-41":
        seed_site = seed_site_2
        seed_r = seed_r_2
    if position in ["right", "dual"]:
        # First define the right-hand A1 position.
                 # 1. Base right-hand seed position 1
                               # 2. Whatever the offset length is X
                                           # 3. Whether or not the site is a
                                           #    dual site
                                                                  # 4. One extra nucleotide
                                                                  #    owing to the bulged
                                                                  #    seed in the left-hand
                                                                  #    lin-41 site.
        seed_0 = seed_pos[1] + len(offset)*(position == "dual") + (seed == "lin-41")
        # Modify the sequence to contain the seed pairing.
        seq = seq[:seed_0 + seed_l] + seed_site + seq[seed_0 + seed_r:]
        # Record the paired positions, adding the offset owing to the right-hand
        # offset length. Recording with the offset reflects the offset yet that
        # has not yet been added, and happens in 9 lines.
        paired_nucs += [
            i for i in range(seed_0 + seed_l + len(offset),
                             seed_0 + seed_r + len(offset))
        ]
        # Add the 3-prime pairing to the sequence and record the positions.
        seq = seq[:seed_0 + thrp_l] + thrp_site + seq[seed_0 + thrp_r:]
        # Record the paired positions with respect to the 3-prime end. Here the
        # right-hand offset doesn't need to be counted, since the 3-prime
        # happens 5-prime of the offset.
        paired_nucs += [i for i in range(seed_0 + thrp_l, seed_0 + thrp_r)]
        # Add the offset nucleotides into the sequence.
        seq = seq[:seed_0 - 9] + offset + seq[seed_0 - 9:]
        # Record the offset sequence.
        paired_nucs += [i for i in range(seed_0 - 9, seed_0 - 9 + len(offset))]
    paired_nucs.sort()
    buffer_seq = [" "]*len(seq)
    for i in paired_nucs:
        buffer_seq[i] = "x"
    buffer_seq = "".join(buffer_seq)
    if print_pairing:
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
    # Make a map of the dinucleotide frequencies from human 3-prime UTRs.
    freq_MrnaDinuc_map = count_utr_sequence(utr_df["sequence"])

    ## 2: LOAD LIN-41 SEQUENCE TO MODEL DUAL-SITE ARCHITECTURE #################
    # Load the lin-41 3-prime utr, taken from wormbase, to parse so as to obtain
    # a background context to use with the other site architectures.
    with open("ThreePrimePaperLet7ReporterAssay/lin41_3pUTR.txt") as file:
        l41_utr = file.read().strip().upper()
    # Using meta-data, get the sequences of the two seed sites and the shared
    # 3p site that is found in the lin41 3 prime utr.
    l41_seed_sites = [_mirna[i] for i in ["8mer-bA5", "7mer-m8w6"]]
    l41_3p_site = _mirna["9mer-m11.19"]
    # Define the left-hand and right-hand limits of the lin41 site.
    l41_site_l = find_all(l41_utr, l41_3p_site)[0]
    l41_site_r = max([l41_utr.find(i) + len(i) for i in l41_seed_sites])
    # Get the length of the lin-41 site.
    l41_site_len = l41_site_r - l41_site_l
    # The added nucleotides are the length of the plasmid fragment minus the
    # the full length of the lin41 dual site.
    N_lr = kLEN + 3 - l41_site_len
    N_l = math.floor(N_lr/2)
    N_r = math.ceil(N_lr/2)
    # Define the lin41 reporter sequence to be used.
    l41_rep_seq = l41_utr[l41_site_l - N_l:l41_site_r + N_r]

    # Calculate the length of the bridge between the two dual sites.
    seed_site_r = [l41_rep_seq.find(i) + len(i) for i in l41_seed_sites]
    thrp_site_l = find_all(l41_rep_seq, l41_3p_site)
    n_bridge = thrp_site_l[1] - seed_site_r[0]
    # Define the limits of the let7a site type:
    # DIST          N_l       11          3        9
    # FEATURE     left flank | 3p site 1 | loop 1 | seed 1 | 
    # DIST            28          11          3       8*        N_r - 1*
    # FEATURE     intersite  | 3p site 2 | loop 2 | seed 2 |  right flank
    # * The reason the N_r is one nucleotide shorter is that the seed2 site
    # is actually a 7mer-m8, not an 8mer, but I am including the non-A
    # nucleotide from the intervening sequence as part of the site.
    # Define the limits of the lin41 sequence:
    l41_l_thrp1 = N_l
    l41_r_thrp1 = l41_l_thrp1 + 9
    l41_l_seed1 = l41_r_thrp1 + 3               # This has an offset of 1
    l41_r_seed1 = l41_l_seed1 + 9               # This has 1 bulged nucleotide
    l41_l_thrp2 = l41_r_seed1 + n_bridge
    l41_r_thrp2 = l41_l_thrp2 + 9
    l41_l_seed2 = l41_r_thrp2 + 3               # This has an offset of 1
    l41_r_seed2 = l41_l_seed2 + 7

    # Get the starting positions for the seed A1 target nucleotides in a
    # non lin-41 site (i.e., one with no offsets and no bulged target site for
    # its first site).
    rep_r_seed1 = l41_r_seed1 - 2 # 2 because of +1 offset and seed bulge.
    rep_r_seed2 = l41_r_seed2 - 2 # 2 = +3 (because of +1 offset, bulge, and +1 offset)
                                  #     -1 (because this is a 7mer-m8-like site, so starts at 2)

    l_seq     = l41_rep_seq[            : l41_l_thrp1]
    thrp1_seq = l41_rep_seq[l41_l_thrp1 : l41_r_thrp1]
    loop1_seq = l41_rep_seq[l41_r_thrp1 : l41_l_seed1]
    seed1_seq = l41_rep_seq[l41_l_seed1 : l41_r_seed1]
    inter_seq = l41_rep_seq[l41_r_seed1 : l41_l_thrp2]
    thrp2_seq = l41_rep_seq[l41_l_thrp2 : l41_r_thrp2]
    loop2_seq = l41_rep_seq[l41_r_thrp2 : l41_l_seed2]
    seed2_seq = l41_rep_seq[l41_l_seed2 : l41_r_seed2]
    r_seq     = l41_rep_seq[l41_r_seed2 :            ]

    # Check that the addition of the parts recapitulates the original site.
    bi1_seq = thrp1_seq + loop1_seq + seed1_seq
    bi2_seq = thrp2_seq + loop2_seq + seed2_seq
    full_seq_check = l_seq + bi1_seq + inter_seq + bi2_seq + r_seq
    print(l41_rep_seq == full_seq_check)

    # Make the lin-41 sequence context, by replacing the entire sequence of the
    # two dual sites with sequence sampled from human 3-prime UTR dinucleotide
    # frequencies.
    # 1. Initialize as left-hand flanking sequence, replacing the first
    #    3-prime site.
    l41_bg_context = add_sequence(l_seq, 9, freq_MrnaDinuc_map)
    # 2. Add the loop sequence.
    l41_bg_context += (loop1_seq[0] + loop1_seq[2])
    # 3. Add sequence replacing the first seed site.
    l41_bg_context = add_sequence(l41_bg_context, 8, freq_MrnaDinuc_map)  
    # 4. Add the inter-sequence sequence.
    l41_bg_context += inter_seq
    # 5. Add sequence replacing the second 3-prime site.
    l41_bg_context = add_sequence(l41_bg_context, 9, freq_MrnaDinuc_map)
    # 6. Add the second loop sequence.
    l41_bg_context += (loop2_seq[0] + loop2_seq[2])
    # 7. Add the second seed sequence.
    l41_bg_context = add_sequence(l41_bg_context, 7, freq_MrnaDinuc_map)
    # 8. Add the right-hand flanking sequence.
    l41_bg_context += r_seq

    print("ACCTCTTTTCCTCAAATTGCACCAACTCAAGTATACCTTTCTCAATT--GTGTGAAAAGACGCGATGTAAATATCGCAATCCCTTTTACAATTC--ATTGCTCCATGAACCATTGAAACCTTCTCCCGTACTCCCACCAATAGATT")

                                             #123456789!1234567#         # Length = 17 nt         
    const_5p   =                             "TCTACAGTCCGACGATC"         # October 11th, has the sequencing primer
    rep_5p_utr =  "GTACAAATAACCACGCTGGGTTCAGAGTTCTACAGTCCGACGATC"

    const_3p    = "TACCAATGCCCTGGCTC"         # October 11th, 17nt with GC changed to TA.
    # PRIMERS ##################################################################
    # Next are the original gibson 5-prime primer, as well as a new gibson
    # 5-prime primer. The new 5-prime primer includes the BsrGI site as well as
    # six nucleotides further 5-prime of the sequence, which is thought to
    # maximize efficient cleavage by the restriction enzyme.
    gibson_pool_fwd_old =    "CTGTACAAATAACCACGCTGGGTTCAGAGTTCTACAGTCCGACGATC"
    # 210630 THESE FOUR PRIMERS ARE THE PRIMERS TO USE FOR THE LIBRARY DESIGN.
    # They allow 36- and 37-nt of overlap 5p and 3p of the pool, when
    # performing Gibson assembly, as well as allowing restriction-based cloning
    # using BsrGI and BstXI to make the restriction fragments at the 5p and 3p
    # ends of the pool.
    gibson_pool_fwd = "ACGAGCTGTACAAATAACCACGCTGGGTTCAGAGTTCTACAGTCCGACGATC"
    gibson_pool_rev = "GGTATTTGTGAGCCAGGGCAGAGCCAGGGCATTGGTA"
    # For the primer to amplify 
    gibson_plasmid_fwd = "TACCAATGCCCTGGCTCTGCCCTGGCTCACAAATACC"
    gibson_plasmid_rev = "AACTCTGAACCCAGCGTGGTTATTTGTACAGCTCGT"
    # THIS CONTENT IS FROM WHEN GENERATING THE ORIGINAL PLASMID LIBRARY, AND SO
    # REFLECTS HOW THE STATE OF THE EXTANT PLASMID LIBRARY BEING STORED BY ASIA.
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

    site_list = [[   "lin-41", "9mer-m11.19",    "T",  "dual"],
                 [     "8mer",            "",     "",  "left"],
                 [     "8mer",            "",     "", "right"],
                 [     "8mer",            "",     "",  "dual"],
                 [  "7mer-m8",            "",     "",  "left"],
                 [  "7mer-m8",            "",     "", "right"],
                 [  "7mer-m8",            "",     "",  "dual"],
                 [  "7mer-A1",            "",     "",  "left"],
                 [  "7mer-A1",            "",     "", "right"],
                 [  "7mer-A1",            "",     "",  "dual"],
                 [     "6mer",            "",     "",  "left"],
                 [     "6mer",            "",     "", "right"],
                 [     "6mer",            "",     "",  "dual"],
                 [  "8mer-w6",            "",     "",  "left"],
                 [  "8mer-w6",            "",     "", "right"],
                 [  "8mer-w6",            "",     "",  "dual"],
                 [  "8mer-w6", "4mer-m13.16",     "",  "left"],
                 [  "8mer-w6", "4mer-m13.16",     "", "right"],
                 [  "8mer-w6", "4mer-m13.16",     "",  "dual"],
                 [  "8mer-w6", "4mer-m13.16",    "A",  "left"],
                 [  "8mer-w6", "4mer-m13.16",    "A", "right"],
                 [  "8mer-w6", "4mer-m13.16",    "A",  "dual"],
                 [  "8mer-w6", "4mer-m13.16", "AAAA",  "left"],
                 [  "8mer-w6", "4mer-m13.16", "AAAA", "right"],
                 [  "8mer-w6", "4mer-m13.16", "AAAA",  "dual"],
                 [  "8mer-w6", "4mer-m13.16",    "T",  "left"],
                 [  "8mer-w6", "4mer-m13.16",    "T", "right"],
                 [  "8mer-w6", "4mer-m13.16",    "T",  "dual"],
                 [  "8mer-w6", "4mer-m13.16", "TTTT",  "left"],
                 [  "8mer-w6", "4mer-m13.16", "TTTT", "right"],
                 [  "8mer-w6", "4mer-m13.16", "TTTT",  "dual"],
                 [  "8mer-w6", "9mer-m11.19",     "",  "left"],
                 [  "8mer-w6", "9mer-m11.19",     "", "right"],
                 [  "8mer-w6", "9mer-m11.19",     "",  "dual"],
                 [  "8mer-w6", "9mer-m11.19",    "A",  "left"],
                 [  "8mer-w6", "9mer-m11.19",    "A", "right"],
                 [  "8mer-w6", "9mer-m11.19",    "A",  "dual"],
                 [  "8mer-w6", "9mer-m11.19", "AAAA",  "left"],
                 [  "8mer-w6", "9mer-m11.19", "AAAA", "right"],
                 [  "8mer-w6", "9mer-m11.19", "AAAA",  "dual"],
                 [  "8mer-w6", "9mer-m11.19",    "T",  "left"],
                 [  "8mer-w6", "9mer-m11.19",    "T", "right"],
                 [  "8mer-w6", "9mer-m11.19",    "T",  "dual"],
                 [  "8mer-w6", "9mer-m11.19", "TTTT",  "left"],
                 [  "8mer-w6", "9mer-m11.19", "TTTT", "right"],
                 [  "8mer-w6", "9mer-m11.19", "TTTT",  "dual"],
                 [  "8mer-w6", "9mer-m13.21",     "",  "left"],
                 [  "8mer-w6", "9mer-m13.21",     "", "right"],
                 [  "8mer-w6", "9mer-m13.21",     "",  "dual"],
                 [  "8mer-w6", "9mer-m13.21",    "A",  "left"],
                 [  "8mer-w6", "9mer-m13.21",    "A", "right"],
                 [  "8mer-w6", "9mer-m13.21",    "A",  "dual"],
                 [  "8mer-w6", "9mer-m13.21", "AAAA",  "left"],
                 [  "8mer-w6", "9mer-m13.21", "AAAA", "right"],
                 [  "8mer-w6", "9mer-m13.21", "AAAA",  "dual"],
                 [  "8mer-w6", "9mer-m13.21",    "T",  "left"],
                 [  "8mer-w6", "9mer-m13.21",    "T", "right"],
                 [  "8mer-w6", "9mer-m13.21",    "T",  "dual"],
                 [  "8mer-w6", "9mer-m13.21", "TTTT",  "left"],
                 [  "8mer-w6", "9mer-m13.21", "TTTT", "right"],
                 [  "8mer-w6", "9mer-m13.21", "TTTT",  "dual"],
                 ["No_site_0",            "",     "",      ""]]

    # site_list = [[   "lin-41", "9mer-m11.19",    "T",  "dual"],
    #              [     "8mer",            "",     "",  "left"]]

    # Modification made on 210726 to make 100 different versions of the library,
    # so that Peter/Thy can check each of them.
    for lib_i in range(100):
        variants_all = []
        for i_context in range(14):
            print("Context %s #######################################" %i_context)
            if i_context == 0:
                seq = l41_bg_context
            else:
                seq = generate_utr_sequence(146, freq_MrnaDinuc_map)

            # Pre-allocate "complete" to make sure each variant gets to the end.
            complete = False
            attempts = 0
            seq_0 = seq
            while complete == False:
                variants = []
                complete = True
                for site_type_i in site_list:
                    site, thrp, offset, position = site_type_i
                    if attempts == 0:
                        print_pairing = True
                    else:
                        print_pairing = False
                        print_pairing = True
                    print("adding site")
                    # Add site architecture into flanking sequence.
                    print_pairing = True
                    variant, allowed_pos = add_sites_to_seq(
                        seq, _mirna, site, thrp, offset, position,
                        [rep_r_seed1, rep_r_seed2], print_pairing=print_pairing
                    )
                    # Define the limits of the let7a site type:
                    # DIST          N_l       9          3        9
                    # FEATURE     left flank | 3p site 1 | loop 1 | seed 1 | 
                    # DIST            28          9          3       8*        N_r - 1*
                    # FEATURE     intersite  | 3p site 2 | loop 2 | seed 2 |  right flank
                    if i_context == 0:
                        left_context = [i for i in range(N_l)]                      # Left Hand
                        loop_1l = [N_l +  9]
                        loop_1r = [N_l + 10 + len(offset)*(position in ["left", "dual"])]       # Loop
                        intersite = [i + len(offset)*(position in ["left", "dual"]) + (site == "lin-41")
                                     for i in range(N_l + 19 , N_l + 19 + 28)] # Intersite
                        loop_2l = [N_l + 56 + len(offset)*(position in ["left", "dual"]) + (site == "lin-41")]
                        loop_2r = [N_l + 57 + len(offset)*(position in ["left", "dual"]) + (site == "lin-41") + len(offset)*(position in ["right", "dual"])]
                        right_context = [i + len(offset)*(position in ["left", "dual"]) + (site == "lin-41") + len(offset)*(position in ["right", "dual"])
                                         for i in range(N_l + 65, len(seq))]
                        extra_pos = left_context + loop_1l + loop_1r + intersite + loop_2l + loop_2r + right_context
                    print("checking site")
                    temp = [" "]*len(variant)
                    for pos_i in extra_pos:
                        temp[pos_i] = "x"
                    print("".join(temp))
                    print(variant)
                    paired_pos = check_sequence_for_sites(const_5p + variant + const_3p, _mirna, site, offset, position, [rep_r_seed1, rep_r_seed2], len(const_5p),
                                                          print_pairing=print_pairing)
                    # Code correcting for the addition of the constant sequence.
                    paired_pos = [i - len(const_5p) for i in paired_pos]
                    paired_pos = [i for i in paired_pos if i >= 0 if i < len(variant)]
                    if i_context == 0:
                        allowed_pos += extra_pos
                        allowed_pos =  list(set(allowed_pos))
                        allowed_pos.sort()
                    union_list = list(set(paired_pos) | set(allowed_pos))
                    union_list.sort()
                    if union_list != allowed_pos:
                        if i_context == 0:
                            print(i_context)
                            print(site_type_i)
                            print("union list")
                            print(union_list)
                            print("paired_pos")
                            print(paired_pos)
                            print("allowed_pos")
                            print(allowed_pos)
                        complete = False
                        if attempts == 10:
                            print("Try a new sequence.")
                            if i_context == 0:
                                l_seq_context = add_sequence(l_seq, 9, freq_MrnaDinuc_map)
                                l_seq_context += (loop1_seq[0] + loop1_seq[2])
                                l_seq_context = add_sequence(l_seq_context, 8, freq_MrnaDinuc_map)  

                                inter_context = add_sequence(inter_seq, 9, freq_MrnaDinuc_map)
                                inter_context += (loop2_seq[0] + loop2_seq[2])
                                inter_context = add_sequence(inter_context, 7, freq_MrnaDinuc_map)

                                l41_bg_context_new = (l_seq_context + inter_context + r_seq)

                                seq = l41_bg_context_new
                            else:
                                seq = generate_utr_sequence(146, freq_MrnaDinuc_map)
                            attempts = 0
                        else:
                            print("Mutate sequence attempt %s; %s failed." %(attempts, site + thrp + offset + position))
                            seq_new = mutate_seq(seq, site, allowed_pos, paired_pos, offset, position, [rep_r_seed1, rep_r_seed2])
                            if i_context == 0 and site == "lin-41":
                                print(seq_0)
                                buffer = [" "]*len(seq)
                                for j in range(len(seq)):
                                    if seq_0[j] != seq_new[j]:
                                        buffer[j] = "x"
                                buffer_str = "".join(buffer)
                                print(buffer_str)
                                print(seq_new)

                            seq = seq_new
                            attempts +=1
                        break
                    variants.append(const_5p + variant + const_3p)
            # Part where 6 more no-site reporters are made per contex.
            for j in range(6):
                seq_l = seq[            : rep_r_seed1 - 20]
                seq_m = seq[rep_r_seed1 : rep_r_seed2 - 20]
                seq_r = seq[rep_r_seed2 : ]
                complete = False
                while complete == False:
                    complete = True
                    seq_l_ext = add_sequence(seq_l, 20, freq_MrnaDinuc_map)
                    seq_m_ext = add_sequence(seq_m, 20, freq_MrnaDinuc_map)
                    variant = seq_l_ext + seq_m_ext + seq_r
                    paired_pos = check_sequence_for_sites(
                        variant, _mirna, "", "", "", [rep_r_seed1, rep_r_seed2], len(const_5p),
                    )
                    if i_context == 0:
                        paired_pos == 0
                    if len(paired_pos) != 0:
                        complete = False
                variants.append(const_5p + variant + const_3p)
            variants_all += variants
        # Finished loop.
        # Add extra constant sequence so that each RNA is the length of the longest
        # current library member.
        lens_all = [len(i) for i in variants_all]
        max_len = max(lens_all)
        for i in range(len(variants_all)):
            variants_all[i] = variants_all[i] + "TGCCCTGGCTCACA"[:max_len - lens_all[i]]
        ########## PRINT COMPARISON OF LIN-41 REPORTER SITE TO THE #################
        ########## ENDOGENOUS SITE #################################################
        # Compare lin-41 site architecture to the lin-41 site being used in the
        # reporter.
        l41_rep_seq_use = const_5p + l41_rep_seq + const_3p
        print(l41_rep_seq_use)
        # Pre-allocate the mutation string
        mut_list = [" "]*len(l41_rep_seq_use)
        for i in range(len(l41_rep_seq_use)):
            if l41_rep_seq_use[i] != variants_all[0][i]:
                mut_list[i] = "x"
        mut_str = "".join(mut_list)
        print(mut_str)
        print(variants_all[0])
        print("\n")
        print("\n".join(variants_all))

        freq_MrnaDinuc_map_final = count_utr_sequence(variants_all)
        print(freq_MrnaDinuc_map)
        print(freq_MrnaDinuc_map_final)


        # Modify site list table
        site_list += [["No_site_%s" %i, "", "", ""] for i in range(1, 7)]
        print(site_list)
        print(len(site_list))

        output_df_cols = ["Seed", "ThreePrime", "Bulge", "Location", "Context",
                          "Sequence"]
        output_df = pd.DataFrame(index=range(952), columns=output_df_cols)
        # Iterate over the variant list to make the final datatable to output.
        for i in range(output_df.shape[0]):
            site_i = i % len(site_list)
            context_num = math.floor(i/len(site_list))
            if context_num == 0:
                context = "lin-41_context"
            else:
                context = "Context_%s" %context_num

            row_add = site_list[site_i] + [context] + [variants_all[i]]
            output_df.iloc[i, :] = row_add

        print(output_df.iloc[:5, :])
        # print(output_df)

        output_df.to_csv("ThreePrimePaperLet7ReporterAssay/libraries_210726/oligo_df_210726_library%s.txt" %lib_i, sep="\t")

        with open("ThreePrimePaperLet7ReporterAssay/sequence_checks_210726/sequences_check_210726_library%s.txt" %lib_i, "w") as file:
            file.write("\n".join(variants_all))


    print_time_elapsed(time_start)
    return()
    


################################################################################




if __name__ == "__main__":
    main()

