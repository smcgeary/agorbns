import imp # Used to import general.py
import time # Used for time
imp.load_source("general",
                ("/lab/bartel1_ata/mcgeary/computation/"
                 "AgoRBNS/general/general.py"
                )
               )
from general import get_rc, seq_mirna_map


def get_site_seq(mirna, start, stop, wobble=False,
                           mismatch=False, bulge=False):
    """Returns the sequence (in DNA form) based on the complementarity of
        the site.

    Args:
        mirna: The miRNA.
        start: The first nucleotide contained in the site (1 = first position).
        stop: The last nucleotide contained in the site.
        wobble: If not False, a number specifying the identify of the wobble
            paired site.

    Returns:
        A string giving the sequence corresponding to the requisite site.
    """
    # Convert miRNA sequence to list (which can be updated)
    rc_list = [i for i in seq_mirna_map[mirna]]
    # Make the sequence consistent with the A1 rule.
    rc_list[0] = "T"
    # Update for wobbles.
    if wobble:
        wobble_nucleotide_map = {"G" : "A", "U" : "C"}
        if type(wobble) == list:
            for wob in wobble:
                rc_list[wob - 1] = wobble_nucleotide_map[rc_list[wob - 1]]
        else:
            wobble_nucleotide_map = {"G" : "A", "U" : "C"}
            rc_list[wobble - 1] = wobble_nucleotide_map[rc_list[wobble - 1]]
    # Update for mismatches.
    if mismatch:
        for mm_ind in range(0, len(mismatch), 2):
            rc_list[mismatch[mm_ind]-1] = get_rc(mismatch[mm_ind + 1])
    # Define the end before bulge introduction.
    rc_list = rc_list[ : stop]
    # Define bulge.
    if bulge:
        rc_list.insert(bulge[0] - 1, get_rc(bulge[1]))
    # Define starting position.
    rc_list = rc_list[start - 1 : ]
    # Get and return the reverce complement of the merged string.
    site_seq = get_rc("".join(rc_list))
    return site_seq


def get_seq_site_map(mirna,sitelist):
    """Returns the sequence (in DNA form) based on the complementarity of
        the site.

    Args:
        mirna_seq: The miRNA sequence.
        name: The site type name

    Returns:
        A string giving the sequence corresponding to the requisite site.
    """

    # Load file with site types for the miRNA.
    sites_file_name = ("/lab/bartel1_ata/mcgeary/computation/"
                       "AgoRBNS/AssignSiteTypes/sites.%s_%s.txt" % (mirna, sitelist))
    # Get sites from file.
    with open(sites_file_name) as file_in:
        sites = file_in.read().split("\n")

    seq_site_map = {
        "4mer-A1" : get_site_seq(mirna, 1, 4),
        "4mer-m2.5" : get_site_seq(mirna, 2, 5),
        "4mer-m3.6" : get_site_seq(mirna, 3, 6),
        "4mer-m4.7" : get_site_seq(mirna, 4, 7),
        "4mer-m8" : get_site_seq(mirna, 5, 8),


        "5mer-A1" : get_site_seq(mirna, 1, 5),
        "5mer-m2.6" : get_site_seq(mirna, 2, 6),
        "5mer-m3.7" : get_site_seq(mirna, 3, 7),
        "5mer-m8" : get_site_seq(mirna, 4, 8),

        "6mer-A1" : get_site_seq(mirna, 1, 6),
        "6mer-A1bA5" : get_site_seq(mirna, 1, 6, bulge = [5, "A"]),
        "6mer-A1bA6" : get_site_seq(mirna, 1, 6, bulge = [6, "A"]),
        "6mer-A1bT6" : get_site_seq(mirna, 1, 6, bulge = [6, "T"]),
        "6mer-A1bG6" : get_site_seq(mirna, 1, 6, bulge = [6, "G"]),
        "6mer-A1mmG5" : get_site_seq(mirna, 1, 6, mismatch = [5, "G"]),
        "6mer-A1mmC5" : get_site_seq(mirna, 1, 6, mismatch = [5, "C"]),
        "6mer-A1mmT6" : get_site_seq(mirna, 1, 6, mismatch = [6, "T"]),
        "6mer-A1mmT5" : get_site_seq(mirna, 1, 6, mismatch = [5, "T"]),
        "6mer-A1mmG6" : get_site_seq(mirna, 1, 6, mismatch = [6, "G"]),
        "6mer-A1bT(4.6)" : get_site_seq(mirna, 1, 6, bulge = [4, "T"]),
        "6mer" : get_site_seq(mirna, 2, 7),
        "6mer-bA5" : get_site_seq(mirna, 2, 7, bulge = [5, "A"]),
        "6mer-mmA5" : get_site_seq(mirna, 2, 7, mismatch = [5, "A"]),
        "6mer-mmC5" : get_site_seq(mirna, 2, 7, mismatch = [5, "C"]),
        "6mer-mmG5" : get_site_seq(mirna, 2, 7, mismatch = [5, "G"]),
        "6mer-mmT5" : get_site_seq(mirna, 2, 7, mismatch = [5, "T"]),
        "6mer-mmA6" : get_site_seq(mirna, 2, 7, mismatch = [6, "A"]),
        "6mer-mmC6" : get_site_seq(mirna, 2, 7, mismatch = [6, "C"]),
        "6mer-mmG6" : get_site_seq(mirna, 2, 7, mismatch = [6, "G"]),
        "6mer-mmT6" : get_site_seq(mirna, 2, 7, mismatch = [6, "T"]),
        "6mer-bT6" : get_site_seq(mirna, 2, 7, bulge = [6, "T"]),
        "6mer-bT7" : get_site_seq(mirna, 2, 7, bulge = [7, "T"]),
        "6mer-bA(6.7)" : get_site_seq(mirna, 2, 7, bulge = [6, "A"]),
        "6mer-bG(6.7)" : get_site_seq(mirna, 2, 7, bulge = [6, "G"]),
        "6mer-bT(4.6)" : get_site_seq(mirna, 2, 7, bulge = [4, "T"]),
        "6mer-m8" : get_site_seq(mirna, 3, 8),
        "6mer-m8bA5" : get_site_seq(mirna, 3, 8, bulge = [5, "A"]),
        "6mer-m8bT6" : get_site_seq(mirna, 3, 8, bulge = [6, "T"]),
        "6mer-m8bA8" : get_site_seq(mirna, 3, 8, bulge = [8, "A"]),
        "6mer-m8bC8" : get_site_seq(mirna, 3, 8, bulge = [8, "C"]),
        "6mer-m8bT(4.6)" : get_site_seq(mirna, 3, 8, bulge = [4, "T"]),
        "6mer-m8bA(6.7)" : get_site_seq(mirna, 3, 8, bulge = [6, "A"]),
        "6mer-m8bG(6.7)" : get_site_seq(mirna, 3, 8, bulge = [6, "G"]),
        "6mer-m8bT(7.8)" : get_site_seq(mirna, 3, 8, bulge = [7, "T"]),
        "6mer-m8mmG5" : get_site_seq(mirna, 3, 8, mismatch = [5, "G"]),
        "6mer-m8mmC5" : get_site_seq(mirna, 3, 8, mismatch = [5, "C"]),
        "6mer-m8mmT6" : get_site_seq(mirna, 3, 8, mismatch = [6, "T"]),
        "6mer-m8mmT5" : get_site_seq(mirna, 3, 8, mismatch = [5, "T"]),
        "6mer-m8mmG6" : get_site_seq(mirna, 3, 8, mismatch = [6, "G"]),


        "6mer-m4.9" : get_site_seq(mirna, 4, 9),
        "6mer-m5.10" : get_site_seq(mirna, 5, 10),
        "6mer-m6.11" : get_site_seq(mirna, 6, 11),
        "6mer-m7.12" : get_site_seq(mirna, 7, 12),
        "6mer-m8.13" : get_site_seq(mirna, 8, 13),
        "6mer-m9.14" : get_site_seq(mirna, 9, 14),
        "6mer-m10.15" : get_site_seq(mirna, 10, 15),
        "6mer-m11.16" : get_site_seq(mirna, 11, 16),
        "6mer-m12.17" : get_site_seq(mirna, 12, 17),
        "6mer-m13.18" : get_site_seq(mirna, 13, 18),
        "6mer-m14.19" : get_site_seq(mirna, 14, 19),
        "6mer-m15.20" : get_site_seq(mirna, 15, 20),
        "6mer-m16.21" : get_site_seq(mirna, 16, 21),
        "6mer-m17.22" : get_site_seq(mirna, 17, 22),


        "7mer-A1" : get_site_seq(mirna, 1, 7),
        "7mer-A1mmA5" : get_site_seq(mirna, 1, 7, mismatch = [5, "A"]),
        "7mer-A1mmC5" : get_site_seq(mirna, 1, 7, mismatch = [5, "C"]),
        "7mer-A1mmG5" : get_site_seq(mirna, 1, 7, mismatch = [5, "G"]),
        "7mer-A1mmT5" : get_site_seq(mirna, 1, 7, mismatch = [5, "T"]),
        "7mer-A1mmG6" : get_site_seq(mirna, 1, 7, mismatch = [6, "G"]),
        "7mer-A1mmT6" : get_site_seq(mirna, 1, 7, mismatch = [6, "T"]),
        "7mer-A1bT2" : get_site_seq(mirna, 1, 7, bulge = [2, "T"]),
        "7mer-A1bT3" : get_site_seq(mirna, 1, 7, bulge = [3, "T"]),
        "7mer-A1bA5" : get_site_seq(mirna, 1, 7, bulge = [5, "A"]),
        "7mer-A1bT6" : get_site_seq(mirna, 1, 7, bulge = [6, "T"]),
        "7mer-A1bT7" : get_site_seq(mirna, 1, 7, bulge = [7, "T"]),
        "7mer-A1bT(4.6)" : get_site_seq(mirna, 1, 7, bulge = [5, "T"]),
        "7mer-A1bA(6.7)" : get_site_seq(mirna, 1, 7, bulge = [6, "A"]),
        "7mer-A1bG(6.7)" : get_site_seq(mirna, 1, 7, bulge = [6, "G"]),
        "7mer-m8" : get_site_seq(mirna, 2, 8),
        "7mer-m8bC(4.6)" : get_site_seq(mirna, 2, 8, bulge = [4, "C"]),
        "7mer-m8bT(4.6)" : get_site_seq(mirna, 2, 8, bulge = [4, "T"]),
        "7mer-m8bA5" : get_site_seq(mirna, 2, 8, bulge = [5, "A"]),
        "7mer-m8bA6" : get_site_seq(mirna, 2, 8, bulge = [6, "A"]),
        "7mer-m8bT6" : get_site_seq(mirna, 2, 8, bulge = [6, "T"]),
        "7mer-m8bT2" : get_site_seq(mirna, 2, 8, bulge = [2, "T"]),
        "7mer-m8bT3" : get_site_seq(mirna, 2, 8, bulge = [3, "T"]),
        "7mer-m8bA(6.7)" : get_site_seq(mirna, 2, 8, bulge = [6, "A"]),
        "7mer-m8bG(6.7)" : get_site_seq(mirna, 2, 8, bulge = [6, "G"]),
        "7mer-m8bG7" : get_site_seq(mirna, 2, 8, bulge = [7, "G"]),
        "7mer-m8bT(7.8)" : get_site_seq(mirna, 2, 8, bulge = [7, "T"]),
        "7mer-m8bA8" : get_site_seq(mirna, 2, 8, bulge = [8, "A"]),
        "7mer-m8bC8" : get_site_seq(mirna, 2, 8, bulge = [8, "C"]),
        "7mer-m8mmA5" : get_site_seq(mirna, 2, 8, mismatch = [5, "A"]),
        "7mer-m8mmC5" : get_site_seq(mirna, 2, 8, mismatch = [5, "C"]),
        "7mer-m8mmG5" : get_site_seq(mirna, 2, 8, mismatch = [5, "G"]),
        "7mer-m8mmT5" : get_site_seq(mirna, 2, 8, mismatch = [5, "T"]),
        "7mer-m8mmG6" : get_site_seq(mirna, 2, 8, mismatch = [6, "G"]),
        "7mer-m8mmC5" : get_site_seq(mirna, 2, 8, mismatch = [5, "C"]),
        "7mer-m8mmT5" : get_site_seq(mirna, 2, 8, mismatch = [5, "T"]),
        "7mer-m8mmT6" : get_site_seq(mirna, 2, 8, mismatch = [6, "T"]),
        "7mer-m8mmA7" : get_site_seq(mirna, 2, 8, mismatch = [7, "A"]),
        "7mer-m8mmC7" : get_site_seq(mirna, 2, 8, mismatch = [7, "C"]),
        "7mer-m8mmG7" : get_site_seq(mirna, 2, 8, mismatch = [7, "G"]),
        "7mer-m8mmT7" : get_site_seq(mirna, 2, 8, mismatch = [7, "T"]),



        "8mer" : get_site_seq(mirna, 1, 8),
        "8mer-bT(4.6)" : get_site_seq(mirna, 1, 8, bulge = [5, "T"]),
        "8mer-bT2" : get_site_seq(mirna, 1, 8, bulge = [2, "T"]),
        "8mer-bT3" : get_site_seq(mirna, 1, 8, bulge = [3, "T"]),
        "8mer-bA5" : get_site_seq(mirna, 1, 8, bulge = [5, "A"]),
        "8mer-bT6" : get_site_seq(mirna, 1, 8, bulge = [6, "T"]),
        "8mer-bA(6.7)" : get_site_seq(mirna, 1, 8, bulge = [6, "A"]),
        "8mer-bG(6.7)" : get_site_seq(mirna, 1, 8, bulge = [6, "G"]),
        "8mer-bG7" : get_site_seq(mirna, 1, 8, bulge = [7, "G"]),
        "8mer-bT(7.8)" : get_site_seq(mirna, 1, 8, bulge = [7, "T"]),
        "8mer-bA8" : get_site_seq(mirna, 1, 8, bulge = [8, "A"]),
        "8mer-bC8" : get_site_seq(mirna, 1, 8, bulge = [8, "C"]),

        "8mer-mmA2" : get_site_seq(mirna, 1, 8, mismatch = [2, "A"]),
        "8mer-mmC2" : get_site_seq(mirna, 1, 8, mismatch = [2, "C"]),
        "8mer-mmG2" : get_site_seq(mirna, 1, 8, mismatch = [2, "G"]),
        "8mer-mmT2" : get_site_seq(mirna, 1, 8, mismatch = [2, "T"]),


        "8mer-mmA3" : get_site_seq(mirna, 1, 8, mismatch = [3, "A"]),
        "8mer-mmC3" : get_site_seq(mirna, 1, 8, mismatch = [3, "C"]),
        "8mer-mmG3" : get_site_seq(mirna, 1, 8, mismatch = [3, "G"]),
        "8mer-mmT3" : get_site_seq(mirna, 1, 8, mismatch = [3, "T"]),

        "8mer-mmA4" : get_site_seq(mirna, 1, 8, mismatch = [4, "A"]),
        "8mer-mmC4" : get_site_seq(mirna, 1, 8, mismatch = [4, "C"]),
        "8mer-mmG4" : get_site_seq(mirna, 1, 8, mismatch = [4, "G"]),
        "8mer-mmT4" : get_site_seq(mirna, 1, 8, mismatch = [4, "T"]),

        "8mer-mmA5" : get_site_seq(mirna, 1, 8, mismatch = [5, "A"]),
        "8mer-mmC5" : get_site_seq(mirna, 1, 8, mismatch = [5, "C"]),
        "8mer-mmG5" : get_site_seq(mirna, 1, 8, mismatch = [5, "G"]),
        "8mer-mmT5" : get_site_seq(mirna, 1, 8, mismatch = [5, "T"]),

        "8mer-mmA6" : get_site_seq(mirna, 1, 8, mismatch = [6, "A"]),
        "8mer-mmC6" : get_site_seq(mirna, 1, 8, mismatch = [6, "C"]),
        "8mer-mmG6" : get_site_seq(mirna, 1, 8, mismatch = [6, "G"]),
        "8mer-mmT6" : get_site_seq(mirna, 1, 8, mismatch = [6, "T"]),

        "8mer-mmA7" : get_site_seq(mirna, 1, 8, mismatch = [7, "A"]),
        "8mer-mmC7" : get_site_seq(mirna, 1, 8, mismatch = [7, "C"]),
        "8mer-mmG7" : get_site_seq(mirna, 1, 8, mismatch = [7, "G"]),
        "8mer-mmT7" : get_site_seq(mirna, 1, 8, mismatch = [7, "T"]),


        "8mer-mmA5mmA6mmA7" : get_site_seq(mirna, 1, 8, mismatch = [5, "A", 6, "A", 7, "A"]),
        "8mer-mmA5mmA6mmC7" : get_site_seq(mirna, 1, 8, mismatch = [5, "A", 6, "A", 7, "C"]),
        "8mer-mmA5mmA6mmG7" : get_site_seq(mirna, 1, 8, mismatch = [5, "A", 6, "A", 7, "G"]),

        "8mer-mmA5mmC6mmA7" : get_site_seq(mirna, 1, 8, mismatch = [5, "A", 6, "C", 7, "A"]),
        "8mer-mmA5mmC6mmC7" : get_site_seq(mirna, 1, 8, mismatch = [5, "A", 6, "C", 7, "C"]),
        "8mer-mmA5mmC6mmG7" : get_site_seq(mirna, 1, 8, mismatch = [5, "A", 6, "C", 7, "G"]),
        "8mer-mmA5mmC6mmT7" : get_site_seq(mirna, 1, 8, mismatch = [5, "A", 6, "C", 7, "T"]),

        "8mer-mmA5mmG6mmA7" : get_site_seq(mirna, 1, 8, mismatch = [5, "A", 6, "G", 7, "A"]),
        "8mer-mmA5mmG6mmC7" : get_site_seq(mirna, 1, 8, mismatch = [5, "A", 6, "G", 7, "C"]),
        "8mer-mmA5mmG6mmG7" : get_site_seq(mirna, 1, 8, mismatch = [5, "A", 6, "G", 7, "G"]),
        "8mer-mmA5mmG6mmT7" : get_site_seq(mirna, 1, 8, mismatch = [5, "A", 6, "G", 7, "T"]),

        "8mer-mmA5mmT6mmA7" : get_site_seq(mirna, 1, 8, mismatch = [5, "A", 6, "T", 7, "A"]),
        "8mer-mmA5mmT6mmC7" : get_site_seq(mirna, 1, 8, mismatch = [5, "A", 6, "T", 7, "C"]),
        "8mer-mmA5mmT6mmG7" : get_site_seq(mirna, 1, 8, mismatch = [5, "A", 6, "T", 7, "G"]),
        "8mer-mmA5mmT6mmT7" : get_site_seq(mirna, 1, 8, mismatch = [5, "A", 6, "T", 7, "T"]),
        
        "8mer-mmC5mmA6mmA7" : get_site_seq(mirna, 1, 8, mismatch = [5, "C", 6, "A", 7, "A"]),
        "8mer-mmC5mmA6mmC7" : get_site_seq(mirna, 1, 8, mismatch = [5, "C", 6, "A", 7, "C"]),
        "8mer-mmC5mmA6mmT7" : get_site_seq(mirna, 1, 8, mismatch = [5, "C", 6, "A", 7, "T"]),

        "8mer-mmC5mmC6mmA7" : get_site_seq(mirna, 1, 8, mismatch = [5, "C", 6, "C", 7, "A"]),
        "8mer-mmC5mmC6mmG7" : get_site_seq(mirna, 1, 8, mismatch = [5, "C", 6, "C", 7, "G"]),
        "8mer-mmC5mmC6mmT7" : get_site_seq(mirna, 1, 8, mismatch = [5, "C", 6, "C", 7, "T"]),

        "8mer-mmC5mmG6mmA7" : get_site_seq(mirna, 1, 8, mismatch = [5, "C", 6, "G", 7, "A"]),
        "8mer-mmC5mmG6mmC7" : get_site_seq(mirna, 1, 8, mismatch = [5, "C", 6, "G", 7, "C"]),
        "8mer-mmC5mmG6mmG7" : get_site_seq(mirna, 1, 8, mismatch = [5, "C", 6, "G", 7, "G"]),
        "8mer-mmC5mmG6mmT7" : get_site_seq(mirna, 1, 8, mismatch = [5, "C", 6, "G", 7, "T"]),

        "8mer-mmC5mmT6mmA7" : get_site_seq(mirna, 1, 8, mismatch = [5, "C", 6, "T", 7, "A"]),
        "8mer-mmC5mmT6mmC7" : get_site_seq(mirna, 1, 8, mismatch = [5, "C", 6, "T", 7, "C"]),
        "8mer-mmC5mmT6mmG7" : get_site_seq(mirna, 1, 8, mismatch = [5, "C", 6, "T", 7, "G"]),
        "8mer-mmC5mmT6mmT7" : get_site_seq(mirna, 1, 8, mismatch = [5, "C", 6, "T", 7, "T"]),

   
        "8mer-mmG5mmA6mmA7" : get_site_seq(mirna, 1, 8, mismatch = [5, "G", 6, "A", 7, "A"]),
        "8mer-mmG5mmA6mmC7" : get_site_seq(mirna, 1, 8, mismatch = [5, "G", 6, "A", 7, "C"]),
        "8mer-mmG5mmA6mmG7" : get_site_seq(mirna, 1, 8, mismatch = [5, "G", 6, "A", 7, "G"]),
        "8mer-mmG5mmA6mmT7" : get_site_seq(mirna, 1, 8, mismatch = [5, "G", 6, "A", 7, "T"]),

        "8mer-mmG5mmC6mmA7" : get_site_seq(mirna, 1, 8, mismatch = [5, "G", 6, "C", 7, "A"]),
        "8mer-mmG5mmC6mmC7" : get_site_seq(mirna, 1, 8, mismatch = [5, "G", 6, "C", 7, "C"]),
        "8mer-mmG5mmC6mmG7" : get_site_seq(mirna, 1, 8, mismatch = [5, "G", 6, "C", 7, "G"]),
        "8mer-mmG5mmC6mmT7" : get_site_seq(mirna, 1, 8, mismatch = [5, "G", 6, "C", 7, "T"]),

        "8mer-mmG5mmG6mmA7" : get_site_seq(mirna, 1, 8, mismatch = [5, "G", 6, "G", 7, "A"]),
        "8mer-mmG5mmG6mmC7" : get_site_seq(mirna, 1, 8, mismatch = [5, "G", 6, "G", 7, "C"]),
        "8mer-mmG5mmG6mmG7" : get_site_seq(mirna, 1, 8, mismatch = [5, "G", 6, "G", 7, "G"]),
        "8mer-mmG5mmG6mmT7" : get_site_seq(mirna, 1, 8, mismatch = [5, "G", 6, "G", 7, "T"]),
 
        "8mer-mmG5mmT6mmA7" : get_site_seq(mirna, 1, 8, mismatch = [5, "G", 6, "T", 7, "A"]),
        "8mer-mmG5mmT6mmC7" : get_site_seq(mirna, 1, 8, mismatch = [5, "G", 6, "T", 7, "C"]),
        "8mer-mmG5mmT6mmG7" : get_site_seq(mirna, 1, 8, mismatch = [5, "G", 6, "T", 7, "G"]),
        "8mer-mmG5mmT6mmT7" : get_site_seq(mirna, 1, 8, mismatch = [5, "G", 6, "T", 7, "T"]),


        "8mer-mmT5mmA6mmA7" : get_site_seq(mirna, 1, 8, mismatch = [5, "T", 6, "A", 7, "A"]),
        "8mer-mmT5mmA6mmC7" : get_site_seq(mirna, 1, 8, mismatch = [5, "T", 6, "A", 7, "C"]),
        "8mer-mmT5mmA6mmG7" : get_site_seq(mirna, 1, 8, mismatch = [5, "T", 6, "A", 7, "G"]),
        "8mer-mmT5mmA6mmT7" : get_site_seq(mirna, 1, 8, mismatch = [5, "T", 6, "A", 7, "T"]),

        "8mer-mmT5mmC6mmA7" : get_site_seq(mirna, 1, 8, mismatch = [5, "T", 6, "C", 7, "A"]),
        "8mer-mmT5mmC6mmC7" : get_site_seq(mirna, 1, 8, mismatch = [5, "T", 6, "C", 7, "C"]),
        "8mer-mmT5mmC6mmG7" : get_site_seq(mirna, 1, 8, mismatch = [5, "T", 6, "C", 7, "G"]),

        "8mer-mmT5mmG6mmA7" : get_site_seq(mirna, 1, 8, mismatch = [5, "T", 6, "G", 7, "A"]),
        "8mer-mmT5mmG6mmC7" : get_site_seq(mirna, 1, 8, mismatch = [5, "T", 6, "G", 7, "C"]),
        "8mer-mmT5mmG6mmG7" : get_site_seq(mirna, 1, 8, mismatch = [5, "T", 6, "G", 7, "G"]),
        "8mer-mmT5mmG6mmT7" : get_site_seq(mirna, 1, 8, mismatch = [5, "T", 6, "G", 7, "T"]),
 
        "8mer-mmT5mmT6mmA7" : get_site_seq(mirna, 1, 8, mismatch = [5, "T", 6, "T", 7, "A"]),
        "8mer-mmT5mmT6mmC7" : get_site_seq(mirna, 1, 8, mismatch = [5, "T", 6, "T", 7, "C"]),
        "8mer-mmT5mmT6mmG7" : get_site_seq(mirna, 1, 8, mismatch = [5, "T", 6, "T", 7, "G"]),
        "8mer-mmT5mmT6mmT7" : get_site_seq(mirna, 1, 8, mismatch = [5, "T", 6, "T", 7, "T"]),









        "6mer-m3.8" : get_site_seq(mirna, 3, 8),
        "6mer-m4.9" : get_site_seq(mirna, 4, 9),
        "6mer-m5.10" : get_site_seq(mirna, 5, 10),
        "6mer-m6.11" : get_site_seq(mirna, 6, 11),
        "6mer-m7.12" : get_site_seq(mirna, 7, 12),
        "6mer-m8.13" : get_site_seq(mirna, 8, 13),
        "6mer-m9.14" : get_site_seq(mirna, 9, 14),
        "6mer-m10.15" : get_site_seq(mirna, 10, 15),
        "6mer-m11.16" : get_site_seq(mirna, 11, 16),
        "6mer-m12.17" : get_site_seq(mirna, 12, 17),
        "6mer-m13.18" : get_site_seq(mirna, 13, 18),
        "6mer-m14.19" : get_site_seq(mirna, 14, 19),
        "6mer-m15.20" : get_site_seq(mirna, 15, 20),
        "6mer-m16.21" : get_site_seq(mirna, 16, 21),
        "6mer-m17.22" : get_site_seq(mirna, 17, 22),



        "7mer-m3.9" : get_site_seq(mirna, 3, 9),

        "8mer-m2.9mmA2mmA4mmA6" : get_site_seq(mirna, 2, 9, mismatch = [2, "A", 4, "A", 6, "A"]),
        "8mer-m2.9mmA2mmA4mmC6" : get_site_seq(mirna, 2, 9, mismatch = [2, "A", 4, "A", 6, "C"]),
        "8mer-m2.9mmA2mmA4mmG6" : get_site_seq(mirna, 2, 9, mismatch = [2, "A", 4, "A", 6, "G"]),
        "8mer-m2.9mmA2mmA4mmT6" : get_site_seq(mirna, 2, 9, mismatch = [2, "A", 4, "A", 6, "T"]),
        "8mer-m2.9mmA2mmC4mmA6" : get_site_seq(mirna, 2, 9, mismatch = [2, "A", 4, "C", 6, "A"]),
        "8mer-m2.9mmA2mmC4mmC6" : get_site_seq(mirna, 2, 9, mismatch = [2, "A", 4, "C", 6, "C"]),
        "8mer-m2.9mmA2mmC4mmG6" : get_site_seq(mirna, 2, 9, mismatch = [2, "A", 4, "C", 6, "G"]),
        "8mer-m2.9mmA2mmC4mmT6" : get_site_seq(mirna, 2, 9, mismatch = [2, "A", 4, "C", 6, "T"]),
        "8mer-m2.9mmA2mmG4mmA6" : get_site_seq(mirna, 2, 9, mismatch = [2, "A", 4, "G", 6, "A"]),
        "8mer-m2.9mmA2mmG4mmC6" : get_site_seq(mirna, 2, 9, mismatch = [2, "A", 4, "G", 6, "C"]),
        "8mer-m2.9mmA2mmG4mmG6" : get_site_seq(mirna, 2, 9, mismatch = [2, "A", 4, "G", 6, "G"]),
        "8mer-m2.9mmA2mmG4mmT6" : get_site_seq(mirna, 2, 9, mismatch = [2, "A", 4, "G", 6, "T"]),
        "8mer-m2.9mmA2mmT4mmA6" : get_site_seq(mirna, 2, 9, mismatch = [2, "A", 4, "T", 6, "A"]),
        "8mer-m2.9mmA2mmT4mmC6" : get_site_seq(mirna, 2, 9, mismatch = [2, "A", 4, "T", 6, "C"]),
        "8mer-m2.9mmA2mmT4mmG6" : get_site_seq(mirna, 2, 9, mismatch = [2, "A", 4, "T", 6, "G"]),
        "8mer-m2.9mmA2mmT4mmT6" : get_site_seq(mirna, 2, 9, mismatch = [2, "A", 4, "T", 6, "T"]),

        "8mer-m2.9mmC2mmA4mmA6" : get_site_seq(mirna, 2, 9, mismatch = [2, "C", 4, "A", 6, "A"]),
        "8mer-m2.9mmC2mmA4mmC6" : get_site_seq(mirna, 2, 9, mismatch = [2, "C", 4, "A", 6, "C"]),
        "8mer-m2.9mmC2mmA4mmG6" : get_site_seq(mirna, 2, 9, mismatch = [2, "C", 4, "A", 6, "G"]),
        "8mer-m2.9mmC2mmA4mmT6" : get_site_seq(mirna, 2, 9, mismatch = [2, "C", 4, "A", 6, "T"]),
        "8mer-m2.9mmC2mmC4mmA6" : get_site_seq(mirna, 2, 9, mismatch = [2, "C", 4, "C", 6, "A"]),
        "8mer-m2.9mmC2mmC4mmC6" : get_site_seq(mirna, 2, 9, mismatch = [2, "C", 4, "C", 6, "C"]),
        "8mer-m2.9mmC2mmC4mmG6" : get_site_seq(mirna, 2, 9, mismatch = [2, "C", 4, "C", 6, "G"]),
        "8mer-m2.9mmC2mmC4mmT6" : get_site_seq(mirna, 2, 9, mismatch = [2, "C", 4, "C", 6, "T"]),
        "8mer-m2.9mmC2mmG4mmA6" : get_site_seq(mirna, 2, 9, mismatch = [2, "C", 4, "G", 6, "A"]),
        "8mer-m2.9mmC2mmG4mmC6" : get_site_seq(mirna, 2, 9, mismatch = [2, "C", 4, "G", 6, "C"]),
        "8mer-m2.9mmC2mmG4mmG6" : get_site_seq(mirna, 2, 9, mismatch = [2, "C", 4, "G", 6, "G"]),
        "8mer-m2.9mmC2mmG4mmT6" : get_site_seq(mirna, 2, 9, mismatch = [2, "C", 4, "G", 6, "T"]),
        "8mer-m2.9mmC2mmT4mmA6" : get_site_seq(mirna, 2, 9, mismatch = [2, "C", 4, "T", 6, "A"]),
        "8mer-m2.9mmC2mmT4mmC6" : get_site_seq(mirna, 2, 9, mismatch = [2, "C", 4, "T", 6, "C"]),
        "8mer-m2.9mmC2mmT4mmG6" : get_site_seq(mirna, 2, 9, mismatch = [2, "C", 4, "T", 6, "G"]),
        "8mer-m2.9mmC2mmT4mmT6" : get_site_seq(mirna, 2, 9, mismatch = [2, "C", 4, "T", 6, "T"]),

        "8mer-m2.9mmG2mmA4mmA6" : get_site_seq(mirna, 2, 9, mismatch = [2, "G", 4, "A", 6, "A"]),
        "8mer-m2.9mmG2mmA4mmC6" : get_site_seq(mirna, 2, 9, mismatch = [2, "G", 4, "A", 6, "C"]),
        "8mer-m2.9mmG2mmA4mmG6" : get_site_seq(mirna, 2, 9, mismatch = [2, "G", 4, "A", 6, "G"]),
        "8mer-m2.9mmG2mmA4mmT6" : get_site_seq(mirna, 2, 9, mismatch = [2, "G", 4, "A", 6, "T"]),
        "8mer-m2.9mmG2mmC4mmA6" : get_site_seq(mirna, 2, 9, mismatch = [2, "G", 4, "C", 6, "A"]),
        "8mer-m2.9mmG2mmC4mmC6" : get_site_seq(mirna, 2, 9, mismatch = [2, "G", 4, "C", 6, "C"]),
        "8mer-m2.9mmG2mmC4mmG6" : get_site_seq(mirna, 2, 9, mismatch = [2, "G", 4, "C", 6, "G"]),
        "8mer-m2.9mmG2mmC4mmT6" : get_site_seq(mirna, 2, 9, mismatch = [2, "G", 4, "C", 6, "T"]),
        "8mer-m2.9mmG2mmG4mmA6" : get_site_seq(mirna, 2, 9, mismatch = [2, "G", 4, "G", 6, "A"]),
        "8mer-m2.9mmG2mmG4mmC6" : get_site_seq(mirna, 2, 9, mismatch = [2, "G", 4, "G", 6, "C"]),
        "8mer-m2.9mmG2mmG4mmG6" : get_site_seq(mirna, 2, 9, mismatch = [2, "G", 4, "G", 6, "G"]),
        "8mer-m2.9mmG2mmG4mmT6" : get_site_seq(mirna, 2, 9, mismatch = [2, "G", 4, "G", 6, "T"]),
        "8mer-m2.9mmG2mmT4mmA6" : get_site_seq(mirna, 2, 9, mismatch = [2, "G", 4, "T", 6, "A"]),
        "8mer-m2.9mmG2mmT4mmC6" : get_site_seq(mirna, 2, 9, mismatch = [2, "G", 4, "T", 6, "C"]),
        "8mer-m2.9mmG2mmT4mmG6" : get_site_seq(mirna, 2, 9, mismatch = [2, "G", 4, "T", 6, "G"]),
        "8mer-m2.9mmG2mmT4mmT6" : get_site_seq(mirna, 2, 9, mismatch = [2, "G", 4, "T", 6, "T"]),

        "8mer-m2.9mmT2mmA4mmA6" : get_site_seq(mirna, 2, 9, mismatch = [2, "T", 4, "A", 6, "A"]),
        "8mer-m2.9mmT2mmA4mmC6" : get_site_seq(mirna, 2, 9, mismatch = [2, "T", 4, "A", 6, "C"]),
        "8mer-m2.9mmT2mmA4mmG6" : get_site_seq(mirna, 2, 9, mismatch = [2, "T", 4, "A", 6, "G"]),
        "8mer-m2.9mmT2mmA4mmT6" : get_site_seq(mirna, 2, 9, mismatch = [2, "T", 4, "A", 6, "T"]),
        "8mer-m2.9mmT2mmC4mmA6" : get_site_seq(mirna, 2, 9, mismatch = [2, "T", 4, "C", 6, "A"]),
        "8mer-m2.9mmT2mmC4mmC6" : get_site_seq(mirna, 2, 9, mismatch = [2, "T", 4, "C", 6, "C"]),
        "8mer-m2.9mmT2mmC4mmG6" : get_site_seq(mirna, 2, 9, mismatch = [2, "T", 4, "C", 6, "G"]),
        "8mer-m2.9mmT2mmC4mmT6" : get_site_seq(mirna, 2, 9, mismatch = [2, "T", 4, "C", 6, "T"]),
        "8mer-m2.9mmT2mmG4mmA6" : get_site_seq(mirna, 2, 9, mismatch = [2, "T", 4, "G", 6, "A"]),
        "8mer-m2.9mmT2mmG4mmC6" : get_site_seq(mirna, 2, 9, mismatch = [2, "T", 4, "G", 6, "C"]),
        "8mer-m2.9mmT2mmG4mmG6" : get_site_seq(mirna, 2, 9, mismatch = [2, "T", 4, "G", 6, "G"]),
        "8mer-m2.9mmT2mmG4mmT6" : get_site_seq(mirna, 2, 9, mismatch = [2, "T", 4, "G", 6, "T"]),
        "8mer-m2.9mmT2mmT4mmA6" : get_site_seq(mirna, 2, 9, mismatch = [2, "T", 4, "T", 6, "A"]),
        "8mer-m2.9mmT2mmT4mmC6" : get_site_seq(mirna, 2, 9, mismatch = [2, "T", 4, "T", 6, "C"]),
        "8mer-m2.9mmT2mmT4mmG6" : get_site_seq(mirna, 2, 9, mismatch = [2, "T", 4, "T", 6, "G"]),
        "8mer-m2.9mmT2mmT4mmT6" : get_site_seq(mirna, 2, 9, mismatch = [2, "T", 4, "T", 6, "T"]),


        "7mer-m4.10" : get_site_seq(mirna, 4, 10),
        "7mer-m5.11" : get_site_seq(mirna, 5, 11),
        "7mer-m6.12" : get_site_seq(mirna, 6, 12),
        "7mer-m7.13" : get_site_seq(mirna, 7, 13),
        "7mer-m8.14" : get_site_seq(mirna, 8, 14),
        "7mer-m9.15" : get_site_seq(mirna, 9, 15),
        "7mer-m10.16" : get_site_seq(mirna, 10, 16),
        "7mer-m11.17" : get_site_seq(mirna, 11, 17),
        "7mer-m12.18" : get_site_seq(mirna, 12, 18),
        "7mer-m13.19" : get_site_seq(mirna, 13, 19),
        "7mer-m14.20" : get_site_seq(mirna, 14, 20),
        "7mer-m15.21" : get_site_seq(mirna, 15, 21),
        "7mer-m16.22" : get_site_seq(mirna, 16, 22),

        "8mer-m2.9" : get_site_seq(mirna, 2, 9),
        "8mer-m3.10" : get_site_seq(mirna, 3, 10),
        "8mer-m4.11" : get_site_seq(mirna, 4, 11),
        "8mer-m5.12" : get_site_seq(mirna, 5, 12),
        "8mer-m6.13" : get_site_seq(mirna, 6, 13),
        "8mer-m7.14" : get_site_seq(mirna, 7, 14),
        "8mer-m8.15" : get_site_seq(mirna, 8, 15),
        "8mer-m9.16" : get_site_seq(mirna, 9, 16),
        "8mer-m10.17" : get_site_seq(mirna, 10, 17),
        "8mer-m11.18" : get_site_seq(mirna, 11, 18),
        "8mer-m12.19" : get_site_seq(mirna, 12, 19),
        "8mer-m13.20" : get_site_seq(mirna, 13, 20),
        "8mer-m14.21" : get_site_seq(mirna, 14, 21),
        "8mer-m15.22" : get_site_seq(mirna, 15, 22),
        "8mer-m16.23" : get_site_seq(mirna, 16, 23),
        "8mer-m17.24" : get_site_seq(mirna, 17, 24),
        "8mer-m18.25" : get_site_seq(mirna, 18, 25),

        "9mer-m1.9" : get_site_seq(mirna, 1, 9),
        "9mer-m2.10" : get_site_seq(mirna, 2, 10),
        "9mer-m3.11" : get_site_seq(mirna, 3, 11),
        "9mer-m4.12" : get_site_seq(mirna, 4, 12),
        "9mer-m5.13" : get_site_seq(mirna, 5, 13),
        "9mer-m6.14" : get_site_seq(mirna, 6, 14),
        "9mer-m7.15" : get_site_seq(mirna, 7, 15),
        "9mer-m8.16" : get_site_seq(mirna, 8, 16),
        "9mer-m9.17" : get_site_seq(mirna, 9, 17),
        "9mer-m10.18" : get_site_seq(mirna, 10, 18),
        "9mer-m11.19" : get_site_seq(mirna, 11, 19),
        "9mer-m12.20" : get_site_seq(mirna, 12, 20),
        "9mer-m13.21" : get_site_seq(mirna, 13, 21),
        "9mer-m14.22" : get_site_seq(mirna, 14, 22),
        "9mer-m15.23" : get_site_seq(mirna, 15, 23),
        "9mer-m16.24" : get_site_seq(mirna, 16, 24),
        "9mer-m17.25" : get_site_seq(mirna, 17, 25),

        "10mer-m1.10" : get_site_seq(mirna, 1, 10),
        "10mer-m2.11" : get_site_seq(mirna, 2, 11),
        "10mer-m3.12" : get_site_seq(mirna, 3, 12),
        "10mer-m4.13" : get_site_seq(mirna, 4, 13),
        "10mer-m5.14" : get_site_seq(mirna, 5, 14),
        "10mer-m6.15" : get_site_seq(mirna, 6, 15),
        "10mer-m7.16" : get_site_seq(mirna, 7, 16),
        "10mer-m8.17" : get_site_seq(mirna, 8, 17),
        "10mer-m9.18" : get_site_seq(mirna, 9, 18),
        "10mer-m10.19" : get_site_seq(mirna, 10, 19),
        "10mer-m11.20" : get_site_seq(mirna, 11, 20),
        "10mer-m12.21" : get_site_seq(mirna, 12, 21),
        "10mer-m13.22" : get_site_seq(mirna, 13, 22),
        "10mer-m14.23" : get_site_seq(mirna, 14, 23),
        "10mer-m15.24" : get_site_seq(mirna, 15, 24),
        "10mer-m16.25" : get_site_seq(mirna, 16, 25),

        "11mer-m1.11" : get_site_seq(mirna, 1, 11),
        "11mer-m2.12" : get_site_seq(mirna, 2, 12),
        "11mer-m3.13" : get_site_seq(mirna, 3, 13),
        "11mer-m4.14" : get_site_seq(mirna, 4, 14),
        "11mer-m5.15" : get_site_seq(mirna, 5, 15),
        "11mer-m6.16" : get_site_seq(mirna, 6, 16),
        "11mer-m7.17" : get_site_seq(mirna, 7, 17),
        "11mer-m8.18" : get_site_seq(mirna, 8, 18),
        "11mer-m9.19" : get_site_seq(mirna, 9, 19),
        "11mer-m10.20" : get_site_seq(mirna, 10, 20),
        "11mer-m11.21" : get_site_seq(mirna, 11, 21),
        "11mer-m12.22" : get_site_seq(mirna, 12, 22),
        "11mer-m13.23" : get_site_seq(mirna, 13, 23),
        "11mer-m14.24" : get_site_seq(mirna, 14, 24),
        "11mer-m15.25" : get_site_seq(mirna, 15, 25),

        "12mer-m1.12" : get_site_seq(mirna, 1, 12),
        "12mer-m2.13" : get_site_seq(mirna, 2, 13),
        "12mer-m3.14" : get_site_seq(mirna, 3, 14),
        "12mer-m4.15" : get_site_seq(mirna, 4, 15),
        "12mer-m5.16" : get_site_seq(mirna, 5, 16),
        "12mer-m6.17" : get_site_seq(mirna, 6, 17),
        "12mer-m7.18" : get_site_seq(mirna, 7, 18),
        "12mer-m8.19" : get_site_seq(mirna, 8, 19),
        "12mer-m9.20" : get_site_seq(mirna, 9, 20),
        "12mer-m10.21" : get_site_seq(mirna, 10, 21),
        "12mer-m11.22" : get_site_seq(mirna, 11, 22),
        "12mer-m12.23" : get_site_seq(mirna, 12, 23),
        "12mer-m13.24" : get_site_seq(mirna, 13, 24),
        "12mer-m14.25" : get_site_seq(mirna, 14, 25),

        "13mer-m1.13" : get_site_seq(mirna, 1, 13),
        "13mer-m2.14" : get_site_seq(mirna, 2, 14),
        "13mer-m3.15" : get_site_seq(mirna, 3, 15),
        "13mer-m4.16" : get_site_seq(mirna, 4, 16),
        "13mer-m5.17" : get_site_seq(mirna, 5, 17),
        "13mer-m6.18" : get_site_seq(mirna, 6, 18),
        "13mer-m7.19" : get_site_seq(mirna, 7, 19),
        "13mer-m8.20" : get_site_seq(mirna, 8, 20),
        "13mer-m9.21" : get_site_seq(mirna, 9, 21),
        "13mer-m10.22" : get_site_seq(mirna, 10, 22),
        "13mer-m11.23" : get_site_seq(mirna, 11, 23),
        "13mer-m12.24" : get_site_seq(mirna, 12, 24),
        "13mer-m13.25" : get_site_seq(mirna, 13, 25),
        }

    # Check miRNA length and sequence for non-uniform site types.
    mirna_seq = seq_mirna_map[mirna]
    if len(mirna_seq) == 23:
        seq_site_map.update(
            {"6mer-m18.23" : get_site_seq(mirna, 18, 23),
             "7mer-m17.23" : get_site_seq(mirna, 17, 23),
             "8mer-m16.23" : get_site_seq(mirna, 16, 23),
             "9mer-m15.23" : get_site_seq(mirna, 15, 23),
             "10mer-m14.23" : get_site_seq(mirna, 14, 23),
             "11mer-m13.23" : get_site_seq(mirna, 13, 23),
             "12mer-m12.23" : get_site_seq(mirna, 12, 23),
             "13mer-m11.23" : get_site_seq(mirna, 11, 23)})
    if mirna_seq[2] in ["G", "U"]:
        seq_site_map.update(
            {"6mer-A1w3" : get_site_seq(mirna, 1, 6, wobble = 3)})
    for i in [1, 2, 3, 4, 5, 6, 7]:
        if mirna_seq[i] in ["G","U"]:
            seq_site_map.update(
                {"8mer-w%s" % (i + 1) : get_site_seq(mirna, 1, 8,
                                                     wobble = i + 1),
                 "7mer-m8w%s" % (i + 1) : get_site_seq(mirna, 2, 8,
                                                       wobble = i + 1)})
# let-7 wobbles
    if mirna_seq == "UGAGGUAGUAGGUUGUAUAGUU":
        seq_site_map.update(
            {"8mer-w2w4" : get_site_seq(mirna, 1, 8, wobble = [2, 4]),
             "8mer-w2w5" : get_site_seq(mirna, 1, 8, wobble = [2, 5]),
             "8mer-w2w6" : get_site_seq(mirna, 1, 8, wobble = [2, 6]),
             "8mer-w2w8" : get_site_seq(mirna, 1, 8, wobble = [2, 8]),
             "8mer-w4w5" : get_site_seq(mirna, 1, 8, wobble = [4, 5]),
             "8mer-w4w6" : get_site_seq(mirna, 1, 8, wobble = [4, 6]),
             "8mer-w4w8" : get_site_seq(mirna, 1, 8, wobble = [4, 8]),
             "8mer-w5w6" : get_site_seq(mirna, 1, 8, wobble = [5, 6]),
             "8mer-w5w8" : get_site_seq(mirna, 1, 8, wobble = [5, 8]),
             "8mer-w6w8" : get_site_seq(mirna, 1, 8, wobble = [6, 8]),
             "8mer-w2w4w5" : get_site_seq(mirna, 1, 8, wobble = [2, 4, 5]),
             "8mer-w2w4w6" : get_site_seq(mirna, 1, 8, wobble = [2, 4, 6]),
             "8mer-w2w4w8" : get_site_seq(mirna, 1, 8, wobble = [2, 4, 8]),
             "8mer-w2w5w6" : get_site_seq(mirna, 1, 8, wobble = [2, 5, 6]),
             "8mer-w2w5w8" : get_site_seq(mirna, 1, 8, wobble = [2, 5, 8]),
             "8mer-w2w6w8" : get_site_seq(mirna, 1, 8, wobble = [2, 6, 8]),
             "8mer-w4w5w6" : get_site_seq(mirna, 1, 8, wobble = [4, 5, 6]),
             "8mer-w4w5w8" : get_site_seq(mirna, 1, 8, wobble = [4, 5, 8]),
             "8mer-w4w6w8" : get_site_seq(mirna, 1, 8, wobble = [4, 6, 8]),
             "8mer-w5w6w8" : get_site_seq(mirna, 1, 8, wobble = [5, 6, 8]),
             "8mer-w2w4w5w6" : get_site_seq(mirna, 1, 8, wobble = [2, 4, 5, 6]),
             "8mer-w2w4w5w8" : get_site_seq(mirna, 1, 8, wobble = [2, 4, 5, 8]),
             "8mer-w2w4w6w8" : get_site_seq(mirna, 1, 8, wobble = [2, 4, 6, 8]),
             "8mer-w2w5w6w8" : get_site_seq(mirna, 1, 8, wobble = [2, 5, 6, 8]),
             "8mer-w4w5w6w8" : get_site_seq(mirna, 1, 8, wobble = [4, 5, 6, 8]),
             "8mer-w2w4w5w6w8" : get_site_seq(mirna, 1, 8, wobble = [2, 4, 5, 6, 8]),

             "7mer-m8w2w4" : get_site_seq(mirna, 2, 8, wobble = [2, 4]),
             "7mer-m8w2w5" : get_site_seq(mirna, 2, 8, wobble = [2, 5]),
             "7mer-m8w2w6" : get_site_seq(mirna, 2, 8, wobble = [2, 6]),
             "7mer-m8w2w8" : get_site_seq(mirna, 2, 8, wobble = [2, 8]),
             "7mer-m8w4w5" : get_site_seq(mirna, 2, 8, wobble = [4, 5]),
             "7mer-m8w4w6" : get_site_seq(mirna, 2, 8, wobble = [4, 6]),
             "7mer-m8w4w8" : get_site_seq(mirna, 2, 8, wobble = [4, 8]),
             "7mer-m8w5w6" : get_site_seq(mirna, 2, 8, wobble = [5, 6]),
             "7mer-m8w5w8" : get_site_seq(mirna, 2, 8, wobble = [5, 8]),
             "7mer-m8w6w8" : get_site_seq(mirna, 2, 8, wobble = [6, 8]),
             "7mer-m8w2w4w5" : get_site_seq(mirna, 2, 8, wobble = [2, 4, 5]),
             "7mer-m8w2w4w6" : get_site_seq(mirna, 2, 8, wobble = [2, 4, 6]),
             "7mer-m8w2w4w8" : get_site_seq(mirna, 2, 8, wobble = [2, 4, 8]),
             "7mer-m8w2w5w6" : get_site_seq(mirna, 2, 8, wobble = [2, 5, 6]),
             "7mer-m8w2w5w8" : get_site_seq(mirna, 2, 8, wobble = [2, 5, 8]),
             "7mer-m8w2w6w8" : get_site_seq(mirna, 2, 8, wobble = [2, 6, 8]),
             "7mer-m8w4w5w6" : get_site_seq(mirna, 2, 8, wobble = [4, 5, 6]),
             "7mer-m8w4w5w8" : get_site_seq(mirna, 2, 8, wobble = [4, 5, 8]),
             "7mer-m8w4w6w8" : get_site_seq(mirna, 2, 8, wobble = [4, 6, 8]),
             "7mer-m8w5w6w8" : get_site_seq(mirna, 2, 8, wobble = [5, 6, 8]),
             "7mer-m8w2w4w5w6" : get_site_seq(mirna, 2, 8, wobble = [2, 4, 5, 6]),
             "7mer-m8w2w4w5w8" : get_site_seq(mirna, 2, 8, wobble = [2, 4, 5, 8]),
             "7mer-m8w2w4w6w8" : get_site_seq(mirna, 2, 8, wobble = [2, 4, 6, 8]),
             "7mer-m8w2w5w6w8" : get_site_seq(mirna, 2, 8, wobble = [2, 5, 6, 8]),
             "7mer-m8w4w5w6w8" : get_site_seq(mirna, 2, 8, wobble = [4, 5, 6, 8]),
             "7mer-m8w2w4w5w6w8" : get_site_seq(mirna, 2, 8, wobble = [2, 4, 5, 6, 8]),

             "7mer-A1w2" : get_site_seq(mirna, 1, 7, wobble = 2),
             "7mer-A1w4" : get_site_seq(mirna, 1, 7, wobble = 4),
             "7mer-A1w5" : get_site_seq(mirna, 1, 7, wobble = 5),
             "7mer-A1w6" : get_site_seq(mirna, 1, 7, wobble = 6),



             "7mer-A1w2w4" : get_site_seq(mirna, 1, 7, wobble = [2, 4]),
             "7mer-A1w2w5" : get_site_seq(mirna, 1, 7, wobble = [2, 5]),
             "7mer-A1w2w6" : get_site_seq(mirna, 1, 7, wobble = [2, 6]),
             "7mer-A1w4w5" : get_site_seq(mirna, 1, 7, wobble = [4, 5]),
             "7mer-A1w4w6" : get_site_seq(mirna, 1, 7, wobble = [4, 6]),
             "7mer-A1w5w6" : get_site_seq(mirna, 1, 7, wobble = [5, 6]),
             "7mer-A1w2w4w5" : get_site_seq(mirna, 1, 7, wobble = [2, 4, 5]),
             "7mer-A1w2w4w6" : get_site_seq(mirna, 1, 7, wobble = [2, 4, 6]),
             "7mer-A1w2w5w6" : get_site_seq(mirna, 1, 7, wobble = [2, 5, 6]),
             "7mer-A1w4w5w6" : get_site_seq(mirna, 1, 7, wobble = [4, 5, 6]),
             "7mer-A1w2w4w5w6" : get_site_seq(mirna, 1, 7, wobble = [2, 4, 5, 6]),


             "6mer-A1w2" : get_site_seq(mirna, 1, 6, wobble = 2),
             "6mer-A1w4" : get_site_seq(mirna, 1, 6, wobble = 4),
             "6mer-A1w5" : get_site_seq(mirna, 1, 6, wobble = 5),
             "6mer-A1w6" : get_site_seq(mirna, 1, 6, wobble = 6),


             "6mer-A1w2w4" : get_site_seq(mirna, 1, 6, wobble = [2, 4]),
             "6mer-A1w2w5" : get_site_seq(mirna, 1, 6, wobble = [2, 5]),
             "6mer-A1w2w6" : get_site_seq(mirna, 1, 6, wobble = [2, 6]),
             "6mer-A1w4w5" : get_site_seq(mirna, 1, 6, wobble = [4, 5]),
             "6mer-A1w4w6" : get_site_seq(mirna, 1, 6, wobble = [4, 6]),
             "6mer-A1w5w6" : get_site_seq(mirna, 1, 6, wobble = [5, 6]),
             "6mer-A1w2w4w5" : get_site_seq(mirna, 1, 6, wobble = [2, 4, 5]),
             "6mer-A1w2w4w6" : get_site_seq(mirna, 1, 6, wobble = [2, 4, 6]),
             "6mer-A1w2w5w6" : get_site_seq(mirna, 1, 6, wobble = [2, 5, 6]),
             "6mer-A1w4w5w6" : get_site_seq(mirna, 1, 6, wobble = [4, 5, 6]),
             "6mer-A1w2w4w5w6" : get_site_seq(mirna, 1, 6, wobble = [2, 4, 5, 6]),


             "6mer-w2" : get_site_seq(mirna, 2, 7, wobble = 2),
             "6mer-w4" : get_site_seq(mirna, 2, 7, wobble = 4),
             "6mer-w5" : get_site_seq(mirna, 2, 7, wobble = 5),
             "6mer-w6" : get_site_seq(mirna, 2, 7, wobble = 6),


             "6mer-w2w4" : get_site_seq(mirna, 2, 7, wobble = [2, 4]),
             "6mer-w2w5" : get_site_seq(mirna, 2, 7, wobble = [2, 5]),
             "6mer-w2w6" : get_site_seq(mirna, 2, 7, wobble = [2, 6]),
             "6mer-w4w5" : get_site_seq(mirna, 2, 7, wobble = [4, 5]),
             "6mer-w4w6" : get_site_seq(mirna, 2, 7, wobble = [4, 6]),
             "6mer-w5w6" : get_site_seq(mirna, 2, 7, wobble = [5, 6]),
             "6mer-w2w4w5" : get_site_seq(mirna, 2, 7, wobble = [2, 4, 5]),
             "6mer-w2w4w6" : get_site_seq(mirna, 2, 7, wobble = [2, 4, 6]),
             "6mer-w2w5w6" : get_site_seq(mirna, 2, 7, wobble = [2, 5, 6]),
             "6mer-w4w5w6" : get_site_seq(mirna, 2, 7, wobble = [4, 5, 6]),
             "6mer-w2w4w5w6" : get_site_seq(mirna, 2, 7, wobble = [2, 4, 5, 6]),

             "6mer-m8w2" : get_site_seq(mirna, 3, 8, wobble = 2),
             "6mer-m8w4" : get_site_seq(mirna, 3, 8, wobble = 4),
             "6mer-m8w5" : get_site_seq(mirna, 3, 8, wobble = 5),
             "6mer-m8w6" : get_site_seq(mirna, 3, 8, wobble = 6),
             "6mer-m8w8" : get_site_seq(mirna, 3, 8, wobble = 8),


             "6mer-m8w4w5" : get_site_seq(mirna, 3, 8, wobble = [4, 5]),
             "6mer-m8w4w6" : get_site_seq(mirna, 3, 8, wobble = [4, 6]),
             "6mer-m8w4w8" : get_site_seq(mirna, 3, 8, wobble = [4, 8]),
             "6mer-m8w5w6" : get_site_seq(mirna, 3, 8, wobble = [5, 6]),
             "6mer-m8w5w8" : get_site_seq(mirna, 3, 8, wobble = [5, 8]),
             "6mer-m8w6w8" : get_site_seq(mirna, 3, 8, wobble = [6, 8]),
             "6mer-m8w4w5w6" : get_site_seq(mirna, 3, 8, wobble = [4, 5, 6]),
             "6mer-m8w4w5w8" : get_site_seq(mirna, 3, 8, wobble = [4, 5, 8]),
             "6mer-m8w4w6w8" : get_site_seq(mirna, 3, 8, wobble = [4, 6, 8]),
             "6mer-m8w5w6w8" : get_site_seq(mirna, 3, 8, wobble = [5, 6, 8]),
             "6mer-m8w4w5w6w8" : get_site_seq(mirna, 3, 8, wobble = [4, 5, 6, 8])})

# miR-155 wobbles
    if mirna_seq == "UUAAUGCUAAUCGUGAUAGGGGU":
        seq_site_map.update(
            {"8mer-w2w5" : get_site_seq(mirna, 1, 8, wobble = [2, 5]),
             "8mer-w2w6" : get_site_seq(mirna, 1, 8, wobble = [2, 6]),
             "8mer-w2w8" : get_site_seq(mirna, 1, 8, wobble = [2, 8]),
             "8mer-w5w6" : get_site_seq(mirna, 1, 8, wobble = [5, 6]),
             "8mer-w5w8" : get_site_seq(mirna, 1, 8, wobble = [5, 8]),
             "8mer-w6w8" : get_site_seq(mirna, 1, 8, wobble = [6, 8]),
             "8mer-w2w5w6" : get_site_seq(mirna, 1, 8, wobble = [2, 5, 6]),
             "8mer-w2w5w8" : get_site_seq(mirna, 1, 8, wobble = [2, 5, 8]),
             "8mer-w2w6w8" : get_site_seq(mirna, 1, 8, wobble = [2, 6, 8]),
             "8mer-w5w6w8" : get_site_seq(mirna, 1, 8, wobble = [5, 6, 8]),
             "8mer-w2w5w6w8" : get_site_seq(mirna, 1, 8, wobble = [2, 5, 6, 8]),

             "7mer-m8w2w5" : get_site_seq(mirna, 2, 8, wobble = [2, 5]),
             "7mer-m8w2w6" : get_site_seq(mirna, 2, 8, wobble = [2, 6]),
             "7mer-m8w2w8" : get_site_seq(mirna, 2, 8, wobble = [2, 8]),
             "7mer-m8w5w6" : get_site_seq(mirna, 2, 8, wobble = [5, 6]),
             "7mer-m8w5w8" : get_site_seq(mirna, 2, 8, wobble = [5, 8]),
             "7mer-m8w6w8" : get_site_seq(mirna, 2, 8, wobble = [6, 8]),
             "7mer-m8w2w5w6" : get_site_seq(mirna, 2, 8, wobble = [2, 5, 6]),
             "7mer-m8w2w5w8" : get_site_seq(mirna, 2, 8, wobble = [2, 5, 8]),
             "7mer-m8w2w6w8" : get_site_seq(mirna, 2, 8, wobble = [2, 6, 8]),
             "7mer-m8w5w6w8" : get_site_seq(mirna, 2, 8, wobble = [5, 6, 8]),
             "7mer-m8w2w5w6w8" : get_site_seq(mirna, 2, 8, wobble = [2, 5, 6, 8]),

             "7mer-A1w2" : get_site_seq(mirna, 1, 7, wobble = 2),
             "7mer-A1w5" : get_site_seq(mirna, 1, 7, wobble = 5),
             "7mer-A1w6" : get_site_seq(mirna, 1, 7, wobble = 6),

             "6mer-A1w2" : get_site_seq(mirna, 1, 6, wobble = 2),
             "6mer-A1w5" : get_site_seq(mirna, 1, 6, wobble = 5),
             "6mer-A1w6" : get_site_seq(mirna, 1, 6, wobble = 6),

             "6mer-w2" : get_site_seq(mirna, 2, 7, wobble = 2),
             "6mer-w5" : get_site_seq(mirna, 2, 7, wobble = 5),
             "6mer-w6" : get_site_seq(mirna, 2, 7, wobble = 6),

             "6mer-m8w5" : get_site_seq(mirna, 3, 8, wobble = 5),
             "6mer-m8w6" : get_site_seq(mirna, 3, 8, wobble = 6),
             "6mer-m8w8" : get_site_seq(mirna, 3, 8, wobble = 8),

             "7mer-A1w2w5" : get_site_seq(mirna, 1, 7, wobble = [2, 5]),
             "7mer-A1w2w6" : get_site_seq(mirna, 1, 7, wobble = [2, 6]),
             "7mer-A1w5w6" : get_site_seq(mirna, 1, 7, wobble = [5, 6]),
             "7mer-A1w2w5w6" : get_site_seq(mirna, 1, 7, wobble = [2, 5, 6]),

             "6mer-A1w2w5" : get_site_seq(mirna, 1, 6, wobble = [2, 5]),
             "6mer-A1w2w6" : get_site_seq(mirna, 1, 6, wobble = [2, 6]),
             "6mer-A1w5w6" : get_site_seq(mirna, 1, 6, wobble = [5, 6]),
             "6mer-A1w2w5w6" : get_site_seq(mirna, 1, 6, wobble = [2, 5, 6]),

             "6mer-w2w5" : get_site_seq(mirna, 2, 7, wobble = [2, 5]),
             "6mer-w2w6" : get_site_seq(mirna, 2, 7, wobble = [2, 6]),
             "6mer-w5w6" : get_site_seq(mirna, 2, 7, wobble = [5, 6]),
             "6mer-w2w5w6" : get_site_seq(mirna, 2, 7, wobble = [2, 5, 6]),

             "6mer-m8w5w6" : get_site_seq(mirna, 3, 8, wobble = [5, 6]),
             "6mer-m8w5w8" : get_site_seq(mirna, 3, 8, wobble = [5, 8]),
             "6mer-m8w6w8" : get_site_seq(mirna, 3, 8, wobble = [6, 8]),
             "6mer-m8w5w6w8" : get_site_seq(mirna, 3, 8, wobble = [5, 6, 8])})


#miR-124
    if mirna_seq == "UAAGGCACGCGGUGAAUGCCAA":
        seq_site_map.update(
            {"8mer-w4w5" : get_site_seq(mirna, 1, 8, wobble = [4, 5]),
             "7mer-m8w4w5" : get_site_seq(mirna, 2, 8, wobble = [4, 5]),
             "7mer-A1w4" : get_site_seq(mirna, 1, 7, wobble = [4]),
             "7mer-A1w5" : get_site_seq(mirna, 1, 7, wobble = [5]),
             "7mer-A1w4w5" : get_site_seq(mirna, 1, 7, wobble = [4, 5]),
             "6mer-w4" : get_site_seq(mirna, 2, 7, wobble = [4]),
             "6mer-w5" : get_site_seq(mirna, 2, 7, wobble = [5]),
             "6mer-w4w5" : get_site_seq(mirna, 2, 7, wobble = [4, 5]),
             "6mer-A1w4" : get_site_seq(mirna, 1, 6, wobble = [4]),
             "6mer-A1w5" : get_site_seq(mirna, 1, 6, wobble = [5]),
             "6mer-A1w4w5" : get_site_seq(mirna, 1, 6, wobble = [4, 5]),
             "6mer-m8w4" : get_site_seq(mirna, 3, 8, wobble = [4]),
             "6mer-m8w5" : get_site_seq(mirna, 3, 8, wobble = [5]),
             "6mer-m8w4w5" : get_site_seq(mirna, 3, 8, wobble = [4, 5])})
#miR-1
    if mirna_seq == "UGGAAUGUAAAGAAGUAUGUAU":
        seq_site_map.update(
            {"8mer-w2w3" : get_site_seq(mirna, 1, 8, wobble = [2, 3]),
             "8mer-w2w6" : get_site_seq(mirna, 1, 8, wobble = [2, 6]),
             "8mer-w2w7" : get_site_seq(mirna, 1, 8, wobble = [2, 7]),
             "8mer-w2w8" : get_site_seq(mirna, 1, 8, wobble = [2, 8]),
             "8mer-w3w6" : get_site_seq(mirna, 1, 8, wobble = [3, 6]),
             "8mer-w3w7" : get_site_seq(mirna, 1, 8, wobble = [3, 7]),
             "8mer-w3w8" : get_site_seq(mirna, 1, 8, wobble = [3, 8]),
             "8mer-w6w7" : get_site_seq(mirna, 1, 8, wobble = [6, 7]),
             "8mer-w6w8" : get_site_seq(mirna, 1, 8, wobble = [6, 8]),
             "8mer-w7w8" : get_site_seq(mirna, 1, 8, wobble = [7, 8]),

             "8mer-w2w3w6" : get_site_seq(mirna, 1, 8, wobble = [2, 3, 6]),
             "8mer-w2w3w7" : get_site_seq(mirna, 1, 8, wobble = [2, 3, 7]),
             "8mer-w2w3w8" : get_site_seq(mirna, 1, 8, wobble = [2, 3, 8]),
             "8mer-w2w6w7" : get_site_seq(mirna, 1, 8, wobble = [2, 6, 7]),
             "8mer-w2w6w8" : get_site_seq(mirna, 1, 8, wobble = [2, 6, 8]),
             "8mer-w2w7w8" : get_site_seq(mirna, 1, 8, wobble = [2, 7, 8]),
             "8mer-w3w6w7" : get_site_seq(mirna, 1, 8, wobble = [3, 6, 7]),
             "8mer-w3w6w8" : get_site_seq(mirna, 1, 8, wobble = [3, 6, 8]),
             "8mer-w3w7w8" : get_site_seq(mirna, 1, 8, wobble = [3, 7, 8]),
             "8mer-w6w7w8" : get_site_seq(mirna, 1, 8, wobble = [6, 7, 8]),

             "8mer-w2w3w6w7" : get_site_seq(mirna, 1, 8, wobble = [2, 3, 6, 7]),
             "8mer-w2w3w6w8" : get_site_seq(mirna, 1, 8, wobble = [2, 3, 6, 8]),
             "8mer-w2w3w7w8" : get_site_seq(mirna, 1, 8, wobble = [2, 3, 7, 8]),
             "8mer-w2w6w7w8" : get_site_seq(mirna, 1, 8, wobble = [2, 6, 7, 8]),
             "8mer-w3w6w7w8" : get_site_seq(mirna, 1, 8, wobble = [3, 6, 7, 8]),

             "8mer-w2w3w6w7w8" : get_site_seq(mirna, 1, 8, wobble = [2, 3, 6, 7, 8]),


             "7mer-m8w2w3" : get_site_seq(mirna, 2, 8, wobble = [2, 3]),
             "7mer-m8w2w6" : get_site_seq(mirna, 2, 8, wobble = [2, 6]),
             "7mer-m8w2w7" : get_site_seq(mirna, 2, 8, wobble = [2, 7]),
             "7mer-m8w2w8" : get_site_seq(mirna, 2, 8, wobble = [2, 8]),
             "7mer-m8w3w6" : get_site_seq(mirna, 2, 8, wobble = [3, 6]),
             "7mer-m8w3w7" : get_site_seq(mirna, 2, 8, wobble = [3, 7]),
             "7mer-m8w3w8" : get_site_seq(mirna, 2, 8, wobble = [3, 8]),
             "7mer-m8w6w7" : get_site_seq(mirna, 2, 8, wobble = [6, 7]),
             "7mer-m8w6w8" : get_site_seq(mirna, 2, 8, wobble = [6, 8]),
             "7mer-m8w7w8" : get_site_seq(mirna, 2, 8, wobble = [7, 8]),

             "7mer-m8w2w3w6" : get_site_seq(mirna, 2, 8, wobble = [2, 3, 6]),
             "7mer-m8w2w3w7" : get_site_seq(mirna, 2, 8, wobble = [2, 3, 7]),
             "7mer-m8w2w3w8" : get_site_seq(mirna, 2, 8, wobble = [2, 3, 8]),
             "7mer-m8w2w6w7" : get_site_seq(mirna, 2, 8, wobble = [2, 6, 7]),
             "7mer-m8w2w6w8" : get_site_seq(mirna, 2, 8, wobble = [2, 6, 8]),
             "7mer-m8w2w7w8" : get_site_seq(mirna, 2, 8, wobble = [2, 7, 8]),
             "7mer-m8w3w6w7" : get_site_seq(mirna, 2, 8, wobble = [3, 6, 7]),
             "7mer-m8w3w6w8" : get_site_seq(mirna, 2, 8, wobble = [3, 6, 8]),
             "7mer-m8w3w7w8" : get_site_seq(mirna, 2, 8, wobble = [3, 7, 8]),
             "7mer-m8w6w7w8" : get_site_seq(mirna, 2, 8, wobble = [6, 7, 8]),

             "7mer-m8w2w3w6w7" : get_site_seq(mirna, 2, 8, wobble = [2, 3, 6, 7]),
             "7mer-m8w2w3w6w8" : get_site_seq(mirna, 2, 8, wobble = [2, 3, 6, 8]),
             "7mer-m8w2w3w7w8" : get_site_seq(mirna, 2, 8, wobble = [2, 3, 7, 8]),
             "7mer-m8w2w6w7w8" : get_site_seq(mirna, 2, 8, wobble = [2, 6, 7, 8]),
             "7mer-m8w3w6w7w8" : get_site_seq(mirna, 2, 8, wobble = [3, 6, 7, 8]),
             
             "7mer-m8w2w3w6w7w8" : get_site_seq(mirna, 2, 8, wobble = [2, 3, 6, 7, 8]),

             "7mer-A1w2" : get_site_seq(mirna, 1, 7, wobble = 2),
             "7mer-A1w3" : get_site_seq(mirna, 1, 7, wobble = 3),
             "7mer-A1w6" : get_site_seq(mirna, 1, 7, wobble = 6),
             "7mer-A1w7" : get_site_seq(mirna, 1, 7, wobble = 7),


             "7mer-A1w2w3" : get_site_seq(mirna, 1, 7, wobble = [2, 3]),
             "7mer-A1w2w6" : get_site_seq(mirna, 1, 7, wobble = [2, 6]),
             "7mer-A1w2w7" : get_site_seq(mirna, 1, 7, wobble = [2, 7]),
             "7mer-A1w3w6" : get_site_seq(mirna, 1, 7, wobble = [3, 6]),
             "7mer-A1w3w7" : get_site_seq(mirna, 1, 7, wobble = [3, 7]),
             "7mer-A1w6w7" : get_site_seq(mirna, 1, 7, wobble = [6, 7]),

             "7mer-A1w2w3w6" : get_site_seq(mirna, 1, 7, wobble = [2, 3, 6]),
             "7mer-A1w2w3w7" : get_site_seq(mirna, 1, 7, wobble = [2, 3, 7]),
             "7mer-A1w2w6w7" : get_site_seq(mirna, 1, 7, wobble = [2, 6, 7]),
             "7mer-A1w3w6w7" : get_site_seq(mirna, 1, 7, wobble = [3, 6, 7]),

             "7mer-A1w2w3w6w7" : get_site_seq(mirna, 1, 7, wobble = [2, 3, 6, 7]),

             "6mer-A1w2" : get_site_seq(mirna, 1, 6, wobble = 2),
             "6mer-A1w3" : get_site_seq(mirna, 1, 6, wobble = 3),
             "6mer-A1w6" : get_site_seq(mirna, 1, 6, wobble = 6),

             "6mer-A1w2w3" : get_site_seq(mirna, 1, 6, wobble = [2, 3]),
             "6mer-A1w2w6" : get_site_seq(mirna, 1, 6, wobble = [2, 6]),
             "6mer-A1w3w6" : get_site_seq(mirna, 1, 6, wobble = [3, 6]),

             "6mer-A1w2w3w6" : get_site_seq(mirna, 1, 6, wobble = [2, 3, 6]),

             "6mer-w2" : get_site_seq(mirna, 2, 7, wobble = 2),
             "6mer-w3" : get_site_seq(mirna, 2, 7, wobble = 3),
             "6mer-w6" : get_site_seq(mirna, 2, 7, wobble = 6),
             "6mer-w7" : get_site_seq(mirna, 2, 7, wobble = 7),


             "6mer-w2w3" : get_site_seq(mirna, 2, 7, wobble = [2, 3]),
             "6mer-w2w6" : get_site_seq(mirna, 2, 7, wobble = [2, 6]),
             "6mer-w2w7" : get_site_seq(mirna, 2, 7, wobble = [2, 7]),
             "6mer-w3w6" : get_site_seq(mirna, 2, 7, wobble = [3, 6]),
             "6mer-w3w7" : get_site_seq(mirna, 2, 7, wobble = [3, 7]),
             "6mer-w6w7" : get_site_seq(mirna, 2, 7, wobble = [6, 7]),

             "6mer-w2w3w6" : get_site_seq(mirna, 2, 7, wobble = [2, 3, 6]),
             "6mer-w2w3w7" : get_site_seq(mirna, 2, 7, wobble = [2, 3, 7]),
             "6mer-w2w6w7" : get_site_seq(mirna, 2, 7, wobble = [2, 6, 7]),
             "6mer-w3w6w7" : get_site_seq(mirna, 2, 7, wobble = [3, 6, 7]),

             "6mer-w2w3w6w7" : get_site_seq(mirna, 2, 7, wobble = [2, 3, 6, 7]),

             "6mer-m8w3" : get_site_seq(mirna, 3, 8, wobble = 3),
             "6mer-m8w6" : get_site_seq(mirna, 3, 8, wobble = 6),
             "6mer-m8w7" : get_site_seq(mirna, 3, 8, wobble = 7),
             "6mer-m8w8" : get_site_seq(mirna, 3, 8, wobble = 8),


             "6mer-m8w3w6" : get_site_seq(mirna, 3, 8, wobble = [3, 6]),
             "6mer-m8w3w7" : get_site_seq(mirna, 3, 8, wobble = [3, 7]),
             "6mer-m8w3w8" : get_site_seq(mirna, 3, 8, wobble = [3, 8]),
             "6mer-m8w6w7" : get_site_seq(mirna, 3, 8, wobble = [6, 7]),
             "6mer-m8w6w8" : get_site_seq(mirna, 3, 8, wobble = [6, 8]),
             "6mer-m8w7w8" : get_site_seq(mirna, 3, 8, wobble = [7, 8]),

             "6mer-m8w3w6w7" : get_site_seq(mirna, 3, 8, wobble = [3, 6, 7]),
             "6mer-m8w3w6w8" : get_site_seq(mirna, 3, 8, wobble = [3, 6, 8]),
             "6mer-m8w3w7w8" : get_site_seq(mirna, 3, 8, wobble = [3, 7, 8]),
             "6mer-m8w6w7w8" : get_site_seq(mirna, 3, 8, wobble = [6, 7, 8]),

             "6mer-m8w3w6w7w8" : get_site_seq(mirna, 3, 8, wobble = [3, 6, 7, 8])})

# Wobbles for lsy-6:

    if mirna_seq == "UUUUGUAUGAGACGCAUUUCGA":
        seq_site_map.update(
            {"8mer-w2w3" : get_site_seq(mirna, 1, 8, wobble = [2, 3]),
             "8mer-w2w4" : get_site_seq(mirna, 1, 8, wobble = [2, 4]),
             "8mer-w2w5" : get_site_seq(mirna, 1, 8, wobble = [2, 5]),
             "8mer-w2w6" : get_site_seq(mirna, 1, 8, wobble = [2, 6]),
             "8mer-w2w8" : get_site_seq(mirna, 1, 8, wobble = [2, 8]),
             "8mer-w3w4" : get_site_seq(mirna, 1, 8, wobble = [3, 4]),
             "8mer-w3w5" : get_site_seq(mirna, 1, 8, wobble = [3, 5]),
             "8mer-w3w6" : get_site_seq(mirna, 1, 8, wobble = [3, 6]),
             "8mer-w3w8" : get_site_seq(mirna, 1, 8, wobble = [3, 8]),
             "8mer-w4w5" : get_site_seq(mirna, 1, 8, wobble = [4, 5]),
             "8mer-w4w6" : get_site_seq(mirna, 1, 8, wobble = [4, 6]),
             "8mer-w4w8" : get_site_seq(mirna, 1, 8, wobble = [4, 8]),
             "8mer-w5w6" : get_site_seq(mirna, 1, 8, wobble = [5, 6]),
             "8mer-w5w8" : get_site_seq(mirna, 1, 8, wobble = [5, 8]),
             "8mer-w6w8" : get_site_seq(mirna, 1, 8, wobble = [6, 8]),

             "8mer-w2w3w4" : get_site_seq(mirna, 1, 8, wobble = [2, 3, 4]),
             "8mer-w2w3w5" : get_site_seq(mirna, 1, 8, wobble = [2, 3, 5]),
             "8mer-w2w3w6" : get_site_seq(mirna, 1, 8, wobble = [2, 3, 6]),
             "8mer-w2w3w8" : get_site_seq(mirna, 1, 8, wobble = [2, 3, 8]),
             "8mer-w2w4w5" : get_site_seq(mirna, 1, 8, wobble = [2, 4, 5]),
             "8mer-w2w4w6" : get_site_seq(mirna, 1, 8, wobble = [2, 4, 6]),
             "8mer-w2w4w8" : get_site_seq(mirna, 1, 8, wobble = [2, 4, 8]),
             "8mer-w2w5w6" : get_site_seq(mirna, 1, 8, wobble = [2, 5, 6]),
             "8mer-w2w5w8" : get_site_seq(mirna, 1, 8, wobble = [2, 5, 8]),
             "8mer-w2w6w8" : get_site_seq(mirna, 1, 8, wobble = [2, 6, 8]),
             "8mer-w3w4w5" : get_site_seq(mirna, 1, 8, wobble = [3, 4, 5]),
             "8mer-w3w4w6" : get_site_seq(mirna, 1, 8, wobble = [3, 4, 6]),
             "8mer-w3w4w8" : get_site_seq(mirna, 1, 8, wobble = [3, 4, 8]),
             "8mer-w3w5w6" : get_site_seq(mirna, 1, 8, wobble = [3, 5, 6]),
             "8mer-w3w5w8" : get_site_seq(mirna, 1, 8, wobble = [3, 5, 8]),
             "8mer-w3w6w8" : get_site_seq(mirna, 1, 8, wobble = [3, 6, 8]),
             "8mer-w4w5w6" : get_site_seq(mirna, 1, 8, wobble = [4, 5, 6]),
             "8mer-w4w5w8" : get_site_seq(mirna, 1, 8, wobble = [4, 5, 8]),
             "8mer-w4w6w8" : get_site_seq(mirna, 1, 8, wobble = [4, 6, 8]),
             "8mer-w5w6w8" : get_site_seq(mirna, 1, 8, wobble = [5, 6, 8]),

             "8mer-w2w3w4w5" : get_site_seq(mirna, 1, 8, wobble = [2, 3, 4, 5]),
             "8mer-w2w3w4w6" : get_site_seq(mirna, 1, 8, wobble = [2, 3, 4, 6]),
             "8mer-w2w3w4w8" : get_site_seq(mirna, 1, 8, wobble = [2, 3, 4, 8]),
             "8mer-w2w3w5w6" : get_site_seq(mirna, 1, 8, wobble = [2, 3, 5, 6]),
             "8mer-w2w3w5w8" : get_site_seq(mirna, 1, 8, wobble = [2, 3, 5, 8]),
             "8mer-w2w3w6w8" : get_site_seq(mirna, 1, 8, wobble = [2, 3, 6, 8]),
             "8mer-w2w4w5w6" : get_site_seq(mirna, 1, 8, wobble = [2, 4, 5, 6]),
             "8mer-w2w4w5w8" : get_site_seq(mirna, 1, 8, wobble = [2, 4, 5, 8]),
             "8mer-w2w4w6w8" : get_site_seq(mirna, 1, 8, wobble = [2, 4, 6, 8]),
             "8mer-w2w5w6w8" : get_site_seq(mirna, 1, 8, wobble = [2, 5, 6, 8]),
             "8mer-w3w4w5w6" : get_site_seq(mirna, 1, 8, wobble = [3, 4, 5, 6]),
             "8mer-w3w4w5w8" : get_site_seq(mirna, 1, 8, wobble = [3, 4, 5, 8]),
             "8mer-w3w4w6w8" : get_site_seq(mirna, 1, 8, wobble = [3, 4, 6, 8]),
             "8mer-w3w5w6w8" : get_site_seq(mirna, 1, 8, wobble = [3, 5, 6, 8]),
             "8mer-w4w5w6w8" : get_site_seq(mirna, 1, 8, wobble = [4, 5, 6, 8]),

             "8mer-w2w3w4w5w6" : get_site_seq(mirna, 1, 8, wobble = [2, 3, 4, 5, 6]),
             "8mer-w2w3w4w5w8" : get_site_seq(mirna, 1, 8, wobble = [2, 3, 4, 5, 8]),
             "8mer-w2w3w4w6w8" : get_site_seq(mirna, 1, 8, wobble = [2, 3, 4, 6, 8]),
             "8mer-w2w3w5w6w8" : get_site_seq(mirna, 1, 8, wobble = [2, 3, 5, 6, 8]),
             "8mer-w2w4w5w6w8" : get_site_seq(mirna, 1, 8, wobble = [2, 4, 5, 6, 8]),
             "8mer-w3w4w5w6w8" : get_site_seq(mirna, 1, 8, wobble = [3, 4, 5, 6, 8]),

             "8mer-w2w3w4w5w6w8" : get_site_seq(mirna, 1, 8, wobble = [2, 3, 4, 5, 6, 8]),


             "7mer-m8w2w3" : get_site_seq(mirna, 2, 8, wobble = [2, 3]),
             "7mer-m8w2w4" : get_site_seq(mirna, 2, 8, wobble = [2, 4]),
             "7mer-m8w2w5" : get_site_seq(mirna, 2, 8, wobble = [2, 5]),
             "7mer-m8w2w6" : get_site_seq(mirna, 2, 8, wobble = [2, 6]),
             "7mer-m8w2w8" : get_site_seq(mirna, 2, 8, wobble = [2, 8]),
             "7mer-m8w3w4" : get_site_seq(mirna, 2, 8, wobble = [3, 4]),
             "7mer-m8w3w5" : get_site_seq(mirna, 2, 8, wobble = [3, 5]),
             "7mer-m8w3w6" : get_site_seq(mirna, 2, 8, wobble = [3, 6]),
             "7mer-m8w3w8" : get_site_seq(mirna, 2, 8, wobble = [3, 8]),
             "7mer-m8w4w5" : get_site_seq(mirna, 2, 8, wobble = [4, 5]),
             "7mer-m8w4w6" : get_site_seq(mirna, 2, 8, wobble = [4, 6]),
             "7mer-m8w4w8" : get_site_seq(mirna, 2, 8, wobble = [4, 8]),
             "7mer-m8w5w6" : get_site_seq(mirna, 2, 8, wobble = [5, 6]),
             "7mer-m8w5w8" : get_site_seq(mirna, 2, 8, wobble = [5, 8]),
             "7mer-m8w6w8" : get_site_seq(mirna, 2, 8, wobble = [6, 8]),

             "7mer-m8w2w3w4" : get_site_seq(mirna, 2, 8, wobble = [2, 3, 4]),
             "7mer-m8w2w3w5" : get_site_seq(mirna, 2, 8, wobble = [2, 3, 5]),
             "7mer-m8w2w3w6" : get_site_seq(mirna, 2, 8, wobble = [2, 3, 6]),
             "7mer-m8w2w3w8" : get_site_seq(mirna, 2, 8, wobble = [2, 3, 8]),
             "7mer-m8w2w4w5" : get_site_seq(mirna, 2, 8, wobble = [2, 4, 5]),
             "7mer-m8w2w4w6" : get_site_seq(mirna, 2, 8, wobble = [2, 4, 6]),
             "7mer-m8w2w4w8" : get_site_seq(mirna, 2, 8, wobble = [2, 4, 8]),
             "7mer-m8w2w5w6" : get_site_seq(mirna, 2, 8, wobble = [2, 5, 6]),
             "7mer-m8w2w5w8" : get_site_seq(mirna, 2, 8, wobble = [2, 5, 8]),
             "7mer-m8w2w6w8" : get_site_seq(mirna, 2, 8, wobble = [2, 6, 8]),
             "7mer-m8w3w4w5" : get_site_seq(mirna, 2, 8, wobble = [3, 4, 5]),
             "7mer-m8w3w4w6" : get_site_seq(mirna, 2, 8, wobble = [3, 4, 6]),
             "7mer-m8w3w4w8" : get_site_seq(mirna, 2, 8, wobble = [3, 4, 8]),
             "7mer-m8w3w5w6" : get_site_seq(mirna, 2, 8, wobble = [3, 5, 6]),
             "7mer-m8w3w5w8" : get_site_seq(mirna, 2, 8, wobble = [3, 5, 8]),
             "7mer-m8w3w6w8" : get_site_seq(mirna, 2, 8, wobble = [3, 6, 8]),
             "7mer-m8w4w5w6" : get_site_seq(mirna, 2, 8, wobble = [4, 5, 6]),
             "7mer-m8w4w5w8" : get_site_seq(mirna, 2, 8, wobble = [4, 5, 8]),
             "7mer-m8w4w6w8" : get_site_seq(mirna, 2, 8, wobble = [4, 6, 8]),
             "7mer-m8w5w6w8" : get_site_seq(mirna, 2, 8, wobble = [5, 6, 8]),

             "7mer-m8w2w3w4w5" : get_site_seq(mirna, 2, 8, wobble = [2, 3, 4, 5]),
             "7mer-m8w2w3w4w6" : get_site_seq(mirna, 2, 8, wobble = [2, 3, 4, 6]),
             "7mer-m8w2w3w4w8" : get_site_seq(mirna, 2, 8, wobble = [2, 3, 4, 8]),
             "7mer-m8w2w3w5w6" : get_site_seq(mirna, 2, 8, wobble = [2, 3, 5, 6]),
             "7mer-m8w2w3w5w8" : get_site_seq(mirna, 2, 8, wobble = [2, 3, 5, 8]),
             "7mer-m8w2w3w6w8" : get_site_seq(mirna, 2, 8, wobble = [2, 3, 6, 8]),
             "7mer-m8w2w4w5w6" : get_site_seq(mirna, 2, 8, wobble = [2, 4, 5, 6]),
             "7mer-m8w2w4w5w8" : get_site_seq(mirna, 2, 8, wobble = [2, 4, 5, 8]),
             "7mer-m8w2w4w6w8" : get_site_seq(mirna, 2, 8, wobble = [2, 4, 6, 8]),
             "7mer-m8w2w5w6w8" : get_site_seq(mirna, 2, 8, wobble = [2, 5, 6, 8]),
             "7mer-m8w3w4w5w6" : get_site_seq(mirna, 2, 8, wobble = [3, 4, 5, 6]),
             "7mer-m8w3w4w5w8" : get_site_seq(mirna, 2, 8, wobble = [3, 4, 5, 8]),
             "7mer-m8w3w4w6w8" : get_site_seq(mirna, 2, 8, wobble = [3, 4, 6, 8]),
             "7mer-m8w3w5w6w8" : get_site_seq(mirna, 2, 8, wobble = [3, 5, 6, 8]),
             "7mer-m8w4w5w6w8" : get_site_seq(mirna, 2, 8, wobble = [4, 5, 6, 8]),

             "7mer-m8w2w3w4w5w6" : get_site_seq(mirna, 2, 8, wobble = [2, 3, 4, 5, 6]),
             "7mer-m8w2w3w4w5w8" : get_site_seq(mirna, 2, 8, wobble = [2, 3, 4, 5, 8]),
             "7mer-m8w2w3w4w6w8" : get_site_seq(mirna, 2, 8, wobble = [2, 3, 4, 6, 8]),
             "7mer-m8w2w3w5w6w8" : get_site_seq(mirna, 2, 8, wobble = [2, 3, 5, 6, 8]),
             "7mer-m8w2w4w5w6w8" : get_site_seq(mirna, 2, 8, wobble = [2, 4, 5, 6, 8]),
             "7mer-m8w3w4w5w6w8" : get_site_seq(mirna, 2, 8, wobble = [3, 4, 5, 6, 8]),

             "7mer-m8w2w3w4w5w6w8" : get_site_seq(mirna, 2, 8, wobble = [2, 3, 4, 5, 6, 8]),

             "7mer-A1w2" : get_site_seq(mirna, 1, 7, wobble = 2),
             "7mer-A1w3" : get_site_seq(mirna, 1, 7, wobble = 3),
             "7mer-A1w4" : get_site_seq(mirna, 1, 7, wobble = 4),
             "7mer-A1w5" : get_site_seq(mirna, 1, 7, wobble = 5),
             "7mer-A1w6" : get_site_seq(mirna, 1, 7, wobble = 6),


             "7mer-A1w2w3" : get_site_seq(mirna, 1, 7, wobble = [2, 3]),
             "7mer-A1w2w4" : get_site_seq(mirna, 1, 7, wobble = [2, 4]),
             "7mer-A1w2w5" : get_site_seq(mirna, 1, 7, wobble = [2, 5]),
             "7mer-A1w2w6" : get_site_seq(mirna, 1, 7, wobble = [2, 6]),
             "7mer-A1w3w4" : get_site_seq(mirna, 1, 7, wobble = [3, 4]),
             "7mer-A1w3w5" : get_site_seq(mirna, 1, 7, wobble = [3, 5]),
             "7mer-A1w3w6" : get_site_seq(mirna, 1, 7, wobble = [3, 6]),
             "7mer-A1w4w5" : get_site_seq(mirna, 1, 7, wobble = [4, 5]),
             "7mer-A1w4w6" : get_site_seq(mirna, 1, 7, wobble = [4, 6]),
             "7mer-A1w5w6" : get_site_seq(mirna, 1, 7, wobble = [5, 6]),

             "7mer-A1w2w3w4" : get_site_seq(mirna, 1, 7, wobble = [2, 3, 4]),
             "7mer-A1w2w3w5" : get_site_seq(mirna, 1, 7, wobble = [2, 3, 5]),
             "7mer-A1w2w3w6" : get_site_seq(mirna, 1, 7, wobble = [2, 3, 6]),
             "7mer-A1w2w4w5" : get_site_seq(mirna, 1, 7, wobble = [2, 4, 5]),
             "7mer-A1w2w4w6" : get_site_seq(mirna, 1, 7, wobble = [2, 4, 6]),
             "7mer-A1w2w5w6" : get_site_seq(mirna, 1, 7, wobble = [2, 5, 6]),
             "7mer-A1w3w4w5" : get_site_seq(mirna, 1, 7, wobble = [3, 4, 5]),
             "7mer-A1w3w4w6" : get_site_seq(mirna, 1, 7, wobble = [3, 4, 6]),
             "7mer-A1w3w5w6" : get_site_seq(mirna, 1, 7, wobble = [3, 5, 6]),
             "7mer-A1w4w5w6" : get_site_seq(mirna, 1, 7, wobble = [4, 5, 6]),

             "7mer-A1w2w3w4w5" : get_site_seq(mirna, 1, 7, wobble = [2, 3, 4, 5]),
             "7mer-A1w2w3w4w6" : get_site_seq(mirna, 1, 7, wobble = [2, 3, 4, 6]),
             "7mer-A1w2w3w5w6" : get_site_seq(mirna, 1, 7, wobble = [2, 3, 5, 6]),
             "7mer-A1w2w4w5w6" : get_site_seq(mirna, 1, 7, wobble = [2, 4, 5, 6]),
             "7mer-A1w3w4w5w6" : get_site_seq(mirna, 1, 7, wobble = [3, 4, 5, 6]),

             "7mer-A1w2w3w4w5w6" : get_site_seq(mirna, 1, 7, wobble = [2, 3, 4, 5, 6]),


             "6mer-A1w2w3" : get_site_seq(mirna, 1, 6, wobble = [2, 3]),
             "6mer-A1w2w4" : get_site_seq(mirna, 1, 6, wobble = [2, 4]),
             "6mer-A1w2w5" : get_site_seq(mirna, 1, 6, wobble = [2, 5]),
             "6mer-A1w2w6" : get_site_seq(mirna, 1, 6, wobble = [2, 6]),
             "6mer-A1w3w4" : get_site_seq(mirna, 1, 6, wobble = [3, 4]),
             "6mer-A1w3w5" : get_site_seq(mirna, 1, 6, wobble = [3, 5]),
             "6mer-A1w3w6" : get_site_seq(mirna, 1, 6, wobble = [3, 6]),
             "6mer-A1w4w5" : get_site_seq(mirna, 1, 6, wobble = [4, 5]),
             "6mer-A1w4w6" : get_site_seq(mirna, 1, 6, wobble = [4, 6]),
             "6mer-A1w5w6" : get_site_seq(mirna, 1, 6, wobble = [5, 6]),

             "6mer-A1w2w3w4" : get_site_seq(mirna, 1, 6, wobble = [2, 3, 4]),
             "6mer-A1w2w3w5" : get_site_seq(mirna, 1, 6, wobble = [2, 3, 5]),
             "6mer-A1w2w3w6" : get_site_seq(mirna, 1, 6, wobble = [2, 3, 6]),
             "6mer-A1w2w4w5" : get_site_seq(mirna, 1, 6, wobble = [2, 4, 5]),
             "6mer-A1w2w4w6" : get_site_seq(mirna, 1, 6, wobble = [2, 4, 6]),
             "6mer-A1w2w5w6" : get_site_seq(mirna, 1, 6, wobble = [2, 5, 6]),
             "6mer-A1w3w4w5" : get_site_seq(mirna, 1, 6, wobble = [3, 4, 5]),
             "6mer-A1w3w4w6" : get_site_seq(mirna, 1, 6, wobble = [3, 4, 6]),
             "6mer-A1w3w5w6" : get_site_seq(mirna, 1, 6, wobble = [3, 5, 6]),
             "6mer-A1w4w5w6" : get_site_seq(mirna, 1, 6, wobble = [4, 5, 6]),

             "6mer-A1w2w3w4w5" : get_site_seq(mirna, 1, 6, wobble = [2, 3, 4, 5]),
             "6mer-A1w2w3w4w6" : get_site_seq(mirna, 1, 6, wobble = [2, 3, 4, 6]),
             "6mer-A1w2w3w5w6" : get_site_seq(mirna, 1, 6, wobble = [2, 3, 5, 6]),
             "6mer-A1w2w4w5w6" : get_site_seq(mirna, 1, 6, wobble = [2, 4, 5, 6]),
             "6mer-A1w3w4w5w6" : get_site_seq(mirna, 1, 6, wobble = [3, 4, 5, 6]),

             "6mer-A1w2w3w4w5w6" : get_site_seq(mirna, 1, 6, wobble = [2, 3, 4, 5, 6]),

             "6mer-m8w3" : get_site_seq(mirna, 3, 8, wobble = 3),
             "6mer-m8w4" : get_site_seq(mirna, 3, 8, wobble = 4),
             "6mer-m8w5" : get_site_seq(mirna, 3, 8, wobble = 5),
             "6mer-m8w6" : get_site_seq(mirna, 3, 8, wobble = 6),
             "6mer-m8w8" : get_site_seq(mirna, 3, 8, wobble = 8),

             "6mer-w2" : get_site_seq(mirna, 2, 7, wobble = 2),
             "6mer-w3" : get_site_seq(mirna, 2, 7, wobble = 3),
             "6mer-w4" : get_site_seq(mirna, 2, 7, wobble = 4),
             "6mer-w5" : get_site_seq(mirna, 2, 7, wobble = 5),
             "6mer-w6" : get_site_seq(mirna, 2, 7, wobble = 6),

             "6mer-A1w2" : get_site_seq(mirna, 1, 6, wobble = 2),
             "6mer-A1w3" : get_site_seq(mirna, 1, 6, wobble = 3),
             "6mer-A1w4" : get_site_seq(mirna, 1, 6, wobble = 4),
             "6mer-A1w5" : get_site_seq(mirna, 1, 6, wobble = 5),
             "6mer-A1w6" : get_site_seq(mirna, 1, 6, wobble = 6),

             "6mer-w2w3" : get_site_seq(mirna, 2, 7, wobble = [2, 3]),
             "6mer-w2w4" : get_site_seq(mirna, 2, 7, wobble = [2, 4]),
             "6mer-w2w5" : get_site_seq(mirna, 2, 7, wobble = [2, 5]),
             "6mer-w2w6" : get_site_seq(mirna, 2, 7, wobble = [2, 6]),
             "6mer-w3w4" : get_site_seq(mirna, 2, 7, wobble = [3, 4]),
             "6mer-w3w5" : get_site_seq(mirna, 2, 7, wobble = [3, 5]),
             "6mer-w3w6" : get_site_seq(mirna, 2, 7, wobble = [3, 6]),
             "6mer-w4w5" : get_site_seq(mirna, 2, 7, wobble = [4, 5]),
             "6mer-w4w6" : get_site_seq(mirna, 2, 7, wobble = [4, 6]),
             "6mer-w5w6" : get_site_seq(mirna, 2, 7, wobble = [5, 6]),

             "6mer-w2w3w4" : get_site_seq(mirna, 2, 7, wobble = [2, 3, 4]),
             "6mer-w2w3w5" : get_site_seq(mirna, 2, 7, wobble = [2, 3, 5]),
             "6mer-w2w3w6" : get_site_seq(mirna, 2, 7, wobble = [2, 3, 6]),
             "6mer-w2w4w5" : get_site_seq(mirna, 2, 7, wobble = [2, 4, 5]),
             "6mer-w2w4w6" : get_site_seq(mirna, 2, 7, wobble = [2, 4, 6]),
             "6mer-w2w5w6" : get_site_seq(mirna, 2, 7, wobble = [2, 5, 6]),
             "6mer-w3w4w5" : get_site_seq(mirna, 2, 7, wobble = [3, 4, 5]),
             "6mer-w3w4w6" : get_site_seq(mirna, 2, 7, wobble = [3, 4, 6]),
             "6mer-w3w5w6" : get_site_seq(mirna, 2, 7, wobble = [3, 5, 6]),
             "6mer-w4w5w6" : get_site_seq(mirna, 2, 7, wobble = [4, 5, 6]),

             "6mer-w2w3w4w5" : get_site_seq(mirna, 2, 7, wobble = [2, 3, 4, 5]),
             "6mer-w2w3w4w6" : get_site_seq(mirna, 2, 7, wobble = [2, 3, 4, 6]),
             "6mer-w2w3w5w6" : get_site_seq(mirna, 2, 7, wobble = [2, 3, 5, 6]),
             "6mer-w2w4w5w6" : get_site_seq(mirna, 2, 7, wobble = [2, 4, 5, 6]),
             "6mer-w3w4w5w6" : get_site_seq(mirna, 2, 7, wobble = [3, 4, 5, 6]),

             "6mer-w2w3w4w5w6" : get_site_seq(mirna, 2, 7, wobble = [2, 3, 4, 5, 6]),


             "6mer-m8w3w4" : get_site_seq(mirna, 3, 8, wobble = [3, 4]),
             "6mer-m8w3w5" : get_site_seq(mirna, 3, 8, wobble = [3, 5]),
             "6mer-m8w3w6" : get_site_seq(mirna, 3, 8, wobble = [3, 6]),
             "6mer-m8w3w8" : get_site_seq(mirna, 3, 8, wobble = [3, 8]),
             "6mer-m8w4w5" : get_site_seq(mirna, 3, 8, wobble = [4, 5]),
             "6mer-m8w4w6" : get_site_seq(mirna, 3, 8, wobble = [4, 6]),
             "6mer-m8w4w8" : get_site_seq(mirna, 3, 8, wobble = [4, 8]),
             "6mer-m8w5w6" : get_site_seq(mirna, 3, 8, wobble = [5, 6]),
             "6mer-m8w5w8" : get_site_seq(mirna, 3, 8, wobble = [5, 8]),
             "6mer-m8w6w8" : get_site_seq(mirna, 3, 8, wobble = [6, 8]),

             "6mer-m8w3w4w5" : get_site_seq(mirna, 3, 8, wobble = [3, 4, 5]),
             "6mer-m8w3w4w6" : get_site_seq(mirna, 3, 8, wobble = [3, 4, 6]),
             "6mer-m8w3w4w8" : get_site_seq(mirna, 3, 8, wobble = [3, 4, 8]),
             "6mer-m8w3w5w6" : get_site_seq(mirna, 3, 8, wobble = [3, 5, 6]),
             "6mer-m8w3w5w8" : get_site_seq(mirna, 3, 8, wobble = [3, 5, 8]),
             "6mer-m8w3w6w8" : get_site_seq(mirna, 3, 8, wobble = [3, 6, 8]),
             "6mer-m8w4w5w6" : get_site_seq(mirna, 3, 8, wobble = [4, 5, 6]),
             "6mer-m8w4w5w8" : get_site_seq(mirna, 3, 8, wobble = [4, 5, 8]),
             "6mer-m8w4w6w8" : get_site_seq(mirna, 3, 8, wobble = [4, 6, 8]),
             "6mer-m8w5w6w8" : get_site_seq(mirna, 3, 8, wobble = [5, 6, 8]),

             "6mer-m8w3w4w5w6" : get_site_seq(mirna, 3, 8, wobble = [3, 4, 5, 6]),
             "6mer-m8w3w4w5w8" : get_site_seq(mirna, 3, 8, wobble = [3, 4, 5, 8]),
             "6mer-m8w3w4w6w8" : get_site_seq(mirna, 3, 8, wobble = [3, 4, 6, 8]),
             "6mer-m8w3w5w6w8" : get_site_seq(mirna, 3, 8, wobble = [3, 5, 6, 8]),
             "6mer-m8w4w5w6w8" : get_site_seq(mirna, 3, 8, wobble = [4, 5, 6, 8]),

             "6mer-m8w3w4w5w6w8" : get_site_seq(mirna, 3, 8, wobble = [3, 4, 5, 6, 8])})


    if mirna_seq[5] in ["G", "U"]:
        seq_site_map.update(
            {"7mer-m8w6bC6" : get_site_seq(mirna, 2, 8, wobble = 6,
                                           bulge = [6, "C"]),
            "8mer-w6bC6" : get_site_seq(mirna, 1, 8, wobble = 6,
                                           bulge = [6, "C"]),
            "7mer-A1w6bC6" : get_site_seq(mirna, 1, 7, wobble = 6,
                                           bulge = [6, "C"]),
            "6mer-w6bC6" : get_site_seq(mirna, 2, 7, wobble = 6,
                                           bulge = [6, "C"]),
            "6mer-A1w6bC6" : get_site_seq(mirna, 1, 6, wobble = 6,
                                           bulge = [6, "C"]),
            "6mer-m8w6bC6" : get_site_seq(mirna, 3, 8, wobble = 6,
                                           bulge = [6, "C"])})
    if mirna_seq[8] in ["G", "U"]:
        seq_site_map.update(
            {"13mer-m8.20w9": get_site_seq(mirna, 8, 20, wobble = 9),
             "13mer-m9.21w9": get_site_seq(mirna, 9, 21, wobble = 9),
             "12mer-m8.19w9": get_site_seq(mirna, 8, 19, wobble = 9),
             "12mer-m9.20w9": get_site_seq(mirna, 9, 20, wobble = 9),
             "11mer-m8.18w9": get_site_seq(mirna, 8, 18, wobble = 9),
             "11mer-m9.19w9": get_site_seq(mirna, 9, 19, wobble = 9),
             "10mer-m8.17w9": get_site_seq(mirna, 8, 17, wobble = 9),
             "10mer-m9.18w9": get_site_seq(mirna, 9, 18, wobble = 9),
             "9mer-m8.16w9":  get_site_seq(mirna, 8, 16, wobble = 9),
             "9mer-m9.17w9":  get_site_seq(mirna, 9, 17, wobble = 9),
             "8mer-m8.15w9":  get_site_seq(mirna, 8, 15, wobble = 9),
             "8mer-m9.16w9":  get_site_seq(mirna, 9, 16, wobble = 9),
             "7mer-m8.14w9":  get_site_seq(mirna, 8, 14, wobble = 9),
             "7mer-m9.15w9":  get_site_seq(mirna, 9, 15, wobble = 9)})

    if mirna_seq[17] in ["G", "U"]:
        seq_site_map.update(
            {"13mer-m8.20w18": get_site_seq(mirna, 8, 20, wobble = 18),
             "13mer-m9.21w18": get_site_seq(mirna, 9, 21, wobble = 18),
             "13mer-m10.22w18": get_site_seq(mirna, 10, 22, wobble = 18),
             "12mer-m8.19w18": get_site_seq(mirna, 8, 19, wobble = 18),
             "12mer-m9.20w18": get_site_seq(mirna, 9, 20, wobble = 18),
             "12mer-m10.21w18": get_site_seq(mirna, 10, 21, wobble = 18),
             "12mer-m11.22w18": get_site_seq(mirna, 11, 22, wobble = 18),
             "11mer-m8.18w18": get_site_seq(mirna, 8, 18, wobble = 18),
             "11mer-m9.19w18": get_site_seq(mirna, 9, 19, wobble = 18),
             "11mer-m10.20w18": get_site_seq(mirna, 10, 20, wobble = 18),
             "11mer-m11.21w18": get_site_seq(mirna, 11, 21, wobble = 18),
             "11mer-m12.22w18": get_site_seq(mirna, 12, 22, wobble = 18),
             "10mer-m9.18w18": get_site_seq(mirna, 9, 18, wobble = 18),
             "10mer-m10.19w18": get_site_seq(mirna, 10, 19, wobble = 18),
             "10mer-m11.20w18": get_site_seq(mirna, 11, 20, wobble = 18),
             "10mer-m12.21w18": get_site_seq(mirna, 12, 21, wobble = 18),
             "10mer-m13.22w18": get_site_seq(mirna, 13, 22, wobble = 18),
             "9mer-m10.18w18": get_site_seq(mirna, 10, 18, wobble = 18),
             "9mer-m11.19w18": get_site_seq(mirna, 11, 19, wobble = 18),
             "9mer-m12.20w18": get_site_seq(mirna, 12, 20, wobble = 18),
             "9mer-m13.21w18": get_site_seq(mirna, 13, 21, wobble = 18),
             "9mer-m14.22w18": get_site_seq(mirna, 14, 22, wobble = 18),
             "8mer-m11.18w18": get_site_seq(mirna, 11, 18, wobble = 18),
             "8mer-m12.19w18": get_site_seq(mirna, 12, 19, wobble = 18),
             "8mer-m13.20w18": get_site_seq(mirna, 13, 20, wobble = 18),
             "8mer-m14.21w18": get_site_seq(mirna, 14, 21, wobble = 18),
             "8mer-m15.22w18": get_site_seq(mirna, 15, 22, wobble = 18),
             "7mer-m12.18w18": get_site_seq(mirna, 12, 18, wobble = 18),
             "7mer-m13.19w18": get_site_seq(mirna, 13, 19, wobble = 18),
             "7mer-m14.20w18": get_site_seq(mirna, 14, 20, wobble = 18),
             "7mer-m15.21w18": get_site_seq(mirna, 15, 21, wobble = 18),
             "7mer-m16.22w18": get_site_seq(mirna, 16, 22, wobble = 18)})

    if mirna_seq[17] in ["G", "U"] and mirna_seq[8] in ["G", "U"]:
        seq_site_map.update(
            {"13mer-m8.20w9w18": get_site_seq(mirna, 8, 20, wobble = [9,18]),
             "13mer-m9.21w9w18": get_site_seq(mirna, 9, 21, wobble = [9,18]),
             "12mer-m8.19w9w18": get_site_seq(mirna, 8, 19, wobble = [9,18]),
             "12mer-m9.20w9w18": get_site_seq(mirna, 9, 20, wobble = [9,18]),
             "11mer-m8.18w9w18": get_site_seq(mirna, 8, 18, wobble = [9,18]),
             "11mer-m9.19w9w18": get_site_seq(mirna, 9, 19, wobble = [9,18]),
             "10mer-m9.18w9w18": get_site_seq(mirna, 9, 18, wobble = [9,18])})


    if mirna_seq[12] in ["G", "U"]:
        seq_site_map.update(
            {"13mer-m8.20w13" : get_site_seq(mirna, 8, 20, wobble = 13),
             "13mer-m9.21w13" : get_site_seq(mirna, 9, 21, wobble = 13),
             "13mer-m10.22w13" : get_site_seq(mirna, 10, 22, wobble = 13),

             "12mer-m8.19w13" : get_site_seq(mirna, 8, 19, wobble = 13),
             "12mer-m9.20w13" : get_site_seq(mirna, 9, 20, wobble = 13),
             "12mer-m10.21w13" : get_site_seq(mirna, 10, 21, wobble = 13),
             "12mer-m11.22w13" : get_site_seq(mirna, 11, 22, wobble = 13),

             "11mer-m8.18w13" : get_site_seq(mirna, 8, 18, wobble = 13),
             "11mer-m9.19w13" : get_site_seq(mirna, 9, 19, wobble = 13),
             "11mer-m10.20w13" : get_site_seq(mirna, 10, 20, wobble = 13),
             "11mer-m11.21w13" : get_site_seq(mirna, 11, 21, wobble = 13),
             "11mer-m12.22w13" : get_site_seq(mirna, 12, 22, wobble = 13),

             "10mer-m8.17w13" : get_site_seq(mirna, 8, 17, wobble = 13),
             "10mer-m9.18w13" : get_site_seq(mirna, 9, 18, wobble = 13),
             "10mer-m10.19w13" : get_site_seq(mirna, 10, 19, wobble = 13),
             "10mer-m11.20w13" : get_site_seq(mirna, 11, 20, wobble = 13),
             "10mer-m12.21w13" : get_site_seq(mirna, 12, 21, wobble = 13),
             "10mer-m13.22w13" : get_site_seq(mirna, 13, 22, wobble = 13),

             "9mer-m8.16w13" : get_site_seq(mirna, 8, 16, wobble = 13),
             "9mer-m9.17w13" : get_site_seq(mirna, 9, 17, wobble = 13),
             "9mer-m10.18w13" : get_site_seq(mirna, 10, 18, wobble = 13),
             "9mer-m11.19w13" : get_site_seq(mirna, 11, 19, wobble = 13),
             "9mer-m12.20w13" : get_site_seq(mirna, 12, 20, wobble = 13),
             "9mer-m13.21w13" : get_site_seq(mirna, 13, 21, wobble = 13),

             "8mer-m8.15w13" : get_site_seq(mirna, 8, 15, wobble = 13),
             "8mer-m9.16w13" : get_site_seq(mirna, 9, 16, wobble = 13),
             "8mer-m10.17w13" : get_site_seq(mirna, 10, 17, wobble = 13),
             "8mer-m11.18w13" : get_site_seq(mirna, 11, 18, wobble = 13),
             "8mer-m12.19w13" : get_site_seq(mirna, 12, 19, wobble = 13),
             "8mer-m13.20w13" : get_site_seq(mirna, 13, 20, wobble = 13),

             "7mer-m8.14w13" : get_site_seq(mirna, 8, 14, wobble = 13),
             "7mer-m9.15w13" : get_site_seq(mirna, 9, 15, wobble = 13),
             "7mer-m10.16w13" : get_site_seq(mirna, 10, 16, wobble = 13),
             "7mer-m11.17w13" : get_site_seq(mirna, 11, 17, wobble = 13),
             "7mer-m12.18w13" : get_site_seq(mirna, 12, 18, wobble = 13),
             "7mer-m13.19w13" : get_site_seq(mirna, 13, 19, wobble = 13)})
        if len(mirna_seq) == 23:
            seq_site_map.update(
                {"13mer-m11.23w13" : get_site_seq(mirna, 11, 23, wobble = 13),
                 "12mer-m12.23w13" : get_site_seq(mirna, 12, 23, wobble = 13),
                 "11mer-m13.23w13" : get_site_seq(mirna, 13, 23, wobble = 13)})
    if mirna_seq[13] in ["G", "U"]:
        seq_site_map.update(
            {"13mer-m8.20w14" : get_site_seq(mirna, 8, 20, wobble = 14),
             "13mer-m9.21w14" : get_site_seq(mirna, 9, 21, wobble = 14),
             "13mer-m10.22w14" : get_site_seq(mirna, 10, 22, wobble = 14),

             "12mer-m8.19w14" : get_site_seq(mirna, 8, 19, wobble = 14),
             "12mer-m9.20w14" : get_site_seq(mirna, 9, 20, wobble = 14),
             "12mer-m10.21w14" : get_site_seq(mirna, 10, 21, wobble = 14),
             "12mer-m11.22w14" : get_site_seq(mirna, 11, 22, wobble = 14),

             "11mer-m8.18w14" : get_site_seq(mirna, 8, 18, wobble = 14),
             "11mer-m9.19w14" : get_site_seq(mirna, 9, 19, wobble = 14),
             "11mer-m10.20w14" : get_site_seq(mirna, 10, 20, wobble = 14),
             "11mer-m11.21w14" : get_site_seq(mirna, 11, 21, wobble = 14),
             "11mer-m12.22w14" : get_site_seq(mirna, 12, 22, wobble = 14),

             "10mer-m8.17w14" : get_site_seq(mirna, 8, 17, wobble = 14),
             "10mer-m9.18w14" : get_site_seq(mirna, 9, 18, wobble = 14),
             "10mer-m10.19w14" : get_site_seq(mirna, 10, 19, wobble = 14),
             "10mer-m11.20w14" : get_site_seq(mirna, 11, 20, wobble = 14),
             "10mer-m12.21w14" : get_site_seq(mirna, 12, 21, wobble = 14),
             "10mer-m13.22w14" : get_site_seq(mirna, 13, 22, wobble = 14),

             "9mer-m8.16w14" : get_site_seq(mirna, 8, 16, wobble = 14),
             "9mer-m9.17w14" : get_site_seq(mirna, 9, 17, wobble = 14),
             "9mer-m10.18w14" : get_site_seq(mirna, 10, 18, wobble = 14),
             "9mer-m11.19w14" : get_site_seq(mirna, 11, 19, wobble = 14),
             "9mer-m12.20w14" : get_site_seq(mirna, 12, 20, wobble = 14),
             "9mer-m13.21w14" : get_site_seq(mirna, 13, 21, wobble = 14),
             "9mer-m14.22w14" : get_site_seq(mirna, 14, 22, wobble = 14),

             "8mer-m8.15w14" : get_site_seq(mirna, 8, 15, wobble = 14),
             "8mer-m9.16w14" : get_site_seq(mirna, 9, 16, wobble = 14),
             "8mer-m10.17w14" : get_site_seq(mirna, 10, 17, wobble = 14),
             "8mer-m11.18w14" : get_site_seq(mirna, 11, 18, wobble = 14),
             "8mer-m12.19w14" : get_site_seq(mirna, 12, 19, wobble = 14),
             "8mer-m13.20w14" : get_site_seq(mirna, 13, 20, wobble = 14),
             "8mer-m14.21w14" : get_site_seq(mirna, 14, 21, wobble = 14),

             "7mer-m8.14w14" : get_site_seq(mirna, 8, 14, wobble = 14),
             "7mer-m9.15w14" : get_site_seq(mirna, 9, 15, wobble = 14),
             "7mer-m10.16w14" : get_site_seq(mirna, 10, 16, wobble = 14),
             "7mer-m11.17w14" : get_site_seq(mirna, 11, 17, wobble = 14),
             "7mer-m12.18w14" : get_site_seq(mirna, 12, 18, wobble = 14),
             "7mer-m13.19w14" : get_site_seq(mirna, 13, 19, wobble = 14),
             "7mer-m14.20w14" : get_site_seq(mirna, 14, 20, wobble = 14)})        
        if len(mirna_seq) == 23:
            seq_site_map.update(
                {"13mer-m11.23w14" : get_site_seq(mirna, 11, 23, wobble = 14),
                 "12mer-m12.23w14" : get_site_seq(mirna, 12, 23, wobble = 14),
                 "11mer-m13.23w14" : get_site_seq(mirna, 13, 23, wobble = 14),
                 "10mer-m14.23w14" : get_site_seq(mirna, 14, 23, wobble = 14)})



    # Restric to just those seen in the data.
    seq_site_map = {
        key: value for key, value in seq_site_map.items() if key in sites}
    # Add those sites for which a site name doesn't exist (ie.e "CGACTTTA").
    seq_site_map.update({i: i for i in sites if i not in seq_site_map})

    return seq_site_map

def assign_site_type_to_read(read_seqs, site_seq_map):
    """Takes a read sequence and identifies the position of all site types'.

    Args:
        read_sequence: The read sequence.
        sites: A list of sequences corresponding to each site type.

    Returns:
        A tuple containing the list of reads and a dictionary of read types.
    """

    # Sub-function to be used on each site:
    def get_unique_sites(read_seq, site_seq_map):
        site_maps = []
        # Iterate over keys in dictionary (site seqs)
        for key, value in site_seq_map.items():
            # Assign left-hand value of site location, if it exists.
            i_l = read_seq.find(key)

            # Check if index is -1
            if i_l != -1:
                i_r = i_l + len(key) - 1
                site = (value, i_l, i_r)
                # if "CCATTCCA" in read_seq:
                #     print(site_seq_map)
                # Check if site_maps is empty (suggesting first map)
                if site_maps == []:

                    site_maps.append(site)

                # Else only add if the sites in the list do not overlap.
                else:
                    add = True
                    site_maps_temp = [i for i in site_maps]
                    for site_map in site_maps:
                        j_l, j_r = site_map[1: ]

                        # Checks if new map is an extension of old map
                        if i_l <= j_l and i_r >= j_r:
                            
                            # Removes map
                            site_maps_temp.remove(site_map)

                        elif ((i_l >= j_l and i_r < j_r)
                              or (i_l > j_l and i_r <= j_r)):
                            add = False
                    if add:
                        site_maps_temp.append((value, i_l, i_r))
                    site_maps = site_maps_temp
        for i, site in enumerate(site_maps):
            if "Centered" in site[0]:
                site_maps[i] = ("Centered", site[1], site[2])
        return site_maps



    # Initialize count dictionary with keys as site names and a value of
    # for each item.
    count_site_map = {value: 0 for value in site_seq_map.values()}
    count_multisite_map = {value: 0 for value in site_seq_map.values()}
    count_site_map.update({"None": 0})
    count_multisite_map.update({"None": 0})
    # Pre-allocate output with "None" for each line.
    output = ["None"]*len(read_seqs)

    for i, read_seq in enumerate(read_seqs):
        # Remove trailing \n
        read_seq = read_seq.strip()[:37]

        site_maps = get_unique_sites(read_seq, site_seq_map)
        if site_maps:  
            output[i] = ", ".join(["%s:%s-%s" % (site, index, end)
                                   for site, index, end in site_maps])
            for site_map in site_maps:
                count_site_map[site_map[0]] += 1
            if len(site_maps) > 1:
                key = ",".join(sorted([j[0] for j in site_maps]))
                if key in count_multisite_map:
                    count_multisite_map[key] += 1
                else:
                    count_multisite_map[key] = 1
            else:
                count_multisite_map[site_maps[0][0]] += 1
        else:
            count_site_map["None"] += 1
            count_multisite_map["None"] += 1


    return output, count_site_map, count_multisite_map

