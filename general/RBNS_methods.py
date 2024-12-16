import imp # Used to import general.py
import time # Used for time
import os
import subprocess
import numpy as np
import pandas as pd
import math
import re
import sys
import itertools as it
import re
from subprocess import PIPE, Popen
pd.set_option('display.width', 1000)
imp.load_source("general",
                ("/lab/bartel1_ata/mcgeary/computation/"
                 "AgoRBNS/general/general.py"
                )
               )

imp.load_source("sitetypes",
                ("/lab/bartel1_ata/mcgeary/computation/"
                 "AgoRBNS/AssignSiteTypes/sitetypes.py"
                )
               )
from general import *
from sitetypes import get_seq_site_map

cyrano = "cgggtgtgcgccggaagagcagcgccgcagtagcagcacacgggaccgaaatggcgtaacgcgcagtcgagcaccgcagcagcgcagaggggccaggctcagcgaacaccacagagcctggcccgaacaacttcaaccagcaaatttgttgaaaattactcttcactccagtccacatgtttccaaaccagcacgagatcaccgctgacgtgcccctcatttgtcatatcaagaagagaaaagaacagtcctgaagggaacacacgaccagatgagaccatcttgtctttgttcttcaaaccgcccttcatcctgataagttgaccttcctccagcgttatttcagggaagtgtttgttcagcataccaagccaagtgttactggtctatgaacatatcatctctgcactttggagaagacaactgaagttttttttaattcacataaactcccccacaaactcttcttgtgttcaaaaaagctaaatgacaacacaattgtgccacaacatggtgccatcgactggtcatattgttcagaataaatcaaatatggttgcagtttttgttaaatgaagataaaaatgcagaaaaaaatgcacatttaaaatctaaatcatgcaagcaagaagacagagctacatccaacagttttattatctgcagattaagtgcctagacgagtctctttttaatagctcctgcctcaatttccatggcttcgctctaaagtaaagcacttctttggacattttagtctttcaaatggctgaaaatgtagtttctgagacgacaaactgtagacagctagctgttctcttaattcacaactcatattggtgctttagcaaaaccatgacccatgctttcatgttttctacagtatatacttgtactgacaaaccaagacaggcagtggcataacaacccaaaaaaaacataagctttataaatcccacagttaaacctttctagcggggtgctattgagttgcaggatgttatgaaaaatagttattttactacctcagaagttaccacagctcagagccacacttggatcaaactataagagcatgcttgtcttttttgtttgcagtaaaataaatatagtcgttatagtctgtctaatggcttggtaatctgatgagaagaggaattaatttttcttttgtatgtgtgaaatatcgtattgagactgtatgttctccagaataagataaacatcatgttcatccattaacgttatgcctgaagaaaggcagtaaaatggttaaatgatcataacatctttaaatagttgctaatagtctaatatgttggacagtaagactggtttctagttatctgctatagagcactgtgaaaataggaatttctcttcatctgttgtgacactacatcttttttgttgttgttggtatagtagttcctatcatacggatggtagatattcaagtgttctgagttttcatgagaggaagagtagcgaggcctcagtgtgggaaatgtgcaatatttgtttctttgtcttgttttgagatgaatatgtatattgtacaaacaagtgacaagttgttcgcaaaaacaacaaaatcaccaatgtcttccattaaatgatgtatagtttcatctgcttgctttagggtggcgggatttgtagtgccagcgtagacatggacatgcatgtacggtagcgatgcttgctgagtttaatagtcgtagtgacattcagtataggctataataaagtcagttagtacataaatataggcaaagatttagttagtaaacaccatagatgcagtatccacaccttggttaatttgattaagaagcgagtagtagggctgcacgattctggctaaaatgaaaatcacgtttttttgttgttgttttttttttaatcaagatcgcgattttctcacaattctatagatgtagaataaaggtttatatgagaaggctattattattattgttaattttcaaatatatggtaaaacaacattagcattaggcctgtgtcattctaacgttttgtataggcttggaatggtggactggttcatgtctttaacaacattattagaaaatttttttcactggtaaaattaaaccagtagtatcgtcattaatggtctatctacaatagataatgttacatgtggaaaattaactaaaatgctcacattcacatgtttgtgcaataatagaccttttgatttgaatagttttgtggagccctttgtgttgaattcatgccttttatggtaacctgttgcatgttttctgacgctttagttcgagtgcatatctgaatcagttcatttttgttcctgtcaatataatttagtagctataacaaccacaacgatctctttaaaagaaaataacaatagtttaaaaaatgagcaacaggaaacattaaaacaattttatatacatcaaatagccaaattttgctggaaaaatgctcttttggtaacatatttcattcggtataatttgtatctttattttaaagtatttttgttggtcagttttctacaaaataaacgaaaagaacgaatgcattttagctgaagtgatgtgcaaatacttcttcaataataagggaattatgaattaaattcaagtaaacctattataacataaaatatagtatttctaataattttatttatttaatagatttgagtttgaattacagtatcacaatactatttgattccacaatacttcaactggtatagtatcttaacgtaaattaatagtatcgcgacaaccctaataaatatagtatgggacacgtaacagcatttggaggtgtccacattgctgagcaacttgcctaaaacaatgtgaagacttgtgcgactgcccttacgcttctctcctgacttaatatacgtaggtgaatcgtacaacagctggattaagatcatgtgtagggtcgaatcaagatcacgttcttttatcgattaatcgtgcagccctagcaagtaagcatttactcttggttttagtttttaagtgacatacgctgtagtttattcactggtaatcactattagttgatgataacgtcatagcatgctgaggattagtaaaatggcggttggtattgatgtgcgactagccatgctgaaatgcatatatatgtttcataagcactgtgtcagtttgtcctagtttgtactgtgcaaaatgaagttctgggctcaaaatagtcgactcagatcaaaacgtctgtcccaaataggccaaaatttggctcagtagcttagaatacgcaggtgttgaagaacttcattcagtacctgtttcaaggttgctaaaactgtactggacagaaagaatagacttctctcttcctaaatagaactgtgtttgttttcatacctctaatattgcattctccatcaggtaatgttgtctcaaactagagctgtgggtagattattcctagtagcagaacatgtagattttaactcaacgtgtagaaaattaagtggcataaagcaattatcaaagcattgtaatgtagataaatacctcaatgactggaatgcaaactatatatttttgaatatactcgaaaaaatatgcaagccatgatgtgtattaagttcatcttttagaagcaataatcagtttatgcatgagttaaagtaccggaacagtttagtggtgtttgtaatatctgtagcatatctgtgacagtgacaagtgtttgatcctgaatatctcattgctcataattaaaagttttccattcctccatttttcatggacaccaaataaatatgttaaggtagcactgggggtctgactaataaagaggattttttgggtgggttgacgaaagttagaagaatcatctacatatgcagggcttttgctcttattaaatttgcatgtttacctaatcattgtggaaagcacataagctgtatgtttattgacagcttattatcattaacagcaatatgagcagtgttctgggtagcaatttgtaatactgatgactcattatcctcatgaagcctggacaaacattttaccactttacatttacaaaccccatttggagtaaagtattttttaaagagtgcaaaaaaaacctgatgcatggtagaaatggataagaaatgttttgtgctgcatgcaataattgagctttgcttgaaaatgtttatgcagaaatgtgatgacactaaatctatgagtggaccagtgttgtgactttctttttaaatgcatggaatcccaaagttgatctctactgtaaattattcacagacgtgttcattgagaaatgtttgattactgtatgatggaattgtaacaagtttgttcttcaagaggcaaagtcagatgcagtgttttagcaaagtcttaacaattgtgagttgaaacgggaaaaactatttgtgttatgttcaatgcaaataatattgcagaatataagttaatatgctttgcagacctgatttgcgttgctgattttgcctctgaatgtgaaactctatgcttagtctgatctatgcactacagaagcgtcccgcatggttttgtacaaaatattaaattttagcaataaaaaaaagcagcatggtgccaatcat"

RC_MAP = {
    "C" : "G",
    "G" : "C",
    "T" : "A",
    "U" : "A",
    "N" : "N",
    "A" : "T", 
}

WOB_MAP = {
    "U" : "C",
    "T" : "C",
    "G" : "A",
}

def order_site_map(mirna, sitelist):
    path = "AssignSiteTypes/sites.%s_%s.txt" %(mirna, sitelist)
    with open(path, "rb") as f:
        list_ = {site.strip(): i for i, site in enumerate(f.readlines())}
    return(list_)

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)  # allows duplicate elements
    return [[i] for i in s] + list(it.chain.from_iterable(it.combinations(s, r)
                                        for r in range(2, len(s)+1)))

def get_all_wobbles(kmer):
    w_map = {"G" : "A", "T" : "C"}
    w_tuples = [(i, w_map[nt]) for i, nt in enumerate(kmer) if nt in ["G", "T"]]
    all_w_tuples = list(powerset(w_tuples))
    w_kmers = []
    w_pos = []
    for w_tuple in all_w_tuples:
        kmer_temp = list(kmer)
        for w in w_tuple:
            kmer_temp[w[0]] = w[1]
        w_kmers.append("".join(kmer_temp))
        w_pos.append([w[0] for w in w_tuple])
    return(list(zip(w_kmers, w_pos)))

def get_all_mismatches(kmer, internal_only=True):
    nts = ["A", "C", "G", "T"]
    # Assign the length of the output and the range of the loop based on whether
    # or not the ends are being considered with respect to the mismatches.
    if internal_only:
        p_l, p_r = (1, len(kmer) - 1)
        len_out = len(kmer) - 2
        i_shift = -1
    else:
        p_l, p_r = (0, len(kmer))
        len_out = len(kmer)
        i_shift = 0
    # Pre-allocate the output lists.
    mm_kmers = [None]*len_out
    mm_pos = [None]*len_out
    mm_nts = [None]*len_out
    mm_mmnts = [None]*len_out
    for i in range(p_l, p_r):
        nt = kmer[i]
        mms = [j for j in nts if j != nt]
        mismatches = [kmer[:i] + mm + kmer[i+1:] for mm in mms]
        mm_kmers[i + i_shift] = [kmer[:i] + mm + kmer[i+1:] for mm in mms]
        mm_pos[i + i_shift] = [i for mm in mms]
        mm_nts[i + i_shift] = [nt for mm in mms]
        mm_mmnts[i + i_shift] = [mm for mm in mms]
    mm_kmers = [j for i in mm_kmers for j in i]
    mm_pos = [j for i in mm_pos for j in i]
    mm_nts = [j for i in mm_nts for j in i]
    mm_mmnts = [j for i in mm_mmnts for j in i]
    return(list(zip(mm_kmers, mm_pos, mm_nts, mm_mmnts)))


def get_all_bulges(kmer):
    nts = ["A", "C", "G", "T"]
    b_kmers = [None]*(len(kmer) - 1)
    b_pos = [None]*(len(kmer) - 1)
    b_nts = [None]*(len(kmer) - 1)
    for i in range(1, len(kmer)):
        b_kmers[i - 1] = [kmer[:i] + nt + kmer[i:] for nt in nts]
        b_pos[i - 1] = [i for nt in nts]
        b_nts[i - 1] = ["b%s" %nt for nt in nts]
    b_kmers = [j for i in b_kmers for j in i]
    b_pos = [j for i in b_pos for j in i]
    b_nts = [j for i in b_nts for j in i]
    return(list(zip(b_kmers, b_pos, b_nts)))

def get_all_deletions(kmer):
    d_kmers = [None]*(len(kmer) - 2)
    d_pos = [None]*(len(kmer) - 2)
    for i in range(1, len(kmer) - 1):
        d_kmers[i - 1] = kmer[:i] + kmer[i+1:]
        d_pos[i - 1] = i
    print(d_kmers)
    print(d_pos)
    sys.stdout.flush()
    return(list(zip(d_kmers, d_pos)))




################################################################################
class Sequence:
    def __init__(self, seq):
        self.seq = seq
    def __repr__(self):
        return(self.seq)
    def __len__(self):
        return(len(self.seq))
    def rev(self, start=0, stop=None):
        return("".join([i for i in self.seq[start:stop]][::-1]))
    def rev_com(self, start=0, stop=None):
        return("".join([RC_MAP[i] for i in self.seq[start:stop]][::-1]))
    def rna(self):
        return(re.sub("T", "U", self.seq))
    def best_match(self, seq2):
        return(LCSubString(self.seq, seq2))
################################################################################
def LowerCaseKmerString(subkmer, kmer, w_kmer=False):
    p = kmer.find(subkmer)
    string = "".join([kmer.lower()[:p], kmer[p:p + len(subkmer)],
                      kmer.lower()[p + len(subkmer):]])
    if w_kmer:
        string = re.sub(subkmer, w_kmer, string)
    return(string)


class Mirna(Sequence):
    SEQ_NAME_MAP = {
        "miR-1"       : "UGGAAUGUAAAGAAGUAUGUAU",
        "miR-155"     : "UUAAUGCUAAUCGUGAUAGGGGU",
        "miR-124"     : "UAAGGCACGCGGUGAAUGCCAA",
        "let-7a-21nt" : "UGAGGUAGUAGGUUGUAUAGU",
        "let-7a"      : "UGAGGUAGUAGGUUGUAUAGUU",
        "lsy-6"       : "UUUUGUAUGAGACGCAUUUCGA",
        "miR-7"       : "UGGAAGACUAGUGAUUUUGUUGU",
        "miR-7-22nt"  : "UGGAAGACUAGUGAUUUUGUUG",
        "miR-7-23nt"  : "UGGAAGACUAGUGAUUUUGUUGU",
        "miR-7-24nt"  : "UGGAAGACUAGUGAUUUUGUUGUU",
        "miR-7-25nt"  : "UGGAAGACUAGUGAUUUUGUUGUUU",
        "miR-671"     : "AGGAAGCCCUGGAGGGGCUGG",
        "miR-122"     : "UGGAGUGUGACAAUGGUGUUUG",
        # Added for HeLa Transfection, from miRBASE #
        "miR-137"  : "UUAUUGCUUAAGAAUACGCGUAG",
        "miR-139"  : "UCUACAGUGCACGUGUCUCCAGU",
        "miR-143"  : "UGAGAUGAAGCACUGUAGCUC",
        "miR-144"  : "UACAGUAUAGAUGAUGUACU",
        "miR-153"  : "UUGCAUAGUCACAAAAGUGAUC",
        "miR-182"  : "UUUGGCAAUGGUAGAACUCACACU",
        "miR-199a" : "CCCAGUGUUCAGACUACCUGUUC",
        "miR-204"  : "UUCCCUUUGUCAUCCUAUGCCU",
        "miR-205"  : "UCCUUCAUUCCACCGGAGUCUG",
        "miR-216b" : "AAAUCUCUGCAGGCAAAUGUGA",
        "miR-223"  : "UGUCAGUUUGUCAAAUACCCCA",
        # Added for Virkam HCT studies, from miRBASE #
        "let-7c"   : "UGAGGUAGUAGGUUGUAUGGUU",
        "miR-16"   : "UAGCAGCACGUAAAUAUUGGCG",
        "miR-103"  : "AGCAGCAUUGUACAGGGCUAUGA",
        "miR-106b" : "UAAAGUGCUGACAGUGCAGAU",
        "miR-200a" : "UAACACUGUCUGGUAACGAUGU",
        "miR-200b" : "UAAUACUGCCUGGUAAUGAUGA",
        "miR-215"  : "AUGACCUAUGAAUUGACAGAC",
        "miR-21"   : "UAGCUUAUCAGACUGAUGUUGA",
        # Synthetic miRNAs used for three-prime paper with N.B. #
        #                   12345678
        "let-7a_plus1"   : "UGAGGUAGUAUGGUUGUAUAG",
        "let-7a_minus1"  : "UGAGGUAGUGGUUGUAUAGUA",
        "let-7a_miR-155" : "UGAGGUAGAAUCGUGAUAGGGGU",
        "miR-155_let-7a" : "UUAAUGCUUAGGUUGUAUAGU",
        # Star sequences for the re-analysis of the threeprime score #
        #                   12345678
        "miR-137*"        : "ACGCGUAUUCUUAAGCAAUAAAU",
        "miR-205*"        : "GACUCCGGUGGAAUGAAGCAAU",
        "miR-155*"        : "CCCUAUCACGAUUAGCAUUAAAU",
        "miR-223*"        : "GGGUAUUUGACAAACUGAUAAU",
        "miR-144*"        : "UACAUCAUCUAUACUCUAAU",
        "miR-143*"        : "GCUACAGUGCUUCAUCUUAAU",
        "miR-153*"        : "UCACUUUUGUGACUAUGUAAAU",
        "miR-216b*"       : "ACAUUUGCCUGCAGAGAUUUAU",
        "miR-199a8"       : "ACAGGUAGUCUGAACACUGCGAU",
        "miR-204*"        : "GCAUAGGAUGACAAAGGCAAAU",
        "miR-139*"        : "UGGAGACACGUGCACUGUACAAU",
        "miR-182*"        : "UGUGAGUUCUACCAUUGCUAAAAU",
        "miR-7*"          : "AACAAAAUCACUAGUCUUCUAAU",
        "miR-1*"          : "ACAUACUUCUUUACAUUCUAAU",
        "miR-124*"        : "GGCAUUCACCGCGUGCUUUAAU",
        "lsy-6*"          : "GAAAUGCGUCUCAUACAAAAAU",
        "let-7a*"         : "CUAUACAACCUACUACCUUAAU"        
    }
    CANONICAL_SITE_MAP = {"8mer-m1.8" : "8mer",
                          "7mer-m2.8" : "7mer-m8",
                          "7mer-m1.7" : "7mer-A1",
                          "6mer-m3.8" : "6mer-m8",
                          "6mer-m2.7" : "6mer",
                          "6mer-m1.6" : "6mer-A1",
                          "5mer-m4.8" : "5mer-m8",
                          "5mer-m1.5" : "5mer-A1"}
    SEED_SITES = ["8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-A1", "6mer-m8"]
    RANGE_SITE_MAP = {
        "8mer"    : (0, 8),
        "7mer-m8" : (1, 8),
        "7mer-A1" : (0, 7),
        "6mer"    : (1, 7),
        "6mer-m8" : (2, 8),
        "6mer-A1" : (0, 6),
        "5mer-A1" : (0, 5),
        "5mer-m8" : (3, 8),
    }
    REGEX_CDNST_MAP = {
        "CDNST.1" : (1, 6,        []),
        "CDNST.2" : (0, 7,       [4]),
        "CDNST.3" : (1, 9, [1, 3, 5]),
        "CDNST.4" : (0, 8, [4, 5, 6])
    }
    MODS = ["d", "w", "b", "mm"]
    def __init__(self, name):
        self.name = name
        self.seq = self.SEQ_NAME_MAP[name]
    def __getitem__(self, site):
        if "CDNST" in site:
            l, r, mm = self.REGEX_CDNST_MAP[site]
            m_list = [RC_MAP[i] for i in list(self.seq)]
            for pos_m in mm:
                m_list[pos_m] = "[^%s]" %(m_list[pos_m])
            out = m_list[l:r][::-1]
            return("".join(out))
        # Checks for non-complementary sites:
        s_split = site.split("mer")
        if site == "7mer-m3.9":
            l, r, mm = (1, 9, [1])
            m_list = [RC_MAP[i] for i in list(self.seq)]
            for pos_m in mm:
                m_list[pos_m] = "[^%s]" %(m_list[pos_m])
            out = m_list[l:r][::-1]
            out2 = [i for i in out]
            out2[-1] = "$"
            out1str = "".join(out)
            out2str = "".join(out2)
            out = "(?:%s|%s)" %(out1str, out2str)
            return(out)
        if site == "11mer-m3.13":
            l, r, mm = (1, 13, [1])
            m_list = [RC_MAP[i] for i in list(self.seq)]
            for pos_m in mm:
                m_list[pos_m] = "[^%s]" %(m_list[pos_m])
            out = m_list[l:r][::-1]
            out2 = [i for i in out]
            out2[-1] = "$"
            out1str = "".join(out)
            out2str = "".join(out2)
            out = "(?:%s|%s)" %(out1str, out2str)
            return(out)

        if "mer" not in site:   # Atypical sites     
            return(site)
        else:
            s_suffix = s_split[1]
            # Get the list of base modifications:
            #           (mm or w) + (A or C or T or G) <- 1 or 0  + (num or (num.num))
            f_l    = re.findall(r'^(?:A|C|G|T)+-', site)
            f_r    = re.findall(r'-(?:A|C|G|T)+$', site)
            if f_l:
                site = re.sub(f_l[0], "", site)
            if f_r:
                site = re.sub(f_r[0], "", site)
            mm_w = re.findall(r'(?:mm|w)(?:A|C|G|T)?(?:\(\d+\.\d+\)|\d+)', site)
            #           (b or d) + (A or C or T or G) <- 1 or 0  + (num or (num.num))
            b_d = re.findall(r'(?:b|d)(?:A|C|G|T)?(?:\(\d+\.\d+\)|\d+)', site)
            if mm_w + b_d:
                site = re.sub("-$", "", site.split((mm_w + b_d)[0])[0])
            # Define left and right indeces of the target site:
            if site in self.RANGE_SITE_MAP:
                (l, r) = self.RANGE_SITE_MAP[site]
            else:   # Converts m3.9 to  [2, 9]
                (l, r) = (int(i) for i in site.split("-m")[1].split("."))
                l -= 1
            # Split the mirna into a list of its positions:
            m_list = list(self.seq)
            for i in mm_w:
                if i[0] == "w":
                    if i[1] == "(":
                        pos = int(i.split(".")[0][2:]) - 1
                    else:
                        pos = int(i[1:]) - 1
                    nt = WOB_MAP[m_list[pos]]
                else:
                    if i[3] == "(":
                        pos = int(i.split(".")[0][4:]) - 1
                    else:
                        pos = int(i[3: ]) - 1
                    nt = RC_MAP[i[2]]
                m_list[pos] = nt
            for i in b_d:
                if i[0] == "b":
                    if i[2] == "(":
                        pos = int(i.split(".")[0][3:]) - 1
                    else:
                        pos = int(i[2:]) - 1
                    nt = RC_MAP[i[1]] + m_list[pos]
                else:
                    if i[1] == "(":
                        pos = int(i.split(".")[0][2:]) - 1
                    else:
                        pos = int(i[1:]) - 1
                    nt = ""
                m_list[pos] = nt
            base_seq = Sequence("".join(m_list[l:r])).rev_com()
            if f_l:
                base_seq = f_l[0][:-1] + base_seq
            if f_r:
                base_seq = base_seq + f_r[0][1:]
            return(base_seq)
    def get_kmer_match(self, kmer, print_site=False):
        kmerlist = list(kmer)
        mir_rc_seq = self.rev_com()
        if kmer in mir_rc_seq:
            k_l =  mir_rc_seq.find(kmer)
            k_r = k_l + len(kmer)
            m_l, m_r = [len(self) - i for i in [k_r, k_l]]
            if m_l <= 4 and m_r > 8:
                kmerlist[:m_r - 8] = [i.lower() for i in kmerlist[:m_r - 8]]
                m_r = 8
            kmer = "".join(kmerlist)
            site_name = "%smer-m%s.%s" %(m_r - m_l, m_l +1, m_r)
            if site_name in self.CANONICAL_SITE_MAP:
                site_name = self.CANONICAL_SITE_MAP[site_name]
            return([[site_name, kmer, k_l, list(range(m_l, m_r))]])
        else:
            return([])
    def get_kmer_bulge_match(self, kmer, print_site=True):
        mir_rc_seq = self.rev_com()
        out = []
        for i_b in range(1, len(kmer) - 2):
            l_side = kmer[:i_b]
            r_side = kmer[i_b + 1:]
            b = kmer[i_b]
            b_kmer = l_side + r_side
            match = self.get_kmer_match(b_kmer, print_site=False)
            if match:
                site_name, kmer_match, k_l, mir_range = match[0]
                m_r = max(mir_range)
                shift_kmer = len([c for c in kmer if c.islower()])
                if site_name[-1] == "r":
                    spacer = "-"
                else:
                    spacer = ""
                m_b = len(self) - k_l - i_b + 1
                if m_r >= m_b - 1:
                    site_name = "%s%sb%s%s" %(site_name, spacer, b, m_b)
                    kmer_new = kmer_match[:i_b]+"b"+ b.lower() + kmer_match[i_b:]
                    out.append([site_name, kmer_new, k_l, mir_range])
        base_sites = [i[0].split("b")[0] for i in out]
        if len(out) > 1 and len(set(base_sites)) == 1:
            pos = [int(i[0].split("b")[1][1:]) for i in out][::-1]
            nt = out[0][0].split("b")[1][0]
            site_name = "%sb%s(%s.%s)" %(base_sites[0], nt, min(pos), max(pos))
            out[0][0] = site_name
            out = [out[0]]
        return(out)
    def get_kmer_wobbles(self, kmer, print_site=False):
        mir_rc_seq = self.rev_com()
        out = []
        wobbles_all = get_all_wobbles(kmer)
        wob_dict = {"A": "H", "C" : "V"}
        for wob_tup, i_w in wobbles_all:
            w = [kmer[i] for i in i_w]
            matches_list = self.get_kmer_match(wob_tup, print_site=False)
            for match_i in matches_list:
                site_name, kmer_match, k_l, mir_range = match_i
                m_r = max(mir_range)
                kmer_list = list(kmer_match)
                lower_capture = r"^((?:a|c|g|t)*)((?:A|C|G|T|V|H)+)((?:a|c|g|t)*)((?:A|C|G|T|V|H)+)((?:a|c|g|t)*)$"
                shift_kmer = sum([len(re.sub(lower_capture, i, kmer_match))
                                  for i in ["\\1", "\\5"]])
                mir_w = [m_r - i + shift_kmer for i in i_w][::-1]
                if max(mir_w) <= m_r and 0 not in mir_w:
                    for i in zip(i_w, mir_w):
                        kmer_list[i[0]] = wob_dict[kmer_list[i[0]]]
                        mir_range.remove(i[1])
                    if len(mir_range) > 0:
                        kmer_match = "".join(kmer_list)
                        if site_name[-1] == "r":
                            spacer = "-"
                        else:
                            spacer = ""
                        string1 = "w" + ",".join([str(i + 1) for i in mir_w])
                        site_name = site_name + spacer + string1
                        out.append([site_name, kmer_match, k_l, mir_range])
            if len(out) > 0:
                break
        return(out)
    def get_kmer_1_mismatch(self, kmer, print_site=False):
        mir_rc_seq = self.rev_com()
        out = []
        for i_mm in range(0, len(kmer) - 1):
            l_side = kmer[:i_mm]
            r_side = kmer[i_mm + 1:]
            mm = kmer[i_mm]
            mm_kmers = [l_side + nt + r_side for nt in ["A", "C", "G", "T"]
                        if nt != mm]
            matches_list = [self.get_kmer_match(mm_kmer, print_site=False)
                            for mm_kmer in mm_kmers]
            matches_list = [i for i in matches_list if i]
            for match_i in matches_list:
                site_name, kmer_match, k_l, mir_range = match_i[0]
                m_r = max(mir_range)
                kmer_match = kmer_match[:i_mm] + mm.lower() + kmer_match[i_mm + 1:]
                lower_capture = r"^((?:a|c|g|t)*)((?:A|C|G|T)+)((?:a|c|g|t)*)((?:A|C|G|T)+)((?:a|c|g|t)*)$"
                shift_kmer = sum([len(re.sub(lower_capture, i, kmer_match))
                                  for i in ["\\1", "\\5"]])
                if site_name[-1] == "r":
                    spacer = "-"
                else:
                    spacer = ""
                mir_mm = m_r - i_mm + shift_kmer
                if ((mir_rc_seq[k_l + i_mm] == "A" and kmer[i_mm] == "G")
                    or (mir_rc_seq[k_l + i_mm] == "C" and kmer[i_mm] == "T")):
                    string1 = "w"
                    wob_dict = {"G": "H", "T" : "V"}

                    kmer_match = kmer_match[:i_mm] + wob_dict[kmer[i_mm]] + kmer_match[i_mm + 1:]
                else:
                    string1 = "mm%s" %(mm)
                if mir_mm < m_r:
                    mir_range.remove(mir_mm)
                    site_name = site_name + spacer + string1 + str(mir_mm + 1)
                    out.append([site_name, kmer_match, k_l, mir_range])
        return(out)
    def get_kmer_bulge_1_mismatch(self, kmer, print_site=True):
        mir_rc_seq = self.rev_com()
        out = []
        for i_b in range(1, len(kmer) - 2):
            l_side = kmer[:i_b]
            r_side = kmer[i_b + 1:]
            b = kmer[i_b]
            b_kmer = l_side + r_side
            match = self.get_kmer_1_mismatch(b_kmer, print_site=True)
            if match:
                site_name, kmer_match, k_l, mir_range = match[0]
                m_r = max(mir_range)
                shift_kmer = len([c for c in kmer if c.islower()])
                lower_capture = r"^((?:a|c|g|t)*)((?:A|C|G|T)+)((?:a|c|g|t)*)((?:A|C|G|T)+)((?:a|c|g|t)*)$"
                shift_kmer = sum([len(re.sub(lower_capture, i, kmer_match))
                                  for i in ["\\1", "\\5"]])
                if site_name[-1] == "r":
                    spacer = "-"
                else:
                    spacer = ""
                m_b = len(self) - k_l - i_b + 1
                if m_r >= m_b:
                    site_name = "%s%sb%s%s" %(site_name, spacer, b, m_b)
                    kmer_new = kmer_match[:i_b]+"b"+ b.lower() + kmer_match[i_b:]
                    out.append([site_name, kmer_new, k_l, mir_range])
        return(out)
    def get_kmer_bulge_wobbles(self, kmer, print_site=True):
        mir_rc_seq = self.rev_com()
        out = []
        for i_b in range(1, len(kmer) - 2):
            l_side = kmer[:i_b]
            r_side = kmer[i_b + 1:]
            b = kmer[i_b]
            b_kmer = l_side + r_side
            match = self.get_kmer_wobbles(b_kmer, print_site=True)
            if match:
                site_name, kmer_match, k_l, mir_range = match[0]
                m_r = max(mir_range)
                shift_kmer = len([c for c in kmer if c.islower()])
                lower_capture = r"^((?:a|c|g|t)*)((?:A|C|G|T)+)((?:a|c|g|t)*)((?:A|C|G|T)+)((?:a|c|g|t)*)$"
                shift_kmer = sum([len(re.sub(lower_capture, i, kmer_match))
                                  for i in ["\\1", "\\5"]])
                if site_name[-1] == "r":
                    spacer = "-"
                else:
                    spacer = ""
                m_b = len(self) - k_l - i_b + 1
                if m_r >= m_b:
                    site_name = "%s%sb%s%s" %(site_name, spacer, b, m_b)
                    kmer_new = kmer_match[:i_b]+"b"+ b.lower() + kmer_match[i_b:]
                    out.append([site_name, kmer_new, k_l, mir_range])
        return(out)
    def get_best_name(self, kmer_sort, suboptimal=False, print_sites=True,
                      output_for_table=False):
        kmer, R, pos_mean, pos_sd, pos_mean_dif, pos_sd_dif = kmer_sort
        starting_pos = len(kmer)
        self.seq = "U" + self.seq[1:]
        mirna_rc_seq = self.rev_com()
        k_len = len(kmer)
        leave = False
        l_mir = len(self)
        match_result = []
        b_i = False
        leave_in_two = False
        leave_in_one = False
        while not leave:
        # Check each independent subkmer:
            subkmers = [kmer[i:i+k_len] for
                        i in range(0, len(kmer) - k_len + 1)]
            s_l = [kmer[:i].lower() for i in range(0, len(kmer) - k_len + 1)]
            s_r = [kmer[i+k_len:].lower() for i in range(0, len(kmer) - k_len + 1)]
            if leave_in_one:
                leave = True
            if leave_in_two:
                leave_in_one = True        
            for i_k, subkmer in enumerate(subkmers):
                m = self.get_kmer_match(subkmer)
                if m:
                    match_result.append(m[0] + [s_l[i_k], s_r[i_k]])
                    leave = True
                b = self.get_kmer_bulge_match(subkmer)
                w = self.get_kmer_wobbles(subkmer)
                mm = self.get_kmer_1_mismatch(subkmer, print_site=False)
                b_w = self.get_kmer_bulge_wobbles(subkmer)
                b_mm = self.get_kmer_bulge_1_mismatch(subkmer)
                # if b_w or b_mm or b or mm or w:
                #     leave_in_two = True
                for m_i in b + mm + w + b_mm + b_w:
                    if m_i[0] not in [i[0] for i in match_result]:
                        match_result.append(m_i + [s_l[i_k], s_r[i_k]])
            k_len -= 1
        mir_space = 30
        sites = [i[0] for i in match_result]
        strings = [i[4] + i[1] + i[5] for i in match_result]
        ranges = [i[3] for i in match_result]
        
        m1A_weight = 0.4
        m2_8_weight = 1
        m9_11_weight = 0.2
        m12_16_weight = 0.7
        m17_end_weight = 0.5
        mir_weights = ([m1A_weight] + [m2_8_weight]*7 + [m9_11_weight]*3 +
                       [m12_16_weight]*5 + [m17_end_weight]*(len(self) - 16))
        
        positional_scores = [[mir_weights[j - 1] for j in i] for i in ranges]
        dinucs = [[[j + ranges[n][0], i[j:j+2]] for j in range(len(i) - 2 + 1)
                               if i[j].isupper() and i[j + 1].isupper()]
                              for n, i in enumerate(strings)]
        dinuc_lens = [len(i)*0.5 for i in dinucs]
        # Compute the raw score based on the positions of pairing and weights:
        raw_scores = [sum([mir_weights[j - 1] for j in i]) for i in ranges]
        # Adds the contiguousness score to the raw score, positionl inspecific
        raw_scores = [sum(i) for i in zip(raw_scores, dinuc_lens)]
        strings_wobble_removed = [re.sub("b.", "", i) for i in strings]
        lower_capture = r"^((?:a|c|g|t)*)((?:A|C|G|T|H|V)+)((?:a|c|g|t)*)((?:A|C|G|T|H|V)+)((?:a|c|g|t)*)$"
        
        # for site, rang, string in zip(sites, ranges, strings):
        #     print(site + " " + string)
        #     print(sum([i in site for i in ["b2", "b3", "b4", "b5"]]))
        #     print(rang)
        #     print(sum([i in rang for i in [2, 3, 4, 5]]))
        m2_5match = [5*(sum([i in site for i in ["b2", "b3", "b4", "b5"]]) == 0
                        and sum([i in rang for i in [1, 2, 3, 4]]) == 4)
                     for site, rang in zip(sites, ranges)]
        raw_scores = [sum(i) for i in zip(raw_scores, m2_5match)]
        # Really awful way to get the number of wobbles, bulges, mismatches,
        # left-hand flanking g/a/t nucleotides,
        #and right-hand flanking g/a/t nucleotides.
        len_w = [len([nt for nt in string if nt in ["H", "V"]])
                   for string in strings]
        len_b = [len([nt for nt in string if nt == "b"])
                   for string in strings]
        len_mm = [len(re.sub(lower_capture, "\\3", re.sub("b.", "", i)))
                  for i in strings_wobble_removed]        
        len_g_l = [len([j for j in re.sub(lower_capture, "\\1", re.sub("b.", "", i))
                        if j == "g"])
                                   for i in strings]
        len_g_r = [len([j for j in re.sub(lower_capture, "\\5", re.sub("b.", "", i))
                        if j == "g"])
                                   for i in strings]
        len_a_l = [len([j for j in re.sub(lower_capture, "\\1", re.sub("b.", "", i))
                        if j == "a"])
                                   for i in strings]
        len_a_r = [len([j for j in re.sub(lower_capture, "\\5", re.sub("b.", "", i))
                        if j == "a"])
                                   for i in strings]
        len_t_l = [len([j for j in re.sub(lower_capture, "\\1", re.sub("b.", "", i))
                        if j == "a"])
                                   for i in strings]
        len_t_r = [len([j for j in re.sub(lower_capture, "\\5", re.sub("b.", "", i))
                        if j == "a"])
                                   for i in strings]
        
        kmer_data = list(zip(strings, raw_scores, len_w, len_b, len_mm, len_g_l, len_g_r,
                   len_a_l, len_a_r, len_t_l, len_t_r))
        w_weight = 0.5
        b_weight = -0.4
        mm_weight = -0.5
        flank_g_weight = -0.5
        flank_at_weight = 0.4
        final_scores = [i[1] + i[2]*w_weight + i[3]*b_weight + i[4]*mm_weight +
                  sum(i[5:7])*flank_g_weight + sum(i[7:])*flank_at_weight
                  for i in kmer_data]
        first = True
        out = []
        for i, match_i in enumerate(match_result):
            if final_scores[i] == max(final_scores):
                site, string, k_l, m_r, s_l, s_r = match_i
                out.append(site)
                bulge = False
                if print_sites :
                    if "b" in string:
                        b_ind = string.find("b")
                        group = r'(?:mm|w)(?:A|C|G|T)?(?:\(\d+\.\d+\)|\d+)'
                        group1 = r'b((?:a|c|g|t))'
                        group2 = r'b(?:a|c|g|t)'
                        b_string = re.search(group1, string).group(1)
                        string = re.sub(group2, "", string)
                        bulge = True
                    if output_for_table:
                        if first:
                            print("  %s" %(kmer) + "  %6.2f" %(R) +
                                  "  %4.1f" %(pos_mean) + "   %4.1f" %(pos_sd))
                    else:
                        if first:

                            print((" "*(10 + k_l - len(s_l)) + s_l + string + s_r +
                                  " "*(mir_space - k_l - len(string) - len(s_r)) +
                                  site +" "*(20 - len(site)) +
                                  "%4.1f" %(final_scores[i]) + "  %s" %(kmer) +
                                  "  %6.2f" %(R) + "  %4.1f" %(pos_mean) +
                                  "   %4.1f" %(pos_sd) + "    %4.1f" %(pos_mean_dif) +
                                  "    %4.1f" %(pos_sd_dif)))
                            first = False
                        else:
                            print((" "*(10 + k_l - len(s_l)) + s_l + string + s_r +
                                  " "*(mir_space - k_l - len(string) - len(s_r)) +
                                  site + " "*(20 - len(site)) +
                                  str(final_scores[i])))
                        if bulge:
                            print((" "*(10 + k_l + b_ind - 1) + b_string + "/"))
        if suboptimal and print_sites and not output_for_table:
            print("Suboptimal sites:")
            for i, match_i in enumerate(match_result):
                if final_scores[i] != max(final_scores):
                    site, string, k_l, m_r, s_l, s_r = match_i
                    bulge = False
                    if "b" in string:
                        b_ind = string.find("b")
                        group = r'(?:mm|w)(?:A|C|G|T)?(?:\(\d+\.\d+\)|\d+)'
                        group1 = r'b((?:a|c|g|t))'
                        group2 = r'b(?:a|c|g|t)'
                        b_string = re.search(group1, string).group(1)
                        string = re.sub(group2, "", string)
                        bulge = True
                    print((" "*(k_l - len(s_l)) + s_l + string + s_r + " "*(mir_space - k_l - len(string) - len(s_r))+ site + " "*(20 - len(site)) + str(final_scores[i])))
                    if bulge:
                        print((" "*(k_l + b_ind - 1) + b_string + "/"))
        return(out)
    def get_all_8mer_mismatch_names(self):
        seq_8mer = self["8mer"]
        seq_8mer_list = list(seq_8mer)
        # Make the list of each of the sequences and names that are all 16 mismatch
        # possibilities.
        seqs_8merMm = [[seq_8mer[:i] + j + seq_8mer[i+1:]
                        for j in ["A", "C", "G", "T"] if j != seq_8mer[i]]
                       for i in range(1, 7)]
        seqs_8merMm = [j for i in seqs_8merMm for j in i]
        names_8merMm = [["8mer-mm%s%s" %(j, 8 - i) for j in ["A", "C", "G", "T"]
                         if j != seq_8mer[i]]
                        for i in range(1, 7)]
        names_8merMm = [j for i in names_8merMm for j in i]
        return(names_8merMm)

################################################################################
LOWERCASE_MAP = {"A": "a", "C" : "c", "G" : "g", "T" : "t"}
def lowercase(matchobj):
    return LOWERCASE_MAP[matchobj.group(0)[2]]


class Site(Sequence):
    def __init__(self, _mirna, name):
        self.name = name
        self.seq = _mirna[name]
    def __str__(self):
        return(re.sub("\[\^.\]", lowercase, self.seq))
    def __len__(self):
        if self.name == "7mer-m3.9":
            return(7)
        if self.name == "CDNST.1":
            return(5)
        if self.name == "CDNST.2":
            return(7)
        if self.name == "CDNST.3":
            return(7)
        if self.name == "CDNST.4":
            return(8)
        if self.name == "11mer-m3.13":
            return(11)
        else:
            return(len(str(self)))
################################################################################
class SiteList:
    BASE_DIR = HOME_DIR + "AssignSiteTypes/site_categories"
    FLANKS = ["%s.%s" %(i[:2], i[2:]) for i in get_kmer_list(4)]
    def __init__(self, _mirna, list_name, read_len):
        mirname = _mirna.name
        print(_mirna.name)
        if list_name == "programmed":
            can_names = ["8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-m8",
                         "6mer-A1"]
            mm_names = _mirna.get_all_8mer_mismatch_names()
            site_names = can_names + mm_names
        elif list_name == "papercutoff3":
            f_name1 = "%s/%s/sites.%s_%s.txt" %(self.BASE_DIR, "papercutoff",
                                     "-".join(mirname.split("-")[:2]),
                                     "papercutoff")
            f_name2 = "%s/%s/sites.%s_%s.txt" %(self.BASE_DIR, "papercutoff2",
                                     "-".join(mirname.split("-")[:2]),
                                     "papercutoff2")
            with open(f_name1) as f:
                site_names1 = f.read().split("\n")
            with open(f_name2) as f:
                site_names2 = f.read().split("\n")
            site_names = list(set(site_names1).intersection(set(site_names2)))

        else:
            if mirname in ["let-7a_miR-155", "miR-155_let-7a"]:
                mirname_base = mirname
            else:
                mirname_base = "-".join(mirname.split("-")[:2])
            f_name = "%s/%s/sites.%s_%s.txt" %(self.BASE_DIR, list_name,
                                     mirname_base, list_name)
            with open(f_name) as f:
                site_names = f.read().split("\n")
        self.name = list_name
        self.mirna = _mirna.name
        self.sites = [Site(_mirna, s) for s in site_names]


        with open('general/site_info.txt',"r+") as site_pos_file:
            site_pos_map = dict()
            for line in site_pos_file:
                line_strip = line.strip().split("\t")
                key = line_strip[0]
                if key in site_names and key != "None":
                    value = int(line_strip[2]) + int(line_strip[3]) - 1
                    site_pos_map[key] = value
        self.m1pos = site_pos_map



        self.container = {"names" : [s.name for s in self.sites],
                          "seqs" : [str(s) for s in self.sites]}
        self.top_counts = {name: {"TGT" : 0, "ACA" : 0} for name in
                           [site.name for site in self.sites] + ["None"]}
        self.multi_counts = {"TGT" : [], "ACA" : []}
        self.single_pos_counts = {name: {"TGT" : [0]*read_len,
                                         "ACA" : [0]*read_len} for name in
                                  [site.name for site in self.sites]}
        self.top_pos_counts = {name: {"TGT" : [0]*read_len,
                                         "ACA" : [0]*read_len} for name in
                                  [site.name for site in self.sites]}


        self.flank_counts = {barcode : {site.name : {flank : 0
                                                for flank in SiteList.FLANKS}
                                        for site in self.sites}
                             for barcode in ["TGT", "ACA"]}
    def __getitem__(self, key):
        return self.container[key]

    def list_site_seqs(self):
        return [s.seq for s in self.sites]

    def reorder(self, kds):
        kds_full = kds.loc[:,"Full"].iloc[:kds.shape[0] - 3]
        sites_and_kds = list(zip(self.sites, kds_full))
        sites_sorted = sorted(sites_and_kds, key=lambda sites_and_kds: sites_and_kds[1])
        sites_new = [i[0] for i in sites_sorted]
        self.sites = sites_new

    def get_best_site_ind(self, site_names):
        # Returns a list of indeces of the best sites:
        site_rank_map = {s.name: i for (i, s) in enumerate(self.sites)}
        ranks = [site_rank_map[s] for s in site_names]
        ind_best = [i for (i, r) in enumerate(ranks) if r == min(ranks)]
        return(ind_best)

    def top_counts_df(self):
        output = pd.DataFrame.from_dict(self.top_counts, orient="index")
        return(output.reindex(self["names"]+["None"], copy=False))

    def flank_counts_df(self, barcode):
        output = pd.DataFrame.from_dict(self.flank_counts[barcode], orient="index")
        return(output.transpose().reindex(SiteList.FLANKS, copy=False))

    def single_pos_counts_df(self, barcode="TGT"):
        output = pd.DataFrame.from_dict({name : self.single_pos_counts[name][barcode] for
                                         name in list(self.single_pos_counts.keys())},
                                        orient="index")
        return(output)

    def top_pos_counts_df(self, barcode="TGT"):
        output = pd.DataFrame.from_dict({name : self.top_pos_counts[name][barcode] for
                                         name in list(self.top_pos_counts.keys())},
                                        orient="index")
        return(output)


    def multi_counts_df(self):
        keys = list(set(self.multi_counts["TGT"] + self.multi_counts["ACA"]))
        counts = {name: {"TGT" : 0, "ACA" : 0} for name in keys + ["None"]}
        for multi in self.multi_counts["TGT"]:
            counts[multi]["TGT"] += 1
        for multi in self.multi_counts["ACA"]:
            counts[multi]["ACA"] += 1
        output = pd.DataFrame.from_dict(counts, orient="index").astype(int)
        return(output)

class SiteListRep:
    BASE_DIR = HOME_DIR + "AssignSiteTypes/site_categories"
    FLANKS = ["%s.%s" %(i[:2], i[2:]) for i in get_kmer_list(4)]
    def __init__(self, _mirna, list_name, exclude=[]):
        mirname = _mirna.name
        if list_name == "topcanonical":
            site_names = ["8mer", "7mer-m8", "7mer-A1"]
        elif list_name == "fourcanonical":
            site_names = ["8mer", "7mer-m8", "7mer-A1", "6mer"]
        else:
            f_name = "%s/%s/sites.%s_%s.txt" %(self.BASE_DIR, list_name,
                                     "-".join(mirname.split("-")[:2]), list_name)
            with open(f_name) as f:
                site_names = f.read().split("\n")
        self.name = list_name
        self.mirna = _mirna.name
        self.sites = [Site(_mirna, s) for s in site_names]
        self.container = {"names" : [s.name for s in self.sites],
                          "seqs" : [str(s) for s in self.sites]}
        self.top_counts = {name: 0 for name in
                           [site.name for site in self.sites] + ["None"]}
        self.multi_counts = []
        self.flank_counts = {site.name : {flank : 0
                                                for flank in SiteList.FLANKS}
                                        for site in self.sites}
    def __getitem__(self, key):
        return self.container[key]

    def list_site_seqs(self):
        return [s.seq for s in self.sites]

    def reorder(self, kds):
        kds_full = kds.loc[:,"Full"].iloc[:kds.shape[0] - 3]
        sites_and_kds = list(zip(self.sites, kds_full))
        sites_sorted = sorted(sites_and_kds, key=lambda sites_and_kds: sites_and_kds[1])
        sites_new = [i[0] for i in sites_sorted]
        self.sites = sites_new

    def get_best_site_ind(self, site_names):
        # Returns a list of indeces of the best sites:
        site_rank_map = {s.name: i for (i, s) in enumerate(self.sites)}
        ranks = [site_rank_map[s] for s in site_names]
        ind_best = [i for (i, r) in enumerate(ranks) if r == min(ranks)]
        return(ind_best)

    def top_counts_df(self):
        output = pd.DataFrame.from_dict(self.top_counts, orient="index")
        return(output.reindex(self["names"]+["None"], copy=False))

    def flank_counts_df(self, barcode):
        output = pd.DataFrame.from_dict(self.flank_counts[barcode], orient="index")
        return(output.transpose().reindex(SiteList.FLANKS, copy=False))

    def single_pos_counts_df(self):
        output = pd.DataFrame.from_dict({name : self.single_pos_counts[name]["TGT"] for
                                         name in list(self.single_pos_counts.keys())},
                                        orient="index")
        return(output)

    def top_pos_counts_df(self):
        output = pd.DataFrame.from_dict({name : self.top_pos_counts[name]["TGT"] for
                                         name in list(self.top_pos_counts.keys())},
                                        orient="index")
        return(output)


    def multi_counts_df(self):
        keys = list(set(self.multi_counts["TGT"] + self.multi_counts["ACA"]))
        counts = {name: {"TGT" : 0, "ACA" : 0} for name in keys + ["None"]}
        for multi in self.multi_counts["TGT"]:
            counts[multi]["TGT"] += 1
        for multi in self.multi_counts["ACA"]:
            counts[multi]["ACA"] += 1
        output = pd.DataFrame.from_dict(counts, orient="index")
        return(output)


class TwelvemerSiteList:
    def __init__(self, _mirna, mir_start):
        self.name = "12mers_m%s.%s" %(mir_start, int(mir_start) + 3)
        self.mirna = _mirna.name
        m_l = int(mir_start) - 1
        m_r = m_l + 4
        self.l = m_l
        self.r = m_r
        self.constant = _mirna.rev_com(start=m_l, stop=m_r)
        self.sites = ["".join([kmer[:10 - m_r], self.constant, kmer[10 - m_r:]])
                      for kmer in get_kmer_list(8)]
        self.top_counts = {name: {"TGT" : 0, "ACA" : 0} for name in
                           self.sites + ["None"]}
        self.multi_counts = {"TGT" : [], "ACA" : []}
        self.container = {"names" : self.sites,
                          "seqs" : self.sites}
    def __getitem__(self, key):
        return self.container[key]
    def top_counts_df(self):
        output = pd.DataFrame.from_dict(self.top_counts, orient="index")
        return(output.reindex(self["names"]+["None"], copy=False))
    def multi_counts_df(self):
        keys = list(set(self.multi_counts["TGT"] + self.multi_counts["ACA"]))
        counts = {name: {"TGT" : 0, "ACA" : 0} for name in keys + ["None"]}
        for multi in self.multi_counts["TGT"]:
            counts[multi]["TGT"] += 1
        for multi in self.multi_counts["ACA"]:
            counts[multi]["ACA"] += 1
        output = pd.DataFrame.from_dict(counts, orient="index")
        return(output)




class SixteenmerSiteList:
    def __init__(self, _mirna, mir_start, split16):
        self.name = "16mers_m%s.%s_%s" %(mir_start, int(mir_start) + 3, split16)
        self.mirna = _mirna.name
        m_l = int(mir_start) - 1
        m_r = m_l + 4
        self.l = m_l
        self.r = m_r
        self.split = split16
        self.constant = _mirna.rev_com(start=m_l, stop=m_r)
        pos = Read.SPLIT16_POS_MAP[self.split]
        self.sites = ["".join([kmer[:pos - m_r], self.constant, kmer[pos - m_r:]])
                      for kmer in get_kmer_list(8)]
        self.top_counts = {name: {"TGT" : 0, "ACA" : 0} for name in
                           self.sites + ["None"]}
        self.multi_counts = {"TGT" : [], "ACA" : []}
        self.container = {"names" : self.sites,
                          "seqs" : self.sites}
    def __getitem__(self, key):
        return self.container[key]
    def top_counts_df(self):
        output = pd.DataFrame.from_dict(self.top_counts, orient="index")
        return(output.reindex(self["names"]+["None"], copy=False))
    def multi_counts_df(self):
        keys = list(set(self.multi_counts["TGT"] + self.multi_counts["ACA"]))
        counts = {name: {"TGT" : 0, "ACA" : 0} for name in keys + ["None"]}
        for multi in self.multi_counts["TGT"]:
            counts[multi]["TGT"] += 1
        for multi in self.multi_counts["ACA"]:
            counts[multi]["ACA"] += 1
        output = pd.DataFrame.from_dict(counts, orient="index")
        return(output)

def get_kmer_index(kmer):
    nt_dict = {"A": 0, "C": 1, "G": 2, "T": 3}
    return(sum([4**(len(kmer) - i - 1)*nt_dict[nt]
              for i, nt in enumerate(kmer)]))
    

class KmerList:
    BCS = ["TGT", "ACA"]
    def __init__(self, kmer_len, n_constant, rand_length, fast=False,
                 buffer_=False):
        self.len = int(kmer_len)
        self.n_constant = n_constant
        self.rand_length= rand_length
        if buffer_:
            self.positions = n_constant + rand_length - 1 - kmer_len + 1
        else:
            self.positions = rand_length + 2*n_constant - kmer_len + 1
        num_kmers = 4**self.len
        self.kmers = {bc : [0]*num_kmers for bc in self.BCS}
        self.pos_mean = {bc : [0]*num_kmers for bc in self.BCS}
        self.pos_sd = {bc : [0]*num_kmers for bc in self.BCS}
        self.pos_mean_dif = {bc : [0]*num_kmers for bc in self.BCS}
        self.pos_sd_dif = {bc : [0]*num_kmers for bc in self.BCS}
        self.sum = {bc : 0 for bc in self.BCS}

    def __getitem__(self, bc):
        return self.kmers[bc]

    def __setitem__(self, bc, val):
        self.kmers[bc] = val
        return self

    def __len__(self):
        return self.len

    def __add__(self, other):
        if other.__class__.__name__ == "KmerList":
            _new = KmerList(len(self), self.n_constant, self.rand_length)
            for bc in self.BCS:
                _new.kmers[bc] = [sum(i) for i in zip(self.kmers[bc], other.kmers[bc])]
                for i, tup in enumerate(zip(self.kmers[bc], self.pos_mean[bc],
                                            self.pos_sd[bc], other.kmers[bc],
                                            other.pos_mean[bc], other.pos_sd[bc])):
                    n_1, mean_1, sd_1, n_2, mean_2, sd_2 = tup
                    if n_1 > 0 or n_2 > 0:
                        mean_x2_1   = sd_1**2 + mean_1**2 
                        mean_x2_2   = sd_2**2 + mean_2**2 
                        mean_x2_new = (mean_x2_1*n_1 + mean_x2_2*n_2)/float(n_1 + n_2)
                        mean_new    = (mean_1*n_1 + mean_2*n_2)/float(n_1 + n_2)
                        sd_new      = (mean_x2_new - mean_new**2)**0.5
                        _new.pos_mean[bc][i] = mean_new
                        _new.pos_sd[bc][i] = sd_new
            _new.sum = {bc : sum(_new.kmers[bc]) for bc in self.BCS}
            return(_new)

    def __radd__(self, other):
        if other == 0:
            return(self)
        elif other.__class__.__name__ == "KmerList":
            return(self.__add__(other))

    def __iadd__(self, other):
        time_start = time.time()
        if other.__class__.__name__ == "KmerList":
            for bc in self.BCS:
                for i, tup in enumerate(zip(self.kmers[bc], self.pos_mean[bc], self.pos_sd[bc],
                                            other.kmers[bc], other.pos_mean[bc],
                                            other.pos_sd[bc])):
                    n_1, mean_1, sd_1, n_2, mean_2, sd_2 = tup
                    if n_1 > 0 or n_2 > 0:
                        mean_x2_1   = sd_1**2 + mean_1**2 
                        mean_x2_2   = sd_2**2 + mean_2**2 
                        mean_x2_new = (mean_x2_1*n_1 + mean_x2_2*n_2)/float(n_1 + n_2)
                        mean_new    = (mean_1*n_1 + mean_2*n_2)/float(n_1 + n_2)
                        sd_new      = (mean_x2_new - mean_new**2)**0.5
                        self.pos_mean[i] = mean_new
                        self.pos_sd[i] = sd_new
            for bc in self.BCS:
                self.kmers[bc] = [sum(i) for i in zip(self.kmers[bc], other.kmers[bc])]
            self.sum = (sum(_new.kmers[bc] for bc in self.BCS))
            return(self)
        elif other.__class__.__name__ == "Read":
            seq = other.seq
            bc = other.barcode
            for pos in range(self.positions):
                # print(pos)
                # time_1 = time.time()

                kmer = seq[pos:pos + len(self)]
                ind = get_kmer_index(kmer)
                n = self[bc][ind]
                # time_2 = time.time()
                # print("Got the new kmers: %f" %(time_2 - time_1))
                mean_old = self.pos_mean[bc][ind]
                sd_old = self.pos_sd[bc][ind]
                mean_x2_old = sd_old**2 + mean_old**2
                mean_x2_new = mean_x2_old*float(n)/(float(n + 1)) + float(pos**2)/(n + 1)
                mean_new = mean_old*float(n)/float(n + 1) + float(pos)/(n + 1)
                sd_new = (mean_x2_new - mean_new**2)**0.5
                # time_3 = time.time()
                # print("Calculated the new summary statistics: %f" %(time_3 - time_2))
                self.pos_mean[bc][ind] = mean_new
                self.pos_sd[bc][ind] = sd_new
                self[bc][ind] += 1
                self.sum[bc] += 1
                # time_4 = time.time()
                # print("Updated the kmer object: %f" %(time_4 - time_3))
            # time_end = time.time()
            # print("About to exit __iadd__: %f" %(time_end - time_start))
            return(self)

    def __truediv__(self, other):
        if other.__class__.__name__ != "KmerList":
            return("ERROR, both objects must be KmerLists")
        else:
            # if self.fast and other.fast:
            _kmerlist_div = KmerList(self.len, self.n_constant,
                                     self.rand_length, fast=True)
            sum_1 = {bc: float(self.sum[bc])
                     for bc in ["TGT", "ACA"]}
            sum_2 = {bc: float(other.sum[bc]) + 4**other.len
                     for bc in ["TGT", "ACA"]}
            kmers = get_kmer_list(self.len)
            for bc in ["TGT", "ACA"]:
                for i, tup in enumerate(zip(self.kmers[bc], self.pos_mean[bc],
                                            self.pos_sd[bc], other.kmers[bc],
                                            other.pos_mean[bc],
                                            other.pos_sd[bc])):
                    n_1, mean_1, sd_1, n_2, mean_2, sd_2 = tup
                    mean_dif = mean_2 - mean_1
                    _kmerlist_div.pos_mean_dif[bc][i] = mean_dif
                    if n_1 > 0 and n_2 > 0:
                        sd_dif = (sd_1**2/float(n_1) + sd_2**2/float(n_2))**0.5
                    else:
                        sd_dif = 0
                    if sum_1[bc] != 0:
                        _kmerlist_div.kmers[bc][i] = (float(n_1)/sum_1[bc])/((float(n_2 + 1)/sum_2[bc]))
                    else:
                        _kmerlist_div.kmers[bc][i] = 0
                _kmerlist_div.pos_mean[bc] = self.pos_mean[bc]
                _kmerlist_div.pos_sd[bc] = self.pos_sd[bc]
            return(_kmerlist_div)

    def write(self, path, bc="TGT"):
        with open(path, "w+") as file_out:
            file_out.write("\t".join(["kmer", "counts_p", "meanpos_p",
                                      "sdpos_p", "counts_c", "meanpos_c",
                                      "sdpos_c"]) + "\n")
            file_out.write("\n".join(
                ["%s\t%s\t%s\t%s\t%s\t%s\t%s" %(kmer, count_p, m_p, sd_p,
                                                count_c, m_c, sd_c)
                 for (kmer, count_p, m_p, sd_p, count_c, m_c, sd_c)
                 in zip(get_kmer_list(self.len), self.kmers["TGT"],
                        self.pos_mean["TGT"], self.pos_sd["TGT"],
                        self.kmers["ACA"], self.pos_mean["ACA"],
                        self.pos_sd["ACA"])]))

    def read(self, mirna, experiment, condition, n_constant, n_ex, uniq=False,
             buffer=False, final=False, cutoffdil=False, temp_path=False,
             resub=False):
        if temp_path:
            kmers_path = "I_combined_temp/%s" %temp_path
        else:
            extension = "_%s_k%s" %(n_constant, self.len)
            if uniq:
                extension = "%s_uniq" %extension
            if buffer:
                extension = "%s_buffer3p" %extension
            if int(n_ex) != 0:
                extension = "%s_ex%s" %(extension, n_ex)
            if final:
                extension = "%s_final" %extension
            if cutoffdil:
                extension = "%s_dil" %extension
            if resub:
                analysis_type = "kmers_cutoff_resub"
            else:
                analysis_type = "kmers_cutoff_final"
            kmers_path = get_analysis_path(mirna.name, experiment, condition,
                                           analysis_type, ext=extension)
        # print((os.system("ls -l %s" %(kmers_path))))
        if not os.path.isfile(kmers_path):
            raise ValueError("%s file not found" %kmers_path)
        with open(kmers_path, "r+") as file_in:
            lines = list(file_in)
            if lines[0].split("\t")[0] == "kmer":
                lines = lines[1:]
            else:
                raise ValueError("Incorrect header on file.")
            line_sum = sum([int(line.split("\t")[1]) for line in lines])
            self.kmers["TGT"] = [int(line.split("\t")[1]) for line in lines]
            self.pos_mean["TGT"] = [float(line.split("\t")[2]) for line in lines]
            self.pos_sd["TGT"] = [float(line.split("\t")[3]) for line in lines]
            self.kmers["ACA"] = [int(line.split("\t")[4]) for line in lines]
            self.pos_mean["ACA"] = [float(line.split("\t")[5]) for line in lines]
            self.pos_sd["ACA"] = [float(line.split("\t")[6]) for line in lines]
            self.sum["TGT"]   = sum(self.kmers["TGT"])
            self.sum["ACA"]   = sum(self.kmers["ACA"])
            return(self)

    def top_and_adjacent(self, bc="TGT"):
        nucs = ["A", "C", "G", "T"]
        kmer_tup = list(zip(get_kmer_list(self.len), self.kmers[bc], self.pos_mean[bc],
                       self.pos_sd[bc], self.pos_mean_dif[bc], self.pos_sd_dif[bc]))
        k_sorted = sorted(kmer_tup, key=lambda kmer: kmer[1], reverse=True)[:10]
        top_kmer = k_sorted[0][0]
        adjacent_kmers = ([top_kmer] +
                          ["%s%s" %(i, top_kmer[:-1]) for i in nucs] +
                          ["%s%s" %(top_kmer[1:], i) for i in nucs])
        inds = [get_kmer_index(i) for i in adjacent_kmers]
        Rs = [kmer_tup[i][1] for i in inds]
        Rs_init = Rs[0]
        Rs = [i/Rs_init for i in Rs]
        spaces = [1, 0, 0, 0, 0, 2, 2, 2, 2]
        adjacent_tups = list(zip(adjacent_kmers, Rs, spaces))
        for tup in adjacent_tups:
            print((" "*tup[2] + "%s\t%6.2f" %(tup[0], tup[1])))
        return

    def top_substitute(self, bc="TGT"):
        nucs = ["A", "C", "G", "T"]
        kmer_tup = list(zip(get_kmer_list(self.len), self.kmers[bc], self.pos_mean[bc],
                       self.pos_sd[bc], self.pos_mean_dif[bc], self.pos_sd_dif[bc]))
        k_sorted = sorted(kmer_tup, key=lambda kmer: kmer[1], reverse=True)[:10]
        top_kmer = k_sorted[0][0]
        adjacent_kmers = ([top_kmer] +
                          ["%s%s" %(i, top_kmer[1:]) for i in nucs] +
                          ["%s%s" %(top_kmer[:-1], i) for i in nucs])
        inds = [get_kmer_index(i) for i in adjacent_kmers]
        Rs = [kmer_tup[i][1] for i in inds]
        Rs_init = Rs[0]
        adjacent_kmers = adjacent_kmers[1:]
        Rs = Rs[1:]
        Rs = [i/Rs_init for i in Rs]
        first_letter = [i[0] for i in adjacent_kmers]
        last_letter = [i[-1] for i in adjacent_kmers]
        for i in range(len(Rs)):
            if first_letter[i] != top_kmer[0]:
                first_letter[i] = first_letter[i].lower()
            if last_letter[i] != top_kmer[-1]:
                last_letter[i] = last_letter[i].lower()
        adjacent_kmer_strings = [first_letter[i] + top_kmer[1:-1] + last_letter[i] for i in range(len(Rs))]
        adjacent_tups = list(zip(adjacent_kmer_strings, Rs))
        for tup in adjacent_tups:
            print(("%s\t%6.2f" %(tup[0], tup[1])))
        return



    def head_with_mirna_sites(self, _mirna, bc="TGT", n=10, sub=False,
                              trim=False, start=False):
        if trim:
            n = 1
        kmer_tup = list(zip(get_kmer_list(self.len), self.kmers[bc], self.pos_mean[bc],
                       self.pos_sd[bc], self.pos_mean_dif[bc], self.pos_sd_dif[bc]))
        k_sorted = sorted(kmer_tup, key=lambda kmer: kmer[1], reverse=True)[:n]
        if not trim or start:
            print((" "*(10 + len(_mirna)-8)+"87654321" + " "*52 + "_"*10 +
                  "POSITION" + "_"* 10))
            print((" "*10 + _mirna.rev() + " "*8 + "site" + " "*16 + "score" +
                  " "*(self.len + 3) + "R" + " "*7 + "p_mean" + " "*2 + "p_sd" +
                  " "*4 + "d_mean" + " "*2 + "d_sd"))
        top = []
        for k_sort in k_sorted:    
            out = _mirna.get_best_name(k_sort, suboptimal=sub)
            top.append(out[0])
        counts = [top.count(i) for i in top]
        # print(counts)
        # print((set(counts)))
        if list(set(counts)) == [1]:
            most_common_site = "NA"
        else: 
            most_common_site = max(set(top), key=top.count)
        top_site = top[0]
        if not trim:
            print(("top site:\t%s" %(top_site)))
            print(("common site:\t%s" %(most_common_site)))

    def output_for_table(self, _mirna, num, bc="TGT", n=1, sub=False, start=False):
        kmer_tup = list(zip(get_kmer_list(self.len), self.kmers[bc], self.pos_mean[bc],
                       self.pos_sd[bc], self.pos_mean_dif[bc], self.pos_sd_dif[bc]))
        k_sorted = sorted(kmer_tup, key=lambda kmer: kmer[1], reverse=True)[:n]
        if start:
            # print((" "*(10 + len(_mirna)-8)+"87654321" + " "*52 + "_"*10 +
            #       "POSITION" + "_"* 10))
            print(("No." + " "*7 + "kmer" + " "*7 + "R" + " "*2 + "p_mean" + " "*4 +
                   "p_sd" + " "*16 + "site"))

        out = _mirna.get_best_name(k_sorted[0], suboptimal=sub,
                                   print_sites=False, output_for_table=True)
        # print("\t".join([str(i) for i in k_sorted[0]] + out))
        kmer, R, pos_mean, pos_sd = k_sorted[0][:4]
        print("%3.0f" %num + " %s" %kmer + "%8.2f" %R + "%8.1f" %pos_mean +
              "%8.1f" %pos_sd + " "*(20 - len(out[0])) + "%s" %out[0])



    def df(self):
        return(pd_dict(self.kmers).transpose()[["TGT", "ACA"]])


class KmerListTest:
    BCS = ["TGT", "ACA"]
    def __init__(self, kmer_len, n_constant, rand_length, fast=False,
                 buffer_=False):
        self.len = int(kmer_len)
        self.n_constant = n_constant
        self.rand_length= rand_length
        if buffer_:
            self.positions = n_constant + rand_length - 1 - kmer_len + 1
        else:
            self.positions = rand_length + 2*n_constant - kmer_len + 1
        num_kmers = 4**self.len
        self.kmers = {bc : [0]*num_kmers for bc in self.BCS}
        self.pos_mean = {bc : [0]*num_kmers for bc in self.BCS}
        self.pos_sd = {bc : [0]*num_kmers for bc in self.BCS}
        self.pos_mean_dif = {bc : [0]*num_kmers for bc in self.BCS}
        self.pos_sd_dif = {bc : [0]*num_kmers for bc in self.BCS}
        self.sum = {bc : 0 for bc in self.BCS}

    def __getitem__(self, bc):
        return self.kmers[bc]

    def __setitem__(self, bc, val):
        self.kmers[bc] = val
        return self

    def __len__(self):
        return self.len

    def __add__(self, other):
        if other.__class__.__name__ in ["KmerList", "KmerListTest"]:
            _new = KmerList(len(self), self.n_constant, self.rand_length)
            for bc in self.BCS:
                _new.kmers[bc] = [sum(i) for i in zip(self.kmers[bc], other.kmers[bc])]
                for i, tup in enumerate(zip(self.kmers[bc], self.pos_mean[bc],
                                            self.pos_sd[bc], other.kmers[bc],
                                            other.pos_mean[bc], other.pos_sd[bc])):
                    n_1, mean_1, sd_1, n_2, mean_2, sd_2 = tup
                    if n_1 > 0 or n_2 > 0:
                        mean_x2_1   = sd_1**2 + mean_1**2 
                        mean_x2_2   = sd_2**2 + mean_2**2 
                        mean_x2_new = (mean_x2_1*n_1 + mean_x2_2*n_2)/float(n_1 + n_2)
                        mean_new    = (mean_1*n_1 + mean_2*n_2)/float(n_1 + n_2)
                        sd_new      = (mean_x2_new - mean_new**2)**0.5
                        _new.pos_mean[bc][i] = mean_new
                        _new.pos_sd[bc][i] = sd_new
            _new.sum = {bc : sum(_new.kmers[bc]) for bc in self.BCS}
            return(_new)

    def __radd__(self, other):
        if other == 0:
            return(self)
        elif other.__class__.__name__ in ["KmerList", "KmerListTest"]:
            return(self.__add__(other))

    def __iadd__(self, other):
        print("in __iadd__")
        time_start = time.time()
        if other.__class__.__name__ in ["KmerList", "KmerListTest"]:
            for bc in self.BCS:
                for i, tup in enumerate(zip(self.kmers[bc], self.pos_mean[bc], self.pos_sd[bc],
                                            other.kmers[bc], other.pos_mean[bc],
                                            other.pos_sd[bc])):
                    n_1, mean_1, sd_1, n_2, mean_2, sd_2 = tup
                    if n_1 > 0 or n_2 > 0:
                        mean_x2_1   = sd_1**2 + mean_1**2 
                        mean_x2_2   = sd_2**2 + mean_2**2 
                        mean_x2_new = (mean_x2_1*n_1 + mean_x2_2*n_2)/float(n_1 + n_2)
                        mean_new    = (mean_1*n_1 + mean_2*n_2)/float(n_1 + n_2)
                        sd_new      = (mean_x2_new - mean_new**2)**0.5
                        self.pos_mean[i] = mean_new
                        self.pos_sd[i] = sd_new
            for bc in self.BCS:
                self.kmers[bc] = [sum(i) for i in zip(self.kmers[bc], other.kmers[bc])]
            self.sum = (sum(_new.kmers[bc] for bc in self.BCS))
            return(self)
        elif other.__class__.__name__ == "Read":
            seq = other.seq
            bc = other.barcode
            time_quick_loop_start = time.time()
            inds = [get_kmer_index(seq[pos:pos + len(self)]) for pos in range(self.positions)]
            for ind in inds:
                self[bc][ind] += 1
            time_quick_loop_end = time.time()
            print("Time for quick loop: %f" %(time_quick_loop_end - time_quick_loop_start))   
            sys.stdout.flush()         
            time_loop_start = time.time()
            for pos in range(self.positions):
                time_1 = time.time()

                kmer = seq[pos:pos + len(self)]
                ind = get_kmer_index(kmer)
                # n = self[bc][ind]
                time_2 = time.time()
                # print("Got the new kmers: %f" %(time_2 - time_1))
                # mean_old = self.pos_mean[bc][ind]
                # sd_old = self.pos_sd[bc][ind]
                # mean_x2_old = sd_old**2 + mean_old**2
                # mean_x2_new = mean_x2_old*float(n)/(float(n + 1)) + float(pos**2)/(n + 1)
                # mean_new = mean_old*float(n)/float(n + 1) + float(pos)/(n + 1)
                # sd_new = (mean_x2_new - mean_new**2)**0.5
                time_3 = time.time()
                # print("Calculated the new summary statistics: %f" %(time_3 - time_2))
                # self.pos_mean[bc][ind] = mean_new
                # self.pos_sd[bc][ind] = sd_new
                self[bc][ind] += 1
                # self.sum[bc] += 1
                time_4 = time.time()
                # print("Updated the kmer object: %f" %(time_4 - time_3))
            time_loop_end = time.time()
            print("Time for loop: %f" %(time_loop_end - time_loop_start))
            time_end = time.time()
            print("About to exit __iadd__: %f" %(time_end - time_start))
            sys.stdout.flush()
            return(self)

    def __truediv__(self, other):
        if other.__class__.__name__ not in ["KmerList", "KmerListTest"]:
            return("ERROR, both objects must be KmerLists")
        else:
            # if self.fast and other.fast:
            _kmerlist_div = KmerList(self.len, self.n_constant,
                                     self.rand_length, fast=True)
            sum_1 = {bc: float(self.sum[bc])
                     for bc in ["TGT", "ACA"]}
            sum_2 = {bc: float(other.sum[bc]) + 4**other.len
                     for bc in ["TGT", "ACA"]}
            kmers = get_kmer_list(self.len)
            for bc in ["TGT", "ACA"]:
                for i, tup in enumerate(zip(self.kmers[bc], self.pos_mean[bc],
                                            self.pos_sd[bc], other.kmers[bc],
                                            other.pos_mean[bc],
                                            other.pos_sd[bc])):
                    n_1, mean_1, sd_1, n_2, mean_2, sd_2 = tup
                    mean_dif = mean_2 - mean_1
                    _kmerlist_div.pos_mean_dif[bc][i] = mean_dif
                    if n_1 > 0 and n_2 > 0:
                        sd_dif = (sd_1**2/float(n_1) + sd_2**2/float(n_2))**0.5
                    else:
                        sd_dif = 0
                    if sum_1[bc] != 0:
                        _kmerlist_div.kmers[bc][i] = (float(n_1)/sum_1[bc])/((float(n_2 + 1)/sum_2[bc]))
                    else:
                        _kmerlist_div.kmers[bc][i] = 0
                _kmerlist_div.pos_mean[bc] = self.pos_mean[bc]
                _kmerlist_div.pos_sd[bc] = self.pos_sd[bc]
            return(_kmerlist_div)

    def write(self, path, bc="TGT"):
        with open(path, "w+") as file_out:
            file_out.write("\t".join(["kmer", "counts_p", "meanpos_p",
                                      "sdpos_p", "counts_c", "meanpos_c",
                                      "sdpos_c"]) + "\n")
            file_out.write("\n".join(
                ["%s\t%s\t%s\t%s\t%s\t%s\t%s" %(kmer, count_p, m_p, sd_p,
                                                count_c, m_c, sd_c)
                 for (kmer, count_p, m_p, sd_p, count_c, m_c, sd_c)
                 in zip(get_kmer_list(self.len), self.kmers["TGT"],
                        self.pos_mean["TGT"], self.pos_sd["TGT"],
                        self.kmers["ACA"], self.pos_mean["ACA"],
                        self.pos_sd["ACA"])]))

    def read(self, mirna, experiment, condition, n_constant, n_ex, uniq=False,
             buffer=False, final=False, temp_path=False):
        if temp_path:
            kmers_path = "I_combined_temp/%s" %temp_path
        else:
            extension = "_%s_k%s" %(n_constant, self.len)
            if uniq:
                extension = "%s_uniq" %(extension)
            if buffer:
                extension = "%s_buffer3p" %(extension)
            if int(n_ex) != 0:
                extension = "%s_ex%s" %(extension, n_ex)
            if final:
                extension = "%s_final" %(extension)
            analysis_type = "kmers_cutoff_final"
            kmers_path = get_analysis_path(mirna.name, experiment, condition,
                                           analysis_type, ext=extension)
        print(kmers_path)
        print((os.system("ls -l %s" %(kmers_path))))
        if not os.path.isfile(kmers_path):
            raise ValueError("%s%s file not found" %(
                condition, extension
            ))
        with open(kmers_path, "r+") as file_in:
            lines = list(file_in)
            if lines[0].split("\t")[0] == "kmer":
                lines = lines[1:]
            else:
                raise ValueError("Incorrect header on file.")
            line_sum = sum([int(line.split("\t")[1]) for line in lines])
            self.kmers["TGT"] = [int(line.split("\t")[1]) for line in lines]
            self.pos_mean["TGT"] = [float(line.split("\t")[2]) for line in lines]
            self.pos_sd["TGT"] = [float(line.split("\t")[3]) for line in lines]
            self.kmers["ACA"] = [int(line.split("\t")[4]) for line in lines]
            self.pos_mean["ACA"] = [float(line.split("\t")[5]) for line in lines]
            self.pos_sd["ACA"] = [float(line.split("\t")[6]) for line in lines]
            self.sum["TGT"]   = sum(self.kmers["TGT"])
            self.sum["ACA"]   = sum(self.kmers["ACA"])
            return(self)

    def top_and_adjacent(self, bc="TGT"):
        nucs = ["A", "C", "G", "T"]
        kmer_tup = list(zip(get_kmer_list(self.len), self.kmers[bc], self.pos_mean[bc],
                       self.pos_sd[bc], self.pos_mean_dif[bc], self.pos_sd_dif[bc]))
        k_sorted = sorted(kmer_tup, key=lambda kmer: kmer[1], reverse=True)[:10]
        top_kmer = k_sorted[0][0]
        adjacent_kmers = ([top_kmer] +
                          ["%s%s" %(i, top_kmer[:-1]) for i in nucs] +
                          ["%s%s" %(top_kmer[1:], i) for i in nucs])
        inds = [get_kmer_index(i) for i in adjacent_kmers]
        Rs = [kmer_tup[i][1] for i in inds]
        Rs_init = Rs[0]
        Rs = [i/Rs_init for i in Rs]
        spaces = [1, 0, 0, 0, 0, 2, 2, 2, 2]
        adjacent_tups = list(zip(adjacent_kmers, Rs, spaces))
        for tup in adjacent_tups:
            print((" "*tup[2] + "%s\t%6.2f" %(tup[0], tup[1])))
        return

    def top_substitute(self, bc="TGT"):
        nucs = ["A", "C", "G", "T"]
        kmer_tup = list(zip(get_kmer_list(self.len), self.kmers[bc], self.pos_mean[bc],
                       self.pos_sd[bc], self.pos_mean_dif[bc], self.pos_sd_dif[bc]))
        k_sorted = sorted(kmer_tup, key=lambda kmer: kmer[1], reverse=True)[:10]
        top_kmer = k_sorted[0][0]
        adjacent_kmers = ([top_kmer] +
                          ["%s%s" %(i, top_kmer[1:]) for i in nucs] +
                          ["%s%s" %(top_kmer[:-1], i) for i in nucs])
        inds = [get_kmer_index(i) for i in adjacent_kmers]
        Rs = [kmer_tup[i][1] for i in inds]
        Rs_init = Rs[0]
        adjacent_kmers = adjacent_kmers[1:]
        Rs = Rs[1:]
        Rs = [i/Rs_init for i in Rs]
        first_letter = [i[0] for i in adjacent_kmers]
        last_letter = [i[-1] for i in adjacent_kmers]
        for i in range(len(Rs)):
            if first_letter[i] != top_kmer[0]:
                first_letter[i] = first_letter[i].lower()
            if last_letter[i] != top_kmer[-1]:
                last_letter[i] = last_letter[i].lower()
        adjacent_kmer_strings = [first_letter[i] + top_kmer[1:-1] + last_letter[i] for i in range(len(Rs))]
        adjacent_tups = list(zip(adjacent_kmer_strings, Rs))
        for tup in adjacent_tups:
            print(("%s\t%6.2f" %(tup[0], tup[1])))
        return



    def head_with_mirna_sites(self, _mirna, bc="TGT", n=10, sub=False):
        kmer_tup = list(zip(get_kmer_list(self.len), self.kmers[bc], self.pos_mean[bc],
                       self.pos_sd[bc], self.pos_mean_dif[bc], self.pos_sd_dif[bc]))
        k_sorted = sorted(kmer_tup, key=lambda kmer: kmer[1], reverse=True)[:n]
        print((" "*(10 + len(_mirna)-8)+"87654321" + " "*52 + "_"*10 +
              "POSITION" + "_"* 10))
        print((" "*10 + _mirna.rev() + " "*8 + "site" + " "*16 + "score" +
              " "*(self.len + 3) + "R" + " "*7 + "p_mean" + " "*2 + "p_sd" +
              " "*4 + "d_mean" + " "*2 + "d_sd"))
        top = []
        for k_sort in k_sorted:    
            out = _mirna.get_best_name(k_sort, suboptimal=sub)
            top.append(out[0])
        counts = [top.count(i) for i in top]
        print(counts)
        print((set(counts)))
        if list(set(counts)) == [1]:
            most_common_site = "NA"
        else: 
            most_common_site = max(set(top), key=top.count)
        top_site = top[0]
        print(("top site:\t%s" %(top_site)))
        print(("common site:\t%s" %(most_common_site)))
    def df(self):
        return(pd_dict(self.kmers).transpose()[["TGT", "ACA"]])



class KmerPositionalList:
    def __init__(self, kmer_len, position_max):
        self.len = int(kmer_len)
        self.kmers = {barcode : {kmer : [0]*position_max for kmer
                                        in get_kmer_list(kmer_len)}
                      for barcode in ["TGT", "ACA"]}
        self.sum = 0
    def __getitem__(self, barcode):
        return self.kmers[barcode]

    def __len__(self):
        return self.len

    def __add__(self, other):
        if other.__class__.__name__ == "KmerPositionalList":
            print("position_max:")
            position_max = len(list(self.kmers["TGT"].items())[0][1])
            print(position_max)
            _new = KmerPositionalList(len(self), position_max)
            _new.kmers = {barcode : {
                kmer : [i_1 + i_2 for (i_1, i_2)
                        in zip(self.kmers[barcode][kmer],
                               other.kmers[barcode][kmer])]
                for kmer in self.kmers["TGT"].keys()}
                          for barcode in ["TGT", "ACA"]}
            return(_new)

    def __radd__(self, other):
        if other == 0:
            return(self)
        elif other.__class__.__name__ == "KmerPositionalList":
            return(self.__add__(other))

    def __iadd__(self, other):
        if other.__class__.__name__ == "KmerPositionalList":

            self.kmers = {bc : dict(list(self.kmers[bc].items()) +
                                list(others.kmers[bc].items()))
                          for bc in ["TGT", "ACA"]}
            return(self)
        elif other.__class__.__name__ == "KmerPositionalList":

            self.kmers = {bc : dict(list(self.kmers[bc].items()) +
                                list(other.kmers[bc].items()))
                          for bc in ["TGT", "ACA"]}
            return(self)
        elif other.__class__.__name__ == "Read":
            seq = other.seq
            len_seq = len(seq)
            bc = other.barcode
            for pos in range(len_seq - len(self) + 1):
                kmer_pos = seq[pos:pos + len(self)]
                self[bc][seq[pos:pos + len(self)]] += 1
                self.sum += 1
            return(self)


    def __truediv__(self, other):
        if other.__class__.__name__ != "KmerList":
            return("ERROR, both objects must be KmerLists")
        else:
            _kmerlist_div = KmerList(self.len)
            _kmerlist_div.kmers = {
            "TGT" : {kmer : ((float(self.kmers["TGT"][kmer]) + 1)/
                        (self.sum +4**self.len)/
                        (float(other["TGT"][kmer]) + 1)*
                        (other.sum + 4**self.len))
                        for kmer in list(self["TGT"].keys())},
            "ACA" : {kmer : ((float(self.kmers["ACA"][kmer]) + 1)/
                        (self.sum +4**self.len)/
                        (float(other["ACA"][kmer]) + 1)*
                        (other.sum + 4**self.len))
                        for kmer in list(self["ACA"].keys())}}
            return(_kmerlist_div)


    def write(self, path, bc="TGT"):
        with open(path, "w+") as file_out:
            file_out.write("\n".join(["%s\t%s" %(kmer, count) for (kmer, count)
                                      in sorted(list(self[bc].items()),
                                                key=lambda item: item[0])]))

    def read(self, mirna, experiment, condition, n_constant, n_ex):
        extension = "_%s_k%s" %(n_constant, self.len)
        kmers_path = get_analysis_path(mirna.name, experiment, condition,
                                       "kmers", ext=extension)
        if n_ex == 0:
            with open(kmers_path, "r+") as file_in:
                self.kmers["TGT"] = {line.split("\t")[0] : int(line.split("\t")[1])
                                     for line in file_in.readlines()}
                return(self)
        else:
            kmer_partial_path = kmers_path.split(".txt")[0]
            call_string = "ls %s*" %(kmer_partial_path)
            files = Popen(call_string, shell=True,
                          stdout=PIPE).communicate()[0].strip().split("\n")
            found_file = False
            for file in files:
                commas = len([x for x in re.finditer(",", file)])
                if commas == int(n_ex) - 1 and "ex:" in file:
                    file_final = file
                    found_file = True
            if not found_file:
                return("File not found!")
            print(file_final)
            with open(file_final, "r+") as file_in:
                self.kmers["TGT"] = {line.split("\t")[0] : int(line.split("\t")[1])
                                     for line in file_in.readlines()}
                return(self)

    def head(self, barcode="TGT"):
        kmer_tup = [i for i in list(self.kmers[barcode].items())]
        k_sorted = sorted(kmer_tup, key=lambda kmer: kmer[1], reverse=True)[:10]
        print(("\n".join(["\t".join([str(j) for j in i]) for i in k_sorted])))
        return

    def head_with_mirna_sites(self, _mirna, barcode="TGT"):
        kmer_tup = [i for i in list(self.kmers[barcode].items())]
        k_sorted = sorted(kmer_tup, key=lambda kmer: kmer[1], reverse=True)[:10]
        print((" "*(len(_mirna)-8)+"87654321"))
        print((_mirna.rev()))
        for k_sort in k_sorted:    
            _mirna.get_best_name(k_sort)
            return
        return

    def df(self):
        return(pd_dict(self.kmers).transpose()[["TGT", "ACA"]])



def merge_data_frames(dfs):
    df_out = dfs[0]
    if len(dfs) > 1:
        for df_i in dfs[1:]:
            df_out = df_out.add(df_i, fill_value=0)
    return(df_out)

def write_bulge_sitelist_file(mirna_name):
    _mirna = Mirna(mirna_name)
    seq = re.sub("U", "T", _mirna.seq)
    print(seq)
    # Make dictionary in order to deal with redundant sequences:
    bulgesite_seq_map = dict()
    for i_p, pos in enumerate(seq):
        print(i_p)
        print(pos)
        # Only considers bulges starting at position 3-7:
        if i_p > 1 and i_p < 7:
            for nt in DNTS:
                site = "8mer-b%s%s" %(nt, i_p + 1)
                print(site)
                target_seq = _mirna[site]
                print(target_seq)
                if target_seq not in bulgesite_seq_map:
                    bulgesite_seq_map[target_seq] = []
                bulgesite_seq_map[target_seq] += [site]
    bulge_list = []
    for key, val in list(bulgesite_seq_map.items()):
        if len(val) > 1:
            print(val)
            pattern = "^%s-b" %("8mer")
            pos_range = [int(re.split(pattern, i)[1][1:]) for i in val]
            nt = re.split(pattern, val[0])[1][0]
            val = ["8mer-b%s(%s.%s)" %(nt, min(pos_range), max(pos_range))]
        bulge_list += val
    # Get the two seed 7mer sequences:
    site_7mers = [_mirna[i] for i in ["7mer-A1", "7mer-m8"]]
    # Take only those bulge sites which are in neither 7mer:
    list_final = [i for i in bulge_list
                  if np.sum([_mirna[i].find(j)
                             for j in site_7mers]) == -2]
    # Write the sitelist to the sitelist file:
    out_file = HOME_DIR + "AssignSiteTypes/site_categories/bulge/sites.%s_bulge.txt" %(mirna_name)
    sites_all = Mirna.SEED_SITES[:6] + list_final
    for i in sites_all:
        print(("%s\t%s" %(i, _mirna[i])))
    with open(out_file, "wb") as f:
        f.write("\n".join(sites_all))

def write_del_sitelist_file(mirna_name):
    _mirna = Mirna(mirna_name)
    seq = _mirna.seq
    print(seq)
    # Make dictionary in order to deal with redundant sequences:
    del_site_seq_map = {}
    seed_site_seq_list = [_mirna[site] for site in Mirna.SEED_SITES]
    target_seq = _mirna["8mer"]
    for i_p in range(1, 7):
        del_site_seq = "%s%s" %(target_seq[:7 - i_p], target_seq[8 - i_p:])
        del_site_name = "8mer-d%s" %(i_p + 1)
        print(("%s:\t%s" %(del_site_name, del_site_seq)))
        print(("5p-%s-3p" %(seq)))
        print(("3p-" +
              (del_site_seq[:7 - i_p] + "." + del_site_seq[7 - i_p:])[::-1] +
              "-5p"))
        if del_site_seq not in seed_site_seq_list:
            if del_site_seq not in del_site_seq_map:
                del_site_seq_map[del_site_seq] = [del_site_name]
            else:
                del_site_seq_map[del_site_seq].append(del_site_name)
    list_final = []
    for key, val in list(del_site_seq_map.items()):
        if len(val) > 1:
            # Make the "8mer-d(2.4) style name:
            val = ["8mer-d(%s)" %(".".join([i.split("d")[1] for i in val]))]
        list_final.append(val[0])
    # Write the sitelist to the sitelist file:
    out_file = HOME_DIR + "AssignSiteTypes/site_categories/del/sites.%s_del.txt" %(mirna_name)
    print(out_file)
    sites_all = Mirna.SEED_SITES[:6] + list_final
    for i in sites_all:
        print(("%s\t%s" %(i, _mirna[i])))
    with open(out_file, "wb") as f:
        f.write("\n".join(sites_all))
    with open(out_file, "r") as f:
        print((f.readlines()))
    print(out_file)



################################################################################
class Read(Sequence):   
    LIB_5p       = "GGGCAGAGTTCTACAGTCCGACGATC"
    LIB_3p       = "TCGTATGCCGTCTTCTGCTTG" # Does not include TGT or ACA
    LIB_3p_miR_7 = "TCGTATGCCGCTGTGTGCTTG"
    SEQ_BARCODE_MAP = {
        "TGT" : "TGT",
        "ACA" : "ACA",
        "TCG" : "TGT",
        "TG"  : "TGT",
        ""    : ""
    }
    SPLIT16_POS_MAP = {
        "left" : 12,
        "right" : 8,
    }
    def __init__(self, seq, rand_len, _mirna, n_constant, experiment, buffer3p=False):
        # Define the left and right hand sides of the barcode region:
        self.barcode = self.SEQ_BARCODE_MAP[seq[rand_len:]]
        self.n_constant = n_constant
        if _mirna.name == "miR-1" and experiment == "equilibrium":
            added_string = ""
        else:
            added_string = self.barcode
        # Define the sequence to only have the number of constant sequences:
        n_3p = n_constant
        self.seq = "".join([self.LIB_5p[::-1][:max(0, n_constant)][::-1],
                            seq[-1*min(0, n_constant):rand_len + min(0, n_3p)],
                            (added_string + self.LIB_3p)[:max(0, n_3p)]])
        if (_mirna.name in ["miR-7-23nt", "miR-7-24nt", "miR-7-25nt"] and
            experiment in ["equilibrium2_nb", "equilibrium_tp"]):
            LIB_3p = self.LIB_3p_miR_7
        else:
            LIB_3p = self.LIB_3p
        self.fullseq = "".join([self.LIB_5p,
                            seq[:rand_len],
                            (added_string + LIB_3p)])
        self.sites = []
        self.topsite = None
    #...........................................................................
    class ReadSite(Site):
        def __init__(self, _sitelist, name, l, r):
            self.name = name
            self.seq = Mirna(_sitelist.mirna)[name]
            self.l = l
            self.r = r
    # Methods that extract information from the top site
    # Using the read and site information:
    def site_overlaps(self):
        #                                                  |**| <- This distance must be positive
        # ACACACA                                    l[0]  r[0]
        #                                            ACACACA
        # 5mer-m2.6                                           l[1]r[1]
        #                                                     ATTCC
        # GGGCAGAGTTCTACAGTCCGACGATCTTTATTTGTCACACAATACACACAAAATTCCCCGATATCGTATGCCGTCTTCTGCTTG
        #                                            -------  -----                           
        # FROM DIAGRAM l[0]
        ls = [s.l for s in self.sites][1:]
        rs = [s.r for s in self.sites][:-1]
        distance = [r - l for r, l in zip(rs, ls)]
        return([max(0, i) for i in distance])
    def all_site_names(self, overlap=True):
        # Prints the list of all sites associated with the read:
        names = [s.name for s in self.sites]
        if overlap and sum(self.site_overlaps()):
            o_i = 0
            for i, ol in enumerate(self.site_overlaps()):
                if ol > 0:
                    names[o_i] = "%s|%s|%s" %(names[o_i], ol, names[i + 1])
                    names[i + 1] = "remove"
                else:
                    o_i = i - 1
            names = [i for i in names if i != "remove"]
        return names
    def site_name(self):
        # Gives the name of the top site:
        s = self.topsite
        if s:
            return s.name
    def site_pos(self, full=False):
        # Gives the position of the top site:
        s = self.topsite
        if s:
            l, r = (s.l, s.r)
            if full:
                l, r = (i - self.n_constant + len(self.LIB_5p) for i in (l, r))
            return l, r

    def all_sites_pos(self, full=False):
        # Gives the position of the top site:
        sites = self.sites
        if sites:         
            # if full:
            return [(site.l - self.n_constant + len(self.LIB_5p),
                     site.r - self.n_constant + len(self.LIB_5p))
                    for site in sites]


    def site_flank(self, length=2):
        # Gives the flanking dinucleotides fo the top site:
        s = self.topsite
        if s:
            read = self.seq
            read_full = self.fullseq
            align = read_full.find(read)
            return "%s.%s" %(read_full[s.l + align - length:s.l + align],
                             read_full[s.r + align: s.r + length + align])
    def get_all_sites_in_read(self, _sitelist, ties=False, buffer_=False):
        # Returns a list of matches with ("name", start, stop)
        if _sitelist.__class__.__name__ == "SiteList":
            read_sites = [(_site.name, k, k + len(_site))
                          for _site in _sitelist.sites
                          for k in find_all(self.seq, _site.seq)]
            read_sites = sorted(read_sites, key=lambda read_sites: read_sites[1])
            if read_sites:
                read_sites_temp = [i for i in read_sites]
                read_sites = [read_sites_temp[0]]
                if len(read_sites_temp) > 1:
                    for t in read_sites_temp:
                        l_t, r_t = t[1:]
                        add = True
                        for s in read_sites:
                            l_s, r_s = s[1:]
                            if l_t >= l_s and r_t <= r_s:
                                add = False
                            elif ((l_t <= l_s and r_t > r_s)
                                  or (l_t < l_s and r_t >= r_s)):
                                read_sites.remove(s)
                        if add == True:
                            read_sites.append(t)
                r_all = [site[2] - self.n_constant for site in read_sites]
                print_later = False
                site_names_pre = [site[0] for site in read_sites]
                if buffer_:
                    read_sites_new = []
                    read_sites = [read_site for (read_site, r_i) in zip(read_sites, r_all) if r_i <= 36]
                if read_sites:
                    self.sites = [Read.ReadSite(_sitelist, s[0], s[1], s[2])
                                  for s in read_sites]
                    site_names = [i[0] for i in read_sites]
                    topsite = [read_sites[i] for i in
                                 _sitelist.get_best_site_ind(site_names)]
                    if not ties:
                        topsite = [topsite[0]]
                    self.topsite = [Read.ReadSite(_sitelist, s[0], s[1], s[2])
                                    for s in topsite][0]
                else:
                    self.sites = [Read.ReadSite(_sitelist, "None", None, None)]
                    self.topsite = Read.ReadSite(_sitelist, "None", None, None)
            else:
                self.sites = [Read.ReadSite(_sitelist, "None", None, None)]
                self.topsite = Read.ReadSite(_sitelist, "None", None, None)
        elif _sitelist.__class__.__name__ == "TwelvemerSiteList":
            read_sites = [(self.seq[k + _sitelist.r - 10:k + _sitelist.r + 2],
                           k + _sitelist.r - 10, k + _sitelist.r + 2) for
                          k in find_all(self.seq, _sitelist.constant) if
                          k > 9 - _sitelist.r and # This is RIGHT
                          k < len(self.seq) - _sitelist.r - 1] # THIS IS RIGHT
            if read_sites:
                self.sites = [Read.ReadSite(_sitelist, s[0], s[1], s[2])
                              for s in read_sites]
                # This makes no sense
                self.topsite = self.sites[0]
            else:
                self.sites = [Read.ReadSite(_sitelist, "None", None, None)]
                self.topsite = Read.ReadSite(_sitelist, "None", None, None)
        elif _sitelist.__class__.__name__ == "SixteenmerSiteList":
            split = _sitelist.split
            pos = Read.SPLIT16_POS_MAP[_sitelist.split]
            read_sites = [(self.seq[k + _sitelist.r - pos:k + _sitelist.r + 12 - pos],
                           k + _sitelist.r - pos, k + _sitelist.r +  12 - pos) for
                          k in find_all(self.seq, _sitelist.constant) if
                          k > pos - 1 - _sitelist.r and # This is RIGHT
                          k < len(self.seq) - _sitelist.r - (12 - pos) + 1] # THIS IS RIGHT
            if read_sites:
                self.sites = [Read.ReadSite(_sitelist, s[0], s[1], s[2])
                              for s in read_sites]
                # This makes no sense
                self.topsite = self.sites[0]
            else:
                self.sites = [Read.ReadSite(_sitelist, "None", None, None)]
                self.topsite = Read.ReadSite(_sitelist, "None", None, None)
                
    def print_sites(self, _sitelist):
        # Prints the whole list of sites:
        if self.sites and self.sites[0].name != "None":
            print("printing sites")
            seq_map = " "*len(self.seq)
            for _rs in self.sites:
                print((_rs.name))
                print((_sitelist.name))
                if "12mers" in _sitelist.name:
                    print((_sitelist.l))
                if _rs.name[:6] == "Random":
                    print((" "*(_rs.l + (10 - int(_rs.name.split(".")[1]))) + str(_rs)))
                else:
                    print((" "*_rs.l + str(_rs)))
                l, r = (_rs.l, _rs.r)
                print((" "*l + self.seq[l:r]))
                for pos in range(l, r):
                    if seq_map[pos] == " ":
                        seq_map = seq_map[:pos] + "-" + seq_map[pos+1:]
                    elif seq_map[pos] == "-":
                        seq_map = seq_map[:pos] + "=" + seq_map[pos+1:]
            print((self.seq))
            print(seq_map)
            # print("012345678901234567890123456789012345678901234567890"[:len(self.seq)])
            # print("          1         2         3         4         5"[:len(self.seq)])

            return

class Utr(Sequence):   
    def __init__(self, seq):
        # Define the left and right hand sides of the barcode region:
        self.seq = seq
        self.sites = []
        self.topsite = None
    #...........................................................................
    class UtrSite(Site):
        def __init__(self, _sitelist, name, l, r):
            self.name = name
            self.seq = Mirna(_sitelist.mirna)[name]
            self.l = l
            self.r = r
    # Methods that extract information from the top site
    # Using the read and site information:
    def site_overlaps(self):
        #                                                  |**| <- This distance must be positive
        # ACACACA                                    l[0]  r[0]
        #                                            ACACACA
        # 5mer-m2.6                                           l[1]r[1]
        #                                                     ATTCC
        # GGGCAGAGTTCTACAGTCCGACGATCTTTATTTGTCACACAATACACACAAAATTCCCCGATATCGTATGCCGTCTTCTGCTTG
        #                                            -------  -----                           
        # FROM DIAGRAM l[0]
        ls = [s.l for s in self.sites][1:]
        rs = [s.r for s in self.sites][:-1]
        distance = [r - l for r, l in zip(rs, ls)]
        return([max(0, i) for i in distance])
    def all_site_names(self, overlap=True):
        # Prints the list of all sites associated with the read:
        names = [s.name for s in self.sites]
        if overlap and sum(self.site_overlaps()):
            o_i = 0
            for i, ol in enumerate(self.site_overlaps()):
                if ol > 0:
                    names[o_i] = "%s|%s|%s" %(names[o_i], ol, names[i + 1])
                    names[i + 1] = "remove"
                else:
                    o_i = i - 1
            names = [i for i in names if i != "remove"]
        return names
    def site_name(self):
        # Gives the name of the top site:
        s = self.topsite
        if s:
            return s.name
    def site_pos(self):
        # Gives the position of the top site:
        s = self.topsite
        if s:
            return s.l, s.r
    def site_flank(self, length=2):
        # Gives the flanking dinucleotides fo the top site:
        s = self.topsite
        if s:
            read = self.seq
            return "%s.%s" %(read[s.l - length:s.l], read[s.r:s.r + length])
    def get_all_sites_in_utr(self, _sitelist, ties=False):
        # Returns a list of matches with ("name", start, stop)
        if _sitelist.__class__.__name__ == "SiteListRep":
            utr_sites = [(_site.name, k, k + len(_site))
                          for _site in _sitelist.sites
                          for k in find_all(self.seq, _site.seq)]
            utr_sites = sorted(utr_sites, key=lambda utr_sites: utr_sites[1])
            if utr_sites:
                utr_sites_temp = [i for i in utr_sites]
                utr_sites = [utr_sites_temp[0]]
                if len(utr_sites_temp) > 1:
                    for t in utr_sites_temp:
                        # print("*********************************************")
                        # print(utr_sites)
                        # print("_____________________________________________")
                        # print(t)
                        l_t, r_t = t[1:]
                        # print(l_t)
                        # print(r_t)
                        add = True
                        # print("......")
                        for s in utr_sites:
                            # print(s)
                            l_s, r_s = s[1:]
                            # print(l_s)
                            # print(r_s)
                            if l_t >= l_s and r_t <= r_s:
                                add = False
                            elif ((l_t <= l_s and r_t > r_s)
                                  or (l_t < l_s and r_t >= r_s)):
                                utr_sites.remove(s)
                        if add == True:
                            utr_sites.append(t)
                self.sites = [Utr.UtrSite(_sitelist, s[0], s[1], s[2])
                              for s in utr_sites]
                site_names = [i[0] for i in utr_sites]
                topsite = [utr_sites[i] for i in
                             _sitelist.get_best_site_ind(site_names)]
                if not ties:
                    topsite = [topsite[0]]
                self.topsite = [Utr.UtrSite(_sitelist, s[0], s[1], s[2])
                                for s in topsite][0]
            else:
                self.sites = [Utr.UtrSite(_sitelist, "None", None, None)]
                self.topsite = Utr.UtrSite(_sitelist, "None", None, None)
        elif _sitelist.__class__.__name__ == "TwelvemerSiteListRep":
            read_sites = [(self.seq[k + _sitelist.r - 10:k + _sitelist.r + 2],
                           k + _sitelist.r - 10, k + _sitelist.r + 2) for
                          k in find_all(self.seq, _sitelist.constant) if
                          k > 9 - _sitelist.r and # This is RIGHT
                          k < len(self.seq) - _sitelist.r - 1] # THIS IS RIGHT
            if read_sites:
                self.sites = [Utr.UtrSite(_sitelist, s[0], s[1], s[2])
                              for s in read_sites]
                # This makes no sense
                self.topsite = self.sites[0]
            else:
                self.sites = [Utr.UtrSite(_sitelist, "None", None, None)]
                self.topsite = Utr.UtrSite(_sitelist, "None", None, None)
        elif _sitelist.__class__.__name__ == "SixteenmerSiteListRep":
            split = _sitelist.split
            pos = Utr.SPLIT16_POS_MAP[_sitelist.split]
            utr_sites = [(self.seq[k + _sitelist.r - pos:k + _sitelist.r + 12 - pos],
                           k + _sitelist.r - pos, k + _sitelist.r +  12 - pos) for
                          k in find_all(self.seq, _sitelist.constant) if
                          k > pos - 1 - _sitelist.r and # This is RIGHT
                          k < len(self.seq) - _sitelist.r - (12 - pos) + 1] # THIS IS RIGHT
            if utr_sites:
                self.sites = [Utr.UtrSite(_sitelist, s[0], s[1], s[2])
                              for s in read_sites]
                # This makes no sense
                self.topsite = self.sites[0]
            else:
                self.sites = [Utr.UtrSite(_sitelist, "None", None, None)]
                self.topsite = Utr.UtrSite(_sitelist, "None", None, None)


class ReporterVariant(Sequence):   
    def __init__(self, seq):
        # Define the left and right hand sides of the barcode region:
        self.seq = seq
        self.sites = []
        self.topsite = None
    #...........................................................................
    class ReporterVariantSite(Site):
        def __init__(self, _sitelist, name, l, r):
            self.name = name
            self.seq = Mirna(_sitelist.mirna)[name]
            self.l = l
            self.r = r
    # Methods that extract information from the top site
    # Using the read and site information:
    def site_overlaps(self):
        #                                                  |**| <- This distance must be positive
        # ACACACA                                    l[0]  r[0]
        #                                            ACACACA
        # 5mer-m2.6                                           l[1]r[1]
        #                                                     ATTCC
        # GGGCAGAGTTCTACAGTCCGACGATCTTTATTTGTCACACAATACACACAAAATTCCCCGATATCGTATGCCGTCTTCTGCTTG
        #                                            -------  -----                           
        # FROM DIAGRAM l[0]
        ls = [s.l for s in self.sites][1:]
        rs = [s.r for s in self.sites][:-1]
        distance = [r - l for r, l in zip(rs, ls)]
        return([max(0, i) for i in distance])
    def all_site_names(self, overlap=True):
        # Prints the list of all sites associated with the read:
        names = [s.name for s in self.sites]
        if overlap and sum(self.site_overlaps()):
            o_i = 0
            for i, ol in enumerate(self.site_overlaps()):
                if ol > 0:
                    names[o_i] = "%s|%s|%s" %(names[o_i], ol, names[i + 1])
                    names[i + 1] = "remove"
                else:
                    o_i = i - 1
            names = [i for i in names if i != "remove"]
        return names
    def site_name(self):
        # Gives the name of the top site:
        s = self.topsite
        if s:
            return s.name
    def site_pos(self):
        # Gives the position of the top site:
        s = self.topsite
        if s:
            return s.l, s.r
    def site_flank(self, length=2):
        # Gives the flanking dinucleotides fo the top site:
        s = self.topsite
        if s:
            read = self.seq
            return "%s.%s" %(read[s.l - length:s.l], read[s.r:s.r + length])
    def get_all_sites_in_reporter_variant(self, _sitelist, ties=False):
        # Returns a list of matches with ("name", start, stop)
        if _sitelist.__class__.__name__ == "SiteList":
            var_sites = [(_site.name, k, k + len(_site))
                          for _site in _sitelist.sites
                          for k in find_all(self.seq, _site.seq)]
            var_sites = sorted(var_sites, key=lambda var_sites: var_sites[1])
            if var_sites:
                var_sites_temp = [i for i in var_sites]
                var_sites = [var_sites_temp[0]]
                if len(var_sites_temp) > 1:
                    for t in var_sites_temp:
                        l_t, r_t = t[1:]
                        add = True
                        for s in var_sites:
                            l_s, r_s = s[1:]
                            if l_t >= l_s and r_t <= r_s:
                                add = False
                            elif ((l_t <= l_s and r_t > r_s)
                                  or (l_t < l_s and r_t >= r_s)):
                                var_sites.remove(s)
                        if add == True:
                            var_sites.append(t)
                self.sites = [ReporterVariant.ReporterVariantSite(_sitelist, s[0], s[1], s[2])
                              for s in var_sites]
                site_names = [i[0] for i in var_sites]
                topsite = [var_sites[i] for i in
                             _sitelist.get_best_site_ind(site_names)]
                if not ties:
                    topsite = [topsite[0]]
                self.topsite = [ReporterVariant.ReporterVariantSite(_sitelist, s[0], s[1], s[2])
                                for s in topsite][0]
            else:
                self.sites = [ReporterVariant.ReporterVariantSite(_sitelist, "None", None, None)]
                self.topsite = ReporterVariant.ReporterVariantSite(_sitelist, "None", None, None)
        elif _sitelist.__class__.__name__ == "TwelvemerSiteListRep":
            read_sites = [(self.seq[k + _sitelist.r - 10:k + _sitelist.r + 2],
                           k + _sitelist.r - 10, k + _sitelist.r + 2) for
                          k in find_all(self.seq, _sitelist.constant) if
                          k > 9 - _sitelist.r and # This is RIGHT
                          k < len(self.seq) - _sitelist.r - 1] # THIS IS RIGHT
            if read_sites:
                self.sites = [ReporterVariant.ReporterVariantSite(_sitelist, s[0], s[1], s[2])
                              for s in read_sites]
                # This makes no sense
                self.topsite = self.sites[0]
            else:
                self.sites = [ReporterVariant.ReporterVariantSite(_sitelist, "None", None, None)]
                self.topsite = ReporterVariant.ReporterVariantSite(_sitelist, "None", None, None)
        elif _sitelist.__class__.__name__ == "SixteenmerSiteListRep":
            split = _sitelist.split
            pos = ReporterVariant.SPLIT16_POS_MAP[_sitelist.split]
            var_sites = [(self.seq[k + _sitelist.r - pos:k + _sitelist.r + 12 - pos],
                           k + _sitelist.r - pos, k + _sitelist.r +  12 - pos) for
                          k in find_all(self.seq, _sitelist.constant) if
                          k > pos - 1 - _sitelist.r and # This is RIGHT
                          k < len(self.seq) - _sitelist.r - (12 - pos) + 1] # THIS IS RIGHT
            if var_sites:
                self.sites = [ReporterVariant.ReporterVariantSite(_sitelist, s[0], s[1], s[2])
                              for s in read_sites]
                # This makes no sense
                self.topsite = self.sites[0]
            else:
                self.sites = [ReporterVariant.ReporterVariantSite(_sitelist, "None", None, None)]
                self.topsite = ReporterVariant.ReporterVariantSite(_sitelist, "None", None, None)



def add_sites_from_read(_sitelist, _read, buffer_=False):
    # Get all sites:
    _read.get_all_sites_in_read(_sitelist, buffer_=buffer_)
    names = [_site.name for _site in _read.sites]
    sequences = [_site.seq for _site in _read.sites]
    if _sitelist.__class__.__name__ == "SiteList":
        top = _read.site_name()
        if top != "None":
            _sitelist.top_counts[top][_read.barcode] += 1
            _sitelist.top_pos_counts[top][_read.barcode][_read.site_pos()[0]] += 1
            distances = ([r - l for l, r in
                         zip([i.r for i in _read.sites][:-1],
                             [i.l for i in _read.sites][1:])])
            distances = [i for i in distances if i >= 0]

            names_raw = _read.all_site_names()
            l = [i.l for i in _read.sites]
            r = [i.r for i in _read.sites]
            if len(names_raw) == 1:
                _sitelist.single_pos_counts[top][_read.barcode][_read.site_pos()[0]] += 1
            for i, i_d in enumerate(distances):
                names_raw[i] = "%s,(%s)," %(names_raw[i], i_d)
                if(names_raw[i] == "6mer,(0),6mer"):
                    _read.print_sites()
            multi = "".join(names_raw)
            # if multi == "9mer-m12.20,(17),8mer-mmT5":
            _sitelist.multi_counts[_read.barcode] += [multi]
        else:
            _sitelist.top_counts["None"][_read.barcode] += 1
            _sitelist.multi_counts[_read.barcode] += ["None"]
    elif (_sitelist.__class__.__name__ in
            ["TwelvemerSiteList", "SixteenmerSiteList"]):
        sites = _read.sites
        if sites[0] != None:
            for site in sites:
                _sitelist.top_counts[str(site)][_read.barcode] += 1.0/len(sites)
            multi = ",".join(_read.all_site_names())
            _sitelist.multi_counts[_read.barcode] += [multi]



def add_sites_from_utr(_sitelist, _utr):
    # Get all sites:
    _utr.get_all_sites_in_utr(_sitelist)
    if _sitelist.__class__.__name__ == "SiteList":
        top = _utr.site_name()
        if top != "None":
            _sitelist.top_counts[top][_utr.barcode] += 1
            _sitelist.top_pos_counts[top][_utr.barcode][_utr.site_pos()[0]] += 1
            distances = ([r - l for l, r in
                         zip([i.r for i in _utr.sites][:-1],
                             [i.l for i in _utr.sites][1:])])
            distances = [i for i in distances if i >= 0]

            names_raw = _utr.all_site_names()
            l = [i.l for i in _utr.sites]
            r = [i.r for i in _utr.sites]
            if len(names_raw) == 1:
                _sitelist.single_pos_counts[top][_utr.barcode][_utr.site_pos()[0]] += 1
            for i, i_d in enumerate(distances):
                names_raw[i] = "%s,(%s)," %(names_raw[i], i_d)
                if(names_raw[i] == "6mer,(0),6mer"):
                    _utr.print_sites()
            multi = "".join(names_raw)
            _sitelist.multi_counts[_utr.barcode] += [multi]
        else:
            _sitelist.top_counts["None"][_utr.barcode] += 1
            _sitelist.multi_counts[_utr.barcode] += ["None"]
    elif (_sitelist.__class__.__name__ in
            ["TwelvemerSiteList", "SixteenmerSiteList"]):
        sites = _utr.sites
        if sites[0] != None:
            for site in sites:
                _sitelist.top_counts[str(site)][_utr.barcode] += 1.0/len(sites)
            multi = ",".join(_utr.all_site_names())
            _sitelist.multi_counts[_utr.barcode] += [multi]
   
def add_sites_from_variant(_sitelist, _variant):
    # Get all sites:
    _variant.get_all_sites_in_reporter_variant(_sitelist)
    if _sitelist.__class__.__name__ == "SiteList":
        top = _variant.site_name()
        if top != "None":
            _sitelist.top_counts[top] += 1
            _sitelist.top_pos_counts[top][_variant.site_pos()[0]] += 1
            distances = ([r - l for l, r in
                         zip([i.r for i in _variant.sites][:-1],
                             [i.l for i in _variant.sites][1:])])
            distances = [i for i in distances if i >= 0]

            names_raw = _variant.all_site_names()
            l = [i.l for i in _variant.sites]
            r = [i.r for i in _variant.sites]
            if len(names_raw) == 1:
                _sitelist.single_pos_counts[top][_variant.site_pos()[0]] += 1
            for i, i_d in enumerate(distances):
                names_raw[i] = "%s,(%s)," %(names_raw[i], i_d)
                if(names_raw[i] == "6mer,(0),6mer"):
                    _variant.print_sites()
            multi = "".join(names_raw)
            _sitelist.multi_counts += [multi]
        else:
            _sitelist.top_counts["None"] += 1
            _sitelist.multi_counts += ["None"]
    elif (_sitelist.__class__.__name__ in
            ["TwelvemerSiteList", "SixteenmerSiteList"]):
        sites = _variant.sites
        if sites[0] != None:
            for site in sites:
                _sitelist.top_counts[str(site)] += 1.0/len(sites)
            multi = ",".join(_variant.all_site_names())
            _sitelist.multi_counts += [multi]




def add_flanks_from_read(_sitelist, _read, buffer_=False):
    # Get all sites:
    _read.get_all_sites_in_read(_sitelist, buffer_=buffer_)
    top = _read.site_name()
    seq = _read.topsite
    if top != "None":
        flank = _read.site_flank()
        _sitelist.flank_counts[_read.barcode][top][flank] += 1




#...............................................................................
class RawRead(Read):
    def __init__(self, four_line):
        self.header1, self.seq, self.header2, self.qc = [i for i in four_line]
        self.N = "N" not in self.seq
        self.bin_score = self.header1[-1] == 1
        self.B = "B" not in self.qc
        self.multi = self.header1.split("#")[1].split("/")[0]
        self.sites = []
    def check_multiplex_barcode(self, multi_seq):
        return(self.multi == multi_seq)
    def check_library_barcode(self, lib_len, lib_bcs):
        if type(lib_bcs) == "list":
            return(self.seq[lib_len:lib_len + 3] in lib_bcs)
        else:
            return(self.seq[lib_len:lib_len + 3] == lib_bcs)
    def check_if_passes(self, multi_seq, lib_len, lib_bcs):
        return([self.bin_score, self.N, self.B,
                self.check_library_barcode(lib_len, lib_bcs),
                self.check_multiplex_barcode(multi_seq)])

###########END CLASSES##########################################################











def RNAplfold(read, win, file_num, logprob=False):
    # This should give the probability of the region in the read complementary
    # to position "mir_start" to "mir_stop" in the miRNA sequence.
    time1 = time.time()
    os.chdir("AnalyzeStructure/RNAplfold_temp")
    time1 = time.time()
    name = randomword(20) + str(file_num)
    # Make the fasta file:
    time1 = time.time()
    with open(name + ".fa", "wb") as fa:
        fa.write(">%s\n%s" %(name, read))
    # call RNAplfold:
    l, w, u = [len(read), len(read), win]
    if logprob:
        o, pl_ext = ["-O ", "_openen"]
    else:
        o, pl_ext = ["", "_lunp"]
    time1 = time.time()
    call = 'RNAplfold -L %s -W %s -u %s %s< %s.fa' %(l, w, u, o, name)
    subprocess.call([call], shell=True, stdout = subprocess.PIPE)
    data = pd.read_csv(name + pl_ext, sep = '\t', header = 1).set_index(' #i$')
    data = data.iloc[:, :u]
    for ext in [".fa", "_dp.ps", pl_ext]:
        os.remove(name + ext)
    os.chdir("../..")
    return data

def extract_raw_read(file):
    """Retrieves the next four lines of the file.

    Args:
        file: The file object being read.

    Returns:
        A list of four reads, corresponding to the header, the sequence,
        the header, and the quality scores corresponding to each base.
        As a consequence, the read file is updated to be four lines ahead.
    """
    out = [file.readline().strip() for i in range(4)]
    if out[3]:
        return RawRead(out)

def check_cd1as(mirna, start, stop, species="mmu"):
    with open("general/Sequences/%s_Cdr1as.fa" %(species), "r+") as file_open:
        cdr1as_seq = "".join([i.strip() for i in list(file_open)[1:]])

    _mirna = Mirna(mirna)
    mir_rc = _mirna.rev_com()
    mir_rc = mir_rc[:-1] + "A"
    mir_site = mir_rc[len(_mirna.seq) - stop:
                              len(_mirna.seq) - start + 1]
    print((_mirna.seq[::-1] + "\t %s" %(mirna)))
    print((str.lower(mir_rc)[:len(_mirna.rev_com()) - stop] + mir_site + "\t nt %s-%s" %(start, stop)))
    match_pos = find_all(cdr1as_seq, mir_site)
    for match_pos_i in match_pos:
        print((cdr1as_seq[match_pos_i - len(mir_rc) + len(mir_site): match_pos_i + len(mir_site) + 10] + "\t nt%s-%s" %(match_pos_i + 1, match_pos_i + len(mir_site))))
        print((" "*(len(mir_rc) - len(mir_site)) + mir_site))

def get_threep_site_deltaG(mirna, start, stop, mm=None, wobble=True, just_complement=False, mir_sub=False):
    _mirna = Mirna(mirna)
    threep_seq = _mirna.seq[8:]
    # Convert standard indeces to python indeces
    start_p = start - 9
    stop_p = stop - 8
    # Get the sequence of the miRNA three-p site
    target_seq = threep_seq[start_p:stop_p]
    if type(mm) == str:
        # Split up the mismatch string into the nucleotide and the position.
        mm_nuc = mm[0]
        start_mm_p = int(mm[1:]) - start
        # Convert the target sequence string into a list, to replace the wt
        # nucleotide with the reverse complement of the mismatch, such that the
        # reverse complement function will generate the intended target site
        # nucleotide.
        target_site_list = list(target_seq)
        target_site_list[start_mm_p] = get_rc(mm_nuc, rna=True)
        target_seq = "".join(target_site_list)
    # If not starting at position 9, get the three left-hand mismatch
    # nucleotides.
    if start_p != 0 and not just_complement:
        mm_5p = [i for i in NTS if i != threep_seq[start_p - 1]]
        print_ind = start_p - 1
    else:
        mm_5p = [""]
        print_ind = start_p
    # If not starting at the last nucleotide of the miRNA, get the three right-
    # hand mismatch nucleotides.
    if stop_p != len(threep_seq) and not just_complement:
        mm_3p = [i for i in NTS if i != threep_seq[stop_p]]
    else:
        mm_3p = [""]
    # Make the full list of all of the target sites with appended 5prime and
    # threeprime nucleotides.
    target_sites = [get_rc("%s%s%s" %(mm_5pi, target_seq, mm_3pi), rna=True)
                    for mm_5pi in mm_5p for mm_3pi in mm_3p]
    if mir_sub:
        threep_seq = threep_seq[start_p:stop_p]
    # print(mirna + " " + "_"*(len(threep_seq) - len(mirna) - 1))
    # print(threep_seq)
    # for target_site in target_sites:
    #     print(" "*max(start - 9 - int(not just_complement), 0) + target_site[::-1])
    # # Calculate the delta G of each of them and take the average to get the
    # final output of the function.
    delta_g = np.mean([get_mfe(threep_seq, target_site, wobble=wobble)
                              for target_site in target_sites])
    return(delta_g)


def main():
    # dG1 = get_threep_site_deltaG("let-7a-21nt", 11, 20)
    # dG2 = get_threep_site_deltaG("let-7a-21nt", 11, 20, wobble=False)
    # dG3 = get_threep_site_deltaG("let-7a-21nt", 11, 20, mm="U11")
    # dG4 = get_threep_site_deltaG("let-7a-21nt", 11, 20, mm="U11", just_complement=True)
    # dG2 = get_threep_site_deltaG("let-7a-21nt", 11, 21, mm="G11")
    # dG3 = get_threep_site_deltaG("let-7a-21nt", 11, 21, mm="G12")
    # dG4 = get_threep_site_deltaG("let-7a-21nt", 11, 21, mm="G13")
    # print("done")
    # print(dG1)
    # print(dG2)
    # print(get_mfe("GGUUG", "CAACC"))
    # print(get_mfe("GGUUG", "CAGCC"))

    # mirna = "let-7a-21nt"
    # start = 9
    # stop = 13
    # # print(LCSuffOld("hello", "hello"))
    # # print(LCSubStringOld("jello", "hellok"))
    # check1 = get_threep_site_deltaG(mirna, 9, 13)
    # check2 = get_threep_site_deltaG(mirna, 18, 21)
    # check3 = get_threep_site_deltaG(mirna, 16, 20)

    # print([check1, check2, check3])

    # _mirna.get_best_name(["TTCACCGCGT", 0, 0, 0, 0, 0])

    # _mirna.get_best_name(["TTCACCGCGTG", 0, 0, 0, 0, 0])

    # temp_list = KmerList(3, 5, 37, fast=False)
    # temp_list.read("miR-1", "equilibrium_2_tp", "I", "5", "1",
    #                temp_path="1562177926.7916143_kpdclpumwpdaaidnsysj.txt")

    # print(temp_list.kmers["TGT"])
    # print(temp_list.pos_mean["TGT"])
    # print(temp_list.pos_sd["TGT"])
    # print(temp_list.sum["TGT"])


    # _mirna = Mirna("miR-124")
    # kmer_temp = ["AAGTGCAATTA", 1, 0, 0, 0, 0]
    # sub = True
    # _mirna.get_best_name(kmer_temp, suboptimal=sub)
    return
    # _mirna.get_kmer_match("AGTCTTCCA")
    # _mirna.get_kmer_match("AGTACTTCCA")
    # _mirna.get_kmer_match("AGTCTCCCA")

    # _mirna.get_kmer_bulge_match("GTCTTCCA")
    # _mirna.get_kmer_bulge_match("GTCTATCCA")
    # # _mirna.get_kmer_bulge_match("AGTCTTCCA")
    # _mirna.get_kmer_bulge_match("AGTACTTCCA")

    # # _mirna.get_kmer_bulge_match("AGTCTCCCA")

    # _mirna.get_kmer_1_mismatch("AGTCTTCCA")
    # _mirna.get_kmer_1_mismatch("AGTACTTCCA")
    # _mirna.get_kmer_1_mismatch("AGTCTCCCA")
    # print(_mirna["8mer"])
    # # _mirna.get_best_name([_mirna["8mer-bA8"], 20], suboptimal=True)
    # _mirna = Mirna("miR-7-23nt")
    # print(" "*(len(_mirna) - 8) + "87654321")
    # print(_mirna.rev())
    # print(_mirna["8mer"])
    # _mirna.get_best_name([_mirna["8mer-mmA7bG7"], 20], suboptimal=True)
    # _mirna.get_best_name([_mirna["7mer-m8mmA7bG7"], 20], suboptimal=True)
    # # return
    # # print(_mirna["7mer-A1mmT5"])
    # for mirna in ["miR-1", "let-7a", "miR-155", "miR-124", "lsy-6", "miR-7-23nt"]:
    #     print(mirna + "_"*(40 - len(mirna)))
    #     _mirna = Mirna(mirna)
    #     _sitelist = SiteList(_mirna, "paper")
    #     sitenames = [i.name for i in _sitelist.sites]
    #     siteseqs = [i.seq for i in _sitelist.sites]
    #     site_test = zip(sitenames, siteseqs)
    #     # print(site_test)
    #     for site_tuple in site_test:
    #         # print(site_tuple[0])
    #         seq_test = "AA" + site_tuple[1] + "AA"
    #         best = _mirna.get_best_name([site_tuple[1], 1], print_sites=False)
    #         print(site_tuple[0] + " "*(20 - len(site_tuple[0])) + best[0] + " "*(20 - len(best[0])) + str(site_tuple[0] == best[0]))
    # print(_mirna)

    # _mirna.get_best_name(["CTTCCTAT", 1, 1, 1, 1, 1])
    # print(get_all_mismatches("ACACACA"))


################################################################################
if __name__ == "__main__":
    main()

