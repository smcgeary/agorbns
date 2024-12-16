import imp # Used to import general.py
import time # Used for time
import os
import subprocess
import numpy as np
import pandas as pd
import math
import sys
import itertools as it
import RNA
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
import random, string


class Mirna:
	SEQ_NAME_MAP = {
	    "miR-1"      : "UGGAAUGUAAAGAAGUAUGUAU",
	    "miR-155"    : "UUAAUGCUAAUCGUGAUAGGGGU",
	    "miR-124"    : "UAAGGCACGCGGUGAAUGCCAA",
	    "let-7a"     : "UGAGGUAGUAGGUUGUAUAGUU",
	    "lsy-6"      : "UUUUGUAUGAGACGCAUUUCGA",
	    "miR-7"      : "UGGAAGACUAGUGAUUUUGUUGUUU",
	    "miR-7-25nt" : "UGGAAGACUAGUGAUUUUGUUGUU",
	}

	def __init__(self, name):
		self.name = name
		self.seq = self.SEQ_NAME_MAP[name]


class ReadSeq:
	FLANK_5P = "GGGCAGAGTTCTACAGTCCGACGATC"
	FLANK_3P = "TGTTCGTATGCCGTCTTCTGCTTG"
	def __init__(self, seq):
		self.seq = seq
	def fullseq(self):
		return "".join([self.FLANK_5P, self.seq, self.FLANK_3P])
	


	# FLANK_5p = "GGGCAGAGTTCTACAGTCCGACGATC"
	# # FLANK_3p = "TGTTCGTATGCCGTCTTCTGCTTG"




readtemp = ReadSeq("ATACACTCTAGTTTCCCACTACTCCA")

print(readtemp.seq)
print(readtemp.fullseq())
Mirna = Mirna("miR-1")

filename = "/lab/solexa_bartel/mcgeary/AgoRBNS/miR-1/equilibrium/reads/I.txt"

with filename.open("r") as file:
	for i in range(10):
		print(file.readline)
