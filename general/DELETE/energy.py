
nn_dH = {
    "init"  :   3.61,
    "AA.UU" :  -6.82,
    "AU.AU" :  -9.38,
    "UA.UA" :  -7.69,
    "CU.AG" : -10.48,
    "CA.UG" : -10.44,
    "GU.AC" : -11.40,
    "GA.UC" : -12.44,
    "CG.CG" : -10.64,
    "GG.CC" : -13.39,
    "GC.GC" : -14.88,
    "tAU"   :   3.72,
    "sym"   :   0.00,
}

nn_dS = {
    "init"  :  -1.5,
    "AA.UU" : -19.0,
    "AU.AU" : -26.7,
    "UA.UA" : -20.5,
    "CU.AG" : -27.1,
    "CA.UG" : -26.9,
    "GU.AC" : -29.5,
    "GA.UC" : -32.5,
    "CG.CG" : -26.7,
    "GG.CC" : -32.7,
    "GC.GC" : -36.9,
    "tAU"   :  10.5,
    "sym"   :  -1.4,
}

nn_dG = {
    "init"  :  4.09,
    "AA" : -0.93,
    "AC" : -2.24,
    "AG" : -2.08,
    "AH": -0.55,
    "AU" : -1.10,
    "AV" : -1.36,

    "CA" : -2.11,
    "CC" : -3.26,
    "CG" : -2.36,
    "CH" : -1.41,
    "CU" : -2.08,
    "CV" : -2.11,

    "GA" : -2.35,
    "GC" : -3.42,
    "GG" : -3.26,
    "GH" : -1.53,
    "GU" : -2.24,
    "GV" : -2.51,

    "HA" : -1.27,
    "HC" : -2.51,
    "HG" : -2.11,
    "HH" : -0.5,
    "HU" : -1.36,
    "HV" : +1.29,

    "UA" : -1.33,
    "UC" : -2.35,
    "UG" : -2.11,
    "UH" : -1.00,
    "UU" : -0.93,
    "UV" : -1.27,

    "VA" : -1.00,
    "VC" : -1.53,
    "VG" : -1.41,
    "VH" : +0.3,
    "VU" : -0.55,
    "VV" : -0.5,

    "GHVC" : -4.12,

    "AU_end" : 0.45,
    "sym" : 0.43,

}

# Order of dictionaries is penultimate pair, then 5' nucleotide, then 3' nucleotide

term_mm_dG = {
    "A" : {"A" : {"A" : -0.8, "C" : -1.0, "G" : -0.8},
           "C" : {"A" : -0.6, "C" : -0.7, "U" : -0.7},
           "G" : {"A" : -0.8, "G" : -0.8},
           "U" : {"C" : -0.8, "U" : -0.8},
           },
    "C" : {"A" : {"A" : -1.5, "C" : -1.5, "G" : -1.4},
           "C" : {"A" : -1.0, "C" : -1.1, "U" : -0.8},
           "G" : {"A" : -1.4, "G" : -1.6},
           "U" : {"C" : -1.4, "U" : -1.2},
           },
    "G" : {"A" : {"A" : -1.1, "C" : -1.5, "G" : -1.3},
           "C" : {"A" : -1.1, "C" : -0.7, "U" : -0.5},
           "G" : {"A" : -1.6, "G" : -1.4},
           "U" : {"C" : -1.0, "U" : -0.7},
           },
    "H" : {"A" : {"A" : -0.3, "C" : -1.0, "G" : -0.8},
           "C" : {"A" : -0.6, "C" : -0.7, "U" : -0.7},
           "G" : {"A" : -0.6, "G" : -0.8},
           "U" : {"C" : -0.8, "U" : -0.6},
           },
    "U" : {"A" : {"A" : -1.0, "C" : -0.8, "G" : -1.1},
           "C" : {"A" : -0.7, "C" : -0.6, "U" : -0.5},
           "G" : {"A" : -1.1, "G" : -1.2},
           "U" : {"C" : -0.6, "U" : -0.5},
           },
    "V" : {"A" : {"A" : -1.0, "C" : -0.8, "G" : -1.1},
           "C" : {"A" : -0.7, "C" : -0.6, "U" : -0.5},
           "G" : {"A" : -0.5, "G" : -0.8},
           "U" : {"C" : -0.6, "U" : -0.5},
           },

}

with open("general/temp_target_site.fa", "rb") as file_in:
    test_seq = "".join([i.strip() for i in file_in.readlines()][1:])

trim = test_seq.find("CAAGATTGTGCCACCACACTCCA")
test_seq = test_seq[trim - 40: trim + 60]




ind_start = test_seq.find("ACACTCC")
ind_stop = ind_start + 8




intern_mm_dG = {i : {j : {k : { l : None for l in nts} for k in nts} for j in wob_nts} for i in wob_nts}
line = mm_file.readline()
for f5p in ["A", "C", "G", "U", "H", "V"]:
    for x in nts:
        if line:
            fields = map(float,line.strip().split("\t"))
            index=0
            for f3p in ["A", "C", "G", "U", "H", "V"]:
                for y in nts:
                    intern_mm_dG[f5p][f3p][x][y] = fields[index]
                    index+=1
        line = mm_file.readline()

mm_file = open("general/mismatch_table.txt")

