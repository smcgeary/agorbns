# Lists used in all AGO-RBNS analysis

# kRomanNumerals <- c("i", "ii", "iii", "iv", "v", "vi", "vii", "viii", "ix", "x",
#                     "xi", "xii", "xiii", "xiv", "xv", "xvi", "xvii", "xviii",
#                     "xix", "xx", "xxi", "xxii", "xxiii", "xxiv", "xxv", "xxvi",
#                     "xxvii", "xxviii", "xxix", "xxx")

# kNextLetter <- c("B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M",
#                  "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y",
#                  "Z")
# names(kNextLetter) <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K",
#                         "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V",
#                         "W", "X", "Y")



kMirnas <- c("miR-1",
             "let-7a",
             "miR-155",
             "miR-124",
             "lsy-6",
             "miR-7-23nt",
             "let-7a-21nt",
             "miR-7")

kMirnasEquil <- kMirnas[1:6]

# kMirnasThy <- c("miR-9",
#                 "miR-21",
#                 "miR-31",
#                 "miR-27a",
#                 "miR-27b",
#                 "miR-122",
#                 "miR-129",
#                 "miR-145",
#                 "miR-503")


# kCenteredSites <- c("11mer-m3.13",
#                     "12mer-m3.14",
#                     "11mer-m4.14",
#                     "12mer-m4.15")

kSeedSites <- c("8mer",
                "7mer-m8",
                "7mer-A1",
                "6mer",
                "6mer-m8",
                "6mer-A1")

kCanonicalSites <- c("8mer",
                     "7mer-m8",
                     "7mer-A1",
                     "6mer")

kPositionalSites <- c("8mer",
                      "6mer",
                      "6mer-m8",
                      "11mer-m3.13",
                      "11mer-m4.14",
                      "11mer-m5.15",
                      "11mer-m6.16",
                      "11mer-m7.17",
                      "11mer-m8.18",
                      "11mer-m9.19",
                      "11mer-m10.20",
                      "11mer-m11.21",
                      "11mer-m12.22",
                      "11mer-m13.23",
                      "None")

## Used for colorization of sites based on broad overall categories ###########
k7mer8merCanSites <- c("8mer", "7mer-m8", "7mer-A1")

k6merCanSites <- c("6mer", "6mer-A1", "6mer-m8")

kE6merSites <- c("8mer-bT(7.8)", "8mer-bA8", "7mer-m8bT(7.8)", "8mer-w8",
                 "7mer-m8w8", "8mer-bT(3.5)", "8mer-bA(6.7)", "8mer-bG(6.7)",
                 "8mer-bG7", "8mer-mmG7bG7", "8mer-mmA7bG7", "8mer-mmC7bG7",
                 "7mer-m8bA8", "AA-8mer-mmC7")

k3PSites <- c("11mer-m9.19", "11mer-m10.20", "11mer-m11.21", "11mer-m12.22",
              "11mer-m13.23", "10mer-m9.18", "10mer-m10.19", "10mer-m11.20",
              "10mer-m12.21", "10mer-m13.22", "10mer-m14.23", "9mer-m9.17",
              "9mer-m10.18", "9mer-m11.19", "9mer-m12.20", "9mer-m13.21",
              "9mer-m14.22", "9mer-m15.23",
              "11mer-m13.23w17", "11mer-m12.22w13", "11mer-m13.23w13",
              "11mer-m12.22w20", "11mer-m12.22w17", "11mer-m12.22w14",
              "10mer-m13.22w13", "10mer-m13.22w17", "10mer-m13.22w14",
              "10mer-m13.22w14", "10mer-m13.22w20", "11mer-m9.19w9",
              "11mer-m9.19w18", "10mer-m9.18w9", "10mer-m7.16",
              "10mer-m10.19w11")

kWeirdSites <- c("ACACACA", "GCTTCCGC", "CACACAC",  # miR-1
                 "GCACTTTA", "CTTCCT", "TCCTGCGC", # let-7a
                 "AACGAGGA", "CTCAGCA", "AATAAAG", # miR-155
                 "ACGACAA", "GACAAC", "AACGAGG",
                 "GACAACA", "ACAACTC",
                 "TCCGCCACA", "TCCGCCAC", "CCGCCACA", # miR-124
                 "CCGCCAC", "TCACCCGC", "CTCTGCCC", "ATCGGCGG",
                 "CCTCCGCC", "CACAGAA", "CTCTGCC", # lsy-6
                 "CGCTTCCG", "CCTCCGCA", "TACGGCTA", # miR-7
                 "AACGACTA", "GCACCA", "CGTTTCCG",
                 "CGCTTTCG", "CCGCTG", "GGAGGGAG",
                 "GCACAC",
                 "CTTCCG",
                 # miR-1 papercutoff
                 "CCACACACA", "CGCTTCCGC",
                 # miR-1 papercutoff2
                 "CACGCTTCC", "AGCTTCACA", "CTTCCGCTG", "CCAGCACGC",
                 "CGCTTCACA", "TACCCATGA", "TGCTTCACA", "GTCTTCACA",
                 # let-7a papercutoff
                 "TGCACTTTA", "AGCACTTTA", "CGCACTTTA",
                 # miR-155 papercutoff
                 "AACGAGGAA", "TAACGAGGA", "AACTCAGCA", "AAATAAAGA",
                 "AACGAGGAG", "AACGAGGAT", "TACTCAGCA", "AGACGACAA",
                 "ATAATAAAG", "AATAAAGAA", "CACTCAGCA", "ATAACGAGG",
                 "AAAATAAAG", "AAACGACAA", "CAATAAAGA", "ATGACAACA",
                 "ACGACAACA", "CTCAGCAAT", "AATAAAGAT", "AATAAAGAC",
                 "TCGACAACA", "CAACTCAGC", "ACGACAACT", "CGACAACTC",
                 # miR-155 papersubmission
                 "AATAAAGA", "ACTCAGCA", "ACGACAAC", "CGACAACA", "ATGACAAC",
                 "GACAACA",
                 # miR-7 papercutoff
                 "CCGCACCAC", "CGCTTCCGT", "TTCCGCTGC", "CGCACCACA",
                 "TCCGCACCA", "ATACGGCTA", "CGTTTCCGC", "AACGACTAA",
                 "TGCACCACA", "CCTCCGCAC", "TCCGCTGCG", "CGCTTTCGC",
                 "GGAGGGAGG", "CACAGCGGC", "CAAGCACAC", "GCACACAGC",
                 # miR-7- papercutoff2
                 "GCTTCCGCT", "TACACGCCA", "AACACGCCA",
                 # miR-7 papersubmission4
                 "CGCTTCCGCTG",
                 "CCACGCCTT", "CCTAAGTGC", "CTTAAGTGC") # miR-124 tp

# kWeirdPositionalSites <- c("TCAA-6mer-m8xG5", "CTAA-5mer-m8", "TCAA-6mer-m8w5")

kSiteCategories <- c(rep("Canonical7-8", length(k7mer8merCanSites)),
                     rep("Canonical6", length(k6merCanSites)),
                     rep("E6merSites", length(kE6merSites)),
                     rep("ThreePSites", length(k3PSites)),
                     rep("WeirdSites", length(kWeirdSites)), "None")
names(kSiteCategories) <- c(k7mer8merCanSites, k6merCanSites, kE6merSites,
                            k3PSites, kWeirdSites, "None")


# kKmerSiteLists <- c("8mers",
#                     "9mers",
#                     "10mers",
#                     "11mers",
#                     "12mers",
#                     "16mers")

# # List of the experiments for which the sitelist "programmed" confers differ
# # meaning and thus there is different behavior of the code. This applies to
# # `MakeSiteCountTable.py` and `FitSiteKds.R`.
# kExpsThreeP <- c("equil_c_nb", "equil_s_nb", "equil_sc_nb", "equil_c2_nb",
#                  "equil_c_alt_nb", "equil_c2_alt_nb", "equil_sc_alt_nb")


                                 #         1         2
                                 #1234567890123456789012345
kMirnaSeqs <- c(`miR-1`         ="UGGAAUGUAAAGAAGUAUGUAU",
                `let-7a-21nt`   ="UGAGGUAGUAGGUUGUAUAGU",
                `let-7a`        ="UGAGGUAGUAGGUUGUAUAGUU",
                `miR-155`       ="UUAAUGCUAAUCGUGAUAGGGGU",
                `miR-124`       ="UAAGGCACGCGGUGAAUGCCAA",
                `lsy-6`         ="UUUUGUAUGAGACGCAUUUCGA",
                `miR-7-22nt`    ="UGGAAGACUAGUGAUUUUGUUG",
                `miR-7-23nt`    ="UGGAAGACUAGUGAUUUUGUUGU",
                `miR-7-24nt`    ="UGGAAGACUAGUGAUUUUGUUGUU",
                # `miR-7-25nt`    ="UGGAAGACUAGUGAUUUUGUUGUUU",
                `let-7a_plus1`  ="UGAGGUAGUAUGGUUGUAUAG",
                `let-7a_minus1` ="UGAGGUAGUGGUUGUAUAGUA",
                `let-7a_miR-155`="UGAGGUAGAAUCGUGAUAGGGGU",
                `miR-155_let-7a`="UUAAUGCUUAGGUUGUAUAGU"
                # `miR-9`         ="UCUUUGGUUAUCUAGCUGUAUGA",
                # `miR-21`        ="UAGCUUAUCAGACUGAUGUUGA",
                # `miR-31`        ="AGGCAAGAUGCUGGCAUAGCU",
                # `miR-27a`       ="UUCACAGUGGCUAAGUUCCGC",
                # `miR-27b`       ="UUCACAGUGGCUAAGUUCUGC",
                # `miR-122`       ="UGGAGUGUGACAAUGGUGUUUG",
                # `miR-129`       ="CUUUUUGCGGUCUGGGCUUGC",
                # `miR-145`       ="GUCCAGUUUUCCCAGGAAUCCCU",
                )

# kMirnaStarSeqs <- c(`miR-1`      ="ACAUACUUCUUUAUAUGCCCAUA",
#                     `let-7a-21nt`="CUAUACAAUCUACUGUCUUUC",
#                     `let-7a`     ="CUAUACAAUCUACUGUCUUUC",
#                     `miR-155`    ="CUCCUACAUAUUAGCAUUAACA",
#                     `miR-124`    ="CGUGUUCACAGCGGACCUUGAU",
#                     `lsy-6`      ="GAAAUGCGUCUAGUAUCAAAAUC",
#                     `miR-7-23nt` ="CAACAAAUCACAGUCUGCCAUA")


# kCaptureOligoSeqs <- c(`miR-1`       ="UCUUCCUCCGCACCACACACAUUCCAACCUUACACAC",
#                        `let-7a`      ="UCUUCCUGCGCACCAAGCCUACCUCAACUUUACACAC",
#                        `let-7a-21nt` ="UCUUCCUCCGCACCACACCUACCUCAACCUUACACAC", 
#                        `miR-155`     ="AAAUAAAGACGACAACUCAGCAUUAAACCUUACACAC",
#                        `miR-124`     ="UCUUCCUCCGCCACAGAAGUGCCUUAACCUUACACAC",
#                        `lsy-6`       ="UCUUCCUCCGCCACAGAAAUACAAAAACCUUACACAC",
#                        `miR-7-23nt`  ="UCUUCCUCCGCACCACACGUCUUCCAACCUUACACAC")

# kCompetitorOligoSeqs <- c(`miR-1`       ="AAGGTTGGAATGTGTGTGGTGCGGAGGAAGA",
#                           `let-7a`      ="AAAGTTGAGGTAGGCTTGGTGCGCAGGAAGA",
#                           `let-7a-21nt` ="AAGGTTGAGGTAGGTGTGGTGCGGAGGAAGA",
#                           `miR-155`     ="AAGGTTTAATGCTGAGTTGTCGTCTTTATTT",
#                           `miR-124`     ="AAGGTTAAGGCACTTCTGTGGCGGAGGAAGA",
#                           `lsy-6`       ="AAGGTTTTTGTATTTCTGTGGCGGAGGAAGA",
#                           `miR-7-23nt`  ="AAGGTTGGAAGACGTGTGGTGCGGAGGAAGA")

# kLibrarySeqs <- c(`miR-1_pilot`  ="GGGUUCAGAGUUCUACAGUCCGACGAUCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNUAUGCCGUCUUCUGCUUG",
#                   `miR-155_pilot`="GGGUUCAGAGUUCUACAGUCCGACGAUCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNUAUGCCGUCUUCUGCUUG",
#                   `miR-1`        ="GGGCAGAGUUCUACAGUCCGACGAUCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNUCGUAUGCCGUCUUCUGCUUG",
#                   `let-7a`       ="GGGCAGAGUUCUACAGUCCGACGAUCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNUGUUCGUAUGCCGUCUUCUGCUUG",
#                   `miR-155`      ="GGGCAGAGUUCUACAGUCCGACGAUCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNUGUUCGUAUGCCGUCUUCUGCUUG",
#                   `miR-124`      ="GGGCAGAGUUCUACAGUCCGACGAUCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNUGUUCGUAUGCCGUCUUCUGCUUG",
#                   `lsy-6`        ="GGGCAGAGUUCUACAGUCCGACGAUCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNUGUUCGUAUGCCGUCUUCUGCUUG",
#                   `miR-7-22nt`   ="GGGCAGAGUUCUACAGUCCGACGAUCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNUGUUCGUAUGCCGCUGUGUGCUUG",
#                   `miR-7-23nt`   ="GGGCAGAGUUCUACAGUCCGACGAUCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNUGUUCGUAUGCCGCUGUGUGCUUG",
#                   `miR-7-24nt`   ="GGGCAGAGUUCUACAGUCCGACGAUCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNUGUUCGUAUGCCGCUGUGUGCUUG",
#                   `miR-7-25nt`   ="GGGCAGAGUUCUACAGUCCGACGAUCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNUGUUCGUAUGCCGCUGUGUGCUUG")

# kLibrarySeqs5P <- sapply(kLibrarySeqs, function(lib) {
#                            strsplit(lib, split="N")[[1]][1]
#                          })

# kLibrarySeqs3P <- sapply(kLibrarySeqs, function(lib) {
#                            rev(strsplit(lib, split="N")[[1]])[1]
#                          })

# SequenceList <- list(`mirna`     =kMirnaSeqs,
#                      `star`      =kMirnaStarSeqs,
#                      `capture`   =kCaptureOligoSeqs,
#                      `competitor`=kCompetitorOligoSeqs,
#                      `lib5p`     =kLibrarySeqs5P,
#                      `lib3p`     =kLibrarySeqs3P)

# SequenceListLabels <- c(`mirna`     ="miRNA sequence",
#                         `star`      ="miRNA* sequence",
#                         `capture`   ="Capture oligo",
#                         `competitor`="Competitor oligo",
#                         `lib5p`     ="Library 5-prime sequence",
#                         `lib3p`     ="Library 3-prime sequence")


kDinucs <- c("AA", "AT", "TA", "TT",
             "AC", "CA", "TC", "CT",
             "AG", "GA", "TG", "GT",
             "CC", "GC", "CG", "GG")

# kNucs  <- c("A", "U", "C", "G")

kDNucs <- c("A", "T", "C", "G")

# kComplements <- c(`A`="T", `C`="G", `G`="C", `T`="A")

kFlanks <- paste0(rep(sort(kDinucs), each=length(kDinucs)), ".", rep(sort(kDinucs), length(kDinucs)))


# # Make the vector of the massively parallel reporter assay sequences.
# reporter_assay_seqs <- t(
#   matrix(unlist(read.delim("ReporterScreen/final_order_twist.fa",
#                            header=FALSE)), nrow=2)
# )
# reporter_assay_seqs[, 1] <- gsub("> ", reporter_assay_seqs[, 1], replace="")

# kMPRA <- reporter_assay_seqs[, 2]
# names(kMPRA) <- reporter_assay_seqs[, 1]


