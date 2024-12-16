# source("general/general.R")
# source("general/ThrPFunctions_temp.R")
# source("general/ThrPFunctions.R")
# source("Repression/repression_workspace.R")


mm_inds <- grep("^8mer-mm[ACTG][2-7]_Kd$", rownames(kds_global), perl=TRUE)
ref_kd <- GeoMean(kds_global[mm_inds, 2])



kds_reduced <- kds_global[!grepl("NA", rownames(kds_global)), ]
kds_reduced <- kds_reduced[grepl("Comp", rownames(kds_reduced)), ]
kds_reduced <- kds_reduced[grepl("mer-m", rownames(kds_reduced)), ]
kds_reduced <- kds_reduced[grepl("\\.", rownames(kds_reduced), perl=TRUE), ]


print(max(kds_reduced[, 2])/min(kds_reduced[, 2]))
print(ref_kd/min(kds_reduced[, 2]))

break


MakeGEOProcessedReporterCountsTable <- function() {
  mirnas <- c("let-7a-21nt", "let-7a-21nt", "miR-1", "miR-155", "let-7a_plus1",
              "let-7a_minus1", "let-7a_miR-155", "miR-155_let-7a")
  experiments <- c("equil_c2_nb", "equil_c_nb", "equil_c_nb", "equil_sc_nb",
                   "equil_c_nb", "equil_c_nb", "equil_c_nb", "equil_c_nb")


  mirnas <- rep(c("miR-1", "miR-1", "let-7a", "let-7a-21nt"), 2)

  conditions <- rep(c("no_duplex_parallel", rep("duplex_parallel", 3)), 2)

  reps <- rep(c(1, 2), each=4)

  table_conditions <- cbind(cbind(mirnas, conditions), reps)

  reporter_cpm_df <- do.call("cbind", apply(table_conditions, 1, function(row) {
    mirna <- row[1]
    condition <- row[2]
    rep <- as.integer(row[3])
    counts_df <- SubfunctionCall(GetThreePrimeReporterCounts,
                                 experiment="twist_reporter_assay_3p_2_tp")
    if (condition == "no_duplex_parallel" & rep == 1) {
      return(counts_df)
    } else {
      return(counts_df[, 7, drop=FALSE])
    }
  }))
  # Format the table to have appropriate column names.
  mirnas[1] <- "mock"
  mirnas[5] <- "mock"
  colnames(reporter_cpm_df)[7:14] <- sprintf("cpm_%s_rep%s", mirnas, reps)
  # Write the output file.
  out_path <- "/lab/solexa_bartel/mcgeary/GEO_McGearyBisaria2021_3/processed_files/"
  write.table(MakeAltMatrix(reporter_cpm_df),
              file=sprintf("%smpra_cpm.txt", out_path), sep="\t", quote=FALSE,
              row.names=FALSE, col.names=FALSE)
}


MakeGEOProcessedReporterCountsTable()






break

graphics.off()


df_1 <- GetThreePrimeReporterCounts(
  "let-7a", "twist_reporter_assay_3p_2_tp", "duplex_parallel", aggregate=TRUE
)

df_2 <- GetThreePrimeReporterCounts(
  "let-7a", "twist_reporter_assay_3p_2_tp", "duplex_parallel", aggregate=TRUE,
  all_lin41_UTR=TRUE
)


print(head(df_1))
print(tail(df_1))
print(head(df_2))
print(tail(df_2))
break







test0 <- SitesXCounts("let-7a-21nt", experiment="equil_c2_nb", n_constant=3,
                      sitelist="progthrp_suppcomp")


test1 <- SitesXCounts("let-7a-21nt", experiment="equil_c2_nb", n_constant=3,
                      sitelist="progthrp_suppcomp", new2=TRUE)

test2 <- SitesXCounts("let-7a-21nt", experiment="equil_c2_nb", n_constant=3,
                      sitelist="progthrp_suppcomp", new2=TRUE, zeros_override=TRUE)

test3 <- SitesXCounts("let-7a-21nt", experiment="equil_c2_nb", n_constant=3,
                      sitelist="progthrp_suppcomp", new2=TRUE, include_zeros=TRUE)


test3_norm <- as.data.frame(t(t(test3)/colSums(test3)))
print(colSums(test3))
print(colSums(test3_norm))
A_ratios <- test3_norm[grep("9mer-m11\\.19\\|15\\|AAAA", perl=TRUE, rownames(test3_norm), value=TRUE), ]
U_ratios <- test3_norm[grep("9mer-m11\\.19\\|15\\|TTTT", perl=TRUE, rownames(test3_norm), value=TRUE), ]

print(A_ratios)
print(U_ratios)

A_ratios_1 <- A_ratios[, c("0.4", "1.265", "4", "12.65", "40")]
U_ratios_1 <- U_ratios[, c("0.4", "1.265", "4", "12.65", "40")]


A_ratios[, c("12.65")] <- rowMeans(A_ratios[, c("12.65", "12.65_2")])
U_ratios[, c("12.65")] <- rowMeans(U_ratios[, c("12.65", "12.65_2")])

A_ratios_2 <- A_ratios[, c("0.4", "1.265", "4", "12.65", "40")]
U_ratios_2 <- U_ratios[, c("0.4", "1.265", "4", "12.65", "40")]

print(sum(A_ratios_1)/sum(U_ratios_1))
print(sum(A_ratios_2)/sum(U_ratios_2))

# print(dim(test1))
# print(dim(test2))
# print(dim(test3))

# print(colSums(test1))
# print(colSums(test2))
# print(colSums(test3))




break




break



# out1 <- MakeLoopNucleotideEnrichmentMatrix(12, 6:9, 4, AUbin=TRUE)
# out2 <- MakeLoopNucleotideEnrichmentMatrix(12, 6:9, 5, AUbin=TRUE)

# print(out1)
# print(out2)
# break





# matrix_test <- matrix(0, nrow=18, ncol=6,
#                       dimnames=list(paste0(rep(11:13, each=6), "_", 0:5),
#                                     c("AU_pred", "A_pred", "U_pred",
#                                       "AU_obs",  "A_obs",  "U_obs")))

# print(matrix_test)

# row_i <- 1
# for (pos_i in 11:13) {
#   for (offset_i in 0:5) {
#     # print(pos_i)
#     # print(offset_i)
#     out_nuc <- MakeLoopNucleotideEnrichmentMatrix(pos_i, 6:9, offset_i)
#     out_nuc_mean <- apply(out_nuc[c(1, 4), ], 2, GeoMean)
#     R_pred_A <- exp(sum(log(out_nuc[1, ])))
#     R_pred_U <- exp(sum(log(c(out_nuc_mean[1], out_nuc[4, -1]))))
#     R_pred_AU <- exp(sum(log(apply(out_nuc[c(1, 4), ], 2, GeoMean))))
#     out_AUbin <- MakeLoopNucleotideEnrichmentMatrix(pos_i, 6:9, offset_i, AUbin=TRUE)
#     matrix_test[row_i, 1] <- R_pred_AU
#     matrix_test[row_i, 2] <- R_pred_A
#     matrix_test[row_i, 3] <- R_pred_U
#     matrix_test[row_i, 4:6] <- out_AUbin
#     row_i <- row_i + 1
#   }
# }

loop_lens <- rep(2:4, each=6) + rep(0:5, times=3)
cols_lens <- rep(rainbow(6), times=3)

pch_use = rep(c(0, 1, 2), each=6)




# print(matrix_test)

dev.new(xpos=20, ypos=20, height=4, width=4)
plot(matrix_test[, 1], matrix_test[, 4], log="xy", xlim=c(1, 4), ylim=c(1, 4), col=cols_lens, pch=pch_use)


dev.new(xpos=420, ypos=20, height=4, width=4)
plot(matrix_test[, 2], matrix_test[, 5], log="xy", xlim=c(1, 4), ylim=c(1, 4), col=cols_lens, pch=pch_use)

dev.new(xpos=820, ypos=20, height=4, width=4)
plot(matrix_test[, 3], matrix_test[, 6], log="xy", xlim=c(1, 4), ylim=c(1, 4), col=cols_lens, pch=pch_use)


out_trinuc <- MakeLoopNucleotideEnrichmentMatrix(11, 6:9, 4, trinuc=TRUE)

out_trinuc_sim <- MakeLoopNucleotideEnrichmentMatrix(11, 6:9, 4, trinuc_sim=TRUE)
out_trinuc_sim_dinuc <- MakeLoopNucleotideEnrichmentMatrix(11, 6:9, 4, dinuc=TRUE, trinuc_sim=TRUE)

dev.new(xpos=20, ypos=420, height=4, width=4)
plot(c(out_trinuc), c(out_trinuc_sim), log="xy", xlim=c(2^-3.5, 2^3.5), ylim=c(2^-3.5, 2^3.5))

dev.new(xpos=420, ypos=420, height=4, width=4)
plot(c(out_trinuc), c(out_trinuc_sim_dinuc), log="xy", xlim=c(2^-3.5, 2^3.5), ylim=c(2^-3.5, 2^3.5))


out_trinuc <- MakeLoopNucleotideEnrichmentMatrix(12, 6:9, 4, trinuc=TRUE)

out_trinuc_sim <- MakeLoopNucleotideEnrichmentMatrix(12, 6:9, 4, trinuc_sim=TRUE)
out_trinuc_sim_dinuc <- MakeLoopNucleotideEnrichmentMatrix(12, 6:9, 4, dinuc=TRUE, trinuc_sim=TRUE)

dev.new(xpos=820, ypos=420, height=4, width=4)
plot(c(out_trinuc), c(out_trinuc_sim), log="xy", xlim=c(2^-3.5, 2^3.5), ylim=c(2^-3.5, 2^3.5))

dev.new(xpos=1220, ypos=420, height=4, width=4)
plot(c(out_trinuc), c(out_trinuc_sim_dinuc), log="xy", xlim=c(2^-3.5, 2^3.5), ylim=c(2^-3.5, 2^3.5))


# print(out_nuc)
# print(out_AUbin)



break



break



CheckThreePSites <- function(mirna, bulge=FALSE) {
  experiments <- c(`let-7a-21nt`="equil_c2_nb",
                   `miR-1`="equil_c_nb",
                   `miR-155`="equil_sc_nb")
  lims <- list(`let-7a-21nt`=matrix(c(11, 14,
                                      11, 15,
                                      11, 16,
                                      11, 17,
                                      11, 18,
                                      9, 17,
                                      9, 18,
                                      10, 20), ncol=2, byrow=TRUE),
               `miR-1`=matrix(c(12, 15,
                                12, 16,
                                12, 17,
                                12, 18,
                                11, 18,
                                11, 19,
                                11, 20,
                                11, 21), ncol=2, byrow=TRUE),
               `miR-155`=matrix(c(13, 16,
                                  13, 17,
                                  13, 18,
                                  15, 21,
                                  15, 22,
                                  13, 21,
                                  13, 22,
                                  13, 23), ncol=2, byrow=TRUE))
  experiment <- experiments[mirna]
  lim_use <- lims[[mirna]]
  # print(lim_use)

  offset_lim_use=c(-2, 6)

  better <- c()
  sig_better <- c()
  for (row_i in seq(nrow(lim_use))) {
    start_mm <- lim_use[row_i, 1]
    stop_mm <- lim_use[row_i, 2]
    # print(sprintf("%s-%s", start_mm, stop_mm))
    kds_mmdb <- GetThreePrimeMmBDKds(mirna, experiment, start_mm, stop_mm, n_constant=3, kd_fc=TRUE, bulge=bulge, win_average=3, best_average=TRUE, offset_lim=offset_lim_use)
    kds_mmdb_ub <- GetThreePrimeMmBDKds(mirna, experiment, start_mm, stop_mm, n_constant=3, kd_fc=TRUE, bulge=bulge, win_average=3, best_average=TRUE, lower_bound=TRUE, offset_lim=offset_lim_use)
    kds_mmdb_lb <- GetThreePrimeMmBDKds(mirna, experiment, start_mm, stop_mm, n_constant=3, kd_fc=TRUE, bulge=bulge, win_average=3, best_average=TRUE, upper_bound=TRUE, offset_lim=offset_lim_use)
    wt <- kds_mmdb["wt", drop=FALSE]
    wt_ub <- kds_mmdb_ub["wt", drop=FALSE]
    mm <- kds_mmdb[-1, drop=FALSE]
    mm_lb <- kds_mmdb_lb[-1, drop=FALSE]

    mm_inds_keep <- !is.na(mm)
    mm <- mm[mm_inds_keep, drop=FALSE]
    mm_lb <- mm_lb[mm_inds_keep, drop=FALSE]
    # print("_________________________")
    # print(lim_use[row_i, ])
    # print(wt)
    # print(wt_ub)
    # print(mm)
    # print(mm_lb)
    if (sum(mm > wt) > 0) {
      # print("better")
      mismatch_sites <- sprintf("%smer-m%s.%smm%s", stop_mm - start_mm + 1, start_mm, stop_mm, names(mm)[which(mm > wt)])
      better <- c(better, mismatch_sites)
      # print(mismatch_sites)
    }
    if (sum(mm_lb > wt_ub) > 0) {
      # print("significantly better")
      mismatch_sites <- sprintf("%smer-m%s.%smm%s", stop_mm - start_mm + 1, start_mm, stop_mm, names(mm)[which(mm > wt)])
      sig_better <- c(sig_better, mismatch_sites)      
      # print(mismatch_sites)
    }
  }
  # Print which misamtches are better than teh wild type sites.
  print(mirna)
  print("better:")
  print(better)
  print("significantly better:")
  print(sig_better)
  return()
}



CheckRandThreePSites <- function(mirna, bulge=FALSE) {
  lims <- list(`let-7a`=matrix(c(12, 15,
                                 11, 15,
                                 11, 16,
                                 11, 17,
                                 10, 17), ncol=2, byrow=TRUE),
               `miR-1`=matrix(c(12, 15,
                                12, 16,
                                13, 18,
                                12, 18,
                                12, 19), ncol=2, byrow=TRUE),
               `miR-155`=matrix(c(13, 16,
                                  13, 17,
                                  13, 18,
                                  15, 21,
                                  15, 22), ncol=2, byrow=TRUE),
               `miR-124`=matrix(c(11, 14,
                                  11, 15,
                                  11, 16,
                                  11, 17,
                                  11, 18), ncol=2, byrow=TRUE),
               `lsy-6`=matrix(c(12, 15,
                                12, 16,
                                12, 17,
                                10, 16,
                                 9, 16), ncol=2, byrow=TRUE),
               `miR-7-23nt`=matrix(c(13, 16,
                                     11, 15,
                                     13, 18,
                                     11, 17,
                                     10, 17), ncol=2, byrow=TRUE))

  if (mirna == "miR-7-23nt")  experiment <- "equilibrium2_nb"
  else                        experiment <- "equilibrium"
  lim_use <- lims[[mirna]]
  # print(lim_use)

  offset_lim_use <- c(-2, 6)

  better <- c()
  sig_better <- c()
  for (row_i in seq(nrow(lim_use))) {
    start_mm <- lim_use[row_i, 1]
    stop_mm <- lim_use[row_i, 2]
    # print(sprintf("%s-%s", start_mm, stop_mm))
    kds_mmdb    <- GetThreePrimeMmBDKds(mirna, experiment, start_mm, stop_mm, n_constant=3, kd_fc=TRUE, bulge=bulge, win_average=3, corrected_kds=FALSE, best_average=TRUE, offset_lim=offset_lim_use)
    kds_mmdb_ub <- GetThreePrimeMmBDKds(mirna, experiment, start_mm, stop_mm, n_constant=3, kd_fc=TRUE, bulge=bulge, win_average=3, corrected_kds=FALSE, best_average=TRUE, offset_lim=offset_lim_use, lower_bound=TRUE)
    kds_mmdb_lb <- GetThreePrimeMmBDKds(mirna, experiment, start_mm, stop_mm, n_constant=3, kd_fc=TRUE, bulge=bulge, win_average=3, corrected_kds=FALSE, best_average=TRUE, offset_lim=offset_lim_use, upper_bound=TRUE)
    wt <- kds_mmdb["wt", drop=FALSE]
    wt_ub <- kds_mmdb_ub["wt", drop=FALSE]
    mm <- kds_mmdb[-1, drop=FALSE]
    mm_lb <- kds_mmdb_lb[-1, drop=FALSE]

    mm_inds_keep <- (!is.na(mm) & !is.na(mm_lb))
    mm <- mm[mm_inds_keep, drop=FALSE]
    mm_lb <- mm_lb[mm_inds_keep, drop=FALSE]
    # print("_________________________")
    # print(lim_use[row_i, ])
    # print(wt)
    # print(wt_ub)
    # print(mm)
    # print(mm_lb)
    if (sum(mm > wt) > 0) {
      # print("better")
      mismatch_sites <- sprintf("%smer-m%s.%smm%s", stop_mm - start_mm + 1, start_mm, stop_mm, names(mm)[which(mm > wt)])
      better <- c(better, mismatch_sites)
      # print(mismatch_sites)
    }
    if (sum(mm_lb > wt_ub) > 0) {
      # print("significantly better")
      mismatch_sites <- sprintf("%smer-m%s.%smm%s", stop_mm - start_mm + 1, start_mm, stop_mm, names(mm)[which(mm > wt)])
      sig_better <- c(sig_better, mismatch_sites)      
      # print(mismatch_sites)
    }
  }
  # Print which misamtches are better than teh wild type sites.
  print(mirna)
  print("better:")
  print(better)
  print("significantly better:")
  print(sig_better)
  return()
}

# CheckThreePSites("let-7a-21nt")

# CheckThreePSites("let-7a-21nt", bulge=TRUE)

# CheckThreePSites("miR-1")

# CheckThreePSites("miR-1", bulge=TRUE)

# CheckThreePSites("miR-155")

# CheckThreePSites("miR-155", bulge=TRUE)


CheckRandThreePSites("let-7a")
CheckRandThreePSites("let-7a", bulge=TRUE)

CheckRandThreePSites("miR-1")
CheckRandThreePSites("miR-1", bulge=TRUE)

CheckRandThreePSites("miR-155")
CheckRandThreePSites("miR-155", bulge=TRUE)

CheckRandThreePSites("miR-124")
CheckRandThreePSites("miR-124", bulge=TRUE)

CheckRandThreePSites("lsy-6")
CheckRandThreePSites("lsy-6", bulge=TRUE)

CheckRandThreePSites("miR-7-23nt")
CheckRandThreePSites("miR-7-23nt", bulge=TRUE)


break


test_data <- LoadBeckerEtAlData("let-7a")
test_data_alt <- LoadBeckerEtAlData("let-7a", alt=TRUE)

test_data_no_mm <- LoadBeckerEtAlData("let-7a", no_mm=TRUE)

print(head(test_data))
print(head(test_data_alt))

inds_diff <- which(test_data[, 1] != test_data_alt[, 1])
print(cbind(test_data[inds_diff, 1], test_data_alt[inds_diff, 1]))

inds_diff_2 <- which(test_data[, 1] != test_data_no_mm[, 1])
print(cbind(test_data[inds_diff_2, 1], test_data_no_mm[inds_diff_2, 1]))

inds_diff_3 <- which(test_data_alt[, 1] != test_data_no_mm[, 1])
print(cbind(test_data_alt[inds_diff_3, 1], test_data_no_mm[inds_diff_3, 1]))


break

sitelist <- "progthrp_suppcomp"
kds_prog <- EquilPars("let-7a-21nt", "equil_c2_nb", n_constant=3,
                      sitelist=sitelist)

kds_prog_c <- SwapProgrammedKds(kds=kds_let7, suppcomp=grepl("suppcomp", sitelist))
kds_rand <- EquilPars("let-7a", "equilibrium", n_constant=3, sitelist="randthrp_suppcomp")

print(head(kds_let7, n=30))
print(head(kds_let7_corrected, n=30))
print(head(kds_rand, n=30))


mm_ratio_prog <- (kds_prog["8mer-mmG6_Kd", ])/(kds_prog["8mer_Kd", ])
mm_ratio_prog_c <- (kds_prog_c["8mer-mmG6_Kd", ])/(kds_prog_c["8mer_Kd", ])
mm_ratio_rand <- (kds_rand["8mer-mmG6_Kd", ])/(kds_rand["8mer_Kd", ])

print(mm_ratio_prog)
print(mm_ratio_prog_c)
print(mm_ratio_rand)

out_mat <- MakePairingMatrix("let-7a-21nt", "equil_c2_nb", 1, n_constant=3,
                             corrected_kds=FALSE, kd_fc=TRUE)

out_mat_c <- MakePairingMatrix("let-7a-21nt", "equil_c2_nb", 1, n_constant=3,
                             corrected_kds=TRUE, kd_fc=TRUE)

out_mat_c_site <- MakePairingMatrix("let-7a-21nt", "equil_c2_nb", 1, n_constant=3, sitelist="progthrp",
                             corrected_kds=TRUE, kd_fc=TRUE, site_base="8mer-mmG6")


print(10^out_mat)
print(10^out_mat_c)
print(10^out_mat_c_site)



break

break



# test_1 <- MakeAllLibraryMismatches(10)

# test_2 <- MakeAllLibraryMismatches(10, bulge=TRUE)

test_1a_mm <- MakePositionalMismatchDf("let-7a-21nt", "equil_c2_nb", 10)

d_all_mm <- d_all

test_1a_bu <- MakePositionalMismatchDf("let-7a-21nt", "equil_c2_nb", 10, bulge=TRUE)

d_all_bu <- d_all

print(head(test_1a_mm))
print(head(d_all_mm))


print(head(test_1a_bu))
print(head(d_all_bu))


break


dev.new(height=4, width=4)
image(pairing_l7$MLE[-1, -1])

dev.new(xpos=400, ypos=0, height=4, width=4)
image(pairing_l7_shift$MLE[-nrow(pairing_l7_shift$MLE), -ncol(pairing_l7_shift$MLE)])

break

preds_orig <- GetTargetScanPreds("HeLa", train_mirnas="seventeen", norm=TRUE)

preds_new <- GetTargetScanPreds("HeLa", train_mirnas="seventeen", norm=TRUE, two_bmodes=TRUE)

preds_tpm <- GetTargetScanPreds("HeLa", train_mirnas="seventeen", tpms=TRUE, norm=TRUE)

features_old <- GetTargetScanFeatureFile("HeLa")
features_tbmodes <- GetTargetScanFeatureFile("HeLa", two_bmodes=TRUE)

feature_pairings <- features_old[, "Threep_pairing"]

pairings <- sapply(feature_pairings, function(item) {
  outs <- unlist(strsplit(item, split="(?:_|\\|)", perl=TRUE))
  return(as.integer(outs[2]) - as.integer(outs[1]) + 1)  
})
pairings[is.na(pairings)] <- 0


lens_uniq <- sort(unique(pairings))
lens_uniq <- seq(0, max(lens_uniq))

matrix_preds <- matrix(0, nrow=length(lens_uniq), ncol=ncol(preds_tpm), dimnames=list(lens_uniq, colnames(preds_tpm)))

matrix_tpms <- matrix(0, nrow=length(lens_uniq), ncol=ncol(preds_tpm), dimnames=list(lens_uniq, colnames(preds_tpm)))

for (len_i in lens_uniq) {
  for (mirna in colnames(preds_tpm)) {
    inds_mir <- which(features_old[, "miRNA.family"] == mirna)
    inds_len <- which(pairings == len_i)
    inds_use <- intersect(inds_len, inds_mir)
    trans <- features_old[inds_use, "Gene.ID"]
    matrix_preds[len_i, mirna] <- mean(preds_orig[trans, mirna], na.rm=TRUE)
    matrix_tpms[len_i, mirna] <- mean(preds_tpm[trans, mirna], na.rm=TRUE)
  }
}

# dev.new(xpos=20, ypos=20, height=4, width=4)

# plot(lens_uniq, matrix_preds[, 1], type="l", ylim=c(-1.5, 0.5))
# cols <- c("blue", "green", "purple", "red", "orange", "darkmagenta", "darkslateblue", "cyan")
# for (col_i in seq(2, ncol(matrix_preds))) {
#   lines(lens_uniq, matrix_preds[, col_i], type="l", col=cols[col_i - 1])
# }

# dev.new(xpos=520, ypos=20, height=4, width=4)

# plot(lens_uniq, matrix_tpms[, 1], type="l", ylim=c(-1.5, 0.5))
# cols <- c("blue", "green", "purple", "red", "orange", "darkmagenta", "darkslateblue", "cyan")
# for (col_i in seq(2, ncol(matrix_preds))) {
#   lines(lens_uniq, matrix_tpms[, col_i], type="l", col=cols[col_i - 1])
# }


inds_m124 <- which(features_old[, "miRNA.family"] == "mir124")


inds_pairing <- grep("^11\\|.*_(?:2|3|4|5)", features_tbmodes[, "Threep_pairing"], perl=TRUE)

inds_look <- intersect(inds_m124, inds_pairing)

print(inds_look)

transcripts_old = features_old[inds_look, "Gene.ID"]
transcripts_new = features_tbmodes[inds_look, "Gene.ID"]

print(transcripts_old)
print(transcripts_new)

dev.new(xpos=20, ypos=20, height=5, width=5)
plot(preds_orig[, "mir124"], preds_tpm[, "mir124"], col=rgb(0, 0, 0, alpha=0.1), xlim=c(-2.5, 1), ylim=c(-2.5, 1))
segments(x0=-2.5, y0=-2.5, x1=1, y1=1, lty=2, col="gray")
points(preds_orig[transcripts_old, "mir124"], preds_tpm[transcripts_old, "mir124"], col="blue")

dev.new(xpos=520, ypos=20, height=5, width=5)
plot(preds_new[, "mir124"], preds_tpm[, "mir124"], col=rgb(0, 0, 0, alpha=0.1), xlim=c(-2.5, 1), ylim=c(-2.5, 1))
segments(x0=-2.5, y0=-2.5, x1=1, y1=1, lty=2, col="gray")
points(preds_new[transcripts_old, "mir124"], preds_tpm[transcripts_old, "mir124"], col="blue")

print("Cor all")
print(cor(unlist(preds_orig), unlist(preds_tpm))^2)
print(cor(unlist(preds_new), unlist(preds_tpm))^2)
print("Cor subset")
print(cor(preds_orig[, "mir124"], preds_tpm[, "mir124"])^2)
print(cor(preds_new[, "mir124"], preds_tpm[, "mir124"])^2)
print("loss subset")
print(sum((preds_orig[transcripts_old, "mir124"] -  preds_tpm[transcripts_old, "mir124"])^2, na.rm=TRUE))
print(sum((preds_new[transcripts_old, "mir124"] -  preds_tpm[transcripts_old, "mir124"])^2, na.rm=TRUE))


inds_let7 <- which(features_old[, "miRNA.family"] == "let7")

inds_look <- intersect(inds_let7, inds_pairing)

print(inds_look)

transcripts_old = features_old[inds_look, "Gene.ID"]
transcripts_new = features_tbmodes[inds_look, "Gene.ID"]

# print(transcripts_old)
# print(transcripts_new)

dev.new(xpos=20, ypos=520, height=5, width=5)
plot(preds_orig[, "let7"], preds_tpm[, "let7"], col=rgb(0, 0, 0, alpha=0.1), xlim=c(-2.5, 1), ylim=c(-2.5, 1))
segments(x0=-2.5, y0=-2.5, x1=1, y1=1, lty=2, col="gray")
points(preds_orig[transcripts_old, "let7"], preds_tpm[transcripts_old, "let7"], col="blue")

dev.new(xpos=520, ypos=520, height=5, width=5)
plot(preds_new[, "let7"], preds_tpm[, "let7"], col=rgb(0, 0, 0, alpha=0.1), xlim=c(-2.5, 1), ylim=c(-2.5, 1))
segments(x0=-2.5, y0=-2.5, x1=1, y1=1, lty=2, col="gray")
points(preds_new[transcripts_old, "let7"], preds_tpm[transcripts_old, "let7"], col="blue")


print("Cor all")
print(cor(unlist(preds_orig), unlist(preds_tpm))^2)
print(cor(unlist(preds_new), unlist(preds_tpm))^2)
print("Cor subset")
print(cor(preds_orig[, "let7"], preds_tpm[, "let7"])^2)
print(cor(preds_new[, "let7"], preds_tpm[, "let7"])^2)
print("loss subset")
print(sum((preds_orig[transcripts_old, "let7"] -  preds_tpm[transcripts_old, "let7"])^2, na.rm=TRUE))
print(sum((preds_new[transcripts_old, "let7"] -  preds_tpm[transcripts_old, "let7"])^2, na.rm=TRUE))



break


graphics.off()


tf_feat1 <- read.table("/lab/solexa_bartel/klin/miRNA_models_data_old/ts7_outputs/refseq_utr90/features.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

tf_feat2 <- read.table("/lab/bartel4_ata/kathyl/RNA_Seq/outputs/ts7/features.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)


CheckForUniqueTranscripts <- function(tf_feat) {
  colnames(tf_feat) <- gsub("\\.", replacement="_", colnames(tf_feat), perl=TRUE)
  print(tf_feat[1:5, 1:4])
  genes <- tf_feat$Gene_ID
  mirnas <- tf_feat$miRNA_family
  print(head(genes))
  print(head(mirnas))
  mirna_gene_pairs <- paste0(mirnas, "_", genes)
  print(head(mirna_gene_pairs))
  uniq_mirna_gene_pairs <- unique(mirna_gene_pairs)
  print(length(uniq_mirna_gene_pairs))
  print(length(mirna_gene_pairs))
}

CheckForUniqueTranscripts(tf_feat1)
CheckForUniqueTranscripts(tf_feat2)
# WritePairingAndOffsetModel("let-7a-21nt", "equil_c2_nb", decomp_pairing=TRUE)
# WritePairingAndOffsetModel("miR-1", "equil_c_nb", decomp_pairing=TRUE)
# WritePairingAndOffsetModel("miR-155", "equil_sc_nb", decomp_pairing=TRUE)

# WritePairingAndOffsetModel("let-7a", "equilibrium", sitelist="randthrp_suppcomp", len_lim=c(4, 8), corrected_kds=FALSE, decomp_pairing=TRUE)
# WritePairingAndOffsetModel("miR-1", "equilibrium", sitelist="randthrp_suppcomp", len_lim=c(4, 8), corrected_kds=FALSE, decomp_pairing=TRUE)
# WritePairingAndOffsetModel("miR-155", "equilibrium", sitelist="randthrp_suppcomp", len_lim=c(4, 8), corrected_kds=FALSE, decomp_pairing=TRUE)
# WritePairingAndOffsetModel("miR-124", "equilibrium", sitelist="randthrp_suppcomp", len_lim=c(4, 8), corrected_kds=FALSE, decomp_pairing=TRUE)
# WritePairingAndOffsetModel("lsy-6", "equilibrium", sitelist="randthrp_suppcomp", len_lim=c(4, 8), corrected_kds=FALSE, decomp_pairing=TRUE)
# WritePairingAndOffsetModel("miR-7-23nt", "equilibrium2_nb", sitelist="randthrp_suppcomp", len_lim=c(4, 8), corrected_kds=FALSE, decomp_pairing=TRUE)



# coefs_all_2 <- coefs_all_global

# inds_shared <- intersect(rownames(coefs_all_1), rownames(coefs_all_2))

# dev.new(xpos=20, ypos=20, height=5, width=5)

# plot(coefs_all_1[inds_shared, 1], coefs_all_2[inds_shared, 1])


break
# out <- c()

# for (i in 9:15) {
#   print("____________")
#   print(9)
#   print(i)
#   difs <- CalculateThreePMismatchRange("miR-155", "equil_sc_nb", 9, i)
#   print(difs)
#   difs_use <- difs[which(difs != 0)]
#   out <- c(out, difs_use)
# }
# for (i in 9:14) {
#   print("____________")
#   print(10)
#   print(i)
#   difs <- CalculateThreePMismatchRange("miR-155", "equil_sc_nb", 10, i)
#   print(difs)
#   difs_use <- difs[which(difs != 0)]
#   out <- c(out, difs_use)
# }

# print(median(out))
# print(GeoMean(out))

# quick test of range for miR-155

inds_use <- which(!is.na(colMeans(out_mat_check_1)))

out_use <- out_mat_check_1[, inds_use]

print(dim(out_use))

means <- colMeans(out_use)

sort_inds <- order(means, decreasing=TRUE)

out_use <- out_use[, sort_inds]

out_diff <- t(t(out_use) - colMeans(out_use))

print(out_diff)
print(colMeans(out_diff))
print(dim(out_diff))
dims_use <- round(54/4)

test_check <- rowMeans(out_diff[, 1:dims_use])

print(test_check)

print(10^(max(test_check) - min(test_check)))

break


dev.new(xpos=20, ypos=20, height=5, width=5)
plot(ecdf(out), log="x", xlim=c(0.1, 200))


alt_test <- apply(out_mat_check_1, 2, range, na.rm=TRUE)
norm_check <- 10^(alt_test[2, ] - alt_test[1, ])
dev.new(xpos=520, ypos=20, height=5, width=5)

plot(ecdf(norm_check), log="x", xlim=c(0.1, 200))

break


# print(CalculateThreePMismatchRange("miR-155", "equil_sc_nb", 9, 9))
# print(CalculateThreePMismatchRange("miR-155", "equil_sc_nb", 9, 10))
# print(CalculateThreePMismatchRange("miR-155", "equil_sc_nb", 9, 11))
# print(CalculateThreePMismatchRange("miR-155", "equil_sc_nb", 9, 12))
# print(CalculateThreePMismatchRange("miR-155", "equil_sc_nb", 9, 13))
# print(CalculateThreePMismatchRange("miR-155", "equil_sc_nb", 9, 14))
# print(CalculateThreePMismatchRange("miR-155", "equil_sc_nb", 9, 15))

# print(CalculateThreePMismatchRange("miR-155", "equil_sc_nb", 10, 9))
# print(CalculateThreePMismatchRange("miR-155", "equil_sc_nb", 10, 10))
# print(CalculateThreePMismatchRange("miR-155", "equil_sc_nb", 10, 11))
# print(CalculateThreePMismatchRange("miR-155", "equil_sc_nb", 10, 12))
# print(CalculateThreePMismatchRange("miR-155", "equil_sc_nb", 10, 13))
print(CalculateThreePMismatchRange("miR-155", "equil_sc_nb", 10, 14))




break


mirna <- "let-7a-21nt"
experiment <- "equil_c2_nb"


check1 <- FitPairingOffsetAndMismatchModelWithError(mirna, experiment)


mirna <- "miR-155"
experiment <- "equil_sc_nb"

check2 <- FitPairingOffsetAndMismatchModelSingle(mirna, experiment)

break

experiments_use <- matrix(c("let-7a-21nt", "equil_c2_nb",
                            "miR-1", "equil_c_nb",
                            "miR-155", "equil_sc_nb",
                            "let-7a_minus1", "equil_c_nb",
                            "let-7a_plus1", "equil_c_nb",
                            "let-7a_miR-155", "equil_c_nb",
                            "miR-155_let-7a", "equil_c_nb"), nrow=7, byrow=TRUE)



apply(experiments_use, 1, function(row_i) {
  mirna <- row_i[1]
  experiment <- row_i[2]
  print(mirna)
  print(experiment)
  time_1 <- proc.time()[3]
  check0 <<- FitPairingAndOffsetModelSingle(mirna, experiment)
  time_2 <- proc.time()[3]
  check1 <<- FitPairingAndOffsetModelSingle(mirna, experiment, additive=TRUE)
  time_3 <- proc.time()[3]
  check2 <<- FitPairingAndOffsetModelSingle(mirna, experiment, intercept=TRUE)
  time_4 <- proc.time()[3]
  print(time_2 - time_1)
  print(time_3 - time_2)
  print(time_4 - time_3)
  print(cor(check0$data$logkd, check0$values)^2)
  print(cor(check1$data$logkd, check1$values)^2)
  print(cor(check2$data$logkd, check2$values)^2)
})

break
# dev.new(xpos=520, ypos=520, height=5, width=5)


# offset_lim_use <- c(4, 5)
# check2 <- FitPairingAndOffsetModel(mirna, experiment, offset_lim=offset_lim_use, fixed_offset=4, par_replacement_value=0)
# rownames(check2$coefs)[nrow(check2$coefs)] <- as.character(offset_lim_use)[2]

# check3 <- FitPairingAndOffsetModel(mirna, experiment, weights_new=TRUE)

dev.new(xpos=20, ypos=20, height=5, width=5)
inds_plot <- which(check1$data$offset == 4)
plot(check0$values[inds_plot], check0$data$logkd[inds_plot], col="blue", xlim=c(0, 3), ylim=c(0, 3))
points(check1$values[inds_plot], check1$data$logkd[inds_plot])
segments(x0=-3, y0=-3, x1=10, y1=10, lty=2)

dev.new(xpos=520, ypos=20, height=5, width=5)

plot(check0$coefs[, 1], check1$coefs[, 1], col="blue", xlim=c(-3, 5), ylim=c(-3, 5))
# points(check1$values[inds_plot], check1$data$logkd[inds_plot])
segments(x0=-3, y0=-3, x1=10, y1=10, lty=2)

break



break


inds_plot2 <- which(check2$data$offset == 4)


dev.new(xpos=520, ypos=20, height=5, width=5)
plot(check0$values[inds_plot], check0$data$logkd[inds_plot], col="blue", xlim=c(0, 3), ylim=c(0, 3))
points(check2$values[inds_plot2], check2$data$logkd[inds_plot2])
# points(check0$values[inds_plot], check1$data$logkd[inds_plot], , col="red")
segments(x0=-3, y0=-3, x1=10, y1=10, lty=2)


# inds_plot <- which(check1$data$offset == offset_lim_use[2])

# inds_plot2 <- which(check2$data$offset == offset_lim_use[2])


dev.new(xpos=1020, ypos=20, height=5, width=5)
plot(check0$values[inds_plot], check0$data$logkd[inds_plot], col="blue", xlim=c(0, 3), ylim=c(0, 3))
points(check3$values[inds_plot], check3$data$logkd[inds_plot])
# points(check0$values[inds_plot], check1$data$logkd[inds_plot], , col="red")
segments(x0=-3, y0=-3, x1=10, y1=10, lty=2)




# dev.new(xpos=20, ypos=520, height=5, width=5)
# plot(check2$values, check2$data$logkd, xlim=c(0, 3), ylim=c(0, 3))
# # points(check0$values[inds_plot], check1$data$logkd[inds_plot], , col="red")
# segments(x0=-3, y0=-3, x1=10, y1=10, lty=2)



break






# print(cor(test_model_1$data$logkd, test_model_1$values)^2)

# coefs_global_1 <- coefs_global

# dmodel.dpars_1 <- dmodel.dpars

# test_model_2 <- FitPairingOffsetAndMismatchModelSingle(mirna, experiment, log_plus_one=TRUE, fixed_offset=1, fixed_mm=1)
# print(cor(test_model_2$data$logkd, test_model_2$values)^2)



data_temp <- test_model_2$data

f1 <- gnls(logkd ~ log10(exp(a + b + c) + 1), data=data_temp,
           params=list(a ~pairing, b ~ offset, c ~ mm), start=c(rep(2, 107)), verbose=TRUE)




# test_model_2 <- FitPairingOffsetAndMismatchModelSingle(mirna, experiment, fixed_offset=1, fixed_mm=1, log_plus_one=TRUE)
# print(cor(test_model_1$data$logkd, test_model_1$values)^2)

break

dev.new(xpos=520, ypos=20, height=4, width=4)
cols <- rep("red", nrow(test_model_1$coefs))
inds_blue <- grep("\\|", rownames(test_model_1$coefs), perl=TRUE)
cols[inds_blue] <- "blue"
plot(test_model_1$coefs[, 1], test_model_2$coefs[, 1], col=cols)
break

# test_model_2 <- FitPairingAndOffsetModelNew(mirna, experiment, F_method=FALSE)


# print(cor(test_model_2$data$logkd, test_model_2$values)^2)


for (i in seq(73)) {
  print(sum(solve(A_mat[-i, -i])))
}

# test_model_mismatch <- FitPairingOffsetAndMismatchModelSingle(
#   "let-7a-21nt", "equil_c2_nb", fixed_offset=0, fixed_mm=1
# )


break






# test_model_2 <- FitPairingAndOffsetModel("let-7a_miR-155", "equil_c_nb", fixed_offset=-4, exponential=FALSE, lower_alt=TRUE, F_method=TRUE)




plot(test_model_1$coefs[, 1], test_model_2$coefs[, 1])

break

# data_temp <- test_model_1$data

# test_1_plot <- data_temp[which(data_temp$pairing == "11|14"), ]

# plot(test_1_plot$offset, test_1_plot$logkd, ylim=c(-1, 1))


# test_1_plot <- data_temp[which(data_temp$pairing == "10|14"), ]

# points(test_1_plot$offset, test_1_plot$logkd, col="blue")

# test_1_plot <- data_temp[which(data_temp$pairing == "10|13"), ]

# points(test_1_plot$offset, test_1_plot$logkd, col="green")

# test_1_plot <- data_temp[which(data_temp$pairing == "11|17"), ]

# points(test_1_plot$offset, test_1_plot$logkd, col="purple")


# break

# print(data_temp[which(data_temp$pairing == "10|13"), ])

# print(data_temp[which(data_temp$pairing == "12|15"), ])





break
test_model <- FitPairingAndOffsetModelForAllOffsets("let-7a_miR-155", "equil_c_nb")

break


check <- CalculateThreePScore(c(9, 18), 0)

print(check)


test_l7 <- GetThreePSiteDeltaG("let-7a", 10, 20)

test_m155 <- GetThreePSiteDeltaG("miR-155", 12, 22)

print(test_l7)
print(test_m155)

break

GetMismatch8merReadFractions <- function(mirna, experiment) {
  # Get the mismatch 8mer sites
  sites_use <- sapply(GetAll8merMmSites(mirna), GetSiteSeq, mirna=mirna)
  grep_strs <- sprintf("^.{25}%s", sites_use)
  # Get the path.
  path_base <- "/lab/solexa_bartel/mcgeary/AgoRBNS/%s/%s/reads/I.txt"
  path_use <- sprintf(path_base, mirna, experiment)
  wc_command <- sprintf("wc -l %s", path_use)
  full_lines <- as.numeric(unlist(strsplit(system(wc_command, intern=TRUE), split=" "))[1])


  out <- sapply(grep_strs, function(grep_str) {

    command <- sprintf("grep -P '%s' %s | wc -l", grep_str, path_use)

    print(command)
    check <- as.numeric(system(command, intern=TRUE))
    print(check)
    return(check)
  })
  print(out)
  return(out/full_lines)

}


percentages_l7_1 <- GetMismatch8merReadFractions("let-7a-21nt", "equil_c2_nb")
percentages_l7_2 <- GetMismatch8merReadFractions("let-7a-21nt", "equil_c_nb")

percentages_m1 <- GetMismatch8merReadFractions("miR-1", "equil_c_nb")

percentages_m155 <- GetMismatch8merReadFractions("miR-155", "equil_sc_nb")

percentages_l7p1 <- GetMismatch8merReadFractions("let-7a_plus1", "equil_c_nb")
percentages_l7m1 <- GetMismatch8merReadFractions("let-7a_minus1", "equil_c_nb")

percentages_l7m155 <- GetMismatch8merReadFractions("let-7a_miR-155", "equil_c_nb")
percentages_m155l7 <- GetMismatch8merReadFractions("miR-155_let-7a", "equil_c_nb")



print(percentages)

break

break



break


PlotKdsAgainstThreePScore <- function(
  mirna, experiment, plotscore=TRUE, xpos=20, ypos=20, height=4, width=4
) {
  model <- SubfunctionCall(FitPairingAndOffsetModelForAllOffsets)
  scores_3p <- SubfunctionCall(GetThreePrimeScoreMatrix)

  val_3p_score <- apply(model$data, 1, function(row) {
    CalculateThreePScore(c(as.integer(row[3]), as.integer(row[4])), as.integer(row[7]))
  })

  dev.new(xpos=xpos, ypos=ypos, height=height, width=height)
  if (plotscore) {
    plot(val_3p_score, model$data$logkd)
  } else {
    plot(model$values, model$data$logkd)
  }
}

PlotKdsAgainstThreePScore("let-7a-21nt", "equil_c2_nb")

break

print(scores_3p)





mirna <- "miR-1"
experiment <- "equil_c_nb"

model <- FitPairingAndOffsetModelForAllOffsets(mirna, experiment)

print(model$pairing$MLE)

scores_3p <- GetThreePrimeScoreMatrix(mirna)

print(scores_3p)

val_3p_score <- apply(model$data, 1, function(row) {
  CalculateThreePScore(c(as.integer(row[3]), as.integer(row[4])), as.integer(row[7]))
})

dev.new(xpos=20, ypos=420, height=4, width=4)

plot(val_3p_score, model$data$logkd)

dev.new(xpos=420, ypos=420, height=4, width=4)

plot(model$values, model$data$logkd)




mirna <- "miR-155"
experiment <- "equil_sc_nb"

model <- FitPairingAndOffsetModelForAllOffsets(mirna, experiment)

print(model$pairing$MLE)

scores_3p <- GetThreePrimeScoreMatrix(mirna)

print(scores_3p)

val_3p_score <- apply(model$data, 1, function(row) {
  CalculateThreePScore(c(as.integer(row[3]), as.integer(row[4])), as.integer(row[7]))
})

dev.new(xpos=20, ypos=820, height=4, width=4)

plot(val_3p_score, model$data$logkd)

dev.new(xpos=420, ypos=820, height=4, width=4)

plot(model$values, model$data$logkd)





break




break


becker_df  <- MakeBeckerThreePDataFrame("let-7a")
counts <- sapply(GetAll8merMmSites("let-7a"), function(site) {
  nrow(subset(becker_df, mm == site & offset >= 0 & offset <= 4))
})

print(counts)
print(names(counts)[which.min(counts)])
print(names(counts)[which.max(counts)])

break


becker_matrix <- MakeBeckerPairingMatrix("let-7a", 0, "6mer")
dev.new(xpos=20, ypos=20, height=3, width=3)
image(becker_matrix)

becker_matrix <- MakeBeckerPairingMatrix("let-7a", 1, "6mer")
dev.new(xpos=20, ypos=20, height=3, width=3)
image(becker_matrix)

becker_matrix <- MakeBeckerPairingMatrix("let-7a", 2, "6mer")
dev.new(xpos=220, ypos=20, height=3, width=3)
image(becker_matrix)

becker_matrix <- MakeBeckerPairingMatrix("let-7a", 3, "6mer")
dev.new(xpos=420, ypos=20, height=3, width=3)
image(becker_matrix)

becker_matrix <- MakeBeckerPairingMatrix("let-7a", 4, "6mer")
dev.new(xpos=620, ypos=20, height=3, width=3)
image(becker_matrix)

becker_matrix <- MakeBeckerPairingMatrix("let-7a", 5, "6mer")
dev.new(xpos=820, ypos=20, height=3, width=3)
image(becker_matrix)

becker_matrix <- MakeBeckerPairingMatrix("let-7a", 6, "6mer")
dev.new(xpos=1020, ypos=20, height=3, width=3)
image(becker_matrix)


break
counts_l7 <- CountAllMismatchSites("let-7a-21nt", "equil_c2_nb")
counts_m1 <- CountAllMismatchSites("miR-1", "equil_c_nb")
counts_m155 <- CountAllMismatchSites("miR-155", "equil_sc_nb")


print(counts_l7)
print(counts_m1)
print(counts_m155)
break



test <- MakeAllMirnaDeltaDeltaGMismatchDf("programmed", 9)
print(test)


break




sXc <- SitesXCounts("miR-155", "equil_sc_nb", 3, "progthrp_suppcomp")
sXc_mm <- SitesXCounts("miR-155", "equil_sc_nb", 3, "progthrp_suppcomp", start_mm=13, stop_mm=23, new=TRUE)
sXc_mm2 <- SitesXCounts("miR-155", "equil_sc_nb", 3, "progthrp_suppcomp", start_mm=13, stop_mm=22, new=TRUE)


print(sXc["11mer-m13.23|13|Comp", ])
print(sXc_mm["11mer-m13.23|13|Comp", ])
print(sXc_mm["11mer-m13.23mmC23|13|Comp", ])
print(sXc_mm["11mer-m13.23mmG23|13|Comp", ])
print(sXc_mm["11mer-m13.23mmT23|13|Comp", ])

print(sXc_mm2["10mer-m13.22|13|Comp", ])

print(sXc["10mer-m13.22|13|Comp", ])

print(sXc_mm["11mer-m13.23bC(19.23)|13|Comp", ])

extras <- matrix(c(0, 0, 31, 35, 68, 23, 22, 2,
                   0, 0, 21, 36, 30, 39, 4, 2,
                   2, 2, 13, 22, 12, 10, 3, 3), nrow=3, byrow=TRUE)
colnames(extras) <- colnames(sXc_mm)

matrix_out <- rbind(sXc_mm[c("11mer-m13.23mmC23|13|Comp", "11mer-m13.23mmG23|13|Comp", "11mer-m13.23mmT23|13|Comp"), ],
                    extras)

print(matrix_out)
print(colSums(matrix_out))




p_wt <- t(t(sXc)/colSums(sXc))
p_mm <- t(t(sXc_mm)/colSums(sXc_mm))
p_mm2 <- t(t(sXc_mm2)/colSums(sXc_mm2))


R_wt <- (p_wt/(p_wt[, 1]))[, 3:7]
R_mm <- (p_mm/(p_mm[, 1]))[, 3:7]
R_mm2 <- (p_mm2/(p_mm2[, 1]))[, 3:7]

dev.new(xpos=20, ypos=20, height=5, width=5)

plot(as.numeric(colnames(R_mm)), R_mm["11mer-m13.23|13|Comp", ], log="xy", xlim=c(0.1, 100), ylim=c(0.1, 1000))

lines(as.numeric(colnames(R_wt)), R_wt["11mer-m13.23|13|Comp", ], col="black", lwd=2)

lines(as.numeric(colnames(R_mm)), R_mm["11mer-m13.23mmC23|13|Comp", ], col="blue")
lines(as.numeric(colnames(R_mm)), R_mm["11mer-m13.23mmG23|13|Comp", ], col="red")
lines(as.numeric(colnames(R_mm)), R_mm["11mer-m13.23mmT23|13|Comp", ], col="forestgreen")

lines(as.numeric(colnames(R_mm2)), R_mm2["10mer-m13.22|13|Comp", ], col="gray", lwd=3, lty=3)











sXc <- SitesXCounts("miR-155", "equil_sc_nb", 3, "progthrp_suppcomp")
sXc_mm <- SitesXCounts("miR-155", "equil_sc_nb", 3, "progthrp_suppcomp", start_mm=13, stop_mm=22, new=TRUE)
sXc_mm2 <- SitesXCounts("miR-155", "equil_sc_nb", 3, "progthrp_suppcomp", start_mm=13, stop_mm=21, new=TRUE)


print(sXc["10mer-m13.22|12|Comp", ])
print(sXc_mm["10mer-m13.22|13|Comp", ])
print(sXc_mm["10mer-m13.22mmA23|13|Comp", ])
print(sXc_mm["10mer-m13.22mmG23|13|Comp", ])
print(sXc_mm["10mer-m13.22mmT23|13|Comp", ])

print(sXc_mm2["9mer-m13.21|13|Comp", ])

print(sXc["9mer-m13.21|13|Comp", ])

# print(sXc_mm["10mer-m13.22bC(19.23)|13|Comp", ])

# extras <- matrix(c(0, 0, 31, 35, 68, 23, 22, 2,
#                    0, 0, 21, 36, 30, 39, 4, 2,
#                    2, 2, 13, 22, 12, 10, 3, 3), nrow=3, byrow=TRUE)
# colnames(extras) <- colnames(sXc_mm)

# matrix_out <- rbind(sXc_mm[c("11mer-m13.23mmC23|13|Comp", "11mer-m13.23mmG23|13|Comp", "11mer-m13.23mmT23|13|Comp"), ],
#                     extras)

# print(matrix_out)
# print(colSums(matrix_out))




p_wt <- t(t(sXc + 1)/colSums(sXc + 1))
p_mm <- t(t(sXc_mm + 1)/colSums(sXc_mm + 1))
p_mm2 <- t(t(sXc_mm2 + 1)/colSums(sXc_mm2 + 1))


R_wt <- (p_wt/(p_wt[, 1]))[, 3:7]
R_mm <- (p_mm/(p_mm[, 1]))[, 3:7]
R_mm2 <- (p_mm2/(p_mm2[, 1]))[, 3:7]

dev.new(xpos=520, ypos=20, height=5, width=5)
plot(as.numeric(colnames(R_mm)), R_mm["10mer-m13.22|13|Comp", ], log="xy", xlim=c(0.1, 100), ylim=c(0.1, 1000))

lines(as.numeric(colnames(R_wt)), R_wt["10mer-m13.22|13|Comp", ], col="black", lwd=2)

lines(as.numeric(colnames(R_mm)), R_mm["10mer-m13.22mmA22|13|Comp", ], col="blue")
lines(as.numeric(colnames(R_mm)), R_mm["10mer-m13.22mmT22|13|Comp", ], col="forestgreen")
lines(as.numeric(colnames(R_mm)), R_mm["10mer-m13.22mmG22|13|Comp", ], col="red")

lines(as.numeric(colnames(R_mm)), R_mm["10mer-m13.22bA22|13|Comp", ], col="blue", lty=2)
lines(as.numeric(colnames(R_mm)), R_mm["10mer-m13.22bT22|13|Comp", ], col="forestgreen", lty=2)
lines(as.numeric(colnames(R_mm)), R_mm["10mer-m13.22bG22|13|Comp", ], col="red", lty=2)

lines(as.numeric(colnames(R_mm)),
      colMeans(R_mm[c("10mer-m13.22mmA22|13|Comp", "10mer-m13.22mmT22|13|Comp",
                      "10mer-m13.22mmG22|13|Comp", "10mer-m13.22bA22|13|Comp",
                      "10mer-m13.22bT22|13|Comp", "10mer-m13.22bG22|13|Comp"), ]), lwd=4, col="gray80")



lines(as.numeric(colnames(R_mm2)), R_mm2["9mer-m13.21|13|Comp", ], lwd=3, lty=3)







check <- EquilPars("miR-155", "equil_sc_nb", 3, "progthrp_suppcomp")
check_mm <- EquilPars("miR-155", "equil_sc_nb", 3, "progthrp_suppcomp", start_mm=13, stop_mm=22, new=TRUE)
check_mm2 <- EquilPars("miR-155", "equil_sc_nb", 3, "progthrp_suppcomp", start_mm=13, stop_mm=21, new=TRUE)


inds_13.22_wt <- grep("^10mer-m13\\.22\\|.*\\|Comp_Kd$", rownames(check), value=TRUE, perl=TRUE)
inds_13.22_mm <- grep("^10mer-m13\\.22\\|.*\\|Comp_Kd$", rownames(check_mm), value=TRUE, perl=TRUE)
inds_13.22mA22_mm <- grep("^10mer-m13\\.22mmA22\\|.*\\|Comp_Kd$", rownames(check_mm), value=TRUE, perl=TRUE)
inds_13.22mG22_mm <- grep("^10mer-m13\\.22mmG22\\|.*\\|Comp_Kd$", rownames(check_mm), value=TRUE, perl=TRUE)
inds_13.22mU22_mm <- grep("^10mer-m13\\.22mmT22\\|.*\\|Comp_Kd$", rownames(check_mm), value=TRUE, perl=TRUE)

inds_13.22bA22_mm <- grep("^10mer-m13\\.22bA22\\|.*\\|Comp_Kd$", rownames(check_mm), value=TRUE, perl=TRUE)
inds_13.22bG22_mm <- grep("^10mer-m13\\.22bG22\\|.*\\|Comp_Kd$", rownames(check_mm), value=TRUE, perl=TRUE)
inds_13.22bU22_mm <- grep("^10mer-m13\\.22bT22\\|.*\\|Comp_Kd$", rownames(check_mm), value=TRUE, perl=TRUE)


inds_13.22_mm2 <- grep("^10mer-m13\\.22\\|.*\\|Comp_Kd$", rownames(check_mm2), value=TRUE, perl=TRUE)
inds_13.21_mm2 <- grep("^9mer-m13\\.21\\|.*\\|Comp_Kd$", rownames(check_mm2), value=TRUE, perl=TRUE)

x_coords <- function(names) {
  as.numeric(gsub("^(.*)\\|(.*)\\|(.*)$", names, replacement="\\2", perl=TRUE))
}

print(check[inds_13.22_wt, ])
print(check_mm[inds_13.22_mm, ])
print(check_mm2[inds_13.22_mm2, ])

dev.new(xpos=820, ypos=20, height=5, width=5)
plot(x_coords(inds_13.22_wt), check[inds_13.22_wt, 2], log="y", ylim=c(1e-4, 1), type="l")



lines(x_coords(inds_13.22_mm), check_mm[inds_13.22_mm, 2], col="gray", lwd=2, lty=2)
lines(x_coords(inds_13.22_mm2), check_mm2[inds_13.22_mm2, 2], col="gray10", lwd=3, lty=3)


lines(x_coords(inds_13.22mA22_mm), check_mm[inds_13.22mA22_mm, 2], col="blue")
lines(x_coords(inds_13.22mU22_mm), check_mm[inds_13.22mU22_mm, 2], col="forestgreen")
lines(x_coords(inds_13.22mG22_mm), check_mm[inds_13.22mG22_mm, 2], col="red")

lines(x_coords(inds_13.22bA22_mm), check_mm[inds_13.22bA22_mm, 2], col="blue", lty=2)
lines(x_coords(inds_13.22bU22_mm), check_mm[inds_13.22bU22_mm, 2], col="forestgreen", lty=2)
lines(x_coords(inds_13.22bG22_mm), check_mm[inds_13.22bG22_mm, 2], col="red", lty=2)



out_mat <- matrix(NA, nrow=17, ncol=6, dimnames=list(9:25, c("mmA", "mmU", "mmG", "bA", "bU", "bG")))

out_mat[as.character(x_coords(inds_13.22mA22_mm)), "mmA"] <- check_mm[inds_13.22mA22_mm, 2]
out_mat[as.character(x_coords(inds_13.22mG22_mm)), "mmG"] <- check_mm[inds_13.22mG22_mm, 2]
out_mat[as.character(x_coords(inds_13.22mU22_mm)), "mmU"] <- check_mm[inds_13.22mU22_mm, 2]

out_mat[as.character(x_coords(inds_13.22bA22_mm)), "bA"] <- check_mm[inds_13.22bA22_mm, 2]
out_mat[as.character(x_coords(inds_13.22bG22_mm)), "bG"] <- check_mm[inds_13.22bG22_mm, 2]
out_mat[as.character(x_coords(inds_13.22bU22_mm)), "bU"] <- check_mm[inds_13.22bU22_mm, 2]


print(out_mat)

lines(rownames(out_mat), 1/apply(1/out_mat[, 1:3], 1, mean, na.rm=TRUE), lwd=4, col="violet")

lines(rownames(out_mat), 1/apply(1/out_mat, 1, mean, na.rm=TRUE), lwd=4, col="purple")




lines(x_coords(inds_13.21_mm2), check_mm2[inds_13.21_mm2, 2], lwd=1)
















mirna <- "miR-155"
experiment <- "equil_sc_nb"
start_mm <- 13
stop_mm <- 22

check_mm <- GetThreePrimeMmBDKds(mirna, experiment, start_mm, stop_mm, bulge=FALSE, n_constant=3, kd_fc=TRUE,
  win_average=1, best_average=FALSE, new=TRUE, corrected_kds=TRUE)

check_bu <- GetThreePrimeMmBDKds(mirna, experiment, start_mm, stop_mm, bulge=TRUE, n_constant=3, kd_fc=TRUE,
  win_average=1, best_average=FALSE, new=TRUE, corrected_kds=TRUE)

stop_mm <- 21


check_wt <- GetThreePrimeMmBDKds(mirna, experiment, start_mm, stop_mm, bulge=FALSE, n_constant=3, kd_fc=TRUE,
  win_average=1, best_average=FALSE, new=TRUE, corrected_kds=TRUE)


dev.new(xpos=20, ypos=20, height=5, width=5)
plot(rownames(check_mm), check_mm[, "wt"], log="y", xlim=c(-4, 16), ylim=c(1e-1, 1e3), type="o")

apply(check_mm[, 2:4], 2, function(row) {
  lines(rownames(check_mm), row, col="blue")
})


apply(check_bu[, 2:4], 2, function(row) {
  lines(rownames(check_mm), row, col="red")
})

check <- cbind(check_mm[, 2:4], check_bu[, 2:4])

lines(rownames(check), rowMeans(check, na.rm=TRUE), col="forestgreen")


lines(rownames(check_wt), check_wt[, "wt"], col="purple")


best_mm <- GetThreePrimeMmBDKds(mirna, experiment, start_mm, 22, bulge=FALSE, n_constant=3, kd_fc=TRUE,
  win_average=3, best_average=TRUE, new=TRUE, corrected_kds=TRUE)

best_bu <- GetThreePrimeMmBDKds(mirna, experiment, start_mm, 22, bulge=TRUE, n_constant=3, kd_fc=TRUE,
  win_average=3, best_average=TRUE, new=TRUE, corrected_kds=TRUE)

best_wt <- GetThreePrimeMmBDKds(mirna, experiment, start_mm, 21, bulge=FALSE, n_constant=3, kd_fc=TRUE,
  win_average=3, best_average=TRUE, new=TRUE, corrected_kds=TRUE)


model_mm <- GetThreePModelMismatchCoefficients(mirna, experiment, start_mm, 22, bulge=FALSE)$pairing

model_bu <- GetThreePModelMismatchCoefficients(mirna, experiment, start_mm, 22, bulge=TRUE)$pairing

model_mmbu_mm <- GetThreePModelMismatchCoefficients(mirna, experiment, start_mm, 22, mm_and_bulge=TRUE)$pairing

model_mmbu_b <- GetThreePModelMismatchCoefficients(mirna, experiment, start_mm, 22, mm_and_bulge=TRUE, bulge=TRUE)$pairing


model_wt <- GetThreePModelMismatchCoefficients(mirna, experiment, start_mm, 21)$pairing
model_mmbu_wt <- GetThreePModelMismatchCoefficients(mirna, experiment, start_mm, 21, mm_and_bulge=TRUE)$pairing



# print(best_mm)
# print(best_bu)
# print(best_wt)

# print(model_mm[1:4])
# print(model_bu[1:4])

print(model_mmbu_mm[1:4])
print(model_mmbu_b[1:4])

# print(model_wt[1:4])
print(model_mmbu_wt[1:4])

break








# check_all <- intersect(rownames(check), rownames(check_mm))
# check_all <- grep("&", check_all, value=TRUE, invert=TRUE)
# check_all <- grep("Comp", check_all, value=TRUE)

dev.new(xpos=20, ypos=20, height=8, width=8)
plot(check[check_all, 2], check_mm[check_all, 2], log="xy")
identify(check[check_all, 2], check_mm[check_all, 2], labels=check_all)


break

m1_11.21m <- MakeSiteMismatchMatrix("miR-1", "equil_c_nb", 11, 21, model_values=TRUE)
m1_11.20m <- MakeSiteMismatchMatrix("miR-1", "equil_c_nb", 11, 20, model_values=TRUE)
m1_11.19m <- MakeSiteMismatchMatrix("miR-1", "equil_c_nb", 11, 19, model_values=TRUE)
m1_11.18m <- MakeSiteMismatchMatrix("miR-1", "equil_c_nb", 11, 18, model_values=TRUE)
m1_12.18m <- MakeSiteMismatchMatrix("miR-1", "equil_c_nb", 12, 18, model_values=TRUE)
m1_12.17m <- MakeSiteMismatchMatrix("miR-1", "equil_c_nb", 12, 17, model_values=TRUE)
m1_12.16m <- MakeSiteMismatchMatrix("miR-1", "equil_c_nb", 12, 16, model_values=TRUE)
m1_12.15m <- MakeSiteMismatchMatrix("miR-1", "equil_c_nb", 12, 15, model_values=TRUE)

m1_11.21 <- MakeSiteMismatchMatrix("miR-1", "equil_c_nb", 11, 21, model_values=FALSE)
m1_11.20 <- MakeSiteMismatchMatrix("miR-1", "equil_c_nb", 11, 20, model_values=FALSE)
m1_11.19 <- MakeSiteMismatchMatrix("miR-1", "equil_c_nb", 11, 19, model_values=FALSE)
m1_11.18 <- MakeSiteMismatchMatrix("miR-1", "equil_c_nb", 11, 18, model_values=FALSE)
m1_12.18 <- MakeSiteMismatchMatrix("miR-1", "equil_c_nb", 12, 18, model_values=FALSE)
m1_12.17 <- MakeSiteMismatchMatrix("miR-1", "equil_c_nb", 12, 17, model_values=FALSE)
m1_12.16 <- MakeSiteMismatchMatrix("miR-1", "equil_c_nb", 12, 16, model_values=FALSE)
m1_12.15 <- MakeSiteMismatchMatrix("miR-1", "equil_c_nb", 12, 15, model_values=FALSE)

break


m155_13.23m <- MakeSiteMismatchMatrix("miR-155", "equil_sc_nb", 13, 23, model_values=TRUE)
m155_13.22m <- MakeSiteMismatchMatrix("miR-155", "equil_sc_nb", 13, 22, model_values=TRUE)
m155_13.21m <- MakeSiteMismatchMatrix("miR-155", "equil_sc_nb", 13, 21, model_values=TRUE)
m155_15.22m <- MakeSiteMismatchMatrix("miR-155", "equil_sc_nb", 15, 22, model_values=TRUE)
m155_15.21m <- MakeSiteMismatchMatrix("miR-155", "equil_sc_nb", 15, 21, model_values=TRUE)



m155_13.23 <- MakeSiteMismatchMatrix("miR-155", "equil_sc_nb", 13, 23, model_values=FALSE)
m155_13.22 <- MakeSiteMismatchMatrix("miR-155", "equil_sc_nb", 13, 22, model_values=FALSE)
m155_13.21 <- MakeSiteMismatchMatrix("miR-155", "equil_sc_nb", 13, 21, model_values=FALSE)
m155_15.22 <- MakeSiteMismatchMatrix("miR-155", "equil_sc_nb", 15, 22, model_values=FALSE)
m155_15.21 <- MakeSiteMismatchMatrix("miR-155", "equil_sc_nb", 15, 21, model_values=FALSE)


break

# opt_out_1 <- GetThreePModelMismatchCoefficients("let-7a-21nt", "equil_c2_nb", 10, 20)
# opt_out_2 <- GetThreePModelMismatchCoefficients("let-7a-21nt", "equil_c2_nb", 9, 18)

# pairing_alt_1 <- log10(GetThreePrimeMmBDKds("let-7a-21nt", "equil_c2_nb", 10, 20, best_average=TRUE, win_average=3, kd_fc=TRUE))
# pairing_alt_2 <- log10(GetThreePrimeMmBDKds("let-7a-21nt", "equil_c2_nb", 9, 18, best_average=TRUE, win_average=3, kd_fc=TRUE))

dG1 <- GetThreePSiteDeltaG("let-7a-21nt", 11, 21)
dG2 <- GetThreePSiteDeltaG("let-7a-21nt", 11, 21, mm="G11")
dG3 <- GetThreePSiteDeltaG("let-7a-21nt", 11, 21, mm="G12")
dG4 <- GetThreePSiteDeltaG("let-7a-21nt", 11, 21, mm="G13")

# dG5 <- GetThreePSiteDeltaG("let-7a-21nt", 11, 21, mm="bG13")


print(dG1)
print(dG2)
print(dG3)
print(dG4)
# print(dG5)
break


pairing_1 <- opt_out_1$pairing
pairing_2 <- opt_out_2$pairing

pairing_1 <- pairing_1[2:length(pairing_1)] - pairing_1["wt"]
pairing_2 <- pairing_2[2:length(pairing_2)] - pairing_2["wt"]

pairing_alt_1 <- pairing_alt_1[2:length(pairing_alt_1)] - pairing_alt_1["wt"]
pairing_alt_2 <- pairing_alt_2[2:length(pairing_alt_2)] - pairing_alt_2["wt"]


# dev.new(xpos=520, ypos=20, height=5, width=5)
# plot(names(opt_out_1$offsets), opt_out_1$offsets)
# points(names(opt_out_2$offsets), opt_out_2$offsets, type="o")

dev.new(xpos=1020, ypos=20, height=5, width=5)
ind_use <- intersect(names(pairing_1), names(pairing_2))
plot(pairing_1[ind_use], pairing_2[ind_use])
text(min(pairing_1[ind_use]), max(pairing_2[ind_use]),
     label=cor(pairing_1[ind_use], pairing_2[ind_use])^2, adj=c(0, 1))



dev.new(xpos=20, ypos=520, height=5, width=5)
ind_use <- intersect(names(pairing_alt_1), names(pairing_alt_2))
plot(pairing_alt_1[ind_use], pairing_alt_2[ind_use])
text(min(pairing_alt_1[ind_use], na.rm=TRUE), max(pairing_alt_2[ind_use], na.rm=TRUE),
     label=cor(pairing_alt_1[ind_use], pairing_alt_2[ind_use], use="pairwise.complete.obs")^2, adj=c(0, 1))


break

coefs_single <- FitPairingAndOffsetModel(
  "let-7a-21nt", "equil_c2_nb", 3, "progthrp_suppcomp", corrected_kds=TRUE,
  additive=TRUE
)

break

print(coefs_single)
coefs_1 <- FitPairingOffsetAndMismatchModelWithError(
  "let-7a-21nt", "equil_c2_nb", 3, "progthrp", corrected_kds=TRUE
)

plot(1:nrow(coefs_1$mm), coefs_1$mm[, 1])
lines(1:nrow(coefs_1$mm), coefs_1$mm[, 2], lty=2)
lines(1:nrow(coefs_1$mm), coefs_1$mm[, 3], lty=2)


break

# coefs_3 <- FitPairingOffsetAndMismatchModelSingle(
#   "let-7a-21nt", "equil_c2_nb", 3, "progthrp", corrected_kds=TRUE,
#   fixed_offset=-3, fixed_mm=1
# )

# coefs_4 <- FitPairingOffsetAndMismatchModelSingle(
#   "let-7a-21nt", "equil_c2_nb", 3, "progthrp", corrected_kds=TRUE,
#   fixed_offset=-3, fixed_mm=2
# )

dev.new(xpos=20, ypos=20, height=5, width=5)

plot(-4:16,  coefs_1$coefs[as.character(-4:16), 1], type="o", ylim=c(0, 1))
lines(-4:16, coefs_1$coefs[as.character(-4:16), 2], lty=2, type="l")
lines(-4:16, coefs_1$coefs[as.character(-4:16), 3], lty=2, type="l")

lines(-4:16, coefs_2$coefs[as.character(-4:16), 1], type="o", col="red")
lines(-4:16, coefs_2$coefs[as.character(-4:16), 2], lty=3, type="l", col="red")
lines(-4:16, coefs_2$coefs[as.character(-4:16), 3], lty=3, type="l", col="red")

lines(-4:16, coefs_3$coefs[as.character(-4:16), 1], type="o", col="purple")
lines(-4:16, coefs_3$coefs[as.character(-4:16), 2], lty=2, type="l", col="purple")
lines(-4:16, coefs_3$coefs[as.character(-4:16), 3], lty=2, type="l", col="purple")

lines(-4:16, coefs_4$coefs[as.character(-4:16), 1], type="o", col="green")
lines(-4:16, coefs_4$coefs[as.character(-4:16), 2], lty=3, type="l", col="green")
lines(-4:16, coefs_4$coefs[as.character(-4:16), 3], lty=3, type="l", col="green")

dev.new(xpos=520, ypos=20, height=5, width=5)
plot(1:18,  coefs_1$coefs[(nrow(coefs_1$coefs) - 17):nrow(coefs_1$coefs), 1], type="o", ylim=c(0.5, 1.5))
lines(1:18, coefs_1$coefs[(nrow(coefs_1$coefs) - 17):nrow(coefs_1$coefs), 2], type="l", lty=2)
lines(1:18, coefs_1$coefs[(nrow(coefs_1$coefs) - 17):nrow(coefs_1$coefs), 3], type="l", lty=2)


lines(1:18, coefs_2$coefs[(nrow(coefs_1$coefs) - 17):nrow(coefs_1$coefs), 1], type="o", col="red")
lines(1:18, coefs_2$coefs[(nrow(coefs_1$coefs) - 17):nrow(coefs_1$coefs), 2], type="l", lty=2, col="red")
lines(1:18, coefs_2$coefs[(nrow(coefs_1$coefs) - 17):nrow(coefs_1$coefs), 3], type="l", lty=2, col="red")

lines(1:18, coefs_3$coefs[(nrow(coefs_1$coefs) - 17):nrow(coefs_1$coefs), 1], type="o", col="purple")
lines(1:18, coefs_3$coefs[(nrow(coefs_1$coefs) - 17):nrow(coefs_1$coefs), 2], type="l", lty=2, col="purple")
lines(1:18, coefs_3$coefs[(nrow(coefs_1$coefs) - 17):nrow(coefs_1$coefs), 3], type="l", lty=2, col="purple")

lines(1:18, coefs_4$coefs[(nrow(coefs_1$coefs) - 17):nrow(coefs_1$coefs), 1], type="o", col="green")
lines(1:18, coefs_4$coefs[(nrow(coefs_1$coefs) - 17):nrow(coefs_1$coefs), 2], type="l", lty=2, col="green")
lines(1:18, coefs_4$coefs[(nrow(coefs_1$coefs) - 17):nrow(coefs_1$coefs), 3], type="l", lty=2, col="green")

break


test_0 <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a-21nt/equil_c2_nb/kds_PAPER/3_progthrp_suppcomp_PAPER_sorted.txt")
test_1 <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a-21nt/equil_c_nb/kds_PAPER/3_progthrp_suppcomp_PAPER_sorted.txt")
test_2 <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a-21nt/equil_c_nb/kds_PAPER/3_progthrp_suppcomp_nbomitc_PAPER_sorted.txt")


dev.new(xpos=20, ypos=20, height=5, width=5)
plot(ecdf(test_1["4mer-m11.14|13|Comp_Kd", ]))
plot(ecdf(test_0["4mer-m11.14|13|Comp_Kd", ]), add=TRUE, col="blue")
plot(ecdf(test_2["4mer-m11.14|13|Comp_Kd", ]), add=TRUE, col="red")


pars_0 <- EquilPars("let-7a-21nt", "equil_c2_nb", 3, "progthrp_suppcomp")
pars_1 <- EquilPars("let-7a-21nt", "equil_c_nb", 3, "progthrp_suppcomp")
pars_2 <- EquilPars("let-7a-21nt", "equil_c_nb", 3, "progthrp_suppcomp", nbomitc=TRUE)
break
plot(pars_0[, 1], pars_1[, 1], log="xy")
break

sXc_1 <- SitesXCounts("let-7a-21nt", "equil_c2_nb", 3, "progthrp_suppcomp")
sXc_2 <- SitesXCounts("let-7a-21nt", "equil_c_nb", 3, "progthrp_suppcomp")

f_1 <- t(t(sXc_1)/colSums(sXc_1))
f_2 <- t(t(sXc_2)/colSums(sXc_2))

R_1 <- f_1[, 3:7]/f_1[, 2]
R_2 <- f_2[, 3:7]/f_2[, 2]

sites_unique <- intersect(rownames(R_1), rownames(R_2))
sites_unique <- grep("&", sites_unique, invert=TRUE, value=TRUE)

col_use <- 6
# plot(R_1[sites_unique, col_use], R_2[sites_unique, col_use])
# identify(R_1[sites_unique, col_use], R_2[sites_unique, col_use], labels=sites_unique)

plot(c(40, 12.65, 4, 1.265, 0.4), R_1["7mer-m11.17|14|Comp", ], type="o", log="xy", ylim=c(0.1, 10))
lines(c(40, 12.65, 4, 1.265, 0.4), R_2["7mer-m11.17|14|Comp", ], col="red", type="o")
lines(c(40, 12.65, 4, 1.265, 0.4), R_1["8mer-mmG7", ], col="black", lty=2, type="o")
lines(c(40, 12.65, 4, 1.265, 0.4), R_2["8mer-mmG7", ], col="red", lty=2, type="o")


lines(c(40, 12.65, 4, 1.265, 0.4), R_1["7mer-m11.17|19|Comp", ], type="o", lty=3)
lines(c(40, 12.65, 4, 1.265, 0.4), R_2["7mer-m11.17|19|Comp", ], col="red", lty=3, type="o")



break

print(head(R_1))
print(head(R_2))
break



mir_prog <- "let-7a-21nt"
exp_prog <- "equil_c2_nb"
sitelist_prog_old <- "programmed_suppcomp"
sitelist_prog_new <- "progthrp"

offset_lim_use <- c(4, 4)
len_lim_use <- c(5, 11)
pos_lim_use <- c(10, 21)
sumseed_use <- FALSE

# out_1 <- MakeKdDataFrame(mir_prog, exp_prog, sitelist_prog_new,
#                          offset_lim=offset_lim_use, len_lim=len_lim_use,
#                          pos_lim=pos_lim_use, n_constant=3, sumseed=sumseed_use)

print(tail(out_1))

len_use <- 6
pos_start_use <- 11


out_2 <- MakeMismatchByOffsetMatrix(mir_prog, exp_prog, sitelist_prog_new,
                                    len=len_use, pos5p=pos_start_use, n_constant=3)


print(out_2)

mir_rand <- "let-7a"
exp_rand <- "equilibrium"
sitelist_rand <- "randthrp_comp"
# sitelist_prog_new <- "progthrp"

offset_lim_use <- c(-2, 12)
len_lim_use <- c(5, 11)
# pos_lim_use <- c(10, 21)
# sumseed_use <- FALSE
out_3 <- MakeMismatchByOffsetMatrix(mir_rand, exp_rand, sitelist_rand, corrected_kds=FALSE,
                                    len=len_use, supp_base=TRUE, pos5p=pos_start_use, n_constant=3)

print(out_3)

# identify(log10(1/out_2[, 4]), out_3[, 3], labels=rownames(out_2))
break


out_2 <- MakeKdDataFrame(mir_prog, exp_prog, sitelist_prog_new,
                         offset_lim=offset_lim_use, len_lim=len_lim_use,
                         pos_lim=pos_lim_use, n_constant=3, sumseed=sumseed_use, corrected_kds=TRUE)

out_3 <- MakeKdDataFrame(mir_prog, exp_prog, sitelist_prog_new,
                         offset_lim=offset_lim_use, len_lim=len_lim_use,
                         pos_lim=pos_lim_use, n_constant=3, sumseed=FALSE)

out_4 <- MakeKdDataFrame(mir_prog, exp_prog, sitelist_prog_new,
                         offset_lim=offset_lim_use, len_lim=len_lim_use,
                         pos_lim=pos_lim_use, n_constant=3, sumseed=FALSE, corrected_kds=TRUE)


print(unique(out_1$len))
print(unique(out_1$pos5p))
print(unique(out_1$pos3p))
print(unique(out_1$offset))

dev.new(xpos=20, ypos=20, height=5, width=5)
plot(out_1$logkd, out_2$logkd, xlim=c(0, 3), ylim=c(0, 3))
segments(x0=0, y0=0, x1=3, y1=3, lty=2)
points(out_1$logkd, out_3$logkd, col="red")
points(out_1$logkd, out_4$logkd, col="blue")


break

out  <- MakePairingMatrix("let-7a-21nt", "equil_c2_nb", 4, corrected_kds=FALSE)
out2 <- MakePairingMatrix("let-7a-21nt", "equil_c2_nb", 4, corrected_kds=FALSE, sumseed=TRUE)
out3 <- MakePairingMatrix("let-7a", "equilibrium", 4, sitelist="randthrp_suppcomp")[1:nrow(out), 1:ncol(out)]
out4 <- MakePairingMatrix("let-7a", "equilibrium", 4, sitelist="randthrp_suppcomp", sumseed=TRUE)[1:nrow(out), 1:ncol(out)]
out5 <- MakePairingMatrix("let-7a", "equilibrium", 0, sitelist="randthrp_comp")[1:nrow(out), 1:ncol(out)]
out6 <- MakePairingMatrix("let-7a", "equilibrium", 0, sitelist="randthrp_comp", sumseed=TRUE)[1:nrow(out), 1:ncol(out)]



plot(out3, out4, xlim=c(0, 2), ylim=c(0, 2))
points(out, out3, col="blue")
points(out, out4, col="red")
points(out, out5, col="purple")
points(out, out6, col="forestgreen")
break



sXc_old <- SitesXCounts(mir_prog, exp_prog, 3, sitelist_prog_old)

sXc_new <- SitesXCounts(mir_prog, exp_prog, 3, sitelist_prog_new)


# print(sXc_old["8mer-mmA2|-4|Supp", , drop=FALSE])
# print(sXc_new["8mer-mmA2|-4|Supp", , drop=FALSE])

kds_old <- EquilPars(mir_prog, exp_prog, 3, sitelist_prog_old, collapsemm=TRUE)

kds_new <- EquilPars(mir_prog, exp_prog, 3, sitelist_prog_new, sumseed=TRUE)

i_x <- intersect(rownames(kds_old), rownames(kds_new))

plot(kds_old[i_x, 2], kds_new[i_x, 2], log="xy")

# identify(kds_old[i_x, 2], kds_new[i_x, 2], labels=i_x)

break




test_mm_supp <- FitPairingOffsetAndMismatchModel(
  "let-7a", "equilibrium", sitelist="bipartite_random", corrected_kds=FALSE,
  intercept=FALSE, exponential=TRUE, kd_fc=TRUE, supp_base=TRUE,
  add_comp_sites=TRUE, collapsemm_addcomp=TRUE
)

break


test_alt <- GetEmpiricalBestSeedMismatchCoefficients(
  "miR-7-23nt", "equilibrium2_nb", sitelist="bipartite_random", corrected_kds=FALSE,
  supp_base=TRUE, kd_fc=TRUE, len_min=4, len_max=5, pos_3p_max=17, offsetmin=0, offsetmax=10, add_comp_sites=TRUE, collapsemm_addcomp=TRUE
)

print(test_alt)
break


test_mm_supp <- FitPairingOffsetAndMismatchModel(
  "miR-7-23nt", "equilibrium2_nb", sitelist="bipartite_random", corrected_kds=FALSE,
  supp_base=TRUE, add_comp_sites=TRUE, collapsemm_addcomp=TRUE, kd_fc=TRUE, len_min=4, len_max=5,
  pos_3p_max=17, offsetmin=0, offsetmax=10,
  intercept=FALSE, justmm=TRUE
)

test_mm <- FitPairingOffsetAndMismatchModel(
  "miR-7-23nt", "equilibrium2_nb", sitelist="bipartite_random", corrected_kds=FALSE,
  supp_base=TRUE, add_comp_sites=TRUE, collapsemm_addcomp=TRUE, kd_fc=TRUE, len_min=4, len_max=5,
  pos_3p_max=17, offsetmin=0, offsetmax=10, intercept=FALSE
)



dev.new(xpos=20, ypos=20, height=4, width=4)
plot(test_mm_supp$pairing, test_mm$pairing)

dev.new(xpos=420, ypos=20, height=4, width=4)
plot(test_mm_supp$offsets, test_mm$offsets)

break


break


df_1 <- df_all

df_use <- which(df_1$mm == "8mer")

print(head(df_1[df_use, ]))
test_new <- FitPairingOffsetAndMismatchModel(
  "miR-7-23nt", "equilibrium2_nb", sitelist="bipartite_random", corrected_kds=FALSE,
  supp_base=TRUE, add_comp_sites=TRUE, collapsemm_addcomp=TRUE, kd_fc=TRUE, len_min=4, len_max=5,
  pos_3p_max=17, offsetmin=0, offsetmax=6
)

df_2 <- df_all

df_use_2 <- which(df_2$mm == "8mer")

print(head(df_2[df_use_2, ]))

plot(df_2[df_use_2, 1], df_1[df_use, 1], xlim=c(-2, 2), ylim=c(-2, 2))


df_use <- which(df_1$mm == "8mer-mm")

print(head(df_1[df_use, ]))


# test_new <- FitPairingOffsetAndMismatchModel(
#   "miR-7-23nt", "equilibrium2_nb", sitelist="bipartite_random", corrected_kds=FALSE,
#   supp_base=TRUE, add_comp_sites=TRUE, kd_fc=TRUE, len_min=4, len_max=5,
#   pos_3p_max=17, offsetmin=0, offsetmax=10
# )

# df_2 <- df_all

df_use_2 <- which(df_2$mm == "8mer-mm")

print(head(df_2[df_use_2, ]))


dev.new(xpos=20, ypos=20, height=4, width=4)
plot(test_new$pairing, test$pairing)
dev.new(xpos=420, ypos=20, height=4, width=4)
plot(test_new$offsets, test$offsets)


break

# test <- FitPairingOffsetAndMismatchModel(
#   "let-7a", "equilibrium", sitelist="bipartite_random", corrected_kds=FALSE,
#   supp_base=TRUE, add_comp_sites=TRUE, collapsemm_addcomp=TRUE
# )


break

# test <- EquilPars("miR-124", "equilibrium", n_constant=3, sitelist="bipartite_random")

# print(head(test))
# break

# model_fit_0 <- FitPairingAndOffsetModel(
#   "let-7a", "equilibrium", kd_fc=TRUE, sitelist="bipartite_random_suppcomp", corrected_kds=FALSE, collapsemm=TRUE, use_global_df=FALSE, fixed_offset=0, F_method=FALSE, weights=TRUE, len_max=8
# )

# model_fit_8mer <- FitPairingAndOffsetModel(
#   "let-7a", "equilibrium", kd_fc=TRUE, sitelist="bipartite_random", site_base="8mer", corrected_kds=FALSE, use_global_df=FALSE, fixed_offset=0, F_method=FALSE, weights=FALSE, len_max=6, pos_3p_max=17, offsetmin=-2, offsetmax=6
# )

# model_fit_7merm8 <- FitPairingAndOffsetModel(
#   "let-7a", "equilibrium", kd_fc=TRUE, sitelist="bipartite_random", site_base="7mer-m8", corrected_kds=FALSE, use_global_df=FALSE, fixed_offset=0, F_method=FALSE, weights=FALSE, len_max=6, pos_3p_max=17, offsetmin=-2, offsetmax=6
# )

# model_fit_7merA1 <- FitPairingAndOffsetModel(
#   "let-7a", "equilibrium", kd_fc=TRUE, sitelist="bipartite_random", site_base="7mer-A1", corrected_kds=FALSE, use_global_df=FALSE, fixed_offset=0, F_method=FALSE, weights=TRUE, len_max=6, pos_3p_max=17, offsetmin=-2, offsetmax=6
# )

# model_fit_6mer <- FitPairingAndOffsetModel(
#   "let-7a", "equilibrium", kd_fc=TRUE, sitelist="bipartite_random", site_base="6mer", corrected_kds=FALSE, use_global_df=FALSE, fixed_offset=0, F_method=FALSE, weights=TRUE, len_max=6, pos_3p_max=17, offsetmin=-2, offsetmax=6
# )

# model_fit_6merm8 <- FitPairingAndOffsetModel(
#   "let-7a", "equilibrium", kd_fc=TRUE, sitelist="bipartite_random", site_base="6mer-m8", corrected_kds=FALSE, use_global_df=FALSE, fixed_offset=0, F_method=FALSE, weights=TRUE, len_max=6, pos_3p_max=17, offsetmin=-2, offsetmax=6
# )

# model_fit_6merA1 <- FitPairingAndOffsetModel(
#   "let-7a", "equilibrium", kd_fc=TRUE, sitelist="bipartite_random", site_base="6mer-A1", corrected_kds=FALSE, use_global_df=FALSE, fixed_offset=0, F_method=FALSE, weights=TRUE, len_max=6, pos_3p_max=17, offsetmin=-2, offsetmax=6
# )


# models <- list(model_fit_8mer, model_fit_7merm8, model_fit_7merA1, model_fit_6mer, model_fit_6merm8, model_fit_6merA1)

# diffs <- unlist(lapply(models, function(model) {
#   model_fit_coefs_shared <- intersect(rownames(model$coefs), rownames(model_fit_0$coefs))
#   coefs_supp <- model_fit_0$coefs[model_fit_coefs_shared, ]
#   coefs_model <- model$coefs[model_fit_coefs_shared, ]
#   rows_pairing <- grep("\\|", rownames(coefs_model), perl=TRUE, value=TRUE)

#   pars_supp <- coefs_supp[rows_pairing, 1]
#   pars_model <- coefs_model[rows_pairing, 1]
#   plot(pars_model, pars_supp)
#   out <- lm(pars_model ~ pars_supp - 1)

#   out$coefficients
# }))

# names(diffs) <- c("8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-m8", "6mer-A1")

# out <- GetEmpiricalBestSeedMismatchCoefficients("let-7a", "equilibrium", sitelist="bipartite_random", kd_fc=TRUE, corrected_kds=FALSE, supp_base=TRUE, len_max=6, offsetmin=-2, offsetmax=6)



coefs_all_supp <- FitPairingOffsetAndMismatchModel(
  "let-7a", "equilibrium", supp_base=TRUE, sitelist="bipartite_random",
  corrected_kds=FALSE, kd_fc=TRUE, use_global_df=FALSE, len_max=6,
  pos_3p_max=17, offsetmin=-2, offsetmax=6
)
kds <- EquilPars("let-7a", "equilibrium", n_constant=3, sitelist="bipartite_random")
# print(head(kds))
kds_use <- kds[sprintf("%s_Kd", names(diffs)), 2]
# print(kds_use)

dev.new(xpos=20, ypos=20, height=4, width=4)
plot(kds_use, exp(diffs), log="xy")

dev.new(xpos=420, ypos=20, height=4, width=4)
plot(kds_use, rowMeans(out[, 1:round(0.25*ncol(out))]), log="x")

dev.new(xpos=820, ypos=20, height=4, width=4)
image(t(out))

dev.new(xpos=20, ypos=420, height=4, width=4)
plot(kds_use, coefs_all_supp$mm[names(diffs)], log="x")


break
coefs_fit_0 <- model_fit_1$coefs[model_fit_coefs_shared, ]
coefs_fit   <- model_fit_2$coefs[model_fit_coefs_shared, ]



rows_pairing <- grep("\\|", rownames(coefs_fit_0), perl=TRUE, value=TRUE)
rows_coefs   <- grep("\\|", rownames(coefs_fit_0), perl=TRUE, value=TRUE, invert=TRUE)
dev.new(xpos=20, ypos=20, height=4, width=4)
par(mar=c(3, 3, 1, 1))
plot(coefs_fit_0[rows_pairing, 1], coefs_fit[rows_pairing, 1])

dev.new(xpos=420, ypos=20, height=4, width=4)
par(mar=c(3, 3, 1, 1))
plot(as.integer(rows_coefs), coefs_fit_0[rows_coefs, 1], type="o", ylim=c(0, 1))
lines(as.integer(rows_coefs), coefs_fit[rows_coefs, 1], type="o", col="red")

break

model_fit_2 <- FitPairingAndOffsetModel(
  "let-7a", "equilibrium", kd_fc=TRUE, sitelist="bipartite_random_suppcomp", supp_base=TRUE, corrected_kds=FALSE, use_global_df=FALSE, fixed_offset=0, F_method=TRUE, weights=TRUE
)

# model_fit_3 <- FitPairingAndOffsetModel(
#   "let-7a", "equilibrium", kd_fc=TRUE, sitelist="bipartite_random_suppcomp", supp_base=TRUE, corrected_kds=FALSE, use_global_df=FALSE, fixed_offset=0, F_method=TRUE, F_exact=FALSE, weights=TRUE
# )



# model_fit_2 <- FitPairingAndOffsetModel(
#   "let-7a-21nt", "equil_c2_nb", kd_fc=TRUE, use_global_df=TRUE, F_method=TRUE, fixed_offset=0
# )

# model_fit_3 <- FitPairingAndOffsetModel(
#   "let-7a-21nt", "equil_c2_nb", kd_fc=TRUE, use_global_df=TRUE, F_method=TRUE, F_exact=FALSE, fixed_offset=0
# )


dev.new(xpos=420, ypos=20, height=4, width=4)
par(mar=c(3, 3, 1, 1))

plot(-4:16, model_fit$coefs[58:78, 1], type="o", ylim=c(0, 1))
lines(-4:16, model_fit$coefs[58:78, 2], lty=2)
lines(-4:16, model_fit$coefs[58:78, 3], lty=2)

lines(-4:16, model_fit_2$coefs[58:78, 1], type="o", col="red")
lines(-4:16, model_fit_2$coefs[58:78, 2], lty=2, col="red")
lines(-4:16, model_fit_2$coefs[58:78, 3], lty=2, col="red")

lines(-4:16, model_fit_3$coefs[58:78, 1], type="o", col="blue")
lines(-4:16, model_fit_3$coefs[58:78, 2], lty=2, col="blue")
lines(-4:16, model_fit_3$coefs[58:78, 3], lty=2, col="blue")



break

segments(x0=model_fit$coefs[, 2], y0=model_fit_2$coefs[, 1], x1=model_fit$coefs[, 3])
segments(x0=model_fit$coefs[, 1], y0=model_fit_2$coefs[, 2], y1=model_fit_2$coefs[, 3])

# dev.new(xpos=420, ypos=20, height=4, width=4)
# par(mar=c(3, 3, 1, 1))
# plot(model_fit$coefs[, 1], model_fit_3$coefs[, 1])

# segments(x0=model_fit$coefs[, 2], y0=model_fit_3$coefs[, 1], x1=model_fit$coefs[, 3])
# segments(x0=model_fit$coefs[, 1], y0=model_fit_3$coefs[, 2], y1=model_fit_3$coefs[, 3])



break
model_fit <- FitPairingAndOffsetModelForAllOffsets(
  "let-7a", "equilibrium", sitelist="bipartite_random", corrected_kds=FALSE,
  kd_fc=TRUE, site_base="8mer"
)
break

model_fit <- FitPairingOffsetAndMismatchModel(
  "let-7a", "equilibrium", sitelist="bipartite_random", corrected_kds=FALSE,
  kd_fc=TRUE, supp_base=TRUE, offsetmin=-2, len_max=6
)

break


break

 # Testing the AIC with the thrp function

mirna <- "let-7a"
experiment <- "equilibrium"
sitelist <- "bipartite_random"


# out_base <- FitPairingAndOffsetModelForAllOffsets(mirna, experiment, sitelist=sitelist, corrected_kds=FALSE, len_max=7)

# out_8mer <- FitPairingAndOffsetModelForAllOffsets(mirna, experiment, sitelist=sitelist, corrected_kds=FALSE, site_base="8mer", len_max=7)
# out_7merm8 <- FitPairingAndOffsetModelForAllOffsets(mirna, experiment, sitelist=sitelist, corrected_kds=FALSE, site_base="7mer-m8", len_max=7)
# out_7merA1 <- FitPairingAndOffsetModelForAllOffsets(mirna, experiment, sitelist=sitelist, corrected_kds=FALSE, site_base="7mer-A1", len_max=8)
# out_6mer <- FitPairingAndOffsetModelForAllOffsets(mirna, experiment, sitelist=sitelist, corrected_kds=FALSE, site_base="6mer", len_max=8)
# out_6merA1 <- FitPairingAndOffsetModelForAllOffsets(mirna, experiment, sitelist=sitelist, corrected_kds=FALSE, site_base="6mer-A1", len_max=8)
# out_6merm8 <- FitPairingAndOffsetModelForAllOffsets(mirna, experiment, sitelist=sitelist, corrected_kds=FALSE, site_base="6mer-m8", len_max=8)



dev.new(xpos=20, ypos=20, height=4, width=4)
plot(c(out_base$pairing$MLE), c(out_6mer$pairing$MLE), col="lightblue", xlim=c(0, 3), ylim=c(0, 2.5))

points(c(out_base$pairing$MLE), c(out_8mer$pairing$MLE), col="purple")
points(c(out_base$pairing$MLE), c(out_7merm8$pairing$MLE), col="red")
points(c(out_base$pairing$MLE), c(out_7merA1$pairing$MLE), col="blue")
points(c(out_base$pairing$MLE), c(out_6mer$pairing$MLE), col="cyan")
points(c(out_base$pairing$MLE), c(out_6merA1$pairing$MLE), col="lightblue")
points(c(out_base$pairing$MLE), c(out_6merm8$pairing$MLE), col="pink")

dev.new(xpos=420, ypos=20, height=4, width=4)
plot(-2:16, out_base$offsets[, 1], type="o", ylim=c(0, 1))
lines(-2:16, out_8mer$offsets[, 1], type="o", col="purple")
lines(-2:16, out_7merm8$offsets[, 1], type="o", col="red")
lines(-2:16, out_7merA1$offsets[, 1], type="o", col="blue")
# lines(-2:16, out_6mer$offsets[, 1], type="o", col="cyan")
lines(-2:16, out_6merA1$offsets[, 1], type="o", col="lightblue")
# lines(-2:16, out_6merm8$offsets[, 1], type="o", col="pink")

break



AIC_1 <- GetModelAIC("let-7a-21nt", "equil_c2_nb")

AIC_2 <- GetModelAIC("let-7a-21nt", "equil_c2_nb", modelnew=TRUE, xpos=420)

AIC_min <- min(AIC_1, AIC_2)

AIC_all <- c(AIC_1, AIC_2)

print(exp((AIC_min - AIC_all)/2))

SSEs <- c(sum( (data_fit - model_sim_old)^2), sum( (data_fit - model_sim_new)^2))
n <- length(data_fit)
sigma2s <- 1/n*SSEs
logLs <- -n/2*( log(2*pi*sigma2s) + 1)
k <- c(139, 478)
AICs <- 2*k - 2*logLs

probs <- exp((min(AICs) - AICs)/2)

print(SSEs)
print(sigma2s)
print(logLs)
print(AICs)
print(AIC_all)
print(probs)

break


out_mat_1 <- GetEmpiricalBestSeedMismatchCoefficients("let-7a-21nt", "equil_c2_nb")


out_mat_2 <- GetEmpiricalBestSeedMismatchCoefficients("let-7a", "equilibrium", sitelist="bipartite_random", corrected_kds=FALSE)


dev.new(xpos=20, ypos=20, height=4, width=4)
image(t(out_mat_1))

dev.new(xpos=420, ypos=20, height=4, width=4)
image(t(out_mat_2))

break

# for (len_i in 4:11) {
#   for (pos_j in seq(9, nchar(kMirnaSeqs["let-7a-21nt"]) - len_i + 1)) {
#     test <- log10(1/MakeMismatchByOffsetMatrix("let-7a-21nt", "equil_c2_nb", len_i, pos_j, kd_fc=TRUE))
#     ind_use <- which(colMeans(test, na.rm=TRUE) == max(colMeans(test, na.rm=TRUE), na.rm=TRUE))
#     test_use <- test[, ind_use]
#     data_out_prog <- cbind(data_out, test_use - mean(test_use, na.rm=TRUE))
#     colnames(data_out_prog)[ncol(data_out_prog)] <- sprintf("l%sp%s", len_i, pos_j)
#   }
# }

data_out_prog_order_cols <- order(apply(data_out_prog[, 3:ncol(data_out_prog)], 2, sd, na.rm=TRUE))

print(data_out_prog_order_cols)

data_out_prog_new <- cbind(data_out_prog[, 1:2], data_out_prog[, rev(data_out_prog_order_cols) + 2])

data_out_prog_new <- data_out_prog_new[, which(!is.na(colMeans(data_out_prog_new)))]
dev.new(xpos=20, ypos=20, height=4, width=4)
image(t(data_out_prog_new))


end_stop <- 10

data_out_prog_new_ave <- rowMeans(data_out_prog_new[, 3:end_stop], na.rm=TRUE)

dev.new(xpos=420, ypos=20, height=4, width=4)
plot(data_out_prog_new[, 1], data_out_prog_new_ave)

dev.new(xpos=820, ypos=20, height=4, width=4)
plot(data_out_prog_new_ave, coefs_p_o$mm - mean(coefs_p_o$mm, na.rm=TRUE))

dev.new(xpos=1220, ypos=20, height=4, width=4)
plot(data_out_prog_new_ave, rowMeans(coefs_p_n$offsetXmm) - mean(rowMeans(coefs_p_n$offsetXmm), na.rm=TRUE))


# data_out_rand <- matrix(log10(c(kds_mm_prog[, 2], kds_mm_rand[, 2])), nrow=nrow(kds_mm_prog), ncol=2,
#   dimnames=list(GetAll8merMmSites("let-7a"), c("prog", "rand")))


# for (len_i in 4:11) {
#   print(len_i)
#   for (pos_j in seq(9, nchar(kMirnaSeqs["let-7a-21nt"]) - len_i + 1)) {
#     print(pos_j)
#     test <- log10(1/MakeMismatchByOffsetMatrix("let-7a", "equilibrium", len_i, pos_j, sitelist="bipartite_random", kd_fc=TRUE, corrected_kds=FALSE))
#     ind_use <- which(colMeans(test, na.rm=TRUE) == max(colMeans(test, na.rm=TRUE), na.rm=TRUE))
#     test_use <- test[, ind_use]
#     data_out_rand <- cbind(data_out_rand, test_use - mean(test_use, na.rm=TRUE))
#     colnames(data_out_rand)[ncol(data_out_rand)] <- sprintf("l%sp%s", len_i, pos_j)
#   }
# }

data_out_rand_order_cols <- order(apply(data_out_rand[, 3:ncol(data_out_rand)], 2, sd, na.rm=TRUE))

print(data_out_rand_order_cols)

data_out_rand_new <- cbind(data_out_rand[, 1:2], data_out_rand[, rev(data_out_rand_order_cols) + 2])

data_out_rand_new <- data_out_rand_new[, which(!is.na(colMeans(data_out_rand_new)))]
dev.new(xpos=20, ypos=420, height=4, width=4)
image(t(data_out_rand_new))


end_stop <- 6

data_out_rand_new_ave <- rowMeans(data_out_rand_new[, 3:end_stop], na.rm=TRUE)



dev.new(xpos=420, ypos=420, height=4, width=4)
plot(data_out_rand_new[, 2], data_out_rand_new_ave)

dev.new(xpos=820, ypos=420, height=4, width=4)
plot(data_out_rand_new_ave, coefs_p_o$mm - mean(coefs_p_o$mm, na.rm=TRUE))

dev.new(xpos=1220, ypos=420, height=4, width=4)
plot(data_out_rand_new_ave, rowMeans(coefs_p_n$offsetXmm) - mean(rowMeans(coefs_p_n$offsetXmm), na.rm=TRUE))

dev.new(xpos=1520, ypos=420, height=4, width=4)
plot(data_out_prog_new_ave, data_out_rand_new_ave)

dev.new(xpos=20, ypos=820, height=4, width=4)
plot(coefs_p_o$mm, coefs_r_o$mm)


break

graphics.off()
dev.new(xpos=20, ypos=20, height=4, width=4)

plot(names(coefs_p_o$offsets), coefs_p_o$offsets)

offXmm_o_mat <- coefs_p_o$offsets %o% coefs_p_o$mm

dev.new(xpos=20, ypos=20, height=4, width=4)
image(offXmm_o_mat)

dev.new(xpos=420, ypos=20, height=4, width=4)
image(t(coefs_p_n$offsetXmm))

dev.new(xpos=820, ypos=20, height=4, width=4)
plot(c(offXmm_o_mat), c(t(coefs_p_n$offsetXmm)))


dev.new(xpos=20, ypos=420, height=4, width=4)
image(coefs_p_o$pairing)

dev.new(xpos=420, ypos=420, height=4, width=4)
image(coefs_p_n$pairing)

dev.new(xpos=820, ypos=420, height=4, width=4)
plot(c(coefs_p_o$pairing), c(coefs_p_n$pairing))



break


# coefs_rand_old <- FitPairingOffsetAndMismatchModel("let-7a", "equilibrium", sitelist="bipartite_random", corrected_kds=FALSE, kd_fc=FALSE)
# coefs_rand_new <- FitPairingOffsetAndMismatchModelNew("let-7a", "equilibrium", sitelist="bipartite_random", corrected_kds=FALSE, kd_fc=FALSE)

# coefs_prog_old_fc <- FitPairingOffsetAndMismatchModel("let-7a-21nt", "equil_c2_nb", sitelist="programmed")
# coefs_prog_new_fc <- FitPairingOffsetAndMismatchModelNew("let-7a-21nt", "equil_c2_nb", sitelist="programmed")

# coefs_rand_old_fc <- FitPairingOffsetAndMismatchModel("let-7a", "equilibrium", sitelist="bipartite_random", corrected_kds=FALSE)
# coefs_rand_new_fc <- FitPairingOffsetAndMismatchModelNew("let-7a", "equilibrium", sitelist="bipartite_random", corrected_kds=FALSE)



# prog_data_fc <- MakeFullNucleotideAndMisMatchContributionDf(
#   "let-7a-21nt", "equil_c2_nb", n_constant=3, sitelist="programmed",
#   offsetmin=-4, offsetmax=16, corrected_kds=TRUE)

# prog_data <- MakeFullNucleotideAndMisMatchContributionDf(
#   "let-7a-21nt", "equil_c2_nb", n_constant=3, sitelist="programmed",
#   offsetmin=-4, offsetmax=16, corrected_kds=TRUE, kd_fc=FALSE)

# rand_data_fc <- MakeFullNucleotideAndMisMatchContributionDf(
#   "let-7a", "equilibrium", n_constant=3, sitelist="bipartite_random",
#   offsetmin=-4, offsetmax=16, corrected_kds=FALSE)

# rand_data <- MakeFullNucleotideAndMisMatchContributionDf(
#   "let-7a", "equilibrium", n_constant=3, sitelist="bipartite_random",
#   offsetmin=-4, offsetmax=16, corrected_kds=FALSE, kd_fc=FALSE)


pairing <- coefs_rand_old_fc$pairing
offsets <- coefs_rand_old_fc$offsets
mm <- coefs_rand_old_fc$mm
base <- coefs_rand_old_fc$base

model_sim <- apply(rand_data_fc, 1, function(row) {
  pairing[as.character(as.integer(row[4])), as.character(as.integer(row[3]))]*mm[row[7]]*offsets[row[6]] + base
})


model_sim_no_mm <- apply(rand_data_fc, 1, function(row) {
  pairing[as.character(as.integer(row[4])), as.character(as.integer(row[3]))]*offsets[row[6]] + base
})




dev.new(xpos=20, ypos=20, height=4, width=4)
plot(rand_data_fc$logkd, model_sim_no_mm, xlim=c(-0.5, 3), ylim=c(-0.5, 3))
points(rand_data_fc$logkd, model_sim, col="blue")
text(-0.5, 2, labels=cor(rand_data_fc$logkd, model_sim_no_mm)^2, adj=0)
text(-0.5, 1.5, labels=cor(rand_data_fc$logkd, model_sim)^2, adj=0)


pairing <- coefs_rand_new_fc$pairing
offsets <- coefs_rand_new_fc$offsets
mm <- coefs_rand_new_fc$mm

model_sim <- apply(rand_data_fc, 1, function(row) {
  pairing[as.character(as.integer(row[4])), as.character(as.integer(row[3]))]*offsets[row[6]] + mm[row[7]]
})

model_sim_no_mm <- apply(rand_data_fc, 1, function(row) {
  pairing[as.character(as.integer(row[4])), as.character(as.integer(row[3]))]*offsets[row[6]]
})




dev.new(xpos=420, ypos=20, height=4, width=4)
plot(rand_data_fc$logkd, model_sim_no_mm, xlim=c(-0.5, 3), ylim=c(-0.5, 3))
points(rand_data_fc$logkd, model_sim, col="blue")
text(-0.5, 1.5, labels=cor(rand_data_fc$logkd, model_sim_no_mm)^2, adj=0)
text(-0.5, 1.25, labels=cor(rand_data_fc$logkd, model_sim)^2, adj=0)

pairing <- coefs_rand_old$pairing
offsets <- coefs_rand_old$offsets
mm <- coefs_rand_old$mm
base <- coefs_rand_old$base

model_sim <- apply(rand_data, 1, function(row) {
  pairing[as.character(as.integer(row[4])), as.character(as.integer(row[3]))]*mm[row[7]]*offsets[row[6]] + base
})


model_sim_no_mm <- apply(rand_data, 1, function(row) {
  pairing[as.character(as.integer(row[4])), as.character(as.integer(row[3]))]*offsets[row[6]] + base
})




dev.new(xpos=820, ypos=20, height=4, width=4)
plot(rand_data$logkd, model_sim_no_mm, xlim=c(0, 4), ylim=c(0, 4))
points(rand_data$logkd, model_sim, col="blue")
text(0, 3, labels=cor(rand_data$logkd, model_sim_no_mm)^2, adj=0)
text(0, 2.5, labels=cor(rand_data$logkd, model_sim)^2, adj=0)



pairing <- coefs_rand_new$pairing
offsets <- coefs_rand_new$offsets
mm <- coefs_rand_new$mm

model_sim <- apply(rand_data, 1, function(row) {
  pairing[as.character(as.integer(row[4])), as.character(as.integer(row[3]))]*offsets[row[6]] + mm[row[7]]
})

model_sim_no_mm <- apply(rand_data, 1, function(row) {
  pairing[as.character(as.integer(row[4])), as.character(as.integer(row[3]))]*offsets[row[6]]
})


dev.new(xpos=1220, ypos=20, height=4, width=4)
plot(rand_data$logkd, model_sim_no_mm, xlim=c(0, 4), yli=c(0, 4))
points(rand_data$logkd, model_sim, col="blue")
text(0, 3, labels=cor(rand_data$logkd, model_sim_no_mm)^2, adj=0)
text(0, 2.5, labels=cor(rand_data$logkd, model_sim)^2, adj=0)



# dev.new(xpos=20, ypos=520, height=4, width=4)
# plot(1/(10^coefs_rand_old$mm), kds_mm_rand[sprintf("%s_Kd", names(coefs_rand_old$mm)), 2], log="xy", xlim=c(1e-3, 1), ylim=c(1e-3, 1))
# segments(x0=1e-3, y0=1e-3, x1=1, y1=1, lty=2)

# dev.new(xpos=420, ypos=520, height=4, width=4)
# plot(1/(10^coefs_rand_new$mm), kds_mm_rand[sprintf("%s_Kd", names(coefs_rand_old$mm)), 2], log="xy", xlim=c(1e-3, 1), ylim=c(1e-3, 1))
# segments(x0=1e-3, y0=1e-3, x1=1, y1=1, lty=2)

# dev.new(xpos=820, ypos=520, height=4, width=4)
# plot(1/(10^coefs_rand_old_fc$mm), kds_mm_rand[sprintf("%s_Kd", names(coefs_rand_old$mm)), 2], log="xy", xlim=c(1e-3, 1), ylim=c(1e-3, 1))
# segments(x0=1e-3, y0=1e-3, x1=1, y1=1, lty=2)

# dev.new(xpos=1220, ypos=520, height=4, width=4)
# plot(1/(10^coefs_rand_new_fc$mm), kds_mm_rand[sprintf("%s_Kd", names(coefs_rand_old$mm)), 2], log="xy", xlim=c(1e-3, 1), ylim=c(1e-3, 1))
# segments(x0=1e-3, y0=1e-3, x1=1, y1=1, lty=2)


# break

dev.new(xpos=20, ypos=420, height=4, width=4)
plot(1/(10^coefs_rand_old_fc$mm), kds_mm_rand[sprintf("%s_Kd", names(coefs_rand_old$mm)), 2], log="xy", xlim=c(1e-3, 1), ylim=c(1e-3, 1))
segments(x0=1e-3, y0=1e-3, x1=1, y1=1, lty=2)

dev.new(xpos=420, ypos=420, height=4, width=4)
plot(1/(10^coefs_rand_new_fc$mm), kds_mm_rand[sprintf("%s_Kd", names(coefs_rand_old$mm)), 2], log="xy", xlim=c(1e-3, 1), ylim=c(1e-3, 1))
segments(x0=1e-3, y0=1e-3, x1=1, y1=1, lty=2)

dev.new(xpos=820, ypos=420, height=4, width=4)
plot(1/(10^coefs_rand_old$mm), kds_mm_rand[sprintf("%s_Kd", names(coefs_rand_old$mm)), 2], log="xy", xlim=c(1e-3, 1), ylim=c(1e-3, 1))
segments(x0=1e-3, y0=1e-3, x1=1, y1=1, lty=2)

dev.new(xpos=1220, ypos=420, height=4, width=4)
plot(1/(10^coefs_rand_new$mm), kds_mm_rand[sprintf("%s_Kd", names(coefs_rand_old$mm)), 2], log="xy", xlim=c(1e-3, 1), ylim=c(1e-3, 1))
segments(x0=1e-3, y0=1e-3, x1=1, y1=1, lty=2)


# dev.new(xpos=20, ypos=520, height=4, width=4)
# plot(1/(10^coefs_rand_old$mm), kds_mm_rand[sprintf("%s_Kd", names(coefs_rand_old$mm)), 2], log="xy", xlim=c(1e-3, 1), ylim=c(1e-3, 1))
# segments(x0=1e-3, y0=1e-3, x1=1, y1=1, lty=2)

# dev.new(xpos=420, ypos=520, height=4, width=4)
# plot(1/(10^coefs_rand_new$mm), kds_mm_rand[sprintf("%s_Kd", names(coefs_rand_old$mm)), 2], log="xy", xlim=c(1e-3, 1), ylim=c(1e-3, 1))
# segments(x0=1e-3, y0=1e-3, x1=1, y1=1, lty=2)

# dev.new(xpos=820, ypos=520, height=4, width=4)
# plot(1/(10^coefs_rand_old_fc$mm), kds_mm_rand[sprintf("%s_Kd", names(coefs_rand_old$mm)), 2], log="xy", xlim=c(1e-3, 1), ylim=c(1e-3, 1))
# segments(x0=1e-3, y0=1e-3, x1=1, y1=1, lty=2)

# dev.new(xpos=1220, ypos=520, height=4, width=4)
# plot(1/(10^coefs_rand_new_fc$mm), kds_mm_rand[sprintf("%s_Kd", names(coefs_rand_old$mm)), 2], log="xy", xlim=c(1e-3, 1), ylim=c(1e-3, 1))
# segments(x0=1e-3, y0=1e-3, x1=1, y1=1, lty=2)


break
dev.new(xpos=1020, ypos=20, height=5, width=5)
plot(coefs_prog$mm, kds_mm_prog[sprintf("%s_Kd", names(coefs_rand$mm)), 2], log="y")



# kds_prog_mm

# coefs_rand_10 <- FitPairingOffsetAndMismatchModel("let-7a", "equilibrium", sitelist="bipartite_random", len_max=10, corrected_kds=FALSE)
# coefs_rand_9 <- FitPairingOffsetAndMismatchModel("let-7a", "equilibrium", sitelist="bipartite_random", len_max=9, corrected_kds=FALSE)
# coefs_rand_8 <- FitPairingOffsetAndMismatchModel("let-7a", "equilibrium", sitelist="bipartite_random", len_max=8, corrected_kds=FALSE)
# coefs_rand_7 <- FitPairingOffsetAndMismatchModel("let-7a", "equilibrium", sitelist="bipartite_random", len_max=7, corrected_kds=FALSE)
# coefs_rand_6 <- FitPairingOffsetAndMismatchModel("let-7a", "equilibrium", sitelist="bipartite_random", len_max=6, corrected_kds=FALSE)


# plot(-4:16, coefs_rand_suppcomp_11$offsets[, 1], type="o")
# lines(-4:16, coefs_rand_suppcomp_10$offsets[, 1], col="red")
# lines(-4:16, coefs_rand_suppcomp_9$offsets[, 1], col="orangered")
# lines(-4:16, coefs_rand_suppcomp_8$offsets[, 1], col="forestg
#   reen")
# lines(-4:16, coefs_rand_suppcomp_7$offsets[, 1], col="blue")
# lines(-4:16, coefs_rand_suppcomp_6$offsets[, 1], col="purple")
# lines(-4:16, coefs_rand_suppcomp$offsets[, 1])

break
plot(coefs_prog$mm, coefs_rand$mm)


break



kds_prog <- ApplyKdCorrection("let-7a-21nt", "equil_c2_nb", prog_n_constant=3, prog_sitelist="programmed",
                              rand_n_constant=3, rand_sitelist="bipartite_random")


kds_rand <- EquilPars("let-7a", "equilibrium", n_constant=3, sitelist="bipartite_random")

kds_mm_prog <- kds_prog[sprintf("%s_Kd", GetAll8merMmSites("let-7a")), ]

kds_mm_rand <- kds_rand[sprintf("%s_Kd", GetAll8merMmSites("let-7a")), ]


site_len <- "5mer-m11.15"


dev.new(xpos=20, ypos=20, height=4, width=4)
plot(kds_mm_prog[, 2], kds_mm_rand[, 2], log="xy", xlim=c(1e-3, 1), ylim=c(1e-3, 1))
dev.new(xpos=420, ypos=20, height=4, width=4)
rows_use <- sprintf("%s|13|%s_Kd", site_len, GetAll8merMmSites("let-7a"))
plot(kds_prog[rows_use, 2]/kds_mm_prog[, 2], kds_mm_prog[, 2], log="xy", xlim=c(1e-3, 1), ylim=c(1e-3, 1))

dev.new(xpos=820, ypos=20, height=4, width=4)
rows_use <- sprintf("%s|14|%s_Kd", site_len, GetAll8merMmSites("let-7a"))
plot(kds_prog[rows_use, 2]/kds_mm_prog[, 2], kds_mm_prog[, 2], log="xy", xlim=c(1e-3, 1), ylim=c(1e-3, 1))

dev.new(xpos=1220, ypos=20, height=4, width=4)
rows_use <- sprintf("%s|15|%s_Kd", site_len, GetAll8merMmSites("let-7a"))
plot(kds_prog[rows_use, 2]/kds_mm_prog[, 2], kds_mm_prog[, 2], log="xy", xlim=c(1e-3, 1), ylim=c(1e-3, 1))

dev.new(xpos=420, ypos=420, height=4, width=4)
rows_use <- sprintf("%s|13|%s_Kd", site_len, GetAll8merMmSites("let-7a"))
plot(kds_rand[rows_use, 2]/kds_mm_rand[, 2], kds_mm_rand[, 2], log="xy", xlim=c(1e-3, 1), ylim=c(1e-3, 1))

dev.new(xpos=820, ypos=420, height=4, width=4)
rows_use <- sprintf("%s|14|%s_Kd", site_len, GetAll8merMmSites("let-7a"))
plot(kds_rand[rows_use, 2]/kds_mm_rand[, 2], kds_mm_rand[, 2], log="xy", xlim=c(1e-3, 1), ylim=c(1e-3, 1))

dev.new(xpos=1220, ypos=420, height=4, width=4)
rows_use <- sprintf("%s|15|%s_Kd", site_len, GetAll8merMmSites("let-7a"))
plot(kds_rand[rows_use, 2]/kds_mm_rand[, 2], kds_mm_rand[, 2], log="xy", xlim=c(1e-3, 1), ylim=c(1e-3, 1))





break


# offsets_16 <- FitPairingAndOffsetModel("miR-155", "equil_sc_nb", fixed_offset=16)$pairing


# plot(c(offsets_4$MLE), c(offsets_16$MLE))

break


plot(rownames(offsets_4), offsets_4[, 1])

segments(x0=as.integer(rownames(offsets_4)), y0=offsets_4[, 2], y1=offsets_4[, 3])
lines(rownames(offsets_4), offsets_all[, 1], col="red")
segments(x0=as.integer(rownames(offsets_all))+ 0.1, y0=offsets_all[, 2], y1=offsets_all[, 3], col="red")


break


mm1_kds <- EquilPars("let-7a-21nt", "equil_c2_nb", 0, "programmed_suppcomp", start_mm=10, stop_mm=20)
# mm2_kds <- EquilPars("miR-155", "equil_sc_nb", 0, "programmed_suppcomp", start_mm=13, stop_mm=22)


kds_9mer <- wt_kds[grep("^10mer-m10\\.19\\|.*\\|Comp_Kd", rownames(wt_kds), perl=TRUE), 2, drop=FALSE]
kds_11mer <- wt_kds[grep("^11mer-m10\\.20\\|.*\\|Comp_Kd", rownames(wt_kds), perl=TRUE), 2, drop=FALSE]

kds_11mermm1 <- mm1_kds[grep("^11mer-m10\\.20\\|.*\\|Comp_Kd", rownames(mm1_kds), perl=TRUE), 2, drop=FALSE]
# kds_11mermm2 <- mm2_kds[grep("^11mer-m13\\.23\\|.*\\|Comp_Kd", rownames(mm2_kds), perl=TRUE), 2, drop=FALSE]

kds_11mermm1_A20 <- mm1_kds[grep("^11mer-m10\\.20mmA20\\|.*\\|Comp_Kd", rownames(mm1_kds), perl=TRUE), 2, drop=FALSE]
kds_11mermm1_G20 <- mm1_kds[grep("^11mer-m10\\.20mmG20\\|.*\\|Comp_Kd", rownames(mm1_kds), perl=TRUE), 2, drop=FALSE]
kds_11mermm1_T20 <- mm1_kds[grep("^11mer-m10\\.20mmT20\\|.*\\|Comp_Kd", rownames(mm1_kds), perl=TRUE), 2, drop=FALSE]

kds_11mermm1_A19 <- mm1_kds[grep("^11mer-m10\\.20mmA19\\|.*\\|Comp_Kd", rownames(mm1_kds), perl=TRUE), 2, drop=FALSE]
kds_11mermm1_C19 <- mm1_kds[grep("^11mer-m10\\.20mmC19\\|.*\\|Comp_Kd", rownames(mm1_kds), perl=TRUE), 2, drop=FALSE]
kds_11mermm1_G19 <- mm1_kds[grep("^11mer-m10\\.20mmG19\\|.*\\|Comp_Kd", rownames(mm1_kds), perl=TRUE), 2, drop=FALSE]

kds_11mermm1_C18 <- mm1_kds[grep("^11mer-m10\\.20mmC18\\|.*\\|Comp_Kd", rownames(mm1_kds), perl=TRUE), 2, drop=FALSE]
kds_11mermm1_G18 <- mm1_kds[grep("^11mer-m10\\.20mmG18\\|.*\\|Comp_Kd", rownames(mm1_kds), perl=TRUE), 2, drop=FALSE]
kds_11mermm1_T18 <- mm1_kds[grep("^11mer-m10\\.20mmT18\\|.*\\|Comp_Kd", rownames(mm1_kds), perl=TRUE), 2, drop=FALSE]

kds_11mermm1_A10 <- mm1_kds[grep("^11mer-m10\\.20mmA10\\|.*\\|Comp_Kd", rownames(mm1_kds), perl=TRUE), 2, drop=FALSE]
kds_11mermm1_C10 <- mm1_kds[grep("^11mer-m10\\.20mmC10\\|.*\\|Comp_Kd", rownames(mm1_kds), perl=TRUE), 2, drop=FALSE]
kds_11mermm1_G10 <- mm1_kds[grep("^11mer-m10\\.20mmG10\\|.*\\|Comp_Kd", rownames(mm1_kds), perl=TRUE), 2, drop=FALSE]


# print(kds_11mer)
# print(kds_10mer)
print(kds_11mermm1)
# print(kds_11mermm2)

# print(out_mat[, 1, drop=FALSE])
# print(out_mat_average[, 1, drop=FALSE])
dev.new(xpos=20, ypos=20, height=4, width=4)
plot(1:nrow(kds_11mer), rev(kds_11mermm1[, 1]), type="l", log="y", ylim=c(1e-4, 1))
lines(1:nrow(kds_11mer), rev(kds_11mermm1_A20[, 1]), col="red")
lines(1:nrow(kds_11mer), rev(kds_11mermm1_G20[, 1]), col="red")
lines(1:nrow(kds_11mer), rev(kds_11mermm1_T20[, 1]), col="red")

lines(1:nrow(kds_11mer), rev(kds_11mermm1_A19[, 1]), col="violet")
lines(1:nrow(kds_11mer), rev(kds_11mermm1_C19[, 1]), col="violet")
lines(1:nrow(kds_11mer), rev(kds_11mermm1_G19[, 1]), col="violet")

# lines(1:nrow(kds_11mer), rev(kds_11mermm1_C18[, 1]), col="green")
# lines(1:nrow(kds_11mer), rev(kds_11mermm1_G18[, 1]), col="green")
# lines(1:nrow(kds_11mer), rev(kds_11mermm1_T18[, 1]), col="green")

lines(1:nrow(kds_11mer), rev(kds_11mermm1_A10[, 1]), col="forestgreen")
lines(1:nrow(kds_11mer), rev(kds_11mermm1_C10[, 1]), col="forestgreen")
lines(1:nrow(kds_11mer), rev(kds_11mermm1_G10[, 1]), col="forestgreen")

break

# lines(1:(nrow(kds_11mer)), rev(out_mat[, 1]), lty=3, lwd=3)

# segments(x0=1, y0=0.003246591, x1=nrow(kds_11mer), col=rgb(0, 0, 0, alpha=0.5))

lines(2:(nrow(kds_11mer) - 1), rev(out_mat_average[, 1]))
lines(2:(nrow(kds_11mer) - 1), rev(out_mat_average[, 2]), lty=3, col="green")
lines(2:(nrow(kds_11mer) - 1), rev(out_mat_average[, 3]), lty=3, col="green")
lines(2:(nrow(kds_11mer) - 1), rev(out_mat_average[, 4]), lty=3, col="green")
# lines(2:(nrow(kds_11mer) - 1), rev(1/rowMeans(1/out_mat_average[, 2:4])), col="green")

lines(2:(nrow(kds_11mer) - 1), rev(out_mat_average[, 5]), lty=3, col="purple")
lines(2:(nrow(kds_11mer) - 1), rev(out_mat_average[, 6]), lty=3, col="purple")
lines(2:(nrow(kds_11mer) - 1), rev(out_mat_average[, 7]), lty=3, col="purple")

lines(2:(nrow(kds_11mer) - 1), rev(out_mat_average[, 8]), lty=3, col="orangered ")
lines(2:(nrow(kds_11mer) - 1), rev(out_mat_average[, 9]), lty=3, col="orangered ")
lines(2:(nrow(kds_11mer) - 1), rev(out_mat_average[, 10]), lty=3, col="orangered  ")

lines(2:(nrow(kds_11mer) - 1), rev(out_mat_average[, 11]), lty=3, col="orangered  ")
lines(2:(nrow(kds_11mer) - 1), rev(out_mat_average[, 12]), lty=3, col="orangered  ")
lines(2:(nrow(kds_11mer) - 1), rev(out_mat_average[, 13]), lty=3, col="orangered  ")
break

lines(2:(nrow(kds_11mer) - 1), rev(1/rowMeans(1/out_mat_average[, 5:7])), col="purple")
lines(2:(nrow(kds_11mer) - 1), rev(1/rowMeans(1/out_mat_average[, 8:10])), col="purple")
lines(2:(nrow(kds_11mer) - 1), rev(1/rowMeans(1/out_mat_average[, 11:13])), col="purple")
lines(2:(nrow(kds_11mer) - 1), rev(1/rowMeans(1/out_mat_average[, 14:16])), col="purple")
lines(2:(nrow(kds_11mer) - 1), rev(1/rowMeans(1/out_mat_average[, 17:19])), col="purple")



check <- GetBestAverageOffsetMismatchKds("miR-155", "equil_sc_nb", 0, 13, 22)
lines(2:(nrow(kds_11mer)), rev(out_mat_average[, 1]), col="forestgreen")


# break
# lines(1:nrow(kds_11mer), rev(kds_10mer[2:nrow(kds_10mer), 1]), col="blue")
# lines(1:nrow(kds_11mer), rev(kds_11mermm1[, 1]), col="red")
# lines(1:nrow(kds_11mer), rev(kds_11mermm2[, 1]), col="purple")

# lines(1:nrow(kds_11mer), rev(kds_11mermm1_C[, 1]), col="green")
# lines(1:nrow(kds_11mer), rev(kds_11mermm1_G[, 1]), col="green")
# lines(1:nrow(kds_11mer), rev(kds_11mermm1_T[, 1]), col="green")



break

message("wt kds:")
print(wt_kds["11mer-m13.23|14|Comp_Kd", ])
print(mm1_kds["11mer-m13.23|14|Comp_Kd", ])
print(mm2_kds["11mer-m13.23|14|Comp_Kd", ])

print(wt_kds["10mer-m13.22|14|Comp_Kd", ])
print(mm1_kds["10mer-m13.22|14|Comp_Kd", ])
print(mm2_kds["10mer-m13.22|14|Comp_Kd", ])

print(mm1_kds["11mer-m13.23mmC23|14|Comp_Kd", ])
print(mm1_kds["11mer-m13.23mmG23|14|Comp_Kd", ])
print(mm1_kds["11mer-m13.23mmT23|14|Comp_Kd", ])


break

print(bg_kds["11mer-m10.20mmA20|14|Comp_Kd", ])
print(test_kds["11mer-m10.20mmA20|14|Comp_Kd", ])

print(bg_kds["11mer-m10.20mmG20|14|Comp_Kd", ])
print(test_kds["11mer-m10.20mmG20|14|Comp_Kd", ])

print(bg_kds["11mer-m10.20mmT20|14|Comp_Kd", ])
print(test_kds["11mer-m10.20mmT20|14|Comp_Kd", ])

print(bg_kds["10mer-m10.19|14|Comp_Kd", ])
print(test_kds["10mer-m10.19|14|Comp_Kd", ])

print(bg_kds["10mer-m10.19|14|Comp_Kd", ])

break
print(test_kds["11mer-m10.20mmA10|14|Comp_Kd", ])
print(test_kds["11mer-m10.20mmC10|14|Comp_Kd", ])
print(test_kds["11mer-m10.20mmG10|14|Comp_Kd", ])
print(test_kds["10mer-m11.20|15|Comp_Kd", ])
print(bg_kds["10mer-m11.20|15|Comp_Kd", ])




break

test_counts_0 <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a-21nt/equil_c2_nb/programmed_site_counts/I_0_m11.16mmsd_suppcomp_test.txt", row.names=1)
test_counts_5 <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a-21nt/equil_c2_nb/programmed_site_counts/I_5_m11.16mmsd_suppcomp_test.txt", row.names=1)

dev.new(xpos=20, ypos=20, height=5, width=5)
sites_use <- intersect(rownames(test_counts_0), rownames(test_counts_5))
plot(test_counts_0[sites_use, 1], test_counts_5[sites_use, 1], log="xy")
identify(test_counts_0[sites_use, 1], test_counts_5[sites_use, 1], labels=sites_use)


# print(head(test_kds))

# print(test_kds["6mer-m11.16mmT16|15|Comp_Kd",])
# print(test_kds["6mer-m11.16mmC16|15|Comp_Kd",])
# print(test_kds["6mer-m11.16mmG16|15|Comp_Kd",])
# print(test_kds["5mer-m11.15|15|Comp_Kd",])

break



library(deSolve)
graphics.off()
dev.new(xpos=20, ypos=20, height=5, width=5)
counts_old <- read.table(file="/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a-21nt/equil_c2_nb/programmed_site_counts/I_5_test.txt")
counts_new <- read.table(file="/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a-21nt/equil_c2_nb/programmed_site_counts/I_5_new_test.txt")

counts_old <- read.table(file="/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/bipartite_site_counts/I_5_test.txt", row.names=1)
counts_new <- read.table(file="/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/bipartite_site_counts/I_5_new_test.txt", row.names=1)

inds_use <- intersect(rownames(counts_old), rownames(counts_new))

plot(counts_old[inds_use, 1], counts_new[inds_use, 1], log="xy")
identify(counts_old[inds_use, 1], counts_new[inds_use, 1], labels=inds_use)


break

DoubleSiteModel <- function(a, kd1, kd2) {
  theta_1 <- a / (a + kd1)
  theta_2 <- a / (a + kd2)
  total_occ <- 1 - (1 - theta_1)*(1 - theta_2)
  return(total_occ)
}

kd1 <- 1
kd2 <- 10

a <- 10^seq(-2, 2, length.out=30)

occs <- DoubleSiteModel(a, kd1, kd2)



occs_alt <- (a^2 + (kd1 + kd2)*a)/(a^2 + (kd1 + kd2)*a + kd1*kd2)
KdApp <- function(kd1, kd2) {
  A <- 1
  B <- (kd1 + kd2)
  C <- -kd1*kd2

  root1 <- (-B + sqrt(B^2 - 4*A*C))/2*A
  root2 <- (-B - sqrt(B^2 - 4*A*C))/2*A
  return(max(root1, root2))  
}

kd1_range <- 10^seq(-4, 0, length.out=30)
kd2 <- 0.1 

kd_app <- sapply(kd1_range, KdApp, kd2=kd2)

dev.new(xpos=20, ypos=20, height=4, width=4)
plot(kd1_range, kd_app, log="xy", xlim=c(1e-4, 1), ylim=c(1e-4, 1))
break
print(root1)
print(root2)




plot(a, occs, log="x")
lines(a, occs_alt, col="red")
segments(x0=root1, y0=0, y1=1, col="blue")
segments(x0=root2, y0=0, y1=1, col="red")
break



graphics.off()

PlotBestThreePrimeSite()

kds_programmed <- EquilPars(mirna="let-7a-21nt", n_constant=0,
                            experiment="equil_c2_nb",
                            sitelist="programmed_suppcomp")
# kds_random <- EquilPars(mirna="let-7a-21nt", n_constant=0,
#                         experiment="equil_c_nb", sitelist="programmed_suppcomp")

kds_random <- EquilPars(mirna="let-7a", n_constant=5,
                        experiment="equilibrium", sitelist="bipartite_random_suppcomp")

kds_random_no_combI <- EquilPars(mirna="let-7a", n_constant=5,
                        experiment="equilibrium", sitelist="bipartite_random_suppcomp",
                        combined=FALSE)

plot(kds_random[, 2], kds_random_no_combI[, 2], log="xy")
break

sXc <- SitesXCounts("let-7a", "equilibrium", n_constant=5,
                    sitelist="bipartite_random_suppcomp")

print(head(sXc))
break


# break
kds_random <- EquilPars(mirna="let-7a", experiment="equilibrium",
                        sitelist="bipartite_random_suppcomp")

# kds_random <- EquilPars(mirna="let-7a-21nt", experiment="equilibrium_nb",
#                         sitelist="bipartite_random_suppcomp")


# names_drop <- sprintf("%s_Kd", rownames(sXc_rand)[which(rowSums(sXc_rand[, 3:8]) == 0)])
# print(names_drop)
# break

# kds_random <- kds_random[setdiff(rownames(kds_random), names_drop), ]

# names_prog <- rownames(kds_random)


# names_ind <- grep("\\|", names_prog, perl=TRUE)
# names_num <- grep("\\|", names_prog, perl=TRUE, value=TRUE)
# print(head(names_num))
# site_l <- gsub("^(.*)\\|(.*)\\|(.*)$", replacement="\\1", names_num, perl=TRUE)
# pos_m <- as.integer(gsub("^(.*)\\|(.*)\\|(.*)$", replacement="\\2", names_num, perl=TRUE)) + 9
# site_r <- gsub("^(.*)\\|(.*)\\|(.*)$", replacement="\\3", names_num, perl=TRUE)

# names_new <- sprintf("%s|%s|%s", site_l, pos_m, site_r)
# print(head(names_new))

# rownames(kds_random)[names_ind] <- names_new






print(head(kds_global))
print(head(kds_programmed))
print(head(kds_random))
# break
names_use <- intersect(rownames(kds_programmed), rownames(kds_random))
names_use <- grep("Comp", names_use, value=TRUE)
# names_use <- grep("Supp", names_use, value=TRUE)

# names_use <- grep("^9mer", names_use, value=TRUE, perl=TRUE)
print(names_use)
# names_use <- grep("Supp", names_use, value=TRUE, invert=TRUE)
tick <- 0
for (kmer_len in 4:11) {
  names_grep <- grep(sprintf("^%smer", kmer_len), names_use, perl=TRUE, value=TRUE)
  xjump <- ((kmer_len - 4)%%4)*300
  yjump <- floor((kmer_len - 4)/4)*300
  dev.new(xpos=20 + xjump, ypos=20 + yjump, height=3, width=3)

  # plot(x_coords, kds_programmed[names_grep, 2], log="y", type="o", ylim=c(1e-3, 1))
  # lines(x_coords, kds_random[names_grep, 2], col="red", pch=20)
  # break


  plot(kds_programmed[names_grep, 2], kds_random[names_grep, 2], log="xy",
       xlim=c(1e-4, 1), ylim=c(1e-4, 1))
  # identify(kds_programmed[names_use, 2], kds_random[names_use, 2], labels=names_use)

  r_squared <- cor(log(kds_programmed[names_grep, 2]), log(kds_random[names_grep, 2]))^2
  text(1e-3, 0.5, label=sprintf("r^2: %3.2f", r_squared), adj=0)  
}

break


# kon is 3.9e8 M-1 s-1

kons <- c(9, 8, 7, 6, 5, 4, 3, 2, 1, 0.5, 0.01, 0.005, 0.001)

cols <- c("red", rep(c("purple", "forestgreen", "blue", "lightblue"), 3))

dev.new(xpos=20, ypos=20, height=5, width=3)

plot(c(0), c(0), type="l", col="white", xlim=c(0, 150), ylim=c(0, 0.5))

k_obs <- c(0, 0, 0, 0)

koff = 0.0035

frac_equil <- rep(0, length(kons))

for (i in 1:length(kons)) {

  parameters <- c(kon = kons[i]*1e9/(1e12)*60, koff=koff)

  state <- c(l=0, a=0, x=0.5)

  Dil <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dx <- l*a*kon - koff*x
      da <- koff*x - l*a*kon
      dl <- koff*x - l*a*kon
      list(c(dl, da, dx))
    })
  }

  times <- seq(0, 2000, by=0.01)
  Exp <- function(times, b1, b0, k) {
    b1*(1 - exp(-1*k*times)) + b0
  }

  Cost <- function(pars, y) {
    b1 <- pars[1]
    b0 <- pars[2]
    k  <- pars[3]
    y_exp <- Exp(times, b1, b0, k)
    sum((y - y_exp)^2)
  }
  out <- ode(y=state, times=times, func=Dil, parms=parameters, maxsteps=10000)

  pars_init <- c(0.25, 0.25, 0.035)
  opt <- optim(pars_init, Cost, y=out[, 4])
  print(opt)
  opt_pars <- opt$par
  b1 <- opt_pars[1]
  b0 <- opt_pars[2]
  k <- opt_pars[3]
  k_obs[i] <- k
  head(out)
  lines(times, out[, 4], col=cols[i])
  frac_equil[i] <- out[nrow(out), 4]/out[1, 4]
  print(length(times))
  y_plot <- Exp(times, b1, b0, k)
  print(head(y_plot))
  lines(times, Exp(times, b1, b0, k), col=cols[i], lty=2)

}

dev.new(height=5, width=5, xpos=320, ypos=20)


plot(k_obs, koff/(1 - frac_equil), col=cols)

koff_calc <- koff/(1 - frac_equil)

lm_fit <- lm(koff_calc ~ k_obs)

coefs <- lm_fit$coefficients
b <- coefs[1]
m <- coefs[2]

segments(x0=0.1, y0=0.1, x1=koff, y1=koff, lty=2)

x_fit <- seq(0, 0.1, length.out=20)
y_fit <- x_fit*m + b

lines(x_fit, y_fit, col="gray")

break



PlotSiteMismatches_test <- function(
  mirna, experiment, n_constant, start_mm, stop_mm, bulge=FALSE, xpos=20, ypos=20,
  height=4, width=4
) {
  # test_kds <- EquilPars("let-7a-21nt", "equil_c2_nb", 0, "programmed_suppcomp")
  kds <- SubfunctionCall(EquilPars, sitelist="programmed_suppcomp")
  print(dim(kds))
  # Build the string to search for the correct sites.
  site_len <- stop_mm - start_mm + 1
  print(site_len)
  if (bulge) {
    str_mmb <- "b"
  } else {
    str_mmb <- "mm"
  }
  # Make vectors of the wt sites, and the mismatch/bulge sites with names.
  str_grep_wt <- sprintf("^%smer-m%s.%s\\|.*\\|Comp_Kd", site_len, start_mm, stop_mm, str_mmb)
  names_wt <- grep(str_grep_wt, rownames(kds), perl=TRUE, value=TRUE)
  kds_wt <- kds[grep(str_grep_wt, rownames(kds)), 2]
  names(kds_wt) <- names_wt
  str_grep_mmb <- sprintf("^%smer-m%s.%s%s", site_len, start_mm, stop_mm, str_mmb)
  names_mmb <- grep(str_grep_mmb, rownames(kds), perl=TRUE, value=TRUE)
  names_mmb <- grep("Comp_Kd", names_mmb, value=TRUE)
  kds_mmb <- kds[names_mmb, 2]
  names(kds_mmb) <- names_mmb

  wt_pos <- gsub("^.*\\|(.*)\\|.*$", replacement="\\1", names(kds_wt),
                 perl=TRUE)
  names(kds_wt) <- wt_pos
  # Make the first matrix, that is all of mismatch types aligned by position,
  # in order to pick the offset range (3 long) that has the Kd in the wild type
  # case.
  # First change the rownames to only have the mismatch/bulge type and position
  # of the site.
  names(kds_mmb) <- gsub(sprintf("^%smer-m%s\\.%s%s(.*)\\|(.*)\\|Comp_Kd$",
                                site_len, start_mm, stop_mm, str_mmb),
                    replacement="\\1\\|\\2", names(kds_mmb))
  mmb_sites <- sapply(names(kds_mmb), function(name_i) {
    unlist(strsplit(name_i, split="\\|"))[1]
  })
  print(mmb_sites)
  # Create the column name vector for the matrix to be made.
  colnames_mat <- c("wt", unique(mmb_sites))
  # Pre-allocate matrix.
  out_mat <- matrix(0, nrow=length(kds_wt), ncol=length(colnames_mat),
                    dimnames=list(names(kds_wt), colnames_mat))
  out_mat[names(kds_wt), 1] <- kds_wt
  # Populate the matrix with all of the values.
  for (name_i in names(kds_mmb)) {
    mmb_pos <- unlist(strsplit(name_i, split="\\|"))
    out_mat[mmb_pos[2], mmb_pos[1]] <- kds_mmb[name_i]
  }
  # Make the matrix where each of three rows are averaged together.
  out_mat_average <- t(sapply(1:(nrow(out_mat) - 2), function(start_row) {
    colMeans(log(out_mat)[start_row:(start_row + 2), ], na.rm=TRUE)
  }))
  # Identify which row to use, and use that to make the vector with only
  # that positional average.
  inds_use <- which.min(out_mat_average[, 1])
  mm_vector <- out_mat_average[inds_use, ]
  names(mm_vector) <- colnames(out_mat_average)
  # make the new matrix that will be used to visualize
  mm_mat <- matrix(mm_vector["wt"],
                    nrow=4,
                    ncol=(length(mm_vector) - 1)/3,
                    dimnames=list(c("A", "C", "G", "T"),
                                  as.character((start_mm + 1):(stop_mm - 1))))
  image(mm_mat)


  out_b <- out_mat_use[grep("b", names(out_mat_use), value=TRUE)]

  print(out_mm)
  print(out_b)


  for (mm_i in names(mm_vector)[-1]) {
    print(mm_i)
    nuc <- substr(mm_i, start=1, stop=1)
    pos <- substr(mm_i, start=2, stop=nchar(mm_i))
    print(nuc)
    print(pos)
    value <- mm_vector[mm_i]
    mm_mat[nuc, pos] <- value
    image(mm_mat)
  }



  break
  kds_comp_mmb <- kds_mmb[grep("Comp", rownames(kds_mmb)), ]
  print(dim(kds_mmb))
  print(dim(kds_comp_mmb))
  inds_base <- grep("^9mer-m11\\.19\\|.*\\|Comp_Kd$", rownames(test_kds_mm),
                    perl=TRUE, value=TRUE)
  print(inds_base)
  pos_base <- gsub("^.*\\|(.*)\\|.*$", replacement="\\1", inds_base, perl=TRUE)
  print(pos_base)

  mmb_sites <- gsub("^9mer-m11\\.19(.*)\\|(.*)\\|Comp_Kd$",
                    replacement="\\1\\|\\2", rownames(kds_comp_mmb))

  print(mmb_sites)

  unique_mmb_sites <- unique(gsub("^(.*)\\|.*$", replacement="\\1", mmb_sites))

  print(unique_mmb_sites)




  out_mat[, 1] <- test_kds_mm[inds_base, 2]

  for (row_i in 1:length(mmb_sites)) {
    print(mmb_sites[row_i])
    row_split <- unlist(strsplit(mmb_sites[row_i], split="\\|"))
    print(row_split)
    out_mat[row_split[2], row_split[1]] <- kds_comp_mmb[row_i, 2]
  }


  out_mat_use <- 10^colMeans(out_mat[c("14", "15", "16"), ])

  mm_mat <- matrix(out_mat_use["wt"],
                    nrow=4,
                    ncol=length(grep("mm", names(out_mat_use)))/3,
                    dimnames=list(c("A", "C", "G", "T"),
                                  as.character(12:18)))

  out_mm <- out_mat_use[grep("mm", names(out_mat_use), value=TRUE)]

  out_b <- out_mat_use[grep("b", names(out_mat_use), value=TRUE)]

  print(out_mm)
  print(out_b)


  for (mm_i in names(out_mm)) {

    nuc <- substr(mm_i, start=3, stop=3)
    pos <- substr(mm_i, start=4, stop=nchar(mm_i))
    print(nuc)
    print(pos)
    value <- out_mm[mm_i]
    mm_mat[nuc, pos] <- value
  }

  image(log10(mm_mat))
}

PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 0, 11, 19)
break


# l7 <- FitOffsetAndPairingModel("let-7a-21nt", "equil_c2_nb")
# out <- FitPairingOffsetAndMismatchModel("let-7a-21nt", "equil_c2_nb")

# m1 <- FitOffsetAndPairingModel("miR-1", "equil_c_nb")
# m1_mm <- FitOffsetPairingAndMismatchModel("miR-1", "equil_c_nb")

# m155 <- FitOffsetAndPairingModel("miR-155", "equil_sc_nb")
# m155_mm <- FitOffsetPairingAndMismatchModel("miR-155", "equil_sc_nb")

# l7m155 <- FitOffsetAndPairingModel("let-7a_miR-155", "equil_c_nb")
# l7m155_mm <- FitOffsetPairingAndMismatchModel("let-7a_miR-155", "equil_c_nb")

# m155l7 <- FitOffsetAndPairingModel("miR-155_let-7a", "equil_c_nb")
# m155l7_mm <- FitOffsetPairingAndMismatchModel("miR-155_let-7a", "equil_c_nb")


# pairing_orig <- out2$pairing
# offsets_orig <- out2$offsets
# base_orig <- out2$base
pairing <- out$pairing
offsets <- out$offsets
mm <- out$mm
# mm_chimera <- out_chimera$mm
# # offsetsmm <- out$offsetsmm
base <- out$base

# dev.new(xpos=20, ypos=20, height=3, width=3)
# par(mar=c(2, 2, 1, 1))


# plot(-4:16, l7$offsets, col="red", xlim=c(-4, 16), ylim=c(0, 1))
# points(-4:16, l7_mm$offsets, type="o", col="blue")

df_model <- out$data

# dev.new(xpos=30, ypos=20, height=3, width=3)
# par(mar=c(2, 2, 1, 1))



# plot(1:18, l7_mm$mm, col="red",
#      ylim=c(0, 2), type="o")
# lines(1:18, l7m155_mm$mm, col="blue",
#      ylim=c(0, 2), type="o", pch=20)
# dev.new(xpos=20, ypos=320, height=3, width=3)
# par(mar=c(2, 2, 1, 1))
# plot(l7_mm$mm, l7m155_mm$mm)
# mtext(text=round(cor(l7_mm$mm, l7m155_mm$mm)^2, 2), side=3, line=0, at=1)



# dev.new(xpos=320, ypos=20, height=3, width=3)
# par(mar=c(2, 2, 1, 1))

# plot(1:18, m155_mm$mm, col="red",
#      ylim=c(0, 2), type="o")
# lines(1:18, m155l7_mm$mm, col="blue",
#      ylim=c(0, 2), type="o", pch=20)

# dev.new(xpos=320, ypos=320, height=3, width=3)
# par(mar=c(2, 2, 1, 1))
# plot(m155_mm$mm, m155l7_mm$mm)
# mtext(text=round(cor(m155_mm$mm, m155l7_mm$mm)^2, 2), side=3, line=0, at=1)

# dev.new(xpos=620, ypos=20, height=3, width=3)
# par(mar=c(2, 2, 1, 1))
# inds <- intersect(names(m155_mm$mm), names(l7_mm$mm))
# plot(m155_mm$mm[inds], l7_mm$mm[inds])
# mtext(text=round(cor(m155_mm$mm[inds], l7_mm$mm[inds])^2, 2), side=3, line=0, at=1)

# dev.new(xpos=620, ypos=320, height=3, width=3)
# par(mar=c(2, 2, 1, 1))
# inds <- intersect(names(m155_mm$mm), names(l7_mm$mm))
# plot(m155l7_mm$mm[inds], l7m155_mm$mm[inds])
# mtext(text=round(cor(m155l7_mm$mm[inds], l7m155_mm$mm[inds])^2, 2), side=3, line=0, at=1)

# dev.new(xpos=920, ypos=20, height=3, width=3)
# par(mar=c(2, 2, 1, 1))
# inds <- intersect(names(m1_mm$mm), names(l7_mm$mm))
# plot(m1_mm$mm[inds], l7_mm$mm[inds])
# mtext(text=round(cor(m1_mm$mm[inds], l7_mm$mm[inds])^2, 2), side=3, line=0, at=1)

# dev.new(xpos=920, ypos=320, height=3, width=3)
# par(mar=c(2, 2, 1, 1))
# inds <- intersect(names(m1_mm$mm), names(m155_mm$mm))
# plot(m1_mm$mm[inds], m155_mm$mm[inds])
# mtext(text=round(cor(m1_mm$mm[inds], m155_mm$mm[inds])^2, 2), side=3, line=0, at=1)




# break

row_ind <- 1
model_sim <- apply(df_model, 1, function(row) {
  # row <- as.character(row)
  # print(row)
  # # row_ind <<- row_ind + 1
  # # print(c(row[4], row[3]))
  # print(as.character(as.integer(row[3])))
  # print(as.character(as.integer(row[4])))
  # length <- as.integer(row[4]) - as.integer(row[3]) + 1
  # if (length > 5) {
    pairing[as.character(as.integer(row[4])), as.character(as.integer(row[3]))]*mm[row[7]]*offsets[row[6]] + base

    # } else {
    # pairing[as.character(as.integer(row[4])), as.character(as.integer(row[3]))]*offsets[row[6]] + base

    # }

})


# model_vals_check <- pairing[*offsets + mm*offsetsmm + base
dev.new(xpos=220, ypos=20, height=3, width=3)
par(mar=c(2, 2, 1, 1))
plot(model_vals, model_sim)
dev.new(xpos=620, ypos=620, height=3, width=3)
par(mar=c(3, 3, 1, 1))
uniq_5ps <- unique(df_model$pos_5p3p)
uniq_offsets <- unique(df_model$offset)

cols <- rainbow(18, start=0, end=0.7, alpha=0.5)
names(cols) <- names(sort(mm))

cols_offsets <- rainbow(length(uniq_offsets), start=0, end=0.7, alpha=0.5)
names(cols_offsets) <- uniq_offsets
plot(-1:1, -1:1, col="white")
title(xlab="model", line=2)
title(ylab="data", line=2)

dev.new(xpos=920, ypos=620, height=3, width=3)
par(mar=c(3, 3, 1, 1))
plot(0:3, 0:3, col="white")
title(xlab="model", line=2)
title(ylab="data", line=2)





for (offset_i in uniq_offsets) {
  for (p5_i in uniq_5ps) {
    start_stop <- as.integer(unlist(strsplit(p5_i, split="\\|")))
    length_site <- start_stop[2] - start_stop[1] + 1
    if (p5_i == "12|19" & offset_i == "14") {
      print(offset_i)
      print(p5_i)
      inds <- which(df_model$pos_5p3p == p5_i & df_model$offset == offset_i)
      if (length(inds) > 0) {
        mm_types <- as.character(df_model$mm[inds])

        x_vals <- model_sim[inds]
        y_vals <- df_model$logkd[inds]
        x_vals_norm <- x_vals - mean(x_vals)
        y_vals_norm <- y_vals - mean(y_vals)
        lm_fit <- lm(y_vals_norm ~ x_vals_norm)
        slope <- lm_fit$coefficients[2]
        print(slope)
        dev.set(3)
        points(x_vals_norm, y_vals_norm, col=cols[mm_types], pch=20)
        x_lines <- seq(-1, 1, length.out=10)
        lines(x_lines, slope*x_lines, col=rgb(0, 0, 0, alpha=0.1))
        dev.set(4)
        # lm_fit <- lm(y_vals ~ x_vals)
        # slope <- lm_fit$coefficients[2]
        # b <- lm_fit$coefficients[1]

        x_lines <- seq(0, 3, length.out=10)
        print(offset_i)
        print(cols_offsets)
        print(cols_offsets[offset_i])
        points(x_vals, y_vals, col=cols_offsets[offset_i], pch=20)
        if (slope < 0) {
          print(offset_i)
          print(p5_i)
          break
        }
        # points(c(mean(y_vals)), c(sd(y_vals)), pch=20)
        # lines(x_lines, slope*x_lines + b, col=rgb(0, 0, 0, alpha=0.1))
      }
    }
  }
}
break
inds <- which(df_model$pos_5p3p == "11|18" & df_model$offset == "4")

plot(df_model$logkd[inds], )

inds <- which(df_model$pos_5p3p == "11|17" & df_model$offset == "4")

dev.new(xpos=1220, ypos=20, height=4, width=4)
plot(df_model$logkd[inds], model_sim[inds])

dev.new(xpos=20, ypos=20, height=4, width=4)




break



p1_lag <- ccf(off_l7$offsets, off_l7_p1$offsets, plot=FALSE, i=seq(-10, 10, length.out=100), lag.max=5)
m1_lag <- ccf(off_l7$offsets, off_l7_m1$offsets, plot=FALSE, i=seq(-10, 10, length.out=100), lag.max=5)

plot(-5:5, p1_lag$acf, xlim=c(-5, 5), ylim=c(-0.5, 1), type="o")
lines(-5:5, m1_lag$acf, col="red", type="o")


# print(test_l7_m1)
# print(test_l7)
# print(test_l7_p1)

break


break

# test <- MakeFullNucleotideAndMisMatchContributionDf(
#   "miR-1", "equil_c_nb"
# )


# check <- MakeNucleotideContributionMatrix("miR-1", "equil_c_nb", n_constant=0, offset=4)
# print(check)

kds <- EquilPars("let-7a-21nt", "equil_c2_nb", n_constant=0, sitelist="programmed_suppcomp")

counts<- SitesXCounts("let-7a-21nt", "equil_c2_nb", n_constant=0, sitelist="programmed")

df1 <- MakeNucleotideContributionMatrix("let-7a-21nt", "equil_c2_nb", 1, n_constant=0)

df2 <- MakeNucleotideContributionMatrix("let-7a-21nt", "equil_c2_nb", 2, n_constant=0)

df3 <- MakeNucleotideContributionMatrix("let-7a-21nt", "equil_c2_nb", 3, n_constant=0)

df4 <- MakeNucleotideContributionMatrix("let-7a-21nt", "equil_c2_nb", 4, n_constant=0)

df5 <- MakeNucleotideContributionMatrix("let-7a-21nt", "equil_c2_nb", 5, n_constant=0)

df6 <- MakeNucleotideContributionMatrix("let-7a-21nt", "equil_c2_nb", 6, n_constant=0)


dfs <- c(df1["20", "10"], df2["20", "10"], df3["20", "10"], df4["20", "10"],
          df5["20", "10"], df6["20", "10"])
names(dfs) <- 1:6


print(dfs)
inds_kd <- grep("^11mer-m10\\.20\\|.*\\|Comp_Kd", rownames(kds), perl=TRUE)

print(kds[inds_kd, ])

break

inds_counts <- grep("^11mer-m10\\.20\\|.*\\|Comp", rownames(counts), perl=TRUE)
print(counts[inds_counts, ])
break


ypos_global <- 20

# check <- FitOffsetAndPairingModel("let-7a-21nt", "equil_c2_nb", n_constant=0)
# graphics.off()
# dev.new(xpos=20, ypos=20, height=4, width=4)
# plot(model$fitted.values, model$model$logkd)
# text(x=0.2, y=1, labels=cor(model$fitted.values, model$model$logkd)^2, adj=0)

# test <- FitOffsetAndPairingNonlinModel("let-7a-21nt", "equil_c2_nb", n_constant=0)


# check <- FitOffsetAndPairingModel("miR-1", "equil_c_nb", n_constant=0)
# dev.new(xpos=20, ypos=420, height=4, width=4)
# plot(model$fitted.values, model$model$logkd)
# text(x=0.2, y=1, labels=cor(model$fitted.values, model$model$logkd)^2, adj=0)

# test <- FitOffsetAndPairingNonlinModel("miR-1", "equil_c_nb", n_constant=0)

# break
# check <- FitOffsetAndPairingModel("miR-155", "equil_sc_nb", n_constant=0)
# dev.new(xpos=20, ypos=820, height=4, width=4)
# plot(model$fitted.values, model$model$logkd)
# text(x=0.2, y=1, labels=cor(model$fitted.values, model$model$logkd)^2, adj=0)

test <- FitOffsetAndPairingNonlinModel("miR-155", "equil_sc_nb", n_constant=0)
break


# test2 <- FitOffsetAndPairingNonlinModel("let-7a-21nt", "equil_c2_nb",
#                                         n_constant=0, mm=TRUE)




break


test1 <- GetPositionalProgKds("let-7a-21nt", "8mer-mmA7", "equil_c2_nb")

test2 <- GetPositionalProgKds("let-7a-21nt", "8mer-mmA7", "equil_c_nb")

test3 <- GetPositionalProgKds("let-7a-21nt", "8mer-mmA7", "equil_c2_nb", n_constant=0, outto11mers=TRUE, downto4mers=TRUE, end3prand=TRUE, suppcomp=TRUE)

test4 <- GetPositionalProgKds("let-7a-21nt", "8mer-mmA7", "equil_c_nb", n_constant=0, outto11mers=TRUE, downto4mers=TRUE, end3prand=TRUE, suppcomp=TRUE)


check1 <- SitesXCounts("let-7a-21nt", "equil_c2_nb", n_constant=5, sitelist="programmed")
check2 <- SitesXCounts("let-7a-21nt", "equil_c_nb", n_constant=5, sitelist="programmed")

kds1 <- GetPositionalProgKds("let-7a-21nt", "9mer-m11.19", "equil_c2_nb", n_constant=5, outto11mers=FALSE, downto4mers=FALSE, end3prand=FALSE, suppcomp=TRUE)
kds2 <- GetPositionalProgKds("let-7a-21nt", "9mer-m11.19", "equil_c_nb", n_constant=5, outto11mers=FALSE, downto4mers=FALSE, end3prand=FALSE, suppcomp=TRUE)

rownames(kds1) <- gsub("_Kd", replacement="", rownames(kds1))
rownames(kds2) <- gsub("_Kd", replacement="", rownames(kds2))

inds1 <- grep("^9mer-m11\\.19\\|", rownames(check1), perl=TRUE)
inds2 <- grep("^9mer-m11\\.19\\|", rownames(check2), perl=TRUE)

l1 <- Norm(check1[, 1])
l2 <- Norm(check2[, 1])

norm1 <- t(t(check1[, 3:8])/colSums(check1[, 3:8]))
norm2 <- t(t(check2[, 3:8])/colSums(check2[, 3:8]))

R1 <- norm1/l1
R2 <- norm2/l2


print(R1[inds1, ])
print(R2[inds2, ])

R_average_1 <- rowMeans(R1[inds1, 1:5], na.rm=TRUE)
R_average_2 <- rowMeans(R2[inds2, 1:5], na.rm=TRUE)

print(head(R_average_1))
print(head(kds1))
dev.new(xpos=20, ypos=20, height=5, width=5)
plot(R_average_1[intersect(names(R_average_1), rownames(kds1))],
     kds1[intersect(names(R_average_1), rownames(kds1)), 2], log="xy",
     xlim=c(0.01, 10), ylim=c(1e-2, 1e4))

dev.new(xpos=520, ypos=20, height=5, width=5)
plot(R_average_2[intersect(names(R_average_2), rownames(kds2))],
     kds2[intersect(names(R_average_2), rownames(kds2)), 2], log="xy",
     xlim=c(0.01, 10), ylim=c(1e-2, 1e4))


break
print(colSums(norm1))
print(colSums(norm2))
print(head(norm1))
print(head(norm2))
break

dev.new(xpos=20, ypos=20, height=5, width=5)
plot(test1[, 2], test2[, 2], log="xy")

dev.new(xpos=520, ypos=20, height=5, width=5)
plot(test3[, 2], test4[, 2], log="xy")

break


kds <- EquilPars("let-7a-21nt", "equil_c2_nb", n_constant=0, sitelist="programmed_suppcomp",
                 outto11mers=TRUE, downto4mers=TRUE, collapsemm=TRUE, end3prand=TRUE)
print(head(kds))

kds_new <- EquilPars("let-7a-21nt", "equil_c2_nb", n_constant=0, sitelist="programmed_suppcomp",
                     outto11mers=TRUE, downto4mers=TRUE, end3prand=TRUE, collapsemm=TRUE, new=TRUE)
print(head(kds_new))
break



test2 <- MakeNucleotideContributionMatrix("let-7a-21nt", "equil_c2_nb", 4,
                                         n_constant=0, outto11mers=TRUE, downto2mers=TRUE,
                                         end3prand=TRUE)

test3 <- MakeNucleotideContributionMatrix("let-7a-21nt", "equil_c2_nb", 4,
                                         n_constant=0, outto11mers=TRUE, downto3mers=TRUE,
                                         end3prand=TRUE)

test4 <- EquilPars("let-7a-21nt", "equil_c2_nb", n_constant=0,
                   sitelist="programmed_suppcomp", outto11mers=TRUE,
                   downto3mers=TRUE, end3prand=TRUE, collapsemm=FALSE)

test4_collapsemm <- EquilPars("let-7a-21nt", "equil_c2_nb", n_constant=0,
                   sitelist="programmed_suppcomp", outto11mers=TRUE,
                   downto3mers=TRUE, end3prand=TRUE, collapsemm=TRUE)



grep_string <- "^3mer-m[0-9]{1,2}\\.[0-9]{1,2}\\|[0-9]{1,2}\\|Comp_Kd$"

inds_4_3mers <- grep(grep_string, rownames(test4), perl=TRUE)
inds_4_3mers_collapsemm <- grep(grep_string, rownames(test4_collapsemm), perl=TRUE)
print(" ")
print(head(test4, n=24))
print(" ")
print(test4[inds_4_3mers, ]/GeoMean(test4[7:24, 2]))
print(" ")
print(head(test4_collapsemm))
print(" ")
print(test4_collapsemm[inds_4_3mers_collapsemm, ])

plot(test4[inds_4_3mers, 2]/GeoMean(test4[7:24, 2]), test4_collapsemm[inds_4_3mers_collapsemm, 2], log="xy")



break


print(average_kmers)
mirnas_use <- c("let-7a-21nt", "miR-1", "miR-155")
experiments_use <- c("equil_c2_nb", "equil_c_nb", "equil_sc_nb")

for (j in 1:3) {


  mirna <- mirnas_use[j]
  experiment <- experiments_use[j]
  n_constant <- 0
  sitelist <- "programmed_suppcomp"


  probs <- sapply(2:11, function(k) {
    prob_k <- ProbKmerInRandom(k, 25)
    prob_k_plus_1 <- ProbKmerInRandom(k + 1, 25)
    prob_3p_k <- 1 - (1 - prob_k)^length(GetMirna3pKmers(mirna, k))
    prob_not_3p_k_plus_1 <- (1 - prob_k_plus_1)^length(GetMirna3pKmers(mirna, k + 1))
    return(prob_3p_k * prob_not_3p_k_plus_1)
  })



  check <- SitesXCounts(mirna, experiment, n_constant, sitelist, downto2mers=TRUE, outto11mers=TRUE, end3prand=TRUE)
  check <- t(t(check)/colSums(check))

  average_kmers <- sapply(2:11, function(k) {
    grep_string <- sprintf("^%smer-m[0-9]{1,2}\\.[0-9]{1,2}\\|[0-9]{1,2}\\|Comp", k)
    print(grep_string)
    kmers <- check[grep(grep_string, rownames(check), perl=TRUE), 1]
    nonzero_kmers <- kmers[which(kmers != 0)]
    return(mean(nonzero_kmers))
    })

  average_kmers_40 <- sapply(2:11, function(k) {
    grep_string <- sprintf("^%smer-m[0-9]{1,2}\\.[0-9]{1,2}\\|[0-9]{1,2}\\|Comp", k)
    print(grep_string)
    kmers <- check[grep(grep_string, rownames(check), perl=TRUE), 3]
    nonzero_kmers <- kmers[which(kmers != 0)]
    return(mean(nonzero_kmers))
  })

  average_kmers_12.6 <- sapply(2:11, function(k) {
    grep_string <- sprintf("^%smer-m[0-9]{1,2}\\.[0-9]{1,2}\\|[0-9]{1,2}\\|Comp", k)
    print(grep_string)
    kmers <- check[grep(grep_string, rownames(check), perl=TRUE), 4]
    nonzero_kmers <- kmers[which(kmers != 0)]
    return(mean(nonzero_kmers))
  })

  average_kmers_4 <- sapply(2:11, function(k) {
    grep_string <- sprintf("^%smer-m[0-9]{1,2}\\.[0-9]{1,2}\\|[0-9]{1,2}\\|Comp", k)
    print(grep_string)
    kmers <- check[grep(grep_string, rownames(check), perl=TRUE), 5]
    nonzero_kmers <- kmers[which(kmers != 0)]
    return(mean(nonzero_kmers))
  })

  average_kmers_1.2 <- sapply(2:11, function(k) {
    grep_string <- sprintf("^%smer-m[0-9]{1,2}\\.[0-9]{1,2}\\|[0-9]{1,2}\\|Comp", k)
    print(grep_string)
    kmers <- check[grep(grep_string, rownames(check), perl=TRUE), 6]
    nonzero_kmers <- kmers[which(kmers != 0)]
    return(mean(nonzero_kmers))
  })

  average_kmers_0.4 <- sapply(2:11, function(k) {
    grep_string <- sprintf("^%smer-m[0-9]{1,2}\\.[0-9]{1,2}\\|[0-9]{1,2}\\|Comp", k)
    print(grep_string)
    kmers <- check[grep(grep_string, rownames(check), perl=TRUE), 7]
    nonzero_kmers <- kmers[which(kmers != 0)]
    return(mean(nonzero_kmers))
  })

  average_kmers_0 <- sapply(2:11, function(k) {
    grep_string <- sprintf("^%smer-m[0-9]{1,2}\\.[0-9]{1,2}\\|[0-9]{1,2}\\|Comp", k)
    print(grep_string)
    kmers <- check[grep(grep_string, rownames(check), perl=TRUE), 8]
    nonzero_kmers <- kmers[which(kmers != 0)]
    return(mean(nonzero_kmers))
  })



  dev.new(xpos=20, ypos=20 + 400*(j - 1), height=4, width=4)

  plot(probs, average_kmers, log="xy", xlim=c(1e-7, 1e0), ylim=c(1e-7, 1e0))


  test_kds <- EquilPars(mirna, experiment, n_constant, sitelist, outto11mers=TRUE, downto2mers=TRUE, end3prand=TRUE)

  average_kmer_kds <- sapply(2:11, function(k) {
    grep_string <- sprintf("^%smer-m[0-9]{1,2}\\.[0-9]{1,2}\\|[0-9]{1,2}\\|Comp_Kd", k)
    print(grep_string)
    kmers <- test_kds[grep(grep_string, rownames(test_kds), perl=TRUE), 2]
    nonzero_kmers <- kmers[which(kmers != 0)]
    return(GeoMean(nonzero_kmers))
    })

  dev.new(xpos=520, ypos=20 + 400*(j - 1), height=4, width=4)
  average_kmer_kds <- average_kmer_kds[1]/average_kmer_kds

  plot(average_kmers_40/average_kmers, average_kmer_kds, log="xy", type="o", lwd=1, col="red", xlim=c(0.1, 30), ylim=c(1, 300))
  points(average_kmers_12.6/average_kmers, average_kmer_kds, type="o", lwd=1, col="orangered")
  points(average_kmers_4/average_kmers, average_kmer_kds, type="o", lwd=1, col="forestgreen")
  points(average_kmers_1.2/average_kmers, average_kmer_kds, type="o", lwd=1, col="blue")
  points(average_kmers_0.4/average_kmers, average_kmer_kds, type="o", lwd=1, col="purple")
  points(average_kmers_0/average_kmers, average_kmer_kds, type="o", lwd=1, col="gray")





}


break
freq_2mers <- check[grep("^2mer-m[0-9]{1,2}\\.[0-9]{1,2}\\|[0-9]{1,2}\\|Comp", rownames(check), perl=TRUE), col_use]
freq_2mers <- freq_2mers[which(freq_2mers != 0)]
print(freq_2mers)

dev.new(xpos=20, ypos=20, height=4, width=4)
plot(1, type="n", log="x", xlim=c(1e2, 1e6), ylim=c(0, 1))

plot(ecdf(freq_2mers), add=TRUE)

freq_3mers <- check[grep("^3mer-m[0-9]{1,2}\\.[0-9]{1,2}\\|[0-9]{1,2}\\|Comp", rownames(check), perl=TRUE), col_use]
freq_3mers <- freq_3mers[which(freq_3mers != 0)]
print(freq_3mers)

plot(ecdf(freq_3mers), add=TRUE, col="red")

freq_4mers <- check[grep("^4mer-m[0-9]{1,2}\\.[0-9]{1,2}\\|[0-9]{1,2}\\|Comp", rownames(check), perl=TRUE), col_use]
freq_4mers <- freq_4mers[which(freq_4mers != 0)]
print(freq_4mers)

plot(ecdf(freq_4mers), add=TRUE, col="purple")

freq_5mers <- check[grep("^5mer-m[0-9]{1,2}\\.[0-9]{1,2}\\|[0-9]{1,2}\\|Comp", rownames(check), perl=TRUE), col_use]
freq_5mers <- freq_5mers[which(freq_5mers != 0)]
print(freq_5mers)

plot(ecdf(freq_5mers), add=TRUE, col="blue")

freq_6mers <- check[grep("^6mer-m[0-9]{1,2}\\.[0-9]{1,2}\\|[0-9]{1,2}\\|Comp", rownames(check), perl=TRUE), col_use]
freq_6mers <- freq_6mers[which(freq_6mers != 0)]
print(freq_6mers)

plot(ecdf(freq_6mers), add=TRUE, col="forestgreen")




break


check1 <- MakeNucleotideContributionMatrix("miR-155", "equil_sc_nb", 3, outto11mers=TRUE, downto2mers=FALSE, end3prand=TRUE, collapsemm=FALSE)
print(check1)
check2_pre <- EquilPars("miR-155", "equil_sc_nb", 5, "programmed_suppcomp", outto11mers=TRUE, downto2mers=TRUE, end3prand=TRUE, collapsemm=FALSE)
print(check2_pre["2mer-m9.10|12|Comp_Kd", ])
check2 <- MakeNucleotideContributionMatrix("miR-155", "equil_sc_nb", 3, outto11mers=TRUE, downto2mers=TRUE, end3prand=TRUE, collapsemm=FALSE)
print(check2)


break


sXc <- SitesXCounts("let-7a-21nt", "equil_c2_nb", n_constant=5,
                    sitelist="programmed", outto11mers=TRUE)
sXc_collapsed <- SitesXCounts("let-7a-21nt", "equil_c2_nb", n_constant=5,
                    sitelist="programmed_collapsed", outto11mers=TRUE)
sXc_suppcomp <- SitesXCounts("let-7a-21nt", "equil_c2_nb", n_constant=5,
                    sitelist="programmed_suppcomp", outto11mers=TRUE)
sXc_suppcomp_collapsemm <- EquilPars("let-7a-21nt", "equil_c2_nb", n_constant=5,
                    sitelist="programmed_suppcomp", outto11mers=TRUE, collapsemm=TRUE)


print(head(sXc, n=50))
print(head(sXc_collapsed, n=50))
print(head(sXc_suppcomp, n=50))
print(head(sXc_suppcomp_collapsemm, n=50))

break




kds_collapsed <- EquilPars("miR-1", "equil_c_nb", n_constant=5,
                       sitelist="programmed_collapsed")


kds_sorted <- EquilPars("miR-1", "equil_c_nb", n_constant=5,
                        sitelist="programmed_collapsed", sorted=TRUE)



print(head(kds))
print(kds_sorted[1:6, 1:5])

kds_means <- do.call("cbind", list(rowMeans(kds_sorted[,  1: 4]),
                                 rowMeans(kds_sorted[,  5: 8]),
                                 rowMeans(kds_sorted[,  9:12]),
                                 rowMeans(kds_sorted[, 13:16]),
                                 rowMeans(kds_sorted[, 17:20])))

print(head(kds_means))



cols <- kSiteColors[gsub("_Kd", replacement="", rownames(kds_sorted)[1:10])]


graphics.off()

dev.new(xpos=20, ypos=20, height=5, width=5)

sds_all <- apply(kds_sorted, 1, sd)
sds_means <- apply(kds_means, 1, sd)

plot(sds_all, sds_means)

fit <- lm(sds_means ~ sds_all)
print(fit)
m <- fit$coefficients[2]
b <- fit$coefficients[1]
abline(a=b, b=m, lty=2)




alt_error <- t(apply(kds_means, 1, function(row) {
  m <- mean(row)
  n <- length(row)
  se <- sd(row)/sqrt(n)
  error <- qt(0.975, df=n - 1)*se
  return(c(m, m - error, median(row), m + error))  
}))
alt_error <- cbind(kds[, 1], 10^alt_error)



# plot(ecdf(kds_sorted[1, ]), xlim=c(-2, 0), col=cols[1])

# segments(x0=log10(kds[1, 3]), y0=0, y1=1, col=cols[1])
# segments(x0=log10(kds[1, 5]), y0=0, y1=1, col=cols[1])

# sapply(2:10, function(row) {
#   plot(ecdf(kds_sorted[row, ]), add=TRUE)
#   plot(ecdf(kds_sorted[row, ]), add=TRUE, col=cols[row])
#   segments(x0=log10(kds[row, 3]), y0=0, y1=1, col=cols[row])
#   segments(x0=log10(kds[row, 5]), y0=0, y1=1, col=cols[row])

#   segments(x0=log10(alt_error[row, 3]), y0=0, y1=1, col=cols[row], lty=2)
#   segments(x0=log10(alt_error[row, 5]), y0=0, y1=1, col=cols[row], lty=2)


# })

se <- log10(kds[, 3]/kds[, 2])
alt_se <- log10(alt_error[, 3]/alt_error[, 2])

fit_se <- lm(alt_se ~ se)
print(fit_se)
break
dev.new(xpos=520, ypos=20, height=5, width=5)

plot(se, alt_se)

l

colnames(alt_error) <- colnames(kds)
print(head(alt_error))
segments()

break





break
# kds_l7_sc <- EquilPars("let-7a-21nt", "equil_sc_nb", n_constant=5,
#                        sitelist="programmed_collapsed")

# kds_l7plus1_c <- EquilPars("let-7a_plus1", "equil_c_nb", n_constant=5,
#                          sitelist="programmed_collapsed")
# kds_l7minus1_c <- EquilPars("let-7a_minus1", "equil_c_nb", n_constant=5,
#                          sitelist="programmed_collapsed")
# kds_l7m155_c <- EquilPars("let-7a_minus1", "equil_c_nb", n_constant=5,
#                          sitelist="programmed_collapsed")

kds_rand <- EquilPars("let-7a", "equilibrium", 5, "resubmissionfinal")

kds_sorted_rand <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/kds_PAPER/5_resubmissionfinal_PAPER_sorted.txt")

kds_sorted <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a-21nt/equil_c2_nb/kds_PAPER/5_programmed_collapsed_PAPER_sorted.txt")




# for (row_use in 1:6) {
  row_use <- 1
  dev.new(xpos=20+400*((row_use - 1)%%3), ypos=20 + 400*floor((row_use - 1)/3),
          height=4, width=4)
  kds_row <- as.numeric(kds_sorted_rand[row_use, ])

  kds_means <- c(mean(kds_row[1:40]), mean(kds_row[41:80]), mean(kds_row[81:120]),
                 mean(kds_row[121:160]), mean(kds_row[161:200]))
  m <- mean(kds_means)
  n <- length(kds_means)
  se <- sd(kds_means)/sqrt(n)
  error <- qt(0.975, df=n - 1)*se
  # print(head(10^kds_sorted))
  plot(ecdf(kds_sorted_rand[row_use, ]))
  segments(x0=log10(kds_rand[row_use, 2]), y0=0, y1=1)
  segments(x0=log10(kds_rand[row_use, 3]), y0=0, y1=1, lty=2)
  segments(x0=log10(kds_rand[row_use, 5]), y0=0, y1=1, lty=2)

  segments(x0=kds_means[1], y0=0, y1=1, lty=3)
  segments(x0=kds_means[2], y0=0, y1=1, lty=3)
  segments(x0=kds_means[3], y0=0, y1=1, lty=3)
  segments(x0=kds_means[4], y0=0, y1=1, lty=3)
  segments(x0=kds_means[5], y0=0, y1=1, lty=3)

  segments(x0=m + error, y0=0, y1=1, lty=2, col="red")
  segments(x0=m - error, y0=0, y1=1, lty=2, col="red")
# }





break

site_8mm <- GetAll8merMmSites("let-7a-21nt")
kd_mm <- kds_l7_c2[sprintf("%s_Kd", site_8mm), ]

kd_mm_alt <- t(sapply(site_8mm, function(site) {
  # site <- gsub("\\.", replacement="\\\\.", site)
  print(site)
  inds <- grep(sprintf("^%s\\|8mer-mm", site), rownames(kds_l7_c2), perl=TRUE)
  # if (length(inds) > 2){
      apply(kds_l7_c2[inds, ], 2, GeoMean)
  # } else {
  #   c(NA, NA, NA, NA, NA)
  # }
}))

plot(kd_mm[, 2], kd_mm_alt[, 2], log='xy')
break


break

# print(dim(kds_l7_c))
# print(dim(kds_l7_c2))
# print(dim(kds_l7_s))
# print(dim(kds_l7_sc))
# print(dim(kds_l7plus1_c))
# print(dim(kds_l7minus1_c))
# print(dim(kds_l7m155_c))

# print(kds_l7_c["5mer-m11.15|8mer-mmA2_Kd", ])
# print(kds_l7_c2["5mer-m11.15|8mer-mmA2_Kd", ])
# print(kds_l7_s["5mer-m11.15|8mer-mmA2_Kd", ])
# print(kds_l7_sc["5mer-m11.15|8mer-mmA2_Kd", ])
# print(kds_l7plus1_c["5mer-m11.15|8mer-mmA2_Kd", ])
# print(kds_l7minus1_c["5mer-m11.15|8mer-mmA2_Kd", ])
# print(kds_l7m155_c["5mer-m11.15|8mer-mmA2_Kd", ])


# inds <- intersect(rownames(kds_l7_c2), rownames(kds_l7_c))

# inds <- grep("\\|8mer-mm[ACTG][2-7]_Kd$", inds, perl=TRUE, value=TRUE)
# inds <- grep("&", inds, value=TRUE, invert=TRUE)

# print(length(inds))

# inds_omit <- intersect(rownames(kds_l7_c2), rownames(kds_l7_c_omit))

# dev.new(xpos=20, ypos=20, height=5, width=5)
# plot(kds_l7_c2[inds, 2], kds_l7_c[inds, 2], log='xy')

# dev.new(xpos=520, ypos=20, height=5, width=5)
# plot(kds_l7_c2[inds, 2], kds_l7_sc[inds, 2], log='xy')

# dev.new(xpos=1020, ypos=20, height=5, width=5)
# plot(kds_l7_c2[inds, 2], kds_l7_s[inds, 2], log='xy')


s_5 <- GetPositionalProgKds("let-7a-21nt", "5mer-m11.15", "equil_c2_nb")
s_6 <- GetPositionalProgKds("let-7a-21nt", "6mer-m11.16", "equil_c2_nb")
s_7 <- GetPositionalProgKds("let-7a-21nt", "7mer-m11.17", "equil_c2_nb")
s_8 <- GetPositionalProgKds("let-7a-21nt", "8mer-m11.18", "equil_c2_nb")
s_9 <- GetPositionalProgKds("let-7a-21nt", "9mer-m11.19", "equil_c2_nb")

list_sites <- list(s_5, s_6, s_7, s_8, s_9)
dev.new(xpos=20, ypos=20, height=5, width=5)
cols <- c("red", "green", "forestgreen", "blue", "purple")
offset <- 0
col_ind <- 1

plot(c(0), c(0), xlim=c(31, 8), ylim=c(1e-3, 1), log='y')
bar_offset <- c(-0.4, -0.2, 0, 0.2, 0.4)

lapply(list_sites, function(df) {
  inds <- grep("8mer-mmC7", rownames(df))
  df_use <- df[inds, ]
  print(df_use)
  lines((30 - length(inds) + 1):(30), df_use[, 2], col=cols[col_ind])
  segments(x0=(30 - length(inds) + 1):(30) + bar_offset[col_ind],
           y0=df_use[, 3], y1=df_use[, 5], col=cols[col_ind], lwd=1)
  # lines((30 - length(inds) + 1):(30), df_use[, 5], col=cols[col_ind], lwd=1, lty=2)
  offset <<- offset - 1
  col_ind <<- col_ind + 1
})

break
plot(kds_l7_c[, 2], kds_l7_c2[, 2], log='xy')
break



# miR-1
kds_m1_c <- EquilPars("miR-1", "equil_c_nb", n_constant=5,
                      sitelist="programmed_collapsed")
# ISSUE: Needs `I` sample currently. DEALT WITH
kds_m1_sc <- EquilPars("miR-1", "equil_sc_nb", n_constant=5,
                       sitelist="programmed_collapsed")

# miR-155
# ISSUE: unknown
kds_m155_c  <- EquilPars("miR-155", "equil_sc_nb", n_constant=5,
                         sitelist="programmed_collapsed")
kds_m155l7_c <- EquilPars("miR-155_let-7a", "equil_c_nb", n_constant=5,
                          sitelist="programmed_collapsed")

PlotRandomRegionSite <- function(mirna, experiment, n_constant=5, site) {
  kds_df <- SubfunctionCall(EquilPars, sitelist="programmed_collapsed")
  inds_site <- grep("^%s|8mer-mm[ACGT][2-7]_Kd", rownames(kds_df), perl=TRUE)
}



inds <- intersect(rownames(kds_l7_c2), rownames(kds_l7_c))

inds <- grep("\\|8mer-mm[ACTG][2-7]_Kd$", inds, perl=TRUE, value=TRUE)
inds <- grep("&", inds, value=TRUE, invert=TRUE)


inds_omit <- intersect(rownames(kds_l7_c2), rownames(kds_l7_c_omit))

dev.new(xpos=20, ypos=20, height=5, width=5)
plot(kds_l7_c2[inds, 2], kds_l7_c[inds, 2], log='xy')

dev.new(xpos=520, ypos=20, height=5, width=5)
plot(kds_l7_c2[inds, 2], kds_l7_c_omit[inds, 2], log='xy')


break






Determine8merReferenceSiteType <- function(mirna, kmer_vec) {
  # Make the vector of letters in the 8mer site type.
  site_8mer <- StrRev(GetSiteSeq(mirna, "8mer"))
  site_8mer_vec <- unlist(strsplit(site_8mer, split=""))
  # Apply a loop over the kmers.
  site_vec <- sapply(kmer_vec, function(kmer) {
    # Condition in case the site is "None"
    if (kmer == "None") {
      return(kmer)
    # Condition in which the kmer is actually a kmer.
    } else {
      suffix <- ""
      kmer_vec <- unlist(strsplit(StrRev(kmer), split=""))
      for (i in 1:8) {
        if (site_8mer_vec[i] != kmer_vec[i]) {
          suffix <- sprintf("%smm%s%s", suffix, kmer_vec[i], i)
        } 
      }
      if (suffix != "") {
        suffix <- sprintf("-%s", suffix)
      }
      out <- sprintf("8mer%s", suffix)
      return(out)
    }
  })
  return(site_vec)
}

GetSiteMutationMatrix <- function(mirna) {
  # First make the final matrix, pre-allocated with all zeros.
  counts <- GetProgrammedKmerFrequencies(mirna, "equilibrium_mmseed_nb",
                                         condition="I")
  counts <- counts[-nrow(counts), , drop=FALSE]

  sites <- Determine8merReferenceSiteType(mirna, rownames(counts)[])
  full_matrix <- matrix(rep("None", (length(sites) + 1)^2),
                        nrow=length(sites) + 1,
                        ncol=length(sites) + 1)
  temp_col <- rep("None", length(sites) + 1)
  # Make the list of miRNA letters.
  site_8mer_vec <- unlist(strsplit(StrRev(GetSiteSeq(mirna, "8mer")), split=""))
  # Make the first column, which describes mutation from the 8mer site
  # to one of the 18 1-nt mismatches.
  sites_single_mismatch <- sites[2:25]
  col_8mer <- temp_col
  s8mer1mm <- gsub("8mer-mm", replacement="", x=sites_single_mismatch)
  pos <- as.integer(substring(s8mer1mm, first=2))
  nuc_start <- site_8mer_vec[pos]
  nuc_end <- substring(s8mer1mm, first=1, last=1)
  nuc_combined <- sprintf("%s%s", nuc_start, nuc_end)
  col_8mer[2:25] <- nuc_combined
  full_matrix[, 1] <- col_8mer
  # Make the next 24 column, which describes mutation from one of the 8mer-1mm
  # sites into either the 8mer (the first entry) or 
  # to one of the 18 1-nt mismatches.
  cols_1mm <- sapply(sites_single_mismatch, function(site) {
    # Get the 1mm suffix, allocate the mismatch nucleotide and pos.
    suffix <- unlist(strsplit(site, split="-"))[2]
    mm1_nuc <- substring(suffix, first=3, last=3)
    mm1_pos <- as.integer(substring(suffix, first=4))
    # Get the indeces of which column entries will used. Note: The first is
    # incorrect and will be changed to `1`, because each 1mm-site can mutate to
    # the 8mer.
    inds <- grep(suffix, sites)
    inds[1] <- 1
    # Series of three string operations intended to determine the identify of
    # the mutated nucleotide.
    mm2_strings <- gsub("8mer-", replacement="",
                        x=grep(suffix, sites, value=TRUE))[-1]
    # Remove the starting mismatch from these strings to yield the second
    # mismatch
    mm2 <- gsub(suffix, replacement="", x=mm2_strings)
    mm2_pos <- as.integer(substring(mm2, first=4, last=4))
    mm2_nuc <- substring(mm2, first=3, last=3)

    nuc_start <- c(mm1_nuc, site_8mer_vec[mm2_pos])
    nuc_end <- c(site_8mer_vec[mm1_pos], mm2_nuc)
    nuc_combined <- sprintf("%s%s", nuc_start, nuc_end)
    temp_col[inds] <- nuc_combined
    # Get the mutations to 1mm sites at the same position (2 possible).
    gsub_string <- sprintf("^8mer-mm.%s$", mm1_pos)
    # The starting nucleotide is the mismatch nucleotide.
    nuc_start <- rep(mm1_nuc, 2)
    # Get the index of the other two mutants at this position.
    inds <- setdiff(grep(gsub_string, sites, perl=TRUE),
                    c(grep(sprintf("^%s$", site), sites, perl=TRUE)))
    nuc_end <- substring(sites[inds], first=8, last=8)
    nuc_combined <- sprintf("%s%s", nuc_start, nuc_end)
    temp_col[inds] <- nuc_combined
    temp_col
  })
  # Allocate these columns
  full_matrix[, 2:25] <- cols_1mm
  # Now determine the indeces for the remaining double sites mutating to single
  # sites.
  double_matrix <- sapply(26:length(sites), function(i_col) {
    # First deal with reversion
    mm2 <- gsub("8mer-", replacement="", x=sites[i_col])
    mm2_vec <- unlist(strsplit(mm2, split="mm"))[-1]
    mm2_pos <- as.integer(rev(substring(mm2_vec, first=2)))
    nuc_start <- rev(substring(mm2_vec, first=1, last=1))
    nuc_end <- site_8mer_vec[mm2_pos]
    nuc_combined <- sprintf("%s%s", nuc_start, nuc_end)
    inds <- sapply(mm2_vec, function(site) {
      which(sites == sprintf("8mer-mm%s", site))
    })
    temp_col[inds] <- nuc_combined
    # Now deal with mutation to one of the other two for each.
    gsub_1 <- sprintf("8mer-mm.%smm%s", mm2_pos[2], mm2_vec[2])
    gsub_2 <- sprintf("8mer-mm%smm.%s", mm2_vec[1], mm2_pos[1])
    nuc_start <- rep(c(substring(mm2, first=3, last=3),
                       substring(mm2, first=7, last=7)), each=2)
    inds <- setdiff(c(grep(gsub_1, sites),
                      grep(gsub_2, sites)), c(i_col))
    other_mm_sites <- gsub("8mer-", replacement="", x=sites[inds])
    nuc_end <- c(substring(other_mm_sites[1:2], first=3, last=3),
                 substring(other_mm_sites[3:4], first=7, last=7))
    nuc_combined <- sprintf("%s%s", nuc_start, nuc_end)
    temp_col[inds] <- nuc_combined
    temp_col
  })
  full_matrix[, 26:length(sites)] <- double_matrix
  print(length(sites))
  print(dim(full_matrix))
  rownames(full_matrix) <- c(sites, "triple_mm")
  colnames(full_matrix) <- c(sites, "triple_mm")
  return(full_matrix)
}

CostFunction <- function(pars, data, mirna, M_strings, col, unmut_inds, pc=1e-8) {
  # Modify the parameter vector so that they are the exponential transform and
  # also so that the "None" parts of the matrix are assigned a value of `0`.
  # print(unmut_inds)
  # print(dim(data))
  # print(data)
  len_unmut <- max(unmut_inds[, 3])
  ncol_counts <- max(unmut_inds[, 2])
  # print(ncol_counts)
  # print(len_unmut)

  pars_unmut <- c(1, exp(pars[1:len_unmut]))
  pars_dinucs  <- exp(pars[(len_unmut + 1): length(pars)])
  pars_extended <- c(pars_dinucs, "None"=0)
  # print(pars_unmut)
  # print(pars_dinucs)
  # print(pars_extended)
  M <- matrix(pars_extended[M_strings], nrow=nrow(M_strings))
  # Correction for mutation away from double mutants.
  M[nrow(M), 26:(nrow(M) - 1)] <- 18*mean(pars_dinucs)
  rownames(M) <- rownames(M_strings)
  colnames(M) <- colnames(M_strings)
  # Define the starting counts of the unmutated, designed library.
  model_unmut <- matrix(0, nrow=nrow(data), ncol=ncol_counts)
  # print(dim(model_unmut))
  # print(dim(data))
  rownames(model_unmut) <- rownames(data)
  # Remove the last value from the data vector
  # Assign the indeces to be used for the model fitting.
  for (row in 1:nrow(unmut_inds)) {
    model_unmut[unmut_inds[row, 1], unmut_inds[row, 2]] <- pars_unmut[unmut_inds[row, 3] + 1]    
  }

  # print(model_unmut)
  # break
  model <- (diag(1 - colSums(M)) + M)%*%(diag(1 - colSums(M)) + M)%*%model_unmut
  data <- data + pc
  # print(colSums(data))
  # print(colSums(model))
  # print(colSums(model)/colSums(data))
  model <- model + pc*colSums(model)/colSums(data)
  # model <- t()
  # model <- model/sum(model)
  model_use <- model[-nrow(M), ]
  model_norm <- t(t(model_use)/colSums(model_use))

  data_use <- data[-nrow(M), ]
  # print(model[80:88, ])
  # print(data[80:88, ])
  # print(data_use[80:88, ])
  # print(c(unlist(data_use))[80:88])
  # print(c(unlist(data_use))[(80+277):(88+277)])
  # print(c(model_norm)[80:88])
  # print(c(model_norm)[(80+277):(88+277)])
  # break


  cost <- sum(-1*data_use[, 1]*log(model_norm[, 1]))
  cost <- sum(-1*data_use*log(model_norm))

  # cost <- sum( (log(counts[-nrow(M), 1]) - log(model[-nrow(M)]))^2 )
  tick <<- tick + 1
  if (tick%%50 == 0) {
    # print(t(cbind(model[-nrow(M)], counts[-nrow(M), 1])))
    print(length(c(model_norm[, 1])))
    print(length(c(unlist(data_use[, 1]))))
    dev.set(2)
    weird_inds <- c(109, 112, 118, 121, 154,
                    157, 160, 163, 166, 169,
                    190, 193, 196, 199, 202,
                    205, 226, 229, 230, 231,
                    232, 239, 240, 241, 257,
                    258, 259)
    col[weird_inds] <- "pink"
    plot(c(model_norm[, 1]), c(unlist(data_use[, 1])), log='xy', col=col,
       xlim=c(1e-8, 1), ylim=c(1e-8, 1), pch=rep(c(1, 19), each=277))
    segments(1e-8, 1e-8, x1=1, y1=1, lty=2)
    # if (tick%%2000 == 0) {
    #   model_norm <<- model_norm
    #   data_use <<- data_use
    #   break
    # identify(c(model_norm[, 1]), c(unlist(data_use[, 1])),
    #          labels=gsub("8mer-mm", replacement="", rownames(model_norm)))
    # }
    dev.set(3)
    plot(c(model_norm[, 2]), c(unlist(data_use[, 2])), log='xy', col=col,
       xlim=c(1e-8, 1), ylim=c(1e-8, 1), pch=rep(c(19, 19), each=277))
    segments(1e-8, 1e-8, x1=1, y1=1, lty=2)
    # if (tick%%2000 == 0) {
    # identify(c(model_norm[, 2]), c(unlist(data_use[, 2])),
    #          labels=gsub("8mer-mm", replacement="", rownames(model_norm)))
    # }


    print(exp(pars))
    print(pars_dinucs)
    print(cost)
  }

  cost
}


FitData <- function(mirna) {
  # Assign the indeces for library counts.
  ind_8       <- 1
  ind_7A1     <- 2:4
  ind_7A1_xl7 <- c(2, 4)
  ind_7m8     <- 23:25
  ind_7m8_xl7 <- c(23, 24)
  ind_6       <- 80:88
  ind_6_xl7   <- c(80, 81, 86, 87)
  ind_comp    <- 5:22
  ind_supp    <- c(ind_8, ind_7m8, ind_7A1, ind_6)

  # Get count data
  if (mirna == "let-7a") {
    mir_exp_matrix <- matrix(c(   "let-7a-21nt",   "equilibrium_mmseed_nb",
                                  "let-7a-21nt",     "equilibrium_seed_nb",
                                       "let-7a",   "equilibrium_mmseed_nb",
                                  "let-7a-21nt", "equilibrium_mmseed_2_nb",
                                 "let-7a_plus1",   "equilibrium_mmseed_nb",
                                "let-7a_minus1",   "equilibrium_mmseed_nb",
                               "let-7a_miR-155",   "equilibrium_mmseed_nb"),
                             ncol=2, byrow=TRUE)
  } else if (mirna == "miR-1") {
    mir_exp_matrix <- matrix(c("miR-1", "equilibrium_mmseed_nb",
                               "miR-1", "equilibrium_mmseed_2_nb"),
                             ncol=2, byrow=TRUE)
  } else if (mirna == "miR-155") {
    mir_exp_matrix <- matrix(c(       "miR-155", "equilibrium_mmseed_nb",
                               "miR-155_let-7a", "equilibrium_mmseed_nb"),
                             ncol=2, byrow=TRUE)
  }
  print(mir_exp_matrix)
  counts <- do.call("cbind", apply(mir_exp_matrix, 1, function(row) {
    GetProgrammedKmerFrequencies(mirna=row[1], experiment=row[2], condition="I")
  }))
  site_names <- Determine8merReferenceSiteType(mirna, rownames(counts))
  # Get the mutation matrix
  M <- GetSiteMutationMatrix(mirna)

  # Make the color vector for plotting
  cols <- rep(
    kSiteColors[c("8mer", "7mer-A1", "8mer-mmT6", "7mer-m8",
                  "11mer-m9.11", "6mer", "11mer-m9.19")],
    times=c(1, 3, 18, 3, 54, 9, 189)
  )


  if (mirna == "miR-155") {
    ind_row <- c(ind_comp, ind_6,
                 ind_comp)
    ind_col <- rep(c(1, 2), times=c(length(c(ind_comp, ind_6)),
                                    length(ind_comp)))
    unmut_inds <- cbind(ind_row, ind_col)
    len_pars <- length(ind_comp) - 1 + length(ind_6)
    pars_unmut_init <- rep(-3, len_pars)
    par_inds <- rep(0, nrow(unmut_inds))
    par_inds[2:(len_pars + 1)] <- 1:len_pars
    par_inds[(len_pars + 3):nrow(unmut_inds)] <- 1:(length(ind_comp) - 1)
    unmut_inds <- cbind(unmut_inds, par_inds)
  } else if (mirna == "miR-1") {
    ind_row <- c(ind_comp,
                 ind_comp, ind_supp)
    ind_col <- rep(c(1, 2), times=c(length(ind_comp),
                                    length(c(ind_comp, ind_supp))))
    unmut_inds <- cbind(ind_row, ind_col)
    len_pars <- length(ind_comp) - 1 + length(ind_supp)
    pars_unmut_init <- rep(-3, len_pars)
    par_inds <- rep(0, nrow(unmut_inds))
    # par_inds[2:(len_pars + 1)] <- 1:len_pars
    par_inds[2:length(ind_comp)] <- 1:(length(ind_comp) - 1)
    # par_inds[(len_pars + 3):nrow(unmut_inds)] <- 1:(length(ind_comp) - 1)
    par_inds[(length(ind_comp) + 2):nrow(unmut_inds)] <- 1:len_pars
    unmut_inds <- cbind(unmut_inds, par_inds)

  }


  # Initialize the parameters for the dinuc mutation frequencies.
  pars_dinuc_init <- rep(-5, 16)
  names(pars_dinuc_init) <- kDinucs
  pars_dinuc_init <- pars_dinuc_init[c(-1, -4, -13, -16)]
  # Combine the two vectors to make the full parameter vector.
  pars_init <- c(pars_unmut_init, pars_dinuc_init)


  opt <- optim(par=pars_init, fn=CostFunction, gr=NULL,
               data=counts, mirna=mirna, M_strings=M, col=cols,
               unmut_inds=unmut_inds,
               method="BFGS", control=list(maxit=10000))
  return(opt$par)


  # ind_8mer <- 1
  # ind_7m8 <- grep("^8mer-mm.1$", Site_Names, perl=TRUE)
  # ind_7A1 <- grep("^8mer-mm.8$", Site_Names, perl=TRUE)
  # ind_6   <- grep("^8mer-mm.1mm.8$", Site_Names, perl=TRUE)
  # # print(ind_6)
  # # break
  # ind_8_1mm <- grep("^8mer-mm.[^18]$", Site_Names, perl=TRUE)
  # inds_change <- c(ind_8mer, ind_7m8, ind_7A1, ind_6, ind_8_1mm)
  # cols <- rep("green", nrow(counts) - 1)
  # cols[inds_change] <- cols_seed_1mm
  # if (mirna == "let-7a-21nt") {
  #   if (experiment == "equilibrium_mmseed_2_nb") {
  #     unmut_inds <- c(1:25, 80:88)    
  #   } else {
  #     unmut_inds <- c(1, 2, 4, 5:22, 23, 24, 80, 81, 86, 87)
  #   }
  # } else if (mirna == "let-7a" && experiment == "equilibrium_mmseed_nb") {
  #   unmut_inds <- c(1, 2, 4, 5:22, 23, 24, 80, 81, 86, 87)    
  # } else if (mirna == "miR-1") {
  #   if (experiment == "equilibrium_mmseed_2_nb") {
  #     unmut_inds <- c(1:25, 80:88)    
  #   } else {
  #     unmut_inds <- c(5:22)
  #   }
  # } else if (mirna == "miR-155") {
  #   unmut_inds <- c(5:22, 80:88)
  # }

  # Initialize the parameters for the unmut counts
}

# pars <- FitData("let-7a")
# pars <- FitData("miR-1")

pars <- FitData("miR-1")

pars_final <- c(1, exp(pars))

pars_comp <- pars_final[1:18]

pars_supp <- pars_final[19:34]

pars_comp_comb <- sapply(1:6, function(ind_start) {
  stop_ind <- 3*ind_start
  start_ind <- stop_ind - 2
  sum(pars_comp[start_ind:stop_ind])
})

pars_supp_comb <- c(pars_supp[1],
                   sum(pars_supp[2:4]),
                   sum(pars_supp[5:7]),
                   sum(pars_supp[8:16]))

print(mean(pars_comp_comb))
print(mean(pars_supp_comb))
print(mean(pars_supp_comb)/mean(pars_comp_comb))
print(GeoMean(pars_supp_comb)/GeoMean(pars_comp_comb))

plot(1:length(pars_supp), pars_supp)

break





pars_extended <- c(exp(opt4$par), "None"=0)
M_use <- M_list[[mirna2]]
M <- matrix(pars_extended[M_use], nrow=nrow(M_use))
M[nrow(M3), 26:(nrow(M3) - 1)] <- 18*mean(exp(pars))
rownames(M) <- rownames(M3)
colnames(M) <- colnames(M3)
start_counts <- rep(0, nrow(counts1))
names(start_counts) <- rownames(counts1)[length(start_counts)]
counts <- counts2[1:length(start_counts), 1, drop=FALSE]
if (mirna2 == "miR-155") {
  inds <- grep("^8mer-mm.1mm.8$", rownames(M), perl=TRUE)
  inds_use <- c(5:22, inds)
} else {
  inds_use <- 5:22
}
start_counts[inds_use] <- counts[inds_use, ]
model <- (diag(1 - colSums(M)) + M)%*%(diag(1 - colSums(M)) + M)%*%start_counts + pc
counts <- counts + pc
model <- model/sum(model)
cost <- sum(-1*counts[-nrow(M), 1]*log(model[-nrow(M)]))
tick <<- tick + 1



ind_8mer <- 1
ind_7m8 <- grep("^8mer-mm.1$", rownames(model), perl=TRUE)
ind_7A1 <- grep("^8mer-mm.8$", rownames(model), perl=TRUE)
ind_6   <- grep("^8mer-mm.1mm.8$", rownames(model), perl=TRUE)
ind_8_1mm <- grep("8mer-mm.[^18]$", rownames(model), perl=TRUE)

inds_change <- c(ind_8mer, ind_7m8, ind_7A1, ind_6, ind_8_1mm)

cols <- rep("green", nrow(model) - 1)

cols[inds_change] <- cols_seed_1mm


  # cols_seed_2mm_and_none <- rep(
  #   c("green", "black"),
  #   times=c(nrow(counts1) - 1 - 3 - 3 - 9 - 18 - 1, 1)
  # )
  


plot(model[-nrow(M)], counts[-nrow(M), 1], log='xy',
   xlim=c(1e-8, 1), ylim=c(1e-8, 1), col=cols, pch=20)
segments(1e-8, 1e-8, x1=1, y1=1, lty=2)
identify(model[-nrow(M)], counts[-nrow(M), 1], labels=rownames(M_use))


break


# mirna1 <- "let-7a"
# mirna2 <- "let-7a"

experiment1 <- "equilibrium_mmseed_nb"
experiment2 <- experiment1

Determine8merReferenceSiteType <- function(mirna, kmer_vec) {
  # Make the vector of letters in the 8mer site type.
  site_8mer <- StrRev(GetSiteSeq(mirna, "8mer"))
  site_8mer_vec <- unlist(strsplit(site_8mer, split=""))
  # Apply a loop over the kmers.
  site_vec <- sapply(kmer_vec, function(kmer) {
    # Condition in case the site is "None"
    if (kmer == "None") {
      return(kmer)
    # Condition in which the kmer is actually a kmer.
    } else {
      suffix <- ""
      kmer_vec <- unlist(strsplit(StrRev(kmer), split=""))
      for (i in 1:8) {
        if (site_8mer_vec[i] != kmer_vec[i]) {
          suffix <- sprintf("%smm%s%s", suffix, kmer_vec[i], i)
        } 
      }
      if (suffix != "") {
        suffix <- sprintf("-%s", suffix)
      }
      out <- sprintf("8mer%s", suffix)
      return(out)
    }
  })
  return(site_vec)
}

counts1 <- SubfunctionCall(GetProgrammedKmerFrequencies, mirna=mirna1,
                           experiment=experiment1, condition="I")
counts2 <- SubfunctionCall(GetProgrammedKmerFrequencies, mirna=mirna2,
                           experiment=experiment2, condition="I")
# Determine the names of each of the possible 8mer sites (1), 7mer-m8 sites 
# (3), 7mer-A1 sites (3), 6mer sites (9), and 8mer-mm sites (18), for both
# miRNAs.



sites_alt <- Determine8merReferenceSiteType(mirna2, rownames(counts1))
sites_alt <- sites_alt[-length(sites_alt)]

GetOverlapMatrix <- function(sites_alt) {
  sites_single_mismatch <- sites_alt[5:(18 + 4)]
  # print(sites_single_mismatch)
  mat_nonzero <- sapply(sites_single_mismatch, function(site) {
    # print(site)
    suffix <- unlist(strsplit(site, split="-"))[2]
    indeces <- grep(suffix, sites_alt)
    # print(indeces)
    indeces[1] <- 1
    row <- rep(0, length(sites_alt))
    row[indeces] <- 1
    row 
  })
  full_matrix <- matrix(rep(0, length(sites_alt)^2), nrow=length(sites_alt),
                        ncol=length(sites_alt))
  full_matrix[, 5:(18 + 4)] <- mat_nonzero
  return(full_matrix)
}

GetOverlapMatrix2 <- function(mirna, sites_alt, nucs_dict) {
  # First make the final matrix, pre-allocated with all zeros.
  full_matrix <- matrix(rep(0, length(sites_alt)^2), nrow=length(sites_alt),
                      ncol=length(sites_alt))
  temp_col <- rep(0, length(sites_alt))
  # Make the list of miRNA letters.
  site_8mer_vec <- unlist(strsplit(StrRev(GetSiteSeq(mirna, "8mer")), split=""))
  sites_single_mismatch <- sites_alt[2:25]
  col_8mer <- temp_col
  sites_8mer1mm_pos <- gsub("8mer-mm", replacement="", x=sites_single_mismatch)
  sites_8mer1mm_pos <- as.integer(substring(sites_8mer1mm_pos, first=2))
  col_8mer[2:25] <- nucs_dict[site_8mer_vec[sites_8mer1mm_pos]]
  # print(col_8mer)
  full_matrix[, 1] <- col_8mer
  mat_nonzero <- sapply(sites_single_mismatch, function(site) {
    suffix <- unlist(strsplit(site, split="-"))[2]
    mm <- substring(suffix, first=3, last=3)
    pos <- as.integer(substring(suffix, first=4))
    indeces <- grep(suffix, sites_alt)
    # Series of three string operations intended to determine the identify of
    # the mutated nucleotide.
    values <- grep(suffix, sites_alt, value=TRUE)
    values <- gsub("8mer-", replacement="", x=values)[-1]
    values <- substring(
      gsub(suffix, replacement="", x=values), first=3, last=3
    )
    values <- c(site_8mer_vec[pos], values)
    # This needs to be pre-appended to the identity of the nucleotide that the
    # mismatch itself is, since mutation back makes it into the 8mer.
    indeces[1] <- 1
    # row <- rep(0, length(sites_alt))
    temp_col[indeces] <- nucs_dict[values]
    temp_col
  })

  full_matrix[, 2:25] <- mat_nonzero

  double_matrix <- sapply(26:length(sites_alt), function(i_col) {
    site <- sites_alt[i_col]
    site <- gsub("8mer-", replacement="", x=site)
    double_site <- unlist(strsplit(site, split="mm"))[-1]
    pos <- as.integer(rev(substring(double_site, first=2)))
    vals <- nucs_dict[site_8mer_vec[pos]]
    inds <- sapply(double_site, function(site) {
      which(sites_alt == sprintf("8mer-mm%s", site))
    })
    temp_col[inds] <- vals
    temp_col
  })
  full_matrix[, 26:length(sites_alt)] <- double_matrix
  return(full_matrix)
}

GetOverlapMatrix3 <- function(mirna, sites_alt) {
  # First make the final matrix, pre-allocated with all zeros.
  full_matrix <- matrix(rep("None", (length(sites_alt) + 1)^2),
                        nrow=length(sites_alt) + 1,
                        ncol=length(sites_alt) + 1)
  temp_col <- rep("None", length(sites_alt) + 1)
  # Make the list of miRNA letters.
  site_8mer_vec <- unlist(strsplit(StrRev(GetSiteSeq(mirna, "8mer")), split=""))
  # Make the first column, which describes mutation from the 8mer site
  # to one of the 18 1-nt mismatches.
  sites_single_mismatch <- sites_alt[2:25]
  col_8mer <- temp_col
  s8mer1mm <- gsub("8mer-mm", replacement="", x=sites_single_mismatch)
  pos <- as.integer(substring(s8mer1mm, first=2))
  nuc_start <- site_8mer_vec[pos]
  nuc_end <- substring(s8mer1mm, first=1, last=1)
  nuc_combined <- sprintf("%s%s", nuc_start, nuc_end)
  col_8mer[2:25] <- nuc_combined
  full_matrix[, 1] <- col_8mer
  # Make the next 24 column, which describes mutation from one of the 8mer-1mm
  # sites into either the 8mer (the first entry) or 
  # to one of the 18 1-nt mismatches.
  cols_1mm <- sapply(sites_single_mismatch, function(site) {
    # Get the 1mm suffix, allocate the mismatch nucleotide and pos.
    suffix <- unlist(strsplit(site, split="-"))[2]
    mm1_nuc <- substring(suffix, first=3, last=3)
    mm1_pos <- as.integer(substring(suffix, first=4))
    # Get the indeces of which column entries will used. Note: The first is
    # incorrect and will be changed to `1`, because each 1mm-site can mutate to
    # the 8mer.
    inds <- grep(suffix, sites_alt)
    inds[1] <- 1
    # Series of three string operations intended to determine the identify of
    # the mutated nucleotide.
    mm2_strings <- gsub("8mer-", replacement="",
                        x=grep(suffix, sites_alt, value=TRUE))[-1]
    # Remove the starting mismatch from these strings to yield the second
    # mismatch
    mm2 <- gsub(suffix, replacement="", x=mm2_strings)
    mm2_pos <- as.integer(substring(mm2, first=4, last=4))
    mm2_nuc <- substring(mm2, first=3, last=3)
    # mm_second_end_nuc <- 
    nuc_start <- c(mm1_nuc, site_8mer_vec[mm2_pos])
    nuc_end <- c(site_8mer_vec[mm1_pos], mm2_nuc)
    nuc_combined <- sprintf("%s%s", nuc_start, nuc_end)
    # temp_col[inds] <- dinucs_dict[nuc_combined]
    temp_col[inds] <- nuc_combined
    temp_col
  })
  # Allocate these columns
  full_matrix[, 2:25] <- cols_1mm
  # Now determine the indeces for the remaining double sites mutating to single
  # sites.
  double_matrix <- sapply(26:length(sites_alt), function(i_col) {
    # First deal with reversion
    mm2 <- gsub("8mer-", replacement="", x=sites_alt[i_col])
    mm2_vec <- unlist(strsplit(mm2, split="mm"))[-1]
    mm2_pos <- as.integer(rev(substring(mm2_vec, first=2)))
    nuc_start <- rev(substring(mm2_vec, first=1, last=1))
    nuc_end <- site_8mer_vec[mm2_pos]
    nuc_combined <- sprintf("%s%s", nuc_start, nuc_end)
    inds <- sapply(mm2_vec, function(site) {
      which(sites_alt == sprintf("8mer-mm%s", site))
    })
    temp_col[inds] <- nuc_combined
    # Now deal with mutation to one of the other two for each.
    gsub_1 <- sprintf("8mer-mm.%smm%s", mm2_pos[2], mm2_vec[2])
    gsub_2 <- sprintf("8mer-mm%smm.%s", mm2_vec[1], mm2_pos[1])
    nuc_start <- rep(c(substring(mm2, first=3, last=3),
                       substring(mm2, first=7, last=7)), each=2)
    inds <- setdiff(c(grep(gsub_1, sites_alt),
                      grep(gsub_2, sites_alt)), c(i_col))
    other_mm_sites <- gsub("8mer-", replacement="", x=sites_alt[inds])
    nuc_end <- c(substring(other_mm_sites[1:2], first=3, last=3),
                 substring(other_mm_sites[3:4], first=7, last=7))
    nuc_combined <- sprintf("%s%s", nuc_start, nuc_end)
    temp_col[inds] <- nuc_combined
    temp_col
  })
  full_matrix[, 26:length(sites_alt)] <- double_matrix
  return(full_matrix)
}


# M1 <- GetOverlapMatrix(sites_alt)

pars_init <- c(`A`=-5, `C`=-5.1, `G`=-5.2, `T`=-5.3)

pc <- 1e-8

M2 <- GetOverlapMatrix2(mirna2, sites_alt, exp(pars_init))

pars_dinuc_init <- rep(-5, 16)
names(pars_dinuc_init) <- kDinucs


pars_dinuc_init <- pars_dinuc_init[c(-1, -4, -13, -16)]
print(pars_dinuc_init)

tick <- 0

M3 <- GetOverlapMatrix3(mirna2, sites_alt)

CostM3 <- function(pars) {
  pars_extended <- c(exp(pars), "None"=0)
  M <- matrix(pars_extended[M3], nrow=nrow(M3))
  M[nrow(M3), 26:(nrow(M3) - 1)] <- 18*mean(exp(pars))
  start_counts <- rep(0, nrow(counts1))
  names(start_counts) <- rownames(counts1)[length(start_counts)]
  counts <- counts2[1:length(start_counts), 1, drop=FALSE]
  start_counts[5:22] <- counts[5:22, ]
  model <- (diag(1 - colSums(M)) + M)%*%(diag(1 - colSums(M)) + M)%*%start_counts + pc
  counts <- counts + pc
  model <- model/sum(model)
  cost <- sum(-1*counts[-nrow(M), 1]*log(model[-nrow(M)]))
  tick <<- tick + 1
  if (tick%%100 == 0) {
    plot(model[-nrow(M)], counts[-nrow(M), 1], log='xy',
       xlim=c(1e-8, 1), ylim=c(1e-8, 1))
    segments(1e-8, 1e-8, x1=1, y1=1, lty=2)
    print(cost)
  }
  cost
}

opt <- optim(par=pars_dinuc_init, fn=CostM3, method="BFGS")
opt2 <- optim(par=opt$par, fn=CostM3, method="BFGS")
opt3 <- optim(par=opt2$par, fn=CostM3, method="BFGS")
opt4 <- optim(par=opt3$par, fn=CostM3, method="BFGS")

pars <- opt4$par

pars_extended <- c(exp(pars), "None"=0)
M <- matrix(pars_extended[M3], nrow=nrow(M3))
M[nrow(M3), 26:(nrow(M3) - 1)] <- 18*mean(exp(pars))
start_counts <- rep(0, nrow(counts1))
names(start_counts) <- rownames(counts1)[length(start_counts)]
counts <- counts2[1:length(start_counts), 1, drop=FALSE]
start_counts[5:22] <- counts[5:22, ]
model <- (diag(1 - colSums(M)) + M)%*%(diag(1 - colSums(M)) + M)%*%start_counts + pc
counts <- counts + pc
model <- model/sum(model)
cost <- sum(-1*counts[-nrow(M), 1]*log(model[-nrow(M)]))
tick <<- tick + 1
plot(model[-nrow(M)], counts[-nrow(M), 1], log='xy',
   xlim=c(1e-8, 1), ylim=c(1e-8, 1))
segments(1e-8, 1e-8, x1=1, y1=1, lty=2)
identify(model[-nrow(M)], counts[-nrow(M), 1], labels=sites_alt[-nrow(M)])



break
M <- GetOverlapMatrix2("let-7a", sites_alt, opt$par)


start_counts <- rep(0, nrow(counts1) - 1)
names(start_counts) <- rownames(counts1)[length(start_counts)]

start_counts[5:22] <- counts[5:22, ]
mut <- 0.001

model <- (diag(nrow(counts)) + M)%*%start_counts + pc
counts <- counts2[1:length(start_counts), 1, drop=FALSE] + pc


graphics.off()
dev.new(xpos=20, ypos=320, height=4, width=4)
print(model)



break


sites1 <- c(GetAll8mersCanonicalSite(mirna1), GetAll8merMmSites(mirna1))
sites2 <- c(GetAll8mersCanonicalSite(mirna2), GetAll8merMmSites(mirna2))

print(cbind(sites_alt[1:10], sites1[1:10]))

break

kmers1 <- sapply(sites1, GetSiteSeq, mirna=mirna1)
kmers2 <- sapply(sites2, GetSiteSeq, mirna=mirna2)
print(counts1[kmers1, ])
print(counts2[kmers2, ])
cols <- sapply(c(rep(kSiteColors[c("8mer", "7mer-m8", "7mer-A1", "6mer", "8mer-mmT6")],
            times=c(1, 3, 3, 9, 18)), "gray"), ConvertRColortoRGB, alpha=0.5)
x <- counts1[c(kmers1, "None"), 1]
y <- counts2[c(kmers2, "None"), 1]

out <- cbind(x, y, c(kmers1, "None"), c(kmers2, "None"))
rownames(out) <- c(sites1, "None")

print(out)

break


# PlotPairwiseInputFrequencies("let-7a-21nt", "equilibrium_mmseed_2_nb",
#                              "miR-1", "equilibrium_mmseed_2_nb", xpos=520)

# PlotPairwiseInputFrequencies("let-7a", "equilibrium_mmseed_nb",
#                              "let-7a_miR-155", "equilibrium_mmseed_nb")




break

cols_use <- kSiteColors[rep(c("8mer", "7mer-m8", "7mer-A1", "6mer"), times=c(1, 3, 3, 9))]
graphics.off()
dev.new(xpos=20, ypos=20, height=5, width=5)
plot(l7_comp_1[kmers_l7, 1], l7_comp_2[kmers_l7, 1], log='xy',
     xlim=c(1e-7, 1e-2), ylim=c(1e-7, 1e-2), pch=1, col=cols_use)
points(l7_comp_1[kmers_l7, 1], l7_comp_supmix[kmers_l7, 1], pch=20, col=cols_use)



break



print(kmer_table)




break


MakeFig0KdTable <- function() {
  base_path <- "/lab/solexa_bartel/mcgeary/AgoRBNS/GEO_McGearyLin2019/processed_files/agorbns_miR-1_nt"
  
  file_name <- "/lab/solexa_bartel/klin/miRNA_models_data/kds/kds.txt"
  file <- fread(file_name, sep="\t")
  print(head(file))
  file_keep <- which(file$mir == "mir1")
  file <- data.frame(file[file_keep, ])
  print(head(file))
  rownames(file) <- unlist(file[, 1])
  print(head(file))
  file_global <<- file
  kd_list <- c("AAACATTCCAGG", "AAACATTCCAGT", "AAACATTCCATA", "AAACATTCCATC", "AAACATTCCATG",
               "ACTCATTCCTCT", "ACTCATTCCTGA", "ACTCATTCCTGC", "ACTCATTCCTGG", "ACTCATTCCTGT",
               "GGACATTCCAGC", "GGACATTCCAGG", "GGACATTCCAGT", "GGACATTCCATA", "GGACATTCCATC")
  
  kds_figure <- file[kd_list,]
  kds_figure
}

out <- MakeFig0KdTable()
break

# MakeGEOMpraProcessedDataFiles()

# counts1 <- MakeReporterMatrix(1)
# counts2 <- MakeReporterMatrix(2)

# print(head(counts1))
# print(head(counts2))

break



MakeGEOMpraProcessedDataFiles()
break
# dir1 <- "/lab/solexa_bartel/mcgeary/AgoRBNS/"
# dir2 <- "/equilibrium/site_counts/4_5_12mers_"

# list_lin_means <- c()
# list_log_means <- c()

# list_max <- c()

# for (mirna in c("let-7a", "miR-124")) {
#   print(mirna)
#   site <- GetSiteSeq(mirna, "8mer")
#   print(site)
#   for (nt_range in c("1-4", "2-5", "3-6", "4-7", "5-8")) {
#     print(nt_range)
#     path <- sprintf("%s%s%s%s.txt", dir1, mirna, dir2, nt_range)
#     print(path)
#     kmers <- read.table(path, row.names=1, header=FALSE, sep="\t")
#     len_use <- nrow(kmers) - 1
#     inds_use <- grep(sprintf("..%s..", site), rownames(kmers), perl=TRUE)
#     list_lin_means <- c(list_lin_means, mean(kmers[inds_use, 1]))
#     list_log_means <- c(list_log_means, exp(mean(log(kmers[inds_use, 1]))))
#     list_max <- c(list_max, max(kmers[1:len_use, 1]))
#   }
# }
# print(range(list_lin_means))
# print(range(list_log_means))

# print(list_max)

# 1/(4^12)*20e6*(37 - 12 + 1)

# break

# read_lines <- c(49268308, 49268308, 10990446, 7883321, 11859706,
# 13707087
# 11859706
# 13230928
#  6321922


sites <- SitesXCounts("let-7a", "equilibrium_mmseed_nb", "0", "programmed")
sites_2 <- SitesXCounts("let-7a-21nt", "equilibrium_mmseed_nb", "0", "programmed")


sites_8mer_mmA7 <- sites[grep("^.*\\|.*\\|8mer-mmA7$", rownames(sites), perl=TRUE), ]


# sites <- t(t(sites)/colSums(sites))
# sites <- sites/sites[, 1]

# sites_2 <- t(t(sites_2)/colSums(sites_2))
# sites_2 <- sites_2/sites_2[, 1]


inds_5mer_1 <- grep("^9mer-m11.19\\|14\\|", rownames(sites), perl=TRUE)
inds_5mer_2 <- grep("^9mer-m11.19\\|14\\|", rownames(sites_2), perl=TRUE)



double_sites_1 <- grep("&", rownames(sites))
double_sites_2 <- grep("&", rownames(sites_2))

single_sites_1 <- setdiff(rownames(sites), rownames(sites)[double_sites_1])
single_sites_2 <- setdiff(rownames(sites_2), rownames(sites_2)[double_sites_2])



print(colSums(sites[double_sites_1, ])/colSums(sites))
print(colSums(sites_2[double_sites_2, ])/colSums(sites_2))



print(colSums(sites))
print(colSums(sites_2))

print(nrow(sites))
print(nrow(sites_2))

print(length(single_sites_1))
print(length(single_sites_2))

print(length(double_sites_1))
print(length(double_sites_2))


temp_kds <- read.table("ThreePrimeTargetingPaper/temp_let-7a_kds.txt", sep="\t",
                       row.names=1, header=TRUE)

temp_kds_2 <- read.table("ThreePrimeTargetingPaper/temp_let-7a-21nt_kds.txt", sep="\t",
                       row.names=1, header=TRUE)


rownames(temp_kds) <- gsub("^6mer\\|28\\|", x=rownames(temp_kds), perl=TRUE, replacement="7mer-m8\\|28\\|")
rownames(temp_kds_2) <- gsub("^6mer\\|28\\|", x=rownames(temp_kds_2), perl=TRUE, replacement="7mer-m8\\|28\\|")

rownames(temp_kds) <- gsub("^7mer-A1\\|27\\|", x=rownames(temp_kds), perl=TRUE, replacement="8mer\\|27\\|")
rownames(temp_kds_2) <- gsub("^7mer-A1\\|27\\|", x=rownames(temp_kds_2), perl=TRUE, replacement="8mer\\|27\\|")



sites_1 <- rownames(temp_kds)
sites_2 <- rownames(temp_kds_2)

sites_shared <- intersect(sites_1, sites_2)
sites_shared_single <- grep("&", sites_shared, invert=TRUE, value=TRUE)

kds_1 <- temp_kds[sites_shared, 1]
kds_2 <- temp_kds_2[sites_shared, 1]

kds_single_1 <- temp_kds[sites_shared_single, 1]
kds_single_2 <- temp_kds_2[sites_shared_single, 1]


kds_paired <- t(rbind(kds_1, kds_2))

kds_paired_single <- t(rbind(kds_single_1, kds_single_2))


print(head(kds_paired_single))
print(cor(kds_paired[, 1], kds_paired[, 2])^2)
print(cor(kds_paired_single[, 1], kds_paired_single[, 2])^2)




rownames(kds_paired) <- sites_shared
rownames(kds_paired_single) <- sites_shared_single
print(head(kds_paired))

inds_8 <- grep("^8mer\\|.*\\|.*_Kd", sites_shared, value=TRUE, perl=TRUE)
inds_7m8 <- grep("^7mer-m8\\|.*\\|.*_Kd", sites_shared, value=TRUE, perl=TRUE)
inds_7a1 <- grep("^7mer-A1\\|.*\\|.*_Kd", sites_shared, value=TRUE, perl=TRUE)
inds_6 <- grep("^6mer\\|.*\\|.*_Kd", sites_shared, value=TRUE, perl=TRUE)
inds_6m8 <- grep("^6mer-m8\\|.*\\|.*_Kd", sites_shared, value=TRUE, perl=TRUE)
inds_6a1 <- grep("^6mer-A1\\|.*\\|.*_Kd", sites_shared, value=TRUE, perl=TRUE)

dev.new(xpos=20, ypos=20, height=5, width=5)
plot(kds_paired_single[, 1], kds_paired_single[, 2], col=rgb(0, 0, 0, alpha=0.2), pch=20)

points(kds_paired_single[1:18, 1], kds_paired_single[1:18, 2], col=ConvertRColortoRGB("forestgreen", alpha=0.2), pch=20)
points(kds_paired_single[inds_8, 1], kds_paired_single[inds_8, 2], col=ConvertRColortoRGB(kSiteColors["8mer"], alpha=0.2), pch=20)
points(kds_paired_single[inds_7m8, 1], kds_paired_single[inds_7m8, 2], col=ConvertRColortoRGB(kSiteColors["7mer-m8"], alpha=0.2), pch=20)
points(kds_paired_single[inds_7a1, 1], kds_paired_single[inds_7a1, 2], col=ConvertRColortoRGB(kSiteColors["7mer-A1"], alpha=0.2), pch=20)
points(kds_paired_single[inds_6, 1], kds_paired_single[inds_6, 2], col=ConvertRColortoRGB(kSiteColors["6mer"], alpha=0.2), pch=20)
points(kds_paired_single[inds_6m8, 1], kds_paired_single[inds_6m8, 2], col=ConvertRColortoRGB(kSiteColors["6mer-m8"], alpha=0.2), pch=20)
points(kds_paired_single[inds_6a1, 1], kds_paired_single[inds_6a1, 2], col=ConvertRColortoRGB(kSiteColors["6mer-A1"], alpha=0.2), pch=20)






GetSeedKdMatrix <- function(site, kd_data, xpos=20, ypos=20) {
  grep_query <- sprintf("^%s\\|.*\\|8mer-mm[ACGT]._Kd$", site)
  message(grep_query)
  print(head(kd_data))
  inds_kd <- grep(grep_query, rownames(kd_data), perl=TRUE)
  print(inds_kd)
  site_names <- rownames(kd_data)[inds_kd]
  mismatch_seeds <- sapply(site_names, function(name) {
    split_name <- unlist(strsplit(name, split="\\|"))[2:3]
  })
  lib_positions <- range(as.integer(mismatch_seeds[1, ]))
  print(lib_positions)
  programmed_mm <- sort(unique(mismatch_seeds[2, ]))
  print(programmed_mm)
  out_matrix <- matrix(NaN, nrow=length(programmed_mm),
                       ncol=lib_positions[2] - lib_positions[1] + 1,
                       dimnames=list(programmed_mm,
                                     as.character(seq(lib_positions[2],
                                                      lib_positions[1]))))
  for (site_name in site_names) {
    col_name <- mismatch_seeds[1, site_name]
    mm_name <- mismatch_seeds[2, site_name]
    print(col_name)
    print(mm_name)
    out_matrix[mm_name, col_name] <- kd_data[site_name, 1]
  }
  # print(out_matrix)
  # dev.new(xpos=xpos, ypos=ypos, height=3, width=3.5)
  # par(mar=c(1, 1, 1, 1))
  # image(t(out_matrix))
  out_matrix
}

mat_8mer_1 <- GetSeedSiteKd("5mer-m8", temp_kds)
mat_8mer_2 <- GetSeedSiteKd("6mer-m8", temp_kds_2)

average_mm_1 <- rowMeans(mat_8mer_1, na.rm=TRUE)
average_mm_2 <- rowMeans(mat_8mer_2, na.rm=TRUE)

average_pos_1 <- colMeans(mat_8mer_1, na.rm=TRUE)
average_pos_2 <- colMeans(mat_8mer_2, na.rm=TRUE)


names(average_mm_1) <- rownames(mat_8mer_1)
names(average_mm_2) <- rownames(mat_8mer_2)



dev.new(xpos=520, ypos=20, height=5, width=5)
plot(average_mm_1, average_mm_2)

dev.new(xpos=1020, ypos=20, height=5, width=5)
plot(average_pos_1, average_pos_2)

seed_sites_1 <- temp_kds[1:18, 1]
seed_sites_2 <- temp_kds_2[1:18, 1]

names(seed_sites_1) <- sapply(rownames(temp_kds)[1:18], function(name) {
  unlist(strsplit(name, split="\\|"))[3]
})

names(seed_sites_2) <- sapply(rownames(temp_kds_2)[1:18], function(name) {
  unlist(strsplit(name, split="\\|"))[3]
})


dev.new(xpos=20, ypos=520, height=5, width=5)
plot(average_mm_1, seed_sites_1[names(average_mm_1)], xlim=c(-3, 0), ylim=c(-3, 0))

dev.new(xpos=520, ypos=520, height=5, width=5)
plot(average_mm_2, seed_sites_2[names(average_mm_2)], xlim=c(-3, 0), ylim=c(-3, 0))

dev.new(xpos=1020, ypos=520, height=5, width=5)
plot(names(average_pos_1), average_pos_1, type="o", xlim=c(1, 30), ylim=c(-3, 0), pch=20)
lines(names(average_pos_1), average_pos_2, col="red", lwd=1, pch=20)





break




split_8mers <- grep("8mer\\|", rownames(temp_kds), perl=TRUE, value=TRUE)

pos_8mers <- sapply(split_8mers, function(name) {
  as.integer(unlist(strsplit(name, split="\\|", perl=TRUE))[2])
  })

split_7merm8s <- grep("7mer-m8\\|", rownames(temp_kds), perl=TRUE, value=TRUE)

pos_7merm8s <- sapply(split_7merm8s, function(name) {
  as.integer(unlist(strsplit(name, split="\\|", perl=TRUE))[2])
  })



kds_8mer <- temp_kds[grep("8mer\\|", rownames(temp_kds), perl=TRUE), 1]

kds_7merm8 <- temp_kds[grep("7mer-m8\\|", rownames(temp_kds), perl=TRUE), 1]

kds_7merA1 <- temp_kds[grep("7mer-A1\\|", rownames(temp_kds), perl=TRUE), 1]

kds_6mer <- temp_kds[grep("6mer\\|", rownames(temp_kds), perl=TRUE), 1]

kds_6merm8 <- temp_kds[grep("6mer-m8\\|", rownames(temp_kds), perl=TRUE), 1]

kds_6merA1 <- temp_kds[grep("6mer-A1\\|", rownames(temp_kds), perl=TRUE), 1]

print(mean(kds_8mer))
print(sd(kds_8mer))
print(mean(kds_7merm8))
print(sd(kds_7merm8))
print(mean(kds_7merA1))
print(sd(kds_7merA1))
print(mean(kds_6mer))
print(sd(kds_6mer))
print(mean(kds_6merm8))
print(sd(kds_6merm8))
print(mean(kds_6merA1))
print(sd(kds_6merA1))


check <- EquilPars("let-7a")[1:6, 2]

kds_equil <- c(rep(check[1], length(kds_8mer)),
               rep(check[2], length(kds_7merm8)),
               rep(check[3], length(kds_7merA1)),
               rep(check[4], length(kds_6mer)),
               rep(check[5], length(kds_6merA1)),
               rep(check[6], length(kds_6merm8)))


cols_plot <- c(rep(kSiteColors["8mer"], length(kds_8mer)),
               rep(kSiteColors["7mer-m8"], length(kds_7merm8)),
               rep(kSiteColors["7mer-A1"], length(kds_7merA1)),
               rep(kSiteColors["6mer"], length(kds_6mer)),
               rep(kSiteColors["6mer-A1"], length(kds_6merA1)),
               rep(kSiteColors["6mer-m8"], length(kds_6merm8)))

print(kds_equil)

dev.new(xpos=20, ypos=20, height=5, width=5)


kds_all_1 <- 10^c(kds_8mer, kds_7merm8, kds_7merA1, kds_6mer, kds_6merA1, kds_6merm8)

plot(kds_equil, 10^c(kds_8mer, kds_7merm8, kds_7merA1, kds_6mer, kds_6merA1, kds_6merm8), col=cols_plot, log='xy', xlim=c(0.0001, 10), ylim=c(0.0001, 10))
segments(0.0001, 0.0001, x1=10, y1=10)

dev.new(xpos=520, ypos=20, height=5, width=5)

temp_kds <- read.table("ThreePrimeTargetingPaper/temp_let7-21nt_kds.txt", sep="\t",
                       row.names=1, header=TRUE)


rownames(temp_kds) <- gsub("^6mer\\|27\\|", x=rownames(temp_kds), perl=TRUE, replacement="7mer-m8\\|27\\|")

rownames(temp_kds) <- gsub("^7mer-A1\\|26\\|", x=rownames(temp_kds), perl=TRUE, replacement="8mer\\|26\\|")

split_8mers <- grep("8mer\\|", rownames(temp_kds), perl=TRUE, value=TRUE)

pos_8mers <- sapply(split_8mers, function(name) {
  as.integer(unlist(strsplit(name, split="\\|", perl=TRUE))[2])
  })

split_7merm8s <- grep("7mer-m8\\|", rownames(temp_kds), perl=TRUE, value=TRUE)

pos_7merm8s <- sapply(split_7merm8s, function(name) {
  as.integer(unlist(strsplit(name, split="\\|", perl=TRUE))[2])
  })



kds_8mer <- temp_kds[grep("8mer\\|", rownames(temp_kds), perl=TRUE), 1]

kds_7merm8 <- temp_kds[grep("7mer-m8\\|", rownames(temp_kds), perl=TRUE), 1]

kds_7merA1 <- temp_kds[grep("7mer-A1\\|", rownames(temp_kds), perl=TRUE), 1]

kds_6mer <- temp_kds[grep("6mer\\|", rownames(temp_kds), perl=TRUE), 1]

kds_6merm8 <- temp_kds[grep("6mer-m8\\|", rownames(temp_kds), perl=TRUE), 1]

kds_6merA1 <- temp_kds[grep("6mer-A1\\|", rownames(temp_kds), perl=TRUE), 1]

print(mean(kds_8mer))
print(sd(kds_8mer))
print(mean(kds_7merm8))
print(sd(kds_7merm8))
print(mean(kds_7merA1))
print(sd(kds_7merA1))
print(mean(kds_6mer))
print(sd(kds_6mer))
print(mean(kds_6merm8))
print(sd(kds_6merm8))
print(mean(kds_6merA1))
print(sd(kds_6merA1))


check <- EquilPars("let-7a")[1:6, 2]

kds_equil <- c(rep(check[1], length(kds_8mer)),
               rep(check[2], length(kds_7merm8)),
               rep(check[3], length(kds_7merA1)),
               rep(check[4], length(kds_6mer)),
               rep(check[5], length(kds_6merA1)),
               rep(check[6], length(kds_6merm8)))


cols_plot <- c(rep(kSiteColors["8mer"], length(kds_8mer)),
               rep(kSiteColors["7mer-m8"], length(kds_7merm8)),
               rep(kSiteColors["7mer-A1"], length(kds_7merA1)),
               rep(kSiteColors["6mer"], length(kds_6mer)),
               rep(kSiteColors["6mer-A1"], length(kds_6merA1)),
               rep(kSiteColors["6mer-m8"], length(kds_6merm8)))

print(kds_equil)


kds_all_2 <- 10^c(kds_8mer, kds_7merm8, kds_7merA1, kds_6mer, kds_6merA1, kds_6merm8)


plot(kds_equil, 10^c(kds_8mer, kds_7merm8, kds_7merA1, kds_6mer, kds_6merA1, kds_6merm8), col=cols_plot, log='xy', xlim=c(0.0001, 10), ylim=c(0.0001, 10))
segments(0.0001, 0.0001, x1=10, y1=10)

dev.new(xpos=1020, ypos=20, height=5, width=5)
plot(kds_all_1, kds_all_2, col=cols_plot, log='xy', xlim=c(0.0001, 10), ylim=c(0.0001, 10))


break


sites_both <- intersect(rownames(temp_kds), rownames(temp_kds_2))

print(head(temp_kds, n=18))
print(head(temp_kds_2, n=18))

plot(temp_kds[sites_both, ], temp_kds_2[sites_both, ])

cor(temp_kds[sites_both, ], temp_kds_2[sites_both, ])


break

temp_kds <- 10^temp_kds

prog_sites <- unique(sapply(rownames(temp_kds), function(site_name) {
  unlist(strsplit(site_name, split="\\|"))[3]
}))[1:18]


print(prog_sites)


output_table <- matrix(0, nrow=18, ncol=19,
                       dimnames=list(prog_sites, as.character(seq(27, 9))))

print(output_table)


for (prog_site in prog_sites) {
  inds <- grep(sprintf("^7mer-m11.17\\|.*\\|%s", prog_site),
               rownames(temp_kds), perl=TRUE)
  kds_prog_site <- temp_kds[inds, , drop=FALSE]

  col_names <- sapply(rownames(kds_prog_site), function(name) {
    unlist(strsplit(name, split="\\|"))[2]
    })
  print(prog_site)
  print(col_names)
  col_names <<- col_names
  output_table[prog_site, col_names] <- kds_prog_site[, 1]
}

output_table_means <- exp(colMeans(log(output_table)))

plot(colnames(output_table), c(output_table_means), log='y')

# print(output_table)

break


inds_1 <- grep("^5mer-m11.15\\|.*\\|8mer-mmC7", rownames(temp_kds), perl=TRUE)
inds_2 <- grep("^5mer-m11.15\\|.*\\|8mer-mmA7", rownames(temp_kds), perl=TRUE)

inds_3 <- grep("^8mer\\|.*\\|8mer-mmC7", rownames(temp_kds), perl=TRUE)
inds_4 <- grep("^8mer\\|.*\\|8mer-mmA7", rownames(temp_kds), perl=TRUE)



kds_check_1 <- temp_kds[inds_1, , drop=FALSE]
kds_check_2 <- temp_kds[inds_2, , drop=FALSE]
kds_check_3 <- temp_kds[inds_3, , drop=FALSE]
kds_check_4 <- temp_kds[inds_4, , drop=FALSE]


print(rownames(kds_check))

order_1 <- sapply(rownames(kds_check_1), function(name) {
  thing_1 <- unlist(strsplit(name, split="\\|"))[2]
  as.integer(thing_1)
})

order_2 <- sapply(rownames(kds_check_2), function(name) {
  thing_1 <- unlist(strsplit(name, split="\\|"))[2]
  as.integer(thing_1)
})

order_3 <- sapply(rownames(kds_check_3), function(name) {
  thing_1 <- unlist(strsplit(name, split="\\|"))[2]
  as.integer(thing_1)
})

order_4 <- sapply(rownames(kds_check_4), function(name) {
  thing_1 <- unlist(strsplit(name, split="\\|"))[2]
  as.integer(thing_1)
})





plot(order_3, 10^kds_check_3[, 1], type="l", xlim = c(30, 1), ylim=c(0.01, 1), log='y')

lines(order_2, 10^kds_check_2[, 1], col="red")

lines(order_1, 10^kds_check_1[, 1], col="purple")

lines(order_4, 10^kds_check_4[, 1], col="pink")





break

PlotPositionalEnrichmentForProgramedLibrary(
  "let-7a", "equilibrium_mmseed_nb", 40, 0, 8, seedex=TRUE, ps=0.000001, toppos=5,
  pdf.plot="2.C",
)

break

# pc_I_nb <- read.table("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/ThreePrimeTargetingPaper/pythonfigures/Figure2/rawdata/I_kmer_counts.txt",
#                       sep="\t", row.names=1, header=TRUE)
# pc_A_nb <- read.table("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/ThreePrimeTargetingPaper/pythonfigures/Figure2/rawdata/A40_kmer_counts.txt",
#                       sep="\t", row.names=1, header=TRUE)




# PlotPositionalEnrichmentForProgramedLibrary(
#   "let-7a", "equilibrium_mmseed_nb", 40, 0, 8, seedex=TRUE, ps=0.01, toppos=3,
# )




# print(image(as.matrix(ps_R_4[names(output_sorted)[1:20], ])))


break

# output <- apply(npc_I, 1, function(row
# norm_counts_I <- total_counts_I / sum(total_counts_I)
# norm_counts_40 <- total_counts_40 / sum(total_counts_40)

R_40 <- norm_counts_40/norm_counts_I

names(R_40) <- rownames(position_counts_I)

print(R_40[1:10, drop=FALSE])

R_40_ordered <- R_40[order(-R_40)]

print(R_40_ordered[1:10])

break



# site_info_pyr_pyr <- GetKdProperties(kds_pyr_pyr, "let-7a")
# site_info_pur_pur <- GetKdProperties(kds_pur_pur, "let-7a")
# site_info_A_C <- GetKdProperties(kds_A_C, "let-7a")
# site_info_wobble <- GetKdProperties(kds_wobble, "let-7a")
site_info_bulge <- GetKdProperties(kds_bulge, "let-7a", bulge=TRUE)


# pyr_pyr_remaining <- grep("mer", site_info_pyr_pyr[, 1], invert=TRUE, value=TRUE)
# pur_pur_remaining <- grep("mer", site_info_pur_pur[, 1], invert=TRUE, value=TRUE)
# A_C_remaining <- grep("mer", site_info_A_C[, 1], invert=TRUE, value=TRUE)
# wobble_remaining <- grep("mer", site_info_wobble[, 1], invert=TRUE, value=TRUE)
bulge_remaining <- grep("mer", site_info_bulge[, 1], invert=TRUE, value=TRUE)




break



shared_1_2 <- intersect(rownames(test1), rownames(test2))
shared_2_3 <- intersect(rownames(test2), rownames(test3))
shared_3_4 <- intersect(rownames(test3), rownames(test4))
shared_4_5 <- intersect(rownames(test4), rownames(test5))
shared_5_1 <- intersect(rownames(test5), rownames(test1))




# col_all <- rgb(0, 0, 0, alpha=0.5)
# dev.new(pos=20, ypos=20, height=5, width=5)
# plot(test1[shared_1_2, ], test2[shared_1_2, ], col=col_all)
# dev.new(xpos=520, ypos=20, height=5, width=5)
# plot(test2[shared_2_3, ], test3[shared_2_3, ], col=col_all)
# dev.new(xpos=1020, ypos=20, height=5, width=5)
# plot(test3[shared_3_4, ], test4[shared_3_4, ], col=col_all)
# dev.new(xpos=20, ypos=520, height=5, width=5)
# plot(test4[shared_4_5, ], test5[shared_4_5, ], col=col_all)
# dev.new(xpos=520, ypos=520, height=5, width=5)
# plot(test5[shared_5_1, ], test1[shared_5_1, ], col=col_all)


break


# test <- EquilPars("miR-1", "equilibrium", sitelist="resubmissionfinal", buffer=TRUE, combined=FALSE, single=TRUE)
# test_nb <- EquilPars("miR-7-23nt", "equilibrium2_nb", sitelist="resubmissionfinal", single=TRUE)
# test_nb_bypass <- EquilPars("miR-7-23nt", "equilibrium2_nb", sitelist="resubmissionfinal", single=TRUE, AGOfixedbypass=TRUE)

# print(test_nb)

# plot(test_nb[, 2], test_nb_bypass[, 2], log='xy')
break

ln_test <- log(test)

low_bound <- ln_test[, 2] - ln_test[, 3]
high_bound <- ln_test[, 5] - ln_test[, 2]
error_average <-  (low_bound + high_bound)/2

print(cbind(test[, 2], exp(error_average)*test[, 2]))
break
# model <- GetFlankLinearModel()
# print(summary(model))
sXc <- SitesXCounts("miR-1", sitelist="resubmissionfinal", buffer=TRUE)

pars_all <- EquilPars("miR-1", sitelist="resubmissionfinal", buffer=TRUE,
                  combined=FALSE, single=TRUE)
par_names <- rownames(pars_all)
pars <- pars_all[, 1]
names(pars) <- par_names
names(pars)[length(pars) - 1] <- "bg"
names(pars)[length(pars)] <- "AGO"


kds <- pars[1:nrow(sXc)]
l <- 100*Norm(sXc[, 1])
L <- sum(l)

# a_R <- FreeAgoR(kds, l, 1)

a <- FreeAgo(kds, l, 1)

kds_check_1 <- kds_check
l_check_1 <- l_check
# a_C_alt <- FreeAgoCAlt(log(pars), l, L, length(l))
# print(a_R)
# print(a)
# print(a_C_alt)



# print(sXc)

# print(pars)

l <- Norm(sXc[, 1])*100
data <- sXc[, 3:7]


print(CostC(pars=c(log(pars)), data=unlist(data), dil=sapply(colnames(data), as.numeric)/100, l=l,
            L=sum(l), Y=colSums(data), n_i=nrow(data), n_j=ncol(data), n_mir=1, zero_grad=nrow(data), n_z=1,
            fixed=FALSE, upper_=rep(10, length(pars)), lower_=rep(-5, length(pars)),
            plot_=FALSE, plotname=NULL, tempname=NULL))



print(CostNew(log10(pars), sXc, combined=FALSE))

theta <- log(pars)

y <- data

print(CostMethodsCheck(theta, l, y))

grad_R <- GradientNew(log10(pars), sXc, combined=FALSE)

grad_C <- GradC(pars=c(log(pars)), data=unlist(data), dil=sapply(colnames(data), as.numeric)/100, l=l,
            L=sum(l), Y=colSums(data), n_i=nrow(data), n_j=ncol(data), n_mir=1, zero_grad=nrow(data), n_z=1,
            fixed=FALSE, upper_=rep(10, length(pars)), lower_=rep(-5, length(pars)),
            plot_=FALSE, plotname=NULL, tempname=NULL)


grad_R_approx <- grad(CostNew, log10(pars), sXc=sXc, combined=FALSE)

grad_C_check <- GradientMethodsCheck(theta, l, y)

plot(grad_C, grad_C_check)


lines(c(min(min(grad_C, min(grad_C_check))), max(max(grad_C), max(grad_C_check))),
      c(min(min(grad_C, min(grad_C_check))), max(max(grad_C), max(grad_C_check))))




break


CheckRegister <- function(pos) {
  command <- sprintf(
    "cut -c %s /lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium_mmseed_nb/reads/I.txt | sort | uniq -c",
  pos
  )
  out <- system(command, intern=TRUE)
  out <- sapply(out, function(str_i) {
    str_i <- strsplit(str_i, split=" ")[[1]]
    str_i <- str_i[which(str_i != "")]
    as.integer(str_i[1])
  })
  names(out) <- sapply(names(out), function(out_i) {
    split <- strsplit(out_i, split=" ")[[1]]
    split[length(split)]
  })
  out_fixed <- rep(0, 4)
  names(out_fixed) <- c("A", "C", "G", "T")
  for (name in names(out)) {
    out_fixed[name] <- out[name]
  }
  out_fixed
}

# window_pos <- cbind(sapply(1:40, CheckRegister))

print(window_pos)

site_8mer <- GetSiteSeq("let-7a", "8mer")

GetPos_mm <- function(site, pos) {
  pos_nuc <- substr(site, start=pos, stop=pos)
  print(pos_nuc)
  other_nucs <- setdiff(c("A", "C", "G", "T"), pos_nuc)
  print(other_nucs)
  vec_site <- strsplit(site, split="")[[1]]
  mm_sites <- sapply(other_nucs, function(mm_nuc) {
    vec_mm <- vec_site
    vec_mm[pos] <- mm_nuc
    paste0(vec_mm, collapse="")
  })
  mm_sites
}

test <- GetPos_mm(site_8mer, 2)

all_mismatches <- c(sapply(2:7, GetPos_mm, site=site_8mer))

print(all_mismatches)

print(test)

counts_out <- rep(0, length(all_mismatches))
counts_out_random <- rep(0, length(all_mismatches))
names(counts_out) <- all_mismatches
names(counts_out_random) <- all_mismatches

total_random <- as.integer(
  strsplit(system("wc -l /lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/reads/I.txt", intern=TRUE),
            split=" ")[[1]][1]
)
total_programmed <- as.integer(
  strsplit(system("wc -l /lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium_mmseed_nb/reads/I.txt", intern=TRUE),
            split=" ")[[1]][1]
)

for (mismatch in all_mismatches) {
  command <- sprintf(
    "grep %s /lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium_mmseed_nb/reads/I.txt | wc -l",
    mismatch
  )
  out <- as.integer(system(command, intern=TRUE))
  counts_out[mismatch] <- out/total_programmed
  command <- sprintf(
    "grep %s /lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/reads/I.txt | wc -l",
    mismatch
  )
  out <- as.integer(system(command, intern=TRUE))
  counts_out_random[mismatch] <- out/total_random
}

print(counts_out)
print(counts_out_random)
print(counts_out/counts_out_random)

break


PrintKdValues <- function(kd.matrix) {
    kd.matrix <- kd.matrix[1:(nrow(kd.matrix) - 3), ]
    kds <- kd.matrix$Mean
    # Part where text is formatted:
    kd_mags <- floor(log10(kds))
    # error_upper <- (kd.matrix$Upper_CI-kds)/(10^kd_mags)
    error_lower <- (kds - kd.matrix$Lower_CI)/(10^kd_mags)
    # error_mags_upper <- floor(log10(error_upper))
    error_mags_lower <- floor(log10(error_lower))
    text_matrix <- cbind(kds/(10^kd_mags), c(error_lower),
                             -error_mags_lower)
    rownames(text_matrix) <- rownames(kd.matrix)
    text_matrix_full <- cbind(text_matrix, kd_mags)

    print(text_matrix_full)
    sapply(1:nrow(text_matrix_full), function(i) {
      site <- rownames(text_matrix_full)[i]
      row <- text_matrix_full[i, ]
      print(sprintf("%s: %.*f +/- %.*f x 10^%s", site, row[3], row[1], row[3],
                    row[2], row[4]))
    })
}


test1 <- EquilPars("miR-124", sitelist="paperfinal", singleonly=TRUE)
test2 <- EquilPars("miR-124", sitelist="paperfinal")
test3 <- EquilPars("miR-124", sitelist="resubmissionfinal", singleonly=TRUE)
test4 <- EquilPars("miR-124", sitelist="resubmissionfinal")


print(test1["7mer-A1_Kd", 2]/test1["7mer-m8_Kd", 2])
print(test2["7mer-A1_Kd", 2]/test2["7mer-m8_Kd", 2])
print(test3["7mer-A1_Kd", 2]/test3["7mer-m8_Kd", 2])
print(test4["7mer-A1_Kd", 2]/test4["7mer-m8_Kd", 2])

break

flanks_test <- GetFlankKds("miR-1", "8mer", combined=TRUE, buffer=TRUE)[1:256, 2, drop=FALSE]

print(flanks_test)

min_test <- min(flanks_test)
max_test <- max(flanks_test)

min_flank <- rownames(flanks_test)[which(flanks_test == min_test)]
max_flank <- rownames(flanks_test)[which(flanks_test == max_test)]

print(min_flank)
print(max_flank)




# break


# max(1/GetFlankKds("miR-1", "8mer", combined=TRUE, buffer=TRUE))
# max(1/GetFlankKds("miR-1", "8mer", combined=TRUE, buffer=TRUE)[1:256, 2])


# PieChartOfColorRamp <- function(start_color, end_color, n_shades) {
#   color_func <- colorRampPalette(c(start_color, end_color))
#   shades <- color_func(n_shades)
#   slices <- rep(1/n_shades, n_shades)
#   dev.new(xpos=20, ypos=20, height=5, width=5)
#   pie(slices, labels = shades, col=shades)
# }



# PieChartOfColorRamp(kSiteColors["9mer-m11.19"], kSiteColors["11mer-m9.19"], 10)

sapply(kMirnas, function(mirna) {
  if (mirna == "miR-7-23nt") {
    experiment <- "equilibrium2_nb"
  } else {
    experiment <- "equilibrium"
  }
  if (mirna == "miR-1") {
    combined <- FALSE
    buffer <- TRUE
  } else {
    combined <- TRUE
    buffer <- FALSE
  }
  occupancies <- SubfunctionCall(GetSiteOccupancy, sitelist="resubmissionfinal",
                                 singleonly=TRUE, remove_sites=FALSE)
  print(occupancies)
  print(colSums(occupancies[kCanonicalSites,]/colSums(occupancies[-nrow(occupancies), ])))
})

break

labels <- c("8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-A1", "8mer-bT(4.6)", "8mer-w6",
  "GCTTCCGC", "8mer-mmC5", "7mer-A1bT(4.6)", "5mer-m2.6", "8mer-mmT6",
  "7mer-m8w6", "6mer-m8", "None")

              
labels_2 <- c("6mer-A1", "8mer-bT(4.6)", "8mer-w6", "8mer-mmC5",
              "7mer-A1bT(4.6)", "5mer-m2.6", "8mer-mmT6", "7mer-m8w6",
              "6mer-m8")

print(colSums(occupancies[labels_2, ]))
break
# test <- EquilPars("let-7a", sitelist="resubmissionfinal", singleonly=TRUE)
# PrintKdValues(test)

break
# test <- EquilPars("miR-155", sitelist="resubmissionfinal", singleonly=TRUE)
# PrintKdValues(test)

# test <- EquilPars("lsy-6", sitelist="resubmissionfinal", singleonly=TRUE)
# PrintKdValues(test)

# test <- EquilPars("miR-7-23nt", exp="equilibrium2_nb",
#                   sitelist="resubmissionfinal", singleonly=TRUE)
# PrintKdValues(test)


GetReporterAssaySites <- function() {
  sites <- names(kMPRA)
  print(head(sites))
  sites_unique <- unique(gsub("^(.*)_(.*)_(.*)$", sites, replace="\\1_\\2"))
  print(head(sites_unique))
  print(length(sites_unique))
}

GetReporterAssaySites()




CompareMpraAndFinalSiteList <- function(mirna, sitelist="resubmissionfinal") {
  sitelist <- SubfunctionCall(GetSiteList, sitelist=sitelist)
  sites_mpra <- SubfunctionCall(GetMpraSites)
  sites_missing_mpra <- setdiff(sitelist, sites_mpra)
  # print(sites_missing_mpra)
  # sites_missing_mpra <<- sites_missing_mpra
  if (length(sites_missing_mpra) > 0) {
    message("Sites in final sitelist but absent from MPRA:")  
  }
  for (site in sites_missing_mpra) {
    site_sequence <- GetSiteSeq(mirna, site)
    message(sprintf("\t%s:\t %s", site, length(grep(site_sequence, kMPRA))))
  }
  sites_not_needed_mpra <- setdiff(sites_mpra, sitelist)
  if (length(sites_not_needed_mpra) > 0) {
    message("Sites in MPRA that aren't needed:")  
  }
  for (site in sites_not_needed_mpra) {
    message(sprintf("\t%s", site))
  }
}

for (mirna in kMirnas) {
  message(sprintf("%s ________________________", mirna))
  if (mirna == "miR-7-23nt") {
    mirna <- "miR-7"
  }
  CompareMpraAndFinalSiteList(mirna)
  # CompareMpraAndFinalSiteList(mirna, sitelist="paperfinal")

}

break
# Input from RNA library 1
I_L1 <- Norm(KmersXCountsVector("miR-1", "equil_pilot", "I", 0, 8, 0,
                                       buffer=TRUE))
# Input from RNA library 2
I_L2 <- Norm(KmersXCountsVector("miR-1", "equilibrium", "I", 0, 8, 0,
                                 buffer=TRUE))
# Input from RNA library 3
I_L3.i <- Norm(KmersXCountsVector("let-7a", "equilibrium", "I", 0, 8, 0))
I_L3.ii <- Norm(KmersXCountsVector("miR-1", "kin_pilot", "I_TGT", 0, 8, 0))
I_L3.iii <- Norm(KmersXCountsVector("miR-124", "equilibrium", "I", 0, 8, 0))

# Input from RNA library 3_nb transcription 1
I_L3.nb.i <- Norm(KmersXCountsVector("let-7a-21nt", "equilibrium_nb", "I", 0, 8,
                                    0))
# Input from RNA library 3_nb transcription 2
I_L3.nb.ii <- Norm(KmersXCountsVector("miR-7-23nt", "equilibrium_nb", "I", 0, 8,
                                    0))
# Input from RNA library 4_nb transcrip
I_L4.nb.i <- Norm(KmersXCountsVector("miR-7-23nt", "equilibrium2_nb", "I", 0,
                                        8, 0))
I_L4.nb.ii <- Norm(KmersXCountsVector("miR-7-24nt", "equilibrium3_nb", "I", 0,
                                      8, 0))
I_L4.tp.i   <- Norm(KmersXCountsVector("miR-7-24nt", "equilibrium_tp", "I", 0, 8,
                                     0))



I_all <- cbind(I_L1, I_L2, I_L3.i, I_L3.ii, I_L3.iii, I_L3.nb.i, I_L3.nb.ii,
               I_L4.nb.i, I_L4.nb.ii, I_L4.tp.i)
colnames(I_all) <- c("I_l1", "I_L2", "I_L3.i", "I_L3.ii", "I_L3.iii",
                     "I_L3.nb.i", "I_L3.nb.ii", "I_L4.nb.i", "I_L4.nb.ii",
                     "I_L4.tp.i")

dev.new(xpos=20, ypos=20, height=5, width=5)
corrplot(cor(log(I_all))^2, method="number")

dev.new(xpos=520, ypos=20, height=4, width=4)
plot(I_L3.nb.ii[, 1], I_l4.nb.i[, 1], log="xy", xlim=c(1e-7, 1e-4),
     ylim=c(1e-7, 1e-4))
text(x=1e-4, y=1e-4, cor(log(I_L3.nb.ii[, 1]), log(I_l4.nb.i[, 1])))

dev.new(xpos=920, ypos=20, height=4, width=4)
plot(I_L4.nb.ii[, 1], I_l4.tp.i[, 1], log="xy", xlim=c(1e-7, 1e-4),
     ylim=c(1e-7, 1e-4))
text(x=1e-4, y=1e-4, cor(log(I_L4.nb.ii[, 1]), log(I_l4.tp.i[, 1])))

dev.new(xpos=1320, ypos=20, height=4, width=4)
plot(I_L4.nb.i[, 1], I_l4.nb.ii[, 1], log="xy", xlim=c(1e-7, 1e-4),
     ylim=c(1e-7, 1e-4))
text(x=1e-4, y=1e-4, cor(log(I_L4.nb.i[, 1]), log(I_l4.nb.ii[, 1])))



print(head(I_all))
break


le7.mi155_A <- Norm(KmersXCountsVector("let-7a", "equilibrium", "40", 0, 8, 0))

R_40_sm <- le7.mi155_A/le7.mi155_I

# mi155_I <- Norm(KmersXCountsVector("miR-155", "equilibrium", "I", 0, 8, 0))
# ls6_I <- Norm(KmersXCountsVector("lsy-6", "equilibrium", "I", 0, 8, 0))
# Initial kinetic input. This should be the same as the prior two inputs.






# Namita Input for first miR-7 experiment. Should be the same as her input
# for let-7a.
le7_A_nb <- Norm(KmersXCountsVector("let-7a-21nt", "equilibrium_nb", "40", 0, 8,
                                    0))

R_40_nb <- le7_A_nb/le7_I_nb

R_40_nb_sm <- le7_A_nb/le7.mi155_I

R_40_sm_nb <- le7.mi155_A/le7_I_nb



mi7_I_nb_2   <- Norm(KmersXCountsVector("miR-7-23nt", "equilibrium2_nb", "I", 0,
                                        8, 0))
mi7_I_nb_3 <- Norm(KmersXCountsVector("miR-7-24nt", "equilibrium3_nb", "I", 0,
                                      8, 0))
m7_I_tp   <- Norm(KmersXCountsVector("miR-7-24nt", "equilibrium_tp", "I", 0, 8,
                                     0))




dev.new(xpos=20, ypos=20, height=5, width=5)
plot(mi[, 1], mi7_I_nb_3[, 1], log="xy", xlim=c(1e-7, 1e-4),
     ylim=c(1e-7, 1e-4))




plot(le7.mi155_A[, 1], le7_A_nb[, 1], log="xy", xlim=c(1e-7, 1e-2),
     ylim=c(1e-7, 1e-2))

# dev.new(xpos=520, ypos=20, height=5, width=5)
# plot(le7.mi155_I[, 1], mi124.ls6_I[, 1], log="xy", xlim=c(1e-7, 1e-4),
#      ylim=c(1e-7, 1e-4))

# dev.new(xpos=1020, ypos=20, height=5, width=5)
# plot(le7.mi155_I[, 1], mi1_kin_I[, 1], log="xy", xlim=c(1e-7, 1e-4),
#      ylim=c(1e-7, 1e-4))

# dev.new(xpos=20, ypos=520, height=5, width=5)
# plot(mi7_I_nb_2[, 1], mi7_I_nb_3[, 1], log="xy", xlim=c(1e-7, 1e-4),
#      ylim=c(1e-7, 1e-4))

dev.new(xpos=20, ypos=520, height=5, width=5)
# plot(le7_I_nb[, 1], mi7_I_nb[, 1], log="xy", xlim=c(1e-7, 1e-4),
#      ylim=c(1e-7, 1e-4))
plot(le7.mi155_I[, 1], le7_I_nb[, 1], log="xy", xlim=c(1e-7, 1e-4),
     ylim=c(1e-7, 1e-4))


dev.new(xpos=520, ypos=520, height=5, width=5)
plot(R_40_sm[, 1], R_40_nb[, 1], log="xy", xlim=c(1e-1, 1e3),
     ylim=c(1e-1, 1e3))

dev.new(xpos=1020, ypos=520, height=5, width=5)
plot(R_40_sm[, 1], R_40_nb_sm[, 1], log="xy", xlim=c(1e-2, 1e3),
     ylim=c(1e-2, 1e3))

dev.new(xpos=1520, ypos=520, height=5, width=5)
plot(R_40_nb_sm[, 1], R_40_nb[, 1], log="xy", xlim=c(1e-2, 1e3),
     ylim=c(1e-2, 1e3))



break

# identify(le7_I_nb[, 1], mi7_I_nb[, 1], labels=rownames(let7_I_nb))

R_l7.m7 <- le7_I_nb/mi7_I_nb



R_sort <- R_l7.m7[order(R_l7.m7[, 1], decreasing=TRUE), 1, drop=FALSE]

print(head(R_sort, n=20))

print(tail(R_sort, n=20)^(-1))



break


# dev.new(xpos=20, ypos=520, height=5, width=5)
# plot(mi124_I[, 1], ls6_I[, 1], log="xy", xlim=c(1e-7, 1e-4),
#      ylim=c(1e-7, 1e-4))

# dev.new(xpos=520, ypos=520, height=5, width=5)
# plot(le7_I[, 1], mi7_I_nb[, 1], log="xy", xlim=c(1e-7, 1e-4),
#      ylim=c(1e-7, 1e-4))

# dev.new(xpos=1020, ypos=520, height=5, width=5)
# plot(mi7_I_nb[, 1], mi7_I_nb_2[, 1], log="xy", xlim=c(1e-7, 1e-4),
#      ylim=c(1e-7, 1e-4))




# identify(mi7_I_nb[, 1], mi7_I_nb_2[, 1], labels=rownames(mi1_I))



break


date <- "190816"

for (mirna in kMirnas) {
  if (mirna == "miR-1") buffer <- TRUE
  else                  buffer <- FALSE
  if (mirna == "miR-7-23nt") experiment <- "equilibrium2_nb"
  else                       experiment <- "equilibrium"
  if (mirna %in% c("miR-1", "miR-7-24nt")) combined <- FALSE
  else                                     combined <- TRUE
  PlotSiteKds(mirna, experiment=experiment, sitelist="papersubmission",
              adjusted_height=TRUE, combined=combined, buffer=buffer,
              pdf.plot=sprintf("%s/%s_papersubmission_kds", date, mirna))
}

break
# CheckReporterAssaySites("miR-1", pdf.plot="190815/miR-1_reporter_ambsites")
# PlotReporterAssayKdVsL2fc("miR-1", "twist_reporter_assay_v2",
#                           combined_reps=TRUE, exclude_comp=FALSE,
#                           pdf.plot="190815/miR-1_reporter_assay")

# CheckReporterAssaySites("let-7a", pdf.plot="190815/let-7a_reporter_ambsites")
# PlotReporterAssayKdVsL2fc("let-7a", "twist_reporter_assay_v2",
#                           combined_reps=TRUE, exclude_comp=FALSE,
#                           pdf.plot="190815/let-7a_reporter_assay")

# CheckReporterAssaySites("miR-155",
#                         pdf.plot="190815/miR-155_reporter_ambsites")
# PlotReporterAssayKdVsL2fc("miR-155", "twist_reporter_assay_v2",
#                           combined_reps=TRUE, exclude_comp=FALSE,
#                           pdf.plot="190815/miR-155_reporter_assay")

# CheckReporterAssaySites("lsy-6", pdf.plot="190815/lsy-6_reporter_ambsites")
# PlotReporterAssayKdVsL2fc("lsy-6", "twist_reporter_assay_v2",
#                           combined_reps=TRUE, exclude_comp=FALSE,
#                           pdf.plot="190815/lsy-6_reporter_assay")

# CheckReporterAssaySites("miR-7", pdf.plot="190815/miR-7_reporter_ambsites")
# PlotReporterAssayKdVsL2fc("miR-7", "twist_reporter_assay_v2",
#                           combined_reps=TRUE, exclude_comp=FALSE,
#                           pdf.plot="190815/miR-7_reporter_assay")


break


print(variant_inds_2)
# print(test[variant_inds_2])

print(variant_inds_3)
# print(test[variant_inds_3])

print(variant_inds_4)
# print(test[variant_inds_4])




break

GetReporterVariantSequences("let-7a")


break





PlotSiteKds("miR-155", sitelist="papercutoff", adjusted_height=TRUE,
            showalignment=TRUE, pdf.plot="190812/miR-155_kds_papercutoff")

PlotSiteKds("miR-155", sitelist="papersubmission", adjusted_height=TRUE,
            pdf.plot="190812/miR-155_kds_papersubmission")

PlotSiteKds("miR-155", sitelist="paperfinal", adjusted_height=TRUE,
            pdf.plot="190812/miR-155_kds_paperfinal")


break


test <- ReporterCounts("miR-155", "twist_reporter_assay", "duplex", rep=1)

all_mirna_sites <- rownames(test)

mirna <- "miR-155"

mirna_sites <- grep(sprintf("%s_", mirna), all_mirna_sites, value=TRUE)

mirna_sites <- gsub(sprintf("%s_", mirna), mirna_sites, replace="")

print(mirna_sites)



break



PlotSiteKds("miR-155", "equilibrium", 5, "papercutoff", single=FALSE,
            adjusted_height=TRUE)


PlotSiteKds("miR-155", "equilibrium", 5, "papersubmission", single=FALSE,
            adjusted_height=TRUE)

break

# PlotComparison("let-7a", pdf.plot="190804/let-7a_detection_becker")

# PlotComparison("miR-21", pdf.plot="190804/miR-21_detection_becker")

# MacRaeA1DeltaGCalculation()

# PlotCompetitorOligoSiteSchematic(sitelist="papercutoff2")
# break


# MakeCompetitorPlot("miR-1", "equilibrium", n_constant=5,
#                    pdf.plot="190805/miR-1_sm_comp_enrich")
# MakeCompetitorPlot("miR-124", "equilibrium", n_constant=5,
#                    pdf.plot="190805/miR-124_sm_comp_enrich")
# MakeCompetitorPlot("miR-7-23nt", "equilibrium2_nb", n_constant=5,
#                    pdf.plot="190805/miR-7-23nt_nb_comp_enrich")

# MakeCompetitorPlot("miR-1", "equilibrium_tp", n_constant=5,
#                    pdf.plot="190805/miR-1_tp_comp_enrich")
# MakeCompetitorPlot("miR-124", "equilibrium_2_tp", n_constant=5,
#                    pdf.plot="190805/miR-124_tp_comp_enrich")
# MakeCompetitorPlot("miR-7-24nt", "equilibrium_tp", n_constant=5,
#                    pdf.plot="190805/miR-7-24nt_tp_comp_enrich")




# break
# PlotSiteKds("miR-1", "equilibrium", sitelist="papercutoff", singleonly=FALSE,
#             buffer=TRUE, added.text=TRUE, adjusted_height=TRUE,
#             pdf.plot="190804/miR-1_sm_papercutoff_kds")
# PlotSiteKds("miR-1", "equilibrium", sitelist="papercutoff2", singleonly=FALSE,
#             adjusted_height=TRUE, buffer=TRUE, added.text=TRUE,
#             pdf.plot="190804/miR-1_sm_papercutoff2_kds")
# PlotSiteKds("miR-1", "equilibrium_tp", sitelist="papercutoff", singleonly=FALSE,
#             adjusted_height=TRUE, added.text=TRUE,
#             pdf.plot="190804/miR-1_tp_papercutoff_kds")
# PlotSiteKds("miR-1", "equilibrium_tp", sitelist="papercutoff2", singleonly=FALSE,
#             adjusted_height=TRUE, added.text=TRUE,
#             pdf.plot="190804/miR-1_tp_papercutoff2_kds")

# PlotSiteKds("miR-1", "equilibrium", sitelist="papercutoff3", singleonly=FALSE,
#             adjusted_height=TRUE, added.text=TRUE, buffer=TRUE,
#             pdf.plot="190804/miR-1_sm_papercutoff3_kds")
# PlotSiteKds("miR-1", "equilibrium_tp", sitelist="papercutoff3", singleonly=FALSE,
#             adjusted_height=TRUE, added.text=TRUE,
#             pdf.plot="190804/miR-1_tp_papercutoff3_kds")







# PlotSiteKds("miR-124", "equilibrium", sitelist="papercutoff", singleonly=FALSE,
#             added.text=TRUE, adjusted_height=TRUE,
#             pdf.plot="190804/miR-124_sm_papercutoff_kds")
# PlotSiteKds("miR-124", "equilibrium", sitelist="papercutoff2", singleonly=FALSE,
#             adjusted_height=TRUE, added.text=TRUE,
#             pdf.plot="190804/miR-124_sm_papercutoff2_kds")
# PlotSiteKds("miR-124", "equilibrium_2_tp", sitelist="papercutoff", singleonly=FALSE,
#             adjusted_height=TRUE, added.text=TRUE,
#             pdf.plot="190804/miR-124_tp_papercutoff_kds")
# PlotSiteKds("miR-124", "equilibrium_2_tp", sitelist="papercutoff2", singleonly=FALSE,
#             adjusted_height=TRUE, added.text=TRUE,
#             pdf.plot="190804/miR-124_tp_papercutoff2_kds")

# PlotSiteKds("miR-124", "equilibrium", sitelist="papercutoff3", singleonly=FALSE,
#             adjusted_height=TRUE, added.text=TRUE,
#             pdf.plot="190804/miR-124_sm_papercutoff3_kds")
# PlotSiteKds("miR-124", "equilibrium_2_tp", sitelist="papercutoff3",
#             singleonly=FALSE, adjusted_height=TRUE, added.text=TRUE,
#             pdf.plot="190804/miR-124_tp_papercutoff3_kds")

# # miR-7#########################################################################
# PlotSiteKds("miR-7-23nt", "equilibrium2_nb", sitelist="papercutoff", singleonly=FALSE,
#             added.text=TRUE, adjusted_height=TRUE,
#             pdf.plot="190804/miR-7_nb_papercutoff_kds")
# PlotSiteKds("miR-7-23nt", "equilibrium2_nb", sitelist="papercutoff2", singleonly=FALSE,
#             added.text=TRUE, adjusted_height=TRUE,
#             pdf.plot="190804/miR-7_nb_papercutoff2_kds")
# PlotSiteKds("miR-7-23nt", "equilibrium2_nb", sitelist="papercutoff3", singleonly=FALSE,
#             added.text=TRUE, adjusted_height=TRUE,
#             pdf.plot="190804/miR-7_nb_papercutoff3_kds")

# PlotSiteKds("miR-7-24nt", "equilibrium_tp", sitelist="papercutoff",
#             singleonly=FALSE, adjusted_height=TRUE, added.text=TRUE,
#             pdf.plot="190804/miR-7_tp_papercutoff_kds")
# PlotSiteKds("miR-7-24nt", "equilibrium_tp", sitelist="papercutoff2",
#             singleonly=FALSE, adjusted_height=TRUE, added.text=TRUE,
#             pdf.plot="190804/miR-7_tp_papercutoff2_kds")
# PlotSiteKds("miR-7-24nt", "equilibrium_tp", sitelist="papercutoff3",
#             singleonly=FALSE, adjusted_height=TRUE, added.text=TRUE,
#             pdf.plot="190804/miR-7_tp_papercutoff3_kds")

PlotPairwiseKds("miR-7-23nt", "equilibrium2_nb", 5, "papercutoff3",
                "miR-7-24nt", "equilibrium_tp", 5, "papercutoff3",
                pdf.plot="190804/miR-7_nb_vs_tp_papercutoff3")

break


output_miR21 <- GetMismatchStretchVariantBecker("miR-21")
output_let7 <- GetMismatchStretchVariantBecker("let-7a")

list_Kds_miR21 <- GetSiteSeqKds("miR-21", plotting=FALSE)
list_Kds_let7a <- GetSiteSeqKds("let-7a", plotting=FALSE)

PlotBeckerFigure3D(output_miR21, output_let7, pdf.plot="190727/Becker_3D")
PlotBeckerFigure3D(output_miR21, output_let7, average_multiple=TRUE,
                   pdf.plot="190727/Becker_3D_dif")

PlotBeckerFigure3DAlt(output_miR21, output_let7, list_Kds_miR21, list_Kds_let7a,
                      pdf.plot="190727/Becker_3D_dif2")



break


# PlotPairwiseReporterCounts("miR-1", "twist_reporter_assay", "duplex",
#                            pdf.plot="190803/v1_pairwise_miR-1")
# PlotPairwiseReporterCounts("let-7a", "twist_reporter_assay", "duplex", 
#                            pdf.plot="190803/v1_pairwise_let-7a")
# PlotPairwiseReporterCounts("miR-155", "twist_reporter_assay", "duplex",
#                            pdf.plot="190803/v1_pairwise_miR-155")
# PlotPairwiseReporterCounts("miR-124", "twist_reporter_assay", "duplex",
#                            pdf.plot="190803/v1_pairwise_miR-124")
# PlotPairwiseReporterCounts("lsy-6", "twist_reporter_assay", "duplex",
#                            pdf.plot="190803/v1_pairwise_lsy-6")
# PlotPairwiseReporterCounts("miR-7", "twist_reporter_assay", "duplex",
#                            pdf.plot="190803/v1_pairwise_miR-7")

# PlotPairwiseReporterCounts("miR-1", "twist_reporter_assay_v2", "duplex", 
#                            pdf.plot="190803/v2_pairwise_miR-1")
# PlotPairwiseReporterCounts("let-7a", "twist_reporter_assay_v2", "duplex",
#                            pdf.plot="190803/v2_pairwise_let-7a")
# PlotPairwiseReporterCounts("miR-155", "twist_reporter_assay_v2", "duplex",
#                            pdf.plot="190803/v2_pairwise_miR-155")
# PlotPairwiseReporterCounts("miR-124", "twist_reporter_assay_v2", "duplex",
#                            pdf.plot="190803/v2_pairwise_miR-124")
# PlotPairwiseReporterCounts("lsy-6", "twist_reporter_assay_v2", "duplex",
#                            pdf.plot="190803/v2_pairwise_lsy-6")
# PlotPairwiseReporterCounts("miR-7", "twist_reporter_assay_v2", "duplex",
#                            pdf.plot="190803/v2_pairwise_miR-7")

# # PlotReporterAssayKdVsL2fc("miR-1", "twist_reporter_assay",
# #                           pdf.plot="190803/miR-1_v1_rep2")
# # PlotReporterAssayKdVsL2fc("let-7a", "twist_reporter_assay",
# #                           pdf.plot="190803/let-7a_v1_rep2")
# # PlotReporterAssayKdVsL2fc("miR-155", "twist_reporter_assay",
# #                           pdf.plot="190803/miR-155_v1_rep2")
# # PlotReporterAssayKdVsL2fc("miR-124", "twist_reporter_assay",
# #                           pdf.plot="190803/miR-124_v1_rep2")
# # PlotReporterAssayKdVsL2fc("lsy-6", "twist_reporter_assay",
# #                           pdf.plot="190803/lsy-6_v1_rep2")
# # PlotReporterAssayKdVsL2fc("miR-7", "twist_reporter_assay",
# #                           pdf.plot="190803/miR-7_v1_rep2")

# # PlotReporterAssayKdVsL2fc("miR-1", "twist_reporter_assay_v2",
# #                           pdf.plot="190803/miR-1_v2_rep2")
# # PlotReporterAssayKdVsL2fc("let-7a", "twist_reporter_assay_v2",
# #                           pdf.plot="190803/let-7a_v2_rep2")
# # PlotReporterAssayKdVsL2fc("miR-155", "twist_reporter_assay_v2",
# #                           pdf.plot="190803/miR-155_v2_rep2")
# # PlotReporterAssayKdVsL2fc("miR-124", "twist_reporter_assay_v2",
# #                           pdf.plot="190803/miR-124_v2_rep2")
# # PlotReporterAssayKdVsL2fc("lsy-6", "twist_reporter_assay_v2",
# #                           pdf.plot="190803/lsy-6_v2_rep2")
# # PlotReporterAssayKdVsL2fc("miR-7", "twist_reporter_assay_v2",
# #                           pdf.plot="190803/miR-7_v2_rep2")


# PlotReporterAssayKdVsL2fc("miR-1", "twist_reporter_assay_v2", rep=1,
#                           pdf.plot="190803/miR-1_v2_rep1")
# PlotReporterAssayKdVsL2fc("let-7a", "twist_reporter_assay_v2", rep=1,
#                           pdf.plot="190803/let-7a_v2_rep1")
# PlotReporterAssayKdVsL2fc("miR-155", "twist_reporter_assay_v2", rep=1,
#                           pdf.plot="190803/miR-155_v2_rep1")
# PlotReporterAssayKdVsL2fc("miR-124", "twist_reporter_assay_v2", rep=1,
#                           pdf.plot="190803/miR-124_v2_rep1")
# PlotReporterAssayKdVsL2fc("lsy-6", "twist_reporter_assay_v2", rep=1,
#                           pdf.plot="190803/lsy-6_v2_rep1")
# PlotReporterAssayKdVsL2fc("miR-7", "twist_reporter_assay_v2", rep=1,
#                           pdf.plot="190803/miR-7_v2_rep1")



# PlotReporterAssayKdVsL2fc("miR-1", "twist_reporter_assay", rep=1,
#                           pdf.plot="190803/miR-1_v1_rep1")
# PlotReporterAssayKdVsL2fc("let-7a", "twist_reporter_assay", rep=1,
#                           pdf.plot="190803/let-7a_v1_rep1")
# PlotReporterAssayKdVsL2fc("miR-155", "twist_reporter_assay", rep=1,
#                           pdf.plot="190803/miR-155_v1_rep1")
# PlotReporterAssayKdVsL2fc("miR-124", "twist_reporter_assay", rep=1,
#                           pdf.plot="190803/miR-124_v1_rep1")
# PlotReporterAssayKdVsL2fc("lsy-6", "twist_reporter_assay", rep=1,
#                           pdf.plot="190803/lsy-6_v1_rep1")
# PlotReporterAssayKdVsL2fc("miR-7", "twist_reporter_assay", rep=1,
#                           pdf.plot="190803/miR-7_v1_rep1")


# Plot of v2 where the reps are combined.
PlotReporterAssayKdVsL2fc("miR-1", "twist_reporter_assay_v2", combined_reps=TRUE,
                          pdf.plot="190803/miR-1_v2")
PlotReporterAssayKdVsL2fc("let-7a", "twist_reporter_assay_v2", combined_reps=TRUE,
                          pdf.plot="190803/let-7a_v2")
PlotReporterAssayKdVsL2fc("miR-155", "twist_reporter_assay_v2", combined_reps=TRUE,
                          pdf.plot="190803/miR-155_v2")
PlotReporterAssayKdVsL2fc("miR-124", "twist_reporter_assay_v2", combined_reps=TRUE,
                          pdf.plot="190803/miR-124_v2")
PlotReporterAssayKdVsL2fc("lsy-6", "twist_reporter_assay_v2", combined_reps=TRUE,
                          pdf.plot="190803/lsy-6_v2")
PlotReporterAssayKdVsL2fc("miR-7", "twist_reporter_assay_v2", combined_reps=TRUE,
                          pdf.plot="190803/miR-7_v2")




# PlotAllPairwiseReporterCounts("miR-1", "twist_reporter_assay", "duplex", "8mer",
#                               pdf.plot="190803/miR-1_v1_8mer_cor")
# PlotAllPairwiseReporterCounts("let-7a", "twist_reporter_assay", "duplex", "8mer",
#                               pdf.plot="190803/let-7a_v1_8mer_cor")
# PlotAllPairwiseReporterCounts("miR-155", "twist_reporter_assay", "duplex",
#                               "8mer", pdf.plot="190803/miR-155_v1_8mer_cor")
# PlotAllPairwiseReporterCounts("miR-124", "twist_reporter_assay", "duplex",
#                               "8mer", pdf.plot="190803/miR-124_v1_8mer_cor")
# PlotAllPairwiseReporterCounts("lsy-6", "twist_reporter_assay", "duplex", "8mer",
#                               pdf.plot="190803/lsy-6_v1_8mer_cor")
# PlotAllPairwiseReporterCounts("miR-7", "twist_reporter_assay", "duplex", "8mer",
#                               pdf.plot="190803/miR-7_v1_8mer_cor")

# PlotAllPairwiseReporterCounts("miR-1", "twist_reporter_assay_v2", "duplex",
#                               "8mer", pdf.plot="190803/miR-1_v2_8mer_cor")
# PlotAllPairwiseReporterCounts("let-7a", "twist_reporter_assay_v2", "duplex",
#                               "8mer", pdf.plot="190803/let-7a_v2_8mer_cor")
# PlotAllPairwiseReporterCounts("miR-155", "twist_reporter_assay_v2", "duplex",
#                               "8mer", pdf.plot="190803/miR-155_v2_8mer_cor")
# PlotAllPairwiseReporterCounts("miR-124", "twist_reporter_assay_v2", "duplex",
#                               "8mer", pdf.plot="190803/miR-124_v2_8mer_cor")
# PlotAllPairwiseReporterCounts("lsy-6", "twist_reporter_assay_v2", "duplex",
#                               "8mer", pdf.plot="190803/lsy-6_v2_8mer_cor")
# PlotAllPairwiseReporterCounts("miR-7", "twist_reporter_assay_v2", "duplex",
#                               "8mer", pdf.plot="190803/miR-7_v2_8mer_cor")





break



GetSiteList <- function(mirna, sitelist) {
  path <- sprintf("AssignSiteTypes/site_categories/%s/sites.%s_%s.txt", 
                     sitelist, mirna, sitelist)
  unlist(read.table(path, stringsAsFactors=FALSE))
}


temp <- GetSiteList("miR-1", "papercutoff")
temp2 <- GetSiteList("miR-1", "papercutoff2")

print(temp)
print(intersect(temp, temp2))
print(setdiff(temp, temp2))

break


# PlotSiteKds_("miR-124", "equilibrium_2_tp", 5, "paperfinal", adjusted_height=TRUE,
#              singleonly=FALSE,
#             pdf.plot="190722/miR-124_kds_paperfinal")



# PlotSiteKds_("miR-124", "equilibrium_2_tp", 5, "papercutoff", adjusted_height=TRUE,
#              singleonly=FALSE,
#             pdf.plot="190722/miR-124_kds_new")

# break



GetPerfectTargetMismatchStretch <- function(mirna, mm_start, mm_stop) {
  # Get the starting string containing the perfectly complementary target site.
  perfect_str <- SubfunctionCall(GetSiteSeq, "21mer-m1.21")
  n <- 21
  # Convert this string to a vector of nucleotide letters.
  perfect_vec <- unlist(strsplit(perfect_str, split=""))
  # Loop over the list of contiguous postions to be changed to their complment
  # (i.e., changed back to the same nucleotide identity as that of the miRNA).
  for (pos in (21 - mm_start:mm_stop + 1)) {
    # kComplements in Lists.R
    perfect_vec[pos] <- kComplements[perfect_vec[pos]]
  }
  # Return the string formed from collapsing the vector of nucleotide strings.
  paste0(perfect_vec, collapse="")
}


GetMismatchStretchVariantBecker <- function(mirna) {
  # Load the data table
  data <- LoadBeckerEtAlData(mirna)
  # Make the vector of variant names.
  variants <- rownames(data)
  # Make the output table.
  output_matrix <- matrix(NaN, nrow=21, ncol=21)
  # Loop over the start and stop positions of the mismatch stretch.
  for (mm_start in 1:21) {
    for (mm_stop in mm_start:21) {
      # Make the string corresponding to the perfect match with the appropriate
      # stretch of mismatches.
      regex_str <- SubfunctionCall(GetPerfectTargetMismatchStretch)
      # Subset the data.
      data_subset <- data[sprintf("AAAAA%sAAAAA", regex_str), 1]

      output_matrix[21 - mm_start + 1, mm_stop] <- data_subset
    }
  }
  output_matrix
}


# output_miR21 <- GetMismatchStretchVariantBecker("miR-21")
# output_let7 <- GetMismatchStretchVariantBecker("let-7a")

PlotBeckerFigure3C <- function(output_miR21, output_let7, xpos=20, ypos=20,
                               height=5, width=5, pdf.plot=FALSE) {
  # Define the bounds for the 
  output_miR21 <- t(output_miR21[, 21:1])[, 21:1]

  output_full <- matrix(NaN, nrow=23, ncol=23)

  for (row_i in 1:21) {
    for (col_i in 1:(21 - row_i + 1)) {
      output_full[row_i, col_i] <- output_miR21[row_i, col_i]
    }
  }

  for (row_i in 1:21) {
    for (col_i in (21 - row_i + 1):21) {
      output_full[row_i + 2, col_i + 2] <- output_let7[row_i, col_i]
    }
  }
  xlefts <- rep(seq(0, 22), each=23) 
  xright <- xlefts + 1
  ybottom <- rep(seq(22, 0), 23)
  ytop <- ybottom + 1
  output <- output_full

  output[which(output >= 5000)] <- 5000
  output <- log10(output)
  output <- output - min(output, na.rm=TRUE)
  output <- output/max(output, na.rm=TRUE)
  print(output[10:21, 10:21])
  col.inds <- round(output*100)
  print(col.inds)

  # Parameters specifying the start and end of the rainbow colorspace usee:
  c_s <- 0.9 # color_start
  c_e <- 0.70  # color_end
  color.dist = plasma(100)
  col.inds <- sapply(col.inds, function(col.ind) {
    min(max(1, col.ind), 100)
    })

  SubfunctionCall(GenericFigureSaveFile)
  xmin <- 0
  xmax <- 23
  ymin <- xmin
  ymax <- xmax
  BlankPlot()

  rect(xlefts, ybottom, xright, ytop, col=color.dist[col.inds], border=NA)
  # rect(xlefts, ybottom, xright, ytop, col=color.dist[col.inds])
  xy <- GetPlotFractionalCoords(0, 1)
  text(0:20 + 0.5, xy[2], labels=1:21, xpd=NA, adj=c(0.5, 0), cex=0.8)

  

  xy <- GetPlotFractionalCoords(0, 1)
  text(xy[1], 0:20 + 2.5, labels=1:21, xpd=NA, adj=c(1, 0.5), cex=0.8)
  xy <- GetPlotFractionalCoords(1, 1)
  text(xy[1], 0:20 + 0.5, labels=1:21, xpd=NA, adj=c(0, 0.5), cex=0.8)
  xy <- GetPlotFractionalCoords(0, 0)
  text(0:20 + 2.5, xy[2], labels=1:21, xpd=NA, adj=c(0.5, 1), cex=0.8)

  # Left side label.
  xy <- GetPlotFractionalCoords(-0.075, 0.55)
  text(xy[1], xy[2], labels="Ending complement match", srt=90, xpd=NA,
       adj=c(0.5, 0))
  # Right site label
  xy <- GetPlotFractionalCoords(1.075, 0.45)
  text(xy[1], xy[2], labels="Beginning complement match", srt=270, xpd=NA,
       adj=c(0.5, 0))
  # Top side label
  xy <- GetPlotFractionalCoords(0.45, 1.075)
  text(xy[1], xy[2], labels="Beginning complement match", xpd=NA,
       adj=c(0.5, 0))
  # Right site label
  xy <- GetPlotFractionalCoords(0.55, -0.075)
  text(xy[1], xy[2], labels="Ending complement match", xpd=NA, adj=c(0.5, 0))

  xy <- GetPlotFractionalCoords(0, -0.07)
  text(xy[1], xy[2], label="let-7a", adj=c(0, 0), xpd=NA)

  xy <- GetPlotFractionalCoords(0.01, 0.025)
  text(xy[1], xy[2], label="miR-21", adj=c(1, 0), xpd=NA)
  segments(-0.5, -0.5, 23, 23, xpd=NA)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}



# list_Kds_let7a <- ComparePositions9Through13("let-7a", pdf.plot="190724/let-7a_9-13_ecdf")
# list_Kds_miR21 <- ComparePositions9Through13("miR-21", pdf.plot="190724/miR-21_9-13_ecdf")




PlotBeckerFigure3C(output_miR21, output_let7, pdf.plot="190727/Becker_3C")
PlotBeckerFigure3D(output_miR21, output_let7, pdf.plot="190727/Becker_3D")
PlotBeckerFigure3D(output_miR21, output_let7, average_multiple=TRUE,
                   pdf.plot="190727/Becker_3D_dif")

PlotBeckerFigure3DAlt(output_miR21, output_let7, list_Kds_miR21, list_Kds_let7a,
                      pdf.plot="190727/Becker_3D_dif2")




# GetSiteSeqKds("let-7a", pdf.plot="190724/let-7a_ecdf")

# PlotComparison("let-7a", pdf.plot="190724/let-7a_vs_let-7a_scatter")


# GetSiteSeqKds("miR-21", pdf.plot="190724/miR-21_ecdf")

# PlotComparison("miR-21", pdf.plot="190724/miR_21_vs_miR-1_scatter")


# list_Kds_let7a_alt <- ComparePositions9Through13(
#   "let-7a", alt=TRUE
# )

# break

# list_Kds_let7a_alt <- ComparePositions9Through13(
#   "let-7a", alt=TRUE, pdf.plot="190724/let-7a_alt_9-13_ecdf"
# )
# list_Kds_miR21_alt <- ComparePositions9Through13(
#   "miR-21", alt=TRUE, pdf.plot="190724/miR-21_alt_9-13_ecdf"
# )



break
t11mer_sequence <- RevComplement(substr(kMirnaSeqs["let-7a"], 1, 11))


seqs <- sapply(rownames(list_Kds_let7a[["11mer-m1.11"]]), function(variant_sequence) {
  unlist(strsplit(variant_sequence, split=t11mer_sequence))[1]
})

target_sequence <- RevComplement(substr(kMirnaSeqs["let-7a"], 12, 21))
target_sequence_alt <- GetSiteSeq("let-7a", "10mer-m12.21")

seq_8mer <- GetSiteSeq("let-7a", "8mer")


print(target_sequence)
print(target_sequence_alt)



for (i in 1:length(seqs)) {
  message(sprintf("====================================%s===========================", i))
  seq <- seqs[i]
  print(seq)
  substring_test <- LCSubString(target_sequence, seq)
  length_sub <- substring_test[[1]]
  starts_sub_1 <- substring_test[[2]][[1]] + 1
  starts_sub_2 <- substring_test[[2]][[2]] + 1
  starts <- cbind(starts_sub_1, starts_sub_2)
  print(substring_test)
  apply(starts, 1, function(row) {
    start_1 <- row[1]
    start_2 <- row[2]
    message(target_sequence)
    message(paste0(rep(" ", start_1 - 1), substr(target_sequence, start_1, start_1 + length_sub - 1), collapse=""))  
    message(seq)
    message(paste0(rep(" ", start_2 - 1), substr(seq, start_2, start_2 + length_sub - 1), collapse=""))  

  })
  break
  seed_test <- LCSubString(seq_8mer, seq)
  print(substring_test)
  print(seed_test)
  print(list_Kds_let7a[["11mer-m1.11"]][i,1])

}

break





# GetSiteSeqKds("miR-21")


break


I_check <- KmersXCountsVector("miR-124", "equilibrium_tp", "I", 5, 6, 0)

I_check_2.1 <- KmersXCountsVector("miR-124", "equilibrium_2_tp", "I,1", 5, 6, 0)
I_check_2.2 <- KmersXCountsVector("miR-124", "equilibrium_2_tp", "I,2", 5, 6, 0)

A40_check <- KmersXCountsVector("miR-124", "equilibrium_tp", "40", 5, 6, 0)

A40_check_2.1 <- KmersXCountsVector("miR-124", "equilibrium_2_tp", "40,1", 5, 6, 0)
A40_check_2.2 <- KmersXCountsVector("miR-124", "equilibrium_2_tp", "40,2", 5, 6, 0)

par <- kPlotParameters
dev.new(xpos=20, ypos=20, height=4, width=4)
plot(I_check[, 1], I_check_2.1[, 1], log='xy')
dev.copy2pdf(file="figures/190711/I_vs_I_2.1.pdf")

dev.new(xpos=420, ypos=20, height=4, width=4)
plot(I_check[, 1], I_check_2.2[, 1], log='xy')
dev.copy2pdf(file="figures/190711/I_vs_I_2.2.pdf")



dev.new(xpos=820, ypos=20, height=4, width=4)
plot(A40_check[, 1], A40_check_2.1[, 1], log='xy')
dev.copy2pdf(file="figures/190711/A40_vs_A40_2.1.pdf")

dev.new(xpos=1220, ypos=20, height=4, width=4)
plot(A40_check[, 1], A40_check_2.2[, 1], log='xy')
dev.copy2pdf(file="figures/190711/A40_vs_A40_2.2.pdf")



R_check <- Norm(A40_check)/Norm(I_check)
R_check_2.1 <- Norm(A40_check_2.1)/Norm(I_check_2.1)
R_check_2.2 <- Norm(A40_check_2.2)/Norm(I_check_2.2)

dev.new(xpos=20, ypos=420, height=4, width=4)
plot(R_check[, 1], R_check_2.1[, 1], log='xy')
dev.copy2pdf(file="figures/190711/R_vs_R_2.1.pdf")

dev.new(xpos=420, ypos=420, height=4, width=4)
plot(R_check[, 1], R_check_2.2[, 1], log='xy')
dev.copy2pdf(file="figures/190711/R_vs_R_2.2.pdf")


Cross_I_check_1 <- Norm(I_check)/Norm(I_check_2.1)
Cross_A40_check_1 <- Norm(A40_check)/Norm(A40_check_2.1)

Cross_I_check_2 <- Norm(I_check)/Norm(I_check_2.2)
Cross_A40_check_2 <- Norm(A40_check)/Norm(A40_check_2.2)

dev.new(xpos=820, ypos=420, height=6, width=6)
plot(Cross_I_check_1[, 1], Cross_A40_check_1[, 1], log='xy')

identify(Cross_I_check_1[, 1], Cross_A40_check_1[, 1],
         labels=rownames(Cross_I_check_1))
dev.copy2pdf(file="figures/190711/Cross_I_40_check.pdf")

# dev.new(xpos=1220, ypos=420, height=4, width=4)
# plot(Cross_I_check_2[, 1], Cross_A40_check_2[, 1], log='xy')

# Cross_R_check_1 <- R_check/R_check_2.1
# Cross_R_check_2 <- R_check/R_check_2.2

# dev.new(xpos=20, ypos=820, height=4, width=4)
# plot(Cross_I_check_1[, 1], Cross_R_check_1[, 1], log='xy')

# dev.new(xpos=420, ypos=820, height=4, width=4)
# plot(Cross_I_check_2[, 1], Cross_R_check_2[, 1], log='xy')





break

PlotPairwiseKds("miR-124", "equilibrium", 5, "paperfinal", "miR-124",
                "equilibrium_2_tp", 5, "paperfinal", xpos=20, ypos=620,
                pdf.plot="190708/miR-124_sm_vs_miR-124_tp")

MakeCompetitorPlot("miR-124", experiment="equilibrium", pdf.plot="190707/miR-124_sm")



break

PlotSiteEnrichments_("miR-7-23nt", experiment="equilibrium2_nb", n_constant=5,
                    sitelist="paperfinal", singleonly=FALSE, combined=FALSE,
                    pdf.plot="190708/miR-7-23nt_2nb_R")

PlotSiteEnrichments_("miR-7-24nt", experiment="equilibrium2_nb", n_constant=5,
                    sitelist="paperfinal", singleonly=FALSE, combined=TRUE,
                    pdf.plot="190708/miR-7-24nt_2nb_R")

PlotSiteEnrichments_("miR-7-25nt", experiment="equilibrium2_nb", n_constant=5,
                    sitelist="paperfinal", singleonly=FALSE, combined=TRUE,
                    pdf.plot="190708/miR-7-25nt_2nb_R")

PlotSiteEnrichments_("miR-7-24nt", experiment="equilibrium_tp", n_constant=5,
                    sitelist="paperfinal", singleonly=FALSE, xpos=620, ypos=20,
                    pdf.plot="190708/miR-7-24nt_tp_R")

PlotSiteEnrichments_("miR-7-24nt", experiment="equilibrium_tp", n_constant=5,
                    sitelist="paperfinal", singleonly=FALSE, minkd=1e-4,
                    xpos=620, ypos=20,
                    pdf.plot="190708/miR-7-24nt_tp_R_minkd1e-4")

PlotSiteEnrichments_("miR-7-24nt", experiment="equilibrium_tp", n_constant=5,
                    sitelist="paperfinal", singleonly=FALSE, minkd=1e-3,
                    xpos=620, ypos=20,
                    pdf.plot="190708/miR-7-24nt_tp_R_minkd1e-3")

PlotPairwiseKds("miR-7-23nt", "equilibrium2_nb", 5, "paperfinal",
                "miR-7-24nt", "equilibrium_tp", 5, "paperfinal",
                combined1=FALSE,
                pdf.plot="190708/miR-7-23nt_2nb_vs_miR-7-24nt_tp")

PlotPairwiseKds("miR-7-23nt", "equilibrium2_nb", 5, "paperfinal",
                "miR-7-24nt", "equilibrium_tp", 5, "paperfinal",
                combined1=FALSE, minkd2=1e-4,
                pdf.plot="190708/miR-7-23nt_2nb_vs_miR-7-24nt_tp_minkd1e-4")

PlotPairwiseKds("miR-7-23nt", "equilibrium2_nb", 5, "paperfinal",
                "miR-7-24nt", "equilibrium_tp", 5, "paperfinal",
                combined1=FALSE, minkd2=1e-3,
                pdf.plot="190708/miR-7-23nt_2nb_vs_miR-7-24nt_tp_minkd1e-3")


break

# PlotSiteEnrichments("miR-1", experiment="equilibrium", n_constant=5, combined=FALSE,
#                     buffer=TRUE, xpos=1220, ypos=20)

PlotPairwiseKds("miR-7-23nt", "equilibrium2_nb", 5, "paperfinal",
                "miR-7-24nt", "equilibrium_tp", 5, "paperfinal",
                combined1=FALSE, xpos=20, ypos=620)

# MakeCompetitorPlot("miR-124", experiment="equilibrium_tp", pdf.plot="190707/miR-124_tp")
# MakeCompetitorPlot("miR-124", experiment="equilibrium_2_tp", pdf.plot="190707/miR-124_2tp")


# MakeCompetitorPlot("miR-7-23nt", experiment="equilibrium2_nb")

# MakeCompetitorPlot("miR-7-24nt", experiment="equilibrium_tp")

# MakeCompetitorPlot("miR-7-23nt", experiment="equilibrium2_nb", pdf.plot="190707/miR-7-23nt_2nb")
# MakeCompetitorPlot("miR-7-24nt", experiment="equilibrium_tp", pdf.plot="190707/miR-7-24nt_tp")


kds_old <- EquilPars("miR-7-23nt", "equilibrium2_nb", singleonly=FALSE, combined=FALSE)

sXc_new <- SitesXCounts("miR-7-24nt", "equilibrium_tp", single=FALSE)

kds_new <- EquilPars("miR-7-24nt", "equilibrium_tp", singleonly=FALSE)


break


n_ex <- 10


temp_data <- SitesXCounts("miR-7-24nt", "equilibrium_tp", n_constant=5)


temp_data_tot <- cbind(rowSums(temp_data[, 1:2]), rowSums(temp_data[, 3:6])) # I and 40
temp_data_tot <- cbind(temp_data_tot, rowSums(temp_data[, 7:8])) # 12.6
temp_data_tot <- cbind(temp_data_tot, rowSums(temp_data[, 9:12])) # 4
temp_data_tot <- cbind(temp_data_tot, rowSums(temp_data[, 13:14])) # 1.26
temp_data_tot <- cbind(temp_data_tot, rowSums(temp_data[, 15:16])) # 0.4
temp_data_tot <- cbind(temp_data_tot, rowSums(temp_data[, 17:18])) # 0

colnames(temp_data_tot) <- c("I", "40", "12.6", "4", "1.26", "0.4", "0")

norm_data <- t(t(temp_data_tot)/colSums(temp_data_tot))
R_data <- norm_data/norm_data[, 1]

sites <- c(rownames(temp_data), "bg", "AGO")


print(R_data)

plot(colnames(R_data)[-1], R_data[1, -1], col=kSiteColors[sites[1]],
     log='xy', xlim=c(1, 100), ylim=c(0.1, 100), pch=20, type="l")
for (row in 2:nrow(temp_data)) {
  lines(colnames(R_data)[-1], R_data[row, -1], pch=20, col=kSiteColors[row])
}




kds_old <- EquilPars("miR-7-23nt", "equilibrium2_nb", combined=FALSE)

kds_new <- EquilPars("miR-7-24nt", "equilibrium_tp")


dev.new(xpos=20, ypos=20, height=5, width=5)
plot(kds_old[, 2], kds_new[, 2], col=kSiteColors[sites], log='xy',
     xlim=c(1e-6, 10), ylim=c(1e-6, 10))

break



for (n_ex in 1:41) {

  test_I <- KmersXCountsVector("miR-124", condition="I_combined",
                               exp="equilibrium_2_tp", n_constant=5, kmer_len=10,
                               n_ex=n_ex)

  test_A <- KmersXCountsVector("miR-124", condition="40",
                               exp="equilibrium_2_tp", n_constant=5, kmer_len=10,
                               n_ex=n_ex)


  norm_A <- as.numeric(test_A[-1, 1])/sum(as.numeric(test_A[-1, 1]))
  norm_I <- as.numeric(test_I[-1, 1])/sum(as.numeric(test_I[-1, 1]))

  names(norm_I) <- rownames(test_I)[2:length(rownames(test_I))]
  names(norm_A) <- names(norm_I)


  print(norm_A["GCATTCACCG"]/norm_I["GCATTCACCG"])

}

break





v1_rep1 <- read.table("ReporterScreen/twist_reporter_assay_rep1_counts.txt")
v1_rep2 <- read.table("ReporterScreen/twist_reporter_assay_rep2_counts.txt")

v2_rep1 <- read.table("ReporterScreen/twist_reporter_assay_v2_rep1_counts.txt")
v2_rep2 <- read.table("ReporterScreen/twist_reporter_assay_v2_rep2_counts.txt")


cols <- rep(rgb(0, 0, 0, alpha=0.5), nrow(v1_rep1))


ind_site <- grep("miR-155_11mer-m13.23_", rownames(v1_rep1))
cols[ind_site] <- "green"

plot(Norm(v2_rep2$miR.155)/Norm(v1_rep2$miR.1),
     Norm(v2_rep1$miR.155)/Norm(v1_rep1$miR.1), log='xy',
     col=rgb(0, 0, 0, alpha=0.2))

points((Norm(v2_rep2$miR.155)/Norm(v1_rep2$miR.1))[ind_site],
     (Norm(v2_rep1$miR.155)/Norm(v1_rep1$miR.1))[ind_site], log='xy',
     col="green", pch=20)


break
# MakeCompetitorPlot("miR-124", exp="equilibrium_2_tp")
# MakeCompetitorPlot("miR-124", exp="equilibrium_tp")

# MakeCompetitorPlot("miR-1", exp="equilibrium_tp")
# MakeCompetitorPlot("miR-1", exp="equilibrium_met_tp")

kds_1 <- EquilPars("miR-124", "equilibrium_2_tp", tp2rep=1)
kds_2 <- EquilPars("miR-124", "equilibrium_2_tp", tp2rep=2)
kds_all <- EquilPars("miR-124", "equilibrium_2_tp")
kds_orig <- EquilPars("miR-124", "equilibrium")

site_names <- gsub("_Kd$", rownames(kds_1), replace="")
site_names <- gsub("_miR-124$", site_names, replace="")

cols <- kSiteColors[site_names]


dev.new(xpos=20, ypos=20, height=5, width=5)
plot(kds_all[, 2], kds_1[, 2], log='xy', col=cols, xlim=c(1e-3, 10), ylim=c(1e-3, 10))
abline(0, 1, lty=2, col=rgb(0, 0, 0, alpha=.5))
dev.new(xpos=520, ypos=20, height=5, width=5)
plot(kds_all[, 2], kds_2[, 2], log='xy', col=cols, xlim=c(1e-3, 10), ylim=c(1e-3, 10))
abline(0, 1, lty=2, col=rgb(0, 0, 0, alpha=.5))
dev.new(xpos=1020, ypos=20, height=5, width=5)
plot(kds_1[, 2], kds_2[, 2], log='xy', col=cols, xlim=c(1e-3, 10), ylim=c(1e-3, 10))
abline(0, 1, lty=2, col=rgb(0, 0, 0, alpha=.5))
dev.new(xpos=20, ypos=520, height=5, width=5)
plot(kds_orig[, 2], kds_all[, 2], log='xy', col=cols, xlim=c(1e-3, 10), ylim=c(1e-3, 10))
abline(0, 1, lty=2, col=rgb(0, 0, 0, alpha=.5))




break

I_m1_1_124sites <- SitesXCountsVector("miR-1", "I", experiment="equilibrium_tp", alt_mirna="miR-124")
I_m1_2_124sites <- SitesXCountsVector("miR-1", "I", experiment="equilibrium_met_tp", alt_mirna="miR-124")

I_m122_124sites <- SitesXCountsVector("miR-122", "I", experiment="equilibrium_tp", alt_mirna="miR-124")

I_m124_1 <- SitesXCountsVector("miR-124", "I", experiment="equilibrium_tp")

I_m124_2.1 <- SitesXCountsVector("miR-124", "I,1", experiment="equilibrium_2_tp")
I_m124_2.2 <- SitesXCountsVector("miR-124", "I,2", experiment="equilibrium_2_tp")

sXc_tp <- SitesXCounts("miR-124", experiment="equilibrium_2_tp", n_constant=5,
                       sitelist="paperfinal")

I_combined <- I_m1_1_124sites + I_m1_2_124sites + I_m122_124sites + I_m124_1 + I_m124_2.1 + I_m124_2.2

plot(sXc_tp[, 3],
     as.numeric(I_combined[, 1]),
     log='xy')
break
print(I_m1_1_124sites)
lims_all <- c(1e-6, 1)
dev.new(xpos=20, ypos=20, height=4, width=4)
plot(Norm(I_m124_2.1[, 1]), Norm(I_m124_2.2[, 1]), log='xy',
     col=kSiteColors[rownames(I_m124_1)], xlim=lims_all, ylim=lims_all)
dev.new(xpos=420, ypos=20, height=4, width=4)
plot(Norm(I_m124_2.1[, 1]), Norm(I_m124_1[, 1]), log='xy',
     col=kSiteColors[rownames(I_m124_1)], xlim=lims_all, ylim=lims_all)
dev.new(xpos=820, ypos=20, height=4, width=4)
plot(Norm(I_m124_2.1[, 1]), Norm(I_m1_1_124sites[, 1]), log='xy',
     col=kSiteColors[rownames(I_m124_1)], xlim=lims_all, ylim=lims_all)
dev.new(xpos=1220, ypos=20, height=4, width=4)
plot(Norm(I_m124_2.1[, 1]), Norm(I_m1_2_124sites[, 1]), log='xy',
     col=kSiteColors[rownames(I_m124_1)], xlim=lims_all, ylim=lims_all)
dev.new(xpos=1620, ypos=20, height=4, width=4)
plot(Norm(I_m124_2.1[, 1]), Norm(I_m122_124sites[, 1]), log='xy',
     col=kSiteColors[rownames(I_m124_1)], xlim=lims_all, ylim=lims_all)


break


I_m122 <- KmersXCountsVector("miR-122", "equilibrium_tp", "I", 5, 6, 0)

I_m1_1 <- KmersXCountsVector("miR-1", "equilibrium_tp", "I", 5, 6, 0)
I_m1_2 <- KmersXCountsVector("miR-1", "equilibrium_met_tp", "I", 5, 6, 0)

I_m124_1 <- KmersXCountsVector("miR-124", "equilibrium_tp", "I", 5, 6, 0)
I_m124_2.1 <- KmersXCountsVector("miR-124", "equilibrium_2_tp", "I,1", 5, 6, 0)
I_m124_2.2 <- KmersXCountsVector("miR-124", "equilibrium_2_tp", "I,2", 5, 6, 0)

A40_m124_1 <- KmersXCountsVector("miR-124", "equilibrium_tp", "40", 5, 6, 0)
A40_m124_2.1 <- KmersXCountsVector("miR-124", "equilibrium_2_tp", "40,1", 5, 6, 0)
A40_m124_2.2 <- KmersXCountsVector("miR-124", "equilibrium_2_tp", "40,2", 5, 6, 0)

R_m124_1 <- Norm(as.numeric(A40_m124_1[-1, 1]))/Norm(as.numeric(I_m124_1[-1, 1]))
R_m124_2 <- Norm(as.numeric(A40_m124_2.1[-1, 1]))/Norm(as.numeric(I_m124_2.1[-1, 1]))


dev.new(xpos=20, ypos=20, height=5, width=5)
plot(as.numeric(I_m124_1[-1, 1]), as.numeric(I_m124_2.2[-1, 1]), log='xy')
points(as.numeric(I_m124_1[-1, 1])[c(1)], as.numeric(I_m124_2.2[-1, 1])[c(1)],
       col="red")
dev.new(xpos=520, ypos=20, height=5, width=5)
plot(as.numeric(A40_m124_1[-1, 1]), as.numeric(A40_m124_2.2[-1, 1]), log='xy')
points(as.numeric(A40_m124_1[-1, 1])[c(1)], as.numeric(A40_m124_2.2[-1, 1])[c(1)], col="red")

dev.new(xpos=1020, ypos=20, height=5, width=5)
plot(R_m124_1, R_m124_2, log='xy')
points(R_m124_1[c(1)], R_m124_2[c(1)], log='xy', col="red")


dev.new(xpos=20, ypos=520, height=5, width=5)
plot(as.numeric(I_m1_1[-1, 1]), as.numeric(I_m124_1[-1, 1]), log='xy')
dev.new(xpos=520, ypos=520, height=5, width=5)
plot(as.numeric(I_m1_2[-1, 1]), as.numeric(I_m124_2.1[-1, 1]), log='xy')
dev.new(xpos=1020, ypos=520, height=5, width=5)
plot(as.numeric(I_m1_1[-1, 1]), as.numeric(I_m1_2[-1, 1]), log='xy')

dev.new(xpos=20, ypos=1020, height=5, width=5)
plot(as.numeric(I_m122[-1, 1]), as.numeric(I_m1_2[-1, 1]), log='xy')

I_all <- do.call("cbind", list(
  as.numeric(I_m1_1[-1, 1]),
  as.numeric(I_m122[-1, 1]),
  as.numeric(I_m1_2[-1, 1]),
  as.numeric(I_m124_1[-1, 1]),
  as.numeric(I_m124_2.1[-1, 1]),
  as.numeric(I_m124_2.2[-1, 1])
))

rownames(I_all) <- rownames(I_m1_1)[-1]
colnames(I_all) <- c("miR-1_1", "miR-122", "miR-1_2", "miR-124_1", "miR-124_2.1", "miR-124_2.2")





break



kds_1 <- EquilPars("miR-124", "equilibrium_2_tp", tp2rep=1)
kds_2 <- EquilPars("miR-124", "equilibrium_2_tp", tp2rep=2)
kds_orig <- EquilPars("miR-124", "equilibrium")

# break
kds_all <- EquilPars("miR-124", "equilibrium_2_tp")

site_names <- gsub("_Kd$", rownames(kds_1), replace="")
site_names <- gsub("_miR-124$", site_names, replace="")

cols <- kSiteColors[site_names]

sXc_original <- SitesXCounts("miR-124", n_constant=5, sitelist="paperfinal")
sXc_tp <- SitesXCounts("miR-124", experiment="equilibrium_tp", n_constant=5,
                       sitelist="paperfinal")
sXc_2_tp <- SitesXCounts("miR-124", experiment="equilibrium_2_tp", n_constant=5,
                         sitelist="paperfinal")

# dev.new(xpos=20, ypos=20, height=5, width=5)
# plot(sXc_original[, 1], sXc_original[, 2], log='xy', col=cols)
# dev.new(xpos=520, ypos=20, height=5, width=5)
# plot(sXc_original[, 2], sXc_tp[, 1], log='xy', col=cols)
# dev.new(xpos=1020, ypos=20, height=5, width=5)
# plot(sXc_original[, 2], sXc_2_tp[, 1], log='xy', col=cols)
# dev.new(xpos=20, ypos=520, height=5, width=5)
# plot(sXc_original[, 2], sXc_2_tp[, 2], log='xy', col=cols)


dev.new(xpos=20, ypos=20, height=5, width=5)
plot(kds_all[, 2], kds_1[, 2], log='xy', col=cols)
dev.new(xpos=520, ypos=20, height=5, width=5)
plot(kds_all[, 2], kds_2[, 2], log='xy', col=cols)
dev.new(xpos=1020, ypos=20, height=5, width=5)
plot(kds_1[, 2], kds_2[, 2], log='xy', col=cols)
dev.new(xpos=20, ypos=520, height=5, width=5)
plot(kds_orig[, 2], kds_all[, 2], log='xy', col=cols)


break








PlotSiteEnrichments <- function(mirna, experiment="equilibrium", n_constant=5,
                                sitelist="paperfinal", plotlist=FALSE,
                                uniq=FALSE, combined=TRUE, L=FALSE,
                                compcorrect=FALSE, wobble=FALSE, tpomit=FALSE,
                                buffer=FALSE, singleonly=FALSE,
                                plot.nconstant=FALSE, flowthrough=FALSE,
                                added.text=FALSE, datalines=FALSE,
                                modellines=TRUE, leg_rescale=1, height=4.5,
                                width=8.1, xpos=20, ypos=20, pdf.plot=FALSE) {
  # Load the count data and the model parameters:
  sXc <- SubfunctionCall(SitesXCounts)
  pars.matrix <- SubfunctionCall(EquilPars)
  # Correct mirna name for miR-7.
  if (mirna %in% c("miR-7-23nt", "miR-7-24nt", "miR-7-25nt")) {
    mirna_temp <- "miR-7"
  }
  # Removed the sites that are ascribed to competitor oligo enrichment.
  removed_sites <- GetRemovedSites(sXc)
  print(removed_sites)
  if (length(removed_sites) != 0) {
    # sXc <- sXc[setdiff(rownames(sXc), removed_sites), , drop=FALSE]
    kd.sites <- gsub("_Kd", replace="", rownames(pars.matrix))
    # pars.matrix <- pars.matrix[setdiff(rownames(pars.matrix), paste0(removed_sites, "_Kd")), ,drop=FALSE]
  }

  pars.model <- log10(pars.matrix$Mean)
  names(pars.model) <- rownames(pars.matrix)
  names(pars.model)[nrow(sXc) + 1] <- "bg"
  names(pars.model)[nrow(sXc) + 2] <- "AGO"
  data <- GetDataEquil(sXc)
  l <- SubfunctionCall(GetInputEquil)
  if (L) {
    l <- l/sum(l) * as.numeric(L)
  }
  # Ago dilutio in the data:
  pars <- pars.model

  A.dil.data <- sapply(gsub(",[12]$", colnames(data), replace="", perl=TRUE),
                       as.numeric)
  A.stock.measured <- kAgoStock[mirna, experiment]
  A.stock.pars <- 10^pars.model["AGO"]
  A.stock.pars <<- A.stock.pars
  pM_from_dil <- A.stock.measured*1000/100
  A.pM.data <- A.dil.data*pM_from_dil
  A.stock.measured <<- A.stock.measured
  xmin <- signif(min(unique(A.pM.data))/(10^0.25), 1)
  xmax <- signif(max(unique(A.pM.data))*(10^0.75), 1)
  print(A.pM.data)
  print(xmin)
  print(xmax)
  A.dil.model <- exp(seq(log(xmin/pM_from_dil),
                         log(xmin/pM_from_dil*10^(length(unique(A.pM.data))/2)),
                     length=100))
  model <- SubfunctionCall(EquilSingleSiteModelFreq, A.dil=A.dil.model)
  A.pM.model <-A.dil.model/100*A.stock.measured *1000
  data.R <- EquilEnrichments(data, l)
  model.R <- EquilEnrichments(model, l)
  # Set up the plotting limits
  if (flowthrough) {
    ymin <- 0.01
    ymax <- 1
  } else {
    ymin <- 0.2
    ymax <- 300
  }
  SubfunctionCall(GenericFigureSaveFile)
  # par(mar= c(3, 0.5, 2, 2))
  BlankPlot(log="xy", adjusted=TRUE)
  # Generate tickmarks for axis.
  xy <- GetPlotFractionalCoords(fx=0.05, fy=0.95, log='xy')
  title.text <- mirna

  if (added.text) {
    title.text <- paste0(mirna, "\n", experiment)
  }

  text(xy[1], xy[2], labels = title.text, adj=c(0, 1))

  mirna.trim <- paste0(strsplit(mirna, split = "-")[[1]][1:2],collapse="-")
  names(pars.model) <- gsub("(.*)_Kd", names(pars.model), replace="\\1")
  if (plot.nconstant) {
    xy <- GetPlotFractionalCoords(0.05, 0.90, log='xy')
    text(xy[1], xy[2], labels=sprintf("%s nt constant region", n_constant), adj=c(0, 1))      
  }

  site_order <- order(pars.matrix$Mean[1:(nrow(data)-1)])
  kd.matrix <- pars.matrix[1:(nrow(data)-1),][site_order, ]
  sites.ordered <- rownames(data)[site_order]
  kds <- kd.matrix$Mean

  # Part where text is formatted:
  kd_mags <- floor(log10(kds))
  # error_upper <- (kd.matrix$Upper_CI-kds)/(10^kd_mags)
  error_lower <- (kds - kd.matrix$Lower_CI)/(10^kd_mags)
  # error_mags_upper <- floor(log10(error_upper))
  error_mags_lower <- floor(log10(error_lower))
  temp_matrix <- rbind(kds/(10^kd_mags), c(error_lower), -error_mags_lower)
  cols <- kSiteColors[sites.ordered]
  names(cols) <- sites.ordered
  xy <- GetPlotFractionalCoords(fx=0.85, fy=1, log='xy')
  temp_legend <- Legend(xy,
                        legend=c(ConvertTtoUandMmtoX(sites.ordered), "None"),
                        col=c(cols, "black"),
                        y.intersp=0.9*leg_rescale)

  if (pdf.plot=="1.D" | pdf.plot=="1.D_uniq") {
    xy <- GetPlotFractionalCoords(fx=0.85, fy=0, log='xy')
    temp_legend <- legend(x=xy[1], y=xy[2], legend=rep("", length(sites.ordered) + 1),
                          col=c(kSiteColors[sites.ordered], "black"), seg.len=1,
                          pch=NA, lwd=1, bty="n",
                          y.intersp=0.86, yjust=0)

    xy <- GetPlotFractionalCoords(fx=0.85, fy=0.46, log='xy')
    text(x=xy[1], y=xy[2], labels=bquote("Relative"~italic(K)[D]*":"),
         adj=c(0, 0))
    pos1 <- 0.15
    text(10^(temp_legend$text$x + pos1),
         10^temp_legend$text$y[1:(length(sites.ordered) + 1)],
         sprintf("%.*f", c(temp_matrix[3, ], 1), c(temp_matrix[1, ], 1)), adj=1)
    pos2 <- pos1 + 0.15
    text(10^(temp_legend$text$x + pos2)[1:length(sites.ordered)],
         10^temp_legend$text$y[1:length(sites.ordered)],
         sapply(kd_mags, function(kd_mag) {
          as.expression(bquote(.("")%+-%.("")))
          }), adj=1)
    pos3 <- pos2 + 0.23
    text(10^(temp_legend$text$x + pos3)[1:length(sites.ordered)],
         10^temp_legend$text$y[1:length(sites.ordered)],
         sprintf("%.*f", temp_matrix[3, ], temp_matrix[2, ]), adj=1)
      pos4 <- pos3 + 0.46
    text(10^(temp_legend$text$x + pos4)[1:length(sites.ordered)],
         10^temp_legend$text$y[1:length(sites.ordered)],
         sapply(kd_mags, function(kd_mag) {
          as.expression(bquote(.("")%*%10^.(kd_mag)))
          }), adj=c(1, 0.4))
  }
  cols <- c(cols, None="black")
  if (flowthrough) {
    segments(A.pM.data, ymin, A.pM.data, ymax, lty=2, col="gray")
  }
  sapply(rownames(data), function(site) {
    if (!flowthrough) {
      Points(A.pM.data, data.R[site, ], col=cols[site], line=datalines)
    }
    if (!(datalines) & modellines) {
      lines(A.pM.model, model.R[site, ], col=cols[site], xpd=NA)        
    }
  })
  AddLogAxis(1, label=AGO_mir_label, adj=TRUE)
  AddLogAxis(2, label="Enrichment")

  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

PlotSiteEnrichments("miR-124", experiment="equilibrium_2_tp")

PlotSiteEnrichments("miR-124")

break

m124 <- SitesXCounts("miR-124")

m124_tp1 <- SitesXCounts("miR-124", experiment="equilibrium_tp")



m124_tp2 <- SitesXCounts("miR-124", experiment="equilibrium_2_tp")



kds_124 <- EquilPars("miR-124")
kds_124_tp1 <- EquilPars("miR-124", experiment="equilibrium_tp", tpomit=TRUE)
kds_124_tp2 <- EquilPars("miR-124", experiment="equilibrium_2_tp")




sites <- gsub("_Kd$", rownames(kds_124), replace="")
sites <- gsub("_miR-124$", sites, replace="")


dev.new(xpos=20, ypos=20, height=5, width=5)

plot(kds_124[, 1], kds_124_tp2[, 1], log='xy', col=kSiteColors[sites],
     xlim=c(7e-4, 10), ylim=c(7e-4, 10))

dev.new(xpos=520, ypos=20, height=5, width=5)

plot(kds_124[, 1], kds_124_tp1[, 1], log='xy', col=kSiteColors[sites],
     xlim=c(7e-4, 10), ylim=c(7e-4, 10))



break


R_m124 <- apply(m124, 2, GetEnrichment, x_I=m124[, 1])

R_m124_tp1 <- apply(m124_tp1, 2, GetEnrichment, x_I=m124_tp1[, 1])
R_m124_tp2 <- apply(m124_tp2, 2, GetEnrichment, x_I=m124_tp2[, 1])

colnames(R_m124_tp1) <- c("I", "I", "40", "12.65", "4", "4", "1.265", "1.265", "0.4", "0.1265", "0.04", "0")

colnames(R_m124_tp2) <- c("I", "I", "40", "40", "12.65", "12.65", "4", "4",
                          "1.265", "1.265", "0.4", "0.4", "0.1265", "0.1265",
                          "0", "0")

PlotEnrichments <- function(R_matrix, pdf.plot=FALSE, xpos=20, ypos=20) {
  R_matrix <- R_matrix[, 3:ncol(R_matrix)]
  xmin <- 0.0004
  xmax <- 1
  ymin <- 0.1
  ymax <- 1000
  SubfunctionCall(GenericFigureSaveFile)
  BlankPlot(log='xy')
  sites <- rownames(R_matrix)
  AddLogAxis(1, label="AGO2-miR-124 (%)", percent=TRUE)
  AddLogAxis(2, label="Enrichment")
  x_vals <- as.numeric(colnames(R_matrix))/100
  for (site in sites) {
    print(site)
    y_vals <- R_matrix[site, ]
    col <- kSiteColors[site]
    print(y_vals)
    print(col)
    lines(x_vals, y_vals, pch=20, col=col)
  }
}

PlotEnrichments(R_m124)


PlotEnrichments(R_m124_tp1, xpos=520)

PlotEnrichments(R_m124_tp2, xpos=1020)


break


m1_1_new <- rowSums(ReporterCounts("miR-1", "duplex_primer1", experiment="twist_reporter_assay_v2_test", mm=FALSE))
m1_2_new <- rowSums(ReporterCounts("miR-1", "duplex_primer2", experiment="twist_reporter_assay_v2_test", mm=FALSE))
m124_1_new <- rowSums(ReporterCounts("miR-124", "duplex_primer1", experiment="twist_reporter_assay_v2_test", mm=FALSE))
m124_2_new <- rowSums(ReporterCounts("miR-124", "duplex_primer2", experiment="twist_reporter_assay_v2_test", mm=FALSE))

m1_1 <- rowSums(ReporterCounts("miR-1", "duplex", rep=1, mm=FALSE))
m1_2 <- rowSums(ReporterCounts("miR-1", "duplex", rep=2, mm=FALSE))

l7_1 <- rowSums(ReporterCounts("let-7a", "duplex", rep=1, mm=FALSE))
l7_2 <- rowSums(ReporterCounts("let-7a", "duplex", rep=2, mm=FALSE))

m155_1 <- rowSums(ReporterCounts("miR-155", "duplex", rep=1, mm=FALSE))
m155_2 <- rowSums(ReporterCounts("miR-155", "duplex", rep=2, mm=FALSE))

m124_1 <- rowSums(ReporterCounts("miR-124", "duplex", rep=1, mm=FALSE))
m124_2 <- rowSums(ReporterCounts("miR-124", "duplex", rep=2, mm=FALSE))

l6_1 <- rowSums(ReporterCounts("lsy-6", "duplex", rep=1, mm=FALSE))
l6_2 <- rowSums(ReporterCounts("lsy-6", "duplex", rep=2, mm=FALSE))

m7_1 <- rowSums(ReporterCounts("miR-7", "duplex", rep=1, mm=FALSE))
m7_2 <- rowSums(ReporterCounts("miR-7", "duplex", rep=2, mm=FALSE))

nd_1 <- rowSums(ReporterCounts("miR-1", "no_duplex", rep=1, mm=FALSE))
nd_2 <- rowSums(ReporterCounts("miR-1", "no_duplex", rep=2, mm=FALSE))


all_old <- m1_1 + m1_2 + l7_1 + l7_2 + m155_1 + m155_2 + m124_1 + m124_2 + l6_1 + l6_2 + m7_1 + m7_2 + nd_1 + nd_2

all_new <- m1_1_new + m1_2_new + m124_1_new + m124_2_new

all_data <- cbind(m1_1, m1_2, l7_1, l7_2, m155_1, m155_2, m124_1, m124_2, l6_1, l6_2, m7_1, m7_2, nd_1, nd_2, m1_1_new, m1_2_new, m124_1_new, m124_2_new)

all_data_norm <- t(t(all_data)/colSums(all_data))


plot(all_data_norm[, 1], all_data_norm[, 15])




break
dev.new(xpos=20, ypos=20, height=5, width=5)
plot(all_old, all_new, log='xy')
dev.new(xpos=20, ypos=20, height=5, width=5)
plot(rowSums(all_old), rowSums(all_new), log='xy')

break

m124_1 <- ReporterCounts("miR-124", "duplex_primer1", experiment="twist_reporter_assay_v2_test", mm=FALSE)
m124_2 <- ReporterCounts("miR-124", "duplex_primer2", experiment="twist_reporter_assay_v2_test", mm=FALSE)


m1_1_nomm <- ReporterCounts("miR-1", "duplex_primer1", experiment="twist_reporter_assay_v2_test", mm=FALSE)





break
PlotReporterTemp <- function(mirna, new=TRUE, pdf.plot=FALSE) {
  mirna_opposites <- c("miR-1", "miR-124")
  print('hi')
  names(mirna_opposites) <- c("miR-124", "miR-1")
  mirna_opposite <- mirna_opposites[mirna]
  print(mirna_opposite)
  if (new) {
    print("true)")
    sig1 <- ReporterCounts(mirna, "duplex_primer1", experiment="twist_reporter_assay_v2_test", mm=TRUE)
    sig2 <- ReporterCounts(mirna, "duplex_primer1", experiment="twist_reporter_assay_v2_test", mm=TRUE)
    bg1 <- ReporterCounts(mirna_opposite, "duplex_primer1", experiment="twist_reporter_assay_v2_test", mm=TRUE)
    bg2 <- ReporterCounts(mirna_opposite, "duplex_primer1", experiment="twist_reporter_assay_v2_test", mm=TRUE)
  } else {
    print("false")
    sig1 <- ReporterCounts(mirna, "duplex", rep="2", experiment="twist_reporter_assay")
    sig2 <- sig1
    bg1 <- ReporterCounts(mirna_opposite, "duplex", rep="2", experiment="twist_reporter_assay")
    bg2 <- bg1
  }
  sig <- Norm(rowSums(sig1)) + Norm(rowSums(sig2))
  bg <- Norm(rowSums(bg1)) + Norm(rowSums(bg2))
  l2fc <- log(sig/bg, base=2)

  cols <- rep("gray", length(sig))

  mirnas <- sapply(rownames(sig1), function(name) {
    strsplit(name, split="_")[[1]][1]
  })

  sites <- sapply(rownames(sig1), function(name) {
    strsplit(name, split="_")[[1]][2]
  })

  inds_change <- which(mirnas == mirna)
  inds_opposite <- which(mirnas == mirna_opposite)

  if (mirna == "miR-1") {
    buffer <- TRUE
    combined <- FALSE
  } else {
    buffer <- FALSE
    combined <- TRUE
  }
  kds <- EquilPars(mirna, buffer=buffer, combined=combined)
  kds_all <- rep(1, length(m1))
  rownames(kds) <- gsub("_Kd", "", rownames(kds))
  names(kds_all) <- sites
  kds_all[inds_change] <- kds[sites[inds_change], 1]

  cols <- rep(ConvertRColortoRGB("gray", alpha=0.1), length(sig))

  cols[inds_change] <- kSiteColors[sites[inds_change]]

  y <- l2fc


  SubfunctionCall(GenericFigureSaveFile)
  xmin <- 0.0001
  xmax <- 2
  ymin <- -0.5
  ymax <- 0.1
  BlankPlot(log='x', inv='x')
  AddLogAxis(1, label="Kd")
  AddLinearAxis(2, label.space=0.1, tick.space=0.02, label="log2(fold change)")
  abline(0, 0, col="gray")
  points(kds_all[-inds_opposite], l2fc[-inds_opposite],
       col=cols[-inds_opposite])

    cutoff <- mean(l2fc[c(-inds_change, -inds_opposite)]) - 2*sd(l2fc[c(-inds_change, -inds_opposite)])
    abline(cutoff, 0, lty=2)

  legend.coords <- GetPlotFractionalCoords(0.0125, 0, log='x', inv='x')

  legend_inds <- which(l2fc[inds_change] <= cutoff)
  inds_final <- inds_change[legend_inds]
  Legend(legend.coords, legend=ConvertTtoUandMmtoX(sites[inds_final]),
         col=cols[inds_final], xjust=0, yjust=0, y.intersp=0.8)



# kd_line <- exp(seq(0, -8, length.out=100))
# l2fc_line <- log(kd_line)*m + b
# lines(kd_line, l2fc_line)


  if (class(pdf.plot) == "character") {
    dev.off()
  }

}

PlotReporterTemp("miR-1", new=FALSE, pdf.plot="190520/miR-1_old")
PlotReporterTemp("miR-124", new=FALSE, pdf.plot="190520/miR-124_old")
PlotReporterTemp("miR-1", new=TRUE, pdf.plot="190520/miR-1_new")
PlotReporterTemp("miR-124", new=TRUE, pdf.plot="190520/miR-124_new")

break

# m1_1 <- ReporterCounts("miR-1", "duplex", rep="1", experiment="twist_reporter_assay")
# m1_2 <- ReporterCounts("miR-1", "duplex", rep="2", experiment="twist_reporter_assay")
# m124_1 <- ReporterCounts("miR-124", "duplex", rep="1", experiment="twist_reporter_assay")
# m124_2 <- ReporterCounts("miR-124", "duplex", rep="2", experiment="twist_reporter_assay")



# dev.new(xpos=20, ypos=20, height=4, width=4)
# plot(m1_1 + 1, m1_2 + 1, log='xy')
# text(2, 500, labels=round(cor(log(c(m1_1 + 1)), log(c(m1_2 + 1)))^2, digits=3))

# dev.new(xpos=420, ypos=20, height=4, width=4)
# plot(m1_1, m124_1, log='xy')
# text(2, 500, labels=round(cor(log(c(m1_1 + 1)), log(c(m124_1 + 1)))^2, digits=3))


# dev.new(xpos=820, ypos=20, height=4, width=4)
# plot(m1_1, m124_2, log='xy')
# text(2, 500, labels=round(cor(log(c(m1_1 + 1)), log(c(m124_2 + 1)))^2, digits=3))


# dev.new(xpos=20, ypos=420, height=4, width=4)
# plot(m1_2, m124_1, log='xy')
# text(2, 500, labels=round(cor(log(c(m1_2 + 1)), log(c(m124_1 + 1)))^2, digits=3))

# dev.new(xpos=420, ypos=420, height=4, width=4)
# plot(m1_2, m124_2, log='xy')
# text(2, 500, labels=round(cor(log(c(m1_2 + 1)), log(c(m124_2 + 1)))^2, digits=3))

# dev.new(xpos=820, ypos=420, height=4, width=4)
# plot(m124_1, m124_2, log='xy')
# text(2, 500, labels=round(cor(log(c(m124_1 + 1)), log(c(m124_2 + 1)))^2, digits=3))




# break
# break
m1 <- Norm(rowSums(m1_2)) + Norm(rowSums(m1_1))

m124 <- Norm(rowSums(m124_1)) + Norm(rowSums(m124_2))


m1_alt <- Norm(m1_1 + 1) + Norm(m1_2 + 1)

m124_alt <- Norm(m124_1 + 1) + Norm(m124_2 + 1)

l2fc <- log(m1/m124, base=2)

# l2fc <- rowMeans(log(m1_alt/m124_alt, base=2))


cols <- rep("gray", length(m1))

mirnas <- sapply(rownames(m1_1), function(name) {
  strsplit(name, split="_")[[1]][1]
})

sites <- sapply(rownames(m1_1), function(name) {
  strsplit(name, split="_")[[1]][2]
})

print(head(mirnas))
print(head(sites))

mirna <- "miR-1"


mirna_opposite <- "miR-124"
inds_change <- which(mirnas == mirna)
inds_opposite <- which(mirnas == mirna_opposite)

if (mirna == "miR-1") {
  buffer <- TRUE
  combined <- FALSE
} else {
  buffer <- FALSE
  combined <- TRUE
}

kds <- EquilPars(mirna, buffer=buffer, combined=combined)

kds_all <- rep(1, length(m1))

rownames(kds) <- gsub("_Kd", "", rownames(kds))

print(kds_all)
names(kds_all) <- sites
kds_all[inds_change] <- kds[sites[inds_change], 1]

kds_all

lm_fit <- lm(l2fc[-inds_opposite] ~ log(kds_all)[-inds_opposite])
print(lm_fit)

b <- lm_fit$coefficients[1]
m <- lm_fit$coefficients[2]

print(b)
print(m)

l2fc_theory <- l


average_matrix <- m1_1 + m1_2 + m124_1 + m124_2


# Make the simulation matrix, which is the average of all of the counts.
frequency_matrix <- Norm(average_matrix)

# Get the log2fc of the count matrices, which uses the linear fit of the
# simulated count matrix.
l2fc_sim <- log(kds_all)*m + b

repression_matrix <- frequency_matrix*2^l2fc_sim

sums_bg <- rowSums(frequency_matrix)
sums_sig <- rowSums(repression_matrix)





cols[inds_change] <- kSiteColors[sites[inds_change]]
dev.new(xpos=20, ypos=20, height=4, width=4)
plot(kds_all[-inds_opposite], l2fc[-inds_opposite],
     col=cols[-inds_opposite], log='x', xlim=c(1, 0.0001), pch=19,
     ylim=c(-0.5, 0.1))

kd_line <- exp(seq(0, -8, length.out=100))
l2fc_line <- log(kd_line)*m + b
lines(kd_line, l2fc_line)


# dev.new(xpos=420, ypos=20, height=4, width=4)
# plot(kds_all[-inds_opposite], log(sums_sig/sums_bg, base=2)[-inds_opposite],
#        log='x', col=cols[-inds_opposite], pch=19, xlim=c(1, 0.0001),
#        ylim=c(-0.5, 0.1))
# lines(kd_line, l2fc_line)


# downsample <- 4e6
# downsample124 <- 1e5
# sampling_bg <- matrix(rmultinom(1, downsample, prob=frequency_matrix), 
#                       nrow=163, ncol=184)

# m1_1_molecules_sim <- matrix(rmultinom(1, downsample, prob=repression_matrix), 
#                       nrow=163, ncol=184)

# m1_2_molecules_sim <- matrix(rmultinom(1, downsample, prob=repression_matrix), 
#                       nrow=163, ncol=184)

# m1_1_sim <- matrix(rmultinom(1, sum(m1_1), prob=m1_1_molecules_sim), 
#                       nrow=163, ncol=184)

# m1_2_sim <- matrix(rmultinom(1, sum(m1_2), prob=m1_2_molecules_sim), 
#                       nrow=163, ncol=184)




# m124_1_molecules_sim <- matrix(rmultinom(1, downsample124, prob=frequency_matrix), 
#                       nrow=163, ncol=184)

# m124_2_molecules_sim <- matrix(rmultinom(1, downsample124, prob=frequency_matrix), 
#                       nrow=163, ncol=184)

# m124_1_sim <- matrix(rmultinom(1, sum(m124_1), prob=m124_1_molecules_sim), 
#                       nrow=163, ncol=184)

# m124_2_sim <- matrix(rmultinom(1, sum(m124_2), prob=m124_2_molecules_sim), 
#                       nrow=163, ncol=184)



# dev.new(xpos=20, ypos=420, height=4, width=4)
# plot(m1_1, m1_2, log='xy', xlim=c(1, 1e5), ylim=c(1, 1e5),
#      col=rgb(0, 0, 0, alpha=0.2))
# text(2, 500, labels=round(cor(log(c(m1_1 + 1)), log(c(m1_2 + 1)))^2, digits=3))

# dev.new(xpos=420, ypos=420, height=4, width=4)
# plot(m1_1_sim, m1_2_sim, log='xy', xlim=c(1, 1e5), ylim=c(1, 1e5),
#      col=rgb(0, 0, 0, alpha=0.2))
# text(2, 500, labels=round(cor(log(c(repression_matrix*sum(m1_1)+1)), log(c(sampling_bg + 1)))^2, digits=3))


# dev.new(xpos=20, ypos=820, height=4, width=4)
# plot(m124_1, m124_2, log='xy', xlim=c(1, 1e5), ylim=c(1, 1e5),
#      col=rgb(0, 0, 0, alpha=0.2))
# text(2, 500, labels=round(cor(log(c(m124_1 + 1)), log(c(m124_2 + 1)))^2, digits=3))

# dev.new(xpos=420, ypos=820, height=4, width=4)
# plot(m124_1_sim, m124_2_sim, log='xy', xlim=c(1, 1e5), ylim=c(1, 1e5),
#      col=rgb(0, 0, 0, alpha=0.2))
# text(2, 500, labels=round(cor(log(c(repression_matrix*sum(m124_1)+1)), log(c(sampling_bg + 1)))^2, digits=3))


# dev.new(xpos=820, ypos=820, height=4, width=4)
# plot(m1_1/m1_2, m124_1/m124_2, log='xy')


# break




# break


# break

mirna <- "miR-124"
mirna_opposite <- "miR-1"
cols <- rep("gray", length(m1))

inds_change <- which(mirnas == mirna)
inds_opposite <- which(mirnas == mirna_opposite)

l2fc <- log(m124/m1, base=2)



# l2fc <- rowMeans(log(m124_alt/m1_alt, base=2), na.rm=TRUE)


if (mirna == "miR-1") {
  buffer <- TRUE
  combined <- FALSE
} else {
  buffer <- FALSE
  combined <- TRUE
}

kds <- EquilPars(mirna, buffer=buffer, combined=combined,)

kds_all <- rep(1, length(m1))

rownames(kds) <- gsub("_Kd", "", rownames(kds))

print(kds_all)
names(kds_all) <- sites
kds_all[inds_change] <- kds[sites[inds_change], 1]

kds_all

cols[inds_change] <- kSiteColors[sites[inds_change]]
dev.new(xpos=520, ypos=20, height=5, width=5)
plot(kds_all[-inds_opposite], l2fc[-inds_opposite],
     col=cols[-inds_opposite], log='x', xlim=c(1, 0.0001), pch=20, ylim=c(-0.5, 0.1))


m1_1 <- ReporterCounts("miR-1", "duplex", rep="2", experiment="twist_reporter_assay")
m1_2 <- ReporterCounts("miR-1", "duplex", rep="2", experiment="twist_reporter_assay")
m124_1 <- ReporterCounts("miR-124", "duplex", rep="2", experiment="twist_reporter_assay")
m124_2 <- ReporterCounts("miR-124", "duplex", rep="2", experiment="twist_reporter_assay")


m1 <- Norm(rowSums(m1_2)) + Norm(rowSums(m1_1))

m124 <- Norm(rowSums(m124_1)) + Norm(rowSums(m124_2))

cols <- rep("gray", length(m1))

mirnas <- sapply(rownames(m1_1), function(name) {
  strsplit(name, split="_")[[1]][1]
})

sites <- sapply(rownames(m1_1), function(name) {
  strsplit(name, split="_")[[1]][2]
})

print(head(mirnas))
print(head(sites))

mirna <- "miR-1"


mirna_opposite <- "miR-124"
inds_change <- which(mirnas == mirna)
inds_opposite <- which(mirnas == mirna_opposite)

if (mirna == "miR-1") {
  buffer <- TRUE
  combined <- FALSE
} else {
  buffer <- FALSE
  combined <- TRUE
}

kds <- EquilPars(mirna, buffer=buffer, combined=combined)

kds_all <- rep(1, length(m1))

rownames(kds) <- gsub("_Kd", "", rownames(kds))

print(kds_all)
names(kds_all) <- sites
kds_all[inds_change] <- kds[sites[inds_change], 1]

kds_all

cols[inds_change] <- kSiteColors[sites[inds_change]]
dev.new(xpos=20, ypos=520, height=5, width=5)
plot(kds_all[-inds_opposite], log(m1/m124, base=2)[-inds_opposite],
     col=cols[-inds_opposite], log='x', xlim=c(1, 0.0001), pch=20,
     ylim=c(-0.5, 0.1))


mirna <- "miR-124"
mirna_opposite <- "miR-1"
cols <- rep("gray", length(m1))

inds_change <- which(mirnas == mirna)
inds_opposite <- which(mirnas == mirna_opposite)


if (mirna == "miR-1") {
  buffer <- TRUE
  combined <- FALSE
} else {
  buffer <- FALSE
  combined <- TRUE
}

kds <- EquilPars(mirna, buffer=buffer, combined=combined,)

kds_all <- rep(1, length(m1))

rownames(kds) <- gsub("_Kd", "", rownames(kds))

print(kds_all)
names(kds_all) <- sites
kds_all[inds_change] <- kds[sites[inds_change], 1]

kds_all

cols[inds_change] <- kSiteColors[sites[inds_change]]
dev.new(xpos=520, ypos=520, height=5, width=5)
plot(kds_all[-inds_opposite], log(m124/m1, base=2)[-inds_opposite],
     col=cols[-inds_opposite], log='x', xlim=c(1, 0.0001), pch=20, ylim=c(-0.5, 0.1))




break


break

MakeCompetitorPlotSingle("miR-124", 13, experiment="equilibrium_tp")
MakeCompetitorPlotSingle("miR-124", 13, experiment="equilibrium", xpos=520)

break

PlotSiteKds("miR-124", "equilibrium", 5, "paperfinal", singleonly=TRUE)

PlotSiteKds("miR-124", "equilibrium_tp", 5, "paperfinal", singleonly=FALSE, tpomit=TRUE)


PlotSiteEnrichments("miR-124", "equilibrium_tp", 5, "paperfinal",
                    singleonly=FALSE, tpomit=TRUE)

PlotSiteEnrichments("miR-124", "equilibrium", 5, "paperfinal",
                    singleonly=FALSE)


PlotPairwiseKds("miR-124", "equilibrium", 5, "paperfinal", "miR-124",
                "equilibrium_tp", 5, "paperfinal", tpomit2=TRUE)


# PlotPairwiseKds("miR-124", "equilibrium", 5, "paperfinal", "miR-124",
#                 "equilibrium_tp", 5, "paperfinal")



MakeCompetitorPlot("miR-124", experiment="equilibrium")

MakeCompetitorPlot("miR-124", experiment="equilibrium_tp")

break


PlotSiteKds <- function(mirna, experiment="equilibrium", n_constant=5,
                        sitelist="paperfinal", combined=TRUE, uniq=FALSE, singleonly=TRUE,
                        papersites=FALSE, collapse_AA=FALSE,
                        adjusted_height=FALSE, added.text=FALSE, L=FALSE,
                        trim_mir_name=TRUE, buffer=FALSE, compcorrect=FALSE, wobble=FALSE,
                        plot.nconstant=FALSE, height=4.5, width=6, xpos=20,
                        ypos=20, pdf.plot=FALSE) {
  # Get kds for all site-types of the mirna.
  sXc <- SubfunctionCall(SitesXCounts)
  pars <- SubfunctionCall(EquilPars)
  rownames(pars) <- gsub("_Kd", replace="", rownames(pars))
  # Assign removed_sites
  removed_sites <- GetRemovedSites(sXc)
  # Check for the AA- dinucleotides conditional
  if (collapse_AA) {
    AA_sites <- grep("^AA-", rownames(sXc), perl=TRUE, value=TRUE)
    # Get the base names of all sites that have both an AA and a non AA site.
    AA_shared_sites <- intersect(gsub("AA-", replace="", AA_sites),
                                 rownames(sXc))
    # Get the separate Kd matrix.
    kds_AAbase <- pars[AA_shared_sites, ]
    # Replace the values of the base Kds with that of the AA-base Kds. 
    pars[AA_shared_sites, ] <- pars[paste0("AA-", AA_shared_sites), ]
    # Add all of the AA sites to the list of sites to be removed.
    removed_sites <- c(AA_sites, removed_sites)
  }
  # Remove the sites to be taken out of the matrix.
  sXc <- sXc[setdiff(rownames(sXc), removed_sites), , drop=FALSE]
  pars <- pars[setdiff(rownames(pars), removed_sites), ,drop=FALSE]
  # Subset the pars matrix to make the kd matrix, order it and define the sites. 
  kds <- pars[1:(nrow(sXc) - 1),]
  kds <- kds[order(kds$Mean),]
  sites <- rownames(kds)

  # Set up the plotting limits.
  xmin <- 1e-4  
  xmax <- 1
  ymin <- 0
  ymax <- nrow(kds) + 0.5
  # Make the y values for the major set of kds to be plotted.
  # Give the parameters allowing for the scaling of the plot window in the event
  # that adjusted height is true. 
  site_line_thickness <- 0.68
  par_margin_bottom <- 3
  par_margin_top <- 2
  if (adjusted_height) {
    height <- (ymax*site_line_thickness + par_margin_bottom + par_margin_top)/4
  }
  # Set up the plot.
  SubfunctionCall(FigureSaveFile)
  # Modify the plot to have different margins.
  par(mar=c(par_margin_bottom, 0.6, par_margin_top, 2))
  BlankPlot(log='x', inv='x')
  # Get the
  x <- kds$Mean
  y <- rev(seq(nrow(kds)))
  names(y) <- sites

  x_error <- list(kds$Upper_CI, kds$Lower_CI)
  cols <- sapply(sites, function(site) {
    print(site)
    if (site %in% names(kSiteCategories)) {
      return(kSiteCategoryColors[kSiteCategories[site]])
    } else {
      return(kSiteCategoryColors["Noncanonical"])
    }
  })
  names(cols) <- sites
  # Check for both conditionals of AA loop:
  if (collapse_AA) {
    # Get the y-values for base sites of the AA-
    y_AAbase <- y[AA_shared_sites]
    x_error_AAbase <- list(kds_AAbase$Lower_CI, kds_AAbase$Upper_CI)
    cols_AAbase <- cols[AA_shared_sites]
    # This makes the lines and labels the base AA sites.
    sapply(rownames(kds_AAbase),function(site) {
      # Get these reference points to plot the function.
      kd <- kds_AAbase[site, ]
      y <- y_AAbase[site]
      label <- ConvertTtoUandMmtoX(site)
      x_seg_l <- kd[2]
      x_seg_r <- kds[site, ][2]

      if (site %in% c("8mer-bT6", "8mer-bC(4.6)")) {
        x <- min(kd[2]/1.6, kd[3]/1.1)
        print(x)
        x <- as.numeric(kd[2]/2)
        print(x)
        adj <- 0
        # x_seg_l <- kd[2]
        segments(unlist(x_seg_l), y, x, y, lty=2)
        str.width <- strwidth(label, units="figure")
        print(str.width)
        break
        pos_start <- (log(x) - log(xmax))/(log(xmin) - log(xmax))
        pos_end <- pos_start + str.width + 0.04
        x_r <- GetPlotFractionalCoords(pos_end, 1, log='x', inv='x')
        segments(x_r[1], y, unlist(x_seg_r), y, lty=2)
      } else {
        x <- kd[5]* 1.1
        adj <- 1
        segments(unlist(x_seg_l), y, unlist(x_seg_r), y, lty=2)
      }
      text_out <- text(x, y, labels=label, adj=adj)
    })
    # Plot the points.  
    Points(kds_AAbase$Mean, y_AAbase, col=cols_AAbase, x_error=x_error_AAbase)
  }
  Points(x, y, col=cols, x_error=x_error)
  AddLogAxis(1, label="Relative Kd")
  title.xy <- GetPlotFractionalCoords(fx=0.05, fy=1, log='x', inv='x')
    if (length(grep("miR-7", mirna)) > 0 & trim_mir_name) {
      title.text <- "miR-7"
    } else {
      title.text <- mirna
    }
  if (added.text) {
    title.text <- paste0(mirna, "\n", experiment)
  }
  if (plot.nconstant) {
    xy <- GetPlotFractionalCoords(0.05, 0.90, inv='x')
    points(c(xy[1]), c(xy[2]), pch=3, col="gray")
    text(xy[1], xy[2], labels=sprintf("%s nt constant region", n_constant), adj=c(0, 1))      
  }
  print(376)
  text(title.xy[1], max(y), labels = title.text, adj=c(0, 0.5), xpd=NA)
  title_pos <- apply(kds, 1, function(kd) {
    min(kd[2]/1.6, kd[3]/1.1)
  })
  if (collapse_AA) {
    if (nrow(kds_AAbase) !=0) {
      AA_inds <- which(sites %in% rownames(kds_AAbase))
      sites[AA_inds] <- paste0("AA-", sites[AA_inds], sep="")  
    }    
  }
  par(xpd=NA) 
  text(title_pos, y, labels=ConvertTtoUandMmtoX(sites), adj=0)
  legend.names <- c("7-8-nt canonical site",
                    "6-nt canonical site",
                    "Enhanced 6mer site",
                    "Noncanonical site",
                    "3'-only site")
  legend.colors <- c(kSiteCategoryColors[c(1, 2, 3, 7, 4)])
  mirna.sites <- list(c(1, 2, 4), # miR-1
                      c(1, 2, 4), # let-7a
                      c(1, 2, 3, 4, 5), # miR-155
                      c(1, 2, 3, 4, 5), # miR-124
                      c(1, 2, 3, 4, 5), # lsy-6
                      c(1, 2, 3, 4), # miR-7-22nt
                      c(1, 2, 3, 4), # miR-7-23nt
                      c(1, 2, 3, 4), # miR-7-24nt
                      c(1, 2, 3, 4)) # miR-7-25nt
  names(mirna.sites) <- c(kMirnas, "miR-7-22nt", "miR-7-24nt", "miR-7-25nt")
  xy <- GetPlotFractionalCoords(fx=0.65, fy=0, log='x', inv='x')      
  if (!(added.text)) {
    Legend(xy, legend=legend.names[mirna.sites[[mirna]]],
           col=legend.colors[mirna.sites[[mirna]]],
           yjust=0)    
  }
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}



PlotSiteKds("miR-124", sitelist="paperfinal", collapse_AA=TRUE, 
            compcorrect=TRUE, adjusted_height=TRUE, xpos=20, ypos=20)

PlotSiteKds("miR-124", sitelist="paperfinal", collapse_AA=TRUE, 
            compcorrect=TRUE, xpos=620, ypos=20)

PlotSiteKds("miR-1", sitelist="paperfinal", buffer=TRUE, combined=FALSE,
            xpos=820, ypos=20, pdf.plot="1.F")







break


CheckTwoSitesDontPass <- function() {

  gray_sites <- c(`miR-155`="AACGAGG_Kd", `miR-124`="TCACCCGC_Kd")
  gray_sites <<- gray_sites
  for (mirna in c("miR-155", "miR-124")) {
    if (mirna == "miR-124") compcorrect <- TRUE
    else                    compcorrect <- FALSE
    pars.matrix <- EquilPars(mirna, experiment="equilibrium", n_constant=5,
                   sitelist="paperfinal", combined=TRUE, singleonly=TRUE,
                   compcorrect=compcorrect)
    kd.matrix <- pars.matrix[1:(nrow(pars.matrix)-3),]
    kds <- kd.matrix$Mean
    names(kds) <- rownames(kd.matrix)
    kd_mags <- floor(log10(kds))
    error_lower <- (kds - kd.matrix$Lower_CI)/(10^kd_mags)
    error_mags_lower <- floor(log10(error_lower))
    temp_matrix <- rbind(kds/(10^kd_mags), c(error_lower), -error_mags_lower)
    colnames(temp_matrix) <- rownames(kd.matrix)
    formatted_kds <- sprintf("%.*f +/- %.*f x 10^%s", c(temp_matrix[3, ]),
                             c(temp_matrix[1, ]), c(temp_matrix[3, ]),
                             c(temp_matrix[2, ]), kd_mags)
    names(formatted_kds) <- rownames(kd.matrix)
    gray_site <- gray_sites[mirna]
    print(formatted_kds[gray_site])
    print(kds[gray_site]/kds["6mer_Kd"])
    print(1/kds[gray_site])
  }
}

CheckTwoSitesDontPass()


break



date <- "190409"
off <- 0
print(off)

sXc <- SitesXCounts("miR-124", sitelist="paperfinal")




sXf <- t(t(sXc)/colSums(sXc))
sXr <- sXf/sXf[, 1]

print(sXr)


conditions <- c("I", "0.4", "1.26", "4", "12.6", "40")

vec_cond_rep <- c("0.4", "1.265", "4", "12.65", "40")


sXr_new_nowob <- sXr
sXr_new_wob <- sXr

for (row_i in grep("AA-", rownames(sXc))) {
  str_site_seq <- GetSiteSeq("miR-124", rownames(sXc)[row_i])
  kmer_len_wob <- GetMaxCompetitorOligoEnrichment("miR-124", str_site_seq,
                                              wobble=TRUE)
  kmer_len_nowob <- GetMaxCompetitorOligoEnrichment("miR-124", str_site_seq,
                                              wobble=FALSE)
  kXc_nowob <- sapply(conditions, function(condition) {
    # Get the total number of kmers.
    counts <- CompetitorKmers("miR-124", condition, kmer_len_nowob, off=off,
                              experiment="equilibrium", wobble=FALSE)
    total_counts <- counts["Total", ]
    kmers <- unlist(counts)
    fracs <- kmers[-length(kmers)]/total_counts
    names(fracs) <- rownames(counts)[-length(kmers)]
    fracs
  })

  kXc_wob <- sapply(conditions, function(condition) {
    # Get the total number of kmers.
    counts <- CompetitorKmers("miR-124", condition, kmer_len_wob, off=off,
                              experiment="equilibrium", wobble=TRUE)
    total_counts <- counts["Total", ]
    kmers <- unlist(counts)
    fracs <- kmers[-length(kmers)]/total_counts
    names(fracs) <- rownames(counts)[-length(kmers)]
    fracs
  })


  kXr_nowob <- kXc_nowob/kXc_nowob[, 1]
  vec_r_nowob <- colMeans(kXr_nowob)
  vec_r_rep_nowob <- vec_r_nowob[2:length(vec_r_nowob)] - sXr["None", vec_cond_rep]
  sXr_new_nowob[row_i, vec_cond_rep] <- sXr_new_nowob[row_i, vec_cond_rep] - vec_r_rep_nowob

  kXr_wob <- kXc_wob/kXc_wob[, 1]
  vec_r_wob <- colMeans(kXr_wob)
  vec_r_rep_wob <- vec_r_wob[2:length(vec_r_wob)] - sXr["None", vec_cond_rep]
  sXr_new_wob[row_i, vec_cond_rep] <- sXr_new_wob[row_i, vec_cond_rep] - vec_r_rep_wob

}

dev.new(xpos=520, ypos=20, height=5, width=5)
plot(sXr[, "40"], sXr_new_nowob[, "40"], log='xy', col=kSiteColors[rownames(sXr)])
# dev.new(xpos=520, ypos=20, height=5, width=5)
# plot(sXr[, "40"], sXr_new_wob[, "40"], log='xy', col=kSiteColors[rownames(sXr)])

break


###################### let-7a ###################################################

mirna1 <- "let-7a"
experiment1 <- "equilibrium"
mirna2 <- "let-7a-21nt"
experiment2 <- "equilibrium_nb"
pdf.plot <- sprintf("%s/%s_%s_vs_%s_%s_kds.pdf", date, mirna1, experiment1,
                    mirna2, experiment2)
PlotPairwiseKds(mirna1=mirna1, experiment1=experiment1, n_constant=5, 
                sitelist1="paperfinal", mirna2=mirna2, experiment2=experiment2,
                n_constant2=5, sitelist2="paperfinal", combined1=TRUE,
                buffer1=FALSE, combined2=TRUE, buffer2=FALSE,
                mirna_start1=FALSE, mirna_start2=FALSE, pdf.plot=pdf.plot)
pdf.plot <- sprintf("%s/%s_%s_competitor_enrichment.pdf", date, mirna1,
                    experiment1)
MakeCompetitorPlot(mirna1, experiment1, pdf.plot=pdf.plot)
pdf.plot <- sprintf("%s/%s_%s_competitor_enrichment_geo.pdf", date, mirna1,
                    experiment1)
MakeCompetitorPlot(mirna1, experiment1, geomean=TRUE, pdf.plot=pdf.plot)
pdf.plot <- sprintf("%s/%s_%s_competitor_enrichment.pdf", date, mirna2,
                    experiment2)
MakeCompetitorPlot(mirna2, experiment2, pdf.plot=pdf.plot)
pdf.plot <- sprintf("%s/%s_%s_competitor_enrichment_geo.pdf", date, mirna2,
                    experiment2)
MakeCompetitorPlot(mirna2, experiment2, geomean=TRUE, pdf.plot=pdf.plot)

###################### miR-1 ###################################################

mirna1 <- "miR-1"
mirna2 <- mirna1
experiment1 <- "equilibrium"
experiment2 <- "equilibrium_met_tp"
pdf.plot <- sprintf("%s/%s_%s_vs_%s_%s_kds.pdf", date, mirna1, experiment1,
                    mirna2, experiment2)
PlotPairwiseKds(mirna1=mirna1, experiment1=experiment1, n_constant=5, 
                sitelist1="paperfinal", mirna2=mirna2, experiment2=experiment2,
                n_constant2=5, sitelist2="paperfinal", combined1=FALSE,
                buffer1=TRUE, combined2=FALSE, buffer2=FALSE,
                mirna_start1=FALSE, mirna_start2=FALSE, pdf.plot=pdf.plot)
pdf.plot <- sprintf("%s/%s_%s_competitor_enrichment.pdf", date, mirna1,
                    experiment1)
MakeCompetitorPlot(mirna1, experiment1, pdf.plot=pdf.plot)
pdf.plot <- sprintf("%s/%s_%s_competitor_enrichment_geo.pdf", date, mirna1,
                    experiment1)
MakeCompetitorPlot(mirna1, experiment1, geomean=TRUE, pdf.plot=pdf.plot)
pdf.plot <- sprintf("%s/%s_%s_competitor_enrichment.pdf", date, mirna2,
                    experiment2)
MakeCompetitorPlot(mirna2, experiment2, pdf.plot=pdf.plot)
pdf.plot <- sprintf("%s/%s_%s_competitor_enrichment_geo.pdf", date, mirna2,
                    experiment2)
MakeCompetitorPlot(mirna2, experiment2, geomean=TRUE, pdf.plot=pdf.plot)

###################### miR-7 ###################################################

# First make all the competitor oligo plots.
experiment <- "equilibrium_nb"
for (mirna in c("miR-7-22nt", "miR-7-23nt", "miR-7-24nt", "miR-7-25nt")) {
  pdf.plot <- sprintf("%s/%s_%s_competitor_enrichment.pdf", date, mirna,
                    experiment)
  MakeCompetitorPlot(mirna, experiment, pdf.plot=pdf.plot)
  pdf.plot <- sprintf("%s/%s_%s_competitor_enrichment_geo.pdf", date, mirna,
                    experiment)
  MakeCompetitorPlot(mirna, experiment, geomean=TRUE, pdf.plot=pdf.plot)
}

experiment <- "equilibrium2_nb"
for (mirna in c("miR-7-23nt", "miR-7-24nt", "miR-7-25nt")) {
  pdf.plot <- sprintf("%s/%s_%s_competitor_enrichment.pdf", date, mirna,
                    experiment)
  MakeCompetitorPlot(mirna, experiment, pdf.plot=pdf.plot)
  pdf.plot <- sprintf("%s/%s_%s_competitor_enrichment_geo.pdf", date, mirna,
                    experiment)
  MakeCompetitorPlot(mirna, experiment, geomean=TRUE, pdf.plot=pdf.plot)
}

mirna <- "miR-7-24nt"
experiment <- "equilibrium3_nb"
pdf.plot <- sprintf("%s/%s_%s_competitor_enrichment.pdf", date, mirna,
                  experiment)
MakeCompetitorPlot(mirna, experiment, pdf.plot=pdf.plot)
pdf.plot <- sprintf("%s/%s_%s_competitor_enrichment_geo.pdf", date, mirna,
                  experiment)
MakeCompetitorPlot(mirna, experiment, geomean=TRUE, pdf.plot=pdf.plot)


MakeCompetitorPlot("miR-7-22nt", "equilibrium_nb", pdf.plot=pdf.plot)
MakeCompetitorPlot(mirna1, experiment1, pdf.plot=pdf.plot)

# Now make the pairwise kd plots.

mirna1 <- "miR-7-22nt"
mirna2 <- "miR-7-23nt"
experiment1 <- "equilibrium_nb"
experiment2 <- experiment1
pdf.plot <- sprintf("%s/%s_%s_vs_%s_%s_kds.pdf", date, mirna1, experiment1,
                    mirna2, experiment2)
PlotPairwiseKds(mirna1=mirna1, experiment1=experiment1, n_constant=5, 
                sitelist1="paperfinal", mirna2=mirna2, experiment2=experiment2,
                n_constant2=5, sitelist2="paperfinal", combined1=TRUE,
                buffer1=FALSE, combined2=TRUE, buffer2=FALSE,
                mirna_start1=FALSE, mirna_start2=FALSE, pdf.plot=pdf.plot)

mirna2 <- "miR-7-24nt"
pdf.plot <- sprintf("%s/%s_%s_vs_%s_%s_kds.pdf", date, mirna1, experiment1,
                    mirna2, experiment2)
PlotPairwiseKds(mirna1=mirna1, experiment1=experiment1, n_constant=5, 
                sitelist1="paperfinal", mirna2=mirna2, experiment2=experiment2,
                n_constant2=5, sitelist2="paperfinal", combined1=TRUE,
                buffer1=FALSE, combined2=TRUE, buffer2=FALSE,
                mirna_start1=FALSE, mirna_start2=FALSE, pdf.plot=pdf.plot)

mirna2 <- "miR-7-25nt"
pdf.plot <- sprintf("%s/%s_%s_vs_%s_%s_kds.pdf", date, mirna1, experiment1,
                    mirna2, experiment2)
PlotPairwiseKds(mirna1=mirna1, experiment1=experiment1, n_constant=5, 
                sitelist1="paperfinal", mirna2=mirna2, experiment2=experiment2,
                n_constant2=5, sitelist2="paperfinal", combined1=TRUE,
                buffer1=FALSE, combined2=TRUE, buffer2=FALSE,
                mirna_start1=FALSE, mirna_start2=FALSE, pdf.plot=pdf.plot)

mirna1 <- "miR-7-23nt"
mirna2 <- "miR-7-23nt"
experiment1 <- "equilibrium2_nb"
experiment2 <- "equilibrium_nb"
pdf.plot <- sprintf("%s/%s_%s_vs_%s_%s_kds.pdf", date, mirna1, experiment1,
                    mirna2, experiment2)
PlotPairwiseKds(mirna1=mirna1, experiment1=experiment1, n_constant=5, 
                sitelist1="paperfinal", mirna2=mirna2, experiment2=experiment2,
                n_constant2=5, sitelist2="paperfinal", combined1=TRUE,
                buffer1=FALSE, combined2=TRUE, buffer2=FALSE,
                mirna_start1=FALSE, mirna_start2=FALSE, pdf.plot=pdf.plot)

mirna2 <- "miR-7-25nt"
pdf.plot <- sprintf("%s/%s_%s_vs_%s_%s_kds.pdf", date, mirna1, experiment1,
                    mirna2, experiment2)
PlotPairwiseKds(mirna1=mirna1, experiment1=experiment1, n_constant=5, 
                sitelist1="paperfinal", mirna2=mirna2, experiment2=experiment2,
                n_constant2=5, sitelist2="paperfinal", combined1=TRUE,
                buffer1=FALSE, combined2=TRUE, buffer2=FALSE,
                mirna_start1=FALSE, mirna_start2=FALSE, pdf.plot=pdf.plot)

mirna1 <- "miR-7-24nt"
mirna2 <- "miR-7-24nt"
experiment1 <- "equilibrium2_nb"
experiment2 <- "equilibrium3_nb"






break





MakeCompSimScatter <- function(mirna, condition, addcomp,
                               experiment="equilibrium", n_constant=5, off=0,
                               kmer_len_start=4, kmer_len_stop=16, geomean=TRUE,
                               xpos=20, ypos=520, height=5, width=5,
                               pdf.plot=FALSE) {
  if (mirna == "miR-7-23nt") {
    conditions <- c("I", "0", "0.4", "1.26", "12.6", "40")
  } else {
    conditions <- c("I", "0", "0.4", "1.26", "4", "12.6", "40")
  }
  conditions <- c("I", condition, condition)
  addcomps <- c(0, 0, addcomp)
  condition_comps <- cbind(conditions, addcomps)
  print(condition_comps)
  Rs <- sapply(kmer_len_start:kmer_len_stop, function(kmer_len) {
    message(kmer_len)
    counts_outer <- apply(condition_comps, 1, function(cond_comps) {
      condition <- cond_comps[1]
      addcomp <- cond_comps[2]
      message(condition)
      print(addcomp)
      # Get the path to file.
      # path <- GetAnalysisPath(mirna, experiment, condition, analysis_type="reads")
      # print(path)
      # Get the total number of kmers.
      # tot_kmers <- CountLines(path)*(37 + 2*n_constant - kmer_len + 1)
      # tot_kmers = tot_kmers*(100 + addcomp)/100
      # print(tot_kmers)
      # Get the

      counts <- CompetitorKmers(mirna, condition, kmer_len, off=off,
                                experiment=experiment, addcomp=addcomp)
      fracs <- counts[1:(nrow(counts) - 1), ]/counts["Total", ]
      print(fracs)
      names(fracs) <- rownames(counts)[1:(nrow(counts) - 1)]
      fracs
    })
    print(counts_outer)
    colnames(counts_outer) <- c("I", condition, sprintf("%s_sim", condition))
    R <- counts_outer/counts_outer[, "I"]
    print(R)
    if (geomean) {
      R_vec <- apply(R, 2, function(column) {
        GeoMean(column[is.finite(column)])
      })
    } else {
      R_vec <- apply(R, 2, function(column) {
        mean(column[is.finite(column)])
      })      
    }
    print(R_vec)
    R_vec
  })

  Rs <<- Rs
  colnames(Rs) <- sprintf("k%s", seq(kmer_len_start, kmer_len_stop))
  xmin <- 3
  xmax <- 17
  ymin <- 0
  if (mirna == "lsy-6") {
    ymax <- 100  
  } else {
    ymax <- 50
  }
  ymax <- 10^ceiling(log10(max(Rs)))
  SubfunctionCall(GenericFigureSaveFile)
  BlankPlot()
  AddLinearAxis(1, tick.space=1, label.space=1,
                label="Length of complementarity")
  AddLinearAxis(2, tick.space=1, label.space=5,
                label="Enrichment")
  kEquilCondColors[sprintf("%s_sim", condition)] <- "grey"
  for (condition_i in rownames(Rs)) {
    print(condition_i)
    print(kEquilCondColors[condition_i])
    lines(kmer_len_start:kmer_len_stop, Rs[condition_i, ], type="o",
          col=kEquilCondColors[condition_i])
    if (condition_i == condition) {
      lines(kmer_len_start:kmer_len_stop, 2*Rs[condition_i, ], type="o",
            col=kEquilCondColors[condition_i], lty=2, lwd=2)
    }
  }
  xy <- GetPlotFractionalCoords(0.05, 0.95)
  legend(x=xy[1], y=xy[2], bty="n", legend=conditions, pch=19,
         col=kEquilCondColors[rownames(Rs)])
  xy <- GetPlotFractionalCoords(0.95, 0.95)
  text(xy[1], xy[2], labels=mirna, adj=1)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

# MakeCompSimScatter("miR-1", 40, "0.000135", ypos=-200)



# break



MakeCompetitorPlot("miR-1", "equilibrium", pdf.plot="190405/miR-1_sm_competitor_enrichment")



MakeCompetitorPlot("miR-1", "equilibrium_tp", pdf.plot="190405/miR-1_tp_competitor_enrichments")

break


date <- "190329"
off <- 1
for (off in seq(1, 3)) {
  print(off)
  for (mirna in kMirnas) {
    print(mirna)
    if (mirna == "miR-7-23nt") {
      experiment <- "equilibrium2_nb"
    } else {
      experiment <- "equilibrium"
    }
    print(experiment)
    pdf.plot <- sprintf("%s/%s_off%s_geomean", date, mirna, off)
    print(pdf.plot)
    MakeCompetitorPlot(mirna, experiment=experiment, off=off,
                       pdf.plot=pdf.plot)
    pdf.plot <- sprintf("%s/%s_off%s_arithmean", date, mirna, off)
    print(pdf.plot)
    MakeCompetitorPlot(mirna, experiment=experiment, off=off, geomean=FALSE,
                       pdf.plot=pdf.plot)
  }  
}


break


GrepCompetitorOligoKmers <- function(mirna, experiment="equilibrium",
                                     n_constant=5, kmer_len_start=4,
                                     kmer_len_stop=12) {
  # Get the competitor oligo sequences, and the subregion for which the averages
  # are desired.
  competitor_seq <- SequenceList[["competitor"]][mirna]
  comp_seq <- substr(competitor_seq, 14, nchar(competitor_seq))
  # Get the sequences of the flanking region to append to each read.
  lib_seq_5p <- kLibrarySeqs5P[mirna]
  lib_seq_3p <- kLibrarySeqs3P[mirna]
  flank_5p <- substr(ConvertUtoT(lib_seq_5p),
                     start=nchar(lib_seq_5p) - n_constant + 1,
                     stop=nchar(lib_seq_5p))
  flank_3p <- substr(ConvertUtoT(lib_seq_3p),
                     start=1,
                     stop=5)
  # Make the empty lists to be added to during the function.
  # Use the "I" and particular condition to get the paths to the reads:
  conditions <- c("I", "0.4", "1.26", "4", "12.6", "40")

  Rs <- sapply(kmer_len_start:kmer_len_stop, function(kmer_len) {
    message(kmer_len)
    competitor_kmers <- GetKmersFromString(comp_seq, len_k=kmer_len)
    rna_kmers <- sapply(competitor_kmers, RevComplement)
    counts_outer <- sapply(conditions, function(condition) {
      message(condition)
      # Get the path to file.
      path <- GetAnalysisPath(mirna, experiment, condition, analysis_type="reads")
      # Get the total number of kmers.
      tot_kmers <- CountLines(path)*(37 + 2*n_constant - kmer_len + 1)
      # Get the 
      counts <- sapply(rna_kmers, GrepShell, path=path,
                            str_left=flank_5p, str_right=flank_3p, stop=37)
      counts/tot_kmers
    })
    R <- counts_outer/counts_outer[, "I"]
    R <- GeoMean(R[is.finite(R)])
    print(R)
    R
  })
  Rs
}


check_12 <- GrepCompetitorOligoKmers("miR-1")
print(check_12)


check_10

print(check_10)

break



GetGeoAverageCompKmers <- function(mirna, experiment="equilibrium",
                                   condition=40, n_ex=0, n_constant=5,
                                   kmer_len_start=4, kmer_len_end=12) {
  # Get the competitor oligo sequences, and the subregion for which the averages
  # are desired.
  competitor_seq <- SequenceList[["competitor"]][mirna]
  comp_seq <- substr(competitor_seq, 14, nchar(competitor_seq))
  # Make the empty lists to be added to during the function.
  fI <- rep(0, kmer_len_end - kmer_len_start + 1)
  fA <- rep(0, kmer_len_end - kmer_len_start + 1)
  Rs <- rep(0, kmer_len_end - kmer_len_start + 1)
  names(Rs) <- sprintf("k%s", kmer_len_start:kmer_len_end)
  names(fI) <- names(Rs)
  names(fA) <- names(Rs)
  # Use the "I" and particular condition to get the paths to the reads:
  conditions <- c("I", condition)
  path_reads <- sapply(conditions, GetAnalysisPath, mirna=mirna,
                         experiment=experiment, analysis_type="reads")
  # Get the number of reads for each condition.
  num_reads <- sapply(path_reads, function(path) {
    command <- sprintf("wc -l %s", path)
    as.numeric(strsplit(system(command, intern=TRUE), split=" ")[[1]][1])
  })
  # Loop through the different kmer lengths to populate the three output lists,
  # and also to get the geometric average enrichment in the subregion of the
  # competitor oligo that does not include the seed ('comp_seq_minus_seed').
  for (kmer_len in kmer_len_start:kmer_len_end) {
    kmer_paths <- sapply(conditions, GetPositionalKmersEquilPath, mirna=mirna,
                         experiment=experiment, kmer_len, n_constant=n_constant)
    mirna_seqs <- kMirnaSeqs[mirna]
    # Get all kmers for plotting purposes.
    competitor_kmers <- GetKmersFromString(comp_seq, len_k = kmer_len)
    rna_kmers <- sapply(competitor_kmers, RevComplement)
    rna_indeces <- sapply(rna_kmers, GetKmerIndex)
    # Get the input kmer frequencies.
    norm_I <- rowSums(t(sapply(rna_indeces, function(index) {
      command <- sprintf("head -%s %s | tail -1", index, kmer_paths[1])
      kmer_counts <- as.numeric(strsplit(system(command, intern=TRUE),
                                         split="\t")[[1]])
      kmer_counts/(num_reads[1]*(37 + 2*n_constant - kmer_len + 1))
    })))
    names(norm_I) <- competitor_kmers
    # Get the Ago-bound kmer frequencies.
    norm_A <- rowSums(t(sapply(rna_indeces, function(index) {
      command <- sprintf("head -%s %s | tail -1", index, kmer_paths[2])
      kmer_counts <- as.numeric(strsplit(system(command, intern=TRUE),
                                         split="\t")[[1]])
      kmer_counts/(num_reads[2]*(37 + 2*n_constant - kmer_len + 1))
    })))
    names(norm_A) <- competitor_kmers

    R <- norm_A/norm_I
    fI[sprintf("k%s", kmer_len)] <- GeoMean(norm_I)
    fA[sprintf("k%s", kmer_len)] <- GeoMean(norm_A)
    Rs[sprintf("k%s", kmer_len)] <- GeoMean(R)
  }
  print(fI)
  print(fA)
  print(Rs)
  print(rbind(fI, fA, Rs))
  out <- cbind(fI, fA, Rs)
  colnames(out) <- c("fI", "fA", "R")
  out
}




Rs_40_new <- GetGeoAverageCompKmers("miR-1")
Rs_12.6_new <- GetGeoAverageCompKmers("miR-1", condition="12.6")
Rs_4_new <- GetGeoAverageCompKmers("miR-1", condition="4")
Rs_1.26_new <- GetGeoAverageCompKmers("miR-1", condition="1.26")
Rs_0.4_new <- GetGeoAverageCompKmers("miR-1", condition="0.4")

dev.new(xpos=20, ypos=20, height=5, width=5)

plot(4:11, Rs_40[,"R"], type="o", ylim=c(0, 6), col=kEquilCondColors["40"])
lines(4:11, Rs_12.6[, "R"], type="o", col=kEquilCondColors["12.6"])
lines(4:11, Rs_4[, "R"], type="o", col=kEquilCondColors["4"])
lines(4:11, Rs_1.26[, "R"], type="o", col=kEquilCondColors["1.26"])
lines(4:11, Rs_0.4[, "R"], type="o", col=kEquilCondColors["0.4"])


break





GetAllCompetitorEnrichments <- function(mirna, experiment="equilibrium",
                                        condition=40, n_ex=0, n_constant=5,
                                        kmer_len_start=4, kmer_len_end=11,
                                        lib=FALSE, xpos=20, ypos=20,
                                        pdf.plot=FALSE) {
  # Get the competitor oligo sequences, and the subregion for which the averages
  # are desired.
  competitor_seq <- SequenceList[["competitor"]][mirna]
  comp_seq_minus_seed <- substr(competitor_seq, 14, nchar(competitor_seq))
  # Make the empty lists to be added to during the function.
  R_list <- list()
  norm_A_list <- list()
  norm_I_list <- list()
  # Use the "I" and particular condition to get the paths to the reads:
  conditions <- c("I", condition)
  path_reads <- sapply(conditions, GetAnalysisPath, mirna=mirna,
                         experiment=experiment, analysis_type="reads")
  # Get the number of reads for each condition.
  num_reads <- sapply(path_reads, function(path) {
    command <- sprintf("wc -l %s", path)
    as.numeric(strsplit(system(command, intern=TRUE), split=" ")[[1]][1])
  })
  # Loop through the different kmer lengths to populate the three output lists,
  # and also to get the geometric average enrichment in the subregion of the
  # competitor oligo that does not include the seed ('comp_seq_minus_seed').
  for (kmer_len in kmer_len_start:kmer_len_end) {
    kmer_paths <- sapply(conditions, GetPositionalKmersEquilPath, mirna=mirna,
                         experiment=experiment, kmer_len, n_constant=n_constant)
    mirna_seqs <- kMirnaSeqs[mirna]
    # Get all kmers for plotting purposes.
    competitor_kmers <- GetKmersFromString(competitor_seq, len_k = kmer_len)
    rna_kmers <- sapply(competitor_kmers, RevComplement)
    rna_indeces <- sapply(rna_kmers, GetKmerIndex)
    # Get the kmers not overlapping the seed for the geometric average.
    competitor_kmers_minus_seed <- GetKmersFromString(comp_seq_minus_seed,
                                                      len_k = kmer_len)    





    norm_I <- rowSums(t(sapply(rna_indeces, function(index) {
      command <- sprintf("head -%s %s | tail -1", index, kmer_paths[1])
      kmer_counts <- as.numeric(strsplit(system(command, intern=TRUE),
                                         split="\t")[[1]])
      kmer_counts/(num_reads[1]*(37 + 2*n_constant - kmer_len + 1))
    })))

    norm_A <- rowSums(t(sapply(rna_indeces, function(index) {
      command <- sprintf("head -%s %s | tail -1", index, kmer_paths[2])
      kmer_counts <- as.numeric(strsplit(system(command, intern=TRUE),
                                         split="\t")[[1]])
      kmer_counts/(num_reads[2]*(37 + 2*n_constant - kmer_len + 1))
    })))
    names(norm_I) <- competitor_kmers
    names(norm_A) <- competitor_kmers

    R <- norm_A/norm_I
    kmer_len <- as.character(kmer_len)
    R_list[[kmer_len]] <- R
    norm_A_list[[kmer_len]] <- norm_A
    norm_I_list[[kmer_len]] <- norm_I
  }
  cols <- rainbow(length(R_list))
  names(cols) <- names(R_list)


  dev.new(xpos=20, ypos=20, height=5, width=5)
  xmin <- 1
  xmax <- nchar(competitor_seq)
  ymin <- 0.1
  ymax <- 100
  par(kPlotParameters)
  BlankPlot(log='y')
  AddLogAxis(2, label="Enrichment")
  sapply(names(R_list), function(name_i) {
    output_i <- R_list[[name_i]]
    lines(seq(length(output_i)) + nchar(names(output_i)[1])/2,
          output_i, col=cols[name_i])
  })
  xy <- GetPlotFractionalCoords(0.5, 0.05, log='y')
  text(1:nchar(competitor_seq), xy[2], unlist(strsplit(competitor_seq, split="")),
     adj=c(0.5, 1), xpd=NA)

  dev.new(xpos=520, ypos=20, height=5, width=5)
  xmin <- 0
  xmax <- nchar(competitor_seq)
  ymin <- 1e-8
  ymax <- 0.1
  par(kPlotParameters)
  BlankPlot(log='y')
  AddLogAxis(2, label="Fraction")
  sapply(names(norm_A_list), function(name_i) {
    output_i <- norm_A_list[[name_i]]
    lines(seq(length(output_i)) + nchar(names(output_i)[1])/2,
          output_i, col=cols[name_i])
  })
  xy <- GetPlotFractionalCoords(0.5, 0.05, log='y')
  text(1:nchar(competitor_seq), xy[2], unlist(strsplit(competitor_seq, split="")),
     adj=c(0.5, 1), xpd=NA)

  dev.new(xpos=1020, ypos=20, height=5, width=5)
  xmin <- 0
  xmax <- nchar(competitor_seq)
  ymin <- 1e-8
  ymax <- 0.1
  par(kPlotParameters)
  BlankPlot(log='y')
  AddLogAxis(2, label="Fraction")
  sapply(names(norm_A_list), function(name_i) {
    output_i <- norm_I_list[[name_i]]
    lines(seq(length(output_i)) + nchar(names(output_i)[1])/2,
          output_i, col=cols[name_i])
  })
  xy <- GetPlotFractionalCoords(0.5, 0.05, log='y')
  text(1:nchar(competitor_seq), xy[2], unlist(strsplit(competitor_seq, split="")),
     adj=c(0.5, 1), xpd=NA)

  if (class(pdf.plot) == "character") {
    dev.off()
  }


}






  # kXp_I <- GetPositionalKmers("miR-1", experiment="equilibrium", n_constant=5,
  #                             condition="I", kmer_len=8, addcomp=FALSE)
  # kXp_1 <- GetPositionalKmers("miR-1", experiment="equilibrium", n_constant=5,
  #                             condition="0.4", kmer_len=8, addcomp=FALSE)
  # kXp_2 <- GetPositionalKmers("miR-1", experiment="equilibrium", n_constant=5,
  #                             condition="1.26", kmer_len=8, addcomp=FALSE)
  # kXp_3 <- GetPositionalKmers("miR-1", experiment="equilibrium", n_constant=5,
  #                             condition="4", kmer_len=8, addcomp=FALSE)
  # kXp_4 <- GetPositionalKmers("miR-1", experiment="equilibrium", n_constant=5,
  #                             condition="12.6", kmer_len=8, addcomp=FALSE)
  # kXp_5 <- GetPositionalKmers("miR-1", experiment="equilibrium", n_constant=5,
  #                             condition="40", kmer_len=8, addcomp=FALSE)


R_1 <- GetEnrichment(rowSums(kXp_1), rowSums(kXp_I))
R_2 <- GetEnrichment(rowSums(kXp_2), rowSums(kXp_I))
R_3 <- GetEnrichment(rowSums(kXp_3), rowSums(kXp_I))
R_4 <- GetEnrichment(rowSums(kXp_4), rowSums(kXp_I))
R_5 <- GetEnrichment(rowSums(kXp_5), rowSums(kXp_I))

full_enrichment_matrix <- cbind(R_1, R_2, R_3, R_4, R_5)

conc <- c(0.4, 1.26, 4, 12.6, 40)

graphics.off()
dev.new(xpos=20, ypos=20, height=5, width=5)
plot(conc, full_enrichment_matrix["ACCACACA",], log='xy')
dev.new(xpos=520, ypos=20, height=5, width=5)
plot(conc, full_enrichment_matrix["ACATTCCA",], log='xy')



PlotSiteOccupancyPerSite <- function(mirna, experiment="equilibrium", n_constant=5,
                              sitelist="paperfinal", plotlist=FALSE, singleonly=TRUE,
                              combined=TRUE, buffer=FALSE, uniq=FALSE,
                              trim=FALSE, bgoff=FALSE, col_name="12.65",
                              pdf.plot=FALSE, xpos=20, ypos=20) {
  print(col_name)
  sXc <- SubfunctionCall(SitesXCounts)
  pars.matrix <- SubfunctionCall(EquilPars)
  num_sites <- nrow(sXc)
  pars.model <- log10(pars.matrix$Mean)
  names(pars.model) <- rownames(pars.matrix)
  names(pars.model)[length(pars.model)] <- "AGO"
  names(pars.model)[length(pars.model) - 1] <- "bg"
  data <- GetDataEquil(sXc)
  data <<- data
  l <- SubfunctionCall(GetInputEquil)
  # Ago dilutio in the data:
  A.stock.measured <- kAgoStock[mirna, "equilibrium"]

  A.dil.data <- sapply(colnames(data), as.numeric)
  A.dil.model <- exp(seq(log(min(A.dil.data)),
                         log(max(A.dil.data)),
                     length=100))

  pM_from_dil <- 10*A.stock.measured
  xmin <- signif(min(A.dil.data)*pM_from_dil, 1)
  xmax <- signif(ceiling(max(A.dil.data)*pM_from_dil), 1)
  A.dil.model <- exp(seq(log(xmin), log(xmax), length=100))/pM_from_dil
  A.pM.data <- A.dil.data*pM_from_dil
  A.pM.model <-A.dil.model*pM_from_dil

  pars <- pars.model
  model <- SubfunctionCall(EquilSingleSiteModelPerSite, A.dil=A.dil.model, addbg=TRUE)

  model <- SubfunctionCall(EquilSingleSiteModelPerSite, A.dil=A.dil.model, addbg=TRUE)


  model.points <- SubfunctionCall(EquilSingleSiteModelPerSite, A.dil=A.dil.data, addbg=TRUE)

  totals <- SubfunctionCall(EquilBoundTotal, A.dil=A.dil.data, addbg=TRUE)

  # xmin <- signif(min(A.dil.data)*pM_from_dil, 1)
  # xmax <- signif(ceiling(max(A.dil.data)*pM_from_dil), 1)
  # ymin <- 1e-3
  # ymax <- 1
  xmin <- 0
  xmax <- 1
  ymin <- 0
  ymax <- 200


  SubfunctionCall(GenericFigureSaveFile)
  BlankPlot()
  cols = kSiteColors[rownames(model)]
  # sapply(rownames(data), function(site) {
  #     lines(A.pM.model, model[site, ], col=cols[site], xpd=NA)        
  # })
  AddLinearAxis(1, tick.space=0.05, label.space=0.2, label="% Occupancy",
                percent=TRUE)
  AddLinearAxis(2, tick.space=5, label.space=20, label="Enrichment")

  # dev.new(xpos=520, ypos=20, height=5, width=5)
  enrich_points <- GetEnrichment(sXc[, col_name], sXc[, 2])

  # row_bool <- which(model.points[, col_name] <= 0.1)


  # model.points_fit <- model.points[row_bool, col_name]
  # enrich_fit <- enrich_points[row_bool]
  # names(model.points_fit) <- rownames(model.points)[row_bool]

  # print(model.points_fit)

  # names(enrich_points) <- names(model.points_fit)

  # print(enrich_fit)

  offset <- 0
  enrich_fit <- enrich_points - offset
  lin_fit <- lm(enrich_fit ~ model.points[, col_name] - 1)

  # # print(A.pM.data)
  # # print(model.points["CACACACA", ])
  # # plot(A.pM.data, model.points["CACACACA", ])


  Points(model.points[, col_name], enrich_points, col=cols,
       xlim=c(0, 0.05), ylim=c(0, 5))
  abline(offset, lin_fit$coefficients)


  print(enrich_fit)
  print(100/totals)
  totals <<- totals
  xy <- GetPlotFractionalCoords(0.05, 0.95)
  text(xy[1], xy[2], label=mirna, adj=0)
  xy <- GetPlotFractionalCoords(0.05, 0.90)
  text(xy[1], xy[2], label=round(lin_fit$coefficients, 1), adj=0)
  xy <- GetPlotFractionalCoords(0.05, 0.85)
  text(xy[1], xy[2], label=round(100/colSums(totals)[col_name], 1), adj=0)
  # xy <- GetPlotFractionalCoords(0.05, 0.85)
  # text(xy[1], xy[2], label=round(lin_fit$coefficients[2], 1), adj=0)


  if (class(pdf.plot) == "character") {
    dev.off()
  }
}


date <- "190328"
mirna <- "miR-124"
for (col_name in c("0.4", "1.265", "4", "12.65", "40")) {
  print(mirna)
  print(col_name)
  pdf.plot <- sprintf("%s/%s_%s_Occ", date, mirna, col_name)
  print(pdf.plot)
  PlotSiteOccupancyPerSite(mirna, combined=TRUE, col_name=col_name,
                           buffer=FALSE, pdf.plot=pdf.plot)
}

break

PlotSiteOccupancyPerSite("miR-1", combined=FALSE, col_name="0.4", buffer=TRUE)
PlotSiteOccupancyPerSite("miR-1", combined=FALSE, col_name="1.265", buffer=TRUE, xpos=520)
PlotSiteOccupancyPerSite("miR-1", combined=FALSE, col_name="4", buffer=TRUE, xpos=1020)
PlotSiteOccupancyPerSite("miR-1", combined=FALSE, col_name="12.65", buffer=TRUE, ypos=520)
PlotSiteOccupancyPerSite("miR-1", combined=FALSE, col_name="40", buffer=TRUE, xpos=520, ypos=520)


break
PlotSiteOccupancyPerSite("miR", combined=TRUE, xpos=520)

PlotSiteOccupancyPerSite("miR-155", combined=TRUE, xpos=1020)

PlotSiteOccupancyPerSite("miR-124", combined=TRUE, ypos=520)

PlotSiteOccupancyPerSite("lsy-6", combined=TRUE, xpos=520, ypos=520)

PlotSiteOccupancyPerSite("miR-7-23nt", experiment="equilibrium2_nb",
                         combined=FALSE, xpos=1020, ypos=520)




break



break

CheckAddedCompetitorSequences <- function(mirna, experiment="equilibrium",
                                          condition=40, n_ex=0, n_constant=5,
                                          addcomp="1", kmer_len=8, height=5,
                                          width=5, xpos=20, ypos=20,
                                          pdf.plot=FALSE) {
  kXp_I <- SubfunctionCall(GetPositionalKmers, condition="I", addcomp=FALSE)
  kXp_A <- SubfunctionCall(GetPositionalKmers, addcomp=FALSE)
  kXp_A_c <- SubfunctionCall(GetPositionalKmers)
  R_A   <- GetEnrichment(rowSums(kXp_A), rowSums(kXp_I))
  R_A_c <- GetEnrichment(rowSums(kXp_A_c), rowSums(kXp_I))

  competitor_seq <- SequenceList[["competitor"]][mirna]
  competitor_kmers <- GetKmersFromString(competitor_seq, len_k = kmer_len)
  rna_kmers <- sapply(competitor_kmers, RevComplement)

  R_A_use <- R_A[rna_kmers]
  R_A_c_use <- R_A_c[rna_kmers]
  names(R_A_use) <- competitor_kmers
  names(R_A_c_use) <- competitor_kmers
  print(R_A_use)
  print(R_A_c_use)
  xmin <- 0
  xmax <- nchar(competitor_seq)
  ymin <- 0.1
  ymax <- 300
  SubfunctionCall(GenericFigureSaveFile)
  BlankPlot(log='y')
  AddLogAxis(2, label="Enrichment")
  x_vals <- 1:nchar(competitor_seq) + kmer_len/2 - 0.5
  print(x_vals)

  # segments(x_vals[1:length(enriched_binding)],
  #          0.04,
  #          x_vals[1:length(enriched_binding)],
  #          apply(cbind(enriched_binding, rel_Kd), 1, max), lty=2, xpd=NA)


  lines(x_vals[1:length(R_A_use)], R_A_use, lwd=1, type="o")

  lines(x_vals[1:length(R_A_use)], R_A_c_use, lwd=1, col="red", xpd=NA, type="o")
  # segments(x_vals[1:length(enriched_binding)],
  #          0.1,
  #          x_vals[1:length(enriched_binding)],
  #          rel_Kd, lwd=1, col="red")
  if (mirna == "miR-7-23nt") {
    mirna <- "miR-7"
  }
  xy <- GetPlotFractionalCoords(0.5, 1, log="y")
  legend(xy[1], xy[2], legend=c("No added reads", sprintf("+ %s%%\ncomp 11mers", addcomp)),
         col=c("black", "red"), bty="n", lwd=1)
  xy <- GetPlotFractionalCoords(0.05, 0.95, log="y")
  text(xy[1], xy[2], labels=mirna, adj=c(0, 1))
  xy <- GetPlotFractionalCoords(0.05, 0.90, log="y")
  text(xy[1], xy[2], labels=sprintf("%s-nt kmers", kmer_len),
       adj=c(0, 1))
  segments(xmin, 1, x1=xmax, lty=2)
  text(1:nchar(competitor_seq), 0.04, unlist(strsplit(competitor_seq, split="")),
       adj=c(0.5, 1), xpd=NA)

  if (class(pdf.plot) == "character") {
    dev.off()
  }


}


graphics.off()

date <- "190327"
mirna <- "miR-1"
addcomp <- "0.3"
for (kmer_len in seq(5, 10)) {
  print(mirna)
  pdf.plot <- sprintf("%s/%s_k%s_comp%s", date, mirna, kmer_len, addcomp)
  print(pdf.plot)
  CheckAddedCompetitorSequences(mirna, kmer_len=kmer_len, addcomp=addcomp,
                                pdf.plot=pdf.plot)
}


# CheckAddedCompetitorSequences("miR-1", kmer_len=6, xpos=420)
# CheckAddedCompetitorSequences("miR-1", kmer_len=7, xpos=820)
# CheckAddedCompetitorSequences("miR-1", kmer_len=8, ypos=420)
# CheckAddedCompetitorSequences("miR-1", kmer_len=9, ypos=420, xpos=420)
# CheckAddedCompetitorSequences("miR-1", kmer_len=10, ypos=420, xpos=820)




break


CompareCompDeltaGWithEnrichment <- function(mirna, experiment="equilibrium",
                                            n_ex=0, n_constant=5, kmer_len=8,
                                            lib=FALSE, xpos=20, ypos=20,
                                            pdf.plot=FALSE) {
  output_I <- SubfunctionCall(GetPositionalKmers, condition="I")
  output_A <- SubfunctionCall(GetPositionalKmers, condition="40")
  output_A <- output_A + SubfunctionCall(GetPositionalKmers, condition="12.6")
  output_R <- GetEnrichment(rowSums(output_A), rowSums(output_I))

  competitor_seq <- SequenceList[["competitor"]][mirna]
  print(competitor_seq)
  competitor_kmers <- GetKmersFromString(competitor_seq, len_k = kmer_len)
  rna_kmers <- sapply(competitor_kmers, RevComplement)
  print(competitor_kmers)

  enriched_binding <- output_R[rna_kmers]
  names(enriched_binding) <- competitor_kmers
  if (class(lib) == "character") {
    path <- file.path("CompetitorOligoMFEs",
                      sprintf("%s_k%s_lib%s.txt", mirna, kmer_len, lib))
  } else {
    path <- file.path("CompetitorOligoMFEs",
                      sprintf("%s_k%s.txt", mirna, kmer_len))

  }
  print(path)
  competitor_dGs <- read.table(path, header=FALSE, stringsAsFactors=FALSE)
  print(competitor_dGs)
  unique_dGs <- unique(competitor_dGs[, 2])
  print(unique_dGs)
  unique_names <- unique(competitor_dGs[, 1])
  print(unique_names)
  unique_dGs <- sapply(unique_names, function(unique_name) {
    competitor_dGs[which(competitor_dGs[, 1] == unique_name)[1], 2]
  })
  names(unique_dGs) <- unique_names
  print(competitor_dGs)
  # print(exp(competitor_dGs))
  print(unique_dGs[competitor_kmers])
  rel_Kd <- exp(mean(unique_dGs))/exp(unique_dGs)[competitor_kmers]
  print(rel_Kd)
  print(head(output_R))
  xmin <- 0
  xmax <- nchar(competitor_seq)
  ymin <- 0.1
  ymax <- 300
  SubfunctionCall(GenericFigureSaveFile)
  BlankPlot(log='y')
  AddLogAxis(2, label="Enrichment")
  x_vals <- 1:nchar(competitor_seq) + kmer_len/2 - 0.5
  print(x_vals)

  segments(x_vals[1:length(enriched_binding)],
           0.04,
           x_vals[1:length(enriched_binding)],
           apply(cbind(enriched_binding, rel_Kd), 1, max), lty=2, xpd=NA)


  lines(x_vals[1:length(enriched_binding)], enriched_binding, lwd=1, type="o")

  lines(x_vals[1:length(enriched_binding)], rel_Kd, lwd=1, col="red", xpd=NA, type="o")
  # segments(x_vals[1:length(enriched_binding)],
  #          0.1,
  #          x_vals[1:length(enriched_binding)],
  #          rel_Kd, lwd=1, col="red")
  if (mirna == "miR-7-23nt") {
    mirna <- "miR-7"
  }
  xy <- GetPlotFractionalCoords(0.5, 1, log="y")
  legend(xy[1], xy[2], legend=c("40% sample", "Ka"),
         col=c("black", "red"), bty="n", lwd=1)
  xy <- GetPlotFractionalCoords(0.05, 0.95, log="y")
  text(xy[1], xy[2], labels=mirna, adj=c(0, 1))
  text(1:nchar(competitor_seq), 0.04, unlist(strsplit(competitor_seq, split="")),
       adj=c(0.5, 1), xpd=NA)

  if (class(pdf.plot) == "character") {
    dev.off()
  }


}


break
# CompareCompDeltaGWithEnrichment("miR-1")
# # CompareCompDeltaGWithEnrichment("miR-1", lib="5p", xpos=520)
# # CompareCompDeltaGWithEnrichment("miR-1", lib="3p", xpos=1020)
# # break
# CompareCompDeltaGWithEnrichment("let-7a", xpos=520)
# CompareCompDeltaGWithEnrichment("miR-155", xpos=1020)
# CompareCompDeltaGWithEnrichment("miR-124", ypos=520)
# CompareCompDeltaGWithEnrichment("lsy-6", xpos=520, ypos=520)
# CompareCompDeltaGWithEnrichment("miR-7-23nt", experiment="equilibrium2_nb",
#                                 xpos=1020, ypos=520)






break



PlotOverlapPositions <- function(mirna, seq_type, condition,
                                 experiment="equilibrium", n_constant=5,
                                 kmer_len=8, n_ex=0, mean_norm=FALSE, height=4,
                                 width=4, xpos=20, ypos=20, pdf.plot=FALSE) {

  output_I <- SubfunctionCall(GetPositionalKmers, condition="I")
  output_A <- SubfunctionCall(GetPositionalKmers)
  output_A <- output_A + SubfunctionCall(GetPositionalKmers, condition="12.6")


  output_R <- GetEnrichment(output_A, output_I)


  kmers <- rownames(output_I)  

  kmers <- rownames(output_R)

  overlap_seq <- SequenceList[[seq_type]][mirna]
  overlap_seq_use <- RevComplement(overlap_seq)

  print(overlap_seq_use)
  range_ks <- nchar(overlap_seq_use) - kmer_len + 1
  xmin <- 0
  xmax <- nchar(overlap_seq_use)
  ymin <- -5
  ymax <- ncol(output_R)
  SubfunctionCall(GenericFigureSaveFile)
  BlankPlot()
  ymin <- 0
  AddLinearAxis(2, tick.space=1, label.space=5, label="Position")
  ymin <- -5
  overlap_seq_use_kmers <- sapply(1:range_ks, function(start) {
      substr(overlap_seq_use, start=start, stop=start + kmer_len - 1)
  })
  overlap_seq_use_kmer_inds <- sapply(overlap_seq_use_kmers, function(kmer) {
    which(rownames(output_R) == kmer)
  })
  R_inds <- output_R[rev(overlap_seq_use_kmer_inds),]
  R_inds <<- R_inds
  # image(as.matrix(log10(R_inds)))

  c_s <- 0.9 # color_start
  c_e <- 0.70  # color_end

  logR_ind <- log10(R_inds)
  kmers <- rownames(logR_ind)
  logR_ind <- do.call(data.frame,lapply(logR_ind, function(x) replace(x, is.infinite(x),NA)))
  rownames(logR_ind) <- kmers
  logR_ind <<- logR_ind

  if (mean_norm) {
    mean_R_vec <- apply(logR_ind, 1, mean, na.rm=TRUE)
    logR_ind <- logR_ind - mean_R_vec
    # max_R_vec <- apply(logR_ind, 1, max, na.rm=TRUE)
    # min_R_vec <- apply(logR_ind, 1, min, na.rm=TRUE)
    # print(logR_ind[1:6, 1:6])
    # print(min_R_vec[1:6])
    # print((logR_ind - min_R_vec)[1:6, 1:6])
    # R_transform <- (logR_ind - min_R_vec)/(max_R_vec - min_R_vec)*50 + 25
  }
  # } else {
    max_R <- max(logR_ind, na.rm=TRUE)
    min_R <- min(logR_ind, na.rm=TRUE)
    R_transform <- (logR_ind - min_R)/(max_R - min_R)*100

  # }

  R_transform <<- R_transform


  R_discrete <- round(R_transform)
  R_discrete <<- R_discrete

  color.dist = rev(rainbow(100, start=c_s, end=c_e))
  color.dist = c("gray", color.dist)
  col.inds <- sapply(unlist(R_discrete), function(col.ind) {
    min(max(1, col.ind), 100)
    })


  col.inds <<- col.inds
  col.inds[is.na(col.inds)] <- 0
  col.inds[which(col.inds == Inf)] <- 0
  col.inds <- col.inds + 1
  print(color.dist[col.inds[1:24]])

  print(col.inds[1:40])



  rect(xleft=rep(0:(nrow(R_inds) - 1) + 4, ncol(R_inds)),
       ybottom=rep(0:(ncol(R_inds) - 1), each=nrow(R_inds)),
       xright=rep(1:nrow(R_inds) + 4, ncol(R_inds)),
       ytop=rep(1:ncol(R_inds), each=nrow(R_inds)),
       col=color.dist[col.inds],
       border=NA)

  rect(9, -1, 10, 40, col=NA, lwd=0.5, xpd=NA)
  if (mirna == "miR-7-23nt") {
    mirna <- "miR-7"
  }
  xy <- GetPlotFractionalCoords(0.05, 0.95)

  text(xy[1], xy[2], labels=mirna, adj=0, xpd=NA)
  xy <- GetPlotFractionalCoords(0.5, 1)

  # legend(xy[1], xy[2], legend=paste0(conditions, "%"), col=kEquilCondColors[conditions],
  #        bty="n", lwd=1)
  text(1:nchar(overlap_seq), -2, unlist(strsplit(overlap_seq, split="")),
       adj=c(0.5, 0.5), xpd=NA)
  xy <- GetPlotFractionalCoords(0.5, 0)
  text(xy[1], xy[2], labels=SequenceListLabels[seq_type], adj=0.5, xpd=NA)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
} 



date <- "190315"
for (mirna in kMirnas) {
  print(mirna)
  if (mirna == "miR-7-23nt") experiment <- "equilibrium2_nb"
  else                       experiment <- "equilibrium"
  seq_type <- "mirna"
  for (seq_type in c("mirna", "capture", "lib5p", "lib3p", "star", "competitor")) {
  #   print(seq_type)
    if (seq_type == "capture")         width <- 5
    else if (seq_type == "competitor") width <- 4.5
    else                               width <- 4
    pdf.plot <- sprintf("%s/%s_%s", date, seq_type, mirna)
    print(pdf.plot)
    PlotOverlapPositions(mirna, seq_type=seq_type, condition=40,
                         experiment=experiment, width=width, pdf.plot=pdf.plot)

    pdf.plot <- sprintf("%s/%s_%s_meannorm", date, seq_type, mirna)
    print(pdf.plot)
    PlotOverlapPositions(mirna, seq_type=seq_type, condition=40,
                         experiment=experiment, mean_norm=TRUE, width=width,
                         pdf.plot=pdf.plot)


  }
}
# PlotOverlapIdentity("miR-1", seq_type="mirna", pdf.plot="190313/mirna_miR-1")
# PlotOverlapIdentity("miR-1", seq_type="star", xpos=420)
# PlotOverlapIdentity("miR-1", seq_type="capture", xpos=820)
# PlotOverlapIdentity("miR-1", seq_type="competitor", xpos=1220)
# PlotOverlapIdentity("miR-1", seq_type="lib5p", ypos=420)
# PlotOverlapIdentity("miR-1", seq_type="lib3p", ypos=420, xpos=420)


break
PlotCaptureIdentity("let-7a", lib_seq=T, xpos=520)
PlotCaptureIdentity("miR-155", lib_seq=T, xpos=1020)
PlotCaptureIdentity("miR-124", lib_seq=T, ypos=520)
PlotCaptureIdentity("lsy-6", lib_seq=T, ypos=520, xpos=520)
PlotCaptureIdentity("miR-7-23nt", experiment="equilibrium2_nb", lib_seq=T, xpos=1020, ypos=520)
                    # pdf.plot="190311/capture_miR-7-23nt")


break






Plot8merFlankVersusCapture <- function(mirna, experiment="equilibrium",
                                       site="6mer-m8", n_constant=5, kmer_len=8,
                                       combined=TRUE, buffer=FALSE,
                                       n_ex=0, rc=FALSE, mir_seq=FALSE,
                                       lib_seq = FALSE, pos="5p", pdf.plot=FALSE,
                                       height=5, width=5, xpos=20, ypos=280) {

  flank_kds <- SubfunctionCall(GetFlankKds)
  flank_kds <- flank_kds[grep("\\|", rownames(flank_kds), perl=TRUE), ]

  rownames(flank_kds) <- gsub("^.*\\|(.*)_Kd", rownames(flank_kds),
                              replace="\\1", perl=TRUE)
  print(flank_kds)



  kds_average <- 10^SubfunctionCall(GetAverageFlanks, sitelist="paperfinal",
                                    buffer=FALSE, combined=TRUE)


  flanks_use <- rownames(flank_kds)

  print(flank_kds)
  print(kds_average)


  cols <- rep("gray", nrow(flank_kds))
  inds_red <- grep("AC$", flanks_use)
  cols[inds_red] <- "red"
  inds_blue <- grep("^AC", flanks_use)
  cols[inds_blue] <- "blue"
  inds_purple <- grep("AC.AC", flanks_use)
  cols[inds_purple] <- "purple"

  print(length(kds_average[flanks_use]))
  print(length(flank_kds[, 2]))

  plot(flank_kds[, 2], kds_average[flanks_use], log='xy', col=cols)
  identify(flank_kds[, 2], kds_average[flanks_use], labels=flanks_use)


  kmers <- rownames(output)[2:nrow(output)]
  output <- as.numeric(unlist(output[2:nrow(output), 1, drop=FALSE]))
  output <- as.matrix(output, nrow=lenth(output), ncol=1)
  rownames(output) <- kmers
  colnames(output) <- "I"
  if (mirna == "miR-7-23nt") {
    conditions <- c("40", "12.6", "1.26", "0.4", "0")
  } else {
    conditions <- c("40", "12.6", "4", "1.26", "0.4", "0")    
  }
  for (condition in conditions) {
    output_i <- SubfunctionCall(KmersXCountsVector)
    output_i <- as.numeric(unlist(output_i[2:nrow(output_i), 1, drop=FALSE]))
    output_i <- as.matrix(output_i, nrow=lenth(output_i), ncol=1)
    rownames(output_i) <- kmers
    colnames(output_i) <- condition
    output <- cbind(output, output_i)
  }
  output <- MatNorm(output)
  output_R <- output/(output[, 1])
  capt_oligo <- kCaptureOligos[mirna]
  if (mir_seq) capt_oligo <- kMirnaSeqs[mirna]
  if (lib_seq) {
    capt_oligo <- kLibMirna[[mirna]]
    capt_oligo <- RevComplement(capt_oligo)
  }
  capt_oligo <- gsub("U", c(capt_oligo), replace="T")
  if (rc) capt_oligo <- RevComplement(capt_oligo)
  print(capt_oligo)
  break
  range_ks <- nchar(capt_oligo) - kmer_len + 1
  xmin <- 0
  xmax <- nchar(capt_oligo)
  ymin <- 0.1
  ymax <- 300
  SubfunctionCall(GenericFigureSaveFile)
  BlankPlot(log='y')
  AddLogAxis(2, label="Enrichment")
  capt_oligo_kmers <- sapply(1:range_ks, function(start) {
      substr(capt_oligo, start=start, stop=start + kmer_len - 1)
  })
  capt_oligo_kmer_inds <- sapply(capt_oligo_kmers, function(kmer) {
    which(rownames(output_R) == kmer)
  })
  R_inds <- output_R[capt_oligo_kmer_inds, -1]
  R_non_inds <- output_R[-c(capt_oligo_kmer_inds), -1]
  R_non_inds_average <- exp(colMeans(log(R_non_inds)))
  R_non_inds_SD <- apply(R_non_inds, 2, function(column) {
                          exp(sd(log(column)))
                         })
  R_inds_x <- 1:nrow(R_inds) + kmer_len/2 - 0.5
  R_noninds_x <- setdiff(1:nchar(capt_oligo) + 0.5, R_inds_x)


  sapply(rev(c(conditions)), function(condition) {
    ys <- R_inds[, condition]
    # plot(ecdf(log10(R_non_inds[, condition])))
    # plot(ecdf(log10(ys)), add=TRUE)
    # return()
    names(ys) <- R_inds_x
    ys_noninds <- rep(R_non_inds_average[condition], length(R_noninds_x))
    names(ys_noninds) <- R_noninds_x
    ys_all <- c(ys, ys_noninds)
    ys_all <- ys_all[order(as.numeric(names(ys_all)))]
    x_all <- as.numeric(names(ys_all))
    x_all <<- x_all
    ys_all <<- ys_all
    lines(x_all, ys_all, col=kEquilCondColors[condition])
    R_noninds_x_1 <- R_noninds_x[which(R_noninds_x < R_inds_x[1])]
    R_noninds_x_2 <- R_noninds_x[which(R_noninds_x > R_inds_x[1])]
    lines(R_noninds_x_1,
          rep(R_non_inds_average[condition]*(R_non_inds_SD[condition]^2),
              length(R_noninds_x_1)), lty=2, col=kEquilCondColors[condition])
    lines(R_noninds_x_1,
          rep(R_non_inds_average[condition]/(R_non_inds_SD[condition]^2),
              length(R_noninds_x_1)), lty=2, col=kEquilCondColors[condition])

    lines(R_noninds_x_2,
          rep(R_non_inds_average[condition]*(R_non_inds_SD[condition]^2),
              length(R_noninds_x_2)), lty=2, col=kEquilCondColors[condition])
    lines(R_noninds_x_2,
          rep(R_non_inds_average[condition]/(R_non_inds_SD[condition]^2),
              length(R_noninds_x_2)), lty=2, col=kEquilCondColors[condition])
    })
  if (mirna == "miR-7-23nt") {
    mirna <- "miR-7"
  }
  xy <- GetPlotFractionalCoords(1, 0.95, log='y')

  text(xy[1], xy[2], labels=mirna, adj=1, xpd=NA)
  xy <- GetPlotFractionalCoords(0.05, 1, log='y')

  legend(xy[1], xy[2], legend=paste0(conditions, "%"), col=kEquilCondColors[conditions],
         bty="n", lwd=1, adj=0)



  text(1:nchar(capt_oligo), 0.2, unlist(strsplit(capt_oligo, split="")),
       adj=c(0.5, 0.5), xpd=NA)
  segments(x0=27.5, y0=ymin, y1=ymax, lty=2)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
} 


# Plot8merFlankVersusCapture("miR-1", combined=FALSE, buffer=TRUE)


# break





PlotLuciferaseExperiment <- function(date, swap=TRUE) {
  path <- file.path("ReporterScreen", "LuciferaseExperiments",
                    sprintf("%s_luciferase_exp.txt", date))
  print(path)

  file <- as.data.frame(read.table(path, header=TRUE, row.names=1,
                                stringsAsFactors=FALSE, sep="\t"), stringsAsFactors=FALSE)
  file$Rel <- file$Renilla / file$Firefly
  print(file)
  if (date == "190227" & swap) {
    for (row_ind in 1:nrow(file)) {
      if (file$Site[row_ind] == "miR-1") {
        file$Sites[row_ind] <- "miR-124"
      } else {
        file$Sites[row_ind] <- "miR-1"
      }
    }
  }
  average_df <- matrix(c(sapply(unique(file[, 1]), function(ind_1) {
    sapply(unique(file[, 2]), function(ind_2) {
      sapply(unique(file[, 3]), function(ind_3) {
        bool_1 <- file[, 1] == ind_1
        bool_2 <- file[, 2] == ind_2
        bool_3 <- file[, 3] == ind_3
        ind_all <- which(bool_1 & bool_2 & bool_3)
        row_mean <- mean(file$Rel[ind_all], na.rm=TRUE)
        row_sd <- sd(file$Rel[ind_all], na.rm=TRUE)
        c(ind_1, ind_2, ind_3, row_mean, row_sd)
        })
      })
    })), byrow=TRUE, nrow=nrow(file)/3, ncol=5)
  average_df <- data.frame(average_df, stringsAsFactors=FALSE)
  average_df <<- average_df
  colnames(average_df) <- c(colnames(file)[1:3], "Mean", "SD")
  average_df$Mean <- as.numeric(average_df$Mean)
  average_df$SD <- as.numeric(average_df$SD)
  print(average_df)
  average_df <<- average_df

  norm_df <- data.frame(t(sapply(1:nrow(average_df), function(row_ind) {
    print(row_ind)
    if (date == "190224") {
      mirna <- average_df$miRNA[row_ind]
      utr <- average_df$UTR[row_ind]
      norm_ind <- which(average_df$miRNA == mirna &
                       average_df$UTR == utr &
                       average_df$Sites != mirna)
    } else if (date == "190227") {
      norm_ind <- ceiling(row_ind/2)*2 - 1
    } else if (date == "190228") {
      norm_ind <- ceiling(row_ind/2)*2
    }
    print(row_ind)
    print(row_ind/2)
    print(round(row_ind/2))
    print(norm_ind)
    mean_row <- average_df$Mean[row_ind]
    sd_row <- average_df$SD[row_ind]
    mean_ctrl <- average_df$Mean[norm_ind]
    sd_ctrl <- average_df$SD[norm_ind]
    mean_norm <- mean_row/mean_ctrl
    if (mean_norm == 1) {
      sd_norm <- 0
    } else {
      sd_norm <- mean_norm*sqrt((sd_row/mean_row)^2 + (sd_ctrl/mean_ctrl)^2)
    }
    out <- average_df[row_ind, ]
    out[4] <- mean_norm
    out[5] <- sd_norm
    out
  })))
  for (ind_col in 1:ncol(norm_df)) {
    norm_df[, ind_col] <- unlist(norm_df[, ind_col])
  }
  # norm_df$miRNA <- unlist(norm_df$miRNA)
  # norm_df$Sites <- unlist(norm_df$Sites)
  # norm_df$UTR <- unlist(norm_df$UTR)
  # norm_df$SD <- unlist(norm_df$SD)
  if (date == "190224") {
    norm_df$miRNA_UTR <- paste(norm_df$miRNA, norm_df$UTR, sep="_")
    norm_df$UTR_site <- paste(norm_df$UTR, ", ", norm_df$Sites, " sites", sep="")    
    list_temp <- list(norm_df$UTR_site, norm_df$miRNA)

  } else if (date %in% c("190227", "190228")) {
    norm_df$div <- paste(norm_df[, 1], "\n", norm_df[, 2], sep=" ")
  list_temp <- list(norm_df$Site, norm_df$div)
  }

  norm_df$Mean <- as.numeric(norm_df$Mean)
  tabbedMeans <- tapply(norm_df$Mean,
                        list_temp,
                        function(x) c(x = x))
  tabbedSDs <- tapply(norm_df$SD,
                        list_temp,
                        function(x) c(x = x))
  print(tabbedMeans)
  print(tabbedSDs)
  print(norm_df$Mean)
  print(norm_df$SD)
  if (date == "190224") {
    tabbedMeans <- tabbedMeans[c(2, 1, 4, 3), ]
    tabbedSDs <- tabbedSDs[c(2, 1, 4, 3), ]    
  } else if (date %in% c("190227", "190228")) {
    tabbedMeans <- tabbedMeans[, c(3, 4, 1, 2)]
    tabbedSDs <- tabbedSDs[, c(3, 4, 1, 2)]
  }

  print(tabbedMeans)


  print(norm_df)
  graphics.off()
  if (date == "190224") {
    legend_title <- "3' UTR:"
    names_arg <- NULL
  } else if (date %in% c("190227", "190228")) {
    legend_title <- "3' UTR sites:"
    names_arg <- c("", "", "", "")
  }
  dev.new(xpos=20, ypos=20, hight=5, width=5)
  barCenters <- barplot(height=tabbedMeans, beside=TRUE, border=NA, names.arg=names_arg,
                        ylim=c(0, 1.5), xlim=c(0, 12), legend.text=TRUE,
                        args.legend=list(title=legend_title, border=NA, box.lty=0))
  arrows(barCenters, tabbedMeans - tabbedSDs,
         barCenters, tabbedMeans + tabbedSDs, lwd=1.5, angle=90, code=3, length=0.05)
  if (date == "190224") {
    text(-1, -0.09, labels="miRNA:", xpd=NA)  
  } else if (date == "190227") {
    text(-1, -0.09, labels="TRL (ng):", xpd=NA)
    text(c(2, 5, 8, 11), -0.09, labels=c(0, 100, 0, 100), xpd=NA)
    segments(x0=c(1, 7), y0=-0.14, x1=c(6, 12), y1=-0.14, xpd=NA)
    text(-1, -0.19, labels="pUC (ug):", xpd=NA)
    text(c(3.5, 9.5), -0.19, labels=c(1, 0), xpd=NA)
  } else if (date == "190228") {
    text(-1, -0.09, labels="TRL (ng):", xpd=NA)
    text(c(2, 5, 8, 11), -0.09, labels=c(0, 100, 0, 100), xpd=NA)
    segments(x0=c(1, 7), y0=-0.14, x1=c(6, 12), y1=-0.14, xpd=NA)
    text(-1, -0.19, labels="miRNA (pmol):", xpd=NA)
    text(c(3.5, 9.5), -0.19, labels=c(50, 15), xpd=NA)
  }
  out_path <- file.path("ReporterScreen", "LuciferaseExperiments",
                    sprintf("%s_luciferase_exp.pdf", date))
  for (i in seq(0, 1) ) {
    x_inds <- seq(1, 3) + 6*i
    y_inds <- seq(4, 6) + 6*i
    x <- file$Rel[x_inds]
    y <- file$Rel[y_inds]
    print(x)
    print(y)
    print(class(x))
    print(class(y))
    print(x[1])
    alt_data <- data.frame(cbind(c("1", "1", "1", "2", "2", "2"), c(x, y)))
    colnames(alt_data) <- c("group", "Rel")
    alt_data$Rel <- as.numeric(alt_data$Rel)
    # print(wilcox.test(x, y=y, alternative="greater", paired=FALSE, exact=TRUE))
    # print(wilcox_test(Rel ~ group, data=alt_data, alternative="greater", distribution="exact"))
  }
  dev.copy2pdf(file=out_path)
}




out <- PlotLuciferaseExperiment("190224")
out <- PlotLuciferaseExperiment("190227")

out <- PlotLuciferaseExperiment("190228")

break




sXc <- SitesXCounts("miR-1", buffer=TRUE)

sXc_flanks <- SiteFlanksXCounts("miR-1", "8mer", experiment="equilibrium", n_constant=5,
                sitelist="paperfinal", buffer=TRUE)
sXc <- rbind(sXc_flanks, sXc[-1, ])


pars <- EquilPars("miR-1", sitelist="paperfinal", combined=FALSE, buffer=TRUE)

flank_pars <- GetFlankKds("miR-1", site="8mer", combined=TRUE, buffer=TRUE)
pars <- flank_pars





struct_8mer <- GetPairingFlankData("miR-1", "8mer", "I_combined", sitelist="paperfinal", buffer=TRUE)
struct_7merm8 <- GetPairingFlankData("miR-1", "7mer-m8", "I_combined", sitelist="paperfinal", buffer=TRUE, old=TRUE)


break
# struct_4 <- GetPairingFlankData("miR-1", "8mer", 4, sitelist="paperfinal", buffer=TRUE)
# struct_I <- GetPairingFlankData("miR-1", "8mer", "I", sitelist="paperfinal", buffer=TRUE)
# struct_I_combined <- GetPairingFlankData("miR-1", "8mer", "I_combined", sitelist="paperfinal", buffer=TRUE)


# dev.new(xpos=20, ypos=20, height=5, width=5)
dev.set(2)
plot(ecdf(struct_I$plfold), xlim=c(0, 1), ylim=c(0, 1))
# plot(ecdf(struct_4$plfold), add=TRUE, col="blue")
# plot(ecdf(struct_40$plfold), add=TRUE, col="red")
p <- seq(0, 1, length.out=100 + 1)

x_p <- (p[1:100] + p[2:(100 + 1)])/2

b_mean <- 0.69
b_var <- 0.027
C <- b_mean*(1 - b_mean)/b_var - 1
shape1_I <- b_mean*C
shape2_I <- (1 - b_mean)*C
lines(x_p, pbeta(x_p, shape1_I, shape2_I), lty=2)

b_mean <- 0.79
b_var <- 0.017
C <- b_mean*(1 - b_mean)/b_var - 1
shape1_4 <- b_mean*C
shape2_4 <- (1 - b_mean)*C

b_mean <- 0.75
b_var <- 0.022
C <- b_mean*(1 - b_mean)/b_var - 1
shape1_40 <- b_mean*C
shape2_40 <- (1 - b_mean)*C



pdf_4 <- dbeta(x_p, shape1_4, shape2_4)/10
pdf_40 <- dbeta(x_p, shape1_40, shape2_40)

pdf_I <- dbeta(x_p, shape1_I, shape2_I)

# lines(x_p, pbeta(x_p, shape1_4, shape2_4), lty=2, col="blue")
# lines(x_p, pbeta(x_p, shape1_40, shape2_40), lty=2, col="red")
lines(x_p, pdf_4, col="blue")
lines(x_p, pdf_40/max(pdf_40)*max(pdf_4), col="red")
lines(x_p, pdf_I/max(pdf_I)*max(pdf_4), col="black")

a_4 <- 0.1
a_40 <- 8*a_4
K <- 1

exp_ <- 2.8

selection <- dbeta(x_p, shape1_I, shape2_I)*a_4/(a_4 + K/(x_p^exp_))
selection <- selection/max(selection)*max(pdf_4)

lines(x_p, selection, lty=2, col="blue")

selection <- dbeta(x_p, shape1_I, shape2_I)*a_40/(a_40 + K/(x_p^exp_))
selection <- selection/max(selection)*max(pdf_4)

lines(x_p, selection, lty=2, col="red")


break




print(log10pars)

l <- Norm(sXc[, 1])*100


kds <- pars[1:nrow(sXc), 2]
names(kds) <- rownames(sXc)

bg <- pars[nrow(sXc) + 1, 2]

A <- pars[nrow(sXc) + 2, 2]

shape1 <- 7
shape2 <- 8

b_var <- c(rep(0.05, nrow(sXc) - 1), 0.1) 
b_mean <- c(0.4, rep(0.5, nrow(sXc) - 2), 0.4)

b_var <- rep(0.03, nrow(sXc))
b_mean <- rep(0.7, nrow(sXc))

pmin <- 1e-4
bins <- 20
# dev.new(xpos=20, ypos=520, height=5, width=5)
dev.set(3)
out <- BetaDist(b_mean, b_var, bins, plot.=TRUE)


A_data <- c(40, 12.6, 4, 1.26, 0.4)/100 * A
A_model <- exp(seq(log(0.01), log(3), length.out=100))


# x_test_0 <- BoundRNA(kds, l, A_data[1])

# x_test_1 <- BoundRNAWithStrucR(kds, l, A_data[1], b_mean[1], b_var[1], bins)
# print("done x_test_1")
# x_test_2 <- BoundRNAWithStrucR(kds, l, A_data[1], rep(b_mean[1], nrow(sXc)), rep(b_var[1], nrow(sXc)), bins)





x_norm <- sapply(A_model, function(A_j) {
	x <- BoundRNA(kds, l, A_j)
	f <- l - x
	x_bg <- Norm(f)*bg
	x + x_bg
})
# kds["8mer"] <- kds["8mer"]/2
x_norm_struc <- sapply(A_model, function(A_j) {
	x <- BoundRNAWithStrucR(kds, l, A_j, b_mean, b_var, bins)
	f <- l - x
	x_bg <- Norm(f)*bg
	x + x_bg
})


colnames(x_norm) <- round(A_model, 2)

R_data <- t(t(sXc[, 3:7])/colSums(sXc[, 3:7]))/Norm(l)

R_model <- t(t(x_norm)/colSums(x_norm))/Norm(l)

R_model_struc <- t(t(x_norm_struc)/colSums(x_norm_struc))/Norm(l)


# dev.new(xpos=520, ypos=20, height=5, width=5)
dev.set(4)
xmin <- min(A_model)
xmax <- max(A_model)
ymin <- 0.3
ymax <- 300
par(kPlotParameters)
BlankPlot(log='xy')
AddLogAxis(1, label="Conc")
AddLogAxis(2, label="Enrichment")

cols <- c(GetColorFunction(kFlanks), kSiteColors[rownames(sXc)[257:nrow(sXc)]])
sapply(1:nrow(sXc), function(row) {
	points(A_data, R_data[row, ], col=cols[row])
	lines(A_model, R_model[row, ], col=cols[row])

})

# dev.new(xpos=1020, ypos=20, height=5, width=5)
dev.set(5)
xmin <- min(A_model)
xmax <- max(A_model)
ymin <- 0.3
ymax <- 300
par(kPlotParameters)
BlankPlot(log='xy')
AddLogAxis(1, label="Conc")
AddLogAxis(2, label="Enrichment")

sapply(1:nrow(sXc), function(row) {
	points(A_data, R_data[row, ], col=cols[row])
	lines(A_model, R_model_struc[row, ], col=cols[row])

})


break
a_struc <- FreeAgoWithStructureR(kds, l, A, shape1, shape2)

print(a)
print(a_struc)

x_struc <- BoundRNAWithStrucR(kds, l, A, shape1, shape2)

print(x)
print(x_struc)
MakeEnrichments <- function(x_b) {
	return(Norm(x_b)/Norm(l))
}
R_un <- MakeEnrichments(x)
R_struc <- MakeEnrichments(x_struc)

plot(R_un, R_struc, log='xy', col=kSiteColors[rownames(sXc)])
abline(0, 1, lty=2)

break
# site_df <- GetRepressionMatrix("miR-1", "paperfinal", bg_method=3, new=TRUE, exrib=TRUE)

mir1_df <- GetRepressionMatrix("miR-1", "paperfinal", flanks=TRUE, bg_method=3, new=TRUE, exrib=TRUE)
let7_df <- GetRepressionMatrix("let-7a", "paperfinal", flanks=TRUE, bg_method=3, new=TRUE, exrib=TRUE)
mir155_df <- GetRepressionMatrix("miR-155", "paperfinal", flanks=TRUE, bg_method=3, new=TRUE, exrib=TRUE)
mir124_df <- GetRepressionMatrix("miR-124", "paperfinal", flanks=TRUE, bg_method=3, new=TRUE, exrib=TRUE)
lsy6_df <- GetRepressionMatrix("lsy-6", "paperfinal", flanks=TRUE, bg_method=3, new=TRUE, exrib=TRUE)
mir7_df <- GetRepressionMatrix("miR-7-23nt", "paperfinal", flanks=TRUE, bg_method=3, new=TRUE, exrib=TRUE)



# GetSiteCoeffieicnts <- function(pars, site_df, n_cutoff=20, kd_cutoff=0.1) {
#   # Set up parameters.
#   num_sites <- ncol(site_df) - 1
#   site_coefs <- pars[1:num_sites]
#   site_occurence <- site_df[, -ncol(site_df)]
#   names(site_coefs) <- colnames(site_occurence)
#   # names(flank_coefs) <- kFlanks
#   fc_model <- rowSums(t(site_coefs*t(site_occurence)))

#   # site_occurence <- site_occurence[, sites]
#   fc <- site_df[["fc"]]/log(2)
#   nosite_fc <- mean(fc[which(rowSums(site_occurence) == 0)])
#   fc <- fc - nosite_fc


#   if (check %% 10 == 0) {
#     plot(fc_model, site_df[["fc"]], xlim=c(2, -5), ylim=c(2, -5))  
#     print(site_coefs)
#   }
#   print(check)
#   print(check %% 100 == 0)
#   check <<- check + 1
#   out <- sum((fc - fc_model)^2)
#   print(out)
#   out
# }

# opt_1 <- optim(pars_sites, GetSiteCoeffieicnts, site_df=site_df, method="BFGS")
# opt_2 <- optim(opt_1$par, GetSiteCoeffieicnts, site_df=site_df)



pars_mir1 <- GetRepressionLinearModel("miR-1", buffer=TRUE, combined=FALSE, bg_method=3, new=TRUE, exrib=TRUE)
pars_let7 <- GetRepressionLinearModel("let-7a", bg_method=3, new=TRUE, exrib=TRUE)
pars_mir155 <- GetRepressionLinearModel("miR-155", bg_method=3, new=TRUE, exrib=TRUE)
pars_mir124 <- GetRepressionLinearModel("miR-124", bg_method=3, new=TRUE, exrib=TRUE)
pars_lsy6 <- GetRepressionLinearModel("lsy-6", bg_method=3, new=TRUE, exrib=TRUE)
pars_mir7 <- GetRepressionLinearModel("miR-7-23nt", experiment="equilibrium2_nb", combined=FALSE, bg_method=3, new=TRUE, exrib=TRUE)

lengths_pars <- c(nrow(pars_mir1), nrow(pars_let7), nrow(pars_mir155),
                  nrow(pars_mir124), nrow(pars_lsy6), nrow(pars_mir7))

print(lengths_pars)

pars_flank_sep <- rep(0, 256)
names(pars_flank_sep) <- kFlanks


flank_coefs_kd <- GetAverageFlanks(sitelist="paperfinal")

pars_flank <- c(pars_mir1[, 1], pars_let7[, 1], pars_mir155[, 1],
                pars_mir124[, 1], pars_lsy6[, 1], pars_mir7[, 1], flank_coefs_kd)


flank_cols <- GetColorFunction(kFlanks)

site_flank_df <- list(mir1_df, let7_df, mir155_df, mir124_df, lsy6_df, mir7_df)

dev.new(xpos=20, ypos=20, height=5, width=5)
dev.new(xpos=520, ypos=20, height=5, width=5)
dev.new(xpos=1020, ypos=20, height=5, width=5)

GetSiteAndFlankCoeffieicnts <- function(pars, site_flank_df, lengths_pars, plotbool=TRUE) {
  # Set up parameters.

  cost <- 0
  start <- 1
  flank_coefs <- pars[(sum(lengths_pars) + 1):(sum(lengths_pars) + 256)]
  fc_all <- c()
  fc_model_all <- c()
  for (i in seq(length(lengths_pars))) {
    num_sites <- lengths_pars[i]
    site_coefs <- pars[start:(start + num_sites - 1)]
    site_occurence <- site_flank_df[[i]][, 1:num_sites]
    flank_occurence <- site_flank_df[[i]][(num_sites + 1):(num_sites + 256)]
    fc_model_sites <- rowSums(t(site_coefs*t(site_occurence)))
    fc_model_flanks <- rowSums(t(flank_coefs*t(flank_occurence)))
    fc_model <- fc_model_sites + fc_model_flanks
    fc <- site_flank_df[[i]][["fc"]]/log(2)
    nosite_fc <- mean(fc[which(rowSums(site_occurence) == 0)])
    fc <- fc - nosite_fc
    fc_all <- c(fc_all, fc)
    fc_model_all <- c(fc_model_all, fc_model)
    start <- start + num_sites
  }
  cost <- sum((fc_all - fc_model_all)^2)
  if (check %% 2 == 0 & plotbool) {
    dev.set(2)
    plot(fc_model_all, fc_all, xlim=c(2, -5), ylim=c(2, -5))
    dev.set(3)

    plot(flank_coefs, flank_coefs_kd, col=flank_cols)
    print(cost)
    print(length(fc_model_all))
    print(length(fc_all))
    dev.set(4)
    plot(pars_flank[1:sum(lengths_pars)], pars[1:sum(lengths_pars)])

  }
  check <<- check + 1
  cost
}

GradFlankFunction <- function(pars, site_flank_df, lengths_pars) {
  grad(GetSiteAndFlankCoeffieicnts, pars, method="simple",
       site_flank_df=site_flank_df, lengths_pars=lengths_pars, plotbool=FALSE)
}

opt <- optim(pars_flank, GetSiteAndFlankCoeffieicnts, gr=GradFlankFunction,
             site_flank_df=site_flank_df, lengths_pars=lengths_pars, method="BFGS")



break







sXc_nb <- read.table(file.path("/lab/solexa_bartel/nbisaria/analysis/miR-1",
                              "equilibrium/seandata",
                              "site_counts_streamline3p_combined/Counts.txt"))

colnames(sXc_nb) <- gsub("^X(.*)$", replace="\\1", colnames(sXc_nb), perl=TRUE)
colnames(sXc_nb) <- gsub("^(.*)6$", replace="\\165", colnames(sXc_nb),
                         perl=TRUE)

sXc_nb <- sXc_nb[, c(7, 1:6)]

sXc <- SitesXCounts("miR-1", buffer=TRUE)
sXc <- sXc[, -1]
colnames(sXc)[1] <- "I"

print(sXc)


print(head(sXc_nb))
dev.new(xpos=20, ypos=20, height=5, width=5)
plot(sXc_nb[, 1], sXc_nb[, 7])


dev.new(xpos=520, ypos=20, height=5, width=5)
plot(sXc_nb[, 1], sXc_nb[, 7], log='xy')

site_8mer <- GetSiteSeq("miR-1", "8mer")
print(site_8mer)

inds_8mer <- grep(paste0("-", site_8mer), rownames(sXc_nb))
print(sXc_nb[inds_8mer, ])

x <- c(colSums(sXc_nb[inds_8mer, ]) + c(unlist(sXc_nb["8mer-1.8",])))
y <- c(unlist(sXc["8mer", ]))

df_lm <- data.frame(x=x, y=y)
plot(x, y, xlim=c(0, 1e6),
     ylim=c(0, 1e6))
abline(0, 1, lty=2)

lm_coefs <- lm(x ~ y, data=df_lm)

print(lm_coefs$coef)

dev.new(xpos=1020, pos=20, height=5, width=5)
plot(c(unlist(sXc_nb["None", ])), c(unlist(sXc["None", ])))
abline(0, 1)


break



# PlotEquilSiteWithInput(sXc, sXc_ft, 1, 1)

# PlotEquilSiteWithInput(sXc, sXc_ft, 1, 2, xpos=520)

# PlotEquilSiteWithInput(sXc, sXc_ft, 1, 3, ypos=520)

# PlotEquilSiteWithInput(sXc, sXc_ft, 1, 4, xpos=520, ypos=520)

PlotEquilSiteWithInput(sXc, sXc_ft, 1, 5, ymin=1e-1, ymax=2, ft_x=TRUE,
                       xpos=1020)


PlotEquilSiteWithInput(sXc, sXc_ft, 1, 6, ymin=1e-1, ymax=2, ft_x=TRUE,
                       ft_axes=TRUE, xpos=1020)


PlotEquilSiteWithInput(sXc, sXc_ft, 1, 7, ymin=1e-1, ymax=2, ft_x=TRUE,
                       ft_axes=TRUE, xpos=1020)

PlotEquilSiteWithInput(sXc, sXc_ft, 1, 8, ymin=1e-1, ymax=2, ft_x=TRUE,
                       ft_axes=TRUE, xpos=1020)

PlotEquilSiteWithInput(sXc, sXc_ft, 1, 9, ymin=1e-1, ymax=2, ft_x=TRUE,
                       ft_axes=TRUE, xpos=1020)

PlotEquilSiteWithInput(sXc, sXc_ft, 1, 10, ymin=1e-1, ymax=2, ft_x=TRUE,
                       ft_axes=TRUE, xpos=1020)


PlotEquilSiteWithInput(sXc_124, sXc_124_ft, 1, 1, ymin=1e-1, ymax=2,
                       ft_axes=TRUE, xpos=1020, ypos=520)

PlotEquilSiteWithInput(sXc_124, sXc_124_ft, 1, 2, ymin=1e-1, ymax=2,
                       ft_axes=TRUE, xpos=1020, ypos=520)

PlotEquilSiteWithInput(sXc_124, sXc_124_ft, 1, 3, ymin=1e-1, ymax=2,
                       ft_axes=TRUE, xpos=1020, ypos=520)

PlotEquilSiteWithInput(sXc_124, sXc_124_ft, 1, 4, ymin=1e-1, ymax=2,
                       ft_axes=TRUE, xpos=1020, ypos=520)

PlotEquilSiteWithInput(sXc_124, sXc_124_ft, 1, 5, ymin=1e-1, ymax=2,
                       ft_axes=TRUE, xpos=1020, ypos=520)

PlotEquilSiteWithInput(sXc_124, sXc_124_ft, 1, 6, ymin=1e-1, ymax=2,
                       ft_axes=TRUE, xpos=1020, ypos=520)




break

Positional12NtKmerList <- function(mirna, pos) {
  mirseq <- MirnaTargetSequence(mirna, pos, pos + 3)
  print(mirseq)
  sapply(GetKmerList(8), function(kmer) {
    pre.str <- substr(kmer, 1, 6 - pos + 1)
    post.str <- substr(kmer, 7 - pos + 1, 8)
    paste0(pre.str, mirseq, post.str)
  })
}




# kmers_1 <- Positional12NtKmerList("miR-1", 1)
# kmers_2 <- Positional12NtKmerList("miR-1", 2)
# kmers_3 <- Positional12NtKmerList("miR-1", 3)
# kmers_4 <- Positional12NtKmerList("miR-1", 4)
# kmers_5 <- Positional12NtKmerList("miR-1", 5)

# kmers_unique <- unique(c(kmers_1, kmers_2, kmers_3, kmers_4, kmers_5))

# print(length(kmers_unique))


rbns_data <- read.table(file.path("/lab/bartel4_ata/kathyl/RNA_Seq/analysis",
                                  "data/no_baseline_analysis/pred_dfs",
                                  "pred_df.txt"), sep="\t", header=TRUE,
                                  stringsAsFactors=FALSE)
mirnas_rbns <- unique(rbns_data$mir)
names(mirnas_rbns) <- gsub("mir", mirnas_rbns, replace="miR-")
names(mirnas_rbns) <- gsub("let7", names(mirnas_rbns), replace="let-7a")
names(mirnas_rbns) <- gsub("lsy6", names(mirnas_rbns), replace="lsy-6")

rbns_canon_data <- read.table(file.path("/lab/bartel4_ata/kathyl/RNA_Seq/analysis",
                                  "data/no_baseline_analysis/pred_dfs",
                                  "pred_df_canon.txt"), sep="\t", header=TRUE,
                                  stringsAsFactors=FALSE)

cnn_data <- read.table(file.path("/lab/bartel4_ata/kathyl/RNA_Seq/analysis",
                                  "data/no_baseline_analysis",
                                  "convnet_preds.txt"), sep="\t", header=TRUE,
                                  stringsAsFactors=FALSE)
colnames(cnn_data)[1] <- "gene"
mirnas_cnn <- unique(cnn_data$mir)
names(mirnas_cnn) <- gsub("mir", mirnas_cnn, replace="miR-")
names(mirnas_cnn) <- gsub("let7", names(mirnas_cnn), replace="let-7a")
names(mirnas_cnn) <- gsub("lsy6", names(mirnas_cnn), replace="lsy-6")

UTR_data <- read.table(file.path("/lab/bartel4_ata/kathyl/RNA_Seq/analysis",
                                  "data/no_baseline_analysis",
                                  "merged.txt"), sep="\t", header=TRUE,
                                  stringsAsFactors=FALSE)


path_ts7 <- file.path("/lab/bartel4_ata/kathyl/RNA_Seq/analysis/data",
                     "no_baseline_analysis/ts7_preds.txt")
ts_data <- read.table(path_ts7, sep="\t", header=TRUE)




col_list <- list("0"=c(red=253.223415, green=196.7886, blue=39.47604)/256,
                 "1"=c(red=249.704415, green=154.92066, blue=60.693315)/256,
                 "2"=c(red=235.438185, green=118.129005, blue=84.864255)/256,
                 "3"=c(red=214.239525, green=85.0629, blue=109.001025)/256,
                 "4"=c(red=187.684845, green=53.406945, blue=134.61654)/256,
                 "5"=c(red=154.398675, green=21.89277, blue=158.78493)/256,
                 "6"=c(red=114.16707, green=0.5304, blue=168.3612)/256,
                 "7"=c(red=69.918705, green=3.087795, blue=158.79411)/256,
                 "8"=c(red=12.847665, green=7.599765, blue=134.633625)/256)
col_cols <- sapply(col_list, function(col_i) {
  rgb(red=col_i["red"], green=col_i["green"], blue=col_i["blue"], alpha=0.8)
  })

PlotBiochemicalModel <- function(mirna, canon=FALSE, height=4, width=4, xpos=20,
                                 ypos=20, pdf.plot=FALSE) {
  print(mirna)
  print(canon)
  xmin <- 0
  ymax <- 2.1

  if (mirna %in% names(mirnas_rbns)) {
    mir_ <- mirnas_rbns[mirna]
    if (canon) {
      rep_data <- rbns_canon_data
    } else {
      rep_data <- rbns_data
    }
    xmax <- 2
    ymin <- -5
  } else if (mirna == "CNN") {
    height <-5
    width <- 5
    rep_data <- cnn_data
    xmax <- 1.5
    ymin <- -4
  } else if (mirna %in% names(mirnas_cnn)) {
    mir_ <- mirnas_cnn[mirna]
    rep_data <- cnn_data
    xmax <- 1.5
    ymin <- -4
  } else if (mirna == "RBNS") {
    height <-5
    width <- 5

    if (canon) {
      rep_data <- rbns_canon_data
    } else {
      rep_data <- rbns_data
    }
    xmax <- 2
    ymin <- -5    
  }
  if (mirna %in% c("CNN", "RBNS")) {
    rep_mirna <- rep_data[,c("gene", "num_canon", "pred", "log2fc_bayes")]
    rep_mirna <- rep_mirna[order(-rep_mirna$num_canon),]
    xy_cor <- GetPlotFractionalCoords(fx=0.80, fy=0.95)
    xy_title <- GetPlotFractionalCoords(fx=0.5, fy=1)
    adj_title <- 0.5
    if (mirna == "RBNS") {
      if (canon) {
        mirna <- "Only 12-nt k-mers with canonical sites"
      } else {
        mirna <- "All six miRNAs"      
      }
    } else {
      mirna <- "All eleven"
    }
  } else {
    rep_mirna <- subset(rep_data, mir==mir_,
                        select=c(gene, num_canon, pred, log2fc_bayes))
    rep_mirna <- rep_mirna[order(-rep_mirna$num_canon),]    
    xy_cor <- GetPlotFractionalCoords(fx=0.80, fy=1)
    xy_title <- GetPlotFractionalCoords(fx=0.05, fy=1)
    adj_title <- 0

  }
  SubfunctionCall(GenericFigureSaveFile)
  BlankPlot()
  ymax <- 2
  cols_i <- rep_mirna$num_canon
  cols_i[which(cols_i >= 8)] <- 8
  cols_i <- as.character(cols_i)
  cols <- col_cols[cols_i]
  points(-rep_mirna$pred, rep_mirna$log2fc_bayes, col=cols, pch=20)
  # xy_cor <- GetPlotFractionalCoords(fx=0.80, fy=1)
  AddCorrelationToPlot(-rep_mirna$pred, rep_mirna$log2fc_bayes,
                       xpos=xy_cor[1], ypos=xy_cor[2], rsquared=TRUE)
  # xy_title <- GetPlotFractionalCoords(fx=0.05, fy=1)
  # if (canon) {
  text(xy_title[1], xy_title[2], labels=mirna, xpd=NA, adj=adj_title)
  # } else {
  #   text(xy[1], xy[2], labels=mirna, xpd=NA, adj=0)  
  # }
  AddLinearAxis(1, label.space=0.5, tick.space=0.5,
                label="Predicted repression (log2)")
  AddLinearAxis(2, label.space=1, tick.space=1,
                label="Measured fold-change (log2")
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

PlotTargetScan <- function(mirna, height=5, width=7, xpos=20, ypos=20,
                           xbayes=FALSE, ybayes=FALSE, pdf.plot=FALSE) {
  path <- file.path("/lab/bartel4_ata/kathyl/RNA_Seq/analysis/data",
                    "no_baseline_analysis/ts7_preds.txt")
  data <- read.table(path, sep="\t", header=TRUE)
  mirnas <- unique(data$mir)
  print(mirnas)
  print(head(data))
  if (mirna %in% names(mirnas_rbns)) {
    mir_ <- mirnas_rbns[mirna]
  } else if (mirna %in% names(mirnas_cnn)) {
    mir_ <- mirnas_cnn[mirna]
  }
  if (mirna %in% c("CNN", "RBNS")) {
    rep_mirna <- subset(data, mir %in% mirnas_rbns, 
                        select=c(Gene_ID, num_canon, pred_bayes, pred_method1,
                                 log2fc_bayes, log2fc_method1)
                        )
    rep_mirna <- rep_mirna[order(-rep_mirna$num_canon),]
    if (mirna == "RBNS") {
      mirna <- "All six, TargetScan7"
    } else {
      mirna <- "All eleven, TargetScan7"
    }
  } else {
    rep_mirna <- subset(data, mir==mir_,
                        select=c(Gene_ID, num_canon, pred_bayes, pred_method1,
                                 log2fc_bayes, log2fc_method1)
                        )
    rep_mirna <- rep_mirna[order(-rep_mirna$num_canon),]    
  }
  SubfunctionCall(GenericFigureSaveFile)
  xmin <- 0
  xmax <- 6
  ymin <- -5
  ymax <- 2.1
  BlankPlot()
  ymax <- 2
  cols_i <- rep_mirna$num_canon
  cols_i[which(cols_i >= 8)] <- 8
  cols_i <- as.character(cols_i)
  cols <- col_cols[cols_i]
  print(length(cols))
  if (xbayes) {
    x <- -rep_mirna$pred_bayes
  } else {
    x <- -rep_mirna$pred_method1
  }
  if (ybayes) {
    y <- rep_mirna$log2fc_bayes
  } else {
    y <- rep_mirna$log2fc_method1
  }
  print(length(x))
  print(length(y))
  # break
  points(x, y, col=cols, pch=20)
  xy <- GetPlotFractionalCoords(fx=0.80, fy=0.95)
  AddCorrelationToPlot(x=x, y=y,
                       xpos=xy[1], ypos=xy[2], rsquared=TRUE)
  xy <- GetPlotFractionalCoords(fx=0.5, fy=1)
  text(xy[1], xy[2], labels=mirna, xpd=NA, adj=0.5)
  AddLinearAxis(1, label.space=1, tick.space=1,
                label="Predicted repression (log2)")
  AddLinearAxis(2, label.space=1, tick.space=1,
                label="Measured fold-change (log2")
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}



# PlotBiochemicalModel(  "miR-1")
# PlotBiochemicalModel( "let-7a", xpos= 420)
# PlotBiochemicalModel("miR-155", xpos= 820)
# PlotBiochemicalModel("miR-124",            ypos=420)
# PlotBiochemicalModel(  "lsy-6", xpos= 420, ypos=420)
# PlotBiochemicalModel  ("miR-7", xpos= 820, ypos=420)
# PlotBiochemicalModel(   "RBNS", xpos=1220, ypos=420)

# PlotBiochemicalModel(  "miR-1", pdf.plot="180927_miR-1_rbns_rep_model")
# PlotBiochemicalModel( "let-7a", pdf.plot="180927_let-7a_rbns_rep_model")
# PlotBiochemicalModel("miR-155", pdf.plot="180927_miR-155_rbns_rep_model")
# PlotBiochemicalModel("miR-124", pdf.plot="180927_miR-124_rbns_rep_model")
# PlotBiochemicalModel(  "lsy-6", pdf.plot="180927_lsy-6_rbns_rep_model")
# PlotBiochemicalModel  ("miR-7", pdf.plot="180927_miR-7_rbns_rep_model")
PlotBiochemicalModel(   "RBNS", pdf.plot="180927_all_rbns_rep_model")
PlotBiochemicalModel(   "RBNS", canon=TRUE,
                     pdf.plot="180927_all_rbns_canon_rep_model")


# PlotTargetScan("RBNS")
# PlotTargetScan("RBNS", xbayes=TRUE)
# PlotTargetScan("RBNS", ybayes=TRUE)
PlotTargetScan("RBNS", xbayes=TRUE, ybayes=TRUE, pdf.plot="180927_all_rbns_ts_rep_model")


# PlotBiochemicalModel("miR-137")                        # 1
# PlotBiochemicalModel( "miR-139", xpos= 420)            # 2
# PlotBiochemicalModel( "miR-143", xpos= 820)            # 3
# PlotBiochemicalModel( "miR-144", xpos=1220)            # 4
# PlotBiochemicalModel( "miR-153",            ypos=420)  # 5
# PlotBiochemicalModel( "miR-182", xpos= 420, ypos=420)  # 6
# PlotBiochemicalModel("miR-199a", xpos= 820, ypos=420)  # 7
# PlotBiochemicalModel( "miR-204", xpos=1220, ypos=420)  # 8
# PlotBiochemicalModel( "miR-205",          , ypos=820)  # 9
# PlotBiochemicalModel("miR-216b", xpos= 420, ypos=820)  #10
# PlotBiochemicalModel( "miR-223", xpos= 820, ypos=820)  #11
# PlotBiochemicalModel(     "CNN", xpos=1220, ypos=820)  #11


PlotBiochemicalModel( "miR-137", pdf.plot="180927_miR-137_cnn_rep_model")
PlotBiochemicalModel( "miR-139", pdf.plot="180927_miR-139_cnn_rep_model")
PlotBiochemicalModel( "miR-143", pdf.plot="180927_miR-143_cnn_rep_model")
PlotBiochemicalModel( "miR-144", pdf.plot="180927_miR-144_cnn_rep_model")
PlotBiochemicalModel( "miR-153", pdf.plot="180927_miR-153_cnn_rep_model")
PlotBiochemicalModel( "miR-182", pdf.plot="180927_miR-182_cnn_rep_model")
PlotBiochemicalModel("miR-199a", pdf.plot="180927_miR-199a_cnn_rep_model")
PlotBiochemicalModel( "miR-204", pdf.plot="180927_miR-204_cnn_rep_model")
PlotBiochemicalModel( "miR-205", pdf.plot="180927_miR-1205_cnn_rep_model")
PlotBiochemicalModel("miR-216b", pdf.plot="180927_miR-216b_cnn_rep_model")
PlotBiochemicalModel( "miR-223", pdf.plot="180927_miR-223_cnn_rep_model")
# PlotBiochemicalModel(     "CNN", pdf.plot="180927_miR-1_rbns_rep_model")




break



Check_mean_and_max_error <- function(mirna="miR-1", experiment="equilibrium",
                                     sitelist="paperfinal", site="8mer", site2="5mer-m2.6",
                                     combined=TRUE, buffer=TRUE) {
  flanks_8mer <- SubfunctionCall(GetFullMirnaSiteFlankKds)
  flanks_8mer_global <<- flanks_8mer
  flanks_8mer_CI <- SubfunctionCall(GetFullMirnaSiteFlankKds_CI)
  flanks_8merxC5 <- SubfunctionCall(GetFullMirnaSiteFlankKds, site=site2)
  flanks_8merxC5_global <<- flanks_8merxC5
  flanks_8mer_CI_global <<- flanks_8mer_CI
  flanks_8merxC5_CI <- SubfunctionCall(GetFullMirnaSiteFlankKds_CI,
                                       site=site2)


  dev.new(xpos=20, ypos=20, height=5, width=5)
  par(kPlotParameters)
  xmin <- 1e-5
  xmax <- 10
  ymin <- xmin
  ymax <- xmax
  BlankPlot(log='xy', inv='xy')
  arrows(flanks_8mer, flanks_8mer_CI[, 1], flanks_8mer, flanks_8mer_CI[, 2],
         length=0.02, lty=2, angle=90, code=3, xpd=NA, col=kSiteColors[site])

  message("8mer")
  message(max(flanks_8mer_CI[, 1]/flanks_8mer_CI[, 2]))
  message(exp(mean(log(flanks_8mer_CI[, 1]/flanks_8mer_CI[, 2]), na.rm=TRUE)))
  message(site2)
  message(max(flanks_8merxC5_CI[, 1]/flanks_8merxC5_CI[, 2], na.rm=TRUE))
  message(exp(mean(log(flanks_8merxC5_CI[, 1]/flanks_8merxC5_CI[, 2]), na.rm=TRUE)))


  arrows(flanks_8merxC5, flanks_8merxC5_CI[, 1], flanks_8merxC5, flanks_8merxC5_CI[, 2],
         length=0.02, lty=2, angle=90, code=3, xpd=NA, col=kSiteColors[site2])

  points(flanks_8mer, flanks_8mer, col=GetColorFunction(names(flanks_8mer)),
         cex=0.1)
  AddLogAxis(1, label="8mer")
  AddLogAxis(2, label="8mer-xC5")  
}

Check_mean_and_max_error()
break



checkm1 <- GetFlankLinearModel(leaveout="miR-1")
checkl7 <- GetFlankLinearModel(leaveout="let-7a")

GetAllFlanksTest <- function(lin_mat, ints=NULL) {
  flanks <- kFlanks
  sapply(kFlanks, function(flank) {
    nucs <- sapply(c(1, 2, 4, 5), function(ind) {
      nuc <- substr(flank, ind, ind)
      print(nuc)
      print(rownames(lin_mat))
      out_ind <- grep(nuc, rownames(lin_mat))
      print(out_ind)
      out_ind
    })
    pairs <- cbind(nucs, 1:4)

    lin_part <- sum(apply(pairs, 1, function(row) {
      lin_mat[row[1], row[2]]
    }))
    if (class(ints) == "list") {
      # print(nucs[1])
      # print(ints[[1]])
    nucs <- sapply(c(1, 2, 4, 5), function(ind) {
      substr(flank, ind, ind)
    })

      if (flank == "CC.CC") {
        print(ints[[1]][nucs[2], nucs[1]])
        print(ints[[2]][nucs[4], nucs[3]])
      }
      int_5psum <- ints[[1]][nucs[2], nucs[1]]
      int_3psum <- ints[[2]][nucs[4], nucs[3]]
      print(flank)
      print(nucs)
      print("in test")
      print(int_5psum)
      print(int_3psum)
      lin_part <- lin_part + int_5psum + int_3psum
    }
    lin_part
  })
}
graphics.off()


out1 <- GetAllFlanksTest(lins_global)
out2 <- GetAllFlanksTest(lins_post_global)

out3 <- GetAllFlanksTest(lins_global, ints=list(int5p, int3p))
out4 <- GetAllFlanksTest(lins_post_global, ints=list(int5p_norm, int3p_norm))

temp_flanks <- GetAverageFlanks(sitelist="paperfinal")
dev.new(xpos=20, ypos=20, height=5, width=5)
plot(temp_flanks, out2, xlim=c(-2, 2), ylim=c(-2, 2), col=GetColorFunction(kFlanks), pch=20)
abline(0, 1, lty=2)

dev.new(xpos=520, ypos=20, height=5, width=5)
plot(temp_flanks, out3, col=GetColorFunction(kFlanks), pch=20, xlim=c(-2, 2),
     ylim=c(-2, 2))
abline(0, 1, lty=2)

dev.new(xpos=1020, ypos=20, height=5, width=5)
plot(temp_flanks, out4, col=GetColorFunction(kFlanks), pch=20, xlim=c(-2, 2),
     ylim=c(-2, 2))
abline(0, 1, lty=2)


break


pars_kds_l7 <- GetFlankKds("let-7a", sitelist="paperfinal", site="8mer")
pars_kds_m1 <- GetFlankKds("miR-1", site="8mer", sitelist="paperfinal", buffer=TRUE, combined=FALSE)
pars_kds_m1_com <- GetFlankKds("miR-1", site="8mer", sitelist="paperfinal", buffer=TRUE, combined=TRUE)

lkds_l7 <- log(pars_kds_l7[1:256, 2])
lkds_m1 <- log(pars_kds_m1[1:256, 2])
lkds_m1_com <- log(pars_kds_m1_com[1:256, 2])

dev.new(xpos=20, ypos=20, height=8, width=8)
# par(mfrow=c(2, 2))
# plot(lkds_l7, lkds_m1)
plot(lkds_l7, lkds_m1_com)
identify(lkds_l7, lkds_m1_com, labels=kFlanks)
plot(lkds_m1, lkds_m1_com)
break

dev.new(xpos=20, ypos=20, height=8, width=8)
par(mfrow=c(2, 2))
x <- avg_I_combined[, 2]
y <- -lkds_I_combined
plot(x, y, col=GetColorFunction(kFlanks), pch=19)
text(0.5, 7, labels=signif(cor(x, y)^2, 2))

x <- avg_I[, 2]
y <- -lkds_I_combined
plot(x, y, col=GetColorFunction(kFlanks), pch=19)
text(0.5, 7, labels=signif(cor(x, y)^2, 2))

x <- avg_I_combined[, 2]
y <- -lkds_I
plot(x, y, col=GetColorFunction(kFlanks), pch=19)
text(0.5, 7, labels=signif(cor(x, y)^2, 2))

x <- avg_I[, 2]
y <- -lkds_I
plot(x, y, col=GetColorFunction(kFlanks), pch=19)
text(0.5, 7, labels=signif(cor(x, y)^2, 2))


break


all_sitelists <- c("canonical", "paperfinal", "paperextendedfinal", "baek", "centered11",
                   "bulge", "del")

for (mirna in mirnas) {
  message("________________________________________________________")
  message(mirna)
  if (mirna == "miR-7-23nt") {
    out <- matrix(NaN, nrow=7, ncol=7, dimnames=list(all_sitelists,
                  c("I", "I_combined", "40", "12.65", "1.265", "0.4", "O")))
  } else {
    out <- matrix(NaN, nrow=7, ncol=8, dimnames=list(all_sitelists,
                  c("I", "I_combined", "40", "12.65", "4", "1.265", "0.4", "O")))

  }
  for (sitelist in all_sitelists) {
    if (mirna == "miR-7-23nt") {
      experiment <- "equilibrium2_nb"
    } else {
      experiment <- "equilibrium"
    }
    if (mirna == "miR-1") {
      buffer <- TRUE
    } else {
      buffer <- FALSE
    }
    sXc <- SubfunctionCall(SitesXCounts)
    out[sitelist,] <- colSums(sXc)
  }
  print(out)
}

break




flanks_8mer_miR1 <- GetFlankKds('miR-1', "equilibrium", site="8mer", buffer=TRUE)
flanks_7merm8_miR1 <- GetFlankKds('miR-1', "equilibrium", site="7mer-m8", buffer=TRUE)

flanks_GCTTCCGC_miR1 <- GetFlankKds('miR-1', "equilibrium", site="GCTTCCGC", buffer=TRUE)

fl_8mer <- flanks_8mer_miR1[grep("\\|", rownames(flanks_8mer_miR1),
                                 perl=TRUE),]
fl_7merm8 <- flanks_7merm8_miR1[grep("\\|", rownames(flanks_7merm8_miR1),
                                 perl=TRUE),]
fl_GCTTCCGC <- flanks_GCTTCCGC_miR1[grep("\\|", rownames(flanks_GCTTCCGC_miR1),
                                         perl=TRUE),]

rownames(fl_8mer) <- gsub("^(.*\\|)(.*)(_Kd$)", rownames(fl_8mer),
                          replacement="\\2", perl=TRUE)

rownames(fl_7merm8) <- gsub("^(.*\\|)(.*)(_Kd$)", rownames(fl_7merm8),
                          replacement="\\2", perl=TRUE)

rownames(fl_GCTTCCGC) <- gsub("^(.*\\|)(.*)(_Kd$)", rownames(fl_GCTTCCGC),
                          replacement="\\2", perl=TRUE)

print(fl_8mer)

pairing_dinucs <- c("AA.TT", "AC.GT", "AG.CT", "AT.AT",
                    "CA.TG", "CC.GG", "CG.CG", "CT.AG",
                    "GA.TC", "GC.GC", "GG.CC", "GT.AC",
                    "TA.TA", "TC.GA", "TG.CA", "TT.AA")

pch1 <- rep(1, nrow(fl_7merm8))
pch2 <- rep(1, nrow(fl_8mer))

names(pch1) <- rownames(fl_7merm8)
names(pch2) <- rownames(fl_8mer)

pch1[pairing_dinucs] <- 19
pch2[pairing_dinucs] <- 19
graphics.off()

dev.new(xpos=20, ypos=20, height=5, width=5)
plot(fl_8mer[intersect(rownames(fl_8mer), rownames(fl_7merm8)), ]$Mean,
     fl_7merm8[intersect(rownames(fl_8mer), rownames(fl_7merm8)), ]$Mean,
     log='xy', col=GetColorFunction(intersect(rownames(fl_8mer),
                                              rownames(fl_7merm8))),
     pch=pch1, lwd=2)
dev.copy2pdf(file="180723 8mer_vs_7merm8_flanks.pdf")


dev.new(xpos=520, ypos=20, height=5, width=5)
plot(fl_8mer$Mean, fl_GCTTCCGC$Mean, log='xy',
     col=GetColorFunction(rownames(fl_GCTTCCGC)),
     pch=pch2, lwd=2)

identify(fl_8mer$Mean, fl_GCTTCCGC$Mean, labels=rownames(fl_8mer))
dev.copy2pdf(file="180723 8mer_vs_GCTTCCGC_flanks.pdf")

break

kds_paperfinal <- EquilPars("miR-1", "equilibrium", sitelist="paperfinal",
                            buffer=TRUE, combined=FALSE)

kds_centered11 <- EquilPars("miR-1", "equilibrium", sitelist="centered11",
                            buffer=TRUE, combined=FALSE)


inds <- intersect(rownames(kds_paperfinal), rownames(kds_centered11))

plot(kds_centered11[inds, 2], kds_paperfinal[inds, 2], pch=20, cex=0.5, log='xy', col=kSiteColors[inds])


segments(0.1, 1e-3, 0.1, 10)
segments(0.2, 1e-3, 0.2, 10)
segments(0.3, 1e-3, 0.3, 10)
segments(0.4, 1e-3, 0.4, 10)

segments(1e-3, 0.1, 10, 0.1)
segments(1e-3, 0.2, 10, 0.2)
segments(1e-3, 0.3, 10, 0.3)
segments(1e-3, 0.4, 10, 0.4)


break


CheckSite <- function(mirna, site, experiment="equilibrium",
                                    n_constant=5, sitelist="paperfinal",
                                    sitelist_ref="paperfinalcheck",
                                    buffer=FALSE, combined=TRUE, fixed=FALSE,
                                    global=FALSE, xpos=20, ypos=20, ident=FALSE,
                                    removedsite=NULL, pdf.plot=FALSE, height=5, width=5) {
  pars.flanks <- SubfunctionCall(GetFlankKds, sitelist=sitelist)
  pars.flanks_ref <- SubfunctionCall(GetFlankKds, sitelist=sitelist_ref)
  flanks_inds <- grep("\\|", rownames(pars.flanks), perl=TRUE)
  flank_kds <- pars.flanks[flanks_inds,]
  flanks_inds_ref <- grep("\\|", rownames(pars.flanks_ref), perl=TRUE)
  flank_kds_ref <- pars.flanks_ref[flanks_inds_ref,]
  rownames(flank_kds) <- gsub("^(.*)\\|(.*)_Kd$", rownames(flank_kds),
                              replace="\\2", perl=TRUE)
  rownames(flank_kds_ref) <- gsub("^(.*)\\|(.*)_Kd$", rownames(flank_kds_ref),
                              replace="\\2", perl=TRUE)

  average_flanks <- exp(GetAverageFlanks())

  SubfunctionCall(FigureSaveFile)
  if (removedsite == "8mer-bA3") {
    different_flanks <- grep("\\.AC", rownames(flank_kds_ref), perl=TRUE)  
  } else if (removedsite == "7mer-m8w7bG7") {
    different_flanks <- grep("^TG", rownames(flank_kds_ref), perl=TRUE)
  } else {
    different_flanks <- setdiff(rownames(flank_kds_ref), rownames(flank_kds))
  }

  pch_all <- rep(1, length(rownames(flank_kds_ref)))
  names(pch_all) <- rownames(flank_kds_ref)
  if (site == "6mer-m8" ) {
    inds_8mer <- grep("\\.TA", rownames(flank_kds_ref), perl=TRUE)
  } else {
    inds_8mer <- different_flanks
  }


  pch_all[different_flanks] <- 19
  pch_all[inds_8mer] <- 19
  xmin <- exp(mean(log(average_flanks)))/30
  xmax <- exp(mean(log(average_flanks)))*30
  ymin <- exp(mean(log(flank_kds_ref[, 2])))/30
  ymax <- exp(mean(log(flank_kds_ref[, 2])))*30
  x_range <- seq(xmin, xmax, length=200)

  BlankPlot(log='xy', inv='xy')
  alphas = rep(0.4, length(rownames(flank_kds_ref)))
  names(alphas) <- names(pch_all)
  alphas[different_flanks] <- 1

  l_average <- log(average_flanks[rownames(flank_kds)])
  l_flanks <- log(flank_kds[, 2])
  model <- lm(l_flanks ~ l_average)
  b <- model$coefficients[1]
  m <- model$coefficients[2]
  lines(x_range, exp(log(x_range)*m + b))

  predict_range <- predict(model, data.frame(l_average=log(x_range)), interval = 'confidence')

  polygon(c(x_range, rev(x_range)), exp(c(predict_range[, 2], rev(predict_range[, 3]))), border=NA, col=ConvertRColortoRGB("red", alpha=0.3))
  lines(x_range, exp(predict_range[, 1]))

  if (removedsite %in% c("8mer-bA3", "7mer-m8w7bG7")) {
    l_average <- log(average_flanks[different_flanks])
    l_flanks <- log(flank_kds[different_flanks, 2])
    model2 <- lm(l_flanks ~ l_average)

    b <- model2$coefficients[1]
    m <- model2$coefficients[2]

    lines(x_range, exp(log(x_range)*m + b))
    predict_range2 <- predict(model2, data.frame(l_average=log(x_range)), interval = 'confidence')
    polygon(c(x_range, rev(x_range)), exp(c(predict_range2[, 2], rev(predict_range2[, 3]))), border=NA, col=ConvertRColortoRGB("blue", alpha=0.3))

    l_flanks <- log(flank_kds_ref[different_flanks, 2])
    model2 <- lm(l_flanks ~ l_average)

    b <- model2$coefficients[1]
    m <- model2$coefficients[2]

    lines(x_range, exp(log(x_range)*m + b), lty=2)



  }



  points(average_flanks[rownames(flank_kds)], flank_kds[, 2],
       col=apply(cbind(rownames(flank_kds_ref), alphas), 1, function(row) {
        GetColorFunction(row[1], alpha=row[2])}), pch=pch_all)
  xy <- GetPlotFractionalCoords(0.1, 0.95, log='xy', inv='xy')
  AddLogAxis(1, label="Average Flank Kd from Linear Model")
  AddLogAxis(2, label="Site Flanking Kd")
  text(xy[1], xy[2], labels=sprintf("%s\n%s", mirna, site), adj=c(0, 1))

  if (ident) {
    identify(average_flanks[rownames(flank_kds)], flank_kds[, 2], labels=rownames(flank_kds))
  }

  if (class(pdf.plot) == "character") {
    dev.off()
    pdf.plot <- paste0(pdf.plot, "_right")
  }

  SubfunctionCall(FigureSaveFile)



  BlankPlot(log='xy', inv='xy')

  l_average <- log(average_flanks[different_flanks])
  l_flanks <- log(flank_kds_ref[different_flanks, 2])
  model2 <- lm(l_flanks ~ l_average)

  b <- model2$coefficients[1]
  m <- model2$coefficients[2]

  lines(x_range, exp(log(x_range)*m + b))
  predict_range2 <- predict(model2, data.frame(l_average=log(x_range)), interval = 'confidence')

  polygon(c(x_range, rev(x_range)), exp(c(predict_range[, 2], rev(predict_range[, 3]))), border=NA, col=ConvertRColortoRGB("red", alpha=0.3))
  lines(x_range, exp(predict_range[, 1]))
  polygon(c(x_range, rev(x_range)), exp(c(predict_range2[, 2], rev(predict_range2[, 3]))), border=NA, col=ConvertRColortoRGB("blue", alpha=0.3))
  points(average_flanks[rownames(flank_kds_ref)], flank_kds_ref[, 2],
       col=apply(cbind(rownames(flank_kds_ref), alphas), 1, function(row) {
        GetColorFunction(row[1], alpha=row[2])}), pch=pch_all)
  xy <- GetPlotFractionalCoords(0.1, 0.95, log='xy', inv='xy')
  AddLogAxis(1, label="Average Flanks from Linear Model")
  AddLogAxis(2, label="Site Flanking Kds")
  text(xy[1], xy[2], labels=sprintf("%s\n%s", mirna, site), adj=c(0, 1))
  if (!is.null(removedsite)) {
    xy <- GetPlotFractionalCoords(0.1, 0.75, log='xy', inv='xy')
    text(xy[1], xy[2], labels=sprintf("%s removed",removedsite), adj=c(0, 1)) 
  }
  if (class(pdf.plot) == "character") {
    dev.off()
  }
  return(flank_kds)
}

graphics.off()
# CheckSite("miR-1", "6mer-m8", buffer=TRUE, combined=FALSE)
# # CheckSite("miR-1", "5mer-m2.6", buffer=TRUE, combined=FALSE)
# # dev.copy2pdf(file="180797 miR-1_6mer-m8_vs_8mer-w2.pdf")
# CheckSite("lsy-6", "7mer-A1")
# # dev.copy2pdf(file="180797 lsy-6_7mer-A1_vs_8mer-w8.pdf")
# CheckSite("lsy-6", "6mer")
# # dev.copy2pdf(file="180797 lsy-6_6mer_vs_7mer-m8w8.pdf")
# CheckSite("miR-7-23nt", "6mer-m8", experiment="equilibrium2_nb", combined=FALSE)
# # dev.copy2pdf(file="180797 miR-7_6mer-m8_vs_8mer-w2.pdf")




# sXc <- SitesXCounts("miR-1", experiment="equilibrium", buffer=TRUE, sitelist="paperfinal")
# sXc_w8 <- SitesXCounts("miR-1", experiment="equilibrium", buffer=TRUE, sitelist="paperfinalcheck")


# print(sXc[rownames(sXc_w8), ] - sXc_w8)


# sXc_7merm8bA8 <- SitesXCounts("miR-7-23nt", experiment="equilibrium2_nb", sitelist="paperfinalcheck7merm8bA8")


### miR-1
# 8mer-w2
# 7mer-m8w2
# CheckSite("miR-1", "6mer-m8", buffer=TRUE, combined=FALSE,
#           removedsite="8mer-w2 and\n7mer-m8w2", pdf.plot="0.180708 miR-1_remove_8mer-w2and7mer-m8w2")



sXc <- SitesXCounts("miR-1", sitelist="paperfinalcheck8merbA3",
                    buffer=TRUE)

sXc_ref <- SitesXCounts("miR-1", buffer=TRUE)

sXc_del <- sXc_ref - sXc[rownames(sXc_ref), ] 

rows_del <- setdiff(rownames(sXc), rownames(sXc_ref))
sXc_ref_del <- sXc[rows_del, ]
print(sXc_del)
print(sXc_ref_del)

print(colSums(sXc_ref_del[2:4, 3, drop=FALSE]))
print(sXc_del["5mer-m2.6",] + sXc_ref_del["7mer-m8w7bG7", ])

CheckSite("miR-1", "6mer-m8", buffer=TRUE, combined=FALSE,
          sitelist="paperfinalcheck8merbA3", removedsite="8mer-bA3",
          sitelist_ref="paperfinal", pdf.plot="0.180709 miR-1_remove_8mer-bA3")

CheckSite("miR-1", "5mer-m2.6", buffer=TRUE, combined=FALSE,
          sitelist="paperfinalcheck7merm8w7bG7", removedsite="7mer-m8w7bG7",
          sitelist_ref="paperfinal", pdf.plot="0.180709 miR-1_remove_7mer-m8w7bG7")



# print(colSums(sXc_del))
# print(sXc[setdiff(rownames(sXc), rownames(sXc_ref)), ])
# print(colSums(sXc[setdiff(rownames(sXc), rownames(sXc_ref)), ]))
# print(colSums(sXc[setdiff(rownames(sXc), rownames(sXc_ref)), ]))

break

CheckSite("miR-1", "5mer-m2.6", buffer=TRUE, sitelist="paperfinalcheck", 
          sitelist_ref="paperfinal", combined=FALSE,
          removedsite="7mer-m8w7bG7", pdf.plot=)

break


### miR-124
# 8mer-bU(7.8)
CheckSite("miR-124", "7mer-A1", combined=FALSE,
          sitelist_ref="paperfinalcheck7merm8bU78",
          removedsite="8mer-bU(7.8)", ypos=920, pdf.plot="0.180708 miR-124_remove_8mer-bU(7,8)")
# 7mer-m8bU(7.8)
CheckSite("miR-124", "6mer",  combined=FALSE,
          sitelist_ref="paperfinalcheck7merm8bU78",
          removedsite="7mer-m8bU(7.8)", ypos=1220, pdf.plot="0.180708 miR-124_remove_7mer-m8bU(7,8)")

### lsy-6
# 8mer-w8
CheckSite("lsy-6", "7mer-A1", combined=FALSE, removedsite="8mer-w8", xpos=720, pdf.plot="0.180708 lsy-6_remove_8mer-w8")
# 7mer-m8w8
CheckSite("lsy-6", "6mer", combined=FALSE, removedsite="7mer-m8w8", xpos=1300, pdf.plot="0.180708 lsy-6_remove_7mer-m8w8")




### miR-7
# 8mer-bU(7.8)
CheckSite("miR-7-23nt", "7mer-A1", sitelist_ref="paperfinalcheck7merm8bU78",
          experiment="equilibrium2_nb", combined=FALSE,
          removedsite="8mer-bU(7.8)", xpos=720, ypos=920, pdf.plot="0.180708 miR-7_remove_8mer-bU(7,8)")
# 7mer-m8bU(7.8)
CheckSite("miR-7-23nt", "6mer", sitelist_ref="paperfinalcheck7merm8bU78",
          experiment="equilibrium2_nb", combined=FALSE,
          removedsite="7merbU(7.8)", pdf.plot="0.180708 miR-7_remove_7mer-m8bU(7,8)")
# 8mer-bA8
CheckSite("miR-7-23nt", "7mer-A1", sitelist_ref="paperfinalcheck7merm8bA8",
          experiment="equilibrium2_nb", combined=FALSE,
          removedsite="8mer-bA8", pdf.plot="0.180708 miR-7_remove_8mer-bA8")
# 7mer-m8bA8
CheckSite("miR-7-23nt", "6mer", sitelist_ref="paperfinalcheck7merm8bA8",
          experiment="equilibrium2_nb", combined=FALSE,
          removedsite="7mer-m8bA8", pdf.plot="0.180708 miR-7_remove_7mer-m8bA8")
# 8mer-m8w2
# 7mer-m8w2
CheckSite("miR-7-23nt", "6mer-m8", experiment="equilibrium2_nb", combined=FALSE,
          removedsite="8mer-w2 and 7mer-m8w2", pdf.plot="0.180708 miR-7_remove_8mer-w2and7mer-m8w2")

break



# # dev.copy2pdf(file="180797 miR-1_6mer-m8_vs_8mer-w2.pdf")
# CheckSite("lsy-6", "7mer-A1")
# # dev.copy2pdf(file="180797 lsy-6_7mer-A1_vs_8mer-w8.pdf")
# CheckSite("lsy-6", "6mer")
# # dev.copy2pdf(file="180797 lsy-6_6mer_vs_7mer-m8w8.pdf")
# CheckSite("miR-7-23nt", "6mer-m8", experiment="equilibrium2_nb", combined=FALSE)
# # dev.copy2pdf(file="180797 miR-7_6mer-m8_vs_8mer-w2.pdf")



CheckSite("miR-7-23nt", "6mer", sitelist_ref="paperfinalcheck7merm8bA8",
          experiment="equilibrium2_nb", combined=FALSE)

dev.copy2pdf(file="180897 miR-7_6mer_vs_7mer-m8bA8.pdf")
CheckSite("miR-7-23nt", "6mer-A1", sitelist_ref="paperfinalcheck7merA1bG7",
          experiment="equilibrium2_nb", ident=TRUE, combined=FALSE)

CheckSite("miR-7-23nt", "6mer-A1", sitelist_ref="paperfinalcheck7merA1bG7",
          experiment="equilibrium2_nb", ident=TRUE, combined=FALSE)


dev.copy2pdf(file="180897 miR-7_6mer-A1_vs_7mer-A1bG7.pdf")



break




string_col <- "4"
string_col <- "40"

# kmers_let7a_I_compare <- data.matrix(data.frame(read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/kmers_cutoff_final/I_5_k10_buffer3p.txt", header=TRUE, stringsAsFactors=FALSE), row.names=1))
# kmers_miR155_I_compare <- data.matrix(data.frame(read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/miR-155/equilibrium/kmers_cutoff_final/I_5_k10_buffer3p.txt", header=TRUE, stringsAsFactors=FALSE), row.names=1))
# kmers_miR124_I_compare <- data.matrix(data.frame(read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/miR-124/equilibrium/kmers_cutoff_final/I_5_k10_buffer3p.txt", header=TRUE, stringsAsFactors=FALSE), row.names=1))
# kmers_lsy6_I_compare <- data.matrix(data.frame(read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/lsy-6/equilibrium/kmers_cutoff_final/I_5_k10_buffer3p.txt", header=TRUE, stringsAsFactors=FALSE), row.names=1))
# kmers_miR7_I_compare <- data.matrix(data.frame(read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/miR-7-23nt/equilibrium2_nb/kmers_cutoff_final/I_5_k10_buffer3p.txt", header=TRUE, stringsAsFactors=FALSE), row.names=1))

kmers_miR1_0_0 <- data.matrix(data.frame(read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/miR-1/equilibrium/kmers_cutoff_final/0_5_k10_buffer3p.txt", header=TRUE, stringsAsFactors=FALSE), row.names=1))
kmers_miR1_I_0 <- data.matrix(data.frame(read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/miR-1/equilibrium/kmers_cutoff_final/I_5_k10_buffer3p.txt", header=TRUE, stringsAsFactors=FALSE), row.names=1))
kmers_miR1_4_0 <- data.matrix(data.frame(read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/miR-1/equilibrium/kmers_cutoff_final/40_5_k10_buffer3p.txt", header=TRUE, stringsAsFactors=FALSE), row.names=1))


n_ex <- 24

kmers_miR1_0_1 <- data.matrix(data.frame(read.table(paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/miR-1/equilibrium/kmers_cutoff_final/0_5_k10_buffer3p_ex", n_ex, ".txt"), header=TRUE, stringsAsFactors=FALSE), row.names=1))
kmers_miR1_I_1 <- data.matrix(data.frame(read.table(paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/miR-1/equilibrium/kmers_cutoff_final/I_5_k10_buffer3p_ex", n_ex, ".txt"), header=TRUE, stringsAsFactors=FALSE), row.names=1))
kmers_miR1_4_1 <- data.matrix(data.frame(read.table(paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/miR-1/equilibrium/kmers_cutoff_final/40_5_k10_buffer3p_ex", n_ex, ".txt"), header=TRUE, stringsAsFactors=FALSE), row.names=1))

div <- 1e8


f_I_0 <- kmers_miR1_I_0[, 1]
f_4_0 <- kmers_miR1_4_0[, 1]
f_0_0 <- kmers_miR1_0_0[, 1]
f_I_1 <- kmers_miR1_I_1[, 1] + median(kmers_miR1_I_1[, 1])
f_4_1 <- kmers_miR1_4_1[, 1] + median(kmers_miR1_4_1[, 1])
f_0_1 <- kmers_miR1_0_1[, 1] + median(kmers_miR1_0_1[, 1])

R_miR1_0 <- f_4_0/f_I_0
R_miR1_0_0 <- f_0_0/f_I_0

R_miR1_1 <- f_4_1/f_I_1
R_miR1_1_0 <- f_0_1/f_I_1


lR_0 <- log10(R_miR1_0)
lR_0_0 <- log10(R_miR1_0_0)

lR_1 <- log10(R_miR1_1)
lR_1_0 <- log10(R_miR1_1_0)


z_0 <- (lR_0 - mean(lR_0))/sd(lR_0)
z_0_0 <- (lR_0_0 - mean(lR_0_0))/sd(lR_0)

z_1 <- (lR_1 - mean(lR_1))/sd(lR_1)
z_1 <- z_1/(max(z_1))*z_0[which(z_1 == max(z_1))]
z_1_0 <- (lR_1_0 - mean(lR_1_0))/sd(lR_1)

dev.new(xpos=20, ypos=20, height=10, width=10)
par(mfrow=c(2, 2))
plot(ecdf(z_0))
plot(ecdf(z_1), add=TRUE, col="blue")
plot(z_0, z_1, col=rgb(0, 0, 0, alpha=0.2), xlim=c(-2, 15), ylim=c(-2, 15))
abline(0, 1, lty=2, col="gray")


plot(ecdf(z_0_0))
plot(ecdf(z_1_0), add=TRUE, col="blue")
plot(z_0_0, z_1_0, col=rgb(0, 0, 0, alpha=0.2), xlim=c(-2, 15), ylim=c(-2, 15))
abline(0, 1, lty=2, col="gray")


break

kmers_let7a_I <- data.matrix(data.frame(read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/kmers_cutoff_final/I_5_k10.txt", header=TRUE, stringsAsFactors=FALSE), row.names=1))
kmers_let7a_4 <- data.matrix(data.frame(read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/kmers_cutoff_final/40_5_k10.txt", header=TRUE, stringsAsFactors=FALSE), row.names=1))

R_let7a <- (kmers_let7a_4[, 1] / sum(kmers_let7a_4[, 1], na.rm=TRUE))/((kmers_let7a_I[, 1] + 1)/sum(kmers_let7a_I[, 1] + 1))

kmers_miR155_I <- data.matrix(data.frame(read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/miR-155/equilibrium/kmers_cutoff_final/I_5_k10.txt", header=TRUE, stringsAsFactors=FALSE), row.names=1))
kmers_miR155_4 <- data.matrix(data.frame(read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/miR-155/equilibrium/kmers_cutoff_final/40_5_k10.txt", header=TRUE, stringsAsFactors=FALSE), row.names=1))

R_miR155 <- (kmers_miR155_4[, 1] / sum(kmers_miR155_4[, 1], na.rm=TRUE))/((kmers_miR155_I[, 1] + 1)/sum(kmers_miR155_I[, 1] + 1))

kmers_miR124_I <- data.matrix(data.frame(read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/miR-124/equilibrium/kmers_cutoff_final/I_5_k10.txt", header=TRUE, stringsAsFactors=FALSE), row.names=1))
kmers_miR124_4 <- data.matrix(data.frame(read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/miR-124/equilibrium/kmers_cutoff_final/40_5_k10.txt", header=TRUE, stringsAsFactors=FALSE), row.names=1))

R_miR124 <- (kmers_miR124_4[, 1] / sum(kmers_miR124_4[, 1], na.rm=TRUE))/((kmers_miR124_I[, 1] + 1)/sum(kmers_miR124_I[, 1] + 1))

kmers_lsy6_I <- data.matrix(data.frame(read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/lsy-6/equilibrium/kmers_cutoff_final/I_5_k10.txt", header=TRUE, stringsAsFactors=FALSE), row.names=1))
kmers_lsy6_4 <- data.matrix(data.frame(read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/lsy-6/equilibrium/kmers_cutoff_final/40_5_k10.txt", header=TRUE, stringsAsFactors=FALSE), row.names=1))

R_lsy6 <- (kmers_lsy6_4[, 1] / sum(kmers_lsy6_4[, 1], na.rm=TRUE))/((kmers_lsy6_I[, 1] + 1)/sum(kmers_lsy6_I[, 1] + 1))

kmers_miR7_I <- data.matrix(data.frame(read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/miR-7-23nt/equilibrium2_nb/kmers_cutoff_final/I_5_k10.txt", header=TRUE, stringsAsFactors=FALSE), row.names=1))
kmers_miR7_4 <- data.matrix(data.frame(read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/miR-7-23nt/equilibrium2_nb/kmers_cutoff_final/40_5_k10.txt", header=TRUE, stringsAsFactors=FALSE), row.names=1))

R_miR7 <- (kmers_miR7_4[, 1] / sum(kmers_miR7_4[, 1], na.rm=TRUE))/((kmers_miR7_I[, 1] + 1)/sum(kmers_miR7_I[, 1] + 1))


R_all_log <- log10(cbind(R_miR1, R_let7a, R_miR155, R_miR124, R_lsy6, R_miR7))
R_all     <- cbind(R_miR1, R_let7a, R_miR155, R_miR124, R_lsy6, R_miR7)





check_8mers <- c("ACATTCCA", "CTACCTCA", "AGCATTAA", "GTGCCTTA", "ATACAAAA", "GTCTTCCA")


ratios <- sapply(check_8mers, function(s_8mer) {
  print(s_8mer)
})

break

check_kmers <- c("ACATTCCA", "CTACCTCA", "AGCATTAA", "GTGCCTTA", "ATACAAAA", "GTCTTCCA",
                 "CGCTTCCG", "GCTTCCGC", "CTTCCGCT", #miR-1
                 "CACACACA", "ACACACAC", # miR-1 1
                 "TGCACTTT", "GCACTTTA", # let-7a 1
                 "TGCGCACC", "CTGCGCAC", # let-7a 2
                 "AACGAGGA", "ACGAGGAA", # miR-155 1
                 "AACTCAGC", "ACTCAGCA", # miR-155 2
                 "TACTCAGC", "ACTCAGCA", # miR-155 3
                 "ATGACAAC", "TGACAACA", # miR-155 4
                 "AAAATAAA", "AAATAAAG", # miR-155 5
                 "CTCAGCAA", "TCAGCAAT", # miR-155 6
                 "CAACTCAG", "AACTCAGC", # miR-155 7
                 "AAACGACA", "AACGACAA", # miR-155 8
                 "TCGACAAC", "CGACAACA", # miR-155 9
                 "CGACAACT", "GACAACTC", # miR-155 10
                 "TAACGAGG", "AACGAGGT", # miR-155 11
                 "ACGACAAC", "CGACAACT", # miR-155 12
                 "AAGTGCAA", "AGTGCAAT", # miR-124 1
                 "CTAAGTGC", "TAAGTGCC", # miR-124 2
                 "CCGCCACA", "CGCCACAG", # miR-124 3
                 "GTCACCCG", "TCACCCGC", # miR-124 4
                 "TCCGCCAC", "CCGCCACA", # lsy-6 1
                 "CCGCCACA", "CGCCACAG", # lsy-6 2
                 "CACAGAAA", "ACAGAAAT", # lsy-6 3
                 "CCTCCGCC", "CTCCGCCA", # lsy-6 4
                 "GCCACAGA", "CCACAGAA", # lsy-6 5
                 "ATCCGCCA", "TCCGCCAC", # lsy-6 6
                 "CCTCTGCC", "CTCTGCCC", # lsy-6 7
                 "AACGAGGA", "ACGAGGAA", # lsy-6 8
                 "TCTTCCTA", "CTTCCTAT", # lsy-6 9
                 "ATTCTTCC", "TTCTTCCT", # lsy-6 10
                 "ATCTTCCT", "TCTTCCTA", # lsy-6 11
                 "TCCTTCCT", "CCTTCCTA", # lsy-6 12
                 "CGCTTCCG", "GCTTCCGC", # miR-7 1
                 "CTTCCGCT", "TTCCGCTG", # miR-7 2
                 "CGCTTCCG", "GCTTCCGC", # miR-7 3
                 "CGCTTCCG", "GCTTCCGT", # miR-7 4
                 "CCGCACCA", "CGCACCAC", # miR-7 5
                 "TTCCGCTG", "TCCGCTGC", # miR-7 6
                 "CGCACCAC", "GCACCACA", # miR-7 7
                 "TCCGCACC", "CCGCACCA", # miR-7 8
                 "CGTTTCCG", "GTTTCCGC", # miR-7 9
                 "TCCGCTGC", "CCGCTGCG", # miR-7 10
                 "CGCTTTCG", "GCTTTCGC", # miR-7 11
                 "CCTCCGCA", "CTCCGCAC", # miR-7 12
                 "AACGACTA", "ACGACTAA", # miR-7 13
                 "ATACGGCT", "TACGGCTA", # miR-7 14
                 "TGCACCAC", "GCACCACA", # miR-7 15
                 "AAGTGCCT",             # miR-124 AA-6mer-m8
                 "AAGTGTCC", "AGTGTCCT", "GTGTCCTT", "TGTCCTTA", # miR-124 AA-8mer-bT6
                 "AAGTGCCC", "AGTGCCCT", "GTGCCCTT", "TGCCCTTA", # miR-124 AA-8mer-C(4.6)
                 "AAGTGACC", "AGTGACCT", "GTGACCTT", "TGACCTTA", # miR-124 AA-8mer-bA5
                 "GGGCTTCC", "GGCTTCCA", # miR-7 8mer-mmG7bG7
                 "GAGCTTCC", "AGCTTCCA", # miR-7 8mer-mmA7bG7
                 "GCGCTTCC", "CGCTTCCA") # miR-7 8mer-mmC7bG7


cols <- c(kMirnaColors, rep.int(kMirnaColors, c(5, 4, 24, 8, 24, 30)), rep(kMirnaColors["miR-124"], 13), rep(kMirnaColors["miR-7-23nt"], 6))
R_means <- rowMeans(R_all)

print(R_all[check_8mers, ])
dev.new(xpos=20, ypos=20, height=5, width=5)
heatmap(log10(R_all[check_8mers, ]), Rowv=NA, Colv=NA, margins=c(8, 8))
print(R_all[check_kmers, ])
dev.new(xpos=520, ypos=20, height=10, width=10)
heatmap(log10(R_all[check_kmers, ]), Rowv=NA, Colv=NA, margins=c(8, 8), scale="none", RowSideColors=cols)

break

kmers_7 <- data.matrix(data.frame(read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/miR-155/equilibrium/kmers_cutoff_final/I_5_k10_ex30.txt", header=TRUE, stringsAsFactors=FALSE), row.names=1))
kmers_8 <- data.matrix(data.frame(read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/miR-155/equilibrium/kmers_cutoff_final/I_5_k10_ex16.txt", header=TRUE, stringsAsFactors=FALSE), row.names=1))

kmers_6_4 <- data.matrix(data.frame(read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/miR-155/equilibrium/kmers_cutoff_final/4_5_k10_ex29.txt", header=TRUE, stringsAsFactors=FALSE), row.names=1))
kmers_7_4 <- data.matrix(data.frame(read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/miR-155/equilibrium/kmers_cutoff_final/4_5_k10_ex30.txt", header=TRUE, stringsAsFactors=FALSE), row.names=1))
kmers_8_4 <- data.matrix(data.frame(read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/miR-155/equilibrium/kmers_cutoff_final/4_5_k10_ex16.txt", header=TRUE, stringsAsFactors=FALSE), row.names=1))

R_6 <- (kmers_6_4[, 1] / sum(kmers_6_4[, 1], na.rm=TRUE))/((kmers_6[, 1] + 1)/sum(kmers_6[, 1] + 1))
R_7 <- (kmers_7_4[, 1] / sum(kmers_7_4[, 1], na.rm=TRUE))/((kmers_7[, 1] + 1)/sum(kmers_7[, 1] + 1))
R_8 <- (kmers_8_4[, 1] / sum(kmers_8_4[, 1], na.rm=TRUE))/((kmers_8[, 1] + 1)/sum(kmers_8[, 1] + 1))

kmers_names <- rownames(kmers_6)

names(R_6) <- kmers_names
names(R_7) <- kmers_names
names(R_8) <- kmers_names


kmers_6_7_del <- kmers_6[, 1] - kmers_7[, 1]
kmers_7_8_del <- kmers_7[, 1] - kmers_8[, 1]

kmers_6_7_del_4 <- kmers_6_4[, 1] - kmers_7_4[, 1]
kmers_7_8_del_4 <- kmers_7_4[, 1] - kmers_8_4[, 1]


names(kmers_6_7_del) <- kmers_names
names(kmers_7_8_del) <- kmers_names
names(kmers_6_7_del_4) <- kmers_names
names(kmers_7_8_del_4) <- kmers_names

print(head(rev(sort(kmers_6_7_del))))
print(head(rev(sort(kmers_6_7_del_4))))


plot(R_6, R_7, log='xy')

break

test_all <- GetPositionalSites("miR-1", "equilibrium", "40", single=TRUE)
test_all_buffer <- GetPositionalSites("miR-1", "equilibrium", "40", single=TRUE, buffer=TRUE)


dev.new(xpos=20, ypos=20, height=5, width=5)

plot(1:ncol(test_all), test_all[1, ], type="o")
points(1:ncol(test_all), test_all_buffer[1, ], type="o", col="blue")


break

######################## KEEP


dev.new(xpos=20, ypos=20, height=12, width=7)

R_40 <- (test_40/test_I)
R_40 <- data.matrix(R_40)
R_40 <- log10(R_40)
R_40_old <- R_40

R_40 <- R_40 - rowMeans(R_40, na.rm=TRUE)
cols_keep <- max(which(colSums(is.na(R_40)) != nrow(R_40)))
R_40 <- R_40[, 1:cols_keep]

plot(ecdf(R_40))
R_40 <- R_40/sd(R_40, na.rm=TRUE)
# plot(1:ncol(R_40_old), R_40_old["6mer-m8",], type="o", ylim=c(-1, 3))
# points(1:ncol(R_40_old), R_40["6mer-m8",], type="o", col="gray")
R_40[which(is.na(R_40))] <- -10

xmin <- -1/ncol(R_40)
xmax <- 1.4
ymin <- -1/nrow(R_40)
ymax <- 1 + 1/nrow(R_40)

par(kPlotParameters)
BlankPlot()

image(t(data.matrix(R_40)), breaks=c(-15, seq(-5, 5)), col=c("gray",topo.colors(15)[3:12]), xpd=NA, add=TRUE)
text(x=1.1, y=seq(0, nrow(R_40) - 1)/(nrow(R_40) - 1), adj=0, labels=rownames(R_40), xpd=NA)

break

PlotSiteEnrichment("7mer-A1")

PlotSiteEnrichment("6mer-m8", ypos=520)







break
########################## END KEEP


# break

# points(1:ncol(test_1.26), test_1.26["6mer-m8", ]/test_I["6mer-m8",], type="o", col="orangered")
# points(1:ncol(test_4), test_4["6mer-m8", ]/test_I["6mer-m8",], type="o", col="forestgreen")
# points(1:ncol(test_12.6), test_12.6["6mer-m8", ]/test_I["6mer-m8",], type="o", col="blue")
# points(1:ncol(test_40), test_40["6mer-m8", ]/test_I["6mer-m8",], type="o", col="purple")



s_6merm8_all <- rbind(test_I["6mer-m8",], test_0.4["6mer-m8", ], test_1.26["6mer-m8", ],
              test_4["6mer-m8", ], test_12.6["6mer-m8", ],
              test_40["6mer-m8", ])


rownames(s_6merm8_all) <- colnames(sXc)[-2][-7]
s_6merm8_all <- t(s_6merm8_all)


ind_n31 <- which(colnames(test_I) == "N31")
ind_n32 <- which(colnames(test_I) == "N32")
ind_n33 <- which(colnames(test_I) == "N33")

s_6merm8_counts <- c(rowSums(test_0.4["6mer-m8", 1:(ind_n31 - 1)]), rowSums(test_1.26["6mer-m8", 1:(ind_n31 - 1)]),
              rowSums(test_4["6mer-m8", 1:(ind_n31 - 1)]), rowSums(test_12.6["6mer-m8", 1:(ind_n31 - 1)]),
              rowSums(test_40["6mer-m8", 1:(ind_n31 - 1)]))
s_6merm8 <- s_6merm8_counts/(rowSums(test_I["6mer-m8", 1:(ind_n31 - 1)]))

s_6merm8_spike_counts <- c(test_0.4["6mer-m8", "N31"], test_1.26["6mer-m8", "N31"],
              test_4["6mer-m8", "N31"], test_12.6["6mer-m8", "N31"],
              test_40["6mer-m8", "N31"])
s_6merm8_spike <- s_6merm8_spike_counts/(test_I["6mer-m8", "N31"])



# s_7merm8w3_counts <- c(test_0.4["7mer-m8w3", "N16"], test_1.26["7mer-m8w3", "N16"],
#               test_4["7mer-m8w3", "N16"], test_12.6["7mer-m8w3", "N16"],
#               test_40["7mer-m8w3", "N16"])
# s_7merm8w3 <- s_7merm8w3_counts/(test_I["7mer-m8w3", "N16"])


s_7merm8w3_counts <- c(rowSums(test_0.4["7mer-m8w3", 1:(ind_n32 - 1)]), rowSums(test_1.26["7mer-m8w3", 1:(ind_n32 - 1)]),
              rowSums(test_4["7mer-m8w3", 1:(ind_n32 - 1)]), rowSums(test_12.6["7mer-m8w3", 1:(ind_n32 - 1)]),
              rowSums(test_40["7mer-m8w3", 1:(ind_n32 - 1)]))
s_7merm8w3 <- s_7merm8w3_counts/(rowSums(test_I["7mer-m8w3", 1:(ind_n32 - 1)]))



s_7merm8w3_spike_counts <- c(test_0.4["7mer-m8w3", "N32"], test_1.26["7mer-m8w3", "N32"],
              test_4["7mer-m8w3", "N32"], test_12.6["7mer-m8w3", "N32"],
              test_40["7mer-m8w3", "N32"])
s_7merm8w3_spike <- s_7merm8w3_spike_counts/(test_I["7mer-m8w3", "N32"])


s_7merm8w3_bump_counts <- c(test_0.4["7mer-m8w3", "N16"], test_1.26["7mer-m8w3", "N16"],
              test_4["7mer-m8w3", "N16"], test_12.6["7mer-m8w3", "N16"],
              test_40["7mer-m8w3", "N16"])
s_7merm8w3_bump <- s_7merm8w3_bump_counts/(test_I["7mer-m8w3", "N16"])




# s_6mermmT3_counts <- c(test_0.4["6mer-mmT3", "N16"], test_1.26["6mer-mmT3", "N16"],
#               test_4["6mer-mmT3", "N16"], test_12.6["6mer-mmT3", "N16"],
#               test_40["6mer-mmT3", "N16"])
# s_6mermmT3 <- s_6mermmT3_counts/(test_I["6mer-mmT3", "N16"])


s_6mermmT3_counts <- c(rowSums(test_0.4["6mer-mmT3", 1:(ind_n33 - 1)]), rowSums(test_1.26["6mer-mmT3", 1:(ind_n33 - 1)]),
              rowSums(test_4["6mer-mmT3", 1:(ind_n33 - 1)]), rowSums(test_12.6["6mer-mmT3", 1:(ind_n33 - 1)]),
              rowSums(test_40["6mer-mmT3", 1:(ind_n33 - 1)]))
s_6mermmT3 <- s_6mermmT3_counts/(rowSums(test_I["6mer-mmT3", 1:(ind_n33 - 1)]))



s_6mermmT3_spike_counts <- c(test_0.4["6mer-mmT3", "N33"], test_1.26["6mer-mmT3", "N33"],
              test_4["6mer-mmT3", "N33"], test_12.6["6mer-mmT3", "N33"],
              test_40["6mer-mmT3", "N33"])
s_6mermmT3_spike <- s_6mermmT3_spike_counts/(test_I["6mer-mmT3", "N33"])

# s_5merm3.7_counts <- c(test_0.4["5mer-m3.7", "N16"], test_1.26["5mer-m3.7", "N16"],
#               test_4["5mer-m3.7", "N16"], test_12.6["5mer-m3.7", "N16"],
#               test_40["5mer-m3.7", "N16"])
# s_5merm3.7 <- s_5merm3.7_counts/(test_I["5mer-m3.7", "N16"])


s_5merm3.7_counts <- c(rowSums(test_0.4["5mer-m3.7", 1:(ind_n32 - 1)]), rowSums(test_1.26["5mer-m3.7", 1:(ind_n32 - 1)]),
              rowSums(test_4["5mer-m3.7", 1:(ind_n32 - 1)]), rowSums(test_12.6["5mer-m3.7", 1:(ind_n32 - 1)]),
              rowSums(test_40["5mer-m3.7", 1:(ind_n32 - 1)]))
s_5merm3.7 <- s_5merm3.7_counts/(rowSums(test_I["5mer-m3.7", 1:(ind_n32 - 1)]))




s_5merm3.7_spike_counts <- c(test_0.4["5mer-m3.7", "N32"], test_1.26["5mer-m3.7", "N32"],
              test_4["5mer-m3.7", "N32"], test_12.6["5mer-m3.7", "N32"],
              test_40["5mer-m3.7", "N32"])
s_5merm3.7_spike <- s_5merm3.7_spike_counts/(test_I["5mer-m3.7", "N32"])



s_6merm8_bump_counts <- c(test_0.4["6mer-m8", "N4"], test_1.26["6mer-m8", "N4"],
              test_4["6mer-m8", "N4"], test_12.6["6mer-m8", "N4"],
              test_40["6mer-m8", "N4"])
s_6merm8_bump <- s_6merm8_bump_counts/(test_I["6mer-m8", "N4"])


# s_6merm8_counts <- c(test_0.4["6mer-m8", "N16"], test_1.26["6mer-m8", "N16"],
#               test_4["6mer-m8", "N16"], test_12.6["6mer-m8", "N16"],
#               test_40["6mer-m8", "N16"])
# s_6merm8 <- s_6merm8_counts/(test_I["6mer-m8", "N16"])

s_8mer_counts <- c(test_0.4["8mer", "N16"], test_1.26["8mer", "N16"],
              test_4["8mer", "N16"], test_12.6["8mer", "N16"],
              test_40["8mer", "N16"])
s_8mer <- s_8mer_counts/(test_I["8mer", "N16"])

s_7merm8_counts <- c(test_0.4["7mer-m8", "N16"], test_1.26["7mer-m8", "N16"],
              test_4["7mer-m8", "N16"], test_12.6["7mer-m8", "N16"],
              test_40["7mer-m8", "N16"])
s_7merm8 <- s_7merm8_counts/(test_I["7mer-m8", "N16"])


s_7merA1_counts <- c(test_0.4["7mer-A1", "N16"], test_1.26["7mer-A1", "N16"],
              test_4["7mer-A1", "N16"], test_12.6["7mer-A1", "N16"],
              test_40["7mer-A1", "N16"])
s_7merA1 <- s_7merA1_counts/(test_I["7mer-A1", "N16"])

s_6mer_counts <- c(test_0.4["6mer", "N16"], test_1.26["6mer", "N16"],
              test_4["6mer", "N16"], test_12.6["6mer", "N16"],
              test_40["6mer", "N16"])
s_6mer <- s_6mer_counts/(test_I["6mer", "N16"])



x_pos <- c(0.4, 1.26, 4, 12.6, 40)


ModelCombination <- function(pars) {
  pars <- exp(pars)
  tots <- c(pars, 1)/(sum(pars) + 1)
  s_6merm8*tots[1] + s_7merm8*tots[2] + s_8mer[3]*tots[3]
}


ResidualFunction <- function(pars) {
  final <- ModelCombination(pars)
  # plot(x_pos, right_spike, type="o", log='xy', ylim=c(0.3, 200))
  # points(x_pos, final, type="o", col="green")
  sum((final - s_6merm8)^2)
}


# dev.new(xpos=20, ypos=20, height=5, width=5)
pars_final <- optim(c(1, 1), ResidualFunction)

print(pars_final$par)
fit_line <- ModelCombination(pars_final$par)

# "AGCAAAGCATGAAAGCGGACAAAATACACTAATACATTCGTAT"
# "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTCGTAT"
# "                               ACATTC"          "6mer-m8"
# "                               ACATTCC"         "7mer-m8"
# "0123456789!123456789@123456789#1"
# "123456789!123456789@123456789#123"

# dev.new(xpos=20, ypos=20, height=5, width=5)
# plot(1:ncol(test_I), test_I["6mer-m8", ], type="o", col="black", log='y')
# points(1:ncol(test_0.4), test_0.4["6mer-m8", ], type="o", col="red")
# points(1:ncol(test_1.26), test_1.26["6mer-m8", ], type="o", col="orangered")
# points(1:ncol(test_4), test_4["6mer-m8", ], type="o", col="forestgreen")
# points(1:ncol(test_12.6), test_12.6["6mer-m8", ], type="o", col="blue")
# points(1:ncol(test_40), test_40["6mer-m8", ], type="o", col="purple")

# dev.new(xpos=520, ypos=20, height=5, width=5)
# par(kPlotParameters)
# xmin <- 5
# xmax <- 40
# ymin <- 0.1
# ymax <- 1e3
# BlankPlot(log='y')
# AddLinearAxis(1, tick.space=1, label.space=5, label="Position")
# AddLogAxis(2, label="Enrichment")
# points(1:ncol(test_0.4), test_0.4["6mer-m8", ]/test_I["6mer-m8",], type="o", col="red")
# points(1:ncol(test_1.26), test_1.26["6mer-m8", ]/test_I["6mer-m8",], type="o", col="orangered")
# points(1:ncol(test_4), test_4["6mer-m8", ]/test_I["6mer-m8",], type="o", col="forestgreen")
# points(1:ncol(test_12.6), test_12.6["6mer-m8", ]/test_I["6mer-m8",], type="o", col="blue")
# points(1:ncol(test_40), test_40["6mer-m8", ]/test_I["6mer-m8",], type="o", col="purple")

# break
PlotSiteEnrichment("7mer-m8")
dev.copy2pdf(file="180527_miR-1_positional_analysis_7mer-m8.pdf", useDingbats=FALSE)
PlotSiteEnrichment("6mer")
dev.copy2pdf(file="180527_miR-1_positional_analysis_6mer.pdf", useDingbats=FALSE)
PlotSiteEnrichment("7mer-m8w3")
dev.copy2pdf(file="180527_miR-1_positional_analysis_7mer-m8w3.pdf", useDingbats=FALSE)
PlotSiteEnrichment("6mer-m8", xpos=520)
dev.copy2pdf(file="180527_miR-1_positional_analysis_6mer-m8.pdf", useDingbats=FALSE)
PlotSiteEnrichment("6mer-mmT3", ypos=520)
dev.copy2pdf(file="180527_miR-1_positional_analysis_6mer-mmT3.pdf", useDingbats=FALSE)
PlotSiteEnrichment("5mer-m3.7", xpos=520, ypos=520)
dev.copy2pdf(file="180527_miR-1_positional_analysis_5mer-m3.7.pdf", useDingbats=FALSE)
# PlotSiteEnrichment("6mer-A1", xpos=520, ypos=520)


# kSiteColors["6mer-m8"] <- "khaki3"
kSiteColors["7mer-m8w3"] <- "forestgreen"
kSiteColors["6mer-mmT3"] <- "brown"
kSiteColors["5mer-m3.7"] <- "goldenrod"
dev.new(xpos=1020, ypos=20, height=5, width=6)
par(kPlotParameters)
xmin <- 0.3
xmax <- 400
ymin <- 0.3
ymax <- 500
BlankPlot(log='xy')
xmax <- 100
AddLogAxis(1, label="Concentration")
AddLogAxis(2, label="Enrichment")
points(x_pos, s_6merm8, type="o", col=kSiteColors["6mer-m8"])
# points(x_pos, s_6merm8_bump, type="o", col=kSiteColors["6mer-m8"], lty=3)
points(x_pos, s_6merm8_spike, type="o", col=kSiteColors["6mer-m8"], lty=2)
points(x_pos, s_8mer, type="o", col=kSiteColors["8mer"])
points(x_pos, s_7merm8, type="o", col=kSiteColors["7mer-m8"])
points(x_pos, s_7merA1, type="o", col=kSiteColors["7mer-A1"])
points(x_pos, s_6mer, type="o", col=kSiteColors["6mer"])
points(x_pos, s_7merm8w3, type="o", col=kSiteColors["7mer-m8w3"])
points(x_pos, s_7merm8w3_spike, type="o", col=kSiteColors["7mer-m8w3"], lty=2)
# points(x_pos, s_7merm8w3_bump, type="o", col=kSiteColors["7mer-m8w3"], lty=3)
points(x_pos, s_6mermmT3, type="o", col=kSiteColors["6mer-mmT3"])
points(x_pos, s_6mermmT3_spike, type="o", col=kSiteColors["6mer-mmT3"], lty=2)
points(x_pos, s_5merm3.7, type="o", col=kSiteColors["5mer-m3.7"])
points(x_pos, s_5merm3.7_spike, type="o", col=kSiteColors["5mer-m3.7"], lty=2)
xy <- GetPlotFractionalCoords(1, 0.5, log='xy')
sites_legend <- c("8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-m8", "7mer-m8w3",
                  "6mer-mmT3", "5mer-m3.7")
print(xy)
legend(xy[1], xy[2], legend=sites_legend, col=kSiteColors[sites_legend], pch=19,
       bty="n", xpd=NA)


# points(x_pos, middle_7merm8*mean(right_spike/middle_7merm8), type="o", col="violet")
# print(mean(right_spike/middle_7merm8))

dev.copy2pdf(file="180527_PositionalEnrichmentplot.pdf")
break


GetMirnaCountData <- function(mirna, condition) {
  path <- GetAnalysisPath(mirna, "AGO_purity", condition,
                           "AGO_pur_counts", ext="_new")
  print(path)
  file_counts <- read.table(path, sep="\t", header=TRUE, row.names=1)
  # colnames(file_counts) <- seq(40)
  file_counts
}




colors_out <- c("gray", # 18nt marker
                "pink", #xtr-mIR-427
                "gray80", # 30nt marker
                "purple", # dme-miR-14
                "green", #5p_adapter
                "black", # None
                "forestgreen", #3p adapter
                "red", # red
                "blue") #blue

barcodes <- c("ATCACG",
              "CGATGT",
              "TTAGGC",
              "TGACCA",
              "ACAGTG",
              "GCCAAT",
              "CAGATC",
              "ACTTGA",
              "GATCAG")

sample_types <- c("S1005_I",
                  "S1005_P",
                  "S1005_P",
                  "S1006_I",
                  "S1006_P",
                  "S1006_P",
                  "S1007_I",
                  "S1007_P",
                  "S1007_P")

mirna_sample_types <- rep(c("Na", "miR-1", "miR-155"), 3)

GetAllMirnaSeqData <- function() {
  out_matrix <- matrix(0)
  out <- do.call("cbind", lapply(c("miR-1", "miR-155"), function(mirna) {
    out <- do.call("cbind", lapply(seq(5, 7), function(i_s) {
      out <- do.call("cbind", lapply(c("I", "P"), function(cond) {
        condition <- paste0("S100", i_s, "_", cond)
        counts <- GetMirnaCountData(mirna, condition)
        count_total <- rowSums(counts)
        names(count_total) <- rownames(counts)
        count_total
      }))

      print("first")
      colnames(out) <- paste0("S100", i_s, "_", c("I", "P"))
      # print(out)
      out
    }))
    print("second")
    # print(out)
    out
  }))

  print("third")
  # print(out)
  dev.new(xpos=20, ypos=20, height=5, width=5)
  exclude_rows <- c("18nt_marker", "Unmapped", "30nt_marker")

  out <- out[, c(1, 2, 8, 3, 4, 10, 5, 6, 12)]
  out <- out[!(rownames(out) %in% exclude_rows),]
  out_global <<- out
  exog_rows <- c("miR-1-3p", "miR-155-5p")
  spike_rows <- c("dme-miR-14-5p", "xtr-miR-427")
  exog_df <- out[exog_rows, ]
  spike_df <- out[spike_rows, ]
  endog_df <- out[setdiff(rownames(out), c(exog_rows, spike_rows)), ]
  input_cols <- grep("I", colnames(out))
  max_val <- 10^ceiling(log10(max(spike_df)))
  max_val <- 6000

  out_old <- read.table("140603_WIGTC-HISEQA_C4W8TACXX_expression_matrix.txt",
                        row.names=1, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  order_out_old <- sapply(barcodes, function(barcode) {
    which(out_old[1,] == barcode)
    })

  out_old <- out_old[, order_out_old]
  # plot(c(spike_df[1, ]), c(out_old[4, ]), xlim=c(0, max_val), ylim=c(0, max_val),)
  # points(spike_df[2, ], out_old[6, ], col="red")


  exog_norm <- exog_df/colSums(spike_df)
  endog_norm <- endog_df/colSums(spike_df)


  total_input_endog <- rowSums(endog_norm[, c(1, 4, 7)])

  print(head(total_input_endog))

  row_order <- order(-total_input_endog)

  endog_norm <- endog_norm[row_order,]
  endog_norm <- rbind(endog_norm[1:10, ], colSums(endog_norm[11:nrow(endog_norm), ]))
  rownames(endog_norm)[nrow(endog_norm)] <- "Remaining miRNAs"
  print(endog_norm)

  final_norm_matrix <- rbind(exog_norm, endog_norm)


  colors_out <- c("red", "purple",
                  sprintf("gray%s",
                          floor(seq(50, 70,
                                    length.out=nrow(final_norm_matrix) - 2))))
  print(exog_norm)

  plot(exog_norm[1, ], exog_df[1,])



  dev.new(xpos=520, ypos=20, height=5, width=8)

  barplot(final_norm_matrix, width=0.5, legend=rownames(final_norm_matrix), col=colors_out)
  return(out)
}


  counts_m1_5_I <- GetMirnaCountData("miR-1", "S1005_I")
  counts_m1_6_I <- GetMirnaCountData("miR-1", "S1006_I")
  counts_m1_7_I <- GetMirnaCountData("miR-1", "S1007_I")
  counts_m1_5_P <- GetMirnaCountData("miR-1", "S1005_P")
  counts_m1_6_P <- GetMirnaCountData("miR-1", "S1006_P")
  counts_m1_7_P <- GetMirnaCountData("miR-1", "S1007_P")

  counts_m155_5_I <- GetMirnaCountData("miR-155", "S1005_I")
  counts_m155_6_I <- GetMirnaCountData("miR-155", "S1006_I")
  counts_m155_7_I <- GetMirnaCountData("miR-155", "S1007_I")
  counts_m155_5_P <- GetMirnaCountData("miR-155", "S1005_P")
  counts_m155_6_P <- GetMirnaCountData("miR-155", "S1006_P")
  counts_m155_7_P <- GetMirnaCountData("miR-155", "S1007_P")
# break
data_out <- GetAllMirnaSeqData()
data_out <- data_out[order(-data_out[, 2]),]





# sapply(kMirnas, function(mirna) {
#   if (mirna == "miR-7-23nt") {
#     experiment <- "equilibrium2_nb"
#   } else {
#     experiment <- "equilibrium"
#   }
#   A <- kAgoStock[mirna, "equilibrium"]
#   sXc <- SitesXCounts(mirna, experiment=experiment)
#   dils <- as.numeric(colnames(sXc)[3:(ncol(sXc) - 1)])
#   A_j <- A * dils
#   ratio <- A_j / 100
#   print(ratio)
# })


break


# graphics.off()

# MakeFigureNow <- function(){

# wee_et_al <- read.table("Wee_et_al_data.txt", sep="\t", header=TRUE)

# salomon_et_al <- read.table("salomon_et_al.txt", sep="\t", header=TRUE)

# salomon_et_al <- salomon_et_al[1:11,]
# salomon_et_al$Sequence <- gsub("(?: |-)", salomon_et_al$Sequence, replace="", perl=TRUE)
# salomon_et_al$Sequence <- gsub("U", salomon_et_al$Sequence, replace="T", perl=TRUE)
# salomon_et_al$Sequence <- sapply(salomon_et_al$Sequence, function(i) {
#   substr(i, 92, 92 + 11)
#   })



# sites <- rownames(SitesXCounts("let-7a"))
# sites <- rev(sites[-length(sites)])
# print(sites)



# for (site in sites) {
#   sequence <- GetSiteSeq("let-7a", site)
#   print(sequence)
#   print(grep(sequence, salomon_et_al$Sequence))
# }

# kds_1.4 <- exp(data.frame(fread(paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a",
#                                    "/equilibrium/kds_PAPER/5_12mers_1-4_PAPER",
#                                    "_logmean.txt"), sep="\t"), row.names=1)/log(10))

# kds_2.5 <- exp(data.frame(fread(paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a",
#                                    "/equilibrium/kds_PAPER/5_12mers_2-5_PAPER",
#                                    "_logmean.txt"), sep="\t"), row.names=1)/log(10))
# kds_3.6 <- exp(data.frame(fread(paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a",
#                                    "/equilibrium/kds_PAPER/5_12mers_3-6_PAPER",
#                                    "_logmean.txt"), sep="\t"), row.names=1)/log(10))

# kds_4.7 <- exp(data.frame(fread(paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a",
#                                    "/equilibrium/kds_PAPER/5_12mers_4-7_PAPER",
#                                    "_logmean.txt"), sep="\t"), row.names=1)/log(10))


# kds_5.8 <- exp(data.frame(fread(paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a",
#                                    "/equilibrium/kds_PAPER/5_12mers_5-8_PAPER",
#                                    "_logmean.txt"), sep="\t"), row.names=1)/log(10))


# all_12mer_kds <- list(kds_1.4, kds_2.5, kds_3.6, kds_4.7, kds_5.8)

# salomon_et_al_seed_only <- salomon_et_al[c(2, 5, 6, 7, 8, 9, 10, 11),]

# final_target_kds <- rep(NaN, nrow(salomon_et_al_seed_only))
# names(final_target_kds) <- salomon_et_al_seed_only$Sequence

# print(final_target_kds)
# for (target in salomon_et_al_seed_only$Sequence) {
#   print(target)
#   target_kds <- c()
#   for (i in seq(5)) {
#     kds_inds <- grep(target, rownames(all_12mer_kds[[i]]))
#     print(kds_inds)
#     print(all_12mer_kds[[i]][kds_inds,])
#     target_kds <- c(target_kds, all_12mer_kds[[i]][kds_inds,])
#     final_target_kds[target] <- exp(mean(log(target_kds)))
#   }
# }

# salomon_kds <- salomon_et_al_seed_only[, 7]/salomon_et_al_seed_only[1, 7]
# final_target_kds <- final_target_kds/final_target_kds[1]

# dev.new(xpos=20, ypos=20, height=5, width=5)

# par(kPlotParameters)

# xmin <- 0.1
# xmax <- 1e5
# ymin <- xmin
# ymax <- xmax
# BlankPlot(log='xy')
# AddLogAxis(1, label="Average 12mer Kd")
# AddLogAxis(2, label="Salomon Kd")
# abline(0, 1, lty=2, col="gray")

# colors_temp <- c("purple", "red", "orangered", "orange", "green", "forestgreen", "cyan", "blue")
# points(final_target_kds, salomon_kds, log='xy', xlim=c(1, 1e5), ylim=c(1, 1e5), pch=19, col=colors_temp)

# legend("topleft", legend=c("8mer", "mm2-3", "mm3-4", "mm4-5", "mm5-6", "mm6-7", "6mer-A1", "7mer-A1"),
#        col=colors_temp, pch=19, bty="n")

# xy <- GetPlotFractionalCoords(0.90, 0.95, log='xy')

# # l_salomon_kds <- log10(salomon_kds)

# # fit <- lm()

# AddCorrelationToPlot(log(salomon_kds), log(final_target_kds), xy[1], xy[2], rsquared=TRUE, adj=1)
# dev.copy2pdf(file="temp_now.pdf")
# }


# MakeFigureNow()
# break

Check12merKdsNew <- function(mirna) {
  dev.new(xpos=20, ypos=20, height=12, width=10)
  par(kPlotParameters)
  par(mfrow=c(7, 5))
  for (mirna in kMirnas) {
    if (mirna %in% c("miR-1", "miR-7-23nt")) {
      str.combined <- "_nocombInput"
    } else {
      str.combined <- ""
    }
    if (mirna == "miR-7-23nt") {
      exp <- "equilibrium2_nb"
      str.global <- "_global"
    } else {
      exp <- "equilibrium"
      str.global <- ""
    }

    for (mir_start in seq(5)) {
      print(mir_start)
      if (mirna == "miR-7-23nt") {
            kds_old <- data.frame(fread(paste0("12merfits/miR-7-23nt_equilibrium2_nb_5_12mers_2__nocombInputglobal_", mir_start,
                            "-", mir_start + 3, ".txt")), row.names=1)

            } else if (mir_start == 1) {
        kds_old <- (t(fread(paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                                  "/", exp, "/kds_PAPER/5_12mers_", mir_start,
                                  "-", mir_start + 3, str.combined, "_PAPER_mean.txt"))))*log(10)

        } else {
        kds_old <- t(fread(paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                                  "/", exp, "/kds_PAPER/5_12mers_", mir_start,
                                  "-", mir_start + 3,  str.combined, "_PAPER_mean.txt")))
        rownames(kds_old) <- kds_old[, 1]
        kds_old <- kds_old[, 2, drop=FALSE]
        kds_old_val <- as.numeric(kds_old)
        kds_old_global <- kds_old
        rownames(kds_old)[nrow(kds_old)] <- "AGO"
        rownames(kds_old)[nrow(kds_old) - 1] <- "bg"
        rownames(kds_old)[1:(nrow(kds_old) - 2)] <- paste0(rownames(kds_old)[1:(nrow(kds_old) - 2)], "_Kd")
        kds_old <- data.frame(log(Logistic(kds_old_val, max=100)), row.names=rownames(kds_old))
      }
      kds_new <- data.frame(fread(paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/", exp, "/kds_PAPER/5_12mers_", mir_start, "-", mir_start + 3, str.combined, str.global, "_PAPER_logmean.txt")), row.names=1)/log(10)
      # kds_old <<- kds_old

      # kds_new <<- kds_new
      plot(c(kds_old[1:65537,]- kds_old[65537,]), c(kds_new[1:65537, ] - kds_new[65537,]), xlim=c(-10, 4), ylim=c(-10, 4))
    }
  }
   mirna <- "miR-7-23nt"
   exp <- "equilibrium2_nb" 
      str.combined <- ""
      str.global <- "_global"

    for (mir_start in seq(5)) {
      print(mir_start)
      if (mirna == "miR-7-23nt") {
            kds_old <- data.frame(fread(paste0("12merfits/miR-7-23nt_equilibrium2_nb_5_12mers_2__nocombInputglobal_", mir_start,
                            "-", mir_start + 3, ".txt")), row.names=1)

            } else if (mir_start == 1) {
        kds_old <- (t(fread(paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                                  "/", exp, "/kds_PAPER/5_12mers_", mir_start,
                                  "-", mir_start + 3, str.combined, "_PAPER_mean.txt"))))*log(10)

        } else {
        kds_old <- t(fread(paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                                  "/", exp, "/kds_PAPER/5_12mers_", mir_start,
                                  "-", mir_start + 3,  str.combined, "_PAPER_mean.txt")))
        rownames(kds_old) <- kds_old[, 1]
        kds_old <- kds_old[, 2, drop=FALSE]
        kds_old_val <- as.numeric(kds_old)
        kds_old_global <- kds_old
        rownames(kds_old)[nrow(kds_old)] <- "AGO"
        rownames(kds_old)[nrow(kds_old) - 1] <- "bg"
        rownames(kds_old)[1:(nrow(kds_old) - 2)] <- paste0(rownames(kds_old)[1:(nrow(kds_old) - 2)], "_Kd")
        kds_old <- data.frame(log(Logistic(kds_old_val, max=100)), row.names=rownames(kds_old))
      }
      kds_new <- data.frame(fread(paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/", exp, "/kds_PAPER/5_12mers_", mir_start, "-", mir_start + 3, str.combined, str.global, "_PAPER_logmean.txt")), row.names=1)/log(10)
      # kds_old <<- kds_old

      # kds_new <<- kds_new
      plot(c(kds_old[1:65537,]- kds_old[65537,]), c(kds_new[1:65537, ] - kds_new[65537,]), xlim=c(-10, 4), ylim=c(-10, 4))
    }



}

Check12merKdsNew()
break

Check12merKds <- function() {
  dev.new(xpos=20, ypos=20, height=10, width=15)
  par(kPlotParameters)
  par(mfcol=c(4, 6))
  for (mirna in kMirnas[1:5]) {
    message("____________________")
    print(mirna)
    exp <- "equilibrium"
    if (mirna %in% c("miR-1", "miR-7-23nt")) {
      str.combined <- "_nocombInput"
    } else {
      str.combined <- ""
    }
    mir_start <- 1
    kds_old <- 10^(t(fread(paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                              "/", exp, "/kds_PAPER/5_12mers_", mir_start,
                              "-", mir_start + 3, str.combined, "_PAPER_mean.txt"))))
    

    sites <- rownames(SitesXCounts(mirna, "equilibrium"))
    sites <- rev(sites[-length(sites)])
    print(sites)
    seqs <- sapply(sites, GetSiteSeq, mirna=mirna)
    print(seqs)
    seq_colors <- kSiteColors[sites]
    # print(seq_colors)

    for (mir_start in seq(2, 5)) {
      kds_new <- t(fread(paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                                "/", exp, "/kds_PAPER/5_12mers_", mir_start,
                                "-", mir_start + 3,  str.combined, "_PAPER_mean.txt")))
      rownames(kds_new) <- kds_new[, 1]
      kds_new <- kds_new[, 2, drop=FALSE]
      kds_new_val <- as.numeric(kds_new)
      kds_new_global <- kds_new
      rownames(kds_new)[nrow(kds_new)] <- "AGO"
      rownames(kds_new)[nrow(kds_new) - 1] <- "bg"
      rownames(kds_new)[1:(nrow(kds_new) - 2)] <- paste0(rownames(kds_new)[1:(nrow(kds_new) - 2)], "_Kd")
      kds_new <- data.frame(Logistic(kds_new_val, max=10), row.names=rownames(kds_new))
      shared_inds <- intersect(rownames(kds_old)[1:4^8], rownames(kds_new)[1:4^8])
      shared_inds_global <<- shared_inds
    
      colors_12mers <- rep(ConvertRColortoRGB("gray", alpha=0.01), nrow(kds_old))
      for (i in 1:length(seq_colors)) {

        seq_i = seqs[i]
        # print(seq_i)
        seq_i_global <<- seq_i
        inds <- grep(seq_i, shared_inds)
        # print(seq_colors[i])
        # print(length(inds))
        colors_12mers[inds] <- ConvertRColortoRGB(seq_colors[i], alpha=0.5)
      }

      xmin <- 1e-6
      xmax <- 30
      ymin <- xmin
      ymax <- xmax
      BlankPlot(log='xy')
      AddLogAxis(1, label=mir_start - 1)
      AddLogAxis(2, label=mir_start)
      xy <- GetPlotFractionalCoords(0.025, 0.95, log='xy')
      text(xy[1], xy[2], mirna, adj=0)
      print(tail(kds_old))
      print(tail(kds_new))
      print(kds_old["None_Kd",])
      print(kds_new["None_Kd",])
      points(kds_old[shared_inds, ]/(kds_old["None_Kd",]),
             kds_new[shared_inds,]/(kds_new["None_Kd",]), col=colors_12mers)
      kds_old <- kds_new
    }
  }

  mirna <- "miR-7-23nt"
  exp <- "equilibrium2_nb"
  str.combined <- "_nocombInput"
      # str.combined <- ""
  mir_start <- 1
  # "12merfits/miR-7-23nt_equilibrium2_nb_5_12mers_2__nocombInputglobal_1-4.txt"
  kds_old <- data.frame(fread(paste0("12merfits/miR-7-23nt_equilibrium2_nb_5_12mers_2__nocombInputglobal_", mir_start,
                            "-", mir_start + 3, ".txt")), row.names=1)
  kds_old_global <<- kds_old
  sites <- rownames(SitesXCounts("miR-7-23nt", "equilibrium2_nb"))
  sites <- rev(sites[-length(sites)])
  seqs <- sapply(sites, GetSiteSeq, mirna="miR-7-23nt")
  print(seqs)
  seq_colors <- kSiteColors[sites]
  print(seq_colors)



  for (mir_start in seq(2, 5)) {
    kds_new <- data.frame(fread(paste0("12merfits/miR-7-23nt_equilibrium2_nb_5_12mers_2__nocombInputglobal_", mir_start,
                            "-", mir_start + 3, ".txt")), row.names=1)
    kds_new_global <<- kds_new
    shared_inds <- intersect(rownames(kds_old)[1:4^8], rownames(kds_new)[1:4^8])
    print(shared_inds)
    colors_12mers <- rep(ConvertRColortoRGB("gray", alpha=0.01), nrow(kds_old))
    for (i in 1:length(seq_colors)) {
      seq_i = seqs[i]
      inds <- grep(seq_i, shared_inds)
      colors_12mers[inds] <- ConvertRColortoRGB(seq_colors[i], alpha=0.5)
    }
    xmin <- 1e-4
    xmax <- 10
    ymin <- xmin
    ymax <- xmax
    BlankPlot(log='xy')
    AddLogAxis(1, label=mir_start - 1)
    AddLogAxis(2, label=mir_start)
    print(tail(kds_old))
    print(tail(kds_new))
    print(kds_old["None_Kd",])
    print(kds_new["None_Kd",])
    points(exp(kds_old[shared_inds, ])/exp(kds_old["None_Kd",]),
           exp(kds_new[shared_inds,])/exp(kds_new["None_Kd",]), col=colors_12mers)
    kds_old <- kds_new
  }


}

CheckMiR712merKds <- function() {
  dev.new(xpos=20, ypos=20, height=8, width=8)
  par(kPlotParameters)
  par(mfrow=c(2, 2))
  mirna <- "miR-7-23nt"
  exp <- "equilibrium2_nb"
  str.combined <- "_nocombInput"
      # str.combined <- ""
  mir_start <- 1
  "12merfits/miR-7-23nt_equilibrium2_nb_5_12mers_2__nocombInputglobal_1-4.txt"
  kds_old <- data.frame(fread(paste0("12merfits/miR-7-23nt_equilibrium2_nb_5_12mers_2__nocombInputglobal_", mir_start,
                            "-", mir_start + 3, ".txt")), row.names=1)
  kds_old_global <<- kds_old
  sites <- rownames(SitesXCounts("miR-7-23nt", "equilibrium2_nb"))
  sites <- rev(sites[-length(sites)])
  seqs <- sapply(sites, GetSiteSeq, mirna="miR-7-23nt")
  print(seqs)
  seq_colors <- kSiteColors[sites]
  print(seq_colors)



  for (mir_start in seq(2, 5)) {
    kds_new <- data.frame(fread(paste0("12merfits/miR-7-23nt_equilibrium2_nb_5_12mers_2__nocombInputglobal_", mir_start,
                            "-", mir_start + 3, ".txt")), row.names=1)
    # rownames(kds_new) <- kds_new[, 1]
    # kds_new <- kds_new[, 2, drop=FALSE]
    # kds_new_val <- as.numeric(kds_new)
    # rownames(kds_new)[nrow(kds_new)] <- "AGO"
    # rownames(kds_new)[nrow(kds_new) - 1] <- "bg"
    # rownames(kds_new)[1:(nrow(kds_new) - 2)] <- paste0(rownames(kds_new)[1:(nrow(kds_new) - 2)], "_Kd")
    # kds_new <- data.frame(exp(kds_new_val), row.names=rownames(kds_new))
    kds_new_global <<- kds_new
    shared_inds <- intersect(rownames(kds_old)[1:4^8], rownames(kds_new)[1:4^8])
    
    colors_12mers <- rep(ConvertRColortoRGB("gray", alpha=0.01), nrow(kds_old))
    for (i in 1:length(seq_colors)) {
      seq_i = seqs[i]
      inds <- grep(seq_i, shared_inds)
      colors_12mers[inds] <- ConvertRColortoRGB(seq_colors[i], alpha=0.5)
    }
    xmin <- 1e-4
    xmax <- 10
    ymin <- xmin
    ymax <- xmax
    BlankPlot(log='xy')
    AddLogAxis(1, label=mir_start - 1)
    AddLogAxis(2, label=mir_start)
    print(tail(kds_old))
    print(tail(kds_new))
    print(kds_old["None_Kd",])
    print(kds_new["None_Kd",])
    points(exp(kds_old[shared_inds, ])/exp(kds_old["None_Kd",]),
           exp(kds_new[shared_inds,])/exp(kds_new["None_Kd",]), col=colors_12mers)
    kds_old <- kds_new
  }
}


Check12merKds()

break
#   kds_1.4 <- t(fread("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/kds_PAPER/5_12mers_1-4_PAPER_mean.txt"))
# kds_2.5 <- t(fread("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/kds_PAPER/5_12mers_2-5_PAPER_mean.txt"))
# rownames(kds_2.5) <- kds_2.5[, 1]
# kds_2.5 <- kds_2.5[, 2, drop=FALSE]
# kds_2.5_val <- as.numeric(kds_2.5)
# rownames(kds_2.5)[nrow(kds_2.5)] <- "AGO"
# rownames(kds_2.5)[nrow(kds_2.5) - 1] <- "bg"
# rownames(kds_2.5)[1:(nrow(kds_2.5) - 2)] <- paste0(rownames(kds_2.5)[1:(nrow(kds_2.5) - 2)], "_Kd")
# kds_2.5 <- data.frame(kds_2.5_val, row.names=rownames(kds_2.5))

# kds_3.6 <- t(fread("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/kds_PAPER/5_12mers_3-6_PAPER_mean.txt"))
# rownames(kds_3.6) <- kds_3.6[, 1]
# kds_3.6 <- kds_3.6[, 2, drop=FALSE]
# kds_3.6_val <- as.numeric(kds_3.6)
# rownames(kds_3.6)[nrow(kds_3.6)] <- "AGO"
# rownames(kds_3.6)[nrow(kds_3.6) - 1] <- "bg"
# rownames(kds_3.6)[1:(nrow(kds_3.6) - 2)] <- paste0(rownames(kds_3.6)[1:(nrow(kds_3.6) - 2)], "_Kd")
# kds_3.6 <- data.frame(kds_3.6_val, row.names=rownames(kds_3.6))


kds_4.7 <- t(fread("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/kds_PAPER/5_12mers_4-7_PAPER_mean.txt"))
rownames(kds_4.7) <- kds_4.7[, 1]
kds_4.7 <- kds_4.7[, 2, drop=FALSE]
kds_4.7_val <- as.numeric(kds_4.7)
rownames(kds_4.7)[nrow(kds_4.7)] <- "AGO"
rownames(kds_4.7)[nrow(kds_4.7) - 1] <- "bg"
rownames(kds_4.7)[1:(nrow(kds_4.7) - 2)] <- paste0(rownames(kds_4.7)[1:(nrow(kds_4.7) - 2)], "_Kd")
kds_4.7 <- data.frame(kds_4.7_val, row.names=rownames(kds_4.7))


kds_5.8 <- t(fread("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/kds_PAPER/5_12mers_5-8_PAPER_mean.txt"))
rownames(kds_5.8) <- kds_5.8[, 1]
kds_5.8 <- kds_5.8[, 2, drop=FALSE]
kds_5.8_val <- as.numeric(kds_5.8)
rownames(kds_5.8)[nrow(kds_5.8)] <- "AGO"
rownames(kds_5.8)[nrow(kds_5.8) - 1] <- "bg"
rownames(kds_5.8)[1:(nrow(kds_5.8) - 2)] <- paste0(rownames(kds_5.8)[1:(nrow(kds_5.8) - 2)], "_Kd")
kds_5.8 <- data.frame(kds_5.8_val, row.names=rownames(kds_5.8))





# break
# kds_3.6 <- data.frame(t(fread("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/kds_PAPER/5_12mers_3-6_PAPER_mean.txt")), row.names=1)
# kds_4.7 <- data.frame(t(fread("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/kds_PAPER/5_12mers_4-7_PAPER_mean.txt")), row.names=1)
# kds_5.8 <- data.frame(t(fread("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/kds_PAPER/5_12mers_5-8_PAPER_mean.txt")), row.names=1)


# kds_1.4 <- data.frame(fread("12merfits/miR-7-23nt_equilibrium2_nb_5_12mers_2__nocombInputglobal_1-4.txt"), row.names=1)
# kds_2.5 <- data.frame(fread("12merfits/miR-7-23nt_equilibrium2_nb_5_12mers_2__nocombInputglobal_2-5.txt"), row.names=1)
# kds_3.6 <- data.frame(fread("12merfits/miR-7-23nt_equilibrium2_nb_5_12mers_2__nocombInputglobal_3-6.txt"), row.names=1)
# kds_4.7 <- data.frame(fread("12merfits/miR-7-23nt_equilibrium2_nb_5_12mers_2__nocombInputglobal_4-7.txt"), row.names=1)
# kds_5.8 <- data.frame(fread("12merfits/miR-7-23nt_equilibrium2_nb_5_12mers_2__nocombInputglobal_5-8.txt"), row.names=1)


share_1.2 <- intersect(rownames(kds_1.4), rownames(kds_2.5))
share_2.3 <- intersect(rownames(kds_2.5), rownames(kds_3.6))
share_3.4 <- intersect(rownames(kds_3.6), rownames(kds_4.7))
share_4.5 <- intersect(rownames(kds_4.7), rownames(kds_5.8))

test_1.2 <- cbind(kds_1.4[share_1.2, ], kds_2.5[share_1.2,])

dev.new(xpos=20, ypos=20, height=5, width=5)
plot(kds_1.4[share_1.2, ], kds_2.5[share_1.2,])
dev.new(xpos=520, ypos=20, height=5, width=5)
plot(kds_2.5[share_2.3, ], kds_3.6[share_2.3,])


break

CheckGlobalFit <- function() {
  par(kPlotParameters)
  par(mfcol=c(3, 5))
  colors_all <- c("purple", "blue", "red")
  names(colors_all) <- c("miR-7-23nt", "miR-7-24nt", "miR-7-25nt")
  color <- ConvertRColortoRGB("black", alpha=0.2)
  sapply(1:5, function(mirna.start) {
    print(mirna.start)
    globalname <- paste0("12merfits/miR-7-23nt_equilibrium2_nb",
                          "_5_12mers_2_global_", mirna.start, "-",
                          as.numeric(mirna.start) + 3, ".txt")
    print(globalname)
    kds_global <- data.frame(fread(globalname), row.names=1)
    print(head(kds_global))
    print(head(sort(exp(kds_global[,]), decreasing=TRUE)))
    sapply(c("miR-7-23nt", "miR-7-24nt", "miR-7-25nt"), function(mirna) {
      print(mirna)
      separatename <- paste0("12merfits/", mirna, "_equilibrium2_nb_5_12mers_", mirna.start,
                             "__nocombInputfull_temp.txt")
      print(separatename)
      separatename2 <- "12merfits/miR-7-23nt_equilibrium2_nb_5_12mers_1__nocombInputfull_temp.txt"
      print(separatename2)
      kds_separate <- data.frame(fread(separatename), row.names=1)
      shared <- intersect(rownames(kds_global), rownames(kds_separate))
      print(exp(head(kds_global[shared, ])))
      print(head(kds_separate[shared, ]))
      plot(exp(kds_global[shared, ]), kds_separate[shared, ],
           log='xy', xlim=c(1e-5, 2e4), ylim=c(1e-5, 2e4), col=ConvertRColortoRGB(colors_all[mirna], alpha=0.2))
    })
  })
  mirna.start <- 2
  globalname <- paste0("12merfits/miR-7-23nt_equilibrium2_nb",
                          "_5_12mers_2_global_", mirna.start, "-",
                          as.numeric(mirna.start) + 3, ".txt")
  kds_global <- data.frame(fread(globalname), row.names=1)

  sapply(3:5, function(mirna.start) {
  globalname <- paste0("12merfits/miR-7-23nt_equilibrium2_nb",
                          "_5_12mers_2_global_", mirna.start, "-",
                          as.numeric(mirna.start) + 3, ".txt")
  kds_global2 <- data.frame(fread(globalname), row.names=1)
  shared <- intersect(rownames(kds_global), rownames(kds_global2))
  xmin <- 1e-5
  xmax <- 2e4
  ymin <- xmin
  ymax <- xmax
  BlankPlot(log='xy')
  AddLogAxis(1, label=mirna.start-1)
  AddLogAxis(2, label=mirna.start)
  abline(0, 1, lty=2)
  points(exp(kds_global[shared, ]), exp(kds_global2[shared, ]), col=color)
  kds_global <<- kds_global2
})

}

CheckGlobalFit()




break
library(numDeriv)
library(Rfast)
library(wCorr)

hela_df <- data.frame(fread("RepressionData/Lin-Shi_transfection_data/counts.txt",
                    sep="\t", header=FALSE), row.names=1)


head(hela_df)



hela_logtpm_df <- data.frame(fread("/lab/bartel4_ata/kathyl/RNA_Seq/analysis/data/final/log_tpm.txt",
                                   sep="\t", header=TRUE), row.names=1)
hela_logtpm_normed_df <- data.frame(fread("/lab/bartel4_ata/kathyl/RNA_Seq/analysis/data/final/log_tpm_normed.txt",
                                   sep="\t", header=TRUE), row.names=1)

hela_merged_df <- data.frame(fread("/lab/bartel4_ata/kathyl/RNA_Seq/analysis/data/final/merged.txt",
                                   sep="\t", header=TRUE), row.names=1)


# hela_tpm_filtered_sm_df <- data.frame(fread("RepressionData/Lin-Shi_transfection_data/tpm_filtered.txt",
#                                       sep="\t", header=TRUE), row.names=1)

# hela_utr_sites_filtered_sm_df <- data.frame(fread("RepressionData/Lin-Shi_transfection_data/utr_sites_filtered.txt",
#                                       sep="\t", header=TRUE), row.names=1)

# logtpm <- log(hela_tpm_filtered_sm_df)


# tpm_kl <- exp(hela_logtpm_df)
# tpm_sm <- hela_tpm_filtered_sm_df
# # My tpms are the same as kathys CHECK

# print(colnames(tpm_kl))
# colnames(tpm_kl) <- sapply(colnames(tpm_kl), function(name) {
#   new_name <- grep(name, colnames(tpm_sm), value=TRUE)
#   return(new_name)
# })
# print(colnames(tpm_kl))
# tpm_kl <- tpm_kl[,colnames(tpm_sm)]
# print(tpm_kl[1:5,1:5])
# print(tpm_sm[1:5,1:5])
# print(dim(tpm_kl))
# print(dim(tpm_sm))
# print((tpm_kl/tpm_sm)[1:5, 1:5])

# batches <- sapply(colnames(tpm_sm), function(name) {
#                   strsplit(name, split="_")[[1]][3]
#                   })
# print(batches)

# batch_matrix <- rep(batches, each=nrow(tpm_sm))
# print(length(batch_matrix))
# mirnas <- sapply(colnames(tpm_sm), function(name) {
#                   print(name)
#                   strsplit(name, split="_")[[1]][1]
#                  })
# print(mirnas)

# mirna_matrix <- matrix("nosite", nrow=nrow(tpm_sm), ncol=ncol(tpm_sm))

# for (i in seq(ncol(tpm_sm))) {
#   # print(i)
#   # print(as.vector(mirnas[i]))
#   # print(mirna_matrix[1:5, 1:3])
#   inds <- which(hela_utr_sites_filtered_sm_df[,i] == 1)
#   mirna_matrix[which(hela_utr_sites_filtered_sm_df[,i] == 1), i] <- as.vector(mirnas[i])
# }

# lm_df <- data.frame(tpm=unlist(tpm_sm), batch=c(batch_matrix), mirna=c(mirna_matrix), gene=rep(rownames(tpm_sm), ncol(tpm_sm)))
# print(lm_df[1:100,])
# # lm_df <- lm_df[1:10,]

# # gene_design <- matrix(0, nrow=nrow(lm_df)*ncol(tpm_sm),ncol=nrow(lm_df))
# # for (i in seq(ncol(lm_df))) {
# #   for (j in seq(nrow(lm_df))) {
# #     gene_design[((i - 1)*nrow(lm_df) + j), j] <- 1
# #   }
# # }

# # mirna_design <- matrix(0, nrow=nrow())
# # break

# # # mirna_design <- matrix(0, nrow=length(mirnas) + 1, ncol=ncol=tpm)

# # break

# fits <- lm(tpm ~ batch + mirna + gene, data=lm_df)

# saveRDS(fits, "new_model.rds")



# tpm_scaled <- t(t(logtpm) - colMeans(logtpm))
# tpm_scaled <- t(t(tpm_scaled) - colMeans(tpm_scaled))
# tpm_scaled <- apply(tpm_scaled, 2, function(col) {
#   col /sd(col)
#   })

# print(colMeans(tpm_scaled))
# print(tpm_scaled[1:5, 1:5])


# tpm.pca <- prcomp(tpm_scaled, retx=FALSE)



# pcs <- tpm.pca$rotation


# mirnas <- matrix(unlist(strsplit(rownames(pcs), split="_")),
#                  nrow=3)
# colors_mirna = rainbow(length(unique(mirnas[1,])))
# colors_batch = rainbow(length(unique(mirnas[3,])))
# names(colors_mirna) <- unique(mirnas[1,])
# names(colors_batch) <- unique(mirnas[3,])

# print(mirnas)
# par(mfrow=c(1, 2))
# plot(pcs[,1], pcs[,2], col=colors_mirna[mirnas[1,]])
# plot(pcs[,1], pcs[,2], col=colors_batch[mirnas[3,]])

# break


PlotBindingMass <- function(mirna, xpos=20, ypos=20) {
  if (mirna %in% c("miR-1", "miR-7-23nt")) {
    combined <- FALSE
  } else {
    combined <- TRUE
  }
  if (mirna == "miR-7-23nt") {
    experiment <- "equilibrium2_nb"
    fixed <- TRUE
  } else {
    fixed <- FALSE
  }
  trial <- SubfunctionCall(GetBindingCapacityKdRepTuples)
  dev.new(xpos = xpos, ypos=ypos, height=4, width=8)


  xmin <- 1e-4
  xmax <- 1
  ymin <- -1.4
  ymax <- 0.4
  par(kPlotParameters)
  par(mfrow=c(1, 2))
  BlankPlot(log='x', inv='x')
  AddLogAxis(1, label="Kd")
  AddLinearAxis(2, tick.space=0.1, label.space=0.2, label="log2fc")
  arrows(trial[, 2], trial[, 3] + trial[, 4], trial[, 2], trial[, 3] - trial[, 4],
         length=0.05*par()$cex, lwd=1.5*par()$cex, angle=90, code=3)

  points(trial[, 2], trial[, 3], col=kSiteColors[rownames(trial)])

  xmin <- 1e-2
  xmax <- 1e2
  ymin <- -1.4
  ymax <- 0.4
  par(kPlotParameters)
  BlankPlot(log='x')
  AddLogAxis(1, label="Binding Mass")
  AddLinearAxis(2, tick.space=0.1, label.space=0.2, label="log2fc")
  arrows(trial[, 1], trial[, 3] + trial[, 4], trial[, 1], trial[, 3] - trial[, 4],
         length=0.05*par()$cex, lwd=1.5*par()$cex, angle=90, code=3)

  points(trial[, 1], trial[, 3], col=kSiteColors[rownames(trial)])

  points(trial[, 1], trial[, 3], col=kSiteColors[rownames(trial)])
}

GetAllBindingPairs <- function() {
  test_m1 <- GetBindingCapacityKdRepTuples("miR-1", combined=TRUE)
  test_l7 <- GetBindingCapacityKdRepTuples("let-7a")
  test_m155 <- GetBindingCapacityKdRepTuples("miR-155")
  test_m124 <- GetBindingCapacityKdRepTuples("miR-124")
  test_ls6 <- GetBindingCapacityKdRepTuples("lsy-6")
  test_m7 <- GetBindingCapacityKdRepTuples("miR-7-23nt",
                                        experiment="equilibrium2_nb",
                                        fixed=TRUE)
  out <- do.call(rbind, list(test_m1, test_l7, test_m155, test_m124, test_ls6,
                            test_m7))
  out
}

pairs_all <- GetAllBindingPairs()
pairs_all <- pairs_all[order(pairs_all[, 1]), ]

data_start <- data.frame(site=rep("bg", nrow(pairs_all)), fc=pairs_all[, 3],
                          sem=pairs_all[, 4], stringsAsFactors=FALSE)
colnames(data_start) <- c("site", "fc", "fc_sem")
out_values <- rep(NA, nrow(data_start) - 1)
l <- sort(pairs_all[, 1])
k <- sort(pairs_all[, 2])
l_cutoffs <- exp(seq(-13, -6, length.out=40))
k_cutoffs <- exp(seq(-7, 0, length.out=40))

cor_mat <- matrix(NaN, nrow=length(l_cutoffs), ncol=length(k_cutoffs),
                  dimnames=list(l_cutoffs, k_cutoffs))
print(cor_mat)


for (i_l in seq(length(l_cutoffs))) {
  print(i_l)
  for (i_k in seq(length(k_cutoffs))) {
    # print(i_k)
    # print("l cutoff:")
    # print(l_cutoffs[i_l])
    # print("k cutoff")
    # print(k_cutoffs[i_k])
    data <- subset(pairs_all, l >= l_cutoffs[i_l] & k <= k_cutoffs[i_k])
    # print(data)
    r_2 <- weightedCorr(data[,3], log(data[, 2]), weights=1/data[, 4])^2
    # print(r_2)
    cor_mat[i_l, i_k] <- r_2
  }
}

xmin <- l_cutoffs[1]
xmax <- l_cutoffs[length(k_cutoffs)]
ymin <- k_cutoffs[1]
ymax <- k_cutoffs[length(k_cutoffs)]

BlankPlot(log='xy')
x_fold <- (l_cutoffs[2]/l_cutoffs[1])^(1/2)
y_fold <- (k_cutoffs[2]/k_cutoffs[1])^(1/2)
print(x_fold)
print(y_fold)

cor_mat_ints <- floor(cor_mat * 100)

colors <- rainbow(100)

cor_mat_min <- min(cor_mat, na.rm=TRUE)
cor_mat_max <- max(cor_mat, na.rm=TRUE)

print(cor_mat_min)
print(cor_mat_max)

colors_vec <- rep("gray", length(c(cor_mat_ints)))

colors_vec[which(!(is.na(cor_mat_ints)))] <- colors[c(cor_mat_ints)[!(is.na(c(cor_mat_ints)))]]
print(colors_vec)

rect(xleft=rep(l_cutoffs/x_fold, each=ncol(cor_mat)),
     ybottom=rep(k_cutoffs/y_fold, nrow(cor_mat)),
     xright=rep(l_cutoffs*x_fold, each=ncol(cor_mat)),
     ytop=rep(k_cutoffs*y_fold, nrow(cor_mat)), col=colors_vec,
     lwd=0)

# image(t(cor_mat))
break

cutoffs <- (pairs_all[-nrow(pairs_all), 1] + pairs_all[-1, 1])/2
  dev.new(xpos = 20, ypos=20, height=4, width=4)


  xmin <- 1e-4
  xmax <- 1
  ymin <- -1.4
  ymax <- 0.4
  par(kPlotParameters)
  par(mfrow=c(1, 2))
  BlankPlot(log='x', inv='x')
  AddLogAxis(1, label="Kd")
  AddLinearAxis(2, tick.space=0.1, label.space=0.2, label="log2fc")


points(pairs_all[-nrow(pairs_all), 1], out_values, type="o")
break

ModelPreds <- function(m, C) {
  kds <- pairs_all[, 2]
  out <- -log(1 + m/(m + kds)*C)
}

m <- 1
C <- 1


repression_matrix <- fread(sprintf("/lab/solexa_bartel/mcgeary/AgoRBNS/%s/repression_hela_cs/lin_model_df/%s.txt",
                                     "miR-1", "paper"), sep="\t")

fc <- repression_matrix[["fc"]]/log(2)

site_occurence <- repression_matrix[,2:(ncol(repression_matrix) - 1), drop=FALSE]

kds <- EquilPars("miR-1", combined=FALSE)
rownames(kds) <- gsub("_Kd", "", rownames(kds))

kds <- kds[1:(nrow(kds) - 3), ]

site_names <- rownames(kds)
kds <- kds$Mean
names(kds) <- site_names

FullModel <- function(m, C) {
  sites <- colnames(site_occurence)
  # print(sites)
  # print(names(kds))
  full_rows <- nrow(site_occurence)
  # full_rows <- 10
  fcs <- -log(1 + C*rowSums(site_occurence*m/(m + kds)), 2)
  fcs
} 


# FullModel(1, 1)

PlotFullModel <- function(m, C, rows=lenth(fc)) {
  xmin <- -5.4
xmax <- 2
ymin <- -5.4
ymax <- 2
par(kPlotParameters)
BlankPlot()
AddLinearAxis(1, tick.space=0.1, label.space=0.2, label="pred. log2fc")
AddLinearAxis(2, tick.space=0.1, label.space=0.2, label="log2fc")

  fc_model<- SubfunctionCall(FullModel)
  points(fc_model, fc, pch=20, col=ConvertRColortoRGB("black", alpha=0.3))
  xy <- GetPlotFractionalCoords(0.5, 0.5)
  text(xy[1], xy[2], cor(fc_model, fc)^2)
}

FullModelResidual <- function(pars) {
  pars <- exp(pars)
  m <- pars[1]
  C <- pars[2]
  fc_model <- SubfunctionCall(FullModel)
  if (tick%%10 == 0) {
    SubfunctionCall(PlotFullModel)  
  }
  tick <<- tick + 1

  sum((fc_model - fc)^2)
}

pars.init <- c(0, 0)
tick <- 0
solution <- optim(pars.init, FullModelResidual)

pars <- exp(solution$par)
m <- pars[1]
C <- pars[2]
PlotFullModel(m, C)
break

PlotModelTrend <- function(m, C) {
  xmin <- 1e-4
xmax <- 1
ymin <- -1.4
ymax <- 0.4
par(kPlotParameters)
BlankPlot(log='x', inv='x')
AddLogAxis(1, label="Kd")
AddLinearAxis(2, tick.space=0.1, label.space=0.2, label="log2fc")

  y <- SubfunctionCall(ModelPreds)
  points(pairs_all[, 2], pairs_all[, 3], col=kSiteColors[rownames(pairs_all)])
  lines(sort(pairs_all[, 2]), y[order(pairs_all[, 2])], col="black")

}

ModelResidual <- function(pars) {
  m <- exp(pars[1])
  C <- exp(pars[2])
  if (tick%%10 == 0) {
    SubfunctionCall(PlotModelTrend)  
  }
  tick <<- tick + 1
  y <- SubfunctionCall(ModelPreds)
  sum((y - pairs_all[,3])^2)
}

tick <- 0
solution <- optim(c(1, 1), ModelResidual)

pars <- exp(solution$par)
m <- pars[1]
C <- pars[2]


PlotModelTrend(m, C)

break
test <- GetBindingCapacityResidual(1, "let-7a", combined=TRUE)

GetFittedValues <- function(pairs) {
  if (nrow(pairs) == 0) {
    pairs_global <<- pairs
    return(cbind(pairs_global, c()))
  }
  # print('full:')
  # print(pairs)
  # print("done full")
  # print(pairs)
  # print("done full")
  linmodel <- lm(pairs[, 2] ~ log(pairs[, 1]))
  vals <- linmodel$fitted.values
  cbind(pairs[, 1], pairs[, 2], vals)
}




GetTotalBindingResidual <- function(cutoff, full_r=FALSE) {
  test_m1 <- GetFittedValues(GetBindingCapacityResidual(cutoff,
                                                        "miR-1", combined=TRUE))
  test_m1 <<- test_m1

  test_l7 <- GetFittedValues(GetBindingCapacityResidual(cutoff, "let-7a"))
  test_m155 <- GetFittedValues(GetBindingCapacityResidual(cutoff, "miR-155"))
  test_m124 <- GetFittedValues(GetBindingCapacityResidual(cutoff, "miR-124"))
  test_ls6 <- GetFittedValues(GetBindingCapacityResidual(cutoff, "lsy-6"))
  test_m7 <- GetFittedValues(GetBindingCapacityResidual(cutoff, "miR-7-23nt",
                                                        experiment="equilibrium2_nb",
                                                        fixed=TRUE))
  all_pairs <- do.call(rbind, list(test_m1, test_l7, test_m155, test_m124, test_ls6,
                            test_m7))

  # all_pairs <<- all_pairs
  # xmin <- 1e-4
  # xmax <- 1
  # ymin <- -1.4
  # ymax <- 0.4

  # dev.set(2)
  # par(kPlotParameters)
  # BlankPlot(log="x", inv="x")

  # AddLogAxis(1, label="Kd")
  # AddLinearAxis(2, tick.space=0.1, label.space=0.2, label="log2fc")

  # points(all_pairs[, 1], all_pairs[, 2], log='x',
  #        col=kSiteColors[rownames(all_pairs)])

  # xmin <- -1.4
  # xmax <- 0.4
  # ymin <- -1.4
  # ymax <- 0.4

  # dev.set(3)
  # par(kPlotParameters)
  # BlankPlot(inv="x")

  # AddLinearAxis(1, tick.space=0.1, label.space=0.2, label="log2fc")
  # AddLinearAxis(2, tick.space=0.1, label.space=0.2, label="log2fc")
  # abline(0, 1, lty=2)
  # points(all_pairs[, 3], all_pairs[, 2],
  #        col=kSiteColors[rownames(all_pairs)])
  r2_all <- cor(log(all_pairs[, 1]), all_pairs[, 2])^2
  r2_each <- cor(all_pairs[, 2], all_pairs[, 3])^2
  # print(r2_all)
  # print(r2_each)
  if (full_r) {
    return(r2_all)
  } else {
    return(r2_each)
  }
}


dev.new(xpos=20, ypos=20, height=5, width=5)
par(kPlotParameters)

xmin <- 1e-2
xmax <- 1e2
ymin <- 0
ymax <- 1
BlankPlot(log="x")
AddLogAxis(1, label="binding mass")
AddLinearAxis(2, tick.space=0.1, label.space=0.2, label=expression(r^2))

course_range <- 10^seq(log10(xmin), log10(xmax), length.out=20)

cutoffs_each <- sapply(course_range, GetTotalBindingResidual)
cutoffs_all <- sapply(course_range, GetTotalBindingResidual,
                      full_r=TRUE)



lines(course_range, cutoffs_each, type="o", xpd=NA)
lines(course_range, cutoffs_all, type = "o", col="red", xpd=NA)

check <- GetTotalBindingResidual(1)
print(check)
break

dev.new(xpos=20, ypos=20, height=5, width=5)
par(kPlotParameters)

segments(xmin, 0, xmax, 0, col="gray")

xy <- GetPlotFractionalCoords(0.05, 0.05, log='x', inv='x')
AddCorrelationToPlot(log(test[, 1]), test[, 2], xy[1], xy[2], rsquared=TRUE)
points(test[,1], test[,2], col=kSiteColors[rownames(test)])

break





# merged <- fread("/lab/bartel4_ata/kathyl/RNA_Seq/analysis/data/final/merged.txt")

# bad <- 0
# for (i in seq(nrow(merged))) {
#   UTR_merged = merged[["sequence"]][i]
#   UTR_alt = UTRs[merged[i][[1]]]
#   if (UTR_merged != UTR_alt) {
#     bad <- bad + 1
#   }
# }
# print(bad)



# GetCanonicalRepression <- function(mirna) {
#   kds <- EquilPars(mirna, "equilibrium", combined=FALSE)
#   rownames(kds) <- gsub("_Kd", "", rownames(kds))
#   print(kds)
#   seqs <- rev(sapply(kSeedSites, GetSiteSeq, mirna=mirna))
#   print(seqs)
#   mirna <- gsub("M", "m", mirna)
#   mirna <- gsub("R", "r", mirna)
#   mirna <- gsub("-", "", mirna)
#   print(mirna)
#   mRNA_sites <- rep("nosite", nrow(merged))
#   for (i in seqs) {
#     inds <- grep(i, merged[["sequence"]])
#     mRNA_sites[inds] <- i
#   }
#   nosite_inds <- which(mRNA_sites == "nosite")
#   log2_fc <- merged[[mirna]]/log(2) - merged[["nosite"]]/log(2)
#   print(head(log2_fc))
#   output <- sapply(c(seqs, "nosite"), function(site) {
#     print(site)
#     data <- log2_fc[which(mRNA_sites == site)]
#     return(mean(data, na.rm=TRUE))
#   })
#   names(output)[length(output)] <- "None"
#   output <- output[c(kSeedSites, "None")]
#   print(output)
#   plot(kds[names(output),]$Mean, output, log='x')
# }

# GetCanonicalRepression("miR-1")


break
# InitializeEquilSitePars <- function(sXc, combined=TRUE) {
#   l <- SubfunctionCall(GetInputEquil)
#   # print(l)
#   data <- SubfunctionCall(GetDataEquil)
#   # print(data)
#   kds <- log10(Norm(l)/Norm(rowSums(data)))/2
#   # kds <- kds*0
#   kds <- kds - kds[length(kds)]
#   names(kds) <- paste0(rownames(sXc), "_Kd")
#   pars <- c(kds, bg=-4, AGO=1)
#   random_pars <- rnorm(length(pars), 0, 3)
#   names(random_pars) <- names(pars)
#   random_pars["None_Kd"] <- 0
#   pars <- random_pars
#   return(pars)
# }

sXc <- SitesXCounts("miR-1", "equilibrium")
sXc_1s <- SingleSitesXCounts("miR-1", "equilibrium")
sXc_2s <- DoubleSitesXCounts("miR-1", "equilibrium")

sXc_pos_4 <- GetPositionalSites("miR-1", "equilibrium", 4, single=TRUE)

sXc_totals <- rowSums(sXc_pos_4)[rownames(sXc)[-nrow(sXc)]]

plot(sXc_totals, sXc_1s[1:(nrow(sXc_1s) - 1), "4"])
print(cbind(sXc_totals, sXc_1s[1:(nrow(sXc_1s) - 1), "4"]))
abline(0, 1, lty=2)

break


temp1_I <- GetPositionalSites("miR-155", "equilibrium", "I")
temp2_I <- GetPositionalSites("miR-155", "equilibrium", "I", single=TRUE)


temp1 <- GetPositionalSites("miR-155", "equilibrium", 4)
temp2 <- GetPositionalSites("miR-155", "equilibrium", 4, single=TRUE)


graphics.off()


break



MakeN_EX_R <- function(mirna, n_ex) {
  if (n_ex == 0) {
    n_ex = ""
  } else {
    n_ex = sprintf("_ex%s", n_ex)
  }
  I4 <- fread(sprintf("/lab/solexa_bartel/mcgeary/AgoRBNS/%s/equilibrium/kmers/I_-3_k8%s.txt", mirna, n_ex))[[2]]
  A4 <- fread(sprintf("/lab/solexa_bartel/mcgeary/AgoRBNS/%s/equilibrium/kmers/4_-3_k8%s.txt", mirna, n_ex))[[2]]
  # cbind(I4, A4)
  cbind(I4, A4, Norm(A4) / Norm(I4))
}


R_l7_0 <- MakeN_EX_R("miR-1", 0)
R_l7_1 <- MakeN_EX_R("miR-1", 1)
R_l7_2 <- MakeN_EX_R("miR-1", 2)
R_l7_3 <- MakeN_EX_R("miR-1", 3)
R_l7_4 <- MakeN_EX_R("miR-1", 4)
R_l7_5 <- MakeN_EX_R("miR-1", 5)

break

A4_test_m1_nex <- fread("/lab/solexa_bartel/mcgeary/AgoRBNS/miR-1/equilibrium/kmers/4_-3_k8_ex6.txt")


# tick <- 0
# sXc <- SitesXCounts("miR-7-23nt", "equilibrium2_nb", sitelist="12mers",
#                     mirna.start=2)
# pars <- InitializeEquilSitePars(sXc)
# pars <- GetSiteKds("miR-7-23nt", "equilibrium2_nb", sitelist="12mers", mirna.start=2)
# model <- ModelC(pars, sXc)
# data <- GetDataEquil(sXc)
# plot(GetDataEquil(sXc), model, log='xy')


# # Compare12merKdsmiR7(type="mean")
# sXc12 <- SitesXCounts("miR-7-23nt", "equilibrium2_nb", 5, "paper")

# kds_23 <- EquilPars("miR-7-23nt", "equilibrium2_nb", 5, "paper")
# kds_24 <- EquilPars("miR-7-24nt", "equilibrium2_nb", 5, "paper")
# kds_25 <- EquilPars("miR-7-25nt", "equilibrium2_nb", 5, "paper")
# kds_all <- EquilPars("miR-7-23nt", "equilibrium2_nb", 5, "paper", global=TRUE)



# state <- c(1, 0, 0, 0, 10)
# names(state) <- c("a_i", "a_ii", "c_i", "c_ii", "l")

# TestModel <- function(t, state, pars) {
#   with(as.list(c(state, pars)), {
#     #dif     bind-birth       conf-birth      bind-birth  conf-death    
#     da_i  <-   k1_r_i *c_i  + alpha_r*a_ii - (l*k1_f_i  + alpha_f)*a_i
#     da_ii <-   k1_r_ii*c_ii + alpha_f*a_i  - (l*k1_f_ii + alpha_r)*a_ii
#     dc_i  <- l*k1_f_i *a_i  + gamma_r*c_ii - (  k1_r_i  + gamma_f)*c_i
#     dc_ii <- l*k1_f_ii*a_ii + gamma_f*c_i  - (  k1_r_ii + gamma_r)*c_ii
#     dl    <- da_i + da_ii 
#     return(list(c(da_i, da_ii, dc_i, dc_ii, dl)))
#   })
# }

# RootFun <- function(t, state, pars) {
#   dstate <- unlist(TestModel(t, state, pars))
#   return(sum(abs(dstate)))
# }

# GenerateParams <- function() {
#   pars <- exp(rnorm(8))
#   names(pars) <- c("k1_f_i", "k1_r_i", "k1_f_ii", "k1_r_ii",
#                    "alpha_f", "alpha_r", "gamma_f", "gamma_r")
#   pars
# }


# SteadyState <- function(pars, state) {
#   tout <- c(0, 1e10)
#   lsodar(state, tout, TestModel, pars, rootfun=RootFun)[2,-1]
# }



# Check_C_II <- function(pars, state) {
#   ode_solution <- SubfunctionCall(SteadyState)
#   with(as.list(c(pars, ode_solution)), {
#     C1 <- 1/(k1_r_i  + gamma_f)
#     C2 <- 1/(k1_r_ii + gamma_r)
#     C3 <- 1/(1 - C1*C2*gamma_f*gamma_r)
#     C4 <- k1_r_ii*C3*C2
#     C5 <- k1_f_i*C4*C1*gamma_f
#     C6 <- k1_f_ii*(1 - C4)
#     a_i_MODEL <- a_i
#     a_ii_MODEL <- a_i_MODEL*(l*C5 + alpha_f)/(l*C6 + alpha_r)
#     c_i_MODEL <-  l*C3*(k1_f_i*C1*a_i_MODEL + gamma_r*k1_f_ii*C1*C2*a_ii_MODEL)
#     c_ii_MODEL <- l*C3*(k1_f_ii*C2*a_ii_MODEL + gamma_f*k1_f_i*C1*C2*a_i_MODEL)
#     print("a_ii:")
#     print(a_ii)
#     print(a_ii_MODEL)
#     print("c_i:")
#     print(c_i)
#     print(c_i_MODEL)
#     print("c_ii:")
#     print(c_ii)
#     print(c_ii_MODEL)
#   })
# }

# pars <- GenerateParams()

# Check_C_II(pars, state)

# break
AggregateAllPaperColors <- function() {
  mirnas_all <- c(kMirnas, "miR-7-23nt")
  exps_all <- c(rep("equilibrium", 5), "equilibrium2_nb")
  mirnas_exps <- cbind(mirnas_all, exps_all)
  print(mirnas_exps)
  sites <- unique(unlist(apply(mirnas_exps, 1, function(row) {
    sites_i <- rownames(SitesXCounts(row[1], experiment=row[2]))
    return(sites_i[-length(sites_i)])
    })))
  print(sites)
  colors <- kSiteColors[sites]
  return(cbind(sites, colors))
}

sites_colors <- AggregateAllPaperColors()
print(sites_colors)
print(sites_colors[is.na(sites_colors[, 2]),])

break



GetVectorPosition <- function(freq_vec, simp="diag") {
  names(freq_vec) <- c("A", "C", "G", "U")
  freq_vec <- Norm(freq_vec)
  matrix_names <- list(names(freq_vec), c("x", "y"))
  PurPymDiagVec <-   matrix(c(-1,  1,
                              -1, -1,
                               1, -1,
                               1,  1), nrow=4, ncol=2, byrow=TRUE,
                            dimnames=matrix_names)
  PurPymEdgeVec <-   matrix(c(-1, -1,
                              -1,  1,
                               1, -1,
                               1,  1), nrow=4, ncol=2, byrow=TRUE,
                            dimnames=matrix_names)
  F <- freq_vec/max(freq_vec)
  mod_F <- sqrt(sum(F^2))
  if (simp == "diag") {
    Mat <- PurPymDiagVec
  } else {
    Mat <- PurPymEdgeVec
  }
  u_v <- c(F %*% Mat)/(sqrt(2)*mod_F)
  names(u_v) <- c("u", "v")
  return(u_v)
}

ConvertToSquare <- function(u_v) {
  names(u_v) <- c("u", "v")
  x_y <- 1/2*c(sqrt(2 + u_v["u"]^2 - u_v["v"]^2 + 2*sqrt(2)*u_v["u"]) -
               sqrt(2 + u_v["u"]^2 - u_v["v"]^2 - 2*sqrt(2)*u_v["u"]),
               sqrt(2 - u_v["u"]^2 + u_v["v"]^2 + 2*sqrt(2)*u_v["v"]) -
               sqrt(2 - u_v["u"]^2 + u_v["v"]^2 - 2*sqrt(2)*u_v["v"]))
  names(x_y) <- c("x", "y")
  return(x_y)
}

xmin <- -1
ymin <- xmin
xmax <- 1
ymax <- xmax

dev.set(2)
par(kPlotParameters)
BlankPlot()
AddLinearAxis(1, tick.space=0.2, label.space=0.2, label="x")
AddLinearAxis(2, tick.space=0.2, label.space=0.2, label="y")

dev.set(3)
par(kPlotParameters)
BlankPlot()
AddLinearAxis(1, tick.space=0.2, label.space=0.2, label="x")
AddLinearAxis(2, tick.space=0.2, label.space=0.2, label="y")


freq_mats <- sapply(1:1000, function(ind) {
  runif(4, 0, 1)
})
u_v_mat <- apply(freq_mats, 2, function(row) {
  out <- GetVectorPosition(row)
  out
})
x_y_mat <- apply(u_v_mat, 2, ConvertToSquare)
dev.set(2)
cols = rainbow(ncol(x_y_mat))
plot(u_v_mat[1,], u_v_mat[2,], col=cols)
dev.set(3)
plot(x_y_mat[1,], x_y_mat[2,], col=cols)
break


# # kmers <- KmersXCounts("miR-7-23nt", "equilibrium2_nb", 40, 5, 8, 0)
# kXc_23 <- KmersXCounts("miR-7-23nt", "equilibrium2_nb", 5, 10, 2)
# kXc_24 <- KmersXCounts("miR-7-24nt", "equilibrium2_nb", 5, 10, 2)
# kXc_25 <- KmersXCounts("miR-7-25nt", "equilibrium2_nb", 5, 10, 2)

graphics.off()
# dev.new(xpos=20, ypos=20, height=12, width=12)
# par(mfrow=c(ncol(kXc), ncol(kXc)))
# par(kPlotParameters)
# xmin <- 10
# xmax <- 1e5
# ymin <- xmin
# ymax <- xmax
# sapply(1:(ncol(kXc)), function(col_ind) {
#   sapply(1:(ncol(kXc)), function(col_ind2) {
#     BlankPlot(log='xy')
#     AddLogAxis(1, label=colnames(kXc)[col_ind])
#     AddLogAxis(2, label=colnames(kXc)[col_ind2])
#     points(kXc[, col_ind], kXc[, col_ind2], pch=20,
#            col=ConvertRColortoRGB("black", alpha=0.2))
#   })
# })

kXc_list <- list(l23=kXc_23, l24=kXc_24, l25=kXc_25)
Compareinput <- function(l1, l2, col_1, col_2, identify=FALSE, col_norm=1) {
  x <- Norm(kXc_list[[sprintf("l%s", l1)]][, col_1]+1)/Norm(kXc_list[[sprintf("l%s", l1)]][, col_norm]+1)
  y <- Norm(kXc_list[[sprintf("l%s", l2)]][, col_2]+1)/Norm(kXc_list[[sprintf("l%s", l2)]][, col_norm]+1)
  plot(x, y, log='xy', xlim=c(1e-2, 1e3), ylim=c(1e-2, 1e3),
       xlab=colnames(kXc_list[[sprintf("l%s", l1)]])[col_1],
       ylab=colnames(kXc_list[[sprintf("l%s", l2)]])[col_2])
  abline(0, 1, lty=2)
  if (identify) {
    identify(x, y, labels=rownames(kXc_23))  
  }
}

Compareinput(24, 25, 2, 2, col_norm=1, identify=TRUE)

break

R_check_0 <- (Norm(kXc[, 2] + 1)/Norm(kXc[, 7] + 1))
R_check_I <- (Norm(kXc[, 2] + 1)/Norm(kXc[, 1] + 1))


names(R_check_0) <- rownames(kXc)
names(R_check_I) <- rownames(kXc)
print(t(head(sort(R_check_I, decreasing=TRUE))))
print(t(head(sort(R_check_0, decreasing=TRUE))))

break

CheckGlobalFit16mers <- function(mirna, n_constant=5) {
  if (mirna == "miR-1") {
    data_i <- 1
    str.nocombI <- "nocombInput_"
  } else {
    data_i <- 2
    str.nocombI <- ""
  }
  if (mirna %in% c("miR-7-23nt", "miR-7-24nt", "miR-7-25t")) {
    experiment <- "equilibrium2_nb"
  } else {
    experiment <- "equilibrium"
  }
  pars_temp_path <- sprintf("16merfits/%s_%s_%s_%s/temp_global.txt", mirna,
                                experiment, n_constant, str.nocombI)
  ind_splits <- c(0, 1)
  names(ind_splits) <- c("left", "right")
  dev.new(xpos=20,ypos=20, height=4, width=10)
  par(mfcol=c(2, 5))
  sapply(1:5, function(mirna.start) {
    sapply(c("left", "right"), function(split16) {
      sXc <- SubfunctionCall(SitesXCounts, sitelist="16mers")
      # Define the number of site types (always the same):
      n_x <- 4^8 + 1
      # This accounts for the fact tha the order is 1_left, 1_right, 2_left, etc.
      kd_start_ind <- (2*(mirna.start - 1) + ind_splits[split16])*n_x
      kds_vec <- exp(data.frame(fread(pars_temp_path, skip=kd_start_ind, nrow=n_x,
                              header=FALSE), row.names=1))
      bg <- data.frame(fread(pars_temp_path, skip=10*n_x, nrow=1,
                              header=FALSE), row.names=1)
      A <- data.frame(fread(pars_temp_path, skip=10*n_x + 1, nrow=1,
                              header=FALSE), row.names=1)

      # Check that the names agree:
      l <- Norm(sXc[, data_i])*100
      bg <- as.numeric(exp(bg[2]))
      A  <- as.numeric(exp(A[2]))
      dils <- as.numeric(colnames(sXc)[3:(ncol(sXc) - 1)])/100
      As <- A*dils
      kds <- unlist(kds_vec)
      names(kds) <- rownames(kds_vec)
      as <- sapply(As, FreeAgo, kds=kds, l=l)
      L <- 100
      totals_all <- colSums(sXc)[3:(ncol(sXc) - 1)]
      n_col <- length(dils)
      model <- matrix(.C("ModelEquil", as.double(kds), as.double(l), as.double(L),
                 as.double(bg), as.double(dils), as.double(As), as.double(as),
                 as.double(totals_all), as.integer(n_x), as.integer(n_col),
                 model=double(n_x*n_col))[["model"]],
            ncol=length(dils), nrow=n_x, byrow=FALSE,
            dimnames=list(rownames(sXc), colnames(sXc)[3:(ncol(sXc) - 1)]))
      data <- sXc[, 3:(ncol(sXc) - 1)]
      par(kPlotParameters)
      xmin <- 0.1
      xmax <- 1e7
      ymin <- xmin
      ymax <- xmax
      BlankPlot(log='xy')
      AddLogAxis(1, label="model")
      AddLogAxis(2, label="data")
      xy <- GetPlotFractionalCoords(fx=0.15, fy=0.95, log='xy')
      text(xy[1], xy[2], sprintf("nt %s-%s", mirna.start,
                                 as.numeric(mirna.start) + 3), adj=0)
      xy <- GetPlotFractionalCoords(fx=0.15, fy=0.85, log='xy')
      text(xy[1], xy[2], split16, adj=0)
      xy <- GetPlotFractionalCoords(fx=0.9, fy=0.1, log='xy')
      text(xy[1], xy[2], mirna, adj=1)

      points(c(model), unlist(data), col=ConvertRColortoRGB("black", alpha=0.2))
    })
  })
}

CheckGlobalFit16mersMiR7 <- function(n_constant=5) {
  experiment <- "equilibrium2_nb"
  data_i <- 2
  pars_temp_path <- sprintf("16merfits/miR-7_equilibrium2_nb_%s_/temp_global.txt", n_constant)
  ind_splits <- c(0, 1)
  names(ind_splits) <- c("left", "right")
  ind_mirs <- 0:2
  names(ind_mirs) <- c("miR-7-23nt", "miR-7-24nt", "miR-7-25nt")

  dev.new(xpos=20,ypos=20, height=12, width=10)
  par(mfcol=c(6, 5))
  sapply(1:5, function(mirna.start) {
    sapply(c("left", "right"), function(split16) {
      # Define the number of site types (always the same):
      n_x <- 4^8 + 1
      # This accounts for the fact tha the order is 1_left, 1_right, 2_left, etc.
      kd_start_ind <- (2*(mirna.start - 1) + ind_splits[split16])*n_x
      kds_vec <- NaN
      while(class(kds_vec) != "data.frame") {
      kds_vec <- try(exp(data.frame(fread(pars_temp_path, skip=kd_start_ind, nrow=n_x,
                              header=FALSE), row.names=1)))
      }
      sapply(c("miR-7-23nt", "miR-7-24nt", "miR-7-25nt"), function(mirna) {
        sXc <- SubfunctionCall(SitesXCounts, sitelist="16mers")
        bg <- NaN
        A <- NaN
        while(class(bg) != "data.frame") {
          bg <- try(data.frame(fread(pars_temp_path, skip=10*n_x + ind_mirs[mirna], nrow=1,
                                  header=FALSE), row.names=1))
        }
        while(class(A) != "data.frame") {
          A <- try(data.frame(fread(pars_temp_path, skip=10*n_x + 3 + ind_mirs[mirna], nrow=1,
                                  header=FALSE), row.names=1))
        }
        # Check that the names agree:
        l <- Norm(sXc[, data_i])*100
        bg <- as.numeric(exp(bg[2]))
        A  <- as.numeric(exp(A[2]))
        dils <- as.numeric(colnames(sXc)[3:(ncol(sXc) - 1)])/100
        As <- A*dils
        kds <- unlist(kds_vec)
        names(kds) <- rownames(kds_vec)
        as <- sapply(As, FreeAgo, kds=kds, l=l)
        L <- 100
        totals_all <- colSums(sXc)[3:(ncol(sXc) - 1)]
        n_col <- length(dils)
        model <- matrix(.C("ModelEquil", as.double(kds), as.double(l), as.double(L),
                   as.double(bg), as.double(dils), as.double(As), as.double(as),
                   as.double(totals_all), as.integer(n_x), as.integer(n_col),
                   model=double(n_x*n_col))[["model"]],
              ncol=length(dils), nrow=n_x, byrow=FALSE,
              dimnames=list(rownames(sXc), colnames(sXc)[3:(ncol(sXc) - 1)]))
        data <- sXc[, 3:(ncol(sXc) - 1)]
        row_inds <- c(sample(nrow(model)-1, 1000, replace=TRUE), nrow(model))
        par(kPlotParameters)
        xmin <- 0.1
        xmax <- 1e8
        ymin <- xmin
        ymax <- xmax
        BlankPlot(log='xy')
        AddLogAxis(1, label="model")
        AddLogAxis(2, label="data")
        xy <- GetPlotFractionalCoords(fx=0.15, fy=0.95, log='xy')
        text(xy[1], xy[2], sprintf("nt %s-%s", mirna.start,
                                   as.numeric(mirna.start) + 3), adj=0)
        xy <- GetPlotFractionalCoords(fx=0.15, fy=0.85, log='xy')
        text(xy[1], xy[2], split16, adj=0)
        xy <- GetPlotFractionalCoords(fx=0.9, fy=0.1, log='xy')
        text(xy[1], xy[2], mirna, adj=1)

        points(c(model[row_inds,]), unlist(data[row_inds,]), col=ConvertRColortoRGB("black", alpha=0.2))
      })
    })
  })
}


CompareGlobal16merKds <- function(mirna, n_constant=5, xpos, ypos) {
  if (mirna == "miR-1") {
    data_i <- 1
    str.nocombI <- "nocombInput_"
  } else {
    data_i <- 2
    str.nocombI <- ""
  }
  if (mirna %in% c("miR-7", "miR-7-23nt", "miR-7-24nt", "miR-7-25nt")) {
    experiment <- "equilibrium2_nb"
  } else {
    experiment <- "equilibrium"
  }
  pars_temp_path <- sprintf("16merfits/%s_%s_%s_%s/temp_global.txt", mirna,
                                experiment, n_constant, str.nocombI)
  print(pars_temp_path)
  ind_splits <- c(0, 1)
  names(ind_splits) <- c("left", "right")
  dev.new(xpos=xpos,ypos=ypos, height=3.5, width=8*3.5/4)
  par(mfrow=c(2, 4))
  n_x <- 4^8 + 1
  kd_start_ind <- 0
  kds_vec_old <- NaN
  while(class(kds_vec_old) != "data.frame") {
  kds_vec_old <- try(exp(data.frame(fread(pars_temp_path, skip=kd_start_ind, nrow=n_x,
                          header=FALSE), row.names=1)))
  }
  xmin <- 1e-5
  xmax <- 1e3
  ymin <- xmin
  ymax <- xmax
  sapply(1:4, function(mirna.ind) {
    kd_start_ind <- 2*mirna.ind*n_x
    kds_vec_new <- NaN
    while(class(kds_vec_new) != "data.frame") {
      kds_vec_new <- try(exp(data.frame(fread(pars_temp_path, skip=kd_start_ind,
                                          nrow=n_x, header=FALSE), row.names=1)))
    }
    overlap <- intersect(rownames(kds_vec_old), rownames(kds_vec_new))
    par(kPlotParameters)
    BlankPlot(log='xy', inv='xy')
    AddLogAxis(1, label=sprintf("nt %s-%s left", mirna.ind, mirna.ind + 3))
    AddLogAxis(2, label=sprintf("nt %s-%s left", mirna.ind + 1, mirna.ind + 4))
    abline(0, 1, lty=2)
    points(kds_vec_old[overlap, 1], kds_vec_new[overlap, 1],
           col=ConvertRColortoRGB("black", alpha=0.2), pch=20)
    kds_vec_old <<- kds_vec_new
    print(length(overlap))
  })
  kd_start_ind <- 1
  kds_vec_old <- exp(data.frame(fread(pars_temp_path, skip=kd_start_ind, nrow=n_x,
                            header=FALSE), row.names=1))
  sapply(1:4, function(mirna.ind) {
    kd_start_ind <- 2*mirna.ind*n_x
    kds_vec_new <- exp(data.frame(fread(pars_temp_path, skip=kd_start_ind,
                                        nrow=n_x, header=FALSE), row.names=1))
    overlap <- intersect(rownames(kds_vec_old), rownames(kds_vec_new))
    par(kPlotParameters)
    BlankPlot(log='xy', inv='xy')
    AddLogAxis(1, label=sprintf("nt %s-%s right", mirna.ind, mirna.ind + 3))
    AddLogAxis(2, label=sprintf("nt %s-%s right", mirna.ind + 1, mirna.ind + 4))
    abline(0, 1, lty=2)
    points(kds_vec_old[overlap, 1], kds_vec_new[overlap, 1],
           col=ConvertRColortoRGB("black", alpha=0.2), pch=20)
    kds_vec_old <<- kds_vec_new
    print(length(overlap))
  })
}



graphics.off()


break


break


temp1 <- read.table("12merfits/miR-7-23nt_equilibrium2_nb_5_12mers_2__nocombInputfull_temp.txt")
temp2 <- read.table("12merfits/miR-7-23nt_equilibrium2_nb_5_12mers_3__nocombInputfull_temp.txt")
temp3 <- read.table("12merfits/miR-7-23nt_equilibrium2_nb_5_12mers_4__nocombInputfull_temp.txt")
temp4 <- read.table("12merfits/miR-7-23nt_equilibrium2_nb_5_12mers_5__nocombInputfull_temp.txt")

temp5 <- read.table("12merfits/miR-7-24nt_equilibrium2_nb_5_12mers_2__nocombInputfull_temp.txt")
temp6 <- read.table("12merfits/miR-7-24nt_equilibrium2_nb_5_12mers_3__nocombInputfull_temp.txt")
temp7 <- read.table("12merfits/miR-7-24nt_equilibrium2_nb_5_12mers_4__nocombInputfull_temp.txt")
temp8 <- read.table("12merfits/miR-7-24nt_equilibrium2_nb_5_12mers_5__nocombInputfull_temp.txt")

temp9  <- read.table("12merfits/miR-7-25nt_equilibrium2_nb_5_12mers_2__nocombInputfull_temp.txt")
temp10 <- read.table("12merfits/miR-7-25nt_equilibrium2_nb_5_12mers_3__nocombInputfull_temp.txt")
temp11 <- read.table("12merfits/miR-7-25nt_equilibrium2_nb_5_12mers_4__nocombInputfull_temp.txt")
temp12 <- read.table("12merfits/miR-7-25nt_equilibrium2_nb_5_12mers_5__nocombInputfull_temp.txt")


rows_i <- intersect(rownames(temp1), rownames(temp2))
xmin <- 1e-6
xmax <- 1e3
ymin <- xmin
ymax <- xmax
dev.new(xpos=20, ypos=20, height=10, width=10)
par(kPlotParameters)
color <- ConvertRColortoRGB("black", alpha=0.3)
par(mfrow=c(5, 3))
BlankPlot(log='xy')
AddLogAxis(1, label="2-5")
AddLogAxis(2, label="3-6")
points(temp1[rows_i, 1], temp2[rows_i, 1], col=color)
abline(0, 1)
rows_i <- intersect(rownames(temp2), rownames(temp3))
BlankPlot(log='xy')
AddLogAxis(1, label="2-5")
AddLogAxis(2, label="3-6")
points(temp2[rows_i, 1], temp3[rows_i, 1], col=color)
abline(0, 1)

rows_i <- intersect(rownames(temp3), rownames(temp4))
BlankPlot(log='xy')
AddLogAxis(1, label="2-5")
AddLogAxis(2, label="3-6")
points(temp3[rows_i, 1], temp4[rows_i, 1], col=color)
abline(0, 1)


rows_i <- intersect(rownames(temp5), rownames(temp6))
BlankPlot(log='xy')
AddLogAxis(1, label="2-5")
AddLogAxis(2, label="3-6")
points(temp5[rows_i, 1], temp6[rows_i, 1], col=color)
abline(0, 1)

rows_i <- intersect(rownames(temp6), rownames(temp7))
BlankPlot(log='xy')
AddLogAxis(1, label="2-5")
AddLogAxis(2, label="3-6")
points(temp6[rows_i, 1], temp7[rows_i, 1], col=color)
abline(0, 1)

rows_i <- intersect(rownames(temp7), rownames(temp8))
BlankPlot(log='xy')
AddLogAxis(1, label="2-5")
AddLogAxis(2, label="3-6")
points(temp7[rows_i, 1], temp8[rows_i, 1], col=color)
abline(0, 1)



rows_i <- intersect(rownames(temp9), rownames(temp10))
BlankPlot(log='xy')
AddLogAxis(1, label="2-5")
AddLogAxis(2, label="3-6")
points(temp9[rows_i, 1], temp10[rows_i, 1], col=color)
abline(0, 1)

rows_i <- intersect(rownames(temp10), rownames(temp11))
BlankPlot(log='xy')
AddLogAxis(1, label="2-5")
AddLogAxis(2, label="3-6")
points(temp10[rows_i, 1], temp11[rows_i, 1], col=color)
abline(0, 1)

rows_i <- intersect(rownames(temp11), rownames(temp12))
BlankPlot(log='xy')
AddLogAxis(1, label="2-5")
AddLogAxis(2, label="3-6")
points(temp11[rows_i, 1], temp12[rows_i, 1], col=color)
abline(0, 1)

BlankPlot(log='xy')
AddLogAxis(1, label="2-5")
AddLogAxis(2, label="3-6")
points(temp1[, 1], temp5[, 1], col=color)
abline(0, 1)

BlankPlot(log='xy')
AddLogAxis(1, label="2-5")
AddLogAxis(2, label="3-6")
points(temp5[, 1], temp9[, 1], col=color)
abline(0, 1)

BlankPlot(log='xy')
AddLogAxis(1, label="2-5")
AddLogAxis(2, label="3-6")
points(temp2[, 1], temp6[, 1], col=color)
abline(0, 1)

BlankPlot(log='xy')
AddLogAxis(1, label="2-5")
AddLogAxis(2, label="3-6")
points(temp6[, 1], temp10[, 1], col=color)
abline(0, 1)



break





InitializeEquilSitePars <- function(sXc, combined=TRUE) {
  l <- SubfunctionCall(GetInputEquil)
  # print(l)
  data <- SubfunctionCall(GetDataEquil)
  # print(data)
  kds <- log10(Norm(l)/Norm(rowSums(data)))/2
  # kds <- kds*0
  kds <- kds - kds[length(kds)]
  names(kds) <- paste0(rownames(sXc), "_Kd")
  pars <- c(kds, bg=-4, AGO=1)
  random_pars <- rnorm(length(pars), 0, 3)
  names(random_pars) <- names(pars)
  random_pars["None_Kd"] <- 0
  pars <- random_pars
  return(pars)
}




sXc <- SitesXCounts("miR-1", "equilibrium")
n_x <- nrow(sXc)

dil <- as.numeric(colnames(sXc)[3:7])
pars <- InitializeEquilSitePars(sXc)
ln_pars <- pars*log(10)	

time1 <- proc.time()[3]
grad_analytical1 <- GradientC(pars, sXc)
time2 <- proc.time()[3]
grad_analytical2 <- GradCNew(ln_pars, sXc, dil, n_x)
time3 <- proc.time()[3]
print(time2 - time1)
print(time3 - time2)
plot(grad_analytical1, grad_analytical2)
abline(0,log(10), col="red")
abline(0, 1/(log(10)), col="blue")

break
grad_numerical <- grad(func=CostCNew, x=ln_pars, sXc=sXc)
grad_analytical <- GradCNew(ln_pars, sXc)
grad_numerical[nrow(sXc)] <- 0
plot(grad_numerical, grad_analytical, xlim = c(-1e7, 1e7), ylim=c(-1e7, 1e7))

sapply(1:10, function(i) {
	pars <- InitializeEquilSitePars(sXc)
	ln_pars <- pars*log(10)	
	grad_numerical <- grad(CostCNew, ln_pars, sXc=sXc)
	grad_analytical <- GradCNew(ln_pars, sXc)
	grad_numerical[nrow(sXc)] <- 0
	points(grad_numerical, grad_analytical,
	       col=c(kSiteColors[rownames(sXc)], "gray", "brown", "black"))
})
abline(0, 1)
break
# cost_old <- CostNew(pars, sXc)
# cost_C <- CostC(pars, sXc)


# cost_C_new <- CostCNew(ln_pars, sXc)
# print(cost_old)
# print(cost_C)
# print(cost_C_new)
grad_numerical <- grad(CostCNew, ln_pars, sXc=sXc)
grad_analytical <- GradCNew(ln_pars, sXc)
grad_analytical[nrow(sXc)] <- 0
plot(grad_numerical, grad_analytical)
abline(0, 1)
break
kds <- 10^pars[1:nrow(sXc)]
A_old <- 10^pars[nrow(sXc) + 2]
l <- GetInputEquil(sXc)
dils <- sapply(colnames(sXc[3:7]), as.numeric)
a_frees <- sapply(A_old*dils/100, FreeAgo, l=l, kds=kds)

model_old <- ModelC(pars, sXc)
model_new <- matrix(cost_newC[[5]], nrow=nrow(sXc), ncol=5,
                    dimnames=list(rownames(sXc), colnames(sXc[, 3:7])), byrow=FALSE)
break


dev.new(xpos=20, ypos=20)
plot(GradientNew(pars, sXc=sXc), GradientC(pars, sXc), xlim = c(-1e7, 1e7), ylim=c(-1e7, 1e7))

sapply(1:10, function(i) {
	pars <- InitializeEquilSitePars(sXc)
	points(GradientNew(pars, sXc=sXc), GradientC(pars, sXc), col=c(kSiteColors[rownames(sXc)], "gray", "brown", "black"))

	})

abline(0, 1)