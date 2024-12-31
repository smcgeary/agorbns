################################################################################
# ThrPFunctions.R
################################################################################
# This script adds functions used for in McGeary, Bisaria, et al. 2022, which
# focuses on analysis of pairing of the miRNA 3-prime end, using AGO-RBNS
# libraries that include a programmed site against a particular miRNA.
#
# It is required for generating the figures in `general/Figures2022.R`, but not
# for generating the figures in `general/Figures2019.R`.
#
# Analyses using the processed data must be performed in order to generat each
# figure, which can be done by copying the commands before each figure subpanel
# into a terminal, assuming the `agorbns` conda environment has been activated.
source("general/general.R")

library(reticulate)
if (!py_available()) {
  use_condaenv(Sys.getenv("CONDA_DEFAULT_ENV"), required = TRUE)
}


# Figure 1DC.
GetPositionalKmersProgrammedLib <- function(
  mirna, experiment, condition, n_constant, kmer_len, seedex=FALSE
)
{
  ext <- sprintf("_%s_k%s", n_constant, kmer_len)
  if (seedex) {
    ext <- sprintf("%s_seedex", ext)
  }
  path <- SubfunctionCall(GetAnalysisPath,
                          analysis="kmers_positional_programmed")
  print(path)
  system(sprintf("ls -l %s", path))
  out <- read.table(path, sep="\t", header=TRUE, row.names=1)
  return(out)
}


LargestCommonSubstring <- function(seq1, seq2) {
  # Get the length of each sequence.
  len_1 <- nchar(seq1)
  len_2 <- nchar(seq2)
  # Get the vector of prefix strings for each of the two sequences.
  # I.E. c("h", "he", "hel", "hell", "hello") for "hello"
  pres1 <- sapply(1:len_1, substr, x=seq1, start=1)
  pres2 <- sapply(1:len_2, substr, x=seq2, start=1)
  # Make a matrix in which each row is a pair of two prefixes, one from each
  # sequence.
  prefix_pairs <- cbind(rep(pres1, each=length(pres2)),
                        rep(pres2, length(pres1)))
  # Find the largest common suffix for each pair of prefixes.
  common_suffixes <- apply(prefix_pairs, 1, function(row) {
    LargestCommonSuffix(row[1], row[2])
  })
  # Return the suffix or suffixes that have the greatest length.
  return(
   unique(
     common_suffixes[
        which(
          nchar(common_suffixes) == max(nchar(common_suffixes))
        )
      ]
    )
  )
}


LargestCommonSuffix <- function(seq1, seq2) {
  # Get the length of each sequence.
  len_1 <- nchar(seq1)
  len_2 <- nchar(seq2)
  # If both sequences are empty strings return an empty string.
  if (len_1 == 0 || len_2 == 0) {
    return("")
  # Otherwise, check if the last character of each sequence is the same.
  # If so, return a string that is the result of this function recursively
  # merged with the the last (shared) character.
  } else if (substr(seq1, len_1, len_1) == substr(seq2, len_2, len_2)) {
    return(
      paste0(
        LargestCommonSuffix(
          substr(seq1, 1, len_1 - 1), substr(seq2, 1, len_2 - 1)
        ),
        substr(seq1, len_1, len_1)
      )
    )
  # Otherwise, return an empty string (there is no common suffix).
  } else {
    return("")
  }
}


GeoMean <- function(vector) {
  exp(mean(log(vector), na.rm=TRUE))
}


ApplyKdCorrection <- function(
  mirna, experiment, prog_n_constant=3, rand_n_constant=3,
  prog_sitelist="progthrp_suppcomp", rand_sitelist="randthrp_comp",
  start_mm=FALSE, stop_mm=FALSE, new=FALSE, new2=FALSE, sumseed=FALSE,
  prog_mm=TRUE, combined_rand=TRUE, buffer_rand=FALSE, span=10,
  return_swapped=TRUE, print_optim=FALSE, xpos=20, ypos=20
) {
  # The `prog_mm` variable is whether or not the kd for the variable site is
  # used in the programmed region, or if the average of it within the random
  # region.
  # 1.A Define the function parameters to load the appropriate random kds. These
  # are the seed-cognate miRNAs, so it will either be let-7a, miR-1, or miR-155.
  # Because miR-1 has the different library, it both doesn't use combined input,
  # and uses the "buffer" flag, which indicates that the last two nucleotides of
  # the data aren't used, due to the enrichment of the 8mer and 8mer variant
  # site types.
  # print("in ApplyKdCorrection")
  if (mirna %in% c("let-7a-21nt", "let-7a_minus1", "let-7a_plus1",
                   "let-7a_miR-155")) {
    mirna_rand <- "let-7a"
  } else if (mirna %in% c("miR-155", "miR-155_let-7a")) {
    mirna_rand <- "miR-155"
  } else if (mirna == "miR-1") {
    mirna_rand <- "miR-1"
    combined_rand <- FALSE
    buffer_rand <- TRUE
  } else {
    message("There is no random data for this miRNA.")
    return()
  }
  if (grepl("suppcomp", prog_sitelist)) grep_base_string <- "Comp"
  else                                  grep_base_string <- "8mer-mm[ACTG][2-7]"
  # 1.B Get the programmed-site kds that are to be corrected.
  kds_p <- SubfunctionCall(EquilPars, sitelist=prog_sitelist,
                           n_constant=prog_n_constant)
  # 1.C Get the random kds for which the comparison is to be made.
  kds_r <- SubfunctionCall(EquilPars, mirna=mirna_rand,
                           experiment="equilibrium", n_constant=rand_n_constant,
                           sitelist=rand_sitelist, combined=combined_rand,
                           buffer=buffer_rand, nbomitc=FALSE, start_mm=FALSE,
                           stop_mm=FALSE, new=FALSE, new2=FALSE)
  # Replace the measured programmed seed mismatches with those of the average
  # of those found in the random portion of the programmed library.
  kds_p_init <- kds_p
  kds_p <- SubfunctionCall(SwapProgrammedKds, kds=kds_p,
                           suppcomp=grepl("suppcomp", prog_sitelist))
  # Keeping this to nail down later why I did this in the first place.
  # if (mirna == "miR-155") {
  #   print(kds_p[c("8mer_Kd", "9mer-m13.21|14|Comp_Kd"), ])
  # }

  if (sumseed) {
    sites_use <- c(sprintf("%s_Kd", kSeedSites), "Comp_Kd")
  } else {
    sites_use <- c(sprintf("%s_Kd", c(kSeedSites, GetAll8merMmSites(mirna))))
  }
  lkds_p_seed <- log(kds_p[sites_use, 2])
  lkds_r_seed <- log(kds_r[sites_use, 2])
  # Define the log(ratio) of the random-library kds to that of the programmed
  # libraries.
  lkds_rp_ratio <- lkds_r_seed - lkds_p_seed
  # Use LOESS in R to get a locally fit regression line between that ratio and
  # log_kd_values of the programmed libraies, in order to get the value that
  # the programmed library kd values should be multiplied by as a function of
  # the programmed kd values themselves.
  # print("before loess_object")
  # print(head(lkds_r_seed))
  # print(head(lkds_rp_ratio))
  # print(head(lkds_p_seed))
  loess_object <- loess(lkds_rp_ratio ~ lkds_p_seed, span=span,
                          control=loess.control(surface="direct"))
  # Use the initial (not swapped) kd values if the flag says so.
  if (!(return_swapped)) kds_p <- kds_p_init
  for (col_i in 1:ncol(kds_p)) {
    # Get the ratios values from the LOESS fit for each kd value.
    data_temp <- log(kds_p[, col_i])
    inds_na <- which(is.na(data_temp))
    if (length(inds_na) != 0) {
      print(kds_p[inds_na, ])
      break
    }
    lkds_rp_ratio <- predict(
      loess_object, data.frame(lkds_p_seed=c(log(kds_p[, col_i])))
    )
    # Update each programmed library kd value using this ratio.
    kds_p[, col_i] <- kds_p[, col_i]*exp(lkds_rp_ratio)
  }
  # Not sure why I was doing this this way, aparently there was something wrong
  # with the 9mer site that I was checking.
  # if (mirna == "miR-155" | mirna == "let-7a-21nt") {
  #   print(kds_p[c("8mer_Kd", "9mer-m13.21|14|Comp_Kd"), ])
  # }
  return(kds_p)
}


# Function to swap Kds between different datasets, used for LOES  correction of
# programmed Kds.
SwapProgrammedKds <- function(kds, prog_mm=TRUE, suppcomp=FALSE, sumseed=FALSE) {
  # Define the base string for the sites that will be used as replacement of the
  # programmed-region site.
  if (suppcomp) grep_base_string <- "Comp"
  else         grep_base_string <- "8mer-mm[ACTG][2-7]"
  # Boolean for wether or not to use the programmed mismatch sites, or to use
  # the geometric mean of that programmed site at all instances within the
  # random region of the programmed library.
  if (prog_mm) {
    # This will just give the 8mer, 7mer-m8, ... ,6mer-A1 sites.
    seed_sites <- sprintf("%s_Kd", kSeedSites)
  } else {
    # This will give the seed sites plus the 18 mismatch sites. This relies on
    # the idea that if the site is not a seed or seed-mismatch site, it will
    # contain the character "|", or either of "None", "AGO", or "bg" in its
    # name.
    seed_sites <- grep("(\\||None|AGO|bg)", rownames(kds), invert=TRUE,
                       perl=TRUE, value=TRUE)     
  }
  # 2.B This relaces the seed site kds of the programmed libraries with the
  # average over all of the positions of the site in random region of the 
  # programmed library.
  if (sumseed) {
    kd_mm <- kds["Comp_Kd", 2]
  } else {
    kd_mm <- GeoMean(kds[7:24, 2])
  }
  for (site in rev(seed_sites)) {
    # Goes from from "8mer_Kd" to "8mer"
    site_trim <- unlist(strsplit(site, split="_"))[1]
    # Make the regex target, and use grep to get the row indeces to be averaged.
    # This will end up being "8mer|.*|Comp_Kd" or
    # "8mer|.*|8mer-mm[ACGT][2–7]_Kd"
    regex_target <- sprintf("%s\\|.*\\|%s_Kd", site_trim, grep_base_string)
    # Use this target string to pull out the correct rows of the kds.
    rows_ave <- grep(regex_target, rownames(kds), perl=TRUE)
    # print(rownames(kds)[rows_ave])
    # print(kds[rows_ave, ])
    # Get the row index for where the averaged Kd value is to be assigned.
    row_replace <- grep(sprintf("^%s$", site), rownames(kds), perl=TRUE)
    if (length(row_replace) == 0) {
      print("length zero")
    }
    # Get the geometric mean of the 18 mismatch sites within the progrmamed
    # region of the library, in order to correct for the fact that the kd being
    # averaged reflects not just the seed site but that it also has one of the
    # 18 possible seed mismatch sites.
    kd_app <- apply(kds[rows_ave,], 2, GeoMean)
    # Perform the correction, which is the agebraic inverse of 
    # kd_app = kd*kd_mm/(kd + kd_mm)
    # kd_app*(kd + kd_mm) = kd*kd_mm
    # kd*kd_app + kd_app*kd_mm = kd*kd_mm
    # kd*(kd_app - kd_mm) = -kd_app*kd_mm
    # kd = -kd_app*kd_mm/(kd_app - kd_mm)
    # kd = kd_app*kd_mm/(kd_mm - kd_app)
    if (length(row_replace) == 0) {
      # kds <- rbind(kd_app*kd_mm/(kd_mm - kd_app), kds)
      kds <- rbind(kd_app, kds)      
      rownames(kds)[1] <- site
    } else {
      # kds[row_replace, ] <- kd_app*kd_mm/(kd_mm - kd_app)
      kds[row_replace, ] <- kd_app
    }
  }
  return(kds)
}



GetAll8merMmSites <- function(mirna) {
  # First get the 8mer site and find the fist nucleotide.
  site_8mer_vec <- strsplit(StrRev(GetSiteSeq(mirna, "8mer")), split="")[[1]]
  # Make each of the internal mismatch names.
  site_names <- sapply(2:7, function(position) {
    # Identify the nucleotide at the current position.
    pos_nuc <- site_8mer_vec[position]
    # Return the other three nucleotides other than the current nucleotide.
    mmnucs <- sort(setdiff(kDNucs, pos_nuc))
    # Return the string that indicates that position and mismatch nucleotide.
    sprintf("8mer-mm%s%s", mmnucs, position)
  })
  # Output `site_names` as a 1-D vector.
  return(as.vector(site_names))
}


GetSiteSeq <- function(mirna, site) {
  source_python("general/RBNS_methods.py")
  main <- py_run_string(sprintf("sequence = Mirna('%s')['%s']", mirna, site))
  main$sequence
}


# GetSiteSeq <- function(mirna, site) {
#   python.load("general/RBNS_methods.py")
#   python.exec(sprintf("_mirna = Mirna('%s')", mirna))
#   python.exec(sprintf("sequence = _mirna['%s']", site))
#   python.get("sequence")
# }


# # String functions that I used to either make 8-nt versions of the four
# # canonical sites (the 8mer, 7mer-m8, 7mer-A1, and 6mer) or the 18
# # 8mer-mismatch sites.



# # This function makes a matrix that has a tab in the first column and the first
# # row, which R doesn't normally do.
# MakeAltMatrix <- function(df) {
#   matrix_new <- matrix("", nrow=nrow(df) + 1, ncol=ncol(df) + 1)
#   for (i in 1:ncol(df)) {
#     matrix_new[2:(nrow(df) + 1), i + 1] <- df[1:nrow(df), i]
#   }
#   matrix_new[1, -1] <- colnames(df)
#   matrix_new[-1, 1] <- rownames(df)
#   return(matrix_new)
# }


# MakeGEOProcessedData <- function() {
#   mirnas <- c("let-7a-21nt", "let-7a-21nt", "miR-1", "miR-155", "let-7a_plus1",
#               "let-7a_minus1", "let-7a_miR-155", "miR-155_let-7a")
#   experiments <- c("equil_c2_nb", "equil_c_nb", "equil_c_nb", "equil_sc_nb",
#                    "equil_c_nb", "equil_c_nb", "equil_c_nb", "equil_c_nb")

#   mirna_labels <- mirnas
#   mirna_labels[5] <- "let-7a(+1)"
#   mirna_labels[6] <- "let-7a(-1)"
#   mirna_labels[7] <- "let-7a-miR-155"
#   mirna_labels[8] <- "miR-155-let-7a"
#   # mir_exp_df <- cbind(mirnas, experiments, mirna_labels)

#   # print(mir_exp_df)
#   out_path <- "/lab/solexa_bartel/mcgeary/GEO_McGearyBisaria2021/processed_files/"
#   for (i in 7:8) {
#     print(i)
#     label <- sprintf("agorbns_%s", mirna_labels[i])
#     if (mirnas[i] == "let-7a-21nt") {
#       if (experiments[i] == "equil_c2_nb") {
#         label <- sprintf("%s_purification1", label)
#       } else {
#         label <- sprintf("%s_purification2", label)
#       }
#     }
#     # 1. Write the sitecount file. #############################################
#     sXc <- SitesXCounts(mirnas[i], experiments[i], n_constant=3,
#                         sitelist="progthrp")
#     sXc <- sXc[, -2]
#     write.table(MakeAltMatrix(sXc), file=sprintf("%s%s_sitecounts.txt",
#                                                  out_path, label),
#                 sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
#     # 2. Write the summed sitecount file. ######################################
#     sXc_combined <- SitesXCounts(mirnas[i], experiments[i], n_constant=3,
#                                  sitelist="progthrp_suppcomp")
#     sXc_combined <- sXc_combined[, -2]
#     print(head(sXc_combined))
#     write.table(MakeAltMatrix(sXc_combined),
#                 file=sprintf("%s%s_summed_8mer_mismatch_sitecounts.txt",
#                              out_path, label),
#                 sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
#     # 3. Write the kd file. ####################################################
#     kds <- ApplyKdCorrection(mirnas[i], experiments[i], rand_n_constant=3,
#                              prog_n_constant=3, rand_sitelist="randthrp",
#                              prog_sitelist="progthrp")
#     kds <- kds[, c(1, 3, 5)]
#     colnames(kds) <- c("MLE", "LowerCI", "UpperCI")
#     print(head(kds))
#     write.table(MakeAltMatrix(kds), file=sprintf("%s%s_kds.txt",
#                                                  out_path, label),
#                 sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
#     # 4. Write the summed kd file. #############################################
#     kds_combined <- ApplyKdCorrection(mirnas[i], experiments[i],
#                                       rand_n_constant=3, prog_n_constant=3,
#                                       rand_sitelist="randthrp_comp",
#                                       prog_sitelist="progthrp_suppcomp")
#     kds_combined <- kds_combined[, c(1, 3, 5)]
#     colnames(kds_combined) <- c("MLE", "LowerCI", "UpperCI")
#     print(head(kds_combined))
#     write.table(MakeAltMatrix(kds_combined),
#                 file=sprintf("%s%s_summed_8mer_mismatch_kds.txt",
#                              out_path, label),
#                 sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
#     if (label != "agorbns_let-7a-21nt_purification2") {
#       if (label == "agorbns_let-7a-21nt_purification1") {
#         label <- "agorbns_let-7a-21nt"
#       }
#       # 5. Write the pairing and offset model coefficients. ####################
#       model1 <- FitPairingAndOffsetModelWithError(mirnas[i], experiments[i])
#       # Re-format the pairing coefficients to be in a table format
#       p_coefs <- model1$pairing
#       p_matrix <- matrix(c(p_coefs$MLE, p_coefs$LowerCI,
#                                         p_coefs$UpperCI), ncol=3)
#       # Name the rows.
#       rownames(p_matrix) <- sprintf(
#         "%s|%s",
#         rep(colnames(p_coefs$MLE), each=nrow(p_coefs$MLE)),
#         rep(rownames(p_coefs$MLE), times=ncol(p_coefs$MLE))
#       )
#       # name the columns.
#       p_matrix <- p_matrix[!is.na(p_matrix[, 1]), , drop=FALSE]
#       colnames(p_matrix) <- c("MLE", "LowerCI", "UpperCI")
#       model_out <- rbind(p_matrix, model1$offset)
#       print(head(model_out))
#       print(tail(model_out))
#       write.table(MakeAltMatrix(model_out),
#                   file=sprintf("%s%s_pairing_and_offset_model_coefficients.txt",
#                                out_path, label),
#                   sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
#       # 6. Write the pairing, offset, and mismatch model coefficients. #########
#       model2 <- FitPairingOffsetAndMismatchModelWithError(mirnas[i],
#                                                           experiments[i])

#       p_coefs <- model2$pairing
#       p_matrix <- matrix(c(p_coefs$MLE, p_coefs$LowerCI,
#                                         p_coefs$UpperCI), ncol=3)
#       # Name the rows.
#       rownames(p_matrix) <- sprintf(
#         "%s|%s",
#         rep(colnames(p_coefs$MLE), each=nrow(p_coefs$MLE)),
#         rep(rownames(p_coefs$MLE), times=ncol(p_coefs$MLE))
#       )
#       # name the columns.
#       p_matrix <- p_matrix[!is.na(p_matrix[, 1]), , drop=FALSE]
#       colnames(p_matrix) <- c("MLE", "LowerCI", "UpperCI")
#       model_out <- rbind(p_matrix, model2$offset, model2$mm)
#       print(head(model_out))
#       print(tail(model_out))
#       write.table(MakeAltMatrix(model_out),
#                   file=sprintf("%s%s_pairing_offset_and_mismatch_model_coefficients.txt",
#                                out_path, label),
#                   sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
#     }
#   }
# }

# # Used in the course of making the processed mismatch tables from Figure 7 and
# # Figure S12. 
# MakeMismatchTable <- function(
#   mirna, experiment, start_mm, stop_mm, bulge=FALSE, n_constant=3,
#   site_base=FALSE, offset_lim=c(-4, 16), kd_fc=FALSE, win_average=1,
#   best_average=FALSE, new=TRUE, corrected_kds=TRUE, upper_bound=FALSE,
#   lower_bound=FALSE
# ) {
#   kdfc_mmdb <<- log10(SubfunctionCall(GetThreePrimeMmBDKds, best_average=TRUE))
#   if (bulge) {
#     ncol_mat <- stop_mm - start_mm
#     colnames_mat <- as.character((start_mm + 1):stop_mm)
#   } else {
#     ncol_mat <- stop_mm - start_mm + 1
#     colnames_mat <- as.character(start_mm:stop_mm)
#   }
#   if (bulge) {
#     nrow_mat <- 5
#     rownames_mat <- c("A", "C", "G", "T", "d")
#   } else {
#     nrow_mat <- 4
#     rownames_mat <- c("A", "C", "G", "T")
#   }
#   # Pre-allocate the output values of the matrix.
#   logkdfc_mat <- matrix(NaN, nrow=nrow_mat, ncol=ncol_mat,
#                         dimnames=list(rownames_mat, colnames_mat))
#   for (mm_i in names(kdfc_mmdb)[-1]) {
#     nuc <- substr(mm_i, start=1, stop=1)
#     if (!(nuc %in% c("A", "C", "G", "T"))) {
#       nuc <- "d"
#       par_pos <- 1
#     } else {
#       par_pos <- 2
#     }
#     if (substr(mm_i, par_pos, par_pos) == "(") {
#       splits <- unlist(strsplit(unlist(strsplit(mm_i, split="\\("))[2],
#                        split="\\)"))[1]
#       range <- unlist(strsplit(splits, split="\\."))
#       pos_list <- as.character(seq(as.integer(range[1]), as.integer(range[2])))
#       for (pos in pos_list) {
#         logkdfc_mat[nuc, pos] <- kdfc_mmdb[mm_i]
#       }
#     } else {
#       pos <- substr(mm_i, start=par_pos, stop=nchar(mm_i))
#       logkdfc_mat[nuc, pos] <- kdfc_mmdb[mm_i]
#     }
#   }
#   if (!bulge) {
#     for (i in colnames(logkdfc_mat)) {
#       mir_nuc <- RevComplement(substr(kMirnaSeqs[mirna], start=i, stop=i))
#       logkdfc_mat[mir_nuc, i] <- kdfc_mmdb["wt"]
#     }
#   }
#   rownames(logkdfc_mat)[4] <- "U"
#   if (bulge) {
#     rownames(logkdfc_mat)[5] <- "Del."
#     logkdfc_mat <- logkdfc_mat[c("Del.", "A", "C", "G", "U"), ]
#   }
#   return(logkdfc_mat)
# }

# # Makes the mismatch and bulge-and-del kd foldchange processed GEO tables.
# MakeGEOProcessedMismatchData <- function() {
#   mirnas <- c("let-7a-21nt", "miR-1", "miR-155")
#   experiments <- c("equil_c2_nb", "equil_c_nb", "equil_sc_nb")
#   mismatch_starts <- matrix(c(11, 11, 11, 11, 11,  9,  9, 10,
#                               12, 12, 12, 12, 11, 11, 11, 11,
#                               13, 13, 13, 15, 15, 13, 13, 13),
#                             ncol=3, byrow=FALSE)

#   mismatch_stops <- matrix(c(14, 15, 16, 17, 18, 17, 18, 20,
#                              15, 16, 17, 18, 18, 19, 20, 21,
#                              16, 17, 18, 21, 22, 21, 22, 23),
#                            ncol=3, byrow=FALSE)
#   out_path <- "/lab/solexa_bartel/mcgeary/GEO_McGearyBisaria2021/processed_files/"
#   for (i in 1:3) {
#     mirna <- mirnas[i]
#     experiment <- experiments[i]
#     mismatch_start_i <- mismatch_starts[, i]
#     mismatch_stop_i <- mismatch_stops[, i]
#     label <- sprintf("agorbns_%s", mirna)
#     if (mirnas[i] == "let-7a-21nt") {
#       if (experiment == "equil_c2_nb") {
#         label <- sprintf("%s_purification1", label)
#       } else {
#         label <- sprintf("%s_purification2", label)
#       }
#     }
#     # Iterate over the number of mismatch tables to be made (8 in total)
#     for (j in 1:nrow(mismatch_starts)) {
#       start_ij <- mismatch_starts[j, i]
#       stop_ij <- mismatch_stops[j, i]
#       # Make the mismatch table
#       mm_df <- 10^MakeMismatchTable(
#         mirna, experiment, start_ij, stop_ij, offset_lim=c(-2, 6), kd_fc=TRUE,
#         win_average=3
#       )
#       bd_df <- 10^MakeMismatchTable(
#         mirna, experiment, start_ij, stop_ij, offset_lim=c(-2, 6), kd_fc=TRUE,
#         win_average=3, bulge=TRUE
#       )
#       colnames(mm_df) <- paste0("Pos. ", colnames(mm_df), sep="")
#       colnames(bd_df) <- paste0("Pos. ", colnames(bd_df), sep="")
#       rownames(bd_df) <- c("Del.", "Bulged A", "Bulged C", "Bulged G", "Bulged U")
#       mm_label <- sprintf("%s_positions%s-%s_mismatch_kd_foldchanges.txt",
#                           label, start_ij, stop_ij)
#       bd_label <- sprintf("%s_positions%s-%s_bulge_and_del_kd_foldchanges.txt",
#                           label, start_ij, stop_ij)
#       write.table(MakeAltMatrix(mm_df), file=sprintf("%s%s", out_path, mm_label),
#                   sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
#       print(mm_label)
#       write.table(MakeAltMatrix(bd_df), file=sprintf("%s%s", out_path, bd_label),
#                   sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
#       print(bd_label)
#     }
#   }
# }











# GetAll8mersCanonicalSite <- function(mirna) {
#   # First get the 8mer site and find the fist nucleotide.
#   m8 <- substr(GetSiteSeq(mirna, "8mer"), 1, 1)
#   mm8 <- sort(setdiff(kDNucs, m8))
#   mmA1 <- sort(setdiff(kDNucs, "A"))
#   sites_7merm8 <- sprintf("7mer-m8-%s", mmA1)
#   sites_7merA1 <- sprintf("%s-7mer-A1", mm8)
#   sites_6mer <- sprintf("%s-6mer-%s", rep(mm8, each=3), rep(mmA1, 3))
#   names_all <- c("8mer", sites_7merm8, sites_7merA1, sites_6mer)
#   names_all
# }



# # This function loads the frequencies of the sites contained within the
# # programmed-region of the three-prime libraries.
# GetProgrammedKmerFrequencies <- function(
#   mirna, experiment, condition
# ) {
#   path <- SubfunctionCall(GetAnalysisPath,
#                           analysis="kmers_positional_programmed_only", ext="_0")
#   kmer_table <- read.table(path, sep="\t", row.names=1, header=FALSE)
#   colnames(kmer_table) <- c(sprintf("%s_%s", mirna, experiment))
#   kmer_table <- kmer_table/sum(kmer_table)
#   kmer_table  
# }





# GetRandomPositionOfKds <- function(kd_data) {
#   # Function that pulls out just the positions within the random region of the
#   # library for each of the Kd values.
#   # Get the names associated with the kds.
#   sites <- rownames(kd_data)
#   # Define the relevant regular expression, that splits up the three parts of
#   # the Kd value (e.g., "5mer-m11.15|9|8mer-mmA7_Kd" will generate three capture
#   # groups, "5mer-m11.15", "9", and "8mer-mmA7_Kd").
#   grep_query <- "^(.*)\\|(.*)\\|(.*)$"
#   positions <- unique(gsub(grep_query, sites, replace="\\2", perl=TRUE))
#   return(positions)
# }


# GetSitePositions <- function(sites) {
#   grep_query <- "^(.*)\\|(.*)\\|(.*)$"
#   as.integer(gsub(grep_query, sites, replace="\\2", perl=TRUE))
# }



# GetProgrammedSiteOfKds <- function(kd_data) {
#   # Function that returns the programmed site associated with each of the Kd
#   # values.
#   # Get the names associated with the kds.
#   sites <- rownames(kd_data)
#   # Define the relevant regular expression, that splits up the three parts of
#   # the Kd value (e.g., "5mer-m11.15|9|8mer-mmA7_Kd" will generate three capture
#   # groups, "5mer-m11.15", "9", and "8mer-mmA7_Kd").
#   grep_query <- "^(.*)\\|(.*)\\|(.*)$"
#   programmed_sites <- unique(gsub(grep_query, sites, replace="\\3", perl=TRUE))
#   return(programmed_sites)
# }


# Makes the dataframe that is used in the linear model.
MakeKdDataFrame <- function(
  mirna, experiment, sitelist, n_constant=0, offset_lim=c(-4, 16),
  pos_lim=c(9, 23), len_lim=c(4, 11), corrected_kds=FALSE, sumseed=FALSE,
  supp_base=FALSE, offset_base=FALSE, site_base=NULL, kd_fc=TRUE
) {
  mm8mer_sites <- GetAll8merMmSites(mirna)
  # Load the Kds ###############################################################
  if (corrected_kds) {
    if (sitelist == "progthrp_suppcomp") rand_sitelist <- "randthrp_comp"
    else                                 rand_sitelist <- "randthrp"
    kds <- SubfunctionCall(
      ApplyKdCorrection, rand_n_constant=n_constant, prog_n_constant=n_constant,
      prog_sitelist=sitelist, rand_sitelist=rand_sitelist
    )
  } else {
    if (mirna == "miR-1" & experiment == "equilibrium") {
      buffer <- TRUE
      combined <- FALSE
    } else if (mirna == "miR-7-23nt" & experiment == "equilibrium2_nb") {
      combined <- FALSE
    }
    kds <- SubfunctionCall(EquilPars)
  }
  kds_global <<- kds
  if (!is.null(site_base) & !grepl("suppcomp", sitelist)) { # For one base site.
    kd_base <- kds[sprintf("%s_Kd", site_base), 2]
    base_str <- site_base   
  } else if (supp_base) {
    print("supp_base is positive")
    if (grepl("suppcomp", sitelist)) {
      if (sumseed) {
        kd_base <- kds["Supp_Kd", 2]
      } else {
        kd_base <- GeoMean(kds[sprintf("%s_Kd", kCanonicalSites), 2])
      }
      base_str <- "Supp"
    } else if (grepl("comp", sitelist)) {
      kd_base_canon <- kds[sprintf("%s_Kd", kSeedSites), 2]
      if (sumseed) {
        kd_comp <- kds["Comp_Kd", 2]
      } else {
        # kd_comp <- kds[""]
        kd_comp <- GeoMean(kds[sprintf("%s_Kd", mm8mer_sites), 2])
      }
      kd_base <- c(kd_base_canon, kd_comp)
      base_str <- c(sprintf("%s", kSeedSites), "Comp")
    }
  } else if (offset_base) {
    if (sumseed) {
      kd_base <- kds["Offset_Kd", 2]
    } else {
      kd_base <- GeoMean(kds[sprintf("%s_Kd", c("6mer-m8", "6mer-A1")), 2])
    }    
    base_str <- "Offset" 
  } else {
    if (grepl("comp", sitelist)) {
      if (sumseed) {
        kd_base <- kds["Comp_Kd", 2]
      } else {
        kd_base <- GeoMean(kds[sprintf("%s_Kd", mm8mer_sites), 2])
      }
      base_str <- "Comp"    
    } else {
      kd_base <- kds[sprintf("%s_Kd", mm8mer_sites), 2]
      base_str <- mm8mer_sites
    }
  }
  names(kd_base) <- base_str
  # if (supp_base & add_comp_sites) {
  #   if (collapsemm_addcomp) {
  #     kd_comp <- kds_comp["Comp_Kd", 2]
  #   } else {
  #     kd_comp <- GeoMean(kds_comp[sprintf("%s_Kd", GetAll8merMmSites(mirna)), 2])            
  #   }
  #   base_kds <- c(base_kds, kd_comp)
  #   names(base_kds)[length(base_kds)] <- "8mer-mm"
  # }
  if (!kd_fc){
    kd_base <- kd_base*0 + 1
  }
  output <- do.call("rbind", lapply(names(kd_base), function(base_site) {
    inds_grep <- grep(
      "\\.", grep(
        "&", grep(
          sprintf("^.*\\|.*\\|%s_Kd$", base_site), rownames(kds), perl=TRUE,
          value=TRUE
        ), invert=TRUE, value=TRUE
      ), value=TRUE
    )
    kds_use <- kds[inds_grep, ]
    s_pos <- as.integer(gsub("^.*\\|(.*)\\|.*_Kd$", replacement="\\1", rownames(kds_use)))
    pos5p <- as.integer(gsub("^.*m(.*)\\..*\\|.*\\|.*_Kd$", replacement="\\1", rownames(kds_use)))
    pos3p <- as.integer(gsub("^.*\\.(.*)\\|.*\\|.*_Kd$", replacement="\\1", rownames(kds_use)))
    len <- pos3p - pos5p + 1
    mm_type <- rep(base_site, length(pos5p))
    offset <- s_pos - pos5p
    logkd <- log10((kd_base[base_site])/(kds_use[, 2]))
    data.frame(`logkd`=logkd, `s_pos`=as.factor(s_pos),
               `pos5p`=pos5p, `pos3p`=pos3p,
               `pairing`=sprintf("%s|%s", pos5p, pos3p), `len`=len,
               `offset`=offset, `mm`=mm_type,
               `offsetXmm`=sprintf("%s|%s", offset, mm_type),
               `pos_5pXmm`=sprintf("%s|%s", pos5p, mm_type),
               `pos_3pXmm`=sprintf("%s|%s", pos3p, mm_type),
               `lenXmm`=sprintf("%s|%s", pos3p - pos5p + 1, mm_type))
  }))

  offset <- output$offset
  pos5p <- output$pos5p
  pos3p <- output$pos3p
  len <- output$len

  inds_keep <- which(offset >= offset_lim[1] & offset <= offset_lim[2] &
                     pos5p >= pos_lim[1] & pos3p <= pos_lim[2] &
                     len >= len_lim[1] & len <= len_lim[2])
  if (mirna == "miR-1") {
    inds_diff <- which(len == 4 & pos5p %in% c(15, 19))
  } else if (mirna == "miR-7-23nt") {
    inds_diff <- which(len == 4 & pos5p %in% c(17, 20))
  } else {
    inds_diff <- c()
  }
  inds_keep <- setdiff(inds_keep, inds_diff)
  output <- output[inds_keep, ]
  # print(head(output))
  # print(dim(output))
  # print("At the End of Make Kd DataFrame")
  return(output)
}

MakePairingMatrix <- function(
  mirna, experiment, offset, n_constant=3, sitelist="progthrp_suppcomp",
  len_lim=c(4, 11), pos_lim=c(9, 23), corrected_kds=FALSE, sumseed=FALSE,
  supp_base=FALSE, offset_base=FALSE, site_base=NULL, error=FALSE, loop=FALSE,
  kd_fc=TRUE
) {
  # Generate the site names.
  kd_data_frame <- SubfunctionCall(MakeKdDataFrame,
                                   offset_lim=c(offset, offset))

  # Define the possible 5-prime starting nucleotides and possible 3-prime
  # starting nucleotides, for overall matrix.
  len_mir <- nchar(kMirnaSeqs[mirna])
  nucs_5p <- 9:(len_mir - len_lim[1] + 1)
  nucs_3p <- (9 + len_lim[1] - 1):len_mir
  # Define the output matrix.
  out_matrix <- matrix(NA, nrow=length(nucs_3p), ncol=length(nucs_5p),
                       dimnames=list(nucs_3p, nucs_5p))
  for (i in 1:nrow(kd_data_frame)) {
    ind_5p <- as.character(kd_data_frame$pos5p[i])
    ind_3p <- as.character(kd_data_frame$pos3p[i])
    out_matrix[ind_3p, ind_5p] <- kd_data_frame$logkd[i]
  }
  return(out_matrix)
}


# MakeMismatchByOffsetMatrix <- function(
#   mirna, experiment, len, pos5p, n_constant=3, sitelist="progthrp",
#   corrected_kds=TRUE, kd_fc=TRUE, rand_base=FALSE, supp_base=FALSE,
#   offset_lim=c(-4, 16), add_comp_sites=FALSE, collapsemm_addcomp=FALSE
# ) {
#   kd_df <- SubfunctionCall(MakeKdDataFrame, len_lim=c(len, len),
#                                    pos_lim=c(pos5p, pos5p + len - 1))
#   if (supp_base) {
#     nrow_mat <- 7
#     row_names <- c(kSeedSites, "Comp")
#   } else {
#     nrow_mat <- 18
#     row_names <- GetAll8merMmSites(mirna)
#   }
#   output_matrix <- matrix(NA, nrow=nrow_mat, ncol=length(-4:16),
#                           dimnames=list(row_names, -4:16))
#   if (nrow(kd_df) != 0) {
#     for (i in 1:nrow(kd_df)) {
#       mm_use <- as.character(kd_df$mm)[i]
#       offset_use <- as.character(kd_df$offset)[i]
#       output_matrix[mm_use, offset_use] <- kd_df$logkd[i]
#     }
#   }  

#   if (rownames(output_matrix)[nrow(output_matrix)] == "Comp") {
#     rownames(output_matrix)[nrow(output_matrix)] <- "Seed mm"
#     output_matrix <- output_matrix[rev(1: nrow(output_matrix)), ]
#   }
#   return(output_matrix)
# }


# CalculateThreePMismatchRange <- function(
#   mirna, experiment, len, pos5p, n_constant=3, sitelist="progthrp",
#   corrected_kds=TRUE, kd_fc=TRUE, rand_base=FALSE, supp_base=FALSE,
#   offset_lim=c(-4, 16), add_comp_sites=FALSE, collapsemm_addcomp=FALSE,
#   win_average=1
# ) {
#   # Get a dataframe of mismatch by offset.
#   logkd_df <- SubfunctionCall(MakeMismatchByOffsetMatrix)
#   if (win_average != 1) {
#     logkd_df_new <- sapply(1:(ncol(logkd_df) - win_average + 1), function(start_col) {
#       rowMeans(logkd_df[, start_col:(start_col + win_average - 1)], na.rm=TRUE)
#     })
#     # Name the columns to reflect their average positions.
#     colnames(logkd_df_new) <- sapply(1:(ncol(logkd_df) - win_average + 1), function(start_col) {
#       mean(as.integer(colnames(logkd_df))[start_col:(start_col + win_average - 1)], na.rm=TRUE)
#     })
#     logkd_df <- logkd_df_new
#   }
#   # Calculate the range of each column, and subtract to get the fold-change
#   # across the column.
#   checks <- colMeans(logkd_df)
#   inds_use <- which(!is.na(checks))
#   if (length(inds_use) == 0) {
#     return(c(1))
#   } else {
#     logkd_df <- logkd_df[, inds_use, drop=FALSE]
#     print("dim:")
#     print(dim(logkd_df))
#     print(inds_use)
#     print(length(inds_use))
#     ranges <- apply(logkd_df, 2, range)
#     if (class(ranges) == "numeric") {
#       return(c(1))
#     } else {
#       ranges_dif <- ranges[2, ] - ranges[1, ]
#       # Return the exponentiated dataframe as the result.
#       return(10^ranges_dif)
#     }
#   }
# }



# CalculateThreePScore <- function(pos_lim, offset) {

#   pos_5p <- pos_lim[1]
#   pos_3p <- pos_lim[2]

#   if (pos_5p > 16 | pos_3p < 13) {
#     return(0)
#   } else {
#     positions <- seq(pos_5p, pos_3p)
#     print(positions)
#   }

# }




# GetEmpiricalBestSeedMismatchCoefficients <- function(
#   mirna, experiment, n_constant=3, offset_lim=c(-4, 16), len_lim=c(4, 11),
#   pos_lim=c(9, 23), pos_5p_lim=FALSE, offset_pick=NA, sitelist="progthrp", corrected_kds=TRUE,
#   supp_base=FALSE, kd_fc=TRUE
# ) {
#   print("GetEmpiricalBestSeedMismatchCoefficients")
#   out_mat <- do.call("cbind", lapply(len_lim[1]:len_lim[2], function(len_i) {
#     if (pos_5p_lim) pos_5p_lim_2 <- min(nchar(kMirnaSeqs[mirna]) - len_i + 1, pos_lim[2])
#     else            pos_5p_lim_2 <- min(nchar(kMirnaSeqs[mirna]), pos_lim[2]) - len_i + 1
#     do.call("cbind", lapply(seq(pos_lim[1], pos_5p_lim_2), function(pos_j) {
#       test <- SubfunctionCall(MakeMismatchByOffsetMatrix, len=len_i,
#                               pos5p=pos_j)
#       offsets_all <- as.integer(colnames(test))
#       test <- test[, which(offsets_all >= offset_lim[1] & offsets_all <= offset_lim[2])]
#       if (!is.na(offset_pick)) {
#         ind_use <- as.character(offset_pick)
#       } else {
#         ind_use <- which.max(colMeans(test, na.rm=TRUE))
#       }
#       out <- test[, ind_use, drop=FALSE]
#       colnames(out)[1] <- sprintf("l%sp%so%s", len_i, pos_j, colnames(out)[1])
#       out
#     }))  
#   }))
#   col_order <- order(apply(out_mat, 2, mean, na.rm=TRUE))
#   out_mat  <- out_mat[, rev(col_order)]
#   out_mat_check_1 <<- out_mat
#   out_mat <- t(t(out_mat) - colMeans(out_mat, na.rm=TRUE))
#   out_mat <- out_mat[, which(!is.na(colMeans(out_mat)))]
#   return(out_mat)
# }


# FitPairingAndOffsetModelSingle <- function(
#   mirna, experiment, n_constant=3, sitelist="progthrp_suppcomp",
#   corrected_kds=TRUE, sumseed=FALSE, offset_lim=c(-4, 16), len_lim=c(4, 11),
#   pos_lim=c(9, 23),  F_method=FALSE, F_exact=TRUE, fixed_offset=0,
#   use_global_df=FALSE, supp_base=FALSE, offset_base=FALSE, site_base=NULL,
#   kd_fc=TRUE, exponential=FALSE, intercept=FALSE, additive=FALSE,
#   lower_alt=TRUE
# ) {
#   # Make a dataframe concatenating all the offset matrices into one.
#   if (!use_global_df) {
#     df_all <<- SubfunctionCall(MakeKdDataFrame)
#   }

#   # Remove any levels that are not within the dataframe    
#   df_all$pairing <- as.character(df_all$pairing)
#   df_all$offset <- as.character(df_all$offset)

#   # This conditional should fix the singular matrix issue that I am experiencing
#   # (200801).
#   for(pairing_i in unique(df_all$pairing)) {
#     inds_i <- which(df_all$pairing == pairing_i)
#     df_sub <- df_all[inds_i, ]
#     df_sub_offsets <- unique(df_sub$offset)
#     if (length(setdiff(df_sub_offsets, c(as.character(fixed_offset)))) == 0) {
#       df_all <- df_all[-inds_i, ]
#     }
#   }
#   df_all <- drop.levels(df_all)

#   # Assign numbers of pairing and offset parameters.
#   n_pairing <- length(unique(df_all$pairing))
#   n_offset <- length(unique(df_all$offset)) - 1
#   # Subfunctions required to fit the model:
#   pars.pairing <- rep(0, n_pairing)
#   pars.offset <- rep(1, n_offset + 1)
#   names(pars.pairing) <- unique(as.character(df_all$pairing))
#   names(pars.offset) <- as.character(sort(as.numeric(unique(as.character(df_all$offset)))))
#   offset_names_use <- setdiff(names(pars.offset), as.character(fixed_offset))
#   AssignPars <- function(pars) {
#     # This function makes the parameters compatible with each of the modeling
#     # and gradient functions.
#     if (exponential) {
#       pars[1:(n_pairing + n_offset)] <- exp(pars[1:(n_pairing + n_offset)])
#     }
#     pars.pairing[names(pars.pairing)] <- pars[names(pars.pairing)]
#     pars.offset[offset_names_use] <- pars[offset_names_use]
#     if (intercept) {
#       pars.b <- pars[length(pars)]
#     } else {
#       pars.b <- 0
#     }
#     return(list(pairing=pars.pairing, offset=pars.offset, b=pars.b)) 
#   }
#   if (additive) {
#     Model <- function(pars) {
#       pars <- AssignPars(pars)
#       return(pars$pairing[df_all$pairing] + pars$offset[df_all$offset] + pars$b)
#     }

#     dModel_dPars <- function(pars) {
#       # Split up the parameters.
#       pars <- AssignPars(pars)
#       # Assign the dG and sigmoids for each row of data.
#       val.pairing <- pars$pairing[df_all$pairing]
#       val.offset <- pars$offset[df_all$offset]
#       #Pre-allocate the parameters.
#       out <- df_all$logkd*0
#       # Calculate the derivative of each of the pairing coefficients.
#       dM.dp_pairing <- sapply(1:n_pairing, function(ind) {
#         inds_par <- which(df_all$pairing == names(pars$pairing)[ind])
#         out[inds_par] <- 1
#         out
#       })
#       # Calculate the derivative of each of the offset coefficients.
#       dM.dp_offset <- sapply(offset_names_use, function(name) {
#         inds_par <- which(df_all$offset == name)
#         out[inds_par] <- 1
#         out
#       })
#       # Calculate the derivative of the base coefficient.
#       out <- cbind(dM.dp_pairing, dM.dp_offset)
#       colnames(out) <- c(names(pars$pairing), offset_names_use)
#       if (exponential) {
#         out <- t(apply(out, 1, function(row) {
#           row * c(pars$pairing, pars$offset[offset_names_use])  
#         }))
#       }
#       if (intercept) {
#         out <- cbind(out, rep(1, nrow(df_all)))
#         colnames(out)[ncol(out)] <- "b"
#       }
#       return(out)
#     }
#   } else {
#     Model <- function(pars) {
#       pars <- AssignPars(pars)
#       val_pairing <- pars$pairing[df_all$pairing]
#       val_offset <- pars$offset[df_all$offset]
#       return(val_pairing*val_offset + pars$b)
#     }

#     dModel_dPars <- function(pars) {
#       # Split up the parameters.
#       pars <- AssignPars(pars)
#       # Assign the dG and sigmoids for each row of data.
#       val.pairing <- pars$pairing[df_all$pairing]
#       val.offset <- pars$offset[df_all$offset]
#       #Pre-allocate the parameters.
#       out <- df_all$logkd*0
#       # Calculate the derivative of each of the pairing coefficients.
#       dM.dp_pairing <- sapply(1:n_pairing, function(ind) {
#         inds_par <- which(df_all$pairing == names(pars$pairing)[ind])
#         out[inds_par] <- val.offset[inds_par]
#         out
#       })
#       # Calculate the derivative of each of the offset coefficients.
#       offset_names_use <- setdiff(names(pars$offset), as.character(fixed_offset))
#       dM.dp_offset <- sapply(offset_names_use, function(name) {
#         inds_par <- which(df_all$offset == name)
#         out[inds_par] <- val.pairing[inds_par]
#         out
#       })
#       # Calculate the derivative of the base coefficient.
#       out <- cbind(dM.dp_pairing, dM.dp_offset)
#       colnames(out) <- c(names(pars$pairing), offset_names_use)
#       if (exponential) {
#         out <- t(t(out)*c(pars$pairing, pars$offset[offset_names_use]))
#       }
#       if (intercept) {
#         out <- cbind(out, rep(1, nrow(df_all)))
#         colnames(out)[ncol(out)] <- "b"
#       }
#       return(out)
#     }

#   }
#   # Loss function calculating the model error
#   tick <<- 0
#   LossFunction <- function(
#     pars, inds_zero_grad=c(), pars_zero_grad=c(), update_tick=TRUE
#   ) {
#     if (length(inds_zero_grad) != 0) {
#       pars[inds_zero_grad] <- pars_zero_grad
#     }
#     # Fit the model using the parameters
#     model_vals <- Model(pars)
#     SSE <- sum((df_all$logkd - model_vals)^2)
#     if (tick %% 10 == 0 & update_tick & additive) {
#       plot(model_vals, df_all$logkd, xlim=c(0, 3), ylim=c(0, 3))
#       segments(x0=0, y0=0, x1=3, y1=3)
#     }
#     if (update_tick) {
#       tick <<- tick + 1
#     }
#     return(SSE)
#   }

#   Gradient <- function(pars, inds_zero_grad=c(), pars_zero_grad=c()) {
#     if (length(inds_zero_grad) != 0) {
#       pars[inds_zero_grad] <- pars_zero_grad
#     }
#     # Split up the parameters.
#     model <- Model(pars)
#     dm_dp <- dModel_dPars(pars)
#     # Calculate the residuals.
#     dL_dM <- 2*(model - df_all$logkd)
#     dL_dp <- t(dL_dM) %*% dm_dp
#     # Determine the L2 regularization term:
#     if (length(pars_zero_grad) != 0) {
#       dL_dp[inds_zero_grad] <- 0
#     }
#     return(dL_dp)
#   }

#   Gradient_numerical <- function(pars, inds_zero_grad=c(), pars_zero_grad=c()) {
#     if (length(pars_zero_grad) != 0) {
#       pars[inds_zero_grad] <- pars_zero_grad
#     }
#     dL_dp <- grad(LossFunction, pars, update_tick=FALSE)
#     if (length(pars_zero_grad) != 0) {
#       dL_dp[inds_zero_grad] <- 0
#     }
#     return(dL_dp)
#   }

#   dModel_dPars_numerical <- function(pars) {
#     return(jacobian(Model, pars))
#   }
#   ################## Initialize the parameters ##################
#   pairing_names <- unique(as.character(df_all$pairing))
#   # Note 200821 - Had to add the "sort" beause now the df_all is made using a
#   # method that doesn't explicitly iterate over all the offset values, so they
#   # are now out of order.
#   offset_names <- sort(setdiff(unique(as.integer(df_all$offset)),
#                        as.character(fixed_offset)))
#   if (exponential) {
#     pars.init <- rep(-2, n_pairing + n_offset)
#   } else {
#     if (mirna == "let-7a" && length(site_base) != 0 && site_base == "7mer-A1") {
#       pars.init <- rep(2, n_pairing + n_offset)
#     } else {
#       pars.init <- rep(0, n_pairing + n_offset)
#     }
#   }

#   pars.init <- pars.init + rnorm(length(pars.init), mean=0, sd=0.01)
#   names(pars.init) <- c(names(pars.pairing), offset_names_use)
#   if (intercept) {
#     pars.init <- c(pars.init, 0)
#     names(pars.init)[length(pars.init)] <- "b"
#   }
#   ##### Get the MLE Estimates of the parameters ################
#   if (lower_alt & !additive) {
#     lower_use <- 0
#   } else {
#     lower_use <- -10
#   }
#   opt <- optim(pars.init, LossFunction, gr=Gradient,  method="L-BFGS-B",
#                lower=rep(lower_use, length(pars.init)),
#                upper=rep(10, length(pars.init)), control=list(maxit=1e7))

#   coefs <- opt$par
#   SSE_global <<- opt$val
#   model_values <- Model(opt$par)

#   ############### Get the error for all but fixed offset ##################
#   # Define degrees of freedom for error estimation:
#   n <- nrow(df_all)
#   p <- length(coefs)
#   if (F_method) {
#     F_test <- qf(0.95, 1, n - p)
#     LossCI <- function(par_val, par_i) {
#       coefs_init <- coefs
#       coefs[par_i] <- par_val
#       if (F_exact) {
#         opt <- optim(coefs, LossFunction, gr=Gradient_numerical, inds_zero_grad=c(par_i),
#                      pars_zero_grad=c(par_val),
#                      lower=rep(lower_use, length(coefs)), upper=rep(10, length(coefs)),
#                      method="L-BFGS-B", control=list(maxit=1e7))
#         SSE <- opt$value
#       } else {
#         SSE <- LossFunction(coefs)
#       }
#       rhs <- (n - p)*(SSE - SSE_global)/SSE_global
#       residual <- (F_test - rhs)^2
#       return(residual)
#     }
#     CIs <- matrix(NA, nrow=length(coefs), ncol=2,
#                   dimnames=list(names(coefs), c("LowerCI", "UpperCI")))

#     for (i in seq(1, length(coefs))) {
#       print(i)
#       coef_low <- optimize(LossCI, par_i=i, lower=c(coefs[i] - 2),
#                            upper=c(coefs[i]))$minimum
#       coef_high <- optimize(LossCI, par_i=i, lower=c(coefs[i]),
#                             upper=c(coefs[i] + 2))$minimum
#       CIs[names(coefs)[i], ] <- c(coef_low, coef_high)
#     }
#   } else {
#     dmodel.dpars <- dModel_dPars(coefs)
#     A_mat <- t(dmodel.dpars) %*% dmodel.dpars
#     sigma_2 <- SSE_global/(n - p)
#     A_mat <- A_mat
#     V_mat <- sigma_2*solve(A_mat)
#     range <- sqrt(diag(V_mat))*qt(0.975, n - p)
#     CIs <- cbind(coefs - range, coefs + range)
#     rownames(CIs) <- names(coefs)
#     colnames(CIs) <- c("LowerCI", "UpperCI")
#   }
#   ######### ########### Format parameters for output ###########################
#   coefs <- cbind(coefs, CIs)
#   colnames(coefs)[1] <- "MLE"
#   ##################### Add in the static parameter ############################
#   row_slot_in <- which(as.numeric(rownames(coefs)[(n_pairing + 1):nrow(coefs)]) == fixed_offset - 1)

#   if (length(row_slot_in) != 0) row_slot_in <- row_slot_in + n_pairing
#   if (exponential) {
#     vals_row <- c(0, 0, 0)
#   } else {
#     vals_row <- c(1, 1, 1)
#   }
#   if (length(row_slot_in) == 0) {
#     coefs <- rbind(coefs[1:n_pairing, ], rbind(vals_row, coefs[(n_pairing + 1):nrow(coefs), ]))
#     rownames(coefs)[n_pairing + 1] <- as.character(fixed_offset)
#   } else {
#       coefs_temp <- rbind(coefs[1:row_slot_in, ], vals_row)
#     if (row_slot_in != nrow(coefs)) {
#       coefs <- rbind(coefs_temp, coefs[(row_slot_in + 1):nrow(coefs), , drop=FALSE])
#     } else {
#       coefs <- coefs_temp
#     }
#     rownames(coefs)[row_slot_in + 1] <- as.character(fixed_offset)
#   }
#   # Extracting the pairing coefficients.
#   if (exponential) {
#     if (intercept) {
#       coefs[1:(nrow(coefs) - 1), ] <- exp(coefs[1:(nrow(coefs) - 1), ])
#     } else {
#       coefs <- exp(coefs)
#     }
#   }
#   coefs_pairing <- coefs[grep("\\|", rownames(coefs), value=TRUE, perl=TRUE), ]
#   # Remove the pairing coefficients.
#   coefs_offset <- coefs[grep("(\\||b)", rownames(coefs), perl=TRUE, invert=TRUE), ]
#   # Order the offset coefficients
#   coefs_offset <- coefs_offset[order(as.integer(rownames(coefs_offset))), ]
#   # This normalizes the parametesr such that the maximum offset coefficient is
#   # 1, while the parameter values can take any value.
#   # if (log_plus_one) {
#   #   coefs_pairing <- coefs_pairing + max(coefs_offset[, 1])
#   #   coefs_offset <- coefs_offset - max(coefs_offset[, 1])
#   # } else {
#     coefs_pairing <- coefs_pairing*max(coefs_offset[, 1])
#     coefs_offset <- coefs_offset/max(coefs_offset[, 1])
#   # }
#   if (intercept) {
#     message("intercept value:")
#     print(coefs[nrow(coefs), , drop=FALSE])
#     coefs <- rbind(coefs_pairing, coefs_offset, coefs[nrow(coefs), , drop=FALSE])
#   } else {
#     coefs <- rbind(coefs_pairing, coefs_offset)
#   }
#   return(list(`coefs`=coefs, `data`=df_all, `values`=model_values))
# }

# FitPairingAndOffsetModelWithError <- function(
#   mirna, experiment, n_constant=3, offset_lim=c(-4, 16), len_lim=c(4, 11),
#   pos_3p_min=c(9, 23), sitelist="progthrp_suppcomp", kd_fc=TRUE,
#   corrected_kds=TRUE, sumseed=FALSE, F_method=FALSE, supp_base=FALSE,
#   offset_base=FALSE, site_base=NULL, exponential=FALSE, intercept=FALSE,
#   additive=FALSE, lower_alt=TRUE
# ) {
#   ############ Solve the first offset and allocate the dataframes ##############
#   # print(mirna)
#   # print(experiment)
#   fit_start <- SubfunctionCall(FitPairingAndOffsetModelSingle,
#                                fixed_offset=offset_lim[1])
#   data <- fit_start$data
#   coefs <- fit_start$coefs
#   values <- fit_start$values
#   # Make matrices used for avereaging the upper- and lower- confidence intervals
#   # for all the parameters.
#   if ("" %in% rownames(coefs)) {
#     ind_change <- which(rownames(coefs) == "")
#     rownames(coefs)[ind_change] <- offset_lim[2]
#   }
#   lowerCI_mat <- matrix(
#     NA, nrow=nrow(coefs), ncol=offset_lim[2] - offset_lim[1] + 1,
#     dimnames=list(rownames(coefs), c(offset_lim[1]:offset_lim[2]))
#   )
#   upperCI_mat <- matrix(
#     NA, nrow=nrow(coefs), ncol=offset_lim[2] - offset_lim[1] + 1,
#     dimnames=list(rownames(coefs), c(offset_lim[1]:offset_lim[2]))
#   )
#   # Identify the rows to use.
#   rows_use <- setdiff(rownames(coefs), as.character(offset_lim[1]))
#   # print(rows_use)
#   # print(coefs)
#   # print(lowerCI_mat)
#   # print(coefs[rows_use, 2])
#   # print(setdiff(rownames(lowerCI_mat), rows_use))
#   # print(setdiff(rows_use, rownames(lowerCI_mat)))
#   # Use these indeces to add each of confidence intervals to.
#   lowerCI_mat[rows_use, as.character(offset_lim[1])] <- coefs[rows_use, 2]
#   upperCI_mat[rows_use, as.character(offset_lim[1])] <- coefs[rows_use, 3]
#   # Repeat by optimizing at each other fixed offset, and assign the values to
#   # the appropriate column of the matrices to be averaged.
#   for (offset_i in (offset_lim[1] + 1):offset_lim[2]) {
#     coefs_new <- SubfunctionCall(FitPairingAndOffsetModelSingle, use_global_df=TRUE,
#                                fixed_offset=offset_i)$coefs
#     coefs_new_global <<- coefs_new
#     rows_use <- setdiff(rownames(coefs_new), as.character(offset_i))
#     if (length(setdiff(rows_use, rownames(lowerCI_mat))) != 0) {
#       row_new <- length(setdiff(rows_use, rownames(lowerCI_mat)))
#       coefs_novel_name <- setdiff(rows_use, rownames(lowerCI_mat))
#       coefs <- rbind(
#         coefs,
#         matrix(NA, nrow=row_new, ncol=ncol(coefs),
#                dimnames=list(setdiff(rows_use, rownames(lowerCI_mat)),
#                              colnames(coefs)))
#       )
#       coefs[coefs_novel_name, 1] <- coefs_new[coefs_novel_name]
#       lowerCI_mat <- rbind(
#         lowerCI_mat,
#         matrix(NA, nrow=row_new, ncol=ncol(lowerCI_mat),
#                dimnames=list(setdiff(rows_use, rownames(lowerCI_mat)),
#                              colnames(lowerCI_mat)))
#       )
#       upperCI_mat <- rbind(
#         upperCI_mat,
#         matrix(NA, nrow=row_new, ncol=ncol(lowerCI_mat),
#                dimnames=list(setdiff(rows_use, rownames(upperCI_mat)),
#                              colnames(upperCI_mat)))
#       )
#     }
#     lowerCI_mat[rows_use, as.character(offset_i)] <- coefs_new[rows_use, 2]
#     upperCI_mat[rows_use, as.character(offset_i)] <- coefs_new[rows_use, 3]
#   }
#   # Assign the average of each of the UpperCI and LowerCI matrices to the
#   # appropriate columns of the starting offset dataframe, and return this
#   # object.

#   coefs[, 2] <- rowMeans(lowerCI_mat, na.rm=TRUE)
#   coefs[, 3] <- rowMeans(upperCI_mat, na.rm=TRUE)

#   ######################### Separate the coefficients ########################
#   coefs_pairing <- coefs[grep("\\|", rownames(coefs), value=TRUE, perl=TRUE), ]
#   coefs_offset <- coefs[grep("(\\||b)", rownames(coefs), perl=TRUE, invert=TRUE), ]
#   ##################### Make pairing coefficient matrices ######################
#   len_mir <- nchar(kMirnaSeqs[mirna])
#   len_k <- 4
#   nucs_5p <- 9:(len_mir - 3)
#   nucs_3p <- nucs_5p + 3
#   pairing_matrix <- matrix(NA, nrow=length(nucs_3p), ncol=length(nucs_5p),
#                        dimnames=list(nucs_3p, nucs_5p))

#   pairing_matrix_lowerCI <- matrix(NA, nrow=length(nucs_3p), ncol=length(nucs_5p),
#                        dimnames=list(nucs_3p, nucs_5p))

#   pairing_matrix_upperCI <- matrix(NA, nrow=length(nucs_3p), ncol=length(nucs_5p),
#                        dimnames=list(nucs_3p, nucs_5p))
#   for (i in 1:nrow(coefs_pairing)) {
#     coef_pairing <- coefs_pairing[i, ]
#     name_pairing <- rownames(coefs_pairing)[i]
#     name_split <- unlist(strsplit(name_pairing, split="\\|"))
#     rowname_val <- name_split[2]
#     colname_val <- name_split[1]
#     pairing_matrix[rowname_val, colname_val] <- coef_pairing[1]
#     pairing_matrix_lowerCI[rowname_val, colname_val] <- coef_pairing[2]
#     pairing_matrix_upperCI[rowname_val, colname_val] <- coef_pairing[3]   
#   }  
#   coefs_pairing <- list(`MLE`=pairing_matrix, `LowerCI`=pairing_matrix_lowerCI,
#                         `UpperCI`=pairing_matrix_upperCI)

#   # # Normalize for the maximum offvalue, such that the maximum ofset is `1`.
#   # coefs_offset <- coefs_offset/max(coefs_offset[, 1])
#   if (intercept) {
#     return(list(`pairing`=coefs_pairing, `offsets`=coefs_offset,
#                 b=coefs[nrow(coefs), , drop=FALSE], `data`=data,
#                 `values`=values))
#   } else {
#     return(list(`pairing`=coefs_pairing, `offsets`=coefs_offset, `data`=data,
#                 `values`=values))
#   }
# }

# DecomposePairingCoefficients <- function(
#   p_coefs, fit_max_len=5, keep_data=TRUE
# ) {
#   # Construct row-wise format matrix for the pairing coefficients, for ease of
#   # export alongside offset coefficients.
#   p_coefs_global <<- p_coefs
#   p_data <- matrix(p_coefs, ncol=1)
#   rownames(p_data) <- sprintf(
#     "%s|%s",
#     rep(colnames(p_coefs), each=nrow(p_coefs)),
#     rep(rownames(p_coefs), times=ncol(p_coefs))
#   )
#   colnames(p_data) <- c("MLE")
#   # Remove the rows with NA from the pairing coefficient matrix.
#   p_data <- p_data[!is.na(p_data[, 1]), , drop=FALSE]
#   # Build the design matrix.
#   # The number of columns in the matrix is the length of the miRNA - 8.
#   len_mir <- as.integer(rownames(p_coefs)[nrow(p_coefs)])
#   M <- matrix(0, nrow(p_data), len_mir - 8)
#   # Populate the design matrix.
#   for (i in 1:nrow(M)) {
#     limits <- as.integer(unlist(strsplit(rownames(p_data)[i], split="\\|"))) - 8
#     start <- limits[1]
#     stop <- limits[2]
#     M[i, limits[1]:limits[2]] <- 1
#   }
#   # Subset which parts of the matrix to use for the fitting, using the
#   # fit_max_len parameter.
#   ranges <- t(sapply(rownames(p_data), function(name_i) {
#     as.integer(unlist(strsplit(name_i, split="\\|"))  )
#   }))
#   lens <- ranges[, 2] - ranges[, 1] + 1
#   inds_use <- which(lens <= fit_max_len)
#   p_data <- p_data[inds_use, ]
#   M <- M[inds_use, ]
#   # Construct the least squares estimate from the matrix M and the y values.
#   Mt <- t(M)
#   inv_MtxM <- solve(Mt%*%M)
#   # Solve for the parameters (algebraically), and replace 0 values with 0.
#   params <- c(inv_MtxM%*%Mt%*%p_data)
#   names(params) <- 9:len_mir
#   # params[which(params < 0)] <- 0
#   # Define the lower and upper bound for the new output matrix.
#   min_nuc <- 9
#   max_nuc <- as.integer(names(params)[length(params)])
#   # Construct the output matrix.
#   p_coefs_pred <- matrix(
#     NA,
#     nrow=length(params), ncol=length(params),
#     dimnames=list(names(params), names(params))
#   )
#   # Iterate over the rows and columns to update the values of the matrix..
#   for (stop in rownames(p_coefs_pred)) {
#     for (start in colnames(p_coefs_pred)) {
#       # Determine the width of the parameter array to use for the summation.
#       len_pairing <- as.integer(stop) - as.integer(start) + 1
#       # Get the indeces of the params vector to sum over.
#       start_ind <- which(names(params) == start)
#       stop_ind <- which(names(params) == stop)
#       # Ensure that the stop index is greater than the start.
#       if (len_pairing >= 1) {
#         # If data exists.
#         if (len_pairing >= 4) {
#           if (!is.na(p_coefs[stop, start])) {
#             # If keep_data flag is true or if data is larger than
#             # what was used to fit.
#             if (keep_data | (len_pairing > fit_max_len)) {
#               p_coefs_pred[stop, start] <- p_coefs[stop, start]
#             } else {
#               p_coefs_pred[stop, start] <- sum(params[start_ind:stop_ind])
#             }
#           }
#         # If data doesn't exist, use model if below lower limit (i.e., 4 nt).
#         } else if (len_pairing < 4) {
#           p_coefs_pred[stop, start] <- sum(params[start_ind:stop_ind])
#         }
#       }
#     }
#   }
#   return(p_coefs_pred)
# }



# WritePairingAndOffsetModel <- function(
#   mirna, experiment, n_constant=3, offset_lim=c(-4, 16), len_lim=c(4, 11),
#   pos_3p_min=c(9, 23), sitelist="progthrp_suppcomp", kd_fc=TRUE,
#   corrected_kds=TRUE, sumseed=FALSE, F_method=FALSE, supp_base=FALSE,
#   offset_base=FALSE, site_base=NULL, exponential=FALSE, intercept=FALSE,
#   additive=FALSE, lower_alt=TRUE, decomp_pairing=FALSE, fit_max_len=5,
#   keep_data=TRUE
# ) {
#   # Fit model and extract the pairing and offset coefficients.
#   model <- SubfunctionCall(FitPairingAndOffsetModelWithError)
#   model_global <<- model
#   model <- model_global
#   p_coefs <- model$pairing$MLE
#   o_coefs <- model$offsets[, 1, drop=FALSE]
#   if (decomp_pairing) {
#     p_coefs <- SubfunctionCall(DecomposePairingCoefficients)
#   }
#   # Construct row-wise format matrix for the pairing coefficients, for ease of
#   # export alongside offset coefficients.
#   p_mat <- matrix(p_coefs, ncol=1)
#   rownames(p_mat) <- sprintf("%s|%s",
#                              rep(colnames(p_coefs), each=nrow(p_coefs)),
#                              rep(rownames(p_coefs), times=ncol(p_coefs)))
#   colnames(p_mat) <- colnames(o_coefs)
#   # Remove the rows with NA from the pairing coefficient matrix.
#   p_mat <- p_mat[!is.na(p_mat[, 1]), , drop=FALSE]
#   # Combine the two coefficient matrices.
#   coefs_all <- rbind(p_mat, o_coefs)
#   print(coefs_all)
#   coefs_all_global <<- coefs_all
#   # Define output path, make sure it exists, and define file name.
#   out_path <- sprintf("/lab/solexa_bartel/mcgeary/AgoRBNS/%s/%s/threep_models/", mirna, experiment)
#   EnsureDirectory(out_path)
#   final_path <- paste0(out_path, "pairing_and_offset")
#   if (decomp_pairing) {
#     if (keep_data)  keep_data.str <- "_keepdata"
#     else            keep_data.str <- ""
#     final_path <- sprintf("%s_decomp_max%s%s", final_path, fit_max_len,
#                           keep_data.str)
#   }
#   final_path <- paste0(final_path, ".txt")
#   print(final_path)
#   # Write coefficient matrix to output file path.
#   write.table(coefs_all, file=final_path, sep="\t", quote=FALSE)
# }


# FitPairingOffsetAndMismatchModelSingle <- function(
#   mirna, experiment, n_constant=3, corrected_kds=TRUE, offsetlim=c(-4, 16),
#   len_lim=c(4, 11), pos_lim=c(9, 23), sitelist="progthrp", sumseed=FALSE,
#   supp_base=FALSE, kd_fc=TRUE, fixed_offset=1, fixed_mm=1, intercept=FALSE,
#   exponential=FALSE, use_global_df=FALSE, F_method=FALSE, F_exact=TRUE,
#   lower_alt=TRUE, log_plus_one=FALSE, lambda_p=0, par_replacement_value=1,
#   pars_init=NULL
# ) {
#   # print("in FitPairingOffsetAndMismatchModel")
#   if (experiment %in% c("equilibrium", "equilibrium2_nb")) {
#     corrected_kds <- FALSE
#     sitelist <- "randthrp_comp"
#   }
#   # Make a dataframe concatenating all the offset matrices into one.
#   if (!use_global_df) {
#     df_all <<- SubfunctionCall(MakeKdDataFrame)
#   }

#   df_all <- drop.levels(df_all)
#   df_all$pairing <- as.character(df_all$pairing)
#   df_all$offset <- as.character(df_all$offset)
#   df_all$mm <- as.character(df_all$mm)

#   n_pairing <- length(unique(df_all$pairing))
#   n_offset <- length(unique(df_all$offset)) - 1
#   n_mm <- length(unique(df_all$mm)) - 1
#   mm_names_all <- unique(df_all$mm)
#   mm_names_all_global <<- mm_names_all
#   t_global <<- 0
#   pars.pairing <- rep(0, n_pairing)
#   names(pars.pairing) <- unique(df_all$pairing)
#   pars.offset <- rep(1, n_offset + 1)
#   names(pars.offset) <- as.character(sort(as.numeric(unique(df_all$offset))))
#   names_offset_use <- setdiff(names(pars.offset), as.character(fixed_offset))
#   pars.mm <- rep(1, n_mm + 1)
#   names(pars.mm) <- unique(df_all$mm)
#   names_mm_use <- names(pars.mm)[-fixed_mm]

#   AssignPars <- function(pars) {
#     # This function makes the parameters compatible with each of the modeling
#     # and gradient functions.
#     if (exponential) {
#       pars[1:(n_pairing + n_offset + n_mm)] <- exp(pars[1:(n_pairing + n_offset + n_mm)])
#     }
#     pars.pairing[names(pars.pairing)] <- pars[names(pars.pairing)]
#     pars.offset[names_offset_use] <- pars[names_offset_use]
#     pars.mm[names_mm_use] <- pars[names_mm_use]
#     # pars.offset_all <- pars[(n_pairing + 1):(n_pairing + n_offset)]
#     # if (n_mm != 0) {
#     #   pars.mm_pre <- pars[(n_pairing + n_offset + 1):(n_pairing + n_offset + n_mm)]
#     # } else {
#     #   pars.mm_pre <- c()
#     # }
#     if (intercept) {
#       pars.b <- pars[length(pars)]
#     } else {
#       pars.b <- 0
#     }
#     # Add in offset value.
#     # offset_vals <- as.numeric(names(pars.offset_all))
#     # pars.offset_min <- pars.offset_all[which(offset_vals < fixed_offset)]
#     # pars.offset_max <- pars.offset_all[which(offset_vals > fixed_offset)]
#     # pars.offset_fixed <- c(par_replacement_value)
#     # names(pars.offset_fixed) <- as.character(fixed_offset)
#     # pars.offset <- c(pars.offset_min, pars.offset_fixed, pars.offset_max)
#     # pars.offset <- pars.offset[order(as.integer(names(pars.offset)))]
#     # Add in mismatch value.
#     # pars.mm <- rep(NA, length(mm_names_all))
#     # names(pars.mm) <- mm_names_all
#     # if (log_plus_one) {
#     #   pars.mm[mm_names_all[fixed_mm]] <- 0
#     # } else {
#     #   pars.mm[mm_names_all[fixed_mm]] <- 1
#     # }
#     # pars.mm[names(pars.mm_pre)] <- pars.mm_pre
#     return(list(pairing=pars.pairing, offset=pars.offset, mm=pars.mm, b=pars.b)) 
#   }


#   if (log_plus_one) {
#     Model <- function(pars) {
#       pars <- AssignPars(pars)
#       # Calculate the dG and offset values for each datapoint.
#       val.pairing <- pars$pairing[df_all$pairing]
#       val.offset <- pars$offset[df_all$offset]
#       val.mm <- pars$mm[df_all$mm]
#       # Calculate the final model.
#       return(log10(exp(val.pairing + val.offset + val.mm) + 1))
#     }

#    dModel_dPars <- function(pars) {
#       # Split up the parameters.
#       pars <- AssignPars(pars)
#       # Assign the dG and sigmoids for each row of data.
#       val.pairing <- pars$pairing[df_all$pairing]
#       val.offset <- pars$offset[df_all$offset]
#       val.mm <- pars$mm[df_all$mm]
#       #Pre-allocate the parameters.
#       out <- df_all$logkd*0
#       Ka <- exp(val.pairing + val.offset + val.mm)
#       # Calculate the derivative of each of the pairing coefficients.
#       dM.dp_pairing <- sapply(1:n_pairing, function(ind) {
#         # print(names(pars$pairing)[ind])
#         inds_par <- which(df_all$pairing == names(pars$pairing)[ind])
#         # print(inds_par)
#         out[inds_par] <- 1/log(10)*(Ka/(Ka + 1))[inds_par]
#         # if (names(pars$pairing)[ind] == "12|21") {
#         #   print((val.offset*val.pairing*sqrt(df_all$weight))[inds_par])
#         #   print(out)
#         #   print(pars)
#         # }
#         return(out)
#       })
#       # Calculate the derivative of each of the offset coefficients.
#       offset_names_use <- setdiff(names(pars$offset), as.character(fixed_offset))
#       dM.dp_offset <- sapply(offset_names_use, function(name) {
#         inds_par <- which(df_all$offset == name)
#         out[inds_par] <- 1/log(10)*(Ka/(Ka + 1))[inds_par]
#         return(out)
#       })
#       mm_names_use <- setdiff(names(pars$mm), mm_names_all[fixed_mm])
#       dM.dp_mm <- sapply(mm_names_use, function(name) {
#         inds_par <- which(df_all$mm == name)
#         out[inds_par] <- 1/log(10)*(Ka/(Ka + 1))[inds_par]
#         return(out)
#       })

#       # Calculate the derivative of the base coefficient.
#       out <- cbind(dM.dp_pairing, dM.dp_offset, dM.dp_mm)
#       colnames(out) <- c(names(pars$pairing), offset_names_use, mm_names_use)
#       if (exponential) {
#         out <- t(apply(out, 1, function(row) {
#           row * c(pars$pairing, pars$offset[offset_names_use], pars$mm[mm_names_use])  
#         }))
#       }
#       if (intercept) {
#         out <- cbind(out, rep(1, nrow(df_all)))
#         colnames(out)[ncol(out)] <- "b"
#       }
#       return(out)
#     }

#   } else {
#     Model <- function(pars) {
#       # Initialize the parameters. Note that the first offset parameter is set to
#       # 1 in this function. This doesn't happen in the cost function, but it does
#       # happen in the gradient function.
#       pars <- AssignPars(pars)
#       # Calculate the dG and offset values for each datapoint.
#       val.pairing <- pars$pairing[df_all$pairing]
#       val.offset <- pars$offset[df_all$offset]
#       val.mm <- pars$mm[df_all$mm]
#       # Calculate the final model.
#       return(val.pairing*val.offset*val.mm + pars$b)
#     }

#     tick <- 0
#     # Loss function calculating the model error

#     dModel_dPars <- function(pars) {
#       # Split up the parameters.
#       pars <- AssignPars(pars)
#       # Assign the dG and sigmoids for each row of data.
#       val.pairing <- pars$pairing[df_all$pairing]
#       val.offset <- pars$offset[df_all$offset]
#       val.mm <- pars$mm[df_all$mm]
#       #Pre-allocate the parameters.
#       out <- df_all$logkd*0
#       # Calculate the derivative of each of the pairing coefficients.
#       dM.dp_pairing <- sapply(1:n_pairing, function(ind) {
#         # print(names(pars$pairing)[ind])
#         inds_par <- which(df_all$pairing == names(pars$pairing)[ind])
#         # print(inds_par)
#         out[inds_par] <- (val.offset*val.mm)[inds_par]
#         # if (names(pars$pairing)[ind] == "12|21") {
#         #   print((val.offset*val.pairing*sqrt(df_all$weight))[inds_par])
#         #   print(out)
#         #   print(pars)
#         # }
#         return(out)
#       })
#       # Calculate the derivative of each of the offset coefficients.
#       offset_names_use <- setdiff(names(pars$offset), as.character(fixed_offset))
#       dM.dp_offset <- sapply(offset_names_use, function(name) {
#         inds_par <- which(df_all$offset == name)
#         out[inds_par] <- (val.pairing*val.mm)[inds_par]
#         return(out)
#       })
#       mm_names_use <- setdiff(names(pars$mm), mm_names_all[fixed_mm])
#       dM.dp_mm <- sapply(mm_names_use, function(name) {
#         inds_par <- which(df_all$mm == name)
#         out[inds_par] <- (val.pairing*val.offset)[inds_par]
#         return(out)
#       })

#       # Calculate the derivative of the base coefficient.
#       out <- cbind(dM.dp_pairing, dM.dp_offset, dM.dp_mm)
#       colnames(out) <- c(names(pars$pairing), offset_names_use, mm_names_use)
#       if (exponential) {
#         out <- t(apply(out, 1, function(row) {
#           row * c(pars$pairing, pars$offset[offset_names_use], pars$mm[mm_names_use])  
#         }))
#       }
#       if (intercept) {
#         out <- cbind(out, rep(1, nrow(df_all)))
#         colnames(out)[ncol(out)] <- "b"
#       }
#       return(out)
#     }
#   }

#   LossFunction <- function(
#     pars, inds_zero_grad=c(), pars_zero_grad=c(), update_tick=FALSE
#   ) {
#     if (length(inds_zero_grad) != 0) {
#       pars[inds_zero_grad] <- pars_zero_grad
#     }
#     # Fit the model using the parameters
#     model_fit <- Model(pars)
#     model_fit <<- model_fit
#     if (tick %% 50 == 0 & update_tick) {
#       cols <- rep("black", nrow(df_all))
#       cols[which(as.character(df_all$mm) == "8mer-mm")] <- "blue"
#       plot(df_all$logkd, model_fit, col=cols)
#       segments(x0=-2, y0=-2, x1=3, y1=3, lty=2)
#     }
#     tick <<- tick + 1

#     SSE <- sum((df_all$logkd - model_fit)^2)

#     # pars <- AssignPars(pars)
#     # L2_error <- lambda_p*sum((pars$pairing - mean(pars$pairing))^2)
#     if (update_tick) {
#       tick <<- tick + 1
#     }
#     return(SSE)
#   }


#   Gradient <- function(
#     pars, inds_zero_grad=c(), pars_zero_grad=c()
#   ) {
#     # if (tick %% 10 == 0) {
#     #     grad_n <- Gradient_numerical(pars)
#     # }
#     if (length(inds_zero_grad) != 0) {
#       pars[inds_zero_grad] <- pars_zero_grad
#     }

#     model <- Model(pars)
#     dm_dp <- dModel_dPars(pars)
#     # Calculate the residuals.
#     dL_dM <- 2*(model - df_all$logkd)
#     dL_dp <- t(dL_dM) %*% dm_dp
#     # Determine the L2 regularization term:
#     # pars <- AssignPars(pars)
#     # par_names <- names(pars$pairing)
#     # dL_dp_L2 <- 2*lambda_p*(pars$pairing - mean(pars$pairing))
#     # dL_dp[1, par_names] <- dL_dp[1, par_names] + dL_dp_L2
 


#     # # Split up the parameters.
#     # model <- SubfunctionCall(Model)
#     # model_global <<- model
#     # pars <- AssignPars(pars)
#     # # pars_base <- exp(pars[length(pars)])
#     # # Assign the dG and sigmoids for each row of data.
#     # val.pairing <- pars$pairing[df_all$pairing]
#     # val.offset <- pars$offset[df_all$offset]
#     # val.mm <- pars$mm[df_all$mm]
#     # # Calculate the residuals.
#     # dL_dx <- 2*(model - df_all$logkd)
#     # # Calculate the derivative of each of the pairing coefficients.
#     # # dSSE.dpars <- matrix(0, nrow=nrow(df_all), ncol=n_pairing + n_offset)
#     # dSSE.dpars_pairing <- sapply(1:n_pairing, function(ind) {
#     #   inds_sum <- which(df_all$pairing == names(pars$pairing)[ind])
#     #   return(sum((dL_dx*val.offset*val.mm)[inds_sum]))
#     # })
#     # # Calculate the derivative of each of the offset coefficients.
#     # offset_names_use <- setdiff(names(pars$offset), as.character(fixed_offset))
#     # dSSE.dpars_offset <- sapply(offset_names_use, function(name) {
#     #   inds_sum <- which(df_all$offset == name)
#     #   return(sum((dL_dx*val.pairing*val.mm)[inds_sum]))
#     # })
#     # mm_names_use <- setdiff(names(pars$mm), mm_names_all[fixed_mm])
#     # dSSE.dpars_mm <- sapply(mm_names_use, function(name) {
#     #   inds_sum <- which(df_all$mm == name)
#     #   return(sum((dL_dx*val.pairing*val.offset)[inds_sum]))
#     # })
#     # # Concatenate the two components to make the entire gradient vector.
#     # dSSE.dpars <- c(dSSE.dpars_pairing, dSSE.dpars_offset, dSSE.dpars_mm)
#     # if (exponential) {
#     #   dSSE.dpars <- dSSE.dpars * c(pars$pairing, pars$offset[offset_names_use], pars$mm[mm_names_use])
#     # }
#     # names(dSSE.dpars) <- c(names(pars$pairing), offset_names_use, mm_names_use)
#     # if (intercept) {
#     #   dSSE.db <- c(`b`=sum(dL_dx))
#     #   dSSE.dpars <- c(dSSE.dpars, dSSE.db)
#     # }
#     if (length(pars_zero_grad) != 0) {
#       dL_dp[inds_zero_grad] <- 0
#     }
#     # if (tick %% 10 == 0) {
#     #     plot(dSSE.dpars, grad_n)
#     #     print("tick")
#     #     print(tick)
#     # }
#     # tick <<- tick + 1
#     return(dL_dp)
#   }

#   Gradient_numerical <- function(
#     pars, inds_zero_grad=c(), pars_zero_grad=c()
#   ) {
#     if (length(pars_zero_grad) != 0) {
#       pars[inds_zero_grad] <- pars_zero_grad
#     }
#     out <- grad(LossFunction, pars, update_tick=FALSE)
#     if (length(pars_zero_grad) != 0) {
#       out[inds_zero_grad] <- 0
#     }
#     return(out)
#   }

#   ################# Initialize the parameters
#   if (exponential) {
#     par_mean <- -2
#   } else {
#     par_mean <- 0.5
#   }
#   if (length(pars_init) != 0) {
#     pars.init <- pars_init[c(unique(df_all$pairing),
#                           sort(as.integer(setdiff(unique(df_all$offset),
#                                        as.character(fixed_offset)))),
#                           mm_names_all[-fixed_mm]), 1]
#     names(pars.init) <- c(unique(df_all$pairing),
#                           sort(as.integer(setdiff(unique(df_all$offset),
#                                        as.character(fixed_offset)))),
#                           mm_names_all[-fixed_mm])
#   } else if (intercept) {
#     pars.init <- c(rnorm(n_pairing + n_offset + n_mm, mean=par_mean, sd=0.1), -2)
#     names(pars.init) <- c(unique(df_all$pairing),
#                           sort(as.integer(setdiff(unique(df_all$offset),
#                                        as.character(fixed_offset)))),
#                           mm_names_all[-fixed_mm], "b")
#   } else {
#     pars.init <- rnorm(n_pairing + n_offset + n_mm, mean=par_mean, sd=0.1)
#     names(pars.init) <- c(unique(df_all$pairing),
#                           sort(as.integer(setdiff(unique(df_all$offset),
#                                        as.character(fixed_offset)))),
#                           mm_names_all[-fixed_mm])
#   }
#   if (exponential) {
#     lower_use <- rep(-10, length(pars.init))
#   } else if (intercept) {
#     lower_use <- rep(c(0, -10), times=c(length(pars.init) - 1, 1))
#   } else {
#     lower_use <- rep(0, length(pars.init))
#   }
#   # lower_use <- rep(-5, length(pars.init))
#   tick <<- 0

#   # grad_a <- Gradient(pars.init)

#   # grad_n <- Gradient_numerical(pars.init)

#   # dev.new(xpos=20, ypos=20, height=5, width=5)
#   # plot(grad_a, grad_n)
#   # dev.new(xpos=520, pos=20, height=5, width=5)
#   # Get the coefficients through the optimization routine.
#   opt <- optim(pars.init, LossFunction, gr=Gradient, method="L-BFGS-B",
#                  lower=lower_use,
#                  upper=rep(10, length(pars.init)), control=list(maxit=1e7))

#   coefs <- opt$par
#   SSE_global <- opt$val
#   model_values <- Model(opt$par)
#   ############### Get the error for all but fixed offset ##################
#   # Define degrees of freedom for error estimation:
#   n <- nrow(df_all)
#   p <- length(coefs)
#   if (F_method) {
#     F_test <- qf(0.95, 1, n - p)
#     LossCI <- function(par_val, par_i) {
#       coefs_init <- coefs
#       coefs[par_i] <- par_val
#       if (F_exact) {
#         opt <- optim(coefs, LossFunction, gr=Gradient, inds_zero_grad=c(par_i),
#                      pars_zero_grad=c(par_val),
#                      lower=lower_use, upper=rep(10, length(coefs)),
#                      method="L-BFGS-B", control=list(maxit=1e7))
#         SSE <- opt$value
#       } else {
#         SSE <- LossFunction(coefs)
#       }
#       rhs <- (n - p)*(SSE - SSE_global)/SSE_global
#       residual <- (F_test - rhs)^2
#       return(residual)
#     }
#     CIs <- matrix(NA, nrow=length(coefs), ncol=2,
#                   dimnames=list(names(coefs), c("LowerCI", "UpperCI")))

#     for (i in seq(1, length(coefs))) {
#       coef_low <- optimize(LossCI, par_i=i, lower=c(coefs[i] - 2),
#                            upper=c(coefs[i]))$minimum
#       coef_high <- optimize(LossCI, par_i=i, lower=c(coefs[i]),
#                             upper=c(coefs[i] + 2))$minimum
#       CIs[names(coefs)[i], ] <- c(coef_low, coef_high)
#     }
#   } else {
#     dmodel.dpars <- dModel_dPars(coefs)
#     dmodel.dpars <- dmodel.dpars
#     A_mat <- t(dmodel.dpars) %*% dmodel.dpars
#     sigma_2 <- SSE_global/(n - p)
#     A_mat <<- A_mat
#     V_mat <- sigma_2*solve(A_mat)
#     range <- sqrt(diag(V_mat))*qt(0.975, n - p)
#     CIs <- cbind(coefs - range, coefs + range)
#     rownames(CIs) <- names(coefs)
#     colnames(CIs) <- c("LowerCI", "UpperCI")
#   }
#   ######### ########### Format parameters for output ###########################
#   coefs <- cbind(coefs, CIs)
#   colnames(coefs)[1] <- "MLE"
#   ##################### Add in the static parameter ############################
#   row_slot_in <- which(as.numeric(rownames(coefs)[(n_pairing + 1):(n_pairing + n_offset)]) == fixed_offset - 1)

#   if (length(row_slot_in) != 0) row_slot_in <- row_slot_in + n_pairing
#   if (exponential) {
#     vals_row <- c(0, 0, 0)
#   } else {
#     vals_row <- c(1, 1, 1)
#   }
#   vals_row <- rep(par_replacement_value, 3)
#   if (length(row_slot_in) == 0) {
#     coefs <- rbind(coefs[1:n_pairing, ], rbind(vals_row, coefs[(n_pairing + 1):nrow(coefs), ]))
#     rownames(coefs)[n_pairing + 1] <- as.character(fixed_offset)
#   } else {
#       coefs_temp <- rbind(coefs[1:row_slot_in, ], vals_row)
#     if (row_slot_in != nrow(coefs)) {
#       coefs <- rbind(coefs_temp, coefs[(row_slot_in + 1):nrow(coefs), , drop=FALSE])
#     } else {
#       coefs <- coefs_temp
#     }
#     rownames(coefs)[row_slot_in + 1] <- as.character(fixed_offset)
#   }


#   row_slot_in_2 <- n_pairing + n_offset + fixed_mm
#   if (exponential) {
#     vals_row <- c(0, 0, 0)
#   } else {
#     vals_row <- c(1, 1, 1)
#   }
#   if (length(row_slot_in_2) == 0) {
#     coefs <- rbind(coefs[1:n_pairing, ], rbind(vals_row, coefs[(n_pairing + 1):nrow(coefs), ]))
#     rownames(coefs)[n_pairing + 1] <- as.character(fixed_offset)
#   } else {
#       coefs_temp <- rbind(coefs[1:row_slot_in_2, ], vals_row)
#     if (row_slot_in_2 != nrow(coefs)) {
#       coefs <- rbind(coefs_temp, coefs[(row_slot_in_2 + 1):nrow(coefs), , drop=FALSE])
#     } else {
#       coefs <- coefs_temp
#     }
#     rownames(coefs)[row_slot_in_2 + 1] <- mm_names_all[fixed_mm]
#   }

#   # Extracting the pairing coefficients.
#   if (exponential) {
#     if (intercept) {
#       coefs[1:(nrow(coefs) - 1), ] <- exp(coefs[1:(nrow(coefs) - 1), ])
#     } else {
#       coefs <- exp(coefs)
#     }
#   }
#   coefs_pairing <- coefs[grep("\\|", rownames(coefs), value=TRUE, perl=TRUE), ]
#   # Remove the pairing coefficients.
#   coefs_offset <- coefs[grep("(\\||b|mm)", rownames(coefs), perl=TRUE, invert=TRUE), ]
#   # Order the offset coefficients
#   coefs_offset <- coefs_offset[order(as.integer(rownames(coefs_offset))), ]
#   coefs_mm <- coefs[grep("mm", rownames(coefs), perl=TRUE), ]
#   # This normalizes the parametesr such that the maximum offset coefficient is
#   # 1, while the parameter values can take any value.

#   if (log_plus_one) {
#     coefs_pairing <- coefs_pairing + max(coefs_offset[, 1])
#     coefs_offset <- coefs_offset - max(coefs_offset[, 1])
#     coefs_pairing <- coefs_pairing + mean(coefs_mm[, 1])
#     coefs_mm <- coefs_mm - mean(coefs_mm[, 1])
#   } else {
#     coefs_pairing <- coefs_pairing*max(coefs_offset[, 1])
#     coefs_offset <- coefs_offset/max(coefs_offset[, 1])
#     coefs_pairing <- coefs_pairing*mean(coefs_mm[, 1])
#     coefs_mm <- coefs_mm/mean(coefs_mm[, 1])
#   }
#   if (intercept) {
#     message("intercept value:")
#     coefs <- rbind(coefs_pairing, coefs_offset, coefs_mm,
#                    coefs[nrow(coefs), , drop=FALSE])
#   } else {
#     coefs <- rbind(coefs_pairing, coefs_offset, coefs_mm)
#   }
#   return(list(`coefs`=coefs, `data`=df_all, `values`=model_values))
# }

# FitPairingOffsetAndMismatchModelWithError <- function(
#   mirna, experiment, n_constant=3, sitelist="progthrp", offset_lim=c(-4, 16),
#   len_lim=c(4, 11), pos_3p_min=c(9, 23), kd_fc=TRUE, corrected_kds=TRUE,
#   sumseed=FALSE, F_method=FALSE, F_exact=TRUE, supp_base=FALSE,
#   offset_base=FALSE, site_base=NULL, loop=FALSE, cutoff=FALSE,
#   exponential=FALSE, intercept=FALSE, log_plus_one=TRUE, lambda_p=0.1,
#   fullpairs=FALSE
# ) {
#   ############ Solve the first offset and allocate the dataframes ##############
#   time_start <- proc.time()[3]
#   fit_start <- SubfunctionCall(FitPairingOffsetAndMismatchModelSingle,
#                                fixed_offset=offset_lim[1], fixed_mm=1)
#   time_end <- proc.time()[3]

#   data <- fit_start$data
#   coefs <- fit_start$coefs
#   values <- fit_start$values
#   # Make matrices used for avereaging the upper- and lower- confidence intervals
#   # for all the parameters.
#   offset_names <- as.character(offset_lim[1]:offset_lim[2])
#   mm_names <- grep("mer", rownames(coefs), value=TRUE)
#   if (fullpairs) {
#     n_cols <- length(offset_names)*length(mm_names)
#     colnames_mat <- sprintf("%s&%s", rep(offset_names, each=length(mm_names)),
#                         rep(mm_names, length(offset_names)))
#   } else {
#     n_cols <- length(offset_names) + length(mm_names)
#     colnames_mat <- c(offset_names, mm_names)
#   }


#   lowerCI_mat <- matrix(
#     NA, nrow=nrow(coefs), ncol=n_cols,
#     dimnames=list(rownames(coefs), colnames_mat)
#   )
#   upperCI_mat <- matrix(
#     NA, nrow=nrow(coefs), ncol=n_cols,
#     dimnames=list(rownames(coefs), colnames_mat)
#   )
#   # Identify the rows to use.
#   rows_use <- setdiff(rownames(coefs), c(as.character(offset_lim[1]), mm_names[1]))
#   col_use <- 1
#   # Use these indeces to add each of confidence intervals to.
#   lowerCI_mat[rows_use, col_use] <- coefs[rows_use, 2]
#   upperCI_mat[rows_use, col_use] <- coefs[rows_use, 3]
#   # Repeat by optimizing at each other fixed offset, and assign the values to
#   # the appropriate column of the matrices to be averaged.
#   if (fullpairs) {
#     mm_start <- 2
#     for (offset_i in (offset_lim[1]):(offset_lim[2])) {
#       offset_time_start <- proc.time()[3]
#       for (mm_i in (mm_start:length(mm_names))) {
#         mm_time_start <- proc.time()[3]
#         coefs_new <- SubfunctionCall(FitPairingOffsetAndMismatchModelSingle, use_global_df=TRUE,
#                                    fixed_offset=offset_i, fixed_mm=mm_i)$coefs
#         # rows_use <- setdiff(rownames(coefs_new), as.character(offset_i))
#         rows_use <- setdiff(rownames(coefs_new), c(as.character(offset_i), mm_names[mm_i]))
#         col_use <- sprintf("%s&%s", offset_i, mm_names[mm_i])
#         # print(col_use)
#         lowerCI_mat[rows_use, col_use] <- coefs_new[rows_use, 2]
#         upperCI_mat[rows_use, col_use] <- coefs_new[rows_use, 3]
#         lowerCI_mat_global <<- lowerCI_mat
#         upperCI_mat_global <<- upperCI_mat
#         mm_time_stop <-proc.time()[3]
#         # message("Time for single optimization:")
#         # print(mm_time_stop - mm_time_start)
#       }
#       offset_time_stop <- proc.time()[3]
#       # message("Time for optimization of one offset with all mismatches:")
#       # print(offset_time_stop - offset_time_start)
#       mm_start <- 1
#     }  
#   } else {
#     mm_i <- 1
#     for (offset_i in ((offset_lim[1] + 1):offset_lim[2])) {
#       # print(offset_i)
#       time_start <- proc.time()[3]
#       coefs_new <- SubfunctionCall(FitPairingOffsetAndMismatchModelSingle, use_global_df=TRUE,
#                                  fixed_offset=offset_i, fixed_mm=mm_i)$coefs
#       time_end <- proc.time()[3]
#       # print(time_end - time_start)
#       # rows_use <- setdiff(rownames(coefs_new), as.character(offset_i))
#       rows_use <- setdiff(rownames(coefs_new), c(as.character(offset_i), mm_names))
#       col_use <- as.character(offset_i)

#       # print(col_use)
#       lowerCI_mat[rows_use, col_use] <- coefs_new[rows_use, 2]
#       upperCI_mat[rows_use, col_use] <- coefs_new[rows_use, 3]
#       lowerCI_mat_global <<- lowerCI_mat
#       upperCI_mat_global <<- upperCI_mat
#     }
#     offset_i <- offset_lim[1]
#     for (mm_i in (2:length(mm_names))) {
#       # print(mm_i)
#       time_start <- proc.time()[3]
#       coefs_new <- SubfunctionCall(FitPairingOffsetAndMismatchModelSingle, use_global_df=TRUE,
#                                  fixed_offset=offset_i, fixed_mm=mm_i, pars_init=coefs)$coefs
#       # rows_use <- setdiff(rownames(coefs_new), as.character(offset_i))
#       time_end <- proc.time()[3]
#       # print(time_end - time_start)
#       rows_use <- setdiff(rownames(coefs_new), c(as.character(offset_lim[1]:offset_lim[2]), mm_names[mm_i]))
#       col_use <- mm_names[mm_i]
#       # print(col_use)
#       lowerCI_mat[rows_use, col_use] <- coefs_new[rows_use, 2]
#       upperCI_mat[rows_use, col_use] <- coefs_new[rows_use, 3]
#       lowerCI_mat_global <<- lowerCI_mat
#       upperCI_mat_global <<- upperCI_mat
#     } 
#   }
#   # Assign the average of each of the UpperCI and LowerCI matrices to the
#   # appropriate columns of the starting offset dataframe, and return this
#   # object.
#   coefs[, 2] <- rowMeans(lowerCI_mat, na.rm=TRUE)
#   coefs[, 3] <- rowMeans(upperCI_mat, na.rm=TRUE)

#   ######################### Separate the coefficients ########################
#   coefs_pairing <- coefs[grep("\\|", rownames(coefs), value=TRUE, perl=TRUE), ]
#   coefs_offset <- coefs[grep("(\\||b|mer)", rownames(coefs), perl=TRUE, invert=TRUE), ]
#   coefs_mm <- coefs[grep("mer", rownames(coefs)), ]
#   ##################### Make pairing coefficient matrices ######################
#   len_mir <- nchar(kMirnaSeqs[mirna])
#   len_k <- 4
#   nucs_5p <- 9:(len_mir - 3)
#   nucs_3p <- nucs_5p + 3
#   pairing_matrix <- matrix(NA, nrow=length(nucs_3p), ncol=length(nucs_5p),
#                        dimnames=list(nucs_3p, nucs_5p))

#   pairing_matrix_lowerCI <- matrix(NA, nrow=length(nucs_3p), ncol=length(nucs_5p),
#                        dimnames=list(nucs_3p, nucs_5p))

#   pairing_matrix_upperCI <- matrix(NA, nrow=length(nucs_3p), ncol=length(nucs_5p),
#                        dimnames=list(nucs_3p, nucs_5p))
#   for (i in 1:nrow(coefs_pairing)) {
#     coef_pairing <- coefs_pairing[i, ]
#     name_pairing <- rownames(coefs_pairing)[i]
#     name_split <- unlist(strsplit(name_pairing, split="\\|"))
#     rowname_val <- name_split[2]
#     colname_val <- name_split[1]
#     pairing_matrix[rowname_val, colname_val] <- coef_pairing[1]
#     pairing_matrix_lowerCI[rowname_val, colname_val] <- coef_pairing[2]
#     pairing_matrix_upperCI[rowname_val, colname_val] <- coef_pairing[3]   
#   }  
#   coefs_pairing <- list(`MLE`=pairing_matrix, `LowerCI`=pairing_matrix_lowerCI,
#                         `UpperCI`=pairing_matrix_upperCI)

#   # # Normalize for the maximum offvalue, such that the maximum ofset is `1`.
#   # coefs_offset <- coefs_offset/max(coefs_offset[, 1])
#   if (intercept) {
#     return(list(`pairing`=coefs_pairing, `offsets`=coefs_offset, `mm`=coefs_mm,
#                 b=coefs[nrow(coefs), , drop=FALSE], `data`=data,
#                 `values`=values))
#   } else {
#     return(list(`pairing`=coefs_pairing, `offsets`=coefs_offset, `mm`=coefs_mm,
#       `data`=data, `values`=values))
#   }
# }





# FitPairingAndOffsetByMismatchModel <- function(
#   mirna, experiment, n_constant=3, offsetmin=-4, offsetmax=16, len_min=4,
#   len_max=11, pos_3p_max=23, pos_3p_min=9, sitelist="programmed",
#   corrected_kds=TRUE, collapsemm=FALSE, supp_base=FALSE, kd_fc=TRUE
# ) {
#   # Make a dataframe concatenating all the offset matrices into one.
#   df_all <- SubfunctionCall(MakeFullNucleotideAndMisMatchContributionDf)
#   lens <- as.integer(df_all$pos_3p) - as.integer(df_all$pos_5p) + 1
#   df_all <- df_all[which(lens >= len_min &
#                          lens <= len_max &
#                          as.integer(df_all$pos_3p) <= pos_3p_max &
#                          as.integer(df_all$pos_3p) >= pos_3p_min), ]
#   df_all <- drop.levels(df_all)
#   n_pars_5p3p <- length(unique(df_all$pos_5p3p))
#   n_pars_offsetXmm <- length(unique(df_all$offsetXmm))
#   LossFunction <- function(pars, printout=TRUE) {
#       # Initialize the parameters.
#     pars.5p3p <<- pars[1:n_pars_5p3p]
#     pars.offsetXmm <<- pars[(n_pars_5p3p + 1):(n_pars_5p3p + n_pars_offsetXmm)]
#     # Get the each of the model terms
#     dG_5p3p <<- pars.5p3p[as.character(df_all$pos_5p3p)]
#     dG_offsetXmm <<- pars.offsetXmm[as.character(df_all$offsetXmm)]
#     # Compute the model.
#     model_vals <<- dG_5p3p*dG_offsetXmm
#     # if (t_global %% 10 == 0) {
#     #   plot(model_vals, df_all$logkd, col=rgb(0, 0, 0, alpha=0.5))
#     #   print(cor(model_vals, df_all$logkd)^2)
#     # }
#     # t_global <<- t_global + 1
#     sum((df_all$logkd - model_vals)^2)
#   }
#   Gradient_numerical <- function(pars) {
#     return(grad(LossFunction, pars, printout=FALSE))
#   }
#   Gradient <- function(pars, printout=FALSE) {
#     # Assign the parameters.
#     pars.5p3p <- pars[1:n_pars_5p3p]
#     pars.offsetXmm <- pars[(n_pars_5p3p +  1):(n_pars_5p3p + n_pars_offsetXmm)]
#     # Assign the delta G values.
#     dG_5p3p <- pars.5p3p[as.character(df_all$pos_5p3p)]
#     dG_offsetXmm <- pars.offsetXmm[as.character(df_all$offsetXmm)]
#     # Calculate the model values.
#     model_vals <- dG_5p3p*dG_offsetXmm
#     dL_dmodel <<- 2*(model_vals - df_all$logkd)
#     p_start <- 1
#     p_stop <- n_pars_5p3p
#     grad.5p3p <- sapply(p_start:p_stop, function(ind) {
#       inds_sum <- which(as.character(df_all$pos_5p3p) == names(pars)[ind])
#       sum((dL_dmodel*dG_offsetXmm)[inds_sum])
#     })
#     p_start <- p_stop + 1
#     p_stop <- p_stop + n_pars_offsetXmm
#     grad.offsetXmm <- sapply(p_start:p_stop, function(ind) {
#       inds_sum <- which(as.character(df_all$offsetXmm) == names(pars)[ind])
#       sum((dL_dmodel*dG_5p3p)[inds_sum])
#     })
#     grad_all <- c(grad.5p3p, grad.offsetXmm)
#     names(grad_all) <- names(pars)
#     grad_all <<- grad_all
#     return(grad_all)
#   }
#   # Initialize the parameters
#   pars.init <- rep(1, n_pars_5p3p + n_pars_offsetXmm)
#   names(pars.init) <- c(unique(as.character(df_all$pos_5p3p)),
#                         unique(as.character(df_all$offsetXmm)))

#   # Get the coefficients through the optimization routine.
#   coefs <- optim(pars.init, LossFunction, gr=Gradient,
#                  method="L-BFGS-B", lower=rep(-10, length(pars.init)),
#                  upper=rep(20, length(pars.init)),
#                  control=list(maxit=1e7))$par
#   # Extract the pairing coefficients.
#   coefs_pairing <- coefs[grep("8mer", names(coefs), value=TRUE, perl=TRUE, invert=TRUE)]
#   # Remove the pairing coefficients from the list.
#   # coefs <- coefs[grep("\\|", names(coefs), perl=TRUE, invert=TRUE)]
#   # Extract the mm coefficients.
#   coefs_offsetXmm <- coefs[grep("8mer", names(coefs), value=TRUE)]

#   # Remove the mm coefficients.
#   # coefs_offset <- coefs[grep("8mer", names(coefs), value=TRUE, invert=TRUE)]
#   # Define the base coefficient.

#   # Extract the offset coefficients.
#   # Reorder the offset coefficients from -4 to +16.
#   # coefs_offset <- coefs_offset[order(as.integer(names(coefs_offset)))]
#   len_mir <- nchar(kMirnaSeqs[mirna])
#   # Define the possible 5-prime starting nucleotides and possible 3-prime
#   # starting nucleotides, for overall matrix.
#   len_k <- 4
#   nucs_5p <- 9:(len_mir - len_k + 1)
#   nucs_3p <- (9 + len_k - 1):len_mir
#   # Define the output matrix.
#   out_matrix <- matrix(NA, nrow=length(nucs_3p), ncol=length(nucs_5p),
#                        dimnames=list(nucs_3p, nucs_5p))

#   for (i in 1:length(coefs_pairing)) {
#     coef_pairing <- coefs_pairing[i]
#     name_pairing <- names(coefs_pairing)[i]
#     name_split <- unlist(strsplit(name_pairing, split="\\|"))
#     rowname_val <- name_split[2]
#     colname_val <- name_split[1]
#     out_matrix[rowname_val, colname_val] <- coef_pairing
#   }
#   coefs_pairing <- out_matrix

#   out_matrix_offsetXmm <- matrix(NA, nrow=18, ncol=length(offsetmin:offsetmax),
#                        dimnames=list(GetAll8merMmSites(mirna),
#                                      as.character(offsetmin:offsetmax)))
#   print(out_matrix_offsetXmm)
#   for (i in 1:length(coefs_offsetXmm)) {
#     coef_offsetXmm <- coefs_offsetXmm[i]
#     name_offsetXmm <- names(coefs_offsetXmm)[i]
#     name_split <- unlist(strsplit(name_offsetXmm, split="\\|"))
#     rowname_val <- name_split[2]
#     colname_val <- name_split[1]
#     print(rowname_val)
#     print(colname_val)
#     out_matrix_offsetXmm[rowname_val, colname_val] <- coef_offsetXmm
#   }
#   coefs_offsetXmm <- out_matrix_offsetXmm

#   # Normalize for the maximum offset coefficient value, such that the maximum
#   # offset is `1`.
#   coefs_pairing <- coefs_pairing*max(colMeans(coefs_offsetXmm))
#   coefs_offsetXmm <- coefs_offsetXmm/max(colMeans(coefs_offsetXmm))
#   # Normalize for the mean mm coefficient value, such that the mean mm
#   # coefficient is `1`.
#   # coefs_pairing <- coefs_pairing*mean(coefs_mm)
#   # coefs_mm <- coefs_mm/mean(coefs_mm)

#   return(list(`pairing`=coefs_pairing, `offsetXmm`=coefs_offsetXmm,
#               `data`=df_all))
# }


# GetSeedSiteInternalDeltaG <- function(
#   mirna, experiment, sitelist="randthrp_suppcomp", n_constant=3, A1=FALSE,
#   m8_sites=FALSE, remove_seed=TRUE
# ) {
#   source_python("general/general.py")
#   if (A1) {
#     start_i <- 1
#   } else {
#     start_i <- 2
#   }
#   # Get the string names of the sites.
#   mismatch_sites <- c("8mer", GetAll8merMmSites(mirna))
#   if (m8_sites) {
#     mismatch_8mer_sites <- sprintf("8mer-mm%s8", ConvertUtoT(setdiff(kNucs, WC_pairs[substr(mirna_seed_seq, nchar(mirna_seed_seq), nchar(mirna_seed_seq))])))
#     mismatch_sites <- c(mismatch_sites, mismatch_8mer_sites)
#   }
#   print(mismatch_sites)
#   # Get the Kds for the comparison.
#   if (mirna == "miR-1") {
#     buffer <- TRUE
#     combined <- FALSE
#   } else if (mirna == "miR-7-23nt") {
#     experiment <- "equilibrium2_nb"
#     combined <- FALSE
#   }
#   kds <- SubfunctionCall(EquilPars)
#   print(head(kds))
#   # Subset the kds to only include those sites that are the 8mer and seed
#   # mismatch sites.
#   inds_use <- sprintf("%s_Kd", mismatch_sites)
#   kds_use <<- kds[inds_use, 2]
#   names(kds_use) <- mismatch_sites
#   print(kds_use)
#   R <- 1.987e-3 # in kcal K-1 mol-1
#   T <- 310.15 # in K
#   dG_obs <- R*T*log(kds_use)
#   ddG_obs <- dG_obs - dG_obs["8mer"]
#   frac_dG_obs  <- dG_obs/dG_obs["8mer"]
#   frac_ddG_obs  <- -ddG_obs/dG_obs["8mer"]

#   WC_pairs <- c(`A`="U", `C`="G", `G`="C", `U`="A")

#   mirna_seed_seq <- substr(kMirnaSeqs[mirna], start_i, 8)
#   print(mirna_seed_seq)
#   mismatch_seqs <- sapply(mismatch_sites, GetSiteSeq, mirna=mirna)
#   if (!(A1)) {
#     mismatch_seqs <- sapply(mismatch_seqs, substr, start=1, stop=7)
#   }
#   print(mismatch_seqs)
#   dG_pred <- sapply(mismatch_seqs, function(mismatch_seq) {
#     main <- py_run_string(sprintf("deltaG = get_mfe('%s', '%s')", mirna_seed_seq, mismatch_seq))
#     main$deltaG
#   })
#   ddG_pred <- dG_pred - dG_pred["8mer"]
#   frac_dG_pred  <- dG_pred/dG_pred["8mer"]
#   frac_ddG_pred  <- -ddG_pred/dG_pred["8mer"]

#   mir_pos <- sapply(names(dG_pred)[-1], function(site_name) {
#     substr(site_name, nchar(site_name), nchar(site_name))
#   })
#   tar_nuc <- sapply(names(dG_pred)[-1], function(site_name) {
#     substr(site_name, nchar(site_name) - 1, nchar(site_name) - 1)
#   })
#   options(warn=2)
#   print("______________")
#   print(kMirnaSeqs[mirna])
#   print(mir_pos)
#   mir_nuc <- sapply(mir_pos, function(mir_pos_i) {
#     substr(kMirnaSeqs[mirna], mir_pos_i, mir_pos_i)  
#   })
#   print(mir_nuc)

#   mir_tar_nuc <- sprintf("m%st%s", mir_nuc, tar_nuc)
#   df_out <- data.frame(
#     `mirna`=rep(mirna, length(dG_pred)), `mir_pos`=c(NA, mir_pos),
#     `mir_nuc`=c(NA, mir_nuc), `tar_nuc`=c(NA, tar_nuc),
#     `mir_tar_nuc`=c(NA, mir_tar_nuc), `dG_obs`=dG_obs, `ddG_obs`=ddG_obs,
#     `frac_dG_obs`=frac_dG_obs, `frac_ddG_obs`=frac_ddG_obs, `dG_pred`=dG_pred,
#     `ddG_pred`=ddG_pred, `frac_dG_pred`=frac_dG_pred,
#     `frac_ddG_pred`=frac_ddG_pred
#   )
#   if (remove_seed) {
#     df_out <- df_out[-1, ]
#   }
#   return(df_out)
# }

# GetAllRandomSeedSiteInternalDeltaG <- function(
#   n_constant, A1=FALSE, m8_sites=FALSE, frac_dG=TRUE
# ) {
#   df_all <- do.call("rbind", lapply(c("let-7a", "miR-1", "miR-155", "miR-124", "lsy-6", "miR-7-23nt"), function(mirna) {
#     SubfunctionCall(GetSeedSiteInternalDeltaG)
#   }))
#   # print(head(df_all))
#   # print(tail(df_all))
#   # print(dim(df_all))
#   if (frac_dG) {
#     cols_use <- which(colnames(df_all) %in% c("frac_ddG_obs", "frac_ddG_pred"))
#   } else {
#     cols_use <- which(colnames(df_all) %in% c("ddG_obs", "ddG_pred"))
#   }
#   out <- aggregate(df_all[, cols_use], list(df_all$mir_tar_nuc), mean, na.rm=TRUE)
#   print(out)
#   return(out)
# }



# GetThreePSiteDeltaG <- function(mirna, start, stop, mm=FALSE, just_complement=TRUE, wobble=TRUE) {
#   source_python("general/RBNS_methods.py")
#   if (just_complement)  jc_string <- ", just_complement=True"
#   else                  jc_string <- ""
#   if (!wobble)  wob_string <- ", wobble=False"
#   else          wob_string <- ""
#   if (class(mm) == "character") {
#     main <- py_run_string(sprintf("deltaG = get_threep_site_deltaG('%s', %s, %s, mm='%s'%s%s)", mirna, start, stop, mm, jc_string, wob_string))
#   } else {
#     main <- py_run_string(sprintf("deltaG = get_threep_site_deltaG('%s', %s, %s%s%s)", mirna, start, stop, jc_string, wob_string))
#   }
#   main$deltaG
# }


# GetDeltaGPredictedMatrix <- function(mirna, just_complement=TRUE, wobble=TRUE) {
#   mirna_seq <- kMirnaSeqs[mirna]
#   # Assign the 5 prime and 3 prime nucleotides within the matrix
#   nucs_5p <- 9:(nchar(mirna_seq) - 4 + 1)
#   nucs_3p <- nucs_5p + 3
#   out_matrix <- matrix(NA, nrow=length(nucs_3p), ncol=length(nucs_5p),
#                        dimnames=list(nucs_3p, nucs_5p))
#   # Iterate over the lengths and starting positions, determine the stopping
#   # nucleotide, and calculate the deltaG.
#   for (len in 4:11) {
#     for (start in 9:(nchar(mirna_seq) - len + 1)) {
#       stop <- start + len - 1
#       delta_g <- SubfunctionCall(GetThreePSiteDeltaG)
#       out_matrix[as.character(stop), as.character(start)] <- delta_g
#     }
#   }
#   print(out_matrix)
#   return(out_matrix)
# }

# ProbKmerInRandom <- function(k, n) {
#   return(1 - (1 - 0.25^(k))^(n - k + 1))
# }

# GetMirna3pKmers <- function(mirna, k) {
#   seq_mirna <- kMirnaSeqs[mirna]
#   seq_mirna <- substr(seq_mirna, start=9, stop=nchar(seq_mirna))
#   mirna_kmers <- sapply(1:(nchar(seq_mirna) - k + 1), function(pos) {
#     substr(seq_mirna, start=pos, stop=pos + k - 1)  
#   })
#   return(unique(mirna_kmers))
# }

# NormalizeByLength <- function(pairing_matrix) {
#   print("in NormalizeByLength")
#   print(pairing_matrix[1:3, 1:3])
#   for (len_i in 4:11) {
#     print(len_i)
#     starts <- colnames(pairing_matrix)
#     print(starts)
#     stops <- as.character(as.integer(starts) + len_i - 1)
#     print(stops)
#     inds_use <- which(stops %in% rownames(pairing_matrix))
#     ind_pairs <- matrix(cbind(starts, stops)[inds_use, ], ncol=2)
#     print(ind_pairs)
#     print(dim(ind_pairs))
#     if (nrow(ind_pairs) != 0) {
#       mean_length <- apply(ind_pairs, 1, function(row) {
#         pairing_matrix[row[2], row[1]]
#       })
#       mean_length <- mean(mean_length)
#       for (row_i in 1:nrow(ind_pairs)) {
#         row <- ind_pairs[row_i, ]
#         pairing_matrix[row[2], row[1]] <- pairing_matrix[row[2], row[1]] - mean_length
#       }
#     }
#   }
#   return(pairing_matrix)
# }

# GetThreePrimeMmBDKds <- function(
#   mirna, experiment, start_mm, stop_mm, bulge=FALSE, n_constant=3,
#   site_base=FALSE, offset_lim=c(-4, 16), kd_fc=FALSE, win_average=1,
#   best_average=FALSE, new=TRUE, corrected_kds=TRUE, upper_bound=FALSE,
#   lower_bound=FALSE
# ) {
#   if (experiment %in% c("equilibrium", "equilibrium2_nb")) {
#     sitelist <- "randthrp_suppcomp"
#     corrected_kds <- FALSE
#     new <- FALSE
#     if (mirna == "miR-1" | mirna == "miR-7-23nt") combined <- FALSE
#     if (mirna == "miR-1") buffer <- TRUE
#   } else {
#     sitelist <- "progthrp_suppcomp"
#   }
#   if (corrected_kds) {
#     kds <- SubfunctionCall(ApplyKdCorrection, prog_n_constant=n_constant,
#                            prog_sitelist="progthrp_suppcomp")
#   } else {
#     kds <- SubfunctionCall(EquilPars)
#   }
#   kds_global <<- kds
#   # Build the string to search for the correct sites.
#   site_len <- stop_mm - start_mm + 1
#   if (bulge) {
#     str_mmbd <- "[bd]"
#   } else {
#     str_mmbd <- "mm"
#   }
#   if (class(site_base) == "character") {
#     base_str <- sprinf("%s_Kd", site_base)
#   } else {
#     base_str <- "Comp_Kd"
#   }
#   # Make vectors of the wt sites, and the mismatch/bulge sites with names.
#   str_grep_wt <- sprintf("^%smer-m%s.%s\\|.*\\|%s", site_len, start_mm,
#                          stop_mm, base_str)
#   if (lower_bound) col_use <- 3
#   else if (upper_bound) col_use <- 5
#   else             col_use <- 2
#   kds_wt <- kds[grep(str_grep_wt, rownames(kds)), col_use]
#   names(kds_wt) <- grep(str_grep_wt, rownames(kds), perl=TRUE, value=TRUE)
#   str_grep_mmb <- sprintf("^%smer-m%s.%s%s", site_len, start_mm, stop_mm, str_mmbd)
#   names_mmb <- grep(str_grep_mmb, rownames(kds), perl=TRUE, value=TRUE)
#   names_mmb <- grep(base_str, names_mmb, value=TRUE)
#   kds_mmb <- kds[names_mmb, col_use]
#   names(kds_mmb) <- names_mmb
#   names(kds_wt) <- gsub("^.*\\|(.*)\\|.*$", replacement="\\1", names(kds_wt),
#                         perl=TRUE)
#   # Make the first matrix, that is all of mismatch types aligned by position,
#   # in order to pick the offset range (3 long) that has the Kd in the wild type
#   # case.
#   # First change the rownames to only have the mismatch/bulge type and position
#   # of the site.
#   names(kds_mmb) <- gsub(sprintf("^%smer-m%s\\.%s%s(.*)\\|(.*)\\|Comp_Kd$",
#                                 site_len, start_mm, stop_mm, str_mmbd),
#                     replacement="\\1\\|\\2", names(kds_mmb))
#   mmb_sites <- sapply(names(kds_mmb), function(name_i) {
#     unlist(strsplit(name_i, split="\\|"))[1]
#   })
#   # Create the column name vector for the matrix to be made.
#   colnames_mat <- c("wt", unique(mmb_sites))
#   # Pre-allocate matrix.
#   wt_pos <- as.numeric(names(kds_wt))

#   matrix_rows <- as.character(seq(min(wt_pos), max(wt_pos)))
#   out_mat <- matrix(NaN, nrow=length(matrix_rows), ncol=length(colnames_mat),
#                     dimnames=list(matrix_rows, colnames_mat))
#   # Assign the wild type kds to the first column of the matrix.
#   out_mat[names(kds_wt), 1] <- kds_wt
#   # Populate the matrix with all of the values. The values being used are from
#   # a vector of al the mismatch-and-position kd values, and the output matrix is
#   # defined as being mismatch-by-position.
#   for (name_i in names(kds_mmb)) {
#     mmb_pos <- unlist(strsplit(name_i, split="\\|"))
#     if (mmb_pos[2] %in% rownames(out_mat)) {
#       out_mat[mmb_pos[2], mmb_pos[1]] <- kds_mmb[name_i]
#     }
#   }
#   # Make the matrix where rows are averaged together.
#   if (win_average != 1) {
#     out_mat_average <- 10^t(sapply(1:(nrow(out_mat) - win_average + 1), function(start_row) {
#       colMeans(log10(out_mat)[start_row:(start_row + win_average - 1), ], na.rm=TRUE)
#     }))
#     rownames(out_mat_average) <- sapply(1:(nrow(out_mat) - win_average + 1), function(start_row) {
#       mean(as.numeric(rownames(out_mat)[start_row:(start_row + win_average - 1)]))
#     })
#   } else {
#     out_mat_average <- out_mat
#   }

#   # Portion where offset limits are taken into account.
#   offsets_all <- as.integer(rownames(out_mat_average)) - start_mm
#   rownames(out_mat_average) <- offsets_all
#   offsets_min <- max(min(offsets_all), offset_lim[1])
#   offsets_max <- min(max(offsets_all), offset_lim[2])
#   offsets_use <- seq(offsets_min, offsets_max)
#   # Added this line in case there are offsets not present in the matrix that are
#   # part of the allowed range.
#   offsets_use <- intersect(offsets_use, rownames(out_mat_average))
#   out_mat_average <- out_mat_average[as.character(offsets_use), , drop=FALSE]
#   # Identify which row to use, and use that to make the vector with only
#   # that positional average.

#   # if (start_mm == 11 & stop_mm == 21) {
#   #   print(out_mat_average)
#   # }



#   if (best_average) {
#     ind_use <- which.min(out_mat_average[, 1])
#     # if (start_mm == 11 & stop_mm == 21) {
#     #   print(ind_use)
#     # }

#     # print("______________")
#     # print(mirna)
#     # print(start_mm)
#     # print(stop_mm)
#     # print(bulge)
#     # print("best offset:")
#     # print(rownames(out_mat_average)[ind_use])
#     names_out <- colnames(out_mat_average)
#     out_mat_average <- out_mat_average[ind_use, ]
#     names(out_mat_average) <- names_out
#   }
#   if (kd_fc) {
#     mm8mer_sites <- GetAll8merMmSites(mirna)
#     kd_ref <- GeoMean(kds[sprintf("%s_Kd", mm8mer_sites), 2])
#     out_mat_average <- kd_ref/out_mat_average    
#   }
#   return(out_mat_average)
# }


# CountAllMismatchSites <- function(mirna, experiment) {
#   len_mir <- nchar(kMirnaSeqs[mirna])
#   counts <- 0
#   for (len in 4:11) {
#     for (start_mm in 9:(len_mir - len + 1)) {
#       stop_mm <- start_mm + len - 1
#       # print(sprintf("%s-%s", start_mm, stop_mm))
#       kd_matrix <- SubfunctionCall(GetThreePrimeMmBDKds)
#       counts <- counts + length(!is.na(kd_matrix))
#       kd_matrix <- SubfunctionCall(GetThreePrimeMmBDKds, bulge=TRUE)
#       counts <- counts + length(!is.na(kd_matrix[, -1]))
#     }
#   }
#   return(counts)
# }




# GetThreePModelMismatchCoefficients <- function(
#   mirna, experiment, start_mm, stop_mm, bulge=FALSE, mm_and_bulge=FALSE,
#   n_constant=3, win_average=1, offset_lim=c(-4, 16), new=TRUE, print=FALSE,
#   corrected_kds=TRUE, kd_fc=FALSE, modelweights=FALSE
# ) {
#   # Pull out the data
#   if (mm_and_bulge) {
#     data_mm <- log10(SubfunctionCall(GetThreePrimeMmBDKds, bulge=FALSE,
#                                      kd_fc=TRUE))
#     data_bd <- log10(SubfunctionCall(GetThreePrimeMmBDKds, bulge=TRUE,
#                                      kd_fc=TRUE))[, -1]
#     inds_b <- grep("[ACGT]", colnames(data_bd), perl=TRUE)
#     inds_d <- grep("[ACGT]", colnames(data_bd), perl=TRUE, invert=TRUE)
#     colnames(data_bd)[inds_b] <- paste0("b", colnames(data_bd)[inds_b])
#     colnames(data_bd)[inds_d] <- paste0("d", colnames(data_bd)[inds_d])
#     data <- cbind(data_mm, data_bd)

#   } else {
#     data <- log10(SubfunctionCall(GetThreePrimeMmBDKds, kd_fc=TRUE))
#   }
#   offsets <- as.integer(rownames(data))
#   inds_keep <- which(offsets >= offset_lim[1] & offsets <= offset_lim[2])
#   data <- data[inds_keep, ]
#   # Define Cost function for minimization.
#   tick <- 0
#   CostFunction <- function(pars, print=FALSE) {
#     offsets <- pars[1:nrow(data)]
#     pairing <- pars[(nrow(data) + 1):length(pars)]
#     # Multiply the column vector by the row vector to get the model matrix.
#     model <- offsets%*%t(pairing)
#     if (tick %% 20 == 0 & print) {
#       plot(c(model), c(data), xlim=c(-1, 3), ylim=c(-1, 3))
#     }
#     tick <<- tick + 1
#     # Return the sum of squared residuals
#     if (modelweights) {
#       return(sum(abs(data)*(data - model)^2, na.rm=TRUE))
#     } else {
#       return(sum((data - model)^2, na.rm=TRUE))
#     }
#   }
#   # Initilize the parameters ###################################################
#   pars_init <- rep(0.2, nrow(data) + ncol(data))
#   names(pars_init) <- c(rownames(data), colnames(data))
#   # Open a plot window if print=TRUE. ##########################################
#   if (print) dev.new(xpos=20, ypos=20, height=5, width=5)
#   lower <- rep(0, length(pars_init))
#   upper <- rep(10, length(pars_init))
#   pars <- optim(par=pars_init, CostFunction, print=print, method="L-BFGS-B",
#                lower=lower, upper=upper, control=list(maxit=10000))$par
#   # Extract the parameters
#   offsets <- pars[1:nrow(data)]
#   pairing <- pars[(nrow(data) + 1):length(pars)]
#   # Renormalize the parameters such that `offsets` is bounded by 1. ############
#   pairing <- pairing*max(offsets, na.rm=TRUE)
#   offsets <- offsets/max(offsets, na.rm=TRUE)
#   # Name the parameters according to the data matrix. ##########################
#   names(offsets) <- rownames(data)
#   names(pairing) <- colnames(data)
#   print(pairing)
#   if (mm_and_bulge) {
#     if (bulge) {
#       inds_keep <- c(1, grep("[bd]", names(pairing), perl=TRUE))
#       pairing <- pairing[inds_keep]
#       names(pairing) <- gsub("b", names(pairing), replacement="")
#       names(pairing) <- gsub("d", names(pairing), replacement="")
#     } else {
#       inds_keep <- grep("[bd]", names(pairing), perl=TRUE, invert=TRUE)
#       pairing <- pairing[inds_keep]
#     }
#   }
#   print(pairing)
#   print(bulge)
#   return(list(`offsets`=offsets, `pairing`=pairing))
# }






# GetModelAIC <- function(
#   mirna, experiment, n_constant=3, sitelist="programmed", modelnew=FALSE, offsetmin=-4, offsetmax=16,
#   xpos=20
# ) {
#   if (modelnew) {
#     model <- SubfunctionCall(FitPairingAndOffsetByMismatchModel)
#     pairing <- model$pairing
#     offsetXmm <- model$offsetXmm
#     num_pars <- length(c(pairing)) + length(c(offsetXmm))
#     data <- model$data
#     model_sim <- apply(data, 1, function(row) {
#       pairing[as.character(as.integer(row[4])), as.character(as.integer(row[3]))]*offsetXmm[row[7], row[6]]
#     })
#     model_sim_new <<- model_sim
#     model_sim
#   } else {
#     model <- SubfunctionCall(FitPairingOffsetAndMismatchModel)  
#     offsets <- model$offsets
#     pairing <- model$pairing
#     mm <- model$mm
#     num_pars <- length(offsets) + length(c(pairing)) + length(mm)
#     data <- model$data
#     model_sim <- apply(data, 1, function(row) {
#       pairing[as.character(as.integer(row[4])), as.character(as.integer(row[3]))]*mm[row[7]]*offsets[row[6]]
#     })
#     model_sim_old <<- model_sim
#     model_sim
#   }
#   SSE <- sum((data$logkd - model_sim)^2)
#   n <- length(data$logkd)
#   logL <- -n/2*(log(2*pi*SSE/n) + 1)
#   AIC_out <- 2*num_pars - 2*logL
#   return(AIC_out)
# }

# MakeSiteMismatchMatrix <- function(
#   mirna, experiment, start_mm, stop_mm, n_constant=3, bulge=FALSE,
#   model_values=FALSE, offset_lim=c(-4, 16), new=TRUE, kd_fc=TRUE, win_average=3,
#   corrected_kds=TRUE, deltaG=FALSE
# ) {
#   if (model_values) {
#     parameters <- SubfunctionCall(GetThreePModelMismatchCoefficients)
#     kdfc_mmdb <- parameters[[2]]
#     if (deltaG) {
#       dG_base <- GetThreePSiteDeltaG(mirna, start_mm, stop_mm)
#       dG_mm <- sapply(names(kdfc_mmdb)[-1], GetThreePSiteDeltaG, mirna=mirna,
#                   start=start_mm, stop=stop_mm)
#       dG_all <- c(dG_base, dG_mm)
#       names(dG_all) <- names(kdfc_mmdb)
#       kdfc_mmdb <- dG_all
#     }
#   } else {
#     kdfc_mmdb <- log10(SubfunctionCall(GetThreePrimeMmBDKds, best_average=TRUE))
#   }

#   # make the new matrix that will be used to visualize the mismatches/bulges,
#   # where the rows are the four possible nucleotides, and the columns are the
#   # positions of the mismatch or bulge.
#   if (bulge) {
#     ncol_mat <- stop_mm - start_mm
#     colnames_mat <- as.character((start_mm + 1):stop_mm)
#   } else {
#     ncol_mat <- stop_mm - start_mm + 1
#     colnames_mat <- as.character(start_mm:stop_mm)
#   }
#   if (bulge) {
#     nrow_mat <- 5
#     rownames_mat <- c("A", "C", "G", "T", "d")
#   } else {
#     nrow_mat <- 4
#     rownames_mat <- c("A", "C", "G", "T")
#   }
#   # Pre-allocate the output values of the matrix.
#   logkdfc_mat <- matrix(NaN, nrow=nrow_mat, ncol=ncol_mat,
#                         dimnames=list(rownames_mat, colnames_mat))
#   print(names(kdfc_mmdb))
#   for (mm_i in names(kdfc_mmdb)[-1]) {
#     nuc <- substr(mm_i, start=1, stop=1)
#     if (!(nuc %in% c("A", "C", "G", "T"))) {
#       nuc <- "d"
#       par_pos <- 1
#     } else {
#       par_pos <- 2
#     }
#     if (substr(mm_i, par_pos, par_pos) == "(") {
#       splits <- unlist(strsplit(unlist(strsplit(mm_i, split="\\("))[2],
#                        split="\\)"))[1]
#       range <- unlist(strsplit(splits, split="\\."))
#       pos_list <- as.character(seq(as.integer(range[1]), as.integer(range[2])))
#       for (pos in pos_list) {
#         logkdfc_mat[nuc, pos] <- kdfc_mmdb[mm_i]
#       }
#     } else {
#       pos <- substr(mm_i, start=par_pos, stop=nchar(mm_i))
#       logkdfc_mat[nuc, pos] <- kdfc_mmdb[mm_i]
#     }
#   }
#   if (!bulge) {
#     for (i in colnames(logkdfc_mat)) {
#       mir_nuc <- RevComplement(substr(kMirnaSeqs[mirna], start=i, stop=i))
#       print(mir_nuc)
#       logkdfc_mat[mir_nuc, i] <- kdfc_mmdb["wt"]
#     }
#   }
#   rownames(logkdfc_mat)[4] <- "U"
#   if (bulge) rownames(logkdfc_mat)[5] <- "Del"
#   return(logkdfc_mat)
# }

# GetDeltaDeltaGResidual <- function(
#   mirna, experiment, start_mm, stop_mm, n_constant=3, bulge=FALSE,
#   model_values=FALSE, new=TRUE, offset_lim=c(-4, 16), kd_fc=TRUE, win_average=3,
#   corrected_kds=TRUE, obs_only=FALSE, obs_and_pred=FALSE, frac_dG=FALSE,
#   error=FALSE
# ) {
#   if (model_values) {
#     parameters <- SubfunctionCall(GetThreePModelMismatchCoefficients)
#     kdfc_mmdb <- parameters[[2]]
#   } else {
#     kdfc_mmdb <- log10(SubfunctionCall(GetThreePrimeMmBDKds, best_average=TRUE))
#     if (error) {
#       kdfc_mmdb_uCI <- log10(SubfunctionCall(GetThreePrimeMmBDKds, best_average=TRUE, lower_bound=TRUE))
#       kdfc_mmdb_lCI <- log10(SubfunctionCall(GetThreePrimeMmBDKds, best_average=TRUE, upper_bound=TRUE))
#       errors <- rowMeans(cbind(kdfc_mmdb_uCI - kdfc_mmdb, kdfc_mmdb - kdfc_mmdb_lCI))
#     }
#   }
#   # Place where ambiguous positions (relevant only to bulges and deletions) are
#   # redundantly represented for later use (i.e., a T(16.17) will be represented
#   # twice in the list, once as a T16, and once as a T17).
#   print(kdfc_mmdb)
#   inds_ambig <- grep("\\(", names(kdfc_mmdb), perl=TRUE)
#   print(inds_ambig)
#   target <- "^(.*)\\((.*)\\.(.*)\\)$"
#   for (i in inds_ambig) {
#     name_i <- names(kdfc_mmdb)[i]
#     nuc <- gsub(target, replacement="\\1", name_i)
#     pos_l <- gsub(target, replacement="\\2", name_i)
#     pos_r <- gsub(target, replacement="\\3", name_i)
#     seqs_use <- as.character(seq(as.integer(pos_l), as.integer(pos_r)))
#     vals_add <- rep(kdfc_mmdb[i], length(seqs_use))
#     names(vals_add) <- sprintf("%s%s", nuc, seqs_use)
#     kdfc_mmdb <- c(kdfc_mmdb, vals_add)
#     if (error) {
#       errors_add <- rep(errors[i], length(seqs_use))
#       names(errors_add) <- names(vals_add)
#       errors <- c(errors, errors_add)
#     }

#   }
#   if (length(inds_ambig) != 0) {
#     kdfc_mmdb <- kdfc_mmdb[-inds_ambig]
#     if (error) {
#       errors <- errors[-inds_ambig]
#     }
#   }
#   R <- 1.987e-3 # in kcal K-1 mol-1
#   T <- 310.15 # in K
#   dG_obs <- -R*T*log(10)*kdfc_mmdb
#   ddG_obs <- dG_obs - dG_obs["wt"]
#   frac_dG_obs  <- dG_obs/dG_obs["wt"]
#   frac_ddG_obs  <- -ddG_obs/dG_obs["wt"]
#   if (error) {
#     dG_obs_errors <- R*T*log(10)*errors
#     ddG_obs_errors <- sqrt(dG_obs_errors[-1]^2 + dG_obs_errors["wt"]^2)
#   }
#   if (!obs_only) {
#     dG_base <- GetThreePSiteDeltaG(mirna, start_mm, stop_mm)
#     dG_mm <- sapply(names(kdfc_mmdb)[-1], GetThreePSiteDeltaG, mirna=mirna,
#                     start=start_mm, stop=stop_mm)
#     dG_pred <- c(dG_base, dG_mm)
#     names(dG_pred)[1] <- "wt"
#     ddG_pred <- dG_pred - dG_pred["wt"]
#     frac_dG_pred <- dG_pred/dG_pred["wt"]
#     frac_ddG_pred <- -ddG_pred/dG_pred["wt"]
#   }
#   # Adding in bulges at multiple positions

#   # Defining the strings associated with the type of mismatch at each position
#   target_nucs <- sapply(names(dG_obs)[-1], substr, start=1, stop=1)
#   target_nucs[!(target_nucs %in% c("A", "C", "G", "T"))] <- "Del."
#   pos <- gsub("(A|C|G|T)", replacement="", names(dG_obs)[-1])
#   miR_seq_list <- unlist(strsplit(substr(kMirnaSeqs[mirna], start_mm, stop_mm),
#                          split=""))
#   names(miR_seq_list) <- start_mm:stop_mm
#   pos_use <- gsub("^\\(*(.*)\\..*$", replacement="\\1", pos)
#   miR_nucs <- miR_seq_list[pos_use]
#   mm_type <- target_nucs
#   if (!bulge) {
#     wG_inds <- which(miR_nucs == "U" & target_nucs == "G")
#     wU_inds <- which(miR_nucs == "G" & target_nucs == "T")
#     mm_type[wG_inds] <- "wG"
#     mm_type[wU_inds] <- "wU"
#   }
#   if (bulge) {    
#     miR_pairing <- sprintf("p%sm%stb%s", pos, miR_nucs,
#                            gsub("T", target_nucs, replacement="U"))    
#     inds_del <- which(target_nucs == "Del.")
#     miR_pairing[inds_del] <- sprintf("p%sm%sdel", pos, miR_nucs)[inds_del]
#   } else {
#     miR_pairing <- sprintf("p%sm%st%s", pos, miR_nucs,
#                            gsub("T", target_nucs, replacement="U"))    
#   }
#   if (obs_only) {
#     if (frac_dG) {
#       output_matrix <- data.frame(`mm_type`=mm_type, `pos`=pos,
#                                   `nuc_pair`=miR_pairing, `ddG_obs`=frac_dG_obs[-1])
#     } else {
#       output_matrix <- data.frame(`mm_type`=mm_type, `pos`=pos,
#                                   `nuc_pair`=miR_pairing, `ddG_obs`=ddG_obs[-1])
#     }
#     if (error) {
#       output_matrix <- cbind(output_matrix, ddG_obs_errors)
#     }
#   } else if (obs_and_pred) {
#     if (frac_dG) {
#     output_matrix <- data.frame(`mm_type`=mm_type, `pos`=pos,
#                                 `nuc_pair`=miR_pairing, `ddG_obs`=frac_dG_obs[-1],
#                                 `ddG_pred`=frac_dG_pred[-1])
#     } else {
#       output_matrix <- data.frame(`mm_type`=mm_type, `pos`=pos,
#                                   `nuc_pair`=miR_pairing, `ddG_obs`=ddG_obs[-1],
#                                   `ddG_pred`=ddG_pred[-1])
#     }
#   } else {
#     frac_ddG_ratio_pre <- frac_ddG_obs[-1]/frac_ddG_pred[-1]
#     frac_ddG_ratio_pre[which(is.infinite(frac_ddG_ratio_pre))] <- NaN
#     output_matrix <- data.frame(`mm_type`=mm_type, `pos`=pos,
#                                 `nuc_pair`=miR_pairing,
#                                 `ddG_pred`=ddG_pred[-1], `ddG_obs`=ddG_obs[-1],
#                                 `ddG_dif` = ddG_obs[-1] - ddG_pred[-1],
#                                 `frac_dG_obs`=frac_dG_obs[-1],
#                                 `frac_dG_pred`=frac_dG_pred[-1],
#                                 `frac_dG_dif`=frac_dG_obs[-1] - frac_dG_pred[-1],
#                                 `frac_ddG_obs`=frac_ddG_obs[-1],
#                                 `frac_ddG_pred`=frac_ddG_pred[-1],
#                                 `frac_ddG_dif`=frac_ddG_obs[-1] - frac_ddG_pred[-1])
#   }
#   return(output_matrix)
# }

# MakeAllDeltaDeltaGMismatchDf <- function(
#   mirna, experiment, len, n_constant=3, bulge=FALSE, model_values=FALSE,
#   mm_and_bulge=TRUE, new=TRUE, kd_fc=TRUE, win_average=3, corrected_kds=TRUE,
#   plot_pred=FALSE, plot_emp=FALSE, ratio=FALSE, obs_and_pred=FALSE,
#   frac_dG=FALSE, frac_ddG=FALSE
# ) {
#   ########## Define the limits of starting positions for that miRNA ############
#   if (mirna == "miR-1") {
#     mm_start_min <- 11
#   } else {
#     mm_start_min <- 9
#   }
#   # mm_start_min <- 10
#   mm_start_max <- nchar(kMirnaSeqs[mirna]) - len + 1 
#   ########################### Make the starting matrix #########################
#   d_all <<- do.call(
#     "rbind", lapply(mm_start_min:mm_start_max, function(start_mm) {
#       GetDeltaDeltaGResidual(mirna, experiment, as.numeric(start_mm),
#         start_mm + len - 1, n_constant=n_constant, bulge=bulge,
#         model_values=model_values, new=new, kd_fc=kd_fc,
#         win_average=win_average, corrected_kds=corrected_kds,
#         obs_and_pred=obs_and_pred, frac_dG=frac_dG)
#   }))
#   ##################### Average each mutation X position #######################
#   if (frac_ddG) {
#     if (plot_pred) {
#       df_pos <- aggregate(d_all$frac_ddG_pred, list(d_all$nuc_pair), mean, na.rm=TRUE)
#     } else if (plot_emp) {
#       df_pos <- aggregate(d_all$frac_ddG_obs, list(d_all$nuc_pair), mean, na.rm=TRUE)
#     } else if (ratio) {
#       df_pos <- aggregate(d_all$frac_ddG_obs, list(d_all$nuc_pair), mean, na.rm=TRUE)
#       df_pos_pred <- aggregate(d_all$frac_ddG_pred, list(d_all$nuc_pair), mean, na.rm=TRUE)
#     } else {
#       df_pos <- aggregate(d_all$frac_ddG_dif, list(d_all$nuc_pair), mean, na.rm=TRUE)
#     }
#   } else if (frac_dG) {
#     if (plot_pred) {
#       df_pos <- aggregate(d_all$frac_dG_pred, list(d_all$nuc_pair), mean, na.rm=TRUE)
#     } else if (plot_emp) {
#       df_pos <- aggregate(d_all$frac_dG_obs, list(d_all$nuc_pair), mean, na.rm=TRUE)
#     } else if (ratio) {
#       df_pos <- aggregate(d_all$frac_dG_obs, list(d_all$nuc_pair), mean, na.rm=TRUE)
#       df_pos_pred <- aggregate(d_all$frac_dG_pred, list(d_all$nuc_pair), mean, na.rm=TRUE)
#     } else if (obs_and_pred) {
#       df_pos <- aggregate(d_all[, 4:5], list(d_all$nuc_pair), mean, na.rm=TRUE)
#     } else {
#       df_pos <- aggregate(d_all$frac_dG_dif, list(d_all$nuc_pair), mean, na.rm=TRUE)
#     }
#   } else {
#     if (plot_pred) {
#       df_pos <- aggregate(d_all$ddG_pred, list(d_all$nuc_pair), mean, na.rm=TRUE)
#     } else if (plot_emp) {
#       df_pos <- aggregate(d_all$ddG_obs, list(d_all$nuc_pair), mean, na.rm=TRUE)
#     } else if (ratio) {
#       df_pos <- aggregate(d_all$ddG_obs, list(d_all$nuc_pair), mean, na.rm=TRUE)
#       df_pos_pred <- aggregate(d_all$ddG_pred, list(d_all$nuc_pair), mean, na.rm=TRUE)
#     } else if (obs_and_pred) {
#       df_pos <- aggregate(d_all[, 4:5], list(d_all$nuc_pair), mean, na.rm=TRUE)
#     } else {
#       df_pos <- aggregate(d_all$ddG_dif, list(d_all$nuc_pair), mean, na.rm=TRUE)
#     }
#   }
#   df_pos_global <<- df_pos
#   pairs <- gsub("^p(.*)(m)(.*)$", df_pos$Group.1, replacement="\\2\\3")
#   if (ratio) {
#     df_pairs_obs <<- data.frame(`pairs`=pairs, `dG_obs`=df_pos$x)
#     df_pairs_pred <<- data.frame(`pairs`=pairs, `dG_pred`=df_pos_pred$x)
#     df_mm_obs <<- aggregate(df_pairs_obs$dG_obs, list(df_pairs_obs$pairs), mean, na.rm=TRUE)
#     df_mm_pred <<- aggregate(df_pairs_pred$dG_pred, list(df_pairs_pred$pairs), mean, na.rm=TRUE)
#     df_mm <- data.frame(`Group.1`=df_mm_obs$Group.1, x = df_mm_obs$x/df_mm_pred$x)
#   } else if (obs_and_pred) {
#     df_pairs <- data.frame(`pairs`=pairs, `ddG_obs`=df_pos$ddG_obs,
#                            `ddG_pred`=df_pos$ddG_pred)
#     df_mm <- aggregate(df_pairs[, 2:3], list(df_pairs$pairs), mean, na.rm=TRUE)
#   } else {
#     df_pairs <- data.frame(`pairs`=pairs, `dG_res`=df_pos$x)
#     df_mm <- aggregate(df_pairs$dG_res, list(df_pairs$pairs), mean, na.rm=TRUE)
#   }
#   df_mm_global <<- df_mm
#   df_mm <- df_mm[order(df_mm[, 1]), ]
#   return(df_mm)
# }

# MakeAllMirnaDeltaDeltaGMismatchDf <- function(
#   library_type="programmed", len, n_constant=3, bulge=FALSE, model_values=FALSE,
#   mm_and_bulge=TRUE, new=TRUE, kd_fc=TRUE, win_average=3, corrected_kds=TRUE,
#   frac_dG=FALSE, frac_ddG=FALSE
# ) {
#   ########## Define the limits of starting positions for that miRNA ############
#   if (library_type == "programmed") {
#     experiment_df <- data.frame(
#       mirna=c("let-7a-21nt", "miR-1", "miR-155"),
#       experiment=c("equil_c2_nb", "equil_c_nb", "equil_sc_nb")
#     )
#   }
#   df_mm_all <- do.call("rbind", apply(experiment_df, 1, function(row) {
#     SubfunctionCall(MakeAllDeltaDeltaGMismatchDf, mirna=row[1],
#                     experiment=row[2], obs_and_pred=TRUE, frac_dG=frac_dG,
#                     frac_ddG=frac_ddG)
#   }))
#   df_mm <- aggregate(df_mm_all[, 2:3], list(df_mm_all$Group.1), mean, na.rm=TRUE)
#   return(df_mm)
# }




# MakePositionalMismatchDf <- function(
#   mirna, experiment, len, n_constant=3, bulge=FALSE, model_values=FALSE,
#   mm_and_bulge=TRUE, new=TRUE, offset_lim=c(-4, 16), kd_fc=TRUE, win_average=3,
#   corrected_kds=TRUE, frac_dG=FALSE
# ) {
#   ########## Define the limits of starting positions for that miRNA ############
#   mm_start_min <- 9
#   mm_start_max <- nchar(kMirnaSeqs[mirna]) - len + 1 
#   ########################### Make the starting matrix #########################
#   d_all <<- do.call(
#     "rbind", lapply(mm_start_min:mm_start_max, function(start_mm) {
#       GetDeltaDeltaGResidual(mirna, experiment, as.numeric(start_mm),
#         start_mm + len - 1, n_constant=n_constant, bulge=bulge,
#         offset_lim=offset_lim, model_values=model_values, new=new, kd_fc=kd_fc,
#         win_average=win_average, corrected_kds=corrected_kds, frac_dG=frac_dG,
#         obs_only=TRUE)
#   }))
#   ##################### Average each mutation X position #######################
#   df_pos <- aggregate(d_all$ddG_obs, list(d_all$nuc_pair), mean, na.rm=TRUE)
#   if (bulge) {
#     pairs <- gsub("^p(.*)(m)(..)(.*)$", df_pos$Group.1, replacement="\\1\\3")
#   } else {
#     pairs <- gsub("^p(.*)(m)(.)(.*)$", df_pos$Group.1, replacement="\\1\\3")    
#   }
#   df_pos_m <<- data.frame(`pairs`=pairs, `dG_obs`=df_pos$x)
#   df_pos <- aggregate(df_pos_m$dG_obs, list(df_pos_m$pairs), mean, na.rm=TRUE)
#   colnames(df_pos)[1] <- "pos_nuc"
#   pos <- gsub("^(.*)([ACUG])$", df_pos$pos_nuc, replacement="\\1")
#   if (bulge) {
#     nuc <- gsub("^(.*)([ACUG])(.)$", df_pos$pos_nuc, replacement="\\3")
#   } else {
#     nuc <- gsub("^(.*)([ACUG])$", df_pos$pos_nuc, replacement="\\2")
#   }
#   df_pos <- cbind(df_pos, pos, nuc, rep(mirna, nrow(df_pos)))
#   colnames(df_pos) <- c("pos_nuc", "dG_obs", "pos", "nuc", "mirna")
#   return(df_pos)
# }


# MakeAllLibraryMismatches <- function(
#   len, n_constant=3, bulge=FALSE, model_values=FALSE, mm_and_bulge=TRUE,
#   new=TRUE, offset_lim=c(-4, 16), kd_fc=TRUE, win_average=3, corrected_kds=TRUE, frac_dG=FALSE,
#   rand_data=FALSE
# ) {
#   if (rand_data) {
#     mirna_and_exp_df <- data.frame(
#       mirnas=c("let-7a", "miR-1", "miR-155", "miR-124", "lsy-6",
#                "miR-7-23nt"),
#       experiments=c("equilibrium", "equilibrium", "equilibrium", "equilibrium",
#                     "equilibrium", "equilibrium2_nb")
#     )
#   } else {
#     mirna_and_exp_df <- data.frame(
#       mirnas=c("let-7a-21nt", "miR-1", "miR-155", "let-7a_plus1",
#                "let-7a_minus1", "let-7a_miR-155", "miR-155_let-7a"),
#       experiments=c("equil_c2_nb", "equil_c_nb", "equil_sc_nb", "equil_c_nb",
#                     "equil_c_nb", "equil_c_nb", "equil_c_nb")
#     )
#     mirna_and_exp_df <- mirna_and_exp_df[1:3, ]
#   }
#   df_all <- do.call(
#     "rbind", apply(mirna_and_exp_df, 1, function(row) {
#       SubfunctionCall(MakePositionalMismatchDf, mirna=row[1], experiment=row[2])
#     })
#   )
#   return(df_all)
# }




# GetBulgeAndMismatchValues <- function(
#   mirna, experiment, len, n_constant=3, model_values=FALSE, modelweights=FALSE,
#   mm_and_bulge=TRUE, new=TRUE, kd_fc=TRUE, win_average=3, corrected_kds=TRUE,
#   threep=TRUE
# ) {
#   ########## Define the limits of starting positions for that miRNA ############
#   mm_start_min <- 9
#   mm_start_max <- nchar(kMirnaSeqs[mirna]) - len + 1 
#   mm_global <- c()
#   bu_global <- c()
#   ratio_global <- c()
#   ########################### Make the starting matrix #########################
#   for (start_mm in mm_start_min:mm_start_max) {
#     if (model_values) {
#       bu <- GetThreePModelMismatchCoefficients(mirna, experiment,
#         as.numeric(start_mm), start_mm + len - 1, n_constant=n_constant,
#         new=new, kd_fc=kd_fc, win_average=win_average, bulge=TRUE,
#         corrected_kds=corrected_kds, modelweights=modelweights)$pairing
#       mm <- GetThreePModelMismatchCoefficients(mirna, experiment,
#         as.numeric(start_mm), start_mm + len - 1, n_constant=n_constant,
#         new=new, kd_fc=kd_fc, win_average=win_average,
#         corrected_kds=corrected_kds, modelweights=modelweights)$pairing

#     } else {
#       bu <- GetThreePrimeMmBDKds(mirna, experiment, as.numeric(start_mm),
#         start_mm + len - 1, n_constant=n_constant, new=new, kd_fc=kd_fc,
#         win_average=win_average, best_average=TRUE, bulge=TRUE,
#         corrected_kds=corrected_kds)
#       mm <- GetThreePrimeMmBDKds(mirna, experiment, as.numeric(start_mm),
#         start_mm + len - 1, n_constant=n_constant, new=new, kd_fc=kd_fc,
#         win_average=win_average, best_average=TRUE,
#         corrected_kds=corrected_kds)

#     }
#     mm_fc <- (mm/mm[1])[-1]
#     bu_fc <- (bu/bu[1])[-1]
#     # If looking at the threeprime end of the site, the mm and bu position
#     # indeces should be the same. If looking at the fiveprime end of the site, 
#     # the bulged position should be greater by 1.
#     print("___________________________")
#     if (threep) {
#       print("threep")
#       target_bu <- sprintf("^[ACGT](\\(.*\\.)*%s\\)*$", start_mm + len - 1)
#       target_bu <- sprintf("^[ACGT]%s$", start_mm + len - 1)
#       target_mm <- sprintf("^[ACGT]%s$", start_mm + len - 1)
#     } else {
#       print("fivep")
#       target_bu <- sprintf("^[ACGT](\\()*%s(\\..*\\))*$", start_mm + 1)
#       target_bu <- sprintf("^[ACGT]%s$", start_mm + 1)
#       target_mm <- sprintf("^[ACGT]%s$", start_mm)
#     }
#     # Get the bulge and mismatch names associated with the position of interest.
#     bu_inds <- grep(target_bu, names(bu), perl=TRUE, value=TRUE)
#     mm_inds <- grep(target_mm, names(mm), perl=TRUE, value=TRUE)
#     out_mm <- mm_fc[mm_inds]
#     out_bu <- bu_fc[bu_inds]
#     # This section renames the bulge indeces to align iwth the mismatches for
#     # mismatches determiend at the 5-prime end.
#     bu_names <- gsub("^([ACGT])(.*)$", replacement="\\1", names(out_bu))
#     if (threep) {
#       bu_names <- sprintf("%s%s", bu_names, start_mm + len - 1)
#     } else {
#       bu_names <- sprintf("%s%s", bu_names, start_mm)
#     }
#     names(out_bu) <- bu_names
#     # Determine which mismatches have a corresponding bulge.
#     inds_shared <- intersect(names(out_bu), names(out_mm))
#     # Add these examples to the three global lists that will be the output.
#     mm_global <- c(mm_global, out_mm[inds_shared])
#     bu_global <- c(bu_global, out_bu[inds_shared])
#     ratio_global <- c(ratio_global, bu_global[inds_shared]/mm_global[inds_shared])
#   }
#   print(mm_global)
#   print(bu_global)
#   print(ratio_global)
#   return(list(`mm`=mm_global, `bu`=bu_global, `ratio`=ratio_global))
# }

# GetAllProgrammedBulgeAndMismatchDifferences <- function(
#   len, n_constant=3, model_values=FALSE, modelweights=FALSE,
#   mm_and_bulge=TRUE, new=TRUE, kd_fc=TRUE, win_average=3, corrected_kds=TRUE,
#   threep=TRUE
# ) {
#   ########## Define the limits of starting positions for that miRNA ############
#   Experiments <- c(`let-7a-21nt`="equil_c2_nb",
#                    `miR-1`="equil_c_nb",
#                    `miR-155`="equil_sc_nb")
#   df_global <- data.frame(mirna=c(), type=c(), end=c(), mm=c(), bu=c(),
#                           ratio=c())
#   print(df_global)
#   for (mirna in c("let-7a-21nt", "miR-1", "miR-155")) {
#     bulge_5p <- SubfunctionCall(GetBulgeAndMismatchValues,
#                                 experiment=Experiments[mirna],
#                                 threep=FALSE)
#     df_5p <- data.frame(mirna=rep(mirna, length(bulge_5p$mm)),
#                         type=names(bulge_5p$mm),
#                         end=rep("5p", length(bulge_5p$mm)),
#                         mm=as.numeric(bulge_5p$mm),
#                         bu=as.numeric(bulge_5p$bu),
#                         ratio=as.numeric(bulge_5p$ratio))
#     df_global <- rbind(df_global, df_5p)

#     bulge_3p <- SubfunctionCall(GetBulgeAndMismatchValues,
#                                 experiment=Experiments[mirna],
#                                 threep=TRUE)
#     df_3p <- data.frame(mirna=rep(mirna, length(bulge_3p$mm)),
#                         type=names(bulge_3p$mm),
#                         end=rep("3p", length(bulge_3p$mm)),
#                         mm=as.numeric(bulge_3p$mm),
#                         bu=as.numeric(bulge_3p$bu),
#                         ratio=as.numeric(bulge_3p$ratio))
#     df_global <- rbind(df_global, df_3p)
#     print(df_global)
#   }
#   return(df_global)
# }


# MakeBeckerThreePDataFrame <- function(
#   mirna, mismatch=FALSE, kd_fc=TRUE, remove_multi=TRUE, remove_CCC=FALSE,
#   no_mm=TRUE
# ) {
#   print("In MakeBeckerThreePDataFrame")
#   data <- LoadBeckerEtAlData(mirna, with_names=TRUE, no_mm=no_mm)
#   if (remove_CCC) {
#     data <- data[grep("CCC", rownames(data), invert=TRUE), ]
#   }
#   all_sites_becker <- as.character(data[, 1])
#   # Deal with duplicate rows
#   inds_multiple <- grep(";", all_sites_becker)
#   if (remove_multi) {
#     data <- data[-inds_multiple, ]    
#   } else {
#     for (row_i in inds_multiple) {
#       row_use_i <- data[row_i, , drop=FALSE]
#       sites_all <- unlist(strsplit(data[row_i, 1], split=";"))
#       fold_split <- length(sites_all)
#       row_use_i[1, 3] <- min(10000, row_use_i[1, 3]*fold_split)
#       row_use_i[1, 4] <- min(10000, row_use_i[1, 4]*fold_split)
#       row_use_i[1, 5] <- min(10000, row_use_i[1, 5]*fold_split)
#       row_use_i[1, 1] <- sites_all[1]
#       data[row_i, ] <- row_use_i
#       for (site_name_i in sites_all[-1]) {
#         row_use_i[1, 1] <- site_name_i
#         data <- rbind(data, row_use_i)
#       }
#     }
#   }

#   # Here aggregate all rows with the same site, getting the Geometric Mean of
#   # their Kd values, and the corresponding confidence intervals.
#   data_global_pre <<- data
#   data <- aggregate(data[, 3:5], list(data[, 1]), GeoMean)
#   data_global <<- data
#   data <- data[-which(data[, 1] == "No seed"), ]
#   data_bipartite <- data[grep("\\|", data[, 1], perl=TRUE), ]
#   data_seed <- data[grep("\\|", data[, 1], perl=TRUE, invert=TRUE), ]


#   data_seed_global <<- data_seed
#   s_pos <- as.integer(gsub("^.*\\|(.*)\\|.*$", replacement="\\1", data_bipartite[, 1]))
#   pos5p <- as.integer(gsub("^.*m(.*)\\..*\\|.*\\|.*$", replacement="\\1", data_bipartite[, 1]))
#   pos3p <- as.integer(gsub("^.*\\.(.*)\\|.*\\|.*$", replacement="\\1", data_bipartite[, 1]))
#   seed_site <- gsub("^.*\\|.*\\|(.*)$", replacement="\\1", data_bipartite[, 1])
#   len <- pos3p - pos5p + 1


#   if (kd_fc) {
#     for (seed_site_i in unique(seed_site)) {
#       inds_i <- which(seed_site == seed_site_i)
#       seed_norm <- data_seed[which(data_seed[, 1] == seed_site_i), ]
#       if (nrow(seed_norm) == 0) {
#         data_bipartite[inds_i, 2] <- 10000/c(data_bipartite[inds_i, 2])
#         data_bipartite[inds_i, 3] <- 10000/c(data_bipartite[inds_i, 3])
#         data_bipartite[inds_i, 4] <- 10000/c(data_bipartite[inds_i, 4])
#       } else {
#         data_bipartite[inds_i, 2] <- as.numeric(seed_norm[2])/c(data_bipartite[inds_i, 2])
#         data_bipartite[inds_i, 3] <- as.numeric(seed_norm[3])/c(data_bipartite[inds_i, 3])
#         data_bipartite[inds_i, 4] <- as.numeric(seed_norm[4])/c(data_bipartite[inds_i, 4])
#       }
#     }
#     logkd <- log10(data_bipartite[, 2])
#   } else {
#     logkd <- log10(100/c(data_bipartite[, 2]))
#   }
#   offset <- s_pos - pos5p
#   return(data.frame(`logkd`=logkd, `s_pos`=as.factor(s_pos),
#              `pos5p`=pos5p, `pos3p`=pos3p,
#              `pairing`=sprintf("%s|%s", pos5p, pos3p), `len`=len,
#              `offset`=offset, `mm`=seed_site,
#              `offsetXmm`=sprintf("%s|%s", offset, seed_site),
#              `pos_5pXmm`=sprintf("%s|%s", pos5p, seed_site),
#              `pos_3pXmm`=sprintf("%s|%s", pos3p, seed_site),
#              `lenXmm`=sprintf("%s|%s", pos3p - pos5p + 1, seed_site),
#              `full_name`=data_bipartite[, 1]))
# }

# MakeBeckerPairingMatrix <- function(
#   mirna, offset_, seed_site_, mismatch=FALSE, kd_fc=TRUE, len_lim=c(4, 11),
#   pos_lim=c(9, 23), get_global_data=TRUE, remove_multi=TRUE, remove_CCC=FALSE,
#   no_mm=FALSE
# ) {
#   print("In MakeBeckerPairingMatrix")
#   # Generate the site names.
#   if (get_global_data) {
#     kd_data_frame_old <<- SubfunctionCall(MakeBeckerThreePDataFrame)
#   }
#   kd_data_frame <- subset(kd_data_frame_old, offset == offset_ & mm==seed_site_)
#   # Define the possible 5-prime starting nucleotides and possible 3-prime
#   # starting nucleotides, for overall matrix.
#   if (mirna == "let-7a")  mirna_use <- "let-7a-21nt"
#   else                    mirna_use <- mirna
#   len_mir <- nchar(kMirnaSeqs[mirna_use])
#   nucs_5p <- 9:(len_mir - len_lim[1] + 1)
#   nucs_3p <- (9 + len_lim[1] - 1):len_mir
#   # Define the output matrix.
#   out_matrix <- matrix(NA, nrow=length(nucs_3p), ncol=length(nucs_5p),
#                        dimnames=list(nucs_3p, nucs_5p))
#   if (nrow(kd_data_frame) != 0) {
#     for (i in 1:nrow(kd_data_frame)) {
#       ind_5p <- as.character(kd_data_frame$pos5p[i])
#       ind_3p <- as.character(kd_data_frame$pos3p[i])
#       out_matrix[ind_3p, ind_5p] <- kd_data_frame$logkd[i]
#     }
#   }
#   return(out_matrix)
# }

# CalculateThreePScore <- function(pos_lim, offset) {
#   pos_5p <- pos_lim[1]
#   pos_3p <- pos_lim[2]
#   region_3p_supp <- seq(13, 16)
#   if (pos_5p > 16 | pos_3p < 13) {
#     return(0)
#   } else {
#     positions <- seq(pos_5p, pos_3p)
#     points_3p_supp <- intersect(positions, region_3p_supp)
#     points_3p_ext <- setdiff(positions, region_3p_supp)
#     points_pairing <- length(points_3p_supp) + 0.5*length(points_3p_ext)
#     offset_diff <- 0.5*(max(abs(offset) - 2, 0))
#     return(points_pairing - offset_diff)
#   }
# }

# GetThreePrimeScoreMatrix <- function(mirna, just_complement=TRUE, wobble=TRUE) {
#   mirna_seq <- kMirnaSeqs[mirna]
#   # Assign the 5 prime and 3 prime nucleotides within the matrix
#   nucs_5p <- 9:(nchar(mirna_seq) - 4 + 1)
#   nucs_3p <- nucs_5p + 3
#   out_matrix <- matrix(NA, nrow=length(nucs_3p), ncol=length(nucs_5p),
#                        dimnames=list(nucs_3p, nucs_5p))
#   # Iterate over the lengths and starting positions, determine the stopping
#   # nucleotide, and calculate the deltaG.
#   for (len in 4:11) {
#     for (start in 9:(nchar(mirna_seq) - len + 1)) {
#       stop <- start + len - 1
#       score <- CalculateThreePScore(c(start, stop), 0)
#       out_matrix[as.character(stop), as.character(start)] <- score
#     }
#   }
#   return(out_matrix)
# }


# MakeLoopNucleotideEnrichmentMatrix <- function(
#   pos, len, offset, conditions="all", n_constant=3, sitelist="progthrp",
#   dinuc=FALSE, dinuc_sim=FALSE, trinuc=FALSE, trinuc_sim=FALSE, AUbin=FALSE,
#   AUbin_sim=FALSE
# ) {
#   # PART 1: Get all of the loop sequences from the input files.
#   # Pre-allocate the loop vectors.
#   loops_I <- c()
#   loops_A <- c()
#   # Define conditions used.
#   if (conditions == "all") {
#     conditions <- c("0.4", "12.6", "4", "12.6", "40")
#   } else {
#     conditions <- c(conditions)
#   }
#   # This list allows the function to pool over multiple 3-prime pairing lengths.
#   for (len_i in len) {
#     # Make the file extension that includes the position of 3-prime pairing, the
#     # positio of 3-prime pairing, and the number of offset nucleotides.
#     ext <- sprintf("_%s_%s_pos%s_len%s_off%s", n_constant, sitelist, pos, len_i,
#                    offset)
#     # Get the input loop sequences.
#     path_I <- GetAnalysisPath("let-7a-21nt", "equil_c2_nb", "I",
#                               "loop_sequences", ext=ext)
#     loops_I <- c(loops_I, unlist(read.table(path_I, stringsAsFactors=FALSE)))
#     # Iterates over the bound conditions.
#     for (condition in conditions) {
#       path <- GetAnalysisPath("let-7a-21nt", "equil_c2_nb", condition,
#                               "loop_sequences", ext=ext)
#       loops_A <- c(loops_A, unlist(read.table(path, stringsAsFactors=FALSE)))
#     }
#   }
#   # PART 2: Iterate over each loop sequence to generate the kmer frequency
#   # tables.
#   # Assign the loop length (repeated variable)
#   len_loop <- nchar(loops_I[1])
#   # Initialize count matrices of both the input and bound loop sequences.
#   # Define the kmers being used, as well as the length of each motif when
#   # correcting for the kmer length being used.
#   if (AUbin) {
#     kmers <- c("A/U", "C/G")
#     len_loop_kmers <- len_loop
#   }
#   if (trinuc) {
#     kmers <- paste0(rep(sort(kDinucs), each=length(kDNucs)), sort(kDNucs))
#     len_loop_kmers <- len_loop - 2
#   } else if (dinuc) {
#     kmers <- sort(kDinucs)
#     len_loop_kmers <- len_loop - 1
#   } else {
#     kmers <- sort(kDNucs)
#     len_loop_kmers <- len_loop
#   }
#   probs_I <- matrix(0, nrow=length(kmers), ncol=len_loop_kmers,
#                     dimnames=list(kmers, seq(len_loop_kmers)))
#   probs_A <- probs_I
#   if (AUbin) {
#     AU_I <- length(grep("^(A|T)*$", loops_I, perl=TRUE))/length(loops_I)
#     AU_A <- length(grep("^(A|T)*$", loops_A, perl=TRUE))/length(loops_A)
#     A_I <- length(grep("^A*$", loops_I, perl=TRUE))/length(loops_I)
#     A_A <- length(grep("^A*$", loops_A, perl=TRUE))/length(loops_A)

#     U_I <- length(grep("^(A|T)T*$", loops_I, perl=TRUE))/length(loops_I)
#     U_A <- length(grep("^(A|T)T*$", loops_A, perl=TRUE))/length(loops_A)


#     return(c(AU_A/AU_I, A_A/A_I, U_A/U_I))
#   } else {
#     loop_list_func <- function(loop) {
#       sapply(seq(len_loop_kmers), function(i) {
#         motif <- substr(loop, start=i, stop=i + len_loop - len_loop_kmers)
#       })
#     }

#     # Iterate over the input loops.
#     for (loop in loops_I) {
#       loop_list <- loop_list_func(loop)
#       for (i in seq(len_loop_kmers)) {
#         nuc <- loop_list[i]
#         probs_I[nuc, i] <- probs_I[nuc, i] + 1
#       }
#     }
#     # Iterate over the bound loops.
#     for (loop in loops_A) {
#       loop_list <- loop_list_func(loop)
#       for (i in seq(len_loop_kmers)) {
#         nuc <- loop_list[i]
#         probs_A[nuc, i] <- probs_A[nuc, i] + 1
#       }
#     }
#     # Normalize both count matrices and take the ratio.
#     probs_I <- probs_I/sum(probs_I)
#     probs_A <- probs_A/sum(probs_A)
#     # Take the ratio of the two normalized matrices as output.
#     R <- probs_A/probs_I
#     if (!dinuc & dinuc_sim) {
#       R_1 <- R[rep(seq(4), each=4), -ncol(R), drop=FALSE]
#       R_2 <- R[rep(seq(4), times=4), -1, drop=FALSE]
#       print(R_1)
#       print(R_2)
#       R <- R_1*R_2
#       rownames(R) <- c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT",
#                        "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")
#       colnames(R) <- seq(len_loop - 1)
#     } else if (!trinuc & trinuc_sim) {
#       if (dinuc) {
#         R_1 <- R[rep(seq(16), each=4), -ncol(R), drop=FALSE]
#         R_2 <- R[, -1, drop=FALSE]
#         print(R_2)
#         # Normalize by each.
#         for (i in seq(0, 3)) {
#           start <- i*4 + 1
#           stop  <- start + 3
#           R_2_means <- apply(R_2[start:stop, ], 2, GeoMean)
#           R_2[start:stop, ] <- t(t(R_2[start:stop, ])/R_2_means)
#         }
#         R_2 <- R_2[rep(seq(16), times=4), ]
#         R <- R_1*R_2
#       } else {
#         R_1 <- R[rep(seq(4), each=16), c(-(ncol(R) - 1), -ncol(R)), drop=FALSE]
#         R_2 <- R[rep(rep(seq(4), each=4), times=4), c(-1, -ncol(R)), drop=FALSE]
#         R_3 <- R[rep(seq(4), times=16), c(-1, -2), drop=FALSE]
#         R <- R_1*R_2*R_3
#       }
#       rownames(R) <- paste0(rep(sort(kDinucs), each=length(kDNucs)), sort(kDNucs))
#       colnames(R) <- seq(len_loop - 2)
#     }
#     return(R)
#   }
# }

# #
# GetThreePrimeReporterCounts <- function(
#   mirna, experiment, condition, analysis_type="counts", rep=1, add_bulges=FALSE,
#   aggregate=FALSE, all_lin41_UTR=FALSE

# ) {
#   file_path <- SubfunctionCall(GetAnalysisPath, ext=sprintf(",%s_120nt", rep))
#   counts_df <- read.table(file_path, stringsAsFactors=FALSE, header=TRUE, sep="\t")
#   # Remove the unmapped counts.
#   counts_df <- counts_df[-nrow(counts_df), ]
#   # Normalize to counts per million mapped reads.
#   counts_df[, "Counts"] <- counts_df[, "Counts"]/sum(counts_df[, "Counts"])*1e6
#   if (add_bulges) {
#     counts_df$Bulge <- gsub("(A|T)", replacement="N", counts_df$Bulge)
#     counts_df_new <- aggregate(
#       x=counts_df$Counts,
#       by=list(counts_df$Seed, counts_df$Location, counts_df$ThreePrime,
#               counts_df$Bulge, counts_df$Context),
#       FUN=sum
#     )
#     colnames(counts_df_new) <- c("Seed", "Location", "ThreePrime", "Bulge", "Context", "Counts")

#     counts_df <- counts_df_new
#   }
#   if (aggregate) {
#     # Convert the `No_site_#` sites to a single `None`.
#     counts_df[, "Seed"] <- gsub("^No_site_[0-9]{1,2}$", "None", counts_df[, "Seed"])
#     # Add up all of the site-type possibilities.
#     sum_df <- aggregate(
#       x=counts_df$Counts,
#       by=list(counts_df$Seed, counts_df$Location, counts_df$ThreePrime,
#               counts_df$Bulge),
#       FUN=sum
#     )
#     counts_df_global <<- counts_df[which(counts_df$Seed == "8mer" &
#                                       counts_df$Location == "left"), ]
#     colnames(sum_df) <- c("Seed", "Location", "ThreePrime", "Bulge", "Counts")
#     # Optional portion in which all of the lin-41 UTR sequence contexts are
#     # added to the dataframe.
#     if (all_lin41_UTR) {
#       # Get all of the lin-41 3-prime UTR contexts.
#       counts_lin41_df <- counts_df[counts_df$Context == "lin-41_context", ]
#       # Make an analogous dataframe to that of the first summed matrix.
#       # Note: the only purpose this serves in practice is to sum the counts for
#       # the seven distinct no-sites corresponding to the lin-41 sequence
#       # context.
#       sum_lin41_df <- aggregate(
#         x=counts_lin41_df$Counts,
#         by=list(counts_lin41_df$Seed, counts_lin41_df$Location,
#                 counts_lin41_df$ThreePrime, counts_lin41_df$Bulge),
#         FUN=sum
#       )
#       colnames(sum_lin41_df) <- c("Seed", "Location", "ThreePrime",
#                                   "Bulge", "Counts")
#       sum_lin41_df$Seed <- paste0(sum_lin41_df$Seed, "_lin41_UTR")
#       sum_df <- rbind(sum_df, sum_lin41_df)
#     }
#     # # Add the lin-41 3'UTR into the dataframe.
#     # lin_41_row <- counts_df[counts_df$Seed == "lin-41" & counts_df$Context == "lin-41_context", ]
#     # lin_41_row[1] <- "lin-41_UTR"
#     # print(lin_41_row)
#     # lin_41_row <- lin_41_row[c("Seed", "Location",
#     #                            "ThreePrime", "Bulge", "Counts")]
#     # sum_df <- rbind(sum_df, lin_41_row)
#     counts_df <- sum_df
#   }
#   return(counts_df)
# }

# GetThreePrimeReporterFoldChanges <- function(
#   mirna, exp_type, experiment="twist_reporter_assay_3p_2_tp", rep=1,
#   add_bulges=FALSE, aggregate=FALSE, all_lin41_UTR=FALSE
# ) {
#   counts_df <<- SubfunctionCall(GetThreePrimeReporterCounts,
#                                condition=sprintf("duplex_%s", exp_type))
#   counts_bg_df <<- SubfunctionCall(GetThreePrimeReporterCounts,
#                                   condition=sprintf("no_duplex_%s", exp_type))
#   counts_global_bg <<- counts_df_global
#   # counts_bg_df <- counts_bg_df[-which(counts_bg_df$Seed == ""), ]
#   log_fc <- log(counts_df[, "Counts"]/counts_bg_df[, "Counts"], 2)
#   fc_df <- cbind(counts_df, log_fc)
#   # print(all_lin41_UTR)
#   # print("in fold changes")
#   # print(head(fc_df, 10))
#   # print(tail(fc_df, 10))
#   return(fc_df)
# }




# GetAllAverages <- function(
#   mirna="let-7", exp_type="parallel", experiment="twist_reporter_assay_3p_2_tp",
#   add_bulges=FALSE, aggregate=TRUE, dual_site=FALSE, all_lin41_UTR=FALSE
# ) {
#   if (mirna == "let-7") {
#     # Add over all four dataframes. ____________________________________________
#     log2fc_df_1 <- SubfunctionCall(
#       GetThreePrimeReporterFoldChanges, mirna="let-7a", rep=1, aggregate=TRUE
#     )
#     log2fc_df_2 <- SubfunctionCall(
#       GetThreePrimeReporterFoldChanges, mirna="let-7a", rep=2, aggregate=TRUE
#     )
#     log2fc_df_3 <- SubfunctionCall(
#       GetThreePrimeReporterFoldChanges, mirna="let-7a-21nt", rep=1,
#       aggregate=TRUE
#     )
#     log2fc_df_4 <- SubfunctionCall(
#       GetThreePrimeReporterFoldChanges, mirna="let-7a-21nt", rep=2,
#       aggregate=TRUE
#     )
#     log2fc_df <- log2fc_df_1[, 1:4]
#     # Add each of the replicates to the dataframe.
#     log2fc_df$log2fc_1 <- log2fc_df_1$log_fc
#     log2fc_df$log2fc_2 <- log2fc_df_2$log_fc
#     log2fc_df$log2fc_3 <- log2fc_df_3$log_fc
#     log2fc_df$log2fc_4 <- log2fc_df_4$log_fc
#   } else if (mirna == "miR-1") {
#     # Only add the two dataframes (since there is only one length of miR-1).
#     log2fc_df_1 <- SubfunctionCall(
#       GetThreePrimeReporterFoldChanges, rep=1, aggregate=TRUE
#     )
#     log2fc_df_2 <- SubfunctionCall(
#       GetThreePrimeReporterFoldChanges, rep=2, aggregate=TRUE
#     )
#     log2fc_df <- log2fc_df_1[, 1:4]
#     # Add each of the replicates to the dataframe.
#     log2fc_df$log2fc_1 <- log2fc_df_1$log_fc
#     log2fc_df$log2fc_2 <- log2fc_df_2$log_fc
#   }
#   # Change the name of the 'No_site' to 'None'.
#   nosite_inds <- grep("No_site", log2fc_df[, 1])
#   log2fc_df[nosite_inds, 1] <- "None"
#   # Make zero-rowed dataframe for output. ______________________________________
#   if (mirna == "let-7")  width <- 4
#   else                    width <- 2
#   if (!dual_site)  width <- 2*width
#   output_df <- data.frame(site=c(), rep=c(), y=c())
#   if (dual_site)  config_str <- "dual"
#   else            config_str <- "single"
#   base_sites <- c("None", "8mer", "7mer-m8", "7mer-A1", "6mer", "8mer-w6")
#   if (dual_site) {
#     log2fc_df <- log2fc_df[which(grepl("None", log2fc_df$Seed) | log2fc_df$Location == "dual"), ]
#   } else {
#     log2fc_df <- log2fc_df[which(grepl("None", log2fc_df$Seed) | log2fc_df$Location != "dual"), ]
#   }

#   for (base_site in base_sites) {
#     inds <- which(log2fc_df$Seed == base_site & log2fc_df$ThreePrime == "")
#     log2fc_df_use <- log2fc_df[inds, 5:ncol(log2fc_df)]
#     if (base_site == "None")  bg_subtract <- mean(unlist(log2fc_df_use))
#     y <- unlist(log2fc_df_use) - bg_subtract
#     for (i in seq(length(y))) {
#       output_df <- rbind(
#         output_df,
#         data.frame(site=sprintf("%s_%s", base_site, config_str), rep=as.character(i), y=y[i])
#       )
#     }
#   }
#   if (all_lin41_UTR) {
#     for (base_site in base_sites) {
#       inds <- which(log2fc_df$Seed == sprintf("%s_lin41_UTR", base_site) & log2fc_df$ThreePrime == "")
#       log2fc_df_global <<- log2fc_df
#       log2fc_df_use <- log2fc_df[inds, 5:ncol(log2fc_df)]
#       if (base_site == "None")  bg_subtract <- mean(unlist(log2fc_df_use))
#       y <- unlist(log2fc_df_use) - bg_subtract
#       for (i in seq(length(y))) {
#         output_df <- rbind(
#           output_df,
#           data.frame(site=sprintf("%s_lin41_UTR_%s", base_site, config_str), rep=as.character(i), y=y[i])
#         )
#       }
#     }
#   }
#   thrp_sites <- c("4mer-m13.16", "9mer-m11.19", "9mer-m13.21")
#   bulges <- c("none", "A", "T", "AAAA", "TTTT")
#   for (thrp_site in thrp_sites) {
#     for (bulge in bulges) {
#       if (bulge == "none")  bulge_use <- ""
#       else                  bulge_use <- bulge
#       inds <- which((log2fc_df[, "Seed"] == "8mer-w6") &
#                     (log2fc_df[, "ThreePrime"] == thrp_site) &
#                     (log2fc_df[, "Bulge"] == bulge_use))
#       log2fc_df_use <- log2fc_df[inds, 5:ncol(log2fc_df)]
#       y <- unlist(log2fc_df_use) - bg_subtract
#       site_label <- sprintf("%s-%s_%s", thrp_site, bulge, config_str)
#       for (i in seq(length(y))) {
#         output_df <- rbind(
#           output_df,
#           data.frame(site=site_label, rep=as.character(i), y=y[i])
#         )
#       }
#     }
#   }
#   if (all_lin41_UTR) {
#     for (thrp_site in thrp_sites) {
#       for (bulge in bulges) {
#         if (bulge == "none")  bulge_use <- ""
#         else                  bulge_use <- bulge
#         inds <- which((log2fc_df[, "Seed"] == "8mer-w6_lin41_UTR") &
#                       (log2fc_df[, "ThreePrime"] == thrp_site) &
#                       (log2fc_df[, "Bulge"] == bulge_use))
#         log2fc_df_use <- log2fc_df[inds, 5:ncol(log2fc_df)]
#         y <- unlist(log2fc_df_use) - bg_subtract
#         site_label <- sprintf("%s-%s_lin41_UTR_%s", thrp_site, bulge, config_str)
#         for (i in seq(length(y))) {
#           output_df <- rbind(
#             output_df,
#             data.frame(site=site_label, rep=as.character(i), y=y[i])
#           )
#         }
#       }
#     }
#   }
#   if (dual_site) {
#     inds <- which(log2fc_df[, "Seed"] == "lin-41")
#     log2fc_df_use <- log2fc_df[inds, 5:ncol(log2fc_df)]
#     y <- unlist(log2fc_df_use) - bg_subtract
#     site_label <- "lin-41"
#     for (i in seq(length(y))) {
#       output_df <- rbind(
#         output_df,
#         data.frame(site=site_label, rep=as.character(i), y=y[i])
#       )
#     }
#     if (all_lin41_UTR) {
#       inds <- which(log2fc_df[, "Seed"] == "lin-41_lin41_UTR")
#       log2fc_df_use <- log2fc_df[inds, 5:ncol(log2fc_df)]
#       y <- unlist(log2fc_df_use) - bg_subtract
#       site_label <- "lin-41_lin41_UTR"
#       for (i in seq(length(y))) {
#         output_df <- rbind(
#           output_df,
#           data.frame(site=site_label, rep=as.character(i), y=y[i])
#         )
#       }
#     }
#   }
#   rownames(output_df) <- 1:nrow(output_df)
#   return(output_df)
# }

# PrintTukeyPValues <- function(all_lin41_UTR=FALSE) {
#   # Get the fold-change values for the single-site counts.
#   single_site_df <- SubfunctionCall(GetAllAverages)
#   # Get the fold-change values for the dual-site counts.
#   dual_site_df <- SubfunctionCall(GetAllAverages, dual_site=TRUE)
#   # Combine the two datasets, removing the redundant "None" sites from the
#   # dual-site dataframe.
#   if (all_lin41_UTR) {
#     all_sites_df <- rbind(single_site_df, dual_site_df[-25:-28, ][-1:-4, ])
#   } else {
#     all_sites_df <- rbind(single_site_df, dual_site_df[-1:-4, ])  
#   }
#   # Generate the Tukey comparison-of-means table, using the TukeyHSD function
#   # on an area of effect linear model.
#   tukey <- TukeyHSD(aov(y ~ site, data=all_sites_df))$site
#   # Define a lower bound for the p-value cut offs, and replace all values below
#   # it with this value.
#   inds_swap <- which(tukey[, 4] < 1e-4)
#   tukey[inds_swap, 4] <- 1e-4

#   print(tukey[c("9mer-m11.19-AAAA_single-9mer-m11.19-A_single",
#                 "9mer-m11.19-AAAA_single-9mer-m11.19-T_single"), , drop=FALSE])
#   print(tukey[c("9mer-m11.19-A_single-9mer-m11.19-none_single",
#                 "9mer-m11.19-T_single-9mer-m11.19-none_single"), , drop=FALSE])
#   print(tukey[c("9mer-m11.19-TTTT_single-9mer-m11.19-none_single",
#                 "9mer-m11.19-TTTT_dual-9mer-m11.19-none_dual"), , drop=FALSE])
#   print(tukey[c("9mer-m13.21-TTTT_single-9mer-m13.21-AAAA_single",
#                 "4mer-m13.16-T_single-4mer-m13.16-A_single"), , drop=FALSE])
#   print(tukey[c("9mer-m13.21-TTTT_dual-9mer-m13.21-AAAA_dual",
#                 "4mer-m13.16-T_dual-4mer-m13.16-A_dual"), , drop=FALSE])
#   print(tukey[c("lin-41-9mer-m11.19-A_dual",
#                 "lin-41-9mer-m11.19-T_dual",
#                 "lin-41-9mer-m11.19-AAAA_dual"), , drop=FALSE])
#   # print(tukey[c("lin-41_UTR-9mer-m11.19-A_dual",
#   #               "lin-41_UTR-9mer-m11.19-T_dual",
#   #               "lin-41_UTR-9mer-m11.19-AAAA_dual"), , drop=FALSE])
#   if (all_lin41_UTR) {
#     seed_sites <- c("8mer", "7mer-m8", "7mer-A1", "6mer", "8mer-w6")
#     sites_3p <- c("4mer-m13.16", "9mer-m11.19", "9mer-m13.21")
#     bulges <- c("none", "A", "T", "AAAA", "TTTT")
#     all_3p_sites <- c(t(outer(sites_3p, bulges, paste, sep="-")))
#     all_sites <- c(seed_sites, all_3p_sites)
#     single_strings <- sapply(all_sites, function(site) {
#       sprintf("%s_lin41_UTR_single-%s_single", site, site)
#     })
#     dual_strings <- sapply(all_sites, function(site) {
#       sprintf("%s_lin41_UTR_dual-%s_dual", site, site)
#     })
#     all_strings <- c(single_strings, dual_strings, "lin-41_lin41_UTR-lin-41")
#     print(tukey[single_strings, , drop=FALSE])

#     print(tukey[dual_strings, , drop=FALSE])
#     print(tukey["lin-41_lin41_UTR-lin-41", , drop=FALSE])

#   }
# }



# CheckLin41EffectSizes <- function(experiment="twist_reporter_assay_3p_2_tp") {
#   # Step 1: Get the aggregate matrix.
#   exp_type <- "parallel"
#   log2fc_df_1 <- SubfunctionCall(
#     GetThreePrimeReporterFoldChanges, mirna="let-7a", rep=1, aggregate=TRUE
#   )
#   log2fc_df_2 <- SubfunctionCall(
#     GetThreePrimeReporterFoldChanges, mirna="let-7a", rep=2, aggregate=TRUE
#   )
#   log2fc_df_3 <- SubfunctionCall(
#     GetThreePrimeReporterFoldChanges, mirna="let-7a-21nt", rep=1, aggregate=TRUE
#   )
#   log2fc_df_4 <- SubfunctionCall(
#     GetThreePrimeReporterFoldChanges, mirna="let-7a-21nt", rep=2, aggregate=TRUE
#   )

#   log2fc_df <- log2fc_df_1[, 1:4]
#   # Add each of the replicates to the dataframe.
#   log2fc_df$log2fc_1 <- log2fc_df_1$log_fc
#   log2fc_df$log2fc_2 <- log2fc_df_2$log_fc
#   log2fc_df$log2fc_3 <- log2fc_df_3$log_fc
#   log2fc_df$log2fc_4 <- log2fc_df_4$log_fc
#   log2fc_df_global <<- log2fc_df
#   nosite_inds <- grep("No_site", log2fc_df[, 1])
#   log2fc_df[nosite_inds, 1] <- "None"
#   log2fc_df_global <<- log2fc_df

#   # Step 2: Get the matrix that only includes the lin-41 context.
#   log2fc_l41_df_1 <- SubfunctionCall(
#     GetThreePrimeReporterFoldChanges, mirna="let-7a", rep=1
#   )
#   log2fc_l41_df_2 <- SubfunctionCall(
#     GetThreePrimeReporterFoldChanges, mirna="let-7a", rep=2
#   )
#   log2fc_l41_df_3 <- SubfunctionCall(
#     GetThreePrimeReporterFoldChanges, mirna="let-7a-21nt", rep=1
#   )
#   log2fc_l41_df_4 <- SubfunctionCall(
#     GetThreePrimeReporterFoldChanges, mirna="let-7a-21nt", rep=2
#   )

#   log2fc_l41_df <- log2fc_l41_df_1[, 1:5]
#   # Add each of the replicates to the dataframe.
#   log2fc_l41_df$log2fc_1 <- log2fc_l41_df_1$log_fc
#   log2fc_l41_df$log2fc_2 <- log2fc_l41_df_2$log_fc
#   log2fc_l41_df$log2fc_3 <- log2fc_l41_df_3$log_fc
#   log2fc_l41_df$log2fc_4 <- log2fc_l41_df_4$log_fc
#   nosite_inds <- grep("No_site", log2fc_l41_df[, 1])
#   log2fc_l41_df[nosite_inds, 1] <- "None"

#   log2fc_all_df <- log2fc_l41_df
#   log2fc_l41_df <- log2fc_l41_df[which(log2fc_l41_df$Context == "lin-41_context"), -5]


#   print(head(log2fc_df))
#   print(head(log2fc_l41_df))
#   print(head(log2fc_all_df))



#   base_sites <- c("None", "8mer", "7mer-m8", "7mer-A1", "6mer", "8mer-w6")
#   for (base_site in base_sites) {
#     if (base_site != "None") {
#       location_bool <- log2fc_df$Location == "dual"
#       location_l41_bool <- log2fc_l41_df$Location == "dual"
#       location_all_bool <- log2fc_all_df$Location == "dual"
#     } else {
#       location_bool <- log2fc_df$Location != "dual"
#       location_l41_bool <- log2fc_l41_df$Location != "dual"
#       location_all_bool <- log2fc_all_df$Location != "dual"
#     }
#     inds <- which((log2fc_df$Seed == base_site) &
#                   (log2fc_df$ThreePrime == "") &
#                   location_bool)
#     inds_context <- which((log2fc_l41_df$Seed == base_site) &
#                   (log2fc_l41_df$ThreePrime == "") &
#                   location_l41_bool)
#     inds_all <- which((log2fc_all_df$Seed == base_site) &
#                   (log2fc_all_df$ThreePrime == "") &
#                   location_all_bool)
#     log2fc_df_use <- log2fc_df[inds, 5:ncol(log2fc_df)]
#     log2fc_l41_df_use <-  log2fc_l41_df[inds_context, 5:ncol(log2fc_l41_df)]
#     log2fc_all_use <- log2fc_all_df[inds_all, 4:ncol(log2fc_all_df)]
#     if (base_site == "None") {
#       bg_subtract <- mean(unlist(log2fc_df_use))
#       bg_subtract_context <- mean(unlist(log2fc_l41_df_use))
#       bg_subtract_context <- bg_subtract
#     }
#     print(sprintf("bg_subtract: %s", bg_subtract))
#     print(sprintf("bg_subtract_context: %s", bg_subtract_context))
#     y <- unlist(log2fc_df_use) - bg_subtract
#     y_context <- unlist(log2fc_l41_df_use) - bg_subtract_context
#     print(base_site)
#     # print(mean(y))
#     # print(mean(y_context))
#     print("division:")
#     print(mean(y_context)/mean(y))
#     print("subtraction:")
#     print(mean(y_context) - mean(y))
#     if (base_site != "None") {
#       # print("here")
#       # print(log2fc_all_use)
#       # print(log2fc_all_use[which.min(rowMeans(log2fc_all_use[, 3:6])), 2])
#       order_use <- order(rowMeans(log2fc_all_use[, 3:6]))
#       # print(order_use)
#       print(log2fc_all_use[order_use, 2][1:5])
#     }
#   }
#   thrp_sites <- c("4mer-m13.16", "9mer-m11.19", "9mer-m13.21")
#   bulges <- c("none", "A", "T", "AAAA", "TTTT")

#   for (thrp_site in thrp_sites) {
#     for (bulge in bulges) {
#       if (bulge == "none") {
#         bulge_use <- ""
#       } else {
#         bulge_use <- bulge
#       }
#       inds <- which((log2fc_df[, "Seed"] == "8mer-w6") &
#                     (log2fc_df[, "ThreePrime"] == thrp_site) &
#                     location_bool &
#                     (log2fc_df[, "Bulge"] == bulge_use))
#       inds_context <- which((log2fc_l41_df[, "Seed"] == "8mer-w6") &
#                     (log2fc_l41_df[, "ThreePrime"] == thrp_site) &
#                     location_l41_bool &
#                     (log2fc_l41_df[, "Bulge"] == bulge_use))
#       inds_all <- which((log2fc_all_df[, "Seed"] == "8mer-w6") &
#                   (log2fc_all_df[, "ThreePrime"] == thrp_site) &
#                   location_all_bool &
#                   (log2fc_all_df[, "Bulge"] == bulge_use))
#       log2fc_df_use <- log2fc_df[inds, 5:ncol(log2fc_df)]
#       log2fc_l41_df_use <-  log2fc_l41_df[inds_context, 5:ncol(log2fc_l41_df)]
#       log2fc_all_use <- log2fc_all_df[inds_all, 4:ncol(log2fc_all_df)]

#       print(sprintf("bg_subtract: %s", bg_subtract))
#       print(sprintf("bg_subtract_context: %s", bg_subtract_context))
#       y <- unlist(log2fc_df_use) - bg_subtract
#       y_context <- unlist(log2fc_l41_df_use) - bg_subtract_context
#       print(thrp_site)
#       print(bulge)
#       print(mean(y))
#       print(mean(y_context))
#       print("division:")
#       print(mean(y_context)/mean(y))
#       print("subtraction:")
#       print(mean(y_context) - mean(y))
#       print(log2fc_all_use)
#       order_use <- order(rowMeans(log2fc_all_use[, 3:6]))
#       # print(order_use)
#       print(log2fc_all_use[order_use, 2][1:5])


#     }
#   }
#   inds <- which(log2fc_df[, "Seed"] == "lin-41")
#   inds_context <- which(log2fc_l41_df[, "Seed"] == "lin-41")
#   log2fc_df_use <- log2fc_df[inds, 5:ncol(log2fc_df)]
#   log2fc_l41_df_use <-  log2fc_l41_df[inds_context, 5:ncol(log2fc_l41_df)]

#   print(sprintf("bg_subtract: %s", bg_subtract))
#   print(sprintf("bg_subtract_context: %s", bg_subtract_context))
#   y <- unlist(log2fc_df_use) - bg_subtract
#   y_context <- unlist(log2fc_l41_df_use) - bg_subtract_context
#   print("lin-41 site:")
#   print(mean(y))
#   print(mean(y_context))
#   print("division:")
#   print(mean(y_context)/mean(y))
#   print("subtraction:")
#   print(mean(y_context) - mean(y))



# }


