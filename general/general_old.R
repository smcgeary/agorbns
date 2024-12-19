# rm(list=ls())
# gc()
options(warn=0)
kSolexaDir <- "/lab/solexa_bartel/mcgeary/AgoRBNS/"
kHomeDir   <- "/lab/bartel1_ata/mcgeary/computation/AgoRBNS/"
# Load relevant libraries:
library(beeswarm)
# library(gplots)
# library(wrswoR)

library(deSolve)
library(data.table)
library(gdata)
# library(corrplot)
# library(deming)
library(numDeriv)
# library(xlsx)
# library(wCorr)
# library(weights)
# library(Matching)
# library(stringr)
library(reticulate)
# library(viridis)

library(rje)
# library(nloptr)
library(inlmisc)
# library(matlib)
library(MASS)
library(scam)
library(shape)

source("general/Lists.R")
SubfunctionCall <- function(f, ...) {
  # Check for the supplied arguments in "...":
  f.args.supplied <- list(...)
  # Remove these arguments from the list of all the arguments of function "f":
  f.arg.names <- setdiff(names(formals(f)), names(f.args.supplied))
  # Check for these arguments in the enviromnent:
  f.args <- lapply(f.arg.names, function(arg) {
    arg.temp <- try(dynGet(arg), silent=TRUE)
    if (class(arg.temp) != "try-error" & class(arg.temp) != "error") arg.temp
    else "remove"
  })
  names(f.args) <- f.arg.names
  # Combine the found arguments with the supplied and perform the function:
  f.args.list <- as.list(c(f.args[f.args != "remove"], f.args.supplied))
  do.call(f, f.args.list)
}

# General metadata##############################################################
kAgoStock <- data.frame(fread(sprintf("%sSolveForKds/k_c_stockago.txt",
                                      kHomeDir), sep="\t"), row.names=1)

kAgoRad <- data.frame(fread(sprintf("%sSolveForKds/equilibrium_radioactivity_quants.txt",
                                    kHomeDir), sep="\t"), row.names=1)
colnames(kAgoRad) <- gsub("^(.*)\\.(.*)$", colnames(kAgoRad), replace="\\1-\\2",
                          perl=TRUE)
colnames(kAgoRad)[ncol(kAgoRad)] <- "miR-7-23nt"



kAgoStockFunc <- function(mirna, logtrans=FALSE) {
  x_vals <- as.numeric(rownames(kAgoRad))
  y_vals <- kAgoRad[, mirna]*100
  # dev.new(xpos=20, ypos=20, height=5, width=5)
  if (logtrans) logstr <- "xy"
  else          logstr <- ""
  # plot(x_vals, y_vals, log=logstr)
  FitAgoData <- function(pars, logtrans=FALSE) {
    pars <- exp(pars)
    A <- pars[1]
    b <- pars[2]
    y_model <- A*x_vals + b
    # lines(x_vals, y_model)
    if (logtrans) {
      y_vals <- log(y_vals)
      y_model <- log(y_model)
    }
    sum((y_vals - y_model)^2)
  }
  out <- exp(optim(c(1, 1), FitAgoData, logtrans=logtrans)$par)
}

kAgoStockNew <- sapply(kMirnas[1:6], kAgoStockFunc)


# new_col <- rep(NA, nrow(kAgoStock))
# names(new_col) <- rownames(kAgoStock)
# kAgoStock <- cbind(cbind(kAgoStock, new_col), new_col)
# kAgoStockOld <- kAgoStock[kMirnas, 1]
# names(kAgoStockOld) <- kMirnas
# colnames(kAgoStock)[5:6] <- c("equilibrium_met_tp", "equilibrium_2_tp")
# kAgoStock["miR-1", "equilibrium_tp"] <- 5
# kAgoStock["miR-1", "equilibrium_met_tp"] <- 5
# kAgoStock["miR-124", "equilibrium_tp"] <- 2
# kAgoStock["miR-124", "equilibrium_2_tp"] <- 2


kAgoStock[kMirnas[1:6],"equilibrium"] <- kAgoStockNew[1, ]

#BASIC MATH FUNCTIONS###########################################################
Norm <- function(vector) {
  # Returns the normalized entries of the vector.
  vector/sum(vector, na.rm=TRUE)
}

MatNorm <- function(matrix) {
  # Returns the normalized entries of the vector.
  out <- apply(matrix, 2, Norm)
  rownames(out) <- rownames(matrix)
  colnames(out) <- colnames(matrix)
  out
}

GetEnrichment <- function(x_A, x_I) {
  Norm(x_A)/Norm(x_I)
}

GeoMean <- function(vector) {
  exp(mean(log(vector), na.rm=TRUE))
}

Logistic <- function(vector, max) {
  return(max/(1 + exp(-vector)))
}

# TODO finish this.
ProbOverlaps <- function(len_site, len_competitor, num_overlaps, len_overlaps,
                         wobble=FALSE) {
  if (wobble) p <- 3/8
  else        p <- 1/4
  prob_1 <- 1 - pbinom(0, size=len_site - len_overlaps + 1, prob=p^len_overlaps)
  print(prob_1)
  print(sprintf("probability of one site register matching one window: %s", prob_1))
  prob_2 <- 1 - pbinom(num_overlaps - 1, size=len_competitor - len_overlaps + 1,
                   prob=prob_1)
  print(sprintf("probability of having at least %s pairing: %s", num_overlaps, prob_2))
  # print(prob_2)
  # print("prob_2 * 4^len_site")
  # print(prob_2*4^len_site)
  # print(4^len_site)

}


#BASIC STRING REPLACEMENT FUNCTIONS#############################################
StrRev <- function(x) {
  # For a single string or a vector of strings, returns each element of the
  # string in reverse order.
  return(sapply(lapply(strsplit(x, NULL), rev), paste, collapse=""))
}

SeqComplement <- function(x, RNA=FALSE) {
  # For a single string or a vector of strings, returns the complementary RNA
  # sequence (not reversed!).
  complement.base <- c(U="A", T="A", G="C", C="G", N="N", V="B", B="V", D="H", H="D")
  if (RNA) complement <- c(complement.base, A="U")
  else complement <- c(complement.base, A="T")
  return(sapply(lapply(strsplit(x, NULL), function(y) complement[y]), paste,
                collapse=""))
}

RevComplement <- function(x, RNA=FALSE) {
  StrRev(SeqComplement(x, RNA=RNA))
}

ConvertTtoU <- function(site) {
  gsub("T", c(site), replace="U")
}

ConvertUtoT <- function(site) {
  gsub("U", c(site), replace="T")
}

ConvertTtoUandMmtoX <- function(site) {
  gsub("mm", c(gsub("T", c(site), replace="U")), replace="x")
}

GetRevMirnaIndex <- function(mirna, ind) {
  return(nchar(kMirnaSeqs[mirna]) - ind + 1)
}

GetSubKmers <- function(sequence, len_sequence) {
  sapply(1:(nchar(sequence) - len_sequence + 1), function(pos) {
    substr(sequence, pos, pos + len_sequence - 1)
  })
}

GetLenComplementaritySeq <- function(str_seq_1, str_seq_2, wobble=FALSE) {
  vec_rc_strings <- c("AT", "TA", "CG", "GC")
  vec_wob_strings <- c("TG", "GT")
  if (wobble) vec_rc_strings <- c(vec_rc_strings, vec_wob_strings)
  vec_seq_1 <- unlist(strsplit(ConvertUtoT(str_seq_1), NULL))
  vec_seq_2 <- rev(unlist(strsplit(ConvertUtoT(str_seq_2), NULL)))
  vec_seq_pairs <- paste(vec_seq_1, vec_seq_2, sep="")
  comp_length <- 0
  comp_start <- 1
  comp_stop <- 0
  comp_wob <- c()
  comp_length_best <- comp_length
  comp_start_best <- comp_start
  comp_stop_best <- comp_stop
  comp_wob_best <- comp_wob
  for (pos in 1:length(vec_seq_pairs)) {
    if (vec_seq_pairs[pos] %in% vec_rc_strings) {
      comp_length <- comp_length + 1
      comp_stop <- pos
      if (vec_seq_pairs[pos] %in% vec_wob_strings) {
        comp_wob <- c(comp_wob, pos)
      }
    } else {
      if (comp_length > comp_length_best) {
        comp_length_best <- comp_length
        comp_start_best <- comp_start
        comp_stop_best <- comp_stop
        comp_wob_best <- comp_wob
      }
      comp_length <- 0
      comp_start <- pos + 1
    }
  }
  if (comp_length > comp_length_best) {
    comp_length_best <- comp_length
    comp_start_best <- comp_start
    comp_stop_best <- comp_stop
    comp_wob_best <- comp_wob
  }
  return(list(`len`=comp_length_best, `start`=comp_start_best,
           `stop`=comp_stop_best, `wob`=comp_wob_best))
}


GetMaxCompetitorOligoEnrichment <- function(mirna, str_site, wobble=FALSE) {
  # Get competitor oligo sequence.
  str_competitor <- SequenceList[["competitor"]][mirna]
  # Make matrix to store each of the greatest complementarity positions when
  # scanning across the matrix.
  mat_comp <- matrix(data=0, nrow=nchar(str_competitor) - nchar(str_site) + 1,
                     ncol=3)
  # Iterate over each possible position of complementarity.
  for (int_start in seq(1, nchar(str_competitor) - nchar(str_site) + 1)) {
    # Get the substring of the competitor for which the complementarity check
    # is to be performed.
    str_comp_region <- substr(str_competitor,
                                  start=int_start,
                                  stop=int_start + nchar(str_site) - 1)
    # Store the output of function in mat_comp in the appropriate row.
    mat_comp[int_start, ] <- unlist(
      GetLenComplementaritySeq(str_site, str_comp_region, wobble=wobble)
    )[1:3]
  }
  # Return the maximum of the first column, which is the length of contiguous
  # complementarity at that starting position within the competitor oligo.
  return(max(mat_comp[, 1]))
}




printq <- function(x) {
  print(noquote(x))
}

GetKmerList <- function(k) {
  if (k == 1) {
    return(sort(kDNucs))
  } else {
    return(c(t(sapply(GetKmerList(k-1), function(kmer) {
      paste0(sort(kDNucs), kmer)
    }))))
  }
}

GetKmersFromString <- function(string, len_k) {
  sapply(1:(nchar(string) - len_k + 1), function(i_start) {
    i_stop <- i_start + len_k - 1
    substring(string, i_start, i_stop)
  })
}

MirnaTargetSequence <- function(mirna, start, stop, RNA=FALSE) {
  # Gives the reverse-compement sequence to the miRNA positions, using R
  # indexing (not pythonic).
  mirna.seq <- substr(kMirnaSeqs[mirna], start, stop)
  return(RevComplement(mirna.seq, RNA=RNA))
}


#I/O functions##################################################################

EnsureDirectory <- function(dir_path) {
  if (!file.exists(dir_path)) {
    dir.create(dir_path)
  }
}

GetAnalysisPath <- function(mirna, experiment, condition, analysis_type, ext="",
                            suffix="txt") {
  sub_dir <- sprintf("%s/%s/%s", mirna, experiment, analysis_type)
  dir_path <- sprintf("%s%s", kSolexaDir, sub_dir)
  EnsureDirectory(dir_path)
  sprintf("%s/%s%s.%s", dir_path, condition, ext, suffix)
}

CountLines <- function(path) {
  command <- sprintf("wc -l %s", path)
  as.numeric(strsplit(system(command, intern=TRUE), split=" ")[[1]][1])
}

GrepShell <- function(key, path, str_left="", str_right="", start=1,
                      stop="length($0)") {
  command_awk <- sprintf("awk '{print \"%s\"substr($0, %s, %s)\"%s\"}' %s",
                         str_left, start, stop, str_right, path)
  command_head <- sprintf("%s | grep %s | head", command_awk, key)
  command_wcl <- sprintf("%s | grep %s | wc -l", command_awk, key)
  system(command_head, intern=TRUE)
  as.numeric(strsplit(system(command_wcl, intern=TRUE), split=" ")[[1]][1])
}

# LCSubString <- function(seq1, seq2) {
#   # python.load("general/RBNS_methods.py")
#   source_python("general/general.py")

#   main <- py_run_string(sprintf("overlap = LCSubString('%s', '%s')", seq1, seq2))
#   main$overlap
# }

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

GetReporterAmbiguousSites <- function(mirna) {
  all_sites <- rownames(ReporterCounts(mirna, "twist_reporter_assay", "duplex",
                                       rep=1))
  # Select all those sites that the correct mirna string.
  sites <- grep(sprintf("%s_", mirna), all_sites, value=TRUE)
  # Select those sites that don't have "mer" in the name (and are therefore
  # constructed from only A T C and G).
  ambiguous_sites <- grep("mer", sites, invert=TRUE, value=TRUE)
  # Remove the miRNA prefix string before output.
  ambiguous_sites <- gsub(sprintf("%s_", mirna), ambiguous_sites, replace="")
  return(ambiguous_sites)
}


GetSiteList <- function(mirna, sitelist) {
  base_path <- file.path("AssignSiteTypes", "site_categories")
  as.vector(
    unlist(
      read.table(
        file.path(
          base_path, sitelist, sprintf("sites.%s_%s.txt", mirna, sitelist)
        )
      )
    )
  )
}

PrintFinalSiteList <- function(mirna, UandX=TRUE) {
  message(mirna)
  message("")
  # Get the site lists from the papercutoff1 and papercutoff2 k-mer analyses.
  sites <- SubfunctionCall(GetSiteList, sitelist="papercutoff")
  if (mirna %in% c("miR-1", "miR-124", "miR-7")) {
    sites2 <- SubfunctionCall(GetSiteList, sitelist="papercutoff2")
    sites <- intersect(sites, sites2)
  } else {
    rc_competitor_seq <- RevComplement(kCompetitorOligoSeqs[mirna])
    inds_keep <- which(
      sapply(sites, function(site) {
        max(nchar(LargestCommonSubstring(rc_competitor_seq, site)))
      }) < 6
    )
    sites <- sites[inds_keep]
  }
  #   for (site in sites) {
  #     print(nchar(LargestCommonSubstring(rc_competitor_seq, site)))
  #   }
  # Print each of the sites found in both site lists, swapping the "mm" naming
  # convention for "x", and swapping "T" nucleotides for "U" nucleotides.
  if (UandX) {
    sites <- gsub("T", gsub("mm", sites, replace="x"), replace="U")
  }
  for (site in sites) {
    message(site)
  }
  message("")
}

GetSiteSeq <- function(mirna, site) {
  # python.load("general/RBNS_methods.py")
  source_python("general/RBNS_methods.py")
  main <- py_run_string(sprintf("sequence = Mirna('%s')['%s']", mirna, site))
  main$sequence
}




MinimalExcludeStringList <- function(list) {
  source_python("general/general.py")

  python_list = sprintf("['%s']", paste0(list, collapse="', '"))
  print(python_list)

  main <- py_run_string(sprintf("result = MinimalExcludeStringList(%s)",
                                python_list))
  main$result
}

GetMirnaCountData <- function(mirna, condition, unique=FALSE, no_adapter=FALSE,
                              no_marker=FALSE) {
  ext <- ""
  if (unique)  ext <- sprintf("%s_unique", ext)
  if (no_adapter) ext <- sprintf("%s_noadapter", ext)
  if (no_marker) ext <- sprintf("%s_nomarker", ext)
  path <- GetAnalysisPath(mirna, "AGO_purity", condition,
                           "AGO_pur_counts", ext=ext)
  file_counts <- read.table(path, sep="\t", header=TRUE, row.names=1)
  file_counts
}

SitesXCountsVector <- function(mirna, condition, experiment="equilibrium", n_constant=5,
                         sitelist="resubmissionfinal", alt_mirna=NULL, mirna.start=NULL,
                         multisite=FALSE, split16=NULL, uniq=FALSE, buffer=FALSE,
                         comp=FALSE, c_k=6, m_k=6) {
  if (multisite) analysis_type <- "multisite_counts"
  else           analysis_type <- "site_counts"
  ext <- sprintf("_%s_%s", n_constant, sitelist)
  if (!is.null(alt_mirna)) {
    ext <- sprintf("%s_alt-%s", ext, alt_mirna)
  }
  if (sitelist %in% kKmerSiteLists) {
    ext <- sprintf("%s_%s-%s", ext, mirna.start, as.integer(mirna.start) + 3)
  }
  if (sitelist == "16mers") {
    ext <- sprintf("%s_%s", ext, split16)
  }
  if (uniq) {
    ext <- sprintf("%s_uniq", ext)
  }
  if (comp) {
    m_k <- min(m_k, c_k)
    ext <- sprintf("%s_comp_c%sm%s", ext, c_k, m_k)
  }
  if (buffer) {
    ext <- sprintf("%s_buffer3p", ext)
  }
  file <- SubfunctionCall(GetAnalysisPath, ext=ext)
  print(file)
  
  sXc <- data.frame(fread(file, sep="\t", header=FALSE, stringsAsFactors=FALSE),
                    row.names=1)
  colnames(sXc) <- c(condition)
  # Remove "Seq" column:
  # Check the two conditions, that the site exists in the input library:
  return(sXc)
}



SitesXCounts <- function(mirna, experiment="equilibrium", n_constant=5,
                         sitelist="resubmissionfinal", mirna.start=NULL,
                         multisite=FALSE, split16=NULL, uniq=FALSE,
                         start_mm=FALSE, stop_mm=FALSE, buffer=FALSE,
                         comp=FALSE, sc_k=6, m_k=6, new=FALSE, new2=FALSE,
                         include_zeros=FALSE, zeros_override=FALSE) {
  prog_sitelists <- c("programmed", "programmed_collapsed",
                      "programmed_suppcomp")
  ext <- sprintf("_%s", n_constant)
  print(sitelist)
  # if (sitelist %in% prog_sitelists & experiment %in% kExpsThreeP) {
  #   analysis_type <- "programmed_site_count_tables"
  #   if (class(start_mm) == "numeric" & class(stop_mm) == "numeric") {
  #     ext <- sprintf("%s_m%s.%smmsd", ext, start_mm, stop_mm)            
  #   }
  #   if (sitelist == "programmed_collapsed") {
  #     ext <- sprintf("%s_collapsed", ext)
  #   } else if (sitelist == "programmed_suppcomp") {
  #     ext <- sprintf("%s_suppcomp", ext)
  #   }
  if (multisite) {
    analysis_type <- "multisite_count_tables"
    ext <- sprintf("%s_%s", ext, sitelist)
  } else {
    analysis_type <- "site_count_tables"
    ext <- sprintf("%s_%s", ext, sitelist)
  }
  if (class(start_mm) == "numeric" & class(stop_mm) == "numeric") {
    ext <- sprintf("%s_m%s.%smmsd", ext, start_mm, stop_mm)            
  }
  condition <- "all_sites"
  if (sitelist %in% kKmerSiteLists) {
    ext <- sprintf("%s_%s-%s", ext, mirna.start, as.integer(mirna.start) + 3)
  }
  if (sitelist == "16mers") {
    ext <- sprintf("%s_%s", ext, split16)
  }
  if (uniq) {
    ext <- sprintf("%s_uniq", ext)
  }
  if (comp) {
    m_k <- min(m_k, c_k)
    ext <- sprintf("%s_comp_c%sm%s", ext, c_k, m_k)
  }
  if (buffer) {
    ext <- sprintf("%s_buffer3p", ext)
  }
  if (new) {
    ext <- sprintf("%s_new", ext)
  }
  if (new2) {
   ext <- sprintf("%s_new2", ext)
  }
  if (include_zeros) {
   ext <- sprintf("%s_zeros", ext)
  }
  file <- SubfunctionCall(GetAnalysisPath, ext=ext)
  print(file)
  sXc <- data.frame(fread(file, sep="\t", header=TRUE, stringsAsFactors=FALSE),
                    row.names=1)
  # Remove "Seq" column, if it exists:
  sXc <- sXc[, grep("Seq", colnames(sXc), inv=TRUE)]
  # Remove the Leading `A` character from the AGO concentration samples.
  colnames(sXc) <- gsub("^A", colnames(sXc), replace="", perl=TRUE)
  # Replaces ".1" with ",1" at the end (this is really bad), because replicates
  # are named with ",1", ",2", ",3", etc.
  colnames(sXc) <- gsub("\\.([1234])$", colnames(sXc), replace=",\\1",
                        perl=TRUE)
  # This is a bandaid for the prior regex conversion because 0.4 is a standard
  # dilution used in the RBNS experiments.
  colnames(sXc)[which(colnames(sXc) == "0,4")] <- "0.4"
  # Check the two conditions, that the site exists in the input library:
  if (!(sitelist %in% c("12mers", "16mers")) & (!include_zeros) & (!zeros_override)) {
    sXc <- sXc[rowSums(sXc[, 3:ncol(sXc)], na.rm=TRUE) > 0,]
    if (!(multisite)) sXc <- sXc[sXc[, 2] > 0,]
  }
  # Conditional that orders the papercutoff3 sitelists based on the papercutoff
  # sitelist.
  if (sitelist == "papercutoff3") {
    sXc_alt <- SubfunctionCall(SitesXCounts, sitelist="papercutoff")
    sites <- rownames(sXc)
    sites_alt <- rownames(sXc_alt)
    sites <- sites_alt[which(sites_alt %in% sites)]
    sXc <- sXc[sites, ]
  }
  if (sitelist == "randthrp") {
    ind_none <- which(rownames(sXc) == "None")
    sXc <- rbind(sXc[-ind_none, ], sXc[ind_none, , drop=FALSE])
  }
  return(sXc)
}

GetRemovedSites <- function(sXc) {
  # Identify unidentified sites ("CACACAC", "GCACTTTA", etc.)
  grep("^(?:A|C|G|T)+$", rownames(sXc), value=TRUE)
}


CompetitorKmers <- function(mirna, condition, kmer_len, off,
                            experiment="equilibrium", n_constant=5, addcomp=0,
                            wobble=FALSE, alt_mirna_comp=FALSE) {
  ext <- sprintf("_%s_k%s_off%s", n_constant, kmer_len, off)
  if (addcomp != "0") ext <- sprintf("%s_addcomp%s", ext, addcomp)
  if (wobble) ext <- sprintf("%s_wob", ext)
  if (class(alt_mirna_comp) == "character") {
    ext <- sprintf("%s_altcomp_%s", ext, alt_mirna_comp)
  }
  path <- SubfunctionCall(GetAnalysisPath, analysis_type="competitor_kmers")
  counts <- read.table(path, sep="\t", row.names=NULL, header=FALSE)
  kmers <- counts[, 1]
  matrix_out <-
  counts_out <- sapply(unique(kmers), function(kmer) {
    row_ind <- which(counts[, 1] == kmer)[1]
    counts[row_ind, 2:ncol(counts)]
  })
  counts_out <- matrix(counts_out, nrow=length(counts_out), ncol=1)
  rownames(counts_out) <- unique(kmers)
  counts_out
}


SingleSitesXCounts <- function(mirna, experiment="equilibrium", n_constant=5,
                         sitelist="resubmissionfinal", uniq=FALSE, buffer=FALSE) {
  sXc <- SitesXCounts(mirna, experiment, n_constant, sitelist, uniq=uniq, buffer=buffer)
  msXc <- SitesXCounts(mirna, experiment, n_constant, sitelist, uniq=uniq, buffer=buffer, multisite=TRUE)
  regex_single <- "(?:,|\\|)"
  inds = grep(regex_single, rownames(msXc), invert=TRUE)
  msXc[inds, ][rownames(sXc), ]
}

DoubleSitesXCounts <- function(mirna, experiment="equilibrium", n_constant=5,
                         sitelist="resubmissionfinal") {
  sXc <- SitesXCounts(mirna, experiment, n_constant, sitelist)
  msXc <- SitesXCounts(mirna, experiment, n_constant, sitelist, multisite=TRUE)

  regex_split         <- ",\\(\\d*\\),"
  regex_split_capture <- ",\\((\\d*)\\),"
  regex_site <- "[^,]*"
  regex_site_capture <- "([^,]*)"
  double_sites = grep(sprintf("^%s%s%s$", regex_site, regex_split,
                              regex_site), rownames(msXc), perl=TRUE)
  nonoverlapping_sites = grep("\\|", rownames(msXc), invert=TRUE, perl=TRUE)
  double_nonoverlapping_sites = intersect(double_sites, nonoverlapping_sites)
  sXc_2s_dist <- msXc[double_nonoverlapping_sites, ]
  site_2s_l <- gsub(sprintf("^%s%s%s$", regex_site_capture,
                           regex_split_capture, regex_site_capture),
                   rownames(sXc_2s_dist), replace="\\1", perl=TRUE)
  site_2s_r <- gsub(sprintf("^%s%s%s$", regex_site_capture,
                           regex_split_capture, regex_site_capture),
                   rownames(sXc_2s_dist), replace="\\3", perl=TRUE)

  double.site.df <- data.frame(left=site_2s_l, right=site_2s_r)
  double.site.df <- cbind(double.site.df, sXc_2s_dist)
  double.pairs <- expand.grid(setdiff(rownames(sXc), "None"),
                              setdiff(rownames(sXc), "None"),
                              stringsAsFactors=FALSE)
  double.pairs.strings <- apply(double.pairs, 1, paste0, collapse=" ")

  double.site.combined <- t(sapply(double.pairs.strings, function(col) {
    splits <- unlist(strsplit(col, split=" "))
    colSums(subset(double.site.df,
                   left==splits[1] & right==splits[2])[, 4:ncol(double.site.df)])
  }))
  double.site.combined
}


extension_list <- c("8mer",
                    "7mer-m8",
                    "10mer-m1.10mmA9mmA10bG7",
                    "9mer-m1.9mmA9bG7")




# Related to flanking nucleotide Kds ###########################################
SiteFlanksXCounts <- function(mirna, site, experiment="equilibrium", n_constant=5,
                         sitelist="resubmissionfinal", uniq=FALSE, buffer=FALSE) {


  analysis_type <- "site_count_tables"
  condition <- paste0(site, "_flanking")
  ext <- sprintf("_%s_%s", n_constant, sitelist)
  if (uniq) {
    ext <- spring("%s_uniq", ext)
  }
  if (buffer) {
    ext <- sprintf("%s_buffer3p", ext)
  }
  file <- SubfunctionCall(GetAnalysisPath, ext=ext)
  sitesXcounts <- read.table(file, stringsAsFactors=FALSE)
  colnames.temp <- colnames(sitesXcounts)[3 : ncol(sitesXcounts)]
  colnames.new <- sapply(colnames.temp, function(name) {
    return(unlist(strsplit(name, split="A", fixed=TRUE))[2])
  })
  colnames(sitesXcounts)[3 : ncol(sitesXcounts)] <- colnames.new
  return(sitesXcounts)
}



SitesAndSingleFlanksXCounts <- function(mirna, site, experiment="equilibrium",
                                        n_constant=5, sitelist="resubmissionfinal",
                                        buffer=FALSE) {
  sXc <- SubfunctionCall(SitesXCounts)
  fXc <- SubfunctionCall(SiteFlanksXCounts)
  site.ind <- which(rownames(sXc) == site)
  left.flanks <- substr(rownames(fXc), 1, 2)
  right.flanks <- substr(rownames(fXc), 4, 5)
  rownames(fXc) <- paste0(site, "|", rownames(fXc))
  rbind(sXc[0:(site.ind - 1), ],
        fXc,
        sXc[(site.ind + 1): nrow(sXc), ])
}

GetFlankKds <- function(mirna, site, experiment="equilibrium", n_constant=5,
                        sitelist="resubmissionfinal", combined=TRUE, buffer=FALSE) {

  analysis_type <- "kds_PAPER"
  condition <- ""
  ext <- sprintf("%s_%s_%s", n_constant, sitelist, site)
  if (buffer) {
    ext <- sprintf("%s_buffer3p", ext)
  }
  if (!combined) {
    ext <- sprintf("%s_nocombInput", ext)
  }
  ext <- sprintf("%s_PAPER", ext)
  params.file <- SubfunctionCall(GetAnalysisPath, ext=ext)
  # system(sprintf("ls -l %s", params.file))
  data <- fread(params.file, fill=TRUE,header=TRUE,
                       stringsAsFactors=FALSE, showProgress=FALSE)
  colnames(data) <- c("", colnames(data)[-ncol(data)])
  out <- data.frame(data, row.names=1)
  out
}

GetPositionalProgKds <- function(mirna, site, experiment, n_constant=0,
                                 lambda=0, suppcomp=FALSE) {
  analysis_type <- "kds_PAPER"
  condition <- ""
  if (suppcomp) {
    ext <- sprintf("%s_programmedmm_%s_bayes%s", n_constant, site, lambda)
  } else {
    ext <- sprintf("%s_programmed_%s", n_constant, site)
  }
  ext <- sprintf("%s_PAPER", ext)
  params.file <- SubfunctionCall(GetAnalysisPath, ext=ext)
  # system(sprintf("ls -l %s", params.file))
  read.table(params.file, header=TRUE, row.names=1,
             stringsAsFactors=FALSE)
}



SitesAndSingleFlanksKds <- function(mirna, site, experiment="equilibrium",
                                    n_constant=5, sitelist="resubmissionfinal",
                                    buffer=FALSE, combined=TRUE) {
  pars.sites <- SubfunctionCall(EquilPars)
  pars.flanks <- SubfunctionCall(GetFlankKds)
  rownames(pars.flanks) <- paste0(site, "|", rownames(pars.flanks))
  site.ind <- which(rownames(pars.sites) == paste0(site, "_Kd"))
  rbind(pars.sites[0:(site.ind - 1), ],
        pars.flanks,
        pars.sites[(site.ind + 1): nrow(pars.sites), ])
}

GetFullMirnaSiteFlankKds <- function(mirna, site, experiment="equilibrium",
                                     n_constant=5, sitelist="resubmissionfinal",
                                     buffer=FALSE, combined=TRUE) {
  tick <<- tick + 1
  out <- setNames(rep(NaN, length(kFlanks)), kFlanks)
  f_kds <- SubfunctionCall(GetFlankKds)
  site <- gsub("\\(", site, replace="\\\\(")
  site <- gsub("\\)", site, replace="\\\\)")
  f_kds <- f_kds[grep(sprintf("^%s\\|", site), rownames(f_kds), perl=TRUE), ]
  rownames(f_kds) <- gsub(sprintf("^%s\\|(.*)_Kd$", site), replace="\\1",
                          rownames(f_kds), perl=TRUE)
  flanks <- rownames(f_kds)
  out[flanks] <- f_kds[flanks, ]$Mean
  return(out)
}

GetFullMirnaSiteFlankKds_CI <- function(mirna, site, experiment="equilibrium",
                                     n_constant=5, sitelist="resubmissionfinal",
                                     buffer=FALSE, combined=TRUE) {
  tick <<- tick + 1
  out <- setNames(rep(NaN, length(kFlanks)), kFlanks)
  out <- cbind(out, out)
  f_kds <- SubfunctionCall(GetFlankKds)
  site <- gsub("\\(", site, replace="\\\\(")
  site <- gsub("\\)", site, replace="\\\\)")
  f_kds <- f_kds[grep(sprintf("^%s\\|", site), rownames(f_kds), perl=TRUE), ]
  rownames(f_kds) <- gsub(sprintf("^%s\\|(.*)_Kd$", site), replace="\\1",
                          rownames(f_kds), perl=TRUE)
  flanks <- rownames(f_kds)
  colnames(out) <- c("Upper_CI", "Lower_CI")
  out[flanks, 1] <- f_kds[flanks, ]$Upper_CI
  out[flanks, 2] <- f_kds[flanks, ]$Lower_CI 
  return(out)
}

CombinedSiteAndFlankColors <- function(sXc) {
  colors <- rep("gray", nrow(sXc))
  flanks.inds <- grep("\\|", rownames(sXc), perl=TRUE)
  flanks     <- gsub("(^.*)\\|(.*$)", rownames(sXc)[flanks.inds], replace="\\2",
                     perl=TRUE)
  colors[flanks.inds] <- GetColorFunction(flanks)
  colors
}

ExtractSingleFlanks <- function(flank.names, side="5p") {
  if (side == "5p") {
    replace.string <- "\\1"
  } else if (side == "3p") {
    replace.string <- "\\2"
  }
  return(gsub("^(.{2,2}).*(.{2,2})$", flank.names, replace=replace.string,
         perl=TRUE))
}

GetFlankDataFrame <- function(experiment="equilibrium", n_constant=5,
                                sitelist="resubmissionfinal", buffer=FALSE,
                                combined=TRUE, leaveout=FALSE) {
  # Output matrix:
  rbindlist(lapply(kMirnas, function(mirna) {
    if (mirna %in% c("miR-1", "miR-7-23nt")) {
      combined <- FALSE
    }
    if (mirna == "miR-7-23nt") {
      experiment <- "equilibrium2_nb"
    }
    if (mirna == "miR-1") {
      buffer <- TRUE
    }
    rbindlist(lapply(kSeedSites, function(site) {
      # Get each site and miRNA Kd table:
      kds <- SubfunctionCall(GetFlankKds)
      kds <- kds[grep(sprintf("^%s\\|", site), rownames(kds), perl=TRUE), ]
      rownames(kds) <- gsub(sprintf("^%s\\|(.*)_Kd$", site), "\\1", rownames(kds), perl=TRUE)
      flank.5p <- ExtractSingleFlanks(rownames(kds))
      flank.3p <- ExtractSingleFlanks(rownames(kds), side="3p")
      data.frame(logkd=log(kds$Mean),
                 mirna=rep(mirna, nrow(kds)),
                 site=rep(site, nrow(kds)),
                 f5p1=substr(flank.5p, 1, 1),
                 f5p2=substr(flank.5p, 2, 2),
                 f3p1=substr(flank.3p, 1, 1),
                 f3p2=substr(flank.3p, 2, 2))
    }))
  }))
}

GetFlankLinearModel <- function(experiment="equilibrium", n_constant=5,
                                sitelist="resubmissionfinal", buffer=FALSE,
                                combined=TRUE, leaveout=FALSE) {
  if (leaveout != FALSE) {
    kMirnas <- setdiff(kMirnas, leaveout)
  }
  kd.data <- rbindlist(lapply(kMirnas, function(mirna) {
    # if (mirna %in% c("miR-7-23nt")) {
    #   combined <- FALSE
    # }
    if (mirna == "miR-7-23nt") {
      experiment <- "equilibrium2_nb"
      sitelist <- "resubmissionfinal6merm8"
    }
    if (mirna == "miR-1") {
      buffer <- TRUE
    }
    rbindlist(lapply(kSeedSites, function(site) {
      # Get each site and miRNA Kd table:
      kds <- SubfunctionCall(GetFlankKds)
      kds <- kds[grep(sprintf("^%s\\|", site), rownames(kds), perl=TRUE), ]
      rownames(kds) <- gsub(sprintf("^%s\\|(.*)_Kd$", site), "\\1", rownames(kds), perl=TRUE)
      flank.5p <- ExtractSingleFlanks(rownames(kds))
      flank.3p <- ExtractSingleFlanks(rownames(kds), side="3p")
      data.frame(logkd=log(kds$Mean),
                 mirna=rep(mirna, nrow(kds)),
                 site=rep(site, nrow(kds)),
                 f5p1=substr(flank.5p, 1, 1),
                 f5p2=substr(flank.5p, 2, 2),
                 f3p1=substr(flank.3p, 1, 1),
                 f3p2=substr(flank.3p, 2, 2))
    }))
  }))
  kd_data_global <<- kd.data
  out <- lm(logkd ~ mirna*site + f5p1*f5p2 + f3p1*f3p2 - 1, data=kd.data)
  out
}

GetSingleMirnaModel <- function(mirna, experiment="equilibrium", n_constant=5,
                                sitelist="resubmissionfinal", buffer=FALSE,
                                combined=TRUE) {
  lm(logkd ~ site - 1, data=SubfunctionCall(GetMirnaFlankDataFrame))
}

GetSingleMirnaModelNoInt <- function(mirna, experiment="equilibrium", n_constant=5,
                                sitelist="resubmissionfinal", buffer=FALSE,
                                combined=TRUE) {
  lm(logkd ~ site - 1, data=SubfunctionCall(GetMirnaFlankDataFrame))
}



# model1 <- GetSingleMirnaModel("miR-1", buffer=TRUE, combined=FALSE)
# model2 <- GetSingleMirnaModelNoInt("miR-1", buffer=TRUE, combined=FALSE)


GetSingleMirnaFlankModel <- function(mirna, experiment="equilibrium", n_constant=5,
                                sitelist="resubmissionfinal", buffer=FALSE,
                                combined=TRUE) {
    kd.data <- rbindlist(lapply(kSeedSites, function(site) {
      # Get each site and miRNA Kd table:
      kds <- SubfunctionCall(GetFlankKds)
      kds <- kds[grep(sprintf("^%s\\|", site), rownames(kds), perl=TRUE), ]
      rownames(kds) <- gsub(sprintf("^%s\\|(.*)_Kd$", site), "\\1", rownames(kds), perl=TRUE)
      flank.5p <- ExtractSingleFlanks(rownames(kds))
      flank.3p <- ExtractSingleFlanks(rownames(kds), side="3p")
      data.frame(logkd=log(kds$Mean),
                 mirna=rep(mirna, nrow(kds)),
                 site=rep(site, nrow(kds)),
                 f5p1=substr(flank.5p, 1, 1),
                 f5p2=substr(flank.5p, 2, 2),
                 f3p1=substr(flank.3p, 1, 1),
                 f3p2=substr(flank.3p, 2, 2))
    }))
  out <- lm(logkd ~ site + f5p1*f5p2 + f3p1*f3p2 - 1, data=kd.data)
}

GetAverageFlanks <- function(experiment="equilibrium", n_constant=5,
                             sitelist="resubmissionfinal", combined=TRUE, buffer=FALSE) {
  # Get the flanking dinucleotide model trained on the other 5 miRNAs:
  R <- 1.987e-3 # in kcal K-1 mol-1
  T <- 310.15 # in K

  flank.model <- SubfunctionCall(GetFlankLinearModel)
  flank.coeffs <- SubfunctionCall(ParseFlankModelCoefs, model=flank.model)
  flank.lin <- flank.coeffs[[1]]
  flank.5int <- flank.coeffs[[2]]
  flank.3int <- flank.coeffs[[3]]

  flank.lin_alt <<- flank.lin*R*T
  flank.5int_alt <<- flank.5int*R*T
  flank.3int_alt <<- flank.3int*R*T
  test.model <- sapply(kFlanks, function(flank) {

    f5p1 <- substr(flank, 1, 1)
    f5p2 <- substr(flank, 2, 2)
    f3p1 <- substr(flank, 4, 4)
    f3p2 <- substr(flank, 5, 5)
    f5p1.val <- flank.lin[f5p1, "f5p1"]
    f5p2.val <- flank.lin[f5p2, "f5p2"]
    f3p1.val <- flank.lin[f3p1, "f3p1"]
    f3p2.val <- flank.lin[f3p2, "f3p2"]
    if (f5p1 != "A" & f5p2 != "A") {
      int5.val <- flank.5int[paste0(2, f5p2), paste0(1, f5p1)]
    } else {
      int5.val <- 0
    }
    if (f3p1 != "A" & f3p2 != "A") {
      int3.val <- flank.3int[paste0(2, f3p2), paste0(1, f3p1)]    

    } else {
      int3.val <- 0
    }
    sum(f5p1.val, f5p2.val, f3p1.val, f3p2.val, int5.val, int3.val)*R*T
  })
  test.model - mean(test.model, na.rm=TRUE)
}

LeaveOneOutFlankModel <- function(mirna, experiment="equilibrium",
                                  n_constant=5, sitelist="resubmissionfinal",
                                  combined=TRUE, xpos=20, ypos=20) {
  # Get the flanking dinucleotide model trained on the other 5 miRNAs:
  flank.model <- SubfunctionCall(GetFlankLinearModel, leaveout=mirna)
  flank.model_temp <<- flank.model$model[1:256, ]
  flank.model <<- flank.model
  if (mirna == "miR-7-23nt") {
    experiment <- "equilibrium2_nb"
    sitelist <- "resubmissionfinal6merm8"
  }
  if (mirna == "miR-1") {
    buffer <- TRUE
  }
  # Get the site coefficients for this miRNA:
  site.mir.model <- SubfunctionCall(GetSingleMirnaModel)
  site.mir.model <<- site.mir.model
  # Get the site coefficients:
  site.coeffs <- site.mir.model$coefficients
  names(site.coeffs) <- gsub("^site(.*)$", names(site.coeffs),
                             replace="\\1", perl=TRUE)
  leave_out_coefficients <- c(site.mir.model$coefficients,
                                   flank.model$coefficients[grep("^f(?:5|3p)",
                                                                 names(flank.model$coefficients),
                                                                 perl=TRUE, value=TRUE)])


  switch_coeffs <- flank.model$coefficients[grep("^f(?:5|3p)",
                                                                 names(flank.model$coefficients),
                                                                 perl=TRUE, value=TRUE)]

  data.flanks <- SubfunctionCall(GetMirnaFlankDataFrame)
  lm_site <- lm(logkd ~ site - 1, data=data.flanks)
  lm_new <- lm(logkd ~ site + f5p1*f5p2 + f3p1*f3p2 - 1, data=data.flanks)
  lm_new$coefficients <- leave_out_coefficients
  lm_flanks_only_new <- lm(logkd ~ f5p1*f5p2 + f3p1*f3p2, data=data.flanks)
  lm_flanks_only_switch <- lm_flanks_only_new
  lm_flanks_only_switch$coefficients <- c(lm_flanks_only_new$coefficients[1], switch_coeffs)
  names(lm_flanks_only_switch$coefficients)[1] <- "(Intercept)"
  flanks_only_pred <<- predict(lm_flanks_only_new, data=data.flanks)
  flanks_only_switch_pred <<- predict(lm_flanks_only_switch, data=data.flanks)
  log_kds_new <- data.flanks$logkd - flanks_only_switch_pred
  log_site_new <- data.flanks$logkd - flanks_only_pred
  data.flanks_alt <- data.flanks
  data.flanks_alt$logkd <- log_kds_new
  lm_site_corrected <- lm(logkd ~ site - 1, data=data.flanks_alt)

  x_vars <- data.flanks$logkd

  y_vars <- predict(lm_new, data=data.flanks)


  lm_new$coefficients <- c(lm_site_corrected$coefficients,
                           lm_flanks_only_switch$coefficients[-1])

  y_vars_new <- predict(lm_new, data=data.flanks)

  # Get the flank coefficients:
  flank.coeffs <- SubfunctionCall(ParseFlankModelCoefs, model=flank.model)
  flank.lin <- flank.coeffs[[1]]
  flank.5int <- flank.coeffs[[2]]
  flank.3int <- flank.coeffs[[3]]
  data.flanks <- SubfunctionCall(GetMirnaFlankDataFrame)
  data.flanks <<- data.flanks
  test.model <- apply(data.flanks, 1, function(row) {
    site <- paste0("site", row["site"])
    f5p1 <- row["f5p1"]
    f5p2 <- row["f5p2"]
    f3p1 <- row["f3p1"]
    f3p2 <- row["f3p2"]
    site.val <- site.coeffs[row["site"]]
    f5p1.val <- flank.lin[f5p1, "f5p1"]
    f5p2.val <- flank.lin[f5p2, "f5p2"]
    f3p1.val <- flank.lin[f3p1, "f3p1"]
    f3p2.val <- flank.lin[f3p2, "f3p2"]
    if (f5p1 != "A" & f5p2 != "A") {
      int5.val <- flank.5int[paste0(2, f5p2), paste0(1, f5p1)]    
    } else {
      int5.val <- 0
    }
    if (f3p1 != "A" & f3p2 != "A") {
      int3.val <- flank.3int[paste0(2, f3p2), paste0(1, f3p1)]    
    } else {
      int3.val <- 0
    }
    sum(site.val, f5p1.val, f5p2.val, f3p1.val, f3p2.val, int5.val, int3.val)
  })
  alt_out <- test.model - mean(test.model, na.rm=TRUE) + mean(data.flanks$logkd, na.rm=TRUE)
  y_vars <- y_vars - mean(y_vars, na.rm=TRUE) + mean(data.flanks$logkd, na.rm=TRUE)
  colors <- kSiteColors[as.character(data.flanks$site)]

  y_vars_new
}


GetMirnaFlankDataFrame <- function(mirna, experiment="equilibrium",
                                   n_constant=5, sitelist="resubmissionfinal",
                                   buffer=FALSE, combined=TRUE) {
  rbindlist(lapply(kSeedSites, function(site) {
      # Get each site and miRNA Kd table:
      kds <- SubfunctionCall(GetFlankKds)
      kds <- kds[grep(sprintf("^%s\\|", site), rownames(kds), perl=TRUE), ]
      rownames(kds) <- gsub(sprintf("^%s\\|(.*)_Kd$", site), "\\1", rownames(kds), perl=TRUE)
      flank.5p <- ExtractSingleFlanks(rownames(kds))
      flank.3p <- ExtractSingleFlanks(rownames(kds), side="3p")
      data.frame(logkd=log(kds$Mean),
                 mirna=rep(mirna, nrow(kds)),
                 site=rep(site, nrow(kds)),
                 f5p1=sapply(flank.5p, substr, start=1, stop=1),
                 f5p2=sapply(flank.5p, substr, start=2, stop=2),
                 f3p1=sapply(flank.3p, substr, start=1, stop=1),
                 f3p2=sapply(flank.3p, substr, start=2, stop=2))
  }))
}

ParseFlankModelCoefs <- function(model) {
  coefs <- c(f5p1A=0, f5p2A=0, f3p1A=0, f3p2A=0, model$coefficients)
  coefs_error <- c(f5p1A=0, f5p2A=0, f3p1A=0, f3p2A=0, summary(model)[[4]][,2])
  CI_lower <- c(f5p1A=0, f5p2A=0, f3p1A=0, f3p2A=0, confint(model)[, 1])
  CI_upper <- c(f5p1A=0, f5p2A=0, f3p1A=0, f3p2A=0, confint(model)[, 2])
  lin <- sapply(c("f5p1", "f5p2", "f3p1", "f3p2"), function(string) {
    grep.key <- paste0("^", string, ".{1}$")
    out <- coefs[grep(grep.key, names(coefs), perl=TRUE)]
    names(out) <- sapply(names(out), substr, start=5, stop=5)
    out
  })

  lin_error_upper <- sapply(c("f5p1", "f5p2", "f3p1", "f3p2"), function(string) {
    grep.key <- paste0("^", string, ".{1}$")
    out <- CI_upper[grep(grep.key, names(CI_upper), perl=TRUE)]
    names(out) <- sapply(names(out), substr, start=5, stop=5)
    out
  })


  lin_error_lower <- sapply(c("f5p1", "f5p2", "f3p1", "f3p2"), function(string) {
    grep.key <- paste0("^", string, ".{1}$")
    out <- CI_lower[grep(grep.key, names(CI_lower), perl=TRUE)]
    names(out) <- sapply(names(out), substr, start=5, stop=5)
    out
  })

  ints <- lapply(c(5, 3), function(flank){
    sapply(c("1C", "1G", "1T"), function(string) {
    grep.key <- paste0("^f", flank, "p", string, ":")
    out <- coefs[grep(grep.key, names(coefs), perl=TRUE)]
    len <- nchar(names(out))
    names(out) <- sapply(names(out), substr, start=len - 1, stop=len)
    out
    })  
  })
  list(lin=lin, int.5p=ints[[1]], int.3p=ints[[2]],
       lin_error_upper=lin_error_upper, lin_error_lower=lin_error_lower)
}


GetFlankLMCoefs <- function(experiment="equilibrium", n_constant=5,
                            sitelist="resubmissionfinal") {
  model <- SubfunctionCall(GetFlankLinearModel)
  SubfunctionCall(ParseFlankModelCoefs)
}


GetPairingFlankData <- function(mirna, site, condition,
                                experiment="equilibrium", n_constant=5,
                                sitelist="resubmissionfinal", buffer=FALSE,
                                mir.start=1, mir.stop=14, noconstant=FALSE,
                                alt_mir_exp_cond=NULL,
                                old=FALSE) {
  ############# THIS FUNCTION TAKES DATA GENERATED BY ##########################
  # folder = "/structural_analysis_PAPER_realfinal/"
  # file = paste0(condition, "_", n_constant, "_", sitelist, "_", mir.start, "-", mir.stop,
  #               collapse="")
  # path <- paste0(kSolexaDir, mirna, "/", experiment, folder, site, "/", file,
  #                ".txt")
  # print(path)
  # if (old) {
  #   analysis_type <- "structural_analysis"  
  # } else {
  analysis_type <- "structural_analysis"  
  # }
  ext <- sprintf("_%s_%s_%s-%s", n_constant, sitelist, mir.start, mir.stop)
  if (!is.null(alt_mir_exp_cond)) {
    ext <- sprintf("%s_alt_%s", ext, alt_mir_exp_cond)
  }
  if (buffer) {
    ext <- sprintf("%s_buffer3p", ext)
  }
  if (noconstant == TRUE) {
    ext <- sprintf("%s_noconstant", ext)
  }
  path <- SubfunctionCall(GetAnalysisPath,
                              condition=sprintf("%s/%s", site, condition),
                              ext=ext)
  system(sprintf("ls -l %s", path))
  # data <- read.table(path, header=TRUE)
  data <- data.frame(read.table(path, header=TRUE))
  return(data)
}

NumFlanks <- function(data) {
  flank_abundances <- aggregate(. ~ flank, data, function(x) {
    length(x)
  })
  flank_num <- flank_abundances[, 2]
  names(flank_num) <- flank_abundances$flank
  return(flank_num)
}



SampleByDinucleotideEnrichment <- function(dist.sample, dist.enriched,
                                           samplesize) {
  weights <- NumFlanks(dist.enriched)/NumFlanks(dist.sample)
  inds <- sample(1:nrow(dist.sample), replace=TRUE, size=samplesize,
                 prob=weights[dist.sample$flank])
  return(inds)
}

GetPlFoldCorrelationMatrix <- function(mirna, site, experiment="equilibrium",
                                       condition="I", n_constant=5,
                                       sitelist="resubmissionfinal", buffer=FALSE,
                                       test=FALSE) {

  kds <- SubfunctionCall(GetFlankKds)
  kds <- kds[grep(sprintf("^%s\\|", site), rownames(kds)), ]
  kd_names <- gsub(sprintf("^%s\\|(.*)_Kd", site), rownames(kds), replace="\\1")
  kds <- kds$Mean
  names(kds) <- kd_names
  ext <- sprintf("_%s_%s_%s", n_constant, sitelist, site)
  if (buffer) {
    ext <- sprintf("%s_buffer3p", ext)
  }
  if (test) {
    ext <- sprintf("%s_test", ext)
  }
  analysis_type <- "plfold_windows_by_flanks"
  path <- SubfunctionCall(GetAnalysisPath, ext=ext)
  data <- data.matrix(data.frame(fread(path, sep="\t"), row.names=1))
  rownames <- unique(sapply(rownames(data), gsub, pattern="^(.*)_w(.*)$",
                        replacement="\\1", perl=TRUE))
  colnames <- unique(sapply(rownames(data), gsub, pattern="^(.*)_w(.*)$",
                        replacement="\\2", perl=TRUE))
  cors_list <- apply(data, 1, cor, y=-log(kds), use="pairwise.complete.obs")
  
  t(matrix(cors_list, nrow=length(rownames), ncol=length(colnames),
           dimnames=list(rownames, colnames)))
}

MakeFlankAveragesNewSingle <- function(mirna, site, win, experiment="equilibrium",
                                    condition="I_combined", n_constant=5,
                                    sitelist="resubmissionfinal") {
  pl_data <- SubfunctionCall(GetPlFoldDataNew)
  flanks <- pl_data[, 2]
  pl_data <- t(t(pl_data[,-c(1, 2)]^(1/win)))
  out <- aggregate(pl_data,by=flanks,FUN=GeoMean)
  flanks.new <- out[,1]
  out <- out[,-1]
  rownames(out) <- flanks.new
  out <- out[order(flanks.new),]
  return(out)
}

# GetPlFoldDataNew <- function(mirna, experiment, condition, n_constant, sitelist,
#                             site, win) {
#   extension <- paste0("_", n_constant, "_", sitelist, "_", site,  "_win", win-1)
#   path <- GetAnalysisPath(mirna, experiment, condition, "plfold_2018_PAPER",
#                           ext=extension)
#   print(path)
#   data <- fread(file=path, header=TRUE, stringsAsFactors=FALSE, sep="\t")
#   return(data)
# }

SampleBySiteAccess <- function(dist.sample, samplesize,
                               prob.factor, exponent=1) {
  inds <- sample_int_rej(nrow(dist.sample), size=samplesize,
                         prob=dist.sample[[prob.factor]]^exponent)
  return(inds)
}



###############################################################################

# KmersXCountsVector <- function(mirna, experiment, condition, n_constant, kmer_len,
#                          n_ex) {
#   ext <- sprintf("_%s_k%s", n_constant, kmer_len)
#   if (n_ex > 0) {
#     ext <- sprintf("%s_ex:%s", ext, extension_list[1])
#   }
#   if (n_ex > 1) {
#     for (i in seq(1,n_ex )[0:(n_ex-1)]) {
#       ext <- sprintf("%s,%s", ext, extension_list[i + 1])
#     }
#   }
#   message(ext)
#   path <- SubfunctionCall(GetAnalysisPath, analysis_type="kmers_cutoff_final", ext=ext)
#   kXc <- data.frame(fread(path, sep="\t", header=FALSE, stringsAsFactors=FALSE),
#                     row.names=1, skip=1)
#   colnames(kXc) <- condition
#   kXc
# }

KmersXCountsVector <- function(mirna, experiment, condition, n_constant,
                               kmer_len, n_ex, buffer=FALSE) {
  ext <- sprintf("_%s_k%s", n_constant, kmer_len)
  if (n_ex > 0) {
    ext <- sprintf("%s_ex%s", ext, n_ex)
  }
  if (buffer) {
    ext <- sprintf("%s_buffer3p", ext)
  }
  path <- SubfunctionCall(GetAnalysisPath, analysis_type="kmers_cutoff_resub", ext=ext)
  message(path)
  kXc <- data.frame(fread(path, sep="\t", header=FALSE, stringsAsFactors=FALSE,
                          skip=1),
                    row.names=1)[, 1, drop=FALSE]
  colnames(kXc) <- condition
  kXc
}

KmersXCounts <- function(mirna, experiment, n_constant, kmer_len, n_ex) {
  if (mirna == "miR-7-23nt" & experiment == "equilibrium2_nb") {
    conditions <- c(40, 12.6, 1.26, 0.4, 0)
  } else {
    conditions <- c(40, 12.6, 4, 1.26, 0.4, 0)
  }
  kXc <- SubfunctionCall(KmersXCountsVector, conditio="I")
  for (condition in conditions) {
    kXc <- cbind(kXc, SubfunctionCall(KmersXCountsVector))
  }
  kXc
}

KmersXCountsAllInput <- function(n_constant, kmer_len) {
  # Loop over the original miRNAs.
  experiment <- "equilibrium"
  condition <- "I"
  n_ex <- 0
  kXc <- do.call("cbind", lapply(kMirnas[-6], function(mirna) {
    SubfunctionCall(KmersXCountsVector)
  }))
  colnames(kXc) <- kMirnas[-6]
  kXc
  # Include Namita's let-7a-21nt experiment:
  mirnas_nb <- c("let-7a-21nt", "miR-7-23nt", "miR-7-24nt", "miR-7-24nt")
  exps_nb <- c("equilibrium_nb", "equilibrium_nb", "equilibrium2_nb",
               "equilibrium3_nb")
  metadata_nb <- cbind(mirnas_nb, exps_nb)
  print(metadata_nb)
  kXc_nb <- do.call("cbind", apply(metadata_nb, 1, function(row) {
    print(row)
    mirna <- row[1]
    experiment <- row[2]
      SubfunctionCall(KmersXCountsVector)
  }))
  colnames(kXc_nb) <- c("let-7a_nb", "miR-7_nb", "miR-7_2nb", "miR-7_3nb")
  cbind(kXc, kXc_nb)
}


# KmersXCountsAllInput <- function(n_constant, kmer_len) {
#   print(n_constant)
#   # Loop over the original miRNAs.
#   experiment <- "equilibrium"
#   condition <- "I"
#   kXc <- do.call("cbind", lapply(kMirnas[-6], function(mirna) {
#     SubfunctionCall(KmersXCountsVector, experiment="equilibrium",
#                     condition="I", n_ex=0)
#   }))
#   colnames(kXc) <- kMirnas[-6]
#   kXc
# }




GetPositionalKmers <- function(mirna, experiment, condition, n_constant,
                               kmer_len, addcomp=FALSE) {
  if (experiment %in% c("kinetics_pilot", "kinetics")) {
    ext_p <- sprintf("_pulse_%s_k%s", n_constant, kmer_len)
    ext_c <- sprintf("_chase_%s_k%s", n_constant, kmer_len)
    paths <- sapply(c(ext_p, ext_c), function(ext) {
      SubfunctionCall(GetAnalysisPath, analysis="positional_kmers")
    })
    out_p <- fread(paths[1], sep="\t", data.table=FALSE)
    out_c <- fread(paths[2], sep="\t", data.table=FALSE)
    rownames(out_p) <- GetKmerList(kmer_len)
    rownames(out_c) <- rownames(out_p)
    return(list(out_p, out_c))
  } else {
    ext <- sprintf("_%s_k%s", n_constant, kmer_len)
    if (class(addcomp) == "character") {
      ext <- sprintf("%s_addcomp%s", ext, addcomp)
    }
    path <- GetAnalysisPath(mirna, experiment, condition, "positional_kmers",
                            ext=ext)
    out <- fread(path, sep="\t", data.table=FALSE)
    rownames(out) <- GetKmerList(kmer_len)
    return(out)
  }
}



GetPositionalKmersEquilPath <- function(mirna, experiment, condition, n_constant,
                               kmer_len) {
  ext <- sprintf("_%s_k%s", n_constant, kmer_len)
  path <- GetAnalysisPath(mirna, experiment, condition, "positional_kmers",
                          ext=ext)
  return(path)
}


# Plots complemtarity to oligo sequences used in the AGO-RBNS experiments.
# Uses the SequenceList objects in Lists.R.
PlotOverlapIdentity <- function(mirna, seq_type, experiment="equilibrium",
                                n_constant=5, kmer_len=8, n_ex=0, height=4,
                                width=4, xpos=20, ypos=20, pdf.plot=FALSE) {

  output <- SubfunctionCall(KmersXCountsVector, condition="I")
  kmers <- rownames(output)[2:nrow(output)]
  output <- as.numeric(unlist(output[2:nrow(output), 1, drop=FALSE]))
  output <- as.matrix(output, nrow=length(output), ncol=1)
  rownames(output) <- kmers
  colnames(output) <- "I"
  if (mirna == "miR-7-23nt") {
    conditions <- c("40", "12.6", "1.26", "0.4", "0")
  } else {
    conditions <- c("40", "12.6", "4", "1.26", "0.4", "0")    
  }
  for (condition in conditions) {
    output_i <- SubfunctionCall(KmersXCountsVector)
    output_i <- as.numeric(unlist(output_i[2:nrow(output_i), 2, drop=FALSE]))
    output_i <- as.matrix(output_i, nrow=length(output_i), ncol=1)
    rownames(output_i) <- kmers
    colnames(output_i) <- condition
    output <- cbind(output, output_i)
  }
  output <- MatNorm(output)
  output_R <- output/(output[, 1])

  overlap_seq <- SequenceList[[seq_type]][mirna]
  overlap_seq_use <- RevComplement(overlap_seq)

  range_ks <- nchar(overlap_seq_use) - kmer_len + 1
  xmin <- 0
  xmax <- nchar(overlap_seq_use)
  ymin <- 0.1
  ymax <- 300
  SubfunctionCall(GenericFigureSaveFile)
  BlankPlot(log='y')
  AddLogAxis(2, label="Enrichment")
  overlap_seq_use_kmers <- sapply(1:range_ks, function(start) {
      substr(overlap_seq_use, start=start, stop=start + kmer_len - 1)
  })
  overlap_seq_use_kmer_inds <- sapply(overlap_seq_use_kmers, function(kmer) {
    which(rownames(output_R) == kmer)
  })
  R_inds <- output_R[rev(overlap_seq_use_kmer_inds), -1]

  R_non_inds <- output_R[-c(overlap_seq_use_kmer_inds), -1]
  R_non_inds_average <- exp(colMeans(log(R_non_inds)))
  R_non_inds_SD <- apply(R_non_inds, 2, function(column) {
                          exp(sd(log(column)))
                         })
  R_inds_x <- 1:nrow(R_inds) + kmer_len/2 - 0.5
  R_noninds_x <- setdiff(1:nchar(overlap_seq_use) + 0.5, R_inds_x)


  sapply(rev(c(conditions)), function(condition) {
    ys <- R_inds[, condition]
    names(ys) <- R_inds_x
    ys_noninds <- rep(R_non_inds_average[condition], length(R_noninds_x))
    names(ys_noninds) <- R_noninds_x
    ys_all <- c(ys, ys_noninds)
    ys_all <- ys_all[order(as.numeric(names(ys_all)))]
    x_all <- as.numeric(names(ys_all))
    x_all <<- x_all
    ys_all <<- ys_all
    lines(x_all, ys_all, col=kEquilCondColors[condition], type="o")
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
  xy <- GetPlotFractionalCoords(0.05, 0.95, log='y')

  text(xy[1], xy[2], labels=mirna, adj=0, xpd=NA)
  xy <- GetPlotFractionalCoords(0.5, 1, log='y')

  legend(xy[1], xy[2], legend=paste0(conditions, "%"), col=kEquilCondColors[conditions],
         bty="n", lwd=1)
  text(1:nchar(overlap_seq), 0.2, unlist(strsplit(overlap_seq, split="")),
       adj=c(0.5, 0.5), xpd=NA)
  xy <- GetPlotFractionalCoords(0.5, 0, log='y')
  text(xy[1], xy[2], labels=SequenceListLabels[seq_type], adj=0.5, xpd=NA)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
} 

GetPositionalSites <- function(mirna, experiment, condition, n_constant=5,
                               sitelist="resubmissionextendedfinal",
                               single=FALSE, buffer=FALSE) {
  ext <- sprintf("_%s_%s", n_constant, sitelist)
  if (buffer) ext <- sprintf("%s_buffer3p", ext)
  if (single) {
    analysis_type <- "single_pos_counts"
  } else {
    analysis_type <- "top_pos_counts"
  }
  path <- SubfunctionCall(GetAnalysisPath)
  print(path)
  system(sprintf("ls -l %s", path))
  file <- data.frame(fread(path), row.names=1)
  file
}

GetKmerIndex <- function(x) {
  nt_inds <- c("A"=0, "C"=1, "G"=2, "T"=3)
  bases <- 4^(seq(nchar(x)-1, 0))
  values <- nt_inds[unlist(strsplit(x, split=""))]
  return(sum(values*bases) + 1)
}

EquilPars <- function(
  mirna, experiment="equilibrium", n_constant=5, sitelist="resubmissionfinal",
  uniq=FALSE, buffer=FALSE, mirna_start=FALSE, combined=TRUE, singleonly=FALSE,
  Xval=FALSE, L=FALSE, compcorrect=FALSE, wobble=FALSE, tpomit=FALSE,
  tpomit2=FALSE, tp2rep=FALSE, minkd=FALSE, AGOfixedbypass=FALSE, nbomitc=FALSE,
  sorted=FALSE, collapsemm=FALSE, sumseed=FALSE, start_mm=FALSE, stop_mm=FALSE,
  new=FALSE, new2=FALSE
) {
  # print("In EquilPars")
  # print(sprintf("start_mm: %s", start_mm))
  # print(sprintf("stop_mm: %s", stop_mm))
  # print(start_mm)
  # print(stop_mm)
  # print(class(start_mm))
  # print(class(stop_mm))
  analysis_type <- "kds_PAPER"
  condition <- ""
  if ((class(start_mm) == "numeric" | class(start_mm) == "integer") &
      (class(stop_mm) == "numeric" | class(stop_mm) == "integer")) {
      str_mm <- sprintf("m%s.%smmsd_", start_mm, stop_mm)
  } else {
    str_mm <- ""
  }
  ext <- sprintf("%s_%s%s", n_constant, str_mm, sitelist)
  if (sitelist %in% kKmerSiteLists) {
    ext <- sprintf("%s_%s-%s", ext, mirna_start, as.integer(mirna_start) + 3)
  }
  if (uniq) {
    ext <- sprintf("%s_uniq", ext)
  }
  if (buffer) {
    ext <- sprintf("%s_buffer3p", ext)
  }
  if (!combined) {
    ext <- sprintf("%s_nocombInput", ext)
  }
  if (singleonly) {
    ext <- sprintf("%s_singleonly", ext)
  }
  if (Xval) {
    ext <- sprintf("%s_Xval", ext)
  }
  if (L) {
    ext <- sprintf("%s_L%s", ext, L)
  }
  if (compcorrect) {
    ext <- sprintf("%s_compcorrect", ext)
    if (wobble) {
      ext <- sprintf("%s_wobble", ext)
    }
  }
  if (tpomit) {
    ext <- sprintf("%s_tpomit", ext)
  }
  if (tpomit2) {
    ext <- sprintf("%s_tpomit2", ext)
  }
  if (tp2rep) {
    ext <- sprintf("%s_rep%s", ext, tp2rep)
  }
  if (minkd) {
    minkd <- gsub("^(.*)e\\-0?(.*)$", sprintf("%.e", minkd),
                  replace="\\1e\\-\\2")
    ext <- sprintf("%s_minkd%s", ext, minkd)
  }
  if (AGOfixedbypass) {
    ext <- sprintf("%s_AGOfixedbypass", ext)
  }
  if (nbomitc) {
    ext <- sprintf("%s_nbomitc", ext)
  }
  if (collapsemm & sitelist %in% c("programmed_suppcomp", "bipartite_random_suppcomp")) {
    ext <- sprintf("%s_collapsemm", ext)
  }
  if (sumseed & sitelist %in% c("progthrp_suppcomp", "randthrp_suppcomp", "randthrp_comp")) {
    ext <- sprintf("%s_sumseed", ext)
  }
  if (new) {
    ext <- sprintf("%s_new", ext)
  }
  if (new2) {
    ext <- sprintf("%s_new2", ext)
  }
  

  # Add 'PAPER' to every file
  # TODO remove this extension from everything.
  ext <- sprintf("%s_PAPER", ext)
  # This is a suffix extension for the programmed libraries, looking at the
  # actual range of values fitted during the optimization to be able to compare
  # to the confidence intervals that I am calculating.
  if (sorted) {
    ext <- sprintf("%s_sorted", ext)
  }
  # Get path to file.
  params.file <- SubfunctionCall(GetAnalysisPath, ext=ext)
  # system(sprintf("ls -l %s", params.file))
  # Get the parameter file.
  params <- data.frame(read.table(params.file, header=TRUE, row.names=1,
                       sep= "\t", stringsAsFactors=FALSE))
  return(params)
}

GetModelOccupancies <- function(mirna, experiment="equilibrium", n_constant=5,
                                sitelist="resubmissionfinal", buffer=FALSE,
                                combined=TRUE) {
  sXc <- SubfunctionCall(SitesXCounts)
  pars.matrix <- SubfunctionCall(EquilPars)
  pars.model <- log10(pars.matrix$Mean)
  names(pars.model) <- rownames(pars.matrix)
  names(pars.model)[length(pars.model)] <- "AGO"
  names(pars.model)[length(pars.model) - 1] <- "bg"
  data <- GetDataEquil(sXc)
  l <- SubfunctionCall(GetInputEquil)
  plot(l, sXc[, 1], log='xy')
  # data <- data[, 1:5]
  # Ago dilutio in the data:
  A.dil <- sapply(colnames(data), as.numeric)
  # Checking using the AGO concnetration from the model:
  A.stock.measured <- pars.matrix[sprintf("AGO_%s", mirna), 1]
  model <- EquilSingleSiteModelFreq(pars.model, sXc, A.dil, addbg=FALSE, combined=combined)
  model.norm <- apply(model, 2, Norm)
  A.pM <- A.dil/100*A.stock.measured*1000
  return(model.norm)
}


# Plots the Line graphs of the read counts for each miRNA-expeiment pair:
PlotReadDistribution <- function(mirna, experiment, makepdf=FALSE) {
  conditions <- c("I", "40", "12.6", "4", "1.26", "0.4", "0")
  cond_colors <- c("black", "violet", "blue", "forestgreen", "orange", "red", "gray")
  names(cond_colors) <- conditions
  xmin <- 1
  xmax <- 200
  ymin <- 1
  ymax <- 1e8
  par(kPlotParameters)
  BlankPlot(log='xy')
  AddLogAxis(1, label="Occurences of read")
  AddLogAxis(2, label="Read frequency")
  if (mirna == "miR-7-23nt" & experiment == "equilibrium2_nb") {
    conditions <- conditions[-4]
  }
  sapply(conditions, function(condition) {
    path <- GetAnalysisPath(mirna, experiment, condition,
                            analysis_type="read_count_distributions")
    data <- try(fread(path, data.table=TRUE))
    if (length(class) >= 1 & class(data) != "try-error") {
      data <- data[order(data[[2]]),]
      lines(data[[2]], data[[1]], col=cond_colors[condition])
    }
  })      
  xy <- GetPlotFractionalCoords(fx=0.7, fy=1, log='xy')
  legend(x=xy[1], y=xy[2], legend=conditions, col=cond_colors[conditions], lwd=2, bty="n")
  xy <- GetPlotFractionalCoords(fx=0.1, fy=0.95, log='xy')
  text(x=xy[1], y=xy[2], labels=mirna, adj=0)
  if (makepdf) {
    dev.copy2pdf(file=sprintf("figures/180327_%s_%s_counthistogram.pdf", mirna,
                              experiment))
  }
}

Check7merm8Against7merA1 <- function(mirna, experiment="equilibrium",
                                       n_constant=5, sitelist="resubmissionfinal",
                                       combined=TRUE, buffer=FALSE, singlyonly=TRUE) {
  kds <- SubfunctionCall(EquilPars)
  kdratio <- kds["7mer-A1_Kd", 2]/kds["7mer-m8_Kd", 2]
  message(mirna)
  message(kdratio)
}

CheckAll7merm8Against7merA11 <- function() {
  for (mirna in kMirnas) {
    if (mirna == "miR-7-23nt") {
      combined <- FALSE
      experiment <- "equilibrium2_nb"
    } else if (mirna == "miR-1") {
      combined <- FALSE
      buffer <- TRUE
      experiment <- "equilibrium"
    } else {
      combined <- TRUE
      buffer <- FALSE
      experiment <- "equilibrium"
    }
    SubfunctionCall(Check7merm8Against7merA1)
  }
}

CheckCenteredAgainst6merm8 <- function(mirna, experiment="equilibrium",
                                       n_constant=5, sitelist="centered11",
                                       buffer=FALSE, combined=TRUE) {
  kds <- SubfunctionCall(EquilPars)
  kd1 <- kds["6mer-m8_Kd", 2]/kds["11mer-m3.13_Kd", 2]
  kd2 <- kds["6mer-m8_Kd", 2]/kds["11mer-m4.14_Kd", 2]
  kd1_low <- kds["6mer-m8_Kd", 2]/kds["11mer-m3.13_Kd", 3]
  kd2_low <- kds["6mer-m8_Kd", 2]/kds["11mer-m4.14_Kd", 3]
  kd1_high <- kds["6mer-m8_Kd", 2]/kds["11mer-m3.13_Kd", 5]
  kd2_high <- kds["6mer-m8_Kd", 2]/kds["11mer-m4.14_Kd", 5]

  message(mirna)
  message(sprintf("m3.13:\t%s\t%s\t%s", kd1_low, kd1, kd1_high))
  message(sprintf("m4.14:\t%s\t%s\t%s", kd2_low, kd2, kd2_high))
}

CheckAllCenteredAgainst6merm8 <- function(sitelist="centered11",
                                          combinedtrue_all=FALSE,
                                          bufferfalse_all=FALSE,
                                          singleonly=TRUE) {
  for (mirna in kMirnas) {
    if (mirna == "miR-7-23nt") {
      combined <- FALSE
      buffer <- FALSE
      experiment <- "equilibrium2_nb"
    } else if (mirna == "miR-1") {
      combined <- FALSE
      buffer <- TRUE
      experiment <- "equilibrium"
    } else {
      combined <- TRUE
      buffer <- FALSE
      experiment <- "equilibrium"
    }
    if (combinedtrue_all) {
      combined <- TRUE
    }
    if (bufferfalse_all) {
      buffer <- FALSE
    }
    SubfunctionCall(CheckCenteredAgainst6merm8)
  }
}

CheckBaek7merm3.9Against6merm8 <- function(mirna, experiment="equilibrium",
                                           n_constant=5, sitelist="baek",
                                           buffer=TRUE, singleonly=FALSE,
                                           combined=TRUE) {
  kds <- SubfunctionCall(EquilPars)
  kd1 <- kds["6mer-m8_Kd", 2]/kds["7mer-m3.9", 2]
  message(sprintf("7mer-m3.8:\t%s", kd1))
}

CheckAllBaek7merm3.9Against6merm8 <- function() {
  for (mirna in kMirnas) {
    if (mirna == "miR-7-23nt") {
      combined <- FALSE
      buffer <- FALSE
      experiment <- "equilibrium2_nb"
    } else if (mirna == "miR-1") {
      combined <- FALSE
      buffer <- TRUE
      experiment <- "equilibrium"
    } else {
      combined <- TRUE
      buffer <- FALSE
      experiment <- "equilibrium"
    }
    SubfunctionCall(CheckBaek7merm3.9Against6merm8)
  }
}


CheckCDNST1 <- function(mirna, experiment="equilibrium", n_constant=5,
                        sitelist="baek", buffer=FALSE, combined=TRUE) {

  kds <- SubfunctionCall(EquilPars)
  kds_norm <- SubfunctionCall(EquilPars, sitelist="resubmissionfinal")

  if (mirna %in% c("miR-1", "miR-7-23nt")){
    kd1 <- 1/kds["CDNST.1_Kd", 2]  
  } else {
    kd1 <- 1/kds["CDNST 1_Kd", 2]  
  }
  if ("5mer-m2.6_Kd" %in% rownames(kds_norm)) {
    kd1_alt <- 1/kds_norm["5mer-m2.6_Kd", 2]
    message(kd1_alt)
  }

}

CheckAllCDNST1 <- function(sitelist="baek") {
  for (mirna in kMirnas) {
    if (mirna == "miR-7-23nt") {
      combined <- FALSE
      experiment <- "equilibrium2_nb"
    } else if (mirna == "miR-1") {
      combined <- FALSE
      buffer <- TRUE
      experiment <- "equilibrium"
    } else {
      combined <- TRUE
      buffer <- FALSE
      experiment <- "equilibrium"
    }
    SubfunctionCall(CheckCDNST1)
  }
}


# GetLinShiRepressionData <- function(best=FALSE, old=FALSE, threePseq=TRUE) {
#   if (best) str1 <- "best_"
#   else str1 <- "all_"
#   if (old) old <- "_old"
#   else old <- ""
#   if (threePseq) threePseq <- "_3pseq"
#   else threePseq <- ""
#   path <- paste0("RepressionData/Lin-Shi_transfection_data/", str1,
#                  "flanking_kds_and_repression", threePseq, old, ".txt")
#   fread(path, data.table=FALSE)
# }

GetSiteOccupancy <- function(mirna, experiment="equilibrium",  n_constant=5,
                             sitelist="resubmissionfinal", buffer=FALSE,
                             combined=TRUE) {
  pars.matrix <- SubfunctionCall(EquilPars)
  sXc <- SubfunctionCall(SitesXCounts)
  pars.model <- log10(pars.matrix$Mean)
  names(pars.model) <- rownames(pars.matrix)
  names(pars.model)[length(pars.model)] <- "AGO"
  names(pars.model)[length(pars.model) - 1] <- "bg"
  data <- GetDataEquil(sXc)
  l <- SubfunctionCall(GetInputEquil)
  A.dil.data <- sapply(colnames(data), as.numeric)
  model.points <- EquilSingleSiteModelFreq(pars.model, sXc, A.dil.data,
                                           combined=combined, addbg=FALSE)
  model.points
}

GetSiteOccupancy <- function(mirna, experiment="equilibrium",  n_constant=5,
                             sitelist="resubmissionfinal", buffer=FALSE,
                             combined=TRUE, singleonly=FALSE,
                             remove_sites=TRUE, addbg=FALSE) {
  pars.matrix <- SubfunctionCall(EquilPars)
  sXc <- SubfunctionCall(SitesXCounts)
  if (remove_sites ) {
    removed_sites <- GetRemovedSites(sXc)  
  } else {
    removed_sites <- c()
  }
  sXc <- sXc[setdiff(rownames(sXc), removed_sites), ]
  pars.matrix <- pars.matrix[setdiff(rownames(pars.matrix),
                                     paste0(removed_sites, "_Kd")), ]
  
  pars.model <- log10(pars.matrix$Mean)
  names(pars.model) <- rownames(pars.matrix)
  names(pars.model)[length(pars.model)] <- "AGO"
  names(pars.model)[length(pars.model) - 1] <- "bg"
  data <- GetDataEquil(sXc)
  l <- SubfunctionCall(GetInputEquil)
  A.dil.data <- sapply(colnames(data), as.numeric)
  model.points <- EquilSingleSiteModelFreq(pars.model, sXc, A.dil.data,
                                           combined=combined, addbg=addbg)
  model.points
}

GetSiteOccupancyPerSite <- function(mirna, experiment="equilibrium",  n_constant=5,
                             sitelist="resubmissionfinal", buffer=FALSE,
                             combined=TRUE) {
  pars.matrix <- SubfunctionCall(EquilPars)
  sXc <- SubfunctionCall(SitesXCounts)
  pars.model <- log10(pars.matrix$Mean)
  names(pars.model) <- rownames(pars.matrix)
  names(pars.model)[length(pars.model)] <- "AGO"
  names(pars.model)[length(pars.model) - 1] <- "bg"
  data <- GetDataEquil(sXc)
  l <- SubfunctionCall(GetInputEquil)
  A.dil.data <- sapply(colnames(data), as.numeric)
  model.points <- EquilSingleSiteModelPerSite(pars.model, sXc, A.dil.data,
                                           combined=combined, addbg=FALSE)
  model.points
}




# GetSiteOccupancy2 <- function(mirna, experiment="equilibrium", n_constant=5,
#                               sitelist="resubmissionfinal", combined=TRUE,
#                               global=FALSE, fixed=FALSE, buffer=FALSE) {
#   message("Site Occupancy")
#   message(mirna)
#   sXc <- SubfunctionCall(SitesXCounts)
#   pars.matrix <- SubfunctionCall(EquilPars)
#   if (global) {
#     pars.matrix <- rbind(pars.matrix[grep("Kd", rownames(pars.matrix)), ],
#                          pars.matrix[grep(mirna, rownames(pars.matrix)), ])
#   }

#   pars.model <- log10(pars.matrix$Mean)
#   names(pars.model) <- rownames(pars.matrix)
#   names(pars.model)[length(pars.model)] <- "AGO"
#   names(pars.model)[length(pars.model) - 1] <- "bg"
#   data <- GetDataEquil(sXc)
#   l <- SubfunctionCall(GetInputEquil)
#   # data <- data[, 1:5]
#   # Ago dilutio in the data:
#   A.dil.data <- sapply(colnames(data), as.numeric)
#   if (mirna %in% c("miR-7-23nt", "miR-7-24nt", "miR-7-25nt")) {
#     pars.matrix_Ago <- SubfunctionCall(EquilPars, fixed=TRUE, global=FALSE)
#     A.stock.measured <- pars.matrix_Ago[sprintf("AGO_%s", mirna), 1]
#   } else {
#     A.stock.measured <- kAgoStock[mirna, "equilibrium"]    
#   }
#   # Checking using the AGO concnetration from the model:
#   A.stock.measured <- pars.matrix[sprintf("AGO_%s", mirna), 1]
#   pars <- pars.model
#   model <- SubfunctionCall(EquilSingleSiteModelFreq, A.dil=A.dil.data, addbg=FALSE)
#   model
# }

GetRepressionMatrix <- function(mirna, sitelist, flanks=FALSE, bg_method=3, exrib=FALSE, new=FALSE, old=FALSE) {
  if (flanks) {
    flank_str = "_flank"
  } else {
    flank_str = ""
  }
  if (new) {
    new_str = "_new"
  } else if (old) {
    new_str = "_old"
  } else {
    new_str = ""
  }
  if (exrib) {
    exrib_str = "_exrib"
  } else {
    exrib_str = ""
  }
  file_path <- sprintf("/lab/solexa_bartel/mcgeary/AgoRBNS/%s/repression_hela_cs/lin%s_model_df/%s_%s%s%s.txt",
                                     mirna, flank_str, sitelist, bg_method, exrib_str, new_str)
  print(file_path)
  repression_matrix <- fread(file_path, sep="\t")
  saved_names <- colnames(repression_matrix)
  repression_matrix <- data.frame(repression_matrix, row.names=1)
  colnames(repression_matrix) <- saved_names[-1]
  repression_matrix
}

GetNumbersOfMrnasWithSites <- function(mirna, sitelist="resubmissionfinal",
                                       bg_method=3, new=FALSE, old=FALSE) {
  data_frame <- SubfunctionCall(GetRepressionMatrix)
  data_frame <- data_frame[, 2:(ncol(data_frame) - 1)]
  print(colSums(data_frame))
}


GetRepressionLinearModel <-  function(mirna, experiment="equilibrium",
                                      n_constant=5, sitelist="resubmissionfinal",
                                      buffer=FALSE, combined=TRUE, 
                                      bg_method=3, exrib=FALSE, new=FALSE, old=FALSE,
                                      intercept=FALSE, best_site=FALSE,
                                      single_site=FALSE, n_cutoff=-Inf,
                                      kd_cutoff=Inf) {
  print(bg_method)
  print(exrib)
  print(new)
  print(old)
  print(intercept)
  print(best_site)
  print(single_site)
  print(n_cutoff)
  print(kd_cutoff)
  repression_matrix <- SubfunctionCall(GetRepressionMatrix)
  # saved_names <- colnames(repression_matrix)
  # repression_matrix <- data.frame(repression_matrix, row.names=1)
  # colnames(repression_matrix) <- saved_names[-1]
  if (mirna == "miR-7-23nt") {
    mirna_trim <- "miR-7"
  } else {
    mirna_trim <- mirna
  }
  pars <- SubfunctionCall(EquilPars)
  if (single_site) {
    repression_matrix <- repression_matrix[which(rowSums(repression_matrix[,-ncol(repression_matrix)]) == 1),]
  }
  rownames(pars) <- gsub("_Kd", "", rownames(pars))
  sites <- grep(sprintf("(?:%s|None)", mirna_trim), rownames(pars), value=TRUE,
                invert=TRUE)
  sites <- colnames(repression_matrix[1:(ncol(repression_matrix) - 1)])
  site_occurence <- repression_matrix[,sites, drop=FALSE]

  n_sites <- colSums(site_occurence)

  kds <- pars[sites, ]$Mean

  sites <- sites[which((n_sites >= n_cutoff & kds <= kd_cutoff) | sites %in% kSeedSites)]
  n_sites <- n_sites[sites]

  if (length(sites) == 0) {
    return()
  }
  order_sites <- seq(length(sites))
  names(order_sites) <- sites
  site_occurence <- site_occurence[, sites]
  fc <- repression_matrix[["fc"]]/log(2)
  nosite_fc <- mean(fc[which(rowSums(site_occurence) == 0)])
  fc <- fc - nosite_fc
  if (intercept) {
    model <- lm(fc ~ ., data=site_occurence)    
  } else {
    model <- lm(fc ~ . - 1, data=site_occurence)

  }
  CIs <- confint(model)
  vals <-  model$fitted.values
  out <- summary(model)$coefficients[, 1:2, drop=FALSE]
  summary_global <<- summary(model)
  CIs <- CIs[which(n_sites > 0),]
  n_sites <- n_sites[which(n_sites > 0)]
  if (intercept) {
    out <- cbind(out[-1, , drop=FALSE], n_sites, CIs[-1, ])    
  } else {
    out <- cbind(out[, , drop=FALSE], n_sites, CIs)    
  }
  rownames(out) <- gsub("`", "", rownames(out))
  return(out)
}

GetRepFlankLinearModel <- function(experiment="equilibrium", n_constant=5,
                                   sitelist="resubmissionfinal", buffer=FALSE,
                                   combined=TRUE, bg_method=3, exrib=FALSE,
                                   new=FALSE, old=FALSE, intercept=FALSE,
                                   best_site=FALSE, single_site=FALSE,
                                   n_cutoff=-Inf, kd_cutoff=Inf) {
  site_occs_all <- matrix(nrow=0, ncol=0)
  flank_occs_all <- matrix(nrow=0, ncol=256)
  fc_all <- c()
  for (mirna in kMirnas) {
    if (mirna == "miR-1") {
      buffer <- TRUE
      combined <- FALSE
      experiment <- "equilibrium"
    } else if (mirna == "miR-7-23nt") {
      buffer <- FALSE
      combined <- FALSE
      experiment <- "equilibrium2_nb"
    } else {
      buffer <- FALSE
      combined <- TRUE
      experiment <- "equilibrium"
    }
    rep_data <- SubfunctionCall(GetRepressionMatrix, flanks=TRUE)
    rownames(rep_data) <- paste0(mirna, "_", rownames(rep_data))

    site_occs <- rep_data[, 1:(ncol(rep_data) - 256 - 1)]
    colnames(site_occs) <- paste0(mirna, "_", colnames(site_occs))

    flank_occs <- rep_data[,(ncol(rep_data) - 256): (ncol(rep_data) - 1)]
    fc <- rep_data[["fc"]]/log(2)
    # Define the rows to keep
    site_inds <- which(rowSums(site_occs) > 0)
    # Define sites to remove due to having zero sites:
    site_inds_col <- which(colSums(site_occs) > 0)
    # Find the average nosite for this miRNA, then substract it.
    nosite_fc <- mean(fc[-site_inds])
    # Remove the rows with no sites.
    site_occs <- site_occs[site_inds, site_inds_col]
    flank_occs <- flank_occs[site_inds,]
    fc <- fc[site_inds] - nosite_fc
    # Update the site_occs_table
    site_occs_all_top_right <- matrix(0, nrow=nrow(site_occs_all),
                                      ncol=ncol(site_occs),
                                      dimnames=list(rownames(site_occs_all),
                                                    colnames(site_occs)))
    site_occs_all_bottom_left <- matrix(0, nrow=nrow(site_occs),
                                      ncol=ncol(site_occs_all),
                                      dimnames=list(rownames(site_occs),
                                                    colnames(site_occs_all)))
    site_occs_all_top <- cbind(site_occs_all, site_occs_all_top_right)
    site_occs_all_bottom <- cbind(site_occs_all_bottom_left, site_occs)
    site_occs_all <- rbind(site_occs_all_top, site_occs_all_bottom)

    fc_all <- c(fc_all, fc)

    site_occs_all_bottom
    # Update the master flank_occs_all
    flank_occs_all <- rbind(flank_occs_all, flank_occs)
  }
  flank_occs_all_image <- as.matrix(flank_occs_all)
  flank_occs_all_image[flank_occs_all_image >=1 ] <- 1
  site_occs_all_image <- as.matrix(site_occs_all)
  site_occs_all_image[site_occs_all_image >=1 ] <- 1
  rep_master <- cbind(site_occs_all, flank_occs_all[, -1])
  rep_master <<- rep_master
  model <- lm(fc_all ~ . - 1, data=rep_master)
  CIs <- confint(model)
  out <- summary(model)$coefficients[, 1:2, drop=FALSE]
  n_sites <- colSums(rep_master)

  out <- cbind(out[, , drop=FALSE], n_sites[which(n_sites > 0)], CIs[which(n_sites > 0), ])
  out <- out[grep("_", rownames(out), invert=TRUE), ]
  out
}

GetAUFreqTable <- function(mirna, sitelist="resubmissionfinal", bg_method=3,
                           exrib=FALSE, new=TRUE) {
  if (new) {
    new_str = "_new"
  } else if (old) {
    new_str = "_old"
  } else {
    new_str = ""
  }
  if (exrib){
    exrib_str = "_exrib"
  } else {
    exrib_str = ""
  }
  extension <- sprintf("_%s%s%s", bg_method, exrib_str, new_str)
  path <- SubfunctionCall(GetAnalysisPath,
                                 experiment="repression_hela_cs",
                                 analysis_type="AUfreq_df", condition=sitelist,
                                 ext=extension)
  read.table(path, sep="\t", header=TRUE, row.names=1)
}

GetFlankTable <- function(mirna, sitelist="resubmissionfinal", bg_method=3, exrib=TRUE, new=TRUE) {
  if (new) {
    new_str = "_new"
  } else if (old) {
    new_str = "_old"
  } else {
    new_str = ""
  }
  if (exrib){
    exrib_str = "_exrib"
  } else {
    exrib_str = ""
  }
  extension <- sprintf("_%s%s%s", bg_method, exrib_str, new_str)
  path <- SubfunctionCall(GetAnalysisPath,
                                 experiment="repression_hela_cs",
                                 analysis_type="site_flank_count_df", condition=sitelist,
                                 ext=extension)
  read.table(path, sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)
}

GetAllFlankTables <- function(sitelist="resubmissionfinal", bg_method=3, exrib=TRUE,
                             new=TRUE, n_cutoff=20, kd_cutoff=0.1) {
  output <- matrix(nrow=0, ncol=258)
  colnames(output) <- c("mirna", "site", GetKmerList(4))
  for (mirna in kMirnas) {
    flank_count_df <- SubfunctionCall(GetFlankTable)
    if (mirna == "miR-1") {
      combined <- FALSE
      buffer <- TRUE
      experiment <- "equilibrium"
    } else if (mirna == "miR-7-23nt") {
      combined <- FALSE
      buffer <- FALSE
      experiment <- "equilibrium2_nb"
    } else {
      combined <- TRUE
      buffer <- FALSE
      experiment <- "equilibrium"
    }
    # Get the Kds for all the sites (not including the 'None' site), order them
    # based on their value, and then remove the "_Kd" at the end of each Kd
    # string.
    # '_Kd' at the end of the string, and 
    kds <- SubfunctionCall(EquilPars)
    kds <- kds[1:(nrow(kds) - 3), ]
    kds <- kds[order(kds$Mean),]
    rownames(kds) <- gsub("_Kd", "", rownames(kds))
    # Order the flank rows according to the Kd values.
    flank_count_df <- flank_count_df[rownames(kds), ]
    # Define the condition for which a site is kept.
    n_bool <- rowSums(flank_count_df) >= n_cutoff
    kd_bool <- kds$Mean <= kd_cutoff
    seed_bool <- rownames(flank_count_df) %in% kSeedSites
    full_bool <- which(seed_bool | (n_bool & kd_bool))
    flank_count_df <- flank_count_df[full_bool, ]
    flank_count_df <- cbind(miRNA=rep(mirna, nrow(flank_count_df)),
                            Site=rownames(flank_count_df),
                            flank_count_df,
                            stringsAsFactors=FALSE)
    rownames(flank_count_df) <- c()
    output <- rbind(output, flank_count_df)
  }
  output
}




MakeFinalSiteFlankTable <- function(sitelist="resubmissionfinal", bg_method=3,
                                    exrib=TRUE, new=TRUE, n_cutoff=20,
                                    kd_cutoff=0.1) {
  # Get the input flanking count dataframe.
  site_flank_df <- SubfunctionCall(GetAllFlankTables)
  # Extract the counts (columnes 3–258, because column 1 is the mirna and
  # column 2 is the site type).
  site_flank_norm_df <- site_flank_df[, 3:258]/rowSums(site_flank_df[, 3:258])
  # Pre-allocate the final matrix.
  out_final <- as.data.frame(
    matrix(0, nrow=nrow(site_flank_norm_df), ncol=32,
           dimnames=list(seq(nrow(out_final)),
                         c(paste0(GetKmerList(2), "_5p"),
                           paste0(GetKmerList(2), "_3p")))
  ))
  # Iterate over the kmers, sum the relevant columns in the normalized flank
  # matrix, and assign it to the correct column in 'out-final'.
  for (kmer in GetKmerList(2)) {
    kmer_5p <- rowSums(
      site_flank_norm_df[, grep(paste0("^", kmer),
                                colnames(site_flank_norm_df),
                                perl=TRUE)]
    )
    kmer_3p <- rowSums(
      site_flank_norm_df[, grep(paste0(kmer, "$"),
                                colnames(site_flank_norm_df),
                                perl=TRUE)]
    )    
    out_final[paste0(kmer, "_5p")] <- kmer_5p
    out_final[paste0(kmer, "_3p")] <- kmer_3p
  }
  # Output the final dataframe with the miRNA and site-type column added.
  output <- cbind(site_flank_df[, 1:2], out_final)
  write.table(output, file="2017_Paper/Referee_Figures/SiteFlankTable.txt",
              sep="\t", row.names=FALSE)
}

GetAllThreePrimeRepUTRs <-  function(experiment="equilbrium", n_constant=5,
                                      sitelist="paper_and_8mer_3p",
                                      combined=TRUE, best_site=FALSE,
                                      single_site=FALSE, n_cutoff=-Inf,
                                      kd_cutoff=Inf, xpos=20, ypos=20) {
  all_reps <- list()
  all_reps_3p <- list()
  all_reps_3p_10_11 <- list()
  ind_out <- 1
  mirnas_use <- c("miR-155", "miR-124", "lsy-6")
  # mirnas_use <- c("miR-124")

  k3PSites <- c(k3PSites, "8mer-m11.18", "8mer-m12.19", "8mer-m13.20",
                "8mer-m14.21", "8mer-m15.22", "8mer-m16.23")
  fc_3p_crude_all <- c()
  fc_3p_10_11_crude_all <- c()
  fc_nosite_all <- c()
  fc_nosite_each <- list()
  for (mirna in mirnas_use) {
    repression_matrix <- fread(sprintf("/lab/solexa_bartel/mcgeary/AgoRBNS/%s/repression_hela_cs/lin_model_df/%s.txt",
                                       mirna, sitelist), sep="\t")
    if (mirna == "miR-124") {
      repression_matrix <<- repression_matrix    
    }
    names_saved <- colnames(repression_matrix)
    colnames(repression_matrix) <- gsub("^X(.*)$", replace="\\1", colnames(repression_matrix), perl=TRUE)
    saved_names <- colnames(repression_matrix)
    repression_matrix <- data.frame(repression_matrix)
    colnames(repression_matrix) <- names_saved
    inds_3p <- which(colnames(repression_matrix) %in% k3PSites)
    inds_3p_10_11 <- which(colnames(repression_matrix) %in% k3PSites[grep("(?:8|9)mer", k3PSites, perl=TRUE, invert=TRUE)])
    row_inds_3p <- which(colSums(repression_matrix[, inds_3p]) > 0)
    row_inds_3p_10_11 <- which(colSums(repression_matrix[, inds_3p_10_11]) > 0)
    all_3p <- rowSums(repression_matrix[, inds_3p])
    all_3p_10_11 <- rowSums(repression_matrix[, inds_3p_10_11])

    repression_matrix_3p <- cbind(repression_matrix[, c(-inds_3p, 
                                                         -ncol(repression_matrix))],
                                   all_3p,
                                   repression_matrix[, ncol(repression_matrix)])
    repression_matrix_3p_10_11 <- cbind(repression_matrix[, c(-inds_3p_10_11, 
                                                         -ncol(repression_matrix))],
                                   all_3p_10_11,
                                   repression_matrix[, ncol(repression_matrix)])

    colnames(repression_matrix_3p)[ncol(repression_matrix_3p)] <- "fc"
    colnames(repression_matrix_3p_10_11)[ncol(repression_matrix_3p_10_11)] <- "fc"
    sites <- colnames(repression_matrix[, 2:(ncol(repression_matrix) - 1)])
    site_occurence <- repression_matrix[,sites]
    site_occurence_3p <- repression_matrix_3p[, 2:(ncol(repression_matrix_3p) - 1)]
    site_occurence_3p_10_11 <- repression_matrix_3p_10_11[, 2:(ncol(repression_matrix_3p_10_11) - 1)]
    fc <- repression_matrix[["fc"]]/log(2)
    names(fc) <- repression_matrix[, 1]
    fc_3p <- repression_matrix_3p[["fc"]]/log(2)
    fc_3p_10_11 <- repression_matrix_3p_10_11[["fc"]]/log(2)
    nosite_fc <- mean(fc[which(rowSums(site_occurence) == 0)])
    nosite_fc_vec <- fc[which(rowSums(site_occurence) == 0)]


    if (mirna == "miR-124") {
      site_occurence <<- site_occurence
    }

    fc_nosite_all <- c(fc_nosite_all, nosite_fc_vec)
    fc_3p_crude <- repression_matrix[["fc"]][row_inds_3p]
    fc_3p_10_11_crude <- repression_matrix[["fc"]][row_inds_3p_10_11]

    fc_3p_crude_all <- c(fc_3p_crude_all, fc_3p_crude)
    fc_3p_10_11_crude_all <- c(fc_3p_10_11_crude_all, fc_3p_10_11_crude)

    num_nosite <- length(which(rowSums(site_occurence) == 0))
    num_nosite_3p <- length(which(rowSums(site_occurence_3p) == 0))
    num_nosite_3p_10_11 <- length(which(rowSums(site_occurence_3p_10_11) == 0))


    fc_nosite_each[[ind_out]] <- nosite_fc_vec
    n_sites <- colSums(site_occurence)
    n_sites_3p <- colSums(site_occurence_3p)
    n_sites_3p_10_11 <- colSums(site_occurence_3p_10_11)
    n_sites <- n_sites[which(n_sites > 0)]
    n_sites_3p <- n_sites_3p[which(n_sites_3p > 0)]
    n_sites_3p_10_11 <- n_sites_3p_10_11[which(n_sites_3p_10_11 > 0)]

    fc <- fc - nosite_fc
    model <- lm(fc ~ ., data=site_occurence)
    model_3p <- lm(fc ~ ., data=site_occurence_3p)
    model_3p_10_11 <- lm(fc ~ ., data=site_occurence_3p_10_11)
    vals <-  model$fitted.values
    vals_3p <-  model_3p$fitted.values
    vals_3p_10_11 <-  model_3p_10_11$fitted.values
    out <- summary(model)$coefficients[, 1:2, drop=FALSE]
    out <- cbind(out, c(num_nosite, n_sites))
    out_3p <- summary(model_3p)$coefficients[, 1:2, drop=FALSE]
    out_3p <- cbind(out_3p, c(num_nosite_3p, n_sites_3p))

    out_3p_10_11 <- summary(model_3p_10_11)$coefficients[, 1:2, drop=FALSE]
    out_3p_10_11 <- cbind(out_3p_10_11, c(num_nosite_3p_10_11, n_sites_3p_10_11))

    out[, 1] <- out[, 1] - out[1, 1]
    out[, 2] <- sqrt(out[, 2]^2 + out[1, 2]^2)
    out <- cbind(out[-1, , drop=FALSE], n_sites)
    out_3p[1,] <- out_3p[1, ] + out_3p[1, 1]
    out_3p[2, ] <- sqrt(out_3p[2, ]^2 + out_3p[1, 2]^2)
    out_3p <- cbind(out_3p[-1, , drop=FALSE], n_sites_3p)

    out_3p_10_11[1,] <- out_3p_10_11[1, ] + out_3p_10_11[1, 1]
    out_3p_10_11[2, ] <- sqrt(out_3p_10_11[2, ]^2 + out_3p_10_11[1, 2]^2)
    out_3p_10_11 <- cbind(out_3p_10_11[-1, , drop=FALSE], n_sites_3p_10_11)

    message("dim of final output vector from linear repression fucntion")
    rownames(out) <- gsub("`", "", rownames(out))
    rownames(out_3p) <- gsub("`", "", rownames(out_3p))
    rownames(out_3p_10_11) <- gsub("`", "", rownames(out_3p_10_11))
    all_reps[[ind_out]] <- out
    all_reps_3p[[ind_out]] <- out_3p
    all_reps_3p_10_11[[ind_out]] <- out_3p_10_11
    ind_out <- ind_out + 1
  }

  ecdf_rep <- ecdf(fc_nosite_all)
  ecdf_rep_3p <- ecdf(fc_3p_crude_all)
  ecdf_rep_3p_10_11 <- ecdf(fc_3p_10_11_crude_all)
  x <- seq(-2, 2, length.out=1000)

  t_test_3p <- t.test(fc_3p_crude_all, fc_nosite_all)

  xmin <- -2
  xmax <- 2
  ymin <- 0
  ymax <- 1

  t_test_3p_10_11 <- t.test(fc_3p_10_11_crude_all, fc_nosite_all)
  t_test_out <<- t_test_3p_10_11
  dev.new(xpos=xpos, ypos=ypos, height=5, width=5)
  par(kPlotParameters)
  BlankPlot()
  AddLinearAxis(1, tick.space=0.5, label.space=1, label="fold-change")
  AddLinearAxis(2, tick.space=0.1, label.space=0.2, label="ecdf")

  xy <- GetPlotFractionalCoords(fx=0.025, fy=0.95)

  text(xy[1], xy[2], labels=sitelist, adj=0)
  xy <- GetPlotFractionalCoords(fx=0.95, fy=0.95)
  AddPValueToPlot(t_test_3p$p.value, xy[1], xy[2], col="red", adj=1)
  xy <- GetPlotFractionalCoords(fx=0.95, fy=0.90)

  AddPValueToPlot(t_test_3p_10_11$p.value, xy[1], xy[2], col="blue", adj=1)


  lines(x, ecdf_rep(x))
  lines(x, ecdf_rep_3p(x), col="red")
  lines(x, ecdf_rep_3p_10_11(x), col="blue")

  nosite_each_m124 <<- fc_nosite_each[[2]]

  return()
}


# Used in "Insights into miRNA targerting".
CheckSiteOccurences <- function(n_constant=5, sitelist="resubmissionfinal") {
  sites_list <- list()
  for (mirna in kMirnas) {
    if (mirna == "miR-1") {
      buffer <- TRUE
    } else {
      buffer <- FALSE
    }
    if (mirna == "miR-7-23nt") {
      experiment <- "equilibrium2_nb"
    } else {
      experiment <- "equilibrium"
    }
    sites_list[[mirna]] <- rownames(SubfunctionCall(SitesXCounts))
  }
  sites_all <- unique(unlist(sites_list))
  site_counts <- rep(0, length(sites_all))
  site_names <- rep("", length(sites_all))

  names(site_counts) <- sites_all
  names(site_names) <- sites_all
  for (mirna in kMirnas) {
    sites_mirna <- sites_list[[mirna]]
    site_counts[sites_mirna] <- site_counts[sites_mirna] + 1
    site_names[sites_mirna] <- sprintf("%s %s", site_names[sites_mirna], mirna)
  }
  out <- data.frame(counts=site_counts, mirnas=site_names)
  out <- out[order(-out$count), ]
  out
}

CheckAllCanonicalOccupancies <- function(experiment="equilibrium",
                                         n_constant=5, sitelist="resubmissionfinal", 
                                         combined=TRUE, buffer=FALSE,
                                         singleonly=FALSE) {
  sapply(kMirnas, function(mirna) {
    message(sprintf("%s__________________________________________", mirna))
    if (mirna %in% c("miR-1", "miR-7-23nt")) {
      combined <- FALSE
    }
    if (mirna == "miR-1") {
      buffer <- TRUE
    }
    if (mirna == "miR-7-23nt") {
      experiment <- "equilibrium2_nb"
    }
    occs <- SubfunctionCall(GetSiteOccupancy)
    c_occs <- colSums(occs[kSeedSites, ])
    c_occs_2 <- colSums(occs[kCanonicalSites, ])

    nc_occs <- colSums(occs[setdiff(rownames(occs), c(kSeedSites, "None")), ])
    message("Canonical site occupancy:")
    print(c_occs)
    message("More stringent canonical site occupancy:")
    print(c_occs_2)
    message("Noncanonical site occupancy:")
    print(nc_occs)
    message("Canonical/noncanonical ratio:")
    print(c_occs/nc_occs)
    print("No-site occupancy:")
    print(colSums(occs["None", , drop=FALSE]))
    print(colSums(occs))
  })
}

CheckAllNonCanonicalLengths <- function(experiment="equilibrium",
                                         n_constant=5, sitelist="resubmissionfinal", 
                                         combined=TRUE, buffer=FALSE,
                                         singleonly=FALSE) {
  sapply(kMirnas, function(mirna) {
    message(sprintf("%s__________________________________________", mirna))
    if (mirna %in% c("miR-1", "miR-7-23nt")) {
      combined <- FALSE
    }
    if (mirna == "miR-1") {
      buffer <- TRUE
    }
    if (mirna == "miR-7-23nt") {
      experiment <- "equilibrium2_nb"
    }
    sXc <- SubfunctionCall(SitesXCounts)
    kds <- SubfunctionCall(EquilPars)
    kds <- kds[1:nrow(sXc), ]
    kds <- kds[paste0(rownames(sXc), "_Kd"), ]$Mean
    names(kds) <- rownames(sXc)
    canon_kds <- kds[c("8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-m8", "6mer-A1")]
    names_nc_kds <- setdiff(names(kds), c("8mer", "7mer-m8", "7mer-A1",
                                                 "6mer", "6mer-m8", "6mer-A1"))
    noncanon_kds <- kds[names_nc_kds]
    for (i in seq(length(noncanon_kds))) {
      for (j in seq(length(canon_kds))) {
        nc_site <- names(noncanon_kds)[i]
        nc_kd <- noncanon_kds[i]
        c_site <- names(canon_kds)[j]
        c_kd <- canon_kds[j]
        nc_seq <- GetSiteSeq(mirna, nc_site)
        c_seq <- GetSiteSeq(mirna, c_site)
        if ((nc_kd <= c_kd) & nchar(nc_seq) <= nchar(c_seq)) {
          print("BETTER")
          print(sprintf("%s is better than %s.", nc_site, c_site))
        }
      }
    }
  })
}


CountAllFlankingSites <- function(mirna, experiment="equilibrium",
                                  n_constant=5, sitelist="resubmissionfinal",
                                  combined=TRUE, buffer=FALSE) {
  print(mirna)
  if (mirna == "miR-1") {
    combined <- FALSE
    buffer <- TRUE
  }
  sXc <- SubfunctionCall(SitesXCounts)
  if (mirna == "miR-1") {
    combined <- TRUE
  }
  sites <- rownames(sXc)
  flank_min <- 256
  for (site in sites[-length(sites)]) {
    fXc <- SubfunctionCall(SiteFlanksXCounts)
    if (nrow(fXc) <= flank_min) {
      flank_min <- nrow(fXc)
      print("new minimum.")
    }
  }
  print(flank_min)
}


CalculateFlowthroughFractions <- function(f_nc1, f_nc2, f_ft, f_i) {
  # Cost function used in this optimization.

  CostFunction <- function(p) {
    # Make vector of weights of each of the three fractions, which corresdpond
    # to the flowthrough, top nitrocellulose and second nitrocellulose,
    # respectively.
    w <- Norm(exp(c(0, p)))
    data <- f_i
    # Matrix multiplication that gives the total input reads.
    model <- w %*% t(cbind(f_ft, f_nc1, f_nc2))
    sum((log(data) - log(model))^2)
  }
  # Optimize using the cost function.
  opt <- optim(c(0, 0), CostFunction)
  # Return normalized vector of parameters from optimization
  Norm(c(1, exp(opt$par)))
}


ReporterCounts <- function(mirna, experiment, condition, rep=NULL, tpm=FALSE,
                           pseudo=FALSE, mm=FALSE) {
  if (mirna == "miR-7-23nt") {
    mirna <- "miR-7"
  }
  analysis_type <- "counts"
  if (length(rep) != 0) {
    ext <- sprintf(",%s", rep)  
  } else {
    ext <- ""
  }
  if (mm) ext <- sprintf("%s_mm", ext)
  path <- SubfunctionCall(GetAnalysisPath)
  print(path)
  data <- read.table(path, stringsAsFactors=FALSE, sep="\t", header=FALSE,
                     colClasses=c("character", "character",
                                  rep("numeric", 184)))
  mirna_site <- paste(data[, 1], data[, 2], sep="_")
  data_full <- as.matrix(data[, 3:186])
  rownames(data_full) <- mirna_site
  if (pseudo) {
    data_full <- data_full + 1
  }
  if (tpm) {
    data_full <- 1e6*Norm(data_full)
  }
  data_full

}

GetBackGroundReporterCounts <- function(mirna, experiment, rep, tpm=TRUE) {
  bg_mirnas <- setdiff(kMirnas, mirna)

  counts_initial <- SubfunctionCall(ReporterCounts, condition="duplex",
                                    tpm=TRUE)

  bg_initial <- SubfunctionCall(ReporterCounts, condition="no_duplex",
                                    tpm=TRUE)

  counts_final <- counts_initial
  # Get the indeces of the mirna for which the background is wanted
  miR_inds <- grep(paste0(mirna, "_"), rownames(counts_initial))
  # Make a list of the five matrices of the expression of these variants from
  # the other five miRNAs.
  bg_sites <- lapply(bg_mirnas, function(mir_bg) {
    counts <- SubfunctionCall(ReporterCounts, condition="duplex", mirna=mir_bg,
                              tpm=TRUE)
    counts[miR_inds, ]
  })
  # Get the average value matrix across the list.
  sums <- Reduce("+", bg_sites)/length(bg_mirnas)
  # Replace those rows in the original counts with the sum from that list.
  counts_final[miR_inds, ] <- sums
  # For the other five mirnas, repeat the process, getting the average value
  # across the four miRNAs that aren't the mirna of interest nor the miRNA for
  # which 
  for (bg_mirna in bg_mirnas) {
    bg_mirnas_2 <- setdiff(bg_mirnas, bg_mirna)
    if (bg_mirna == "miR-7-23nt") {
      bg_mirna <- "miR-7"
    }

    miR_inds <- grep(paste0(bg_mirna, "_"), rownames(counts_initial))
    bg_sites <- lapply(bg_mirnas_2, function(mir_bg) {
      counts <- SubfunctionCall(ReporterCounts, condition="duplex",
                                mirna=mir_bg, tpm=TRUE)
      counts[miR_inds, ]
    })
    # Get the average value matrix across the list.
    sums <- Reduce("+", bg_sites)/length(bg_mirnas_2)
    # Replace those rows in the original counts with the sum from that list.
    counts_final[miR_inds, ] <- sums
  }
  counts_final
}

MakeReporterMatrix <- function(rep, experiment="twist_reporter_assay_v2", mm=FALSE,
                               output=FALSE) {
  counts <- SubfunctionCall(ReporterCounts, mirna="miR-1", condition="duplex")
  rnames_1 <- rep(rownames(counts), each=ncol(counts))
  rnames_2 <- rep(1:ncol(counts), nrow(counts))
  rnames <- paste(rnames_1, rnames_2, sep="_")
  counts_out <- matrix(0, nrow=nrow(counts)*ncol(counts), ncol=7)
  rownames(counts_out) <- rnames
  colnames(counts_out) <- c(kMirnas, "no_duplex")
  colnames(counts_out)[6] <- "miR-7"
  counts_out[, 1] <- c(t(counts))
  for (i_col in 2:6) {
    mirna <- colnames(counts_out)[i_col]
    counts <- SubfunctionCall(ReporterCounts, mirna=mirna, condition="duplex")
    counts_out[, i_col] <- c(t(counts))
  }
  counts_out[, 7] <- c(SubfunctionCall(ReporterCounts, mirna="miR-1", condition="no_duplex"))
  print(head(counts_out))
  if (output) {
    if (mm) mm_str <- "_mm"
    else    mm_str <- ""
    path <- sprintf("ReporterScreen/%s_rep%s_counts%s.txt", experiment, rep, mm_str)
    print(path)
    write.table(counts_out, file=path, row.names=TRUE)
    print("wrote file")
  } else {
    counts_out
  }
}


GetReporterMatrix <- function(experiment, rep, mm) {
  if (mm)  mm_str <- "_mm"
  else          mm_str <- ""

  read.table(sprintf("ReporterScreen/%s_rep%s_counts%s.txt",
                                    experiment, rep, mm_str))
}



MakeReporterMatrixV2test <- function(experiment="twist_reporter_assay_v2_test",
                                     mm=FALSE, output=FALSE) {
  counts <- SubfunctionCall(ReporterCounts, mirna="miR-1", condition="duplex_primer1")
  rnames_1 <- rep(rownames(counts), each=ncol(counts))
  rnames_2 <- rep(1:ncol(counts), nrow(counts))
  rnames <- paste(rnames_1, rnames_2, sep="_")
  counts_out <- matrix(0, nrow=nrow(counts)*ncol(counts), ncol=4)
  rownames(counts_out) <- rnames
  colnames(counts_out) <- c("miR-1_1", "miR-124_1", "miR-1_2", "miR-124_2")
  counts_out[, 1] <- c(t(counts))
  for (i_col in 2:4) {
    col_name <- colnames(counts_out)[i_col]
    col_name_split <- unlist(strsplit(col_name, split="_"))
    mirna <- col_name_split[1]
    condition <- sprintf("duplex_primer%s", col_name_split[2])
    counts <- SubfunctionCall(ReporterCounts, mirna=mirna, condition=condition)
    counts_out[, i_col] <- c(t(counts))
  }
  print(head(counts_out))
  if (output) {
    if (mm) mm_str <- "_mm"
    else    mm_str <- ""
    path <- sprintf("ReporterScreen/v2test_counts%s.txt", mm_str)
    print(path)
    write.table(counts_out, file=path, row.names=TRUE)
    print("wrote file")
  } else {
    counts_out
  }
}



# Used to get the numbers reported in the legend of fig. S2, A.
CheckTwoSitesDontPass <- function() {

  gray_sites <- c(`miR-155`="AACGAGG_Kd", `miR-124`="TCACCCGC_Kd")
  gray_sites <<- gray_sites
  for (mirna in c("miR-155", "miR-124")) {
    if (mirna == "miR-124") compcorrect <- TRUE
    else                    compcorrect <- FALSE
    pars.matrix <- EquilPars(mirna, experiment="equilibrium", n_constant=5,
                   sitelist="resubmissionfinal", combined=TRUE, singleonly=TRUE,
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

CountAllNonCompetitorSites <- function() {
  int_sites <- 0
  int_sites_corrected <- 0
  for (mirna in kMirnas) {
    if (mirna == "miR-7-23nt") experiment <- "equilibrium2_nb"
    else                       experiment <- "equilibrium"
    if (mirna == "miR-1") buffer <- TRUE
    else                  buffer <- FALSE
    sXc <- SubfunctionCall(SitesXCounts)
    int_exsites <- length(SubfunctionCall(GetRemovedSites))
    int_sites <- int_sites + nrow(sXc) - 1
    int_sites_corrected <- int_sites_corrected + nrow(sXc) - 1 - int_exsites
  }
  print(int_sites)
  print(int_sites_corrected)
}

LoadBeckerEtAlData <- function(mirna, with_names=TRUE, alt=FALSE, no_mm=FALSE) {
  # It seems that the miR-21 data has 
  if (with_names) str.with_names <- "_with_names"
  else            str.with_names <- ""
  if (alt) str.alt <- "_alt"
  else     str.alt <- ""
  if (no_mm)  str.no_mm <- "_no_mm"
  else        str.no_mm <- ""
  kPath <- sprintf("Becker,Ober-Reynolds_et_al_data/%s_v2%s%s%s.txt", mirna,
                   str.with_names, str.alt, str.no_mm)
  data <- unique(read.table(kPath, header=TRUE, stringsAsFactors=FALSE,
                     sep="\t"))
  if (with_names) {
    colnames(data) <- c("Sequence", "Site names", "Mismatch site names",
                        "Kd_median (pM)", "Lower_CI (pM)", "Upper_CI",
                        "rel_kon", "R2_kon", "n", "k_cleave (1/s)",
                        "R2_kcleave", "dG (kcal/mol)")
    start_col <- 4
  } else {
    colnames(data) <- c("Sequence", "Kd_median (pM)", "Lower_CI (pM)",
                        "Upper_CI", "rel_kon", "R2_kon", "n", "k_cleave (1/s)",
                        "R2_kcleave", "dG (kcal/mol)")
    start_col <- 2
  }
  sequences <- data[, 1]
  sequences_unique <- unique(sequences)

  if (mirna == "let-7a" & !with_names) {
    data <- data[1:(nrow(data) - 1), ]
  }
  dups <- which(duplicated(sequences))
  data_full <- data
  dup_matrix <- data[which(data[, 1] %in% unique(sequences[dups])), ]
  print(unique(sequences[dups])[1:10])
  print(sprintf("The %s data table has %s sequences repeated", mirna,
               length(unique(sequences[dups]))))
  message("rows with blank Kd values:")
  print(which(data[, start_col] == ""))

  for (dup_seq in unique(sequences[dups])) {
    data_rows <- which(data[, 1] == dup_seq)
    data_temp <- data[data_rows, ]
    if (data_temp[1, start_col] != "<10" & mirna == "miR-21") {
      print("different duplicate issue")
      print(data_temp)
      break
    }
    if (data_temp[1, start_col] == "<10" & data_temp[2, start_col] == "") {
      print(data_temp)
    }
    for (col_i in start_col:ncol(data_temp)) {
      data_temp_temp <- data_temp[, col_i]
      if ((is.na(data_temp_temp) ||
           data_temp_temp[1] == "<10" ||
           data_temp_temp[1] == "<-15.609") &&
          (!is.na(data_temp_temp[2]) &&
           data_temp_temp[2] != "" &&
           data_temp_temp[2] != "<10")) {
        data[data_rows[1], col_i] <- data_temp_temp[2]
      }
    }
    # Remove the duplcate rows
    data <- data[-data_rows[-1], ]
  }
  sequences <- data[, 1]
  rownames(data) <- sequences
  data <- data[, -1]
  for (col_i in (start_col - 1):(start_col + 1)) {
    row_i <- which(data[, col_i] == "<10" | data[, col_i] == " <10")
    data[row_i, col_i] <- 10
    row_i <- which(data[, col_i] == ">10000")
    data[row_i, col_i] <- 10000
    row_i <- which(data[, col_i] == "")
    data[row_i, col_i] <- 10000
    data[row_i, col_i] <- NaN
  }
  row_i <- which(data[, start_col + 5] == "<0.0001")
  data[row_i, start_col + 5] <- 0.0001
  row_i <- which(data[, start_col + 7] == "<-15.609")
  data[row_i, start_col + 7] <- -15.609
  data$n <- trimws(data$n)
  data$`k_cleave (1/s)` <- trimws(data$`k_cleave (1/s)`)
  message("rows with blank Kd values:")
  print(length(which(data[, 1] == "")))
  print(length(which(data[, 1] == " ")))
  for (col_i in (start_col - 1):ncol(data)) {
    data[, col_i] <- as.numeric(data[, col_i])
  }
  return(data)
}

MacRaeA1DeltaGCalculation <- function() {
  R <- 1.987e-3 # in kcal K-1 mol-1
  # T <- 310.15 # in K, they used room temperature
  T <- 298.15 # in K, they used room temperature

  Kds <- c(`A`=0.75, `C`=1.8, `G`=1.9, `U`=1.9)
  CI <- c(`A`=0.04, `C`=0.12, `G`=0.09, `U`=0.1)
  RelKds <- Kds/Kds["A"]
  DelDelG <- -R*T*log(RelKds)
  print(DelDelG)
}

GetMpraSites <- function(mirna) {
  all_sites <- grep(sprintf("^%s_", mirna), names(kMPRA), perl=TRUE, value=TRUE)
  sites <- unique(gsub("(.+)_(.+)_(.+)", all_sites, replace="\\2"))
  return(sites)
}


# FILES THAT MAKE THE PROCESSED DATA FILES OF MCGEARY,LIN,ETAL.,2019 ###########
MakeGEOsmallRNASeqSummaryTable <- function() {
  experiment <- "AGO_purity"
  out_matrix <- matrix(0)
  # Gets the expression data within the experiment:
  out <- do.call("cbind", lapply(c("miR-1", "miR-155"), function(mirna) {
    out <- do.call("cbind", lapply(seq(5, 7), function(i_s) {
      out <- do.call("cbind", lapply(c("I", "P"), function(cond) {
        condition <- paste0("S100", i_s, "_", cond)
        counts <- GetMirnaCountData(mirna, condition, unique=FALSE,
                                    no_marker=FALSE, no_adapter=FALSE)
        count_total <- rowSums(counts)
        names(count_total) <- rownames(counts)
        count_total
      }))
      colnames(out) <- paste0("S100", i_s, "_", c("I", "P"))
      out
    }))
    out
  }))
  exclude_rows <- c("Unmapped", "18nt_marker", "30nt_marker", "5p_adapter", "3p_adapter", "dme-miR-14-5p", "xtr-miR-427")
  out <- out[, c(4, 10)]
  colnames(out) <- c("AGO2-miR-1", "AGO2-miR-155")
  write.table(out,
              file="/lab/solexa_bartel/mcgeary/AgoRBNS/GEO_McGearyLin2019/processed_files/smallRNA-seq_counts.txt",
              sep="\t", row.names=TRUE, quote=FALSE, col.names=NA)
  out <- out[!(rownames(out) %in% exclude_rows),]
  out <- t(t(out)/colSums(out))*1e6
  write.table(out,
              file="/lab/solexa_bartel/mcgeary/AgoRBNS/GEO_McGearyLin2019/processed_files/smallRNA-seq_CPM.txt",
              sep="\t", row.names=TRUE, quote=FALSE, col.names=NA)
}

# break
k12merPathDict <- c(`miR-1`="/lab/solexa_bartel/mcgeary/AgoRBNS/miR-1/equilibrium/kds_PAPER/5_12mers_%s-%s_nocombInput_PAPER_logmean.txt",
`miR-124`="/lab/solexa_bartel/mcgeary/AgoRBNS/miR-124/equilibrium/kds_PAPER/5_12mers_%s-%s_PAPER_logmean.txt",
`miR-155`="/lab/solexa_bartel/mcgeary/AgoRBNS/miR-155/equilibrium/kds_PAPER/5_12mers_%s-%s_PAPER_logmean.txt",
`miR-7`="/lab/solexa_bartel/mcgeary/AgoRBNS/miR-7-23nt/equilibrium2_nb/kds_PAPER/5_12mers_%s-%s_global_PAPER_logmean.txt",
`lsy-6`="/lab/solexa_bartel/mcgeary/AgoRBNS/lsy-6/equilibrium/kds_PAPER/5_12mers_%s-%s_PAPER_logmean.txt",
`let-7a`="/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/kds_PAPER/5_12mers_%s-%s_PAPER_logmean.txt")


MakeGEOsmallAgoRBNSProcessedData <- function(experiment="equilibrium", n_constant=5,
                                             sitelist="resubmissionfinal",
                                             multi=TRUE,
                                             buffer=FALSE) {
  sapply(kMirnas, function(mirna) {
    print(mirna)
    # Check for the buffer condition for miR-1 and the `equilibrium2_nb`
    # condition for miR-7-23nt
    if (mirna == "miR-1") buffer <- TRUE
    if (mirna == "miR-7-23nt") experiment <- "equilibrium2_nb"
    # Load the data.
    msXc <- SubfunctionCall(SitesXCounts)
    # Convert the Ts into Us.
    rownames(msXc) <- gsub("T", "U", rownames(msXc))
    # Remove the "-23nt" suffix from miR-7-23nt.
    if (mirna == "miR-7-23nt") mirna <- "miR-7"
    # Assign the output path and write the file.
    geo_path <- file.path("/lab/solexa_bartel/mcgeary/AgoRBNS",
                          "GEO_McGearyLin2019/processed_files/")
    file_path <- sprintf("%sagorbns_%s_sitecounts.txt", geo_path, mirna)
    print(file_path)
    write.table(msXc, file=file_path, sep="\t", row.names=TRUE, quote=FALSE,
                col.names=NA)
    # Second loop where the kd files are moved:
    kd_string <- k12merPathDict[mirna]
    for (p_s in seq(5)) {
      p_e <- p_s + 3
      kd_in_path <- sprintf(kd_string, p_s, p_e)
      print(kd_in_path)
      kds <- read.table(kd_in_path, header=FALSE, sep="\t")
      print(head(kds))
      rownames(kds) <- kds[, 1]
      rownames(kds) <- gsub("_Kd", "", rownames(kds))
      kds <- kds[, 2, drop=FALSE]
      colnames(kds) <- "ln(Kd)"
      print(head(kds))
      kd_out_path <- sprintf("%sagorbns_%s_nt%s-%s_12mer_kds.txt", geo_path,
                             mirna, p_s, p_e)
      write.table(kds, file=kd_out_path, sep="\t", row.names=TRUE, quote=FALSE,
                  col.names=NA)
    }
  })

}

MakeGEOMpraProcessedDataFiles <- function() {
  in_path_1 <- sprintf(
    "ReporterScreen/twist_reporter_assay_v2_rep1_counts.txt", rep
  )
  file1 <- read.table(in_path_1, row.names=1, header=TRUE)[, 1:6]
  print(head(file))
  
  in_path_2 <- sprintf(
    "ReporterScreen/twist_reporter_assay_v2_rep2_counts.txt", rep
  )
  file2 <- read.table(in_path_2, row.names=1, header=TRUE)[, 1:6]
  print(head(file1))
  colnames(file1) <- gsub("\\.", replace="-", colnames(file1), perl=TRUE)
  colnames(file2) <- gsub("\\.", replace="-", colnames(file2), perl=TRUE)
  colnames(file1) <- paste0(colnames(file1), "_rep1")
  print(head(file1))
  colnames(file2) <- paste0(colnames(file2), "_rep2")
  print(head(file2))

  output_file <- cbind(file1, file2)
  print(head(output_file))

  output_file <- output_file[, c(1, 7, 2, 8, 3, 9, 4, 10, 5, 11, 6, 12)]
  print(head(output_file))
  
  # Convert from counts to TPM:
  output_file <- t(t(output_file)/colSums(output_file))*1e6
  print(head(output_file))
  print(colSums(output_file))

  output_path <- "/lab/solexa_bartel/mcgeary/AgoRBNS/GEO_McGearyLin2019/processed_files/mpra_cpm.txt"
  write.table(output_file, file=output_path, sep="\t", row.names=TRUE, 
              quote=FALSE, col.names=NA)
}

# Defines essential lists used throughout scripts:
source("general/Lists.R")
source("general/PlotFunctions.R")
source("general/ModelingFunctions.R")
print("did modeling Functions")