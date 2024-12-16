################################################################################
#GenerateSiteKds_PAPER.py
################################################################################

source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")

# Initial parameters and constants.
args       <- commandArgs(trailingOnly=TRUE)
# args <- c("let-7a", "equilibrium_mmseed_nb", "5", "programmed_collapsed")
mirna      <- args[1]
experiment <- args[2]
n_constant <- args[3]
sitelist   <- args[4]


kExpsThreeP <- c("equil_c_nb", "equil_s_nb", "equil_sc_nb", "equil_c2_nb",
                 "equil_c_alt_nb", "equil_c2_alt_nb", "equil_sc_alt_nb")

if ("-reps" %in% args) {
  reps <- as.integer(args[which(args == "-reps") + 1])
# } else if (sitelist == "programmed_collapsed" || sitelist == "programmed_suppcomp") {
#   reps <- 4
} else {
  reps <- 40
}



if (sitelist == "12mers") {
  mirna.start <- as.integer(args[which(args == "-mir_start") + 1])
  str.mir_start <- sprintf("_%s-%s", mirna.start, mirna.start + 3)
  programmed <- FALSE
  collapsed <- FALSE
  suppcomp <- FALSE
} else if (sitelist == "programmed_suppcomp" & experiment %in% kExpsThreeP) {
  mirna.start <- FALSE
  str.mir_start <- ""
  programmed <- TRUE
  suppcomp <- TRUE
  collapsed <- FALSE
} else if (sitelist == "programmed_collapsed" & experiment %in% kExpsThreeP) {
  mirna.start <- FALSE
  str.mir_start <- ""
  programmed <- TRUE
  suppcomp <- FALSE
  collapsed <- TRUE
} else if (sitelist == "programmed" & experiment %in% kExpsThreeP) {
  mirna.start <- FALSE
  str.mir_start <- ""
  programmed <- TRUE
  suppcomp <- FALSE
  collapsed <- FALSE
} else {
  mirna.start <- FALSE
  str.mir_start <- ""
  programmed <- FALSE
  suppcomp <- FALSE
  collapsed <- FALSE
}
## I DON'T THINK I EVER USE THIS #
if ("-uniq" %in% args) {
  uniq <- TRUE
  str.unique <- "_uniq"
} else {
  uniq <- FALSE
  str.unique <- ""
}
# This is important for the fitting of the miR-1 equilibrium library that has
# an incorrect 3' tag.
if ("-buffer" %in% args) {
  buffer <- TRUE
  str.buffer <- "_buffer3p"
} else {
  buffer <- FALSE
  str.buffer <- ""
}
# This parameter causes the fitting to occur without combining the input
# libraries.
if ("-nocombI" %in% args) {
  combined <- FALSE
  str.combined <- "_nocombInput"
} else {
  combined <- TRUE
  str.combined <- ""
}
# This is used to fix parameters when fitting the random-library equilibrium
# experiments performed by Namita with miR-7 of a range of lengths (22–25 nt).
# The 
if ("-fixed" %in% args & mirna %in% c("miR-7-23nt", "miR-7-24nt", "miR-7-25nt")
    & experiment == "equilibrium2_nb") {
  fixed <- TRUE
  altfixed <- FALSE
  str.fixed <- "_fixed"
} else if ("-altfixed" %in% args & mirna %in% c("miR-7-23nt", "miR-7-24nt",
                                                "miR-7-25nt")
    & experiment == "equilibrium2_nb") {
  fixed <- TRUE
  altfixed <- TRUE
  str.fixed <- "_altfixed"
} else {
  altfixed <- FALSE
  fixed <- FALSE
  str.fixed <- ""
}

# This argument allows for the total library concentration in the modeled
# binding experiment to varied; the default is 100 nM.
if ("-L" %in% args) {
  L <- as.integer(args[which(args == "-L") + 1])
  str.L <- sprintf("_L%s", L)
} else {
  L <- 100
  str.L <- ""
}

# This argument causes the equilibrium binding model to be fit using the site
# count tables that only contain the counts for single instances of each site in
# the site list, with the read counts associated with reads containing multiple
# sites removed.
if ("-single" %in% args) {
  single <- TRUE
  str.single <- "_singleonly"
} else {
  single <- FALSE
  str.single <- ""
}

# This parameter causes the model to be fit using site count tables for which
# reads that contain at least ??? nucleotides of complementarity to the
# competitor oligo were removed, in order to check whether or not their
# inclusion causes a major difference in the fitted Kd values by the binding
# model.
if ("-compcorrect" %in% args) {
  compcorrect <- TRUE
  str.compcorrect <- "_compcorrect"
} else {
  compcorrect <- FALSE
  str.compcorrect <- ""
}

# This parameter causes the model to be fit using site count tables in which
# wobbles are tolerated. In truth I do not remember ever using this parameter
# at the moment have no idea when or why I did this.
if ("-wobble" %in% args) {
  wobble <- TRUE
  str.wobble <- "_wobble"
} else {
  wobble <- FALSE
  str.wobble <- ""
}

# I have no idea what this does at the moment.
if ("-tpomit" %in% args) {
  tpomit <- TRUE
  str.tpomit <- "_tpomit"
} else {
  tpomit <- FALSE
  str.tpomit <- ""
}

# I have no idea what this does at the moment.
if ("-AGOfixed_bypass" %in% args) {
  AGOfixed_bypass <- TRUE
  str.AGOfixedbypass <- "_AGOfixedbypass"
} else {
  AGOfixed_bypass <- FALSE
  str.AGOfixedbypass <- ""
}

# The flag `tpomit2` indicates whether or not the lowest concentration in Thy
# Pham's miR-1 data should be omitted when fittig the Kds, in order to check
# how the Kds look incomparison to the original miR-1 data. The reason for doing
# this is to check whether or not it is a good idea to use her data, with the
# benfit being that the library used is consistent with all the other non-miR-1
# data sets in the paper other than miR-7 (which requires a different library
# due to the presence of a miR-7 6mer-m8 site in the constant sequence of the
# library).
if ("-tpomit2" %in% args) {
  tpomit2 <- TRUE
  str.tpomit2 <- "_tpomit2"
} else {
  tpomit2 <- FALSE
  str.tpomit2 <- ""
}


if (mirna == "miR-7-24nt" && experiment == "equilibrium_tp") {
  fixed8mer <- TRUE
} else {
  fixed8mer <- FALSE
}

# Flag argument that allows for changing the minimum Kd value during the
# optimization routine. This was implemented for looking at Thy's miR-7-24nt
# data, for which there is not saturation. By changing the minimum Kd value to
# 1e-3 (which is close to the Kd value obtained from the combined optimization
# over Namita's data sets), the Kd values agree readily with Namita's data
# sets. Allowing the Kd value to be the initially decided (by me) value of
# 1e-6 causes Thy's miR-7-24nt data to take on the minimum value of 1e-6.
if ("-minkd" %in% args) {
  minkd <- as.double(args[which(args == "-minkd") + 1])
  str.minkd <- sprintf("_minkd%s", args[which(args == "-minkd") + 1])
} else {
  minkd <- 1e-6
  str.minkd <- ""
}

if ("-nbomitc" %in% args) {
  nbomitc <- TRUE
  str.nbomitc <- "_nbomitc"
} else {
  nbomitc <- FALSE
  str.nbomitc <- ""
}
# This is for the programmed libraries, and allows for the splitting of greater
# than 9-nt sites into all of its consituent 9-nt sites.
# if ("-split9mers" %in% args) {
#   split9mers <- TRUE
#   str.split9mers <- "_split9mers"
# } else {
#   split9mers <- FALSE
#   str.split9mers <- ""
# }
# This is for the programmed libraries, and allows for the extension out to
# 11-nt sites.
# if ("-outto11mers" %in% args) {
#   outto11mers <- TRUE
#   str.outto11mers <- "_outto11mers"
# } else {
#   outto11mers <- FALSE
#   str.outto11mers <- ""
# }

# if ("-downto2mers" %in% args) {
#   downto2mers <- TRUE
#   str.downto2mers <- "_downto2mers"
# } else {
#   downto2mers <- FALSE
#   str.downto2mers <- ""
# }

# if ("-downto3mers" %in% args) {
#   downto3mers <- TRUE
#   str.downto3mers <- "_downto3mers"
# } else {
#   downto3mers <- FALSE
#   str.downto3mers <- ""
# }

# if ("-downto4mers" %in% args) {
#   downto4mers <- TRUE
#   str.downto4mers <- "_downto4mers"
# } else {
#   downto4mers <- FALSE
#   str.downto4mers <- ""
# }

# if ("-end3prand" %in% args) {
#   end3prand <- TRUE
#   str.end3prand <- "_end3prand"
# } else {
#   end3prand <- FALSE
#   str.end3prand <- ""
# }

# 
if ("-collapsemm" %in% args && sitelist == "programmed_suppcomp") {
  collapsemm <- TRUE
  str.collapsemm <- "_collapsemm"
} else {
  collapsemm <- FALSE
  str.collapsemm <- ""
}

# if ("-new" %in% args) {
#   new <- TRUE
#   str.new <- "_new"
# } else {
#   new <- FALSE
#   str.new <- ""
# }

# MAIN #########################################################################
# 1. Get data table:
mir_list <- unlist(strsplit(mirna, ","))
exp_list  <- unlist(strsplit(experiment, ","))
mir_exp_list <- data.frame(mirna=rep(mir_list, length(exp_list)),
                           experiment=rep(exp_list, each=length(mir_list)))

sXc <- apply(mir_exp_list, 1, function(row) {
  mirna <- row[1]
  sXc_i <- SitesXCounts(row[1], experiment=row[2], n_constant=n_constant,
                        sitelist=sitelist, mirna.start=mirna.start,
                        uniq=uniq, buffer=buffer)
  print(head(sXc_i))
  if (sitelist == "programmed") {
    sXc_i_collapsed <- SitesXCounts(
      row[1], experiment=row[2], n_constant=n_constant, sitelist="programmed_suppcomp",
      mirna.start=mirna.start, uniq=uniq, buffer=buffer)
    sXc_i <- rbind(sXc_i[grep("&", rownames(sXc_i), invert=TRUE), ],
                   sXc_i_collapsed[grep("&", rownames(sXc_i_collapsed)), ])
  }
  print(sXc_i[grep("^5mer-m9\\.13\\|", rownames(sXc_i), perl=TRUE), ])
  if (single) {
    msXc <- SitesXCounts(row[1], experiment=row[2], n_constant=n_constant,
                         sitelist=sitelist, mirna.start=mirna.start,
                         uniq=uniq, buffer=buffer, multisite=TRUE)
    regex_split         <- ",\\(\\d*\\),"
    regex_split_capture <- ",\\((\\d*)\\),"
    regex_site <- "[^,]*"
    regex_site_capture <- "([^,]*)"

    single_sites = grep(regex_split, rownames(msXc), invert=TRUE)
    double_sites = grep(sprintf("^%s%s%s$", regex_site, regex_split,
                                regex_site), rownames(msXc), perl=TRUE)
    nonoverlapping_sites = grep("\\|", rownames(msXc), invert=TRUE, perl=TRUE)
    single_nonoverlapping_sites = intersect(single_sites, nonoverlapping_sites)
    sXc_out <- msXc[single_nonoverlapping_sites, ][rownames(sXc_i), ]

  } else {
    sXc_out <- sXc_i
  }
  # Portion of code where miR-124 sites are corrected. #########################
  if (compcorrect) {
    sXf <- t(t(sXc_out)/colSums(sXc_out))
    sXr <- sXf/sXf[, 1]
    conditions <- c("I", "0.4", "1.26", "4", "12.6", "40")
    vec_cond_rep <- c("0.4", "1.265", "4", "12.65", "40")
    sXr_new <- sXr
    for (row_i in grep("AA-", rownames(sXc_out))) {
      str_site_seq <- GetSiteSeq(mirna, rownames(sXc_out)[row_i])
      kmer_len <- GetMaxCompetitorOligoEnrichment(mirna, str_site_seq,
                                                  wobble=wobble)
      kXc <- sapply(conditions, function(condition) {
        # Get the total number of kmers.
        counts <- CompetitorKmers(mirna, condition, kmer_len, off=0,
                                  experiment=experiment, wobble=wobble)
        total_counts <- counts["Total", ]
        kmers <- unlist(counts)
        fracs <- kmers[-length(kmers)]/total_counts
        names(fracs) <- rownames(counts)[-length(kmers)]
        fracs
      })
      kXr <- kXc/kXc[, 1]
      vec_r <- colMeans(kXr)
      vec_r_rep <- vec_r[2:length(vec_r)] - sXr["None", vec_cond_rep]
      sXr_new[row_i, vec_cond_rep] <- sXr_new[row_i, vec_cond_rep] - vec_r_rep
    }
    sXr <<- sXr
    sXr_new <<- sXr_new

    sXc_out_new <- sXc_out
    sXc_out_new[, vec_cond_rep] <- (
      sXc_out[, vec_cond_rep]/sXr[, vec_cond_rep]*sXr_new[, vec_cond_rep]
    )
    sXc_out <- sXc_out_new
  }
  ######## END COMP CORRECT PART OF PARAMETER PARSING ##########################
  # PART OF SCRIPT DEALING WITH THE TP EXPERIMENTS (HOW TO DEAL WITH REPLICATES)
  if (mirna == "miR-124" & experiment == "equilibrium_tp" & tpomit) {
    sXc_out <- sXc_out[, c(-5, -6, -7, -8)]
  }
  # Removes the 0.1265 concentration for Thy Pham's miR-1 data, to make it
  # consistent with the rest of the original experiments.
  if (mirna == "miR-1" & experiment == "equilibrium_tp" & tpomit2) {
    sXc_out <- sXc_out[, c(-8)]
  }

  if (mirna == "miR-124" & experiment == "equilibrium_2_tp") {
  	if ("-rep1" %in% args) {
  		sXc_out <- sXc_out[, c(1, 3, 4, 6, 8, 10, 12, 14, 16)]
  		str.tp2rep <<- "_rep1"
  	} else if ("-rep2" %in% args) {
  		sXc_out <- sXc_out[, c(2, 3, 5, 7, 9, 11, 13, 15, 17)]
		str.tp2rep <<- "_rep2"
  	} else {
  		str.tp2rep <<- ""
  	}
    colnames(sXc_out) <- gsub("\\.([12])$", colnames(sXc_out), replace=",\\1")
  } else if (mirna == "miR-7-24nt" & experiment == "equilibrium_tp") {
    colnames(sXc_out)[2] <- "I_combined"
    str.tp2rep <<- ""
  } else {
  	str.tp2rep <<- ""
  }
  # RELAVENT TO THE THREEP LIBRARIES ###########################################
  if (mirna == "let-7a-21nt") {
    if (experiment == "equil_c2_nb" || experiment == "equil_c2_alt_nb") {
      sXc_out[, "12.65"] <- sXc_out[, "12.65"] + sXc_out[, "12.65_2"]
      sXc_out <- sXc_out[, -which(colnames(sXc_out) == "12.65_2")]
    } else if (experiment == "equil_s_nb") {
      sXc_out[, "12.65"] <- sXc_out[, "12.65"] + sXc_out[, "12.65_2"]
      sXc_out[,     "4"] <- sXc_out[,     "4"] + sXc_out[,     "4_2"]
      sXc_out <- sXc_out[, -which(colnames(sXc_out) %in% c("4_2", "12.65_2"))]      
    }
    if (experiment == "equil_c_nb" & nbomitc) {
      sXc_out <- sXc_out[, -7]
    }
  }
  # This conditional deals with potentially having different amounts of
  # compensatory and supplementary sites in the programmed region.
  if (collapsemm) {
    supp_inds <- intersect(c("8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-m8", "6mer-A1"), rownames(sXc_out))
    comp_inds <- grep("^8mer-mm[ACTG][2-7]$", rownames(sXc_out), perl=TRUE)
    # If there are any supplemntary sites, then assign the first row to 1,
    # and if not, assign the first row to 0 (so that it isn't used).
    if (length(supp_inds) != 0) {
      supp_final_ind <- 1
    } else {
      supp_final_ind <- 0
    }
    sXc_out[supp_final_ind, ] <- colSums(sXc_out[supp_inds, ])
    sXc_out[supp_final_ind + 1, ] <- colSums(sXc_out[comp_inds, ])
    sXc_out <- sXc_out[-1*((supp_final_ind + 2):max(comp_inds)), ]
    rownames(sXc_out)[supp_final_ind] <- "Supp"
    rownames(sXc_out)[supp_final_ind + 1] <- "Comp"
  }
  sXc_out
})

# Assign the strings for the output files.
if (length(sXc) > 1) {
  str.global <- "_global"
} else {
  str.global <- ""
}
str_condition <- sprintf(
  "%s_%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s_parallel_PAPER", n_constant,
  sitelist, str.unique, str.buffer, str.mir_start, str.combined, str.global,
  str.fixed, str.L, str.single, str.compcorrect, str.wobble, str.tpomit,
  str.tpomit2, str.tp2rep, str.minkd, str.AGOfixedbypass, str.nbomitc,
  str.collapsemm
)

kOutputFileMean <- GetAnalysisPath(
  mir_list[[1]], experiment, condition=sprintf("%s_logmean", str_condition),
  analysis_type="kds_PAPER"
)

kOutputFileProg <- GetAnalysisPath(
  mir_list[[1]], experiment, condition=sprintf("%s_sorted", str_condition),
  analysis_type="kds_PAPER"
)

kOutputFileProg_temp <- GetAnalysisPath(
  mir_list[[1]], experiment, condition=str_condition, analysis_type="kds_PAPER"
)

# Assign the standard deviation, and the plot names for the optimization
# routine.
if (sitelist == "12mers") {
  pars_sd <- 0.01
  plot_ <- TRUE
  plotname <- paste0("12merfits/", paste(mirna, experiment, n_constant,
                                         sitelist, mirna.start,
                     sep="_"), ".pdf")
} else if (sitelist == "programmed" || sitelist == "programmed_collapsed") {
  pars_sd <- 0.01
  plot_ <- FALSE
  plotname <- NULL
  factr <- 1e7
} else {
  pars_sd <- 0.1
  plot_ <- FALSE
  plotname <- NULL
  factr <- 1e7
}



InitializeEquilSitePars <- function(sXc, combined=TRUE, fixed=FALSE,
                                    fixed8mer=FALSE, programmed=FALSE) {
  message("In the initialization function.")
  n_mir = length(sXc)
  if (programmed) ps <- 1
  else            ps <- 0
  # ps <- 0.001
  # ps <- 0
  l <- SubfunctionCall(GetInputEquil, sXc=sXc[[1]] + ps)
  data <- SubfunctionCall(GetDataEquil, sXc=sXc[[1]] + ps)
  if (programmed) {
    kds <- log(Norm(l)/Norm(rowSums(data)))
  } else {
    kds <- log(Norm(l)/Norm(rowSums(data)))
  }
  if (programmed) {
    if (collapsemm) {
      inds_mm <- grep("^Comp$", rownames(data), perl=TRUE)
    } else if (collapsed || suppcomp) {
      inds_mm <- grep("^8mer-mm[ACGT][2-7]$", names(kds), perl=TRUE)    
    } else {
      inds_mm <- grep("^NA\\|NA\\|8mer-mm[ACGT][2-7]$", names(kds), perl=TRUE)
    }
    max_mm <- inds_mm[which.max(kds[inds_mm])]
    kds <- kds - kds[max_mm]
  } else {
    kds <- kds - kds[length(kds)]    
  }
  names(kds) <- paste0(rownames(sXc[[1]]), "_Kd")
  if (fixed) {
    if (altfixed) {
      combined_fixed <- TRUE
    } else {
      combined_fixed <- FALSE
    }
    pars <- EquilPars(mirna, experiment=experiment, n_constant=n_constant,
                      sitelist=sitelist, uniq=uniq, combined=combined_fixed,
                      global=TRUE)
    bgs <- log(pars[sprintf("bg_%s", mirna),]$Mean)
    As <- log(pars[sprintf("AGO_%s", mirna),]$Mean)
  } else if (fixed8mer) {
    message("About to print the parameters:")
    pars <- EquilPars("miR-7-23nt", experiment="equilibrium2_nb",
                      n_constant=n_constant, sitelist=sitelist, uniq=uniq,
                      combined=combined)
    message("Namita parameters:")
    kds["8mer_Kd"] <- log(pars["8mer_Kd", "Mean"])
    bgs <- rep(log(0.1), n_mir)
    As <- rep(log(20), n_mir)    
  } else {
    bgs <- rep(log(0.1), n_mir)
    As <- rep(log(20), n_mir)    
  }
  names(bgs) <- sprintf("bg_%s", mir_list)
  names(As) <- sprintf("AGO_%s", mir_list)
  pars <- c(kds, bgs, As)
  random_pars <- rnorm(length(pars), 0, pars_sd)
  if (fixed) {
    random_pars[length(pars) - 1] <- 0
    random_pars[length(pars)    ] <- 0
  } else if (fixed8mer) {
    random_pars[1] <- 0
  } else if (programmed) {
    random_pars[max_mm] <- 0
  }
  names(random_pars) <- names(pars)
  pars <- pars + random_pars
  if (!programmed) {
    pars["None_Kd"] <- 0  
  }
  return(pars)
}

tick <- 0
OptimizeEquilSitePars <- function(sXc, pars=NULL, fixed=fixed, AGOfixed=FALSE,
                                  fixed8mer=FALSE, programmed=FALSE,
                                  plotname_=plotname) {
  message("In the optimization function.")
  time_start <- proc.time()[3]
  tick <- 0
  if (is.null(pars) == TRUE) {
    message("About to initialize parameters:")
    initial.pars <- InitializeEquilSitePars(sXc, combined=combined, fixed=fixed,
                                            fixed8mer=fixed8mer,
                                            programmed=programmed)
  } else {
    initial.pars <- pars
  }
  message("Parameters are initialized.")
  L_ <- L
	start_i <- 3
	end_i <- ncol(sXc[[1]]) - 1
  if (length(sXc) > 1) {
    data <- do.call(cbind, sapply(sXc, function(sXc_i) {
      as.matrix(sXc_i[, start_i:end_i])
    }))
  } else {
    data <- sXc[[1]][, start_i:end_i]
  }
  data <- data + 1
	if ((mirna == "miR-124" & experiment == "equilibrium_2_tp") &
	    (!("-rep1" %in% args) & !("-rep2" %in% args))) {
  	n_j <- sapply(sXc, function(sXc_i) {
    	ncol(sXc_i) - start_i - 1
  	})
  } else if ((mirna == "miR-7-24nt" & experiment == "equilibrium_tp")) {
    n_j <- sapply(sXc, function(sXc_i) {
      ncol(sXc_i) - start_i - 1
    })
  } else {
  	n_j <- sapply(sXc, function(sXc_i) {
    	ncol(sXc_i) - start_i
  	})
  }
  ind_l <- which(colnames(sXc[[1]]) == "I_combined")
  l <- Norm(sXc[[1]][, ind_l] + 1)*L_
  Y <- colSums(data)
  dil <- as.numeric(gsub(",[1234]$", colnames(data), replace=""))/100.0
  n_i <- nrow(data)
  n_mir <- length(sXc)
  data_vec <- as.numeric(as.matrix(data))
  lower <- log(rep(c(minkd, 0.01, 0.1), c(n_i, n_mir, n_mir)))
  upper <- log(rep(c(1e1, 10, 10), c(n_i, n_mir, n_mir)))

  # A normalizing constant to divide the cost and gradients by, such that it is
  # a per-site-and-condition negative-log-likelihood, rather than an overall
  # negative log likelihood (which will continue to grow with more sites and
  # conditions, because with each new site and new condition the likelihood of
  # the data must necessarily decrease).
  norm_constant <- n_i*sum(n_j)
  if (fixed) {
    zero_grad <- n_i + c(0, 1, 2)
  } else if (AGOfixed) {
    zero_grad <- n_i + c(0, 1, 2)
  } else if (programmed) {
    # Find the base kd indeces.
    if (collapsemm) {
      zero_grad <- c(which(names(initial.pars) == "Comp_Kd"))
    } else {
      if (suppcomp) {
        base_inds <- grep("^8mer-mm[ACGT][2-7]_Kd$", names(initial.pars), perl=TRUE)
      } else {
        base_inds <- grep("^NA\\|NA\\|8mer-mm[ACGT][2-7]_Kd$", names(initial.pars), perl=TRUE)
      }
      # Find which of these indeces corresponds to the maximum compensatory base site.
      max_base_ind <- which(initial.pars[base_inds] == max(initial.pars[base_inds]))
      # Assign that index as the zero_grad_kd_ind.
      zero_grad_kd_ind <- base_inds[max_base_ind]
      initial.pars[zero_grad_kd_ind] <- 0
      zero_grad <- c(zero_grad_kd_ind)
    }
  } else if (fixed8mer) {
    ind_8mer <- which(names(initial.pars) == "8mer_Kd")
    zero_grad <- c(ind_8mer, n_i)
  } else {
    zero_grad <- c(n_i)
  }
  n_z <- length(zero_grad)
  # Initialize global constants that can be tracked within the `CostC` function,
  # that 
  t_global <<- 0
  time_interval <<- proc.time()[3]
  message("About to enter the optimization routine.")
  print(dim(data))
  print(length(data_vec))
  print(length(initial.pars))
  print(n_z)
  print(initial.pars[zero_grad])
  solution <- optim(initial.pars,
                    CostC,
                    gr = GradC,
                    data = data_vec,
                    dil  = dil,
                    l = l,
                    L = L_,
                    Y = Y,
                    n_i = n_i,
                    n_j = n_j,
                    n_mir = n_mir,
                    zero_grad = zero_grad,
                    n_z = n_z,
                    fixed = fixed,
                    norm_constant=norm_constant,
                    upper_ = upper,
                    lower_ = lower,
                    plot_ = plot_,
                    plotname = plotname_,
                    method = "L-BFGS-B",
                    lower=lower,
                    upper=upper,
                    # control = list(maxit=1000000, fnscale=1),
                    control = list(maxit=1000000, factr=factr, fnscale=1))

  # Print the total amount of time the optimization took.
  time_elapsed <- proc.time()[3] - time_start
  print(time_elapsed)
  # Return the optimized parameters.
  pars <- solution$par
  # if (programmed) {
  #   max_par <- max(pars[7:24])
  #   pars[1:(length(pars) - 2)] <- pars[1:(length(pars) - 2)] - max_par
  # }
  pars/log(10)
}

message("About to do the optimization:")
pars.MLE <- OptimizeEquilSitePars(sXc, fixed=fixed, fixed8mer=fixed8mer,
                                  programmed=programmed)

print(head(pars.MLE))

# write.table(file=sprintf("general/temp_%s_kds.txt", mirna),
#             10^pars.MLE, sep="\t", quote=FALSE)


ncol_span = ncol(sXc[[1]]) - 3

resampled.pars <- matrix(NaN, nrow=length(pars.MLE), ncol=ncol_span*reps,
                         dimnames=list(names(pars.MLE), 1:(ncol_span*reps)))
print("RESAMPLING")

if ((mirna == "miR-124" & experiment == "equilibrium_2_tp") &
		(!("-rep1" %in% args) & !("-rep2" %in% args))) {
	max_j <- ncol(sXc[[1]]) - 6
} else if (mirna == "miR-7-24nt" & experiment == "equilibrium_tp") {
  max_j <- ncol(sXc[[1]]) - 6
} else {
	max_j <- ncol(sXc[[1]]) - 4
}

full_number_of_reps <- 0
for (j in 0:max_j) {
  for (i in 1:reps) {
    i_full <- j*reps + i
    tick <- 0
    # Resample the data:
    sXc.resample <- lapply(sXc, function(sXc_i) {
      sXc_new <- apply(sXc_i, 2, function(col) {rmultinom(1, sum(col), col)})
      rownames(sXc_new) <-rownames(sXc_i)
      return(sXc_new)
    })

    # If not 12mers, remove a column:
    if (sitelist != "12mers") {
      # if (mirna == "miR-7-23nt" & experiment == "equilibrium2_nb") {
      #   min_col <- 4
      # } else {
      #   min_col <- 3
      # }
      sXc.resample.withhold <- lapply(sXc.resample, function(sXc_i) {
    	  if ((mirna == "miR-124" & experiment == "equilibrium_2_tp") &
			  		(!("-rep1" %in% args) & !("-rep2" %in% args))) {
    	  	 leave_out <- j + 4
  	  	 } else {
    	  	 	leave_out <- j + 3
  	  	 }
        sXc_i[, -leave_out]
      })
    } else {
      sXc.resample.withhold <- sXc.resample
      plotname <- paste0("12merfits/", paste(mirna, experiment, n_constant, sitelist, mir_start,
                     sep="_"), "_rep", j*reps + i, ".pdf")
    }
    if (length(sXc.resample.withhold) != 1) {
      sXc.resample.withhold <- sXc.resample
    }
    if (mirna == "miR-7-23nt" & experiment == "equilibrium2_nb") {
      AGOfixed <- TRUE
    } else {
      AGOfixed <- FALSE
    }
    if (AGOfixed_bypass) {
      AGOfixed <- FALSE    
    }
    resampled.pars[, i_full] <- OptimizeEquilSitePars(sXc.resample.withhold,
                                                      pars = log(10)*pars.MLE,
                                                      fixed=fixed,
                                                      AGOfixed=AGOfixed,
                                                      fixed8mer=fixed8mer,
                                                      plotname_=plotname,
                                                      programmed=programmed)
    resampled.pars.sort <- t(apply(resampled.pars[, 1:i_full, drop=FALSE], 1,
                                   sort))
    full_number_of_reps <- full_number_of_reps + 1
    write.table(file=kOutputFileProg, resampled.pars.sort, sep="\t",
                quote=FALSE)
    output <- 10^cbind(pars.MLE,
                       rowMeans(resampled.pars, na.rm=TRUE),
                       resampled.pars.sort[, ceiling(0.025*i_full)],
                       resampled.pars.sort[, ceiling(0.5*i_full)],
                       resampled.pars.sort[, ceiling(0.975*i_full)])
    colnames(output) <- c("Full","Mean", "Lower_CI", "Median", "Upper_CI")
    sapply(mir_list, function(mirna) {
      kOutputFile <- GetAnalysisPath(mirna, experiment, condition=str_condition,
                               analysis_type="kds_PAPER")
      print(kOutputFile)
      write.table(file=kOutputFile, output, sep="\t", quote = FALSE)
    })
    if (sitelist == "12mers") {
      kds.mean <- log(10^rowMeans(resampled.pars, na.rm=TRUE))
      # print(kds.mean[1:5])
      write.table(file=kOutputFileMean, kds.mean, sep="\t", quote=FALSE,
                  row.names=FALSE)
    }
  write.table(file=kOutputFileProg, resampled.pars.sort, sep="\t", quote = FALSE)
  }
}

write.table(file=sprintf("general/temp_%s_final_new_kds.txt", mirna),
            output, sep="\t", quote=FALSE)

print("done")
warnings()