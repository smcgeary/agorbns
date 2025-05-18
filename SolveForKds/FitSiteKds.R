################################################################################
#FitSiteKds.R
################################################################################

source("general/general.R")
source("general/ModelingFunctions.R")

# Initial parameters and constants.
args       <- commandArgs(trailingOnly=TRUE)
mirna      <- args[1]
experiment <- args[2]
n_constant <- args[3]
sitelist   <- args[4]


kExpsThreeP <- c("equil_c_nb", "equil_s_nb", "equil_sc_nb", "equil_c2_nb",
                 "equil_c_alt_nb", "equil_c2_alt_nb", "equil_sc_alt_nb")

if ("-reps" %in% args) {
  reps <- as.integer(args[which(args == "-reps") + 1])
} else {
  reps <- 40
}



if (sitelist == "12mers") {
  mirna.start <- as.integer(args[which(args == "-mir_start") + 1])
  str.mir_start <- sprintf("_%s-%s", mirna.start, mirna.start + 3)
  # programmed <- FALSE
  # collapsed <- FALSE
  # suppcomp <- FALSE
} else if (sitelist == "programmed_suppcomp" & experiment %in% kExpsThreeP) {
  mirna.start <- FALSE
  str.mir_start <- ""
  # programmed <- TRUE
  # suppcomp <- TRUE
  # collapsed <- FALSE
} else if (sitelist == "programmed_collapsed" & experiment %in% kExpsThreeP) {
  mirna.start <- FALSE
  str.mir_start <- ""
  # programmed <- TRUE
  # suppcomp <- FALSE
  # collapsed <- TRUE
} else if (sitelist == "programmed" & experiment %in% kExpsThreeP) {
  mirna.start <- FALSE
  str.mir_start <- ""
  # programmed <- TRUE
  # suppcomp <- FALSE
  # collapsed <- FALSE
} else {
  mirna.start <- FALSE
  str.mir_start <- ""
  # programmed <- FALSE
  # suppcomp <- FALSE
  # collapsed <- FALSE
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
if ("-buffer3p" %in% args) {
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

if ("-sumseed" %in% args && (sitelist == "progthrp_suppcomp" |
                             sitelist == "randthrp_suppcomp" |
                             sitelist == "randthrp_comp")) {
  sumseed <- TRUE
  str.sumseed <- "_sumseed"
} else {
  sumseed <- FALSE
  str.sumseed <- ""
}

if ("-start_mm" %in% args & "-stop_mm" %in% args) {
  start_mm <- as.numeric(args[which(args == "-start_mm") + 1])
  stop_mm <- as.numeric(args[which(args == "-stop_mm") + 1])
  str.full_mm <- sprintf("m%s.%smmsd_", start_mm, stop_mm)
} else {
  start_mm <- FALSE
  stop_mm <- FALSE
  str.full_mm <- ""
}

if ("-new" %in% args) {
  new <- TRUE
  str.new <- "_new"
} else {
  new <- FALSE
  str.new <- ""
}

if ("-new2" %in% args) {
  new2 <- TRUE
  str.new2 <- "_new2"
} else {
  new2 <- FALSE
  str.new2 <- ""
}

if ("-original" %in% args) {
  original <- TRUE
  str.original <- "_original"
} else {
  original <- FALSE
  str.original <- ""
}

if ("-filt_by_original" %in% args) {
  filt_by_original <- TRUE
  str.filt_by_original <- "_filtbyoriginal"
} else {
  filt_by_original <- FALSE
  str.filt_by_original <- ""
}

# MAIN #########################################################################
# 1. Get data table:
mir_list <- unlist(strsplit(mirna, ","))
exp_list  <- unlist(strsplit(experiment, ","))
mir_exp_list <- data.frame(mirna=rep(mir_list, length(exp_list)),
                           experiment=rep(exp_list, each=length(mir_list)))

sXc <- apply(mir_exp_list, 1, function(row) {
  mirna <- row[1]
  sXc_i <- SitesXCounts(
    row[1], experiment=row[2], n_constant=n_constant, sitelist=sitelist,
    mirna.start=mirna.start, uniq=uniq, buffer=buffer, start_mm=start_mm,
    stop_mm=stop_mm, new=new, new2=new2, original=original,
    filt_by_original=filt_by_original
  )
  # This swaps the "8mer&4mer-m9.12|NA|8mer" sites in the progthrp counts for
  # the "8mer&4mer-m9.12|NA|Supp" site types, to cut down on the number of sites
  # that need to be fit by the optimization routine.
  if (sitelist == "progthrp") {
    # Get the suppcomp version of the programmed library
    sXc_i_collapsed <- SitesXCounts(
      row[1], experiment=row[2], n_constant=n_constant,
      sitelist="progthrp_suppcomp", mirna.start=mirna.start, uniq=uniq,
      buffer=buffer, start_mm=start_mm, stop_mm=stop_mm, new=new, new2=new2
    )
    # Remove all of the rows that contain an "&" character from the original
    # site count table, and replace with those of the suppcomp version of the
    # table.
    sXc_i <- rbind(sXc_i[grep("&", rownames(sXc_i), invert=TRUE), ],
                   sXc_i_collapsed[grep("&", rownames(sXc_i_collapsed)), ])
  }
  if (single) {
    msXc <- SitesXCounts(
      row[1], experiment=row[2], n_constant=n_constant, sitelist=sitelist,
      mirna.start=mirna.start, uniq=uniq, buffer=buffer, multisite=TRUE,
      original=original, filt_by_original=filt_by_original
    )
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
  # RELEVANT TO THE THREEP LIBRARIES ###########################################
  if (mirna == "let-7a-21nt") {
    if (experiment == "equil_c2_nb" || experiment == "equil_c2_alt_nb") {
      sXc_out[, "12.65"] <- sXc_out[, "12.65"] + sXc_out[, "12.65_2"]
      sXc_out <- sXc_out[, -which(colnames(sXc_out) == "12.65_2")]
    } else if (experiment == "equil_s_nb") {
      sXc_out[, "12.65"] <- sXc_out[, "12.65"] + sXc_out[, "12.65_2"]
      sXc_out[,     "4"] <- sXc_out[,     "4"] + sXc_out[,     "4_2"]
      sXc_out <- sXc_out[, -which(colnames(sXc_out) %in% c("4_2", "12.65_2"))]      
    }
    # if (experiment == "equil_c_nb" & nbomitc) {
    #   sXc_out <- sXc_out[, -7]
    # }
  }
  # This conditional deals with potentially having different amounts of
  # compensatory and supplementary sites in the programmed region.
  if (sumseed) {
    # Get the indeces for each of the summations.
    if (grepl("suppcomp", sitelist)) {
      supp_inds <- grep("^(8mer|7mer-m8|7mer-A1|6mer)$", rownames(sXc_out),
                        perl=TRUE)

      offset_inds <- grep("^(6mer-m8|6mer-A1)$", rownames(sXc_out), perl=TRUE)
    } else {
      supp_inds <- c()
      offset_inds <- c()
    }
    comp_inds <- grep("^8mer-mm[ACTG][2-7]$", rownames(sXc_out), perl=TRUE)
    # Make the new top of the dataframe.
    # Name the rows.
    if (grepl("suppcomp", sitelist)) {
      sXc_base <- rbind(colSums(sXc_out[comp_inds, ]), 
                        colSums(sXc_out[offset_inds, ]),
                        colSums(sXc_out[supp_inds, ]))
      rownames(sXc_base) <- c("Comp", "Offset", "Supp")
    } else {
      sXc_base <- t(as.matrix(colSums(sXc_out[comp_inds, , drop=FALSE])))
      print(sXc_base)
      rownames(sXc_base) <- c("Comp")
    }
    print(sXc_base)
    # Concatenate the new base with the original dataframe minus the separated
    # seed-site counts.
    sXc_out <- rbind(sXc_base, sXc_out[-c(supp_inds, offset_inds, comp_inds), ])
  }
  cols_data <- setdiff(colnames(sXc_out), c("I", "I_combined", "0"))
  inds_drop <- which(rowSums(sXc_out[, cols_data]) == 0)
  # sXc_out <- sXc_out[inds_drop]
  print(head(sXc_out))
  return(sXc_out)
})

# Assign the strings for the output files.
if (length(sXc) > 1) {
  str.global <- "_global"
} else {
  str.global <- ""
}
str_condition <- sprintf(
  "%s_%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s_PAPER", n_constant, str.full_mm,
  sitelist, str.unique, str.buffer, str.mir_start, str.combined, str.global,
  str.fixed, str.L, str.single, str.compcorrect, str.wobble, str.tpomit,
  str.tpomit2, str.tp2rep, str.minkd, str.AGOfixedbypass, str.nbomitc,
  str.sumseed, str.new, str.new2
)


kOutputFileMean <- GetAnalysisPath(
  mir_list[[1]], experiment, condition=sprintf(
    "%s_logmean%s%s", str_condition, str.original, str.filt_by_original
  ),
  analysis_type="kds_PAPER"
)
kOutputFileProg <- GetAnalysisPath(
  mir_list[[1]], experiment, condition=sprintf(
    "%s_sorted%s%s", str_condition, str.original, str.filt_by_original
  ),
  analysis_type="kds_PAPER"
)

kOutputFileProg_temp <- GetAnalysisPath(
  mir_list[[1]], experiment, condition=sprintf(
    "%s%s%s", str_condition, str.original, str.filt_by_original
  ), analysis_type="kds_PAPER"
)

print(kOutputFileMean)
print(kOutputFileProg)
print(kOutputFileProg_temp)



# Assign the standard deviation, and the plot names for the optimization
# routine.
if (sitelist == "12mers") {
  pars_sd <- 0.01
  plot_ <- TRUE
  plotname <- paste0("12merfits/", paste(mirna, experiment, n_constant,
                                         sitelist, mirna.start,
                     sep="_"), ".pdf")
  factr <- 1e7
} else if (sitelist == "progthrp") {
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
                                    fixed8mer=FALSE) {
  message("In the initialization function.")
  n_mir = length(sXc)
  # TO DO ##### FIGURE OUT THIS PSEUDOCOUNT FOR PROGRAMMED LIBRARIES
  if (sitelist %in% c("progthrp", "progthrp_suppcomp", "randthrp", "randthrp_suppcomp", "randthrp_comp")) ps <- 1
  else            ps <- 0
  # ps <- 0
  print("476")
  print(combined)
  l <- SubfunctionCall(GetInputEquil, sXc=sXc[[1]] + ps)
  print("478")
  data <- SubfunctionCall(GetDataEquil, sXc=sXc[[1]] + ps)
  print("478")
  # Initialize the Kd values.
  kds <- log(Norm(l)/Norm(rowSums(data)))
  # Assign the reference Kd. If there is no "None" site, need to use either the
  # worst mismatchsite, or the site "Comp" that refers to the summed counts of
  # all the mismatched compensatory base seed sites. The index is explicitly put
  # here because it is used later for the randomization and the index for which
  # the derivative should be zero.
  if (!("None" %in% rownames(data))) {
    if (sumseed) {
      ind_ref <<- grep("^Comp$", rownames(data), perl=TRUE)
    } else {
      inds_mm <- grep("^8mer-mm[ACGT][2-7]$", names(kds), perl=TRUE)
      ind_ref <<- inds_mm[which(kds[inds_mm] == max(kds[inds_mm]))]
    }
  } else {
    ind_ref <<- grep("^None$", rownames(data), perl=TRUE)
  } 
  kds <- kds - kds[ind_ref]
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
  names(random_pars) <- names(pars)
  if (fixed) {
    random_pars[length(pars) - 1] <- 0
    random_pars[length(pars)    ] <- 0
  } else if (fixed8mer) {
    random_pars[1] <- 0
  }
  random_pars[ind_ref] <- 0
  pars <- pars + random_pars
  return(pars)
}

tick <- 0
OptimizeEquilSitePars <- function(sXc, pars=NULL, fixed=fixed, AGOfixed=FALSE,
                                  fixed8mer=FALSE, plotname_=plotname) {
  message("In the optimization function.")
  time_start <- proc.time()[3]
  tick <- 0
  if (is.null(pars) == TRUE) {
    message("About to initialize parameters:")
    print("combined:")
    print(combined)
    print(sXc)
    initial.pars <- InitializeEquilSitePars(sXc, combined=combined, fixed=fixed,
                                            fixed8mer=fixed8mer)
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
  # ind_l <- which(colnames(sXc[[1]]) == "I_combined")
  ind_l <- 1 + combined
  l <- Norm(sXc[[1]][, ind_l] + 1)*L_
  # l <- Norm(sXc[[1]][, ind_l])*L_
  # print("here)")
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
  } else if (sitelist %in% c("progthrp", "progthrp_comp", "progthrp_suppcomp")) {
    zero_grad <- c(ind_ref)
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
  # message("About to enter the optimization routine.")
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
                    control = list(maxit=1000000, factr=factr, fnscale=1))
  # Print the total amount of time the optimization took.
  time_elapsed <- proc.time()[3] - time_start
  # print(time_elapsed)
  # Return the optimized parameters.
  pars <- solution$par
  pars/log(10)
}

message("About to do the optimization:")
# PERFORM THE MAIN OPTIMIZATION.
pars.MLE <- OptimizeEquilSitePars(sXc, fixed=fixed, fixed8mer=fixed8mer)


# Pre-allocate the output matrix that will be use to calculate the confidence
# intervals of the Kds.
ncol_span <- ncol(sXc[[1]]) - 3
resampled.pars <- matrix(NaN, nrow=length(pars.MLE), ncol=ncol_span*reps,
                         dimnames=list(names(pars.MLE), 1:(ncol_span*reps)))
# Define the number of differe columns that are removed during the downsampling.
if ((mirna == "miR-124" & experiment == "equilibrium_2_tp") &
		(!("-rep1" %in% args) & !("-rep2" %in% args))) {
	max_j <- ncol(sXc[[1]]) - 6
} else if (mirna == "miR-7-24nt" & experiment == "equilibrium_tp") {
  max_j <- ncol(sXc[[1]]) - 6
} else if (mirna == "let-7a-21nt" & experiment == "equil_c_nb" & nbomitc) {
  max_j <- ncol(sXc[[1]]) - 5
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
         print(leave_out)
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
                                                      plotname_=plotname)
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
      kOutputFile <- GetAnalysisPath(
        mirna, experiment, condition=sprintf(
          "%s%s%s", str_condition, str.original, str.filt_by_original
        ), analysis_type="kds_PAPER"
      )
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

print("final means:")
print(rowMeans(resampled.pars, na.rm=TRUE))

# write.table(file=sprintf("general/temp_%s_final_new_kds.txt", mirna),
            # output, sep="\t", quote=FALSE)

print("done")
warnings()