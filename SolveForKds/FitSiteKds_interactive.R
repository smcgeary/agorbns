################################################################################
#GenerateSiteKds_PAPER.py
################################################################################

source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")
print("Out of general")

# Initial parameters and constants.
# args       <- commandArgs(trailingOnly=TRUE)
args <- c("miR-1", "equilibrium", 5, "paperfinal", "-buffer", "-nocombI")
# args <- c("miR-124", "equilibrium", 5, "paperfinal", "-compcorrect")
# args <- c("miR-124", "equilibrium", 5, "paperfinal")

mirna      <- args[1]
experiment <- args[2]
n_constant <- args[3]
sitelist   <- args[4]

print(args)



if ("-reps" %in% args) {
  reps <- as.integer(args[which(args == "-reps") + 1])
} else {
  reps <- 40
}

if (sitelist == "12mers") {
  mirna.start = as.integer(args[which(args == "-mir_start") + 1])
  str.mir_start = sprintf("_%s-%s", mirna.start, mirna.start + 3)
} else {
  mirna.start = FALSE
  str.mir_start = ""
}

if ("-uniq" %in% args) {
  uniq <- TRUE
  str.unique <- "_uniq"
} else {
  uniq <- FALSE
  str.unique <- ""
}

if ("-buffer" %in% args) {
  buffer <- TRUE
  str.buffer <- "_buffer3p"
} else {
  buffer <- FALSE
  str.buffer <- ""
}

if ("-nocombI" %in% args) {
  combined <- FALSE
  str.combined <- "_nocombInput"
} else {
  combined <- TRUE
  str.combined <- ""
}

if ("-fixed" %in% args & mirna %in% c("miR-7-23nt", "miR-7-24nt", "miR-7-25nt")
    & experiment == "equilibrium2_nb") {
  fixed <- TRUE
  altfixed <- FALSE
  str.fixed <- "_fixed"
} else if ("-altfixed" %in% args &
           mirna %in% c("miR-7-23nt", "miR-7-24nt", "miR-7-25nt") &
           experiment == "equilibrium2_nb") {
  fixed <- TRUE
  altfixed <- TRUE
  str.fixed <- "_altfixed"
} else {
  altfixed <- FALSE
  fixed <- FALSE
  str.fixed <- ""
}

if ("-L" %in% args) {
  L <- as.integer(args[which(args == "-L") + 1])
  str.L <- sprintf("_L%s", L)
} else {
  L <- 100
  str.L <- ""
}

if ("-single" %in% args) {
  single <- TRUE
  str.single <- "_singleonly"
} else {
  single <- FALSE
  str.single <- ""
}

if ("-compcorrect" %in% args) {
  compcorrect <- TRUE
  str.compcorrect <- "_compcorrect"
} else {
  compcorrect <- FALSE
  str.compcorrect <- ""
}

if ("-wobble" %in% args) {
  wobble <- TRUE
  str.wobble <- "_wobble"
} else {
  wobble <- FALSE
  str.wobble <- ""
}

if (mirna == "miR-7-23nt" & experiment == "equilibrium2_nb") {
  AGOfixed <- TRUE
} else {
  AGOfixed <- FALSE
}


# MAIN #########################################################################
# 1. Get data table:
# Make two vectors, one of miRNA names and one of experiment names, and combine
# them into one data frame.
mir_list <- unlist(strsplit(mirna, ","))
exp_list  <- unlist(strsplit(experiment, ","))
mir_exp_list <- data.frame(mirna=rep(mir_list, length(exp_list)),
                           experiment=rep(exp_list,each=length(mir_list)))

# Make the list of site-by-count matricies corresponding each row in the 
# `mir_exp_list` dataframe.
sXc <- apply(mir_exp_list, 1, function(row) {
  mirna <- row[1]
  sXc_i <- SitesXCounts(row[1], experiment=row[2], n_constant=n_constant,
                        sitelist=sitelist, mirna.start=mirna.start,
                        uniq=uniq, buffer=buffer)
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
  # Portion of code where miR-124 sites are corrected.
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
    # sXr <<- sXr
    # sXr_new <<- sXr_new

    sXc_out_new <- sXc_out
    sXc_out_new[, vec_cond_rep] <- (
      sXc_out[, vec_cond_rep]/sXr[, vec_cond_rep]*sXr_new[, vec_cond_rep]
    )
    sXc_out <- sXc_out_new
  }
  sXc_out
})



if (length(sXc) > 1) {
  str.global <- "_global"
} else {
  str.global <- ""
}
str_condition <- sprintf("%s_%s%s%s%s%s%s%s%s%s%s%s_PAPER", n_constant, sitelist,
                         str.unique, str.buffer, str.mir_start, str.combined,
                         str.global, str.fixed, str.L, str.single,
                         str.compcorrect, str.wobble)

kOutputFileMean <- GetAnalysisPath(
  mir_list[[1]], experiment, condition=sprintf("%s_logmean", str_condition),
  analysis_type="kds_PAPER"
)

if (sitelist == "12mers") {
  pars_sd <- 0.01
  plot_ <- TRUE
  plotname <- paste0("12merfits/", paste(mirna, experiment, n_constant,
                                         sitelist, mirna.start,
                     sep="_"), ".pdf")
} else {
  pars_sd <- 0.1
  plot_ <- FALSE
  plotname <- NULL
}


InitializeEquilSitePars <- function(sXc, combined=TRUE, fixed=FALSE) {
  n_mir = length(sXc)
  l <- SubfunctionCall(GetInputEquil, sXc=sXc[[1]])
  data <- SubfunctionCall(GetDataEquil, sXc=sXc[[1]])
  kds <- log(Norm(l)/Norm(rowSums(data)))
  kds <- kds - kds[length(kds)]
  kds <- kds/2
  names(kds) <- paste0(rownames(sXc[[1]]), "_Kd")
  
  if (fixed) {
    if (altfixed) {
      combined_fixed <- TRUE
    } else {
      combined_fixed <- FALSE
    }
    pars <- EquilPars(mirna, experiment=experiment, n_constant=n_constant,
                      sitelist=sitelist, uniq=uniq, combined=combined_fixed, global=TRUE)
    bgs <- log(pars[sprintf("bg_%s", mirna),]$Mean)
    As <- log(pars[sprintf("AGO_%s", mirna),]$Mean)
  } else {
    bgs <- rep(log(0.1), n_mir)
    As <- rep(log(20), n_mir)    
  }
  names(bgs) <- sprintf("bg_%s", mir_list)
  names(As) <- sprintf("AGO_%s", mir_list)
  pars <- c(kds, bgs, As)
  random_pars <- rnorm(length(pars), 0, pars_sd)
  if (fixed) {
    random_pars[length(pars) - 1] = 0
    random_pars[length(pars)    ] = 0
  }
  names(random_pars) <- names(pars)
  pars <- pars + random_pars
  pars["None_Kd"] <- 0
  return(pars)
}

tick <- 0
OptimizeEquilSitePars <- function(sXc, pars=NULL, fixed=fixed, AGOfixed=FALSE,
                                  plotname_=plotname) {
  time_start <- proc.time()[3]
  tick <- 0
  if (is.null(pars) == TRUE) {
    initial.pars <- InitializeEquilSitePars(sXc, combined=combined, fixed=fixed)
  } else {
    initial.pars <- pars
  }
    L_ <- L
  if (length(sXc) > 1) {
    data <- do.call(cbind, sapply(sXc, function(sXc_i) {
      as.matrix(sXc_i[, 3:(ncol(sXc_i) - 1)])
    }))
  } else {
    data <- sXc[[1]][, 3:(ncol(sXc[[1]]) - 1)]
  }
  n_j <- sapply(sXc, function(sXc_i) {
    ncol(sXc_i) - 3
  })
  l <- Norm(sXc[[1]][, 1 + combined])*L_
  Y <- colSums(data)
  dil <- as.numeric(colnames(data))/100.0
  n_i <- nrow(data)
  n_mir <- length(sXc)
  data_vec <- as.numeric(as.matrix(data))
  lower <- log(rep(c(1e-6, 0.001, 0.1), c(n_i, n_mir, n_mir)))
  upper <- log(rep(c(1e4, 10, 10), c(n_i, n_mir, n_mir)))

  if (fixed) {
    zero_grad <- n_i + c(0, 1, 2)
  } else if (AGOfixed) {
    zero_grad <- n_i + c(0, 1, 2)
  } else {
    zero_grad <- c(n_i)
  }
  n_z <- length(zero_grad)
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
                    upper_ = upper,
                    lower_ = lower,
                    plot_ = plot_,
                    plotname = plotname_,
                    method = "L-BFGS-B",
                    lower=lower,
                    upper=upper,
                    control = list(maxit=10000000, factr=10000, fnscale=1))
  output.pars <- solution$par
  time_elapsed <- time_start - proc.time()[3]
  # print(time_elapsed)
  output.pars/log(10)
}

pars.MLE <- OptimizeEquilSitePars(sXc, fixed=fixed)

print(pars.MLE)

break

# print(sprintf("%.16f", 10^pars.MLE))

ncol_span = ncol(sXc[[1]]) - 2

resampled.pars <- matrix(NaN, nrow=length(pars.MLE), ncol=ncol_span*reps,
                         dimnames=list(names(pars.MLE), 1:(ncol_span*reps)))
print("HERE")

for (j in 0:(ncol(sXc[[1]]) - 4)) {
  message("column")
  message(j)
  if (sitelist != "12mers") {
    # if (mirna == "miR-7-23nt" & experiment == "equilibrium2_nb") {
    #   min_col <- 4
    # } else {
    #   min_col <- 3
    # }
    sXc.withhold <- lapply(sXc, function(sXc_i) {
      leave_out <- j + 3
      print(leave_out)
      sXc_i[, -leave_out]
    })
  } else {
    sXc.withhold <- sXc
    plotname <- paste0("12merfits/", paste(mirna, experiment, n_constant, sitelist, mir_start,
                   sep="_"), "_rep", j*reps + i, ".pdf")
  }
  if (length(sXc.withhold) != 1) {
    sXc.withhold <- sXc
  }
  for (i in 1:reps) {
    i_full <- j*reps + i
    tick <- 0
    # Resample the data:
    sXc.resample.withhold <- lapply(sXc.withhold, function(sXc_i) {
      sXc_new <- apply(sXc_i, 2, function(col) {rmultinom(1, sum(col), col)})
      rownames(sXc_new) <-rownames(sXc_i)
      return(sXc_new)
    })
    # If not 12mers, remove a column:
    resampled.pars[, i_full] <- OptimizeEquilSitePars(sXc.resample.withhold,
                                                      pars = log(10)*pars.MLE,
                                                      fixed=fixed, AGOfixed=AGOfixed,
                                                      plotname_=plotname)
    resampled.pars.sort <- t(apply(resampled.pars[, 1:i_full, drop=FALSE], 1,
                                   sort))
    output <- 10^cbind(pars.MLE,
                       rowMeans(resampled.pars, na.rm=TRUE),
                       resampled.pars.sort[, ceiling(0.025*i_full)],
                       resampled.pars.sort[, ceiling(0.5*i_full)],
                       resampled.pars.sort[, ceiling(0.975*i_full)])
    colnames(output) <- c("Full","Mean", "Lower_CI", "Median", "Upper_CI")
    sapply(mir_list, function(mirna) {
      kOutputFile <- GetAnalysisPath(mirna, experiment, condition=str_condition,
                               analysis_type="kds_PAPER")
      # write.table(file=kOutputFile, output, sep="\t", quote = FALSE)
    })
    if (sitelist == "12mers") {
      kds.mean <- log(10^rowMeans(resampled.pars, na.rm=TRUE))
      # print(kds.mean[1:5])
      # write.table(file=kOutputFileMean, kds.mean, sep="\t", quote=FALSE,
      #             row.names=FALSE)
    }
  }
}

temp <- EquilPars("miR-124", "equilibrium", 5, "paperfinal")

plot(unlist(temp), unlist(output), log='xy')


print("done")
print(ceiling(0.025*i_full))
print(ceiling(0.975*i_full))
# write.table(file=sprintf("%s_temp_sampleout.txt", mirna), cbind(pars.MLE, resampled.pars))
print(sXc)
print(output)
warnings()





