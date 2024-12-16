################################################################################
#GenerateSiteTypeKds.py
################################################################################

source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")
print("Out of general")

# Initial parameters and constants.
args       <- commandArgs(trailingOnly=TRUE)
# args <- c("let-7a-21nt", "equil_c2_nb", "5", "5mer-m11.15")

mirna      <- args[1]
experiment <- args[2]
n_constant <- args[3]
site       <- args[4]



kExpsThreeP <- c("equil_c_nb", "equil_s_nb", "equil_sc_nb", "equil_c2_nb",
                 "equil_c_alt_nb", "equil_c2_alt_nb", "equil_sc_alt_nb")

sitelist <- "programmed"

if ("-reps" %in% args) {
  reps <- as.integer(args[which(args == "-reps") + 1])
} else {
  reps <- 20
}

if ("-nocombI" %in% args) {
  combined <- FALSE
  str.combined <- "_nocombInput"
} else {
  combined <- TRUE
  str.combined <- ""
}

if ("-buffer" %in% args) {
  buffer <- TRUE
  str.buffer <- "_buffer3p"
} else {
  buffer <- FALSE
  str.buffer <- ""
}
  fixed <- FALSE
  str.fixed <- ""
# }

if ("-L" %in% args) {
  L <- as.integer(args[which(args == "-L") + 1])
  str.L <- sprintf("_L%s", L)
} else {
  L <- 100
  str.L <- ""
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
if ("-split9mers" %in% args) {
  split9mers <- TRUE
  str.split9mers <- "_split9mers"
} else {
  split9mers <- FALSE
  str.split9mers <- ""
}
# This is for the programmed libraries, and allows for the extension out to
# 11-nt sites.
if ("-outto11mers" %in% args) {
  outto11mers <- TRUE
  str.outto11mers <- "_outto11mers"
} else {
  outto11mers <- FALSE
  str.outto11mers <- ""
}

if ("-downto2mers" %in% args) {
  downto2mers <- TRUE
  str.downto2mers <- "_downto2mers"
} else {
  downto2mers <- FALSE
  str.downto2mers <- ""
}

if ("-downto3mers" %in% args) {
  downto3mers <- TRUE
  str.downto3mers <- "_downto3mers"
} else {
  downto3mers <- FALSE
  str.downto3mers <- ""
}

if ("-downto4mers" %in% args) {
  downto4mers <- TRUE
  str.downto4mers <- "_downto4mers"
} else {
  downto4mers <- FALSE
  str.downto4mers <- ""
}


if ("-end3prand" %in% args) {
  end3prand <- TRUE
  str.end3prand <- "_end3prand"
} else {
  end3prand <- FALSE
  str.end3prand <- ""
}

if ("-new" %in% args) {
  new <- TRUE
  str.new <- "_new"
} else {
  new <- FALSE
  str.new <- ""
}



# MAIN #########################################################################
mir_list <- unlist(strsplit(mirna, ","))
exp_list  <- unlist(strsplit(experiment, ","))
mir_exp_list <- data.frame(mirna=rep(mir_list, length(exp_list)),
                           experiment=rep(exp_list,each=length(mir_list)))

kGlobalSaveRows <- c()

sXc <- apply(mir_exp_list, 1, function(row) {
  sXc_collapsed <- SitesXCounts(
    mirna=row[1], experiment=row[2], n_constant=n_constant,
    sitelist="programmed_collapsed", split9mers=split9mers,
    outto11mers=outto11mers, downto2mers=downto2mers, downto3mers=downto3mers,
    downto4mers=downto4mers, end3prand=end3prand, new=new
  )
  sXc_full <- SitesXCounts(
    mirna=row[1], experiment=row[2], n_constant=n_constant,
    sitelist="programmed", split9mers=split9mers, outto11mers=outto11mers,
    downto2mers=downto2mers, downto3mers=downto3mers, downto4mers=downto4mers,
    end3prand=end3prand, new=new
  )

  sXc_out <- sXc_collapsed

  prog_sites <- unique(gsub("^.*\\|(.*)$", replacement="\\1",
                            rownames(sXc_collapsed)))
  prog_sites <- grep("8mer-mm", prog_sites, value=TRUE)
  print(prog_sites)
  # break
  site_rows <- sprintf("^%s\\|%s$", site, prog_sites, perl=TRUE)
  for (site_row in site_rows) {
    site_row <- gsub("\\.", replacement="\\\\.", site_row)
    site_base <- unlist(strsplit(site_row, split="\\|"))[2]
    print(site_base)
    print(site_row)
    site_ind <- grep(site_row, rownames(sXc_collapsed), perl=TRUE)
    print(site_ind)
    if (length(site_ind) != 0) {
      regex_pos_sites <- sprintf("^%s\\|[0-9]*\\|%s$", site, site_base,
                                 perl=TRUE)
      site_pos_inds <- grep(regex_pos_sites, rownames(sXc_full))
      site_pos_values <- grep(regex_pos_sites, rownames(sXc_full), value=TRUE)
      sXc_temp <- rbind(sXc_out[1:(site_ind - 1), ],
                         sXc_full[site_pos_inds, ])
      sXc_out <- rbind(sXc_temp, sXc_out[(site_ind + 1):nrow(sXc_out), ])
      kGlobalSaveRows <<- c(kGlobalSaveRows, site_pos_values)      
    } 
  }
  # First remove duplicate columns depeneding on which experiment is being
  # modeled in the experiment.
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
  sXc_out
})

kGlobalInds <- sapply(kGlobalSaveRows, function(row) {
  which(rownames(sXc[[1]]) == row)
})

print(sXc[[1]][kGlobalInds, ])

if (length(sXc) > 1) {
  str.global <- "_global"
} else {
  str.global <- ""
}
str_condition <- sprintf("%s_programmed_%s", n_constant, site)
str_condition <- paste0(str_condition, str.buffer, str.combined, str.global, str.fixed,
                         str.L, str.split9mers, str.outto11mers,
                         str.downto2mers, str.downto3mers, str.downto4mers,
                         str.end3prand, str.new, "_PAPER")


print(str.end3prand)
print(str.downto4mers)
print(str.outto11mers)
print(str_condition)


# Separate site sequences from data file.
kOutputFile <- GetAnalysisPath(mirna, experiment, condition=str_condition,
                               analysis_type="kds_PAPER")
print(kOutputFile)
kSitePars <- EquilPars(mir_list[1], exp_list[1], n_constant=n_constant,
                       sitelist="programmed_collapsed", buffer=buffer,
                       combined=combined, split9mers=split9mers,
                       outto11mers=outto11mers, downto2mers=downto2mers,
                       downto3mers=downto3mers, downto4mers=downto4mers,
                       end3prand=end3prand, new=new)

print(kOutputFile)
print(kOutputFile)

InitializeFlankParsEquil <- function(sXc, combined=TRUE) {
  # l <- SubfunctionCall(GetInputEquil, sXc=sXc[[1]] + 1)
  # data <- SubfunctionCall(GetDataEquil, sXc=sXc[[1]] + 1)
  # kds <- log(Norm(l)/Norm(rowSums(data)))/2
  # kds <- kds - kds[length(kds)]
  kd_names <- paste0(rownames(sXc[[1]]), "_Kd")
  par_names <- c(kd_names, sprintf("bg_%s", mirna), sprintf("AGO_%s", mirna))
  pars <- rep(1, length(par_names))
  names(pars) <- par_names
  # Get the names of the shared parameters
  pars_shared <- intersect(par_names, rownames(kSitePars))
  pars[pars_shared] <- log(kSitePars[pars_shared, 2])
  # Get the names of the parameters that are different between the collapsed
  # list and the full list, and assign the shared parameters to the full ones.
  pars_dif_collapsed <- setdiff(rownames(kSitePars), par_names)
  pars_dif_prog <- setdiff(par_names, rownames(kSitePars))
  pars[pars_dif_prog] <- mean(log(kSitePars[pars_dif_collapsed, 2]), na.rm=TRUE)
  return(pars)
}
OptimizeEquilFlankPars <- function(sXc, pars=NULL, fixed=fixed, AGOfixed=FALSE) {
    time_start <- proc.time()[3]
    tick <<- 0
    if (is.null(pars))
      initial.pars <- SubfunctionCall(InitializeFlankParsEquil)
    else
      initial.pars <- pars
    global_pars <<- initial.pars
    L_ <- L
    if (length(sXc) > 1) {
      data <- do.call(cbind, sapply(sXc, function(sXc_i) {
        as.matrix(sXc_i[, 3:(ncol(sXc_i) - 1)] + 1)
      }))
    } else {
      data <- sXc[[1]][, 3:(ncol(sXc[[1]]) - 1)]
    }
    n_j <- sapply(sXc, function(sXc_i) {
      ncol(sXc_i) - 3
    })
    l <- Norm(sXc[[1]][, 1 + combined] + 1)*L_
    Y <- colSums(data)
    dil <- as.numeric(colnames(data))/100.0
    n_i <- nrow(data)
    # n_i <<- n_i
    n_mir <- length(sXc)
    # n_mir <<- n_mir
    data_vec <- as.numeric(as.matrix(data))
    # data_vec <<- data_vec
    lower <- log(rep(c(1e-6, 0.001, 0.1), c(n_i, n_mir, n_mir)))
    upper <- log(rep(c(1e4, 10, 10), c(n_i, n_mir, n_mir)))
    zero_grad <- setdiff(1:(n_i + 2), kGlobalInds)
    n_z <- length(zero_grad)
    norm_constant <- n_i*sum(n_j)
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
                    control = list(maxit=10000000, factr=10000, fnscale=1))

  output.pars <- solution$par
  print(proc.time()[3] - time_start)

  output.pars/log(10)
}

pars.MLE <- OptimizeEquilFlankPars(sXc, fixed=fixed)
print(10^pars.MLE[kGlobalInds])

# inds.pars <- grep("\\|", names(pars.MLE))
resampled.pars <- matrix(NaN, nrow=length(pars.MLE), ncol=5*reps,
                         dimnames=list(names(pars.MLE), 1:(5*reps)))





# Dealing with the miR-7-23nt issue of no saturation:
if (mirna == "miR-7-23nt" & experiment == "equilibrium2_nb") {
  min_col <- 4
} else {
  min_col <- 3
}

list_sXc_resamples <- list()

# leave_outs <- rep(NaN, 5*reps)
# # Iterate over the four columns
# for (j in 0:4) {
#   for (i in 1:reps) {
#     i_full <- j*reps + i
#     print(i_full)
#     tick <- 0
#     sXc.resample <- lapply(sXc, function(sXc_i) {
#       sXc_new <- apply(sXc_i, 2, function(col) {rmultinom(1, sum(col), col)})
#       rownames(sXc_new) <-rownames(sXc_i)
#       return(sXc_new)
#     })

#     print("here")
#     sXc.resample.withhold <- lapply(sXc.resample, function(sXc_i) {
#       leave_out <- sample(min_col:(ncol(sXc_i) - 1), 1)
#       leave_outs[i_full] <<- colnames(sXc_i)[leave_out]
#       sXc_i[, -leave_out]
#     })
#     list_sXc_resamples[[i_full]] <- sXc.resample.withhold[[1]]
#     resampled.pars[, i_full] <- OptimizeEquilFlankPars(sXc.resample.withhold, pars=pars.MLE)
#     resampled.pars.sort <- t(apply(resampled.pars[, 1:i_full, drop=FALSE], 1, sort))
#     output <- 10^cbind(pars.MLE,
#                        rowMeans(resampled.pars, na.rm=TRUE),
#                        resampled.pars.sort[, ceiling(0.025*i_full)],
#                        resampled.pars.sort[, ceiling(0.5*i_full)],
#                        resampled.pars.sort[, ceiling(0.975*i_full)])
#     colnames(output) <- c("Full","Mean", "Lower_CI", "Median", "Upper_CI")
#     rownames(output) <- rownames(resampled.pars)
#   }
# }
# print(kOutputFile)
# write.table(file=kOutputFile, output, sep="\t", quote = FALSE)



for (j in 0:(ncol(sXc[[1]]) - 4)) {
  message("column")
  message(j)
  for (i in 1:reps) {
    i_full <- j*reps + i
    tick <- 0
    # Resample the data:
    sXc.resample <- lapply(sXc, function(sXc_i) {
      sXc_new <- apply(sXc_i, 2, function(col) {rmultinom(1, sum(col), col)})
      rownames(sXc_new) <-rownames(sXc_i)
      return(sXc_new)
    })
    # Remove a column:
    sXc.resample.withhold <- lapply(sXc.resample, function(sXc_i) {
      leave_out <- j + 3
      sXc_i[, -leave_out]
    })
    if (length(sXc.resample.withhold) != 1) {
      sXc.resample.withhold <- sXc.resample
    }
    if (mirna == "miR-7-23nt" & experiment == "equilibrium2_nb") {
      AGOfixed <- TRUE
    } else {
      AGOfixed <- FALSE
    }
    resampled.pars[, i_full] <- OptimizeEquilFlankPars(
      sXc.resample.withhold, pars = log(10)*pars.MLE, fixed=fixed,
      AGOfixed=AGOfixed
    )
    resampled.pars.sort <- t(apply(resampled.pars[, 1:i_full, drop=FALSE], 1,
                                   sort))
    output <- 10^cbind(pars.MLE,
                       rowMeans(resampled.pars, na.rm=TRUE),
                       resampled.pars.sort[, ceiling(0.025*i_full)],
                       resampled.pars.sort[, ceiling(0.5*i_full)],
                       resampled.pars.sort[, ceiling(0.975*i_full)])
    colnames(output) <- c("Full","Mean", "Lower_CI", "Median", "Upper_CI")
    sapply(mir_list, function(mirna) {
      write.table(file=kOutputFile, output[kGlobalInds, ], sep="\t",
                  quote=FALSE)
    })
  }
}


