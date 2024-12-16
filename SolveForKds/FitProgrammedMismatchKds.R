################################################################################
#GenerateSiteTypeKds.py
################################################################################

source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")
print("Out of general")

# Initial parameters and constants.
args       <- commandArgs(trailingOnly=TRUE)
# args <- c("let-7a-21nt", "equil_c2_nb", "5", "8mer-mmA2", "0")

mirna      <- args[1]
experiment <- args[2]
n_constant <- args[3]
site       <- args[4]
bayes_iter <- args[5]
# lambda0    <- args[6]



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

if ("-collapsemm" %in% args) {
  collapsemm <- TRUE
  str.collapsemm <- "_collapsemm"
} else {
  collapsemm <- FALSE
  str.collapsemm <- ""
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
    sitelist="programmed_suppcomp"
  )
  sXc_full <- SitesXCounts(
    mirna=row[1], experiment=row[2], n_constant=n_constant,
    sitelist="programmed"
  )
  sXc_collapsed <<- sXc_collapsed
  sXc_full <<- sXc_full
  sXc_out <- sXc_collapsed

  site_rows <- grep(sprintf("^%s\\|.*\\|Comp$", site),
                    rownames(sXc_collapsed), perl=TRUE, value=TRUE)
  for (site_row in site_rows) {
    site_ind <- which(rownames(sXc_out) == site_row)
    site_row <- gsub("\\.", replacement="\\\\.", site_row)
    site_row <- gsub("\\|", replacement="\\\\|", site_row)
    site_row <- gsub("Comp", replacement="8mer-mm[ACTG][2-7]", site_row)
    site_mm_inds <- grep(site_row, rownames(sXc_full), perl=TRUE)
    site_mm_values <- grep(site_row, rownames(sXc_full), perl=TRUE, value=TRUE)
    sXc_temp <- rbind(sXc_out[1:(site_ind - 1), ],
                       sXc_full[site_mm_inds, ])
    sXc_out <- rbind(sXc_temp, sXc_out[(site_ind + 1):nrow(sXc_out), ])
    kGlobalSaveRows <<- c(kGlobalSaveRows, site_mm_values)
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

# Define the global indeces that will be used within the C script to subtract
# The relevant values of the bayesian prior.
kGlobalInds <- sapply(kGlobalSaveRows, function(row) {
  which(rownames(sXc[[1]]) == row)
})

if (length(sXc) > 1) {
  str.global <- "_global"
} else {
  str.global <- ""
}

str_condition <- sprintf("%s_%s_%s_bayes%s%s%s%s%s_PAPER", n_constant, "programmedmm", 
                         site, bayes_iter, str.buffer, str.combined, str.global, str.fixed,
                         str.L, str.collapsemm)

# Allocate the output file.
kOutputFile <- GetAnalysisPath(mirna, experiment, condition=str_condition,
                               analysis_type="kds_PAPER")

# Load the parameters fit with the optimization scheme in which the mismatched
# sites are all collapsed together.
kSitePars <- EquilPars(mir_list[1], exp_list[1], n_constant=n_constant,
                       sitelist="programmed_suppcomp", buffer=buffer,
                       combined=combined, collapsemm=collapsemm)



InitializeFlankParsEquil <- function(sXc, combined=TRUE) {
  l <- SubfunctionCall(GetInputEquil, sXc=sXc[[1]] + 1)
  data <- SubfunctionCall(GetDataEquil, sXc=sXc[[1]] + 1)
  kds_init <- log(Norm(l)/Norm(rowSums(data)))
  # kds <- kds - kds[length(kds)]
  kd_names <- paste0(rownames(sXc[[1]]), "_Kd")
  names(kds_init) <- kd_names
  par_names <- c(kd_names, sprintf("bg_%s", mirna), sprintf("AGO_%s", mirna))
  # Pre-assign a 1 to each parameter
  pars <- rep(1, length(par_names))
  names(pars) <- par_names
  # Get the names of the shared parameters
  pars_shared <- intersect(par_names, rownames(kSitePars))
  pars[pars_shared] <- log(kSitePars[pars_shared, 2])
  kds_shared <- pars_shared[1:(length(pars_shared) - 2)]
  kds_not_shared <- setdiff(par_names, rownames(kSitePars))

  # Fit a linear model to the relationship between the initialized kd values
  # and the parameters fit in the fitting with the collapsed mismatch sites.
  lmfit <- lm(pars[kds_shared] ~ kds_init[kds_shared])
  m <- lmfit$coefficients[2]
  b <- lmfit$coefficients[1]
  # Initialize the separated site Kds according to this relationship.
  pars[kds_not_shared] <- kds_init[kds_not_shared]*m + b
  pars_dif_collapsed <- setdiff(rownames(kSitePars), par_names)

  # for (site_dif in pars_dif_collapsed) {
  #   site_use <- gsub("\\.", replacement="\\\\.", site_dif)
  #   site_use <- gsub("\\|", replacement="\\\\|", site_use)
  #   site_use <- gsub("Comp", replacement="8mer-mm[ACTG][2-7]", site_use)
  #   pars[grep(site_use, names(pars), perl=TRUE)] <- log(kSitePars[site_dif, 2])
  # }
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
  zero_grad <- setdiff(1:(n_i), kGlobalInds)
  # zero_grad <- c()
  reg_inds <- kGlobalInds
  reg_vals <- initial.pars[kGlobalInds]
  # reg_vals <- rep(0, length(reg_vals))
  reg_vals <<- reg_vals
  n_z <- length(zero_grad)
  n_reg <- length(reg_inds)
  norm_constant <- n_i*sum(n_j)
  # reg_sd <- 0
  reg_mean <- 0
  lambda <- 0
  lambda0 <- 0

  solution <- optim(initial.pars,
                  CostCPrior,
                  gr = GradCPrior,
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
                  reg_mean=reg_mean,
                  lambda=lambda,
                  lambda0=lambda0,
                  reg_inds=reg_inds,
                  reg_vals=reg_vals,
                  n_reg=n_reg,
                  upper_ = upper,
                  lower_ = lower,
                  plot_ = plot_,
                  plotname = plotname_,
                  method = "L-BFGS-B",
                  lower=lower,
                  upper=upper,
                  control = list(maxit=10000000, fnscale=1))
  output.pars <- solution$par
  reg_mean <- mean(output.pars[kGlobalInds] - reg_vals)
  reg_sd <- sd(output.pars[kGlobalInds] - reg_vals)
  if (bayes_iter != 0) {
    for (rep in 1:bayes_iter) {
      print(rep)
      # Update lambda and lambda0.
      lambda <- 1/(2*reg_sd^2)
      lambda0 <- log(reg_sd*sqrt(2*pi))
      # Perform the new optimization.
      solution <- optim(output.pars,
                        CostCPrior,
                        gr = GradCPrior,
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
                        reg_mean=reg_mean,
                        lambda=lambda,
                        lambda0=lambda0,
                        reg_inds=reg_inds,
                        reg_vals=reg_vals,
                        n_reg=n_reg,
                        upper_ = upper,
                        lower_ = lower,
                        plot_ = plot_,
                        plotname = plotname_,
                        method = "L-BFGS-B",
                        lower=lower,
                        upper=upper,
                        control = list(maxit=10000000, fnscale=1))
      # Update the parameters, and the regularization mean and standard
      # deviation.
      output.pars <- solution$par
      reg_mean <- mean(output.pars[kGlobalInds] - reg_vals)
      reg_sd <- sd(output.pars[kGlobalInds] - reg_vals)
    }
  }
  message("out of optimization")
  print(proc.time()[3] - time_start)
  # Return the parameters.
  output.pars/log(10)
}

pars.MLE <- OptimizeEquilFlankPars(sXc, fixed=fixed)

# Pre-allocate the output matrix.
resampled.pars <- matrix(NaN, nrow=length(pars.MLE), ncol=5*reps,
                         dimnames=list(names(pars.MLE), 1:(5*reps)))

# Dealing with the miR-7-23nt issue of no saturation:
if (mirna == "miR-7-23nt" & experiment == "equilibrium2_nb") {
  min_col <- 4
} else {
  min_col <- 3
}

list_sXc_resamples <- list()


for (j in 0:(ncol(sXc[[1]]) - 4)) {
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


