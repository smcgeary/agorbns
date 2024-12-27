################################################################################
#FitFlankKds.R
################################################################################

source("general/general.R")
source("general/ModelingFunctions.R")
print("Out of general")

# Initial parameters and constants.
args       <- commandArgs(trailingOnly=TRUE)
mirna      <- args[1]
experiment <- args[2]
n_constant <- args[3]
sitelist   <- args[4]
site       <- args[5]

# mirna <- "miR-124"
# experiment <- "equilibrium"
# n_constant <- 5
# sitelist <- "paper"
# site <- "AA-8mer-bT6"

# reps <- 20
# combined <- TRUE
# str.combined <- ""
# fixed <- FALSE
# str.fixed <- ""
# L <- 100
# str.L <- ""


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

if ("-buffer3p" %in% args) {
  buffer <- TRUE
  str.buffer <- "_buffer3p"
} else {
  buffer <- FALSE
  str.buffer <- ""
}

# if ("-fixed" %in% args & mirna %in% c("miR-7-23nt", "miR-7-24nt", "miR-7-25nt")
#     & experiment == "equilibrium2_nb") {
#   fixed <- TRUE
#   str.fixed <- "_fixed"
# } else {
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

# MAIN #########################################################################
mir_list <- unlist(strsplit(mirna, ","))
exp_list  <- unlist(strsplit(experiment, ","))
mir_exp_list <- data.frame(mirna=rep(mir_list, length(exp_list)),
                           experiment=rep(exp_list,each=length(mir_list)))


sXc <- apply(mir_exp_list, 1, function(row) {
  SitesAndSingleFlanksXCounts(row[1], site, experiment=row[2], n_constant=n_constant,
                              sitelist=sitelist, buffer=buffer)
})

if (length(sXc) > 1) {
  str.global <- "_global"
} else {
  str.global <- ""
}
str_condition <- sprintf("%s_%s_%s%s%s%s%s%s_PAPER", n_constant, sitelist, 
                         site, str.buffer, str.combined, str.global, str.fixed,
                         str.L)
print(str_condition)


# Separate site sequences from data file.
kOutputFile <- GetAnalysisPath(mirna, experiment,
                               condition=sprintf("%s_%s_%s%s%s_PAPER",
                                                 n_constant, sitelist, site,
                                                 str.buffer, str.combined),
                               analysis_type="kds_PAPER")
print(kOutputFile)
site.pars <- EquilPars(mir_list[1], exp_list[1], n_constant=n_constant,
                       sitelist=sitelist, buffer=buffer, combined=combined)



InitializeFlankParsEquil <- function(sXc, combined=TRUE) {
  n_mir = length(sXc)
  l <- SubfunctionCall(GetInputEquil, sXc=sXc[[1]] + 1)
  data <- SubfunctionCall(GetDataEquil, sXc=sXc[[1]] + 1)
  kds <- log(Norm(l)/Norm(rowSums(data)))/2
  kds <- kds - kds[length(kds)]
  names(kds) <- paste0(rownames(sXc[[1]]), "_Kd")
  pars <- c(kds, bg=-4, AGO=1)
  names(pars)[length(pars)-1] <- sprintf('bg_%s', mirna)
  names(pars)[length(pars)] <- sprintf('AGO_%s', mirna)
  inds_shared <- intersect(rownames(site.pars), names(pars))
  pars[inds_shared] <- log(site.pars[inds_shared,]$Mean)
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
    # zero_grad <- which(names(initial.pars) %in% rownames(site.pars))
    # zero_grad <- c(n_i)
    norm_constant <- n_i*sum(n_j)

    if (fixed) {
      zero_grad <- n_i + c(0, 1, 2)
    } else if (AGOfixed) {
      zero_grad <- n_i + c(0, 1, 2)
    } else {
      zero_grad <- c(n_i)
    }
    n_z <- length(zero_grad)
    # pars_global <<- initial.pars
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

print(sXc)
pars.MLE <- OptimizeEquilFlankPars(sXc, fixed=fixed)
print(10^pars.MLE)

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
    # If not 12mers, remove a column:
    if (sitelist != "12mers") {
      # if (mirna == "miR-7-23nt" & experiment == "equilibrium2_nb") {
      #   min_col <- 4
      # } else {
      #   min_col <- 3
      # }
      sXc.resample.withhold <- lapply(sXc.resample, function(sXc_i) {
        leave_out <- j + 3
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
    if (i == 1) {
      print(colnames(sXc.resample.withhold[[1]]))    
    }
    if (mirna == "miR-7-23nt" & experiment == "equilibrium2_nb") {
      AGOfixed <- TRUE
    } else {
      AGOfixed <- FALSE
    }
    print(sXc.resample.withhold)
    resampled.pars[, i_full] <- OptimizeEquilFlankPars(sXc.resample.withhold,
                                                      pars = log(10)*pars.MLE,
                                                      fixed=fixed, AGOfixed=AGOfixed)
    resampled.pars.sort <- t(apply(resampled.pars[, 1:i_full, drop=FALSE], 1,
                                   sort))
    output <- 10^cbind(pars.MLE,
                       rowMeans(resampled.pars, na.rm=TRUE),
                       resampled.pars.sort[, ceiling(0.025*i_full)],
                       resampled.pars.sort[, ceiling(0.5*i_full)],
                       resampled.pars.sort[, ceiling(0.975*i_full)])
    colnames(output) <- c("Full","Mean", "Lower_CI", "Median", "Upper_CI")
    sapply(mir_list, function(mirna) {
      # kOutputFile <- GetAnalysisPath(mirna, experiment, condition=str_condition,
      #                          analysis_type="kds_PAPER")
      print(kOutputFile)
      write.table(file=kOutputFile, output, sep="\t", quote = FALSE)
    })
    if (sitelist == "12mers") {
      kds.mean <- log(10^rowMeans(resampled.pars, na.rm=TRUE))
      # print(kds.mean[1:5])
      write.table(file=kOutputFileMean, kds.mean, sep="\t", quote=FALSE,
                  row.names=FALSE)
    }
  }
}


