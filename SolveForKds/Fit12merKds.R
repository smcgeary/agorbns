################################################################################
#GenerateSiteKds_PAPER.py
################################################################################

## THIS IS THE OFFICIAL SCRIPT FOR FINALIZING THE MANUSCRIPT.
## THIS SCRIPT HAS THESE DECISIONS MADE:
## 1.) PROTEIN IS FIT AS A PARAMETER
## 2.) A SINGLE BACKGROUND PARAMETER IS GIVEN TO ALL FIVE CONCENTRATIONS.
## 3.) LIST, START, AND STOP POSITION ARE MANDATORY PARAMETERS.
# ## 4.) ALL INPUTS ARE COMBINED ACROSS FOR THESE ANALYSES.
source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")
print("Out of general")
source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/ModelingFunctions.R")
print("out of modeling functions")
library(numDeriv)

# Initial parameters and constants.
args       <- commandArgs(trailingOnly=TRUE)
mirna      <- args[1]
experiment <- args[2]
n_constant <- args[3]
sitelist   <- args[4]

# mirna <- "miR-1"
# experiment <- "equilibrium"
# n_constant <- 5

# sitelist <- "12mers"



if ("-reps" %in% args) {
  reps <- as.integer(args[which(args == "-reps") + 1])
} else {
  reps <- 20
}
if (sitelist == "12mers") {
  mir_start = as.integer(args[which(args == "-mir_start") + 1])
  mir_start <- 1
  str.mir_start = sprintf("_%s-%s", mir_start, mir_start + 3)
} else {
  mir_start = FALSE
  str.mir_start = ""
}

str.mir_start

if ("-fixed" %in% args & mirna %in% c("miR-7-23nt", "miR-7-24nt", "miR-7-25nt")
    & experiment == "equilibrium2_nb") {
  fixed <- TRUE
  str.fixed <- "_fixed"
} else {
  fixed <- FALSE
  str.fixed <- ""
}


if ("-pyclasses" %in% args) {
  pyclasses <- TRUE
  str.classes <- "_classes"
} else {
  pyclasses <- FALSE
  str.classes <- ""    
}
if ("-nocombI" %in% args) {
  combined <- FALSE
  str.combined <- "_nocombInput"
} else {
  combined <- TRUE
  str.combined <- ""
}

print(mirna)
# MAIN #########################################################################
# 1. Get data table:
sXc <- SitesXCounts(mirna, experiment, n_constant, sitelist,
                    mirna.start=mir_start)
print(head(sXc))
# sXc_names_collapsed <-  unique(substr(rownames(sXc)[-nrow(sXc)], 3, 10))
# sXc_alt <- matrix(NaN, nrow=length(sXc_names_collapsed) + 1, ncol=ncol(sXc),
#                   dimnames=list(c(sXc_names_collapsed, "None"), colnames(sXc)))
# for (name in sXc_names_collapsed) {
#   sXc_alt[name,] <- colSums(sXc[grep(paste0("^.{2}",
#                                                       name,
#                                                       ".{2}$"),rownames(sXc)),])
# }

# sXc_alt["None", ] <- unlist(sXc["None", ])
# print(head(sXc_alt))
# print(tail(sXc_alt))
# print(colSums(sXc))
# print(colSums(sXc_alt))
# print(dim(sXc_alt))

kOutputFile <- GetAnalysisPath(mirna, experiment,
                               condition=sprintf("%s_%s%s%s_PAPER", n_constant,
                                                 sitelist, str.mir_start, 
                                                 str.combined),
                               analysis_type="kds_PAPER")
kOutputFileFull <- GetAnalysisPath(mirna, experiment,
                               condition=sprintf("%s_%s%s%s_PAPER_full", n_constant,
                                                 sitelist, str.mir_start, 
                                                 str.combined),
                               analysis_type="kds_PAPER")
kOutputFileMean <- GetAnalysisPath(mirna, experiment,
                               condition=sprintf("%s_%s%s%s_PAPER_mean", n_constant,
                                                 sitelist, str.mir_start, 
                                                 str.combined),
                               analysis_type="kds_PAPER")



print(kOutputFile)
print(kOutputFileFull)
print(kOutputFileMean)
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
if (sitelist == "12mers") {
  pars_sd <- 0.01
  plot_ <- TRUE
  plot_dir <- paste0("12merfits/", paste(mirna, experiment, n_constant,
                                       sitelist, mir_start, str.combined,
                                       sep="_"))
} else {
  pars_sd <- 3
  plot_ <- TRUE
  plot_dir <- paste0("12merfits/", paste(mirna, experiment, n_constant,
                                       sitelist, str.combined, sep="_"))
}
plotname <- paste0(plot_dir, "/full.pdf")
tempname <- paste0(plot_dir, "full_temp.txt")
if (!file.exists(plot_dir)) {
  dir.create(plot_dir)
}


InitializeEquilSitePars <- function(sXc, combined=TRUE) {
  if (sitelist=="12mers" & experiment=="equilibrium2_nb") {
    L_ <- 300
  } else {
    L_ <- 100
  }
  l <- Norm(sXc[, 1 + combined] + 1)*L_
  data <- SubfunctionCall(GetDataEquil)
  kds <- log10(Norm(l)/Norm(rowSums(data)))
  kds <- kds - kds[length(kds)]
  kds <- kds
  names(kds) <- paste0(rownames(sXc), "_Kd")
  pars <- c(kds, bg=-1, AGO=1)
  random_pars <- rnorm(length(pars), 0, pars_sd)
  names(random_pars) <- names(pars)
  pars <- pars + random_pars
  pars["None_Kd"] <- 0
  return(pars)
}



tick <- 0
# pars <- InitializeEquilSitePars(sXc, combined=combined)*log(10)
OptimizeEquilSitePars <- function(sXc, pars=NULL, plotname_=plotname,
                                  tempname_=tempname) {
  time_start <- proc.time()[3]
  tick <- 0
  if (is.null(pars))
    initial.pars <- InitializeEquilSitePars(sXc, combined=combined)*log(10)
  else
    initial.pars <- pars*log(10)
  if (sitelist=="12mers" & experiment=="equilibrium2_nb") {
    L_ <- 300
  } else {
    L_ <- 100
  }
  data <- as.matrix(sXc[, 3:(ncol(sXc) - 1)])
  n_j <- ncol(sXc) - 3
  l <- Norm(sXc[, 1 + combined] + 1)*L_
  Y <- colSums(data)
  dil <- as.numeric(colnames(data))/100.0
  n_i <- nrow(data)
  n_mir <- 1
  data_vec <- as.numeric(as.matrix(data))
  lower <- log(rep(c(1e-6, 0.001, 0.1), c(n_i, n_mir, n_mir)))
  upper <- log(rep(c(1e4, 10, 10), c(n_i, n_mir, n_mir)))
  zero_grad <- c(n_i)
  n_z <- length(zero_grad)
  # print(head(data))
  # print(n_j)
  # print(Y)
  # print(dil)
  # print(n_i)
  # print(n_mir)
  # print(zero_grad)
  # print(n_z)
  # print(fixed)
  # print(head(l))
  # print(sum(l))
  # print(L_)
  # print(print(initial.pars[is.na(initial.pars)]))
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


  output.pars <- solution$par/log(10)
  print(proc.time()[3] - time_start)
  output.pars
}
pars.MLE <- OptimizeEquilSitePars(sXc)
print(min(pars.MLE))
print(max(pars.MLE))
if (sitelist == "12mers") {
  kds.full <- pars.MLE
  print(kds.full[1:5])
  names(kds.full) <- rownames(sXc)
  write.table(file=kOutputFileFull, t(kds.full), sep="\t", quote=FALSE,
              row.names=FALSE)
}
resampled.pars <- matrix(NaN, nrow=length(pars.MLE), ncol=5*reps,
                         dimnames=list(names(pars.MLE), 1:(5*reps)))

for (j in 0:(ncol(sXc) - 4)) {
  print(j)
  for (i in 1:reps) {
    i_full <- j*reps + i
    # print(i_full)
    tick <- 0
    sXc.resample <- apply(sXc, 2, function(col) {rmultinom(1, sum(col), col)})+1
    rownames(sXc.resample) <- rownames(sXc)
    if (sitelist != "12mers") {
      sXc.resample.withhold <- sXc.resample[, -(j + 3)]

    } else {
      sXc.resample.withhold <- sXc.resample
      plotname <- paste0(plot_dir, "/", j*reps + i, ".pdf")
      tempname <- paste0(plot_dir, j*reps + i, "_temp.txt")

    }
    resampled.pars[, i_full] <- OptimizeEquilSitePars(sXc.resample.withhold, pars.MLE,
                                                      plotname_=plotname, tempname_=tempname)
    resampled.pars.sort <- t(apply(resampled.pars[, 1:i_full, drop=FALSE], 1, sort))
    output <- 10^cbind(pars.MLE,
                     rowMeans(resampled.pars, na.rm=TRUE),
                     resampled.pars.sort[, ceiling(0.025*i_full)],
                     resampled.pars.sort[, ceiling(0.5*i_full)],
                     resampled.pars.sort[, ceiling(0.975*i_full)])
    colnames(output) <- c("Full","Mean", "Lower_CI", "Median", "Upper_CI")
    write.table(file=kOutputFile, output, sep="\t", quote = FALSE)
    if (sitelist == "12mers") {
      kds.mean <- rowMeans(resampled.pars, na.rm=TRUE)
      # print(kds.mean[1:5])
      write.table(file=kOutputFileMean, t(kds.mean), sep="\t", quote=FALSE,
                 row.names=FALSE)
    }
  }
}
print(output)
print(kOutputFile)
warnings()





