################################################################################
#GenerateSiteKds_PAPER.py
################################################################################

## THIS IS THE OFFICIAL SCRIPT FOR FINALIZING THE MANUSCRIPT.
## THIS SCRIPT HAS THESE DECISIONS MADE:
## 1.) PROTEIN IS FIT AS A PARAMETER
## 2.) A SINGLE BACKGROUND PARAMETER IS GIVEN TO ALL FIVE CONCENTRATIONS.
## 3.) LIST, START, AND STOP POSITION ARE MANDATORY PARAMETERS.
# ## 4.) ALL INPUTS ARE CO5BINED ACROSS FOR THESE ANALYSES.
source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")
print("Out of general")
source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/ModelingFunctions.R")
print("out of modeling functions")
library(numDeriv)

# Initial parameters and constants.
mirna      <- "miR-7-23nt"
experiment <- "equilibrium2_nb"
n_constant <- 5
sitelist   <- "12mers"
mir_start <- 1
reps <- 20
str.mir_start = sprintf("_%s-%s", mir_start, mir_start + 3)

combined <- FALSE
if (combined) {
  str.combined <- ""
} else {
  str.combined <- "_nocombInput"
}

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

if (!file.exists(plot_dir)) {
  dir.create(plot_dir)
}

plot_dir <- paste0("12merfits/", paste(mirna, experiment, n_constant,
                                       sitelist, mir_start, str.combined,
                                       sep="_"))
plotname <- paste0(plot_dir, "/full.pdf")

tempname <- paste0(plot_dir, "full_temp.txt")


InitializeEquilSitePars <- function(sXc, combined=TRUE) {
  l <- SubfunctionCall(GetInputEquil)
  data <- SubfunctionCall(GetDataEquil)
  kds <- log10(Norm(l)/Norm(rowSums(data)))
  kds <- kds - kds[length(kds)]
  kds <- kds/2
  names(kds) <- paste0(rownames(sXc), "_Kd")
  pars <- c(kds, bg=-4, AGO=1)
  random_pars <- rnorm(length(pars), 0, pars_sd)
  names(random_pars) <- names(pars)
  pars <- pars + random_pars
  pars["None_Kd"] <- 0
  return(pars)
}



tick <- 0
# pars <- InitializeEquilSitePars(sXc, combined=combined)*log(10)
OptimizeEquilSitePars <- function(sXc, pars=NULL, plotname_=plotname) {
  time_start <<- proc.time()[3]
  tick <- 0
  if (is.null(pars))
    initial.pars <- InitializeEquilSitePars(sXc)*log(10)
  else
    initial.pars <- pars*log(10)
  if (sitelist=="12mers" & experiment=="equilibrium2_nb") {
    L_ <- 300
  } else {
    L_ <- 100
  }
  data <- sXc[, 3:(ncol(sXc) - 1)]

  Y <- colSums(data)

  l <- Norm(sXc[, 1 + combined])*L_
  y <- as.double(as.numeric(as.matrix(sXc)))
  Y <- colSums(data)

  dil = sapply(colnames(sXc)[3:(ncol(sXc) - 1)], as.numeric)
  n_i <- nrow(sXc)
  n_j <- 3:(ncol(sXc) - 1)
  solution <- optim(initial.pars,
                    CostCNewInteractive,
                    gr = GradCNew,
                    sXc = sXc,
                    dil=dil,
                    l=l,
                    L=L_,
                    Y=Y,
                    n_i=n_i,
                    n_j=n_j,
                    plot_ = plot_,
                    plotname = plotname_,
                    tempname=tempname,
                    method = "L-BFGS-B",
                    lower=log(10)*c(rep(-6, length=length(initial.pars)-2), -4, -0.5),
                    upper=log(10)*c(rep(4, length=length(initial.pars)-2), 1, 1),
                    control = list(maxit=10000000, factr=1000, fnscale=1))
  output.pars <- solution$par
  print(proc.time()[3] - time_start)
  output.pars
}
pars.MLE <- OptimizeEquilSitePars(sXc)/log(10)
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
    sXc.resample <- apply(sXc, 2, function(col) {rmultinom(1, sum(col), col)})
    rownames(sXc.resample) <- rownames(sXc)
    if (sitelist != "12mers") {
      sXc.resample.withhold <- sXc.resample[, -(j + 3)]

    } else {
      sXc.resample.withhold <- sXc.resample
      plotname <- paste0(plot_dir, "/", j*reps + i, ".pdf")
    }
    resampled.pars[, i_full] <- OptimizeEquilSitePars(sXc.resample.withhold,
                                                      plotname_=plotname)
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




