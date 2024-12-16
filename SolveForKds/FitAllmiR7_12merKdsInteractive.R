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
mirna      <- "miR-7-23nt"
experiment <- "equilibrium2_nb"
n_constant <- 5
sitelist   <- "12mers"
reps <- 20
str.global = "_globalfit"



combined <- FALSE
if (combined) {
  str.combined <- ""
} else {
  str.combined <- "_nocombInput"
}

mir_start <- 2
# MAIN #########################################################################
# 1. Get data table:
mirnas <- c("miR-7-23nt", "miR-7-24nt", "miR-7-25nt")
sXc <- lapply(mirnas, function(mirna) {
                    out <- lapply(1:5, function(mir_start) {
                              SitesXCounts(mirna, experiment, n_constant,
                                  sitelist, mirna.start=mir_start)})
                    names(out) <- sprintf("nt%s-%s", 1:5, 4:8)
                    out
                  })
names(sXc) <- mirnas
kOutputFile <- GetAnalysisPath(mirna, experiment,
                               condition=sprintf("%s_%s%s%s_allcombined_PAPER", n_constant,
                                                 sitelist, str.global, 
                                                 str.combined),
                               analysis_type="kds_PAPER")
kOutputFileFull <- GetAnalysisPath(mirna, experiment,
                               condition=sprintf("%s_%s%s%s_allcombined_PAPER_full", n_constant,
                                                 sitelist, str.global, 
                                                 str.combined),
                               analysis_type="kds_PAPER")
kOutputFileMean <- GetAnalysisPath(mirna, experiment,
                               condition=sprintf("%s_%s%s%s_allcombined_PAPER_mean", n_constant,
                                                 sitelist, str.global, 
                                                 str.combined),
                               analysis_type="kds_PAPER")
print(kOutputFile)
print(kOutputFileFull)
print(kOutputFileMean)
plot_ <- TRUE
plot_dir <- paste0("12merfits/", paste(mirna, experiment, n_constant,
                                       sitelist, mir_start, str.combined,
                                       sep="_"))
plotname <- paste0(plot_dir, "/full.pdf")

tempname <- paste0(plot_dir, "global")
if (!file.exists(plot_dir)) {
  dir.create(plot_dir)
}

InitializeEquilSitePars <- function(sXc, combined=TRUE) {
  if (combined) {
    input_i <- 2
  } else {
    input_i <- 1
  }
  kds1_4 <- log(Norm(sXc[["miR-7-23nt"]][["nt1-4"]][, input_i] + 1)/
             Norm(rowSums(sXc[["miR-7-23nt"]][["nt1-4"]][,3:(ncol(sXc[["miR-7-23nt"]][["nt1-4"]]))]) + 1))
  kds2_5 <- log(Norm(sXc[["miR-7-23nt"]][["nt2-5"]][, input_i] + 1)/
             Norm(rowSums(sXc[["miR-7-23nt"]][["nt2-5"]][,3:(ncol(sXc[["miR-7-23nt"]][["nt2-5"]]))]) + 1))
  kds3_6 <- log(Norm(sXc[["miR-7-23nt"]][["nt3-6"]][, input_i] + 1)/
             Norm(rowSums(sXc[["miR-7-23nt"]][["nt3-6"]][,3:(ncol(sXc[["miR-7-23nt"]][["nt3-6"]]))]) + 1))
  kds4_7 <- log(Norm(sXc[["miR-7-23nt"]][["nt4-7"]][, input_i] + 1)/
             Norm(rowSums(sXc[["miR-7-23nt"]][["nt4-7"]][,3:(ncol(sXc[["miR-7-23nt"]][["nt4-7"]]))]) + 1))
  kds5_8 <- log(Norm(sXc[["miR-7-23nt"]][["nt5-8"]][, input_i] + 1)/
             Norm(rowSums(sXc[["miR-7-23nt"]][["nt5-8"]][,3:(ncol(sXc[["miR-7-23nt"]][["nt5-8"]]))]) + 1))
  kds1_4 <- kds1_4 - kds1_4[length(kds1_4)]
  kds2_5 <- kds2_5 - kds2_5[length(kds2_5)]
  kds3_6 <- kds3_6 - kds3_6[length(kds3_6)]
  kds4_7 <- kds4_7 - kds4_7[length(kds4_7)]
  kds5_8 <- kds5_8 - kds5_8[length(kds5_8)]

  kds <- c(kds1_4, kds2_5, kds3_6, kds4_7, kds5_8)
  kds <- kds/2
  names(kds) <- c(rownames(sXc[["miR-7-23nt"]][["nt1-4"]]),
                  rownames(sXc[["miR-7-23nt"]][["nt2-5"]]),
                  rownames(sXc[["miR-7-23nt"]][["nt3-6"]]),
                  rownames(sXc[["miR-7-23nt"]][["nt4-7"]]),
                  rownames(sXc[["miR-7-23nt"]][["nt5-8"]]))
  names(kds) <- paste0(names(kds), "_Kd")
  pars <- c(kds, b23=log(1), b24=log(1), b25=log(1), AGO23=log(4), AGO24=log(4), AGO25=log(4))
  return(pars)
}

tick <- 0
pars <- InitializeEquilSitePars(sXc, combined=combined)

OptimizeEquilSitePars <- function(sXc, pars=NULL, plotname_=plotname,
                                  tempname_=tempname) {
  time_start <- proc.time()[3]
  tick <- 0
  if (is.null(pars))
    initial.pars <- InitializeEquilSitePars(sXc)
  else
    initial.pars <- pars*log(10)
  print(length(initial.pars))
  L_ <- 100
  dils_all <- c(rep(as.numeric(colnames(sXc[["miR-7-23nt"]][["nt2-5"]])[3:6]), 5),
                rep(as.numeric(colnames(sXc[["miR-7-24nt"]][["nt2-5"]])[3:7]), 5),
                rep(as.numeric(colnames(sXc[["miR-7-25nt"]][["nt2-5"]])[3:7]), 5))
  dils_all <<- dils_all
  print(length(dils_all))
  n_x <- nrow(sXc[["miR-7-23nt"]][["nt1-4"]])
  print(n_x)
  print(n_x*5)
  L <- L_
  sXc_vec <- unlist(sapply(sXc, function(x) {
    sapply(x, function(y) as.double(as.numeric(as.matrix(y))))
  }))
  sXc_vec_1 <<- sXc_vec
  print(length(sXc_vec))
  print(n_x*length(dils_all))
  initial.pars <<- initial.pars
  lower <- log(10)*c(rep(-6, length=length(initial.pars)-6), rep(-1, 3), rep(-0.5, 3))
  upper <- log(10)*c(rep(2, length=length(initial.pars)-6), rep(1, 6))

  lower <<- lower
  upper <<- upper

  n_x <<- n_x
  solution <- optim(initial.pars,
                    CostCMir7All,
                    gr = GradCMir7All,
                    sXc = sXc_vec,
                    dils=dils_all,
                    n_x=n_x,
                    plot_ = plot_,
                    plotname = plotname_,
                    tempname=tempname_,
                    L = L_,
                    method = "L-BFGS-B",
                    lower=lower,
                    upper=upper,
                    control = list(maxit=1000000, factr=1000, fnscale=1))
  output.pars <- solution$par/log(10)
  print(proc.time()[3] - time_start)
  output.pars
}
pars.MLE <- OptimizeEquilSitePars(sXc)
check_data <- unlist(sXc[[1]][[1]][, 3:(ncol(sXc[[1]][[1]]) - 1)])
data_internal <- data_global[1:(n_x*4)]
# plot(data_internal, check_data)
break
kds.full <- pars.MLE
print(kds.full[1:5])
names(kds.full) <- names(pars.MLE)
write.table(file=kOutputFileFull, t(kds.full), sep="\t", quote=FALSE,
            row.names=FALSE)

# break
# resampled.pars <- matrix(NaN, nrow=length(pars.MLE), ncol=5*reps,
#                          dimnames=list(names(pars.MLE), 1:(5*reps)))

# for (j in 0:(ncol(sXc) - 4)) {
#   print(j)
#   for (i in 1:reps) {
#     i_full <- j*reps + i
#     # print(i_full)
#     tick <- 0
#     sXc.resample <- apply(sXc, 2, function(col) {rmultinom(1, sum(col), col)})
#     rownames(sXc.resample) <- rownames(sXc)
#     if (sitelist != "12mers") {
#       sXc.resample.withhold <- sXc.resample[, -(j + 3)]

#     } else {
#       sXc.resample.withhold <- sXc.resample
#       plotname <- paste0(plot_dir, "/", j*reps + i, ".pdf")
#     }
#     resampled.pars[, i_full] <- OptimizeEquilSitePars(sXc.resample.withhold,
#                                                       plotname_=plotname)
#     resampled.pars.sort <- t(apply(resampled.pars[, 1:i_full, drop=FALSE], 1, sort))
#     output <- 10^cbind(pars.MLE,
#                      rowMeans(resampled.pars, na.rm=TRUE),
#                      resampled.pars.sort[, ceiling(0.025*i_full)],
#                      resampled.pars.sort[, ceiling(0.5*i_full)],
#                      resampled.pars.sort[, ceiling(0.975*i_full)])
#     colnames(output) <- c("Full","Mean", "Lower_CI", "Median", "Upper_CI")
#     write.table(file=kOutputFile, output, sep="\t", quote = FALSE)
#     if (sitelist == "12mers") {
#       kds.mean <- rowMeans(resampled.pars, na.rm=TRUE)
#       # print(kds.mean[1:5])
#       write.table(file=kOutputFileMean, t(kds.mean), sep="\t", quote=FALSE,
#                  row.names=FALSE)
#   }
# }
# }
# print(output)
# print(kOutputFile)
# warnings()




