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
if (sitelist == "12mers") {
  mir_start = as.integer(args[which(args == "-mir_start") + 1])
  str.mir_start = sprintf("_%s-%s", mir_start, mir_start + 3)
} else {
  mir_start = FALSE
  str.mir_start = ""
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


if ("-fixed" %in% args & mirna %in% c("miR-7-23nt", "miR-7-24nt", "miR-7-25nt")
    & experiment == "equilibrium2_nb") {
  fixed <- TRUE
  str.fixed <- "_fixed"
} else {
  fixed <- FALSE
  str.fixed <- ""
}



# MAIN #########################################################################
# 1. Get data table:

mir_list <- unlist(strsplit(mirna, ","))
exp_list  <- unlist(strsplit(experiment, ","))
mir_exp_list <- data.frame(mirna=rep(mir_list, length(exp_list)),
                           experiment=rep(exp_list,each=length(mir_list)))

print("Mirna experiment list:")
print(mir_exp_list)
sXc <- apply(mir_exp_list, 1, function(row) {
  SitesXCounts(row[1], experiment=row[2], n_constant=n_constant,
               sitelist=sitelist, buffer=buffer, mirna.start=mirna.start)
})

if (length(sXc) > 1) {
  str.global <- "_global"
} else {
  str.global <- ""
}
str_condition <- sprintf("%s_%s%s%s%s%s%s_Xval_PAPER", n_constant, sitelist,
                         str.mir_start, str.buffer, str.combined, str.global, str.fixed)
print(str_condition)

if (sitelist == "12mers") {
  pars_sd <- 0.01
  plot_ <- TRUE
  plotname <- paste0("12merfits/", paste(mirna, experiment, n_constant,
                                         sitelist, mir_start,
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
  
  if (fixed == TRUE) {
    pars <- EquilPars(mirna, experiment=experiment, n_constant=n_constant,
                      sitelist=sitelist, global=TRUE)
    bgs <- log(pars[sprintf("bg_%s", mirna),]$Mean)
    As <- log(pars[sprintf("AGO_%s", mirna),]$Mean)
  } else {
    bgs <- rep(log(0.1), n_mir)
    As <- rep(log(10), n_mir)    
  }
  names(bgs) <- sprintf("bg_%s", mir_list)
  names(As) <- sprintf("AGO_%s", mir_list)
  pars <- c(kds, bgs, As)
  random_pars <- rnorm(length(pars), 0, pars_sd)
  if (fixed == TRUE) {
    random_pars[length(pars) - 1] = 0
    random_pars[length(pars)    ] = 0
  }
  names(random_pars) <- names(pars)
  pars <- pars + random_pars
  pars["None_Kd"] <- 0
  return(pars)
}


tick <- 0
OptimizeEquilSitePars <- function(sXc, pars=NULL, fixed=fixed,
                                  plotname_=plotname, combined=FALSE) {
  time_start <- proc.time()[3]
  tick <- 0
  if (is.null(pars) == TRUE) {
    initial.pars <- InitializeEquilSitePars(sXc, combined=combined, fixed=fixed)
  } else {
    initial.pars <- pars
  }
  if (sitelist=="12mers" & experiment=="equilibrium2_nb") {
    L_ <- 100
  } else {
    L_ <- 100
  }
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
  zero_grad <- c(n_i)
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
  # print(proc.time()[3] - time_start)
  # print(tick)
  output.pars/log(10)
}
pars.MLE <- OptimizeEquilSitePars(sXc, fixed=fixed, combined=combined)
print(sprintf("%.16f", 10^pars.MLE))
resampled.pars <- matrix(NaN, nrow=length(pars.MLE), ncol=5,
                         dimnames=list(names(pars.MLE), 1:5))
print("HERE")

for (j in 0:4) {
    i_full <- j + 1
    tick <- 0
    sXc.withhold <- lapply(sXc, function(sXc_i) {
      leave_out <- j + 3
      sXc_i[, -leave_out]
    })
    resampled.pars[, i_full] <- OptimizeEquilSitePars(sXc.withhold,
                                                      fixed=fixed,
                                                      plotname_=plotname)
    # resampled.pars.sort <- t(apply(resampled.pars[, 1:i_full, drop=FALSE], 1, sort))
    output <- 10^cbind(pars.MLE, resampled.pars)
    colnames(output) <- c("Full","-40", "-12.65", "-4", "-1.265", "-0.4")
    sapply(mir_list, function(mirna) {
      kOutputFile <- GetAnalysisPath(mirna, experiment, condition=str_condition,
                               analysis_type="kds_PAPER")
      print(kOutputFile)
      write.table(file=kOutputFile, output, sep="\t", quote = FALSE)

    })
    if (sitelist == "12mers") {
      kds.mean <- Logit(10^rowMeans(resampled.pars, na.rm=TRUE), 100)
      # print(kds.mean[1:5])
      # write.table(file=kOutputFileMean, t(kds.mean), sep="\t", quote=FALSE,
      #            row.names=FALSE)
  }
}
print(output)
warnings()





