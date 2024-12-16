################################################################################
#GenerateSiteKds_PAPER.py
################################################################################

# THIS IS THE OFFICIAL SCRIPT FOR FINALIZING THE MANUSCRIPT.
# THIS SCRIPT HAS THESE DECISIONS MADE:
# 1.) PROTEIN IS FIT AS A PARAMETER
# 2.) A SINGLE BACKGROUND PARAMETER IS GIVEN TO ALL FIVE CONCENTRATIONS.
# 3.) LIST, START, AND STOP POSITION ARE MANDATORY PARAMETERS.
## 4.) ALL INPUTS ARE COMBINED ACROSS FOR THESE ANALYSES.
library(colorspace)
library(multicore)
library(data.table)
library(numDeriv)
graphics.off()
source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")
# Initial parameters and constants.
args       <- commandArgs(trailingOnly=TRUE)
mirna      <- args[1]
experiment <- args[2]
n_constant <- args[3]
sitelist   <- args[4]


if ("-reps" %in% args) {
  reps <- as.integer(args[which(args == "-reps") + 1])
} else {
  reps <- 20
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
if ("-fixed" %in% args) {
  fixed <- TRUE
  str.fixed <- "_fixed"
} else {
  fixed <- FALSE
  str.fixed <- ""
}

# mirna <- "miR-155"
# experiment <- "equilibrium"
# n_constant <- 5
# sitelist <- "mismatch_and_threeprime"
# reps <- 50
# combined <- FALSE
# str.combined <- ""
# fixed <- FALSE
# str.fixed <- ""
# Loads general functions used in AgoRBNS analysis.

# TODO make a new table of miRNA concentrations or just convert within the
# script.
stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
                       "/SolveForKds/k_c_stockago.txt"), row.names=1,
                       header=TRUE, sep="\t")
k.c.stockago <- stockago[mirna,experiment]
k.c.lib <- 100

# MAIN #########################################################################
# 1. Get data table:

mir_list <- unlist(strsplit(mirna, ","))
exp_list  <- unlist(strsplit(experiment, ","))
mir_exp_list <- data.frame(mirna=rep(mir_list, length(exp_list)),
                           experiment=rep(exp_list,each=length(mir_list)))

str_condition <- sprintf("%s_%s%s_singleonly_PAPER", n_constant, sitelist,
                         str.combined)
print(str_condition)


print("Mirna experiment list:")
print(mir_exp_list)
sXc <- apply(mir_exp_list, 1, function(row) {
  sXc_1s <- SitesXCounts(row[1], experiment=row[2], n_constant=n_constant,
                         sitelist=sitelist)
  msXc <- SitesXCounts(row[1], experiment=row[2], n_constant=n_constant,
                       sitelist=sitelist, multisite=TRUE)



  regex_split         <- ",\\(\\d*\\),"
  regex_split_capture <- ",\\((\\d*)\\),"
  regex_site <- "[^,]*"
  regex_site_capture <- "([^,]*)"

  single_sites = grep(regex_split, rownames(msXc), invert=TRUE)
  double_sites = grep(sprintf("^%s%s%s$", regex_site, regex_split,
                              regex_site), rownames(msXc), perl=TRUE)
  nonoverlapping_sites = grep("\\|", rownames(msXc), invert=TRUE, perl=TRUE)
  single_nonoverlapping_sites = intersect(single_sites, nonoverlapping_sites)
  # double_nonoverlapping_sites = intersect(double_sites, nonoverlapping_sites)
  sXc_out <- msXc[single_nonoverlapping_sites, ][rownames(sXc_1s), ]
  sXc_out
})

pars_sd <- 0.1
plot_ <- FALSE
plotname <- NULL
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

InitializeEquilSitePars(sXc, combined=combined)


tick <- 0
OptimizeEquilSitePars <- function(sXc, pars=NULL, fixed=fixed,
                                  plotname_=plotname) {
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
  print(proc.time()[3] - time_start)
  # print(tick)
  output.pars/log(10)
}
pars.MLE <- OptimizeEquilSitePars(sXc, fixed=fixed)
# print(sprintf("%.16f", 10^pars.MLE))
resampled.pars <- matrix(NaN, nrow=length(pars.MLE), ncol=5*reps,
                         dimnames=list(names(pars.MLE), 1:(5*reps)))



# dev.new(xpos=20, ypos=20, height=4, width=4)
# par(kPlotParameters)
# kds_paper <- EquilPars(mirna, sitelist="paper")
# sites_paper <- rownames(kds_paper)
# BlankPlot(log="xy", inv="x")
# AddLogAxis(1, label="Relative Kd; single site fit")
# AddLogAxis(2, label="single-site/paper fit; relative Kd")
# kds_paper_del <- kds_paper$Full/c((10^pars.MLE[sites_paper]))
# points(10^pars.MLE[sites_paper], kds_paper_del,
#        col=kSiteColors[c(gsub("_Kd", "", sites_paper), "bg", "AGO")],)



for (j in 0:4) {
  for (i in 1:reps) {
    i_full <- j*reps + i
    print(i_full)
    tick <- 0
    sXc.resample <- lapply(sXc, function(sXc_i) {
      sXc_new <- apply(sXc_i, 2, function(col) {rmultinom(1, sum(col), col)})
      rownames(sXc_new) <-rownames(sXc_i)
      return(sXc_new)
    })
    if (sitelist != "12mers") {
      sXc.resample.withhold <- lapply(sXc.resample, function(sXc_i) {
        leave_out <- sample(3:(ncol(sXc_i) - 1), 1)
        leave_out <- j + 3
        sXc_i[, -leave_out]
      })
    } else {
      sXc.resample.withhold <- sXc.resample
      plotname <- paste0("12merfits/", paste(mirna, experiment, n_constant, sitelist, mir_start,
                     sep="_"), "_rep", j*reps + i, ".pdf")
    }
    resampled.pars[, i_full] <- OptimizeEquilSitePars(sXc.resample.withhold,
                                                      fixed=fixed,
                                                      plotname_=plotname)
    resampled.pars.sort <- t(apply(resampled.pars[, 1:i_full, drop=FALSE], 1, sort))
    output <- 10^cbind(pars.MLE,
                     rowMeans(resampled.pars, na.rm=TRUE),
                     resampled.pars.sort[, ceiling(0.025*i_full)],
                     resampled.pars.sort[, ceiling(0.5*i_full)],
                     resampled.pars.sort[, ceiling(0.975*i_full)])
    colnames(output) <- c("Full","Mean", "Lower_CI", "Median", "Upper_CI")
    sapply(mir_list, function(mirna) {
      kOutputFile <- GetAnalysisPath(mirna, experiment, condition=str_condition,
                               analysis_type="kds_PAPER")
      write.table(file=kOutputFile, output, sep="\t", quote = FALSE)
    })
  }
}

print(output)
warnings()

