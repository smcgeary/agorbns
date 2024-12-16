################################################################################
#GenerateSiteKds_PAPER.py
################################################################################

source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")
print("Out of general")
graphics.off()
library(numDeriv)
# Initial parameters and constants.
args       <- commandArgs(trailingOnly=TRUE)
mirna      <- args[1]
experiment <- args[2]
n_constant <- args[3]


# mirna <- "miR-1"
# experiment <- "equilibrium"
# n_constant <- 5


# mirna <- "miR-7-23nt,miR-7-24nt,miR-7-25nt"
# experiment <- "equilibrium2_nb"
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


if ("-L" %in% args) {
  L <- as.integer(args[which(args == "-L") + 1])
  str.L <- sprintf("_L%s", L)
} else {
  L <- 100
  str.L <- ""
}

# MAIN #########################################################################
# 1. Get data table:
mir_list <- unlist(strsplit(mirna, ","))

sXc <- lapply(mir_list, function(mirna) {
  out <- lapply(1:5, function(mir_start) {
    sXc_i <- SitesXCounts(mirna, experiment, n_constant, "12mers",
                 mirna.start=mir_start)
    # rbind(head(sXc_i, 4), tail(sXc_i, 1))
    sXc_i
  })
  names(out) <- sprintf("nt%s-%s", 1:5, 4:8)
  out
})

if (length(sXc) > 1) {
  str.global <- "_global"
} else {
  str.global <- ""
}
kOutputFileMeanList <- sapply(seq(5), function(mir_start) {
  GetAnalysisPath(mir_list[[1]], experiment,
                  condition=sprintf("%s_12mers_%s-%s%s%s_PAPER_logmean",
                                    n_constant, mir_start, mir_start + 3,
                                    str.combined, str.global),
                                    analysis_type="kds_PAPER")
  })

print(kOutputFileMeanList)

print(head(sXc[[1]][[1]]))



pars_sd <- 0.01
plot_ <- TRUE
# plotname <- paste0("12merfits/", paste(mirna, experiment, n_constant,
#                                        "12mers", mirna.start,
#                    sep="_"), ".pdf")

InitializeEquilSitePars <- function(sXc, combined=TRUE) {
  n_mir = length(sXc)
  print(n_mir)
  kds <- unlist(lapply(sXc[[1]], function(sXc_i) {
    frac_I <- Norm(sXc_i[, 1 + combined] + 1)
    frac_data <- Norm(rowSums(sXc_i[, 3:(ncol(sXc_i) - 1)] + 1))
    kds <- log(frac_I/frac_data)
    kds <- kds - kds[length(kds)]
    kds
  }))
  names(kds) <- paste0(names(kds), "_Kd")
  bgs <- rep(log(0.1), n_mir)
  As <- rep(log(10), n_mir)
  print(bgs)

  names(bgs) <- sprintf("bg_%s", mir_list)
  names(As) <- sprintf("AGO_%s", mir_list)
  pars <- c(kds, bgs, As)
  return(pars)
}

tick <- 0
OptimizeEquilSitePars <- function(sXc, pars=NULL, combined=TRUE,
                                  plotname_=plotname) {
  print("in optimization")
  time_start <- proc.time()[3]
  tick <- 0
  if (is.null(pars) == TRUE) {
    initial.pars <- InitializeEquilSitePars(sXc, combined=combined)
  } else {
    initial.pars <- pars
  }
  L_ <- L
  n_j <- c(sapply(sXc, function(sXc_i) {
    sapply(sXc_i, function(sXc_ij) {
      ncol(sXc_ij) - 3
    })
  }))
  print(n_j)


  data_mat <- do.call(cbind, sapply(sXc, function(sXc_i) {
    lapply(sXc_i, function(sXc_ij) {
      print(dim(sXc_ij))
      as.matrix(sXc_ij[, 3:(ncol(sXc_ij) - 1)])
    })
  }))
  # l <- Norm(sXc[[1]][, 1 + combined])*L_
  l <- do.call(c, lapply(sXc[[1]], function(sXc_i) {
    out <- Norm(sXc_i[, 1 + combined] + 1)*L_
    names(out) <- rownames(sXc_i)
    out
  }))
  # print(l)
  dil <- as.numeric(colnames(data_mat))/100.0
  print(dil)
  n_i <- nrow(data_mat)
  print(n_i)
  n_mir <- length(sXc)
  print(n_mir)
  n_pos <- 5
  print(n_pos)

  Y <- colSums(data_mat)
  # print(Y)
  data_vec <- as.numeric(as.matrix(data_mat))
  # print(data_vec)
  lower <- log(rep(c(1e-6, 0.001, 0.1), c(n_i*5, n_mir, n_mir)))
  upper <- log(rep(c(1e4, 10, 10), c(n_i*5, n_mir, n_mir)))
  zero_grad <- seq(5)*n_i
  n_z <- length(zero_grad)
  print("here")
  print(n_z)
  print(dil)
  solution <- optim(initial.pars,
                    CostC12mers,
                    gr = GradC12mers,
                    data = data_vec,
                    dil  = dil,
                    l = l,
                    L = L_,
                    Y = Y,
                    n_i = n_i,
                    n_j = n_j,
                    n_mir = n_mir,
                    n_pos = n_pos,
                    zero_grad = zero_grad,
                    n_z = n_z,
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
  print(time_elapsed)
  output.pars
}
pars.MLE <- OptimizeEquilSitePars(sXc, combined=combined)
# print(sprintf("%.16f", 10^pars.MLE))
resampled.pars <- matrix(NaN, nrow=length(pars.MLE), ncol=5*reps,
                         dimnames=list(names(pars.MLE), 1:(5*reps)))
print("HERE")

for (j in 0:4) {
  for (i in 1:reps) {
    i_full <- j*reps + i
    tick <- 0
    # Resample the data:
    sXc.resample <- lapply(sXc, function(sXc_i) {
      out <- lapply(sXc_i, function(sXc_ij) {
        sXc_new <- apply(sXc_ij, 2, function(col) {rmultinom(1, sum(col), col)})
        rownames(sXc_new) <-rownames(sXc_ij)
        return(sXc_new)        
        })
    })
    # sXc <- lapply(mir_list, function(mirna) {
    #   out <- lapply(1:5, function(mir_start) {
    #     sXc_i <- SitesXCounts(mirna, experiment, n_constant, "12mers",
    #                  mirna.start=mir_start)
    #     rbind(head(sXc_i, 4), tail(sXc_i, 1))
    #     # sXc_i
    #   })
    #   names(out) <- sprintf("nt%s-%s", 1:5, 4:8)
    #   out
    # })




    sXc.resample.withhold <- sXc.resample
    #   plotname <- paste0("12merfits/", paste(mirna, experiment, n_constant, sitelist, mir_start,
    #                  sep="_"), "_rep", j*reps + i, ".pdf")
    # }
    resampled.pars[, i_full] <- OptimizeEquilSitePars(sXc.resample.withhold,
                                                      combined=combined,
                                                      plotname_=plotname)
    resampled.pars.sort <- t(apply(resampled.pars[, 1:i_full, drop=FALSE], 1,
                                   sort))
    output <- 10^cbind(pars.MLE,
                       rowMeans(resampled.pars, na.rm=TRUE),
                       resampled.pars.sort[, ceiling(0.025*i_full)],
                       resampled.pars.sort[, ceiling(0.5*i_full)],
                       resampled.pars.sort[, ceiling(0.975*i_full)])
    colnames(output) <- c("Full","Mean", "Lower_CI", "Median", "Upper_CI")
    # sapply(mir_list, function(mirna) {
    #   kOutputFile <- GetAnalysisPath(mirna, experiment, condition=str_condition,
    #                            analysis_type="kds_PAPER")
    #   print(kOutputFile)
    #   write.table(file=kOutputFile, output, sep="\t", quote = FALSE)
    # })
    kds.mean <- log(10^rowMeans(resampled.pars, na.rm=TRUE))
    n_i <- nrow(sXc[[1]][[1]])
    kds.mean_1.4 <- kds.mean[1:(n_i)]
    kds.mean_2.5 <- kds.mean[(n_i + 1):(2*n_i)]
    kds.mean_3.6 <- kds.mean[(2*n_i + 1):(3*n_i)]
    kds.mean_4.7 <- kds.mean[(3*n_i + 1):(4*n_i)]
    kds.mean_5.8 <- kds.mean[(4*n_i + 1):(5*n_i)]
    kds_means <- list(kds.mean_1.4, kds.mean_2.5, kds.mean_3.6,
                      kds.mean_4.7, kds.mean_5.8)
    for (i in seq(5)) {
      outputfile <- kOutputFileMeanList[i]
      print(outputfile)
      kds_out <- kds_means[[i]]
      write.table(file=outputfile, kds_out, sep="\t", quote=FALSE,
                  row.names=TRUE, col.names=FALSE)

    }
  }
}


print(output)
warnings()





