################################################################################
#GenerateSiteKds_PAPER.py
################################################################################

source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")
print("Out of general")

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

if (sitelist == "12mers") {
  mirna.start = as.integer(args[which(args == "-mir_start") + 1])
  str.mir_start = sprintf("_%s-%s", mirna.start, mirna.start + 3)
} else {
  mirna.start = FALSE
  str.mir_start = ""
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
  str.fixed <- "_fixed"
} else {
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
               sitelist=sitelist, mirna.start=mirna.start, buffer=buffer)
})

if (length(sXc) > 1) {
  str.global <- "_global"
} else {
  str.global <- ""
}
str_condition <- sprintf("%s_%s%s%s%s%s%s_6merm8cor_PAPER", n_constant, sitelist,
                         str.buffer, str.mir_start, str.combined, str.global,
                         str.fixed, str.L)
print(str_condition)
# kOutputFileFull <- GetAnalysisPath(mirna, experiment,
#                                condition=sprintf("%s_%s%s%s%s_PAPER_full", n_constant,
#                                                  sitelist, str.mir_start, 
#                                                  str.combined),
#                                analysis_type="kds_PAPER")
kOutputFileMean <- GetAnalysisPath(mir_list[[1]], experiment,
                                   condition=sprintf("%s_%s%s%s%s_6merm8cor_PAPER_logmean",
                                                     n_constant, sitelist,
                                                     str.mir_start,
                                                     str.combined, str.global),
                                   analysis_type="kds_PAPER")
# print(kOutputFileMean)
# print("here")
# break

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
  
  if (fixed == TRUE) {
    pars <- EquilPars(mirna, experiment=experiment, n_constant=n_constant,
                      sitelist=sitelist, combined=combined, global=TRUE)
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
  zero_grad <- c(n_i)
  n_z <- length(zero_grad)
  i_6m8 <- which(rownames(sXc[[1]]) == "6mer-m8")
  l_adj <- l
  l_adj[i_6m8] <- l_adj[i_6m8] + L_
  print(i_6m8)

  solution <- optim(initial.pars,
                    CostC_8mer,
                    gr = GradC_8mer,
                    data = data_vec,
                    dil  = dil,
                    l = l,
                    l_adj = l,
                    L = L_,
                    Y = Y,
                    n_i = n_i,
                    n_j = n_j,
                    n_mir = n_mir,
                    zero_grad = zero_grad,
                    n_z = n_z,
                    i_6m8 = i_6m8,
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
  print(time_elapsed)
  output.pars/log(10)
}
pars.MLE <- OptimizeEquilSitePars(sXc, fixed=fixed)
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
      sXc_new <- apply(sXc_i, 2, function(col) {rmultinom(1, sum(col), col)})
      rownames(sXc_new) <-rownames(sXc_i)
      return(sXc_new)
    })
    # If not 12mers, remove a column:
    if (sitelist != "12mers") {
      if (mirna == "miR-7-23nt" & experiment == "equilibrium2_nb") {
        min_col <- 4
      } else {
        min_col <- 3
      }
      sXc.resample.withhold <- lapply(sXc.resample, function(sXc_i) {
        leave_out <- sample(min_col:(ncol(sXc_i) - 1), 1)
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

print(output)
warnings()





