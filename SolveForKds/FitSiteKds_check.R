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
# args       <- commandArgs(trailingOnly=TRUE)
# mirna      <- args[1]
# experiment <- args[2]
# n_constant <- args[3]
# sitelist   <- args[4]

mirna <- "miR-1"
experiment <- "equilibrium"
n_constant <- 5
sitelist <- "paper"


# if ("-reps" %in% args) {
#   reps <- as.integer(args[which(args == "-reps") + 1])
# } else {
  reps <- 20
# }
# if (sitelist == "12mers") {
#   mir_start = as.integer(args[which(args == "-mir_start") + 1])
#   str.mir_start = sprintf("_%s-%s", mir_start, mir_start + 3)
# } else {
  mir_start = FALSE
  str.mir_start = ""
# }

# if ("-nocombI" %in% args) {
#   combined <- FALSE
#   str.combined <- "_nocombInput"
# } else {
  combined <- FALSE
  str.combined <- ""
# }

# if ("-fixed" %in% args & mirna %in% c("miR-7-23nt", "miR-7-24nt", "miR-7-25nt")
#     & experiment == "equilibrium2_nb") {
#   fixed <- TRUE
#   str.fixed <- "_fixed"
# } else {
  fixed <- FALSE
  str.fixed <- ""
# }

L <- 100

# MAIN #########################################################################
# 1. Get data table:

mir_list <- unlist(strsplit(mirna, ","))
exp_list  <- unlist(strsplit(experiment, ","))
mir_exp_list <- data.frame(mirna=rep(mir_list, length(exp_list)),
                           experiment=rep(exp_list,each=length(mir_list)))

print("Mirna experiment list:")
print(mir_exp_list)
sXc <- apply(mir_exp_list, 1, function(row) {
  sXc_sites <- SitesXCounts(row[1], experiment=row[2], n_constant=n_constant,
               sitelist=sitelist, mirna.start=mirna.start)
  sXc_flanks <- SiteFlanksXCounts(row[1], "8mer", experiment=row[2], n_constant=5,
                    sitelist="paperfinal", buffer=TRUE)
  out <- rbind(sXc_flanks, sXc_sites[-1, ])
})

if (length(sXc) > 1) {
  str.global <- "_global"
} else {
  str.global <- ""
}
str_condition <- sprintf("%s_%s%s%s%s%s_PAPER", n_constant, sitelist,
                         str.mir_start, str.combined, str.global, str.fixed)
print(str_condition)
# kOutputFileFull <- GetAnalysisPath(mirna, experiment,
#                                condition=sprintf("%s_%s%s%s%s_PAPER_full", n_constant,
#                                                  sitelist, str.mir_start, 
#                                                  str.combined),
#                                analysis_type="kds_PAPER")
# kOutputFileMean <- GetAnalysisPath(mirna, experiment,
#                                condition=sprintf("%s_%s%s%s%s_PAPER_mean", n_constant,
#                                                  sitelist, str.mir_start, 
#                                                  str.combined),
#                                analysis_type="kds_PAPER")
# print(kOutputFile)

if (sitelist == "12mers") {
  pars_sd <- 0.01
  plot_ <- TRUE
  plotname <- paste0("12merfits/", paste(mirna, experiment, n_constant,
                                         sitelist, mir_start,
                     sep="_"), ".pdf")
} else {
  pars_sd <- 2
  plot_ <- FALSE
  plotname <- NULL
}
stretch <- 1.5

print(mirna)
print(experiment)
InitializeEquilSitePars <- function(sXc, combined=TRUE, fixed=FALSE) {
  n_mir = length(sXc)
  l <- SubfunctionCall(GetInputEquil, sXc=sXc[[1]])
  data <- SubfunctionCall(GetDataEquil, sXc=sXc[[1]])
  kds <- log(Norm(l)/Norm(rowSums(data)))
  kds <- kds - kds[length(kds)]
  kds <- kds/2
  names(kds) <- paste0(rownames(sXc[[1]]), "_Kd")
  
  if (fixed) {
    if (altfixed) {
      combined_fixed <- TRUE
    } else {
      combined_fixed <- FALSE
    }
    pars <- EquilPars(mirna, experiment=experiment, n_constant=n_constant,
                      sitelist=sitelist, uniq=uniq, combined=combined_fixed, global=TRUE)
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
  if (fixed) {
    random_pars[length(pars) - 1] = 0
    random_pars[length(pars)    ] = 0
  }
  names(random_pars) <- names(pars)
  pars <- pars + random_pars
  pars["None_Kd"] <- 0
  return(pars)
}

tick <- 0
OptimizeEquilSitePars <- function(sXc, pars=NULL, fixed=fixed, AGOfixed=FALSE,
                                  plotname_=plotname) {
  print(combined)
  print(fixed)
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
  dil <- dil * 1/(stretch^c(0, 1, 2, 3, 4))

  print(dil)
  n_i <- nrow(data)
  n_mir <- length(sXc)
  data_vec <- as.numeric(as.matrix(data))
  lower <- log(rep(c(1e-6, 0.001, 0.1), c(n_i, n_mir, n_mir)))
  upper <- log(rep(c(1e4, 10, 10), c(n_i, n_mir, n_mir)))

  if (fixed) {
    zero_grad <- n_i + c(0, 1, 2)
  } else if (AGOfixed) {
    zero_grad <- n_i + c(0, 1, 2)
  } else {
    zero_grad <- c(n_i)
  }
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
  time_elapsed <- time_start - proc.time()[3]
  # print(time_elapsed)
  output.pars/log(10)
}

ResampleData <- function(sXc) {
  sXc_new <- apply(sXc, 2, function(col) {rmultinom(1, sum(col), col)})
  rownames(sXc_new) <-rownames(sXc)
  return(sXc_new)
}


# sXc_minus_1 <- list(ResampleData(sXc[[1]][,c(-3,-4, -5)]))
# sXc_minus_2 <- list(ResampleData(sXc[[1]][,c(-4,-5, -6)]))
# sXc_minus_3 <- list(ResampleData(sXc[[1]][,c(-5,-6, -7)]))
# sXc_minus_4 <- list(ResampleData(sXc[[1]][,c(-6, -7, -3)]))
# sXc_minus_5 <- list(ResampleData(sXc[[1]][,c(-7, -3, -4)]))

# sXc_minus_1 <- list(ResampleData(sXc[[1]][,c(-3,-4)]))
# sXc_minus_2 <- list(ResampleData(sXc[[1]][,c(-4,-5)]))
# sXc_minus_3 <- list(ResampleData(sXc[[1]][,c(-5,-6)]))
# sXc_minus_4 <- list(ResampleData(sXc[[1]][,c(-6, -7)]))
# sXc_minus_5 <- list(ResampleData(sXc[[1]][,c(-7, -3)]))

# sXc_minus_1 <- list(ResampleData(sXc[[1]][,c(-3)]))
# sXc_minus_2 <- list(ResampleData(sXc[[1]][,c(-4)]))
# sXc_minus_3 <- list(ResampleData(sXc[[1]][,c(-5)]))
# sXc_minus_4 <- list(ResampleData(sXc[[1]][,c(-6)]))
# sXc_minus_5 <- list(ResampleData(sXc[[1]][,c(-7)]))


pars.MLE <- OptimizeEquilSitePars(sXc, fixed=fixed)
# pars.MLE_1 <- OptimizeEquilSitePars(sXc_minus_1, fixed=fixed)
# pars.MLE_2 <- OptimizeEquilSitePars(sXc_minus_2, fixed=fixed)
# pars.MLE_3 <- OptimizeEquilSitePars(sXc_minus_3, fixed=fixed)
# pars.MLE_4 <- OptimizeEquilSitePars(sXc_minus_4, fixed=fixed)
# pars.MLE_5 <- OptimizeEquilSitePars(sXc_minus_5, fixed=fixed)


print('got here')


PlotSiteEnrichment <- function(sXc, pars) {
  sXc_norm <- t(t(sXc)/colSums(sXc))

  R_plots <- sXc_norm/sXc_norm[, 2]
  col_inds <- 3:(ncol(R_plots) - 1)
  l = Norm(sXc_norm[,2])*100

  dils <- as.numeric(colnames(sXc)[col_inds])
  dils <- dils * 1/(stretch^c(0, 1, 2, 3, 4))
  print(dils)
  dils_model <- 10^seq(log10(dils[5]/3), log10(dils[1]*3), length.out=20)
  print(dils_model)
  dev.new(xpos=20, ypos=20, height=5, width=5)
  par(kPlotParameters)
  xmin <- min(dils_model)
  xmax <- max(dils_model)
  ymin <- 0.1
  ymax <- 300
  BlankPlot(log='xy')
  AddLogAxis(1, label="Concentration")
  AddLogAxis(2, label="R")

  A <- 10^pars["AGO_miR-1"]*dils_model/100
  kds <- 10^pars[1:nrow(sXc)]

  a_free <- sapply(A, function(A_i) {
    FreeAgo(kds=kds, l=l, A=A_i)
    })
  x_mat <- sapply(A, function(A_i) {
    a <- FreeAgo(kds=kds, l=l, A=A_i)
    x <- l*a/(a + kds)
    return(x)
    })
  f <- apply(x_mat, 2, function(x_col) {
    l - x_col
    })
  bg <- 10^pars["bg_miR-1"]
  x_tots <- x_mat + t(t(f)/colSums(f))*bg
  x_norm <- t(t(x_tots)/colSums(x_tots))
  R_model <- x_norm/Norm(l)
  cols <- c(GetColorFunction(flanks=kFlanks),
            kSiteColors[rownames(sXc)[257:nrow(sXc)]])
  xy <- GetPlotFractionalCoords(0.05, 0.95, log='xy')
  text(xy[1], xy[2], labels=stretch)

  sapply(1:nrow(sXc_norm), function(row_i) {
    points(dils,
           R_plots[row_i, col_inds], col=cols[row_i])
    lines(dils_model,
           R_model[row_i, ], col=cols[row_i])

    })
}

PlotSiteEnrichment(sXc[[1]], pars.MLE)

break
dev.set(2)
PlotSiteEnrichment(sXc_minus_1[[1]], pars.MLE_1)
dev.set(3)
PlotSiteEnrichment(sXc_minus_2[[1]], pars.MLE_2)
dev.set(4)
PlotSiteEnrichment(sXc_minus_3[[1]], pars.MLE_3)
dev.set(5)
PlotSiteEnrichment(sXc_minus_4[[1]], pars.MLE_4)
dev.set(6)
PlotSiteEnrichment(sXc_minus_5[[1]], pars.MLE_5)

break


break

# print(sprintf("%.16f", 10^pars.MLE))
# resampled.pars <- matrix(NaN, nrow=length(pars.MLE), ncol=5*reps,
#                          dimnames=list(names(pars.MLE), 1:(5*reps)))
# print("HERE")

for (j in 0:4) {
  # print(j)
  for (i in 1:reps) {
    i_full <- j*reps + i
    # print(i_full)
    tick <- 0
    sXc.resample <- lapply(sXc, function(sXc_i) {
      sXc_new <- apply(sXc_i, 2, function(col) {rmultinom(1, sum(col), col)})
      rownames(sXc_new) <-rownames(sXc_i)
      return(sXc_new)
    })
    print("original:")
    print(sXc)
    print("resample:")
    print(sXc.resample)
    if (sitelist != "12mers") {
      sXc.resample.withold <- lapply(sXc.resample, function(sXc_i) {
        leave_out <- sample(3:(ncol(sXc_i) - 1), 1)
        leave_out <- j + 3
        print(leave_out)
        sXc_i[, -leave_out]
      })
      print("resample and dropout:")
      print(sXc.resample.withold)
      print("i_full:")
      print(i_full)
      # sXc.resample.withhold <- sXc.resample
    } else {
      sXc.resample.withhold <- sXc.resample
      plotname <- paste0("12merfits/", paste(mirna, experiment, n_constant, sitelist, mir_start,
                     sep="_"), "_rep", j*reps + i, ".pdf")
    }
    resampled.pars[, i_full] <- OptimizeEquilSitePars(sXc.resample.withhold,
                                                      fixed=fixed,
                                                      plotname_=plotname)
    print(resampled.pars)
    resampled.pars.sort <- t(apply(resampled.pars[, 1:i_full, drop=FALSE], 1, sort))
    print(resampled.pars.sort)
    output <- 10^cbind(pars.MLE,
                     rowMeans(resampled.pars, na.rm=TRUE),
                     resampled.pars.sort[, ceiling(0.025*i_full)],
                     resampled.pars.sort[, ceiling(0.5*i_full)],
                     resampled.pars.sort[, ceiling(0.975*i_full)])
    colnames(output) <- c("Full","Mean", "Lower_CI", "Median", "Upper_CI")
    sapply(mir_list, function(mirna) {
      kOutputFile <- GetAnalysisPath(mirna, experiment, condition=str_condition,
                               analysis_type="kds_PAPER")
      # write.table(file=kOutputFile, output, sep="\t", quote = FALSE)

    })
    if (sitelist == "12mers") {
      kds.mean <- Logit(10^rowMeans(resampled.pars, na.rm=TRUE), 100)
      # print(kds.mean[1:5])
      # write.table(file=kOutputFileMean, t(kds.mean), sep="\t", quote=FALSE,
      #            row.names=FALSE)
  }
}
}

print(output)
warnings()





