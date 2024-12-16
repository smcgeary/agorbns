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
if ("-reps" %in% args) {
  reps <- as.integer(args[which(args == "-reps") + 1])
} else {
  reps <- 1
}
if ("-nocombI" %in% args) {
  combined <- FALSE
  str.combined <- "nocombInput_"
} else {
  combined <- TRUE
  str.combined <- ""
}

sitelist <- "16mers"
# MAIN #########################################################################
# 1. Get data table:
sXc <- lapply(1:5, function(mirna_start) {
  out <- lapply(c("left", "right"), function(split) {
    SitesXCounts(mirna, experiment, n_constant, sitelist,
                 mirna.start=mirna_start,
                 split16=split) + 1
  })
  names(out) <- c("left", "right")
  out
})
names(sXc) <- sprintf("nt%s-%s", 1:4, 4:7)

kOutputFile <- GetAnalysisPath(mirna, experiment,
                               condition=sprintf("%s_16mers_global_%sPAPER", n_constant, 
                                                 str.combined),
                               analysis_type="kds_PAPER")
kOutputFileFull <- GetAnalysisPath(mirna, experiment,
                               condition=sprintf("%s_16mers_global_%sPAPER_full", n_constant, 
                                                 str.combined),
                               analysis_type="kds_PAPER")
kOutputFileMean <- GetAnalysisPath(mirna, experiment,
                               condition=sprintf("%s_16mers_global_%sPAPER_mean", n_constant, 
                                                 str.combined),
                               analysis_type="kds_PAPER")
print(kOutputFile)
print(kOutputFileFull)
print(kOutputFileMean)

print(sapply(sXc, function(x) {
  sapply(x, dim)
}))

plot_ <- TRUE
plot_dir <- paste0("16merfits/", paste(mirna, experiment, n_constant, str.combined,
                                       sep="_"))

tempname <- paste0(plot_dir, "/temp")

if (!file.exists(plot_dir)) {
  dir.create(plot_dir)
}



InitializeEquilSitePars <- function(sXc, combined=FALSE) {
  if (combined) {
    input_i <- 2
  } else {
    input_i <- 1
  }
  kds <- log(Norm(sXc[, input_i])/
             Norm(rowSums(sXc[, 3:(ncol(sXc)-1)])))
  kds <- kds - kds[length(kds)]
  names(kds) <- paste0(rownames(sXc), "_Kd")
  pars <- c(kds, bg=log(1), AGO=log(1))
  return(pars)
}

# # TESTING THE INDIVIDUAL COST FUNCTIONS #
# costs_individual <- sapply(sXc, function(x) {
#   sapply(x, function(sXc_i) {
#     tick <<- 0
#     pars <- InitializeEquilSitePars(sXc_i)
#     print(pars[1])
#     CostCNew(pars, sXc_i, dil=as.numeric(colnames(sXc_i)[3:7]),
#              n_x=nrow(sXc_i), input_i=1)
#   })
# })

# tick <<- 0
# time.past <<- proc.time()[3]
# sXc_1 <- sXc[[1]][[1]]
# pars_1 <- InitializeEquilSitePars(sXc_1)

# grad_analytical_1 <- GradCNew(pars_1, sXc_1,
#                               dil=as.numeric(colnames(sXc_1)[3:7]),
#                               n_x=nrow(sXc_1), input_i=1)


# grad_numerical_1 <- grad(CostCNew, pars_1, sXc=sXc_1,
#                          dil=as.numeric(colnames(sXc_1)[3:7]),
#                          n_x=nrow(sXc_1), input_i=1)

# plot(grad_analytical_1, grad_numerical_1)
# break
# print(costs_individual)
InitializeEquilSiteParsAll <- function(sXc, combined=FALSE) {
  if (combined) {
    input_i <- 2
  } else {
    input_i <- 1
  }
  kds <- c()
  total_exps <- 0
  sapply(sXc, function(x) {
    sapply(x, function(sXc_i) {
      kds_i <- log(Norm(sXc_i[, input_i])/
                   Norm(rowSums(sXc_i[, 3:(ncol(sXc_i) - 1)])))
      kds_i <- kds_i - kds_i[length(kds_i)]
      names(kds_i) <- paste0(rownames(sXc_i), "_Kd")
      print(length(kds_i))
      kds <<- c(kds, kds_i)
      total_exps <<- total_exps + 1
    })
  })
  c(kds, bg=log(1), AGO=log(1))
}

# pars <- InitializeEquilSiteParsAll(sXc)


i_counter <- 0

n_xs <- c()
dils <- c()
data_ls <- c()
data_rs <- c()
l_is <- c()
num_exps <- 0
sapply(sXc, function(x) {
  sapply(x, function(sXc_i) {
    n_xs <<- c(n_xs, nrow(sXc_i))
    cond <- colnames(sXc_i)
    dils <<- c(dils, as.numeric(cond[!(cond %in% c("I", "I_combined", "0"))]))
    i_temp <- i_counter
    i_counter <<- i_counter + ncol(sXc_i)
    data_ls <<- c(data_ls, i_temp + which(colnames(sXc_i)=="I_combined"))
    data_rs <<- c(data_rs, i_temp + which(colnames(sXc_i)=="0") - 1)
    l_is <<- c(l_is, i_temp + which(colnames(sXc_i)=="I") - 1)
    num_exps <<- num_exps + 1
  })
})

print(n_xs)
# print(dils)
# print(data_ls)
# print(data_rs)
# print(l_is)

# print(tail(pars, n=30))

# sXc_vec <- unlist(sapply(sXc, function(x) {
#     sapply(x, function(y) as.double(as.numeric(as.matrix(y))))
# }))




# grad_numerical <- grad(CostC_Global_16mers, pars, sXc=sXc_vec, dils=dils,
#                        n_x=n_xs[1], data_ls=data_ls, data_rs=data_rs, l_is=l_is,
#                        num_exps=num_exps, L=100, plot_=plot_,
#                                 plotname=plotname, tempname=tempname)

# grad_simple <- grad(CostC_Global_16mers, pars, sXc=sXc_vec, dils=dils,
#                        n_x=n_xs[1], data_ls=data_ls, data_rs=data_rs, l_is=l_is,
#                        num_exps=num_exps, L=100, plot_=plot_,
#                                 plotname=plotname, tempname=tempname, method="simple")

# for (i in 1:8) {
#   print((n_xs[1])*i)
#   grad_numerical[(n_xs[1])*i] <- 0
#   grad_simple[(n_xs[1])*i] <- 0
#   print(grad_numerical)
# }
# time0 <- proc.time()[3]
# for (i in 1:10) {
#   gradient_data <- GradC_Global_16mers(pars, sXc_vec, dils, n_xs[1], data_ls,
#                                        data_rs, l_is, num_exps, L=100,
#                                        plot_=plot_, plotname=plotname,
#                                        tempname=tempname)  
# }
# time1 <- proc.time()[3]
# print(time1 - time0)
# time2 <- proc.time()[3]
# for (i in 1:10) {
#   gradient_efficient <- GradC_Global_16mers_Parallel(pars, sXc_vec, dils, n_xs[1], data_ls,
#                                        data_rs, l_is, num_exps, L=100,
#                                        plot_=plot_, plotname=plotname,
#                                        tempname=tempname)
# }
# time3 <- proc.time()[3]
# print(time3 - time2)
# # par(mfrow=c(1, 3))
# plot(gradient_data, gradient_efficient)
# break
# abline(0, 1)
# plot(gradient_data, grad_simple)
# abline(0, 1)
# plot(grad_numerical, grad_simple)
# abline(0, 1)

# print(cbind(gradient_data, grad_simple, grad_numerical))
# break

# tick <- 0
# sXc <- sXc_l

#  print(CostCNew(pars_l, sXc_l, dil=as.numeric(colnames(sXc)[3:7]), n_x=nrow(sXc),
#                 input_i=1))
#  print(CostCNew(pars_r, sXc_r, dil=as.numeric(colnames(sXc)[3:7]), n_x=nrow(sXc),
#                 input_i=1))

#  break
# # pars <- InitializeEquilSitePars(sXc, combined=combined)*log(10)
OptimizeEquilSitePars <- function(sXc, pars=NULL, plotname_=plotname,
                                  tempname_=tempname) {
  time_start <- proc.time()[3]
  time.past <<- time_start
  tick <<- 0
  if (is.null(pars))
    initial.pars <- InitializeEquilSiteParsAll(sXc)
  else
    initial.pars <- pars
  if (experiment=="equilibrium2_nb") {
    L_ <- 300
  } else {
    L_ <- 100
  }
  if (combined) {
    input_i_ <- 2
  } else {
    input_i_ <- 1
  }
  i_counter <- 0
  n_xs <- c()
  dils <- c()
  data_ls <- c()
  data_rs <- c()
  l_is <- c()
  num_exps <- 0
  sapply(sXc, function(x) {
    sapply(x, function(sXc_i) {
      n_xs <<- c(n_xs, nrow(sXc_i))
      cond <- colnames(sXc_i)
      dils <<- c(dils, as.numeric(cond[!(cond %in% c("I", "I_combined", "0"))]))
      i_temp <- i_counter
      i_counter <<- i_counter + ncol(sXc_i)
      data_ls <<- c(data_ls, i_temp + which(colnames(sXc_i)=="I_combined"))
      data_rs <<- c(data_rs, i_temp + which(colnames(sXc_i)=="0") - 1)
      l_is <<- c(l_is, i_temp + which(colnames(sXc_i)=="I_combined") - 1)
      num_exps <<- num_exps + 1
    })
  })
  print("colsums external:")
  print(colSums(sXc[[1]][[1]]))
  sXc_vec <- unlist(sapply(sXc, function(x) {
      sapply(x, function(y) as.double(as.numeric(as.matrix(y))))
  }))
  time_0 <- proc.time()[3]
  solution <- optim(initial.pars,
                    CostC_Global_16mers,
                    gr = GradC_Global_16mers,
                    sXc = sXc_vec,
                    dils=dils,
                    n_x=n_xs[1],
                    num_exps=num_exps,
                    data_ls = data_ls,
                    data_rs = data_rs,
                    l_is=l_is,
                    plot_ = plot_,
                    plotname = plotname_,
                    tempname = tempname_,
                    L = L_,
                    method = "L-BFGS-B",
                    lower=log(10)*c(rep(-8, length=length(initial.pars)-2), -2, -0.5),
                    upper=log(10)*c(rep(4, length=length(initial.pars)-2), 1, 1),
                    control = list(maxit=1e8, factr=1e3, fnscale=1))
  # time_1 <- proc.time()[3]
  # print(time_1 - time_0)
  # solution2 <- optim(initial.pars,
  #                   CostC_Global_16mers,
  #                   gr = GradC_Global_16mers_Parallel,
  #                   sXc = sXc_vec,
  #                   dils=dils,
  #                   n_x=n_xs[1],
  #                   num_exps=num_exps,
  #                   data_ls = data_ls,
  #                   data_rs = data_rs,
  #                   l_is=l_is,
  #                   plot_ = plot_,
  #                   plotname = plotname_,
  #                   tempname = tempname_,
  #                   L = L_,
  #                   method = "L-BFGS-B",
  #                   lower=log(10)*c(rep(-8, length=length(initial.pars)-2), 0, -0.5),
  #                   upper=log(10)*c(rep(4, length=length(initial.pars)-2), 1, 1),
  #                   control = list(maxit=10, factr=1e3, fnscale=1))
  # time_2 <- proc.time()[3]
  # print(time_2 - time_1)


  output.pars <- solution$par/log(10)
  print(proc.time()[3] - time_start)
  output.pars
}

pars.MLE <- OptimizeEquilSitePars(sXc)
write.table(file=kOutputFileFull, pars.MLE, sep="\t", quote=FALSE,
              row.names=FALSE)

break
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





