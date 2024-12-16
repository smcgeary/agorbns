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
n_constant <- 5
if ("-reps" %in% args) {
  reps <- as.integer(args[which(args == "-reps") + 1])
} else {
  reps <- 1
}
  combined <- TRUE
  str.combined <- ""

sitelist <- "16mers"
# MAIN #########################################################################
# 1. Get data table:
mirnas <- c("miR-7-23nt", "miR-7-24nt", "miR-7-25nt")
experiment <- "equilibrium2_nb"
sXc <- lapply(1:5, function(mirna_start) {
  out <- lapply(c("left", "right"), function(split) {
    out <- lapply(mirnas, function(mirna) {
      SitesXCounts(mirna, experiment, n_constant, sitelist,
                   mirna.start=mirna_start, split16=split) + 1
    })
    names(out) <- mirnas
    out
  })
  names(out) <- c("left", "right")
  out
})
names(sXc) <- sprintf("nt%s-%s", 1:4, 4:7)

plot_ <- TRUE
plot_dir <- paste0("16merfits/", paste("miR-7", experiment, n_constant,
                                       str.combined, sep="_"))
tempname <- paste0(plot_dir, "/temp")
if (!file.exists(plot_dir)) {
  dir.create(plot_dir)
}


InitializeEquilSiteParsAll <- function(sXc, combined=FALSE) {
  if (combined) {
    input_i <- 2
  } else {
    input_i <- 1
  }
  kds <- c()
  total_exps <- 0
  sapply(sXc, function(x) {
    sapply(x, function(y) {
      sXc_i <- y[[1]]
      kds_i <- log(Norm(sXc_i[, input_i])/
                   Norm(rowSums(sXc_i[, 3:(ncol(sXc_i) - 1)]))) + 1
      kds_i <- kds_i - kds_i[length(kds_i)]
      names(kds_i) <- paste0(rownames(sXc_i), "_Kd")
      kds <<- c(kds, kds_i)
      total_exps <<- total_exps + 1
    })
  })
  c(kds, b23=log(1), b24=log(1), b25=log(1),
    AGO23=log(1), AGO24=log(1), AGO25=log(1))
}

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
  dils <- c()
  data_ls <- c()
  data_rs <- c()
  l_is <- c()
  sapply(sXc, function(x) {
    sapply(x, function(y) {
      l_is <<- c(l_is, i_counter + which(colnames(y[[1]])=="I_combined") - 1)
      sapply(y, function(sXc_i) {
        cond <- colnames(sXc_i)
        dils <<- c(dils, as.numeric(cond[!(cond %in% c("I", "I_combined", "0"))]))
        i_temp <- i_counter
        i_counter <<- i_counter + ncol(sXc_i)
        data_ls <<- c(data_ls, i_temp + which(colnames(sXc_i)=="I_combined"))
        data_rs <<- c(data_rs, i_temp + which(colnames(sXc_i)=="0") - 1)
      })
    })
  })
  sXc_vec <- unlist(sapply(sXc, function(x) {
    sapply(x, function(y) {
      sapply(y, function(z) as.double(as.numeric(as.matrix(z))))
    })
  }))
  n_exp <- 10
  n_x <- nrow(sXc[[1]][[1]][[1]])
  solution <- optim(initial.pars,
                    CostEquilGlobal16mersMir7_C,
                    gr = GradEquilGlobal16mersMir7_C,
                    sXc = sXc_vec,
                    dils=dils,
                    data_ls = data_ls,
                    data_rs = data_rs,
                    l_is=l_is,
                    n_x=n_x,
                    n_exp=n_exp,
                    plot_ = plot_,
                    plotname = plotname_,
                    tempname = tempname_,
                    L = L_,
                    method = "L-BFGS-B",
                    lower=log(10)*c(rep(-8, length=length(initial.pars)-6), rep(-1, 6)),
                    upper=log(10)*c(rep(4, length=length(initial.pars)-6), rep(1, 6)),
                    control = list(maxit=1e8, factr=1e3, fnscale=1))
  output.pars <- solution$par/log(10)
  print(proc.time()[3] - time_start)
  output.pars
}

pars.MLE <- OptimizeEquilSitePars(sXc)
write.table(file=kOutputFileFull, pars.MLE, sep="\t", quote=FALSE,
              row.names=FALSE)
