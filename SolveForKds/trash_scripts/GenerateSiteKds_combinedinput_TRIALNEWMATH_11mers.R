################################################################################
#GenerateSiteKds_combinedinput_PAPER.py
################################################################################

## THIS IS THE OFFICIAL SCRIPT FOR FINALIZING THE MANUSCRIPT.
## THIS SCRIPT HAS THESE DECISIONS MADE:
## 1.) PROTEIN IS FIT AS A PARAMETER
## 2.) A SINGLE BACKGROUND PARAMETER IS GIVEN TO ALL FIVE CONCENTRATIONS.
## 3.) SITELIST, START, AND STOP POSITION ARE MANDATORY PARAMETERS.
# ## 4.) ALL INPUTS ARE COMBINED ACROSS FOR THESE ANALYSES.
library(colorspace)
library(multicore)
library(data.table)
# # Initial parameters and constants.
# args  <- commandArgs(trailingOnly=TRUE)
# mirna <- args[1]
# mirna <- "let-7a"
# # Region within random sequence from which site types orginiates,
# # going from position [26 - "start" : 26 + 37 + "stop"]
# start <- as.integer(args[2])
# stop  <- as.integer(args[3])

start <- 5
stop <- 5


# # Parameter specifying which list of sites is used for the analysis.
# I.E "Current" includes all the current site types, "canonical" is just
# the six seed-based 6mers, 7mers, and 8mers, along with no site., etc.
# sitelist <- args[4]
sitelist <- "11mers"
# if (sitelist == "10mers") {
#   mirna.start <- as.integer(args[5])
#   mirna.stop <- as.integer(args[6])
# } else {
#   mirna.start <- NULL
#   mirna.stop <- NULL
# }

mirna.start <- 2
mirna.stop <- 5
# # Experiment name
experiment <- "equilibrium"

# Loads general functions used in AgoRBNS analysis.
source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")

# Loads the colors associated with each site type, for plotting purposes.

# Loads the table of Ago–miRNA concentrations, for the purposes of modeling
# them into the structure.
# NOTE These are actually higher than the real concentrations, because before
# I began the structural analysis I had to use 1.5-2X the amount of AGO.
# TODO make a new table of miRNA concentrations or just convert within the
# script.
stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
                       "/SolveForKds/k_c_stockago.txt"), row.names=1,
                       header=TRUE, sep="\t")
print(stockago)
k.c.stockago <- stockago[mirna,experiment]
print(k.c.stockago)
k.c.lib = 100


# MAIN #########################################################################
# 1. Get data table:
sitesXcounts <- GetSitesXCounts(mirna,
                                        experiment,
                                        start,
                                        stop,
                                        sitelist,
                                        mirna.start = mirna.start,
                                        mirna.stop = mirna.stop)
# Separate site sequences from data file.
if (sitelist == "11mers") {
  seqs <- rownames(sitesXcounts)
} else {
  seqs <- sitesXcounts[,1]
  sitesXcounts <- sitesXcounts[,-1]
}
# if (sitelist == "10mers") {
#   sitesXcounts_new <- t(sapply(1:((nrow(sitesXcounts)-1)/256),function(x){colSums(sitesXcounts[1:256+(x-1)*256,])}))
#   sitesXcounts_new <- rbind(sitesXcounts_new,sitesXcounts[nrow(sitesXcounts),,drop=FALSE])
#   sitesXcounts <- sitesXcounts_new
# }

# Remove any rows with no reads whoatsoever.
sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
# Remove any rows for which there are no input reads:
sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]

data <- sitesXcounts[,2:6]
data <- data.matrix(data)

# data_new <- t(sapply(1:((nrow(data)-1)/4),function(x){colSums(data[1:4+(x-1)*4,])}))
# rownames(data_new) <- rownames(data)[seq(1,nrow(data)-1,by=4)]
rownames(data)[nrow(data)] <- "None"
# data <- rbind(data_new,data[nrow(data),,drop=FALSE])

num.kds <- nrow(data)-1
num.bgs <- 1


# Define the vector of total site type concentrations:
c.I.tots <- Norm(sitesXcounts[, 1]) * k.c.lib
c.totals <- matrix(c.I.tots, nrow=length(c.I.tots), ncol=5, byrow=FALSE)
colnames(c.totals) <- colnames(data)
rownames(c.totals) <- rownames(data)

# Define the vector of AGO concentrations:
c.agos <- sapply(colnames(data), function(x) {
   as.numeric(x) / 100
  })
print(c.agos)
# PARAMETER INITIALIZATION:
# Initialize starting Kds, which are set to 1/the enrichment of each site type
# in the A12.6 experiment.
r_numerator <- Norm(data[,2]+.001)
r_denomenator <- Norm(c.I.tots+.0001)
pars_init <- c(Logit(r_denomenator/r_numerator,10)[-length(r_numerator)],0,1)
pars <- pars_init
pars[which(is.na(pars))] <- 0

names(pars) <- c(rownames(data)[-nrow(data)], "bg", "AGO")


tick <- 1
print("up to model")

if (sitelist == "11mers") {
  specific_figure_string <- paste0(start, "-", stop, "_", sitelist, "_",
    mirna.start, "-", mirna.stop, "_singlebg_combinedinput_PAPER_final")
} else {
  specific_figure_string <- paste0(start, "-", stop, "_", sitelist, "_",
    "singlebg_combinedinput_PAPER_final")
}



ModelLikelihood <- function(pars) {
  time_prior <- proc.time()

  # Split up the parameters into the kd and background parameters.
  kds  <- c(Logistic(pars[1 : num.kds], 10),1)
  names(kds) <- rownames(data)
  bgs  <- rep(10^pars[(num.kds + 1)], 5)
  B <- bgs[1]
  stock.ago <- 10^pars[num.kds + 2]
  print(stock.ago)
  print(B)
  # Get the bound counts in each experiment:


  f.mat <- matrix(
      rep(sapply(c.agos, function(ago.percent) {
           return(GetFreeAgo(kds, c.I.tots, ago.percent * stock.ago))
         }
         ), nrow(data)), nrow=nrow(data), ncol=ncol(data), byrow=TRUE)

  print("hi")
  A.mat <- matrix(c.agos*stock.ago, nrow=nrow(data), ncol = ncol(data), byrow=TRUE)

  K.mat <- matrix(kds, nrow=nrow(data), ncol = ncol(data), byrow=FALSE)

  l.mat <- matrix(c.I.tots, nrow=nrow(data), ncol = ncol(data), byrow=FALSE)

  R.mat <- matrix(colSums(data), nrow = nrow(data), ncol = ncol(data), byrow=TRUE)

  L <- 100
  # Get the amount of background binding by subtracting the bound from the
  # total sites in each experiment, normalizing. Must transpose to multiply
  # each column.

  c.final <- (-R.mat*l.mat * (f.mat^2 + (L - A.mat)*f.mat + K.mat*B) *
                  (f.mat^3 - (2*A.mat + B - L - K.mat)*f.mat^2 +
                    (A.mat^2 + (B - L)*A.mat - L*B - 
                     (2*A.mat + B - L)*K.mat)*f.mat +
                    (A.mat^2 + (B - L)*A.mat - L*B)*K.mat)^-1)
  colnames(c.final) <- colnames(data)
  # lc.final <- log(c.final + 1)
  # ldata <- log(data + 1)
  sumoflogsquares <- sum((c.final - data)^2)

  if (tick > 1) {
  out <<- rbind(out, c(pars, sumoflogsquares))

    if (sitelist == "11mers") {
      out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                       experiment, "/kds_PAPER/", start, "-", stop, "_", 
                       sitelist, "_", mirna.start, "-", mirna.stop,
                       "_singlebg_combinedinput_PAPER_final5.txt")
      PlotSiteKdOptimization(out[,c(seq(1,num.kds,length=256),(ncol(out)-2):ncol(out))], specific_figure_string, mirna, 256, num.bgs, colors=colors)

    } else {
      out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                       experiment, "/kds_PAPER/", start, "-", stop, "_", 
                       sitelist, 
                       "_singlebg_combinedinput_PAPER_final5.txt")
            PlotSiteKdOptimization(out, specific_figure_string, mirna, num.kds, num.bgs, colors=colors)

    }
    write.table(file=out_file, out, sep="\t", quote=FALSE, row.names=FALSE)
    print(sumoflogsquares)
    plot(c(c.final),c(data), col = rep(c(rgb(0,0,0,alpha =0.4),rgb(1,0,0,alpha=0.4),rgb(0,1,0,alpha=0.4),
                                      rgb(1,1,0,alpha=0.4), rgb(1,0,1, alpha = 0.4)),each=length(c(c.final))/5),log='xy')

  }
  tick <<- tick + 1
  return(sumoflogsquares)
}




Gradient <- function(pars) {
  time_prior <- proc.time()
  # Split up the parameters into the kd and background parameters.
  kds  <- c(Logistic(pars[1 : num.kds], 10),1)
  names(kds) <- rownames(data)
  bgs  <- rep(10^pars[(num.kds + 1)], 5)
  B <- bgs[1]
  stock.ago <- 10^pars[num.kds + 2]

  f.mat <- matrix(
      rep(sapply(c.agos, function(ago.percent) {
           return(GetFreeAgo(kds, c.I.tots, ago.percent * stock.ago))
         }
         ), nrow(data)), nrow=nrow(data), ncol=ncol(data), byrow=TRUE)
  A.mat <- matrix(c.agos*stock.ago, nrow=nrow(data), ncol = ncol(data), byrow=TRUE)
  K.mat <- matrix(kds, nrow=nrow(data), ncol = ncol(data), byrow=FALSE)
  l.mat <- matrix(c.I.tots, nrow=nrow(data), ncol = ncol(data), byrow=FALSE)
  R.mat <- matrix(colSums(data), nrow = nrow(data), ncol = ncol(data), byrow=TRUE)
  L <- 100
  # Get the amount of background binding by subtracting the bound from the
  # total sites in each experiment, normalizing. Must transpose to multiply
  # each column.
  k_8 <- kds["8mer"]
  c.final <- (-R.mat*l.mat * (f.mat^2 + (L - A.mat)*f.mat + K.mat*B) *
                  (f.mat^3 - (2*A.mat + B - L - K.mat)*f.mat^2 +
                    (A.mat^2 + (B - L)*A.mat - L*B - 
                     (2*A.mat + B - L)*K.mat)*f.mat +
                    (A.mat^2 + (B - L)*A.mat - L*B)*K.mat)^-1)
  dF.base <- (1 + colSums(l.mat * K.mat * (f.mat + K.mat)^(-2)))^(-1)

  dF.dK.mat <- sweep(f.mat * l.mat * (f.mat + K.mat)^(-2),MARGIN=2,dF.base, "*")


  C1.mat <- L - B - 2*A.mat


  C2.mat <- A.mat^2 + (B - L)*A.mat - L*B


  # The d (each model point) d Free derivative:)
  dc.ai.dF <- R.mat * l.mat * (
    (
      f.mat^4
    ) + (
      2 * (L - A.mat) * f.mat^3
    ) + (
      ((L - A.mat)^2 + K.mat * (4 * B + A.mat)) * f.mat^2
    ) + (
      2 * K.mat * (A.mat * (L - B) + B * (K.mat - 2 * L) - (A.mat + B)^2) * f.mat
    ) + (
      K.mat * (A.mat^3 + 2 * A.mat^2 * (B - L) - B * (B - L) * (K.mat + L) +
               A.mat * (B^2 - 2 * B * K.mat - 3 * B * L + L^2))
    )
  ) * ((K.mat + f.mat)^2 * (C2.mat + C1.mat * f.mat + f.mat^2)^2)^(-1)




  dc.ai.dKj_specific <- R.mat * l.mat * (
    (
      f.mat^4
    ) + (
      (2 * C1.mat + A.mat) * f.mat^3
    ) + (
      (C2.mat + (C1.mat + A.mat) * C1.mat) * f.mat^2
    ) + (
      (C2.mat * (C1.mat + A.mat)) * f.mat
    )
  ) * ((K.mat + f.mat)^2 * (C2.mat + C1.mat * f.mat + f.mat^2)^2)^(-1)

  dc.ai.dA_specific <- R.mat * l.mat * (
    (
      -f.mat^3
    ) + (
      (2 * (A.mat - L)) * f.mat^2
    ) + (
      ((A.mat - L) * C1.mat - 2 * B * K.mat + C2.mat) * f.mat
    ) + (
      - C1.mat * B * K.mat
    )
  ) * ((K.mat + f.mat) * (C2.mat + C1.mat * f.mat + f.mat^2)^2)^(-1)



  dc.ai.dB <- R.mat * l.mat * (
    (
      -f.mat^3
    ) + (
      (2 * (A.mat - L) - K.mat) * f.mat^2
    ) + (
      ((L - A.mat) * (A.mat - L) - K.mat * (B + C1.mat)) * f.mat
    ) + (
      + A.mat * B * K.mat - B * K.mat * L - K.mat * C2.mat
    )
  ) * ((K.mat + f.mat) * (C2.mat + C1.mat * f.mat + f.mat^2)^2)^(-1)

  lc.final <- log(c.final + 1)
  ldata <- log(data + 1)

  time_new <- proc.time()
  time_prior <- time_new


  residuals <- c.final - data
  # residuals_norm <- residuals * (c.final+1)^(-1)
  grad_derivs <- unlist(mclapply(seq(num.kds), function(index) {
    base <- sweep(dc.ai.dF,MARGIN=2,c(dF.dK.mat[index,]),"*")
    base[index,] <- base[index,] + dc.ai.dKj_specific[index,]
    return(10*exp(pars[index]) * (exp(pars[index]) + 1)^(-2) * 2 * sum(residuals*base))
    }, mc.cores=16))
  time_new <- proc.time()
  print("parallel")
  print(time_new - time_prior)

  dc.ai.dA <- sweep(sweep(dc.ai.dF,MARGIN=2,dF.base,"*") + dc.ai.dA_specific,
                    MARGIN=2,c.agos, "*")

  colnames(c.final) <- colnames(data)


  gradient_with_Ago <- log(10) * stock.ago * 2 * sum(residuals * dc.ai.dA)


  gradient_with_bg <- log(10) * B * 2 * sum(residuals * dc.ai.dB)

  return(c(grad_derivs, gradient_with_bg, gradient_with_Ago))
}

out <- matrix(c(pars, ModelLikelihood(pars)),nrow=1)
pars_update <- pars
colnames(out) <- c(rownames(data)[-length(rownames(data))],
                   "bg", "AGO", "-logp")


  if (sitelist == "11mers") {
    n <- 256
    colors <- c(rainbow_hcl(n, c=50, l=70, start = 0, end = 360*(n-1)/n),"red", "black")
    out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                       experiment, "/kds_PAPER/", start, "-", stop, "_", 
                       sitelist, "_", mirna.start, "-", mirna.stop,
                       "_singlebg_combinedinput_PAPER_final5.txt")

  } else {
    colors = FALSE
    out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                       experiment, "/kds_PAPER/", start, "-", stop, "_", 
                       sitelist, 
                       "_singlebg_combinedinput_PAPER_final5.txt")

  }



print(length(pars))

    solution <- optim(pars,
                    ModelLikelihood,
                    gr = Gradient,
                    method = "L-BFGS-B",
                    lower=c(rep(-20, nrow(data)-1), -3,-0.5),
                    upper=c(rep(4, nrow(data)-1), 0, 1),
                    control = list(maxit=10000))



  pars <- solution$par
  sumoflogsquares_current <- solution$value
  print(sumoflogsquares_current)
  out <- rbind(out, c(pars, sumoflogsquares_current))
  print(out[,(ncol(out)-1)])
warnings()

  write.table(file=out_file, out, sep="\t", quote=FALSE, row.names=FALSE)
# }



