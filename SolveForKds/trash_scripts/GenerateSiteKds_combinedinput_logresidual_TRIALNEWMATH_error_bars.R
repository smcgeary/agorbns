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
# # Initial parameters and constants.
args  <- commandArgs(trailingOnly=TRUE)
mirna <- args[1]
# # Region within random sequence from which site types orginiates,
# # going from position [26 - "start" : 26 + 37 + "stop"]
start <- as.integer(args[2])
stop  <- as.integer(args[3])
# # Parameter specifying which list of sites is used for the analysis.
# I.E "Current" includes all the current site types, "canonical" is just
# the six seed-based 6mers, 7mers, and 8mers, along with no site., etc.
sitelist <- args[4]
if (sitelist == "12mers") {
  mirna.start <- as.integer(args[5])
  mirna.stop <- as.integer(args[6])
} else {
  mirna.start <- NULL
  mirna.stop <- NULL
}


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
k.c.stockago = stockago[mirna,experiment]
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
if (sitelist == "12mers") {
  seqs <- rownames(sitesXcounts)
} else {
  seqs <- sitesXcounts[,1]
  sitesXcounts <- sitesXcounts[,-1]
}

# sitesXcounts_new <- t(sapply(1:((nrow(sitesXcounts)-1)/64),function(x){colSums(sitesXcounts[1:64+(x-1)*64,])}))
# sitesXcounts_new <- rbind(sitesXcounts_new,sitesXcounts[nrow(sitesXcounts),,drop=FALSE])
# sitesXcounts <- sitesXcounts_new

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


num.kds <- nrow(data)
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

# PARAMETER INITIALIZATION:
# Initialize starting Kds, which are set to 1/the enrichment of each site type
# in the A12.6 experiment.
r_numerator <- Norm(data[,2])
r_denomenator <- Norm(c.I.tots)
kds.init <- ((r_numerator / r_denomenator) + 1)^(-1)
pars.init <- c(Logit(kds.init, max = 10), -1, log10(k.c.stockago))
names(pars.init) <- c(rownames(data), "bg", "AGO")
pars <- pars.init

tick <- 1
print("up to model")

if (sitelist == "12mers") {
  specific_figure_string <- paste0(start, "-", stop, "_", sitelist, "_",
    mirna.start, "-", mirna.stop, "_singlebg_combinedinput_logresidual_PAPER.txt")
} else {
  specific_figure_string <- paste0(start, "-", stop, "_", sitelist, "_",
    "gradient_logresidual_PAPER")
}



ModelLikelihood <- function(pars, out_model, plot = FALSE) {
  time_prior <- proc.time()
  pars_update <<- rbind(pars_update, pars)

  # Split up the parameters into the kd and background parameters.
  kds  <- Logistic(pars[1 : num.kds], 10)
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

  c.final <- (-R.mat*l.mat * (f.mat^2 + (L - A.mat)*f.mat + K.mat*B) *
                  (f.mat^3 - (2*A.mat + B - L - K.mat)*f.mat^2 +
                    (A.mat^2 + (B - L)*A.mat - L*B - 
                     (2*A.mat + B - L)*K.mat)*f.mat +
                    (A.mat^2 + (B - L)*A.mat - L*B)*K.mat)^-1)

  colnames(c.final) <- colnames(data)
  lc.final <- log(c.final)
  ldata <- log(data + 1)
  sumoflogsquares <<- sum((lc.final - ldata)^2)

  out <<- rbind(out, c(pars, sumoflogsquares))

  tick <<- tick + 1
  return(sumoflogsquares)
}



ModelFunction <- function(pars, plot = FALSE) {
  time_prior <- proc.time()
  # Split up the parameters into the kd and background parameters.
  kds  <- Logistic(pars[1 : num.kds], 10)
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
  return(c.final)
}


DerivativeManualMOdelFunction <- function(pars) {
    c.final <- ModelFunction(pars)
    return(c.final[1,1])
}

f.MAT <- function(pars){
    kds  <- Logistic(pars[1 : num.kds], 10)
  names(kds) <- rownames(data)
  bgs  <- rep(10^pars[(num.kds + 1)], 5)
  B <- bgs[1]
  stock.ago <- 10^pars[num.kds + 2]
  f.mat <- matrix(
      rep(sapply(c.agos, function(ago.percent) {
           return(GetFreeAgo(kds, c.I.tots, ago.percent * stock.ago))
         }
         ), nrow(data)), nrow=nrow(data), ncol=ncol(data), byrow=TRUE)
  return(f.mat)
}

A.MAT <- function(pars){
    kds  <- Logistic(pars[1 : num.kds], 10)
  names(kds) <- rownames(data)
  bgs  <- rep(10^pars[(num.kds + 1)], 5)
  B <- bgs[1]
  stock.ago <- 10^pars[num.kds + 2]
  A.mat <- matrix(c.agos*stock.ago, nrow=nrow(data), ncol = ncol(data), byrow=TRUE)
  return(A.mat)
}

A.MAT <- function(pars){
    kds  <- Logistic(pars[1 : num.kds], 10)
  names(kds) <- rownames(data)
  bgs  <- rep(10^pars[(num.kds + 1)], 5)
  B <- bgs[1]
  stock.ago <- 10^pars[num.kds + 2]
  A.mat <- matrix(c.agos*stock.ago, nrow=nrow(data), ncol = ncol(data), byrow=TRUE)
  return(A.mat)
}

K.MAT <- function(pars){
    kds  <- Logistic(pars[1 : num.kds], 10)
  names(kds) <- rownames(data)
  bgs  <- rep(10^pars[(num.kds + 1)], 5)
  B <- bgs[1]
  stock.ago <- 10^pars[num.kds + 2]
  K.mat <- matrix(kds, nrow=nrow(data), ncol = ncol(data), byrow=FALSE)
  return(K.mat)
}

l.MAT <- function(pars){
    kds  <- Logistic(pars[1 : num.kds], 10)
  names(kds) <- rownames(data)
  bgs  <- rep(10^pars[(num.kds + 1)], 5)
  B <- bgs[1]
  stock.ago <- 10^pars[num.kds + 2]
  l.mat <- matrix(c.I.tots, nrow=nrow(data), ncol = ncol(data), byrow=FALSE)
  return(l.mat)
}

R.MAT <- function(pars){
    kds  <- Logistic(pars[1 : num.kds], 10)
  names(kds) <- rownames(data)
  bgs  <- rep(10^pars[(num.kds + 1)], 5)
  B <- bgs[1]
  stock.ago <- 10^pars[num.kds + 2]
  R.mat <- matrix(colSums(data), nrow = nrow(data), ncol = ncol(data), byrow=TRUE)
  return(R.mat)
}

c.final.MAT <- function(pars) {
  f.mat <- f.MAT(pars)
  A.mat <- A.MAT(pars)
  K.mat <- K.MAT(pars)
  l.mat <- l.MAT(pars)
  R.mat <- R.MAT(pars)
  L=100
    B <- 10^pars[(num.kds + 1)]

    c.final <- (-R.mat*l.mat * (f.mat^2 + (L - A.mat)*f.mat + K.mat*B) *
                  (f.mat^3 - (2*A.mat + B - L - K.mat)*f.mat^2 +
                    (A.mat^2 + (B - L)*A.mat - L*B - 
                     (2*A.mat + B - L)*K.mat)*f.mat +
                    (A.mat^2 + (B - L)*A.mat - L*B)*K.mat)^-1)
    return(c.final)
}


df.base.MAT <- function(pars){
    f.mat <- f.MAT(pars)
  K.mat <- K.MAT(pars)
  l.mat <- l.MAT(pars)

  dF.base <- (1 + colSums(l.mat * K.mat * (f.mat + K.mat)^(-2)))^(-1)
  return(dF.base)
}





Gradient <- function(pars, out_model, plot = FALSE) {
  time_prior <- proc.time()
  # Split up the parameters into the kd and background parameters.
  kds  <- Logistic(pars[1 : num.kds], 10)
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
  c.final <<- c.final
  # print("up to c.final")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new


  dF.base <- (1 + colSums(l.mat * K.mat * (f.mat + K.mat)^(-2)))^(-1)

  # print("dF.base")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new



  dF.dK.mat <- sweep(f.mat * l.mat * (f.mat + K.mat)^(-2),MARGIN=2,dF.base, "*")

  # print("dF.dK.mat")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new


  C1.mat <- L - B - 2*A.mat

  # print("C1.mat")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new


  C2.mat <- A.mat^2 + (B - L)*A.mat - L*B

  # print("C2.mat")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new


  # The d (each model point) d Free derivative:)
  # print("dc.ai.dF")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new

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

  # print("dc.ai.dKj_specific")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new



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

  # print("dc.ai.dA_specific")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new


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

  # print("dc.ai.dB")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new



  # grad_derivs <- sapply(seq(length(kds)), function(index) {
  #   base <- sweep(dc.ai.dF,MARGIN=2,c(dF.dK.mat[index,]),"*")
  #   base[index,] <- base[index,] + dc.ai.dKj_specific[index,]
  #   return(10*exp(pars[index]) * (exp(pars[index]) + 1)^(-2) * 2 * sum((log(c.final) - log(data+1)*base*c.final^(-1))))
  #   })
  lc.final <- log(c.final)
  ldata <- log(data + 1)
  residuals <- log((c.final) / (data + 1)) / (c.final)
  # grad_derivs_original <- sapply(seq(length(kds)), function(index) {
  #   base <- sweep(dc.ai.dF,MARGIN=2,c(dF.dK.mat[index,]),"*")
  #   base[index,] <- base[index,] + dc.ai.dKj_specific[index,]
  #   return(10*exp(pars[index]) * (exp(pars[index]) + 1)^(-2) * 2 * sum((lc.final - ldata)*base*c.final^(-1)))
    # })

    grad_derivs <- 2*(10*exp(pars[1:num.kds]) * (exp(pars[1:num.kds]) + 1)^(-2) *
                  (colSums(residuals*dc.ai.dF) %*%
                   t(dF.dK.mat
                  ) + rowSums(dc.ai.dKj_specific*residuals)))





    # print("grad deriv")

    # time_new <- proc.time()
    # print(time_new - time_prior)
    # time_prior <- time_new


  dc.ai.dA <- sweep(sweep(dc.ai.dF,MARGIN=2,dF.base,"*") + dc.ai.dA_specific,
                    MARGIN=2,c.agos, "*")



  # print("AGO deriv")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new



  colnames(c.final) <- colnames(data)
  use <- which(c.final > 0 & data > 0)

  # print("use")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new

  # sumoflogsquares <<- sum((log(c.final) - log(data + 1))^2)
  # sumoflogsquares <<- sum((unlist(c.final) - unlist(data))^2)

  # print("logsquares")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new


  gradient_with_Ago <- log(10)*stock.ago*2*sum((log(c.final) - log(data + 1))*dc.ai.dA * c.final^(-1))
  # gradient_with_Ago <- log(10)*stock.ago*2*sum((unlist(c.final) - unlist(data))*dc.ai.dA)

  # print("gradient with Ago")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new


  gradient_with_bg <- log(10)*B*2*sum((log(c.final) - log(data + 1))*dc.ai.dB * c.final^(-1))
  # gradient_with_bg <- log(10)*B*2*sum((unlist(c.final) - unlist(data))*dc.ai.dB)

  # print("gradient with bg")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new


  return(c(grad_derivs, gradient_with_bg, gradient_with_Ago))
}

print("529")

out <- rbind(c(pars, 10000000), c(pars, 10000000))
pars_update <- pars
colnames(out) <- c(rownames(data),
                   "bg", "AGO", "-logp")



  if (sitelist == "12mers") {
    n <- num.kds
    colors <- c(rainbow_hcl(n, c=50, l=70, start = 0, end = 360*(n-1)/n), "red", "gray")
    out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                       experiment, "/kds_PAPER/", start, "-", stop, "_", 
                       sitelist, "_", mirna.start, "-", mirna.stop,
                       "_singlebg_combinedinput_logresidual_PAPER.txt")

  } else {
    colors = FALSE
    out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                       experiment, "/kds_PAPER/", start, "-", stop, "_", 
                       sitelist, 
                       "_gradient_PAPER.txt")

  }

print(pars)
print("841")

    solution <- optim(pars.init,
                    ModelLikelihood,
                    gr = Gradient,
                    out_model = out,
                    plot = TRUE,
                    method = "L-BFGS-B",
                    lower=rep(-20,length(pars)),
                    upper=rep(2, length(pars)),
                    control = list(maxit=10000))

  pars.final <- solution$par
pars_bootstrap <- pars.final
print(pars_bootstrap)
  time_start <- proc.time()
data_original <- data
c.I.tots_orig <- c.I.tots

for (i in seq(500)) {
  print(i)

c.I.tots <- Norm(rmultinom(1,sum(sitesXcounts[, 1]),prob=sitesXcounts[,1])+1) * k.c.lib
# print(c.I.tots)
# print(c.I.tots_orig)
c.totals <- matrix(c.I.tots, nrow=length(c.I.tots), ncol=5, byrow=FALSE)
colnames(c.totals) <- colnames(data)
rownames(c.totals) <- rownames(data)


data_new <- matrix(apply(data_original,2, function(column) {
  rmultinom(1, sum(column),prob=column)
  }), nrow=nrow(data_original), ncol = ncol (data_original), byrow=FALSE)
rownames(data_new) <- rownames(data)
colnames(data_new) <- colnames(data)
data <- data_new

    solution <- optim(pars.init,
                    ModelLikelihood,
                    gr = Gradient,
                    out_model = out,
                    plot = FALSE,
                    method = "L-BFGS-B",
                    lower=rep(-20,length(pars)),
                    upper=rep(2, length(pars)),
                    control = list(maxit=10000))

    pars_bootstrap <- rbind(pars_bootstrap, solution$par)
}

  pars_sorted <- apply(pars_bootstrap,2,sort)
  pars_upper <- pars_sorted[25,]
  pars_lower <- pars_sorted[475,]
  pars_final <- cbind(pars.final, pars_upper, pars_lower)
  sumoflogsquares_current <- solution$value
  print("final sum of squares")
  print(sumoflogsquares_current)

  out_file_sorted <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                         experiment, "/kds_PAPER/", start, "-", stop, "_", 
                         sitelist,
                         "_errorall_PAPER.txt")


  out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                         experiment, "/kds_PAPER/", start, "-", stop, "_", 
                         sitelist,
                         "_error_PAPER.txt")
  print(out_file)
    write.table(file=out_file, pars_final, sep="\t", quote=FALSE, row.names=TRUE)
    write.table(file=out_file_sorted, t(pars_sorted), sep="\t", quote=FALSE, row.names=TRUE)

