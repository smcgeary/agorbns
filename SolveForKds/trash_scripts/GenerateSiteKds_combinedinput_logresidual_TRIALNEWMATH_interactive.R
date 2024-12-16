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
# # Initial parameters and constants.
# args  <- commandArgs(trailingOnly=TRUE)
# mirna <- args[1]
# # # Region within random sequence from which site types orginiates,
# # # going from position [26 - "start" : 26 + 37 + "stop"]
# start <- as.integer(args[2])
# stop  <- as.integer(args[3])
# # # Parameter specifying which list of sites is used for the analysis.
# # I.E "Current" includes all the current site types, "canonical" is just
# # the six seed-based 6mers, 7mers, and 8mers, along with no site., etc.
# sitelist <- args[4]
# if (sitelist == "12mers") {
#   mirna.start <- as.integer(args[5])
#   mirna.stop <- as.integer(args[6])
# } else {
#   mirna.start <- NULL
#   mirna.stop <- NULL
# }


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


# if (sitelist == "12mers") {
#   sitesXcounts_new <- t(sapply(1:((nrow(sitesXcounts)-1)/4),function(x){colSums(sitesXcounts[1:4+(x-1)*4,])}))
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
kds.init <- ((r_numerator / r_denomenator) + 1)^(-3)
pars.init <- c(Logit(kds.init, max = 10), -0.5, log10(k.c.stockago)+2)


names(pars.init) <- c(rownames(data), "bg", "AGO")
pars <- pars.init

tick <- 1
print("up to model")

if (sitelist == "12mers") {
  specific_figure_string <- paste0(start, "-", stop, "_", sitelist, "_",
    mirna.start, "-", mirna.stop, "_singlebg_combinedinput_logresidual_PAPER_test_test.txt")
} else {
  specific_figure_string <- paste0(start, "-", stop, "_", sitelist, "_",
    "gradient_logresidual_PAPER_test_test")
}



ModelLikelihood <- function(pars, out_model) {
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
  plot(c(lc.final),c(ldata))

  sumoflogsquares <<- sum((lc.final - ldata)^2)

  out <<- rbind(out, c(pars, sumoflogsquares))
  if (tick%%1 == 0) {

    # if (sitelist == "12mers") {
    #   out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
    #                    experiment, "/kds_PAPER/", start, "-", stop, "_", 
    #                    sitelist, "_", mirna.start, "-", mirna.stop,
    #                    "_singlebg_combinedinput_logresidual_PAPER.txt")
    # } else {
    #   out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
    #                    experiment, "/kds_PAPER/", start, "-", stop, "_", 
    #                    sitelist, 
    #                    "_gradient_PAPER.txt")
    # }
    # print(out_file)
    # PlotSiteKdOptimization(out, specific_figure_string, mirna, num.kds, num.bgs, colors=colors)
    # write.table(file=out_file, out, sep="\t", quote=FALSE, row.names=FALSE)
    # print(sumoflogsquares)
  }
  tick <<- tick + 1
  return(sumoflogsquares)
}



ModelFunction <- function(pars) {
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





Gradient <- function(pars, out_model) {
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

  time_new <- proc.time()
  time_prior <- time_new


  # grad_derivs <- sapply(seq(length(kds)), function(index) {
  #   base <- sweep(dc.ai.dF,MARGIN=2,c(dF.dK.mat[index,]),"*")
  #   base[index,] <- base[index,] + dc.ai.dKj_specific[index,]
  #   return(10*exp(pars[index]) * (exp(pars[index]) + 1)^(-2) * 2 * sum((lc.final - ldata)*base*c.final^(-1)))
  #   })

  # time_new <- proc.time()
  # print('no parallel')
  # print(time_new - time_prior)
  # time_prior <- time_new

  grad_derivs <- unlist(mclapply(seq(length(kds)), function(index) {
    base <- sweep(dc.ai.dF,MARGIN=2,c(dF.dK.mat[index,]),"*")
    base[index,] <- base[index,] + dc.ai.dKj_specific[index,]
    return(10*exp(pars[index]) * (exp(pars[index]) + 1)^(-2) * 2 * sum((lc.final - ldata)*base*c.final^(-1)))
    }, mc.cores=16))
  time_new <- proc.time()
  print("parallel")
  print(time_new - time_prior)


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
   Ago_check <- log(10)*stock.ago*2*colSums((log(c.final) - log(data + 1))*dc.ai.dA * c.final^(-1))

  # gradient_with_Ago <- log(10)*stock.ago*2*sum((unlist(c.final) - unlist(data))*dc.ai.dA)

  # print("gradient with Ago")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new


  gradient_with_bg <- log(10)*B*2*sum((log(c.final) - log(data + 1))*dc.ai.dB * c.final^(-1))
  bg_check <- log(10)*B*2*colSums((log(c.final) - log(data + 1))*dc.ai.dB * c.final^(-1))
  # gradient_with_bg <- log(10)*B*2*sum((unlist(c.final) - unlist(data))*dc.ai.dB)

  # print("gradient with bg")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new
  # par(mfrow=c(1,1))
  # plot(pars,c(grad_derivs, gradient_with_bg, gradient_with_Ago),ylim=c(-10,10),col = colors)
  # plot(rep(c(gradient_with_bg, gradient_with_Ago), each = 5), c(bg_check, Ago_check), xlim = c(-10,10),col=c("red","orange","yellow","green","blue"))

  return(c(grad_derivs, gradient_with_bg, gradient_with_Ago))
}

print("529")
# CovarianceMatrix <- function(pars, out_model) {
#   time_prior <- proc.time()
#   # Split up the parameters into the kd and background parameters.
#   kds  <- Logistic(pars[1 : num.kds], 10)
#   names(kds) <- rownames(data)
#   bgs  <- rep(10^pars[(num.kds + 1)], 5)
#   B <- bgs[1]
#   stock.ago <- 10^pars[num.kds + 2]

#   f.mat <- matrix(
#       rep(sapply(c.agos, function(ago.percent) {
#            return(GetFreeAgo(kds, c.I.tots, ago.percent * stock.ago))
#          }
#          ), nrow(data)), nrow=nrow(data), ncol=ncol(data), byrow=TRUE)
#   A.mat <- matrix(c.agos*stock.ago, nrow=nrow(data), ncol = ncol(data), byrow=TRUE)
#   K.mat <- matrix(kds, nrow=nrow(data), ncol = ncol(data), byrow=FALSE)
#   l.mat <- matrix(c.I.tots, nrow=nrow(data), ncol = ncol(data), byrow=FALSE)
#   R.mat <- matrix(colSums(data), nrow = nrow(data), ncol = ncol(data), byrow=TRUE)
#   L <- 100
#   # Get the amount of background binding by subtracting the bound from the
#   # total sites in each experiment, normalizing. Must transpose to multiply
#   # each column.
#   k_8 <- kds["8mer"]
#   c.final <- (-R.mat*l.mat * (f.mat^2 + (L - A.mat)*f.mat + K.mat*B) *
#                   (f.mat^3 - (2*A.mat + B - L - K.mat)*f.mat^2 +
#                     (A.mat^2 + (B - L)*A.mat - L*B - 
#                      (2*A.mat + B - L)*K.mat)*f.mat +
#                     (A.mat^2 + (B - L)*A.mat - L*B)*K.mat)^-1)
#   c.bounds <- as.matrix(
#     sapply(c.agos, function(ago.percent) {
#              return(GetBoundRNA(kds, c.I.tots, ago.percent * stock.ago))
#            }
#            ))


#   c.frees <- c.totals - c.bounds
#   c.bgs <- t(t(c.frees) * B / colSums(c.frees))
#   c.all <- c.bounds + c.bgs
#   c.final.script <- data.frame(t(t(c.all) / colSums(c.all) * colSums(data)))
#   colnames(c.final) <- colnames(data)



#   B_orig <- B
#   B <- B + 0.01

#   c.new.equation <- (-R.mat*l.mat * (f.mat^2 + (L - A.mat)*f.mat + K.mat*B) *
#                   (f.mat^3 - (2*A.mat + B - L - K.mat)*f.mat^2 +
#                     (A.mat^2 + (B - L)*A.mat - L*B - 
#                      (2*A.mat + B - L)*K.mat)*f.mat +
#                     (A.mat^2 + (B - L)*A.mat - L*B)*K.mat)^-1)

#   c.frees <- c.totals - c.bounds
#   c.bgs <- t(t(c.frees) * B / colSums(c.frees))
#   c.all <- c.bounds + c.bgs
#   c.new.script <- data.frame(t(t(c.all) / colSums(c.all) * colSums(data)))


#   c.final_diff.eqn <<- (c.new.equation - c.final)/0.01
#   c.final_diff.script <<- (c.new.script - c.final.script)/0.01

#   c.final <<- c.final
#   # print("up to c.final")
#   # time_new <- proc.time()
#   # print(time_new - time_prior)
#   # time_prior <- time_new

#   B <- B_orig
#   dF.base <- (1 + colSums(l.mat * K.mat * (f.mat + K.mat)^(-2)))^(-1)

#   # print("dF.base")
#   # time_new <- proc.time()
#   # print(time_new - time_prior)
#   # time_prior <- time_new



#   dF.dK.mat <- sweep(f.mat * l.mat * (f.mat + K.mat)^(-2),MARGIN=2,dF.base, "*")

#   # print("dF.dK.mat")
#   # time_new <- proc.time()
#   # print(time_new - time_prior)
#   # time_prior <- time_new


#   C1.mat <- L - B - 2*A.mat

#   # print("C1.mat")
#   # time_new <- proc.time()
#   # print(time_new - time_prior)
#   # time_prior <- time_new


#   C2.mat <- A.mat^2 + (B - L)*A.mat - L*B

#   # print("C2.mat")
#   # time_new <- proc.time()
#   # print(time_new - time_prior)
#   # time_prior <- time_new


#   # The d (each model point) d Free derivative:)
#   # print("dc.ai.dF")
#   # time_new <- proc.time()
#   # print(time_new - time_prior)
#   # time_prior <- time_new

#   dc.ai.dF <- R.mat * l.mat * (
#     (
#       f.mat^4
#     ) + (
#       2 * (L - A.mat) * f.mat^3
#     ) + (
#       ((L - A.mat)^2 + K.mat * (4 * B + A.mat)) * f.mat^2
#     ) + (
#       2 * K.mat * (A.mat * (L - B) + B * (K.mat - 2 * L) - (A.mat + B)^2) * f.mat
#     ) + (
#       K.mat * (A.mat^3 + 2 * A.mat^2 * (B - L) - B * (B - L) * (K.mat + L) +
#                A.mat * (B^2 - 2 * B * K.mat - 3 * B * L + L^2))
#     )
#   ) * ((K.mat + f.mat)^2 * (C2.mat + C1.mat * f.mat + f.mat^2)^2)^(-1)




#   dc.ai.dKj_specific <- R.mat * l.mat * (
#     (
#       f.mat^4
#     ) + (
#       (2 * C1.mat + A.mat) * f.mat^3
#     ) + (
#       (C2.mat + (C1.mat + A.mat) * C1.mat) * f.mat^2
#     ) + (
#       (C2.mat * (C1.mat + A.mat)) * f.mat
#     )
#   ) * ((K.mat + f.mat)^2 * (C2.mat + C1.mat * f.mat + f.mat^2)^2)^(-1)

#   # print("dc.ai.dKj_specific")
#   # time_new <- proc.time()
#   # print(time_new - time_prior)
#   # time_prior <- time_new



#   dc.ai.dA_specific <- R.mat * l.mat * (
#     (
#       -f.mat^3
#     ) + (
#       (2 * (A.mat - L)) * f.mat^2
#     ) + (
#       ((A.mat - L) * C1.mat - 2 * B * K.mat + C2.mat) * f.mat
#     ) + (
#       - C1.mat * B * K.mat
#     )
#   ) * ((K.mat + f.mat) * (C2.mat + C1.mat * f.mat + f.mat^2)^2)^(-1)

#   # print("dc.ai.dA_specific")
#   # time_new <- proc.time()
#   # print(time_new - time_prior)
#   # time_prior <- time_new


#   dc.ai.dB <- R.mat * l.mat * (
#     (
#       -f.mat^3
#     ) + (
#       (2 * (A.mat - L) - K.mat) * f.mat^2
#     ) + (
#       ((L - A.mat) * (A.mat - L) - K.mat * (B + C1.mat)) * f.mat
#     ) + (
#       + A.mat * B * K.mat - B * K.mat * L - K.mat * C2.mat
#     )
#   ) * ((K.mat + f.mat) * (C2.mat + C1.mat * f.mat + f.mat^2)^2)^(-1)
#   dc.ai.dB <<- dc.ai.dB


#   lc.final <- log(c.final)
#   ldata <- log(data + 1)
#   grad_derivs <- matrix(sapply(seq(length(kds)), function(index) {
#     base <- sweep(dc.ai.dF,MARGIN=2,c(dF.dK.mat[index,]),"*")
#     base[index,] <- base[index,] + dc.ai.dKj_specific[index,]
#     return( c(10*exp(-pars[index]) * (exp(-pars[index]) + 1)^(-2) * base * c.final^(-1)))
#     }), nrow = length(data), ncol = length(kds), byrow = FALSE)



#   grad_derivs <<- grad_derivs
#   dc.ai.dA <- sweep(sweep(dc.ai.dF,MARGIN=2,dF.base,"*") + dc.ai.dA_specific,
#                     MARGIN=2,c.agos, "*")



#   colnames(c.final) <- colnames(data)
#   use <- which(c.final > 0 & data > 0)

#   # print("use")
#   # time_new <- proc.time()
#   # print(time_new - time_prior)
#   # time_prior <- time_new

#   # sumoflogsquares <<- sum((log(c.final) - log(data + 1))^2)
#   # sumoflogsquares <<- sum((unlist(c.final) - unlist(data))^2)

#   # print("logsquares")
#   # time_new <- proc.time()
#   # print(time_new - time_prior)
#   # time_prior <- time_new


#   gradient_with_Ago <- c(log(10)*stock.ago*dc.ai.dA * c.final^(-1))

#   gradient_with_bg <- log(10)*B*2*dc.ai.dB * c.final^(-1)
#   gradient_with_bg <<- gradient_with_bg
#   gradient_with_bg <- c(gradient_with_bg)
#   print(gradient_with_bg[1:10])
#   Matrix_all <- cbind(grad_derivs, gradient_with_bg, gradient_with_Ago)
#   Matrix_all <<- Matrix_all
#   sumoflogsquares_norm <<- sum((lc.final - ldata)^2) * (length(data) - length(par))^(-1)
#   V_matrix <- (t(Matrix_all)%*%Matrix_all)^(-1)*sumoflogsquares_norm
#   return(list(V_matrix,lc.final[2,1],grad_derivs[2,1]))
# }

# check <- Gradient(pars)




# # out_AGO <- c(0,0,0)
# #  for (i in seq(pars["AGO"]-0.1,pars["AGO"]+0.1,length=100)) {
# #   pars_temp <- pars
# #   pars_temp["AGO"] <- i
# #   out_AGO <- rbind(out_AGO,ModelLikelihood(pars_temp, out_AGO))
# #  }

#  out_bg <- c(0,0,0)
#  for (i in seq(pars["bg"]-0.1,pars["bg"]+0.1,length=300)) {
#   pars_temp <- pars
#   pars_temp["bg"] <- i
#   out_bg <- rbind(out_bg,ModelLikelihood(pars_temp, out_bg))
#  }
# out_bg <- out_bg[-1, ]

# bg_ <- out_bg[, 1]
# Grad_bg <- out_bg[-nrow(out_bg), 3]
# dOBJ_bg <- out_bg[-1, 2] - out_bg[-nrow(out_bg), 2]
# dbg <- bg_[-1] - bg_[-length(bg_)]
# dOBJ.dbg <- dOBJ_bg / dbg
# # dFree2.dA <- dFree2 / dA
# graphics.off()
# dev.new(xpos = 20, ypos = 20, width = 15)
# par(mfrow=c(1, 3))

# plot(K_8mer[-length(K_8mer)],Grad_8mer,col="gray",type="l",lwd = 5)
# lines(K_8mer[-length(K_8mer)], dOBJ.dK_8mer,lwd=1, lty = 2)

# plot(A_[-length(A_)], Grad_AGO, col="gray",type="l",lwd = 5)
# lines(A_[-length(A_)], dOBJ.dA,lwd=1, lty = 2)


# plot(bg_[-length(bg_)], Grad_bg, col="gray",type="l",lwd = 5,ylim=c(0, max(Grad_bg,dOBJ.dbg)))
# lines(bg_[-length(bg_)], dOBJ.dbg,lwd=1, lty = 2)

out <- rbind(c(pars, 10000000), c(pars, 10000000))
pars_update <- pars
colnames(out) <- c(rownames(data),
                   "bg", "AGO", "-logp")


# ## Optimization of parameters:
# for (i in seq(1, 2000)) {
#   print(i)
#   # Define the "scale" vector to be used as the parscale vector.
#   sumoflogsquares_round <- ModelLikelihood(pars, out)
#   scale <- sapply(seq(1, length(pars)), function(par_i) {
#     if (pars[par_i] > 10) {
#       pars[par_i] <<- 10
#     }
#     temp <- pars
#     temp[par_i] <- temp[par_i] - 1
#     func_out <- abs(ModelLikelihood(temp, out) - sumoflogsquares_round)
#     return(func_out)
#   })
#   scale[is.na(scale)] <- 0
#   scale <- scale + 0.001
#   scale <- abs(scale)
#   names(scale) <- names(pars)

  # solution <- optim(pars,
  #                   ModelLikelihood,
  #                   out_model = out,
  #                   method = "Nelder-Mead",
  #                   control = list(maxit=10000))

  if (sitelist == "12mers") {
    n <- num.kds - 1
    colors <- c(rainbow_hcl(n, c=50, l=70, start = 0, end = 360*(n-1)/n), "black", "gray", "red")
    out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                       experiment, "/kds_PAPER/", start, "-", stop, "_", 
                       sitelist, "_", mirna.start, "-", mirna.stop,
                       "_gradient_PAPER_final.txt")

  } else {
    colors = FALSE
    out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                       experiment, "/kds_PAPER/", start, "-", stop, "_", 
                       sitelist, 
                       "_gradient_PAPER_final.txt")

  }

print("841")



    solution <- optim(pars,
                    ModelLikelihood,
                    gr = Gradient,
                    out_model = out,
                    method = "L-BFGS-B",
                    lower=rep(-40,length(pars)),
                    upper=rep(10, length(pars)),
                    control = list(maxit=10000))
  pars.final <- solution$par




pars_bootstrap <- pars.final
break
print(pars_bootstrap)
  time_start <- proc.time()
data_original <- data

  Model_error <- (diag(Vmat.out)/length(data))^(1/2)
for (i in seq(100)) {




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
                    method = "L-BFGS-B",
                    lower=rep(-20,length(pars)),
                    upper=rep(4, length(pars)),
                    control = list(maxit=10000))

    pars_bootstrap <- rbind(pars_bootstrap, solution$par)
    print(apply(pars_bootstrap,2,sd))
    plot(apply(pars_bootstrap,2,sd),Model_error, xlim = c(0, 0.3), ylim = c(0, 0.3))
}

  time_end <- proc.time()
  Vmat.out <- CovarianceMatrix(pars)[[1]]
  Model_error <- (diag(Vmat.out)/length(data))^(1/2)
  pars <- solution$par
  sumoflogsquares_current <- solution$value
  print("final sum of squares")
  print(sumoflogsquares_current)

# parm_8mer_input <- seq(pars["8mer"] - 1, pars["8mer"] + 0.5, length= 300)
# out_8mer <- c()
# deriv_8mer <- c()

# for (i in parm_8mer_input) {
#   pars_temp <- pars
#   pars_temp["8mer"] <- i
#   print(pars_temp["8mer"])
#   out_8mer <- c(out_8mer,CovarianceMatrix(pars_temp)[[2]])
#   deriv_8mer <- c(deriv_8mer,CovarianceMatrix(pars_temp)[[3]])

#  }

# d_out_8mer <- c(out_8mer)[-1] - c(out_8mer)[-length(out_8mer)]
# d_parm <- parm_8mer_input[-1] - parm_8mer_input[-length(parm_8mer_input)]
# d_manual <- d_out_8mer/ d_parm

# plot(d_manual,deriv_8mer[-length(deriv_8mer)],type="l",lwd=5,col="red")
# segments(min(deriv_8mer),min(deriv_8mer),max(deriv_8mer), max(deriv_8mer))

#   out <- rbind(out, c(pars, sumoflogsquares_current))

#   print("time taken:")
#   print(time_end[3] - time_start[3])

#   write.table(file=out_file, out, sep="\t", quote=FALSE, row.names=FALSE)
# # }



