################################################################################
#GenerateSiteKds_combinedinput_PAPER.py
################################################################################

## THIS IS THE OFFICIAL SCRIPT FOR FINALIZING THE MANUSCRIPT.
## THIS SCRIPT HAS THESE DECISIONS MADE:
## 1.) PROTEIN IS FIT AS A PARAMETER
## 2.) A SINGLE BACKGROUND PARAMETER IS GIVEN TO ALL FIVE CONCENTRATIONS.
## 3.) SITELIST, START, AND STOP POSITION ARE MANDATORY PARAMETERS.
# ## 4.) ALL INPUTS ARE COMBINED ACROSS FOR THESE ANALYSES.
library(numDeriv)
# # Initial parameters and constants.
# args  <- commandArgs(trailingOnly=TRUE)
# mirna <- args[1]
# # Region within random sequence from which site types orginiates,
# # going from position [26 - "start" : 26 + 37 + "stop"]
# start <- as.integer(args[2])
# stop  <- as.integer(args[3])
mirna.start <- 2
mirna.stop <- 5
# # # Parameter specifying which list of sites is used for the analysis.
# # # I.E "Current" includes all the current site types, "canonical" is just
# # # the six seed-based 6mers, 7mers, and 8mers, along with no site., etc.
# sitelist <- args[4]
# # Experiment name
# experiment <- "equilibrium"

# Loads general functions used in AgoRBNS analysis.
source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")

# Loads the colors associated with each site type, for plotting purposes.
site_cols <- read.table(
  "/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/site_info.txt",
  row.names=1, header=FALSE, sep="\t", stringsAsFactors=FALSE)[, 1, drop=FALSE]

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

# sitesXcounts_new <- t(sapply(1:((nrow(sitesXcounts)-1)/4),function(x){colSums(sitesXcounts[1:4+(x-1)*4,])}))
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
ModelLikelihood <- function(pars, out_model) {
  time_prior <- proc.time()
  # pars_update <<- rbind(pars_update, pars)

  # Split up the parameters into the kd and background parameters.
  kds  <- Logistic(pars[1 : num.kds], 10)
  names(kds) <- rownames(data)
  bgs  <- rep(10^pars[(num.kds + 1)], 5)
  B <- bgs[1]
  stock.ago <- 10^pars[num.kds + 2]
  # print(pars[num.kds + 2])
  # print(stock.ago)
  # print(sapply(c.agos, function(ago.percent) {
  #          return(GetFreeAgo(kds, c.I.tots, ago.percent * stock.ago))}))
  # Get the bound counts in each experiment:

  # print("first time")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new

  f.mat <- matrix(
      rep(sapply(c.agos, function(ago.percent) {
           return(GetFreeAgo(kds, c.I.tots, ago.percent * stock.ago))
         }
         ), nrow(data)), nrow=nrow(data), ncol=ncol(data), byrow=TRUE)

  # print("f.mat")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new

  A.mat <- matrix(c.agos*stock.ago, nrow=nrow(data), ncol = ncol(data), byrow=TRUE)

  # print("A.mat")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new



  K.mat <- matrix(kds, nrow=nrow(data), ncol = ncol(data), byrow=FALSE)

  #   print("K.mat")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new


  l.mat <- matrix(c.I.tots, nrow=nrow(data), ncol = ncol(data), byrow=FALSE)

  # print("l.mat")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new

  R.mat <- matrix(colSums(data), nrow = nrow(data), ncol = ncol(data), byrow=TRUE)

  # print("R.mat")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new



  L <- 100
  # Get the amount of background binding by subtracting the bound from the
  # total sites in each experiment, normalizing. Must transpose to multiply
  # each column.






  # c.final <<- c.final
  c.final <- (-R.mat*l.mat * (f.mat^2 + (L - A.mat)*f.mat + K.mat*B) *
                  (f.mat^3 - (2*A.mat + B - L - K.mat)*f.mat^2 +
                    (A.mat^2 + (B - L)*A.mat - L*B - 
                     (2*A.mat + B - L)*K.mat)*f.mat +
                    (A.mat^2 + (B - L)*A.mat - L*B)*K.mat)^-1)
  print("make c.final_alt")
  time_new <- proc.time()
  print(time_new - time_prior)
  time_prior <- time_new


  colnames(c.final) <- colnames(data)
  sumoflogsquares <<- sum((c.final - data)^2)

  # if (tick%%1 == 0) {
  #   # print(pars)
  #   # print(kds)
  #   # print(stock.ago)
  #   # print(bgs[1])
  #   print(sumoflogsquares)
  #   dev.set(5)
  #   plot(c(0.01, 10000000), c(0.01, 10000000), log = 'xy')
  #   sapply(1:ncol(c.final), function(column) {
  #     points(c.final[,column], data[,column])
  #   })
  # }
  tick <<- tick + 1
  # plot(seq(length(pars)),pars)
  return(sumoflogsquares)
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
    c.bounds <- as.matrix(
    sapply(c.agos, function(ago.percent) {
             return(GetBoundRNA(kds, c.I.tots, ago.percent * stock.ago))
           }
           ))

  c.final <- (-R.mat*l.mat * (f.mat^2 + (L - A.mat)*f.mat + K.mat*B) *
                  (f.mat^3 - ((L - A.mat) - (A.mat + B) + K.mat)*f.mat^2 +
                    (A.mat^2 + (B - L)*A.mat - L*B + 
                     (L - 2*A.mat - B)*K.mat)*f.mat +
                    (A.mat^2 + (B - L)*A.mat - L*B)*K.mat)^-1)

  c.final <- -R.mat * l.mat * (
    (
      f.mat^2
    ) + (
      (L - A.mat) * f.mat
    ) + (
      K.mat * B
    )
  ) * (
    (
      f.mat^3
    ) + (
      ((L - A.mat) - (A.mat + B) + K.mat) * f.mat^2
    ) + (
      (K.mat * ((L - A.mat) - (A.mat + B)) - (L - A.mat) * (A.mat + B)) * f.mat 
    ) + (
      -K.mat * (L - A.mat) * (A.mat + B)
    )
    )^(-1)

  c.frees <- c.totals - c.bounds
  c.bgs <- t(t(c.frees) * bgs / colSums(c.frees))
  c.all <- c.bounds + c.bgs
  c.final <- data.frame(t(t(c.all) / colSums(c.all) * colSums(data)))


  c.final <<- c.final
  print("up to c.final")
  time_new <- proc.time()
  print(time_new - time_prior)
  time_prior <- time_new


  dF.base <- (1 + colSums(l.mat * K.mat * (f.mat + K.mat)^(-2)))^(-1)

  print("dF.base")
  time_new <- proc.time()
  print(time_new - time_prior)
  time_prior <- time_new



  dF.dK.mat <- sweep(f.mat * l.mat * (f.mat + K.mat)^(-2),MARGIN=2,dF.base, "*")

  print("dF.dK.mat")
  time_new <- proc.time()
  print(time_new - time_prior)
  time_prior <- time_new


  C1.mat <- L - B - 2*A.mat

  print("C1.mat")
  time_new <- proc.time()
  print(time_new - time_prior)
  time_prior <- time_new


  C2.mat <- A.mat^2 + (B - L)*A.mat - L*B

  print("C2.mat")
  time_new <- proc.time()
  print(time_new - time_prior)
  time_prior <- time_new


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

  check_function <- function(f.mat) { 
               out <- (-R.mat*l.mat * (f.mat^2 + (L - A.mat)*f.mat + K.mat*B) *
                  (f.mat^3 - (2*A.mat + B - L - K.mat)*f.mat^2 +
                    (A.mat^2 + (B - L)*A.mat - L*B - 
                     (2*A.mat + B - L)*K.mat)*f.mat +
                    (A.mat^2 + (B - L)*A.mat - L*B)*K.mat)^-1)
               return(out[1,1])}

  check_dFdparm <- function(parm){
    kds[1] <- parm
    out <- sapply(c.agos, function(ago.percent) {
           return(GetFreeAgo(kds, c.I.tots, ago.percent * stock.ago))
           })
    return(out[1])
  }
    check_dFdAGO <- function(stock.ago){
    out <- sapply(c.agos, function(ago.percent) {
           return(GetFreeAgo(kds, c.I.tots, ago.percent * stock.ago))
           })
    return(out[1])
  }


  f_4 <- (R.mat * l.mat *(f.mat^4)* ((K.mat + f.mat)^2 * (C2.mat + C1.mat * f.mat + f.mat^2)^2)^(-1))[1,1]
  f_3 <- (R.mat * l.mat *(2 * (L - A.mat) * f.mat^3)* ((K.mat + f.mat)^2 * (C2.mat + C1.mat * f.mat + f.mat^2)^2)^(-1))[1,1]
  f_2 <- (R.mat * l.mat *(((L - A.mat)^2 + K.mat * (2 * (L - B) - A.mat)) * f.mat^2)* ((K.mat + f.mat)^2 * (C2.mat + C1.mat * f.mat + f.mat^2)^2)^(-1))[1,1]
  f_1 <- (R.mat * l.mat *(2 * K.mat * (L * A.mat + B * K.mat - 2 * A.mat * B - (A.mat + B)^2) * f.mat)* ((K.mat + f.mat)^2 * (C2.mat + C1.mat * f.mat + f.mat^2)^2)^(-1))[1,1]
  f_0 <- (R.mat * l.mat *(  K.mat * (B * (C2.mat + C1.mat * K.mat) - C2.mat * (L - A.mat) ))* ((K.mat + f.mat)^2 * (C2.mat + C1.mat * f.mat + f.mat^2)^2)^(-1))[1,1]

  f_grad_numeric <- grad(check_dFdAGO, stock.ago)
  print(c.agos * stock.ago)
  print(sapply(c.agos, function(ago.percent) {
           return(GetFreeAgo(kds, c.I.tots, ago.percent * stock.ago))}))
  # original
  # dc.ai.dF <- R.mat * l.mat * (
  #   (
  #     f.mat^4
  #   ) + (
  #     2 * (L - A.mat) * f.mat^3
  #   ) + (
  #     ((L - A.mat) * (C1.mat + K.mat) - (C2.mat + C1.mat * K.mat) + 3 * B * K.mat) * f.mat^2
  #   ) + (
  #     2 * K.mat * (B * (C1.mat + K.mat) - C2.mat) * f.mat
  #   ) + (
  #     K.mat * (B * (C2.mat + C1.mat * K.mat) - C2.mat * (L - A.mat)) 
  #   )
  # ) * ((K.mat + f.mat)^2 * (C2.mat + C1.mat * f.mat + f.mat^2)^2)^(-1)







  print("dc.ai.dF")
  time_new <- proc.time()
  print(time_new - time_prior)
  time_prior <- time_new



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

  print("dc.ai.dKj_specific")
  time_new <- proc.time()
  print(time_new - time_prior)
  time_prior <- time_new



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

  print("dc.ai.dA_specific")
  time_new <- proc.time()
  print(time_new - time_prior)
  time_prior <- time_new


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

  print("dc.ai.dB")
  time_new <- proc.time()
  print(time_new - time_prior)
  time_prior <- time_new


  f_grad <<- c.agos[1]*dF.base[1]
  grad_derivs <- sapply(seq(length(kds)), function(index) {

    base <- sweep(dc.ai.dF,MARGIN=2,c(dF.dK.mat[index,]),"*")

    if(index==1)
    base[index,] <- base[index,] + dc.ai.dKj_specific[index,]



    return(10*exp(pars[index]) * (exp(pars[index]) + 1)^(-2) * 2 * sum((c.final - data)*base))
    })

    print("grad deriv")

    time_new <- proc.time()
    print(time_new - time_prior)
    time_prior <- time_new


  dc.ai.dA <- sweep(sweep(dc.ai.dF,MARGIN=2,dF.base,"*") + dc.ai.dA_specific,
                    MARGIN=2,c.agos, "*")

  print("AGO deriv")
  time_new <- proc.time()
  print(time_new - time_prior)
  time_prior <- time_new



  colnames(c.final) <- colnames(data)
  use <- which(c.final > 0 & data > 0)

  print("use")
  time_new <- proc.time()
  print(time_new - time_prior)
  time_prior <- time_new

  sumoflogsquares <<- sum((unlist(c.final) - unlist(data))^2)

  print("logsquares")
  time_new <- proc.time()
  print(time_new - time_prior)
  time_prior <- time_new


  gradient_with_Ago <- log(10)*stock.ago*2*sum((unlist(c.final) - unlist(data))*dc.ai.dA)

  print("gradient with Ago")
  time_new <- proc.time()
  print(time_new - time_prior)
  time_prior <- time_new


  gradient_with_bg <- log(10)*B*2*sum((unlist(c.final) - unlist(data))*dc.ai.dB)

  print("gradient with bg")
  time_new <- proc.time()
  print(time_new - time_prior)
  time_prior <- time_new


  return(list(c(f_grad, f_grad_numeric, f_4, f_3, f_2, f_1, f_0, f.mat[1,1], c.final[1,1]),c(grad_derivs, gradient_with_bg, gradient_with_Ago)))
}

# pars_test <- pars
# pars_test["AGO"] <- 4
# grad_analytic <- Gradient(pars_test)
# break
# grad_numeric <- grad(ModelLikelihood,pars_test)
# plot(grad_analytic, grad_numeric, col = c(site_cols[names(pars),]))
# segments(min(grad_numeric),min(grad_numeric),max(grad_numeric), max(grad_numeric), lty = 2)
# break


# out_8mer <- c(0,0,0,0,0,0,0)
#  for (i in seq(pars["8mer"]-0.1,pars["8mer"]+0.1,length=100)) {
#   pars_temp <- pars
#   pars_temp["8mer"] <- i
#   print(pars_temp["8mer"])
#   out_8mer <- rbind(out_8mer,ModelLikelihood(pars_temp, out_8mer))
#  }

out_AGO <- c(0,0,0,0,0,0,0,0,0)
 for (i in seq(-5,6,length=50)) {
  pars_temp <- pars
  pars_temp["8mer"] <- i
  out_AGO <- rbind(out_AGO,Gradient(pars_temp)[[1]])
 }
out_AGO <-out_AGO[-1,]
plot(out_AGO[,8],out_AGO[,1],type="l",log='x',col="red")
lines(out_AGO[,8],out_AGO[,2],col="blue", lty=2)
# lines(out_AGO[,8],out_AGO[,3],col="orange")
# lines(out_AGO[,8],out_AGO[,4],col="blue")
# lines(out_AGO[,8],out_AGO[,5],col="purple")
# lines(out_AGO[,8],out_AGO[,6],col="green")
# lines(out_AGO[,8],out_AGO[,7],col="red")
# lines(out_AGO[,8],out_AGO[,9],col="red",lwd=2)

break
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

# specific_figure_string <- paste0(start, "-", stop, "_", sitelist, "_",
#   "singlebg_combinedinput_logresidual_PAPER.txt")

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

  solution <- optim(pars,
                    ModelLikelihood,
                    gr = Gradient,
                    out_model = out,
                    method = "L-BFGS-B",
                    lower = rep(-20, length(pars)),
                    upper = rep(7, length(pars)),
                    control = list(maxit=10000))


  pars <- solution$par
  sumoflogsquares_current <- solution$value
  print(sumoflogsquares_current)
#   out <- rbind(out, c(pars, sumoflogsquares_current))

#   PlotSiteKdOptimization(out, specific_figure_string, mirna, num.kds, num.bgs)
#   out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
#                      experiment, "/kds_PAPER/", start, "-", stop, "_", 
#                      sitelist, 
#                      "_singlebg_combinedinput_logresidual_PAPER.txt")
  
#   write.table(file=out_file, out, sep="\t", quote=FALSE, row.names=FALSE)
# }



