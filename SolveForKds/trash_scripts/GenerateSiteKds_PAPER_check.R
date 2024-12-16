################################################################################
#GenerateSiteKds_PAPER.py
################################################################################

## THIS IS THE OFFICIAL SCRIPT FOR FINALIZING THE MANUSCRIPT.
## THIS SCRIPT HAS THESE DECISIONS MADE:
## 1.) PROTEIN IS FIT AS A PARAMETER
## 2.) A SINGLE BACKGROUND PARAMETER IS GIVEN TO ALL FIVE CONCENTRATIONS.
## 3.) LIST, START, AND STOP POSITION ARE MANDATORY PARAMETERS.
# ## 4.) ALL INPUTS ARE COMBINED ACROSS FOR THESE ANALYSES.
library(colorspace)
library(multicore)
library(data.table)
source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")
source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/ModelingFunctions.R")
source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/PlotFunctions.R")
# Initial parameters and constants.
# args       <- commandArgs(trailingOnly=TRUE)
# mirna      <- args[1]
# exp        <- args[2]
# n_constant <- args[3]
# sitelist   <- args[4]
# if (length(args) == 5) {
#   reps = as.integer(args[5])
# } else {
#   reps = 200
# }
# Loads general functions used in AgoRBNS analysis.

mirna <- "miR-1"
experiment <- "equilibrium"
n_constant <- 5
sitelist <- "canonical"
k.c.stockago <- stockago[mirna, exper]
k.c.lib <- 100
#23456789!123456789@123456789#123456789$123456789
# MAIN #########################################################################
# 1. Get data table:
sXc <- SitesXCounts(mirna, experiment, n_constant, sitelist)
# Separate site sequences from data file.
if (sitelist %in% kKmerSiteLists) {
  seqs <- rownames(sXc)
} else {
  seqs <- sXc[,1]
  sXc <- sXc[,-1]
}
data <- sXc[, 3:7]
data <- data.matrix(data)
data.all <- data

# Define all constants for script:
kAgoStock <- stockago[mirna, experiment]
kInput <- Norm(sXc[, 2])*100
kInputMatrix <- matrix(kInput, nrow=length(kInput), ncol=5, byrow=FALSE)
colnames(kInputMatrix) <- colnames(data)
rownames(kInputMatrix) <- rownames(data)
kAgoDils <- sapply(colnames(data), as.numeric)/100 # Ago concentration
kOutputFile <- paste0(kSolexaDir, mirna, "/", experiment, "/kds_PAPER/", n_constant, "_",
                   sitelist, "_PAPER.txt")

InitializePars <- function(freqs=data, l=kInput) {
  # Function to initialize the parameters:
  kds <- log10(1/rowMeans(EquilEnrichments(freqs, l)))
  names(kds) <- rownames(data)
  pars <- c(kds, bg=-1, AGO=log10(kAgoStock))
  return(pars)
}
pars <- InitializePars()


model.p <- EquilSingleSiteModel(pars, sXc)
print(MultinomialCost(model.p, data.all))
print(MultinomialCost2(model.p, data.all))

break
f.x_M <- function(p) {
  # Split up the parameters into the kd and background parameters.
  P   <- 10^p
  kds_c <- P[kIndsKds]
  names(kds_c) <- sitenames
  l_c   <- P[kIndsFrI]*kL
  A.  <- P[kIndsAgo]*kADil.
  x_M <- sapply(A., GetBoundRNA, l=l_c, kds=kds_c)
  return(x_M)  
}

f.g_M <- function(x_M, p, d=kData) {
  # Split up the parameters into the kd and background parameters.
  P   <- 10^p
  G   <- P[kIndsBgs]
  l_M <- P[kIndsFrI]%o%kL.
  rownames(l_M) <- sitenames
  X_c <- colSums(x_M)
  L_c <- colSums(l_M)
  g_M <- t(t(l_M - x_M)/(L_c - X_c))
  return(g_M) 
}



x_M <- f.x_M(pars)
g_M <- f.g_M(x_M, pars)

break

MakeModelPrediction <- function(x, bgs, data = kData, l = kInputMatrix) {
  # Combine the ago-bound and contaminant-bound site types:
  # Written out the numerator for maximum accuracy:
  # Form of equation is :        x(L - X) + B(l - x)
  #                          D * –––––––––––––––––––
  #                                (L - X)(X + B)
  x.ago <- x[1:(2 * kNumSites), ] 
  x.con <- x[(2 * kNumSites + 1):(4 * kNumSites), ]
  x     <- x.ago + x.con
  X     <- colSums(x)
  L     <- colSums(l)
  D     <- colSums(kData)
  B     <- bgs
  # Transpose D-multiply x and l for row-multiplication:
  Dx.     <- D * t(x)
  Dl.     <- D * t(l)
  pred.  <- (Dx.*L - Dx.*X + Dl.*B - Dx.*B) /
                (L*B + L*X -  B*X - X^2)
  pred   <- t(pred.)
  return(pred)
}


ModelFunction <- function(p) {
  # print(c.totals[1,1])
  # print(c.I.tots[1])
  # Split up the parameters into the kd and background parameters.
  p. <- 10^pars
  kds     <- 10^p.[kIndsKds]
  l       <- 10^p.[kIndsIFrI]
  G       <- 10^p.[kIndsBgs]
  A.      <- 10^p.[kIndsAgo]

  # Get the bound counts in each exp:
  A. <- sapply(colnames(data), function(x) {
    as.numeric(x) / 100
  })
  x <- as.matrix(
    sapply(A., function(percent) {
      return(GetBoundRNA(kds, frI*100, percent * stock.ago))
    }
  ))
  # Get the amount of background binding by subtracting the bound from the
  # total sites in each exp, normalizing. Must transpose to multiply
  # each column.
  c.frees <- frI%o%columnI - c.bounds
  c.bgs <- t(t(c.frees) * bgs / colSums(c.frees))
  c.all <- c.bounds + c.bgs
  c.final <- data.frame(t(t(c.all) / colSums(c.all) * colSums(kData)))
  return(c.final)
}



ModelLikelihood <- function(p){
  model <- ModelFunction(p)
  model_norm <- t(t(model) / colSums(model))
  loglikelihood <- -sum(data*log(model_norm))
  return(loglikelihood)
}

CostLogRes <- function(pred, column = FALSE) {
  if (column != FALSE) {
    pred <- pred[, column, drop=FALSE]
    kData <- kData[, column, drop=FALSE]
  }
  cost <- log((pred + 1) / (kData + 1))^2
  cost.total <- sum(cost)
  return(cost.total)
}

CostMultinom <- function(pred, column = FALSE) {
  # In this function the model is of the read
  probs <- t(t(pred + 1) / colSums(pred + 1))
  if (column != FALSE) {
    probs <- probs[, column, drop=FALSE]
    kData <- kData[, column, drop=FALSE]
  }
  # print(dim(probs))
  # print(dim(kData))
  bestprobs <- t(t(kData + 1) / colSums(kData + 1))
  cost <- -(kData + 1) * log(probs) + (kData + 1) * log(bestprobs)
  # cost[(nrow(cost)/2+1):nrow(cost)] <- cost[(nrow(cost)/2+1):nrow(cost)] * 10
  cost.total <- sum(cost)
  return(cost.total)
}

DCostDPredLogRes <- function(pred) {
  d.cost.d.pred <- 2 * log((pred + 1) / (kData + 1)) / (pred + 1)
  return(d.cost.d.pred)
}

DCostDPredMultinom <- function(pred) {
  probs <- (pred + 1) %*% (1 / colSums(pred + 1))
  d.cost.d.pred <-  (1 - (kData + 1) / (pred + 1))
  # d.cost.d.pred[(nrow(d.cost.d.pred)/2+1):nrow(d.cost.d.pred)] <- d.cost.d.pred[(nrow(d.cost.d.pred)/2+1):nrow(d.cost.d.pred)] * 10

  return(d.cost.d.pred/1000)
}


dGdp <- function(p) {
  # Split up the parameters into the kd and background parameters.
  p. <- 10^p
  kds       <- 10^p.[kIndsKds]
  frI       <- 10^p.[kIndsFrI]
  bg        <- 10^p.[kIndsBgs]
  stock.ago <- 10^p.[kIndsAgo]

  c.agos <- sapply(colnames(data), function(x) {
    as.numeric(x) / 100
  })

   f.mat <- matrix(
      rep(sapply(c.agos, function(percent) {
           return(GetFreeAgo(kds, frI*100, percent * stock.ago))
         }
         ), nrow(data)), nrow=nrow(data), ncol=ncol(data), byrow=TRUE)

  f.jvec <- sapply(c.agos, function(ago.percent) {
           return(GetFreeAgo(kds, frI*100, ago.percent * stock.ago))})

  A.mat <- matrix(c.agos*stock.ago, nrow=nrow(data), ncol=ncol(data),
                  byrow=TRUE)

  A.jvec <- c.agos * stock.ago

  K.mat <- matrix(kds, nrow=nrow(data), ncol=ncol(data), byrow=FALSE)
  K.ivec <- kds

  l.mat <- matrix(frI*100, nrow=nrow(data), ncol=ncol(data), byrow=FALSE)
  l.ivec <- frI*100

  R.mat <- matrix(colSums(data), nrow=nrow(data), ncol=ncol(data), byrow=TRUE)
  R.jvec <- colSums(data) 
  L <- 100
  B <- bg[1]
  # Get the amount of background binding by subtracting the bound from the
  # total sites in each exp, normalizing. Must transpose to multiply
  # each column.

  time.init <- proc.time()
  c.vec_num <- -(
                 (l.ivec %*% t(R.jvec * f.jvec * (f.jvec + L - A.jvec)) +
                 (l.ivec * K.ivec * B)  %*% t(R.jvec))
                ) 
  C1.jvec <- L - B - 2 * A.jvec
  C2.jvec <- A.jvec^2 + (B - L) * A.jvec - L * B


  c.vec_dem <- t(
                (f.jvec^3 + f.jvec^2 * C1.jvec + f.jvec * C2.jvec) +
                t(K.ivec %*% t(f.jvec^2 + f.jvec * C1.jvec + C2.jvec))
                )^-1
  c.final <- c.vec_num * c.vec_dem

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

  residuals <- -data/c.final

  grad_derivs <- (log(10)*kds *
                  (colSums(residuals*dc.ai.dF) %*%
                   t(dF.dK.mat
                  ) + rowSums(dc.ai.dKj_specific*residuals)))
  
  dc.ai.dA <- sweep(sweep(dc.ai.dF,MARGIN=2,dF.base,"*") + dc.ai.dA_specific,
                    MARGIN=2,c.agos, "*")

  colnames(c.final) <- colnames(data)

  gradient_with_Ago <- log(10)*stock.ago*sum(residuals * dc.ai.dA)

  gradient_with_bg <- log(10)*B*sum(residuals * dc.ai.dB)
  gradient_all <- c(grad_derivs, gradient_with_bg, gradient_with_Ago)
  names(gradient_all) <- names(pars[c(kIndsKds, kIndsBgs, kIndsAgo)])
  return(gradient_all)
}

break

NewCostFunction <- function(pars) {

}


print(data)
solution <- optim(pars,
                  ModelLikelihood,
                  gr = Gradient,
                  method = "L-BFGS-B",
                  lower=c(rep(-13, length=num.kds), -2, -1),
                  upper=c(rep(10, length=num.kds), 1, 2),
                  control = list(maxit=100000, factr=1e2))



break
pars <- solution$par

model <- ModelFunction(pars)

plot(colnames(data),rep(1,5),type="l",lty=2,xlim=c(0.1,100),ylim=c(0.3,1000),log='xy')
sapply(seq(nrow(model)), function(row){
  print(row)
  points(colnames(data),data[row,]/colSums(data)*sum(c.I.tots)/c.I.tots[row],col=kSiteColors[rownames(data)[row],])
  lines(colnames(data),model[row,]/colSums(model)*sum(c.I.tots)/c.I.tots[row],col=kSiteColors[rownames(data)[row],])
  })
print(out_file)
print(pars)
pars.save <- pars
print(reps)

pars_loocv <- matrix(NaN,nrow=length(pars),ncol=reps)
rownames(pars_loocv) <- names(pars)
colnames(pars_loocv) <- seq(reps)
for (i_trial in seq(reps)) {
  print(i_trial)
  tick <- 0
  i_col <- sample(1:ncol(data.all),1)
  data.temp <- data.all[,-i_col]
  data <- apply(data.temp, 2, function(col) {rmultinom(1,size=sum(col), prob=col)})
  rownames(data) <- rownames(data.temp)
  colnames(data) <- colnames(data.temp)
  input_resample <- rmultinom(1,size=sum(sXc[,2]),prob=sXc[,2])
  c.I.tots <- Norm(input_resample+1)*k.c.lib
  c.totals <- matrix(
                rep(c.I.tots,ncol(data)), nrow=length(c.I.tots), ncol=ncol(data), byrow=FALSE)
  colnames(c.totals) <- colnames(data)
  rownames(c.totals) <- rownames(data)
  pars <- pars.save
  solution <- optim(pars,
                    ModelLikelihood,
                    gr = Gradient,
                    method = "L-BFGS-B",
                    lower=c(rep(-13, length=num.kds), -2, -1),
                    upper=c(rep(10, length=num.kds), 1, 2),
                    control = list(maxit=100000, factr=1e2))

  pars <- solution$par

  pars_loocv[,i_trial] <- pars

}

pars_loocv_sort <- t(apply(pars_loocv,1,sort))
print(pars_loocv_sort)
print(pars.save)
print(10^pars.save)
print(10^cbind(pars.save,rowMeans(pars_loocv)))
output <- 10^cbind(pars.save,rowMeans(pars_loocv),
                   pars_loocv_sort[,ceiling(0.025*reps)],
                   pars_loocv_sort[,ceiling(0.5*reps)],
                   pars_loocv_sort[,ceiling(0.975*reps)])
colnames(output) <- c("Full","Mean", "Lower_CI", "Median", "Upper_CI")
write.table(file=out_file, output, sep="\t", quote = FALSE)
print(output)
print(out_file)


warnings()





