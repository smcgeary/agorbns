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
# Initial parameters and constants.
args       <- commandArgs(trailingOnly=TRUE)
mirna      <- args[1]
experiment        <- args[2]
n_constant <- args[3]
sitelist   <- args[4]
if (length(args) == 5) {
  reps = as.integer(args[5])
} else {
  reps = 200
}
# mirna <- "miR-1"
# experiment <- "equilibrium"
# n_constant <- 5
# sitelist <- "bulge"
# MAIN #########################################################################
# 1. Get data table:
sXc <- SitesXCounts(mirna, experiment, n_constant, sitelist, pyclasses=TRUE)
# Separate site sequences from data file.
data <- sXc[, 3:7]
data <- data.matrix(data)
data.saved <- data
# Define all constants for script:
kAgoStock <- stockago[mirna, experiment]
kInput <- Norm(sXc[, 2])*100
kAgoDils <- sapply(colnames(data), as.numeric)/100 # Ago concentration
kOutputFile <- paste0(kSolexaDir, mirna, "/", experiment, "/kds_PAPER/", n_constant, "_",
                   sitelist, "_classes_PAPER.txt")

kIndsKds <- seq(nrow(data)-2)
InitializePars <- function(data) {
  # Function to initialize the parameters:
  kds <- log10(Norm(kInput)/Norm(rowSums(data)))
  names(kds) <- rownames(data)
  pars <- c(kds, bg=-1, AGO=log10(kAgoStock))
  return(pars)
}
# PARAMETER INITIALIZATION:
# Initialize starting Kds, which are set to 1/the enrichment of each site type
# in the A12.6 experiment.
r_numerator <- Norm(rowSums(data))
r_denomenator <- Norm(kInput)
kds.init <- c(r_denomenator) / c(r_numerator)
kds.init <- kds.init/max(kds.init)
pars.init <- c(log10(kds.init), -1, 0)

names(pars.init) <- c(rownames(data), "bg", "AGO")

pars <- pars.init
tick <- 1

print("up to model")


kOutputFile <- paste0(kSolexaDir, mirna, "/", experiment, "/kds_PAPER/",
                   n_constant, "_", sitelist, "_classes_PAPER.txt")



Enrichent <- function(matrix) {
  sapply(matrix, 2, function(column) {
    Norm(column)/kInput/100
  })
}

ModelFunction <- function(pars, data) {
  # Split up the parameters into the kd and background parameters.
  kds  <- 10^pars[1 : nrow(data)]
  names(kds) <- rownames(data)
  bgs  <- 10^pars[nrow(data) + 1]
  A.stock <- 10^pars[nrow(data) + 1 + 1]
  # print(A.stock*c.agos)
  # Get the bound counts in each experiment:
  c.agos <- sapply(colnames(data), function(x) {
  as.numeric(x)/100
  })

  c.bounds <- as.matrix(
    sapply(c.agos, function(ago.percent) {
      return(BoundRNA(kds, kInput, ago.percent * A.stock))
    }
  ))
  # Get the amount of background binding by subtracting the bound from the
  # total sites in each experiment, normalizing. Must transpose to multiply
  # each column.
  c.frees <- apply(c.bounds, 2, function(c.col) kInput - c.col)
  c.bgs <- t(t(c.frees) * bgs / colSums(c.frees))
  c.all <- c.bounds + c.bgs
  c.final <- data.frame(t(t(c.all) / colSums(c.all) * colSums(data)))
  c.final_new <<- c.final
  return(c.final)
}


Enrichment <- function(matrix) {
  apply(matrix, 2, function(column) {
    Norm(column)/kInput*100
  })
}

plotmodel <- function(model, data) {
  xrow <- c(40, 12.6, 4, 1.26, 0.4)
  Rmodel <- Enrichment(model)
  Rdata <- Enrichment(data)
  for (row in seq(nrow(Rmodel))) {
    lines(xrow, Rmodel[row,])
    points(xrow, Rdata[row,])
  }
}



ModelLikelihood <- function(pars, data){

  model <- ModelFunction(pars, data)
    lc.final <- log(model+1)
  ldata <- log(data+1)
  model_norm <- t(t(model) / colSums(model))
  # sumofsquares <- sum((lc.final - ldata)^2)
  sumofsquares <- -sum(data*log(model_norm))
  # if (tick %% 1000 == 0 & tick > 1) {
  #   dev.set(2)
  #   plot(c(0.0001,0.0002),
  #        c(0.0001,0.0002),
  #        col="white",
  #        xlim=c(0.1,100),
  #        ylim=c(0.5,500),
  #        log='xy')
  #   title(main=mirna, font.main = 1)
  #   plotmodel(model, data)
  #   print(sumofsquares)
  #   print(sort(pars[1:nrow(data)]))
  # }
  # tick <<- tick + 1
  return(sumofsquares)
}

Gradient <- function(pars, data) {
  # Split up the parameters into the kd and background parameters.
  kds  <- 10^pars[1 : nrow(data)]
  names(kds) <- rownames(data)
  c.agos <- sapply(colnames(data), function(x) {
  as.numeric(x)/100
})

  B  <- 10^pars[nrow(data) + 1]
  A.stock <- 10^pars[nrow(data) + 2]
   f.mat <- matrix(
      rep(sapply(c.agos, function(ago.percent) {
           return(FreeAgo(kds, kInput, ago.percent * A.stock))
         }
         ), nrow(data)), nrow=nrow(data), ncol=ncol(data), byrow=TRUE)

  f.jvec <- sapply(c.agos, function(ago.percent) {
           return(FreeAgo(kds, kInput, ago.percent * A.stock))})

  A.mat <- matrix(c.agos*A.stock, nrow=nrow(data), ncol=ncol(data),
                  byrow=TRUE)

  A.jvec <- c.agos * A.stock

  K.mat <- matrix(kds, nrow=nrow(data), ncol=ncol(data), byrow=FALSE)
  K.ivec <- kds

  l.mat <- matrix(kInput, nrow=nrow(data), ncol=ncol(data), byrow=FALSE)
  l.ivec <- kInput

    B.mat  <- matrix(B, nrow=nrow(data), ncol =ncol(data), byrow=FALSE)


  R.mat <- matrix(colSums(data), nrow=nrow(data), ncol=ncol(data), byrow=TRUE)
  R.jvec <- colSums(data) 
  L <- 100

  # Get the amount of background binding by subtracting the bound from the
  # total sites in each experiment, normalizing. Must transpose to multiply
  # each column.

  time.init <- proc.time()
  c.vec_num <- -(
                 (l.ivec %*% t(R.jvec * f.jvec * (f.jvec + L - A.jvec)) +
                 (l.ivec * K.ivec)  %*% t(B* R.jvec))
                ) 
  C1.jvec <- L - B - 2 * A.jvec
  C2.jvec <- A.jvec^2 + (B - L) * A.jvec - L * B


  c.vec_dem <- t(
                (f.jvec^3 + f.jvec^2 * C1.jvec + f.jvec * C2.jvec) +
                t(K.ivec %*% t(f.jvec^2 + f.jvec * C1.jvec + C2.jvec))
                )^-1
  c.final <- c.vec_num * c.vec_dem

  time.init <- proc.time()
  dF.base <- (1 + colSums(l.mat * K.mat * (f.mat + K.mat)^(-2)))^(-1)
  time.new <- proc.time()

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

  time.init <- proc.time()


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

  gradient_with_Ago <- log(10)*A.stock*sum(residuals * dc.ai.dA)

  gradient_with_bg <- log(10)*B*sum(residuals * dc.ai.dB)
  gradient_all <- c(grad_derivs, gradient_with_bg, gradient_with_Ago)
  names(gradient_all) <- names(pars)
  gradient_all
}


print(data)
solution <- optim(pars,
                  ModelLikelihood,
                  gr = Gradient,
                  data= data,
                  method = "L-BFGS-B",
                  lower=rep(-3, length(pars)),
                  upper=rep(3, length(pars)),
                  control = list(maxit=100000, factr=1e2))
pars.save <- solution$par

print(10^solution$par)
print("done first one")
pars_loocv <- matrix(NaN,nrow=length(pars),ncol=reps)
rownames(pars_loocv) <- names(pars)
colnames(pars_loocv) <- seq(reps)

for (i_trial in seq(reps)) {
  print(i_trial)
  tick <- 0
  i_col <- sample(1:ncol(data.saved), 1)
  data.temp <- data.saved[, -i_col]
  data <- apply(data.temp, 2, function(col) {rmultinom(1,size=sum(col), prob=col)})
  rownames(data) <- rownames(data.temp)
  colnames(data) <- colnames(data.temp)
  input_resample <- rmultinom(1, size=sum(sXc[, 2]), prob=sXc[, 2])
  kInput <- Norm(input_resample+1)*100
  names(kInput) <- rownames(data)
  pars <- pars.init
  solution <- optim(pars,
                    ModelLikelihood,
                    data = data,
                    gr = Gradient,
                    method = "L-BFGS-B",
                    lower=c(rep(-13, length=nrow(data)), -2, -1),
                    upper=c(rep(10, length=nrow(data)), 1, 2),
                    control = list(maxit=100000, factr=1e2))

  pars <- solution$par
  print(10^solution$par)

  pars_loocv[,i_trial] <- pars
  
  if (i_trial%%10 == 0) {
    pars_loocv_sort <- t(apply(pars_loocv[!(is.na(pars_loocv[, 1])), ], 1, sort))
    output <- 10^cbind(pars.save,rowMeans(pars_loocv, na.rm=TRUE),
                       pars_loocv_sort[,ceiling(0.025*i_trial)],
                       pars_loocv_sort[,ceiling(0.5*i_trial)],
                       pars_loocv_sort[,ceiling(0.975*i_trial)])
    colnames(output) <- c("Full","Mean", "Lower_CI", "Median", "Upper_CI")
    write.table(file=kOutputFile, output, sep="\t", quote = FALSE)
  }
}

warnings()