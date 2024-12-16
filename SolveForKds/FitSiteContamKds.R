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
# Initial parameters and constants.
# args       <- commandArgs(trailingOnly=TRUE)
# mirna      <- args[1]
# experiment <- args[2]
# n_constant <- args[3]
# sitelist   <- args[4]
# if (length(args) == 5) {
#   reps = as.integer(args[5])
# } else {
#   reps = 200
# }
# # Region within random sequence from which site types orginiates,
# # going from position [26 - "start" : 26 + 37 + "stop"]stop  <- as.integer(args[3])
# # Parameter specifying which list of sites is used for the analysis.
# # I.E "Current" includes all the current site types, "canonical" is just
# # the six seed-based 6mers, 7mers, and 8mers, along with no site., etc.

# Loads general functions used in AgoRBNS analysis.

# TODO make a new table of miRNA concentrations or just convert within the
# script.
stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
                       "/SolveForKds/k_c_stockago.txt"), row.names=1,
                       header=TRUE, sep="\t")
k.c.stockago <- stockago[mirna,experiment]
k.c.lib <- 100

# MAIN #########################################################################
# 1. Get data table:
sitesXcounts <- GetSitesXCounts(mirna,
                                experiment,
                                n_constant,
                                sitelist,
                                mirna.start = mirna.start,
                                mirna.stop = mirna.stop)
print(sitesXcounts)
# Separate site sequences from data file.
if (sitelist %in% kmer_list_names) {
  seqs <- rownames(sitesXcounts)
} else {
  seqs <- sitesXcounts[,1]
  sitesXcounts <- sitesXcounts[,-1]
}

# Remove any rows with no reads whoatsoever.
sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
# Remove any rows for which there are no input reads:
sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]

data <- sitesXcounts[,3:7]
data <- data.matrix(data)

data.all <- data
rownames(data)[nrow(data)] <- "None"
# data <- rbind(data_new,data[nrow(data),,drop=FALSE])

kNumSites <- nrow(data)
kNumBgs <- 1

# Define the vector of total site type concentrations:
c.I.tots <- Norm(sitesXcounts[, 2]+1) * k.c.lib
names(c.I.tots) <- rownames(data)
c.total <- matrix(c.I.tots, nrow=length(c.I.tots), ncol=5, byrow=FALSE)

colnames(c.total) <- colnames(data)
rownames(c.total) <- rownames(data)


# Define the vector of AGO concentrations:
print("hi")
# PARAMETER INITIALIZATION:
# Initialize starting Kds, which are set to 1/the enrichment of each site type
# in the A12.6 exp.
r_numerator <- Norm(rowSums(data))
r_denomenator <- Norm(c.I.tots)
kds.init <- c(r_denomenator) / c(r_numerator)
kds.init <- kds.init/max(kds.init)
pars.init <- c(log10(kds.init), 1, 1, -1, -1)


names(pars.init) <- c(rownames(data), "kd.con", "AGO", "CON", "bg")

pars <- pars.init

print("up to model")

colors <- FALSE

out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                   experiment, "/kds_PAPER/", n_constant, "_", 
                   sitelist, "_contam_equil_PAPER.txt")


ModelFunction <- function(pars, l=c.I.tots) {
  # print(c.totals[1,1])
  # print(c.I.tots[1])
  # Split up the parameters into the kd and background parameters.
  kds.a  <- 10^pars[1:num.kds]
  kd.b   <- 10^pars[num.kds + 1]
  A      <- 10^pars[num.kds + 2]
  B      <- 10^pars[num.kds + 3]
  bg     <- 10^pars[num.kds + 4]
  # Get the dilution of protein in each experiment:
  dils <- sapply(colnames(data), function(x) {
    as.numeric(x) / 100
  })
  # Calculate the bound RNA:
  c.bound <- as.matrix(sapply(dils, function(dil) {
    a.b. <- GetFreeAgoAndContam(kds.a, kd.b, l, A, B)
    a <- a.b.[1]
    b <- a.b.[2]
    a.b.occ.list <- GetOccupanciesContaminant(a, b, kds.a, kd.b)
    xa <- l * a.b.occ.list[[1]]
    xb <- l * a.b.occ.list[[2]]
    x <- xa + xb
    return(x)
  }))
  # Get the amount of background binding by subtracting the bound from the
  # total sites in each exp, normalizing. Must transpose to multiply
  # each column.
  c.free  <- c.total - c.bound
  c.bg    <- bg * t(t(c.free) / colSums(c.free))
  c.all   <- c.bound + c.bg
  c.final <- t(t(c.all) / colSums(c.all) * colSums(data))
  colnames(c.final) <- colnames(data)
  return(c.final)
}
MakeX <- function(pars, l=c.I.tots, column = FALSE) {
  # print(c.totals[1,1])
  # print(c.I.tots[1])
  # Split up the parameters into the kd and background parameters.
  kds.a  <- 10^pars[1:kNumSites]
  kd.b   <- 10^pars[kNumSites + 1]
  A      <- 10^pars[kNumSites + 2]
  B      <- 10^pars[kNumSites + 3]
  bg     <- 10^pars[kNumSites + 4]

  # Get the dilution of protein in each experiment:
  dils <- sapply(colnames(data), function(x) {
    as.numeric(x) / 100
  })
  # Calculate the bound RNA:
  x. <- as.matrix(sapply(dils, function(dil) {
    a.b. <- GetFreeAgoAndContam(kds.a, kd.b, l, A, B)
    a <- a.b.[1]
    b <- a.b.[2]
    a.b.occ.list <- GetOccupanciesContaminant(a, b, kds.a, kd.b)
    xa <- l * a.b.occ.list[[1]]
    xb <- l * a.b.occ.list[[2]]
    x <- c(xa, xb)
    return(x)
  }))
  if (column != FALSE) {
    x. <- x.[,column]
  }
  return(x.)
}

ExtractBgs <- function(pars) {
    bg     <- 10^pars[kNumSites + 4]
    return(bg)
}


MakeModelPrediction <- function(x, bg, l = c.I.tots) {
  # Combine the ago-bound and contaminant-bound site types:
  # Written out the numerator for maximum accuracy:
  # Form of equation is :        x(L - X) + B(l - x)
  #                          D * –––––––––––––––––––
  #                                (L - X)(X + B)
  print(x)
  xa    <- x[1:kNumSites, ] 
  xb    <- x[(kNumSites + 1):(2*kNumSites), ]
  x     <- xa + xb
  X     <- colSums(x)
  L     <- sum(l)
  D     <- colSums(data)
  B     <- bg
  # Transpose D-multiply x and l for row-multiplication:
  Dx.     <- D * t(x)
  print(Dx.)
  Dl.     <- D %o% l
  pred.  <- (Dx.*L - Dx.*X + Dl.*B - Dx.*B) /
                (L*B + L*X -  B*X - X^2)
  pred   <- t(pred.)
  print(pred)
  return(pred)
}


pred.1 <- ModelFunction(pars.init)

x <- MakeX(pars.init)
bg <- ExtractBgs(pars.init)
print(bg)
pred.2 <- MakeModelPrediction(x, bg)

plot(c(pred.1), c(pred.2), log = 'xy')
abline(0, 1, lty = 2)

ModelLikelihood <- function(pars){
  model <- ModelFunction(pars)
  model_norm <- t(t(model) / colSums(model))
  loglikelihood <- -sum(data*log(model_norm))
  return(loglikelihood)
}


DXDParams <-function(pars, dil, l = c.I.tots) {
  # The site type on and off rates, the background terms for each column,
  # the and the total concentration of Ago (a) and the contaminant (b):
  kds.a  <- 10^pars[1:kNumSites]
  kd.b   <- 10^pars[kNumSites + 1]
  A      <- 10^pars[kNumSites + 2] * dil
  B      <- 10^pars[kNumSites + 3] * dil

  # Assignment of the Kds (being the ratio of the on and off rates)
  if (print.diag == TRUE) print("solving for free ago and contam:")
  a.b. <- GetFreeAgoAndContam(kds.a, kd.b, l, A, B)
  a <- a.b.[1]
  b <- a.b.[2]
  a.b.occ.list <- GetOccupanciesContaminant(a, b, kds.a, kd.b)
  xa <- l * a.b.occ.list[[1]]
  xb <- l * a.b.occ.list[[2]]
  x <- xa + xb


  s0 <- (l * kds.a * kd.b)/(a*kd.b + b*kds.a + kds.a*kd.b)^2
  sx <- (b + kd.b )*(l * kds.a * kd.b)/(a*kd.b + b*kds.a + kds.a*kd.b)^2
  sy <- (a + kds.a)*(l * kds.a * kd.b)/(a*kd.b + b*kds.a + kds.a*kd.b)^2

  S0 <- sum(s0)
  Sx <- sum(sx)
  Sy <- sum(sy)
  # Define constant denomenator term:
  C.den <- (1 + Sy)*(1 + Sx) - a*b*S0*S0

  d.x.d.A <- a *  (c(   sx, -b*s0, 1, 0)*(1 + Sy) +
                 c(-a*s0,    sy, 0, 1)*(b * S0)   )/C.den
  d.x.d.B <- b *  (c(   sx, -b*s0, 1, 0)*(a * S0) +
                 c(-a*s0,    sy, 0, 1)*(1 + Sx)   )/C.den
  print(d.x.d.A)
  print(d.x.d.B)
  d.x.d.Kb <- kd.b * b/kd.b*((c(   sx, -b*s0, 1, 0)*(-a*S0                 ) +
                          c(-a*s0,    sy, 0, 1)*(Sy + Sx*Sy - a*b*S0*S0)  )/C.den + 
                          c( a*s0,   -sy, 0, 0)            )
  d.x.d.Ks <-  t(kds.a*t((c(  sx,  -b*s0, 1, 0)%o%(a/kds.a*s0*((kd.b + b)*(1 + Sy) - a*b*S0)/C.den) +
                  c(-a*s0,    sy, 0, 1)%o%(-a*b/kds.a*s0/C.den)                                   ) +
                  a*rbind(diag(-sx/kds.a), diag(b*s0/kds.a))))[, 1:kNumSites]
  # Final calculation of the d.x0.d.pars matrix:

  d.x.d.pars <- log(10) * cbind(d.x.d.Ks, d.x0.d.Kb, d.x0.d.A, d.x0.d.B)
  return(d.x0.d.pars)
}

d.x1.d.pars.n <- jacobian(MakeX, pars.init, column=1)
d.x1.d.pars.a <- DXDParams(pars.init,as.numeric(colnames(x)[1])/100)

plot(c(d.x1.d.pars.n), c(d.x1.d.pars.a), log = 'xy')
abline(0, 1, lty = 2)
break

DPredDXColumn <- function(x, bgs, col) {
  # Combine the ago-bound and contaminant-bound site types:
  R <- kDataTotals[col]
  l <- kInputMatrix[, col]
  L <- sum(l)
  B <- bgs[col]
  x <- x[1:(2 * kNumSites), col] + x[(2 * kNumSites + 1):(4 * kNumSites), col]
  X <- sum(x)
  Cl <- L - X
  Cb <- X + B

  dpidxi.submat <- diag(rep((Cl - B) / (Cl * Cb), length = 2*kNumSites))
  dpidxi.mat <- cbind(dpidxi.submat, dpidxi.submat, 0, 0)
  dpxidj.vector <- (B * (Cb - Cl) * (l - x) - Cl^2 * x) / (Cl * Cb)^2
  dpidxj.submat <- matrix(dpxidj.vector, nrow = 2 * kNumSites,
                          ncol = 4 * kNumSites)
  dpidxj.mat    <- cbind(dpidxj.submat, 0, 0)
  dpdx.t <- R * (dpidxi.mat + dpidxj.mat)
  return(dpdx.t)
}

DCostDPredMultinom <- function(pred) {
  probs <- (pred + 1) %*% (1 / colSums(pred + 1))
  d.cost.d.pred <-  (1 - (kData + 1) / (pred + 1))
  d.cost.d.pred[(nrow(d.cost.d.pred)/2+1):nrow(d.cost.d.pred)] <- d.cost.d.pred[(nrow(d.cost.d.pred)/2+1):nrow(d.cost.d.pred)] * 10

  return(d.cost.d.pred/1000)
}


Gradient <- function(pars) {
  # Split up the parameters into the kd and background parameters.
  kds  <- 10^pars[1 : kNumSites]
  B  <- 10^pars[(kNumSites + 1)]
  stock.ago <- 10^pars[kNumSites + 2]
  c.agos <- sapply(colnames(data), function(x) {
    as.numeric(x) / 100
  })

   f.mat <- matrix(
      rep(sapply(c.agos, function(percent) {
           return(GetFreeAgo(kds, c.I.tots, percent * stock.ago))
         }
         ), nrow(data)), nrow=nrow(data), ncol=ncol(data), byrow=TRUE)

  f.jvec <- sapply(c.agos, function(ago.percent) {
           return(GetFreeAgo(kds, c.I.tots, ago.percent * stock.ago))})

  A.mat <- matrix(c.agos*stock.ago, nrow=nrow(data), ncol=ncol(data),
                  byrow=TRUE)

  A.jvec <- c.agos * stock.ago

  K.mat <- matrix(kds, nrow=nrow(data), ncol=ncol(data), byrow=FALSE)
  K.ivec <- kds

  l.mat <- matrix(c.I.tots, nrow=nrow(data), ncol=ncol(data), byrow=FALSE)
  l.ivec <- c.I.tots

  R.mat <- matrix(colSums(data), nrow=nrow(data), ncol=ncol(data), byrow=TRUE)
  R.jvec <- colSums(data) 
  L <- 100

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
  names(gradient_all) <- names(pars)
  return(gradient_all)
}






print(data)
solution <- optim(pars,
                  ModelLikelihood,
                  gr = Gradient,
                  method = "L-BFGS-B",
                  lower=c(rep(-13, length=kNumSites), -2, -1),
                  upper=c(rep(10, length=kNumSites), 1, 2),
                  control = list(maxit=100000, factr=1e2))

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
  input_resample <- rmultinom(1,size=sum(sitesXcounts[,2]),prob=sitesXcounts[,2])
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
                    lower=c(rep(-13, length=kNumSites), -2, -1),
                    upper=c(rep(10, length=kNumSites), 1, 2),
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





