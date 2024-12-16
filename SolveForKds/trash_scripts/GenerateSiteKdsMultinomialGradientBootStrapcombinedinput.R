################################################################################
#GenerateSiteKds_PAPER.py
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
library('Matrix')

source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")
# Initial parameters and constants.
args  <- commandArgs(trailingOnly=TRUE)
mirna <- args[1]

# Region within random sequence from which site types orginiates,
# going from position [26 - "start" : 26 + 37 + "stop"]
start <- as.integer(args[2])
stop  <- as.integer(args[3])
# Parameter specifying which list of sites is used for the analysis.
# I.E "Current" includes all the current site types, "canonical" is just
# the six seed-based 6mers, 7mers, and 8mers, along with no site., etc.
sitelist <- args[4]
if (sitelist %in% kmer_list_names) {
  mirna.start <- as.integer(args[5])
  mirna.stop <- as.integer(args[6])
} else {
  mirna.start <- NULL
  mirna.stop <- NULL
}

# # Experiment name
experiment <- "equilibrium"
# print(experiment)

# Loads general functions used in AgoRBNS analysis.

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
k.c.stockago <- stockago[mirna,experiment]
k.c.lib <- 100

# MAIN #########################################################################
# 1. Get data table:
sitesXcounts <- GetSitesXCountsCombined(mirna,
                                experiment,
                                start,
                                stop,
                                sitelist,
                                mirna.start = mirna.start,
                                mirna.stop = mirna.stop)
# Separate site sequences from data file.
if (sitelist %in% kmer_list_names) {
  seqs <- rownames(sitesXcounts)
} else {
  seqs <- sitesXcounts[,1]
  sitesXcounts <- sitesXcounts[,-1]
}
sitesXcounts_original <- sitesXcounts
if (sitelist %in% kmer_list_names) {
  n <- 256
  colors <- c(rainbow_hcl(n, c=50, l=70, start=0, end=360 * (n - 1) / n),"red",
              "black")

  out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                     experiment, "/kds_PAPER/", start, "-", stop, "_", 
                     sitelist, "_", mirna.start, "-", mirna.stop,
                     "_singlebg_combinedinput_multinomial_bootstraps_PAPER.txt")

    out_ci90_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                     experiment, "/kds_PAPER/", start, "-", stop, "_", 
                     sitelist, "_", mirna.start, "-", mirna.stop,
                     "_singlebg_combinedinput_multinomial_ci90_PAPER.txt")

    out_ci95_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                     experiment, "/kds_PAPER/", start, "-", stop, "_", 
                     sitelist, "_", mirna.start, "-", mirna.stop,
                     "_singlebg_combinedinput_multinomial_ci95_PAPER.txt")



} else {
  colors <- FALSE

  out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                     experiment, "/kds_PAPER/", start, "-", stop, "_", 
                     sitelist, "_singlebg_combinedinput_multinomial_bootstraps_PAPER.txt")

  out_ci90_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                     experiment, "/kds_PAPER/", start, "-", stop, "_", 
                     sitelist, "_singlebg_combinedinput_multinomial_ci90_PAPER.txt")

  out_ci95_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                     experiment, "/kds_PAPER/", start, "-", stop, "_", 
                     sitelist, "_singlebg_combinedinput_multinomial_ci95_PAPER.txt")


}
startingkds <- GetKds(mirna, experiment, start, stop, sitelist,scaled=FALSE)


out_boots <- mclapply(1:500, function(index) {
print(index)

ModelFunction <- function(pars) {
  # Split up the parameters into the kd and background parameters.
  kds  <- c(Logistic(pars[1 : num.kds], 10),1)
  names(kds) <- rownames(data)
  bgs  <- rep(10^pars[(num.kds + 1)], 5)
  stock.ago <- 10^pars[num.kds + 2]
  # Get the bound counts in each experiment:
  c.bounds <- as.matrix(
    sapply(c.agos, function(ago.percent) {
      return(GetBoundRNA(kds, c.I.tots, ago.percent * stock.ago))
    }
  ))

  # Get the amount of background binding by subtracting the bound from the
  # total sites in each experiment, normalizing. Must transpose to multiply
  # each column.
  c.frees <- c.totals - c.bounds
  c.bgs <- t(t(c.frees) * bgs / colSums(c.frees))
  c.all <- c.bounds + c.bgs
  c.final <- data.frame(t(t(c.all) / colSums(c.all) * colSums(data)))
  c.final_new <<- c.final
  return(c.final)
}
tick <- 0
ModelLikelihood <- function(pars){

  model <- ModelFunction(pars)
  model_norm <- t(t(model) / colSums(model))
  sumofsquares <- -sum(data*log(model_norm))
  if (tick%%1000000 == 0) {
    print(sumofsquares)
  }
  tick <<- tick + 1
  return(sumofsquares)
}

Gradient <- function(pars) {
  # Split up the parameters into the kd and background parameters.
  kds  <- c(Logistic(pars[1 : num.kds], 10),1)
  names(kds) <- rownames(data)
  B  <- 10^pars[(num.kds + 1)]
  stock.ago <- 10^pars[num.kds + 2]

   f.mat <- matrix(
      rep(sapply(c.agos, function(ago.percent) {
           return(GetFreeAgo(kds, c.I.tots, ago.percent * stock.ago))
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
  # total sites in each experiment, normalizing. Must transpose to multiply
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

  grad_derivs <- (10*exp(pars[1:num.kds]) * (exp(pars[1:num.kds]) + 1)^(-2) *
                  (colSums(residuals*dc.ai.dF) %*%
                   t(dF.dK.mat[-nrow(dF.dK.mat),]
                  ) + rowSums(dc.ai.dKj_specific*residuals)[-nrow(dF.dK.mat)]))
  
  dc.ai.dA <- sweep(sweep(dc.ai.dF,MARGIN=2,dF.base,"*") + dc.ai.dA_specific,
                    MARGIN=2,c.agos, "*")

  colnames(c.final) <- colnames(data)

  gradient_with_Ago <- log(10)*stock.ago*sum(residuals * dc.ai.dA)
  gradient_with_bg <- log(10)*B*sum(residuals * dc.ai.dB)
  gradient_all <- c(grad_derivs, gradient_with_bg, gradient_with_Ago)
  names(gradient_all) <- names(pars)
  return(gradient_all)
}
  tick <- 1

  sitesXcounts_original <- sitesXcounts_original[rowSums(sitesXcounts_original)>0,]
  # Remove any rows for which there are no input reads:
  sitesXcounts_original <- sitesXcounts_original[sitesXcounts_original[,1]>0,]

  sitesXcounts <- apply(sitesXcounts_original, 2, function(col) {rmultinom(1,size=sum(col), prob=col)})
  rownames(sitesXcounts) <- rownames(sitesXcounts_original)

  # Remove any rows with no reads whoatsoever.
  data <- sitesXcounts[,2:6]
  data <- data.matrix(data)
  pars <- startingkds
  rownames(data)[nrow(data)] <- "None"
  data_inds <- which(rownames(data) %in% c(names(startingkds),"None"))
  data <- data[data_inds,]
  # names(pars) <- c(rownames(data)[-nrow(data)], "bg", "AGO")
  num.kds <- nrow(data)-1
  num.bgs <- 1

  # Define the vector of total site type concentrations:
  c.I.tots <- Norm(sitesXcounts[, 1]+1) * k.c.lib
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
  solution <- optim(pars,
                    ModelLikelihood,
                    gr = Gradient,
                    method = "L-BFGS-B",
                    lower=c(rep(-13, length=num.kds), -2, -1),
                    upper=c(rep(10, length=num.kds), 1, 2),
                    control = list(maxit=100000, factr=1e2))

  pars <- solution$par
  sumofsquares <- solution$value
  output_vector <- c(pars, sumofsquares)
  names(output_vector)[length(output_vector)] <- "'-log(p)'"
  return(output_vector)

},mc.cores=20)






out <- matrix(unlist(out_boots), byrow=TRUE, nrow=length(out_boots))
colnames(out) <- names(out_boots[[1]])
write.table(file=out_file, t(out), sep="\t", quote=FALSE, row.names=FALSE)

bounds90 <- apply(out, 2, function(col) {
  col <- sort(col)
  return(c(mean(col), col[c(25, 250, 475)]))
  })
print(bounds90)
bounds90 <- rbind(startingkds,bounds90)
print(bounds90)
rownames(bounds90) <- c("original","mean", "lower", "median", "upper")
write.table(file=out_ci90_file, t(bounds90), sep="\t", quote = FALSE)

bounds95 <- apply(out, 2, function(col) {
  col <- sort(col)
  return(c(mean(col), col[c(12, 250, 487)]))
  })
bounds95 <- rbind(startingkds,bounds95)
rownames(bounds95) <- rownames(bounds90)
write.table(file=out_ci95_file, t(bounds95), sep="\t", quote = FALSE)



warnings()





