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
print(experiment)

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
sitesXcounts <- GetSitesXCounts(mirna,
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

# Remove any rows with no reads whoatsoever.
sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
# Remove any rows for which there are no input reads:
sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]

data <- sitesXcounts[,2:6]
data <- data.matrix(data)

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

# PARAMETER INITIALIZATION:
# Initialize starting Kds, which are set to 1/the enrichment of each site type
# in the A12.6 experiment.



freqs_input <- 
r_numerator <- Norm(rowSums(data))
r_denomenator <- Norm(c.I.tots)
kds.init <- c(r_denomenator) / c(r_numerator)
kds.init <- kds.init/max(kds.init)
pars.init <- c(Logit(kds.init[-length(kds.init)], 10), -1, 1)

rand_length <- nchar(rownames(data)[1]) + mirna.start - mirna.stop - 1

print(rand_length)

pars_freq <- c(1, 2, 3, 4, 5, 6)

AssignFrequencies

if (sitelist %in% kmer_list_names[-1]) {
  list_ind <- which(kmer_list_names == sitelist) - 1
  sitelist_prior <- kmer_list_names[list_ind]
  in_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                       experiment, "/kds_PAPER/", start, "-", stop, "_", 
                       sitelist_prior, "_", mirna.start, "-", mirna.stop,
                       "_singlebg_multinomial2_PAPER_last.txt")
  pars_raw <- fread(in_file, sep="\t", header=TRUE, colClasses="numeric", 
                    stringsAsFactors=FALSE)
  pars_table <- matrix(unlist(pars_raw), nrow=nrow(pars_raw),
                       ncol=ncol(pars_raw), byrow=FALSE)
  pars_old <- unlist(pars_table[nrow(pars_table),])
  names(pars_old) <- colnames(pars_raw)

  kds.init <- sapply(rownames(data)[-nrow(data)], function(name) {
    if (list_ind %% 2 == 1) {
          smaller_name <- substr(name, 2, nchar(name))
    } else {
          smaller_name <- substr(name, 1, nchar(name)-1)
    }
    sub_name <- paste0(smaller_name)
    ind <- which(names(pars_old) == sub_name)
    return(pars_old[ind])  
  })

  pars.init <- c(kds.init, -1, 1)
}

names(pars.init) <- c(rownames(data)[-nrow(data)], "bg", "AGO")

pars <- pars.init
tick <- 1

print("up to model")

if (sitelist %in% kmer_list_names) {
  n <- 256
  colors <- c(rainbow_hcl(n, c=50, l=70, start=0, end=360 * (n - 1) / n),"red",
              "black")

  out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                     experiment, "/kds_PAPER/", start, "-", stop, "_", 
                     sitelist, "_", mirna.start, "-", mirna.stop,
                     "_singlebg_multinomial2_PAPER.txt")
  out_last_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                          experiment, "/kds_PAPER/", start, "-", stop, "_", 
                          sitelist, "_", mirna.start, "-", mirna.stop,
                          "_singlebg_multinomial2_PAPER_last.txt")

  specific_figure_string <- paste0(start, "-", stop, "_", sitelist, "_",
                                   mirna.start, "-", mirna.stop,
                                   "_singlebg_multinomial2_PAPER")
} else {
  colors <- FALSE

  out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                     experiment, "/kds_PAPER/", start, "-", stop, "_", 
                     sitelist, "_singlebg_multinomial2_PAPER.txt")

  out_last_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                          experiment, "/kds_PAPER/", start, "-", stop, "_", 
                          sitelist, "_singlebg_multinomial2_PAPER_last.txt")

  specific_figure_string <- paste0(start, "-", stop, "_", sitelist,
                                   "_singlebg_multinomial2_PAPER")
}


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

ModelLikelihood <- function(pars){

  model <- ModelFunction(pars)
  model_norm <- t(t(model) / colSums(model))
  sumofsquares <- -sum(data*log(model_norm))

  if (tick %% 20 == 0 & tick > 1) {
  out <<- rbind(out, c(pars, sumofsquares))
    if (sitelist %in% kmer_list_names) {
      PlotSiteKdOptimization(out[,c(seq(1,num.kds,length=256),(ncol(out)-2):ncol(out))], specific_figure_string, mirna, 256, num.bgs,
        colors=colors)

    } else {
      PlotSiteKdOptimization(out, specific_figure_string, mirna, num.kds, num.bgs, colors=colors)
    }
    write.table(file=out_file, out, sep="\t", quote=FALSE, row.names=FALSE)
    write.table(file=out_last_file, out[nrow(out),,drop=FALSE], sep="\t",
                quote=FALSE, row.names=FALSE)
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



out <- matrix(c(pars, ModelLikelihood(pars)), nrow=1)
out_names <- c(rownames(data)[-length(rownames(data))],
               "bg", "AGO", "-logp")
colnames(out) <- out_names

QuickDeriv <-function(params,width) {
  params_orig <- params
  sapply(seq(length(params)), function(i) {

    params <- params_orig
    params[i] <- params[i] + width
    return((ModelLikelihood(params) - ModelLikelihood(params_orig))/(width))
  })
}



solution <- optim(pars,
                  ModelLikelihood,
                  gr = Gradient,
                  method = "L-BFGS-B",
                  lower=c(rep(-10, length=num.kds), -2, -1),
                  upper=c(rep(10, length=num.kds), 1, 2),
                  control = list(maxit=100000, factr=1e2))

pars <- solution$par
sumoflogsquares_current <- solution$value
out <- rbind(out, c(pars, sumofsquares))

print(sumoflogsquares_current)
write.table(file=out_file, out, sep="\t", quote=FALSE, row.names=FALSE)
write.table(file=out_last_file, out[nrow(out),,drop=FALSE], sep="\t",
                quote=FALSE, row.names=FALSE)

warnings()





