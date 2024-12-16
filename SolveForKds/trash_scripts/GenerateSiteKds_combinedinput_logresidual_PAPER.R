################################################################################
#GenerateSiteKds_combinedinput_PAPER.py
################################################################################

## THIS IS THE OFFICIAL SCRIPT FOR FINALIZING THE MANUSCRIPT.
## THIS SCRIPT HAS THESE DECISIONS MADE:
## 1.) PROTEIN IS FIT AS A PARAMETER
## 2.) A SINGLE BACKGROUND PARAMETER IS GIVEN TO ALL FIVE CONCENTRATIONS.
## 3.) SITELIST, START, AND STOP POSITION ARE MANDATORY PARAMETERS.
## 4.) ALL INPUTS ARE COMBINED ACROSS FOR THESE ANALYSES.

# Initial parameters and constants.
args  <- commandArgs(trailingOnly=TRUE)
mirna <- args[1]
# Region within random sequence from which site types orginiates,
# going from position [26 - "start" : 26 + 37 + "stop"]
start <- as.integer(args[2])
stop  <- as.integer(args[3])

# # Parameter specifying which list of sites is used for the analysis.
# # I.E "Current" includes all the current site types, "canonical" is just
# # the six seed-based 6mers, 7mers, and 8mers, along with no site., etc.
sitelist <- args[4]
# Experiment name
experiment <- "equilibrium"

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
sitesXcounts <- GetSitesXCountsCombined(mirna,
                                        experiment,
                                        start,
                                        stop,
                                        sitelist)
# Separate site sequences from data file.
seqs <- sitesXcounts[,1]
sitesXcounts <- sitesXcounts[,-1]

# Remove any rows with no reads whoatsoever.
sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
# Remove any rows for which there are no input reads:
sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]

data <- sitesXcounts[,2:6]
rownames(data)[nrow(data)] <- "None"

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

ModelLikelihood <- function(pars, out_model) {
  # Split up the parameters into the kd and background parameters.
  kds  <- Logistic(pars[1 : num.kds], 10)
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
  c.final <<- c.final
  colnames(c.final) <- colnames(data)
  use <- which(c.final > 0 & data > 0)
  sumoflogsquares <<- sum((log(unlist(c.final)[use]) - log(unlist(data)[use]))^2)
  if (tick%%10000 == 0) {
    print(sumoflogsquares)
  }
  tick <<- tick + 1
  return(sumoflogsquares)
}

 
sumoflogsquares_initial <- ModelLikelihood(pars, out)
out <- matrix(c(pars, sumoflogsquares_initial), nrow=1, byrow=TRUE)

colnames(out) <- c(rownames(data),
                   "bg", "AGO", "-logp")

specific_figure_string <- paste0(start, "-", stop, "_", sitelist, "_",
  "singlebg_combinedinput_logresidual_PAPER.txt")

## Optimization of parameters:
for (i in seq(1, 2000)) {
  print(i)
  # Define the "scale" vector to be used as the parscale vector.
  sumoflogsquares_round <- ModelLikelihood(pars, out)
  scale <- sapply(seq(1, length(pars)), function(par_i) {
    if (pars[par_i] > 10) {
      pars[par_i] <<- 10
    }
    temp <- pars
    temp[par_i] <- temp[par_i] - 1
    func_out <- abs(ModelLikelihood(temp, out) - sumoflogsquares_round)
    return(func_out)
  })
  scale[is.na(scale)] <- 0
  scale <- scale + 0.001
  scale <- abs(scale)
  names(scale) <- names(pars)
  solution <- optim(pars,
                    ModelLikelihood,
                    out_model = out,
                    method = "Nelder-Mead",
                    control = list(maxit = 5000, alpha = 1, beta = 0.5, gamma = 2, parscale = scale^(-1/6)))
  pars <- solution$par
  sumoflogsquares_current <- solution$value
  print(sumoflogsquares_current)
  out <- rbind(out, c(pars, sumoflogsquares_current))

  PlotSiteKdOptimization(out, specific_figure_string, mirna, num.kds, num.bgs)
  out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                     experiment, "/kds_PAPER/", start, "-", stop, "_", 
                     sitelist, 
                     "_singlebg_combinedinput_logresidual_PAPER.txt")
  
  write.table(file=out_file, out, sep="\t", quote=FALSE, row.names=FALSE)
}



