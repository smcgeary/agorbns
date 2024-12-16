################################################################################
#GenerateSiteKds_combinedinput_PAPER.py
################################################################################

## THIS IS THE OFFICIAL SCRIPT FOR FINALIZING THE MANUSCRIPT.
## THIS SCRIPT HAS THESE DECISIONS MADE:
## 1.) PROTEIN IS FIT AS A PARAMETER
## 2.) A SINGLE BACKGROUND PARAMETER IS GIVEN TO ALL FIVE CONCENTRATIONS.
## 3.) SITELIST, START, AND STOP POSITION ARE MANDATORY PARAMETERS.
## 4.) ALL INPUT IS NOT COMBINED FOR THESE ANALYSES.

# Initial parameters and constants.
# args  <- commandArgs(trailingOnly=TRUE)
# mirna <- args[1]
# # Region within random sequence from which site types orginiates,
# # going from position [26 - "start" : 26 + 37 + "stop"]
# start <- as.integer(args[2])
# stop  <- as.integer(args[3])

# # # Parameter specifying which list of sites is used for the analysis.
# # # I.E "Current" includes all the current site types, "canonical" is just
# # # the six seed-based 6mers, 7mers, and 8mers, along with no site., etc.
# sitelist <- args[4]
# if (sitelist == "12mers") {
#   mirna.start <- as.integer(args[5])
#   mirna.stop <- as.integer(args[6])
# } else {
#   mirna.start <- NULL
#   mirna.stop <- NULL
# }
# Experiment name

mirna <- "let-7a"
start <- 5
stop <- 5
sitelist <- "12mers"
mirna.start <- 2
mirna.stop <- 5
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
if (sitelist == "12mers") {
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
print(c.agos)
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
      time_1 <- proc.time()

  # Split up the parameters into the kd and background parameters.
  kds  <- Logistic(pars[1 : num.kds], 10)
  names(kds) <- rownames(data)
  bgs  <- rep(10^pars[(num.kds + 1)], 5)
  stock.ago <- 10^pars[num.kds + 2]

  # Get the bound counts in each experiment:
      time_2 <- proc.time()

  c.bounds <- as.matrix(
    sapply(c.agos, function(ago.percent) {
             return(GetBoundRNA(kds, c.I.tots, ago.percent * stock.ago))
           }
           ))
      time_3 <- proc.time()

  # Get the amount of background binding by subtracting the bound from the
  # total sites in each experiment, normalizing. Must transpose to multiply
  # each column.
  c.frees <- c.totals - c.bounds
  c.bgs <- t(t(c.frees) * bgs / colSums(c.frees))
  c.all <- c.bounds + c.bgs
  c.final <- data.frame(t(t(c.all) / colSums(c.all) * colSums(data)))
  c.final <<- c.final
  colnames(c.final) <- colnames(data)
  loglikelihood <- sum((unlist(c.final) - unlist(data))^2)

  if (tick%%1 == 0) {
          time_4 <- proc.time()

    print(loglikelihood)
    print(dim(c.final))
    print(dim(data))
    print(time_4 - time_1)
    print(time_3 - time_2)
    plot(c(0.01, 10000000), c(0.01, 10000000), log = 'xy')
    sapply(1:ncol(c.final), function(column) {
      points(c.final[,column], data[,column])
    })
  }
  tick <<- tick + 1
  return(loglikelihood)
}

ObjectiveFunctionGradient <- function(pars) {
  kds  <- Logistic(pars[1 : num.kds], 10)
  names(kds) <- rownames(data)
  bgs  <- rep(10^pars[(num.kds + 1)], 5)
  stock.ago <- 10^pars[num.kds + 2]
  c.bounds <- as.matrix(
    sapply(c.agos, function(ago.percent) {
             return(GetBoundRNA(kds, c.I.tots, ago.percent * stock.ago))
           }
           ))
      time_3 <- proc.time()
  c.free.agos <- sapply(c.agos, function(ago.percent) {
    GetFreeAgo(kds, c.I.tots, ago.percent * stock.ago)
  })
  c.free.agos <<- c.free.agos
  weird_gradient_thing <- sapply(1:length(c.free.agos), function(col_ind){
    (1 + sum(kds*c.I.tots * (c.free.agos[col_ind] + kds)^(-2)))^-1
    })
  print(weird_gradient_thing)
  print(c.free.agos)
  print(c.free.agos + kds)
  # Get the amount of background binding by subtracting the bound from the
  # total sites in each experiment, normalizing. Must transpose to multiply
  # each column.
  c.frees <- c.totals - c.bounds
  c.bgs <- t(t(c.frees) * bgs / colSums(c.frees))
  c.all <- c.bounds + c.bgs
  c.final <- data.frame(t(t(c.all) / colSums(c.all) * colSums(data)))



}
# ObjectiveFunctionGradient(pars)
# break
loglikelihood_initial <- ModelLikelihood(pars, out)
out <- matrix(c(pars, loglikelihood_initial), nrow=1, byrow=TRUE)

colnames(out) <- c(rownames(data),
                   "bg", "AGO", "-logp")
print(out)

specific_figure_string <- paste0(start, "-", stop, "_", sitelist, "_",
  "singlebg_multinomial_PAPER.txt")

if (sitelist == "12mers") {
  specific_figure_string <- paste0(start, "-", stop, "_", sitelist, "_",
                                   mirna.start, "-", mirna.stop,
                                   "_singlebg_multinomial_PAPER.txt")

} else {
  specific_figure_string <- paste0(start, "-", stop, "_", sitelist, "_",
                                   "singlebg_multinomial_PAPER_test.txt")
}

## Optimization of parameters:
for (i in seq(1, 2000)) {
  print(i)
  # Define the "scale" vector to be used as the parscale vector.
  loglikelihood_round <- ModelLikelihood(pars, out)
  scale <- sapply(seq(1, length(pars)), function(par_i) {
    if (pars[par_i] > 10) {
      pars[par_i] <<- 10
    }
    temp <- pars
    temp[par_i] <- temp[par_i] - 1
    func_out <- abs(ModelLikelihood(temp, out) - loglikelihood_round)
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
  loglikelihood_current <- solution$value
  print(loglikelihood_current)

  out <- rbind(out, c(pars, loglikelihood_current))

  PlotSiteKdOptimization(out, specific_figure_string, mirna, num.kds, num.bgs)
  # if (sitelist == "12mers") {
  #   out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
  #                      experiment, "/kds_PAPER/", start, "-", stop, "_", 
  #                      sitelist, "_", mirna.start, "-", mirna.stop,
  #                      "_singlebg_multinomial_PAPER_test.txt")

  # } else {
  #   out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
  #                      experiment, "/kds_PAPER/", start, "-", stop, "_", 
  #                      sitelist, 
  #                      "_singlebg_multinomial_PAPER.txt")
  # }
  # write.table(file=out_file, out, sep="\t", quote=FALSE, row.names=FALSE)
}



