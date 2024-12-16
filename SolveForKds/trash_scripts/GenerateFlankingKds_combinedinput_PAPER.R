################################################################################
#GenerateSiteTypeKds.py
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
site <- args[5]
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
k.c.lib = 100


tick <- 0

# MAIN #########################################################################

# Get siteXcounts, called sitesXcounts:
sitesXcounts <- GetSitesXCountsCombined(mirna, experiment, start, stop, sitelist)

seqs <- sitesXcounts[,1]
sitesXcounts <- sitesXcounts[,-1]

# Remove any rows with no reads whoatsoever.
sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
# Remove any rows for which there are no input reads:
sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]

params <- GetKds(mirna, experiment, start, stop, sitelist, scaled = FALSE)

num.kds <- nrow(sitesXcounts)
num.bgs <- 1

# Get vector of single site counts s.c.
s.c <- as.numeric(sitesXcounts[site, ])

sfXc <- GetSiteFlanksXCounts(mirna, experiment, site, start, stop, sitelist)
original_flanks <- sfXc[,1]
sfXc <- sfXc[,-1]
sfXc <- sfXc[rowSums(sfXc[, 2:6]) > 0,]
colnames(sfXc)[1] <- "I"
sfXc <- round(t((t(sfXc) / colSums(sfXc)) * s.c))
# Assign parameters for Kds and background, subracting 1 from rows due to "None"
# Kd assigned to 1 nM, and the fact that the I and A0 columns are not modeled
# in the experiment, so there's no background term assigned ot them.

# Get the number of flanking sites
num.sf <- nrow(sfXc)
# Get the number of parameters (kds + bg + 1 log(prob))
bgs <- rep(10^params[(num.kds + 1): (num.kds + num.bgs)], 5)
stock.ago <- 10^params[num.kds + num.bgs + 1]

# Get the site kds
kds.s <- params[1:num.kds]
# Omit the site kd for which the flanking sites are being fit.
kds.s <- kds.s[names(kds.s) != site]


# Assign the total site concentration in each experiment, and initialize the
# data column to be used.
# k.c.lib should be 100, for 100 nM.
print(kds.s)
colors_sites <- site_cols[rownames(sitesXcounts)[rownames(sitesXcounts) != site],]

sitesXcounts <- rbind(sitesXcounts, sfXc)
sitesXcounts <- sitesXcounts[rownames(sitesXcounts) != site, ]

c.I.tots <- Norm(sitesXcounts[, 1]) * k.c.lib
names(c.I.tots) <- rownames(sitesXcounts)
# Define the vector of total site type concentrations:
c.totals <- matrix(c.I.tots, nrow=length(c.I.tots), ncol=5, byrow=FALSE)
# Remove the I and A0 columns from the data to be fit to the model. 
data <- sitesXcounts[, 2:6]
rownames(c.totals) <- rownames(data)
colnames(c.totals) <- colnames(data)

colors_flanks <- sapply(rownames(sfXc), GetColorFunction)
colors_all <- c(colors_sites, colors_flanks)
pch_all <- c(rep(1, length(colors_sites)), rep(19, nrow(sfXc)))

c.agos <- sapply(colnames(data), function(x) {
   as.numeric(x) / 100
  })

# Initialize starting Kds, which are set to 1/the enrichment of each site type
# in the A12.6 experiment.
site_colors = sample(colors(),nrow(data))
enrichments.init <- (Norm(c.I.tots)/Norm(rowSums(data)+1))[(nrow(data) - nrow(sfXc) + 1):nrow(data)]
pars.init <- Logit(enrichments.init, max = max(enrichments.init+1))
pars.init <- pars.init - mean(pars.init) + params[site]
names(pars.init) <- rownames(sfXc)
pars <- pars.init
# Define function of just kds and pars.
ModelLikelihood <- function(pars, out_model) {
  kds <- Logistic(c(kds.s, pars), max = 10) 
  # Get the bound counts in each experiment:
  c.bounds <- as.matrix(
    sapply(c.agos, function(ago.percent) {
             return(GetBoundRNA(kds, c.I.tots, ago.percent * stock.ago))
           }
           ))
  c.bounds <<- c.bounds

  # Get the amount of background binding by subtracting the bound from the
  # total sites in each experiment, normalizing. Must transpose to multiply
  # each column.
  c.frees <- c.totals - c.bounds
  c.bgs <- t(t(c.frees) * bgs / colSums(c.frees))
  c.all <- c.bounds + c.bgs
  c.final <- data.frame(t(t(c.all) / colSums(c.all) * colSums(data)))
  c.final <<- c.final
  colnames(c.final) <- colnames(data)
  loglikelihood <<- -1*sum(sapply(colnames(c.totals), function(col_name) {
    dmultinom(x=data[which(c.final[, col_name] > 0), col_name], 
              prob=c.final[which(c.final[, col_name] > 0),col_name],log = TRUE)
    }))

  if (tick%%10000 == 0) {
    print(loglikelihood)
  }
  tick <<- tick + 1
  return(loglikelihood)
}

loglikelihood_initial <- ModelLikelihood(pars, out)
out <- matrix(c(pars, loglikelihood_initial), nrow=1, byrow=TRUE)

colnames(out) <- c(names(pars), "-logp")

specific_figure_string <- paste0(site, "_", start, "-", stop, "_", sitelist, "_",
  "flanking_singlebg_combinedinput_multinomial_PAPER.txt")
## Optimization of parameters:
for (i in seq(1, 2000)) {
  print(i)
  # Define the "scale" vector to be used as the parscale vector.
  loglikelihood_round <- ModelLikelihood(pars, out)
  scale <- sapply(seq(1, length(pars)), function(par_i) {
    if (pars[par_i] > 5) {
      pars[par_i] <<- 5
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

  PlotSiteFlankKdOptimization(out, specific_figure_string, mirna, num.kds, num.bgs, colors = colors_flanks)
  out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                     experiment, "/kds_PAPER/", site, "_", start, "-", stop, "_", 
                     sitelist, 
                     "_flanking_singlebg_combinedinput_multinomial_PAPER.txt")
  
  write.table(file=out_file, out, sep="\t", quote=FALSE, row.names=FALSE)
}


