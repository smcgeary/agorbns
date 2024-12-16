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

# Parameter specifying which list of sites is used for the analysis.
# I.E "Current" includes all the current site types, "canonical" is just
# the six seed-based 6mers, 7mers, and 8mers, along with no site., etc.
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


## Functions ###################################################################
## I/O functions


# MAIN #########################################################################
# 1. Get data table:
sitesXcounts <- GetSitesXcountsCombined(experiment, 
                                        mirna,
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


data <- sitesXcounts[,3:8]
rownames(data)[nrow(data)] <- "None"

num.kds <- nrow(data)
num.bgs <- 1

c.I.tots <- Norm(sitesXcounts[, 2]) * k.c.lib
c.totals <- matrix(c.I.tots, nrow=length(c.I.tots), ncol=5, byrow=FALSE)


data <- data[, 2 : 6]
colnames(c.I.tots) <- colnames(data)
rownames(c.I.tots) <- rownames(data)

# PARAMETER INITIALIZATION:
# Initialize starting Kds, which are set to 1/the enrichment of each site type
# in the A12.6 experiment.
r_numerator <- Norm(data[,2])
r_denomenator <- Norm(c.I.tots[,1])
kds.init <- ((r_numerator / r_denomenator) + 1)^(-1)
pars.init <- c(Logit(kds.init, max = 1), -1, log10(k.c.stockago))
names(pars.init) <- c(rownames(data), "bg", "AGO")

pars <- pars.init

tick <- 0

# Define the vector of total site type concentrations:
# Define the vector of AGO concentrations:
c.agos <- sapply(colnames(data), function(x) {
   as.numeric(x) / 100
  })

time_all <- proc.time()


print(data)
print(c.totals)
GetModelFrequencies <- function(pars, out_model) {
  time_0 <- proc.time()
  # Split up the parameters into the kd and background parameters.
  kds  <- Logistic(pars[1 : num.kds], 10)
  names(kds) <- rownames(data)
  bgs  <- rep(10^pars[(num.kds + 1)],5)
  stock.ago <- 10^pars[num.kds + 2]
  time_1 <- proc.time()

  # Get the bound counts in each experiment:
  c.bounds <- as.matrix(
    sapply(c.agos, function(x) {
             return(GetBoundRNA(kds, c.totals[, 1], x*stock.ago))
           }
           ))

  time_2 <- proc.time()

  # Use the free Ago concentrations to get the amount of each complex bound
  # to Ago.

  # Get the amount of background binding by subtracting the bound from the
  # total sites in each experiment, normalizing. Must transpose to multiply
  # each column.
  c.frees <- c.totals - c.bounds
  c.bgs <- t(t(c.frees) * bgs / colSums(c.frees))
  c.all <- c.bounds + c.bgs
  c.final <- data.frame(t(t(c.all) / colSums(c.all) * colSums(data)))
  # print(c.final)
  colnames(c.final) <- colnames(data)
  time_4 <- proc.time()
  prob_pois <<- sum(sapply(colnames(c.totals), function(col_name) {
    dmultinom(x=data[which(c.final[, col_name] > 0),col_name], prob=c.final[which(c.final[, col_name] > 0),col_name],log = TRUE)
    }))
  time_5 <- proc.time()

  if (tick%%10000 == 0) {
    print(-1*prob_pois)


    # MakeSiteIterationPlot(out,paste0(experiment,"_", start, "-", stop, "_", win_left, "-", win_right, "_", sitelist,"_Basic_notfixedproteinonebackground_combinedinput"),mirna, num.kds, num.bgs)  


  }

  tick <<- tick + 1
  return(-1*prob_pois)
}



out <- rbind(c(pars, 100000),c(pars, 100000))
colnames(out) <- c(rownames(data),
                   "bg", "AGO", "-logp")
 
initial_prob <- GetModelFrequencies(pars, out)
out <- rbind(c(pars, initial_prob),c(pars, initial_prob))

print(out)
colnames(out) <- c(rownames(data),
                   "bg", "AGO", "-logp")
print(out)

print("NOW IT'S OPTIMIZING!!")

time_all <- proc.time()


# # # Solve the first run of the function, and create the output matrix.
for (i in seq(1, 2000)) {
  initial_prob <- GetModelFrequencies(pars, out)
  scale <- sapply(seq(1, length(pars)), function(par_i) {
    if (pars[par_i] > 10) {
      pars[par_i] <<- 10
    }
    temp <- pars
    temp[par_i] <- temp[par_i] - 1
    func_out <- abs(GetModelFrequencies(temp, out) - initial_prob)
    return(func_out)
  })
  scale[is.na(scale)] <- 0
  scale <- scale + 0.001
  scale <- abs(scale)
  names(scale) <- names(pars)
  pars_new <- optim(pars,
                    GetModelFrequencies,
                    out_model = out,
                    method = "Nelder-Mead",
                    control = list(maxit = 5000, alpha = 1, beta = 0.5, gamma = 2, parscale = scale^(-1/6)))
print("ADD")

# # Initialize starting Kds, which are set to 1/the enrichment of each site type

pars <- pars_new$par
    out<<- rbind(out, c(pars, -1*prob_pois))
          out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/", experiment, "/kds_with_structure/",
                   start, "-", stop, "_", win_left, "-", win_right, "_", sitelist, "_Basic_singlebg_combinedinput_PAPER.txt")
      write.table(file=out_file, out, sep="\t", quote=FALSE, row.names=FALSE)


}


# Assign output file for the entire sequnce of the optimization and write
# to it.

