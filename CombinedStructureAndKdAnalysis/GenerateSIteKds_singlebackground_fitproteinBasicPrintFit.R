################################################################################
#GenerateSiteTypeKds.py
################################################################################

# Initial parameters and constants.
# args = commandArgs(trailingOnly=TRUE)
# mirna = args[1]
# # # Region within random sequence from which site types orginiates,
# # # going from position [26 - "start" : 26 + 37 + "stop"]
# start = as.integer(args[2])
# stop = as.integer(args[3])
# # # Region within the miRNA sequence for which the structural accessibility
# # # is being measured. [win_left : win_right], where win_left is the 5' most
# # # position withinthe miRNA, and win_right is the 3' most nucleotide position
# # #  within the miRNA.
# win_left = as.integer(args[4])
# win_right = as.integer(args[5])
experiment <- "equilibrium"
# mirna <- "let-7a"
# experiment <- "equilibrium_nb"
start <- 5
stop <- 5
read_len <- 37 + start + stop
# Loads general functions used in AgoRBNS analysis.
source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")



# Loads the colors associated with each site type, for plotting purposes.
site_cols <- read.table(
  "/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/site_info.txt",
  row.names=1, header=FALSE, sep="\t", stringsAsFactors=FALSE)[, 1, drop = FALSE]
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

#______________
GetSitesXcounts <- function(experiment,mirna) {
  sites_file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
                       mirna, "/",experiment,"/full_site_count_tables/all_sites_",
                       start, "-", stop, ".txt")
  print(sites_file_name)
  sitesXcounts <- read.table(sites_file_name)
  return(sitesXcounts)
}



GetEquilParameters <- function(experiment,mirna) {
  file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/equilibrium/kds_with_structure/",
                   start, "-", stop, "_1-15_Basic_singlebg.txt")
  return(read.table(file_name,header = TRUE))
}


# Modeling Functions
#___________
GetOccupancy <- function(c.freeago, kds) {
  # Generates occupancy matrix for the entire input matrix which has rows
  # that are all possible flnaking nucleotide combinations, and columns that
  # are the probabilities of being unpaired in the window.
  #
  # Args:
  # c.freeago: The concentration of free AGO in the binding reaction
  # kds: a list of all kds corresponding to the individual site types including
  #   a non-specific binding kd.
  # IDs: The matrix which gives the site types in column 1 and the flanking
  #    nucleotide identity in column 2.
  # freq.flankXp: The matrix which splits the input into all site and flanking
  #   identities, as well as along the possible average windows for being
  #   unpaired.
  # p.exponent: The exponent used for the weight of geometric average
  #   probability of being unpaired.
  #
  # Returns:
  # A matrix giving the coresponding occuapncy (between 0 and 1) for each
  #   position in the matrix, where 

  return(c.freeago * (c.freeago + kds)^(-1))
}

#______________
GetFreeResidual <- function(c.freeago, kds, c.tots, c.ago) {
  occs_I <- GetOccupancy(c.freeago, kds)
  residual <- (c.ago - c.freeago - sum(occs_I*c.tots))
  if (is.na(residual)) {
    print(kds)
    print(c.ago)
  }
  return(residual)
}

#_________
GetBoundRNA <- function(kds, c.tots, c.ago) {
  if (c.ago > 0) {
    c.free <- NaN
    try(c.free <- uniroot(GetFreeResidual,
                       c(0, c.ago),
                       kds = kds,
                       c.tots = c.tots,
                       c.ago = c.ago,
                       tol = 0.001*.Machine$double.eps^0.25)$root)
    if (is.na(c.free)) {
      c.free <- optimize(GetFreeResidual,
                       c(0, c.ago),
                       kds = kds,
                       c.tots = c.tots,
                       c.ago = c.ago,
                       tol = 0.001*.Machine$double.eps^0.25)$minimum
    }

    
    c.bound <- GetOccupancy(c.free, kds)*c.tots
    return(c.bound)
    } else {
      return(rep(0,length(c.tots)))
    }
}

# Output Functions:




# MAIN #########################################################################
# 1. Get data table:
sitesXcounts <- GetSitesXcounts(experiment,mirna)
colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
  colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
    return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
  }
)
seqs <- sitesXcounts[,1]
sitesXcounts <- sitesXcounts[,-1]
sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]


data <- sitesXcounts[,1:7]

data_full <- data
# ind_current <- seq(1,10)

# data <- data_full[ind_current,]
# data <- rbind(data,rowSums(data_full[-ind_current,]))
rownames(data)[nrow(data)] <- "None"

num.kds <- dim(data)[1]
num.bgs <- 1


c.I.tots <- matrix(unlist(rep(data[,1]/sum(data[,1])*k.c.lib,
                  5)), ncol = 5, byrow = FALSE)
print(c.I.tots)
data <- data[,2:6]
colnames(c.I.tots) <- colnames(data)
# colnames(c.I.tots)[1] <- "0"
# colnames(data)[1] <- "0"
rownames(c.I.tots) <- rownames(data)
# Initialize starting Kds, which are set to 1/the enrichment of each site type
# in the A12.6 experiment.

params <- GetEquilParameters(experiment, mirna)
pars.init <- as.numeric(params[nrow(params), -ncol(params)])
trial <- pars.init
names(trial) <- c(rownames(data), "bg")
# trial[trial == Inf] <- 0
# trial[is.na(trial)] <- 0
# trial[trial == -1*Inf] <- 0

tick <- 0
# Define function of just kds and pars.
c.totals <- c.I.tots
c.agos <- sapply(colnames(data), function(x) {
   as.numeric(x) * k.c.stockago / 100
  })

time_all <- proc.time()
print("hi")
GetModelFrequencies <- function(pars, out_model, blank = FALSE) {
  # print(pars)

  time_0 <- proc.time()
  # Split up the parameters into the kd and background parameters.
  kds  <- c(Logistic(pars[1 : num.kds], 10))
  if (blank == TRUE) {
    kds <- rep(1,num.kds)
  }
  # print(kds)
  # coccur <- outer(freqs,freqs,FUN = "*")*double.windows + GetOverlapProbs(co.matrix.new,weights)
  # print(rowSums(coccur))
  bgs  <- rep(10^pars[(num.kds + 1)],5)
  # print(kds)
  names(kds) <- rownames(data)
  time_1 <- proc.time()
  # Solve for the free Ago concentration in each experiment.

  # Initialize a matrix with the same total concentration of each site type
  # for each experiment.
  c.bounds <- as.matrix(
    sapply(c.agos, function(x) {
             return(GetBoundRNA(kds, c.totals[, 1], x))
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
  print(c.final)
  colnames(c.final) <- colnames(data)
  prob_pois <<- sum(sapply(colnames(c.totals), function(col_name) {
    dmultinom(x=data[which(c.final[, col_name] > 0),col_name], prob=c.final[which(c.final[, col_name] > 0),col_name],log = TRUE)
    }))

  plot(c(0.0001,0.0002),
         c(0.0001,0.0002),
         col="white",
         xlim=10^c(-3,18),
         ylim=10^c(0,10),
         log='xy',
         ann = FALSE,
         axes = FALSE)
    title(main=mirna, font.main = 1)
    title(xlab = "Predicted counts", adj = 0.25, line = 1)
    title(ylab = "Observed counts", adj = 0.4, line = 1)
    segments(1,1, x1=1,10^7, lty = 2)
    segments(1, 1, x1=10^7, 10^7, lty = 3)

    axis(1,sapply(10^(-3:8), function(x) {x * c(1,2,3,4,5,6,7,8,9)}), labels = FALSE, pos = 1, lwd = 2)
    axis(2,sapply(10^(0:8), function(x) {x * c(1,2,3,4,5,6,7,8,9)}),labels = FALSE, pos = 10^-3,lwd = 2)
    axis(1,10^(-3:8), pos = 1, lwd = 0)
    axis(2,10^(0:8), pos = 10^-3,lwd = 0)


    sapply(1:ncol(data), function(col) {
      x_model <- c.all[,col]
      y_data <- data[, col]
      x_plot <- x_model/sum(x_model)*sum(y_data)
      y_plot <- y_data
       points(x_plot, y_plot,col=site_cols[rownames(data),], pch=20)
    })

    legend(x = 10^9, y = 10^9, legend = rownames(data), col = site_cols[rownames(data),], pch = 20, cex = 0.95, bty = "n", ncol = 3)



  return(-1*prob_pois)
}
dev.new(xpos = 20, ypos = 20, width = 12, height = 8)

GetModelFrequencies(trial)
dev.copy2pdf(file = paste0("170329 ", mirna, "_equilibrium_1bg_fixedprotein_fit.pdf"))
# break
GetModelFrequencies(trial, blank = TRUE)
dev.copy2pdf(file = paste0("170329 ", mirna, "_equilibrium_1bg_fixedprotein_null.pdf"))


break

prob_init <- -GetModelFrequencies(trial)

GetConfidenceParam <- function(ind){
  pars <- trial
  par_ind <- pars[ind]
  range <- c()
  width <- 1.0
  while (length(range)<= 300){
    width <- width / 10
  limits <- sapply(seq(par_ind - width, par_ind + width, length = 1000), function(par) {
    par_use <- pars
    par_use[ind] <- par
    return(c(par, 2*(prob_init + GetModelFrequencies(par_use))))
    })
  range <- which(limits[2,] <= qchisq(0.95,1))
}
  return(c(min(limits[1,range]),max(limits[1,range])))
}

# dev.copy2pdf(file = paste0("170329 ", mirna, "_equilibrium_null.pdf"))


error_bars <- sapply(1:length(trial),GetConfidenceParam)


# out <- read.table(file = paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
#            mirna, "/",experiment,"/kds_with_structure/", k.c.stockago,
#            "_Namita_all_sites"),stringsAsFactors = FALSE)


# trial <- out[nrow(out),-ncol(out)]

# out <- rbind(c(trial, 100000),c(trial, 100000))
# colnames(out) <- c(rownames(data)[-nrow(data)],
#                    "40", "12.6", "4", "12.6", "0.4","AGO","-logp")

# initial_prob <- GetModelFrequencies(trial, out)
# out <- rbind(c(trial, initial_prob),c(trial, initial_prob))

# print(out)
# colnames(out) <- c(rownames(data)[-nrow(data)],
#                    "40", "12.6", "4", "12.6", "0.4","AGO","-logp")



# Assign output file for the entire sequnce of the optimization and write
# to it.

