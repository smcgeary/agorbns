################################################################################
#GenerateSiteTypeKds.py
################################################################################

# Initial parameters and constants.
# args = commandArgs(trailingOnly=TRUE)
# mirna = args[1]
# # Region within random sequence from which site types orginiates,
# # going from position [26 - "start" : 26 + 37 + "stop"]
# start = args[2]
# stop = args[3]
# # Region within the miRNA sequence for which the structural accessibility
# # is being measured. [win_left : win_right], where win_left is the 5' most
# # position withinthe miRNA, and win_right is the 3' most nucleotide position
# #  within the miRNA.
# win_left = args[4]
# win_right = args[5]
# mirna <- "let-7a"
start <- -3
stop <- -3
win_left <- 1
win_right <- 15

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
                       header=FALSE, sep="\t")
k.c.stockago = stockago[mirna,1]
k.c.lib = 100

GetKmerOverlap <- function(seq1, seq2) {
  m <- nchar(seq1)
  n <- nchar(seq2)
  if (n > m) {
    temp <- seq1
    temp_i <- m
    seq1 <- seq2
    m <- n
    seq2 <- temp
    n <- temp_i
  }
  overlaps = c()
  for (i in seq(1, m + n - 1)) {
    m_l <- max(m + 1 - i, 1)
    m_r <- min(m + n - i, m)
    n_l <- max(1, 1 - m + i)
    n_r <- min(i, n)
    if (substr(seq1, m_l, m_r) == substr(seq2, n_l, n_r) &
        m_l != n_l & m_r != n_r) {
      if (m_l == 1) {
        combined_seq <- paste0(substr(seq2, 1, n_r - 1), seq1)
      } else if (m_r == m) {
        combined_seq <- paste0(substr(seq1, 1, m_r - 1), seq2)
      } else {
        combined_seq <- seq1
      }
      if (nchar(combined_seq) <= read_len) {
        overlaps <- c(overlaps, combined_seq)

        }
    }
  }
  return(overlaps)
}



## Functions ###################################################################
## I/O functions

#______________
GetSitesXcounts <- function(mirna) {
  sites_file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
                       mirna, "/equilibrium/full_site_count_tables/all_sites_",
                       start, "-", stop, ".txt")
  sitesXcounts <- read.table(sites_file_name)
  return(sitesXcounts)
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
  residual <- (c.ago - c.freeago - sum(occs_I*c.tots))^2
  return(residual)
}

#_________
GetBoundRNA <- function(kds, c.tots, c.ago) {
  if (c.ago > 0) {
    c.free <- optimize(GetFreeResidual,
                       c(0, c.ago),
                       kds = kds,
                       c.tots = c.tots,
                       c.ago = c.ago)$minimum
    c.bound <- GetOccupancy(c.free, kds)*c.tots
    return(c.bound)
    } else {
      return(rep(0,length(c.tots)))
    }
}# Output Functions:
CheckMaxDifference <- function(out) {
  row.last <- dim(out)[1]
  row.secondtolast <- row.last-1
  diffs <- sapply(1:dim(out)[2],function(x) {
    abs((out[row.last,x]-out[row.secondtolast,x])/out[row.secondtolast,x])
  })
  return(max(diffs))
}





# MAIN #########################################################################
# 1. Get data table:
sitesXcounts <- GetSitesXcounts(mirna)

# Fix column names, with respect to the "A" in front of it.
colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
  colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
    return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
  }
)
print(sitesXcounts)
# print(setdiff(rownames(sitesXcounts),rownames(site_cols)))
# Get data from input library with flanking counts, the average probability of
# pairing, and the variance for this distribution.
sitesXcounts <- sitesXcounts[c(50:10, nrow(sitesXcounts)),]



# Assign parameters for Kds and background, subracting 1 from rows due to "None"
# Kd assigned to 1 nM, and the fact that the I and A0 columns are not modeled
# in the experiment, so there's no background term assigned to them.

# Assign the number of kd parameters and background parameters:
num.kds <- dim(sitesXcounts)[1]-1
num.bgs <- dim(sitesXcounts)[2]-3
# Assign the total site concentration in each experiment, and initialize the
# data column to be used.
# k.c.lib should be 100, for 100 nM.
# Remove the I and A0 columns from the data to be fit to the model. 
data <- sitesXcounts[,2:8]
print('HI')
site_seqs <- as.character(sitesXcounts[-nrow(sitesXcounts),1])
site.indeces <- lapply(site_seqs,GetNucleotideProbIndeces)
print('HI')

site.windows <- read_len - nchar(site_seqs) + 1

double.windows <- outer(nchar(site_seqs),nchar(site_seqs),FUN="+")
GetDoubleWindows <- function(sum) {
  return(2 * sum(seq(1,read_len - sum + 1)))
}

double.windows[] <- vapply(double.windows,FUN=GetDoubleWindows,numeric(1))

co.matrix <- sapply(site_seqs,function(site_1) {
                    sapply(site_seqs, function(site_2) {
                           return(GetKmerOverlap(site_1, site_2))
                           })
  })
print("hi")
test.weights <- rep(1, 4^3)

co.matrix.new <- GetMatrixWeights(co.matrix)
co.matrix.probs <- GetOverlapProbs(co.matrix.new,test.weights)

c.I.tots <- matrix(unlist(rep(data[,1]/sum(data[,1])*k.c.lib,
                  7)), ncol = 7, byrow = FALSE)
# print(c.I.tots)
colnames(c.I.tots) <- colnames(data)
colnames(c.I.tots)[1] <- "0"
colnames(data)[1] <- "0"
rownames(c.I.tots) <- rownames(sitesXcounts)
# Initialize starting Kds, which are set to 1/the enrichment of each site type
# in the A12.6 experiment.
r_numerator <- data[, 2] / sum(data[, 2])
r_denomenator <- sitesXcounts["I"] / sum(sitesXcounts["I"])
kds.init <- (r_numerator / r_denomenator)[1 : num.kds, 1]
weights.init <- rep(0,length(trinucs)-1)
pars.init <- c(-log10(kds.init) - 1, weights.init, -1, -1, -1, -1, -1)

trial <- pars.init
trial[trial == Inf] <- 0
trial[is.na(trial)] <- 0
trial[trial == -1*Inf] <- 0

tick <- 0
# Define function of just kds and pars.
c.totals <- c.I.tots
c.agos <- sapply(colnames(data), function(x) {
   as.numeric(x) / 100 * k.c.stockago
  })
GetModelFrequencies <- function(pars, out_model) {
  # Split up the parameters into the kd and background parameters.
  kds  <- 10^c(pars[1 : num.kds], 1)
  weights <- c(10^c(pars[(num.kds+1):(num.kds+length(trinucs) - 1)]), 1)
  freqs <- c(sapply(site.indeces,GetKmerProb,probs_freqs = weights))*site.windows
  coccur <- outer(freqs,freqs,FUN = "*")*double.windows + GetOverlapProbs(co.matrix.new,weights)
  names(freqs) <- rownames(data)[-nrow(data)]
  rownames(coccur) <- names(freqs)
  colnames(coccur) <- names(freqs)
  # print(rowSums(coccur))
  freqs <- freqs - rowSums(coccur)
  freqs <- c(freqs, max(1 - sum(freqs), 0))
  sapply(1:length(c.agos), function(col) {
    c.totals[,col] <<- freqs*k.c.lib
    })
  bgs  <- c(1,10^pars[(num.kds + length(trinucs)) : (num.kds + length(trinucs) + num.bgs - 1)],1)
  # print(kds)
  # print(bgs)
  # print(c.totals[(nrow(c.totals)-10):nrow(c.totals),])
  names(kds) <- rownames(sitesXcounts)
  # Solve for the free Ago concentration in each experiment.

  # Initialize a matrix with the same total concentration of each site type
  # for each experiment.
  c.bounds <- as.matrix(
    sapply(c.agos, function(x) {
             return(GetBoundRNA(kds, c.totals[, 1], x))
           }
           ))
  # Use the free Ago concentrations to get the amount of each complex bound
  # to Ago.

  # Get the amount of background binding by subtracting the bound from the
  # total sites in each experiment, normalizing. Must transpose to multiply
  # each column.
  c.frees <- c.totals - c.bounds
  c.bgs <- t(t(c.frees) * bgs / colSums(c.frees))
  c.all <- c.bounds + c.bgs
  c.final <- data.frame(t(t(c.all) / colSums(c.all)))
  colnames(c.final) <- colnames(data)
  prob_pois <- sum(sapply(colnames(c.totals), function(col_name) {
    dmultinom(x=data[which(c.final[, col_name] > 0),col_name], prob=c.final[which(c.final[, col_name] > 0),col_name],log = TRUE)
    }))
  if (tick%%10 == 0) {
    out<<- rbind(out, c(pars, -1*prob_pois))
    print(-1*prob_pois)
    print(c.final[nrow(c.final),]/colSums(c.final))
    print(data["None",]/colSums(data))
    # setEPS()
    # postscript(file=paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/figures/kds/",
    # mirna,"/", "structure_fit_", start, "-", stop, "_", win_left, "-", win_right, "_preliminary.eps"))
    plot(c(0.0001,0.0002),
         c(0.0001,0.0002),
         col="white",
         xlim=10^c(-3,10),
         ylim=10^c(0,10),
         log='xy',
         ann = FALSE,
         axes = FALSE)
    title(main=mirna, font.main = 1)
    segments(1,1, x1=1,10^7, lty = 2)
    axis(1,sapply(10^(-3:8), function(x) {x * c(1,2,3,4,5,6,7,8,9)}), labels = FALSE, pos = 1, lwd = 2)
    axis(2,sapply(10^(0:8), function(x) {x * c(1,2,3,4,5,6,7,8,9)}),labels = FALSE, pos = 10^-3,lwd = 2)
    axis(1,10^(-3:8), pos = 1, lwd = 0)
    axis(2,10^(0:8), pos = 10^-3,lwd = 0)


    sapply(1:ncol(data), function(col) {
      x_model <- c.final[, col]
      y_data <- data[, col]
      x_plot <- x_model/sum(x_model)*sum(y_data)
      y_plot <- y_data
      points(x_plot, y_plot,col=c(site_cols[rownames(data),],"black"),lwd=2,pch=20)
    })
    # dev.off()
      out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/equilibrium/kds_with_structure/", k.c.stockago, "_",
                   start, "-", stop, "_", win_left, "-", win_right, "_preliminary.txt")
      write.table(file=out_file, out, sep="\t", quote=FALSE, row.names=FALSE)

    MakeSiteIterationPlot(out,paste0("structure", start, "-", stop, "_", win_left, "-", win_right,"_preliminary"),mirna)  
 
  }
  tick <<- tick + 1
  return(-1*prob_pois)
}
# trial <- out[nrow(out),-ncol(out)]

out <- rbind(c(trial, 100000),c(trial, 100000))
# out <- out[(nrow(out)-2):nrow(out),]

# print(out)
colnames(out) <- c(rownames(sitesXcounts)[-nrow(sitesXcounts)],trinucs[1:(length(trinucs) - 1)],
                   "40", "12.6", "4", "12.6", "0.4","-logp")
# print(out)

print("NOW IT'S OPTIMIZING!!")

Outputmodel <- function(pars) {
  # Split up the parameters into the kd and background parameters.
  kds  <- 10^c(pars[1 : num.kds], 1)
  freqs <- c(10^c(pars[(num.kds+1):(2*num.kds)]), 1)
  freqs <- freqs/sum(freqs)
  sapply(1:length(c.agos), function(col) {
    c.totals[,col] <- freqs*k.c.lib
    })
  bgs  <- c(1,10^pars[(2*num.kds + 1) : (2*num.kds + num.bgs)], 1)
  # print(kds)
  # print(bgs)
  names(kds) <- rownames(sitesXcounts)
  # Solve for the free Ago concentration in each experiment.

  # Initialize a matrix with the same total concentration of each site type
  # for each experiment.
  c.bounds <- as.matrix(
    sapply(c.agos, function(x) {
             return(GetBoundRNA(kds, c.totals[, 1], x))
           }
           ))
  
  # Use the free Ago concentrations to get the amount of each complex bound
  # to Ago.

  # Get the amount of background binding by subtracting the bound from the
  # total sites in each experiment, normalizing. Must transpose to multiply
  # each column.
  c.frees <- c.totals - c.bounds
  c.bgs <- t(t(c.frees) * bgs / colSums(c.frees))
  c.all <- c.bounds + c.bgs
  c.final <- data.frame(t(t(c.all) / colSums(c.all)))
  colnames(c.final) <- colnames(data)
  sapply(ncol(data), function(col) {
    c.final[,col] <- c.final[,col] * sum(data[,col])
    })
  return(c.final)
}
# print(trial)
# trial <- pars_new$par
# # # Solve the first run of the function, and create the output matrix.
pars_new <- optim(trial,
                  GetModelFrequencies,
                  out_model = out,
                  method = "Nelder-Mead",
                  control = list(maxit = 200000))



# Assign output file for the entire sequnce of the optimization and write
# to it.

