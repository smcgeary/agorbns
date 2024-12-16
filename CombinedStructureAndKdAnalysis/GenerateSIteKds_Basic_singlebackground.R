################################################################################
#GenerateSiteTypeKds.py
################################################################################

# Initial parameters and constants.
args = commandArgs(trailingOnly=TRUE)
mirna = args[1]
# # Region within random sequence from which site types orginiates,
# # going from position [26 - "start" : 26 + 37 + "stop"]
start = as.integer(args[2])
stop = as.integer(args[3])
sitelist = args[4]
print(sitelist)
experiment <- "equilibrium"
# # Region within the miRNA sequence for which the structural accessibility
# # is being measured. [win_left : win_right], where win_left is the 5' most
# # position withinthe miRNA, and win_right is the 3' most nucleotide position
# #  within the miRNA.
# win_left = args[4]
# win_right = args[5]
# mirna <- "let-7a"
# experiment <- "equilibrium_nb"
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
                       header=TRUE, sep="\t")
k.c.stockago = stockago[mirna,experiment]
k.c.lib = 100

## Functions ###################################################################
## I/O functions

#______________
GetSitesXcounts <- function(experiment,mirna) {
  sites_file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
                       mirna, "/",experiment,"/full_site_count_tables/all_sites_",
                       start, "-", stop, "_", sitelist, ".txt")
  print(sites_file_name)
  sitesXcounts <- read.table(sites_file_name)
  return(sitesXcounts)
}


ConvertParameters <- function(out, row) {
  pars_row <- out[row,-ncol(out)]
  return(c(Logistic(pars_row[1:num.kds], max = 10), 10^(pars_row[(num.kds + 1) : length(pars_row)])))
}

GetBackground <- function(out, row) {
  pars <- ConvertParameters(out, row)
  return(pars[c("0.4", "1.26", "4", "12.6", "40")])
}

# break
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
sitesXcounts <- GetSitesXcounts(experiment,mirna)
colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
  colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
    return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
  }
)
print(sitesXcounts)
seqs <- sitesXcounts[,1]
sitesXcounts <- sitesXcounts[,-1]
sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]

cols_new = c("magenta","purple","blue","green","red","orange","yellow","black")
names(cols_new) <- c(4, 5, 6, 7, 8, 9, 10,"None")

cols_final <- cols_new[sapply(rownames(sitesXcounts), function(name) { unlist(strsplit(name, "mer"))[1]})]

data <- sitesXcounts[,1:7]

data_full <- data
# ind_current <- seq(1,10)

# data <- data_full[ind_current,]
# data <- rbind(data,rowSums(data_full[-ind_current,]))
rownames(data)[nrow(data)] <- "None"

num.kds <- dim(data)[1]
num.bgs <- dim(data)[2]-2


c.I.tots <- matrix(unlist(rep(data[,1]/sum(data[,1])*k.c.lib,
                  5)), ncol = 5, byrow = FALSE)
# print(c.I.tots)
data <- data[,2:6]
colnames(c.I.tots) <- colnames(data)
# colnames(c.I.tots)[1] <- "0"
# colnames(data)[1] <- "0"
rownames(c.I.tots) <- rownames(data)
# Initialize starting Kds, which are set to 1/the enrichment of each site type
# in the A12.6 experiment.
r_numerator <- data[, 2] / sum(data[, 2])
r_denomenator <- c.I.tots[,1] / sum(c.I.tots[,1])
kds.init <- ((r_numerator / r_denomenator)[1 : num.kds] + 1)^(-1)
pars.init <- c(Logit(kds.init, max = 1), -1)
names(pars.init) <- c(rownames(data), "bg")
trial <- pars.init
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

GetModelFrequencies <- function(pars, out_model) {
  # print(pars)

  time_0 <- proc.time()
  # Split up the parameters into the kd and background parameters.
  kds  <- Logistic(pars[1 : num.kds], 10)
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
  # print(c.final)
  colnames(c.final) <- colnames(data)
  time_4 <- proc.time()
  prob_pois <<- sum(sapply(colnames(c.totals), function(col_name) {
    dmultinom(x=data[which(c.final[, col_name] > 0),col_name], prob=c.final[which(c.final[, col_name] > 0),col_name],log = TRUE)
    }))
  time_5 <- proc.time()

  if (tick%%10000 == 0) {
    print(-1*prob_pois)
    # out<<- rbind(out, c(pars, -1*prob_pois))

    # print(tick / (proc.time() - time_all)[3])
    # print(out[nrow(out),ncol(out)] - out[nrow(out)-1,ncol(out)])
    # print(bgs)
    # print(10^pars[length(pars)])
    # setEPS()
    # postscript(file=paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/figures/kds/",
    # mirna,"/", "structure_fit_", start, "-", stop, "_", win_left, "-", win_right, "_preliminary.eps"))
    # plot(c(0.0001,0.0002),
    #      c(0.0001,0.0002),
    #      col="white",
    #      xlim=10^c(-3,10),
    #      ylim=10^c(0,10),
    #      log='xy',
    #      ann = FALSE,
    #      axes = FALSE)
    # title(main=mirna, font.main = 1)
    # segments(1,1, x1=1,10^7, lty = 2)
    # segments(0.01,0.01,10000000,10000000, lty = 2)

    # axis(1,sapply(10^(-3:8), function(x) {x * c(1,2,3,4,5,6,7,8,9)}), labels = FALSE, pos = 1, lwd = 2)
    # axis(2,sapply(10^(0:8), function(x) {x * c(1,2,3,4,5,6,7,8,9)}),labels = FALSE, pos = 10^-3,lwd = 2)
    # axis(1,10^(-3:8), pos = 1, lwd = 0)
    # axis(2,10^(0:8), pos = 10^-3,lwd = 0)


    # sapply(1:5, function(col) {
    #   x_model <- c.final[, col]
    #   y_data <- data[, col]
    #   x_plot <- x_model/sum(x_model)*sum(y_data)
    #   y_plot <- y_data
    #    points(x_plot, y_plot,col=site_cols[rownames(c.final),], lwd=2,pch=20)
    # })

    # plot(c(0.0001,0.0002),
    #      c(0.0001,0.0002),
    #      col="white",
    #      xlim=10^c(-3,1),
    #      ylim=10^c(0,10),
    #      log='xy',
    #      ann = FALSE,
    #      axes = FALSE)
    # title(main=mirna, font.main = 1)

    # axis(1,sapply(10^(-3:1), function(x) {x * c(1,2,3,4,5,6,7,8,9)}), labels = FALSE, pos = 1, lwd = 2)
    # axis(2,sapply(10^(0:8), function(x) {x * c(1,2,3,4,5,6,7,8,9)}),labels = FALSE, pos = 10^-3,lwd = 2)
    # axis(1,10^(-3:1), pos = 1, lwd = 0)
    # axis(2,10^(0:8), pos = 10^-3,lwd = 0)


    # sapply(1:nrow(data), function(row) {
    #   # print("hi")
    #   x_model <- c.final[row,]
    #   y_data <- data[row,]
    #    lines(c(0.4, 1.26, 4, 12.6, 40)/100, x_model,col=site_cols[rownames(c.final)[row],], lwd = 2)
    #    points(c(0.4, 1.26, 4, 12.6, 40)/100, y_data,col=site_cols[rownames(c.final)[row],],pch=20)

    # })



    # dev.off()

    MakeSiteIterationPlot(out,paste0(experiment,"_", start, "-", stop, "_", win_left, "-", win_right, "_", sitelist, "_Basic_fixedproteinonebackground"),mirna, num.kds, num.bgs)  
  # time_kds <- time_1 - time_0
  # print("time_kds")
  # print(time_kds)

  # time_bound <- time_2 - time_1
  # print("time to solve free ago and get bound RNA:")
  # print(time_bound)


  #  time_get_probabilities <- time_4 - time_2
  # print("Time to add background, and generate multinomial probs:")
  # print(time_get_probabilities)


  # time_fit <- time_5 - time_4
  # print("time to get -logp for fit with data:")
  # print(time_fit)


  }

  tick <<- tick + 1
  return(-1*prob_pois)
}

# out <- read.table(file = paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
#            mirna, "/",experiment,"/kds_with_structure/", k.c.stockago,
#            "_Namita_all_sites"),stringsAsFactors = FALSE)


# trial <- out[nrow(out),-ncol(out)]

out <- rbind(c(trial, 100000),c(trial, 100000))
colnames(out) <- c(rownames(data),
                   "bg","-logp")
 
initial_prob <- GetModelFrequencies(trial, out)
out <- rbind(c(trial, initial_prob),c(trial, initial_prob))

print(out)
colnames(out) <- c(rownames(data),
                   "bg","-logp")
print(out)

print("NOW IT'S OPTIMIZING!!")

time_all <- proc.time()

# print(trial)
# trial <- pars_new$par
# # # Solve the first run of the function, and create the output matrix.
for (i in seq(1, 2000)) {
  initial_prob <- GetModelFrequencies(trial, out)
  scale <- sapply(seq(1, length(trial)), function(par_i) {
    if (trial[par_i] > 10) {
      trial[par_i] <<- 10
    }
    temp <- trial
    temp[par_i] <- temp[par_i] - 1
    func_out <- abs(GetModelFrequencies(temp, out) - initial_prob)
    return(func_out)
  })
  scale[is.na(scale)] <- 0
  scale <- scale + 0.001
  scale <- abs(scale)
  names(scale) <- names(trial)
  pars_new <- optim(trial,
                    GetModelFrequencies,
                    out_model = out,
                    method = "Nelder-Mead",
                    control = list(maxit = 5000, alpha = 1, beta = 0.5, gamma = 2, parscale = scale^(-1/6)))
print("ADD")
# ind_new <- seq(1,max(ind_current)+10)

# data <- data_full[ind_new,]
# data <- rbind(data,rowSums(data_full[-ind_new,]))
# rownames(data)[nrow(data)] <- "None"
# num.kds <- dim(data)[1]-1
# num.bgs <- dim(data)[2]-2

# c.I.tots <- matrix(unlist(rep(data[,1]/sum(data[,1])*k.c.lib,
#                   5)), ncol = 5, byrow = FALSE)
# # print(c.I.tots)
# data <- data[,2:6]
# colnames(c.I.tots) <- colnames(data)
# # colnames(c.I.tots)[1] <- "0"
# # colnames(data)[1] <- "0"
# rownames(c.I.tots) <- rownames(data)
# c.totals <- c.I.tots
# print(trial[1:num.kds])
# # Initialize starting Kds, which are set to 1/the enrichment of each site type

trial <- pars_new$par
    out<<- rbind(out, c(trial, -1*prob_pois))
          out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/", experiment, "/kds_with_structure/",
                   start, "-", stop, "_", win_left, "-", win_right, "_", sitelist, "_Basic_singlebg.txt")
      write.table(file=out_file, out, sep="\t", quote=FALSE, row.names=FALSE)


}


# Assign output file for the entire sequnce of the optimization and write
# to it.

