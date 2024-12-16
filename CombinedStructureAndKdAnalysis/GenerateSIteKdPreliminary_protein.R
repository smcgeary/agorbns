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
start <- 5
stop <- 5
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
k.c.stockago = stockago[mirna, experiment]
k.c.lib = 100



## Functions ###################################################################
## I/O functions

#______________
GetSitesXcounts <- function(mirna, experiment) {
  sites_file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
                       mirna, "/",experiment,"/full_site_count_tables/all_sites_",
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
sitesXcounts <- GetSitesXcounts(mirna, experiment)

# Fix column names, with respect to the "A" in front of it.
colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
  colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
    return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
  }
)
# print(sitesXcounts)
# print(setdiff(rownames(sitesXcounts),rownames(site_cols)))
# Get data from input library with flanking counts, the average probability of
# pairing, and the variance for this distribution.



# Assign parameters for Kds and background, subracting 1 from rows due to "None"
# Kd assigned to 1 nM, and the fact that the I and A0 columns are not modeled
# in the experiment, so there's no background term assigned to them.

# Assign the number of kd parameters and background parameters:
num.bgs <- dim(sitesXcounts)[2]-3
# Assign the total site concentration in each experiment, and initialize the
# data column to be used.
# k.c.lib should be 100, for 100 nM.
# Remove the I and A0 columns from the data to be fit to the model. 
data <- sitesXcounts[,2:8]
site_seqs <- as.character(sitesXcounts[-nrow(sitesXcounts),1])


remove <- which(!(rownames(data) %in% c("8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-m8", "6mer-A1",
                 "5mer-A1", "5mer-m2.6", "5mer-m3.7", "5mer-m8", "4mer-A1","4mer-m2.5",
                 "4mer-m3.6", "4mer-m4.7","4mer-m8", "None")))
data[nrow(data),] <- data[nrow(data),] + colSums(data[remove,])
data <- data[-remove,]
print('HI')
num.kds <- dim(data)[1]-1


c.I.tots <- matrix(unlist(rep(data[,1]/sum(data[,1])*k.c.lib,
                  5)), ncol = 5, byrow = FALSE)
r_denomenator <- data[,1] / sum(data[,1])

data <- data[,2:6]
# print(c.I.tots)
colnames(c.I.tots) <- colnames(data)
rownames(c.I.tots) <- rownames(data)
print("hi")
# Initialize starting Kds, which are set to 1/the enrichment of each site type
# in the A12.6 experiment.
r_numerator <- data[, 3] / sum(data[, 3])
kds.init <- ((r_numerator / r_denomenator)[1 : num.kds])^(1/2)
pars.init <- c(-log10(kds.init), -2, -2, -2, -2, -2, log10(k.c.stockago))
print('hi')
trial <- pars.init
trial[trial == Inf] <- 0
trial[is.na(trial)] <- 0
trial[trial == -1*Inf] <- 0

tick <- 0
# Define function of just kds and pars.
c.totals <- c.I.tots
c.agos <- sapply(colnames(data), function(x) {
   as.numeric(unlist(strsplit(x,split="_"))[[1]]) / 100
  })
print("hi")
GetModelFrequencies <- function(pars, out_model = out) {
  # Split up the parameters into the kd and background parameters.
  kds  <- 10^c(pars[1 : num.kds], 1)
  bgs  <- 10^pars[(num.kds + 1) : (num.kds + num.bgs)]
  # print(bgs)
  # print(c.totals[(nrow(c.totals)-10):nrow(c.totals),])
  names(kds) <- rownames(data)
  # Solve for the free Ago concentration in each experiment.

  # Initialize a matrix with the same total concentration of each site type
  # for each experiment.
  c.bounds <- as.matrix(
    sapply((10^pars[length(pars)])*c.agos, function(x) {
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
  if (tick%%1000 == 0) {
    out<<- rbind(out, c(pars, -1*prob_pois))
    print(-1*prob_pois)
    print("AGO")
    print(10^pars[length(pars)])
    print("background ratio")
    print(colSums(c.bgs) / colSums(c.all))
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
    # axis(1,sapply(10^(-3:8), function(x) {x * c(1,2,3,4,5,6,7,8,9)}), labels = FALSE, pos = 1, lwd = 2)
    # axis(2,sapply(10^(0:8), function(x) {x * c(1,2,3,4,5,6,7,8,9)}),labels = FALSE, pos = 10^-3,lwd = 2)
    # axis(1,10^(-3:8), pos = 1, lwd = 0)
    # axis(2,10^(0:8), pos = 10^-3,lwd = 0)


    # sapply(1:ncol(data), function(col) {
    #   x_model <- c.final[, col]
    #   y_data <- data[, col]
    #   x_plot <- x_model/sum(x_model)*sum(y_data)
    #   y_plot <- y_data
    #   points(x_plot, y_plot,col=c(site_cols[rownames(data),],"black"),lwd=2,pch=20)
    # })
        x <- c(40,12.65,4,1.265,0.4)*k.c.stockago/100*1000
    y <- c(1,1,1,1,1)
    sites.norm <- c.totals[,1] / sum(c.totals[,1])

    data.norm <- t(t(data)/colSums(data))

    data.R <- data.norm/(sites.norm)

    model.R <- c.final/(sites.norm)


    xmin <- floor(0.5*min(x))
    xmax <- ceiling(2*max(x))

    ymin <- 0.2
    ymax <- ceiling(max(data.R)*2)
    yextension <- (ymax/ymin)
    xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
    xs <- xs[xs >= xmin & xs <= xmax]
    xmin <- min(xs)
    xmax <- max(xs)
    xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
    ys <- c(sapply(seq(floor(log10(ymin)), ceiling(log10(ymax))), function(x) seq(10)*10^x))
    ys <- ys[ys >= ymin & ys <= ymax]
    ymin <- min(ys)
    ymax <- max(ys)
    yl <- 10^seq(ceiling(log10(min(ys))), floor(log10(max(ys))))

    plot(x, y,log='xy', xlim=c(xmin, xmax*(xmax/xmin)^0.55), ylim=c(ymin, ymax), type="l",
       col="white", axes=FALSE, ann=FALSE)        
    # Generate tickmarks for axis.

    axis(1, at=xl,
         labels=sapply(xl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=ymin, lwd=0, cex.axis=1.5, padj=0.2)
    axis(1, at=xs, labels=FALSE,
         pos=ymin, lwd=2)
    # Label the axis at each order of magnitude.

    axis(2, at=yl,
         labels=sapply(yl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
    axis(2, at=ys, labels=FALSE,
         pos=xmin, lwd=2)

    title(main = mirna, font.main=1, cex.main=1.5, line=-2, adj=0.1)
    title(xlab = "[AGO2-miRNA] (pM)", cex.lab=1.5, line=2, adj=0.3)
    title(ylab = "Enrichment", cex.lab=1.5, line=2)
    legend_names <- rownames(data)
    # centered_sites <- c("11mer-m3.13", "12mer-m3.14",
    #                 "11mer-m4.14", "12mer-m4.15")
    # centered_index <- sapply(centered_sites, function(site){
    #                        which(legend_names == site)
    #                      }
    #                      )
    # legend_names[centered_index] <- "Centered"
    # legend_names <- unique(legend_names)
    legend(x=xmax,y=ymax,legend=legend_names, pch=19, col=site_cols[legend_names, ], cex=1.2, bty="n")
    for (name in rownames(data)) {
      points(x, data.R[name, ], col=site_cols[name, ], pch=19, lwd=3)
      lines(x, model.R[name, ], col=site_cols[name, ], lwd=2)
      
    }


    # dev.off()
      out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/", experiment, "/kds_with_structure/",
                   start, "-", stop, "_", win_left, "-", win_right, "_preliminary_protein_temp.txt")
      write.table(file=out_file, out, sep="\t", quote=FALSE, row.names=FALSE)

    MakeSiteIterationPlot(out,paste0(experiment, start, "-", stop, "_", win_left, "-", win_right,"_preliminary_protein_temp"),mirna)  
 
  }
  tick <<- tick + 1
  return(-1*prob_pois)
}
# trial <- out[nrow(out),-ncol(out)]

# out <- rbind(c(trial, 100000),c(trial, 100000))
# out <- out[(nrow(out)-2):nrow(out),]

# print(out)
colnames(out) <- c(rownames(data)[-nrow(data)],
                   "40", "12.6", "4", "12.6", "0.4", "AGO","-logp")
# print(out)

print("NOW IT'S OPTIMIZING!!")

GetModel <- function(pars, out_model) {
  # Split up the parameters into the kd and background parameters.
  kds  <- 10^c(pars[1 : num.kds], 1)
  bgs  <- 10^pars[(num.kds + 1) : (num.kds + num.bgs)]
  # print(bgs)
  # print(c.totals[(nrow(c.totals)-10):nrow(c.totals),])
  names(kds) <- rownames(data)
  # Solve for the free Ago concentration in each experiment.

  # Initialize a matrix with the same total concentration of each site type
  # for each experiment.
  c.bounds <- as.matrix(
    sapply((10^pars[length(pars)])*c.agos, function(x) {
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

  return(t(t(c.final)*colSums(data)))
}
# print(trial)
trial <- pars_new$par
# # # Solve the first run of the function, and create the output matrix.

par_scale <- sapply(1:length(trial), function(i) {
  pars_temp <- trial
  i_0 <- GetModelFrequencies(pars_temp)
  pars_temp[i] <- pars_temp[i] + 0.1 * pars_temp[i]
  i_1 <- GetModelFrequencies(pars_temp)
  return(abs(i_1 - i_0))


})
print(par_scale)
par_scale[par_scale == 0] <- 1
for (rep in seq(1,1000)) { 
  pars_new <- optim(trial,
                    GetModelFrequencies,
                    out_model = out,
                    method = "Nelder-Mead",
                    control = list(maxit = 20000, parscale = 1/par_scale))
  print(pars_new$value)
  trial <- pars_new$par
}



# Assign output file for the entire sequnce of the optimization and write
# to it.

