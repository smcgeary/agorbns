################################################################################
#GenerateSiteTypeKds.py
################################################################################

# Initial parameters and constants.
args = commandArgs(trailingOnly=TRUE)
mirna = args[1]
# Region within random sequence from which site types orginiates,
# going from position [26 - "start" : 26 + 37 + "stop"]
start = args[2]
stop = args[3]
# Region within the miRNA sequence for which the structural accessibility
# is being measured. [win_left : win_right], where win_left is the 5' most
# position withinthe miRNA, and win_right is the 3' most nucleotide position
#  within the miRNA.
win_left = args[4]
win_right = args[5]
# mirna <- "miR-155"
# start <- -3
# stop <- -3
# win_left <- 1
# win_right <- 15
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
print(site_cols)
k.c.stockago = stockago[mirna,1]
k.c.lib = 100

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

#_____________________
GetSitesXFlanks_counts <- function(mirna,site) {
  sites_file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
                       mirna, "/equilibrium/full_site_count_tables/",
                       site, "_flanking_",
                       start, "-", stop, ".txt")
  sitesXcounts <- read.table(sites_file_name)
  return(sitesXcounts)
}

GetSiteKds <- function(mirna, k.c.stockago, start, stop, win_left, win_right) {
  sites_file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/equilibrium/kds_with_structure/", k.c.stockago, "_",
                   start, "-", stop, "_", win_left, "-", win_right, "_stoch_new.txt")
  site.kds <- t(read.table(sites_file_name,header=FALSE))
  return(as.numeric(site.kds[,ncol(site.kds)]))
}


#__________________
GetSitesXStructures <- function(mirna, condition, start, stop, win_left,
                                win_right) {
  # Generates structural matrix for assigning beta parameters in the model.
  #
  # Args:
  # mirna: The mirna in the experiment
  # condition: The experimental condition being studied.
  # start: The left most position in the random sequence, where a positive
  #   number denotes including constant sequence, e.g. 3, includes 3
  #   nucleotides of constant sequence, -3 removes the first 3 nucleotides
  #   of random sequence.
  # stop: Same as start but for the 3' end of the random sequence.
  # win_left
  # Returns:
  # A list of three items: one containing the counts of each site-by-flank
  # combination, a second with the geometric-mean-normalized probability
  # of being unpaired with the specificied window
  sites_file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
                            mirna, "/equilibrium/site_flank_unpaired_data/",
                            condition, "_", start, "-", stop, "_", win_left,
                            "-", win_right, "_noconstant.txt")

  sitesXcounts <- read.table(sites_file_name, stringsAsFactors=FALSE)
  # Removes the "X" in the column names of the header.
  # (I.e. 'X7mer.m8' becomes '7mer.m8')
  colnames(sitesXcounts) <- sapply(colnames(sitesXcounts), function(name) {
    return(paste0(rev(strsplit(noquote(name), split="X")[[1]])[1]))
    })

  # Converts the "." between site type and designation into a dash.
  # (I.e. '7mer.m8' becomes '7mer-m8')
  colnames(sitesXcounts) <- sapply(colnames(sitesXcounts), function(name) {
    return(paste0(strsplit(name,split="mer\\.")[[1]],collapse="mer-"))
    })

  # Converts the "." between BN and a number into an open parenthesis.
  # (I.e. '7mer-m8bT.4.6.' becomes '7mer-m8bT(4.6.')
  colnames(sitesXcounts) <- sapply(colnames(sitesXcounts), function(name) {
    return(paste0(strsplit(name,split="bT\\.")[[1]],collapse="bT("))
    })
  colnames(sitesXcounts) <- sapply(colnames(sitesXcounts), function(name) {
    return(paste0(strsplit(name,split="bG\\.")[[1]],collapse="bG("))
    })
  colnames(sitesXcounts) <- sapply(colnames(sitesXcounts), function(name) {
    return(paste0(strsplit(name,split="bA\\.")[[1]],collapse="bA("))
    })

  # Converts the "." at the end to an open parenthesis.
  # (I.e. '7mer-m8bT.4.6.' becomes '7mer-m8bT(4.6.')
  colnames(sitesXcounts) <- sapply(colnames(sitesXcounts), function(name) {
    if (substr(name,nchar(name),nchar(name))==".") {
      return(paste0(substr(name,1,nchar(name)-1),")"))
    } else {
      return(name)
    }
  })

  # Generate site & flank columns for output dataframe.
  site <- rep(colnames(sitesXcounts), each=nrow(sitesXcounts))
  flank <- rep(rownames(sitesXcounts), ncol(sitesXcounts))

  # Generate the counts ("n"), mean structural accessibility ("u"), and
  # the variance ("rho") for a particular site-flanking dinucleotide
  # combination.
   n <- c(apply(sitesXcounts,2,function(col){
    sapply(col,function(x){
      if (x == "0") {
        return(0)
      } else {
        col = strsplit(x,split=",")[[1]]
        return(as.numeric(col[1]))
      }
      })
    }))

  u <- c(apply(sitesXcounts,2,function(col){
    sapply(col,function(x){
      if (x == "0") {
        return(NaN)
      } else {
        col = strsplit(x,split=",")[[1]]
        return(as.numeric(col[2]))
      }
      })
    }))

  rho <- c(apply(sitesXcounts,2,function(col){
    sapply(col,function(x){
      if (x == "0") {
        return(NaN)
      } else {
        col = strsplit(x,split=",")[[1]]
        return(as.numeric(col[3]))
      }
      })
    }))
  output <- data.frame(site, flank, n, u, rho)
  return(output)
}


# TODO: do i need this?
#______

GetLinearCoefficients <- function(linear_model, pop=TRUE) {
  coeffs <- summary(linear_model)$coefficients
  B <- coeffs[1, 1]
  site  <- coeffs[which(substr(rownames(coeffs), 1, 4)=="site"), 1]
  names(site) <- sapply(names(site), function(name) {
    return(substr(name, 5, nchar(name)))
  })
  flank <- coeffs[which(substr(rownames(coeffs), 1, 5)=="flank"), 1]
  names(flank) <- sapply(names(flank), function(name) {
    return(substr(name, 6, nchar(name)))
  })
  return(list(B, site, flank))

}

GetFitValues <- function(data, coeffs) {
    B <- coeffs[[1]]
    site <- coeffs[[2]]
    flank <- coeffs[[3]]

    flank.total <- unique(as.character(data$flank))
    site.total <- unique(as.character(data$site))
    site <- c(0, site)
    flank <- c(0, flank)

    names(site)[1] <- setdiff(site.total, names(site))
    names(flank)[1] <- setdiff(flank.total, names(flank))

    model.values <- (B + site[as.character(data$site)] + 
                     flank[as.character(data$flank)])

    model.values <- (1 + exp(-1 * model.values))^(-1)
    return(model.values)
}

GetFitValues <- function(data,coeffs,variance=FALSE) {


    B <- coeffs[[1]]
    site <- coeffs[[2]]
    flank <- coeffs[[3]]
    flank_total <- unique(as.character(data$flank))
    if (variance == TRUE) {
      site_total <- as.character(unique(subset(data, n > 1, select = c(site)))$site)
    } else {
    site_total <- unique(as.character(data$site))

    }
    site <- c(0, site)

    flank <- c(0, flank)
    names(site)[1] <- setdiff(site_total,names(site))
    names(flank)[1] <- setdiff(flank_total,names(flank))
    if (length(site) < length(unique(data$site))) {
      remaining <- setdiff(unique(as.character(data$site)), names(site))
      vars <- rep(0, length(remaining))
      names(vars) <- remaining
      site <- c(site,vars)
    }
    print(site)
    # if (variance == TRUE) {

    # }
    model.values <- (B + site[as.character(data$site)]
      + flank[as.character(data$flank)])
    model.values <- (1 + exp(-1*model.values))^(-1)
    return(model.values)
}


GetBetaParams <- function(u,rho) {
  alpha <- ((1 - u)/ rho - 1 / u) * u^2
  beta <- alpha * (1 / u - 1)

  if (is.infinite(alpha)) {
    alpha <- 1
  }
  if (is.infinite(beta)) {
    beta <- 1
  }

  return(c(alpha,beta))

}

# Modeling Functions
#___________
GetOccupancy <- function(c.freeago, kds, IDs, freq.flankXp, p.exponent) {
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
  x <- as.numeric(colnames(freq.flankXp))

  out <- sapply(1 : nrow(freq.flankXp), function(row) {
    kd <- kds[IDs[row,1]]
    return(c.freeago * (c.freeago + kd * x^(-p.exponent))^(-1))
    })
  return(t(out))
}

#______________
GetFreeResidual <- function(c.freeago, kds, IDs, freq.flankXp, p.exponent,
                            c.ago) {
  occs_I <- GetOccupancy(c.freeago, kds, IDs, freq.flankXp, p.exponent)
  residual <- (c.ago - c.freeago - sum(occs_I*freq.flankXp))^2
  return(residual)
}

#_________
GetFreeAgo <- function(kds, IDs, freq.flankXp, p.exponent, c.ago) {
  c.free <- optimize(GetFreeResidual,
                     c(0, c.ago),
                     kds = kds,
                     IDs = IDs,
                     freq.flankXp = freq.flankXp,
                     p.exponent = p.exponent,
                     c.ago = c.ago)$minimum
  return(c.free)
}

# Output Functions:
CheckMaxDifference <- function(out) {
  row.last <- dim(out)[1]
  row.secondtolast <- row.last-1
  diffs <- sapply(1:dim(out)[2],function(x) {
    if (out[row.last,x]-out[row.secondtolast,x] != 0) {
      return(abs((out[row.last,x]-out[row.secondtolast,x])/out[row.secondtolast,x]))
    } else {
      return(0)
    }
  })
  return(max(diffs))
}




# MAIN #########################################################################

# 1. Get data table:
sitesXcounts <- GetSitesXcounts(mirna)

# Fix column names, with respect to the "A" in front of it.
colnames(sitesXcounts)[2:dim(sitesXcounts)[2]] <- sapply(
  colnames(sitesXcounts)[2:dim(sitesXcounts)[2]], function(x) {
    return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
  }
)
print(sitesXcounts)
# Get data from input library with flanking counts, the average probability of
# pairing, and the variance for this distribution.
I.data <- GetSitesXStructures(mirna, "I", start, stop, win_left, win_right)

# 2. Get data for the five time points 
A.data <- sapply(c("40", "12.6", "4", "1.26", "0.4"), function(condition) {
                   as.numeric(GetSitesXStructures(mirna, condition, start,
                                                  stop, win_left, win_right)$n)
                 })
I.data$rho[which(is.nan(I.data$rho))] <- NA
I.data$rho[which(I.data$rho == Inf)] <- NA
A.data <- rbind(A.data,c(sitesXcounts["None",2:6,drop=FALSE]))
A.data <- matrix(as.numeric(A.data),byrow=FALSE,ncol=5)
colnames(A.data) <- c("40","12.6","4","1.26","0.4")
print(unique(I.data$site))
print(dim(A.data))
print(dim(I.data))
# 3. Generate structural matrix
mean.fit <- glm(u ~ site + flank, quasi(link = "logit"), data=I.data)

var.fit <- glm(rho ~ site + flank, quasi(link = "logit"),
               data=I.data[which(I.data$rho > 0), ], na.action = na.exclude)

mean.coeffs <- GetLinearCoefficients(mean.fit)
var.coeffs <- GetLinearCoefficients(var.fit)

I.u.model <- GetFitValues(I.data, mean.coeffs)
rho_model <- lm(log10(I.data$rho[which(I.data$n > 0 & I.data$rho > 0)]) ~ log10(I.data$u[which(I.data$n > 0 & I.data$rho > 0)]))
rho_factor_1 <- rho_model$coefficient[1]
rho_factor_2 <- rho_model$coefficient[2]
I.rho.model <- 10^(rho_factor_1 + log10(I.u.model)*rho_factor_2)
I.beta.params <- apply(cbind(I.u.model, I.rho.model), 1, function(row) {
  GetBetaParams(row[1], row[2])
  })
# The range of site accesibilities from 0 to 1
unpaired.range <- seq(0.000001,0.99999,length=101)
sites.unpaired.dist <- t(apply(I.beta.params, 2, function(col) {
  dbeta(unpaired.range, col[1], col[2])
}))

rownames(sites.unpaired.dist) <- c()

# Generates the input matrix that includes the no site distribution
I.dist <- rbind(sites.unpaired.dist, dbeta(unpaired.range, 1, 1))

# Takes the input distribution and divides each row by the sum of that row
I.dist.norm <- t(apply(cbind(I.dist,c(I.data$n,sitesXcounts["None", 1])), 1, function(row){
  data <- row[1:(length(row)-1)]
  tot <- row[length(row)]
  return(data / sum(data) * tot)
}))

c.I.matrix <- I.dist.norm / sum(I.dist.norm) * 100
c.I.all <- rowSums(c.I.matrix)

I.IDs <- rbind(cbind(as.character(I.data$site), as.character(I.data$flank)),c("None","None"))
colnames(c.I.matrix) <- unpaired.range

temp_kds <- rep(1,nrow(sitesXcounts))
names(temp_kds) <- rownames(sitesXcounts)



# Assign parameters for Kds and background, subracting 1 from rows due to "None"
# Kd assigned to 1 nM, and the fact that the I and A0 columns are not modeled
# in the experiment, so there's no background term assigned to them.

# Assign the number of kd parameters and background parameters:
num.kds <- dim(sitesXcounts)[1]-1
num.bgs <- dim(sitesXcounts)[2]-2
# Assign the total site concentration in each experiment, and initialize the
# data column to be used.
# k.c.lib should be 100, for 100 nM.

# Remove the I and A0 columns from the data to be fit to the model. 
data <- A.data

# Initialize starting Kds, which are set to 1/the enrichment of each site type
# in the A12.6 experiment.
pars.init <- GetSiteKds(mirna, k.c.stockago, start, stop, win_left, win_right)

trial <- pars.init[-length(pars.init)]
val <- pars.init[length(pars.init)]
tick <- 0
# Define function of just kds and pars.
print(dim(c.I.matrix))

GetModelFrequencies <- function(pars, c.I.model = c.I.matrix, I.IDs.model = I.IDs, out_model = out) {
  # Split up the parameters into the kd and background parameters.
  kds  <- 10^c(pars[1 : num.kds], 1)
  bgs  <- 10^pars[(num.kds + 1) : (num.kds + num.bgs)]
  k    <- pars[num.kds + num.bgs + 1]
  bg.a <- 10^pars[num.kds + num.bgs + 2]
  bg.b <- 10^pars[num.kds + num.bgs + 3]
  I.beta.nosite <- dbeta(unpaired.range, 1 + bg.a, 1 + bg.b)
  I.beta.nosite <- I.beta.nosite / sum(I.beta.nosite) * rowSums(c.I.model)[nrow(c.I.model)]
  c.I.model[nrow(c.I.model), ] <- I.beta.nosite
  names(kds) <- rownames(sitesXcounts)
  # Solve for the free Ago concentration in each experiment.
  c.freeagos <<- sapply(colnames(data), function(x) {
    c.ago <- as.numeric(x) / 100 * k.c.stockago*0.6
      return(GetFreeAgo(kds, I.IDs, c.I.model, k, c.ago))
    }
  )

  # Initialize a matrix with the same total concentration of each site type
  # for each experiment.
  c.totals <- as.matrix(
    sapply(colnames(data), function(x) {
             return(rowSums(c.I.model))
           }
           ))
  
  # Use the free Ago concentrations to get the amount of each complex bound
  # to Ago.
  c.bounds <- as.matrix(
    sapply(c.freeagos, function(x) {
             return(rowSums((GetOccupancy(x, kds, I.IDs.model, c.I.model, k) * c.I.model)))
           }
           ))

  # Get the amount of background binding by subtracting the bound from the
  # total sites in each experiment, normalizing. Must transpose to multiply
  # each column.
  c.frees <- c.totals - c.bounds

  c.bgs <- t(t(c.frees) * bgs / colSums(c.frees))

  c.all <- c.bounds + c.bgs
  c.final <- data.frame(t(t(c.all) / colSums(c.all)))
  colnames(c.final) <- colnames(data)
  sites_unique <- as.character(rownames(sitesXcounts)[1:(nrow(sitesXcounts)-1)])
  dinucs <- unique(substr(as.character(I.IDs.model[1:(nrow(I.IDs.model)-1),2]),2,3))
  paired_types <- expand.grid(dinucs,sites_unique)
  c.final_collapse <- t(apply(paired_types,1,function(type) {
        sums <- colSums(c.final[which(I.IDs.model[,1] == type[2] & substr(I.IDs.model[,2],2,3) == type[1]),])
        return(sums)
        }))
  Centered_ind <- which(paired_types[, 2] == "Centered")
  c.final_collapse <- rbind(c.final_collapse[1:(min(Centered_ind) - 1), ],
                            colSums(c.final_collapse[Centered_ind, ]),
                            c.final_collapse[(max(Centered_ind) + 1):nrow(c.final_collapse), ],
                            c.final[nrow(c.final), , drop=FALSE])

    data_collapse <- t(apply(paired_types,1,function(type) {
        sums <- colSums(data[which(I.IDs.model[,1] == type[2] & substr(I.IDs.model[,2],2,3) == type[1]),])
        return(sums)
        }))
    data_collapse <- rbind(data_collapse[1:(min(Centered_ind) - 1), ],
                            colSums(data_collapse[Centered_ind, ]),
                            data_collapse[(max(Centered_ind) + 1):nrow(data_collapse), ],
                            data[nrow(data), , drop=FALSE])
  IDs_collapse <- matrix(apply(paired_types,1, function(type){
      c(as.character(type[2]),as.character(type[1]))
      }), byrow=TRUE,ncol=2)
  IDs_collapse <- rbind(IDs_collapse[1:(min(Centered_ind) - 1), ],
                            c("Centered", "NA"),
                            IDs_collapse[(max(Centered_ind) + 1):nrow(IDs_collapse), ])

  prob_pois <- sum(sapply(colnames(c.final_collapse), function(col_name) {
    dmultinom(x=data_collapse[which(c.final_collapse[, col_name] > 0),col_name], prob=c.final_collapse[which(c.final_collapse[, col_name] > 0),col_name],log = TRUE)
    }))
  print(-1*prob_pois)
  if (tick%%10 == 0) {
    out<<- rbind(out, c(pars, -1*prob_pois))
    setEPS()
    postscript(file=paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/figures/kds/",
    mirna,"/", "structure_fit_", start, "-", stop, "_", win_left, "-", win_right, "_round2_new.eps"))
    plot(c(0.0001,0.0002),c(0.0001,0.0002),col="white",xlim=c(0.00000005,.02),ylim=c(0.00000005,.02), log='xy')
    for (name in colnames(c.final)) {
      x_model <- c.final_collapse[, name]
      y_data <- data_collapse[, name]
      x_plot <- x_model/sum(x_model)
      y_plot <- y_data/sum(y_data)
      points(x_plot, y_plot,col=c(site_cols[IDs_collapse[,1],],"black"),lwd=2,pch=20)
    }
    dev.off()
      out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/equilibrium/kds_with_structure/", k.c.stockago, "_",
                   start, "-", stop, "_", win_left, "-", win_right, "_round2_new.txt")
      write.table(file=out_file, out, sep="\t", quote=FALSE, row.names=FALSE)

    MakeSiteIterationPlot(out,paste0("structure", start, "-", stop, "_", win_left, "-", win_right,"_round2_new"))  
 
  }
  tick <<- tick + 1
  return(-1*prob_pois)
}

out <- rbind(pars.init,pars.init)
colnames(out) <- c(rownames(sitesXcounts)[-nrow(sitesXcounts)],
                   "40", "12.6", "4", "12.6", "0.4",
                   "unpaired_k", "bg_alpha", "bg_beta","-logp")
print(out)

print("NOW IT'S OPTIMIZING!!")
print(pars.init)

# # # Solve the first run of the function, and create the output matrix.
pars_new <- optim(trial,
                  GetModelFrequencies,
                  c.I.model = c.I.matrix,
                  I.IDs.model = I.IDs,
                  out_model = out,
                  control = list(maxit = 200000))


# Assign output file for the entire sequnce of the optimization and write
# to it.

