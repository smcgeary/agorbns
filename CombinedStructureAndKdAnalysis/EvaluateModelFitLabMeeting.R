source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")


experiment <- "equilibrium"
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

k.c.stockago <- stockago[mirna, experiment]
k.c.lib = 100





GetSitesXcounts <- function(mirna) {
  sites_file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
                       mirna, "/equilibrium/full_site_count_tables/all_sites_",
                       start, "-", stop, ".txt")
  sitesXcounts <- read.table(sites_file_name)
  return(sitesXcounts)
}

# # #_____________________
# # GetSitesXFlanks_counts <- function(mirna,site) {
# #   sites_file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
# #                        mirna, "/equilibrium/full_site_count_tables/",
# #                        site, "_flanking_",
# #                        start, "-", stop, ".txt")
# #   sitesXcounts <- read.table(sites_file_name)
# #   return(sitesXcounts)
# # }

# # #__________________
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
  print(rownames(sitesXcounts)[1:10])
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


# # # TODO: do i need this?
# # #______

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

# # GetEquilParameters <- function(experiment,mirna) {
# #   file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
# #                    "/equilibrium/kds_with_structure/",
# #                    start, "-", stop, "_1-15_Basic_singlebg.txt")
# #   return(read.table(file_name,header = TRUE))
# # }
# # # Modeling Functions
# # #___________
GetOccupancy <- function(c.freeago, kds, x, p.exponent) {
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
    out <- c.freeago * (c.freeago + kds %*% t(x^(-p.exponent)))^(-1)
  return(out)
}

# # #______________
GetFreeResidual <- function(c.freeago, kds, freq.flankXp, p.exponent,
                            c.ago) {
  occs_I <- GetOccupancy(c.freeago, kds, as.numeric(colnames(freq.flankXp)), p.exponent)
  residual <- (c.ago - c.freeago - sum(occs_I*freq.flankXp))
  return(residual)
}

# # #_________
GetFreeAgo <- function(kds, freq.flankXp, p.exponent, c.ago) {
  c.free <- uniroot(GetFreeResidual,
                     c(0, c.ago),
                     kds = kds,
                     freq.flankXp = freq.flankXp,
                     p.exponent = p.exponent,
                     c.ago = c.ago)$root
  return(c.free)
}

# # # Output Functions:
# # CheckMaxDifference <- function(out) {
# #   row.last <- dim(out)[1]
# #   row.secondtolast <- row.last-1
# #   diffs <- sapply(1:dim(out)[2],function(x) {
# #     if (out[row.last,x]-out[row.secondtolast,x] != 0) {
# #       return(abs((out[row.last,x]-out[row.secondtolast,x])/out[row.secondtolast,x]))
# #     } else {
# #       return(0)
# #     }
# #   })
# #   return(max(diffs))
# # }




# # # MAIN #########################################################################

# # # 1. Get data table:
# sitesXcounts <- GetSitesXcounts(mirna)
# # sites_norm <- sitesXcounts/c(sitesXcounts["None",])
# # sites_norm_norm <- sites_norm[,-1]/(sites_norm[, 2])
# # print(round(sites_norm_norm, 4))
# # # # Fix column names, with respect to the "A" in front of it.
# colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
#   colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
#     return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
#   }
# )
# print(sitesXcounts)
# # # # Get data from input library with flanking counts, the average probability of
# # # # pairing, and the variance for this distribution.
# I.data <- GetSitesXStructures(mirna, "I", start, stop, win_left, win_right)
# inds_for_counts <- sapply(I.data$site, function(s) {
#   which(rownames(sitesXcounts) == s)
# })
# I.data <- I.data[order(inds_for_counts),]
# # # print(dim(I.data))
# # # print("hi")
# # # # 2. Get data for the five time points 
# A.data <- sapply(c("40", "12.6", "4", "1.26", "0.4"), function(condition) {
#                    as.numeric(GetSitesXStructures(mirna, condition, start,
#                                                   stop, win_left, win_right)$n)
#                  })
# A.data <- A.data[order(inds_for_counts),]
# # print(dim(A.data))
# # # print((A.data))
# I.data$rho[which(is.nan(I.data$rho))] <- NA
# I.data$rho[which(I.data$rho == Inf)] <- NA
# A.data <- rbind(A.data,c(sitesXcounts["None",2:6,drop=FALSE]))
# A.data <- matrix(as.numeric(A.data),byrow=FALSE,ncol=5)
# colnames(A.data) <- c("40","12.6","4","1.26","0.4")
# # print(unique(I.data$site))
# # # print(dim(A.data))
# # # print(dim(I.data))
# # # # 3. Generate structural matrix
# mean.fit <- glm(u ~ site + flank, quasi(link = "logit"), data=I.data)

# var.fit <- glm(rho ~ site + flank, quasi(link = "logit"),
#                data=I.data[which(I.data$rho > 0), ], na.action = na.exclude)

# mean.coeffs <- GetLinearCoefficients(mean.fit)
# var.coeffs <- GetLinearCoefficients(var.fit)

# I.u.model <- GetFitValues(I.data, mean.coeffs)
# rho_model <- lm(log10(I.data$rho[which(I.data$n > 0 & I.data$rho > 0)]) ~ log10(I.data$u[which(I.data$n > 0 & I.data$rho > 0)]))
# rho_factor_1 <- rho_model$coefficient[1]
# rho_factor_2 <- rho_model$coefficient[2]
# I.rho.model <- 10^(rho_factor_1 + log10(I.u.model)*rho_factor_2)
# I.beta.params <- apply(cbind(I.u.model, I.rho.model), 1, function(row) {
#   GetBetaParams(row[1], row[2])
#   })
# # # The range of site accesibilities from 0 to 1
# unpaired.range <- seq(0.000001,0.99999,length=21)
# sites.unpaired.dist <- t(apply(I.beta.params, 2, function(col) {
#   dbeta(unpaired.range, col[1], col[2])
# }))
# num.kds <- nrow(sitesXcounts)
# rownames(sites.unpaired.dist) <- c()

# # # # Generates the input matrix that includes the no site distribution
# I.dist <- rbind(sites.unpaired.dist, dbeta(unpaired.range, 1, 1))

# # # # Takes the input distribution and divides each row by the sum of that row
# I.dist.norm <- t(apply(cbind(I.dist,c(I.data$n,sitesXcounts["None", 2])), 1, function(row){
#   data <- row[1:(length(row)-1)]
#   tot <- row[length(row)]
#   return(data / sum(data) * tot)
# }))

# c.I.matrix <- I.dist.norm / sum(I.dist.norm) * 100
# c.I.all <- rowSums(c.I.matrix)

# I.IDs <- rbind(cbind(as.character(I.data$site), as.character(I.data$flank)),c("None","None"))
# colnames(c.I.matrix) <- unpaired.range

# # temp_kds <- rep(1,nrow(sitesXcounts))
# # # names(temp_kds) <- rownames(sitesXcounts)


# # print(413)
# # # # Assign parameters for Kds and background, subracting 1 from rows due to "None"
# # # # Kd assigned to 1 nM, and the fact that the I and A0 columns are not modeled
# # # # in the experiment, so there's no background term assigned to them.

# # # # Assign the number of kd parameters and background parameters:
# # num.kds <- dim(sitesXcounts)[1]
# # num.bgs <- 1
# # # # Assign the total site concentration in each experiment, and initialize the
# # # # data column to be used.
# # # # k.c.lib should be 100, for 100 nM.

# # # # Remove the I and A0 columns from the data to be fit to the model. 
# data <- A.data

GetEquilParameters <- function(experiment,mirna) {
  file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/equilibrium/kds_with_structure/",
                   start, "-", stop, "_1-15_Basic_singlebg.txt")
  return(read.table(file_name,header = TRUE))
}

GetStructureParameters <- function(experiment, mirna) {
      file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                     "/equilibrium/kds_with_structure/2.8_",
                       start, "-", stop, "_", win_left, "-", win_right, "_stoch_new_LABMEETING.txt")
  return(read.table(file_name,header = TRUE))

}


# # # Initialize starting Kds, which are set to 1/the enrichment of each site type
# # # in the A12.6 experiment.
# # r_numerator <- sitesXcounts[, 2] / sum(sitesXcounts[, 2])
# # r_denomenator <- sitesXcounts["I"] / sum(sitesXcounts["I"])
# # kds.init <- (r_numerator / r_denomenator)[1 : num.kds, 1]
# # # pars.init <- c(-2*log10(kds.init), -1, -1, -1, -1, -1, 0.5, 0, 0)
pars.init <- GetStructureParameters("equilibrium",mirna)
pars.null <- GetEquilParameters("equilibrium", mirna)

pars.null <- unlist(c(pars.null[nrow(pars.null),1 : (num.kds+1)], 0, 1, 1))
# # print(pars.init)
trial <- unlist(c(pars.init[nrow(pars.init), -ncol(pars.init)]))
names(trial) <- c(rownames(sitesXcounts), "bg","k", "bg1", "bg2")
names(pars.null) <- names(trial)
# tick <- 0
# # Define function of just kds and pars.
# print("hi")

GetModelFrequencies <- function(pars, out_model, c.I.model = c.I.matrix, I.IDs.model = I.IDs) {

  # Split up the parameters into the kd and background parameters.
  kds  <- Logistic(pars[I.IDs.model[,1]], 10)
  # c.I.model <- c.I.model[1:10,]
  # kds <- kds[1:10]
  # I.IDs.model <- I.IDs.model[1:10,]
  # print(kds)
  bgs  <- rep(10^pars[num.kds+1],5)
  k    <- pars[num.kds+2]
  bg.a <- 10^pars[num.kds + 3]
  bg.b <- 10^pars[num.kds + 4]
  # print(k)
  # print(bg.a)
  # print(bg.b)
  I.beta.nosite <- dbeta(unpaired.range, 1 + bg.a, 1 + bg.b)
  I.beta.nosite <- I.beta.nosite / sum(I.beta.nosite) * rowSums(c.I.model)[nrow(c.I.model)]
  c.I.model[nrow(c.I.model), ] <- I.beta.nosite
  # Solve for the free Ago concentration in each experiment.
  c.freeagos <- sapply(colnames(data), function(x) {
    c.ago <- as.numeric(x) / 100 * k.c.stockago
      return(GetFreeAgo(kds, c.I.model, k, c.ago))
    }
  )
  # Initialize a matrix with the same total concentration of each site type
  # for each experiment.
  c.totals <<- as.matrix(
    sapply(colnames(data), function(x) {
             return(rowSums(c.I.model))
           }
           ))

  # Use the free Ago concentrations to get the amount of each complex bound
  # to Ago.
  print("hello")
  c.bounds <<- as.matrix(
    sapply(c.freeagos, function(x) {
             return(rowSums((GetOccupancy(x, kds, as.numeric(colnames(c.I.model)), k) * c.I.model)))
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
  c.final_collapse <<- rbind(t(apply(paired_types,1,function(type) {
        sums <- colSums(c.final[which(I.IDs.model[,1] == type[2] & substr(I.IDs.model[,2],2,3) == type[1]),])
        return(sums)
        })), c.final[nrow(c.final), ])
  # Centered_ind <- which(paired_types[, 2] == "Centered")
  # c.final_collapse <- rbind(c.final_collapse[1:(min(Centered_ind) - 1), ],
  #                           colSums(c.final_collapse[Centered_ind, ]),
  #                           c.final_collapse[(max(Centered_ind) + 1):nrow(c.final_collapse), ],
  #                           c.final[nrow(c.final), , drop=FALSE])

    data_collapse <<- rbind(t(apply(paired_types,1,function(type) {
        sums <- colSums(data[which(I.IDs.model[,1] == type[2] & substr(I.IDs.model[,2],2,3) == type[1]),])
        return(sums)
        })), data[nrow(data), ])
    # data_collapse <- rbind(data_collapse[1:(min(Centered_ind) - 1), ],
    #                         colSums(data_collapse[Centered_ind, ]),
    #                         data_collapse[(max(Centered_ind) + 1):nrow(data_collapse), ],
    #                         data[nrow(data), , drop=FALSE])
  IDs_collapse <<- rbind(matrix(apply(paired_types,1, function(type){
      c(as.character(type[2]),as.character(type[1]))
      }), byrow=TRUE,ncol=2), c("None", "NA"))
  # IDs_collapse <- rbind(IDs_collapse[1:(min(Centered_ind) - 1), ],
  #                           c("Centered", "NA"),
  #                           IDs_collapse[(max(Centered_ind) + 1):nrow(IDs_collapse), ])

  prob_pois <- sum(sapply(colnames(c.final_collapse), function(col_name) {
    dmultinom(x=data[which(c.final[, col_name] > 0),col_name], prob=c.final[which(c.final[, col_name] > 0),col_name],log = TRUE)
    }))
  # print(-1*prob_pois)

    # setEPS()
    # postscript(file=paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/figures/kds/",
    # mirna,"/", "structure_fit_", start, "-", stop, "_", win_left, "-", win_right, "_stoch_new_LABMEETING_fit.eps"))
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
    segments(0.01,0.01,10000000,10000000, lty = 2)

    axis(1,sapply(10^(-3:8), function(x) {x * c(1,2,3,4,5,6,7,8,9)}), labels = FALSE, pos = 1, lwd = 2)
    axis(2,sapply(10^(0:8), function(x) {x * c(1,2,3,4,5,6,7,8,9)}),labels = FALSE, pos = 10^-3,lwd = 2)
    axis(1,10^(-3:8), pos = 1, lwd = 0)
    axis(2,10^(0:8), pos = 10^-3,lwd = 0)

    for (name in colnames(c.final)) {
      # print(c.final[1:10,name])
      # print(data[1:10,name])
      x_model <- c.final_collapse[, name]
      # print(rownames(x_model))
      y_data <- data_collapse[, name]
      x_plot <- x_model/sum(x_model)*sum(y_data)
      y_plot <- y_data
      # print(c(site_cols[IDs_collapse[,1],],"black"))
      points(x_plot, y_plot,col=site_cols[IDs_collapse[,1],],lwd=2,pch=20)
    }

  
  tick <<- tick + 1
  return(-1*prob_pois)
}


SimulateModel <- function(pars, c.I.model = c.I.matrix, I.IDs.model = I.IDs) {
  # Split up the parameters into the kd and background parameters.

  # Split up the parameters into the kd and background parameters.
  kds  <- as.vector(Logistic(pars[I.IDs.model[,1]], 10), mode = 'numeric')
  # print(kds)
  # c.I.model <- c.I.model[1:10,]
  # kds <- kds[1:10]
  # I.IDs.model <- I.IDs.model[1:10,]
  # print(kds)
  bgs  <- rep(10^pars[num.kds+1],5)
  k    <- pars[num.kds+2]
  bg.a <- 10^pars[num.kds + 3]
  bg.b <- 10^pars[num.kds + 4]
  print(k)
  print(bg.a)
  print(bg.b)
  I.beta.nosite <- dbeta(unpaired.range, 1 + bg.a, 1 + bg.b)
  I.beta.nosite <- I.beta.nosite / sum(I.beta.nosite) * rowSums(c.I.model)[nrow(c.I.model)]
  c.I.model[nrow(c.I.model), ] <- I.beta.nosite
  # Solve for the free Ago concentration in each experiment.
  c.freeagos <- sapply(colnames(data), function(x) {
    c.ago <- as.numeric(x) / 100 * k.c.stockago
      return(GetFreeAgo(kds, c.I.model, k, c.ago))
    }
  )
  # Initialize a matrix with the same total concentration of each site type
  # for each experiment.
  c.totals <<- as.matrix(
    sapply(colnames(data), function(x) {
             return(rowSums(c.I.model))
           }
           ))

  # Use the free Ago concentrations to get the amount of each complex bound
  # to Ago.
  # print("hello")
  c.bounds <<- as.matrix(
    sapply(c.freeagos, function(x) {
             return(rowSums((GetOccupancy(x, kds, as.numeric(colnames(c.I.model)), k) * c.I.model)))
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
  c.final_collapse <<- rbind(t(apply(paired_types,1,function(type) {
        sums <- colSums(c.final[which(I.IDs.model[,1] == type[2] & substr(I.IDs.model[,2],2,3) == type[1]),])
        return(sums)
        })), c.final[nrow(c.final), ])
  # Centered_ind <- which(paired_types[, 2] == "Centered")
  # c.final_collapse <- rbind(c.final_collapse[1:(min(Centered_ind) - 1), ],
  #                           colSums(c.final_collapse[Centered_ind, ]),
  #                           c.final_collapse[(max(Centered_ind) + 1):nrow(c.final_collapse), ],
  #                           c.final[nrow(c.final), , drop=FALSE])

    data_collapse <<- rbind(t(apply(paired_types,1,function(type) {
        sums <- colSums(data[which(I.IDs.model[,1] == type[2] & substr(I.IDs.model[,2],2,3) == type[1]),])
        return(sums)
        })), data[nrow(data), ])
    # data_collapse <- rbind(data_collapse[1:(min(Centered_ind) - 1), ],
    #                         colSums(data_collapse[Centered_ind, ]),
    #                         data_collapse[(max(Centered_ind) + 1):nrow(data_collapse), ],
    #                         data[nrow(data), , drop=FALSE])
  IDs_collapse <<- rbind(matrix(apply(paired_types,1, function(type){
      c(as.character(type[2]),as.character(type[1]))
      }), byrow=TRUE,ncol=2), c("None", "NA"))
  # IDs_collapse <- rbind(IDs_collapse[1:(min(Centered_ind) - 1), ],
  #                           c("Centered", "NA"),
  #                           IDs_collapse[(max(Centered_ind) + 1):nrow(IDs_collapse), ])

  out_final <- apply(rbind(c.final,colSums(data)), 2, function(col) {
    return(col[1:(length(col)-1)]*col[length(col)])
  })
  return(list(out_final, c.final_collapse, data_collapse, IDs_collapse))
}

ExtractSubsample <- function(pars,site,condition) {
  pars <- as.numeric(pars[-nrow(pars),ncol(pars)])
  print(pars)
  model <- SimulateModel(pars)[[1]]
  extract_model <- model[which(I.IDs[-nrow(data), 1] == as.character(site)), which(colnames(model) == condition)]
  extract_data <- data[which(I.IDs[-nrow(data), 1] == as.character(site)), which(colnames(data) == condition)]
  out <- cbind(extract_model, extract_data)
  rownames(out) <- I.IDs[which(I.IDs[-nrow(data), 1] == as.character(site)), 2]
  return(out)
}



GetLeftFlanks <- function(extract_data,col) {
  flanks_all <- substr(rownames(extract_data),1,2)
  out <- sapply(unique(flanks_all), function(flank) {
    return(sum(extract_data[which(flanks_all == flank),col]))
  })
  names(out) <- unique(flanks_all)
  return(out)
}

GetRightFlanks <- function(extract_data,col) {
  flanks_all <- substr(rownames(extract_data),3,4)
  out <- sapply(unique(flanks_all), function(flank) {
    return(sum(extract_data[which(flanks_all == flank),col]))
  })
  names(out) <- unique(flanks_all)
  return(out)
}

CheckFlanks <- function(site, condition) {
  extract_data <- data[which(I.IDs[-nrow(data), 1] %in% site), which(colnames(data) == condition)]
  extract_input <- I.data$n[which(I.IDs[-nrow(data), 1] %in% site)]
  flanks_right <- substr(I.IDs[which(I.IDs[-nrow(data), 1] == as.character(site)), 2], 3, 4)
  flanks_left <- substr(I.IDs[which(I.IDs[-nrow(data), 1] == as.character(site)), 2], 1, 2)
  sums_right <- sapply(unique(flanks_right), function(flank) {
    return(sum(extract_data[which(flanks_right == flank)])/sum(extract_data))
  })
  sums_left <- sapply(unique(flanks_right), function(flank) {
    return(sum(extract_data[which(flanks_left == flank)])/sum(extract_data))
  })

  sums_right_input <- sapply(unique(flanks_right), function(flank) {
    return(sum(extract_input[which(flanks_right == flank)])/sum(extract_input))
  })
  sums_left_input <- sapply(unique(flanks_right), function(flank) {
    return(sum(extract_input[which(flanks_left == flank)])/sum(extract_input))
  })
  out <- cbind(sums_left/sums_left_input, sums_right/sums_right_input)
  names(out) <- unique(flanks_right)
  cols <- c("blue", "purple", "red", "green")
  names(cols) <- c("A", "C", "G", "T")
  out[is.na(out)] <- 0
  plot(out[, 1], out[, 2], col = cols[substr(unique(flanks_right),1,1)], cex = 4, pch = 19, main = paste0(site,"\n",condition))
  points(out[, 1], out[, 2], col = cols[substr(unique(flanks_right),2,2)], cex = 2, pch = 19)

  legend("topleft", legend = unique(flanks_right), col = rainbow(16, v = 0.6)[inds], cex = 1, pch = 19)
  return(out)
}

CheckFlanksModel <- function(pars,site, condition) {
  pars <- as.numeric(pars[-nrow(pars),ncol(pars)])
  model <- SimulateModel(pars)[[1]]
  extract_data <- model[which(I.IDs[-nrow(data), 1] %in% site), which(colnames(data) == condition)]
  extract_input <- I.data$n[which(I.IDs[-nrow(data), 1] %in% site)]
  flanks_right <- substr(I.IDs[which(I.IDs[-nrow(data), 1] == as.character(site)), 2], 3, 4)
  flanks_left <- substr(I.IDs[which(I.IDs[-nrow(data), 1] == as.character(site)), 2], 1, 2)
  sums_right <- sapply(unique(flanks_right), function(flank) {
    return(sum(extract_data[which(flanks_right == flank)])/sum(extract_data))
  })
  sums_left <- sapply(unique(flanks_right), function(flank) {
    return(sum(extract_data[which(flanks_left == flank)])/sum(extract_data))
  })

  sums_right_input <- sapply(unique(flanks_right), function(flank) {
    return(sum(extract_input[which(flanks_right == flank)])/sum(extract_input))
  })
  sums_left_input <- sapply(unique(flanks_right), function(flank) {
    return(sum(extract_input[which(flanks_left == flank)])/sum(extract_input))
  })
  out <- cbind(sums_left/sums_left_input, sums_right/sums_right_input)
  names(out) <- unique(flanks_right)
  cols <- c("blue", "purple", "red", "green")
  names(cols) <- c("A", "C", "G", "T")
  inds <- which(out[,1] > 0 & out[, 2] > 0)
  plot(out[inds, 1], out[inds, 2], col = cols[substr(unique(flanks_right),1,1)][inds], cex = 4, pch = 19, main = paste0(site,"\n",condition))
  points(out[inds, 1], out[inds, 2], col = cols[substr(unique(flanks_right),2,2)][inds], cex = 2, pch = 19)

  legend("topleft", legend = unique(flanks_right)[inds], col = rainbow(16, v = 0.6)[inds], cex = 1, pch = 19)
  return(out)
}

CheckFlanksSingle <- function(site, condition,left,right) {
  left_ind <- c(1,2)[left]
  right_ind <- c(3,4)[right]
  extract_data <- data[which(I.IDs[-nrow(data), 1] == as.character(site)), which(colnames(data) == condition)]
  extract_input <- I.data$n[which(I.IDs[-nrow(data), 1] == as.character(site))]
  flanks_right <- substr(I.IDs[which(I.IDs[-nrow(data), 1] == as.character(site)), 2], right_ind, right_ind)
  flanks_left <- substr(I.IDs[which(I.IDs[-nrow(data), 1] == as.character(site)), 2], left_ind, left_ind)
  sums_right <- sapply(unique(flanks_right), function(flank) {
    return(sum(extract_data[which(flanks_right == flank)])/sum(extract_data))
  })
  sums_left <- sapply(unique(flanks_right), function(flank) {
    return(sum(extract_data[which(flanks_left == flank)])/sum(extract_data))
  })

  sums_right_input <- sapply(unique(flanks_right), function(flank) {
    return(sum(extract_input[which(flanks_right == flank)])/sum(extract_input))
  })
  sums_left_input <- sapply(unique(flanks_right), function(flank) {
    return(sum(extract_input[which(flanks_left == flank)])/sum(extract_input))
  })
  out <- cbind(sums_left/sums_left_input, sums_right/sums_right_input)
  names(out) <- unique(flanks_right)
  cols <- c("blue", "purple", "red", "green")
  names(cols) <- c("A", "C", "G", "T")
  inds <- which(out[,1] > 0 & out[, 2] > 0)
  plot(out[inds, 1], out[inds, 2], col = cols[unique(flanks_right)][inds], cex = 4, pch = 19)

  legend("topleft", legend = unique(right)[inds], col = rainbow(16, v = 0.6)[inds], cex = 1, pch = 19)
  return(out)
}


TestModel <- function(pars) {
  print(pars)
  model <- SimulateModel(pars)[[1]]
  out <- matrix(sapply(unique(I.data$site), function(site_type) {
                         sapply(seq(1, ncol(model)), function(col_num) {
                          inds <- which(I.IDs[-nrow(data), 1] == as.character(site_type) & I.data$n > 0 & data[-nrow(data), col_num] > 0)
                           cor(log10(model[inds, col_num]/I.data$n[inds]),
                               log10(data[inds, col_num]/I.data$n[inds]))
                            })
                          }), byrow=TRUE, ncol=5)
  rownames(out) <- unique(I.data$site)
  colnames(out) <- colnames(data)
  return(out)
}

PrintModel <- function(pars) {
    c.final <- SimulateModel(pars)[[1]]
    data_sites <- t(sapply(rownames(sitesXcounts), function(site_name) {
    colSums(data[which(I.IDs[,1] == site_name),,drop=FALSE])/colSums(data)
    }))

  model_sites <- t(sapply(rownames(sitesXcounts), function(site_name) {
    colSums(c.final[which(I.IDs[,1] == site_name),,drop=FALSE])/colSums(c.final)
    }))

  data.R <- apply(data_sites,2,function(col){
    col * (sitesXcounts[,"I"]/sum(sitesXcounts[,"I"]))^(-1)
    })
  model.R <- apply(model_sites,2,function(col){
    col * (sitesXcounts[,"I"]/sum(sitesXcounts[,"I"]))^(-1)
    })

    x <- c(40,12.65,4,1.265,0.4)*k.c.stockago/100*1000
    y <- c(1,1,1,1,1)

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
    # dev.set(5)
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
    legend_names <- rownames(sitesXcounts)
    legend(x=xmax,y=ymax,legend=legend_names, pch=19, col=site_cols[legend_names, ], cex=1.2, bty="n")
    for (name in rownames(sitesXcounts)) {
      cs <- col2rgb(site_cols[name, ])/255
      points(x, data.R[name, ], col=site_cols[name, ], pch=19, lwd=3)
      lines(x, model.R[name, ], col=rgb(cs["red", ], cs["green", ], cs["blue", ], alpha=1), lwd=2)
      }
}

RandomModel <- function(pars) {
  model <- SimulateModel(pars)[[1]]
  out <- matrix(sapply(unique(I.data$site), function(site_type) {

                         sapply(seq(1, ncol(model)), function(col_num) {
                          input_dist <- subset(I.data,site == site_type & I.data$n > 0, select = c(n))[,1]
                          # print("input_dist")
                          # print(input_dist)
                          distribution <- model[which(I.IDs[-nrow(data), 1] == as.character(site_type) & I.data$n > 0), col_num]
                          R <- distribution / input_dist
                          # print(R)
                          input_sim <- R * (rmultinom(1, size= sum(input_dist),prob = input_dist)+1)
                          # print(input_sim)
                          random <- rmultinom(1,size = sum(distribution),prob = input_sim)
                          # print(cbind(distribution,random))
                          print(cor(log10(distribution[which(random >0)]),log10(random[which(random > 0)])))
                          out <- cor(log10((distribution/input_dist)[which(random >0)]),log10((random/input_dist)[which(random > 0)]))
                          print(out)
                            if (is.na(out)) {
                              plot(log10(random),log10(distribution))
                            }
                           return(cor(log10((distribution/input_dist)[which(random >0)]),log10((random/input_dist)[which(random > 0)])))
                            })
                          }), byrow=TRUE, ncol=5)
  rownames(out) <- unique(I.data$site)
  colnames(out) <- colnames(data)
  return(out)
}
print("hi")
temp.model <- SimulateModel(trial)
# print(temp.model)
diagnostic <- TestModel(trial)
diagnostic_random <- RandomModel(trial)
dev.new(xpos = 20, ypos = 20, height = 10, width = 20)
par(par_plots)
layout.matrix <- matrix(c(1, 2, 3, 4, 5, 6, 7, 4), nrow = 2, ncol = 4, byrow = TRUE)

layout(mat = layout.matrix, heights = c(0.5, 0.5))
# plot(colnames(diagnostic), rep(0,5), xlim=c(0.3, 100), ylim=c(0, 1), log="x", col = "white")
# sapply(seq(1, nrow(diagnostic)), function(row_num) {
#   lines(colnames(diagnostic), diagnostic[row_num, ], lwd=2, col=site_cols[rownames(diagnostic)[row_num], ], type="o")
#   })
column_colors <- c("gray10", "gray25", "gray30", "gray45", "gray60")
names(column_colors) <- colnames(diagnostic)
plot(colnames(diagnostic), rep(0,5),
  xlim=c(0.3, 50),
  ylim=c(-1, 1),
  log="x",
  col = "white")
title(main = "Actual data:",font = 1)
title(xlab = "% AGO (v/v)", line = 3)
title(ylab = "Correlation of log(counts) model versus data", line = 3)

sapply(seq(1, nrow(diagnostic)), function(row_num) {
  lines(colnames(diagnostic), diagnostic[row_num, ], lwd=2, col=site_cols[rownames(diagnostic)[row_num], ], type="o", pch = 20)
  })

x_y_pairs <- apply(expand.grid(rownames(diagnostic), colnames(diagnostic)), 1, function(row){
  c(diagnostic[row[1], row[2]], (sum(data[which(I.IDs[,1] == row[1]), row[2]])*sum(I.data$n[which(I.IDs[,1]==row[1])]))^(1/2), site_cols[row[1],], column_colors[row[2]])
  })

plot(x_y_pairs[2, ], x_y_pairs[1, ],
  col=x_y_pairs[3, ],
  log='x',
  xlim=c(10,10^6),
  ylim = c(-1, 1),
  cex=2,
  pch = 20)

title(main=mirna)
title(xlab = "[log(input reads) + log(selection reads)]/2", line = 3)
title(ylab = "R (log(model) vs log(data))", line = 3)

plot(x_y_pairs[2, ], x_y_pairs[1, ],
  col=x_y_pairs[4, ],
  log='x',
  xlim=c(10,10^6),
  ylim = c(-1, 1),
  cex=2,
  pch = 20)

title(xlab = "[log(input reads) + log(selection reads)]/2", line = 3)
title(ylab = "R (log(model) vs log(data))", line = 3)

plot(1, type = "n", axes=FALSE, xlab="", ylab="")

legend(x = "top", inset = 0, legend = unique(rownames(diagnostic)), col = site_cols[unique(rownames(diagnostic)),], cex = 1.5, pch = 20, bty = "n", ncol = 2)

plot(colnames(diagnostic_random), rep(0,5),
  xlim=c(0.3, 50),
  ylim=c(-1, 1),
  log="x",
  col = "white")
title(main = "Data simulated from model:",font = 1)
title(xlab = "% AGO (v/v)", line = 3)
title(ylab = "R (log(model) vs log(data))", line = 3)

sapply(seq(1, nrow(diagnostic_random)), function(row_num) {
  lines(colnames(diagnostic_random), diagnostic_random[row_num, ], lwd=2, col=site_cols[rownames(diagnostic_random)[row_num], ], type="o", pch = 20)
  })

x_y_pairs <- apply(expand.grid(rownames(diagnostic_random), colnames(diagnostic_random)), 1, function(row){
  c(diagnostic_random[row[1], row[2]], (sum(data[which(I.IDs[,1] == row[1]), row[2]])*sum(I.data$n[which(I.IDs[,1]==row[1])]))^(1/2), site_cols[row[1],], column_colors[row[2]])
  })

plot(x_y_pairs[2, ], x_y_pairs[1, ],
  col=x_y_pairs[3, ],
  log='x',
  xlim = c(10, 10^6),
  ylim = c(-1, 1),
  cex=2,
  pch = 20)
title(xlab = "[log(input reads) + log(selection reads)]/2", line = 3)
title(ylab = "R (log(model) vs log(data))", line = 3)

plot(x_y_pairs[2, ], x_y_pairs[1, ],
  col=x_y_pairs[4, ],
  log='x',
  xlim = c(10, 10^6),
  ylim = c(-1, 1),
  cex=2,
  pch = 20)
title(xlab = "[log(input reads) + log(selection reads)]/2", line = 3)
title(ylab = "R (log(model) vs. log(data)", line = 3)





