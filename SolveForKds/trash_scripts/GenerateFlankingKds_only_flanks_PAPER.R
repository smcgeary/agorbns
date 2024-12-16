################################################################################
#GenerateSiteTypeKds.py
################################################################################

## THIS IS THE OFFICIAL SCRIPT FOR FINALIZING THE MANUSCRIPT.
## THIS SCRIPT HAS THESE DECISIONS MADE:
## 1.) PROTEIN IS FIT AS A PARAMETER
## 2.) A SINGLE BACKGROUND PARAMETER IS GIVEN TO ALL FIVE CONCENTRATIONS.
## 3.) LIST, START, AND STOP POSITION ARE MANDATORY PARAMETERS.
# ## 4.) ALL INPUTS ARE COMBINED ACROSS FOR THESE ANALYSES.
library(colorspace)
library(multicore)
library(data.table)
library(numDeriv)
source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")
# Initial parameters and constants.
# args       <- commandArgs(trailingOnly=TRUE)
# mirna      <- args[1]
# experiment <- args[2]
# n_constant <- args[3]
# sitelist   <- args[4]
# site       <- args[5]
# if (length(args) == 6) {
#   reps = as.integer(args[6])
# } else {
#   reps = 200
# }
# # Region within random sequence from which site types orginiates,
# # going from position [26 - "start" : 26 + 37 + "stop"]stop  <- as.integer(args[3])
# # Parameter specifying which list of sites is used for the analysis.
# # I.E "Current" includes all the current site types, "canonical" is just
# # the six seed-based 6mers, 7mers, and 8mers, along with no site., etc.

# Loads general functions used in AgoRBNS analysis.

# TODO make a new table of miRNA concentrations or just convert within the
# script.
stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
                       "/SolveForKds/k_c_stockago.txt"), row.names=1,
                       header=TRUE, sep="\t")
k.c.stockago <- stockago[mirna,experiment]
k.c.lib <- 100


# MAIN #########################################################################

# Get siteXcounts, called sitesXcounts:
sitesXcounts <- GetSitesXCounts(mirna,
                                experiment,
                                n_constant,
                                sitelist,
                                mirna.start = mirna.start,
                                mirna.stop = mirna.stop)
print(sitesXcounts)
# Separate site sequences from data file.
if (sitelist %in% kmer_list_names) {
  seqs <- rownames(sitesXcounts)
} else {
  seqs <- sitesXcounts[,1]
  sitesXcounts <- sitesXcounts[,-1]
}

# Remove any rows with no reads whoatsoever.
sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
# Remove any rows for which there are no input reads:
sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]

sitesXcounts_parent <- sitesXcounts

params <- GetSiteKds(mirna, experiment, n_constant, sitelist)
  tick <- 0

  sitesXcounts <- sitesXcounts_parent
  sitesXcounts.all <- sitesXcounts_parent["None", , drop=FALSE]
  # Get vector of single site counts s.c.
    site_colors = c()
  print(sitesXcounts.all)
  for (ind in seq(nrow(sitesXcounts)-1)) {
    site <- rownames(sitesXcounts)[ind]
  sfXc <- GetSiteFlanksXCounts(mirna, experiment, n_constant, sitelist, site)
  colnames(sfXc)[1] <- "I"
  sfXc <- sfXc[rowSums(sfXc[, 2:6]) > 0,]
  sfXc <- sfXc[which(sfXc[,1] > 0),,drop=FALSE]

  rownames(sfXc) <- paste0(site,"|",rownames(sfXc))

  sitesXcounts.all <- rbind(sitesXcounts.all[0:(nrow(sitesXcounts.all)-1),],sfXc, sitesXcounts.all[nrow(sitesXcounts.all),])  # Assign parameters for Kds and background, subracting 1 from rows due to "None"
  site_colors <- c(site_colors, rep(kSiteColors[site, ], nrow(sfXc)))
  }
  sitesXcounts <- sitesXcounts.all+1
  site_colors <- c(site_colors, "black")
  # Kd assigned to 1 nM, and the fact that the I and A0 columns are not modeled
  # in the experiment, so there's no background term assigned ot them.

  # Get the number of flanking sites
  num.bgs <- 1

  # Get the number of parameters (kds + bg + 1 log(prob))
  bgs <- rep(params$Mean[(num.kds + 1): (num.kds + num.bgs)], 5)
  stock.ago <- params$Mean[num.kds + num.bgs + 1]

  # Get the site kds

  # Assign the total site concentration in each experiment, and initialize the
  # data column to be used.
  # k.c.lib should be 100, for 100 nM.
  # sitesXcounts <- sitesXcounts[rownames(sitesXcounts) != site, ]
  # print(sitesXcounts)

  rownames(sitesXcounts)[nrow(sitesXcounts)] <- "None"


  c.I.tots <- Norm(sitesXcounts[, 2]) * k.c.lib
  names(c.I.tots) <- rownames(sitesXcounts)
  # Define the vector of total site type concentrations:
  c.totals <- matrix(c.I.tots, nrow=length(c.I.tots), ncol=5, byrow=FALSE)
  # Remove the I and A0 columns from the data to be fit to the model. 
  data <- sitesXcounts[, 3:7]
  data <- data.matrix(data)
  num.kds <- nrow(data)
  data.all <- data
  rownames(c.totals) <- rownames(data)
  colnames(c.totals) <- colnames(data)
  colors_flanks <- sapply(rownames(sfXc), GetColorFunction)
  colors_all <- c(colors_flanks, kSiteColors[rownames(data)[257:nrow(data)], ])
  colors_all <- site_colors
  pch_all <- 1

  c.agos <- sapply(colnames(data), function(x) {
     as.numeric(x) / 100
    })
  # Initialize starting Kds, which are set to 1/the enrichment of each site type
  # in the A12.6 experiment.
  enrichments.init <- (Norm(c.I.tots)/Norm(rowSums(data)+1))
  pars.init <- log10(c(enrichments.init, 0.1, 10))
  names(pars.init) <- c(rownames(data), "bg", "AGO")
  pars <- pars.init
  # Define function of just kds and pars.


  out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                     experiment, "/kds_PAPER/", n_constant, "_", 
                     sitelist, "_", site,"_onlyflanks_PAPER.txt")

GetSiteKmers <- function(site, cond) {
  data <- read.table(paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                            "/equilibrium/sitekmer_counts/", cond,
                            "_", n_constant,"_", sitelist, "_k8.txt"))
  site_kmers <- data[,site]
  names(site_kmers) <- rownames(data)
  return(site_kmers)
}

GetTopKmers <- function(site, cond) {
  kmers <- GetSiteKmers(site, cond)
  kmers_I <- GetSiteKmers(site, "I_combined")
    print(rev(sort((kmers/sum(as.numeric(kmers))/
      ((kmers_I/sum(as.numeric(kmers_I)))))))[1:20])
}

break

ModelFunction <- function(pars) {
  # print(c.totals[1,1])
  # print(c.I.tots[1])
  # Split up the parameters into the kd and background parameters.
  print(length(pars))
  kds  <- 10^pars[1 : num.kds]
  bgs  <- rep(10^pars[(num.kds + 1)], ncol(data))
  stock.ago <- 10^pars[num.kds + 2]
  # Get the bound counts in each exp:
  c.agos <- sapply(colnames(data), function(x) {
    as.numeric(x) / 100
  })
  c.bounds <- as.matrix(
    sapply(c.agos, function(percent) {
      return(GetBoundRNA(kds, c.I.tots, percent * stock.ago))
    }
  ))

  # Get the amount of background binding by subtracting the bound from the
  # total sites in each exp, normalizing. Must transpose to multiply
  # each column.
  c.frees <- c.totals - c.bounds
  c.bgs <- t(t(c.frees) * bgs / colSums(c.frees))
  c.all <- c.bounds + c.bgs
  c.final <- data.frame(t(t(c.all) / colSums(c.all) * colSums(data)))
  return(c.final)
}
tick <- 0
ModelLikelihood <- function(pars, print.=FALSE){
  model <- ModelFunction(pars)
  print(dim(model))
  model_norm <- t(t(model) / colSums(model))
  loglikelihood <- -sum(data*log(model_norm))
  data_norm <- t(t(data) / colSums(data))
  model_R <- model_norm / Norm(c.I.tots)
  data_R <- data_norm / Norm(c.I.tots)
  if (print.==TRUE) {
  plot(1, "n", log = 'xy', xlim = c(0.1, 100), ylim = c(0.1, 1000))
  sapply(seq(nrow(data)), function(row){
    lines(colnames(data), model_R[row,], col = colors_all[row])
    points(colnames(data), data_R[row, ], col = colors_all[row])
    })
}
  if (tick %% 10 == 0) {
  plot(c(model_norm), c(data_norm), col = colors_all, log = 'xy')
  abline(0, 1, lty = 2)
}
tick <<- tick + 1
  return(loglikelihood)
}

Gradient <- function(pars) {
  # Split up the parameters into the kd and background parameters.
  kds  <- 10^pars[1 : num.kds]
  B  <- 10^pars[(num.kds + 1)]
  stock.ago <- 10^pars[num.kds + 2]
  c.agos <- sapply(colnames(data), function(x) {
    as.numeric(x) / 100
  })

   f.mat <- matrix(
      rep(sapply(c.agos, function(percent) {
           return(GetFreeAgo(kds, c.I.tots, percent * stock.ago))
         }
         ), nrow(data)), nrow=nrow(data), ncol=ncol(data), byrow=TRUE)

  f.jvec <- sapply(c.agos, function(ago.percent) {
           return(GetFreeAgo(kds, c.I.tots, ago.percent * stock.ago))})

  A.mat <- matrix(c.agos*stock.ago, nrow=nrow(data), ncol=ncol(data),
                  byrow=TRUE)

  A.jvec <- c.agos * stock.ago

  K.mat <- matrix(kds, nrow=nrow(data), ncol=ncol(data), byrow=FALSE)
  K.ivec <- kds

  l.mat <- matrix(c.I.tots, nrow=nrow(data), ncol=ncol(data), byrow=FALSE)
  l.ivec <- c.I.tots

  R.mat <- matrix(colSums(data), nrow=nrow(data), ncol=ncol(data), byrow=TRUE)
  R.jvec <- colSums(data) 
  L <- 100

  # Get the amount of background binding by subtracting the bound from the
  # total sites in each exp, normalizing. Must transpose to multiply
  # each column.

  time.init <- proc.time()
  c.vec_num <- -(
                 (l.ivec %*% t(R.jvec * f.jvec * (f.jvec + L - A.jvec)) +
                 (l.ivec * K.ivec * B)  %*% t(R.jvec))
                ) 
  C1.jvec <- L - B - 2 * A.jvec
  C2.jvec <- A.jvec^2 + (B - L) * A.jvec - L * B


  c.vec_dem <- t(
                (f.jvec^3 + f.jvec^2 * C1.jvec + f.jvec * C2.jvec) +
                t(K.ivec %*% t(f.jvec^2 + f.jvec * C1.jvec + C2.jvec))
                )^-1
  c.final <- c.vec_num * c.vec_dem

  dF.base <- (1 + colSums(l.mat * K.mat * (f.mat + K.mat)^(-2)))^(-1)

  dF.dK.mat <- sweep(f.mat * l.mat * (f.mat + K.mat)^(-2),MARGIN=2,dF.base, "*")

  C1.mat <- L - B - 2*A.mat
  C2.mat <- A.mat^2 + (B - L)*A.mat - L*B

  # The d (each model point) d Free derivative:)
  dc.ai.dF <- R.mat * l.mat * (
    (
      f.mat^4
    ) + (
      2 * (L - A.mat) * f.mat^3
    ) + (
      ((L - A.mat)^2 + K.mat * (4 * B + A.mat)) * f.mat^2
    ) + (
      2 * K.mat * (A.mat * (L - B) + B * (K.mat - 2 * L) - (A.mat + B)^2) * f.mat
    ) + (
      K.mat * (A.mat^3 + 2 * A.mat^2 * (B - L) - B * (B - L) * (K.mat + L) +
               A.mat * (B^2 - 2 * B * K.mat - 3 * B * L + L^2))
    )
  ) * ((K.mat + f.mat)^2 * (C2.mat + C1.mat * f.mat + f.mat^2)^2)^(-1)



  dc.ai.dKj_specific <- R.mat * l.mat * (
    (
      f.mat^4
    ) + (
      (2 * C1.mat + A.mat) * f.mat^3
    ) + (
      (C2.mat + (C1.mat + A.mat) * C1.mat) * f.mat^2
    ) + (
      (C2.mat * (C1.mat + A.mat)) * f.mat
    )
  ) * ((K.mat + f.mat)^2 * (C2.mat + C1.mat * f.mat + f.mat^2)^2)^(-1)

  dc.ai.dA_specific <- R.mat * l.mat * (
    (
      -f.mat^3
    ) + (
      (2 * (A.mat - L)) * f.mat^2
    ) + (
      ((A.mat - L) * C1.mat - 2 * B * K.mat + C2.mat) * f.mat
    ) + (
      - C1.mat * B * K.mat
    )
  ) * ((K.mat + f.mat) * (C2.mat + C1.mat * f.mat + f.mat^2)^2)^(-1)

  dc.ai.dB <- R.mat * l.mat * (
    (
      -f.mat^3
    ) + (
      (2 * (A.mat - L) - K.mat) * f.mat^2
    ) + (
      ((L - A.mat) * (A.mat - L) - K.mat * (B + C1.mat)) * f.mat
    ) + (
      + A.mat * B * K.mat - B * K.mat * L - K.mat * C2.mat
    )
  ) * ((K.mat + f.mat) * (C2.mat + C1.mat * f.mat + f.mat^2)^2)^(-1)

  residuals <- -data/c.final

  grad_derivs <- (log(10)*kds *
                  (colSums(residuals*dc.ai.dF) %*%
                   t(dF.dK.mat
                  ) + rowSums(dc.ai.dKj_specific*residuals)))
  
  dc.ai.dA <- sweep(sweep(dc.ai.dF,MARGIN=2,dF.base,"*") + dc.ai.dA_specific,
                    MARGIN=2,c.agos, "*")

  colnames(c.final) <- colnames(data)

  gradient_with_Ago <- log(10)*stock.ago*sum(residuals * dc.ai.dA)

  gradient_with_bg <- log(10)*B*sum(residuals * dc.ai.dB)
  gradient_all <- c(grad_derivs, gradient_with_bg, gradient_with_Ago)
  names(gradient_all) <- names(pars)
  # print(length(gradient_all))
  # print(length(pars))
  return(gradient_all)
}







  print(data)
  solution <- optim(pars,
                    ModelLikelihood,
                    gr = Gradient,
                    method = "L-BFGS-B",
                    lower=c(rep(-13, length=length(pars))),
                    upper=c(rep(10, length=length(pars))),
                    control = list(maxit=100000, factr=1e2))

  pars <- solution$par

print(pars)

ModelLikelihood(pars, print.=TRUE)
  # plot(colnames(data),rep(1,5),type="l",lty=2,xlim=c(0.1,100),ylim=c(0.3,1000),log='xy')
  # sapply(seq(nrow(model)), function(row){
  #   print(row)
  #   points(colnames(data),data[row,]/colSums(data)*sum(c.I.tots)/c.I.tots[row],col=kSiteColors[rownames(data)[row],])
  #   lines(colnames(data),model[row,]/colSums(model)*sum(c.I.tots)/c.I.tots[row],col=kSiteColors[rownames(data)[row],])
  #   })
  # print(out_file)
  # print(pars)
  pars.save <- pars
  # print(reps)

  pars_loocv <- matrix(NaN,nrow=length(pars),ncol=reps)
  rownames(pars_loocv) <- names(pars)
  colnames(pars_loocv) <- seq(reps)


break
  ## Optimization of parameters:
  for (i_trial in seq(reps)) {
    print(i_trial)
    tick <- 0
    i_col <- sample(1:ncol(data.all),1)
    data.temp <- data.all[,-i_col]
    data <- apply(data.temp, 2, function(col) {rmultinom(1,size=sum(col), prob=col)})
    rownames(data) <- rownames(data.temp)
    colnames(data) <- colnames(data.temp)
    input_resample <- rmultinom(1,size=sum(sitesXcounts[,2]),prob=sitesXcounts[,2])
    c.I.tots <- Norm(input_resample+1)*k.c.lib
    c.totals <- matrix(
                  rep(c.I.tots,ncol(data)), nrow=length(c.I.tots), ncol=ncol(data), byrow=FALSE)
    colnames(c.totals) <- colnames(data)
    rownames(c.totals) <- rownames(data)
    pars <- pars.save
    ind_print <- which(rownames(data) == "AC.CA")
    solution <- optim(pars,
                      ModelLikelihood,
                      gr = Gradient,
                      method = "L-BFGS-B",
                      lower=c(rep(-13, length=length(pars))),
                      upper=c(rep(10, length=length(pars))),
                      control = list(maxit=100000, factr=1e2))

    pars <- solution$par

    pars_loocv[,i_trial] <- pars
    
  }


  pars_loocv_sort <- t(apply(pars_loocv,1,sort))
  output <- 10^cbind(pars.save,rowMeans(pars_loocv),
                     pars_loocv_sort[,ceiling(0.025*reps)],
                     pars_loocv_sort[,ceiling(0.5*reps)],
                     pars_loocv_sort[,ceiling(0.975*reps)])
  colnames(output) <- c("Full","Mean", "Lower_CI", "Median", "Upper_CI")
  write.table(file=out_file, output, sep="\t", quote = FALSE)
  print(out_file)


  warnings()

