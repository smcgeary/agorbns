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
library(wrswoR)
source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")
mirna <- "miR-1"
experiment <- "equilibrium"
n_constant <- 5
sitelist <- "paper"
site <- "8mer"
mir.start <- 1
mir.stop <- 15
reps = 100

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

  NumFlanks <- function(data) {
    flank_abundances <- aggregate(. ~ Flank, data, function(x) {
      length(x)
    })
    flank_num <- flank_abundances[,2]
    names(flank_num) <- flank_abundances$Flank
    names(flank_num) <- sapply(names(flank_num), function(name) {
      temp <- unlist(strsplit(name, ""))
      return(paste0(temp[1:2], temp[3:4], collapse = "."))
      })
    print(flank_num[1:10])
    return(flank_num)
  }

sample_distribution <- function(condition, score="plfold") {
  dist.I <- data.frame(GetPairingFlankData(mirna, experiment, "I_combined",
                                   n_constant, sitelist, site, mir.start, mir.stop))
  dist <- data.frame(GetPairingFlankData(mirna, experiment, condition,
                                   n_constant, sitelist, site, mir.start, mir.stop))

  # Make vector to sample from dists 1 and 2
  input <- dist.I[[score]]
  len <- length(input)
  len2 <- nrow(dist)
  # Make the target mean:
  target <- mean(log10(dist[[score]]^15))
  target2 <- sd(log10(dist[[score]]^15))
  # target <- mean(dist2[[score]])
  # Define Cost function:
  tick <- 1
  costfunction_plfold <- function(pars) {
    tick <<- tick + 1
    inds <- c()
    while (length(inds) < len/10){
      inds <- c(inds,
                sample_int_rej(len,
                             size = as.integer(Logistic(pars[1],1)*len),
                             prob = input^pars[2] + 10^pars[3]))

    }
    residual <- (
      (mean(log10(input[inds]^15)) - target)^2 + 
      (sd(log10(input[inds]^15))   - target2)^2
      )
    # residual <- (mean(input[inds]) - target)^2
    return(residual)
  }
  pars <- optim(c(0,1,-10),costfunction_plfold)$par
  inds <- c()
  while (length(inds) < len2) {
    inds <- c(inds,
              sample_int_rej(len,
                           size = as.integer(Logistic(pars[1],1)*len),
                           prob = input^pars[2] + 10^pars[3]))
    }
    NumFlanks(dist.I[inds,])
  }
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
  # Get vector of single site counts s.c.
  s.c <- as.numeric(sitesXcounts[site, ])

  sfXc <- GetSiteFlanksXCounts(mirna, experiment, n_constant, sitelist, site)
  colnames(sfXc)[1] <- "I"
  sfXc <- sfXc[rowSums(sfXc[, 2:6]) > 0,]
  sfXc <- sfXc[which(sfXc[,1] > 0),,drop=FALSE]

  sfXc <- round(t((t(sfXc) / colSums(sfXc)) * s.c))
  sfXc[is.na(sfXc)] <- 0

  #   siteflanksXprob.I <- data.frame(GetPairingFlankData(mirna, experiment, "I_combined",
  #                                  n_constant, sitelist, site, mir.start, mir.stop))
  #   max_prob <- max(siteflanksXprob.I$plfold^15)
  #   min_prob <- min(siteflanksXprob.I$plfold^15)

  # prob_limits <- c(0,10^seq(log10(min_prob*10),log10(max_prob/10),length=15),1)
  # prob_strings <- as.character(prob_limits[1:(length(prob_limits)-1)])

  # sfXc_prob <- matrix(0, nrow=nrow(sfXc)*length(prob_strings), ncol=ncol(sfXc))
  # colnames(sfXc_prob) <- c("I", "I_combined", 40, 12.6, 4, 1.26, 0.4, 0)
  # rownames(sfXc_prob) <- sapply(rownames(sfXc), function(flank) {
  #   sapply(prob_strings, function(prob_string) {
  #     paste0(flank, "_", prob_string, collapse="")
  #   })
  #   })

  # print(sfXc_prob[1:10,])
  sfXc_new <- sfXc
  print(sfXc[1:10,])


  colnames(sfXc_new) <- c("I", "I_combined", 40, 12.6, 4, 1.26, 0.4, 0)
  for (condition in c(40, 12.6, 4, 1.26, 0.4)) {
    print(condition)
      print(sfXc_new[1:10,])
      new_column <- sample_distribution(condition)
      print(new_column[1:10])
      print(sfXc_new[1:10, as.character(condition)])

    sfXc_new[,as.character(condition)] <- new_column
          print(sfXc_new[1:10,])

  }

  # Assign parameters for Kds and background, subracting 1 from rows due to "None"
  # Kd assigned to 1 nM, and the fact that the I and A0 columns are not modeled
  # in the experiment, so there's no background term assigned ot them.

  # Get the number of flanking sites

  num.kds <- length(params$Mean) - 2
  num.bgs <- 1

  num.sf <- nrow(sfXc_new)
  # Get the number of parameters (kds + bg + 1 log(prob))
  bgs <- rep(params$Mean[(num.kds + 1): (num.kds + num.bgs)], 5)
  stock.ago <- params$Mean[num.kds + num.bgs + 1]

  # Get the site kds
  kds.s <- params$Mean[1:num.kds]
  names(kds.s) <- rownames(params)[1:num.kds]
  kd.site <- kds.s[names(kds.s)==site]
  # Omit the site kd for which the flanking sites are being fit.
  kds.s <- kds.s[names(kds.s) != site]

  # Assign the total site concentration in each experiment, and initialize the
  # data column to be used.
  # k.c.lib should be 100, for 100 nM.
  colnames(sfXc_new) <- colnames(sitesXcounts)
  sitesXcounts <- rbind(sitesXcounts, sfXc_new)

  sitesXcounts <- sitesXcounts[rownames(sitesXcounts) != site, ]

  c.I.tots <- Norm(sitesXcounts[, 2]+1) * k.c.lib
  names(c.I.tots) <- rownames(sitesXcounts)
  # Define the vector of total site type concentrations:
  c.totals <- matrix(c.I.tots, nrow=length(c.I.tots), ncol=5, byrow=FALSE)
  # Remove the I and A0 columns from the data to be fit to the model. 
  data <- sitesXcounts[, 3:7]
  data <- data.matrix(data)
  data.all <- data
  rownames(c.totals) <- rownames(data)
  colnames(c.totals) <- colnames(data)
  colors_flanks <- rep(sapply(rownames(sfXc), GetColorFunction), each = 10)
  colors_all <- c(kSiteColors[rownames(data)[1:(num.kds-1)],], colors_flanks)
  pch_all <- c(rep(1, num.kds), rep(19, nrow(sfXc_new)))
  c.agos <- sapply(colnames(data), function(x) {
     as.numeric(x) / 100
    })

  # Initialize starting Kds, which are set to 1/the enrichment of each site type
  # in the A12.6 experiment.
  site_colors = sample(colors(),nrow(data),replace=TRUE)
  enrichments.init <- (Norm(c.I.tots)/Norm(rowSums(data)+1))[(nrow(data) - nrow(sfXc_new) + 1):nrow(data)]
  pars.init <- log10(enrichments.init)
  pars.init <- rep(log10(kd.site),nrow(sfXc_new))
  names(pars.init) <- rownames(sfXc_new)
  pars <- pars.init
  # Define function of just kds and pars.
  print("hi")

  out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                     experiment, "/kds_PAPER/", n_constant, "_", 
                     sitelist, "_", site,"_controlplfold_PAPER.txt")



  ModelFunction <- function(pars) {
    # Split up the parameters into the kd and background parameters.
    kds  <- c(kds.s, 10^pars)
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
    c.bgs <- t(t(c.frees) * bgs[1:ncol(data)] / colSums(c.frees))
    c.all <- c.bounds + c.bgs
    c.final <- data.frame(t(t(c.all) / colSums(c.all) * colSums(data)))
    return(c.final)
  }

  ModelLikelihood <- function(pars,print.=FALSE){
    model <- ModelFunction(pars)
    model_norm <- t(t(model) / colSums(model))
    loglikelihood <- -sum(data*log(model_norm))
    tick <<- tick+1
    if (print. == TRUE) {
    # if (tick %% 20 == 0) {
      plot(unlist(c(model)),c(data), log = 'xy', col = colors_all, xlim = c(1, 1e7), ylim = c(1, 1e7))
      # identify(unlist(c(model)),c(data), labels=rownames(data))

    }
    return(loglikelihood)
  }


  Gradient <- function(pars,print.=FALSE) {
    # Split up the parameters into the kd and background parameters.
    kds  <- c(kds.s, 10^pars)
    B  <- bgs[1]
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

    gradient_all <- grad_derivs[num.kds:length(grad_derivs)]
    names(gradient_all) <- names(pars)
    return(gradient_all)
  }

  solution <- optim(pars,
                    ModelLikelihood,
                    gr = Gradient,
                    method = "L-BFGS-B",
                    lower=c(rep(-13, length=length(pars))),
                    upper=c(rep(10, length=length(pars))),
                    control = list(maxit=100000, factr=1e2))

  pars <- solution$par


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

  pars_loocv <- matrix(NaN,nrow=length(pars),ncol=1)
  rownames(pars_loocv) <- names(pars)
  pars_loocv[,1] <- pars

  ModelLikelihood(pars,print.=TRUE)


  ## Optimization of parameters:
  for (i_trial in seq(reps)) {
    print(i_trial)
    tick <- 0
    i_col <- sample(1:ncol(data.all),1)
    data.temp <- data.all[,-i_col]
    data <- apply(data.temp, 2, function(col) {rmultinom(1,size=sum(col), prob=col)})
    rownames(data) <- rownames(data.temp)
    colnames(data) <- colnames(data.temp)
    print("hi")
    input_resample <- rmultinom(1,size=sum(sitesXcounts[,2]),prob=sitesXcounts[,2])+1
    c.I.tots <- Norm(input_resample)*k.c.lib
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
    ModelLikelihood(pars,print.=TRUE)
    pars_loocv <- cbind(pars_loocv,pars)
      pars_loocv_sort <- t(apply(pars_loocv,1,sort))
    output <- 10^cbind(pars.save,rowMeans(pars_loocv),
                       pars_loocv_sort[,ceiling(0.025*ncol(pars_loocv_sort))],
                       pars_loocv_sort[,ceiling(0.5*ncol(pars_loocv_sort))],
                       pars_loocv_sort[,ceiling(0.975*ncol(pars_loocv_sort))])
    colnames(output) <- c("Full","Mean", "Lower_CI", "Median", "Upper_CI")
    write.table(file=out_file, output, sep="\t", quote = FALSE)

  }




  warnings()

