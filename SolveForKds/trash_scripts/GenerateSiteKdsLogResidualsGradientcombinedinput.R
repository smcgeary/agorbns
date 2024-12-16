################################################################################
#GenerateSiteKds_PAPER.py
################################################################################

## THIS IS THE OFFICIAL SCRIPT FOR FINALIZING THE MANUSCRIPT.
## THIS SCRIPT HAS THESE DECISIONS MADE:
## 1.) PROTEIN IS FIT AS A PARAMETER
## 2.) A SINGLE BACKGROUND PARAMETER IS GIVEN TO ALL FIVE CONCENTRATIONS.
## 3.) SITELIST, START, AND STOP POSITION ARE MANDATORY PARAMETERS.
# ## 4.) ALL INPUTS ARE COMBINED ACROSS FOR THESE ANALYSES.
print("HI")
library(colorspace)
library(multicore)
library(data.table)
source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")
# Initial parameters and constants.
args  <- commandArgs(trailingOnly=TRUE)
mirna <- args[1]
experiment <- args[2]

# # Region within random sequence from which site types orginiates,
# # going from position [26 - "start" : 26 + 37 + "stop"]
start <- as.integer(args[3])
stop  <- as.integer(args[4])
# # Parameter specifying which list of sites is used for the analysis.
# # I.E "Current" includes all the current site types, "canonical" is just
# # the six seed-based 6mers, 7mers, and 8mers, along with no site., etc.
sitelist <- args[5]
if (sitelist %in% kmer_list_names) {
  mirna.start <- as.integer(args[6])
  mirna.stop <- as.integer(args[7])
} else {
  mirna.start <- NULL
  mirna.stop <- NULL
}


# Loads general functions used in AgoRBNS analysis.

# Loads the colors associated with each site type, for plotting purposes.
print("hi")
# Loads the table of Ago–miRNA concentrations, for the purposes of modeling
# them into the structure.
# NOTE These are actually higher than the real concentrations, because before
# I began the structural analysis I had to use 1.5-2X the amount of AGO.
# TODO make a new table of miRNA concentrations or just convert within the
# script.
stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
                       "/SolveForKds/k_c_stockago.txt"), row.names=1,
                       header=TRUE, sep="\t")
k.c.stockago <- stockago[mirna,experiment]
k.c.lib <- 100

# MAIN #########################################################################
# 1. Get data table:
sitesXcounts <- GetSitesXCountsCombined(mirna,
                                experiment,
                                start,
                                stop,
                                sitelist,
                                mirna.start = mirna.start,
                                mirna.stop = mirna.stop)
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

data <- sitesXcounts[,2:6]
data <- data.matrix(data)

rownames(data)[nrow(data)] <- "None"
# data <- rbind(data_new,data[nrow(data),,drop=FALSE])

num.kds <- nrow(data)-1
num.bgs <- 1

# Define the vector of total site type concentrations:
c.I.tots <- Norm(sitesXcounts[, 1]) * k.c.lib
c.totals <- matrix(c.I.tots, nrow=length(c.I.tots), ncol=5, byrow=FALSE)

colnames(c.totals) <- colnames(data)
rownames(c.totals) <- rownames(data)

# Define the vector of AGO concentrations:
c.agos <- sapply(colnames(data), function(x) {
  as.numeric(x) / 100
})
print("hi")
# PARAMETER INITIALIZATION:
# Initialize starting Kds, which are set to 1/the enrichment of each site type
# in the A12.6 experiment.
r_numerator <- Norm(rowSums(data))
r_denomenator <- Norm(c.I.tots)
kds.init <- c(r_denomenator) / c(r_numerator)
kds.init <- kds.init/max(kds.init)
pars.init <- c(Logit(kds.init[-length(kds.init)], 10), -1, 1)

if (sitelist %in% kmer_list_names[-1]) {
  list_ind <- which(kmer_list_names == sitelist) - 1
  sitelist_prior <- kmer_list_names[list_ind]
  in_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                       experiment, "/kds_PAPER/", start, "-", stop, "_", 
                       sitelist_prior, "_", mirna.start, "-", mirna.stop,
                       "_singlebg_combinedinput_logresiduals_PAPER_last.txt")
  pars_raw <- fread(in_file, sep="\t", header=TRUE, colClasses="numeric", 
                    stringsAsFactors=FALSE)
  pars_table <- matrix(unlist(pars_raw), nrow=nrow(pars_raw),
                       ncol=ncol(pars_raw), byrow=FALSE)
  pars_old <- unlist(pars_table[nrow(pars_table),])
  names(pars_old) <- colnames(pars_raw)

  kds.init <- sapply(rownames(data)[-nrow(data)], function(name) {
    if (list_ind %% 2 == 1) {
          smaller_name <- substr(name, 2, nchar(name))
    } else {
          smaller_name <- substr(name, 1, nchar(name)-1)
    }
    sub_name <- paste0(smaller_name)
    ind <- which(names(pars_old) == sub_name)
    return(pars_old[ind])  
  })

  pars.init <- c(kds.init, -1, 1)
}

names(pars.init) <- c(rownames(data)[-nrow(data)], "bg", "AGO")

pars <- pars.init
tick <- 1

print("up to model")

if (sitelist %in% kmer_list_names) {
  n <- 256
  colors <- c(rainbow_hcl(n, c=50, l=70, start=0, end=360 * (n - 1) / n),"red",
              "black")

  out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                     experiment, "/kds_PAPER/", start, "-", stop, "_", 
                     sitelist, "_", mirna.start, "-", mirna.stop,
                     "_singlebg_combinedinput_logresiduals_PAPER.txt")
  out_last_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                          experiment, "/kds_PAPER/", start, "-", stop, "_", 
                          sitelist, "_", mirna.start, "-", mirna.stop,
                          "_singlebg_combinedinput_logresiduals_PAPER_last.txt")

  specific_figure_string <- paste0(start, "-", stop, "_", sitelist, "_",
                                   mirna.start, "-", mirna.stop,
                                   "_singlebg_combinedinput_logresiduals_PAPER")
} else {
  colors <- FALSE

  out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                     experiment, "/kds_PAPER/", start, "-", stop, "_", 
                     sitelist, "_singlebg_combinedinput_logresiduals_PAPER.txt")

  out_last_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                          experiment, "/kds_PAPER/", start, "-", stop, "_", 
                          sitelist, "_singlebg_combinedinput_logresiduals_PAPER_last.txt")

  specific_figure_string <- paste0(start, "-", stop, "_", sitelist,
                                   "_singlebg_combinedinput_logresiduals_PAPER")
}


ModelFunction <- function(pars) {
  # Split up the parameters into the kd and background parameters.
  kds  <- c(Logistic(pars[1 : num.kds], 10),1)
  names(kds) <- rownames(data)
  bgs  <- rep(10^pars[(num.kds + 1)], 5)
  stock.ago <- 10^pars[num.kds + 2]
  # Get the bound counts in each experiment:
  c.bounds <- as.matrix(
    sapply(c.agos, function(ago.percent) {
      return(GetBoundRNA(kds, c.I.tots, ago.percent * stock.ago))
    }
  ))

  # Get the amount of background binding by subtracting the bound from the
  # total sites in each experiment, normalizing. Must transpose to multiply
  # each column.
  c.frees <- c.totals - c.bounds
  c.bgs <- t(t(c.frees) * bgs / colSums(c.frees))
  c.all <- c.bounds + c.bgs
  c.final <- data.frame(t(t(c.all) / colSums(c.all) * colSums(data)))
  c.final_new <<- c.final
  return(c.final)
}

ModelLikelihood <- function(pars, plot.=FALSE){

  model <- ModelFunction(pars)
    lc.final <- log(model+1)
  ldata <- log(data + 1)

  sumofsquares <<- sum((lc.final - ldata)^2)

  if (tick %% 1000 == 0 & tick > 1 & plot.==TRUE) {
  out <<- rbind(out, c(pars, sumofsquares))
    if (sitelist %in% kmer_list_names) {
      PlotSiteKdOptimization(out[,c(seq(1,num.kds,length=256),(ncol(out)-2):ncol(out))], specific_figure_string, mirna, 256, num.bgs,
        colors=colors)

    } else {
      PlotSiteKdOptimization(out, specific_figure_string, mirna, num.kds, num.bgs, colors=colors)
    }
    write.table(file=out_file, out, sep="\t", quote=FALSE, row.names=FALSE)
    write.table(file=out_last_file, out[nrow(out),,drop=FALSE], sep="\t",
                quote=FALSE, row.names=FALSE)

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
    # segments(1, 1, x1=10^7, 10^7, lty = 3)

    # axis(1,sapply(10^(-3:8), function(x) {x * c(1,2,3,4,5,6,7,8,9)}), labels = FALSE, pos = 1, lwd = 2)
    # axis(2,sapply(10^(0:8), function(x) {x * c(1,2,3,4,5,6,7,8,9)}),labels = FALSE, pos = 10^-3,lwd = 2)
    # axis(1,10^(-3:8), pos = 1, lwd = 0)
    # axis(2,10^(0:8), pos = 10^-3,lwd = 0)


    # sapply(1:ncol(data), function(col) {
    #   x_model <- model[,col]
    #   y_data <- data[, col]
    #   x_plot <- x_model/sum(x_model)*sum(y_data)
    #   y_plot <- y_data
    #    points(x_plot, y_plot,col=site_cols[rownames(data),], pch=20)
    # })
    par(mfrow=c(1,2))
    plot(grad(ModelLikelihood,pars),Gradient(pars))
    segments(-10, -10, 10, 10, lty = 2)
    plot(grad(ModelLikelihood,pars,method="simple"),Gradient(pars))
    segments(-10, -10, 10, 10, lty = 2)

    print(sumofsquares)
  }
  tick <<- tick + 1
  return(sumofsquares)
}
Gradient <- function(pars, plot.=FALSE) {
  time_prior <- proc.time()
  # Split up the parameters into the kd and background parameters.
  kds  <- c(Logistic(pars[1 : num.kds], 10),1)
  names(kds) <- rownames(data)
  bgs  <- rep(10^pars[(num.kds + 1)], 5)
  B <- bgs[1]
  stock.ago <- 10^pars[num.kds + 2]

  f.mat <- matrix(
      rep(sapply(c.agos, function(ago.percent) {
           return(GetFreeAgo(kds, c.I.tots, ago.percent * stock.ago))
         }
         ), nrow(data)), nrow=nrow(data), ncol=ncol(data), byrow=TRUE)
  A.mat <- matrix(c.agos*stock.ago, nrow=nrow(data), ncol = ncol(data), byrow=TRUE)
  K.mat <- matrix(kds, nrow=nrow(data), ncol = ncol(data), byrow=FALSE)
  l.mat <- matrix(c.I.tots, nrow=nrow(data), ncol = ncol(data), byrow=FALSE)
  R.mat <- matrix(colSums(data), nrow = nrow(data), ncol = ncol(data), byrow=TRUE)
  L <- 100
  # Get the amount of background binding by subtracting the bound from the
  # total sites in each experiment, normalizing. Must transpose to multiply
  # each column.
  k_8 <- kds["8mer"]
  c.final <- (-R.mat*l.mat * (f.mat^2 + (L - A.mat)*f.mat + K.mat*B) *
                  (f.mat^3 - (2*A.mat + B - L - K.mat)*f.mat^2 +
                    (A.mat^2 + (B - L)*A.mat - L*B - 
                     (2*A.mat + B - L)*K.mat)*f.mat +
                    (A.mat^2 + (B - L)*A.mat - L*B)*K.mat)^-1)
  c.final <<- c.final
  # print("up to c.final")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new


  dF.base <- (1 + colSums(l.mat * K.mat * (f.mat + K.mat)^(-2)))^(-1)

  # print("dF.base")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new



  dF.dK.mat <- sweep(f.mat * l.mat * (f.mat + K.mat)^(-2),MARGIN=2,dF.base, "*")

  # print("dF.dK.mat")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new


  C1.mat <- L - B - 2*A.mat

  # print("C1.mat")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new


  C2.mat <- A.mat^2 + (B - L)*A.mat - L*B

  # print("C2.mat")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new


  # The d (each model point) d Free derivative:)
  # print("dc.ai.dF")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new

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

  # print("dc.ai.dKj_specific")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new



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

  # print("dc.ai.dA_specific")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new


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

  # print("dc.ai.dB")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new



  # grad_derivs <- sapply(seq(length(kds)), function(index) {
  #   base <- sweep(dc.ai.dF,MARGIN=2,c(dF.dK.mat[index,]),"*")
  #   base[index,] <- base[index,] + dc.ai.dKj_specific[index,]
  #   return(10*exp(pars[index]) * (exp(pars[index]) + 1)^(-2) * 2 * sum((log(c.final) - log(data+1)*base*c.final^(-1))))
  #   })
  lc.final <- log(c.final+1)
  ldata <- log(data + 1)
  residuals <- log((c.final + 1) / (data + 1))/(c.final + 1)
  grad_derivs_original <- sapply(seq(length(kds)-1), function(index) {
    base <- sweep(dc.ai.dF,MARGIN=2,c(dF.dK.mat[index,]),"*")
    base[index,] <- base[index,] + dc.ai.dKj_specific[index,]
    return(10*exp(pars[index]) * (exp(pars[index]) + 1)^(-2) * 2 * sum((lc.final - ldata)*base*(c.final+1)^(-1)))
    })

    grad_derivs <- 2*(10*exp(pars[1:num.kds]) * (exp(pars[1:num.kds]) + 1)^(-2) *
                  (colSums(residuals*dc.ai.dF) %*%
                   t(dF.dK.mat[-nrow(dF.dK.mat),]
                  ) + rowSums(dc.ai.dKj_specific*residuals)[-nrow(dF.dK.mat)]))

    # print("grad deriv")

    # time_new <- proc.time()
    # print(time_new - time_prior)
    # time_prior <- time_new


  dc.ai.dA <- sweep(sweep(dc.ai.dF,MARGIN=2,dF.base,"*") + dc.ai.dA_specific,
                    MARGIN=2,c.agos, "*")



  # print("AGO deriv")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new



  colnames(c.final) <- colnames(data)
  use <- which(c.final > 0 & data > 0)

  # print("use")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new

  # sumoflogsquares <<- sum((log(c.final) - log(data + 1))^2)
  # sumoflogsquares <<- sum((unlist(c.final) - unlist(data))^2)

  # print("logsquares")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new


  gradient_with_Ago <- log(10)*stock.ago*2*sum((log(c.final+1) - log(data + 1))*dc.ai.dA * (c.final+1)^(-1))
  # gradient_with_Ago <- log(10)*stock.ago*2*sum((unlist(c.final) - unlist(data))*dc.ai.dA)

  # print("gradient with Ago")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new


  gradient_with_bg <- log(10)*B*2*sum((log(c.final + 1) - log(data + 1))*dc.ai.dB * (c.final+1)^(-1))
  # gradient_with_bg <- log(10)*B*2*sum((unlist(c.final) - unlist(data))*dc.ai.dB)

  # print("gradient with bg")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new


  return(c(grad_derivs, gradient_with_bg, gradient_with_Ago))
}



out <- matrix(c(pars, ModelLikelihood(pars)), nrow=1)
out_names <- c(rownames(data)[-length(rownames(data))],
               "bg", "AGO", "-logp")
colnames(out) <- out_names

QuickDeriv <-function(params,width) {
  params_orig <- params
  sapply(seq(length(params)), function(i) {

    params <- params_orig
    params[i] <- params[i] + width
    return((ModelLikelihood(params) - ModelLikelihood(params_orig))/(width))
  })
}


print(data)
solution <- optim(pars,
                  ModelLikelihood,
                  gr = Gradient,
                  plot.=TRUE,
                  method = "L-BFGS-B",
                  lower=c(rep(-13, length=num.kds), -2, -1),
                  upper=c(rep(10, length=num.kds), 1, 2),
                  control = list(maxit=100000, factr=1e2))

pars <- solution$par
sumofsquares <- solution$value
out <- rbind(out, c(pars, sumofsquares))

print(sumofsquares)
write.table(file=out_file, out, sep="\t", quote=FALSE, row.names=FALSE)
write.table(file=out_last_file, out[nrow(out),,drop=FALSE], sep="\t",
                quote=FALSE, row.names=FALSE)

warnings()





