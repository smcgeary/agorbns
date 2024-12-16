################################################################################
#GenerateSiteKds_combinedinput_PAPER.py
################################################################################

## THIS IS THE OFFICIAL SCRIPT FOR FINALIZING THE MANUSCRIPT.
## THIS SCRIPT HAS THESE DECISIONS MADE:
## 1.) PROTEIN IS FIT AS A PARAMETER
## 2.) A SINGLE BACKGROUND PARAMETER IS GIVEN TO ALL FIVE CONCENTRATIONS.
## 3.) SITELIST, START, AND STOP POSITION ARE MANDATORY PARAMETERS.
# ## 4.) ALL INPUTS ARE COMBINED ACROSS FOR THESE ANALYSES.
library(colorspace)
library(multicore)
library(data.table)
# # Initial parameters and constants.
args  <- commandArgs(trailingOnly=TRUE)
mirna <- args[1]
# # Region within random sequence from which site types orginiates,
# # going from position [26 - "start" : 26 + 37 + "stop"]
start <- as.integer(args[2])
stop  <- as.integer(args[3])
# # Parameter specifying which list of sites is used for the analysis.
# I.E "Current" includes all the current site types, "canonical" is just
# the six seed-based 6mers, 7mers, and 8mers, along with no site., etc.
sitelist <- args[4]
if (sitelist %in% c("12mers", "10mers")) {
  mirna.start <- as.integer(args[5])
  mirna.stop <- as.integer(args[6])
} else {
  mirna.start <- NULL
  mirna.stop <- NULL
}


# # Experiment name
experiment <- "equilibrium"

# Loads general functions used in AgoRBNS analysis.
source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")

# Loads the colors associated with each site type, for plotting purposes.

# Loads the table of Ago–miRNA concentrations, for the purposes of modeling
# them into the structure.
# NOTE These are actually higher than the real concentrations, because before
# I began the structural analysis I had to use 1.5-2X the amount of AGO.
# TODO make a new table of miRNA concentrations or just convert within the
# script.
stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
                       "/SolveForKds/k_c_stockago.txt"), row.names=1,
                       header=TRUE, sep="\t")
print(stockago)
k.c.stockago <- stockago[mirna,experiment]
print(k.c.stockago)
k.c.lib = 100


# MAIN #########################################################################
# 1. Get data table:
sitesXcounts <- GetSitesXCounts(mirna,
                                        experiment,
                                        start,
                                        stop,
                                        sitelist,
                                        mirna.start = mirna.start,
                                        mirna.stop = mirna.stop)
# Separate site sequences from data file.
if (sitelist %in% c("12mers", "10mers")) {
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

# data_new <- t(sapply(1:((nrow(data)-1)/4),function(x){colSums(data[1:4+(x-1)*4,])}))
# rownames(data_new) <- rownames(data)[seq(1,nrow(data)-1,by=4)]
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
print(c.agos)
# PARAMETER INITIALIZATION:
# Initialize starting Kds, which are set to 1/the enrichment of each site type
# in the A12.6 experiment.
r_numerator <- Norm(rowSums(data))
r_denomenator <- Norm(c.I.tots)
# print(pars_10mer)
kds.init <- c(r_denomenator + 0.000001)/c(r_numerator + 0.000001)
kds.init <- kds.init/max(kds.init)
print(range(kds.init))
pars.init <- c(Logit(kds.init[-length(kds.init)], 10), -1,0)
names(pars.init) <- c(rownames(data)[-nrow(data)], "bg", "AGO")
pars <- pars.init
# print(pars)
tick <- 1
print("up to model")

if (sitelist %in% c("12mers", "10mers")) {
  specific_figure_string <- paste0(start, "-", stop, "_", sitelist, "_",
    mirna.start, "-", mirna.stop, "_singlebg_combinedinput_PAPER_cgfinal")
} else {
  specific_figure_string <- paste0(start, "-", stop, "_", sitelist, "_",
    "singlebg_combinedinput_PAPER_cgfinal")
}



ModelLikelihood <- function(pars) {
  time_prior <- proc.time()
  print(pars)
  # Split up the parameters into the kd and background parameters.
  kds  <- c(Logistic(pars[1 : num.kds], 10),1)
  names(kds) <- rownames(data)
  bgs  <- rep(10^pars[(num.kds + 1)], 5)
  B <- bgs[1]
  stock.ago <- 10^pars[num.kds + 2]
  # Get the bound counts in each experiment:


  f.mat <- matrix(
      rep(sapply(c.agos, function(ago.percent) {
           return(GetFreeAgo(kds, c.I.tots, ago.percent * stock.ago))
         }
         ), nrow(data)), nrow=nrow(data), ncol=ncol(data), byrow=TRUE) 
  f.jvec <- sapply(c.agos, function(ago.percent) {
           return(GetFreeAgo(kds, c.I.tots, ago.percent * stock.ago))})

  A.mat <- matrix(c.agos*stock.ago, nrow=nrow(data), ncol = ncol(data), byrow=TRUE)
  A.jvec <- c.agos * stock.ago

  K.mat <- matrix(kds, nrow=nrow(data), ncol = ncol(data), byrow=FALSE)
  K.ivec <- kds

  l.mat <- matrix(c.I.tots, nrow=nrow(data), ncol = ncol(data), byrow=FALSE)
  l.ivec <- c.I.tots

  R.mat <- matrix(colSums(data), nrow = nrow(data), ncol = ncol(data), byrow=TRUE)
  R.jvec <- colSums(data)

  L <- 100
  # Get the amount of background binding by subtracting the bound from the
  # total sites in each experiment, normalizing. Must transpose to multiply
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




  colnames(c.final) <- colnames(data)
  sumofsquares <- sum((c.final - data)^2)

  if (tick %% 1 == 0 & tick > 1) {
  out <<- rbind(out, c(pars, sumofsquares))

if (sitelist %in% c("12mers", "10mers")) {
      out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                       experiment, "/kds_PAPER/", start, "-", stop, "_", 
                       sitelist, "_", mirna.start, "-", mirna.stop,
                       "_singlebg_combinedinput_PAPER_cgfinal.txt")

      out_last_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                 experiment, "/kds_PAPER/", start, "-", stop, "_", 
                 sitelist, "_", mirna.start, "-", mirna.stop,
                 "_singlebg_combinedinput_PAPER_cgfinal_last.txt")

      PlotSiteKdOptimization(out[,c(seq(1,num.kds,length=256),(ncol(out)-2):ncol(out))], specific_figure_string, mirna, 256, num.bgs, colors=colors)

    } else {
      out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                       experiment, "/kds_PAPER/", start, "-", stop, "_", 
                       sitelist, 
                       "_singlebg_combinedinput_PAPER_cgfinal.txt")

      out_last_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                 experiment, "/kds_PAPER/", start, "-", stop, "_", 
                 sitelist,
                 "_singlebg_combinedinput_PAPER_cgfinal_last.txt")



            PlotSiteKdOptimization(out, specific_figure_string, mirna, num.kds, num.bgs, colors=colors)

    }
    write.table(file=out_file, out, sep="\t", quote=FALSE, row.names=FALSE)
    write.table(file=out_last_file, out[nrow(out),,drop=FALSE], sep="\t", quote=FALSE, row.names=FALSE)

    print(sumofsquares)
  }

  tick <<- tick + 1
  return(sumofsquares)
}




Gradient <- function(pars) {
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
  f.jvec <- sapply(c.agos, function(ago.percent) {
           return(GetFreeAgo(kds, c.I.tots, ago.percent * stock.ago))})

  A.mat <- matrix(c.agos*stock.ago, nrow=nrow(data), ncol = ncol(data), byrow=TRUE)
  A.jvec <- c.agos * stock.ago

  K.mat <- matrix(kds, nrow=nrow(data), ncol = ncol(data), byrow=FALSE)
  K.ivec <- kds

  l.mat <- matrix(c.I.tots, nrow=nrow(data), ncol = ncol(data), byrow=FALSE)
  l.ivec <- c.I.tots

  R.mat <- matrix(colSums(data), nrow = nrow(data), ncol = ncol(data), byrow=TRUE)
  R.jvec <- colSums(data)

  L <- 100
  # Get the amount of background binding by subtracting the bound from the
  # total sites in each experiment, normalizing. Must transpose to multiply
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




  time.init <- proc.time()
  dF.base <- (1 + colSums(l.mat * K.mat * (f.mat + K.mat)^(-2)))^(-1)
  time.new <- proc.time()
  dF.base.jvec <- (1 + colSums(l.ivec * K.ivec * (t(f.jvec^2 + t(2 * K.ivec %*% t(f.jvec) + K.ivec^2)))^-1))^-1
  time.new2 <- proc.time()


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

  time.init <- proc.time()

  # dc.ai.dF_vec_num <- R.mat * l.mat * (
  #   (
  #     f.mat^4
  #   ) + (
  #     2 * (L - A.mat) * f.mat^3
  #   ) + (
  #     ((L - A.mat)^2 + K.mat * (4 * B + A.mat)) * f.mat^2
  #   ) + (
  #     2 * K.mat * (A.mat * (L - B) + B * (K.mat - 2 * L) - (A.mat + B)^2) * f.mat
  #   ) + (
  #     K.mat * (A.mat^3 + 2 * A.mat^2 * (B - L) - B * (B - L) * (K.mat + L) +
  #              A.mat * (B^2 - 2 * B * K.mat - 3 * B * L + L^2))
  #   )
  # )

  # time.new <- proc.time()
  # print(dc.ai.dF_vec_num[1:2, 1:2])
  # print(time.new - time.init)

  # time.init <- time.new

  # dc.ai.dF.vec <- (
  #                   (l.ivec %*% (R.jvec * (f.jvec^4 + 2 * (L - A.jvec) * f.jvec^4) + (L - A.mat)^2)) +

  #                   ((l.jvec * ))


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


  residuals <- c.final - data
  

  grad_derivs <- 10*exp(pars[1:num.kds]) * (exp(pars[1:num.kds]) + 1)^(-2) * 2* (colSums(residuals*dc.ai.dF) %*% t(dF.dK.mat[-nrow(dF.dK.mat),]) + rowSums(dc.ai.dKj_specific*residuals)[-nrow(dF.dK.mat)])
  

  dc.ai.dA <- sweep(sweep(dc.ai.dF,MARGIN=2,dF.base,"*") + dc.ai.dA_specific,
                    MARGIN=2,c.agos, "*")

  colnames(c.final) <- colnames(data)


  gradient_with_Ago <- log(10)*stock.ago*2*sum(residuals * dc.ai.dA)


  gradient_with_bg <- log(10)*B*2*sum(residuals * dc.ai.dB)

  return(c(grad_derivs, gradient_with_bg, gradient_with_Ago))
}

out <- matrix(c(pars, ModelLikelihood(pars)),nrow=1)
pars_update <- pars
colnames(out) <- c(rownames(data)[-length(rownames(data))],
                   "bg", "AGO", "-logp")


if (sitelist %in% c("12mers", "10mers")) {
    n <- 256
    colors <- c(rainbow_hcl(n, c=50, l=70, start = 0, end = 360*(n-1)/n),"red", "black")
    out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                       experiment, "/kds_PAPER/", start, "-", stop, "_", 
                       sitelist, "_", mirna.start, "-", mirna.stop,
                       "_singlebg_combinedinput_PAPER_cgfinal.txt")

  } else {
    colors = FALSE
    out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                       experiment, "/kds_PAPER/", start, "-", stop, "_", 
                       sitelist, 
                       "_singlebg_combinedinput_PAPER_cgfinal.txt")

  }



print(length(pars))

    solution <- optim(pars,
                    ModelLikelihood,
                    gr = Gradient,
                    method = "BFGS",
                    control = list(maxit=10000, parscale=100*Gradient(pars)^(-1)))



  pars <- solution$par
  sumoflogsquares_current <- solution$value
  print(sumoflogsquares_current)
  out <- rbind(out, c(pars, sumoflogsquares_current))
  print(out[,(ncol(out)-1)])
warnings()

  write.table(file=out_file, out, sep="\t", quote=FALSE, row.names=FALSE)
# }



