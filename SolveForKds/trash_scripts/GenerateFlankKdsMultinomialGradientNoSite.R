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

# Region within random sequence from which site types orginiates,
# going from position [26 - "start" : 26 + 37 + "stop"]
start <- as.integer(args[2])
stop  <- as.integer(args[3])
# Parameter specifying which list of sites is used for the analysis.
# I.E "Current" includes all the current site types, "canonical" is just
# the six seed-based 6mers, 7mers, and 8mers, along with no site., etc.

site <- args[4]

sitelist <- args[5]
if (sitelist %in% kmer_list_names) {
  mirna.start <- as.integer(args[6])
  mirna.stop <- as.integer(args[7])
} else {
  mirna.start <- NULL
  mirna.stop <- NULL
}

# # # Experiment name
experiment <- "equilibrium"
# print(experiment)

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

# Get siteXcounts, called sitesXcounts:
sitesXcounts <- GetSitesXCounts(mirna,
                                        experiment,
                                        start,
                                        stop,
                                        sitelist,
                                        mirna.start = mirna.start,
                                        mirna.stop = mirna.stop)

seqs <- sitesXcounts[, 1]
names(seqs) <- rownames(sitesXcounts)
sitesXcounts <- sitesXcounts[,-1]

# Remove any rows with no reads whoatsoever.
sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
# Remove any rows for which there are no input reads:
sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]

params <- GetKds(mirna, experiment, start, stop, sitelist, scaled = FALSE,nosite=TRUE)

num.kds <- nrow(sitesXcounts)
num.bgs <- 1

# Get vector of single site counts s.c.
s.c <- as.numeric(sitesXcounts[site, ])

sfXc_start <- GetSiteFlanksXCounts(mirna, experiment, site, start, stop, sitelist)
sfXc <- sfXc_start[,-2]

# print(sfXc)
colnames(sfXc)[1] <- "I"
sfXc <- sfXc[rowSums(sfXc[, 2:6]) > 0,]
sfXc <- sfXc[which(sfXc[,1] > 0),,drop=FALSE]

sfXc <- round(t((t(sfXc) / colSums(sfXc)) * s.c))
sfXc[is.na(sfXc)] <- 0
# Assign parameters for Kds and background, subracting 1 from rows due to "None"
# Kd assigned to 1 nM, and the fact that the I and A0 columns are not modeled
# in the experiment, so there's no background term assigned ot them.

# Get the number of flanking sites
num.sf <- nrow(sfXc)
# Get the site kds
kds.s <- params[1:num.kds]
kd.site <- kds.s[site]
# Omit the site kd for which the flanking sites are being fit.
kds.s <- kds.s[names(kds.s) != site]

# Assign the total site concentration in each experiment, and initialize the
# data.all column to be used.
# k.c.lib should be 100, for 100 nM.
# colors_sites <- site_cols[rownames(sitesXcounts)[rownames(sitesXcounts) != site],]

sitesandflanksXcounts <- rbind(sitesXcounts, sfXc)
sitesandflanksXcounts <- sitesandflanksXcounts[rownames(sitesandflanksXcounts) != site, ]



c.I.tots <- Norm(sitesandflanksXcounts[, 1]) * k.c.lib
names(c.I.tots) <- rownames(sitesandflanksXcounts)
# Define the vector of total site type concentrations:
c.totals <- matrix(c.I.tots, nrow=length(c.I.tots), ncol=5, byrow=FALSE)
# Remove the I and A0 columns from the data.all to be fit to the model. 
data.all <- sitesandflanksXcounts[, 2:6]
rownames(c.totals) <- rownames(data.all)
colnames(c.totals) <- colnames(data.all)
# colors_flanks <- sapply(rownames(sfXc), GetColorFunction)
# colors_all <- c(colors_sites, colors_flanks)
# pch_all <- c(rep(1, length(colors_sites)), rep(19, nrow(sfXc)))

c.agos <- sapply(colnames(data.all), function(x) {
   as.numeric(x) / 100
  })
# Initialize starting Kds, which are set to 1/the enrichment of each site type
# in the A12.6 experiment.
site_colors = sample(colors(),nrow(data.all))
pars.init <- c(kds.s,rep(kd.site,nrow(sfXc)),params[(num.kds + 1): (num.kds + num.bgs + 1)])
# pars.init <- rep(0,length(pars.init))
bgs <- 10^params[num.kds+1]
stock.ago <- 10^params[num.kds+2]
B <- bgs

pars.init <- rep(kd.site,nrow(sfXc))

names(pars.init) <- rownames(sfXc)
pars <- pars.init
print(pars.init)




# Define function of just kds and pars.

num.kds.all <- num.kds + nrow(sfXc) - 1

out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                   experiment, "/kds_PAPER/", site, "_", start, "-", stop, "_", 
                   sitelist, "_flanking_singlebg_nosite_multinomial_PAPER.txt")

out_last_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                        experiment, "/kds_PAPER/", site, "_", start, "-", stop, "_", 
                        sitelist, "_flanking_singlebg_nosite_multinomial_PAPER_last.txt")

specific_figure_string <- paste0(site, "_", start, "-", stop, "_", sitelist,
                                 "_flanking_singlebg_nosite_multinomial_PAPER")





ModelFunction <- function(pars) {
  # Split up the parameters into the kd and background parameters.
  kds  <- Logistic(c(kds.s,pars),10)

  names(kds) <- rownames(data.all)
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
  c.final <- data.frame(t(t(c.all) / colSums(c.all) * colSums(data.all)))
  return(c.final)
}

tick <- 0
ModelLikelihood <- function(pars, print.=FALSE){

  model <- ModelFunction(pars)
  model_norm <- t(t(model) / colSums(model))
  sumofsquares <- -sum(data.all*log(model_norm))
  if (tick > 0) {
    out <<- rbind(out, c(pars, sumofsquares))  
  }

  # if (tick %% 10 == 0 & tick > 1 & print.==TRUE) {
  #   dev.set(3)
  #   plot(c(0.0001,0.0002),
  #        c(0.0001,0.0002),
  #        col="white",
  #        xlim=10^c(-3,10),
  #        ylim=10^c(0,10),
  #        log='xy',
  #        ann = FALSE,
  #        axes = FALSE)
  #   title(main=mirna, font.main = 1)
  #   segments(1,1, x1=1,10^7, lty = 2)
  #   segments(1, 1, x1=10^7, 10^7, lty = 3)

  #   axis(1,sapply(10^(-3:8), function(x) {x * c(1,2,3,4,5,6,7,8,9)}), labels = FALSE, pos = 1, lwd = 2)
  #   axis(2,sapply(10^(0:8), function(x) {x * c(1,2,3,4,5,6,7,8,9)}),labels = FALSE, pos = 10^-3,lwd = 2)
  #   axis(1,10^(-3:8), pos = 1, lwd = 0)
  #   axis(2,10^(0:8), pos = 10^-3,lwd = 0)

  #   sapply(1:ncol(data.all), function(col) {
  #     x_model <- model[,col]
  #     y_data.all <- data.all[, col]
  #     x_plot <- x_model/sum(x_model)*sum(y_data.all)
  #     y_plot <- y_data.all
  #      points(x_plot, y_plot,col=c(site_cols[rownames(data.all)[1:length(kds.s)],], colors_flanks), pch=20)
  #   })
  #   print(sumofsquares)
  # }
  tick <<- tick + 1
  model <<- model
  return(sumofsquares)
}

Gradient <- function(pars,print.=FALSE) {
  # Split up the parameters into the kd and background parameters.
  kds  <- Logistic(c(kds.s, pars), 10)
  names(kds) <- rownames(data.all)
  B <- bgs[1]
   f.mat <- matrix(
      rep(sapply(c.agos, function(ago.percent) {
           return(GetFreeAgo(kds, c.I.tots, ago.percent * stock.ago))
         }
         ), nrow(data.all)), nrow=nrow(data.all), ncol=ncol(data.all), byrow=TRUE)

  f.jvec <- sapply(c.agos, function(ago.percent) {
           return(GetFreeAgo(kds, c.I.tots, ago.percent * stock.ago))})

  A.mat <- matrix(c.agos*stock.ago, nrow=nrow(data.all), ncol=ncol(data.all),
                  byrow=TRUE)

  A.jvec <- c.agos * stock.ago

  K.mat <- matrix(kds, nrow=nrow(data.all), ncol=ncol(data.all), byrow=FALSE)
  K.ivec <- kds

  l.mat <- matrix(c.I.tots, nrow=nrow(data.all), ncol=ncol(data.all), byrow=FALSE)
  l.ivec <- c.I.tots

  R.mat <- matrix(colSums(data.all), nrow=nrow(data.all), ncol=ncol(data.all), byrow=TRUE)
  R.jvec <- colSums(data.all) 
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

  residuals <- -data.all/c.final
  grad_derivs <- (10*exp(pars) * (exp(pars) + 1)^(-2) *
                  (colSums(residuals*dc.ai.dF) %*%
                   t(dF.dK.mat[(length(kds.s)+1):num.kds.all,]
                  ) + rowSums((dc.ai.dKj_specific*residuals)[(length(kds.s)+1):num.kds.all,])))
  
  # dc.ai.dA <- sweep(sweep(dc.ai.dF,MARGIN=2,dF.base,"*") + dc.ai.dA_specific,
  #                   MARGIN=2,c.agos, "*")

  colnames(c.final) <- colnames(data.all)

  # gradient_with_Ago <- log(10)*stock.ago*sum(residuals * dc.ai.dA)

  # gradient_with_bg <- log(10)*B*sum(residuals * dc.ai.dB)
  gradient_all <- c(grad_derivs)
  names(gradient_all) <- names(pars)
  return(gradient_all)
}


out <- matrix(c(pars, ModelLikelihood(pars)), nrow=1)
out_names <- c(names(pars), "-logp")
colnames(out) <- out_names



solution <- optim(pars,
                  ModelLikelihood,
                  gr = Gradient,
                  print.=TRUE,
                  method = "L-BFGS-B",
                  lower=rep(-13, length=length(pars)),
                  upper=rep(10, length=length(pars)),
                  control = list(maxit=100000, factr=1e2))

pars <- solution$par
sumofsquares <- solution$value
out <- rbind(out, c(pars, sumofsquares))

# params_site_new <- Logit(10^mean(log(Logistic(pars[rownames(sfXc)],10))),10)
# params_new <- params
# params_new[names(params)] <- pars[names(params)]

# params_new[site] <- params_site_new
# kds_old <- Logistic(params[1:num.kds],10)
# bgs_old <- 10^params[num.kds+1]
# ago_old <- 10^params[num.kds+2]
# final_old <- c(kds_old, bgs_old, ago_old)
# kds_new <- Logistic(params_new[1:num.kds],10)
# bgs_new <- 10^params_new[num.kds+1]
# ago_new <- 10^params_new[num.kds+2]
# final_new <- c(kds_new, bgs_new, ago_new)

# plot(final_old, final_new,col=c(site_cols[names(final_old[-length(final_old)]),],"brown"),log='xy')

# segments(1e-10, 1e-10, 100, 100, lty=2)
print(out_file)
print(out_last_file)

write.table(file=out_file, out, sep="\t", quote=FALSE, row.names=FALSE)
write.table(file=out_last_file, out[nrow(out),,drop=FALSE], sep="\t",
                quote=FALSE, row.names=FALSE)


