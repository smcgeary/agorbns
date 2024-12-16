source("general/general.R")




GetKineticsData <- function(mirna, experiment="kinetics", n_constant=5,
                            sitelist="paper", kDilRatio=10, kSubSet=FALSE) {
  sitesXcounts <- SitesXCounts(mirna, "equilibrium", n_constant, sitelist)
  print(sitesXcounts)
  sitesXcounts.kinetics <- SitesXCountsKinetics(mirna, experiment, n_constant,
                                              sitelist)
  print(sitesXcounts.kinetics)
  sitesXcounts.p <- sitesXcounts.kinetics[[1]]
  sitesXcounts.c <- sitesXcounts.kinetics[[2]]
  # 2. Load the Kds from the equilibrium fits
  params.e <- GetSiteKds(mirna, "equilibrium", n_constant, sitelist)
  # I never use the sequences, realistically
  seqs <- sitesXcounts[,1]
  # Subset the equilibrium and kinetic data matrices, such that the equilibrium
  # Matrix includes the combined input through to the zero protein, and the
  # kinetics matrices include the pulse through to the equilibrium sample.
  sitesXcounts <- sitesXcounts[, c(-1, -2)]
  indeces.to.use <- which(rowSums(sitesXcounts) > 0 & sitesXcounts[, 1] > 0)
  print(indeces.to.use)
  kData.e <- sitesXcounts[indeces.to.use, 1 : 7]
  kData.p <- sitesXcounts.p[indeces.to.use, 4 : (ncol(sitesXcounts.p) - 2)]
  kData.c <- sitesXcounts.c[indeces.to.use, 4 : (ncol(sitesXcounts.c) - 2)]
  kData.c[, 1] <- 0
  # Conditionally subset the data if flag is not FALSE:
  if (kSubSet != FALSE) {
    print(kSubSet)
    kData.e.sites <- kData.e[1:kSubSet, ]
    kData.e.none  <- colSums(kData.e[(kSubSet + 1):nrow(kData.e), ])
    kData.e       <- rbind(kData.e.sites, kData.e.none)
    rownames(kData.e)[kSubSet + 1] <- "None"
    kData.p.sites <- kData.p[1:kSubSet, ]
    kData.p.none  <- colSums(kData.p[(kSubSet + 1):nrow(kData.p), ])
    kData.p       <- rbind(kData.p.sites, kData.p.none)
    rownames(kData.p)[kSubSet + 1] <- "None"
    kData.c.sites <- kData.c[1:kSubSet, ]
    kData.c.none  <- colSums(kData.c[(kSubSet + 1):nrow(kData.c), ])
    kData.c       <- rbind(kData.c.sites, kData.c.none)
    rownames(kData.c)[kSubSet + 1] <- "None"
  }
  if (mirna == "let-7a") {
  colnames(kData.p)[2:5] <- c("2.1", "5.1", "2.2", "5.2")
  colnames(kData.c) <- colnames(kData.p)

  }
  times <- as.integer(floor(as.numeric(colnames(kData.p))))
  kData <- sapply(unique(times), GetAverageOfReplicates,
                   times=times, data=rbind(kData.p, kData.c))
  print("Kdata:")
  print(kData)
  print("done kData")
  rownames(kData) <- c(sapply(c(".p", ".c"), function(i) {
    sapply(rownames(kData.p), function(site) {
      paste0(site, i , collapse = "")
      })
    }))
  kTimes <- unique(times)/60
  colnames(kData) <- kTimes
  kInputPulseReads <- sitesXcounts.p[indeces.to.use, c("I_combined")]
kInputChaseReads <- sitesXcounts.c[indeces.to.use, c("I_combined")]
if (kSubSet != FALSE) {
  kInputPulseReads.sites <- kInputPulseReads[1:kSubSet]
  kInputPulseReads.none  <- sum(kInputPulseReads[(kSubSet + 1):length(kInputPulseReads)])
  kInputPulseReads       <- c(kInputPulseReads.sites, kInputPulseReads.none)
  kInputChaseReads.sites <- kInputChaseReads[1:kSubSet]
  kInputChaseReads.none  <- sum(kInputChaseReads[(kSubSet + 1):length(kInputChaseReads)])
  kInputChaseReads       <- c(kInputChaseReads.sites, kInputChaseReads.none)
}
# Fit the ratio of pulse to chase, in the native library,
# to correct for the actual difference in concentration from 1:1 in the
# experiment.
kPulseChaseRatio <- (sum(sitesXcounts.c[indeces.to.use,c("I")]) /
                     sum(sitesXcounts.p[indeces.to.use,c("I")]))
kInputPulseConc <- Norm(kInputPulseReads) * kLibraryConcInRxn
kInputChaseConc <- Norm(kInputChaseReads) * kLibraryConcInRxn * kPulseChaseRatio
names(kInputPulseConc) <- rownames(kData.p)
names(kInputChaseConc) <- rownames(kData.c)
# Constants throughout the optimization routine:
kAgoStockConc <- params.e["AGO","Mean"]
kIndsKDs        <<- 1:kNumSites
kIndsKoffs      <<- (kNumSites + 1):(2 * kNumSites)
kIndsContamKD   <<- 2 * kNumSites + 1
kIndsContamKoff <<- 2 * kNumSites + 2
kIndsAgo        <<- 2 * kNumSites + 3
kIndsContam     <<- 2 * kNumSites + 4
kIndsBgs        <<- (2 * kNumSites + 5):(2 * kNumSites + 16)
# Make the matrix to multiply through the adjoint equations:
kFDotXMatrix <- diag(x=1,nrow=4*kNumSites+2)
kFDotXMatrix[nrow(kFDotXMatrix) -1, 1:(2 * kNumSites)] <- 1
kFDotXMatrix[nrow(kFDotXMatrix), (2 * kNumSites + 1):(4 * kNumSites)] <- 1
# Define the total concentration in the INITIAL BINDING (kInputInitial), and in the
# chase (L).
kInputInitial <- c(kInputPulseConc, kInputChaseConc * 0)
kInput <- ((kInputInitial + c(kInputPulseConc * 0, kInputChaseConc) * kDilRatio)
           / (kDilRatio + 1))
kInputMatrix <- matrix(c(kInputInitial / (1 + kDilRatio), 
                         rep(kInput, kNumBgs - 1)),
                       nrow=length(kInput),
                       ncol=kNumBgs,
                       byrow=FALSE)
rownames(kInputMatrix) <- rownames(kData)
colnames(kInputMatrix) <- colnames(kData)
kInputTotals <- colSums(kInputMatrix)
  return(list(kData,kInputMatrix))
}

GetBgs <- function(pars, kNumSites) {
  kIndsBgs        <- (2 * kNumSites + 5):(2 * kNumSites + 16)
  return(10^pars[kIndsBgs])
}

MakeODEPars <- function(pars, kNumSites, l. = kInput) {
  pars <- 10^pars
  kIndsKDs        <- 1:kNumSites
  kIndsKoffs      <- (kNumSites + 1):(2 * kNumSites)
  kIndsContamKD   <- 2 * kNumSites + 1
  kIndsContamKoff <- 2 * kNumSites + 2
  # The site type on and off rates, the background terms for each column,
  # the and the total concentration of Ago (a) and the contaminant (b):
  print(kIndsKDs)
  print(kIndsKoffs)
  kds.a.   <- pars[kIndsKDs]
  koffs.a. <- pars[kIndsKoffs]
  kd.b     <- pars[kIndsContamKD]
  koff.b   <- pars[kIndsContamKoff]
  # Assignment of the Kds (being the ratio of the on and off rates)
  kons.a. <- koffs.a. / kds.a.
  kon.b   <- koff.b / kd.b
  # Removal of Kd parameter as it is no longer useful:
  # Calculate the vector of diluted, bound RNA:
  # Assign the parameter vector and initial conditions vector for the ODE solver:
  parms.k <- c(koffs.a., kons.a., l., koff.b, kon.b)
  return(parms.k)
}

MakeX0 <-function(pars, kNumSites, l0. = kInputInitial) {
  pars <- 10^pars
  kIndsKDs        <- 1:kNumSites
  kIndsKoffs      <- (kNumSites + 1):(2 * kNumSites)
  kIndsContamKD   <- 2 * kNumSites + 1
  kIndsContamKoff <- 2 * kNumSites + 2
  kIndsAgo        <- 2 * kNumSites + 3
  kIndsContam     <- 2 * kNumSites + 4

  # The site type on and off rates, the background terms for each column,
  # the and the total concentration of Ago (a) and the contaminant (b):
  kds.a.   <- pars[kIndsKDs]
  koffs.a. <- pars[kIndsKoffs]
  kd.b     <- pars[kIndsContamKD]
  koff.b   <- pars[kIndsContamKoff]
  A0       <- 0.4 * pars[kIndsAgo]
  B0       <- 0.4 * pars[kIndsContam]
  # Assignment of the Kds (being the ratio of the on and off rates)
  kons.a. <- koffs.a. / kds.a.
  kon.b   <- koff.b / kd.b
  # print(kons.a.)
  # print("kons.a.")
  # print(kon.b)
  # print("kon.b")
  # print(kds.a.)
  # print("kds.a.")
  # print(koffs.a.)
  # print("koffs.a.")
  # print(kd.b)
  # print("kd.b")
  # print(A0)
  # print("A0")
  # print(B0)
  # print("B0")
  # print(l0.)
  # print("l0.")
  # Get the free amount of Ago and contaminant:
  a0.b0. <- GetFreeAgoAndContam(kds.a., kd.b, l0., A0, B0)
  # print("a0.b0.")
  # print(a0.b0.)
  a0 <- a0.b0.[1]
  b0 <- a0.b0.[2]
  a0 <<- a0
  b0 <<- b0
  l0. <<- l0.
  # print("a0")
  # print(a0)
  # print("b0")
  # print(b0)

  # Calculate the initial Ago and contaminant occupancies:
  a.occs0.b.occs0 <- GetOccupanciesContaminant(a0, b0, kds.a., kd.b)
  xap0.xac0. <- l0. * a.occs0.b.occs0[[1]]
  xbp0.xbc0. <- l0. * a.occs0.b.occs0[[2]]
  # print("xap0.xac0.")
  # print(xap0.xac0.)
  # print("xbp0.xbc0.")
  # print(xbp0.xbc0.)
  # Removal of Kd parameter as it is no longer useful:
  # Calculate the vector of diluted, bound RNA:
  # Assign the parameter vector and initial conditions vector for the ODE solver:
  x0. <- c(xap0.xac0., xbp0.xbc0., a0.b0.) / (kDilRatio + 1)
  return(x0.)
}

MakeXTimeCourse <- function(x0, parms.k, kNumSites, kCScriptDir, kCScriptName, kTimesAll, kTimes, verbose.=FALSE) {
  # Load and run the ODE:
  dyn.load(kCScriptDir)
  # print(x0)
  # print(parms.k)
  x.t <- t(ode(y          = x0,
               times      = kTimesAll,
               func       = "ode_deriv",
               parms      = parms.k,
               method     = "lsodes",
               dllname    = kCScriptName,
               sparsetype = "sparseint",
               initfunc   = "ode_p_init",
               nout       = 1,
               verbose    = verbose.)[,2:(4*kNumSites+3)])
  dyn.unload(kCScriptDir)
  x <- x.t[, which(kTimesAll %in% kTimes)]
  print(x)
  print(kTimes)
  colnames(x) <- kTimes
  # Calculate total bound pulse and chase sites:
  return(x)
}


MakeModelPrediction <- function(x, bgs, data = kData, l = kInputMatrix) {
  # Combine the ago-bound and contaminant-bound site types:
  # Written out the numerator for maximum accuracy:
  # Form of equation is :        x(L - X) + B(l - x)
  #                          D * –––––––––––––––––––
  #                                (L - X)(X + B)
  kNumSites = nrow(kData)/2
  # print(dim(x))
  # print(dim(l))
  # print(length(bgs))
  # print(dim(data))
  x.ago <- x[1:(2 * kNumSites), ] 
  x.con <- x[(2 * kNumSites + 1):(4 * kNumSites), ]
  x     <- x.ago + x.con
  X     <- colSums(x)
  L     <- colSums(l)
  D     <- colSums(kData)
  B     <- bgs
  B <<- B
  L <<- L
  X <<- X
  # Transpose D-multiply x and l for row-multiplication:
  Dx.     <- D * t(x)
  Dl.     <- D * t(l)
  Dl. <<- Dl.
  Dx. <<- Dx.
  # print(Dx.)
  # print(Dl.)
  # print(B)
  pred.  <- (Dx.*L - Dx.*X + Dl.*B - Dx.*B) /
                (L*B + L*X -  B*X - X^2)
  pred   <- t(pred.)
  return(pred)
}


GetKineticsModelGeneral <- function(mirna, experiment, n_constant, sitelist, costfunction,dil=FALSE,pars=FALSE) {
  if (pars == FALSE) {
    pars <- GetSiteKoffs(mirna, experiment, n_constant, sitelist, costfunction, dil=dil)
  }
  if (dil == TRUE) {
    kDilRatio <- 10^pars[length(pars)-12]
    print(kDilRatio)
  } else {
    kDilRatio <- 11
  }
  kData_kInput <- GetKineticsData(mirna, experiment, n_constant, sitelist, kDilRatio=11)
  kData <- kData_kInput[[1]]
  kNumSites <- nrow(kData)/2
  print(kNumSites)
  # system(paste0("python SolveForOffRates/",
  #             "MakeCScriptSingleExponential.py ",
  #             kNumSites))
  # # Performa a system command to compile both of these scripts.
  # system(paste0("R CMD SHLIB SolveForOffRates/SingleODEContam_",
  #               kNumSites, ".c"))

  kInputMatrix <- kData_kInput[[2]]

# Assign names to their name and directory to be used within the optimization
# routine.
  kCScriptName <- paste0("SingleODEContam_", kNumSites)
  kCScriptDir <- paste0("SolveForOffRates/", kCScriptName, ".so")
  try(dyn.unload(kCScriptDir))
  print("hi")
  print(names(pars))
  temp_names <- names(pars)
  pars <- as.numeric(pars)

  names(pars) <- temp_names
  if (dil != FALSE) {
    kDilRatio <- 10^pars[2 * kNumSites + 5]
    pars <- pars[-(2*kNumSites + 5)]
  } else {
    kDilRatio <- 11
  }
  x0      <- MakeX0(pars, kNumSites, l0. = (kDilRatio+1) * kInputMatrix[,1])
  x0 <<- x0
  print("done x0")
  # print(x0)
  parms.k <- MakeODEPars(pars, kNumSites, l. = kInputMatrix[, 2])
  print("done parms.k")
  parms.k <<- parms.k
  # print(parms.k)
  x       <- MakeXTimeCourse(x0, parms.k, kNumSites, kCScriptDir, kCScriptName,kTimesAll, kTimes)
  print("done x")
  x <<- x
  kData <<- kData
  bgs     <- GetBgs(pars, kNumSites)
  bgs <<- bgs
  model     <- MakeModelPrediction(x, bgs, data = kData, l = kInputMatrix)
  print("done model")
  return(model)
}

PlotKineticsFit <- function(mirna, num_sites, sitelist, costfunc, dil=FALSE, col=FALSE) {
  kinetics.data <- GetKineticsData(mirna, "kinetics", num_sites, sitelist,kDilRatio=11)
  kData <- kinetics.data[[1]]
  kInputMatrix <- kinetics.data[[2]]
  data.e <- SitesXCounts(mirna, "equilibrium", 5, "paper")
  colors <- kSiteColors[rownames(data.e), ]
  model <- GetKineticsModelGeneral(mirna, "kinetics", 5, "paper", costfunc, dil=dil)
  model <<- model
  par(kPlotParameters)
  model.norm <- t(t(model + 1) / colSums(model + 1))
  data.norm <- t(t(kData + 1) / colSums(kData + 1))
  if (col != FALSE) {
    model.norm <- model.norm[, col, drop=FALSE]
    data.norm <- data.norm[, col, drop=FALSE]
  }
  plot(c(model.norm), c(data.norm),
       xlim = c(1e-8, 0.8),
       ylim = c(1e-8, 0.8),
       log  = 'xy',
       col  = colors,
       pch  = rep(c(20,1),
       each = kNumSites),
       cex = 1) 
    abline(0, 1, lty = 3)
}

PlotKineticsFitInidividual <- function(mirna, n_constant, sitelist, costfunc, log. = "", dil=FALSE) {
  site_list <- c(read.table(file = paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/AssignSiteTypes/sites.", mirna, "_", sitelist,".txt"),
    stringsAsFactors=FALSE)[,1], "None")
  print(site_list)
  kinetics.data <- GetKineticsData(mirna, "kinetics", n_constant, sitelist, kDilRatio=11)
  kData <- kinetics.data[[1]]
  kInputMatrix <- kinetics.data[[2]]
  data.e <- SitesXCounts(mirna, "equilibrium", 5, sitelist)
  model <- GetKineticsModelGeneral(mirna, "kinetics", 5, sitelist, costfunc, dil = dil)
  model <<- model
  koffs <- GetSiteKoffs(mirna, "kinetics", n_constant, sitelist, costfunc, dil=dil)
  print(koffs)
  koffs <- 10^koffs[(nrow(kData)/2+1):(nrow(kData))]
  names(koffs) <- rownames(data.e)
  print(koffs)

  times <- as.numeric(colnames(model))
  print(times)
  if (sitelist == "canonical") {
    dev.new(xpos = 20, ypos = 20, height = 7, width = 10)
    par(mfrow = c(2, 4))
  } else if (mirna == "miR-124") {

    dev.new(xpos = 20, ypos = 20, height = 8, width = 14)
    par(mfrow = c(4, 6))  
  } else {
    dev.new(xpos = 20, ypos = 20, height = 8, width = 12)
    par(mfrow = c(4, 5))  
  }
  par(kPlotParameters)

  kInput.Pulse <- Norm(kInputMatrix[1:(nrow(kData)/2),1])
  kInput.Chase <- Norm(kInputMatrix[(nrow(kData)/2 + 1):nrow(kData),2])

  kInput.Pulse <<- kInput.Pulse
  kInput.Chase <<- kInput.Chase
  kInput.Pulse <- 1
  kInput.Chase <- 1
  for (site in site_list) {
    print(site)
    color <- kSiteColors[site, ]
    times[1] <- 0.5/60
    site.inds <- grep(paste0(site,"."), rownames(model),fixed = TRUE)
    print(rownames(model))
    print(site.inds)
    model.p <- model[site.inds[1], ]/colSums(model)
    model.c <- model[site.inds[2], ]/colSums(model)
    data.p <- kData[site.inds[1], ]/colSums(kData)
    data.c <- kData[site.inds[2], ]/colSums(kData)
    # model.c[1] <- 1
    # data.c[1] <- 1

    print(model.p)
    print(model.c)
    print(data.p)
    print(data.c)

    ylim. = c(min(model.p, model.c[-1], data.p, data.c[-1]),
             max(model.p, model.c, data.p, data.c))
    # ylim. <- c(0.1, max(model.p, model.c, data.p, data.c)*2)
    print(ylim.)
    plot(times, data.p,
           xlim = c(0.3, 20000)/60,
           ylim = ylim.,
           log  = log.,
           col  = color,
           pch  = 20) 
    points(times, data.c, pch = 1,col=color)
    lines(times, model.p, col=color)
    lines(times, model.c, lty = 2, col=color)
    mtext(site, 3, line = 0.5, adj=0.4)
    print(grep(paste0(site),names(koffs),fixed=TRUE))
    dwell = 1/(koffs[site])
    segments(dwell,ylim.[1],dwell,ylim.[2], lty=2,col=color)
  }
    plot(1, type = "n", ann=FALSE, axes=FALSE)
  title(main = costfunc, cex = 2)

}