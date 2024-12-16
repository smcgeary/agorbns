
graphics.off()
source("general/general.R")
source("general/GenericFigures.R")

MakeTargetPredictionFigure <- function(pdf.plot=FALSE, width=8, height=7) {
	data_pred <- read.table("EMBO_talk_data/TargetPredictionAlgorithms.txt")
  # colnames(dG.df) <- gsub("(^.*)\\.(.*$)", colnames(dG.df), replace="\\1-\\2")
  # R <- 1.987e-3 # in kcal K-1 mol-1
  # T <- 310.15 # in K
  # total_x <- c(1)
  # total_y <- c(1)
  SubfunctionCall(FigureSaveFile)
  xmin <- -0.05
  xmax <- 0.15
  ymin <- 0
  ymax <- nrow(data_pred) + 1
  BlankPlot(inv="y")
  xmin <- 0
  # ymin <- -2.5
  # mirna_labs <- kMirnas
  # mirna_labs[length(mirna_labs)] <- "miR-7"
  print(seq(nrow(data_pred)))
  AddLinearAxis(2, alt_lab=data_pred[, 1], alt_lab_pos=seq(nrow(data_pred)),
                alt_lab_y_dist=-1, line=1,
                label="", angled=FALSE, noline=TRUE, xpd=NA)

  AddLinearAxis(1, tick.space=0.05, label.space=0.05,
                label=expression(italic(r)^2))

  x_start <- 0
  rect(xleft=rep(0, nrow(data_pred)), ybottom=seq(nrow(data_pred))+0.4,
					 xright=data_pred[, 3], ytop=seq(nrow(data_pred))-0.4,
           lwd=0.25, xpd=NA, col="black")

  break
  text(rep(length(kMirnas)*4 + 0.75, 3), -R*T*log(c(2, 10, 50)),
       labels=c("2-fold greater\nbinding affinity",
                "10-fold greater\nbinding affinity",
                "50-fold greater\nbinding affinity"),
       cex=0.8, adj=c(0, 0.5), xpd=NA)
  kds.all <- sapply(kMirnas, function(mirna) {
    if (mirna == "miR-7-23nt") {
      combined <- FALSE
      experiment <- "equilibrium2_nb"
    }
    if (mirna == "miR-1") {
      combined <- FALSE
      buffer <- TRUE
    }
    fixed <- FALSE
    altfixed <- FALSE
    kds <- SubfunctionCall(EquilPars)[sprintf("%s_Kd", kSeedSites), ]$Mean
    names(kds) <- kSeedSites
    rel_kds <- R*T*log(kds[c("7mer-A1", "8mer", "7mer-m8", "8mer")]/
                kds[c("6mer", "7mer-m8", "6mer", "7mer-A1")])
    x_l <- seq(4) - 1 + x_start
    x_r <- x_l + 1
    y_b <- c(0, 0, 0, 0)
    y_t <- rel_kds
    cols_full <- kSiteColors[c("7mer-A1", "7mer-m8")]
    cols_alpha <- ConvertRColortoRGB(kSiteColors[c("7mer-A1", "7mer-m8")],
                                     alpha=0.3)
    cols_alpha <- c("turquoise1", "plum1")
    cols <- c(cols_full[1], cols_alpha[1], cols_full[2], cols_alpha[2])
    rect(x_l, y_b, x_r, y_t, col=cols, border=NA)    
    x_start <<- x_start + 4
    rel_kds
  })
  # segments(x0=seq(from=0.5, to=length(kMirnas) - 0.5, by=0.5)*4,
  #          y0=c(rep(c(0, 0.1), length(kMirnas) - 1), 0),
  #          x1=seq(from=0.5, to=length(kMirnas) - 0.5, by=0.5)*4, y1=ymin,
  #          lty=c(2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1), lwd=0.5, xpd=NA)
  x <- kds.all[c(1, 3),]
  y <- kds.all[c(2, 4), ]
  xy <- GetPlotFractionalCoords(fx=0.85, fy=1.1, inv='y')
  AddCorrelationToPlot(c(x), c(y), xpos=xy[1], ypos=xy[2], adj=0,
                       rsquared=TRUE)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

MakeTargetScanFigure <- function(pdf.plot=FALSE) {
	path <- file.path("/lab/bartel4_ata/kathyl/RNA_Seq/analysis/data",
	                  "TS7_output_files/recreate_ts7_predictions.txt")
	data_pred <- read.table(path, header=TRUE)
  SubfunctionCall(FigureSaveFile)
  xmin <- -2
  xmax <- 0
  ymin <- -2
  ymax <- 2
  BlankPlot(inv="xy")
  # ymin <- -2.5
  # mirna_labs <- kMirnas
  # mirna_labs[length(mirna_labs)] <- "miR-7"
  print(seq(nrow(data_pred)))
  AddLinearAxis(1, tick.space=1, label.space=1, label="Predicted log2fc")
	AddLinearAxis(2, tick.space=1, label.space=1, label="Observed log2fc")

  x <- data_pred[, 3]
  y <- data_pred[, 4]

	Points(x, y, col=ConvertRColortoRGB("black", alpha=0.2))
  xy <- GetPlotFractionalCoords(fx=0.95, fy=0.95, inv='xy')
  AddCorrelationToPlot(c(x), c(y), xpos=xy[1], ypos=xy[2], adj=1,
                       rsquared=TRUE)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

# MakeTargetPredictionFigure()

path_kl <- file.path("/lab/bartel4_ata/kathyl/RNA_Seq/analysis/data",
	                   "TS7_output_files/recreate_ts7_predictions.txt")
path_va <- file.path("EMBO_talk_data/elife-05005-supp1-v1.txt")
data_kl <- read.table(path_kl, header=TRUE, sep="\t", stringsAsFactors=FALSE)
seeds <- unique(data_kl$Seed)
print(seeds)

sapply(seeds, function(seed) {
	message(seed)
	mRNAs <- length(which(data_kl$Seed == seed))
	message(sprintf("Number of genes: %s", round(mRNAs, 1)))
})



data_va <- read.table(path_va, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# MakeTargetScanFigure()

