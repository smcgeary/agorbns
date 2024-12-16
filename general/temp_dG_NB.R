PlotSiteKdsVsSPS <- function(experiment="equilibrium", n_constant=5,
                             sitelist="paperfinal", combined=TRUE, buffer=FALSE,
                             dG.table=3, height=8.5, pdf.plot=FALSE) {
  dG.df <- read.table(file=paste0("canonical_sites_mfe_", dG.table, ".txt"))
  colnames(dG.df) <- gsub("(^.*)\\.(.*$)", colnames(dG.df), replace="\\1-\\2")
  R <- 1.987e-3 # in kcal K-1 mol-1
  T <- 310.15 # in K
  SubfunctionCall(FigureSaveFile)
  ymin <- 1e-3
  ymax <- 1e6
  xmin <- -12
  xmax <- -2
  kSeedSites <- c("6mer", "7mer-m8")
  BlankPlot(log="y", inv="xy")
  x.line <- seq(xmin, xmax, length=20)
  xmax.alt <- log(1e6/1e9)*R*T
  x.line.alt <- seq(xmin, xmax.alt, length=20)
  y.line <- exp(x.line.alt/(R*T))*1e9
  lines(x.line.alt, y.line, col="gray80")

  AddLinearAxis(1, tick.space=1, label.space=2,
                label=expression("Predicted"~Delta*italic(G)~"(kcal/mol)"))
  AddLogAxis(2, label="Kd (nM)")
  Kd.matrix <- sapply(kMirnas, function(mirna) {
    if (mirna == "miR-7-23nt") {
      combined <- FALSE
      experiment <- "equilibrium2_nb"
    }
    if (mirna == "miR-1") {
      combined <- FALSE
      buffer <- TRUE
    }
    kd <- SubfunctionCall(EquilPars)
    out <- kd[paste0(kSeedSites, "_Kd"), ]$Mean/kd["None_Kd", ]$Mean*10
    out
  })
  dG.matrix <- sapply(kMirnas, function(mirna) {
    if (mirna == "miR-7-23nt") {
      mirna <- "miR-7"
    }
    as.numeric(dG.df[kSeedSites, mirna])
  })
  rownames(dG.matrix) <- kSeedSites
  rownames(Kd.matrix) <- kSeedSites
  cols.mirnas <- ConvertRColortoRGB(kMirnaColors[kMirnas], alpha=1)
  cols.sites <- kSiteColors[kSeedSites]
  x1_out <- c()
  x2_out <- c()
  y1_out <- c()
  y2_out <- c()

  sapply(1:ncol(dG.matrix), function(i_c) {
    x <- dG.matrix[, i_c]
    y <- Kd.matrix[, i_c]
    x1_out <<- c(x1_out, x[1])
    x2_out <<- c(x2_out, x[2])
    y1_out <<- c(y1_out, y[1])
    y2_out <<- c(y2_out, y[2])
  })
  linmodel1 <- lm(log10(y1_out)~x1_out)
  linmodel2 <- lm(log10(y2_out)~x2_out)
  b1 <- linmodel1$coefficients[1]
  m1 <- linmodel1$coefficients[2]
  b2 <- linmodel2$coefficients[1]
  m2 <- linmodel2$coefficients[2]
  ypoints1 <- 10^(m1*x.line + b1)
  ypoints2 <- 10^(m2*x.line + b2)

  # global_df <<- data.frame(Kd_6=log10(y1_out), dG_6=x1_out,
  #                         Kd_7=log10(y2_out), dG_7=x2_out)
  print("6mer test")
  print(summary(linmodel1))

  coefs1 <- summary(linmodel1)$coefficients[2, 1:2]
  coefs2 <- summary(linmodel2)$coefficients[2, 1:2]


  print("7mer-m8 test")
  print(summary(linmodel2))

  RT_coef <- 1/(R*T)/log(10)

  prob_RT_1 <- 2*pt((coefs1[1] - RT_coef)/coefs1[2], df=length(kMirnas) - 2)
  prob_RT_2 <- 2*pt((coefs2[1] - RT_coef)/coefs2[2], df=length(kMirnas) - 2)

  message("probability 6mer is RTlnk:")
  message(prob_RT_1)
  message("probability 7mer-m8 is RTlnk:")
  message(prob_RT_2)


  predict1 <- predict(linmodel1, data.frame(x1_out=x.line), interval = 'confidence')
  polygon(c(x.line, rev(x.line)), c(10^predict1[,2], rev(10^predict1[, 3])), col=rgb(0, 0, 0, alpha=0.1), border=NA)
  predict2 <- predict(linmodel2, data.frame(x2_out=x.line), interval = 'confidence')
  polygon(c(x.line, rev(x.line)), c(10^predict2[,2], rev(10^predict2[, 3])), col=rgb(0, 0, 0, alpha=0.05), border=NA)
  lines(x.line, ypoints1, lty=1)
  lines(x.line, ypoints2, lty=2)

  sapply(1:ncol(dG.matrix), function(i_c) {
    x <- dG.matrix[, i_c]
    y <- Kd.matrix[, i_c]
    Points(x[order(x)], y[order(x)], pch=c(0, legend_pch),
           pt.lwd=c(1, 0), col=cols.mirnas[i_c])
  })

  kds_global <<- Kd.matrix
  sps_global <<- dG.matrix

  kMirnas[length(kMirnas)] <- "miR-7"
  leg.xy1 <- GetPlotFractionalCoords(fx=1, fy=0.125, log='y', inv='xy')
  leg.xy1alt <- GetPlotFractionalCoords(fx=0.8, fy=0.125, log='y', inv='xy')

  par(lend=1)
  legend_text <- expression(italic(K)[D]==italic(e)^{-Delta*italic(G)/italic(RT)})      
  lty=c("1", line_dash_length, "1")
  lty=c(1, 2, 1)
  print(lty)
  legend(x=leg.xy1[1], y=leg.xy1[2], legend=c(kSeedSites, legend_text), bty="n",
         pch=c(19, 0, 1), pt.lwd=c(0, 1, 0), pt.cex=pt_cex_final, lty=c(1, 2, 1), lwd=c(1, 1, 1)*par()$lwd,
         seg.len=2.7,
         col=c("black", "black", "gray"), xjust=1, xpd=NA)
  # legend(x=leg.xy1[1], y=leg.xy1[2], legend=c(kSeedSites, legend_text), bty="n",
  #        lty="23", lwd=c(0, 1, 0)*par()$lwd, seg.len=3,
  #        col=c("black", "black", "gray"), pt.cex=pt_cex_final, xjust=1, xpd=NA)
  leg.xy2 <- GetPlotFractionalCoords(fx=0.025, fy=1, log='y', inv='xy')
  Legend(leg.xy2, legend=kMirnas, col=cols.mirnas)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}
