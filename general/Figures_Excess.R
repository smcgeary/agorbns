################################################################################
# FIGURE 
################################################################################

PlotAgoPrepPurityOld <- function(experiment="AGO_purity", unique=FALSE,
                              no_marker=FALSE, no_adapter=FALSE, height=5,
                              width=5, format="built-in", pdf.plot=FALSE) {
  out_matrix <- matrix(0)
  # Gets the expression data within the experiment:
  out <- do.call("cbind", lapply(c("miR-1", "miR-155"), function(mirna) {
    out <- do.call("cbind", lapply(seq(5, 7), function(i_s) {
      out <- do.call("cbind", lapply(c("I", "P"), function(cond) {
        condition <- paste0("S100", i_s, "_", cond)
        counts <- GetMirnaCountData(mirna, condition, unique=unique,
                                    no_marker=no_marker, no_adapter=no_adapter)
        count_total <- rowSums(counts)
        names(count_total) <- rownames(counts)
        count_total
      }))
      colnames(out) <- paste0("S100", i_s, "_", c("I", "P"))
      out
    }))
    out
  }))
  exclude_rows <- c("Unmapped")
  if (!(no_marker)) exclude_rows <- c(exclude_rows, "18nt_marker", "30nt_marker")
  if (!(no_adapter)) exclude_rows <- c(exclude_rows, "5p_adapter", "3p_adapter")
  # Get the rows similar to how I used to do it, without the redundant inputs
  out <- out[, c(1, 2, 8, 3, 4, 10, 5, 6, 12)]
  out <- out[, 4:9]
  # Takes out the markers and unpammped categories from the table.
  out_new <<- out
  print(exclude_rows)
  unmapped  <- out["Unmapped", ]
  print(unmapped)
  out <- out[!(rownames(out) %in% exclude_rows),]
  # Isolate the exogenously loaded miRNA in the S100 extract
  exog_rows <- c("miR-1", "miR-155")
  spike_rows <- c("dme-miR-14-5p", "xtr-miR-427")
  spike_df <- out[spike_rows, ]
  exog_df <- out[exog_rows, ]/colSums(spike_df)
  endog_df <- out[setdiff(rownames(out), c(exog_rows, spike_rows)), ]/colSums(spike_df)
  row_order <- order(-rowSums(endog_df[, c(1, 4)]))
  endog_df <- endog_df[row_order,]
  endog_df <- rbind(endog_df[1:10, ], colSums(endog_df[11:nrow(endog_df), ]))
  rownames(endog_df)[nrow(endog_df)] <- "Remaining miRNAs"
  norm_df <- rbind(exog_df, endog_df)

  cols <- c("red", "purple",
            sprintf("gray%s", floor(seq(20, 80, length.out=nrow(norm_df) - 2))))
  SubfunctionCall(FigureSaveFile)
  xmin <- 0
  xmax <- 9
  ymin <- 0
  ymax <- max(norm_df)*1.2
  BlankPlot()

  percentage_rows <- t(t(exog_df)/(colSums(endog_df) + colSums(exog_df)))
  print(percentage_rows)

  cumul_df <- matrix(0, nrow=nrow(norm_df), ncol=ncol(norm_df))


  for (i_row in 1:nrow(norm_df)) {
    for (j_row in i_row:nrow(norm_df)) {
      cumul_df[j_row,] <- cumul_df[j_row, ] + norm_df[i_row, ]
    }
  }

  rownames(cumul_df) <- rownames(norm_df)
  colnames(cumul_df) <- colnames(norm_df)
  cumul_df_lower <- rbind(c(rep(0, ncol(cumul_df))), cumul_df[-nrow(cumul_df), ])
  rownames(cumul_df_lower) <- rownames(cumul_df)
  conds_r_matrix <- matrix(rep(seq(ncol(cumul_df)), nrow(cumul_df)),
                           byrow=TRUE,
                           nrow=nrow(cumul_df), ncol=ncol(cumul_df))
  conds_l_matrix <- conds_r_matrix - 1
  x_l <- c(conds_l_matrix)
  x_r <- c(conds_r_matrix)
  y_b <- c(cumul_df_lower)
  y_t <- c(cumul_df)
  segments(x0=3, y0=0,y1=250, lty="23")
  rect(x_l, y_b, x_r, y_t, col=cols, border=NA)
  AddLinearAxis(2, tick.space=50, label.space=50, label="Spike\u2212normalized counts")
  AddLinearAxis(1, alt_lab=c("\u2212", "miR-1", "miR-155", "\u2212", "miR-1", "miR-155"), alt_lab_pos=seq(6), label="",
                angled=TRUE, noline=TRUE, alt_lab_y_dist=-5)
  text(x=c(1.5, 4.5), y=c(250, 250), labels=c("+AGO2", "\u2212 AGO2"), xpd=NA)

  percent_nums <- c(percentage_rows[1, 2], percentage_rows[2, 3])
  percent_text <- format(round(c(percentage_rows[1, 2], percentage_rows[2, 3]), 2), nsmall=2)
  print(percent_text)
  text(x=c(1.5, 2.5), y=colSums(norm_df)[c(2, 3)] + 10,
       labels=sprintf("%1.2f%%", percent_nums), cex=0.75)
  xy <- GetPlotFractionalCoords(0.7, 0.9)

  legend_labels <- rownames(norm_df)
  for (i in seq(length(legend_labels))) {
    label_i <- legend_labels[i]
    if (grepl("/", label_i)) {
      labels <- unlist(strsplit(label_i, split="/"))
      label_suffix <- sapply(labels, function(label) {
        unlist(strsplit(label, split="miR-"))[2]
        })
      merged_miRNA <- paste(label_suffix, collapse="/")
      merged_miRNA <- paste0("miR-", merged_miRNA)
      legend_labels[i] <- merged_miRNA
    }
  }
  legend_labels[length(legend_labels)] <- "Other miRNAs"
  legend(xy[1], xy[2], legend=legend_labels, col=cols, pch=15, bty="n",
         y.intersp=0.9, xpd=NA, pt.cex=1.5)


  norm_unmapped <- unmapped/colSums(spike_df)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}




PlotRelKdsVsDDG <- function(experiment="equilibrium", n_constant=5,
                                  sitelist="paper", dG.table=3,
                                  square=FALSE, format="built-in", pdf.plot=FALSE) {
  dG.df <- read.table(file=paste0("canonical_sites_mfe_", dG.table, ".txt"))
  colnames(dG.df) <- gsub("(^.*)\\.(.*$)", colnames(dG.df), replace="\\1-\\2")   
  R <- 1.987e-3 # in kcal K-1 mol-1
  T <- 310.15 # in K
  total_x <- c()
  total_y <- c()
  SubfunctionCall(FigureSaveFile)
  ymin <- 0.005
  ymax <- 2
  xmin <- -4
  xmax <- 0.5
  if (square) {
    adjusted = TRUE
  } else {
    adjusted = FALSE
  }
  BlankPlot(log="y", adjusted=adjusted, inv='xy')
  AddLogAxis(2, label="Kd relative to 6mer")
  if (format == "cairo") {
    AddLinearAxis(1, tick.space=0.5,
                  label=expression(Delta~Delta~italic(G)*"; predicted"))
  } else {
    AddLinearAxis(1, tick.space=0.5,
                  label=expression(Delta*Delta*italic(G)*"; predicted"))

  }
  x.line <- seq(-8, 6, length=20)
  y.line <- exp(x.line/(R*T))
  lines(x.line, y.line, col="gray80", lwd=4)
  for (mirna in kMirnas) {
    color = kMirnaColors[mirna]
    pars.matrix <- SubfunctionCall(EquilPars)
    kds   <- pars.matrix[paste0(kCanonicalSites, "_Kd"), ]$Mean
    names(kds) <- kCanonicalSites
    dG    <- dG.df[kCanonicalSites, mirna]
    names(dG) <- kCanonicalSites
    y <- kds/kds["6mer"]
    x  <- dG - dG["6mer"]
    total_x <- c(total_x, x)
    total_y <- c(total_y, y)
    lmodel <- lm(log(y) ~ x)
    # print(mirna)
    # print(sprintf("r^2: %.4f", cor.test(log(y), x)$estimate^2))
    # print(sprintf("p_value: %.4f", cor.test(log(y), x)$p.value))
    m <- lmodel$coefficients[2]
    # print(m*R*T)
    b <- lmodel$coefficients[1]
    # print(b)
    y.line <- exp(m*x.line + b)
    lines(x.line, y.line, col=color, lty=2)
    points(x, y, col=c(rep(color, 3), "gray30"), cex=1.5, lwd=1.2,
           pch=c(1, 2, 0, 20)) 
  }
  text(3, -3.5, round(cor(total_x, total_y)^2, 3))
  xy <- GetPlotFractionalCoords(fx=0.4, fy=0.3, log='y', inv='xy')
  legend(x=xy[1], y=xy[2], seg.len=1, legend=kMirnas, col=kMirnaColors[kMirnas],
         bty="n", lty=1)
  xy <- GetPlotFractionalCoords(fx=0.7, fy=0.3, log='y', inv='xy')
  legend(x=xy[1], y=xy[2], legend = kCanonicalSites, bty="n",pch=c(1, 2, 0, 20),
         col=c(rep("black", 3), "gray30"))
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}


PlotPositionalEnrichment <- function(mirna, experiment="equilibrium",
                                     condition=40, n_constant=5,
                                     sitelist="paperextendedfinal",
                                     sites=kSeedSites, showdil=FALSE, xpos=20,
                                     ypos=20, width=5, height=5, combI=FALSE,
                                     buffer=FALSE, format="built-in",
                                     pdf.plot=FALSE) {
  sXc <- SubfunctionCall(SitesXCounts)
  # sXc <- SubfunctionCall(SingleSitesXCounts)
  if (combI) {
    I <- SubfunctionCall(GetPositionalSites, condition="I_combined", single=TRUE,
                         buffer=buffer)
  } else {
    I <- SubfunctionCall(GetPositionalSites, condition="I", single=TRUE,
                         buffer=buffer)
  }
  A <- SubfunctionCall(GetPositionalSites,
                       single=TRUE, buffer=buffer)
  if (condition == 12.6) {
    condition <- 12.65
  } else if (condition == 1.26) {
    condition <- 1.265
  }
  condition <- as.character(condition)
  R_all <- A/sum(sXc[, condition])/(I/sum(sXc[, 1]))
  R_overall_2 <- rowSums(A)/sum(sXc[, condition])/(rowSums(I)/sum(sXc[, 1]))
  R_overall <- Norm(sXc[, condition])/Norm(sXc[, 1])
  R_all <- data.matrix(R_all)
  # R_all[which(is.na(R_all))] <- 1
  names(R_overall) <- rownames(sXc)
  if (length(sites)==0) {
    sites <- rownames(sXc)
  }
  sites_use <- intersect(sites, rownames(sXc))
  print(sites_use)
  SubfunctionCall(FigureSaveFile)
  xmin <- 2
  if (identical(sites, kSeedSites)) {
    xmax <- 40
    offsets <- c(2, 2, 1, 1, 2, 0)
  } else {
    # Regex to remove the region after the period.
    ending_positions <- as.integer(gsub("^(.*)\\.(.*)$", sites_use, replace="\\2", perl=TRUE))
    offsets <- -ending_positions - min(-ending_positions)
  }
  xmax <- 39

  names(offsets) <- sites_use

  ymin <- 1e-1
  ymax <- 1e3
  BlankPlot(log='y')
  print("1194")
  for (site in sites_use) {
    segments(xmin, R_overall_2[site], xmax, R_overall_2[site],
             col=ConvertRColortoRGB(kSiteColors[site], alpha=0.5))

    frac_A <- A[site, ]/sum(sXc[, condition])
    frac_I <- I[site, ]/sum(sXc[, 1])
    # Portion of code that only takes the "N0–N36" portion of the enrichment profiles
    col_inds <- grep("N", names(frac_A))
    R <- R_all[site,]
    R_names <- names(R)
    R <- unlist(c(rep(NaN, offsets[site]), R[1:(length(R) - offsets[site])]))
    names(R) <- R_names
    R <- R[ceiling(xmin):floor(xmax)]
    lines(1:length(R), R, type="o", col= kSiteColors[site])        
  }
  AddLinearAxis(1, tick.space=1, label.space=1, label="Position",
                alt_lab=    c(1,  5, 10, 15, 20, 25, 30, 35, 37),
                alt_lab_pos=c(6, 10, 15, 20, 25, 30, 35, 40, 42))
  AddLogAxis(2, label="Enrichment")
  xy <- GetPlotFractionalCoords(0.05, 1, log='y')
  mirna.split <- paste0(strsplit(mirna, split="-")[[1]][1:2], collapse="-")
  text(xy[1], xy[2], labels=mirna.split, adj=c(0, 1))  
  if (showdil) {
    xy <- GetPlotFractionalCoords(0.05, 0.95, log='y')
    text(xy[1], xy[2], labels=condition, adj=c(0, 1))    
  }
  xy <- GetPlotFractionalCoords(0.97, 0, log='y')
  legend(x=xy[1], y=xy[2], bty="n", ncol=2, legend=sites_use, pch=19,
         col=kSiteColors[sites_use], xjust=1, yjust=0)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}



PlotHeatMap <- function(mirna, condition, experiment="equilibrium", n_constant=5,
                        sitelist="mismatch_and_threeprime", buffer=FALSE, xpos=20, ypos=20,
                        height=12, width=7, combI=FALSE, format="built-in",
                        pdf.plot=FALSE) {

  # sXc <- SubfunctionCall(SitesXCounts)
  # sXc_tot <- colSums(sXc)
  # names(sXc_tot) <- colnames(sXc)
  test_I <- SubfunctionCall(GetPositionalSites, condition="I") + 1
  test_A <- SubfunctionCall(GetPositionalSites) + 1
  print(test_I)
  test_I <- Norm(test_I + 1)
  test_A <- Norm(test_A)


  print(test_I)
  # test_I_ps <- Norm(test_I_ps)
  # test_A_ps <- Norm(test_A_ps)

  # test_I <- t(t(test_I)/rowMeans(test_I))
  # test_A <- t(t(test_A)/rowMeans(test_A))

  # print(test_I)
  # print(test_I_ps)

  # test_I_ps_global <<- test_I_ps
  # test_I_global <<- test_I
  sites <- rownames(test_I)
  ymax <- length(sites) + 0.5

  y_dist_124_norm <- 20/(nrow(SitesXCounts("miR-124")) + 0.5)
  if (sitelist == "paper") {
    height <- nrow(sXc) + 1
  }


R_40 <- (test_A/test_I)
R_40 <- data.matrix(R_40)

R_40_global <<- R_40
R_40 <- log10(R_40)
R_40_old <- R_40


row_means <- rowMeans(R_40, na.rm=TRUE)
R_40 <- R_40[order(row_means),]
R_40 <- R_40 - rowMeans(R_40, na.rm=TRUE)
row_sds <- apply(R_40, 1, sd, na.rm=TRUE)
R_40 <<- R_40
# R_40 <- t(t(R_40)/row_sds)
R_40[which(R_40 == Inf)] <- NaN
R_40[which(R_40 == -Inf)] <- NaN

cols_keep <- max(which(colSums(is.na(R_40)) != nrow(R_40)))
# R_40 <- R_40[, 1:cols_keep]

# R_40 <- R_40/sd(R_40, na.rm=TRUE)
R_min <- min(R_40, na.rm=TRUE)
R_max <- max(R_40, na.rm=TRUE)

R_min <- -1
R_max <- 1
# print(R_min)
# print(R_max)

dens <- density(c(R_40), na.rm=TRUE, n=100, from=-1, to=1)
out_global <<- R_40


R_40[which(is.na(R_40))] <- -10

SubfunctionCall(FigureSaveFile)


xmin <- -1/ncol(R_40)
xmax <- 1.4
ymin <- -1/nrow(R_40)
ymax <- 1 + 1/nrow(R_40)

par(kPlotParameters)
BlankPlot()

breaks <- c(-15, seq(R_min, R_max, length.out=21))
break_cols <- c("gray",rev(rainbow(20, start=0.9, end=0.8)))
image(t(data.matrix(R_40)), breaks=breaks, col=break_cols, xpd=NA, add=TRUE)
xy <- GetPlotFractionalCoords(0.05, 1)
text(xy[1], xy[2], labels=mirna, xpd=NA)
text(x=1.1, y=seq(0, nrow(R_40) - 1)/(nrow(R_40) - 1), adj=0, labels=rownames(R_40), xpd=NA)

x <- dens$x
min_lab_x <- min(x)
max_lab_x <- max(x)
ind_global <- 2
colors_points <- sapply(x, function(x_i) {
  # print(x_i)
  # print(breaks)
  inds <- which(breaks >= x_i)
  # print(inds)
  if (length(inds) == 0) {
    ind <- length(break_cols) - 1
  } else {
    ind <- min(inds)
  }
  # print(ind)
  break_cols[ind]
  })
# x <- rep(x, each=2)


x <- x - min(x)
x <- x/(max(x))
y <- dens$y/max(dens$y)
if (sitelist == "paper") {
  offset1 <- GetPlotFractionalCoords(0, -0.05)[2]
} else {
  offset1 <- GetPlotFractionalCoords(0, -0.025)[2]
}
text(0, offset1, labels=10^min_lab_x, xpd=NA)
text(1, offset1, labels=10^max_lab_x, xpd=NA)
# x <- x
print(colors_points)
x.mat <- matrix(x, nrow=2)
y.mat <- matrix(y, nrow=2)
if (sitelist == "paper") {
  offset2 <- GetPlotFractionalCoords(0, -0.1)[2]
} else {
  offset2 <- GetPlotFractionalCoords(0, -0.05)[2]  
}
y <- y*(offset1 - offset2)

trend_out <<- cbind(x, y)
# print(x[1:4])
# print(x.mat[1:2,])
# print(y.mat)
x.mat <- rbind(x.mat[1,],
               x.mat[1,],
               x.mat[2,], 
               x.mat[2,],
               x.mat[2,])
x.mat_alt <- c(rep(x[1], 2), rep(x[2:(length(x) - 1)], each=5), rep(x[length(x)], 3))
y.mat_alt <- c(y[1], rep(y[2:(length(y) - 1)], each=2), y[length(y)])
y.mat_alt <- matrix(y.mat_alt, nrow=2)
y.mat_alt <- rbind(rep(offset2, ncol(y.mat_alt)),
                   y.mat_alt + offset2,
                   rep(offset2, ncol(y.mat_alt)),
                   rep(NA, ncol(y.mat_alt)))
y.mat_alt_global <<- y.mat_alt
x.mat_alt_global <<- x.mat_alt
# print(x.mat_alt)
# print(x.mat)
#   x.mat.use <- x.mat*29 + 1
  # y.mat <- rbind(rep(offset, ncol(y.mat)),
  #                y.mat + offset,
  #                rep(offset, ncol(y.mat)),
  #                rep(NA, ncol(y.mat)))
  # print(y.mat)
  x_out <<- x.mat
  y_out <<- y.mat
   polygon(x.mat_alt, c(y.mat_alt), border=colors_points, col=colors_points, xpd=NA)
 
#   x.col.inds <- sapply(round(x.mat[1,]*m + b), function(i) {min(max(i, 1), 100)})
#   x.cols <- color.dist[x.col.inds]

# print(x.mat)
# print(y.mat)

# x_global <<- x
# y_global <<- y
# lines(x, y, xpd=NA)
# x_out <<- dens$x
# y_out <<- dens$y

  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

check_paperfinal <- GetPositionalSites("miR-1", "equilibrium", 40, sitelist="centered11")
check_paperfinal_I <- GetPositionalSites("miR-1", "equilibrium", "I", sitelist="centered11")
check_centered11 <- GetPositionalSites("miR-1", "equilibrium", 40, sitelist="centered11", buffer=TRUE)
check_centered11_I <- GetPositionalSites("miR-1", "equilibrium", "I", sitelist="centered11", buffer=TRUE)



graphics.off()
dev.new(xpos=20, ypos=20, height=5, width=5)

plot(1:ncol(check_paperfinal),
     (check_paperfinal["6mer-m8", ])/(check_paperfinal_I["6mer-m8",]), log='y',
     type="o", ylim=c(1e-1, 1e3))
lines(1:ncol(check_paperfinal), (check_centered11["6mer-m8", ])/(check_centered11_I["6mer-m8", ]), type="o", col="black")

lines(1:ncol(check_paperfinal), (check_centered11["11mer-m3.13", ])/(check_centered11_I["11mer-m3.13", ]), type="o", col="red")
lines(1:ncol(check_paperfinal), (check_centered11["11mer-m4.14", ])/(check_centered11_I["11mer-m4.14", ]), type="o", col="orange")
lines(1:ncol(check_paperfinal), (check_centered11["11mer-m5.15", ])/(check_centered11_I["11mer-m5.15", ]), type="o", col="forestgreen")
lines(1:ncol(check_paperfinal), (check_centered11["11mer-m6.16", ])/(check_centered11_I["11mer-m6.16", ]), type="o", col="blue")
lines(1:ncol(check_paperfinal), (check_centered11["11mer-m7.17", ])/(check_centered11_I["11mer-m7.17", ]), type="o", col="purple")


dev.new(xpos=520, ypos=20, height=5, width=5)

plot(1:ncol(check_paperfinal),
     (check_paperfinal["6mer-m8", ])/(check_paperfinal_I["6mer-m8",]), log='y',
     type="o", ylim=c(1e-1, 1e3))
lines(1:ncol(check_paperfinal), (check_centered11[ "11mer-m8.18", ])/(check_centered11_I[ "11mer-m8.18", ]), type="o", col="red")
lines(1:ncol(check_paperfinal), (check_centered11[ "11mer-m9.19", ])/(check_centered11_I[ "11mer-m9.19", ]), type="o", col="orange")
lines(1:ncol(check_paperfinal), (check_centered11["11mer-m10.20", ])/(check_centered11_I["11mer-m10.20", ]), type="o", col="forestgreen")
lines(1:ncol(check_paperfinal), (check_centered11["11mer-m11.21", ])/(check_centered11_I["11mer-m11.21", ]), type="o", col="blue")
lines(1:ncol(check_paperfinal), (check_centered11["11mer-m12.22", ])/(check_centered11_I["11mer-m12.22", ]), type="o", col="purple")


break

PlotHeatMap("miR-1", 40, sitelist="centered11", buffer=TRUE)

break



test_canonical_I <- GetPositionalSites("miR-1", experiment="equilibrium", condition="I",
                                                sitelist="canonical", buffer=TRUE)


test_canonical <- GetPositionalSites("miR-1", experiment="equilibrium", condition=40,
                                                sitelist="canonical", buffer=TRUE)

test_baek_I <- GetPositionalSites("miR-1", "equilibrium", condition="I",
                                                sitelist="baek", buffer=TRUE)

test_baekalt_I <- GetPositionalSites("miR-1", "equilibrium", condition="I",
                                                sitelist="baekalt", buffer=TRUE)


test_baek <- GetPositionalSites("miR-1", "equilibrium", condition=40,
                                                sitelist="baek", buffer=TRUE)


test_baekalt <- GetPositionalSites("miR-1", "equilibrium", condition=40,
                                                sitelist="baekalt", buffer=TRUE)


graphics.off()
# PlotHeatMap("miR-1", "40", sitelist="paperfinal", buffer=TRUE)
# PlotHeatMap("miR-1", "40", sitelist="paperextendedfinal", buffer=TRUE)
PlotHeatMap("miR-1", "40", sitelist="baek", buffer=TRUE)

PlotHeatMap("miR-1", "40", sitelist="baekalt", buffer=TRUE)
PlotHeatMap("let-7a", "40", sitelist="baek")
# PlotHeatMap("miR-1", "40", sitelist="centered11", buffer=TRUE)
# PlotHeatMap("miR-1", "40", sitelist="bulge", buffer=TRUE)
# PlotHeatMap("miR-1", "40", sitelist="del", buffer=TRUE)

break


PlotDist <- function() {
    x.mat.use <- x.mat*29 + 1
  y.mat <- rbind(rep(offset, ncol(y.mat)),
                 offset - y.mat,
                 rep(offset, ncol(y.mat)),
                 rep(NA, ncol(y.mat)))
  x.col.inds <- sapply(round(x.mat[1,]*m + b), function(i) {min(max(i, 1), 100)})
  x.cols <- color.dist[x.col.inds]
  polygon(c(y.mat), c(x.mat.use), border=x.cols, col=x.cols)

}
# PlotHeatMap("miR-1", "40", pdf.plot="S1.D_heatmap_miR-1")
# PlotHeatMap("let-7a", "40", pdf.plot="S1.D_heatmap_let-7a")
# PlotHeatMap("miR-155", "40", pdf.plot="S1.D_heatmap_miR-155")
# PlotHeatMap("miR-124", "40", pdf.plot="S1.D_heatmap_miR-124")
# PlotHeatMap("lsy-6", "40", pdf.plot="S1.D_heatmap_lsy-6")
# PlotHeatMap("miR-7-23nt", "40", exp="equilibrium2_nb", pdf.plot="S1.D-heatmap_miR-7")

PlotHeatMap("miR-1", "4", sitelist="paper", pdf.plot="S1.D_heatmap_miR-1_paper")
PlotHeatMap("let-7a", "4", sitelist="paper", pdf.plot="S1.D_heatmap_let-7a_paper")
PlotHeatMap("miR-155", "4", sitelist="paper", pdf.plot="S1.D_heatmap_miR-155_paper")
PlotHeatMap("miR-124", "4", sitelist="paper", pdf.plot="S1.D_heatmap_miR-124_paper")
PlotHeatMap("lsy-6", "4", sitelist="paper", pdf.plot="S1.D_heatmap_lsy-6_paper")
PlotHeatMap("miR-7-23nt", "12.6", sitelist="paper", exp="equilibrium2_nb", pdf.plot="S1.D-heatmap_miR-7_paper", height=7)


break


PlotSiteEnrichment <- function(mirna, site, experiment="equilibrium",
                               condition="40",
                               sitelist="mismatch_and_threeprime", xpos=20,
                               ypos=20, height=5, width=5, format="built-in",
                               pdf.plot=FALSE) {
  sXc <- SubfunctionCall(SitesXCounts)


  sXc_tot <- colSums(sXc)

  test_I <- GetPositionalSites(mirna, experiment, condition="I", sitelist=sitelist, single=TRUE)/(sXc_tot["I"])
  test_A <- GetPositionalSites(mirna, experiment, condition=condition, sitelist=sitelist, single=TRUE)/(sXc_tot[condition])
  sequence <- GetSiteSeq("miR-1", site)
  sequence_length <- nchar(sequence)

  # First plot
  pdf.plot_orig <- pdf.plot
  pdf.plot <- paste0(pdf.plot, "_5p")
  SubfunctionCall(FigureSaveFile)
  xmin <- 0
  xmax <- 47
  ymin <- 0.1
  ymax <- 1e3
  BlankPlot(log='y')

  # Make the left-hand background rectangles:

  rect(0.5, 0.1, 5.5, 1000, border=NA, col="steelblue1", xpd=NA)
  rect(5.5, 0.1, 37.5 + 5, 1000, col="pink", border=NA)
  rect(37.5 + 5, 0.1, 40.5 + 5, 1000, border=NA, col="greenyellow", xpd=NA)
  rect(40.5 + 5, 0.1, 42.5 + 5, 1000, border=NA, col="steelblue1", xpd=NA)

  AddLinearAxis(1, alt_lab_pos=c(1, 6, 10, 15, 20, 25, 30, 35, 42, 47), 
                alt_lab=c("-5", "1", "5", "10", "15", "20", "25", "30", "37", "C5"),tick.space=1, label.space=5, label="Position")
  AddLogAxis(2, label="Enrichment")
  if (site %in% kSiteColors) {
    color <- kSiteColors[site,]
  } else {
    color <- "black"
  }
  R_vals <- test_A[site, 1:ncol(test_I)]/test_I[site, 1:ncol(test_I)]
  R_vals_5p <- R_vals[1:5]

  R_vals_random <- R_vals[(1 + 5):(37 + 5)]
  R_vals_3p <- R_vals[(38 + 5):(40 + 5)]
  R_vals_3p2 <- R_vals[41:42 + 5]

  cex_text <- height/5
  print(R_vals_random)
  points(1:length(R_vals), R_vals, type="o", col=color)
  text(1:5, 1e3, c("C", "G", "A", "T", "C"), adj=c(0.5, 0), xpd=NA)
  text((5 + 1):(5 + 37), 1e3,
       rep("N", 37), adj=c(0.5, 0), xpd=NA, col="gray")
  text((1 + 5):(5 + 37), 1.5e3,
       unlist(strsplit("123456789!123456789@123456789#1234567", split="")),
       adj=c(0.5, 0), xpd=NA, col="gray")
  text(38:40 + 5, 1e3, c("T", "C", "G"), adj=c(0.5, 0), xpd=NA, col="forestgreen")
  text(41:42 + 5, 1e3, c("T", "A"), adj=c(0.5, 0), xpd=NA)


  xy <- GetPlotFractionalCoords(0.05, 0.1, log="y")
  text(xy[1], xy[2], site, adj=0)
  xy <- GetPlotFractionalCoords(0.05, 0.05, log="y")
  text(xy[1], xy[2], sequence, adj=0)
  new_frac <- 1
  R_vals_random_log <- data.matrix(log(R_vals_random))

  R_vals_random_log <<-R_vals_random_log
  R_vals_random_z <- (R_vals_random_log
                      - mean(R_vals_random_log[-length(R_vals_random_log)], na.rm=TRUE))/(2*sd(R_vals_random_log[-length(R_vals_random_log)], na.rm=TRUE))
  message("R_vals_random_z")
  print(R_vals_random_z)
  print(R_vals_random_z[1] > 1)
  if ((!is.na(R_vals_random_z[1]) & R_vals_random_z[1] > 2) | site %in% c("6mer", "7mer-m8")) {
    ind <- length(R_vals_random_z)
    print(ind)  
    xy <- GetPlotFractionalCoords(1, new_frac, log="y")[2]
    seqs_list <- unlist(strsplit(sequence, split=""))
    text(1:sequence_length + ind, seqs_list, adj=c(0.5, 1))
    new_frac <- new_frac - 0.05
    xy <- GetPlotFractionalCoords(1, new_frac, log="y")[2]
    seqs_list_8mer <- unlist(strsplit(GetSiteSeq("miR-1", "8mer"), split=""))
    seqs_list_7merA1 <- unlist(strsplit(GetSiteSeq("miR-1", "7mer-A1"), split=""))
    if (seqs_list[1] == seqs_list_8mer[1]) {
      seqs_list <- seqs_list_8mer
    } else if (seqs_list[1] == seqs_list_7merA1[1]) {
      seqs_list <- seqs_list_7merA1
    }
    # text(1:length(seqs_list) + ind + 5 + 37, xy, seqs_list, adj=c(0.5, 1))
    new_frac <- new_frac - 0.05


    } else if (sum(is.na(R_vals_5p)) != length(R_vals_5p)) {
    for (ind in which(!is.na(R_vals_5p))) {
      print("ind")
      print(ind)
      xy <- GetPlotFractionalCoords(1, new_frac, log="y")[2]
      seqs_list <- unlist(strsplit(sequence, split=""))
      text(1:sequence_length + ind - 1, xy, seqs_list, adj=c(0.5, 1))
      new_frac <- new_frac - 0.05
      xy <- GetPlotFractionalCoords(1, new_frac, log="y")[2]
      seqs_list_8mer <- unlist(strsplit(GetSiteSeq("miR-1", "8mer"), split=""))
      seqs_list_7merA1 <- unlist(strsplit(GetSiteSeq("miR-1", "7mer-A1"), split=""))
      if (seqs_list[1] == seqs_list_8mer[1]) {
        seqs_list <- seqs_list_8mer
      } else if (seqs_list[1] == seqs_list_7merA1[1]) {
        seqs_list <- seqs_list_7merA1
      }
      text(1:length(seqs_list) + ind - 1, xy, seqs_list, adj=c(0.5, 1))
      new_frac <- new_frac - 0.05

    }
  }

  if (class(pdf.plot) == "character") {
    dev.off()
  }


  pdf.plot <- paste0(pdf.plot_orig, "_3p")
  SubfunctionCall(FigureSaveFile)


  xmin <- 0
  xmax <- 47
  ymin <- 0.1
  ymax <- 1e3
  BlankPlot(log='y')
  min_value <- max(0, sequence_length - 5)
  rect(min_value + 0.5, 0.1, 37.5 - min_value + 1, 1000, col="pink", border=NA)
  rect(37.5 - sequence_length + 1 + 5, 0.1, 40.5 - sequence_length + 1 + 5, 1000, border=NA, col="greenyellow", xpd=NA)
  rect(40.5 - sequence_length + 1 + 5, 0.1, 47.5 - sequence_length + 1 + 5, 1000, border=NA, col="steelblue1", xpd=NA)

  AddLinearAxis(1, tick.space=1, label.space=5, label="Position")
  AddLogAxis(2, label="Enrichment")
  if (site %in% kSiteColors) {
    color <- kSiteColors[site,]
  } else {
    color <- "black"
  }
  R_vals <- test_A[site, 1:ncol(test_I)]/test_I[site, 1:ncol(test_I)]
  R_vals_5constant <- R_vals[1:5]

  R_vals_constant <- R_vals[(5 + 38 - sequence_length + 1):length(R_vals)]
  R_vals_random <- R_vals[6:(38 - sequence_length)]
  points(1:length(R_vals), R_vals, type="o", col=color)
  text(1:5, 3e3, c("C", "G", "A", "T", "C"), adj=c(0.5, 0), xpd=NA,
       cex=cex_text)
  text(1:37, 1e3, rep("N", 37), adj=c(0.5, 0), xpd=NA, col="gray")
  text(1:37, 1.5e3,
       unlist(strsplit("123456789!123456789@123456789#1234567", split="")),
       adj=c(0.5, 0), xpd=NA, col="gray")

  # text(sequence_length:37 - sequence_length + 1, 1e3, rep("N", 37), adj=c(0.5, 0), xpd=NA, col="gray")
  text(38:40 - sequence_length + 1 + 5, 1e3,
       c("T", "C", "G"), adj=c(0.5, 0), xpd=NA, col="forestgreen")
  text(41:47 - sequence_length + 1 + 5, 1e3,
       c("T", "A", "T", "G", "C", "C", "G"), adj=c(0.5, 0), xpd=NA)
  xy <- GetPlotFractionalCoords(0.05, 0.1, log="y")
  text(xy[1], xy[2], site, adj=0)
  xy <- GetPlotFractionalCoords(0.05, 0.05, log="y")
  text(xy[1], xy[2], sequence, adj=0)
  new_frac <- 1
  R_vals_random_log <- data.matrix(log(R_vals_random))
  R_vals_random_log <<-R_vals_random_log
  R_vals_random_z <- ((R_vals_random_log
                      - mean(R_vals_random_log[-length(R_vals_random_log)],
                             na.rm=TRUE))/
                      (2*sd(R_vals_random_log[-length(R_vals_random_log)],
                            na.rm=TRUE)))
  message("R_vals_random_z")
  print(R_vals_random_z)
  print(R_vals_random_z[length(R_vals_random_z)] > 1)
  if (R_vals_random_z[length(R_vals_random_z)] > 2 | site %in% c("6mer", "7mer-m8")) {
    ind <- length(R_vals_random_z)  
    xy <- GetPlotFractionalCoords(1, new_frac, log="y")[2]
    seqs_list <- unlist(strsplit(sequence, split=""))
    text(1:sequence_length + ind - sequence_length + 5, xy, seqs_list, adj=c(0.5, 1))
    new_frac <- new_frac - 0.05
    xy <- GetPlotFractionalCoords(1, new_frac, log="y")[2]
    seqs_list_8mer <- unlist(strsplit(GetSiteSeq("miR-1", "8mer"), split=""))
    seqs_list_7merA1 <- unlist(strsplit(GetSiteSeq("miR-1", "7mer-A1"), split=""))
    if (seqs_list[1] == seqs_list_8mer[1]) {
      seqs_list <- seqs_list_8mer
    } else if (seqs_list[1] == seqs_list_7merA1[1]) {
      seqs_list <- seqs_list_7merA1
    }
    text(1:length(seqs_list) + ind - sequence_length + 5, xy, seqs_list, adj=c(0.5, 1))
    new_frac <- new_frac - 0.05


    } else if (sum(is.na(R_vals_constant)) != length(R_vals_constant)) {
    for (ind in which(!is.na(R_vals_constant))) {
      print(ind)
      xy <- GetPlotFractionalCoords(1, new_frac, log="y")[2]
      seqs_list <- unlist(strsplit(sequence, split=""))
      text(1:sequence_length + 37 + ind - 1- sequence_length - sequence_length  + 5, xy, seqs_list, adj=c(0.5, 1))
      new_frac <- new_frac - 0.05
      xy <- GetPlotFractionalCoords(1, new_frac, log="y")[2]
      seqs_list_8mer <- unlist(strsplit(GetSiteSeq("miR-1", "8mer"), split=""))
      seqs_list_7merA1 <- unlist(strsplit(GetSiteSeq("miR-1", "7mer-A1"), split=""))
      if (seqs_list[1] == seqs_list_8mer[1]) {
        seqs_list <- seqs_list_8mer
      } else if (seqs_list[1] == seqs_list_7merA1[1]) {
        seqs_list <- seqs_list_7merA1
      }
      text(1:length(seqs_list) + 37 + ind - 1 - sequence_length - sequence_length + 5, xy, seqs_list, adj=c(0.5, 1))
      new_frac <- new_frac - 0.05

    }
  }
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}


PlotKineticData <- function(mirna, site, experiment="kinetics", n_constant=5, sitelist="paper",
                                   height=5, width=5, xpos=520, ypos=420,
                                   pdf.plot=FALSE) {
  args <- as.list(match.call())
  site.koffs <- SubfunctionCall(GetSiteKoffs, args, experiment="kinetics",
                           costfunc="multinom", subset=FALSE, dil=FALSE)
  names(site.koffs) <- gsub("^(.*)mer\\.(.*)$", names(site.koffs), replace="\\1-\\2", perl=TRUE)
  print(names(site.koffs))
  site.koff <- 10^site.koffs[paste0(site, "_koff")]
  print(1/site.koff)
  sXc <- SitesXCountsKinetics(mirna, experiment, n_constant, sitelist)

  print(sXc)
  sXc.p <- sXc[[1]]
  sXc.c <- sXc[[2]]
  sXc.p <- sXc.p[, 3:ncol(sXc.p)]
  sXc.c <- sXc.c[, 3:ncol(sXc.c)]
  totals <- colSums(sXc.p) + colSums(sXc.c)
  sXc.p <- t(t(sXc.p) / totals)
  sXc.c <- t(t(sXc.c) / totals)
  R.p <- sXc.p/sXc.p[, 1]
  R.c <- sXc.c/sXc.c[, 1]
  # sXc.c <- GetAverageOfReplicates(unique(colnamesapply(sXc.c, 2, Norm)
  if (class(pdf.plot) == "character") {
    pdf(file=paste0("2017_Paper/", pdf.plot, ".pdf"), height=2, width=2)
    par(kPDFParameters)
  } else {
    dev.new(xpos=xpos, ypos=ypos, height=height, width=width)
    par(kPlotParameters)    
  }
  temp <<- rownames(sitesXcounts)
  site.colors <- kSiteColors[rownames(sitesXcounts)]
  # Make dummy plot:
  x <- as.numeric(colnames(sXc.c)[2:(ncol(sXc.p)-2)])
  x[1] <- 1
  y.p <- R.p[site, 2:(ncol(sXc.p)-2), drop=FALSE]
  y.c <- R.c[site, 2:(ncol(sXc.p)-2), drop=FALSE]
  print(x)
  print(y.p)

  xmin <- 1
  xmax <- 100000
  ymin <- 0.03
  ymax <- 1000
  BlankPlot(log="xy")
  # Make x=y line:
  # Make the lines connecting the points to the x = y line:
  # Make axes:
  AddLogAxis(1, label="Time (s)")
  AddLogAxis(2, label="Enrichment")


  color = kSiteColors[site]
  # Add the points to the plot:
  points(x, y.p, col=color)
  points(x, y.c, col=color, pch=1)
  abline(v=1/site.koff, col=kSiteColors[site], lty=2)
  times_unique <- unique(round(x, digits=0))
  print(times_unique)

  lines(times_unique, sapply(times_unique, GetAverageOfReplicates, data=y.p, times=round(x, digits=0)), col=color)
  lines(times_unique, sapply(times_unique, GetAverageOfReplicates, data=y.c, times=round(x, digits=0)), lty=2, col=color)
  xy <- GetPlotFractionalCoords(fx=0.9, fy=0.85, log='xy')
  text(xy[1], xy[2], mirna, adj=1)
  xy <- GetPlotFractionalCoords(fx=0.9, fy=0.8, log='xy')
  text(xy[1], xy[2], site, adj=1, col=kSiteColors[site])

  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

PlotKoffsVsKds <- function(mirna, n_constant=5, sitelist="paper", ident=FALSE, 
                           costfunc="multinom", subset=FALSE, dil=FALSE,
                           alpha=0.5, pdf.plot=FALSE) {
    width <- 5
    height <- 5
    kds <- SubfunctionCall(EquilPars, experiment="equilibrium")
    koffs <- SubfunctionCall(GetSiteKoffs, experiment="kinetics")
    kinetics.koffs <- 10^koffs[grep("koff$", names(koffs), perl=TRUE)]

    xmin <- 5e-4
    xmax <- 10
    ymin <- 5e-4
    ymax <- 10
    if (class(pdf.plot) == "character") {
      pdf(file=paste0("2017_Paper/", pdf.plot, ".pdf"), height=height*2/5,
          width=width*2/5)
      par(kPDFParameters)
    } else {
      dev.new(xpos=xpos, ypos=ypos, height=height, width=width)
      par(kPlotParameters)    
    }
    BlankPlot(log='xy', inv='xy')  
    AddLogAxis(1, label="Kd")
    AddLogAxis(2, label="koff")  
    kds.plot <- kds[1:length(kinetics.koffs),1]
    print(kinetics.koffs)
    colors <- kSiteColors[rownames(kds)][1:length(kinetics.koffs)]
    points(kds.plot, kinetics.koffs, col=colors)

    xy <- GetPlotFractionalCoords(fx=0.9, fy=0.86, log='xy', inv='xy')
    text(xy[1], xy[2], mirna, adj=1)


  if (class(pdf.plot) == "character") {
    dev.off()
  }

}



PlotAllKdAndNumSiteDistribution <- function(experiment="equilibrium",
                                    n_constant=5, sitelist="paper",
                                    sites_use="all", n_cutoff=1,
                                    kd_cutoff=Inf, single_site=FALSE,
                                    best_site=FALSE, combined=TRUE,
                                    global=FALSE, fixed=FALSE, bulk=FALSE,
                                    cat_colors=FALSE, noncanon=FALSE,
                                    merge=FALSE, best=TRUE, old=FALSE,
                                    threePseq=TRUE, height=5, width=5, xpos=20,
                                    ypos=20, format="built-in", pdf.plot=FALSE) {
  xmin <- 1e-2
  xmax <- 1e4
  ymin <- 0
  ymax <- 0.1
  log_all <- 'x'
  SubfunctionCall(FigureSaveFile)
  BlankPlot(log=log_all)  
  AddLogAxis(1, label="Binding Mass (Relative [Input] / Relative Kd)", adj=TRUE)
  AddLinearAxis(2, tick.space=0.05, label.space=0.1, label="(repress/kd res.)^2
                + SE repression model coefficient") 
  n_sites_all <- c()
  y_all <- c()
  kds_all <- c()
  fc_all <- c()
  fc_sem_all <- c()
  cols_all <- c()
  fc_fitted_all <- c()
  mirnas_all <- c()
  sites_all <- c()
  total_sites <- 0

  len_kd <- rep(0, length(kMirnas))
  names(len_kd) <- kMirnas
  len_fc  <- len_kd
  len_sites <- len_kd



  sapply(kMirnas, function(mirna) {
    if (mirna == "miR-7-23nt") {
      experiment <- "equilibrium2_nb"
      fixed <- TRUE
    }
    if (mirna %in% c("miR-1", "miR-7-23nt")) {
      combined <- FALSE
    }
    kds <- SubfunctionCall(EquilPars)
    rownames(kds) <- gsub("_Kd", "", rownames(kds))
    sites <- grep(sprintf("(?:%s|None)", mirna), rownames(kds), perl=TRUE,
                  invert=TRUE, value=TRUE)
    adjusted <- FALSE
    l.xy <- GetPlotFractionalCoords(fx=0.05, fy=0.05, log=log_all)
    yjust <- 0
    l.cex <- 1
    fc <- SubfunctionCall(GetRepressionLinearModel)
    # print(mirna)
    # print(dim(fc))
    # print(nrow(kds))

    len_fc[mirna] <<- nrow(fc)
    len_kd[mirna] <<- nrow(kds[sites, ])
    total_sites <<- total_sites + nrow(fc)
    fc_sites <- fc[, 1]
    fc_sem <- fc[, 2]
    len_sites[mirna] <<- nrow(fc)

    segments(xmin, 0, xmax, 0, lwd=0.5, col="gray")
    kds <- kds[rownames(fc), ]
    sites <- rownames(kds)
    n_sites <- fc[, 3]
    # print(n_sites)

    colors.sites <- kSiteColors[sites]
    kds_sites <- kds$Mean
    names(kds_sites) <- rownames(fc)
    # print(kds_sites)
    kds_sites_lci <- kds$Lower_CI
    kds_sites_uci <- kds$Upper_CI
    lkds <- log(kds_sites)
    lkds_sem <- (log(kds_sites_uci) - log(kds_sites_lci))/2
    lm_fc_kd <- lm(fc_sites ~ lkds)
    lm_fc_kd <<- lm_fc_kd
    lm_fc_kd_deming <- deming(fc_sites ~ lkds, xstd=lkds_sem, ystd=fc_sem)
    m <- lm_fc_kd$coefficients[2]
    b <- lm_fc_kd$coefficients[1]
    m_d <- lm_fc_kd_deming$coefficients[2]
    b_d <- lm_fc_kd_deming$coefficients[1]
    x.line <- exp(seq(-10, 10, by=0.5))
    y.line <- m*log(x.line) + b
    y.line2 <- m_d*log(x.line) + b_d
    y <- sqrt(abs(fc_sites-lm_fc_kd$fitted.values)^2 + fc_sem^2)
    mirnas_all <<- c(mirnas_all, rep(mirna, length(n_sites)))
    sites_all <<- c(sites_all, names(fc_sites))
    cols_all <<- c(cols_all, rep(kMirnaColors[mirna], length(n_sites)))
    fc_all <<- c(fc_all, fc_sites)
    n_sites_all <<- c(n_sites_all, n_sites)
    y_all <<- c(y_all, y)
    kds_all <<- c(kds_all, lkds)
    fc_sem_all <<- c(fc_sem_all, fc_sem)
    fc_fitted_all <<- c(fc_fitted_all, lm_fc_kd$fitted.values)
})

  everything.df <<- data.frame(mirna=mirnas_all, site=sites_all, kd=kds_all, fc=fc_all, fc_sem=fc_sem_all, n_sites=n_sites_all)
  write.table(everything.df, file="everything_df.txt", quote=FALSE, sep="\t")
  lm_fc_kd_all <- lm(fc_all ~ kds_all)
  lm_fc_kd_fitted_all_new <- lm_fc_kd_all$fitted.values
  lm_fc_kd_fitted_all_new <<- lm_fc_kd_fitted_all_new
  y_new <- (fc_all-lm_fc_kd_fitted_all_new)^2 + fc_sem_all^2
  cols_all <- cols_all[order(n_sites_all)]
  y_new <- y_new[order(n_sites_all)]
  y_all <- y_all[order(n_sites_all)]
  n_sites_all <- sort(n_sites_all)
  points(n_sites_all,  y_new, col=cols_all)
  sapply(10^seq(0, 3, length.out=5), function(f) {
    smoothened <- lowess(log(n_sites_all), y_new, f=f, delta=0.00001)  
    lines(exp(smoothened$x), smoothened$y, col=ConvertRColortoRGB("black", alpha=0.5))
  })

  if (class(pdf.plot) == "character") {
    dev.off()
  }
  print("finished plot")
}



MakeSupplementalFigure1_AllPositions <- function() {
  n_constant <- 5
  sitelist <- "mismatch_and_threeprime"
  print("S1_2")
  print("miR-1")
  PlotPositionalEnrichment("miR-1", "equilibrium", 0.4, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.A1")
  PlotPositionalEnrichment("miR-1", "equilibrium", 1.26, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.A2")
  PlotPositionalEnrichment("miR-1", "equilibrium", 4, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.A3")
  PlotPositionalEnrichment("miR-1", "equilibrium", 12.6, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.A4")
  PlotPositionalEnrichment("miR-1", "equilibrium", 40, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.A5")
  PlotPositionalEnrichment("miR-1", "equilibrium", 0, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.A6")
  print("let-7a")
  PlotPositionalEnrichment("let-7a", "equilibrium", 0.4, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.B1")
  PlotPositionalEnrichment("let-7a", "equilibrium", 1.26, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.B2")
  PlotPositionalEnrichment("let-7a", "equilibrium", 4, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.B3")
  PlotPositionalEnrichment("let-7a", "equilibrium", 12.6, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.B4")
  PlotPositionalEnrichment("let-7a", "equilibrium", 40, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.B5")
  PlotPositionalEnrichment("let-7a", "equilibrium", 0, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.B6")
  print("S1_3")
  print("miR-155")
  PlotPositionalEnrichment("miR-155", "equilibrium", 0.4, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.C1")
  PlotPositionalEnrichment("miR-155", "equilibrium", 1.26, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.C2")
  PlotPositionalEnrichment("miR-155", "equilibrium", 4, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.C3")
  PlotPositionalEnrichment("miR-155", "equilibrium", 12.6, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.C4")
  PlotPositionalEnrichment("miR-155", "equilibrium", 40, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.C5")
  PlotPositionalEnrichment("miR-155", "equilibrium", 0, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.C6")

  PlotPositionalEnrichment("miR-155", "equilibrium", 0.4, sitelist=sitelist, showdil=TRUE, sites=k3PSites, n_constant=n_constant, pdf.plot="S1-2.D1")
  PlotPositionalEnrichment("miR-155", "equilibrium", 1.26, sitelist=sitelist, showdil=TRUE, sites=k3PSites, n_constant=n_constant, pdf.plot="S1-2.D2")
  PlotPositionalEnrichment("miR-155", "equilibrium", 4, sitelist=sitelist, showdil=TRUE, sites=k3PSites, n_constant=n_constant, pdf.plot="S1-2.D3")
  PlotPositionalEnrichment("miR-155", "equilibrium", 12.6, sitelist=sitelist, showdil=TRUE, sites=k3PSites, n_constant=n_constant, pdf.plot="S1-2.D4")
  PlotPositionalEnrichment("miR-155", "equilibrium", 40, sitelist=sitelist, showdil=TRUE, sites=k3PSites, n_constant=n_constant, pdf.plot="S1-2.D5")
  PlotPositionalEnrichment("miR-155", "equilibrium", 0, sitelist=sitelist, showdil=TRUE, sites=k3PSites, n_constant=n_constant, pdf.plot="S1-2.D6")
  print("S1_4")
  print("miR-124")
  PlotPositionalEnrichment("miR-124", "equilibrium", 0.4, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.E1")
  PlotPositionalEnrichment("miR-124", "equilibrium", 1.26, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.E2")
  PlotPositionalEnrichment("miR-124", "equilibrium", 4, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.E3")
  PlotPositionalEnrichment("miR-124", "equilibrium", 12.6, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.E4")
  PlotPositionalEnrichment("miR-124", "equilibrium", 40, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.E5")
  PlotPositionalEnrichment("miR-124", "equilibrium", 0, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.E6")

  PlotPositionalEnrichment("miR-124", "equilibrium", 0.4, sitelist=sitelist, showdil=TRUE, sites=k3PSites, n_constant=n_constant, pdf.plot="S1-2.F1")
  PlotPositionalEnrichment("miR-124", "equilibrium", 1.26, sitelist=sitelist, showdil=TRUE, sites=k3PSites, n_constant=n_constant, pdf.plot="S1-2.F2")
  PlotPositionalEnrichment("miR-124", "equilibrium", 4, sitelist=sitelist, showdil=TRUE, sites=k3PSites, n_constant=n_constant, pdf.plot="S1-2.F3")
  PlotPositionalEnrichment("miR-124", "equilibrium", 12.6, sitelist=sitelist, showdil=TRUE, sites=k3PSites, n_constant=n_constant, pdf.plot="S1-2.F4")
  PlotPositionalEnrichment("miR-124", "equilibrium", 40, sitelist=sitelist, showdil=TRUE, sites=k3PSites, n_constant=n_constant, pdf.plot="S1-2.F5")
  PlotPositionalEnrichment("miR-124", "equilibrium", 0, sitelist=sitelist, showdil=TRUE, sites=k3PSites, n_constant=n_constant, pdf.plot="S1-2.F6")
  print("S1_5")
  print("lsy-6")
  PlotPositionalEnrichment("lsy-6", "equilibrium", 0.4, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.G1")
  PlotPositionalEnrichment("lsy-6", "equilibrium", 1.26, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.G2")
  PlotPositionalEnrichment("lsy-6", "equilibrium", 4, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.G3")
  PlotPositionalEnrichment("lsy-6", "equilibrium", 12.6, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.G4")
  PlotPositionalEnrichment("lsy-6", "equilibrium", 40, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.G5")
  PlotPositionalEnrichment("lsy-6", "equilibrium", 0, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.G6")

  PlotPositionalEnrichment("lsy-6", "equilibrium", 0.4, sitelist=sitelist, showdil=TRUE, sites=k3PSites, n_constant=n_constant, pdf.plot="S1-2.H1")
  PlotPositionalEnrichment("lsy-6", "equilibrium", 1.26, sitelist=sitelist, showdil=TRUE, sites=k3PSites, n_constant=n_constant, pdf.plot="S1-2.H2")
  PlotPositionalEnrichment("lsy-6", "equilibrium", 4, sitelist=sitelist, showdil=TRUE, sites=k3PSites, n_constant=n_constant, pdf.plot="S1-2.H3")
  PlotPositionalEnrichment("lsy-6", "equilibrium", 12.6, sitelist=sitelist, showdil=TRUE, sites=k3PSites, n_constant=n_constant, pdf.plot="S1-2.H4")
  PlotPositionalEnrichment("lsy-6", "equilibrium", 40, sitelist=sitelist, showdil=TRUE, sites=k3PSites, n_constant=n_constant, pdf.plot="S1-2.H5")
  PlotPositionalEnrichment("lsy-6", "equilibrium", 0, sitelist=sitelist, showdil=TRUE, sites=k3PSites, n_constant=n_constant, pdf.plot="S1-2.H6")
  print("S1_6")
  print("miR-7-23nt")
  PlotPositionalEnrichment("miR-7-23nt", "equilibrium_nb", 0.4, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.I1")
  PlotPositionalEnrichment("miR-7-23nt", "equilibrium_nb", 1.26, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.I2")
  PlotPositionalEnrichment("miR-7-23nt", "equilibrium_nb", 4, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.I3")
  PlotPositionalEnrichment("miR-7-23nt", "equilibrium_nb", 12.6, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.I4")
  PlotPositionalEnrichment("miR-7-23nt", "equilibrium_nb", 40, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.I5")
  PlotPositionalEnrichment("miR-7-23nt", "equilibrium_nb", 0, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.I6")
  PlotPositionalEnrichment("miR-7-23nt", "equilibrium2_nb", 0.4, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.J1")
  PlotPositionalEnrichment("miR-7-23nt", "equilibrium2_nb", 1.26, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.J2")
  PlotPositionalEnrichment("miR-7-23nt", "equilibrium2_nb", 12.6, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.J4")
  PlotPositionalEnrichment("miR-7-23nt", "equilibrium2_nb", 40, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.J5")
  PlotPositionalEnrichment("miR-7-23nt", "equilibrium2_nb", 0, sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.J6")
  print("S1_7")
  print("equil_pilot")
  PlotPositionalEnrichment("miR-1", "equil_pilot", "L100A10", sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.K1")
  PlotPositionalEnrichment("miR-1", "equil_pilot", "L10A10", sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.K2")
  PlotPositionalEnrichment("miR-1", "equil_pilot", "L100A0", sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.K3")
  PlotPositionalEnrichment("miR-1", "equil_pilot", "L10A0", sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.K4")
  PlotPositionalEnrichment("miR-155", "equil_pilot", "L100A10", sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.K5")
  PlotPositionalEnrichment("miR-155", "equil_pilot", "L10A10", sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.K6")
  PlotPositionalEnrichment("miR-155", "equil_pilot", "L100A0", sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.L1")
  PlotPositionalEnrichment("miR-155", "equil_pilot", "L10A0", sitelist=sitelist, showdil=TRUE, n_constant=n_constant, pdf.plot="S1-2.L2")
  PlotPositionalEnrichment("miR-155", "equil_pilot", "L100A10", sitelist=sitelist, showdil=TRUE, sites=k3PSites, n_constant=n_constant, pdf.plot="S1-2.L3")
  PlotPositionalEnrichment("miR-155", "equil_pilot", "L10A10", sitelist=sitelist, showdil=TRUE, sites=k3PSites, n_constant=n_constant, pdf.plot="S1-2.L4")
  PlotPositionalEnrichment("miR-155", "equil_pilot", "L100A0", sitelist=sitelist, showdil=TRUE, sites=k3PSites, n_constant=n_constant, pdf.plot="S1-2.L5")
  PlotPositionalEnrichment("miR-155", "equil_pilot", "L10A0", sitelist=sitelist, showdil=TRUE, sites=k3PSites, n_constant=n_constant, pdf.plot="S1-2.L6")
}

PlotAllMir7Figures <- function() {
  PlotSiteKds("miR-7-22nt", exp="equilibrium_nb", adjusted_height=TRUE, combined=FALSE, added.text=TRUE,
              pdf.plot="0.miR-7-22nt,equilibrium_nb")
  PlotSiteKds("miR-7-23nt", exp="equilibrium_nb", adjusted_height=TRUE, combined=FALSE, added.text=TRUE,
              pdf.plot="0.miR-7-23nt,equilibrium_nb")
  PlotSiteKds("miR-7-24nt", exp="equilibrium_nb", adjusted_height=TRUE, combined=FALSE, added.text=TRUE,
              pdf.plot="0.miR-7-24nt,equilibrium_nb")
  PlotSiteKds("miR-7-25nt", exp="equilibrium_nb", adjusted_height=TRUE, combined=FALSE, added.text=TRUE,
              pdf.plot="0.miR-7-25nt,equilibrium_nb")
  PlotSiteKds("miR-7-23nt", exp="equilibrium2_nb", adjusted_height=TRUE, combined=FALSE, added.text=TRUE,
              pdf.plot="0.miR-7-23nt,equilibrium2_nb")
  PlotSiteKds("miR-7-24nt", exp="equilibrium2_nb", adjusted_height=TRUE, combined=FALSE, added.text=TRUE,
              pdf.plot="0.miR-7-24nt,equilibrium2_nb")
  PlotSiteKds("miR-7-25nt", exp="equilibrium2_nb", adjusted_height=TRUE, combined=FALSE, added.text=TRUE,
              pdf.plot="0.miR-7-25nt,equilibrium2_nb")

}

PlotAllMir7Figures2 <- function() {
  PlotSiteKds("miR-7-22nt", exp="equilibrium_nb", sitelist="canonical", height=3.5, combined=FALSE, added.text=TRUE,
              pdf.plot="0.miR-7-22nt,equilibrium_nb_canonical")
  PlotSiteKds("miR-7-23nt", exp="equilibrium_nb", sitelist="canonical", height=3.5, combined=FALSE, added.text=TRUE,
              pdf.plot="0.miR-7-23nt,equilibrium_nb_canonical")
  PlotSiteKds("miR-7-24nt", exp="equilibrium_nb", sitelist="canonical", height=3.5, combined=FALSE, added.text=TRUE,
              pdf.plot="0.miR-7-24nt,equilibrium_nb_canonical")
  PlotSiteKds("miR-7-25nt", exp="equilibrium_nb", sitelist="canonical", height=3.5, combined=FALSE, added.text=TRUE,
              pdf.plot="0.miR-7-25nt,equilibrium_nb_canonical")
  PlotSiteKds("miR-7-23nt", exp="equilibrium2_nb", sitelist="canonical", height=3.5, combined=FALSE, added.text=TRUE,
              pdf.plot="0.miR-7-23nt,equilibrium2_nb_canonical")
  PlotSiteKds("miR-7-24nt", exp="equilibrium2_nb", sitelist="canonical", height=3.5, combined=FALSE, added.text=TRUE,
              pdf.plot="0.miR-7-24nt,equilibrium2_nb_canonical")
  PlotSiteKds("miR-7-25nt", exp="equilibrium2_nb", sitelist="canonical", height=3.5, combined=FALSE, added.text=TRUE,
              pdf.plot="0.miR-7-25nt,equilibrium2_nb_canonical")

}

MakeOldFigure2 <- function() {
  PlotSiteKdsVsSPS(pdf.plot="2.K")
  Plot8merVs7merCombined(pdf.plot="2.L")
  PlotRelKdsVsDDG(pdf.plot="2.M")
  print("Done Figure 2")
}

MakeOldFigure3 <- function() {
  PlotSiteFlankEnrichments("miR-1", "8mer", pdf.plot="3.A")
  print("Done Figure 3A")
  PlotSiteFlankKds("miR-1", pdf.plot="3.B")
  print("Done Figure 3B")
  PlotStructureVsFlankingKds("miR-1", "8mer", pdf.plot="3.C")
  print("Done Figure 3C")
  PlotAllSamplePlFlanks("miR-1", "8mer", "0.4", pdf.plot="3.D")
  print("Done Figure 3D")
  print("Done Figure 3")
}

MakeOldSupplementalFigure3 <- function() {
  PlotSiteFlankKds("let-7a", pdf.plot="S3.A")
  print("Done figure S3A")
  PlotSiteFlankKds("miR-155", pdf.plot="S3.B")
  print("Done figure S3B")
  PlotSiteFlankKds("miR-124", pdf.plot="S3.C")
  print("Done figure S3C")
  PlotSiteFlankKds("lsy-6", pdf.plot="S3.D")
  print("Done figure S3D")
  PlotFlankLinModel(pdf.plot="S3.E")
  print("Done figure S3E")
  PlotFlankLinModelCoefficients(pdf.plot="S3.F")
  print("Done figure S3F")
  PlotFlankPlCorrelationWindowsNew("miR-1", "8mer", condition="I_combined",
                                   pdf.plot="S3.G")
  print("Done figure S3G")
  PlotFlanksSamplePl("miR-1", "8mer", 0.4, pdf.plot="S3.H")
  print("Done figure S3H")
  PlotFlanksSamplePl("miR-1", "8mer", 0.4, matchdist=TRUE, pdf.plot="S3.I")
  print("Done figure S3I")
}

MakeOldFigure4 <- function() {
  PlotSiteKdsVsRepression("miR-1", pdf.plot="4.A")
  print("Done figure 4A")
  PlotSiteKdsVsRepression("let-7a", pdf.plot="4.B")
  print("Done figure 4B")
  PlotSiteKdsVsRepression("miR-155", pdf.plot="4.C")
  print("Done figure 4C")
  PlotSiteKdsVsRepression("miR-124", pdf.plot="4.D")
  print("Done figure 4D")
  PlotSiteKdsVsRepression("lsy-6", pdf.plot="4.E")
  print("Done figure 4E")
  PlotAllSeedKdsVsRepression(pdf.plot="4.F")
  print("Done figure 4F") 
}

MakeOldSupplementalFigure4 <- function() {
  PlotSiteKdsVsRepressionNonCanon("miR-1", pdf.plot="S4.A")
  print("Done figure S4A")
  PlotSiteKdsVsRepressionNonCanon("let-7a", pdf.plot="S4.B")
  print("Done figure S4B")
  PlotSiteKdsVsRepressionNonCanon("miR-155", pdf.plot="S4.C")
  print("Done figure S4C")
  PlotSiteKdsVsRepressionNonCanon("miR-124", pdf.plot="S4.D")
  print("Done figure S4D")
  PlotSiteKdsVsRepressionNonCanon("lsy-6", pdf.plot="S4.E")
  print("Done figure S4E")
}

MakeFigure5 <- function() {
  PlotKineticData("miR-1", "8mer", pdf.plot="5.A1")
  PlotKineticData("miR-1", "7mer-m8", pdf.plot="5.A2")
  PlotKineticData("miR-1", "7mer-A1", pdf.plot="5.A3")
  PlotKineticData("miR-1", "6mer", pdf.plot="5.A4")

  PlotKineticData("let-7a", "8mer", pdf.plot="5.B1")
  PlotKineticData("let-7a", "7mer-m8", pdf.plot="5.B2")
  PlotKineticData("let-7a", "7mer-A1", pdf.plot="5.B3")
  PlotKineticData("let-7a", "6mer", pdf.plot="5.B4")

  PlotKineticData("miR-124", "8mer", pdf.plot="5.C1")
  PlotKineticData("miR-124", "7mer-m8", pdf.plot="5.C2")
  PlotKineticData("miR-124", "7mer-A1", pdf.plot="5.C3")
  PlotKineticData("miR-124", "6mer", pdf.plot="5.C4")

  PlotKineticData("lsy-6", "8mer", pdf.plot="5.D1")
  PlotKineticData("lsy-6", "7mer-m8", pdf.plot="5.D2")
  PlotKineticData("lsy-6", "7mer-A1", pdf.plot="5.D3")
  PlotKineticData("lsy-6", "6mer", pdf.plot="5.D4")


  PlotKineticData("miR-124", "11mer-m9.19", pdf.plot="5.E1")
  PlotKineticData("miR-124", "10mer-m9.18", pdf.plot="5.E2")
  PlotKineticData("miR-124", "10mer-m10.19", pdf.plot="5.E3")
  PlotKineticData("miR-124", "9mer-m9.17", pdf.plot="5.E4")
  PlotKineticData("miR-124", "9mer-m10.18", pdf.plot="5.E5")
  PlotKineticData("miR-124", "9mer-m11.19", pdf.plot="5.E6")

  PlotKineticData("lsy-6", "11mer-m9.19", pdf.plot="5.F1")
  PlotKineticData("lsy-6", "10mer-m9.18", pdf.plot="5.F2")
  PlotKineticData("lsy-6", "10mer-m10.19", pdf.plot="5.F3")
  PlotKineticData("lsy-6", "9mer-m9.17", pdf.plot="5.F4")
  PlotKineticData("lsy-6", "9mer-m10.18", pdf.plot="5.F5")
  PlotKineticData("lsy-6", "9mer-m11.19", pdf.plot="5.F6")
  PlotKoffsVsKds("miR-1", pdf.plot="5.G1")
  PlotKoffsVsKds("let-7a", pdf.plot="5.G2")
  PlotKoffsVsKds("miR-124", pdf.plot="5.G3")
  PlotKoffsVsKds("lsy-6", pdf.plot="5.G4")
}

MakeFigure1_Original <- function() {
  message("Making Figure 1")
  PlotEquilSiteWithInput("miR-1", 7, sitelist="canonical", combined=FALSE,
                         buffer=TRUE, pdf.plot="1.C")
  PlotSiteEnrichments("miR-1", sitelist="canonical", combined=FALSE,
                      buffer=TRUE, pdf.plot="1.D")
  PlotSiteEnrichments("miR-1", sitelist="paper", combined=FALSE, buffer=TRUE, pdf.plot="1.E_original")
  PlotSiteKds("miR-1", sitelist="paper", combined=FALSE, buffer=TRUE, pdf.plot="1.F_original")
  PlotSiteOccupancy("miR-1", sitelist="paper", combined=FALSE, buffer=TRUE, pdf.plot="1.G_original")
  message("Done Figure 1")
}

MakeFigure1_Cutoff <- function() {
  message("Making Figure 1")
  PlotEquilSiteWithInput("miR-1", 7, sitelist="canonical", combined=FALSE,
                         buffer=TRUE, pdf.plot="1.C")
  PlotSiteEnrichments("miR-1", sitelist="canonical", combined=FALSE,
                      buffer=TRUE, pdf.plot="1.D")
  PlotSiteEnrichments("miR-1", sitelist="papercutoff", combined=FALSE, 
                      buffer=TRUE, pdf.plot="1.E_cutoff")
  PlotSiteKds("miR-1", sitelist="papercutoff", combined=FALSE, buffer=TRUE,
              pdf.plot="1.F_cutoff")
  PlotSiteOccupancy("miR-1", sitelist="papercutoff", combined=FALSE, buffer=TRUE,
                    pdf.plot="1.G_cutoff")
  message("Done Figure 1")
}


MakeFigure2_Original <- function() {
  message("Making Figure 2")
  PlotSiteKds("let-7a", adjusted_height=TRUE, pdf.plot="2.Ai_original")
  PlotSiteOccupancy("let-7a", adjusted_height=TRUE, pdf.plot="2.Aii_original")
  PlotSiteKds("miR-155", adjusted_height=TRUE, pdf.plot="2.Bi_original")
  PlotSiteOccupancy("miR-155", adjusted_height=TRUE, pdf.plot="2.Bii_original")
  PlotSiteKds("miR-124", adjusted_height=TRUE, pdf.plot="2.Ci_original")
  PlotSiteOccupancy("miR-124", adjusted_height=TRUE, pdf.plot="2.Cii_original")
  PlotSiteKds("lsy-6", adjusted_height=TRUE, pdf.plot="2.Di_original")
  PlotSiteOccupancy("lsy-6", adjusted_height=TRUE, pdf.plot="2.Dii_original")
  PlotSiteKds("miR-7-23nt", exp="equilibrium2_nb", fixed=TRUE,
              adjusted_height=TRUE, pdf.plot="2.Ei_original")
  PlotSiteOccupancy("miR-7-23nt", exp="equilibrium2_nb",
                    fixed=TRUE, adjusted_height=TRUE, pdf.plot="2.Eii_original")
  message("Done Figure 2")
}


MakeFigure2_Cutoff <- function() {
  message("Making Figure 2")
  PlotSiteKds("let-7a", sitelist="papercutoff", adjusted_height=TRUE, pdf.plot="2.Ai_cutoff")
  PlotSiteOccupancy("let-7a", sitelist="papercutoff", adjusted_height=TRUE, pdf.plot="2.Aii_cutoff")
  PlotSiteKds("miR-155", sitelist="papercutoff", adjusted_height=TRUE, pdf.plot="2.Bi_cutoff")
  PlotSiteOccupancy("miR-155", sitelist="papercutoff", adjusted_height=TRUE, pdf.plot="2.Bii_cutoff")
  PlotSiteKds("miR-124", sitelist="papercutoff", adjusted_height=TRUE, pdf.plot="2.Ci_cutoff")
  PlotSiteOccupancy("miR-124", sitelist="papercutoff", adjusted_height=TRUE, pdf.plot="2.Cii_cutoff")
  PlotSiteKds("lsy-6", sitelist="papercutoff", adjusted_height=TRUE, pdf.plot="2.Di_cutoff")
  PlotSiteOccupancy("lsy-6", sitelist="papercutoff", adjusted_height=TRUE, pdf.plot="2.Dii_cutoff")
  PlotSiteKds("miR-7-23nt", exp="equilibrium2_nb", fixed=TRUE,
              sitelist="papercutoff", adjusted_height=TRUE, pdf.plot="2.Ei_cutoff")
  PlotSiteOccupancy("miR-7-23nt", exp="equilibrium2_nb",
                    fixed=TRUE, sitelist="papercutoff", adjusted_height=TRUE, pdf.plot="2.Eii_cutoff")
  message("Done Figure 2")
}

