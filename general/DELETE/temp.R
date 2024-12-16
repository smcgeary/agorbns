  source("general/general.R")

  mirna <- "miR-1"
  experiment <- "equilibrium"
  start <- 5
  stop <- 5
  site <- "8mer"
  sitelist <- "current"
  num.sites <- 20
  log.residual <- FALSE
  combined.input <- TRUE
  bgoff <- FALSE
  model <- TRUE
  site_list <- 10
  connected.points <- FALSE
  params <- GetKds(mirna, experiment, start, stop, sitelist, combined.input=combined.input, log.residual=log.residual, scaled=FALSE)
  params.flanks <- GetFlankKds(mirna, experiment, start, stop, site, sitelist, combined.input=combined.input, log.residual=log.residual,scaled=FALSE)
  if (combined.input == TRUE) {
    sitesXcounts <- GetSitesXCountsCombined(mirna, experiment, start, stop, sitelist)
  } else {
    sitesXcounts <- GetSitesXCounts(mirna, experiment, start, stop, sitelist)
  }
  sitesXcounts <- sitesXcounts[,-1]

  # Remove any rows with no reads whoatsoever.
  sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
  # Remove any rows for which there are no input reads:
  sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]

  kds.s <- unlist(Logistic(params[1:nrow(sitesXcounts)], 10))
  names(kds.s) <- rownames(sitesXcounts)
  s.c <- as.numeric(sitesXcounts[site, ])
  sfXc <- GetSiteFlanksXCounts(mirna, experiment, site, start, stop, sitelist)
  sfXc <- sfXc[,-1]

  colnames(sfXc)[1] <- "I"
  sfXc <- sfXc[rowSums(sfXc[, 2:6]) > 0,]
  sfXc <- sfXc[which(sfXc[,1] > 0),,drop=FALSE]

  sfXc <- round(t((t(sfXc) / colSums(sfXc)) * s.c))
  sfXc[is.na(sfXc)] <- 0
  sitesXcounts.sites <- sitesXcounts[rownames(sitesXcounts) != site,]
  sitesXcounts <- rbind(sitesXcounts, sfXc)
  sitesXcounts <- sitesXcounts[rownames(sitesXcounts) != site, ]
  print(s.c)

# Omit the site kd for which the flanking sites are being fit.
  kds.s <- kds.s[names(kds.s) != site]
  kds <- c(kds.s, Logistic(params.flanks, max=10))

  bgs <- rep(10^params["bg"], 5)
  k.c.stockago <- 10^params["Ago"]

  c.I.tots <- Norm(sitesXcounts[,1])*100
  names(c.I.tots) <- rownames(sitesXcounts)
  data <- sitesXcounts[,2:6]
  c.totals <- matrix(c.I.tots, nrow=nrow(data), ncol = 5, byrow=FALSE)
  colnames(c.totals) <- colnames(data)
  rownames(c.totals) <- rownames(data)
  colors_sites <- kSiteColors[rownames(sitesXcounts)[rownames(sitesXcounts) != site],]
  colors_flanks <- sapply(rownames(sfXc), GetColorFunction)
  colors_all <- c(rep("grey", length(kds.s)), colors_flanks)
  names(colors_all) <- rownames(sitesXcounts)
  c.agos <- sapply(colnames(data), function(x) {
      as.numeric(x) * k.c.stockago / 100
    }
  )

  c.bounds <- as.matrix(
    sapply(c.agos, function(x) {
        return(GetBoundRNA(kds, c.I.tots, x))
      }
    )
  )
  print("percent free")
  print((c.agos - colSums(c.bounds) )/ c.agos)
  c.frees <- c.totals - c.bounds
  c.bgs <- t(t(c.frees) * bgs / colSums(c.frees))

  if (bgoff == TRUE) {
    c.all <- c.bounds
  } else {
    c.all <- c.bounds + c.bgs
  }

  c.final <- data.frame(t(t(c.all) / colSums(c.all)))
  rownames(c.final) <- rownames(data)
  x <- c(40,12.65,4,1.265,0.4)*k.c.stockago/100*1000
  y <- c(1,1,1,1,1)

  sites.norm <- Norm(c.I.tots)
  data.norm <- t(t(data)/colSums(data))

  data.R <- data.norm/(sites.norm)
  model.R <- c.final/(sites.norm)

  xmin <- 0.5*min(x)
  xmax <- 2*max(x)
  ymin <- 0.2
  ymax <- ceiling(max(data.R)*3)
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
  print(site_list)
  if (length(site_list) == 0) {
    site_list_real <- rownames(data)
  } else if (class(site_list) == "character") {
    site_list_real <- c(read.table(file = paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/AssignSiteTypes/sites.", mirna, "_", site_list,".txt"),
  stringsAsFactors=FALSE)[,1], "None",rownames(sfXc))
  } else {
    site_list_real <- c(rownames(sitesXcounts.sites)[order(kds.s)][1:site_list], "None",rownames(sfXc))
  }
  print(site_list_real)
  legend.names <- rownames(data)[order(kds)]

  legend.names <- legend.names[which(legend.names %in% site_list_real)]
  ordered_list <- legend.names
  print(ordered_list)
  centered_sites <- c("11mer-m3.13", "12mer-m3.14",
                  "11mer-m4.14", "12mer-m4.15", "13mer-m3.15", "11mer-m5.15")
  centered_rename <- c("Centered-m3.13", "Centered-m3.14", "Centered-m4.14", "Centered.m4-15", "Centered-m3.15", "Centered-m5.15")
  names(centered_rename) <- centered_sites
  print(kSiteColors[site_list, ])
  print(kSiteColors[ordered_list, ])

  legend.names <- sapply(legend.names, function(site){
                        if (site %in% centered_sites) {
                          return(centered_rename[site])
                        } else {
                          return(site)
                        }
                       }
                       )
  # legend.names <- unique(legend.names)
  print(legend.names)
  legend(x=xmax,y=ymax,legend=legend.names, pch=19, col=kSiteColors[ordered_list, ], cex=0.9, bty="n", ncol = 1)
  for (name in c(site_list_real, rownames(sfXc))) {
    if (connected.points == TRUE) {
      type = "o"
    } else {
      type = "p"
    }
    points(x, data.R[name, ], col=colors_all[name], type = type, pch=19, cex=1.5, lwd=3)
    if (model == TRUE) {
    lines(x, model.R[name, ], col=colors_all[name], lwd=2)      
    }
  }

