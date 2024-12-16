
pt_cex_final <- 1.35
# pch21_lwd.cex <- 0.7 * par()$cex
legend_pch <- 19
line_dash_length <- 2
occupancy_schematic_bg_cex <- 0.45
occupancy_schematic_bg_col <- "#989898"
occupancy_mirna_seq_col <- "#958872"

AGO_mir_label <- "[AGO2-miRNA] (pM)"
AGO_mir1_label <-"[AGO2-miR-1]"
AGO_mir1_label_no_conc <- "AGO2-miR-1 library"



FigureSaveFile <- function(pdf.plot, height=5, width=5, xpos=20, ypos=20) {
  if (class(pdf.plot) == "character") {
    print(pdf.plot)
    # Checks if the first six characters are a date, which means it should go in
    # the "figures" folder.
    if (grepl("^[0-9]*$", substr(pdf.plot, 1, 6))) {
      path_split <- unlist(strsplit(pdf.plot, split="/"))
      date <- path_split[1]
      figure <- path_split[2]
      dir_path <- file.path("figures", date)
      if (!file.exists(dir_path)) {
        dir.create(dir_path)
      }
      path <- paste0(dir_path, "/", figure, ".pdf")
    # Alternative, where the figure goes into the 2017_Paper folder.
    } else {
      figure.split <- unlist(strsplit(pdf.plot, split="\\.", perl=TRUE))
      figure <- figure.split[1]
      if (length(figure.split) > 2) {
        panel <- paste0(figure.split[2:length(figure.split)], collapse=".")
      } else {
        panel <- figure.split[2]    
      }
      figurename <- paste0(figure, panel)
      path <- paste0("2017_Paper/Figure_", figure, "/", figurename, "_raw.pdf")
    }
    pdf(file=path, height=2.2/5*height, width=2.2/5*width, useDingbats=FALSE)
    par(kPDFParameters)
  } else {
    dev.new(xpos=xpos, ypos=ypos, height=height, width=width)
    par(kPlotParameters)    
  }  
}

FigureSaveFile2 <- function(pdf.plot, height=5, width=5, xpos=20, ypos=20) {
  if (class(pdf.plot) == "character") {
    print(pdf.plot)
    # Checks if the first six characters are a date, which means it should go in
    # the "figures" folder.
    if (grepl("^[0-9]*$", substr(pdf.plot, 1, 6))) {
      path_split <- unlist(strsplit(pdf.plot, split="/"))
      date <- path_split[1]
      figure <- path_split[2]
      dir_path <- file.path("figures", date)
      if (!file.exists(dir_path)) {
        dir.create(dir_path)
      }
      path <- paste0(dir_path, "/", figure, ".pdf")
    # Alternative, where the figure goes into the 2017_Paper folder.
    } else {
      figure.split <- unlist(strsplit(pdf.plot, split="\\.", perl=TRUE))
      figure <- figure.split[1]
      if (length(figure.split) > 2) {
        panel <- paste0(figure.split[2:length(figure.split)], collapse=".")
      } else {
        panel <- figure.split[2]    
      }
      figurename <- paste0(figure, panel)
      path <- paste0("ThreePrimeTargetingPaper/Figure_", figure, "/", figurename, "_raw.pdf")
    }
    pdf(file=path, height=2.2/5*height, width=2.2/5*width, useDingbats=FALSE)
    par(kPDFParameters)
  } else {
    dev.new(xpos=xpos, ypos=ypos, height=height, width=width)
    par(kPlotParameters)    
  }  
}

VikramFigureSaveFile <- function(pdf.plot, height=5, width=5, xpos=20, ypos=20) {
  if (class(pdf.plot) == "character") {
    print(pdf.plot)
    path_split <- unlist(strsplit(pdf.plot, split="/"))
    date <- path_split[1]
    figure <- path_split[2]
    dir_path <- file.path("figures", date)
    if (!file.exists(dir_path)) {
      dir.create(dir_path)
    }
    path <- paste0(dir_path, "/", figure, ".pdf")
    pdf(file=path, height=2.2/5*height, width=2.2/5*width, useDingbats=FALSE)
    par(kPDFParameters)
  } else {
    dev.new(xpos=xpos, ypos=ypos, height=height, width=width)
    par(kPlotParameters)    
  }  
}


KineticsScrapFigureSaveFile <- function(pdf.plot, height=5, width=5, xpos=20, ypos=20) {
  if (class(pdf.plot) == "character") {
    print(pdf.plot)
    path <- paste0("KineticsPaper/Scrap/", pdf.plot, ".pdf")
    pdf(file=path, height=2.2/5*height, width=2.2/5*width, useDingbats=FALSE)
    par(kPDFParameters)
  } else {
    dev.new(xpos=xpos, ypos=ypos, height=height, width=width)
    par(kPlotParameters)    
  }  
}

PlotSiteEnrichments_ <- function(mirna, experiment="equilibrium", n_constant=5,
                                sitelist="paperfinal", plotlist=FALSE,
                                uniq=FALSE, combined=TRUE, L=FALSE,
                                compcorrect=FALSE, wobble=FALSE, tpomit=FALSE,
                                buffer=FALSE, singleonly=TRUE, minkd=FALSE,
                                plot.nconstant=FALSE, flowthrough=FALSE,
                                added.text=FALSE, datalines=FALSE,
                                modellines=TRUE, leg_rescale=1, height=4.5,
                                dil_AGO=1, dil_bg=1,
                                width=8.1, xpos=20, ypos=20, pdf.plot=FALSE) {
  # Load the count data and the model parameters:
  sXc <- SubfunctionCall(SitesXCounts)
  pars.matrix <- SubfunctionCall(EquilPars)
  # Correct mirna name for miR-7.
  if (mirna %in% c("miR-7-23nt", "miR-7-24nt", "miR-7-25nt")) {
    mirna_temp <- "miR-7"
  }
  # Removed the sites that are ascribed to competitor oligo enrichment.
  removed_sites <- GetRemovedSites(sXc)
  print(removed_sites)
  if (length(removed_sites) != 0) {
    # sXc <- sXc[setdiff(rownames(sXc), removed_sites), , drop=FALSE]
    kd.sites <- gsub("_Kd", replace="", rownames(pars.matrix))
    print(kd.sites)
    # pars.matrix <- pars.matrix[setdiff(rownames(pars.matrix), paste0(removed_sites, "_Kd")), ,drop=FALSE]
  }

  pars.model <- log10(pars.matrix$Mean)
  names(pars.model) <- rownames(pars.matrix)
  names(pars.model)[nrow(sXc) + 1] <- "bg"
  names(pars.model)[nrow(sXc) + 2] <- "AGO"



  data <- GetDataEquil(sXc, combine_reps=TRUE)
  l <- SubfunctionCall(GetInputEquil)
  if (L) {
    l <- l/sum(l) * as.numeric(L)
  }

  print(data)

  # Ago dilutio in the data:
  pars <- pars.model
  A.dil.data <- sapply(colnames(data), as.numeric)
  # A.stock.measured <- kAgoStock[mirna, experiment]
  # print(A.stock.measured)
  A.stock.pars <- 10^pars.model["AGO"]
  A.stock.pars <- A.stock.pars*dil_AGO
  pars.model["AGO"] <- log10(A.stock.pars)
  pars.model["bg"] <- log10(10^pars.model["bg"]*dil_bg)
  pars <- pars.model
  # pM_from_dil <- A.stock.measured*1000/100
  pM_from_dil <- A.stock.pars*1000/100
  A.pM.data <- A.dil.data*pM_from_dil
  # A.stock.measured <<- A.stock.measured
  xmin <- signif(min(A.pM.data)/(10^0.25), 1)
  xmax <- signif(max(A.pM.data)*(10^0.5), 1)
  print(A.pM.data)
  print(xmin)
  print(xmax)
  print(pM_from_dil)
  A.dil.model <- exp(seq(log(xmin/pM_from_dil),
                         log(xmax/pM_from_dil),
                     length=100))
  model <- SubfunctionCall(EquilSingleSiteModelFreq, A.dil=A.dil.model)

  A.pM.model <-A.dil.model*A.stock.pars *1000/100
  data.R <- EquilEnrichments(data, l)
  model.R <- EquilEnrichments(model, l)
  # Set up the plotting limits
  if (flowthrough) {
    ymin <- 0.01
    ymax <- 1
  } else {
    ymin <- 0.2
    ymax <- 300
  }
  SubfunctionCall(GenericFigureSaveFile)
  # par(mar= c(3, 0.5, 2, 2))

  BlankPlot(log="xy", adjusted=TRUE)
  # Generate tickmarks for axis.
  xy <- GetPlotFractionalCoords(fx=0.05, fy=0.95, log='xy')
  title.text <- mirna

  if (added.text) {
    title.text <- paste0(mirna, "\n", experiment)
  }

  text(xy[1], xy[2], labels = title.text, adj=c(0, 1))

  mirna.trim <- paste0(strsplit(mirna, split = "-")[[1]][1:2],collapse="-")
  names(pars.model) <- gsub("(.*)_Kd", names(pars.model), replace="\\1")

  if (plot.nconstant) {
    xy <- GetPlotFractionalCoords(0.05, 0.90, log='xy')
    text(xy[1], xy[2], labels=sprintf("%s nt constant region", n_constant), adj=c(0, 1))      
  }

  xy <- GetPlotFractionalCoords(0.05, 0.87, log='xy')
  text(xy[1], xy[2], labels=experiment, adj=c(0, 1))      


  if (minkd) {
    xy <- GetPlotFractionalCoords(0.05, 0.80, log='xy')
    text(xy[1], xy[2], labels=sprintf("Min Kd: %.0e", minkd), adj=c(0, 1))      
  }



  site_order <- order(pars.matrix$Mean[1:(nrow(data)-1)])
  kd.matrix <- pars.matrix[1:(nrow(data)-1),][site_order, ]
  sites.ordered <- rownames(data)[site_order]
  kds <- kd.matrix$Mean

  # Part where text is formatted:
  kd_mags <- floor(log10(kds))
  # error_upper <- (kd.matrix$Upper_CI-kds)/(10^kd_mags)
  error_lower <- (kds - kd.matrix$Lower_CI)/(10^kd_mags)
  # error_mags_upper <- floor(log10(error_upper))
  error_mags_lower <- floor(log10(error_lower))
  temp_matrix <- rbind(kds/(10^kd_mags), c(error_lower), -error_mags_lower)
  print(temp_matrix)
  cols <- kSiteColors[sites.ordered]
  names(cols) <- sites.ordered

  xy <- GetPlotFractionalCoords(fx=0.85, fy=1, log='xy')
  xy <- GetPlotFractionalCoords(fx=1, fy=1, log='xy')

  temp_legend <- Legend(xy,
                        legend=c(ConvertTtoUandMmtoX(sites.ordered), "None"),
                        col=c(cols, "black"),
                        y.intersp=0.9*leg_rescale)

  if (pdf.plot=="1.D" | pdf.plot=="1.D_uniq") {
    xy <- GetPlotFractionalCoords(fx=0.85, fy=0, log='xy')
    temp_legend <- legend(x=xy[1], y=xy[2], legend=rep("", length(sites.ordered) + 1),
                          col=c(kSiteColors[sites.ordered], "black"), seg.len=1,
                          pch=NA, lwd=1, bty="n",
                          y.intersp=0.86, yjust=0)

    xy <- GetPlotFractionalCoords(fx=0.85, fy=0.46, log='xy')
    text(x=xy[1], y=xy[2], labels=bquote("Relative"~italic(K)[D]*":"),
         adj=c(0, 0))
    pos1 <- 0.15
    text(10^(temp_legend$text$x + pos1),
         10^temp_legend$text$y[1:(length(sites.ordered) + 1)],
         sprintf("%.*f", c(temp_matrix[3, ], 1), c(temp_matrix[1, ], 1)), adj=1)
    pos2 <- pos1 + 0.15
    text(10^(temp_legend$text$x + pos2)[1:length(sites.ordered)],
         10^temp_legend$text$y[1:length(sites.ordered)],
         sapply(kd_mags, function(kd_mag) {
          as.expression(bquote(.("")%+-%.("")))
          }), adj=1)
    pos3 <- pos2 + 0.23
    text(10^(temp_legend$text$x + pos3)[1:length(sites.ordered)],
         10^temp_legend$text$y[1:length(sites.ordered)],
         sprintf("%.*f", temp_matrix[3, ], temp_matrix[2, ]), adj=1)
      pos4 <- pos3 + 0.46
    text(10^(temp_legend$text$x + pos4)[1:length(sites.ordered)],
         10^temp_legend$text$y[1:length(sites.ordered)],
         sapply(kd_mags, function(kd_mag) {
          as.expression(bquote(.("")%*%10^.(kd_mag)))
          }), adj=c(1, 0.4))
  }
  print(cols)
  cols <- c(cols, None="black")
  print(cols)
  if (flowthrough) {
    segments(A.pM.data, ymin, A.pM.data, ymax, lty=2, col="gray")
  }
  sapply(rownames(data), function(site) {
    if (!flowthrough) {
      Points(A.pM.data, data.R[site, ], col=cols[site], line=datalines)
    }
    if (!(datalines) & modellines) {
      lines(A.pM.model, model.R[site, ], col=cols[site], xpd=NA)        
    }
  })
  AddLogAxis(1, label=AGO_mir_label, adj=TRUE)
  AddLogAxis(2, label="Enrichment")

  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

PlotSiteKds_ <- function(mirna, experiment="equilibrium", n_constant=5,
                        sitelist="paper", height=5, width=6,
                        pdf.plot=FALSE) {
  # Get kds for all site-types of the mirna.
  sXc <- SubfunctionCall(SitesXCounts)
  pars.matrix <- SubfunctionCall(GetSiteKds)
  kds <- pars.matrix[1:nrow(sXc),]
  rownames(kds) <- rownames(sXc)
  kds <- kds[order(kds$Mean),]
  sites <- rownames(kds)
  SubfunctionCall(GenericFigureSaveFile)
  xmin <- 3e-5
  xmax <- 3
  ymin <- 0
  ymax <- length(sites) + 0.5
  BlankPlot(log='x', inv='x')
  y <- nrow(kds) - seq(nrow(kds)) + 1
  arrows(kds$Upper_CI, y, kds$Lower_CI, y, length=0.05*par()$cex, angle=90,
         code=3)
  bg.colors <- sapply(sites, function(site) {
    if (site%in% names(kSiteCategories)) {
      return(kSiteCategoryColors[kSiteCategories[site]])
    } else {
      return(kSiteCategoryColors["Noncanonical"])
    }
  })
  points(kds$Mean,y, pch=21, bg=c(bg.colors), cex=1.2,
         lwd=par()$cex)
  AddLogAxis(1, label="Kd")
  title.xy <- GetPlotFractionalCoords(fx=0.05, fy=0.95, log='x', inv='x')
  text(title.xy[1], title.xy[2], labels=mirna, adj=0)
  text((kds$Lower_CI)/1.1, y, labels=ConvertTtoU(sites), adj=0)
  legend.xy <- GetPlotFractionalCoords(fx=0.60, fy=0.40, log='x', inv='x')
  legend.names <- c("7-8-nt canonical site",
                    "6-nt canonical site",
                    "Enhanced 6mer-containing sites",
                    "Noncanonical sites",
                    "3' sites",
                    "?")
  legend.colors <- c(kSiteCategoryColors[c(1, 2, 3, 7, 4, 5)])
  mirna.sites <- list(c(1, 2, 4, 6),
                      c(1, 2, 4, 6),
                      c(1, 2, 4, 5, 6),
                      c(1, 2, 3, 4, 5),
                      c(1, 2, 4, 5, 6),
                      c(1, 2, 3, 4))

  names(mirna.sites) <- kMirnas
  # legend(x=legend.xy[1], y=legend.xy[2],
  #        legend=legend.names[mirna.sites[[mirna]]],
  #        bty="n",
  #        col="black",
  #        pt.bg=legend.colors[mirna.sites[[mirna]]],
  #        pch=21,
  #        pt.cex=1.5,
  #        pt.lwd=par()$cex)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}


PlotSiteKds_ <- function(mirna, experiment="equilibrium", n_constant=5,
                        sitelist="paperfinal", combined=TRUE, uniq=FALSE,
                        singleonly=TRUE, papersites=FALSE, collapse_AA=TRUE,
                        adjusted_height=FALSE, added.text=FALSE, L=FALSE,
                        trim_mir_name=TRUE, buffer=FALSE, compcorrect=FALSE,
                        tpomit=FALSE, wobble=FALSE, plot.nconstant=FALSE,
                        height=5, width=6, xpos=20, ypos=20, pdf.plot=FALSE) {
  # Get kds for all site-types of the mirna.
  sXc <- SubfunctionCall(SitesXCounts)
  pars.matrix <- SubfunctionCall(EquilPars)
  kd.sites <- gsub("_Kd", replace="", rownames(pars.matrix))

  # Assign removed_sites
  removed_sites <- GetRemovedSites(sXc)
  removed_sites <- c()
  if (collapse_AA) {
    # removed_sites <- GetRemovedSites(sXc)
    AA_sites <- grep("^AA-", rownames(sXc), perl=TRUE, value=TRUE)
    if (length(AA_sites) != 0) {
      # Get the inds o
      AA_inds <- grep("^AA-", rownames(sXc), perl=TRUE)
      # Get the base names of all AA sites.
      AA_sites_base <- gsub("^(AA-)(.*)$", AA_sites, replace="\\2", perl=TRUE)
      # Get the AA sites that are present as base sites in rownames(sXc)
      AA_sites_keep <- AA_sites[AA_sites_base %in% rownames(sXc)]
      # get the base sites that are present as base sites
      sites_base_keep <- AA_sites_base[AA_sites_base %in% rownames(sXc)]
      inds_keep_base <- sapply(sites_base_keep, function(site) {
        which(rownames(sXc) == site)
      })
      AA_inds_keep <- sapply(AA_sites_keep, function(site) {
        which(rownames(sXc) == site)
      })      
    } else {
      AA_sites_keep <- c()
      inds_keep_base <- c()
      AA_inds_keep <- c()
    }
    sXc_AA <- sXc[AA_sites_keep, ]
    removed_sites <- c(AA_sites_keep, removed_sites)
    # sXc_AA <- sXc[AA_sites_keep, ]
    pars.matrix_AA <- pars.matrix[inds_keep_base, ]
    print(pars.matrix_AA)
    if (length(AA_sites) !=0) {
      rownames(sXc_AA) <- sites_base_keep
      pars.matrix[inds_keep_base, 1:5] <- pars.matrix[AA_inds_keep, 1:5]
      rownames(pars.matrix_AA) <- sites_base_keep      
    }
  } else {
    removed_sites <- GetRemovedSites(sXc)
  }
  if (length(removed_sites) != 0) {
    sXc <- sXc[setdiff(rownames(sXc), removed_sites), , drop=FALSE]
    rownames(pars.matrix) <- kd.sites
    pars.matrix <- pars.matrix[setdiff(kd.sites, removed_sites), ,drop=FALSE]
    if (sitelist %in% c("paperfinal", "paper")) {
    }
  }
  kds <- pars.matrix[1:(nrow(sXc) - 1),]
  rownames(kds) <- rownames(sXc)[-nrow(sXc)]
  order_all <- order(kds$Mean)
  kds <- kds[order_all,]
  sites <- rownames(kds)
  xmin <- 8e-5  
  xmax <- 1
  ymin <- 0
  ymax <- length(sites) + 0.5
  y <- nrow(kds) - seq(nrow(kds)) + 1
  names(y) <- sites
  if (adjusted_height) {
    height <- (ymax*0.68 + 5)/4
  }
  SubfunctionCall(GenericFigureSaveFile)
  par(mar=c(3, 0.5, 2, 2))
  BlankPlot(log='x', inv='x')
  segments(1, ymin, 1, max(y) - 1.5, xpd=NA)
  xmin <- 1e-4  
  ind_none <- which(sites == "None")
  x0 <- kds$Upper_CI
  x0[ind_none] <- NA
  x1 <- kds$Lower_CI
  x1[ind_none] <- NA
  x_error <- list(x0, x1)
  cols <- sapply(sites, function(site) {
    print(site)
    if (site %in% names(kSiteCategories)) {
      return(kSiteCategoryColors[kSiteCategories[site]])
    } else {
      return(kSiteCategoryColors["Noncanonical"])
    }
  })
  if (collapse_AA) {
    if (nrow(pars.matrix_AA) != 0) {
    y_AA <- y[rownames(sXc_AA)]
    kds_error <- list(pars.matrix_AA$Lower_CI, pars.matrix_AA$Upper_CI)
      colors_AA <- sapply(rownames(pars.matrix_AA), function(site) {
        if (site %in% names(kSiteCategories)) {
          return(kSiteCategoryColors[kSiteCategories[site]])
        } else {
          return(kSiteCategoryColors["Noncanonical"])
        }
      })
    }
  }
  if (collapse_AA) {
    if(nrow(pars.matrix_AA) != 0) {
      sapply(rownames(pars.matrix_AA),function(site) {
        kd <- pars.matrix_AA[site, ]
        y <- y_AA[site]
        label <- ConvertTtoUandMmtoX(site)
        x_seg_l <- kd[2]
        x_seg_r <- pars.matrix[site, ][2]

        if (site %in% c("8mer-bT6", "8mer-bC(4.6)")) {
          x <- min(kd[2]/1.6, kd[3]/1.1)
          adj <- 0
          x_seg_l <- kd[2]
          segments(unlist(x_seg_l), y, x, y, lty=2, col="blue")
          str.width <- strwidth(label, units="figure")
          pos_start <- (log(x) - log(xmax))/(log(xmin) - log(xmax))
          pos_end <- pos_start + str.width + 0.04
          x_r <- GetPlotFractionalCoords(pos_end, 1, log='x', inv='x')
          segments(x_r[1], y, unlist(x_seg_r), y, lty=2)
        } else {
          x <- kd[5]* 1.1
          adj <- 1
          segments(unlist(x_seg_l), y, unlist(x_seg_r), y, lty=2)
        }
        text_out <- text(x, y, labels=label, adj=adj)
      })  
      Points(pars.matrix_AA$Mean, y_AA, col=colors_AA, x_error=kds_error)
    }
  }
  Points(kds$Mean, y, col=cols, x_error=x_error)
  AddLogAxis(1, label="Relative Kd")
  title.xy <- GetPlotFractionalCoords(fx=0.05, fy=1, log='x', inv='x')
    if (length(grep("miR-7", mirna)) > 0 & trim_mir_name) {
      title.text <- "miR-7"
    } else {
      title.text <- mirna
    }
  if (added.text) {
    title.text <- paste0(mirna, "\n", experiment)
  }
  if (plot.nconstant) {
    xy <- GetPlotFractionalCoords(0.05, 0.90, inv='x')
    points(c(xy[1]), c(xy[2]), pch=3, col="gray")
    text(xy[1], xy[2], labels=sprintf("%s nt constant region", n_constant), adj=c(0, 1))      
  }
  print(376)
  text(title.xy[1], max(y), labels = title.text, adj=c(0, 0.5), xpd=NA)
  title_pos <- apply(kds, 1, function(kd) {
    min(kd[2]/1.6, kd[3]/1.1)
  })
  if (collapse_AA) {
    if (nrow(pars.matrix_AA) !=0) {
      AA_inds <- which(sites %in% rownames(pars.matrix_AA))
      sites[AA_inds] <- paste0("AA-", sites[AA_inds], sep="")  
    }    
  }
  par(xpd=NA) 
  text(title_pos, y, labels=ConvertTtoUandMmtoX(sites), adj=0)
  legend.names <- c("7-8-nt canonical site",
                    "6-nt canonical site",
                    "Enhanced 6mer site",
                    "Noncanonical site",
                    "3'-only site")
  legend.colors <- c(kSiteCategoryColors[c(1, 2, 3, 7, 4)])
  mirna.sites <- list(c(1, 2, 4), # miR-1
                      c(1, 2, 4), # let-7a
                      c(1, 2, 3, 4, 5), # miR-155
                      c(1, 2, 3, 4, 5), # miR-124
                      c(1, 2, 3, 4, 5), # lsy-6
                      c(1, 2, 3, 4), # miR-7-22nt
                      c(1, 2, 3, 4), # miR-7-23nt
                      c(1, 2, 3, 4), # miR-7-24nt
                      c(1, 2, 3, 4)) # miR-7-25nt
  names(mirna.sites) <- c(kMirnas, "miR-7-22nt", "miR-7-24nt", "miR-7-25nt")
  xy <- GetPlotFractionalCoords(fx=0.65, fy=0, log='x', inv='x')      
  if (!(added.text)) {
    Legend(xy, legend=legend.names[mirna.sites[[mirna]]],
           col=legend.colors[mirna.sites[[mirna]]],
           yjust=0)    
  }
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}




PlotPairwiseKds <- function(mirna1, experiment1, n_constant1, sitelist1, 
                            mirna2, experiment2, n_constant2, sitelist2,
                            combined1=TRUE, buffer1=FALSE, combined2=TRUE, 
                            buffer2=FALSE, mirna_start1=FALSE,
                            mirna_start2=FALSE, tpomit1=FALSE, tpomit2=FALSE,
                            minkd1=FALSE, minkd2=FALSE,
                            ident=FALSE, xmin=8e-4, pdf.plot=FALSE,
                            xpos=20, ypos=20) {
  # Get kds for all site-types of the mirna.
  mirna <- mirna1
  experiment <- experiment1
  sitelist <- sitelist1
  n_constant <- n_constant1
  combined <- combined1
  buffer <- buffer1
  mirna_start <- mirna_start1
  tpomit <- tpomit1
  minkd <- minkd1
  kds1 <- SubfunctionCall(EquilPars)
  if (grepl("miR-7", mirna1)) {
    rownames(kds1)[(nrow(kds1) - 1):nrow(kds1)] <- c("bg_miR-7", "AGO_miR-7")
  }
  mirna <- mirna2
  experiment <- experiment2
  sitelist <- sitelist2
  n_constant <- n_constant2
  combined <- combined2
  buffer <- buffer2
  mirna_start <- mirna_start2
  tpomit <- tpomit2
  minkd <- minkd2
  kds2 <- SubfunctionCall(EquilPars)
  print(kds2)
  if (grepl("miR-7", mirna2)) {
    rownames(kds2)[(nrow(kds2) - 1):nrow(kds2)] <- c("bg_miR-7", "AGO_miR-7")
    xmin <- 2e-4
  }

  kds2 <- kds2[rownames(kds1), ]

  print(kds1)
  print(kds2)
  sites <- rownames(kds1)
  sites[length(sites) - 1] <- "bg"
  sites[length(sites)] <- "AGO"
  SubfunctionCall(GenericFigureSaveFile)
  kSiteColors["bg"] <- "brown"
  kSiteColors["AGO"] <- "green"

  # xmin <- 3e-7
  xmax <- 10
  ymin <- xmin
  ymax <- xmax
  BlankPlot(log="xy", inv="xy", adjusted=TRUE)

  x <- kds1$Mean
  y <- kds2$Mean
  arrows(kds1$Upper_CI, y, kds1$Lower_CI, y, length=0.05*par()$cex, angle=90,
         code=3)
  arrows(x, kds2$Upper_CI, x, kds2$Lower_CI, length=0.05*par()$cex, angle=90,
         code=3)
  bg.colors <- kSiteColors[sites]
  points(x=x, y=y, cex=1.2,
         lwd=par()$cex, col=kSiteColors[sites])
  print(sites)
  print(kSiteColors[sites])
  abline(0, 1, lty=2)
  AddLogAxis(1, label=paste(mirna1, experiment1, n_constant1, sitelist1,
                            collapse=" "), adj=TRUE)
  AddLogAxis(2, label=paste(mirna2, experiment2, n_constant2, sitelist2,
                            collapse=" "), adj=TRUE)
  if (ident) {
    identify(x, y, labels=ConvertTtoUandMmtoX(sites))
  }
  title.xy <- GetPlotFractionalCoords(fx=0.05, fy=0.95, log='x', inv='x')
  xy <- GetPlotFractionalCoords(fx=0.05, fy=0.95, log='xy', inv='xy')
  legend(x=xy[1], y=xy[2], legend=c("[AGO2-miRNA]", "[Background]"), 
       pch=19, col=c("green", "brown"), bty="n")
  xy <- GetPlotFractionalCoords(fx=0.8, fy=0.05, log='xy', inv='xy')
  AddCorrelationToPlot(x[1:(length(x) - 2)], y[1:(length(x) - 2)], xy[1], xy[2],
                       rsquared=TRUE)
  # legend(x=xy[1], y=xy[2], legend=c("23", "24", "25"), col = c("blue", "red", "green"),
  #        pch=19, bty="n", cex=0.9)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

PlotPairwiseInput <- function(mirna1, experiment1, n_constant1, sitelist1, 
                            mirna2, experiment2, n_constant2, sitelist2,
                            combined1=TRUE, buffer1=FALSE, combined2=TRUE, 
                            buffer2=FALSE, mirna_start1=FALSE, mirna_start2=FALSE,
                            ident=FALSE, pdf.plot=FALSE, ...) {
  # Get kds for all site-types of the mirna.
  mirna <- mirna1
  experiment <- experiment1
  sitelist <- sitelist1
  n_constant <- n_constant1
  combined <- combined1
  buffer <- buffer1
  mirna_start <- mirna_start1
  kds1 <- SubfunctionCall(EquilPars)
  sXc1 <- SubfunctionCall(SitesXCounts)
  sXc1 <<- sXc1

  mirna <- mirna2
  experiment <- experiment2
  sitelist <- sitelist2
  n_constant <- n_constant2
  combined <- combined2
  buffer <- buffer2
  mirna_start <- mirna_start2
  kds2 <- SubfunctionCall(EquilPars)
  sXc2 <- SubfunctionCall(SitesXCounts)
  sXc2 <<- sXc2


  sites <- rownames(kds1)
  sites[length(sites) - 1] <- "bg"
  sites[length(sites)] <- "AGO"
  SubfunctionCall(FigureSaveFile)


  kSiteColors["bg"] <- "brown"
  kSiteColors["AGO"] <- "green"


  xmin <- 1e-7
  xmax <- 1
  ymin <- xmin
  ymax <- xmax
  BlankPlot(log="xy", inv="xy", adjusted=TRUE)

  print(kds1)
  x <- Norm(sXc1[, 1])
  y <- Norm(sXc2[, 2])
  print(x)
  print(y)

  # arrows(kds1$Upper_CI, y, kds1$Lower_CI, y, length=0.05*par()$cex, angle=90,
  #        code=3)
  # arrows(x, kds2$Upper_CI, x, kds2$Lower_CI, length=0.05*par()$cex, angle=90,
  #        code=3)
  bg.colors <- kSiteColors[sites]
  points(x=x, y=y, cex=1.2,
         lwd=par()$cex, col=kSiteColors[sites])
  abline(0, 1, lty=2)
  AddLogAxis(1, label=paste(mirna1, experiment1, n_constant1, sitelist1, collapse=" "), adj=TRUE)
  AddLogAxis(2, label=paste(mirna2, experiment2, n_constant2, sitelist2, collapse=" "), adj=TRUE)
  title.xy <- GetPlotFractionalCoords(fx=0.05, fy=0.95, log='x', inv='x')
  # text(title.xy[1], title.xy[2], labels=mirna, adj=0)
  # text((kds$Full)/1.6, y, labels=ConvertTtoU(sites), adj=0)
  xy <- GetPlotFractionalCoords(fx=0.5, fy=0.95, log='xy', inv='xy')
  legend(x=xy[1], y=xy[2], legend=ConvertTtoU(sites),
         pch=19, col=bg.colors, bty="n", cex=0.9)
  legend(x=xy[1], y=xy[2], legend=c("23", "24", "25"), col = c("blue", "red", "green"),
         pch=19, bty="n", cex=0.9)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

MakeCompetitorPlot <- function(mirna, experiment="equilibrium", n_constant=5,
                               kmer_len_start=4, kmer_len_stop=16, off=0,
                               alt_mirna_comp=FALSE, geomean=FALSE, xpos=20,
                               ypos=20, height=5, width=5, pdf.plot=FALSE) {
  mirna_base <- paste0(unlist(strsplit(mirna, split="-"))[1:2], collapse="-")
  if (mirna == "miR-7-23nt" & experiment == "equilibrium2_nb") {
    conditions <- c("I", "0", "0.4", "1.26", "12.6", "40")
  } else if (mirna == "miR-124" & experiment == "equilibrium_tp") {
    conditions <- c("I", "0", "0.04", "0.126", "0.4", "1.26", "4", "12.6", "40")
  } else if (mirna == "miR-124" & experiment == "equilibrium_2_tp") {
    conditions <- c("I,1", "I,2", "0,1", "0,2", "0.126,1", "0.126,2",
                    "0.4,1", "0.4,2", "1.26,1", "1.26,2", "4,1", "4,2", "12.6,1",
                    "12.6,2", "40,1", "40,2")
    conditions <- c("I", "0", "0.4", "1.26", "4", "12.6", "40")
  } else if (mirna == "miR-1" & (experiment == "equilibrium_tp" ||
                                 experiment == "equilibrium_met_tp")) {
    conditions <- c("I", "0", "0.126", "0.4", "1.26", "4", "12.6", "40")
  } else {
    conditions <- c("I", "0", "0.4", "1.26", "4", "12.6", "40")
  }
  print("In function")
  Rs <- sapply(kmer_len_start:kmer_len_stop, function(kmer_len) {
    counts_outer <- sapply(conditions, function(condition) {
      # Get the path to file.
      # path <- GetAnalysisPath(mirna, experiment, condition, analysis_type="reads")
      # Get the total number of kmers.
      counts <- CompetitorKmers(mirna, condition, kmer_len, off=off,
                                experiment=experiment,
                                alt_mirna_comp=alt_mirna_comp)
      total_counts <- counts["Total", ]
      kmers <- unlist(counts)
      fracs <- (kmers[-length(kmers)] + 1)/(total_counts + 4**kmer_len)

      names(fracs) <- rownames(counts)[-length(kmers)]
      fracs
    })
    # if (mirna == "miR-124" & experiment == "equilibrium_2_tp") {
    #   counts_outer <- (
    #     counts_outer[, c(1, 3, 5, 7, 9, 11, 13, 15)] +
    #     counts_outer[, c(2, 4, 6, 8, 10, 12, 14, 16)]
    #   )/2
    #   colnames(counts_outer) <- c("I", "0", "0.126", "0.4", "1.26", "4", "12.6",
    #                               "40")
    # }
    R <- counts_outer/counts_outer[, "I"]
    if (geomean) {
      R_vec <- apply(R, 2, function(column) {
        GeoMean(column[is.finite(column)])
      })
    } else {
      R_vec <- apply(R, 2, function(column) {
        mean(column[is.finite(column)])
      })      
    }
    R_vec
  })
  colnames(Rs) <- sprintf("k%s", seq(kmer_len_start, kmer_len_stop))
  SEs <- sapply(kmer_len_start:kmer_len_stop, function(kmer_len) {
    counts_outer <- sapply(conditions, function(condition) {
      # Get the path to file.
      path <- GetAnalysisPath(mirna, experiment, condition, analysis_type="reads")
      # Get the total number of kmers.
      counts <- CompetitorKmers(mirna, condition, kmer_len, off=off,
                                experiment=experiment)
      total_counts <- counts["Total", ]
      kmers <- unlist(counts)
      fracs <- kmers[-length(kmers)]/total_counts
      names(fracs) <- rownames(counts)[-length(kmers)]
      fracs
    })
    # if (mirna == "miR-124" & experiment == "equilibrium_2_tp") {
    # counts_outer <- (
    #   counts_outer[, c(1, 3, 5, 7, 9, 11, 13, 15)] +
    #   counts_outer[, c(2, 4, 6, 8, 10, 12, 14, 16)]
    # )/2
    # colnames(counts_outer) <- c("I", "0", "0.126", "0.4", "1.26", "4", "12.6",
    #                             "40")
    # }
    R <- counts_outer/counts_outer[, "I"]
    if (geomean) {
      SE_vec <- apply(R, 2, function(column) {
        exp(sd(log(column[is.finite(column)]))/
            sqrt(length(column[is.finite(column)])))
      })
    } else {
      SE_vec <- apply(R, 2, function(column) {
        sd(column[is.finite(column)])/sqrt(length(column[is.finite(column)]))
      })      
    }
    SE_vec
  })

  # if (mirna == "miR-124" & experiment == "equilibrium_2_tp") {
  #   conditions <- c("I", "0", "0.126", "0.4", "1.26", "4", "12.6", "40")
  # }
  xmin <- 3
  xmax <- 17
  ymin <- 0
  if (mirna_base == "miR-7") {
    ymax <- 120
    y.tick.space <- 5
    y.label.space <- 20
  } else {
    ymax <- 50
    y.tick.space <- 1
    y.label.space <- 5
  }
  ymax <- 50
  y.tick.space <- 1
  y.label.space <- 5

  SubfunctionCall(FigureSaveFile)
  BlankPlot()
  conditions <- conditions[-2]
  Rs["I", ] <- 0
  ys <- c(Rs[conditions, ])
  xs <- seq(kmer_len_start - 0.5,
            to=kmer_len_stop + 0.5 - 1/length(conditions),
            length.out=length(ys))
  errors <- c(SEs[conditions, ])
  if (geomean) {
    y_top <- ys*errors
    y_bottom <- ys/errors
  } else {
    y_top <- ys + errors
    y_bottom <- ys - errors
  }
  arrows(x0=xs + 1/(2*length(conditions)), y0=y_bottom, y1=y_top,
         length=error_height_final*par()$cex/4, lwd=par()$cex, angle=90, code=3)

  rect(xleft=xs, ybottom=0, xright=xs + 1/length(conditions), ytop=ys,
       border=NA, col=kEquilCondColors[conditions])
  AddLinearAxis(1, tick.space=1, label.space=1,
                label="Length of complementarity")
  AddLinearAxis(2, tick.space=y.tick.space, label.space=y.label.space,
                label="Average enrichment within positions 14-31")
  # for (condition in conditions) {
  #   lines(kmer_len_start:kmer_len_stop, Rs[condition, ], type="o",
  #         col=kEquilCondColors[condition])
  # }
  conditions <- conditions[-1]
  xy <- GetPlotFractionalCoords(0.05, 0.95)
  legend(x=xy[1], y=xy[2], bty="n", legend=conditions, pch=19,
         col=kEquilCondColors[conditions])
  xy <- GetPlotFractionalCoords(0.95, 0.95)
  text(xy[1], xy[2], labels=mirna, adj=1)
  xy <- GetPlotFractionalCoords(0.95, 0.90)
  text(xy[1], xy[2], labels=experiment, adj=1)
  xy <- GetPlotFractionalCoords(0.95, 0.85)
  if (geomean) string <- "Geometric mean"
  else         string <- "Arithmetic mean"
  text(xy[1], xy[2], labels=string, adj=1)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}


MakeCompetitorPlotSingle <- function(mirna, kmer_len, experiment="equilibrium",
                                     n_constant=5, kmer_len_stop=16, off=0,
                                     geomean=FALSE, xpos=20, ypos=20,
                                     height=5, width=5, pdf.plot=FALSE) {
  mirna_base <- paste0(unlist(strsplit(mirna, split="-"))[1:2], collapse="-")
  if (mirna == "miR-7-23nt" & experiment == "equilibrium2_nb") {
    conditions <- c("I", "0", "0.4", "1.26", "12.6", "40")
  } else if (mirna == "miR-124" & experiment == "equilibrium_tp") {
    conditions <- c("I", "0", "0.04", "0.126", "0.4", "1.26", "4", "12.6", "40")
  } else if (mirna == "miR-124" & experiment == "equilibrium_2_tp") {
    conditions <- c("I,1", "I,2", "0,1", "0,2", "0.126,1", "0.126,2",
                    "0.4,1", "0.4,2", "1.26,1", "1.26,2", "4,1", "4,2", "12.6,1",
                    "12.6,2", "40,1", "40,2")
    conditions <- c("I", "0", "0.4", "1.26", "4", "12.6", "40")
  } else {
    conditions <- c("I", "0", "0.4", "1.26", "4", "12.6", "40")
  }
  counts_outer <- sapply(conditions, function(condition) {
    # Get the path to file.
    path <- GetAnalysisPath(mirna, experiment, condition, analysis_type="reads")
    # Get the total number of kmers.
    counts <- CompetitorKmers(mirna, condition, kmer_len, off=off,
                              experiment=experiment)
    total_counts <- counts["Total", ]
    kmers <- unlist(counts)
    fracs <- kmers[-length(kmers)]/total_counts
    names(fracs) <- rownames(counts)[-length(kmers)]
    fracs
  })
  if (mirna == "miR-124" & experiment == "equilibrium_2_tp") {
    counts_outer <- (
      counts_outer[, c(1, 3, 5, 7, 9, 11, 13, 15)] +
      counts_outer[, c(2, 4, 6, 8, 10, 12, 14, 16)]
    )/2
    colnames(counts_outer) <- c("I", "0", "0.126", "0.4", "1.26", "4", "12.6",
                                "40")
  }
  R <- counts_outer/counts_outer[, "I"]
  print(R)


  xmin <- 1
  xmax <- nrow(R)
  ymin <- 0
  if (mirna_base == "miR-7") {
    ymax <- 120
    y.tick.space <- 5
    y.label.space <- 20
  } else {
    ymax <- 50
    y.tick.space <- 1
    y.label.space <- 5
  }
  SubfunctionCall(GenericFigureSaveFile)
  BlankPlot()
  conditions <- conditions[-2]
  xs <- seq(xmin, xmax)
  sapply(conditions, function(condition) {
    lines(xs, R[, condition], col=kEquilCondColors[condition])
    })

  AddLinearAxis(1, tick.space=1, label.space=1,
                label="Length of complementarity")
  AddLinearAxis(2, tick.space=y.tick.space, label.space=y.label.space,
                label="Average enrichment within positions 13-31")
  # for (condition in conditions) {
  #   lines(kmer_len_start:kmer_len_stop, Rs[condition, ], type="o",
  #         col=kEquilCondColors[condition])
  # }
  conditions <- conditions[-1]
  xy <- GetPlotFractionalCoords(0.05, 0.95)
  legend(x=xy[1], y=xy[2], bty="n", legend=conditions, pch=19,
         col=kEquilCondColors[conditions])
  xy <- GetPlotFractionalCoords(0.95, 0.95)
  text(xy[1], xy[2], labels=mirna, adj=1)
  xy <- GetPlotFractionalCoords(0.95, 0.90)
  text(xy[1], xy[2], labels=experiment, adj=1)
  xy <- GetPlotFractionalCoords(0.95, 0.85)
  if (geomean) string <- "Geometric mean"
  else         string <- "Arithmetic mean"
  text(xy[1], xy[2], labels=string, adj=1)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

GetSiteSeqKds <- function(mirna, xpos=20, ypos=20, height=5, width=8,
                          plotting=TRUE, adjust=TRUE, pdf.plot=FALSE) {
  data <- SubfunctionCall(LoadBeckerEtAlData)

  # Get the fully complementary sequence.
  seq_perfect <- SubfunctionCall(GetSiteSeq, "21mer-m1.21")
  # Make the single-mismatch regexpression
  vec_1mm_strs <- sapply(1:nchar(seq_perfect), function(i) {
    sprintf("%s[^%s]%s", substr(seq_perfect, 0, i - 1),
    substr(seq_perfect, i, i), substr(seq_perfect, i + 1, nchar(seq_perfect)))  
  })
  # Make the double-mismatch regexpressions
  vec_2mm_strs <- unlist(sapply(1:(nchar(seq_perfect) - 1), function(i) {
    sapply((i + 1):(nchar(seq_perfect)), function(j) {
      sprintf("%s[^%s]%s[^%s]%s", substr(seq_perfect, 0, i - 1),
              substr(seq_perfect, i, i), substr(seq_perfect, i + 1, j - 1),
              substr(seq_perfect, j, j), substr(seq_perfect, j + 1,
                                                nchar(seq_perfect)))
    })
  }))
  # Make the triple-mismatch regexpressions
  vec_3mm_strs <- unlist(sapply(1:(nchar(seq_perfect) - 2), function(i) {
    sapply((i + 1):(nchar(seq_perfect) - 1), function(j) {
      sapply((j + 1):nchar(seq_perfect), function(k) {
        sprintf("%s[^%s]%s[^%s]%s[^%s]%s", substr(seq_perfect, 0, i - 1),
                substr(seq_perfect, i, i), substr(seq_perfect, i + 1, j - 1),
                substr(seq_perfect, j, j), substr(seq_perfect, j + 1, k - 1),
                substr(seq_perfect, k, k), substr(seq_perfect, k + 1,
                                                  nchar(seq_perfect)))
      })
    })
  }))
  # Make the single deletion regexpressions
  vec_1d <- sapply(1:nchar(seq_perfect), function(i) {
    sprintf("%s%s", substr(seq_perfect, 0, i - 1),
            substr(seq_perfect, i + 1, nchar(seq_perfect)))
  })
  # Make the double deletion regexpressions.
  vec_2d <- unlist(sapply(1:(nchar(seq_perfect) - 1), function(i) {
    sapply((i + 1):nchar(seq_perfect), function(j) {
      sprintf("%s%s%s", substr(seq_perfect, 0, i - 1),
              substr(seq_perfect, i + 1, j - 1),
              substr(seq_perfect, j + 1, nchar(seq_perfect)))

    })
  }))
  # make the 1-7 nt bulges.
  vec_1.7b <- as.vector(sapply(1:(nchar(seq_perfect) - 1), function(i) {
    sapply(c("A", "C", "G", "T"), function(nuc) {
        sprintf("%s%s{1,7}%s", substr(seq_perfect, 0, i), nuc,
                substr(seq_perfect, i + 1, nchar(seq_perfect)))
    })
  }))
  # make the combinatorial dels and bulges.
  vec_1d1b <- as.vector(sapply(vec_1d, function(sequence) {
    sapply(1:nchar(sequence) - 1, function(i) {
      sapply(c("A", "C", "G", "T"), function(nuc) {
          sprintf("%s%s%s", substr(sequence, 0, i), nuc,
                  substr(sequence, i + 1, nchar(sequence)))
      })
    })
  }))
  # Remove the sequences that re-establish perfect complementarity
  vec_1d1b <- unique(grep(seq_perfect, vec_1d1b, value=TRUE, invert=TRUE))
  # Make three-prime regular expressions. Removes all sequences with 6-nt of
  # pairing to the three-prime end.
  vec_3p <- sapply(1:8, function(start) {
    substr(seq_perfect, start=start, stop=start + 5)
  })
  # Make the more stringent 11mer site regexpressions.
  vec_13_2 <- sapply(seq(15, 21), function(stop) {
    message(stop)
    n <- nchar(seq_perfect)
    seq_end <- substr(seq_perfect, 1, n - stop)
    seq_mm <- substr(seq_perfect, n - stop + 1, n - 13)
    seq_seed <- substr(seq_perfect, n - 13 + 1, n)
    sprintf("%s%s%s", seq_end, SeqComplement(seq_mm), seq_seed)
  })
  # Make the more stringent 11mer site regexpressions.
  vec_12_2 <- sapply(seq(15, 21), function(stop) {
    message(stop)
    n <- nchar(seq_perfect)
    seq_end <- substr(seq_perfect, 1, n - stop)
    seq_mm <- substr(seq_perfect, n - stop + 1, n - 12)
    seq_seed <- substr(seq_perfect, n - 12 + 1, n)
    sprintf("%s%s%s", seq_end, SeqComplement(seq_mm), seq_seed)
  })
  # Make the more stringent 11mer site regexpressions.
  vec_11_2 <- sapply(seq(15, 21), function(stop) {
    message(stop)
    n <- nchar(seq_perfect)
    seq_end <- substr(seq_perfect, 1, n - stop)
    seq_mm <- substr(seq_perfect, n - stop + 1, n - 11)
    seq_seed <- substr(seq_perfect, n - 11 + 1, n)
    sprintf("%s%s%s", seq_end, SeqComplement(seq_mm), seq_seed)
  })
  #Make the more stringent 10mer site regexpressions.
  vec_10_2 <- sapply(seq(15, 21), function(stop) {
    message(stop)
    n <- nchar(seq_perfect)
    seq_end <- substr(seq_perfect, 1, n - stop)
    seq_mm <- substr(seq_perfect, n - stop + 1, n - 10)
    seq_seed <- substr(seq_perfect, n - 10 + 1, n)
    sprintf("%s%s%s", seq_end, SeqComplement(seq_mm), seq_seed)
  })
  #Make the more stringent 10mer site regexpressions.
  vec_9_2 <- sapply(seq(15, 21), function(stop) {
    message(stop)
    n <- nchar(seq_perfect)
    seq_end <- substr(seq_perfect, 1, n - stop)
    seq_mm <- substr(seq_perfect, n - stop + 1, n - 9)
    seq_seed <- substr(seq_perfect, n - 9 + 1, n)
    sprintf("%s%s%s", seq_end, SeqComplement(seq_mm), seq_seed)
  })
  # Make the more stringent 8mer site regexpressions.
  vec_8_2 <- sapply(seq(15, 21), function(stop) {
    n <- nchar(seq_perfect)
    seq_end <- substr(seq_perfect, 1, n - stop)
    seq_mm <- substr(seq_perfect, n - stop + 1, n - 8)
    seq_seed <- substr(seq_perfect, n - 8 + 1, n)
    sprintf("%s%s%s", seq_end, SeqComplement(seq_mm), seq_seed)
  })
  # Make the more stringent 7mer-A1 site regexpressions.
  vec_7m8_2 <- sapply(seq(15, 21), function(stop) {
    n <- nchar(seq_perfect)
    seq_end <- substr(seq_perfect, 1, n - stop)
    seq_mm <- substr(seq_perfect, n - stop + 1, n - 8)
    seq_seed <- substr(seq_perfect, n - 8 + 1, n - 1)
    sprintf("%s%s%s", seq_end, SeqComplement(seq_mm), seq_seed)
  })
  # Make the more stringent 7mer-A1 site regexpressions.
  vec_7A1_2 <- sapply(seq(15, 21), function(stop) {
    n <- nchar(seq_perfect)
    seq_end <- substr(seq_perfect, 1, n - stop)
    seq_mm <- substr(seq_perfect, n - stop + 1, n - 7)
    seq_seed <- substr(seq_perfect, n - 7 + 1, n)
    sprintf("%s%s%s", seq_end, SeqComplement(seq_mm), seq_seed)
  })
  # Make the more stringent 7mer-A1 site regexpressions.
  vec_6_2 <- sapply(seq(15, 21), function(stop) {
    n <- nchar(seq_perfect)
    seq_end <- substr(seq_perfect, 1, n - stop)
    seq_mm <- substr(seq_perfect, n - stop + 1, n - 7)
    seq_seed <- substr(seq_perfect, n - 7 + 1, n - 1)
    sprintf("%s%s%s", seq_end, SeqComplement(seq_mm), seq_seed)
  })
  # Make the more stringent 6mer-A1 site regexpressions.
  vec_6A1_2 <- sapply(seq(15, 21), function(stop) {
    message(stop)
    n <- nchar(seq_perfect)
    seq_end <- substr(seq_perfect, 1, n - stop)
    seq_mm <- substr(seq_perfect, n - stop + 1, n - 6)
    seq_seed <- substr(seq_perfect, n - 6 + 1, n)
    sprintf("%s%s%s", seq_end, SeqComplement(seq_mm), seq_seed)
  })
  # Make the more stringent 6mer-m8 site regexpressions.
  vec_6m8_2 <- sapply(seq(15, 21), function(stop) {
    n <- nchar(seq_perfect)
    seq_end <- substr(seq_perfect, 1, n - stop)
    seq_mm <- substr(seq_perfect, n - stop + 1, n - 8)
    seq_seed <- substr(seq_perfect, n - 8 + 1, n - 2)
    sprintf("%s%s%s", seq_end, SeqComplement(seq_mm), seq_seed)
  })


  regex_1mm <- paste0(vec_1mm_strs, collapse="|")
  regex_2mm <- paste0(vec_2mm_strs, collapse="|")
  regex_3mm <- paste0(vec_3mm_strs, collapse="|")

  regex_1d <- paste0(vec_1d, collapse="|")
  regex_2d <- paste0(vec_2d, collapse="|")
  regex_1.7b <- paste0(vec_1.7b, collapse="|")
  regex_1d1b <- paste0(vec_1d1b, collapse="|")
  regex_3p <- paste0(vec_3p, collapse="|")
  regex_13_2 <- paste0(vec_13_2, collapse="|")
  regex_12_2 <- paste0(vec_12_2, collapse="|")
  regex_11_2 <- paste0(vec_11_2, collapse="|")
  regex_10_2 <- paste0(vec_10_2, collapse="|")
  regex_9_2 <- paste0(vec_9_2, collapse="|")
  regex_8_2 <- paste0(vec_8_2, collapse="|")
  regex_7m8_2 <- paste0(vec_7m8_2, collapse="|")
  regex_7A1_2 <- paste0(vec_7A1_2, collapse="|")
  regex_6_2 <- paste0(vec_6_2, collapse="|")
  regex_6m8_2 <- paste0(vec_6m8_2, collapse="|")
  regex_6A1_2 <- paste0(vec_6A1_2, collapse="|")

  regex_13 <- GetSiteSeq(mirna, "13mer-m1.13")
  regex_12 <- GetSiteSeq(mirna, "12mer-m1.12")
  regex_11 <- GetSiteSeq(mirna, "11mer-m1.11")
  regex_10 <- GetSiteSeq(mirna, "10mer-m1.10")
  regex_9 <- GetSiteSeq(mirna, "9mer-m1.9")
  regex_8 <- GetSiteSeq(mirna, "8mer")
  regex_7m8 <- GetSiteSeq(mirna, "7mer-m8")
  regex_7A1 <- GetSiteSeq(mirna, "7mer-A1")
  regex_6 <- GetSiteSeq(mirna, "6mer")
  regex_6m8 <- GetSiteSeq(mirna, "6mer-m8")
  regex_6A1 <- GetSiteSeq(mirna, "6mer-A1")


  list_regex <- list(`Perfect`=seq_perfect, `mm1`=regex_1mm, `mm2`=regex_2mm,
                     `mm3`=regex_3mm, `1d`=regex_1d, `2d`=regex_2d,
                     `1.7b`=regex_1.7b, `1d1b`=regex_1d1b, `6mer-3p`=regex_3p,
                     `13mer-m1.13_alt`=regex_13_2, `12mer-m1.12_alt`=regex_12_2,
                     `11mer-m1.11_alt`=regex_11_2, `10mer-m1.10_alt`=regex_10_2,
                     `9mer-m1.9_alt`=regex_9_2, `8mer_alt`=regex_8_2,
                     `7mer-m8_alt`=regex_7m8_2, `7mer-A1_alt`=regex_7A1_2,
                     `6mer_alt`=regex_6_2, `6mer-m8_alt`=regex_6m8_2,
                     `6mer-A1_alt`=regex_6A1_2, `13mer-m1.13`=regex_13,
                     `12mer-m1.12`=regex_12, `11mer-m1.11`=regex_11,
                     `10mer-m1.10`=regex_10, `9mer-m1.9`=regex_9,
                     `8mer`=regex_8, `7mer-m8`=regex_7m8, `7mer-A1`=regex_7A1,
                     `6mer`=regex_6, `6mer-m8`=regex_6m8, `6mer-A1`=regex_6A1)

  cols <- c(`Perfect`="black", `mm1`="gray50", `mm2`="gray70", `mm3`="gray90", 
            `1d`="orange", `2d`="orangered", `1.7b`="goldenrod",
            `1d1b`="brown", `6mer-3p`="forestgreen")
  cols <- c(cols, `13mer-m1.13_alt`="darkgreen", `12mer-m1.12_alt`="green4",
            `11mer-m1.11_alt`="seagreen3", `10mer-m1.10_alt`="lawngreen",
            `9mer-m1.9_alt`="greenyellow", `8mer_alt`=kSiteColors["8mer"],
            `7mer-m8_alt`=kSiteColors["7mer-m8"],
            `7mer-A1_alt`=kSiteColors["7mer-A1"],
            `6mer_alt`=kSiteColors["6mer"],
            `6mer-m8_alt`=kSiteColors["6mer-m8"],
            `6mer-A1_alt`=kSiteColors["6mer-A1"],
            `13mer-m1.13`="darkgreen", `12mer-m1.12`="green4",
            `11mer-m1.11`="seagreen3", `10mer-m1.10`="lawngreen",
            `9mer-m1.9`="greenyellow",
            kSiteColors[kSeedSites])

  names(cols)[10:20] <- c("13mer-m1.13_alt", "12mer-m1.12_alt",
                          "11mer-m1.11_alt", "10mer-m1.10_alt", "9mer-m1.9_alt",
                          "8mer_alt", "7mer-m8_alt", "7mer-A1_alt",
                          "6mer_alt", "6mer-m8_alt", "6mer-A1_alt")
  list_data <- lapply(names(list_regex), function(name_i) {
    regex_i <- list_regex[[name_i]]
    data_i <- data[grep(regex_i, rownames(data), perl=TRUE), ]
    data <<- data[grep(regex_i, rownames(data), perl=TRUE, invert=TRUE), ]
    data_i  
  })

  names(list_data) <- names(list_regex)
  list_data_global <<- list_data
  if (plotting) {
    SubfunctionCall(GenericFigureSaveFile)
    xmin <- 10
    xmax <- 10000
    ymin <- 0
    ymax <- 1
    
    BlankPlot(log="x", adjusted=TRUE)

    AddLogAxis(1, label="Kd (pM)", adj=TRUE)
    AddLinearAxis(2, tick.space=0.05, label.space=0.2, percent=TRUE,
                  label="Cumulative fraction (%)")
    x_seq <- 10^seq(1, 4, length.out=100)
    ltys <- as.integer(grepl("_2", names(list_data))) + 1
    names(ltys) <- names(list_data)
    lapply(names(list_data), function(name_i) {
      print(name_i)
      data_i <- list_data[[name_i]]
      print(dim(data_i))
      ECDF_func <- ecdf(data_i[, 1])
      if (grepl("_2", list_data)) {
        lty <- 2
      } else {
        lty <- 1
      }
      lines(x_seq, ECDF_func(x_seq), col=cols[name_i], lty=ltys[name_i])
      message(sprintf("Median range: %s - %s pM", min(data_i[, 1], na.rm=TRUE), max(data_i[, 1], na.rm=TRUE)))
      message(sprintf("Min range: %s-%s pM", min(data_i[, 2], na.rm=TRUE), max(data_i[, 2], na.rm=TRUE)))
      message(sprintf("Max range: %s-%s pM", min(data_i[, 3], na.rm=TRUE), max(data_i[, 3], na.rm=TRUE)))
      print(max(data_i[, 1]))
    })
    xy <- GetPlotFractionalCoords(1.05, 0.95, log='x')
    text(xy[1], xy[2], adj=c(0, 0), label=mirna, xpd=NA)
    xy <- GetPlotFractionalCoords(1.05, 0.95, log='x')
    legend(xy[1], xy[2], legend=names(list_data), col=cols, bty="n", lwd=1,
           y.intersp=0.9, xpd=NA, lty=ltys[names(list_data)])
    if (class(pdf.plot) == "character") {
      dev.off()
    }
  } else {
    list_data
  }
}


# ComparePositions9Through13 <- function(mirna, alt=FALSE, xpos=20, ypos=20,
#                                        height=5, width=8, plotting=TRUE,
#                                        adjust=TRUE, pdf.plot=FALSE)
# {
#   data <- SubfunctionCall(LoadBeckerEtAlData)
#   # Get the fully complementary sequence.
#   seq_perfect <- SubfunctionCall(GetSiteSeq, "21mer-m1.21")
#   # Make the single-mismatch regexpression
#   vec_1mm_strs <- sapply(1:nchar(seq_perfect), function(i) {
#     sprintf("%s[^%s]%s", substr(seq_perfect, 0, i - 1),
#     substr(seq_perfect, i, i), substr(seq_perfect, i + 1, nchar(seq_perfect)))  
#   })
#   # Make the double-mismatch regexpressions
#   vec_2mm_strs <- unlist(sapply(1:(nchar(seq_perfect) - 1), function(i) {
#     sapply((i + 1):(nchar(seq_perfect)), function(j) {
#       sprintf("%s[^%s]%s[^%s]%s", substr(seq_perfect, 0, i - 1),
#               substr(seq_perfect, i, i), substr(seq_perfect, i + 1, j - 1),
#               substr(seq_perfect, j, j), substr(seq_perfect, j + 1,
#                                                 nchar(seq_perfect)))
#     })
#   }))
#   # Make the triple-mismatch regexpressions
#   vec_3mm_strs <- unlist(sapply(1:(nchar(seq_perfect) - 2), function(i) {
#     sapply((i + 1):(nchar(seq_perfect) - 1), function(j) {
#       sapply((j + 1):nchar(seq_perfect), function(k) {
#         sprintf("%s[^%s]%s[^%s]%s[^%s]%s", substr(seq_perfect, 0, i - 1),
#                 substr(seq_perfect, i, i), substr(seq_perfect, i + 1, j - 1),
#                 substr(seq_perfect, j, j), substr(seq_perfect, j + 1, k - 1),
#                 substr(seq_perfect, k, k), substr(seq_perfect, k + 1,
#                                                   nchar(seq_perfect)))
#       })
#     })
#   }))
#   # Make the single deletion regexpressions
#   vec_1d <- sapply(1:nchar(seq_perfect), function(i) {
#     sprintf("%s%s", substr(seq_perfect, 0, i - 1),
#             substr(seq_perfect, i + 1, nchar(seq_perfect)))
#   })
#   # Make the double deletion regexpressions.
#   vec_2d <- unlist(sapply(1:(nchar(seq_perfect) - 1), function(i) {
#     sapply((i + 1):nchar(seq_perfect), function(j) {
#       sprintf("%s%s%s", substr(seq_perfect, 0, i - 1),
#               substr(seq_perfect, i + 1, j - 1),
#               substr(seq_perfect, j + 1, nchar(seq_perfect)))

#     })
#   }))
#   # make the 1-7 nt bulges.
#   vec_1.7b <- as.vector(sapply(1:(nchar(seq_perfect) - 1), function(i) {
#     sapply(c("A", "C", "G", "T"), function(nuc) {
#         sprintf("%s%s{1,7}%s", substr(seq_perfect, 0, i), nuc,
#                 substr(seq_perfect, i + 1, nchar(seq_perfect)))
#     })
#   }))
#   # make the combinatorial dels and bulges.
#   vec_1d1b <- as.vector(sapply(vec_1d, function(sequence) {
#     sapply(1:nchar(sequence) - 1, function(i) {
#       sapply(c("A", "C", "G", "T"), function(nuc) {
#           sprintf("%s%s%s", substr(sequence, 0, i), nuc,
#                   substr(sequence, i + 1, nchar(sequence)))
#       })
#     })
#   }))
#   # Remove the sequences that re-establish perfect complementarity
#   vec_1d1b <- unique(grep(seq_perfect, vec_1d1b, value=TRUE, invert=TRUE))
#   # Make three-prime regular expressions. Removes all sequences with 6-nt of
#   # pairing to the three-prime end.
#   vec_3p <- sapply(1:8, function(start) {
#     substr(seq_perfect, start=start, stop=start + 5)
#   })
#   # Make the more stringent 11mer site regexpressions.
#   vec_13_2 <- sapply(seq(15, 21), function(stop) {
#     message(stop)
#     n <- nchar(seq_perfect)
#     seq_end <- substr(seq_perfect, 1, n - stop)
#     seq_mm <- substr(seq_perfect, n - stop + 1, n - 13)
#     seq_seed <- substr(seq_perfect, n - 13 + 1, n)
#     sprintf("%s%s%s", seq_end, SeqComplement(seq_mm), seq_seed)
#   })
#   # Make the more stringent 11mer site regexpressions.
#   vec_12_2 <- sapply(seq(15, 21), function(stop) {
#     message(stop)
#     n <- nchar(seq_perfect)
#     seq_end <- substr(seq_perfect, 1, n - stop)
#     seq_mm <- substr(seq_perfect, n - stop + 1, n - 12)
#     seq_seed <- substr(seq_perfect, n - 12 + 1, n)
#     sprintf("%s%s%s", seq_end, SeqComplement(seq_mm), seq_seed)
#   })
#   # Make the more stringent 11mer site regexpressions.
#   vec_11_2 <- sapply(seq(15, 21), function(stop) {
#     message(stop)
#     n <- nchar(seq_perfect)
#     seq_end <- substr(seq_perfect, 1, n - stop)
#     seq_mm <- substr(seq_perfect, n - stop + 1, n - 11)
#     seq_seed <- substr(seq_perfect, n - 11 + 1, n)
#     sprintf("%s%s%s", seq_end, SeqComplement(seq_mm), seq_seed)
#   })
#   #Make the more stringent 10mer site regexpressions.
#   vec_10_2 <- sapply(seq(15, 21), function(stop) {
#     message(stop)
#     n <- nchar(seq_perfect)
#     seq_end <- substr(seq_perfect, 1, n - stop)
#     seq_mm <- substr(seq_perfect, n - stop + 1, n - 10)
#     seq_seed <- substr(seq_perfect, n - 10 + 1, n)
#     sprintf("%s%s%s", seq_end, SeqComplement(seq_mm), seq_seed)
#   })
#   #Make the more stringent 10mer site regexpressions.
#   vec_9_2 <- sapply(seq(15, 21), function(stop) {
#     message(stop)
#     n <- nchar(seq_perfect)
#     seq_end <- substr(seq_perfect, 1, n - stop)
#     seq_mm <- substr(seq_perfect, n - stop + 1, n - 9)
#     seq_seed <- substr(seq_perfect, n - 9 + 1, n)
#     sprintf("%s%s%s", seq_end, SeqComplement(seq_mm), seq_seed)
#   })
#   # Make the more stringent 8mer site regexpressions.
#   vec_8_2 <- sapply(seq(15, 21), function(stop) {
#     n <- nchar(seq_perfect)
#     seq_end <- substr(seq_perfect, 1, n - stop)
#     seq_mm <- substr(seq_perfect, n - stop + 1, n - 8)
#     seq_seed <- substr(seq_perfect, n - 8 + 1, n)
#     sprintf("%s%s%s", seq_end, SeqComplement(seq_mm), seq_seed)
#   })
#   # Make the more stringent 7mer-A1 site regexpressions.
#   vec_7m8_2 <- sapply(seq(15, 21), function(stop) {
#     n <- nchar(seq_perfect)
#     seq_end <- substr(seq_perfect, 1, n - stop)
#     seq_mm <- substr(seq_perfect, n - stop + 1, n - 8)
#     seq_seed <- substr(seq_perfect, n - 8 + 1, n - 1)
#     sprintf("%s%s%s", seq_end, SeqComplement(seq_mm), seq_seed)
#   })
#   # Make the more stringent 7mer-A1 site regexpressions.
#   vec_7A1_2 <- sapply(seq(15, 21), function(stop) {
#     n <- nchar(seq_perfect)
#     seq_end <- substr(seq_perfect, 1, n - stop)
#     seq_mm <- substr(seq_perfect, n - stop + 1, n - 7)
#     seq_seed <- substr(seq_perfect, n - 7 + 1, n)
#     sprintf("%s%s%s", seq_end, SeqComplement(seq_mm), seq_seed)
#   })
#   # Make the more stringent 7mer-A1 site regexpressions.
#   vec_6_2 <- sapply(seq(15, 21), function(stop) {
#     n <- nchar(seq_perfect)
#     seq_end <- substr(seq_perfect, 1, n - stop)
#     seq_mm <- substr(seq_perfect, n - stop + 1, n - 7)
#     seq_seed <- substr(seq_perfect, n - 7 + 1, n - 1)
#     sprintf("%s%s%s", seq_end, SeqComplement(seq_mm), seq_seed)
#   })
#   # Make the more stringent 6mer-A1 site regexpressions.
#   vec_6A1_2 <- sapply(seq(15, 21), function(stop) {
#     message(stop)
#     n <- nchar(seq_perfect)
#     seq_end <- substr(seq_perfect, 1, n - stop)
#     seq_mm <- substr(seq_perfect, n - stop + 1, n - 6)
#     seq_seed <- substr(seq_perfect, n - 6 + 1, n)
#     sprintf("%s%s%s", seq_end, SeqComplement(seq_mm), seq_seed)
#   })
#   # Make the more stringent 6mer-m8 site regexpressions.
#   vec_6m8_2 <- sapply(seq(15, 21), function(stop) {
#     n <- nchar(seq_perfect)
#     seq_end <- substr(seq_perfect, 1, n - stop)
#     seq_mm <- substr(seq_perfect, n - stop + 1, n - 8)
#     seq_seed <- substr(seq_perfect, n - 8 + 1, n - 2)
#     sprintf("%s%s%s", seq_end, SeqComplement(seq_mm), seq_seed)
#   })


#   regex_1mm <- paste0(vec_1mm_strs, collapse="|")
#   regex_2mm <- paste0(vec_2mm_strs, collapse="|")
#   regex_3mm <- paste0(vec_3mm_strs, collapse="|")

#   regex_1d <- paste0(vec_1d, collapse="|")
#   regex_2d <- paste0(vec_2d, collapse="|")
#   regex_1.7b <- paste0(vec_1.7b, collapse="|")
#   regex_1d1b <- paste0(vec_1d1b, collapse="|")
#   regex_3p <- paste0(vec_3p, collapse="|")
#   regex_13_2 <- paste0(vec_13_2, collapse="|")
#   regex_12_2 <- paste0(vec_12_2, collapse="|")
#   regex_11_2 <- paste0(vec_11_2, collapse="|")
#   regex_10_2 <- paste0(vec_10_2, collapse="|")
#   regex_9_2 <- paste0(vec_9_2, collapse="|")
#   regex_8_2 <- paste0(vec_8_2, collapse="|")
#   regex_7m8_2 <- paste0(vec_7m8_2, collapse="|")
#   regex_7A1_2 <- paste0(vec_7A1_2, collapse="|")
#   regex_6_2 <- paste0(vec_6_2, collapse="|")
#   regex_6A1_2 <- paste0(vec_6A1_2, collapse="|")
#   regex_6m8_2 <- paste0(vec_6m8_2, collapse="|")

#   regex_13 <- GetSiteSeq(mirna, "13mer-m1.13")
#   regex_12 <- GetSiteSeq(mirna, "12mer-m1.12")
#   regex_11 <- GetSiteSeq(mirna, "11mer-m1.11")
#   regex_10 <- GetSiteSeq(mirna, "10mer-m1.10")
#   regex_9 <- GetSiteSeq(mirna, "9mer-m1.9")
#   regex_8 <- GetSiteSeq(mirna, "8mer")
#   regex_7m8 <- GetSiteSeq(mirna, "7mer-m8")
#   regex_7A1 <- GetSiteSeq(mirna, "7mer-A1")
#   regex_6 <- GetSiteSeq(mirna, "6mer")
#   regex_6m8 <- GetSiteSeq(mirna, "6mer-m8")
#   regex_6A1 <- GetSiteSeq(mirna, "6mer-A1")


#   list_regex <- list(`Perfect`=seq_perfect, `mm1`=regex_1mm, `mm2`=regex_2mm,
#                      `mm3`=regex_3mm, `1d`=regex_1d, `2d`=regex_2d,
#                      `1.7b`=regex_1.7b, `1d1b`=regex_1d1b, `6mer-3p`=regex_3p,
#                      `13mer-m1.13_alt`=regex_13_2, `12mer-m1.12_alt`=regex_12_2,
#                      `11mer-m1.11_alt`=regex_11_2, `10mer-m1.10_alt`=regex_10_2,
#                      `9mer-m1.9_alt`=regex_9_2, `8mer_alt`=regex_8_2,
#                      `7mer-m8_alt`=regex_7m8_2, `7mer-A1_alt`=regex_7A1_2,
#                      `6mer_alt`=regex_6_2, `6mer-m8_alt`=regex_6m8_2,
#                      `6mer-A1_alt`=regex_6A1_2, `13mer-m1.13`=regex_13,
#                      `12mer-m1.12`=regex_12, `11mer-m1.11`=regex_11,
#                      `10mer-m1.10`=regex_10, `9mer-m1.9`=regex_9,
#                      `8mer`=regex_8, `7mer-m8`=regex_7m8, `7mer-A1`=regex_7A1,
#                      `6mer`=regex_6, `6mer-m8`=regex_6m8, `6mer-A1`=regex_6A1)

#   cols <- c(`Perfect`="black", `mm1`="gray50", `mm2`="gray70", `mm3`="gray90", 
#             `1d`="orange", `2d`="orangered", `1.7b`="goldenrod",
#             `1d1b`="brown", `6mer-3p`="forestgreen")
#   cols <- c(cols, `13mer-m1.13_alt`="darkgreen", `12mer-m1.12_alt`="green4",
#             `11mer-m1.11_alt`="seagreen3", `10mer-m1.10_alt`="lawngreen",
#             `9mer-m1.9_alt`="greenyellow", `8mer_alt`=kSiteColors["8mer"],
#             `7mer-m8_alt`=kSiteColors["7mer-m8"],
#             `7mer-A1_alt`=kSiteColors["7mer-A1"],
#             `6mer_alt`=kSiteColors["6mer"],
#             `6mer-m8_alt`=kSiteColors["6mer-m8"],
#             `6mer-A1_alt`=kSiteColors["6mer-A1"],
#             `13mer-m1.13`="darkgreen", `12mer-m1.12`="green4",
#             `11mer-m1.11`="seagreen3", `10mer-m1.10`="lawngreen",
#             `9mer-m1.9`="greenyellow",
#             kSiteColors[kSeedSites])

#   names(cols)[10:20] <- c("13mer-m1.13_alt", "12mer-m1.12_alt",
#                           "11mer-m1.11_alt", "10mer-m1.10_alt", "9mer-m1.9_alt",
#                           "8mer_alt", "7mer-m8_alt", "7mer-A1_alt",
#                           "6mer_alt", "6mer-m8_alt", "6mer-A1_alt")
#   list_data <- lapply(names(list_regex), function(name_i) {
#     regex_i <- list_regex[[name_i]]
#     data_i <- data[grep(regex_i, rownames(data), perl=TRUE), ]
#     data <<- data[grep(regex_i, rownames(data), perl=TRUE, invert=TRUE), ]
#     data_i  
#   })
#   names(list_data) <- names(list_regex)
#   list_data_global <<- list_data
#   if (alt) {
#     inds <- 10:20
#   } else {
#     inds <- 21:30
#   }
#   cols <- cols[inds]

#   if (plotting) {
#     SubfunctionCall(GenericFigureSaveFile)
#     xmin <- 10
#     xmax <- 10000
#     ymin <- 0
#     ymax <- 1
    
#     BlankPlot(log="x", adjusted=TRUE)

#     AddLogAxis(1, label="Kd (pM)", adj=TRUE)
#     AddLinearAxis(2, tick.space=0.05, label.space=0.2, percent=TRUE,
#                   label="Cumulative fraction (%)")
#     x_seq <- 10^seq(1, 4, length.out=100)
#     ltys <- as.integer(grepl("_2", names(list_data))) + 1
#     names(ltys) <- names(list_data)
#     ltys <- ltys[inds]
#     lapply(names(list_data)[inds], function(name_i) {
#       print(name_i)
#       data_i <- list_data[[name_i]]
#       print(dim(data_i))
#       ECDF_func <- ecdf(data_i[, 1])
#       if (grepl("_2", list_data)) {
#         lty <- 2
#       } else {
#         lty <- 1
#       }
#       print(ECDF_func(x_seq))
#       print(cols[name_i])
#       lines(x_seq, ECDF_func(x_seq), col=cols[name_i], lty=ltys[name_i])
#       message(sprintf("Median range: %s - %s pM", min(data_i[, 1], na.rm=TRUE), max(data_i[, 1], na.rm=TRUE)))
#       message(sprintf("Min range: %s-%s pM", min(data_i[, 2], na.rm=TRUE), max(data_i[, 2], na.rm=TRUE)))
#       message(sprintf("Max range: %s-%s pM", min(data_i[, 3], na.rm=TRUE), max(data_i[, 3], na.rm=TRUE)))
#       print(max(data_i[, 1]))
#     })
#     xy <- GetPlotFractionalCoords(1.05, 0.95, log='x')
#     text(xy[1], xy[2], adj=c(0, 0), label=mirna, xpd=NA)
#     xy <- GetPlotFractionalCoords(1.05, 0.95, log='x')
#     legend(xy[1], xy[2], legend=names(list_data)[inds], col=cols, bty="n", lwd=1,
#            y.intersp=0.9, xpd=NA, lty=ltys[names(list_data)[inds]])
#     if (class(pdf.plot) == "character") {
#       dev.off()
#     }
#   }
#   list_data
# }





PlotComparison <- function(mirna, xpos=20, ypos=20, height=5, width=5,
                           pdf.plot=FALSE)
{
  kd_Becker <- GetSiteSeqKds(mirna, plotting=FALSE)
  kd_median <- lapply(kd_Becker, function(data) {
    median(data[, 1], na.rm=TRUE)
  })
  if (mirna == "miR-21") {
    mirna_rbns <- "miR-1"
    buffer <- TRUE
    xlab <- "miR-1 "
    ylab <- "miR-21 "
    text_append <- " vs. miR-1"
  } else {
    mirna_rbns <- mirna
    buffer <- FALSE
    xlab <- ""
    ylab <- ""
    text_append <- ""
  }
  kds_RBNS <- SubfunctionCall(EquilPars, mirna=mirna_rbns,
                              experiment="equilibrium")

  SubfunctionCall(FigureSaveFile)
  ymin <- 10
  ymax <- 10000
  xmin <- 1e-4
  xmax <- 1
  
  BlankPlot(log="xy", inv='xy')
  AddLogAxis(1, label=sprintf("%s Relative Kd; McGeary, Lin, et al.", xlab))
  AddLogAxis(2, label=sprintf("%s Kd; Becker et al. (pM)", ylab))
  xy <- GetPlotFractionalCoords(0.95, 0.05, log="xy", inv="xy")
  text(xy[1], xy[2], adj=c(1, 0), label=sprintf("%s%s", mirna, text_append))
  points(kds_RBNS[rep(paste0(kSeedSites, "_Kd"), 2), 1],
         kd_median[c(kSeedSites, sprintf("%s_alt", kSeedSites))],
         col=kSiteColors[rep(kSeedSites, 2)], pch=rep(c(19, 1), each=6),
         cex=rep(c(1, 2), each=6))

  if (mirna == "let-7a") {
    Becker_et_al_delG <- c(`8mer`=-13.70166171, `7mer-m8`=-12.89980191,
                          `7mer-A1`=-12.65070981, `6mer`=-12.10140861)
    R <- 1.987e-3 # in kcal K-1 mol-1
    T <- 310.15 # in K
    Becker_published_Kds_let7a <- exp(Becker_et_al_delG/(R*T))*1e12
    points(kds_RBNS[paste0(kCanonicalSites, "_Kd"), 1],
           Becker_published_Kds_let7a,
           col=kSiteColors[kCanonicalSites], pch=17)
  }
  xy <- GetPlotFractionalCoords(0.05, 0.95, log='xy', inv='xy')
  if (mirna == "let-7a") {
    legend(xy[1], xy[2], legend=c("Median from reanalysis",
                                  "Median of cumulative mismatch targets",
                                  "Published value"), pt.cex=c(1, 2, 1),
           pch=c(19, 1, 17), bty="n")    
  } else {
    legend(xy[1], xy[2], legend=c("Median from reanalysis",
                                  "Median of cumulative mismatch targets"),
           pt.cex=c(1, 2), pch=c(19, 1), bty="n")       
  }
  xy <- GetPlotFractionalCoords(0.05, 0.7, log='xy', inv='xy')
  legend(xy[1], xy[2], legend=kSeedSites, pch=19,
         col=kSiteColors[kSeedSites], bty="n")

  if (class(pdf.plot) == "character") {
      dev.off()    
  }
}





GetPerfectTargetMismatchStretch <- function(mirna, mm_start, mm_stop) {
  # Get the starting string containing the perfectly complementary target site.
  perfect_str <- SubfunctionCall(GetSiteSeq, "21mer-m1.21")
  n <- 21
  # Convert this string to a vector of nucleotide letters.
  perfect_vec <- unlist(strsplit(perfect_str, split=""))
  # Loop over the list of contiguous postions to be changed to their complment
  # (i.e., changed back to the same nucleotide identity as that of the miRNA).
  for (pos in (21 - mm_start:mm_stop + 1)) {
    # kComplements in Lists.R
    perfect_vec[pos] <- kComplements[perfect_vec[pos]]
  }
  # Return the string formed from collapsing the vector of nucleotide strings.
  paste0(perfect_vec, collapse="")
}


GetMismatchStretchVariantBecker <- function(mirna) {
  # Load the data table
  data <- LoadBeckerEtAlData(mirna)
  # Make the vector of variant names.
  variants <- rownames(data)
  # Make the output table.
  output_matrix <- matrix(NaN, nrow=21, ncol=21)
  lowCI_matrix <- output_matrix
  highCI_matrix <- output_matrix

  # Loop over the start and stop positions of the mismatch stretch.
  for (mm_start in 1:21) {
    for (mm_stop in mm_start:21) {
      # Make the string corresponding to the perfect match with the appropriate
      # stretch of mismatches.
      regex_str <- SubfunctionCall(GetPerfectTargetMismatchStretch)
      # Subset the data.

      data_subset <- data[sprintf("AAAAA%sAAAAA", regex_str), 1:3]
      output_matrix[21 - mm_start + 1, mm_stop] <- as.numeric(data_subset[1])
      lowCI_matrix[21 - mm_start + 1, mm_stop] <- as.numeric(data_subset[2])
      highCI_matrix[21 - mm_start + 1, mm_stop] <- as.numeric(data_subset[3])
    }
  }
  list(output_matrix, lowCI_matrix, highCI_matrix)
}


PlotBeckerFigure3C <- function(xpos=20, ypos=20, height=5, width=5,
                               pdf.plot=FALSE) {
  # Define the bounds for the
  output_miR21 <- GetMismatchStretchVariantBecker("miR-21")[[1]]
  output_let7 <- GetMismatchStretchVariantBecker("let-7a")[[1]]
  output_miR21 <- t(output_miR21[, 21:1])[, 21:1]

  output_full <- matrix(NaN, nrow=23, ncol=23)

  for (row_i in 1:21) {
    for (col_i in 1:(21 - row_i + 1)) {
      output_full[row_i, col_i] <- output_miR21[row_i, col_i]
    }
  }

  for (row_i in 1:21) {
    for (col_i in (21 - row_i + 1):21) {
      output_full[row_i + 2, col_i + 2] <- output_let7[row_i, col_i]
    }
  }
  xlefts <- rep(seq(0, 22), each=23) 
  xright <- xlefts + 1
  ybottom <- rep(seq(22, 0), 23)
  ytop <- ybottom + 1
  output <- output_full

  output[which(output >= 5000)] <- 5000
  output <- log10(output)
  output <- output - min(output, na.rm=TRUE)
  output <- output/max(output, na.rm=TRUE)
  print(output[10:21, 10:21])
  col.inds <- round(output*100)
  print(col.inds)

  # Parameters specifying the start and end of the rainbow colorspace usee:
  c_s <- 0.9 # color_start
  c_e <- 0.70  # color_end
  color.dist = plasma(100)
  col.inds <- sapply(col.inds, function(col.ind) {
    min(max(1, col.ind), 100)
    })

  SubfunctionCall(GenericFigureSaveFile)
  xmin <- 0
  xmax <- 23
  ymin <- xmin
  ymax <- xmax
  BlankPlot()

  rect(xlefts, ybottom, xright, ytop, col=color.dist[col.inds], border=NA)
  # rect(xlefts, ybottom, xright, ytop, col=color.dist[col.inds])
  xy <- GetPlotFractionalCoords(0, 1)
  text(0:20 + 0.5, xy[2], labels=1:21, xpd=NA, adj=c(0.5, 0), cex=0.8)

  

  xy <- GetPlotFractionalCoords(0, 1)
  text(xy[1], 0:20 + 2.5, labels=1:21, xpd=NA, adj=c(1, 0.5), cex=0.8)
  xy <- GetPlotFractionalCoords(1, 1)
  text(xy[1], 0:20 + 0.5, labels=1:21, xpd=NA, adj=c(0, 0.5), cex=0.8)
  xy <- GetPlotFractionalCoords(0, 0)
  text(0:20 + 2.5, xy[2], labels=1:21, xpd=NA, adj=c(0.5, 1), cex=0.8)

  # Left side label.
  xy <- GetPlotFractionalCoords(-0.075, 0.55)
  text(xy[1], xy[2], labels="Ending complement match", srt=90, xpd=NA,
       adj=c(0.5, 0))
  # Right site label
  xy <- GetPlotFractionalCoords(1.075, 0.45)
  text(xy[1], xy[2], labels="Beginning complement match", srt=270, xpd=NA,
       adj=c(0.5, 0))
  # Top side label
  xy <- GetPlotFractionalCoords(0.45, 1.075)
  text(xy[1], xy[2], labels="Beginning complement match", xpd=NA,
       adj=c(0.5, 0))
  # Right site label
  xy <- GetPlotFractionalCoords(0.55, -0.075)
  text(xy[1], xy[2], labels="Ending complement match", xpd=NA, adj=c(0.5, 0))

  xy <- GetPlotFractionalCoords(0, -0.07)
  text(xy[1], xy[2], label="let-7a", adj=c(0, 0), xpd=NA)

  xy <- GetPlotFractionalCoords(0.01, 0.025)
  text(xy[1], xy[2], label="miR-21", adj=c(1, 0), xpd=NA)
  segments(-0.5, -0.5, 23, 23, xpd=NA)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

PlotBeckerFigure3D <- function(xpos=20, ypos=20, height=4, width=7,
                               average_multiple=FALSE,
                               pdf.plot=FALSE) {
  # Load the median, lower CI, and upper CI limit matrices
  output_miR21 <- GetMismatchStretchVariantBecker("miR-21")
  message("Made miR-21 data matrices.")
  median_miR21 <- output_miR21[[1]]
  lowCI_miR21 <- output_miR21[[2]]
  highCI_miR21 <- output_miR21[[3]]

  output_let7 <- GetMismatchStretchVariantBecker("let-7a")
  message("Made let-7a data matrices.")
  median_let7 <- output_let7[[1]]
  lowCI_let7 <- output_let7[[2]]
  highCI_let7 <- output_let7[[3]]



  SubfunctionCall(GenericFigureSaveFile)
  xmin <- 0
  xmax <- 21
  ymin <- -17
  ymax <- -10
  BlankPlot()

  AddLinearAxis(1, tick.space=1, label.space=1, label="Final complementary base")
  AddLinearAxis(2, tick.space=1, label.space=1,
                label=expression(Delta*italic(G)~"(kcal/mol)"))

  R <- 1.987e-3 # in kcal K-1 mol-1
  T <- 310.15 # in K
  if (average_multiple) {
    data_let7 <- R*T*c(
      rowMeans(log(median_let7[20:1, 16:21]*10^-12), na.rm=TRUE), log(10*10^-12)
    )
    data_miR21 <- R*T*c(
      rowMeans(log(median_miR21[20:1, 16:21]*10^-12), na.rm=TRUE), log(10*10^-12)
    )
    text_str <- "Average of non-complementarity\nthroughs positions 15-21"
  } else {
    data_let7 <- R*T*log(c(median_let7[20:1, 21], 10)*10^-12)
    data_let7_lowCI <- R*T*log(c(lowCI_let7[20:1, 21], 10)*10^-12)
    data_let7_highCI <- R*T*log(c(highCI_let7[20:1, 21], 10)*10^-12)

    data_miR21 <- R*T*log(c(median_miR21[20:1, 21], 10)*10^-12)
    data_miR21_lowCI <- R*T*log(c(lowCI_miR21[20:1, 21], 10)*10^-12)
    data_miR21_highCI <- R*T*log(c(highCI_miR21[20:1, 21], 10)*10^-12)
    text_str <- "Non-complementarity\nthrough position 21"
  }
  rect(xleft=9.5, ybottom=-17, xright=13.5, ytop=-10, border=NA,
       col=ConvertRColortoRGB("gray", alpha=0.2))
  Points(1:21, data_let7, col="goldenrod", line=TRUE)
  Points(1:21, data_let7_lowCI, col="goldenrod", line=TRUE, lty=2)
  Points(1:21, data_let7_highCI, col="goldenrod", line=TRUE, lty=2)
  Points(1:21, data_miR21, col="forestgreen", line=TRUE)
  Points(1:21, data_miR21_lowCI, col="forestgreen", line=TRUE, lty=2)
  Points(1:21, data_miR21_highCI, col="forestgreen", line=TRUE, lty=2)

  xy <- GetPlotFractionalCoords(0.95, 0.95)
  text(xy[1], xy[2], label=text_str, adj=c(1, 1))
  xy <- GetPlotFractionalCoords(0.025, 0.3)
  legend(xy[1], xy[2], legend=c("miR-21", "let-7a"), bty="n", lwd=1,
         pch=19, xpd=NA, col=c("forestgreen", "goldenrod"))

  if (class(pdf.plot) == "character") {
    dev.off()
  }
}



PlotBeckerFigure3DAlt <- function(output_miR21, output_let7, list_Kds_miR21, list_Kds_let7, xpos=20, ypos=20,
                               height=4, width=7, pdf.plot=FALSE) {
  # Define the bounds for the 


  SubfunctionCall(GenericFigureSaveFile)
  xmin <- 0
  xmax <- 21
  ymin <- -17
  ymax <- -10
  BlankPlot()

  AddLinearAxis(1, tick.space=1, label.space=1, label="Final complementary base")
  AddLinearAxis(2, tick.space=1, label.space=1,
                label=expression(Delta*italic(G)~"(kcal/mol)"))

  R <- 1.987e-3 # in kcal K-1 mol-1
  T <- 310.15 # in K

  Kds_let7a <- sapply(list_Kds_let7, function(Kd_data) {
    GeoMean(Kd_data[, 1])
  })

  Kds_miR21 <- sapply(list_Kds_miR21, function(Kd_data) {
    GeoMean(Kd_data[, 1])
  })

  names_normal <- c("13mer-m1.13", "12mer-m1.12", "11mer-m1.11", "10mer-m1.10",
                    "9mer-m1.9", "8mer", "7mer-A1", "6mer-A1")

  names_alt <- paste0(names_normal, "_alt", sep="")

  print(names_normal)
  print(names_alt)

  print(Kds_let7a)
  data_let7 <- R*T*log(Kds_let7a[names_normal]*10^-12)
  # data_let7_alt <- R*T*log(Kds_let7a[names_alt]*10^-12)

  data_miR21 <- R*T*log(Kds_miR21[names_normal]*10^-12)
  # data_miR21_alt <- R*T*log(Kds_miR21[names_alt]*10^-12)


  data_let7_paper <- R*T*log(c(output_let7[20:1, 21], 10)*10^-12)
  data_miR21_paper <- R*T*log(c(output_miR21[20:1, 21], 10)*10^-12)
  rect(xleft=9.5, ybottom=-17, xright=13.5, ytop=-10, border=NA,
       col=ConvertRColortoRGB("gray", alpha=0.2))

  Points(1:21, data_let7_paper, col=ConvertRColortoRGB("goldenrod", alpha=0.2),
         line=TRUE)
  Points(1:21, data_miR21_paper,
         col=ConvertRColortoRGB("forestgreen", alpha=0.2), line=TRUE)
  print(data_let7)
  Points(13:6, data_let7, col="goldenrod", line=TRUE)
  Points(13:6, data_miR21, col="forestgreen", line=TRUE)

  xy <- GetPlotFractionalCoords(0.95, 0.95)
  text(xy[1], xy[2], label="Average of all sequences\nsatisfying criteria",
       adj=c(1, 1))
  xy <- GetPlotFractionalCoords(0.025, 0.3)
  legend(xy[1], xy[2], legend=c("miR-21", "let-7a"), bty="n", lwd=1, pch=19,
         col=c("forestgreen", "goldenrod"), xpd=NA)


  if (class(pdf.plot) == "character") {
    dev.off()
  }
}


PlotPairwiseReporterCounts <- function(mirna, experiment, condition, xpos=xpos,
                                       ypos=ypos, combined=TRUE, buffer=FALSE,
                                       pdf.plot=FALSE) {
  # Get the two signal data sets:
  if (experiment == "twist_reporter_assay") {
    exp_label <- "V1"
  } else {
    exp_label <- "V2"
  }  
  signal1 <- ReporterCounts(mirna, experiment, condition, rep=1, tpm=TRUE)
  signal2 <- ReporterCounts(mirna, experiment, condition, rep=2, tpm=TRUE)
  # Get the two background data sets:
  bg1 <- GetBackGroundReporterCounts(mirna, experiment, rep=1,
                                     tpm=TRUE)
  bg2 <- GetBackGroundReporterCounts(mirna, experiment, rep=2,
                                     tpm=TRUE)
  # Set up the arguments to get the correct Kd values:
  if (mirna == "miR-1" || mirna == "miR-7") {
    combined <- FALSE
  }
  if (mirna == "miR-1") {
    buffer <- TRUE
  }
  if (mirna == "miR-7") {
    mirna <- "miR-7-23nt"
    experiment <- "equilibrium2_nb"
  } else {
    experiment <- "equilibrium"
  }
  kds <- SubfunctionCall(EquilPars)
  # Relabel miR-7:
  if (mirna == "miR-7-23nt") {
    mirna <- "miR-7"
  }
  # Get the indeces for the miRNA in question:
  mirna_site_inds <- grep(paste0(mirna, "_"), rownames(signal1))
  # Strip the site names off of the sites:
  site_names <- sapply(rownames(bg1), function(mir_site) {
    unlist(strsplit(mir_site, split="_"))[2]
  })
  # Get the colors for the points:
  cols <- rep(ConvertRColortoRGB("gray", alpha=0.1), nrow(bg1))
  cols[mirna_site_inds] <- kSiteColors[site_names[mirna_site_inds]]
  # Make the log2fold-change for each rep, and average them:
  l2fc1 <- log(rowSums(signal1)/rowSums(bg1), 2)
  l2fc2 <- log(rowSums(signal2)/rowSums(bg2), 2)
  # Make a Kd vector, with 1 as the default value for each:
  SubfunctionCall(FigureSaveFile)
  xmin <- -0.5
  xmax <- 0.1
  ymin <- -0.5
  ymax <- 0.1
  BlankPlot()
  AddLinearAxis(1, label.space=0.1, tick.space=0.02, label="log2(fold change); rep 1")
  AddLinearAxis(2, label.space=0.1, tick.space=0.02, label="log2(fold change); rep 2")
  Points(l2fc1, l2fc2, col=cols)

  xy <- GetPlotFractionalCoords(0.9, 0.95)
  text(xy[1], xy[2], label=mirna)
  xy <- GetPlotFractionalCoords(0.05, 0.95)
  text(xy[1], xy[2], label=exp_label, adj=0)

  xy <- GetPlotFractionalCoords(0.95, 0.05)
  AddCorrelationToPlot(l2fc1[mirna_site_inds], l2fc2[mirna_site_inds], xy[1],
                       xy[2], rsquared=TRUE, adj=1)


  if (class(pdf.plot) == "character") {
    dev.off()
  }
}


