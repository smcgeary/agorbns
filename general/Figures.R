source("general/general.R")
library(wCorr)

pt_cex_final <- 1.35
legend_pch <- 19
line_dash_length <- 2
occupancy_schematic_bg_cex <- 0.45
occupancy_schematic_bg_col <- "#989898"
occupancy_mirna_seq_col <- "#958872"

AGO_mir_label <- "[AGO2-miRNA] (pM)"
AGO_mir1_label <-"[AGO2-miR-1]"
AGO_mir1_label_no_conc <- "AGO2-miR-1 library"


slopes_lm <- rep(NA, 7)

names(slopes_lm) <- c(kMirnas, "flanks")




source("general/GenericFigures.R")

################################################################################
# FIGURE 1
################################################################################

# 1C____________________________________________________________________________
PlotEquilSiteWithInput <- function(mirna, column, experiment="equilibrium",
                                   n_constant=5, sitelist="resubmissionfinal",
                                   uniq=FALSE, combined=TRUE, singleonly=TRUE,
                                   height=4.5, width=4.5, buffer=FALSE,
                                   pdf.plot=FALSE) {
  sXc <- SubfunctionCall(SitesXCounts)
  x <- Norm(sXc[,1 + combined])
  y <- Norm(sXc[,column])
  site.colors <- kSiteColors[rownames(sXc)]
  SubfunctionCall(FigureSaveFile)
  xmin <- 5e-4
  xmax <- 1
  ymin <- xmin
  ymax <- xmax
  BlankPlot(log="xy")
  xmax <- 1.00
  # Make x=y line:
  segments(xmin, ymin, xmax, ymax, lty=line_dash_length)
  # Make the lines connecting the points to the x = y line:
  segments(x, x, x, y, lty=line_dash_length, col=site.colors)
  # Make axes:
  AddLogAxis(1, label="Input library (%)", percent=TRUE)
  AddLogAxis(2, label="AGO-bound library (%)", percent=TRUE)
  # Add the points to the plot:
  Points(x, y, col=site.colors)
  R <- y/x
  names(R) <- rownames(sXc)
  A.stock <- kAgoStock[mirna, "equilibrium"] 
  A <- as.numeric(colnames(sXc)[column])/100*A.stock
  A.pM <- round(A*1000, 1)
  xy <- GetPlotFractionalCoords(0.05, 0.95, log='xy')
  ago_text <- " pM AGO2-"
  text(xy[1], xy[2], labels=paste0(A.pM,ago_text, mirna),
       adj=c(0, 1))
  legend.coords <- GetPlotFractionalCoords(0.95, 0.025, log='xy')

  Legend(legend.coords, legend=ConvertTtoUandMmtoX(rownames(sXc)),
         col=kSiteColors[rownames(sXc)], xjust=1, yjust=0, y.intersp=0.9)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

# 1D&E__________________________________________________________________________
PlotSiteEnrichments <- function(mirna, experiment="equilibrium", n_constant=5,
                                sitelist="resubmissionfinal", uniq=FALSE,
                                combined=TRUE, L=FALSE, remove_sites=TRUE,
                                compcorrect=FALSE, wobble=FALSE, tpomit=FALSE,
                                tpomit2=FALSE, buffer=FALSE, singleonly=TRUE,
                                plot.nconstant=FALSE, flowthrough=FALSE,
                                added.text=FALSE, datalines=FALSE,
                                modellines=TRUE, leg_rescale=1, write_kds=FALSE,
                                height=4.5, width=8.1, xpos=20, ypos=20,
                                pdf.plot=FALSE) {
  ################## LOAD AND PROCESS COUNT DATA ###############################
  sXc <- SubfunctionCall(SitesXCounts)
  removed_sites <- GetRemovedSites(sXc)
  if (length(removed_sites) != 0 & remove_sites) {
    sites_keep <- setdiff(rownames(sXc), removed_sites)
    sXc <- sXc[sites_keep, , drop=FALSE]
  }
  data <- GetDataEquil(sXc, combine_reps=TRUE)
  l <- SubfunctionCall(GetInputEquil)
  if (L) {
    l <- l/sum(l) * as.numeric(L)
  }
  data.R <- EquilEnrichments(data, l)
  #################### LOAD AND PROCESS KD MEAUREMENTS #########################
  pars.matrix <- SubfunctionCall(EquilPars)
  # Portion where the `0.1265` concentration column is removed from the sXc
  # matrix if the `tpomit2` flag is true.
  if (mirna == "miR-1" & experiment == "equilibrium_tp" & tpomit2) {
    sXc <- sXc[, -8]
  }
  # Correct mirna name for miR-7.
  if (mirna %in% c("miR-7-23nt", "miR-7-24nt", "miR-7-25nt")) {
    mirna_temp <- "miR-7"
  }
  # Removed the sites that are ascribed to competitor oligo enrichment.
  if (length(removed_sites) != 0 & remove_sites) {
    # Get the names of the last two rows of the parameter matrix.
    non_kd_pars <- grep(mirna, rownames(pars.matrix), value=TRUE)
    # Assemble the vector of all of the names to keep, both the sites that
    # are ambiguous sites and the two non-kD parameters.
    names_keep <- c(paste0(sites_keep, "_Kd"), non_kd_pars)
    # Subset the matrix
    pars.matrix <- pars.matrix[names_keep, , drop=FALSE]
  }
  # Log10-transform the data, and remove the `_Kd` and `_$(mirna)` suffix.
  pars.model <- log10(pars.matrix$Mean)
  names(pars.model) <- rownames(pars.matrix)
  names(pars.model)[nrow(sXc) + 1] <- "bg"
  names(pars.model)[nrow(sXc) + 2] <- "AGO"
  # Get the model curves.
  model <- SubfunctionCall(EquilSingleSiteModelFreq, pars=pars.model,
                           A.dil=A.dil.model)
  model.R <- EquilEnrichments(model, l)

  # Ago dilution in the data:
  A.dil.data <- sapply(colnames(data), as.numeric)
  # Get the 
  A.stock.measured <- kAgoStock[mirna, experiment]
  A.stock.pars <- 10^pars.model["AGO"]
  pM_from_dil <- A.stock.measured*1000/100
  A.pM.data <- A.dil.data*pM_from_dil

  xmin <- signif(min(A.pM.data)/(10^0.25), 1)
  xmax <- signif(max(A.pM.data)*(10^0.75), 1)
  A.dil.model <- exp(seq(log(min(A.pM.data)/(10^0.25)/pM_from_dil),
                         log(max(A.pM.data)*(10^0.25)/pM_from_dil),
                     length=100))

  A.pM.model <-A.dil.model/100*A.stock.measured *1000
  # Set up the plotting limits
  if (flowthrough) {
  	ymin <- 0.01
  	ymax <- 1
  } else {
	  ymin <- 0.2
	  ymax <- 300
  }
  SubfunctionCall(FigureSaveFile)
  BlankPlot(log="xy", adjusted=TRUE)
  ######### Make the title ####################################################
  xy <- GetPlotFractionalCoords(fx=0.05, fy=0.95, log='xy')
  if (added.text) {
    title.text <- paste0(mirna, "\n", experiment)
  } else {
    title.text <- mirna
  }
  text(xy[1], xy[2], labels = title.text, adj=c(0, 1))
  if (plot.nconstant) {
    xy <- GetPlotFractionalCoords(0.05, 0.90, log='xy')
    text(xy[1], xy[2], labels=sprintf("%s nt constant region", n_constant),
         adj=c(0, 1))      
  }
  ########## Make legend with dots #############################################
  rank_by_kd <- order(pars.matrix$Mean[1:(nrow(data)-1)])
  sites.ordered <- rownames(data)[rank_by_kd]
  cols <- kSiteColors[sites.ordered]
  xy <- GetPlotFractionalCoords(fx=0.85, fy=1, log='xy')
  Legend(xy, legend=c(ConvertTtoUandMmtoX(sites.ordered), "None"),
          col=c(cols, "black"), y.intersp=0.9*leg_rescale)
  # Add he lines that contain the Kd values for Figure 1D.
  if (write_kds) {
    ######################## ADD THE KD TEXT VALUES ############################
    xy <- GetPlotFractionalCoords(fx=0.85, fy=0, log='xy')
    temp_legend <- legend(x=xy[1], y=xy[2],
                          legend=rep("", length(sites.ordered) + 1),
                          col=c(kSiteColors[sites.ordered], "black"), seg.len=1,
                          pch=NA, lwd=1, bty="n",
                          y.intersp=0.86, yjust=0)
    kd.matrix <- pars.matrix[1:(nrow(data)-1),][rank_by_kd, ]
    kds <- kd.matrix$Mean
    # Part where text is formatted:
    kd_mags <- floor(log10(kds))
    # error_upper <- (kd.matrix$Upper_CI-kds)/(10^kd_mags)
    error_lower <- (kds - kd.matrix$Lower_CI)/(10^kd_mags)
    # error_mags_upper <- floor(log10(error_upper))
    error_mags_lower <- floor(log10(error_lower))
    text_matrix <- cbind(kds/(10^kd_mags), c(error_lower),
                             -error_mags_lower)
    # Plot the title to the legend ("Relative Kd:").
    x_text <- 10^temp_legend$text$x
    y_text <- 10^temp_legend$text$y
    n_k <- length(x_text)
    xy_title <- GetPlotFractionalCoords(fx=0.85, fy=0.46, log='xy')
    text(x=xy_title[1], y=xy_title[2],
         labels=bquote("Relative"~italic(K)[D]*":"),
         adj=c(0, 0))
    # Plot the Kd value with significant digit.
    pos1 <- 0.15
    text(x=x_text*10^pos1, y=y_text,
         labels=sprintf("%.*f", c(text_matrix[, 3], 1), c(text_matrix[, 1], 1)),
         adj=1)
    # Plot the +/- symbols.
    pos2 <- pos1 + 0.15
    text(x=x_text*10^pos2[-n_k], y=y_text[-n_k],
         labels=sapply(kd_mags, function(kd_mag) {
                       as.expression(bquote(.("")%+-%.("")))
                       }), adj=1)
    # Plot the error value with one significant digit.
    pos3 <- pos2 + 0.23
    text(x=x_text*10^pos3[-n_k], y=y_text[-n_k],
         labels=sprintf("%.*f", text_matrix[, 3], text_matrix[, 2]), adj=1)
    # Plot the order of magnitude (e.g., "10^-1").
    pos4 <- pos3 + 0.46
    text(x=x_text*10^pos4[-n_k], y=y_text[-n_k],
         labels=sapply(kd_mags, function(kd_mag) {
                       as.expression(bquote(.("")%*%10^.(kd_mag)))
                       }), adj=c(1, 0.4))
  }
  cols <- c(cols, None="black")
  if (flowthrough) {
  	segments(A.pM.data, ymin, A.pM.data, ymax, lty=2, col="gray")
  }
  # Add the points and lines:
  sapply(rownames(data), function(site) {
  	if (!flowthrough) {
	    Points(A.pM.data, data.R[site, ], col=cols[site], line=datalines)
		}
		if (!(datalines) & modellines) {
	    lines(A.pM.model, model.R[site, ], col=cols[site], xpd=NA)    		
		}
  })
  # Add axes to the plot.
  AddLogAxis(1, label=AGO_mir_label, adj=TRUE)
  AddLogAxis(2, label="Enrichment")
  # Close the plotting device and exit the function.
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

# 1F, 2A-D______________________________________________________________________
PlotSiteKds <- function(mirna, experiment="equilibrium", n_constant=5,
                        sitelist="resubmissionfinal", combined=TRUE, uniq=FALSE,
                        singleonly=TRUE, papersites=FALSE, remove_sites=TRUE,
                        collapse_AA=TRUE, adjusted_height=FALSE,
                        added.text=FALSE, L=FALSE, trim_mir_name=TRUE,
                        buffer=FALSE, compcorrect=FALSE, tpomit=FALSE,
                        wobble=FALSE, plot.nconstant=FALSE, showalignment=FALSE,
                        height=5, width=6, xpos=20, ypos=20, pdf.plot=FALSE) {
  # Get kds for all site-types of the mirna.
  sXc <- SubfunctionCall(SitesXCounts)
  print(head(sXc))
  pars.matrix <- SubfunctionCall(EquilPars)
  kd.sites <- gsub("_Kd", replace="", rownames(pars.matrix))
  # Assign removed_sites
  if (remove_sites) {
    removed_sites <- GetRemovedSites(sXc)  
  } else {
    removed_sites <- c()
  }
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
    if (length(AA_sites) !=0) {
      rownames(sXc_AA) <- sites_base_keep
      pars.matrix[inds_keep_base, 1:5] <- pars.matrix[AA_inds_keep, 1:5]
      rownames(pars.matrix_AA) <- sites_base_keep      
    }
  } else {
    removed_sites <- GetRemovedSites(sXc)
  }
  # removed_sites <- c()
  if (length(removed_sites) != 0) {
    sXc <- sXc[setdiff(rownames(sXc), removed_sites), , drop=FALSE]
    rownames(pars.matrix) <- kd.sites
    pars.matrix <- pars.matrix[setdiff(kd.sites, removed_sites), ,drop=FALSE]
  }
  kds <- pars.matrix[1:(nrow(sXc) - 1),]
  rownames(kds) <- rownames(sXc)[-nrow(sXc)]
  order_all <- order(kds$Mean)
  print(order_all)
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
  } else {
    height <- 4.5
  }
  SubfunctionCall(FigureSaveFile)
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
          # x_seg_l <- kd[2]
          # segments(unlist(x_seg_l), y, x, y, lty=2, col="blue")
          # str.width <- strwidth(label, units="figure")
          # pos_start <- (log(x) - log(xmax))/(log(xmin) - log(xmax))
          # pos_end <- pos_start + str.width + 0.04
          # x_r <- GetPlotFractionalCoords(pos_end, 1, log='x', inv='x')
          # segments(x_r[1], y, unlist(x_seg_r), y, lty=2)
        } else {
          x <- kd[5]* 1.1
          adj <- 1
          # segments(unlist(x_seg_l), y, unlist(x_seg_r), y, lty=2)
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
    xy_title <- GetPlotFractionalCoords(0.05, 0.90, inv='x')
    text(x=xy_title[1], y=xy_title[2],
         labels=sprintf("%s nt constant region", n_constant),
         adj=c(0, 1))      
  }
  text(title.xy[1], max(y), labels = title.text, adj=c(0, 0.5), xpd=NA)
  title_pos_1 <- apply(kds, 1, function(kd) {
    min(kd[2]/1.6, kd[3]/1.1)
  })

  title_pos_2 <- apply(kds, 1, function(kd) {
    max(kd[2]*1.6, kd[3]*1.1)
  })
  if (mirna == "miR-7-23nt" & experiment == "equilibrium2_nb") {
    pos_switch <- which(y <= 6)
  } else {
    pos_switch <- NULL
  }
  if (mirna == "miR-155" & sitelist == "papercutoff" & showalignment == TRUE) {
    sites_1 <- c("AACGAGGAA", "TAACGAGGA", "AACGAGGAT", "AACGAGGAG",
                 "ATAACGAGG")
    sites_2 <- c("AACTCAGCA", "TACTCAGCA", "CACTCAGCA", "CTCAGCAAT",
                 "CAACTCAGC")
    sites_3 <- c("AAATAAAGA", "ATAATAAAG", "AATAAAGAA", "AAAATAAAG",
                 "AATAAAGAT", "CAATAAAGA", "AATAAAGAC")
    sites_4 <- c("AGACGACAA", "AAACGACAA", "CGACAACTC", "ACGACAACA",
                 "ACGACAACT")
    sites_5 <- c("ATGACAACA", "TCGACAACA")
    weird_sites <- list(sites_1, sites_2, sites_3, sites_4, sites_5)
    cols_line <- c("red", "orange", "forestgreen", "blue", "purple2")
    x_right <- 0.3
    col_ind <- 1
    lapply(weird_sites, function(sites) {
      y_vals <- y[sites]
      x_lefts <- title_pos_1[sites]*2
      x_rights <- x_right
      segments(x0=x_lefts, y0=y[sites], x1=x_right, col=cols_line[col_ind])
      # print(x_right)
      # print(y[sites])
      # print(min(y[sites]))
      # print(max(y[sites]))
      segments(x0=x_right, y0=min(y[sites]), y1=max(y[sites]),
               col=cols_line[col_ind])
      print(x_right)
      x_right <<- x_right*1.2
      col_ind <<- col_ind + 1
    })
  }
  if ("CGCTTCCGCTG" %in% names(title_pos_1)) {
    title_pos_1["CGCTTCCGCTG"] <- kds["CGCTTCCGCTG", 5]*16
  }
  if (collapse_AA) {
    if (nrow(pars.matrix_AA) !=0) {
      AA_inds <- which(sites %in% rownames(pars.matrix_AA))
      sites[AA_inds] <- paste0("AA-", sites[AA_inds], sep="")  
    }    
  }
  par(xpd=NA)

  text(title_pos_1[setdiff(1:length(y), pos_switch)],
       y[setdiff(1:length(y), pos_switch)],
       labels=ConvertTtoUandMmtoX(sites[setdiff(1:length(y), pos_switch)]),
       adj=0)
  if (!is.null(pos_switch)){
    text(title_pos_2[pos_switch], y[pos_switch],
         labels=ConvertTtoUandMmtoX(sites[pos_switch]), adj=1)  
  }

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

# 1G____________________________________________________________________________
PlotSiteOccupancy <- function(mirna, experiment="equilibrium", n_constant=5,
                              sitelist="resubmissionfinal", plotlist=FALSE,
                              singleonly=TRUE, combined=TRUE, width=10,
                              buffer=FALSE, uniq=FALSE, adjusted_height=FALSE,
                              bgoff=FALSE, collapse_AA=FALSE, L=FALSE,
                              compcorrect=FALSE, wobble=FALSE,
                              pdf.plot=FALSE, xpos=20, ypos=20) {
  # Get the site count matrix and the parameter matrix.
  sXc <- SubfunctionCall(SitesXCounts)
  pars.matrix <- SubfunctionCall(EquilPars)
  # Replace the names of the sites with that of the matrix.
  kd.sites <- gsub("_Kd", replace="", rownames(pars.matrix))
  rownames(pars.matrix) <- kd.sites
  # Get the removed sites (complementary to the competitor) and removed their
  # site counts and their Kd values from the parameter matrix.
  removed_sites <- GetRemovedSites(sXc)
  removed_sites <- c()
  sXc <- sXc[setdiff(rownames(sXc), removed_sites), , drop=FALSE]
  pars.matrix <- pars.matrix[setdiff(kd.sites, removed_sites), ,drop=FALSE]

  # Check for the collapsed AA sites conditional.
  if (collapse_AA) {
    AA_sites <- grep("^AA-", rownames(sXc), perl=TRUE, value=TRUE)
    if (length(AA_sites) != 0) {
      # Get the inds o
      AA_inds <- grep("^AA-", rownames(sXc), perl=TRUE)
      # Get the base names of all AA sites.
      AA_sites_base <- gsub("^(AA-)(.*)$", AA_sites, replace="\\2", perl=TRUE)
      # Get the AA sites that are present as base sites in rownames(sXc)
      AA_sites_removed <- AA_sites[!(AA_sites_base %in% rownames(sXc))]
      # get the base sites that are present as base sites
      AA_sites_removed <- c()

    }
  } else {
    AA_sites_removed <- c()
  }
  removed_sites <- c(removed_sites, AA_sites_removed)

  num_sites <- nrow(sXc) - length(AA_sites_removed) - 1
  pars.model <- log10(pars.matrix$Mean)
  names(pars.model) <- rownames(pars.matrix)
  names(pars.model)[length(pars.model)] <- "AGO"
  names(pars.model)[length(pars.model) - 1] <- "bg"
  data <- GetDataEquil(sXc)
  l <- SubfunctionCall(GetInputEquil)
  # Ago dilutio in the data:
  A.stock.measured <- kAgoStock[mirna, "equilibrium"]
  A.dil.data <- sapply(colnames(data), as.numeric)
  A.dil.model <- exp(seq(log(min(A.dil.data)),
                         log(max(A.dil.data)),
                     length=100))

  pM_from_dil <- 10*A.stock.measured
  xmin <- signif(min(A.dil.data)*pM_from_dil, 1)
  xmax <- signif(ceiling(max(A.dil.data)*pM_from_dil), 1)
  A.dil.model <- exp(seq(log(xmin), log(xmax), length=100))/pM_from_dil
  A.pM.data <- A.dil.data*pM_from_dil
  A.pM.model <-A.dil.model*pM_from_dil
  pars <- pars.model
  model <- SubfunctionCall(EquilSingleSiteModelFreq, A.dil=A.dil.model,
                           addbg=FALSE)
  model <<- model
  model.points <- SubfunctionCall(EquilSingleSiteModelFreq, A.dil=A.dil.data,
                                  addbg=FALSE)
  model.points <<- model.points
  model.norm <- apply(model, 2, Norm)
  model.points.norm <- apply(model.points, 2, Norm)
  model.norm <<- model.norm
  model.points.norm <<- model.points.norm
  if (adjusted_height) {
    sXc_lsy6 <- SitesXCounts("lsy-6", sitelist=sitelist, uniq=uniq)
    sites_exclude_lsy6 <- GetRemovedSites(sXc_lsy6)
    num_sites_lsy6 <- nrow(sXc_lsy6) - length(sites_exclude_lsy6) - 1
    height_ratio <- (num_sites + 0.5)/(num_sites_lsy6 + 0.5)
    height <-  ((num_sites - 0.5) * 0.68 + 5)/4
    if (mirna == "let-7a") {
      real_height <- 1.6528
    } else if (mirna == "miR-155") {
      real_height <- 2.2083
    } else {
      real_height <- 1.5694
    }
    height_dist <- (num_sites) * 0.68
    height_dist_ref <- (num_sites_lsy6) * 0.68
    width <- 10
  } else {
    height_ratio <- 1
    height <- 4.5
  }
  # xmin <- ceiling(min(A.pM.model))
  xmax <- floor(max(A.pM.model))
  xmax <- (xmax/xmin)^height_ratio*xmin
  xmax_orig <- xmax
  ymin <- 0
  ymax <- height_ratio
  SubfunctionCall(FigureSaveFile)
  BlankPlot(log="x", adjusted=TRUE)
  if (length(grep("miR-7", mirna)) > 0) {
      title.text <- "miR-7"
    } else {
      title.text <- mirna
    }
  mirna.trim = paste0(strsplit(mirna, split = "-")[[1]][1:2],collapse="-")
  sites <- rownames(data)
  # if ("11mer-m9.19w9w11" %in% sites) {
  #   ind <- which(sites == "11mer-m9.19")
  #   sites[ind] <- "11mer-m9.19w9,11"
  # }
  baseline <- rep(0, length(A.pM.model))
  pars.model <- sort(pars.model[1:nrow(sXc)])
  sites_kd_order <- rev(gsub("(^.*)_Kd$", names(pars.model),
                             replacement="\\1"))

  sites <- sites_kd_order
  nc_sites <- setdiff(sites, c(kCanonicalSites, "None"))
  num_nc_sites <- length(nc_sites)

  if (mirna == "miR-1") {
    y_positions <- seq(length(nc_sites))*0.095  
    y_positions <- y_positions - min(y_positions) + 0.075
  } else {
    y_positions <- seq(length(nc_sites))*0.079  
    y_positions <- y_positions - min(y_positions) + 0.05
  }

  names(y_positions) <- setdiff(nc_sites, removed_sites)

  site_schematics <- c()

  # Printing of the occupancy plots:
  df.global <- data.frame(sitetype=c(), color=c())
  sapply(sites, function(site) {
    kSiteColors["None"] <- "gray70"
    x.polygon <- c(A.pM.model, rev(A.pM.model))
    y.polygon <- c(baseline, rev(baseline + model.norm[site, ]))
    df.site <- data.frame(sitetype=c(site),
                          color=c(ConvertRColortoRGB(kSiteColors[site],
                                                     alpha=0.3)))
    df.global <<- rbind(df.global, df.site)
    polygon(x=x.polygon, y=y.polygon,
            col=ConvertRColortoRGB(kSiteColors[site], alpha=0.3),
            border=FALSE, xpd=NA)    
    kSiteColors["None"] <- "black"
    if (site %in% c(kCanonicalSites, "None")) {
      xy <- GetPlotFractionalCoords(fx=0.4/height_ratio, fy=0.5/height_ratio,
                                    log='x')
      dist <- (A.pM.model - xy[1])^2
      x_ind <- which(dist == min(dist))
      x <- A.pM.model[x_ind]
      y <- baseline[x_ind] + model.norm[site, x_ind]/2 - strheight("A")/2
      col <- kSiteColors[site]
      adj <- c(0.5, 0)
      text(x, y, ConvertTtoUandMmtoX(site), col=col, adj=adj)
    } else if (!(site %in% removed_sites)) {
      if (pdf.plot == "1.G") {
        xy <- GetPlotFractionalCoords(fx=1.07/height_ratio,
                                      fy=y_positions[site]/height_ratio,
                                      log='x')
      } else {
        xy <- GetPlotFractionalCoords(fx=1.1/height_ratio,
                                      fy=y_positions[site]/height_ratio,
                                      log='x')
      }
      x_ind <- length(A.pM.model)
      y_average_right <- baseline[x_ind] + model.norm[site, x_ind]/2
      par(xpd=NA)
      if (pdf.plot == "1.G") {
        xy_int <- max(A.pM.model)*1.05
      } else {
        xy_int <- max(A.pM.model)*1.1
      }
      segments(x0=max(A.pM.model), y0=y_average_right, x1=xy_int,
               y1=y_average_right, lwd=0.5)
      segments(x0=xy_int, y0=y_average_right, x1=xy[1],
               y1=xy[2], lwd=0.5)
      if (pdf.plot == "1.G") {
        xy <- GetPlotFractionalCoords(fx=1.08/height_ratio,
                                      fy=y_positions[site]/height_ratio,
                                      log='x')
      } else {
        xy <- GetPlotFractionalCoords(fx=1.12/height_ratio,
                                      fy=y_positions[site]/height_ratio,
                                      log='x')
      }
      x <- xy[1]
      y <- xy[2] - strheight("A")/2
      col <- "black"
      output <- xy[2]
      names(output) <- site
      site_schematics <<- c(site_schematics, output)
      adj <- c(0, 0)
      text(x, y, ConvertTtoUandMmtoX(site), col=col, adj=adj)   
  }
    par(xpd=FALSE) 
    baseline <<- baseline + model.norm[site, ]
  })
  if (mirna == "miR-7-23nt") {
    mirna.str <- "miR-7"
  } else {
    mirna.str <- mirna
  }
  file_out <- sprintf("aesthetic_tables_for_kathy/%s_site_pastel_colors.txt",
                      mirna.str)
  write.table(file_out, x=df.global, quote=FALSE, row.names=FALSE, sep="\t",
              col.names=TRUE)
  par(xpd=NA)
  # Add Mirna schematic:
  mir_seq_rev <- StrRev(ConvertTtoUandMmtoX(kMirnaSeqs[mirna]))
  mir_seq_rev_list <- unlist(strsplit(mir_seq_rev, split=""))
  mir_len <- length(mir_seq_rev_list)
  num_3pterm <- mir_len - 16
  num_3p <- 4
  num_mid <- 4
  col_bg <- kSiteColors["mirna_bg"]
  col_seed <- kSiteColors["mirna_seed"]
  col_bg <- occupancy_mirna_seq_col

  if (adjusted_height) {
    scale_C <- 1.5  
  } else {
    scale_C <- 1
  }
  mirna_nt_col <- c(rep(col_bg, num_3pterm), rep(col_bg, num_3p),
                  rep(col_bg, num_mid), "black", rep(col_seed, 6), col_bg)
  seed_offset <- 0.02
  mirna_y_pos <- c(rep(0, mir_len - 8), rep(seed_offset, 7), 0)*scale_C
  seed_inds <- seq(mir_len - 8 + 1, mir_len)
  threep_inds <- seq(mir_len - 16 + 1, mir_len - 13 + 1)
  if (mirna == "miR-1") {
    mir_start_x <- 1.55
    y_dip <- 0.041
    y_dip1 <- -0.07
  } else {
    mir_start_x <- 1.82  
    y_dip <- 0.06
    y_dip1 <- -0.12/height_ratio  
  }
  nt_del_x <- 0.038*scale_C
  A.pM.model <<- A.pM.model
  x_fold <- max(A.pM.model)/min(A.pM.model)
  dist <- x_fold^(1/14.5)
  temp_stop <- min(A.pM.model)*(x_fold^3.2)
  temp_stop <- 10^par("usr")[2]
  temp_start <- 10^par("usr")[1]
  if (mirna == "miR-1") {
    dist <- (temp_stop/temp_start)^(1/60)  
    temp_stop <- temp_start*(temp_stop/temp_start)^(0.975)
  } else {
    dist <- (temp_stop/temp_start)^(1/60)
    temp_stop <- temp_start*(temp_stop/temp_start)^(0.956)
  }

  xpos_alt <- exp(log(temp_stop) - log(dist)*seq(length(mir_seq_rev_list) + 2))
  xmax_temp <- xmax
  xmax <- xmax_orig
  xy <- GetPlotFractionalCoords(fx=xmax, fy=y_dip1, log='x')
  x_dist <- rev(xpos_alt)
  x_dist <- x_dist[1:length(mir_seq_rev_list)]
  x_seed <- x_dist[seed_inds]
  text(x=x_dist,
       y=xy[2] + mirna_y_pos - strheight("A")/2, mir_seq_rev_list,
       adj=c(0.5, 0),
       col=mirna_nt_col)
  text(x=x_dist[seed_inds],
       y=xy[2] + mirna_y_pos[seed_inds] - y_dip - strheight("A")/2,
       seq(8, 1, by=-1), adj=c(0.5, 0),
       col=c(mirna_nt_col[seed_inds][-8], "black"))
  segments(x0=x_dist[seed_inds[-length(seed_inds)]],
           y0=xy[2] + seed_offset + 0.8*y_dip,
           y1=xy[2] + seed_offset + 1.8*y_dip)
  target_seq <- RevComplement(kMirnaSeqs[mirna], RNA=TRUE)
  target_seq_list <- unlist(strsplit(target_seq, split=""))
  x_dist_new <- x_dist
  x_dist_space <- x_dist[2]/x_dist[1]
  x_dist_new <- c(x_dist_new,
                  max(x_dist_new)*x_dist_space,
                  max(x_dist_new)*x_dist_space*x_dist_space)
  target_col <- rep(occupancy_schematic_bg_col, length(x_dist_new))
  # Map giving the start and stop coordinates of the seed sites:
  lr_seedsites <- matrix(c(1, 8,
                           2, 8,
                           1, 7,
                           2, 7,
                           3, 8,
                           1, 6,
                           1, 5,
                           4, 8), nrow=8, ncol=2, byrow=TRUE,
                          dimnames=list(c("8mer", "7mer-m8", "7mer-A1", "6mer",
                                          "6mer-m8", "6mer-A1", "5mer-A1",
                                          "5mer-m8"), c("start", "stop")))
  # This matrix defines the lengths of the ends of the three-prime sites.
  lr_3psites <- matrix(c(8, 18,   # 1
                         9, 19,  # 2
                         10, 20, # 3
                         11, 21, # 4
                         12, 22, # 5
                         13, 23, # 6
                         7, 17,  # 7
                         9, 18,  # 8
                         10, 19, # 9
                         11, 20, # 10
                         12, 21, # 11
                         13, 22, # 12
                         14, 23, # 13
                         9, 17,  # 14
                         10, 18, # 15
                         11, 19, # 16
                         12, 20, # 17
                         13, 21, # 18
                         14, 22, # 19
                         15, 23), nrow=20, ncol=2, byrow=TRUE,
                          dimnames=list(c("11mer-m8.18", "11mer-m9.19", "11mer-m10.20",
                                          "11mer-m11.21", "11mer-m12.22",
                                          "11mer-m13.23", "10mer-m7.17",
                                          "10mer-m9.18", "10mer-m10.19",
                                          "10mer-m11.20", "10mer-m12.21",
                                          "10mer-m13.22", "10mer-m14.22",
                                          "9mer-m9.17", "9mer-m10.18",
                                          "9mer-m11.19", "9mer-m12.20",
                                          "9mer-m13.21", "9mer-m14.22",
                                          "9mer-m15.23"), c("start", "stop")))

  lr_all <- rbind(lr_seedsites, lr_3psites)

  wobbles <- c(A="G",C="U")
  ends <- c(A="B", C="D", G="H", T="V", U="V")
  # FIGURE OUT THE SITE SCHEMATICS IN R
  l_seed_tar <- GetRevMirnaIndex(mirna, 8)
  key <- "%mer-m[[:digit:]]+\\.[[:digit:]]+$"
  threep_sites <- grep("mer-m(?:7|8|9|10|11|12|13|14|15|16)+\\.[[:digit:]]+.*$",
                       names(site_schematics), perl=TRUE, value=TRUE)
  l_3p_mir <- 0
  r_3p_mir <- 0
  l_3p_range <- try(
    as.integer(gsub("^(.*mer-m)([[:digit:]]+)\\.([[:digit:]]+)(w*.*)$",
                    threep_sites, replace="\\2"))
  )
  r_3p_range <- try(
    as.integer(gsub("^(.*mer-m)([[:digit:]]+)\\.([[:digit:]]+)(w*.*)$",
                    threep_sites, replace="\\3"))
  )
  if (length(l_3p_range) > 0) {
    positional_sites <- cbind(l_3p_range, r_3p_range)
    positional_sites <- positional_sites[which(rowSums(positional_sites) > 16),]
    if (nrow(positional_sites) > 0) {
      l_3p_mir <- min(positional_sites[, 1])
      r_3p_mir <- max(positional_sites[, 2])  
    }    
  }
  r_3p_tar <- nchar(target_seq) - l_3p_mir + 1
  l_3p_tar <- nchar(target_seq) - r_3p_mir + 1

  # Printing of the site schematics:
  sapply(1:length(site_schematics), function(ind) {
    target_seq_temp <- target_seq_list
    target_nt_col <- c(rep("black", length(mir_seq_rev_list) - 7),
                       rep("blue", 6), "black")
    site <- names(site_schematics)[ind]
    # Remove left-purple flankin sequence sites, like for AA-8mer.
    if (length(grep("^[[:upper:]]+-", site, perl=TRUE)) != 0) {
      prim_left_flank <- gsub("^([[:upper:]]+)(-.*)$", site, replace="\\1")
      site <- gsub("^([[:upper:]]+-)(.*)$", site, replace="\\2")
    } else {
      prim_left_flank <- ""
    }
    # Remove right-hand flanking sequence in site, like for 8mer-AA.
    if (length(grep("^.*-[[:upper:]]+$", site, perl=TRUE)) != 0) {
      prim_right_flank <- gsub("^(.*-)([[:upper:]]+)$", site, replace="\\2")
      site <- gsub("^(.*)(-[[:upper:]]+)$", site, replace="\\1")
    } else {
      prim_right_flank <- ""
    }
    y_dist <- site_schematics[ind]
    l_tar <- 0
    r_tar <- 0
    site_match <- ""
    # Attempt to define what the site type is:
    sapply(rownames(lr_seedsites), function(site_type) {
      if (length(grep(site_type, site, perl=TRUE)) != 0
          & nchar(site_type) > nchar(site_match)) {
        site_match <<- site_type
      }
    })
    if (nchar(site_match) == 0) {
      sapply(rownames(lr_3psites), function(site_type) {
        if (length(grep(site_type, site, perl=TRUE)) != 0
            & nchar(site_type) > nchar(site_match)) {
          site_match <<- site_type
        }
      })
    }
    # Check if a site as been determined yet:
    if (nchar(site_match) > 0) {
      lr_site_mir <- lr_all[site_match,]
      l_mir <- lr_all[site_match, 1]
      r_mir <- lr_all[site_match, 2]
      r_tar <- GetRevMirnaIndex(mirna, l_mir)
      l_tar <- GetRevMirnaIndex(mirna, r_mir)
    }
    # Pull the site out of the site name:
    remainder <- gsub(site_match, site, replace="")
    if (nchar(remainder) > 1 & substr(remainder, 1, 1) == "-") {
      remainder <- substr(remainder, 2, nchar(remainder))
    }
    check_all <- unlist(sapply(c("mm", "w", "b"), grep, x=remainder))
    pos_all <- unlist(sapply(c("mm", "w", "b"), gregexpr, text=remainder))
    if (nchar(remainder) >= 1 & sum(check_all) == 0) {
      if (length(grep("mer", remainder)) != 0) {
        site_type <- strsplit(remainder, split="-")[[1]][2]
        lr_mir <- as.integer(strsplit(substr(site_type, 2, nchar(site_type)),
                                      split="\\.")[[1]])
        l_mir <- lr_mir[1]
        r_mir <- lr_mir[2]
        r_tar <- GetRevMirnaIndex(mirna, l_mir)
        l_tar <- GetRevMirnaIndex(mirna, r_mir)
        remainder <- ""
      }
    }
    # Assign wobble colors:
    if (pos_all[2] != -1) {
      pos_temp <- pos_all
      pos_temp[which(pos_temp == -1)] <- nchar(remainder) + 1
      start_pos <- pos_all[2] + 1
      end_pos <- min(pos_temp[-2]) - 1
      w_pos <- GetRevMirnaIndex(mirna, as.integer(substr(remainder,
                                                         start_pos,
                                                         end_pos)))
      w_nt <- wobbles[target_seq_temp[w_pos]]
      target_seq_temp[w_pos] <- w_nt
      target_nt_col[w_pos] <- "cyan1"
    }
    # Assign mismatch colors:
    if (pos_all[1] != -1) {
      pos_temp <- pos_all
      pos_temp[which(pos_temp == -1)] <- nchar(remainder) + 1
      end_pos <- min(pos_temp[-1]) - 1
      start_pos <- pos_all[1] + 3
      mm_pos <- GetRevMirnaIndex(mirna, as.integer(substr(remainder,
                                                          start_pos,
                                                          end_pos)))
      mm_nt <- ConvertTtoUandMmtoX(substr(remainder, pos_all[1] + 2, pos_all[1] + 2))
      target_seq_temp[mm_pos] <- mm_nt
      target_nt_col[mm_pos] <- "orangered"
    }
    # Right gray
    r_tar_start = r_tar
    l_tar_start = l_tar
    r_seed_bool <- (r_tar_start != 0 &
                    r_tar_start < length(target_seq_temp) &
                    l_tar_start >= l_seed_tar)
    r_3p_bool <- (r_tar_start < r_3p_tar & l_tar_start >= l_3p_tar)
    if ((r_seed_bool | r_3p_bool) & prim_right_flank == "" & !grepl("11mer", site)) {
      target_seq_temp[r_tar_start + 1] <- ends[target_seq_temp[r_tar_start + 1]]
      target_nt_col[r_tar_start + 1] <- target_col[1]
      r_tar <- r_tar + 1    
    }
    # Left gray
    l_seed_bool <- (l_tar_start > l_seed_tar & l_tar_start != 0)
    l_3p_bool <- (l_tar_start > l_3p_tar & r_tar_start <= r_3p_tar)
    if ((l_seed_bool | l_3p_bool) & prim_left_flank == "" & !grepl("11mer", site)) {
      target_seq_temp[l_tar_start - 1] <- ends[target_seq_temp[l_tar_start - 1]]
      target_nt_col[l_tar_start - 1] <- target_col[1]
      l_tar <- l_tar - 1
    }
    if (prim_left_flank != "") {
      l_range <- (l_tar - nchar(prim_left_flank)):(l_tar - 1)
      target_seq_temp[l_range] <- ConvertTtoUandMmtoX(strsplit(prim_left_flank,
                                                               split="")[[1]])
      target_nt_col[l_range] <- "purple"
      l_tar <- l_tar - nchar(prim_left_flank)
    }
    if (prim_right_flank != "") {
      r_range <- (r_tar + 1):(r_tar + nchar(prim_right_flank))
      target_seq_temp[r_range] <- ConvertTtoUandMmtoX(strsplit(prim_right_flank,
                                                               split="")[[1]])
      target_nt_col[r_range] <- "purple"
      r_tar <- r_tar + nchar(prim_right_flank)
    }
    t_inds <- seq(r_tar, l_tar)
    t_include <- seq(1, length(x_dist_new))
    if (length(t_inds) == 1 && t_inds == 0) {
      t_inds <- c()
    } else {
      t_include <- t_include[-t_inds]
    }
    # if (pdf.plot == "1.G") {
    #   t_include <- t_include[-c(1, 2)]
    # }
    points(x=x_dist_new[t_include],
           y=rep(y_dist, length(x_dist_new[t_include])),
           pch=20, lwd=0, cex=occupancy_schematic_bg_cex,
           col=target_col[t_include])
    print(y_dist)
    if (length(t_inds) != 0) { 
      text(x=x_dist_new[t_inds], y=y_dist - strheight("A")/2,
           target_seq_temp[t_inds], adj=c(0.5, 0),
           col=target_nt_col[t_inds])
    }
  })
  # xmax <- (xmax/xmin)^height_ratio*xmin
  AddLogAxis(1, label=AGO_mir_label, adj=TRUE,
             xmax=signif(floor(max(A.pM.data)), 1))
  if (mirna == "miR-1") {
  AddLinearAxis(2, tick.space=0.1, label.space=0.1,
                label="Fraction of AGO-bound RNA", adj=TRUE, line=0.2,
                blank_lab=TRUE, ymax=1)

  } else {
  AddLinearAxis(2, tick.space=0.1, label.space=0.1,
                label="Fraction of AGO-bound RNA", adj=TRUE, line=0.2,
                blank_lab=TRUE, ymax=1)
  }
  title.xy <- GetPlotFractionalCoords(fx=0.95/height_ratio,
                                      fy=0.05/height_ratio, log='x')
  text(title.xy[1], title.xy[2], labels=mirna.trim, adj=1)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

################################################################################
# FIGURE S1
################################################################################

# S1A___________________________________________________________________________
PlotAgoPrepPurity <- function(experiment="AGO_purity", unique=FALSE,
                              no_marker=FALSE, no_adapter=FALSE, geo_test=FALSE,
                              height=4.5, width=4.5, xpos=20, ypos=20,
                              pdf.plot=FALSE) {
  out_matrix <- matrix(0)
  if (geo_test) str_geotest <- "_geotest"
  else          str_geotest <- ""
  # Gets the expression data within the experiment:
  out <- do.call("cbind", lapply(c("miR-1", "miR-155"), function(mirna) {
    out <- do.call("cbind", lapply(seq(5, 7), function(i_s) {
      out <- do.call("cbind", lapply(c("I", "P"), function(cond) {
        condition <- paste0("S100", i_s, "_", cond, str_geotest)
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

  out <<- out

  # Takes out the markers and unpammped categories from the table.
  out_new <<- out
  unmapped  <- out["Unmapped", ]
  out <- out[!(rownames(out) %in% exclude_rows),]
  # Isolate the exogenously loaded miRNA in the S100 extract
  exog_rows <- c("miR-1", "miR-155")
  spike_rows <- c("dme-miR-14-5p", "xtr-miR-427")
  spike_df <- out[spike_rows, ]
  exog_df <- t(t(out[exog_rows, ])/colSums(spike_df))
  endog_df <- t(t(out[setdiff(rownames(out), c(exog_rows, spike_rows)), ])/colSums(spike_df))
  row_order <- order(-rowSums(endog_df[,c(5, 6)]))
  endog_df <- endog_df[row_order,]
  ind_mir4521 <- grep("miR-4521", rownames(endog_df))
  endog_df <- rbind(endog_df[-ind_mir4521,], endog_df[ind_mir4521,])
  endog_df <- rbind(endog_df[1:10, ], colSums(endog_df[11:nrow(endog_df), ]))
  rownames(endog_df)[nrow(endog_df)] <- "Remaining miRNAs"
  norm_df <- rbind(exog_df, endog_df)
  norm_df <- round(t(t(norm_df)/colSums(norm_df))*1e6, 1)
  SubfunctionCall(FigureSaveFile)
  xmin <- 0
  xmax <- 1
  ymin <- 0
  ymax <- 1
  BlankPlot()
  legend_labels <- rownames(norm_df)
  for (i in seq(length(legend_labels))) {
    label_i <- legend_labels[i]
    if (grepl("/", label_i)) {
      labels <- unlist(strsplit(label_i, split="/"))
      label_suffix <- sapply(labels, function(label) {
        if (substr(label, 1, 3) == "let") {
          base_string <<- "let-"
          unlist(strsplit(label, split="let-"))[2]
        } else {
          base_string <<- "miR-"
          unlist(strsplit(label, split="miR-"))[2]
        }
      })
      regex_mir <- "([0-9]+)"
      regex_letter <- "([a-z]?)"
      regex_dash_number <- "((-[1-9])?)"
      regex_3.5p <- "((-(?:3|5)p)?$)"
      regex_all <- paste0(regex_mir, regex_letter, regex_dash_number, regex_3.5p)
      mir_number <- gsub("^([0-9]+)([^0-9]+)", label_suffix, replace="\\1", perl=TRUE)
      mir_number <- gsub(regex_all, label_suffix, replace="\\1", perl=TRUE)
      mir_letter <- gsub(regex_all, label_suffix, replace="\\2", perl=TRUE)
      mir_hyphen_number <- gsub(regex_all, label_suffix, replace="\\3", perl=TRUE)
      mir_p <- gsub(regex_all, label_suffix, replace="\\4", perl=TRUE)
      mir_letter <- mir_letter[which(nchar(mir_letter) > 0)]
      suffix <- mir_number[1]
      if (length(mir_letter) > 0) {
        suffix <- paste0(suffix, mir_letter[1], "-", mir_letter[length(mir_letter)])
      }
      merged_miRNA <- paste0(base_string, suffix)
      legend_labels[i] <- merged_miRNA
    }
  }
  legend_labels[length(legend_labels)] <- "Other miRNAs"
  y_pos_text <- seq(13,1,-1)/16
  y_pos_dist <- y_pos_text[1] - y_pos_text[2]
  y_pos_col_labels <-  y_pos_text[1] + y_pos_dist
  text_cex <- 1
  verts <- c(-0.5, 0.3, 0.65, 1.15)

  column_pos <- (verts[-1] + verts[-length(verts)])/2
  y_pos_lines <- c(y_pos_text, y_pos_text[length(y_pos_text)] - y_pos_dist)

  # Add  column labels:
  text(x=column_pos,
       y=y_pos_col_labels + 0.5*y_pos_dist,
       labels=c("miRNA", "AGO2-miR-1", "AGO2-miR-155"), cex=text_cex, xpd=NA)
  # add row names:
  text(-0.15, y_pos_text, adj=0, labels=legend_labels, cex=text_cex, xpd=NA)
  text(verts[3] - 0.025, y_pos_text, adj=1, format(round(norm_df[, 5], 0), nsmall=0, big.mark=","), cex=text_cex, xpd=NA)
  text(verts[4] - 0.08, y_pos_text, adj=1, format(round(norm_df[, 6], 0), nsmall=0, big.mark=","), cex=text_cex, xpd=NA)
  # Horizontal lines
  segments(verts[1], y_pos_col_labels, verts[length(verts)], y_pos_col_labels, xpd=NA)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

# S1B___________________________________________________________________________
PlotMiR1_KmersCorrelation <- function(n_constant=5, kmer_len=9, buffer=TRUE,
                                      pdf.plot=FALSE, height=4.5, width=4.5) {
  mirna <- "miR-1"
  exps <- c("equil_pilot", "equil_pilot", "equilibrium", "equilibrium")
  conditions <- c("I", "L100A10", "I", "40") # I know this is the best cor.
  if (buffer) {
    str.buffer <- "_buffer3p"
  } else {
    str.buffer <- ""
  }
  counts <- apply(cbind(exps, conditions), 1, function(row) {
    path <- GetAnalysisPath(mirna, row[1], row[2], analysis_type="kmers_cutoff_resub",
                                 ext=sprintf("_%s_k%s%s", n_constant, kmer_len,
                                             str.buffer))
    print(path)
    fread(path)[[2]]
  })

  # counts <- as.matrix(do.call(cbind, counts))
  colnames(counts) <- conditions
  kmers <- GetKmerList(kmer_len)
  R_1 <- Norm(counts[, 2])/Norm(counts[, 1])
  R_2 <- Norm(counts[, 4])/Norm(counts[, 3])
  bg_col <- "gray70"
  cols <- rep(bg_col, length(R_1))
  order_cols <- seq(length(kSeedSites))
  names(order_cols) <- kSiteColors[rev(kSeedSites)]
  cols_per_site <- sapply(rev(kSeedSites), function(site) {
    kSiteColors[site]
    })
  names(cols_per_site) <- rev(kSeedSites)
  names(order_cols) <- kSeedSites
  for (site in rev(kSeedSites)) {
    seq <- GetSiteSeq(mirna, site)
    cols[grep(seq, kmers)] <- kSiteColors[site]
  }
  SubfunctionCall(FigureSaveFile)
  xmin <- 1e-1
  xmax <- 1e3
  ymin <- xmin
  ymax <- xmax
  BlankPlot(log="xy")
  x <- R_1
  y <- R_2
  xy <- GetPlotFractionalCoords(0.05, 0.95, log='xy')
  AddCorrelationToPlot(log(x), log(y), xy[1], xy[2],
                       rsquared=TRUE, adj=c(0, 1))
  inds_bg <- which(cols == bg_col)
  inds_bg <- inds_bg[seq(1, length(inds_bg), by=floor(length(inds_bg)/10000))]
  inds_bg <- inds_bg[1:1e4]
  inds_color <- unlist(sapply(cols_per_site, function(color) {
    which(cols == color)
  }))
  message("number of background points:")
  message(length(inds_bg))
  inds <- c(inds_bg, inds_color)
  x <- x[inds]
  y <- y[inds]
  cols <- cols[inds]

  # Make x=y line:
  abline(0, 1, lty=2)
  # Make axes:
  AddLogAxis(1,
             label=expression("9-nt"~italic(k)*"-mer enrichment; Replicate 1"))
  AddLogAxis(2,
             label=expression(text="9-nt"~italic(k)*"-mer enrichment; Replicate 2"),
             col="blue")
  # Add the points to the plot:
  Points(x, y, col=ConvertRColortoRGB(cols, alpha=0.5))        
  legend.coords <- GetPlotFractionalCoords(1, 0.025, log='xy')
  Legend(legend.coords,
         legend=c(ConvertTtoUandMmtoX(kSeedSites), "None"),
         col=c(kSiteColors[kSeedSites], "gray80"), xjust=1, yjust=0)

  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

PlotEnrichmentsAgainstKds <- function(mirna, experiment="equilibrium",
                                      n_constant=5,
                                      sitelist="resubmissionfinal",
                                      plotlist=FALSE, uniq=FALSE, combined=TRUE,
                                      L=FALSE, buffer=FALSE, singleonly=FALSE,
                                      leg_rescale=1, width=4.5, height=4.5,
                                      pdf.plot=FALSE) {
  # Load the data:
  sXc <- SubfunctionCall(SitesXCounts)
  if (mirna %in% c("miR-7-23nt", "miR-7-24nt", "miR-7-25nt")) {
    mirna_temp <- "miR-7"
  }
  pars.matrix <- SubfunctionCall(EquilPars)
  pars.model <- log10(pars.matrix$Mean)
  names(pars.model) <- rownames(pars.matrix)
  names(pars.model)[nrow(sXc) + 1] <- "bg"
  names(pars.model)[nrow(sXc) + 2] <- "AGO"
  data <- GetDataEquil(sXc)
  l <- SubfunctionCall(GetInputEquil)
  if (L) {
    l <- l/sum(l) * as.numeric(L)
  }
  # Ago dilution in the data:
  pars <- pars.model
  A.dil.data <- sapply(colnames(data), as.numeric)
  A.stock.measured <- kAgoStock[mirna, "equilibrium"]
  A.stock.pars <- 10^pars.model["AGO"]
  pM_from_dil <- A.stock.measured*1000/100
  A.pM.data <- A.dil.data*pM_from_dil
  xmin <- 0.001
  xmax <- 3
  data.R <- EquilEnrichments(data, l)
  # Set up the plotting limits
  ymin <- 4e-1
  ymax <- 1e3
  SubfunctionCall(FigureSaveFile)
  BlankPlot(log="xy", inv="x", adjusted=TRUE)
  # Generate tickmarks for axis.
  xy <- GetPlotFractionalCoords(fx=0.05, fy=0.95, log='xy', inv='x')
  title.text <- mirna
  text(xy[1], xy[2], labels = title.text, adj=c(0, 1))
  names(pars.model) <- gsub("(.*)_Kd", names(pars.model), replace="\\1")

  site_order <- order(pars.matrix$Mean[1:(nrow(data)-1)])
  kd.matrix <- pars.matrix[1:(nrow(data)-1),][site_order, ]
  sites.ordered <- rownames(data)[site_order]
  kds <- c(kd.matrix$Mean, 1)
  names(kds) <- c(rownames(kd.matrix), "None")
  cols <- rep("black", nrow(sXc))
  names(cols) <- sites.ordered
  cols <- c(kSiteColors[sites.ordered], None="black")

  sample_cols <- sapply(seq(25, 75, length.out=5), function(num) {
    sprintf("gray%s", ceiling(num))
    })
  names(sample_cols) <- colnames(data)

  xy <- GetPlotFractionalCoords(fx=0.05, fy=0.90, log='xy', inv='x')
  Legend(xy, legend=c(ConvertTtoUandMmtoX(sites.ordered), "None"),
         col=c(cols, "black"), y.intersp=0.9*leg_rescale)

  A_stock <- kAgoStock[mirna, "equilibrium"]
  dils <- as.numeric(colnames(data.R))
  A_j <- A_stock*dils/100*1000
  mags <- floor(log10(A_j))

  # vals <- round(A_j, sapply(mags, function(mag_i) {
  #   print(max(-mag_i + 1, 0))
  #   max(-mag_i + 1, 0)
  # }))
  vals <- signif(A_j, 2)
  vals_formated <- sapply(1:length(vals), function(i) {
    sprintf("%s pM ",format(vals[i], nsmall=max(-mags[i] + 1, 0)))
  })
  colnames(pars.matrix) <- vals_formated


  xy <- GetPlotFractionalCoords(fx=0.90, fy=0.375, log='xy', inv='x')
  Legend(xy, legend=vals_formated, col=rep(NA, length(vals_formated)), lty=0,
         y.intersp=0.9*leg_rescale, adj=1)

  names(kds) <- sapply(names(kds), function(name_kd) {
    substr(name_kd, 1, nchar(name_kd) - 3)
    })

  xy <- GetPlotFractionalCoords(fx=0.675, fy=0.375, log='xy', inv='x')
  legend(x=xy[1], y=xy[2], legend=rep("", ncol(data.R)), col=sample_cols, lty=1,
         seg.len=0.8, lwd=1.5, y.intersp=0.9*leg_rescale, bty="n")



  names(kds)[length(kds)] <- c("None")

  segments(xmax, ymin, xmin, ymax, lty=2)
  sapply(colnames(data.R), function(col) {
    lines(kds/kds["None"], data.R[names(kds),col]/data.R["None",col],
          col=sample_cols[col], xpd=NA, lwd=1.5)      
    Points(kds/kds["None"], data.R[names(kds),col]/data.R["None",col], col=cols)
  })
  AddLogAxis(1, label="Relative Kd", adj=TRUE)
  AddLogAxis(2, label="Enrichment")

  if (class(pdf.plot) == "character") {
    dev.off()
  }
}



# S1C___________________________________________________________________________
PlotPositionalEnrichment <- function(mirna, experiment="equilibrium",
                                     condition=40, n_constant=5,
                                     sitelist="resubmissionextendedfinal",
                                     sites=kSeedSites, showdil=FALSE, xpos=20,
                                     ypos=20, width=4.5, height=4.5, combI=FALSE,
                                     buffer=FALSE, pdf.plot=FALSE) {
  sXc <- SubfunctionCall(SingleSitesXCounts)
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

  SubfunctionCall(FigureSaveFile)
  xmin <- 2
  if (identical(sites, kSeedSites)) {
    xmax <- 40
    offsets <- c(2, 2, 1, 1, 2, 0)
  } else {
    # Regex to remove the region after the period.
    # Need to get rid of wobbles
    sites_use <- grep("w", sites_use, invert=TRUE, value=TRUE)
    print(sites_use)
    ending_positions <- as.integer(gsub("^(.*)\\.(.*)$", sites_use,
                                        replace="\\2", perl=TRUE))
    print(ending_positions)
    offsets <- -ending_positions - min(-ending_positions)
  }
  xmax <- 39

  names(offsets) <- sites_use

  ymin <- 1e-1
  ymax <- 1e3
  BlankPlot(log='y')
  print(sites_use)
  if (sites == k3PSites) {
    print("THREE PEE SITES")
    color_function <- colorRampPalette(c(kSiteColors["11mer-m9.19"],
                                     kSiteColors["9mer-m11.19"]))
    cols <- color_function(length(sites_use))
    names(cols) <- sites_use
  } else {
    cols <- kSiteColors[sites_use]
    names(cols) <- sites_use
  }
  for (site in sites_use) {
    segments(xmin, R_overall_2[site], xmax, R_overall_2[site],
             col=ConvertRColortoRGB(cols[site], alpha=0.5))

    frac_A <- A[site, ]/sum(sXc[, condition])
    frac_I <- I[site, ]/sum(sXc[, 1])
    # Portion of code that only takes the "N0–N36" portion of the enrichment profiles
    col_inds <- grep("N", names(frac_A))
    R <- R_all[site,]
    R_names <- names(R)
    R <- unlist(c(rep(NaN, offsets[site]), R[1:(length(R) - offsets[site])]))
    names(R) <- R_names
    R <- R[ceiling(xmin):floor(xmax)]
    Points(1:length(R), R, col=cols[site], line=TRUE)        
  }
  AddLinearAxis(1, tick.space=1, label.space=1, label="Position",
                alt_lab=    c(1,  5, 10, 15, 20, 25, 30, 35, 37),
                alt_lab_pos=c(6, 10, 15, 20, 25, 30, 35, 40, 42))
  AddLogAxis(2, label="Enrichment")
  xy <- GetPlotFractionalCoords(0.05, 1, log='y')
  mirna.split <- paste0(strsplit(mirna, split="-")[[1]][1:2], collapse="-")
  text(xy[1], xy[2], labels=mirna.split, adj=c(0, 1))  
  if (showdil) {
    xy <- GetPlotFractionalCoords(0.05, 1, log='y')
    text(xy[1], xy[2], labels=condition, adj=c(0, 1))    
  }
  xy <- GetPlotFractionalCoords(0, 0, log='y')
  Legend(xy, legend=sites_use[1:(length(sites_use)/2)],
         col=cols[1:(length(sites_use)/2)],
         xjust=0, yjust=0)
  xy <- GetPlotFractionalCoords(0.5, 0, log='y')
  Legend(xy, legend=sites_use[(length(sites_use)/2 + 1):length(sites_use)],
         col=cols[(length(sites_use)/2 + 1):length(sites_use)],
         xjust=0, yjust=0)

  if (class(pdf.plot) == "character") {
    dev.off()
  }
}



# S1C___________________________________________________________________________
PlotWorstSiteKdCrossValScatter <- function(mirna, experiment="equilibrium",
                                           n_constant=5, sitelist="resubmissionfinal",
                                           buffer=FALSE, combined=TRUE,
                                           singleonly=FALSE, papersites=FALSE,
                                           height=4.5, width=4.5,
                                           adjusted_height=FALSE, xpos=20,
                                           ypos=20, pdf.plot=FALSE) {
  # Get kds for all site-types of the mirna.
  sXc <- SubfunctionCall(SitesXCounts)
  pars.matrix <- SubfunctionCall(EquilPars, Xval=TRUE)

  if (papersites) {
    pars.matrix_paper <- SubfunctionCall(EquilPars, sitelist="paper")
    pars.matrix <- pars.matrix[intersect(rownames(pars.matrix_paper),
                                         rownames(pars.matrix)), ]
    sites <- gsub("_Kd", replace="", grep("_Kd", rownames(pars.matrix), value=TRUE))
    sXc <- sXc[sites, ]
  }
  xmin <- 1e-4
  xmax <- 10
  ymin <- xmin
  ymax <- xmax
  SubfunctionCall(FigureSaveFile)
  pars.matrix <- pars.matrix[, -1]
  A_stock <- kAgoStock[mirna, "equilibrium"]
  dils <- as.numeric(gsub("X.", replace="", colnames(pars.matrix)))
  A_j <- A_stock*dils/100*1000
  colnames(pars.matrix) <- A_j
  cor_mat <- cor(log(pars.matrix))
  min_cor_mat <- min(cor_mat)
  leave_bool <- FALSE
  BlankPlot(log='xy', inv='xy')
  for (i in seq(1, 5)) {
    for (j in seq(1, 5)) {
      if (cor_mat[i, j] == min_cor_mat) {
        kds_x <- pars.matrix[, i]
        kds_y <- pars.matrix[, j]
        AddLogAxis(1,
                   label=sprintf("Fitted parameters; %s pM sample removed",
                                 signif(A_j[i], 2)))
        AddLogAxis(2,
                   label=sprintf("Fitted parameters; %s pM sample removed",
                                 signif(A_j[j], 2)))
        abline(0, 1, lty=2)
        Points(kds_x, kds_y, col=kSiteColors[c(rownames(sXc), "bg", "AGO")])
        xy <- GetPlotFractionalCoords(0.05, 0.95, log='xy', inv='xy')
        AddCorrelationToPlot(log(kds_x), log(kds_y), xy[1], xy[2],
                             rsquared=TRUE, adj=c(0, 1))
        leave_bool <- TRUE
        break
      }
      if (leave_bool) {
        break
      }
    }
    if (leave_bool) {
      break
    }
  }
  ################################### LEGEND ###################################
  xy <- GetPlotFractionalCoords(1.12, 0, log='xy', inv='xy')
  legend_text <- c(sapply(kSeedSites, function(site) {
    as.expression(bquote(.(site)~italic(K)[D]))
    }), AGO_mir1_label, "[Nonspecific RNA]") 
  Legend(xy, legend=legend_text, col=kSiteColors[c(kSeedSites, "bg", "AGO")],
         xjust=1, yjust=0, y.intersp=0.9)
  #################################### EXIT ####################################
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

# S1D___________________________________________________________________________
PlotSiteKdCrossValMatrix <- function(mirna, experiment="equilibrium",
                                     n_constant=5, sitelist="resubmissionfinal",
                                     buffer=FALSE, combined=TRUE,
                                     singleonly=FALSE, papersites=FALSE,
                                     height=4.5, width=4.5, xpos=20, ypos=20,
                                     pdf.plot=FALSE) {
  # Get kds for all site-types of the mirna.
  sXc <- SubfunctionCall(SitesXCounts)
  pars.matrix <- SubfunctionCall(EquilPars, Xval=TRUE)
  if (papersites) {
    pars.matrix_paper <- SubfunctionCall(EquilPars, sitelist="paper")
    pars.matrix <- pars.matrix[intersect(rownames(pars.matrix_paper),
                                         rownames(pars.matrix)), ]
    sites <- gsub("_Kd", replace="", grep("_Kd", rownames(pars.matrix),
                                          value=TRUE))
    sXc <- sXc[sites, ]
  }
  SubfunctionCall(FigureSaveFile)
  pars.matrix <- pars.matrix[, -1]
  cor_mat <- cor(log(pars.matrix))

  A_stock <- kAgoStock[mirna, "equilibrium"]
  dils <- as.numeric(gsub("X.", replace="", colnames(pars.matrix)))
  A_j <- A_stock*dils/100*1000
  colnames(pars.matrix) <- sprintf("%s pM", signif(A_j, 2))
  cor_mat <- cor(log(pars.matrix))
  cor_title <- expression("Pairwise"~italic(r)^2~" when removing each sample:")
  corrplot(cor_mat^2, title=cor_title, type="lower", number.digits=3,
           method="number", diag=FALSE, addgrid.col="black", shade.col="white",
           addCoef.col="black", col="black", cl.pos="n", bg=NULL,
           tl.pos="lt", tl.col="black", mar=c(1, 0, 3, 1),
           number.font=1, cex.main=1)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

# S1E___________________________________________________________________________
PlotSalomonComparison <- function(main_text_seq=TRUE, pdf.plot=FALSE,
                                  width=4.5, height=4.5){
  # wee_et_al <- read.table("Wee_et_al_data.txt", sep="\t", header=TRUE)
  # The Kds in the salomon et al dataset are in nM.
  salomon_et_al <- read.table("salomon_et_al.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
  salomon_et_al <- salomon_et_al[1:11,]
  salomon_et_al$Sequence <- gsub("(?: |-)", salomon_et_al$Sequence, replace="", perl=TRUE)
  salomon_et_al$Sequence <- gsub("U", salomon_et_al$Sequence, replace="T", perl=TRUE)
  salomon_et_al$Sequence <- sapply(salomon_et_al$Sequence, function(i) {
    substr(i, 92, 92 + 11)
    })

  kds_1.4 <- exp(data.frame(fread(paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a",
                                     "/equilibrium/kds_PAPER/5_12mers_1-4_PAPER",
                                     "_logmean.txt"), sep="\t"), row.names=1)/log(10))

  kds_2.5 <- exp(data.frame(fread(paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a",
                                     "/equilibrium/kds_PAPER/5_12mers_2-5_PAPER",
                                     "_logmean.txt"), sep="\t"), row.names=1)/log(10))
  kds_3.6 <- exp(data.frame(fread(paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a",
                                     "/equilibrium/kds_PAPER/5_12mers_3-6_PAPER",
                                     "_logmean.txt"), sep="\t"), row.names=1)/log(10))

  kds_4.7 <- exp(data.frame(fread(paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a",
                                     "/equilibrium/kds_PAPER/5_12mers_4-7_PAPER",
                                     "_logmean.txt"), sep="\t"), row.names=1)/log(10))


  kds_5.8 <- exp(data.frame(fread(paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a",
                                     "/equilibrium/kds_PAPER/5_12mers_5-8_PAPER",
                                     "_logmean.txt"), sep="\t"), row.names=1)/log(10))
  all_12mer_kds <- list(kds_1.4, kds_2.5, kds_3.6, kds_4.7, kds_5.8)
  salomon_et_al_seed_only <- salomon_et_al[c(2, 5, 6, 7, 8, 9, 10, 11),]
  final_target_kds <- rep(NaN, nrow(salomon_et_al_seed_only))
  names(final_target_kds) <- salomon_et_al_seed_only$Sequence
  # Why is this here?
  for (target in salomon_et_al_seed_only$Sequence) {
    if (main_text_seq & target == "ATCTACAGCAAC") {
      target_use <- "ATCTACGACAAC"
    } else if (main_text_seq & target == "ATCTACCGAAAC") {
      target_use <- "ATCTACCAGAAC"
    } else {
      target_use <- target
    }
    target_kds <- c()
    for (i in seq(5)) {
      kds_inds <- grep(target_use, rownames(all_12mer_kds[[i]]))
      target_kds <- c(target_kds, all_12mer_kds[[i]][kds_inds,])
      final_target_kds[target] <- exp(mean(log(target_kds)))
    }
  }
  salomon_kds <- salomon_et_al_seed_only[, 7]/salomon_et_al_seed_only[1, 7]
  salomon_kds <- salomon_et_al_seed_only[, 7]

  print(salomon_kds)
  # final_target_kds <- final_target_kds*10

  # Make the second comparison:
  Becker_et_al_delG <- c(`8mer`=-13.70166171, `7mer-m8`=-12.89980191,
                        `7mer-A1`=-12.65070981, `6mer`=-12.10140861)
  R <- 1.987e-3 # in kcal K-1 mol-1
  T <- 310.15 # in K

  #Convert Becker ∆G to Kd in nM.
  Becker_et_al_Kds <- exp(Becker_et_al_delG/(R*T))*1e9
  # These Kds are relative Kds.
  McGeary_Lin_Kds <- SubfunctionCall(EquilPars, mirna="let-7a",
                                     experiment="equilibrium")
  McGeary_Lin_Kds <- McGeary_Lin_Kds[paste0(kCanonicalSites, "_Kd", sep=""), 1]
  SubfunctionCall(FigureSaveFile)
  # par(mar= c(3, 3, 2, 4))

  xmin <- 1e-3
  xmax <- 1e3
  ymin <- 1e-4
  ymax <- 1e1
  BlankPlot(log='xy', inv='xy')
  AddLogAxis(1, label="Kd (nM)")
  AddLogAxis(2, label="Relative Kd; AGO-RBNS")
  # AddLogAxis(4, label="Kd; Becker, et al. (pM)")
  data <- data.frame(salomon_name = salomon_et_al_seed_only$Site,
                     salomon_kds = salomon_kds,
                     rbns_name = c("8mer", "mm2-3", "mm3-4", "mm4-5", "mm5-6",
                                   "mm6-7", "6mer-A1", "7mer-A1"),
                     rbns_kds = final_target_kds,
                     cols = c("purple", "red", "orange", "orange", "green",
                              "forestgreen", "cyan", "blue"),
                              stringsAsFactors=FALSE)
  data <- data[!is.na(data$rbns_kds), ]
  # Make the linear model fits for each of the two comparisons:
  lm_becker <- lm(log10(Becker_et_al_Kds) ~ log10(McGeary_Lin_Kds))
  lm_becker <<- lm_becker
  lm_salomon <- lm(log10(data$salomon_kds) ~ log10(data$rbns_kds))
  x_lm_seq <- seq(-4, 1)
  b <- lm_salomon$coefficients[1]
  m <- lm_salomon$coefficients[2]
  y_salomon <- m*x_lm_seq + b
  x_becker <- seq(-4, -1.5, length.out=10)
  b <- lm_becker$coefficients[1]
  m <- lm_becker$coefficients[2]
  y_becker <- m*x_becker + b
  # Plot the two trend lines.
  lines(10^y_becker, 10^x_becker, lwd=0.5, lty=2)
  lines(10^y_salomon, 10^x_lm_seq, lwd=0.5, lty=2)
  # Plot the two data sets:
  # Points(McGeary_Lin_Kds, Becker_et_al_Kds, col=kSiteColors[kCanonicalSites],
  #        pch=1)
  Points(data$salomon_kds, data$rbns_kds, col=data$cols)
  Points(Becker_et_al_Kds, McGeary_Lin_Kds, col=kSiteColors[kCanonicalSites],
         pch=1, pt.lwd=1)

  # Salomon et al text
  # Title
  xy <- GetPlotFractionalCoords(0.05, 0.95, log='xy', inv='xy')
  # text(xy[1], xy[2], label="Becker et al., 2019", adj=c(0, 1))
  text(xy[1], xy[2], label=expression("Becker"~italic("et al")*"., 2019"),
       adj=c(0, 1))
  # Legend
  xy <- GetPlotFractionalCoords(0.025, 0.93, log='xy', inv='xy')
  legend(xy[1], xy[2], legend=kCanonicalSites, col=kSiteColors[kCanonicalSites],
         xjust=0, yjust=1, y.intersp=0.85, pt.lwd=1, pt.cex=1.2, pch=1, cex=1,
         bty="n")
  # Correlation
  xy <- GetPlotFractionalCoords(0.35, 0.65, log='xy', inv='xy')
  AddCorrelationToPlot(log(McGeary_Lin_Kds), log(Becker_et_al_Kds), xy[1],
                       xy[2], rsquared=TRUE, adj=c(1, 1))

  # Becker et al., 2019
  # Title
  xy <- GetPlotFractionalCoords(1.05, 0.53, log='xy', inv='xy')
  text(xy[1], xy[2], label=expression("Salomon"~italic("et al")*"., 2015"), adj=c(1, 1), xpd=NA)
  # Legend
  xy <- GetPlotFractionalCoords(1.05, 0.51, log='xy', inv='xy')
  Legend(xy, legend=data$rbns_name, col=data$cols, xjust=1, yjust=1,
         y.intersp=0.85)
  # Correlation
  xy <- GetPlotFractionalCoords(1.05, 0.85, log='xy', inv='xy')
  AddCorrelationToPlot(log(salomon_kds), log(final_target_kds), xy[1], xy[2],
                       rsquared=TRUE, adj=c(1, 1))
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}



################################################################################
# FIGURE 2
################################################################################

################################################################################
# FIGURE S2
################################################################################
PlotCompetitorOligoSiteSchematic <- function(mirna=NA,
                                             sitelist="papercutoff",
                                             buffer=FALSE, pdf.plot=FALSE,
                                             height=6, width=14.788) {
  SubfunctionCall(FigureSaveFile)
  xmin <- 0
  xmax <- 1
  ymin <- 0
  ymax <- 1
  par(mar=c(0, 0, 0, 0))
  BlankPlot()

  # Define the relevant x-coordinates of this schematic.
  x_mirna <- 0.06
  x_mirna_line <- x_mirna + 0.0075
  x_sites <- x_mirna_line + 0.04
  x_site_line <- x_sites + 0.0075

  x_l <- GetPlotFractionalCoords(0.11, 0)[1]
  x_r <- GetPlotFractionalCoords(0.41, 0)[1]
  x_m <- (x_l + x_r)/2
  x_del <- (x_r - x_l)/30

  x_comp <- 0.475

  y_pos_start <- 0.90
  y_pos <- y_pos_start
  x_offset <- 0
  x_offset_set <- 0.49

  y_comp_heading_pos <- 0.93
  text(x=x_comp, y=y_comp_heading_pos, labels="Length of\npairing (nt)",
       adj=c(0.5, 0))
  text(x=x_comp + x_offset_set, y=y_comp_heading_pos,
       labels="Length of\npairing (nt)", adj=c(0.5, 0))


  y_del_site <- 0.020
  y_del_mirna <- 0.05


  cols_comp <- rep(c("black", "red", "black"), times=c(31 - 8 - 6 + 1, 8, 6))
  for (mirna in c("miR-155", "lsy-6")) {
    y_mirna_top <- y_pos - 0.7*y_del_site
    # Get the competitor sequence
    str_competitor <- SequenceList[["competitor"]][mirna]
    # Get the vector of nucleotides.
    vec_competitor_rev <- rev(GetKmersFromString(str_competitor, 1))
    # Get the x positions for the competitor nucleotides.
    x_pos <- seq(x_l, x_r, length.out=length(vec_competitor_rev))
    # Label the competitor oligo on the plot.
    text(x_pos + x_offset, y=y_pos, labels=vec_competitor_rev, adj=c(0.5, 0),
         col=cols_comp)
    # Label the 3' end.
    text(x_l - x_del + x_offset, y=y_pos, labels=expression(paste(3*minute,"-")),
         adj=c(0.75, 0))
    # Label the 5' end.
    text(x_r + x_del + x_offset, y=y_pos, labels=expression(paste("-", 5*minute)),
         adj=c(0.25, 0))
    # Get the list of sites for that miRNA.
    sites <- unlist(read.table(
      file.path("AssignSiteTypes", "site_categories", sitelist,
                sprintf("sites.%s_%s.txt", mirna, sitelist)),
      stringsAsFactors=FALSE
    ))
    # Pull out the noncognate sites.
    noncognate_sites <- grep("mer", sites, invert=TRUE, value=TRUE)
    # Remove those sites that have fewer than 6 nucleotides of perfect
    # complementarity to the competitor oligo.
    rc_competitor_seq <- RevComplement(kCompetitorOligoSeqs[mirna])
    print(rc_competitor_seq)
    inds_keep <- which(
      sapply(noncognate_sites, function(site) {
        max(nchar(LargestCommonSubstring(rc_competitor_seq, site)))
      }) >= 6
    )
    noncognate_sites <- noncognate_sites[inds_keep]
    # Loop over the noncognate sites to print them.
    site_num <- 1
    y_pos <- y_pos - 0.2*y_del_site
    for (str_site in noncognate_sites) {
      # Define the vertical line giving the sites.
      y_pos <- y_pos - 0.2*y_del_site
      y_site_top <- y_pos - 0.5*y_del_site
      y_pos <- y_pos - y_del_site

      # str_site <- ConvertTtoU(str_site)
      # print(str_site)
      c_competitor_seq <- StrRev(rc_competitor_seq)
      rev_str_site <- StrRev(str_site)
      substr <- LargestCommonSubstring(rc_competitor_seq, str_site)
      # print(substr)
      substrs_competitor_oligo <- GetKmersFromString(rc_competitor_seq,
                                                     nchar(substr))
      int_start <- which(substrs_competitor_oligo == substr)
      substrs_site <- GetKmersFromString(str_site, nchar(substr))
      pos_site <- which(c(substrs_site) == substr)
      x_start <- x_l + x_del*(int_start - 1 - (pos_site - 1))
      vec_site <- GetKmersFromString(ConvertTtoU(str_site), 1)
      x_stop <- x_start + x_del*(length(vec_site) - 1)
      vec_xpos <- seq(x_start, x_stop, length.out=length(vec_site))

      vec_xpos <- seq(x_start, x_stop,
                           length.out=length(vec_site))
      cols <- rep("gray", nchar(str_site))
      cols[pos_site:(pos_site + nchar(substr) - 1)] <- "black"
      # Label the weird sites at the right position.
      text(vec_xpos + x_offset, y=y_pos, col=cols, labels=vec_site,
           adj=c(0.5, 0))
      text(x_comp + x_offset, y=y_pos, labels=sprintf("%s", nchar(substr)),
           adj=c(0.5, 0))
      y_site_bottom <- y_pos - 0.5*y_del_site
      y_mirna_bottom <- y_pos - 0.5*y_del_site
      y_pos <- y_pos - y_del_site
      y_site_med <- (y_site_bottom + y_site_top)/2
      # segments(x0=x_site_line + x_offset, y0=y_site_bottom, y1=y_site_top)
      # text(x_sites + x_offset, y=y_site_med, labels=sprintf("Site %s", site_num), adj=1)
      site_num <- site_num + 1

    }
    y_mirna_med <- (y_mirna_top + y_mirna_bottom)/2
    segments(x0=x_mirna_line + x_offset, y0=y_mirna_bottom - 0.01,
             y1=y_mirna_top)
    text(x_mirna + x_offset, y=y_mirna_med, labels=mirna, adj=c(1, 0))
    y_pos <- y_pos_start
    x_offset <- x_offset_set
  }
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

################################################################################
# FIGURE S3
################################################################################

# S2A-F_________________________________________________________________________
PlotBaekKds <- function(mirna, experiment="equilibrium", n_constant=5,
                        combined=TRUE, singleonly=TRUE, buffer=FALSE, width=7,
                        height=3.5, pdf.plot=FALSE) {
  pars.matrix <- SubfunctionCall(EquilPars, sitelist="baek")
  kds <- pars.matrix[grep(mirna, rownames(pars.matrix), invert=TRUE), ]
  # Make a data frame where each kd is associated with a sitetype
  # category.
  kds <- kds[order(kds$Mean),]
  kds <- kds[setdiff(rownames(kds), "None_Kd"), ]
  message(mirna)
  message(kds["6mer-m8", 2] / kds["7mer-m3.9_Kd", 2])
  rownames(kds) <- gsub("(.*)_Kd", rownames(kds), replace="\\1", perl=TRUE)
  rownames(kds) <- gsub("CDNST\\.", rownames(kds), replace="CDNST ", perl=TRUE)
  sites <- rownames(kds)
  SubfunctionCall(FigureSaveFile)
  xmin <- 1e-4
  xmax <- 3
  ymin <- 0
  ymax <- length(sites) + 2
  BlankPlot(log='x', inv='x')
  segments(1, ymin, 1, ymax, xpd=NA)
  y <- nrow(kds) - seq(nrow(kds)) + 1
  x_error <- list(kds$Upper_CI, kds$Lower_CI)
  Points(kds$Mean, y, x_error=x_error, col=kBaekColors[rownames(kds)])
  AddLogAxis(1, label="Relative Kd")
  title.xy <- GetPlotFractionalCoords(fx=0.95, fy=0.1, log='x', inv='x')
  mirna.split <- paste0(strsplit(mirna, split="-")[[1]][1:2], collapse="-")
  text(title.xy[1], title.xy[2], labels=mirna.split, adj=1)
  text((kds$Lower_CI)/1.3, y, labels=sites, adj=0)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

################################################################################
# FIGURE 3
################################################################################

# 3A____________________________________________________________________________
PlotPositionalKds <- function(experiment="equilibrium", n_constant=5,
                              buffer=FALSE, sitelist="centered11",
                              combined=TRUE, sigleonly=TRUE,
                              width=10, height=5, pdf.plot=FALSE) {
  # Get kds for all site-types of the mirna.
  SubfunctionCall(FigureSaveFile)
  xmin <- 1e-4
  xmax <- 10
  ymin <- 0
  ymax <- length(kPositionalSites)
  BlankPlot(log='x', inv='x')
  # xmax <- 10
  AddLogAxis(1, label="Relative Kd")
  segments(1, ymin, 1, ymax, xpd=NA)

  kMirnas <- kMirnas
  tick <- 0
  max_kds <- 0
  none_kd_offset_start <- 0.8
  none_kd_offset <- none_kd_offset_start
  df.global <- data.frame(mirnas=c(), color=c())
  sapply(kMirnas, function(mirna) {
    if (mirna == "miR-1") {
      buffer <- TRUE
      combined <- FALSE
    } else if (mirna == "miR-7-23nt") {
      experiment <- "equilibrium2_nb"
    }

    singleonly <- TRUE
    kds <- SubfunctionCall(EquilPars)[paste0(kPositionalSites, "_Kd"), ]
    rownames(kds) <- kPositionalSites
    kds <- kds[-nrow(kds), ]
    # message(mirna)
    # message("11mer-m3.13")
    # message(kds["6mer-m8", 2]/kds["11mer-m3.13", 2])
    # message("11mer-m4.14")
    # message(kds["6mer-m8", 2]/kds["11mer-m4.14", 2])
    if (sum(is.na(kds)) == 0) {
      max_kds <<- nrow(kds)
      max_kd_labels <<- rownames(kds)
    }
    color <- kMirnaColors[mirna]
    if (mirna == "miR-7-23nt") {
      mirna.str <- "miR-7"
    } else {
      mirna.str <- mirna
    }
    df.single <- data.frame(mirnas=c(mirna.str), colors=c(ConvertRColortoRGB(color)))
    df.global <<- rbind(df.global, df.single)
    ind_p <- c(1, 2, 3)
    y_p <- (nrow(kds) - seq(nrow(kds)) + 1)[ind_p]
    # y_p[length(y_p)] <- y_p[length(y_p)] + none_kd_offset
    ind_l <- seq(nrow(kds))[-ind_p]
    y_l <- (nrow(kds) - seq(nrow(kds)) + 1)[ind_l]
    Points(x=kds$Mean[ind_p], y=y_p, col=color)
    Points(x=kds$Mean[ind_l], y=y_l, col=color, ln.lwd=1,
          line=TRUE)
    lines(x=kds$Lower_CI[ind_l], y=y_l, col=color, lwd=0.5, lty=2)    
    lines(x=kds$Upper_CI[ind_l], y=y_l, col=color, lwd=0.5, lty=2)
    # Updated counter to move the None over
    none_kd_offset <<- none_kd_offset - 2*none_kd_offset_start/length(kMirnas)
  })
  file_out <- sprintf("aesthetic_tables_for_kathy/all_mirna_colors.txt")
  write.table(file_out, x=df.global, quote=FALSE, row.names=FALSE, col.names=TRUE,
              sep="\t")
  # Add text to plot:
  par(xpd=NA)
  xy_alt <- GetPlotFractionalCoords(fx=-0.01, fy=0.4, log='x', inv='x')
  text(x=xy_alt[1], y=max_kds - seq(max_kds) + 1, labels=max_kd_labels, adj=0,
       col="black")
  # Add legend to plot:
  xy <- GetPlotFractionalCoords(fx=0.75, fy=0.4, log='x', inv='x')
  kMirnas[length(kMirnas)] <- "miR-7"
  Legend(xy, legend=kMirnas, col=kMirnaColors)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

# 3B____________________________________________________________________________
Plot8merVs7merCombined <- function(experiment="equilibrium", n_constant=5,
                                   sitelist="resubmissionfinal", buffer=FALSE,
                                   singleonly=TRUE, dG.table=3, height=4,
                                   width=5, pdf.plot=FALSE) {
  dG.df <- read.table(file=paste0("canonical_sites_mfe_", dG.table, ".txt"))
  colnames(dG.df) <- gsub("(^.*)\\.(.*$)", colnames(dG.df), replace="\\1-\\2")
  R <- 1.987e-3 # in kcal K-1 mol-1
  T <- 310.15 # in K
  total_x <- c(1)
  total_y <- c(1)
  SubfunctionCall(FigureSaveFile)
  xmin <- 0
  xmax <- length(kMirnas)*4 + 7
  ymin <- -2.7
  ymax <- 0
  BlankPlot(inv="y")
  ymin <- -2.5
  mirna_labs <- kMirnas
  mirna_labs[length(mirna_labs)] <- "miR-7"
  AddLinearAxis(1, alt_lab=mirna_labs, alt_lab_pos=(seq(length(kMirnas)) - 1)*4 + 2, label="",
                angled=TRUE, noline=TRUE)

  AddLinearAxis(2, tick.space=0.5, label.space=1,
                label=expression(Delta*Delta*italic(G)~"(kcal/mol)"))

  x_start <- 0
  segments(x0=xmin, y0=-R*T*log(c(2, 10, 50)), x1=length(kMirnas)*4 + 0.5,
           lwd=0.25, xpd=NA)
  text(rep(length(kMirnas)*4 + 0.75, 3), -R*T*log(c(2, 10, 50)),
       labels=c("2-fold greater\nbinding affinity",
                "10-fold greater\nbinding affinity",
                "50-fold greater\nbinding affinity"),
       cex=0.8, adj=c(0, 0.5), xpd=NA)
  kds.all <- sapply(kMirnas, function(mirna) {
    if (mirna == "miR-7-23nt") {
      experiment <- "equilibrium2_nb"
      combined <- TRUE
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

# 3C____________________________________________________________________________
PlotSiteKdsVsSPS <- function(experiment="equilibrium", n_constant=5,
                             sitelist="resubmissionfinal", combined=TRUE, buffer=FALSE,
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


  print(summary(linmodel1))
  print(summary(linmodel2))

  coefs1 <- summary(linmodel1)$coefficients[2, 1:2]
  coefs2 <- summary(linmodel2)$coefficients[2, 1:2]
  print(coefs1)
  print(coefs2)


  RT_coef <- 1/(R*T)/log(10)

  prob_RT_1 <- 2*pt((coefs1[1] - RT_coef)/coefs1[2], df=length(kMirnas) - 2)
  prob_RT_2 <- 2*pt((coefs2[1] - RT_coef)/coefs2[2], df=length(kMirnas) - 2)

  message("probability 6mer is RTlnk:")
  message(prob_RT_1)
  message("probability 7mer-m8 is RTlnk:")
  message(prob_RT_2)


  predict1 <- predict(linmodel1, data.frame(x1_out=x.line), interval = 'confidence')
  polygon(c(x.line, rev(x.line)), c(10^predict1[,2], rev(10^predict1[, 3])), col=ConvertRColortoRGB("gray", alpha=0.3), border=NA)
  predict2 <- predict(linmodel2, data.frame(x2_out=x.line), interval = 'confidence')
  polygon(c(x.line, rev(x.line)), c(10^predict2[,2], rev(10^predict2[, 3])), col=ConvertRColortoRGB("gray", alpha=0.3), border=NA)
  lines(x.line, ypoints1)
  lines(x.line, ypoints2)

  sapply(1:ncol(dG.matrix), function(i_c) {
    x <- dG.matrix[, i_c]
    y <- Kd.matrix[, i_c]
    Points(x[order(x)], y[order(x)], pch=c(1, legend_pch),
           pt.lwd=c(1, 1), col=cols.mirnas[i_c])
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
         pch=c(19, 1, NA), pt.lwd=1, pt.cex=pt_cex_final, lty=c(1, 2, 1), lwd=c(1, 1, 1)*par()$lwd,
         seg.len=2.7,
         col=c("black", "black", "gray"), xjust=1, xpd=NA)
  leg.xy2 <- GetPlotFractionalCoords(fx=0.02, fy=1, log='y', inv='xy')
  legend(x=leg.xy2[1], y=leg.xy2[2], legend=rep("", length(kMirnas)), pch=1,
         pt.lwd=1, pt.cex=pt_cex_final, col=cols.mirnas, xpd=NA, bty="n")
  leg.xy2 <- GetPlotFractionalCoords(fx=0.06, fy=1, log='y', inv='xy')
  legend(x=leg.xy2[1], y=leg.xy2[2], legend=kMirnas, pch=19,
         pt.lwd=1, pt.cex=pt_cex_final, col=cols.mirnas, xpd=NA, bty="n")
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

# 3D-I__________________________________________________________________________
PlotSiteKdsVsRepression <- function(mirna, experiment="equilibrium",
                                    n_constant=5, sitelist="resubmissionfinal",
                                    single_site=FALSE, best_site=FALSE,
                                    bg_method=3, exrib=FALSE, combined=TRUE,
                                    buffer=FALSE, compcorrect=FALSE,
                                    n_cutoff=20, kd_cutoff=0.10, bulk=FALSE,
                                    new=TRUE, old=FALSE, cat_colors=FALSE,
                                    noncanon=FALSE, merge=FALSE, best=TRUE,
                                    threePseq=TRUE, height=5, width=5, xpos=20,
                                    ypos=20,
                                    pdf.plot=FALSE) {
  xmin <- 1e-4
  xmax <- 1
  ymin <- -1.4
  ymax <- 0.4

  kds <- SubfunctionCall(EquilPars)
  rownames(kds) <- gsub("_Kd", "", rownames(kds))
  removed_sites <- GetRemovedSites(kds)
  kds <- kds[setdiff(rownames(kds), removed_sites), ]
  sites <- grep(sprintf("(?:%s|None)", mirna), rownames(kds), perl=TRUE,
                  invert=TRUE, value=TRUE)
  adjusted <- FALSE
  l.xy <- GetPlotFractionalCoords(fx=-0.005, fy=-0.01, log='x', inv='x')
  yjust <- 0
  l.cex <- 1
  SubfunctionCall(FigureSaveFile)
  BlankPlot(log='x', inv='x', adjusted=adjusted)  
  AddLogAxis(1, label="Relative Kd", adj=TRUE)
  AddLinearAxis(2, tick.space=0.1, label.space=0.2, label="Fold change (log2)")
  fc <- SubfunctionCall(GetRepressionLinearModel)
  fc <- fc[setdiff(rownames(fc), removed_sites), ]
  fc_sites <- fc[, 1]
  fc_sem <- fc[, 2]
  fc_nsites <- fc[, 3]
  gray_col <- "gray50"
  segments(xmin, 0, xmax, 0, lwd=0.5, col=gray_col)
  kds <- kds[rownames(fc), ]

  sites <- intersect(rownames(kds), rownames(fc))
  kds_sites <- kds$Mean
  kds_sites_lci <- kds$Lower_CI
  kds_sites_uci <- kds$Upper_CI

  lm_fc_kd <- lm(fc_sites ~ log(kds_sites), weights=1/(fc[, 5] - fc[, 4])^2)
  lm_fc_kd <<- lm_fc_kd

  log_fc_kd <- lm(fc_sites ~ log10(kds_sites), weights=1/(fc[, 5] - fc[, 4])^2)
  log_fc_kd <<- log_fc_kd
  m_slope <- log_fc_kd$coefficients[2]

  slopes_lm[mirna] <<- m_slope


  m <- lm_fc_kd$coefficients[2]
  b <- lm_fc_kd$coefficients[1]

  # kSlopes[mirna] <<- m

  x.line <- exp(seq(log(xmin), -b/m, length.out=20))
  y.line <- m*log(x.line) + b
  lines(x.line, y.line, lwd=0.5)

  # put back cex=l.cex
  l.cex <- 1
  y.intersp <- 0.9
  colors.sites <- kSiteColors[sites]

  
  arrows(kds_sites, fc[, 4], kds_sites, fc[, 5],
         length=error_height_final*par()$cex, lwd=par()$cex, angle=90, code=3)

  arrows(kds_sites_lci, fc_sites, kds_sites_uci, fc_sites,
         length=error_height_final*par()$cex, lwd=par()$cex, angle=90, code=3)

  y_error <- list(fc[, 4], fc[, 5])
  x_error <- list(kds_sites_lci, kds_sites_uci)

  Points(kds_sites, fc_sites, x_error=x_error, y_error=y_error, col=colors.sites)
  sites <- sites[order(kds_sites)]
  colors.sites <- kSiteColors[sites]
  sites <- ConvertTtoUandMmtoX(sites)

  rownames(fc) <- ConvertTtoUandMmtoX(rownames(fc))
  legend_text_table <- cbind(sites, fc[sites, ][, 3])
  legend_text <- apply(legend_text_table, 1, function(row) {
    paste0(row[1], " (", row[2], ")")
  })
  if (mirna == "miR-7-23nt") {
    l1 <- Legend(l.xy, legend=legend_text[1:12],
                 cex=l.cex, col=colors.sites[1:12],
                 yjust=yjust, y.intersp=y.intersp, x.intersp=0.7)
    xy <- GetPlotFractionalCoords(fx=0.55, fy=-0.01, log='x', inv='x')
    l2 <- legend(x=xy[1], y=xy[2],
                 legend=legend_text[13:length(legend_text)],
                 bty="n", pch=legend_pch, cex=l.cex, pt.cex=pt_cex_final,
                 col=colors.sites[13:length(legend_text)], y.intersp=y.intersp,
                 x.intersp=0.7, xjust=0, yjust=yjust, xpd=NA)
  } else {
    l <- legend(x=l.xy[1], y=l.xy[2], legend=legend_text,
                bty="n", cex=l.cex, pch=legend_pch, col=colors.sites,
                pt.cex=pt_cex_final, xjust=0, yjust=yjust, xpd=NA,
                y.intersp=y.intersp, x.intersp=0.7)
  }

  xy <- GetPlotFractionalCoords(fx=0.9, fy=0.95, log='x', inv='x')
  mirna.split <- paste0(strsplit(mirna, split="-")[[1]][1:2], collapse="-")
  text(xy[1], xy[2], mirna.split, adj=1)
  # Tested the sem and the 95% confidence interval, and these work exactly the same.
  xy <- GetPlotFractionalCoords(fx=0.9, fy=0.90, log='x', inv='x')
  AddCorrelationToPlot(log(kds_sites), fc_sites, xpos=xy[1], ypos=xy[2],
                       rsquared=TRUE, adj=1, weights = 1/fc_sem^2)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

################################################################################
# FIGURE S3
################################################################################

# S3A,C,G,I,K___________________________________________________________________
PlotBulgeKds <- function(mirna, experiment="equilibrium", n_constant=5, 
                         combined=TRUE, singleonly=TRUE, buffer=FALSE,
                         height=3.5, width=7, pdf.plot=FALSE) {
  # Get kds for all site-types of the mirna.
  kds <- SubfunctionCall(EquilPars, sitelist="bulge")
  rownames(kds) <- ConvertTtoUandMmtoX(rownames(kds))
  kds <- kds[1:(nrow(kds) - 2), ]
  kds <- kds[setdiff(rownames(kds), "None_Kd"), ]
  rownames(kds) <- gsub("^(.*)_Kd$", rownames(kds), replace="\\1", perl=TRUE)
  bulge_matrix <- matrix(NaN, nrow=5, ncol=4,
                         dimnames=list(paste("8mer-b", seq(3, 7), sep=""), kNucs))
  bulge_matrix_full <- matrix(NaN, nrow=5, ncol=4,
                         dimnames=list(paste("8mer-b", seq(3, 7), sep=""), kNucs))
  SubfunctionCall(FigureSaveFile)
  xmin <- 1e-4
  xmax <- 7
  ymin <- 0
  ymax <- 11.5
  BlankPlot(log='x', inv='x')
  segments(1, ymin, 1, ymax, xpd=NA)

  xmin <- 1e-4
  # Now add the bulge kds to the matrix, appropriately:
  bulge.pivot <- MirnaTargetSequence(mirna, 6, 6)
  nuc.offset <- -3*c(0.15, 0.05, -0.05, -0.15)
  names(nuc.offset) <- kNucs
  N_seed <- length(kSeedSites)
  N_bulge <- nrow(bulge_matrix)
  y_pos_max <- N_seed + N_bulge + 1
  y_pos_seed <- y_pos_max - c(seq(N_seed-1), N_seed + N_bulge)
  y_pos_bulge <- y_pos_max - N_seed + 1 - seq(N_bulge)
  # y_pos_none <- 1
  y_pos <- c(y_pos_seed, rep(y_pos_bulge, each=4))
  eq.blg <- data.frame(kd=c(), min_pos=c(), max_pos=c(), nuc=c(),
                             alpha=c(), lwd=c(), stringsAsFactors=FALSE)
  # Make the eq.blg dataframe and bulge_matrix:
  sapply(grep("-b", rownames(kds)), function(index) {
    # Get the kd:
    kd <- kds$Mean[index]
    # Get the region after "-b" in the sitename:
    bulge_sub <- gsub("^.*-b(.*)$", rownames(kds)[index], replace="\\1",
                      perl=TRUE)
    # Get the nucleotide:
    nuc_ <- substr(bulge_sub, 1, 1)

    # Check for "(m.n)":
    if (substr(bulge_sub, 2, 2) == "(") {
      regex.key <- "^.*\\((.*)\\.(.*)\\).*$"
      min_pos <- as.integer(gsub(regex.key, bulge_sub, replace="\\1", perl=TRUE))
      max_pos <- as.integer(gsub(regex.key, bulge_sub, replace="\\2", perl=TRUE))
      # Add line connecting adjacent bulge sites:
      if (min_pos == 6 & nuc_ == bulge.pivot) {
        lwd_ <- 3
        alpha_ <- 1
      } else {
        lwd_ <- 1
        alpha_ <- 0.4
      }
      eq.blg <<- rbind(eq.blg, data.frame(kd=kd, min_pos=min_pos, max_pos=max_pos,
                                       nuc=nuc_, alpha=alpha_, lwd=lwd_,
                                       stringsAsFactors=FALSE))
      colnames(eq.blg) <<- c("kd", "min_pos", "max_pos", "nuc", "alpha", "lwd")
      # Add the entire range to the Kd matrix:
      sapply(seq(min_pos, max_pos), function(position) {
        bulge_matrix_full[paste0("8mer-b", position), nuc_] <<- kd
      })
    # If no range, add the single position:
    } else {
      position <- paste0("8mer-b", substr(bulge_sub, 2, nchar(bulge_sub)))
      bulge_matrix[position, nuc_] <<- kd
      bulge_matrix_full[position, nuc_] <<- kd
    }
  })
  bulge_matrix <<- bulge_matrix
  offset_matrix <- t(apply(bulge_matrix_full, 1, function(row) {
    out <- rep(NaN, 4)
    not.na <- row[!is.na(row)]
    offsets <- 0.2*seq(-1, 1, length=length(not.na))[order(not.na)]
    out[!is.na(row)] <- offsets
    out
  }))
  eq.blg <<- eq.blg
  colnames(offset_matrix) <- colnames(bulge_matrix)
  # Horizontal lines:
  segments(x0=apply(bulge_matrix_full, 1, min, na.rm=TRUE), y0=y_pos_bulge,
           x1=apply(bulge_matrix_full, 1, max, na.rm=TRUE), y1=y_pos_bulge,
           lty=2, lwd=0.5)
  GetOffset <- function(i_r, i_c) {
    offset_matrix[paste0("8mer-b", i_r), i_c]
  }
  offsets.min <- mapply(GetOffset, eq.blg$min_pos, eq.blg$nuc)*0
  offsets.max <- mapply(GetOffset, eq.blg$max_pos, eq.blg$nuc)*0
  segments(x0=eq.blg$kd, y0=y_pos_max - (eq.blg$min_pos + 3) + offsets.min,
           x1=eq.blg$kd, y1=y_pos_max - (eq.blg$max_pos + 3) + offsets.max,
           lwd=eq.blg$lwd, col=ConvertRColortoRGB(kRNucleotideColors[eq.blg$nuc]))
  y_offset_seed <- rep(0, N_seed)
  y_offset_bulge <- c(t(offset_matrix))
  y_offset_none <- 0
  y_offset <- c(y_offset_seed, y_offset_bulge, y_offset_none)*0

  seed_kds <- kds[kSeedSites, ]$Mean
  names(seed_kds) <- kSeedSites
  kds_list_final <- c(seed_kds, c(t(bulge_matrix)),
                      kds["None",]$Mean)
  pchs <- c(rep(19, N_seed), rep(1, N_bulge*length(kNucs)), 19)
  xmax <- 3
  AddLogAxis(1, label="Relative Kd")

  seed_kds <<- seed_kds
  points(kds_list_final, y_pos + y_offset,
         col=c(rep("gray80", 6), rep(kNucleotideColors, 5), "black"),
         pch=pchs, cex=1.2)
  joined <- rbindlist(lapply(1:nrow(eq.blg), function(index) {
    y <- y_pos_max - 3 - seq(eq.blg$min_pos[index], eq.blg$max_pos[index])
    x <- rep(eq.blg$kd[index], length(y))
    color <- kRNucleotideColors[as.character(eq.blg$nuc[index])]
    if (eq.blg$alpha[index]==1) {
      pch=19
    } else {
      pch=1
    }
    return(data.frame(x=x, y=y, color=color, pch=pch, stringsAsFactors=FALSE))
  }))
  par(xpd=NA)

  points(joined$x, joined$y, col=joined$color, cex=1.2, pch=joined$pch)
  # xy = GetPlotFractionalCoords(fx=1, fy=0.9, log='x', inv='x')
  xmax <- 3
  xy <- GetPlotFractionalCoords(fx=0.025, fy=1.025, log='x', inv='x')

  mirna.split <- paste0(strsplit(mirna, split="-")[[1]][1:2], collapse="-")
  text(xy[1], xy[2], mirna.split, adj=c(0, 0))
  text(20, unique(y_pos),
       labels=c(names(seed_kds), rownames(bulge_matrix)), adj=0)
  # Legend baseline:
  if (pdf.plot == "S4.B") {
    f.l <- 0.535
    f.y <- 0.56
    leg.xy <- GetPlotFractionalCoords(fx=f.l, fy=f.y, log='x', inv='x')
    legend(x=leg.xy[1], y=leg.xy[2], legend=kNucs, col=kNucleotideColors, bty="n",
           pch = 19, xjust=0, yjust=1)
    adj1 <- 0.115
    leg.xy <- GetPlotFractionalCoords(fx=f.l + adj1, fy=f.y + 0.065, log='x', inv='x')
    leg_adj <- 0.8
    par(lheight = 0.8*par()$lheight)
    legend(x=leg.xy[1], y=leg.xy[2],
           legend=c("Unambiguous\nbulged position", "Ambiguous\nbulged position",
                    "Pivot site\n"),
           adj=c(0, leg_adj), xjust=0, yjust=1, lwd=c(NA, 1, 3), pt.lwd=1,
           pt.cex=1.2, pch=c(1, 1, 19), bty="n", y.intersp=1.6, x.intersp=0.4)
    par(lheight = par()$lheight/0.8)
  }


  if (class(pdf.plot) == "character") {
    dev.off()
  }
}


# S3B,D,F,H,J,L_________________________________________________________________
PlotDelKds <- function(mirna, experiment="equilibrium", n_constant=5,
                       combined=TRUE, singleonly=TRUE, buffer=FALSE,
                       sitelist="del", height=3.5, width=7, pdf.plot=FALSE) {
  # Get kds for all site-types of the mirna.
  kds <- SubfunctionCall(EquilPars)
  # Convert the Kds to non Kd names:
  rownames(kds) <- gsub("^(.*)_Kd$", rownames(kds), replace="\\1", perl=TRUE)
  # kds <- kds[setdiff(rownames(kds), "None_Kd"), ]
  # Make matrix with each deleted position (i.e. with redundant Kds):
  kds_full <- matrix(NaN, nrow=12, ncol=ncol(kds),
                         dimnames=list(c(kSeedSites[-length(kSeedSites)],
                                         paste("8mer-d", seq(2, 7), sep=""),
                                         "6mer-A1"), colnames(kds)))
  eq_del <- data.frame(kds=c(), min_pos=c(), max_pos=c(), stringsAsFactors=FALSE)
  inds_fill <- intersect(rownames(kds), rownames(kds_full))
  for (name in colnames(kds)) {
    kds_full[inds_fill, name] <- kds[inds_fill, name]  
  }
    # Split up the "b(N.M)" type sites:
  sapply(grep("-d\\(.*\\..*\\)$", rownames(kds), perl=TRUE), function(index) {
    regex <- "^8mer\\-d\\((.*)\\.(.*)\\)$"
    start_pos <- gsub(regex, rownames(kds)[index], replace="\\1", perl=TRUE)
    end_pos <- gsub(regex, rownames(kds)[index], replace="\\2", perl=TRUE)
    row_labs <- paste0("8mer-d", seq(start_pos, end_pos))
    # Add this information to global variable, to make the thick lines.
    eq_del <<- rbind(eq_del,
                     data.frame(kds=kds$Mean[index],
                                min_pos=which(rownames(kds_full)==row_labs[1]),
                                max_pos=which(rownames(kds_full)==row_labs[length(row_labs)]),
                                stringsAsFactors=FALSE))
    # Update the global matrices:
    for (name in colnames(kds)) {
      kds_full[row_labs, name] <<- kds[index, name]
    }
  })
  SubfunctionCall(FigureSaveFile)
  xmin <- 1e-4
  xmax <- 7
  ymin <- 0
  ymax <- 12.5
  BlankPlot(log='x', inv='x')
  xmax <- 3
  xmin <- 1e-4
  segments(1, ymin, 1, ymax, xpd=NA)

  AddLogAxis(1, label="Relative Kd")

  y_vals = nrow(kds_full) + 1 - seq(nrow(kds_full)) # The y values:
  colors <- rep("gray80", nrow(kds_full))
  pchs <- rep(19, nrow(kds_full))
  colors[grep("8mer-d", rownames(kds_full))] <- "black"
  pchs[grep("8mer-d", rownames(kds_full))] <- 1
  # colors[nrow(kds_full)] <- "black"
  apply(eq_del, 1, function(row) {
    kd <- as.numeric(row[1])
    ind_l <- as.numeric(row[2])
    ind_h <- as.numeric(row[3])
    lines(x=rep(kd, 2), y=y_vals[c(ind_l, ind_h)],
          col="black", lwd=1)
  })
  points(x=kds_full[, "Mean"], y=y_vals, col=colors, pch=pchs, cex=1.2)
  par(xpd=NA)

  text(x=20, y=y_vals, labels=rownames(kds_full), adj=0,
       col="black")
  xy <- GetPlotFractionalCoords(fx=0.025, fy=1.025, log='x', inv='x')

  mirna.split <- paste0(strsplit(mirna, split="-")[[1]][1:2], collapse="-")
  text(xy[1], xy[2], mirna.split, adj=c(0, 0))

  # Add legend to plot:
  if (pdf.plot == "S4.C") {
    # Legend baseline:
    # f.l <- 0.5
    # f.y <- 0.55
    # adj1 <- 0.070
    # leg.xy <- GetPlotFractionalCoords(fx=f.l + adj1, fy=f.y, log='x', inv='x')
    # leg_adj <- 0.5
    # legend(x=leg.xy[1], y=leg.xy[2],
    #        legend=c("Unambiguous", "deleted position", "Ambiguous",
    #                 "deleted position"),
    #        adj=c(0, leg_adj), xjust=0, yjust=1, lwd=c(0, 0, 1, 0),
    #        col=c(NA, NA, "black", 
    #              NA),
    #        bty="n")
    # adj2 <- 0.0925 - 0.060
    # leg.xy <- GetPlotFractionalCoords(fx=f.l + adj1 + adj2, fy=f.y, log='x', inv='x')
    # legend(x=leg.xy[1], y=leg.xy[2], legend=rep("", 4),
    #        col=c("black", NA, "black", NA), bty="n",
    #        xjust=0, yjust=1, pt.cex=1.2, adj=c(0, leg_adj),
    #        pch = 1)    

    f.l <- 0.535
    f.y <- 0.56
    # leg.xy <- GetPlotFractionalCoords(fx=f.l, fy=f.y, log='x', inv='x')
    # legend(x=leg.xy[1], y=leg.xy[2], legend=kNucs, col=kNucleotideColors, bty="n",
    #        pch = 19, xjust=0, yjust=1)
    adj1 <- 0.115
    leg.xy <- GetPlotFractionalCoords(fx=f.l + adj1, fy=f.y + 0.065, log='x', inv='x')
    leg_adj <- 0.8
    par(lheight = 0.8*par()$lheight)
    legend(x=leg.xy[1], y=leg.xy[2],
           legend=c("Unambiguous\ndeleted position", "Ambiguous\ndeleted position"),
           adj=c(0, leg_adj), xjust=0, yjust=1, lwd=c(NA, 1), pt.lwd=1,
           pt.cex=1.2, pch=1, bty="n", y.intersp=1.6, x.intersp=0.4)
    par(lheight = par()$lheight/0.8)
  }
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}


################################################################################
# FIGURE 4
################################################################################
# 4A____________________________________________________________________________
PlotSiteFlankEnrichments <- function(mirna,  site, experiment="equilibrium",
                                     n_constant=5, sitelist="resubmissionfinal",
                                     combined=TRUE, combined_site=TRUE,
                                     buffer=FALSE, bgoff=FALSE,
                                     L=FALSE, height=5, width=6,
                                     pdf.plot=FALSE) {
  sXc <- SubfunctionCall(SitesAndSingleFlanksXCounts)
  pars.matrix <- SubfunctionCall(GetFlankKds)
  pars.matrix <- pars.matrix[grep("mer|None|^AGO|^bg", rownames(pars.matrix)), ]
  pars <- log10(pars.matrix$Mean)
  names(pars) <- rownames(pars.matrix)
  sXc <- sXc[grep("mer|None", rownames(sXc), perl=TRUE), ]
  names(pars)[nrow(sXc) + 1] <- "bg"
  names(pars)[nrow(sXc) + 2] <- "AGO"
  data <- GetDataEquil(sXc)
  l <- SubfunctionCall(GetInputEquil)
  if (L) {
    l <- l/sum(l) * as.numeric(L)
  }
  # Ago dilutio in the data:
  A.dil.data <- sapply(colnames(data), as.numeric)
  A.stock.measured <- kAgoStock[mirna, "equilibrium"]
  A.stock.pars <- 10^pars["AGO"]
  pM_from_dil <- A.stock.measured*1000/100
  A.pM.data <- A.dil.data*pM_from_dil
  xmin <- signif(min(A.pM.data)/(10^0.25), 1)
  xmax <- signif(max(A.pM.data)*(10^0.75), 1)
  A.dil.model <- exp(seq(log(xmin/pM_from_dil),
                         log(xmin/pM_from_dil*10^(5/2)),
                     length=10))

  model <- SubfunctionCall(EquilSingleSiteModelFreq, A.dil=A.dil.model)
  A.pM.model <-A.dil.model*pM_from_dil
  data.R <- EquilEnrichments(data, l)
  data.R <<- data.R
  model.R <- EquilEnrichments(model, l)
  # Set up the plotting limits
  ymin <- 0.2
  ymax <- 500
  SubfunctionCall(FigureSaveFile)
  BlankPlot(log="xy", adjusted=TRUE)
  # Generate tickmarks for axis.
  xy <- GetPlotFractionalCoords(fx=0.05, fy=0.95, log='xy')
  mirna.trim <- paste0(strsplit(mirna, split = "-")[[1]][1:2], collapse="-")
  text(xy[1], xy[2], labels=mirna.trim, adj=c(0, 1))
  kds <- pars[paste0(rownames(data), "_Kd")]
  names(kds) <- rownames(data)
  colors <- CombinedSiteAndFlankColors(sXc)
  names(colors) <- names(kds)
  bg_colors <- which(colors == "gray")
  kds <- c(kds[bg_colors], kds[-bg_colors])
  colors_global <<- colors
  sapply(names(kds), function(site) {
    Points(A.pM.data, data.R[site, ], col=colors[site])
    lines(A.pM.model, model.R[site, ], col=colors[site])      
  })
  AddLogAxis(1, label=AGO_mir_label, adj=TRUE)
  AddLogAxis(2, label="Enrichment")

  Dinucs <- c("CCCC", "CCCA", "AACC", "AAAC", "AAAA")
  Labels <- sapply(seq(0, 4), function(i) {
    eval(substitute(expression((A/U)[x](G/C)[y]), list(x=i, y=4-i)))})
  xy <- GetPlotFractionalCoords(0.85, 0.5, log='xy')
  legend(x=xy[1], y=xy[2], legend=rev(Labels),
         col=GetColorFunction(rev(Dinucs)), bty="n", pch=15, pt.cex=2, xpd=NA)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}


# 4B____________________________________________________________________________
PlotSiteFlankKds <- function(mirna, experiment="equilibrium", n_constant=5,
                             sitelist="resubmissionfinal", combined=TRUE,
                             combined_sites=TRUE, singleonly=TRUE, buffer=FALSE,
                             xmin=log10(1e-7),
                             adjusted_height=FALSE, plot.nconstant=FALSE,
                             height=5, width=8, xpos=20, ypos=20, pdf.plot=FALSE) {

  pars.matrix <- SubfunctionCall(EquilPars, combined=combined_sites)
  site.kds <- pars.matrix[1:(nrow(pars.matrix) - 3), ]
  sites <- gsub("^(.*)_Kd", rownames(site.kds), replace="\\1")
  sXc <- SubfunctionCall(SitesXCounts)
  rownames(site.kds) <- sites
  removed_sites <- GetRemovedSites(sXc)
  print(removed_sites)
  sites_AA <- grep("^AA-", rownames(sXc), perl=TRUE, value=TRUE)
  # if (length(sites_AA) != 0) {
  #   inds_AA <- grep("^AA-", rownames(sXc), perl=TRUE)
  #   sites_AA_base <- gsub("^(AA-)(.*)$", sites_AA, replace="\\2", perl=TRUE)
  #   sites_AA_keep <- sites_AA[sites_AA_base %in% rownames(sXc)]
  #   sites_AA_remove <- setdiff(sites_AA, sites_AA_keep)
  #   removed_sites <- c(sites_AA_remove, removed_sites)
  # }
  sites <- setdiff(sites, removed_sites)
  site.kds <- site.kds[sites, ]
  sites <- sites[order(site.kds$Mean)]
  site.kds <- site.kds[order(site.kds$Mean), ]
  # print(order(site.kds$Mean))
  # print(site.kds)
  # print(sites[order(site.kds$Mean)])
  # break
  # Get kds for all site-types of the mirna.
  # # Use the 8mer flanks to just get the flank strings.
  # Pre-allocate the matrix with the flanking kds.
  tick <<- 1
  flank.kds <- sapply(sites, GetFullMirnaSiteFlankKds, mirna=mirna,
                      experiment=experiment, n_constant=n_constant,
                      sitelist=sitelist, combined=combined, buffer=buffer)
  flank.kds_CI <- sapply(sites, GetFullMirnaSiteFlankKds_CI, mirna=mirna,
                      experiment=experiment, n_constant=n_constant,
                      sitelist=sitelist, combined=combined, buffer=buffer)
  flank.kds_uCI <- flank.kds_CI[1:nrow(flank.kds), ]
  flank.kds_lCI <- flank.kds_CI[(nrow(flank.kds) + 1):(2*nrow(flank.kds)), ]
  rownames(flank.kds_uCI) <- rownames(flank.kds)
  rownames(flank.kds_lCI) <- rownames(flank.kds)
  flank.kds_CI_range <- flank.kds_uCI/flank.kds_lCI
  flank.kds_median_range <- apply(log10(flank.kds_CI_range), 2, median, na.rm=TRUE)
  # Removes site types for which there are literally zero flank Kds.


  print(site.kds)
  site.kds <- site.kds[which(colSums(is.na(flank.kds)) != 256),]
  flank.kds <- flank.kds[, rownames(site.kds)]
  # flank.kds <- flank.kds[, order(site.kds$Mean)]
  data.sites <- rep(colnames(flank.kds), each=nrow(flank.kds))
  data.ranks <- rep(ncol(flank.kds) - seq(ncol(flank.kds)) + 1, each=nrow(flank.kds))
  data.colors <- rep(GetColorFunction(kFlanks, alpha=1), ncol(flank.kds))
  flanks.df <- data.frame(kds=log10(c(flank.kds)),
                          rank=as.numeric(data.ranks),
                          sites=data.sites,
                          cols=data.colors,
                          stringsAsFactors=FALSE)
  if (pdf.plot == "4.B") {
    xmin <- log10(1e-5)
  }
  xmax <- log10(7)
  ymin <- 0
  ymax <- length(sites) + 0.5
  if (adjusted_height) {
    height <- (ymax*0.8 + 5)/4
  } else {
    height <- 5
  }
  SubfunctionCall(FigureSaveFile)
  xpd=NA
  ranks <- unique(flanks.df$rank)
  names(ranks) <- unique(flanks.df$site)
  boxplot(kds ~ rank,
          data       = flanks.df,
          axes       = FALSE,
          xaxt       = "n",
          horizontal = TRUE,
          outline    = FALSE,
          xlim       = c(ymin, ymax),
          ylim       = rev(c(xmin, xmax)),
          ann        = FALSE)
  y <- nrow(site.kds) - seq(nrow(site.kds)) + 1
  segments(0, ymin, 0, max(data.ranks) - 0.5, xpd=NA)
  # xmin <- log10(1e-5)
  AddLogAxis(1, label="Relative Kd", boxplot=TRUE)
  mirna.split <- paste0(strsplit(mirna, split="-")[[1]][1:2], collapse="-")
  xy <- GetPlotFractionalCoords(0.025, 0.95, inv='x')
  text(xy[1], max(data.ranks), labels=mirna.split, adj=c(0, 0), xpd=NA)  
  par(xpd=NA)
  beeswarm(kds ~ rank,
           data        = flanks.df,
           add         = TRUE,
           method      = "swarm",
           corral      = "random",
           corralWidth = 0.5,
           pch         = 1,
           lwd         =1.2,
           cex         = 0.8,
           horizontal  = TRUE,  
           pwcol       = cols,
           axes        = FALSE,
           xpd=NA)
  ymin=0.0001
  flank_mins <- apply(flank.kds, 2 ,function(col) {min(col[!is.na(col)])})
  text(x=log10(flank_mins) - 0.1,
       y=unique(flanks.df$rank),
       labels=ConvertTtoUandMmtoX(unique(flanks.df$sites)),
       adj=0, col= "black", xpd=NA)
  if (pdf.plot == "4.B") {
    bar_pos <- -4.5
  } else {
    bar_pos <- -6.5
  }
  arrows(bar_pos - flank.kds_median_range[names(ranks)]/2, ranks,
         bar_pos + flank.kds_median_range[names(ranks)]/2, ranks,
         length=0.02, angle=90, code=3, xpd=NA)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

# 4C-left_______________________________________________________________________
PlotFlankLinModel <- function (experiment="equilibrium", n_constant=5,
                               sitelist="resubmissionfinal", leaveoneout=TRUE,
                               new_way=TRUE, height=5,
                               width=5, pdf.plot=FALSE) {
  # Extract the flanking dinucleotides:
  flank.lm <- SubfunctionCall(GetFlankLinearModel)
  flank.lm <<- flank.lm
  # Data analysis for the 1st plot:
  # First linear model; splitting up all flanks by 5p and 3p sequence.
  SubfunctionCall(FigureSaveFile)
  # 1st Plot:
  if (new_way) {
    xmin <- 4e-2
    xmax <- 70    
  } else {
    xmin <- 1e-5
    xmax <- 10        
  }
  ymin <- xmin
  ymax <- xmax
  site_mod <- lm(logkd ~ site*mirna, data=flank.lm$model)
  data_temp <- flank.lm$model
  data <- exp(flank.lm$model$logkd)
  data_new <- data/exp(predict(site_mod, data_temp))
  data_out <<- data
  list_out <<- list()
  if (leaveoneout) {
    model <- exp(unlist(sapply(kMirnas, function(mirna) {
      out <- SubfunctionCall(LeaveOneOutFlankModel)
      list_out[[mirna]] <<- out
      out
    })))
  } else {
    model <- exp(flank.lm$fitted.values)  
  }

  flank.lm$model$logkd_mod <- log(model)
  flank.lm_corective <- lm(logkd_mod ~ mirna*site, data=flank.lm$model)

  model_new <- model/exp(flank.lm_corective$fitted.values)
  flanks <- paste(as.character(flank.lm$model$f5p1),
                  as.character(flank.lm$model$f5p2),
                  as.character(flank.lm$model$f3p1),
                  as.character(flank.lm$model$f3p2), sep = "")
  BlankPlot(log='xy', inv='xy')
  AddLogAxis(1, label="Predicted relative Kd")
  AddLogAxis(2, label="Observed relative Kd")
  # inds <- grep("miR-23nt", names(model), invert=TRUE)
  xy <- GetPlotFractionalCoords(fx=0.95, fy=0.05, log='xy')
  abline(0, 1, lty=2)



  if (new_way) {
    model_new <<- model_new
    data_new <<- data_new
    points(x=model_new, y=data_new, col=GetColorFunction(flanks, alpha=0.4))  
    AddCorrelationToPlot(x=log(model_new), y=log(data_new), xpos=xy[1],
                         ypos=xy[2], rsquared=TRUE)
  } else {
    model <<- model
    data <<- data
    points(x=model, y=data, col=GetColorFunction(flanks, alpha=0.2))  
    AddCorrelationToPlot(x=log(model), y=log(data), xpos=xy[1],
                           ypos=xy[2], rsquared=TRUE)
  }
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}


# 4C-right______________________________________________________________________
PlotFlankLinModelCoefficients <- function(experiment="equilibrium",
                                          n_constant=5,
                                          sitelist="resubmissionfinal", width=4,
                                          height=5, xpos=20, ypos=20,
                                          pdf.plot=FALSE) {
  coefs <- SubfunctionCall(GetFlankLMCoefs)
  coefs <<- coefs
  ymin <- -0.5
  ymax <- 0.5
  xmin <- 0
  xmax <- 24
  R <- 1.987e-3 # in kcal K-1 mol-1
  T <- 310.15 # in K

  SubfunctionCall(FigureSaveFile)
  BlankPlot(inv='y')
  axis(1, at=c(4, 8, 12, 16)-2, labels=c("5p1", "5p2", "3p1", "3p2"),
      pos = ymax, lwd=0)

  # Add y-axis:
  AddLinearAxis(2, tick.space=0.1, label.space= 0.1,
                label=expression(Delta*Delta*italic(G)*" (kcal/mol)"))

  x <- seq(16) - 0.5
  lins <- coefs$lin
  lins_error_upper <- coefs$lin_error_upper
  lins_error_lower <- coefs$lin_error_lower
  int5p <- cbind(c(0, 0, 0, 0), rbind(c(0, 0, 0), coefs$int.5p))
  int3p <- cbind(c(0, 0, 0, 0), rbind(c(0, 0, 0), coefs$int.3p))
  nucs <- c("A", "C", "G", "T")
  rownames(int5p) <- nucs
  colnames(int5p) <- nucs
  rownames(int3p) <- nucs
  colnames(int3p) <- nucs
  nt <- names(kNucleotideColors)
  y_out1  <<- c(t(t(lins) - colMeans(lins))[nt, ], int5p[nt, nt] - mean(int5p),
         int3p[nt, nt] - mean(int3p))*R*T

  y_out2  <<- c(t(t(lins) - colMeans(lins))[nt, ], int5p[nt, nt],
         int3p[nt, nt])*R*T
  lins_global <<- lins*R*T
  y <- R*T*t(t(lins) - colMeans(lins))[nt, ]
  lins_post_global <<- y

  int5p <- R*T*int5p
  int3p <- R*T*int3p

  # int5p <<- int5p
  # int3p <<- int3p

  int5p_norm <<- int5p - mean(int5p)
  int3p_norm <<- int3p - mean(int3p)
  y_error_upper <- R*T*t(t(lins_error_upper) - colMeans(lins))[nt, ]
  y_error_lower <- R*T*t(t(lins_error_lower) - colMeans(lins))[nt, ]
  # Define boundaries for plotting:
  y.b <- sapply(y, max, i=0)
  y.t <- sapply(y, min, i=0)
  segments(x0=xmin,y0=R*T*log(c(1/2, 1, 2)), x1=16, lwd=0.25)
  text(x=c(16.5, 16.5), y=R*T*log(c(1/2, 2)), cex=1,
       labels=c("2-fold greater\nbinding affinity",
                "2-fold weaker\nbinding affinity"), adj=c(0, 0.5), xpd=NA)

  abline(v=seq(3)*4, lty=2, lwd=0.5)
  arrows(x0=x[-c(1, 5, 9, 13)], y0=y_error_upper[-1, ], x1=x[-c(1, 5, 9, 13)],
         y1=y_error_lower[-1, ], length=0.05*par()$cex, angle=90, code=3, xpd=NA)
  rect(xleft=x - 0.5, ybottom=y.b, xright=x + 0.5, ytop=y.t,
       col=kNucleotideColors, border=NA)

  xy <- GetPlotFractionalCoords(fx=0.7, fy=0.8, inv="y")
  legend(x=xy[1], y=xy[2], legend=ConvertTtoU(nt), col=kNucleotideColors,
         pt.cex=2, pch=15, bty="n", xpd=NA)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}


# 4D____________________________________________________________________________
PlotStructureVsFlankingKds <- function(mirna, site, condition="I_combined",
                                       combined=TRUE, buffer=FALSE,
                                       experiment="equilibrium", n_constant=5,
                                       sitelist="resubmissionfinal", mir.start=1, mir.stop=14,
                                       absolute=TRUE, noconstant=FALSE,
                                       height=5, width=5, xpos=20, ypos=20,
                                       pdf.plot=FALSE) {
  # Window size for normalizing the pl_fold with the accessibility:
  win <- 1/(mir.stop - mir.start + 1)
  # Get the flanking kds:
  flank.pars <- SubfunctionCall(GetFlankKds)
  flank.pars <- flank.pars[grep(sprintf("^%s\\|", site), rownames(flank.pars)), ]
  kds <- flank.pars$Mean

  # Remove the "." in the flank names:
  names(kds) <- gsub(sprintf("^%s\\|(.*)\\.(.*)_Kd$", site), rownames(flank.pars), replace="\\1\\2")
  kds_4D <<- kds
  # Get the flank structural dataframe:
  if (condition == "I_combined") {
    data_1 <- SubfunctionCall(GetPairingFlankData, condition="I", alt_mir_exp_cond="let-7a_equilibrium_I")
    data_2 <- SubfunctionCall(GetPairingFlankData, condition="I", alt_mir_exp_cond="miR-124_equilibrium_I")
    data_3 <- SubfunctionCall(GetPairingFlankData, condition="I", alt_mir_exp_cond="miR-1_kin_pilot_I_TGT")
    data_4 <- SubfunctionCall(GetPairingFlankData, condition="I")
    data <- rbind(data_1, data_2, data_3, data_4)
  } else {
    data <- SubfunctionCall(GetPairingFlankData)  
  }
  data <<- data
  p_access <- c(by(data, data$flank, function(x) exp(mean(log(x[["plfold"]]))*15)))
  # p_access <<- p_access
  xmin <- 1e-6
  xmax <- 1e-1
  ymin <- 0.0003
  ymax <- 0.2
  SubfunctionCall(FigureSaveFile)
  BlankPlot(log='xy', inv='y')
  # Plot the line showing the linear fit of the relationshp:
  x <- p_access
  y <- kds
  fit <- lm(log10(y) ~ log10(x))
  m <- fit$coefficients[2]
  b <- fit$coefficients[1]
  x_line <- 10^seq(log10(xmin), log10(xmax), length=20)
  y_line <- 10^(m*log10(x_line) + b)
  lines(x_line, y_line, lty = 2,lwd = 0.5)
  # Plot the points:
  points(x=x, y=y[names(x)],
         col=GetColorFunction(names(x)))
  AddLogAxis(1, label="Mean accessibility score")
  AddLogAxis(2, label="Relative Kd", adj=TRUE)
  print("FLANKING DINUCLEOTIDE KD P VALUE")
  print(sprintf("r^2: %.4f", cor.test(log(x), log(y))$estimate^2))
  print(sprintf("p_value: %.6f", cor.test(log(x), log(y))$p.value))
  xy <- GetPlotFractionalCoords(fx=0.05, fy=0.95, log='xy', inv='y')
  AddCorrelationToPlot(x=log(y[names(y)]), y=-log(x), xpos=xy[1],
                       ypos=xy[2], rsquared=TRUE)
  if (condition == "I_combined") {
    condition <- "I"
  }
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

# 4E____________________________________________________________________________
PlotAllSamplePlFlanks <- function(mirna, site, condition, mir.start=1,
                                  mir.stop=14, experiment="equilibrium",
                                  n_constant=5, sitelist="resubmissionfinal",
                                  combined=TRUE, buffer=TRUE, noconstant=FALSE,
                                  absolute=TRUE, depth=5000, p_score="plfold",
                                  pdf.plot=FALSE) {
  height <- 5
  width <- 5
  # Window size for normalizing the pl_fold with the accessibility:
  if (combined) {
    data_1 <- SubfunctionCall(GetPairingFlankData, condition="I", alt_mir_exp_cond="let-7a_equilibrium_I")
    data_2 <- SubfunctionCall(GetPairingFlankData, condition="I", alt_mir_exp_cond="miR-124_equilibrium_I")
    data_3 <- SubfunctionCall(GetPairingFlankData, condition="I", alt_mir_exp_cond="miR-1_kin_pilot_I_TGT")
    data_4 <- SubfunctionCall(GetPairingFlankData, condition="I")
    data.I <- rbind(data_1, data_2, data_3, data_4)
  } else {
    data.I <- SubfunctionCall(GetPairingFlankData)  
  }
  data.A <- SubfunctionCall(GetPairingFlankData)
  inds_AU <- c(sapply(seq(50), function(i) {
    SampleByDinucleotideEnrichment(data.I, data.A, nrow(data.A))
  }))
  GetPercentTAEffect <- function(condition) {
    data.A <- SubfunctionCall(GetPairingFlankData)
    mean.I <- mean(log(data.I[[p_score]]))
    mean.A <- mean(log(data.A[[p_score]]))
    inds <- c(sapply(seq(50), function(i) {
      SampleByDinucleotideEnrichment(data.I, data.A, nrow(data.A))
    }))
    mean.I.sample <- mean(log(data.I[[p_score]][inds]))
    # test_inside <<- c(mean.I, mean.A, mean.I.sample)
    out <- (mean.I.sample - mean.I)/(mean.A - mean.I)
    message("GetPercentTAEFfect:")
    print(out)
    return(out)
  }
  #Plot the probabilities:
  xmin <- 1e-8
  xmax <- 1
  ymin <- 0
  ymax <- 1
  # if (class(pdf.plot) == "character") pdf.plot <- paste0(pdf.plot, "_left")
  SubfunctionCall(FigureSaveFile)
  BlankPlot(log='x')
  # Make the ECDF for the input, 0.4 Ago, and the resampled:
  x_ecdf <- 10^seq(-8, 0, by=0.01)

  ecdf.I <- ecdf(data.I[[p_score]]^(mir.stop - mir.start + 1))
  ecdf.A <- ecdf(data.A[[p_score]]^(mir.stop - mir.start + 1))
  ecdf.I.sample <- ecdf(data.I[[p_score]][inds_AU]^(mir.stop - mir.start + 1))
  lines(x_ecdf, ecdf.I(x_ecdf), lwd=1)
  lines(x_ecdf, ecdf.A(x_ecdf),
        col=kEquilCondColors[as.character(condition)])
  lines(x_ecdf, ecdf.I.sample(x_ecdf),
        col=kEquilCondColors[as.character(condition)], lty="23")
  # Add the points (3 in total) showing the mean value:
  x_GeoMeans <- exp(c(
    mean(log(data.I[[p_score]])),
    mean(log(data.A[[p_score]])),
    mean(log(data.I[[p_score]][inds_AU]))
   ))
  test <- log(x_GeoMeans)
  test_out <<- test
  alt_check <- (test[3] - test[1])/(test[2] - test[1])

  y_GeoMeans <- c(ecdf.I(x_GeoMeans[1]^(mir.stop - mir.start + 1)),
                  ecdf.A(x_GeoMeans[2]^(mir.stop - mir.start + 1)),
                  ecdf.I.sample(x_GeoMeans[3]^(mir.stop - mir.start + 1)))
  points(x=x_GeoMeans^(mir.stop - mir.start + 1), y=y_GeoMeans, col=c("black", "red", "red"),
         pch=c(19, 19, 20))
  AddLogAxis(1, label="Accessibility score", maglabel=2)
  AddLinearAxis(2, tick.space=0.2, label.space=0.2, label="Cumulative fraction")
  xy <- GetPlotFractionalCoords(fx=0.95, fy=0.05, log='x')
  text(xy[1], xy[2],
       label=paste0("Effect: ",round(GetPercentTAEffect(condition)*100,1),"%"),
       adj=1)
  # Add a legend to the plot:
  xy <- GetPlotFractionalCoords(fx=0, fy=1, log='x')
  legend(x=xy[1], y=xy[2],
         legend=c("Input library", AGO_mir1_label_no_conc,
                  "Input matched for flanking\ndinucleotide composition"),
         col=c("black", kEquilCondColors[as.character(condition)],
                 kEquilCondColors[as.character(condition)]),
         lty=c(1, 1, 2), bty="n")
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

################################################################################
# SUPPLEMENTARY FIGURE 5
################################################################################

# S5A___________________________________________________________________________
PlotReporterAssaySchematicLetters <- function(buffer=FALSE, pdf.plot=FALSE,
                                              height=1, width=8.1) {
  # Get the two signal data sets:
  counts_full_1 <- read.table("ReporterScreen/final_order_twist.fa", sep="\t",
                              stringsAsFactors=FALSE)
  variants <- counts_full_1[c(2, 4, 368),]
  variants <- c(variants[1:2], "...", variants[3])
  names(variants) <- c("Variant 1", "Variant 2", "Space", "Variant 184")

  SubfunctionCall(FigureSaveFile)
  xmin <- 0
  xmax <- 1
  ymin <- 0
  ymax <- 1
  par(mar=c(0, 0, 0, 0))
  BlankPlot()
  # Define the left and right hand sides of the text string
  x_l <- GetPlotFractionalCoords(0.15, 0)[1]
  x_r <- GetPlotFractionalCoords(0.93, 0)[1]

  y_t <- GetPlotFractionalCoords(0.05, 0.80)[2]
  y_b <- GetPlotFractionalCoords(0.05, 0.04)[2]

  y_list <- seq(y_t, y_b, length.out=4)
  x_list <-  seq(x_l, x_r)

  flank_num <- 20
  seed_pos <- 87
  seed_stop <- 94
  x_list <- seq(x_l, x_r, length.out=8 + 2*flank_num)
  
  y_t_header <- GetPlotFractionalCoords(0, 1.1)[2]

  cols_list <- c(rep("#00ADEF", flank_num), rep("#ED1C24", 8),
                 rep("#00ADEF", flank_num))
  digit_list <- c(flank_num:1, 1:8, 1: flank_num)
  text(x_list, y_t_header, labels=digit_list, xpd=NA, col=cols_list)

  for (i in seq(length(variants))) {
    variant <- variants[i]
    if (i != 3) {
      variant_text <- ConvertTtoUandMmtoX(substr(variant, seed_pos - 20,
                                                 seed_stop + 20))
      variant_list <- unlist(strsplit(variant_text, split=""))
      text(0, y_list[i], labels=names(variants)[i], adj=c(0, 0))
      text(x_list, y_list[i], labels=variant_list, col=cols_list, adj=c(0.5, 0))
    } else {
      text(mean(c(x_l, x_r)), y_list[i], labels=variant, srt=90, adj=c(0, 0.5))
    }
  }
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

# S5B___________________________________________________________________________
PlotReporterAssayKdVsL2fc <- function(mirna,
                                      experiment="twist_reporter_assay_v2",
                                      rep=1, mm=FALSE, old=FALSE,
                                      sum_variants=TRUE, cutoff=10,
                                      alt_way=FALSE, geomean_bg=FALSE,
                                      no_duplex_bg=FALSE, combined=TRUE,
                                      exclude_comp=TRUE, combined_reps=TRUE,
                                      buffer=FALSE, ident=FALSE, pdf.plot=FALSE,
                                      xpos=xpos, ypos=ypos, height=5, width=5) {
  if (mm & old) mm_str <- "_mm_old"
  else if (mm)  mm_str <- "_mm"
  else          mm_str <- ""
  counts_full <- read.table(sprintf("ReporterScreen/%s_rep%s_counts%s.txt",
                                    experiment, rep, mm_str))
  counts_full <<- counts_full
  colnames(counts_full) <- gsub("\\.", "-", colnames(counts_full))
  if (exclude_comp) {
    counts_full <- counts_full[grep("mer", rownames(counts_full)), ]
  }
  # Remove the variants for which there are no counts in all seven samples:
  inds_zeros <- which(rowSums(counts_full) == 0)
  counts <- counts_full[-inds_zeros, ]
  # Convert both matrices to TPM.
  tpm <- 1e6*t(t(counts)/colSums(counts))
  # Define the two signal vectors
  signal <- tpm[, mirna]
  names(signal) <- rownames(counts)
  # Pre-allocate the signal vector, name it and use it to mmake bg2.
  bg <- rep(NA, length(signal))
  names(bg) <- names(signal)
  # Define the list of the miRNAs in the columns of the experiment.
  mirnas_bg <- c("miR-1", "let-7a", "miR-155", "miR-124", "lsy-6", "miR-7") 
  # Loop over the mirnas.
  for (mirna_2 in mirnas_bg) {
    inds_col <- which(!(colnames(counts) %in% c(unique(c(mirna, mirna_2)),
                                               "no_duplex")))
    inds_row <- grep(sprintf("%s_", mirna_2), names(bg))
    bg[inds_row] <- rowMeans(tpm[inds_row, inds_col])
  }
  if (no_duplex_bg) {
    bg <- tpm[, "no_duplex"]  
  }
  # Renormalize the TPMs for the background, as it is now the average of six. 
  bg <- 1e6*bg/sum(bg)
  rownames_matrix <- matrix(unlist(strsplit(names(signal), split="_")),
                        byrow=TRUE, ncol=3)
  mirna_site_strings <- apply(rownames_matrix, 1, function(row) {
    paste(row[1], row[2], sep="_")
  })
  # Make the log2(fold change) vector:
  l2fc <- log(signal/bg, 2)
  # Replace the Inf and -Inf values with NA, so that they can be excluded from
  # the mean calculation.
  l2fc[which(l2fc == Inf)] <- NA
  l2fc[which(l2fc == -Inf)] <- NA
  l2fc_mean <- rep(NA, length(unique(mirna_site_strings)))
  names(l2fc_mean) <- unique(mirna_site_strings)
  # Allocate the standard error vector from the mean vector.
  l2fc_sem <- l2fc_mean
  # Allocate the sum of the signal and background variants from the mean vector.
  signal_sum <- l2fc_mean
  bg_sum <- l2fc_mean
  # Loop over the 163 mirna-site types in the experiment.
  for (str in unique(mirna_site_strings)) {
    # Find the indeces of all the rows for that miRNA-site type.
    inds_grep <- grep(paste0(str, "_"), names(signal), fixed=TRUE)
    # Use those indeces to add to the signal and background vectors.
    signal_sum[str] <- sum(signal[inds_grep])
    bg_sum[str] <- sum(bg[inds_grep])
    # Restrict the indeces to only be those which are above the cutoff.
    inds_all <- intersect(inds_grep, which(bg >= cutoff))
    # Use these indeces to calculate the mean and standard error of the
    # log2(fold-change) values.
    l2fc_mean[str] <- mean(l2fc[inds_all], na.rm=TRUE)
    l2fc_sem[str] <- sd(l2fc[inds_all], na.rm=TRUE)/sqrt(length(inds_all))
  }
  # Make the log2(fold-change) for the summed variants.
  l2fc_sum <- log(signal_sum/bg_sum, 2)
  if (combined_reps) {
    l2fc_1 <- l2fc
    l2fc_sum_1 <- l2fc_sum

    rep <- setdiff(c(1, 2), rep)
    counts_full <- read.table(sprintf("ReporterScreen/%s_rep%s_counts%s.txt",
                                      experiment, rep, mm_str))
    colnames(counts_full) <- gsub("\\.", "-", colnames(counts_full))
    if (exclude_comp) {
      counts_full <- counts_full[grep("mer", rownames(counts_full)), ]
    }
    # Remove the variants for which there are no counts in all seven samples:
    inds_zeros <- which(rowSums(counts_full) == 0)
    counts <- counts_full[-inds_zeros, ]
    # Convert both matrices to TPM.
    tpm <- 1e6*t(t(counts)/colSums(counts))
    # Define the two signal vectors
    signal <- tpm[, mirna]
    names(signal) <- rownames(counts)
    # Pre-allocate the signal vector, name it and use it to mmake bg2.
    bg <- rep(NA, length(signal))
    names(bg) <- names(signal)
    # Define the list of the miRNAs in the columns of the experiment.
    mirnas_bg <- c("miR-1", "let-7a", "miR-155", "miR-124", "lsy-6", "miR-7") 
    # Loop over the mirnas.
    for (mirna_2 in mirnas_bg) {
      inds_col <- which(!(colnames(counts) %in% c(unique(c(mirna, mirna_2)),
                                               "no_duplex")))
      inds_row <- grep(sprintf("%s_", mirna_2), names(bg))
      bg[inds_row] <- rowMeans(tpm[inds_row, inds_col])
    }
    # Renormalize the TPMs for the background, as it is now the average of six. 
    bg <- 1e6*bg/sum(bg)
    rownames_matrix <- matrix(unlist(strsplit(names(signal), split="_")),
                          byrow=TRUE, ncol=3)
    mirna_site_strings <- apply(rownames_matrix, 1, function(row) {
      paste(row[1], row[2], sep="_")
    })
    # Make the log2(fold change) vector:
    l2fc <- log(signal/bg, 2)
    # Replace the Inf and -Inf values with NA, so that they can be excluded from
    # the mean calculation.
    l2fc[which(l2fc == Inf)] <- NA
    l2fc[which(l2fc == -Inf)] <- NA
    l2fc_mean <- rep(NA, length(unique(mirna_site_strings)))
    names(l2fc_mean) <- unique(mirna_site_strings)
    # Allocate the standard error vector from the mean vector.
    l2fc_sem <- l2fc_mean
    # Allocate the sum of the signal and background variants from the mean vector.
    signal_sum <- l2fc_mean
    bg_sum <- l2fc_mean
    # Loop over the 163 mirna-site types in the experiment.
    for (str in unique(mirna_site_strings)) {
      # Find the indeces of all the rows for that miRNA-site type.
      inds_grep <- grep(paste0(str, "_"), names(signal), fixed=TRUE)
      # Use those indeces to add to the signal and background vectors.
      signal_sum[str] <- sum(signal[inds_grep])
      bg_sum[str] <- sum(bg[inds_grep])
      # Restrict the indeces to only be those which are above the cutoff.
      inds_all <- intersect(inds_grep, which(bg >= cutoff))
      # Use these indeces to calculate the mean and standard error of the
      # log2(fold-change) values.
      l2fc_mean[str] <- mean(l2fc[inds_all], na.rm=TRUE)
      l2fc_sem[str] <- sd(l2fc[inds_all], na.rm=TRUE)/sqrt(length(inds_all))
    }
    # Make the log2(fold-change) for the summed variants.
    l2fc_sum <- (l2fc_sum_1 + log(signal_sum/bg_sum, 2))/2
  }
  # Set up the arguments to get the correct Kd values:
  if (mirna == "miR-1") {
    combined <- FALSE
    buffer <- TRUE
  }
  if (mirna == "miR-1") buffer <- TRUE
  if (mirna == "miR-7") {
    mirna <- "miR-7-23nt"
    experiment_kd <- "equilibrium2_nb"
  } else {
    experiment_kd <- "equilibrium"
  }
  if (mirna == "miR-124") compcorrect <- TRUE
  else                    compcorrect <- FALSE

  kds <- SubfunctionCall(EquilPars, sitelist="paperfinal",
                         experiment=experiment_kd)
  # Relabel miR-7:
  if (mirna == "miR-7-23nt") {
    mirna <- "miR-7"
  }
  mirna_site_inds <- grep(paste0(mirna, "_"), names(l2fc_mean))
  site_names <- sapply(names(l2fc_mean), function(mirna_site) {
    unlist(strsplit(mirna_site, split="_"))[2]
  })
  # Get the colors for the points:
  cols <- rep(ConvertRColortoRGB("gray", alpha=0.1), length(l2fc_mean))
  cols[mirna_site_inds] <- kSiteColors[site_names[mirna_site_inds]]
  # black_inds <- which(names(l2fc_mean) %in% c("miR-155_AACGAGG", "lsy-6_AACGAGGA"))
  # cols[black_inds] <- "black"

  # Make a Kd vector, with 1 as the default value for each:
  kds_all <- rep(1, length(l2fc_mean))
  kds_all[mirna_site_inds] <- kds[paste0(site_names[mirna_site_inds], "_Kd"), 2]
  SubfunctionCall(FigureSaveFile)
  xmin <- 0.0001
  xmax <- 2
  ymin <- -0.5
  ymax <- 0.1
  BlankPlot(log='x', inv='x')
  AddLogAxis(1, label="Relative Kd")
  AddLinearAxis(2, label.space=0.1, tick.space=0.02, label="Fold change (log2)")
  abline(0, 0, col="gray")
  if (sum_variants) {
    y <- l2fc_sum
    cutoff <- mean(l2fc_sum[-mirna_site_inds]) - 2*sd(l2fc_sum[-mirna_site_inds])
    abline(cutoff, 0, lty=2)
  } else {
    y <- l2fc_mean
    arrows(kds_all[mirna_site_inds],
           l2fc_mean[mirna_site_inds] - 1.96*l2fc_sem[mirna_site_inds],
           kds_all[mirna_site_inds],
           l2fc_mean[mirna_site_inds] + 1.96*l2fc_sem[mirna_site_inds],
           length=error_height_final*par()$cex, lwd=par()$cex, angle=90, code=3)
  }
  Points(kds_all, y, col=cols)
  if (ident) {
    identify(kds_all, y, labels=names(y))
  }
  if (sum_variants) {
    legend_inds <- which(l2fc_sum[mirna_site_inds] <= cutoff)
  } else {
    legend_inds <- which(l2fc_mean[mirna_site_inds] + 1.96*l2fc_sem[mirna_site_inds] <= 0)
  }
  xy <- GetPlotFractionalCoords(0.975, 0.975, log='x', inv='x')
  text(xy[1], xy[2], label=mirna, adj=c(1, 1))
  # if (mirna == "miR-1") {
  xy <- GetPlotFractionalCoords(0.025, 0.975, log='x', inv='x')
  if (experiment == "twist_reporter_assay_v2") {
    exp_text <- "V.2"
  } else {
    exp_text <- "V.1"
  }
  # if (combined_reps) {
  #   text(xy[1], xy[2], label=sprintf("%s, Combined reps", exp_text, rep), adj=c(0, 1))    
  # } else {
  #   text(xy[1], xy[2], label=sprintf("%s, Rep. %s", exp_text, rep), adj=c(0, 1))    
  # }
  # }
  # LEGEND.
  legend.coords <- GetPlotFractionalCoords(0.0125, 0, log='x', inv='x')
  # Variables for assigning how many sites are in the left-hand legend.

  # List of sites for the legend, in the order of their Kd values:
  sites_cutoff <- site_names[mirna_site_inds[legend_inds]]
  kds_cutoff <- kds_all[mirna_site_inds[legend_inds]]
  sites_cutoff_order <- sites_cutoff[order(kds_all[mirna_site_inds[legend_inds]])]

  legend_cutoff <- 15
  legend_inds_1 <- 1:min(length(legend_inds), legend_cutoff)
  sites_legend <- sites_cutoff_order[legend_inds_1]

  Legend(legend.coords, legend=ConvertTtoUandMmtoX(sites_legend),
         col=kSiteColors[sites_legend], xjust=0, yjust=0, y.intersp=0.8)
  if (length(legend_inds) > legend_cutoff) {
    legend.coords <- GetPlotFractionalCoords(0.475, 0, log='x', inv='x')
    sites_legend <- sites_cutoff_order[(legend_cutoff + 1):length(sites_cutoff_order)]
    Legend(legend.coords, legend=ConvertTtoUandMmtoX(sites_legend),
           col=kSiteColors[sites_legend], xjust=0, yjust=0, y.intersp=0.8)
  }
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}





################################################################################
# SUPPLEMENTARY FIGURE 6
################################################################################

# S6F___________________________________________________________________________
PlotFlankKdsVsRepression <- function(experiment="equilibrium",
                                     sitelist="resubmissionfinal", new=TRUE,
                                     bg_method=3, exrib=TRUE, ident=FALSE,
                                     xpos=20, ypos=20, height=4, width=4,
                                     pdf.plot=FALSE) {
  rep_df <- SubfunctionCall(GetRepFlankLinearModel)

  rep_df_global <<- rep_df
  break
  rep_df <- rbind(c(0, NA, NA, 0, 0), rep_df)
  rownames(rep_df)[1] <- "AA.AA"
  kds <- 10^SubfunctionCall(GetAverageFlanks, sitelist="resubmissionfinal")
  if (ident) {
    height <- 8
    width <- 8
  }
  SubfunctionCall(FigureSaveFile)
  xmin <- 0.1
  xmax <- 100
  ymin <- -1
  ymax <- 1
  BlankPlot(log='x', inv='x')
  AddLogAxis(1, label="Relative Kd")
  AddLinearAxis(2, tick.space=0.1, label.space=0.5, label="Fold change (log2)")
  names_intersect <- intersect(names(kds), rownames(rep_df))
  x <- kds[names_intersect]
  fc_df <- rep_df[names_intersect, ]
  y <- fc_df[, 1]
  y_error <- list(fc_df[, 4], fc_df[, 5])
  y_weights <- 1/(fc_df[, 4] - fc_df[, 5])^2
  log_kds_fit <- log10(x[-1])
  fc_fit <- y[-1]
  weights_fit <- y_weights[-1]
  x.line <- seq(-1, 2, length.out=400)
  lm_line <- lm(fc_fit ~ log_kds_fit, weights=weights_fit)
  coefs <- lm_line$coefficients
  print(summary(lm_line))
  print(confint(lm_line))
  b <- coefs[1]
  m <- coefs[2]

  kSLopeFlanks <<- m

  slopes_lm["flanks"] <<- m
  predict1 <- predict(lm_line, data.frame(log_kds_fit=x.line), interval='confidence')
  polygon(c(10^x.line, rev(10^x.line)), c(predict1[,2], rev(predict1[, 3])),
          col=ConvertRColortoRGB("gray", alpha=0.3), border=NA)
  abline(b, m)
  xy <- GetPlotFractionalCoords(0.95, 0.95, log='x', inv='x')
  AddCorrelationToPlot(log(x[-1]), y[-1], xpos=xy[1], ypos=xy[2],
                       weights=y_weights[-1], adj=1, rsquared=TRUE)
  Points(x, y, col=GetColorFunction(names(kds)))
  if (ident & (class(pdf.plot) != "character")) {
    identify(x, y, labels=names_intersect)
  }
  if (class(pdf.plot) == "character") {
    dev.off()
  }
  print("finished plot")
}

# S6G___________________________________________________________________________
PlotAllSiteInputDistribution <- function(mirna, site, experiment="equilibrium",
                                         n_constant=5,
                                         sitelist="resubmissionfinal",
                                         buffer=FALSE, combined=FALSE,
                                         mir.start=1, mir.stop=14,
                                         noconstant=FALSE, xpos=20, ypos=20,
                                         height=4, width=4, pdf.plot=FALSE) {

  kConditionsColor <- c("black", "deeppink", "darkmagenta", "blue")
  names(kConditionsColor) <- c("I_combined", "0.4", "4", "40")
  sXc <- SubfunctionCall(SitesXCounts)
  # conditions <- c(colnames(sXc)[1 + combined], colnames(sXc)[3:(ncol(sXc) - 1)])
  # conditions <- gsub("12.65", replace="12.6", conditions)
  # conditions <- gsub("1.265", replace="1.26", conditions)
  SubfunctionCall(FigureSaveFile)
  xmin <- 1e-8
  xmax <- 1
  ymin <- 0
  ymax <- 1
  BlankPlot(log='x')
  x_ecdf <- 10^seq(-8, 0, by=0.01)
  sapply(names(kConditionsColor), function(condition) {
    # if (condition == "I_combined") {
    #   data_1 <- SubfunctionCall(GetPairingFlankData, condition="I", alt_mir_exp_cond="let-7a_equilibrium_I")
    #   data_2 <- SubfunctionCall(GetPairingFlankData, condition="I", alt_mir_exp_cond="miR-124_equilibrium_I")
    #   data_3 <- SubfunctionCall(GetPairingFlankData, condition="I", alt_mir_exp_cond="miR-1_kin_pilot_I_TGT")
    #   data_4 <- SubfunctionCall(GetPairingFlankData, condition="I")
    #   data <- rbind(data_1, data_2, data_3, data_4)
    # } else {
    #   data <- SubfunctionCall(GetPairingFlankData)  
    # }
    data <- SubfunctionCall(GetPairingFlankData)
    log_mean <- exp(mean((mir.stop - mir.start + 1)*log(data[["plfold"]])))
    message(condition)
    message("geometric mean:")
    print(log_mean)
    ecdf <- ecdf(data[["plfold"]]^(mir.stop - mir.start + 1))
    lines(x_ecdf, ecdf(x_ecdf), col=kConditionsColor[condition])
  })
  message("got here")
  # Make the ECDF for the input, 0.4 Ago, and the resampled:
  AddLogAxis(1,
             label=sprintf("Accessibility score",
                           mir.start, mir.stop),
             maglabel=2)
  AddLinearAxis(2, tick.space=0.2, label.space=0.2, label="Cumulative fraction")
  A.stock <- kAgoStock[mirna, "equilibrium"] 
  ago_dils <- A.stock*c(0.004, 0.04, 0.4)*1000
  ago_formats <- sprintf("AGO2-miR-1; %s pM", signif(ago_dils, 2))
  xy <- GetPlotFractionalCoords(0.005, 1.05, log='x')
  # legend(xy[1], xy[2], rep("", 4), lty=1, seg.len=1,
  #        col=kConditionsColor,
  #        bty="n",
  #        xpd=NA)
  xy <- GetPlotFractionalCoords(0, 1.05, log='x')
  legend(xy[1], xy[2], c("Input", ago_formats), col=kConditionsColor, lty=1,
         lwd=1, seg.len=1, y.intersp=0.9, bty="n", xpd=NA)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

# S6H___________________________________________________________________________
PlotFlankPlCorrelationWindowsNew <- function(mirna, site,
                                             experiment="equilibrium",
                                             condition="I_combined", n_constant=5,
                                             sitelist="resubmissionfinal", test=FALSE,
                                             buffer=TRUE,
                                             xpos=20, ypos=20, m=49, b=50,
                                             c_s=0.4, c_e=1, height=4, width=7,
                                             pdf.plot=FALSE) {
  # Load the correlation matrix, and square it for r^2.
  cormatrix <- SubfunctionCall(GetPlFoldCorrelationMatrix)
  cormatrix <- cormatrix^2
  xmin=0
  xmax=36
  ymin=1
  ymax=20
  SubfunctionCall(FigureSaveFile)
  BlankPlot()
  ymin <- 1
  ymax <- 20
  
  AddLinearAxis(2, alt_lab=c(1, 5, 10, 15, 20),
                alt_lab_pos=c(1, 5, 10, 15, 20),
                label="Window length (nt)", line=0.8, adj_pos=0.25)

  cormatrix <- cormatrix[1:ymax, 3:28]
  xlefts <- (
    rep(seq(1, ncol(cormatrix)), ymax) + 
    rep(rep(c(-0.5, 0), each=ncol(cormatrix)), ymax/2)
  )
  xright <- xlefts+1
  ybottom <- rep(seq(1, ymax), each=ncol(cormatrix)) - 0.5
  ytop <- ybottom + 1

  b <- 0
  m <- 100
  col.inds <- round(cormatrix*100)

  # Parameters specifying the start and end of the rainbow colorspace usee:
  c_s <- 0.9 # color_start
  c_e <- 0.70  # color_end
  color.dist = rev(rainbow(100, start=c_s, end=c_e))
  col.inds <- sapply(t(col.inds), function(col.ind) {
    min(max(1, col.ind), 100)
    })
  seedpos = grep("sp", colnames(cormatrix))
  par(xpd=NA)
  
  mir_seq_rev <- StrRev(ConvertTtoUandMmtoX(kMirnaSeqs[mirna]))
  mir_seq_rev_list <- unlist(strsplit(mir_seq_rev, split=""))
  mir_len <- length(mir_seq_rev_list)

  mirpos <- seq(grep(sprintf("^f5p%s", mir_len - 8), colnames(cormatrix),
                     perl=TRUE),
                grep("^sp1$", colnames(cormatrix), perl=TRUE))
  num_3pterm <- mir_len - 16
  num_3p <- 4
  num_mid <- 4
  col_bg <- kSiteColors["mirna_bg"]
  col_seed <- kSiteColors["mirna_seed"]
  mirna_nt_col <- c(rep(col_bg, num_3pterm), rep(col_bg, num_3p),
                  rep(col_bg, num_mid), "black", rep(col_seed, 6), col_bg)
  seed_offset <- 0.5
  mirna_y_pos <- c(rep(0, mir_len - 8), rep(seed_offset, 7), 0) - 1
  # Makes the ticks for nt 2–7.
  text(x=seedpos[-length(seedpos)], 0.5, rep("|", length(seedpos) - 1))
  # Makes the ticks for 
  text(x=min(seedpos) - c(12, 7, 2, 1, -7), 0, rep("|", 4), col="black")

  text(x=mirpos, y=mirna_y_pos, mir_seq_rev_list, adj=0.5,
       col=mirna_nt_col)
  number_inds <- c(min(seedpos) - c(12, 7, 2), seedpos)
  mirna_nt_col[which(mirna_nt_col == col_bg)] <- "black"
  text(x=number_inds, y=mirna_y_pos[number_inds - mirpos[1] + 1] - 1,
       c(20, 15, 10, seq(8, 1, by=-1)), adj=0.5, col=mirna_nt_col[number_inds - mirpos[1] + 1])

  # This find s the best mirna correlation in the matrix:
  best_ind <- which(c(t(cormatrix)) == max(c(t(cormatrix))))
  # This finds the 1-15 window used
  used_window = ncol(cormatrix)*13 +seedpos[1]

  rect(xlefts, ybottom, xright, ytop, col=color.dist[col.inds], border=NA)
  rect(xlefts[c(used_window, best_ind)], ybottom[c(used_window, best_ind)],
       xright[c(used_window, best_ind)], ytop[c(used_window, best_ind)],
       col=color.dist[col.inds][c(used_window, best_ind)], lwd=1,
       border=c("gray", "black"))

  # par(kPlotParameters)
  dens <- density(c(cormatrix), n=300, from=0, to=1)
  x.mat <- matrix(dens$x, nrow=2)
  y.mat <- matrix(dens$y, nrow=2)*1.5
  offset <- ncol(cormatrix) + 4
  x.mat <- rbind(x.mat[1,],
                 x.mat[1,],
                 x.mat[2,], 
                 x.mat[2,],
                 x.mat[2,])
  # Convert x.mat to plot coordinates, starting at 1 and going to 30:
  offset <- 33.5
  x_m <- 19
  x_b <- 1
  x.mat.use <- x.mat*x_m + x_b
  y.mat <- rbind(rep(offset, ncol(y.mat)),
                 offset + y.mat,
                 rep(offset, ncol(y.mat)),
                 rep(NA, ncol(y.mat)))
  x.col.inds <- sapply(round(x.mat[1,]*m + b), function(i) {min(max(i, 1), 100)})
  x.cols <- color.dist[x.col.inds]
  polygon(c(y.mat), c(x.mat.use), border=x.cols, col=x.cols)
  # axis(1, at=seq(0, 1, by=0.1)*x_m + x_b, labels=seq(0, 1, by=0.1), lwd=par()$lwd,
  #      pos=offset,)
  axis(2, at=seq(0, 1, by=0.2)*x_m + x_b, labels=seq(0, 1, by=0.2), lwd=par()$lwd,
       pos=offset,)

  xy <- GetPlotFractionalCoords(1, 0.5)
  text(30, xy[2],
       labels=bquote(italic(r)^2~.(" with flanking dinucleotide")~italic(K)[D]*.("s")),
       xpd=NA, srt=90)


  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

# S6I___________________________________________________________________________
PlotFlanksSamplePl <- function(mirna, site, condition, experiment="equilibrium",
                               n_constant=5, sitelist="resubmissionfinal",
                               buffer=TRUE, mir.start=1, mir.stop=14,
                               noconstant=FALSE, height=4, width=4, depth=20000,
                               p_score="plfold", matchdist=FALSE, xpos=20,
                               ypos=20, pdf.plot=FALSE) {
  kConditionsColor <- c("black", "deeppink", "darkmagenta", "blue")

  # Get data:
  data.I <- SubfunctionCall(GetPairingFlankData, condition="I_combined")
  data.A <- SubfunctionCall(GetPairingFlankData)

  inds_AU <- c(sapply(seq(50), function(i) {
    SampleByDinucleotideEnrichment(data.I, data.A, nrow(data.A))
  }))

  GetPercentTAEffect <- function(condition) {
    data.A <- SubfunctionCall(GetPairingFlankData)
    mean.I <- mean(log(data.I[[p_score]]))
    mean.A <- mean(log(data.A[[p_score]]))
    inds <- c(sapply(seq(50), function(i) {
      SampleByDinucleotideEnrichment(data.I, data.A, nrow(data.A))
    }))
    mean.I.sample <- mean(log(data.I[[p_score]][inds]))
    # test_inside <<- c(mean.I, mean.A, mean.I.sample)
    out <- (mean.I.sample - mean.I)/(mean.A - mean.I)
    message("GetPercentTAEFfect:")
    print(out)
    return(out)
  }



  if (matchdist) {
    # Calculate the mean and sd of the target distribution:
    dist.A.mean <- mean((mir.stop - mir.start + 1)*log(data.A[["plfold"]]))
    dist.A.sd <- sd((mir.stop - mir.start + 1)*log(data.A[["plfold"]]))
    SampleCostFunction <- function(pars) {
      tick <<- tick + 1
      par.exponent <- pars[1]
      par.size <- ceiling(Logistic(pars[2], nrow(data.I)))
      inds <- SampleBySiteAccess(data.I, par.size, "plfold",
                                 exponent=(mir.stop - mir.start + 1)*par.exponent)
      residual.mean <- (mean((mir.stop - mir.start + 1)*log(data.I[["plfold"]][inds])) - dist.A.mean)^2
      residual.sd <- (sd((mir.stop - mir.start + 1)*log(data.I[["plfold"]][inds])) - dist.A.sd)^2
      residual <- residual.mean + residual.sd
      return(residual)
    }
    tick <- 1
    exp.par <- optim(c(0, 0),SampleCostFunction)$par
    par.exponent <- exp.par[1]
    par.size <- ceiling(Logistic(exp.par[2], nrow(data.I)))
  } else {
    par.exponent <- 1
    par.size <- 20000
  }
  inds <- SampleBySiteAccess(data.I, par.size, "plfold",
                             exponent=(mir.stop - mir.start + 1)*par.exponent)
  while(min(NumFlanks(data.I[inds, ])) < 50) {
    inds <- c(inds, SampleBySiteAccess(data.I, par.size, "plfold",
                                       exponent=(mir.stop - mir.start + 1)*par.exponent))
  }
  # ECDF plot:

  xmin <- 1e-8
  xmax <- 1
  ymin <- 0
  ymax <- 1
  if (class(pdf.plot) == "character") pdf.plot <- paste0(pdf.plot, "_left")
  SubfunctionCall(FigureSaveFile)
  BlankPlot(log='x')
  # Make the ECDF for the input, 0.4 Ago, and the input resampled both for
  # structural accessibility and AU-content.
  x_ecdf <- 10^seq(-8, 0, by=0.01)
  ecdf.I <- ecdf(data.I[["plfold"]]^(mir.stop - mir.start + 1))
  ecdf.A <- ecdf(data.A[["plfold"]]^(mir.stop - mir.start + 1))
  ecdf.I.sample <- ecdf(data.I[["plfold"]][inds]^(mir.stop - mir.start + 1))
  ecdf.I.sampleAU <- ecdf(data.I[[p_score]][inds_AU]^(mir.stop - mir.start + 1))
  #.................................. Lines ....................................
  lines(x_ecdf, ecdf.I(x_ecdf), col=kConditionsColor[1])
  lines(x_ecdf, ecdf.A(x_ecdf), col=kConditionsColor[2])
  lines(x_ecdf, ecdf.I.sampleAU(x_ecdf), col=kConditionsColor[3])
  lines(x_ecdf, ecdf.I.sample(x_ecdf), col=kConditionsColor[4])
  #.................................. Axes .....................................
  AddLogAxis(1, label= "Accessibility score", maglabel=2)
  AddLinearAxis(2, tick.space=0.2, label.space=0.2, label="Cumulative fraction")
  #................................. Legend ....................................
  xy <- GetPlotFractionalCoords(fx=0, fy=1.05, log="x")
  A.stock <- kAgoStock[mirna, "equilibrium"] 
  ago_dils <- A.stock*c(as.numeric(condition)/100)*1000
  ago_formats <- sprintf("AGO2-miR-1; %s pM", signif(ago_dils, 2))
  legend(x=xy[1], y=xy[2],
         legend=c("Input ", ago_formats, "Input matched for",
                  "flanking dinucleotide", "composition", "Input matched for",
                  "accessibility score"),
         col=c(kConditionsColor[1:3], NA, NA, kConditionsColor[4]),
         lty=c(1, 1, 1, 0, 0, 1, 0), lwd=1, bty="n", seg.len=1, y.intersp=0.9)
  if (class(pdf.plot) == "character") {
    dev.off()
    pdf.plot <- paste0(unlist(strsplit(pdf.plot, split="_"))[1], "_right")
  }
  ############################## PLOT 2 ########################################
  SubfunctionCall(FigureSaveFile)
  xmin <- 1e-5
  xmax <- 1e-1
  ymin <- xmin
  ymax <- xmax
  BlankPlot(log='xy')
  AddLogAxis(1, label="Sampled flanking dinucleotide frequencies")
  AddLogAxis(2, label="Observed flanking dinucleotide frequencies")
  #................................. Points ....................................
  points(x=Norm(NumFlanks(data.I[inds, ])), y=Norm(NumFlanks(data.A)),
         col=GetColorFunction(names(NumFlanks(data.A))))
  xy <- GetPlotFractionalCoords(fx=0.05, fy=0.95, log='xy')
  #.................................. r^2 ......................................
  AddCorrelationToPlot(x=log(Norm(NumFlanks(data.A))),
                       y=log(Norm(NumFlanks(data.I[inds, ]))),
                       xy[1], xy[2], rsquared=TRUE, adj=c(0, 1))
  #............................... y = x line ..................................
  abline(0, 1, lty = 2)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

PlotFlowthrough <- function(mirna, experiment="equil_flowthrough", n_constant=5,
                            sitelist="resubmissionfinal", combined=FALSE,
                            enrich_against_enrich=FALSE, log_y=TRUE,
                            pdf.plot=FALSE) {
  # Get the data required for the figure.
  sXc <- SubfunctionCall(SitesXCounts)
  # Normalize the data, and calculate enrichments:
  sXf <- t(t(sXc)/colSums(sXc))
  sXR <- sXf/sXf[, "I"]
  # Check if enrich against enrich.
  if (enrich_against_enrich) {
    x <- sXR[, "40_2h"]
    y <- sXR[, "40_2h_nc1"]
  } else {
    x <- sXR[, "40_2h_nc1"]
    y <- sXR[, "40_2h_ft"]    
  }
  # Conditions for giving the fracion remaining of the flowthrough experiments.
  # if (log_y == FALSE & !(enrich_against_enrich)) {
  #   f_nc1 <- sXf[, "40_2h_nc1"]
  #   f_nc2 <- sXf[, "40_2h_nc2"]
  #   f_ft  <- sXf[, "40_2h_ft"]
  #   f_i   <- sXf[, "I"]
  #   p <- CalculateFlowthroughFractions(f_nc1, f_nc2, f_ft, f_i)
  #   y <- y*p[1]
  # }
  if (!(enrich_against_enrich)) {
    f_nc1 <- sXf[, "40_2h_nc1"]
    f_nc2 <- sXf[, "40_2h_nc2"]
    f_ft  <- sXf[, "40_2h_ft"]
    f_i   <- sXf[, "I"]
    p <- CalculateFlowthroughFractions(f_nc1, f_nc2, f_ft, f_i)
    y <- y*p[1]
    y <- 1 - y  
  }
  # Generate plot.
  SubfunctionCall(FigureSaveFile, height=5, width=5)
  xmin <- 1e-1
  xmax <- 2e2
  if (enrich_against_enrich) {
    ymin <- xmin
    ymax <- xmax
    log_str <- "xy"
  } else if (log_y) {
    ymin <- 1e-3  
    ymax <- 1
    log_str <- "xy"

  } else {
    ymin <- 0
    ymax <- 0.8
    log_str <- "x"
  }
  BlankPlot(log=log_str)
  site.colors <- kSiteColors[rownames(sXc)]
  # Make x=y line and the lines connecting the points to the x = y line.
  if (enrich_against_enrich) {
    segments(xmin, ymin, xmax, ymax, lty=line_dash_length)    
    segments(x, x, x, y, lty=line_dash_length, col=site.colors)    
    legend.coords <- GetPlotFractionalCoords(1.12, 0.01, log=log_str)
    y.intersp <- 0.75
  } else if (log_y) {
    # Linear model to plot trend line.
    lm_df <- data.frame(ly=log(y), lx=log(x))
    lm_model <- lm(ly ~ lx, data=lm_df)
    lm_b <- lm_model$coefficients[1]
    lm_m <- lm_model$coefficients[2]
    x_line <- c(log(xmin), log(xmax))
    y_line <- exp(lm_m*x_line + lm_b)
    y_0 <- exp(lm_m*log(x) + lm_b)
    segments(x, y_0, x, y, lty=line_dash_length, col=site.colors)    

    lines(exp(x_line), y_line, lty=2)
    legend.coords <- GetPlotFractionalCoords(1.12, 0.01, log=log_str)
    y.intersp <- 0.75
  } else {
    # segments(xmin, 1, xmax, 1, lty=line_dash_length)    
    segments(x, 0, x, y, lty=line_dash_length, col=site.colors)    
    legend.coords <- GetPlotFractionalCoords(0.4, 0.025, log=log_str)
    y.intersp <- 0.9
  }
  # AXES.........
  if (enrich_against_enrich) {
    AddLogAxis(1, label="Enrichment with NC over HB")
    AddLogAxis(2, label="Enrichment with NC over NC")

  } else {
    AddLogAxis(1, label="Enrichment in nitrocellulose-bound RNA")
    if (log_y) {
      AddLogAxis(2, label="% depleted in flowthrough RNA", percent=TRUE)    
    } else {
      AddLinearAxis(2, tick.space=0.1, label.space=0.2,
                    label="% depleted in flowthrough RNA", percent=TRUE)  
    }    
  }
  # POINTS.......
  Points(x, y, col=site.colors)
  # LEGEND......
  Legend(legend.coords, legend=ConvertTtoUandMmtoX(rownames(sXc)),
         col=kSiteColors[rownames(sXc)], xjust=1, yjust=0, y.intersp=y.intersp)
  # CORRELATION.......
  if (enrich_against_enrich | log_y) {
    xy <- GetPlotFractionalCoords(0.05, 0.975, log='xy')
    AddCorrelationToPlot(log(x), log(y), xy[1], xy[2], rsquared=TRUE)
  # } else if (log_y) {
  #   xy <- GetPlotFractionalCoords(0.95, 0.05, log=log_str)
  #   AddCorrelationToPlot(log(x), log(y), xy[1], xy[2], method="spearman", adj=1)
  } else {
    xy <- GetPlotFractionalCoords(0.05, 0.975, log=log_str)
    AddCorrelationToPlot(log(x), log(y), xy[1], xy[2], method="spearman", adj=0)    
  }
  # DONE.......
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}


PlotInVivoDinucAUfreqs <- function(mirna, experiment="repression_hela_cs",
                                   sitelist="resubmissionfinal", adjusted_height=TRUE,
                                   trim_mir_name=TRUE, height=5, width=6,
                                   new=TRUE, bg_method=3, xpos=20, ypos=20, pdf.plot=FALSE) {
  # Get kds for all site-types of the mirna.
  AUfreq_df <- SubfunctionCall(GetAUFreqTable)
  AUfreq_df[nrow(AUfreq_df), ncol(AUfreq_df)] <- 20
  AUfreq_df <- AUfreq_df[-(nrow(AUfreq_df) - 1), ]
  AUfreq_df <- AUfreq_df[which(AUfreq_df[, 3] >= 20), ]
  sites <- rownames(AUfreq_df)
  xmin <- -0.2  
  xmax <- 1
  ymin <- 0
  ymax <- nrow(AUfreq_df) + 0.5
  y <- nrow(AUfreq_df) - seq(nrow(AUfreq_df)) + 1
  names(y) <- sites
  if (adjusted_height) {
    height <- (ymax*0.80 + 5)/4
  } else {
    height <- 4.5
  }
  SubfunctionCall(FigureSaveFile)
  BlankPlot()
  xmin <- 0
  sem <- AUfreq_df$std/sqrt(AUfreq_df$n)

  x0 <- AUfreq_df$mean - sem
  x1 <- AUfreq_df$mean + sem
  x_error <- list(x0, x1)
  cols <- kSiteColors[sites]
  cols[is.na(cols)] <- "black"
  Points(AUfreq_df$mean, y, col=cols, x_error=x_error)
  AddLinearAxis(1, percent=TRUE, label.space=0.2, tick.space=1, adj_pos=0.9,
                label="Average flanking dinucleotide AU content (%)")
  title.xy <- GetPlotFractionalCoords(fx=0.975, fy=1)
  xmin <- -0.2
    if (length(grep("miR-7", mirna)) > 0 & trim_mir_name) {
      title.text <- "miR-7"
    } else {
      title.text <- mirna
    }
  text(title.xy[1], title.xy[2], labels = title.text, adj=c(1, 0), xpd=NA)
  title_pos <- GetPlotFractionalCoords(-0.1, 0)[1]
    site_labels_table <- cbind(sites, AUfreq_df[, 3])
  site_labels <- apply(site_labels_table, 1, function(row) {
    paste0(row[1], " (", row[2], ")")
  })

  # site_labels <- sites
  site_labels[length(site_labels)] <- "3' UTR average"
 
  text(title_pos, y, labels=site_labels, adj=0, xpd=NA)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

PlotSiteKdRepResidualVsAUFreq <- function(mirna, experiment="equilibrium",
                                    n_constant=5, sitelist="resubmissionfinal",
                                    sites_use="seed",
                                    single_site=FALSE, best_site=FALSE,
                                    combined=TRUE, buffer=FALSE, n_cutoff=20,
                                    kd_cutoff=0.10, bulk=FALSE, new=TRUE, old=FALSE,
                                    cat_colors=FALSE, noncanon=FALSE,
                                    merge=FALSE, best=TRUE, bg_method=3,
                                    threePseq=TRUE, height=5, width=5, xpos=20,
                                    ypos=20,
                                    pdf.plot=FALSE) {
  xmin <- 0
  xmax <- 1
  ymin <- -0.3
  ymax <- 0.3

  kds <- SubfunctionCall(EquilPars)
  rownames(kds) <- gsub("_Kd", "", rownames(kds))
  sites <- grep(sprintf("(?:%s|None)", mirna), rownames(kds), perl=TRUE,
                  invert=TRUE, value=TRUE)
  if (sites_use == "all") {
    n_cutoff <- -Inf
    kd_cutoff <- Inf
    width <- 7
    adjusted <- TRUE
    l.xy <- GetPlotFractionalCoords(fx=1.2, fy=1.1)
    yjust <- 1
    l.cex <- 0.65
  } else {
    adjusted <- FALSE
    l.xy <- GetPlotFractionalCoords(fx=-0.005, fy=-0.01)
    yjust <- 0
    l.cex <- 1
  }
  SubfunctionCall(FigureSaveFile)
  BlankPlot(adjusted=adjusted)  
  AddLinearAxis(2, tic.space=0.01, label.space=0.05, label="Residual", adj=TRUE)
  AddLinearAxis(1, tick.space=0.1, label.space=0.2, label="Average flankig dinucleotid AU content(%)", percent=TRUE)

  fc <- SubfunctionCall(GetRepressionLinearModel)
  fc_sites <- fc[, 1]
  fc_sem <- fc[, 2]
  fc_nsites <- fc[, 3]
  gray_col <- "gray50"
  segments(xmin, 0, xmax, 0, lwd=0.5, col=gray_col)

  AUfreq_df <- SubfunctionCall(GetAUFreqTable)
  kds <- kds[rownames(fc), ]
  AUfreq_df <- AUfreq_df[rownames(fc), ]
  sites <- rownames(kds)
  kds_sites <- kds$Mean
  kds_sites_lci <- kds$Lower_CI
  kds_sites_uci <- kds$Upper_CI

  lm_fc_kd <- lm(fc_sites ~ log(kds_sites), weights=1/(fc[, 5] - fc[, 4])^2)

  m <- lm_fc_kd$coefficients[2]
  b <- lm_fc_kd$coefficients[1]
  fc_sites_pred <- m*log(kds_sites) + b
  fc_sites_residuals <- fc_sites - fc_sites_pred


  y.intersp <- 0.7
  l.cex <- 0.7
  colors.sites <- kSiteColors[sites]

  y_error <- list(fc[, 4], fc[, 5])
  x_error <- list(kds_sites_lci, kds_sites_uci)
  x_error = list(AUfreq_df[, 1] - (AUfreq_df[, 2])/sqrt(AUfreq_df[, 3]),
                 AUfreq_df[, 1] + (AUfreq_df[, 2])/sqrt(AUfreq_df[, 3]))

  Points(AUfreq_df[, 1], fc_sites_residuals, x_error=x_error, col=colors.sites)
  if (sites_use == "all" ) {
    segments(kds_sites, fc_sites, 10^l$rect$l, l$text$y, lwd=0.5,
             col=colors.sites, xpd=NA)
  }

  sites <- sites[order(kds_sites)]
  colors.sites <- kSiteColors[sites]
  sites <- ConvertTtoUandMmtoX(sites)

  rownames(fc) <- ConvertTtoUandMmtoX(rownames(fc))
  rownames(AUfreq_df) <- ConvertTtoUandMmtoX(rownames(AUfreq_df))
  legend_text_table <- cbind(sites, fc[sites, ][, 3], AUfreq_df[sites, 3])
  legend_text <- apply(legend_text_table, 1, function(row) {
    paste0(row[1], " (", row[2], ") (", row[3], ")")
  })

  if (mirna == "miR-7-23nt") {
    l1 <- Legend(l.xy, legend=legend_text[1:12],
                 cex=l.cex, col=colors.sites[1:12],
                 yjust=yjust, y.intersp=y.intersp, x.intersp=0.7)
    xy <- GetPlotFractionalCoords(fx=0.55, fy=-0.01)
    l2 <- legend(x=xy[1], y=xy[2],
                 legend=legend_text[13:length(legend_text)],
                 bty="n", pch=legend_pch, cex=l.cex, pt.cex=pt_cex_final,
                 col=colors.sites[13:length(legend_text)], y.intersp=y.intersp,
                 x.intersp=0.7, xjust=0, yjust=yjust, xpd=NA)
  } else {
    l <- legend(x=l.xy[1], y=l.xy[2], legend=legend_text,
                bty="n", cex=l.cex, pch=legend_pch, col=colors.sites,
                pt.cex=pt_cex_final, xjust=0, yjust=yjust, xpd=NA,
                y.intersp=y.intersp, x.intersp=0.7)
  }

  xy <- GetPlotFractionalCoords(fx=0.9, fy=0.95)
  mirna.split <- paste0(strsplit(mirna, split="-")[[1]][1:2], collapse="-")
  text(xy[1], xy[2], mirna.split, adj=1)
  # Tested the sem and the 95% confidence interval, and these work exactly the same.
  xy <- GetPlotFractionalCoords(fx=0.9, fy=0.90)
  fc_sites_global <<- fc_sites
  kds_sites_global <<- kds_sites
  fc_sem_global <<- fc_sem
  AddCorrelationToPlot(AUfreq_df, fc_sites_residuals, xpos=xy[1], ypos=xy[2],
                       rsquared=TRUE, adj=1)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
  print("finished plot")
}



PlotAllPairwiseReporterCounts <- function(mirna, experiment, condition, site,
                                          combined=TRUE, buffer=FALSE,
                                          xpos=xpos, ypos=ypos, height=4.5,
                                          width=4.5, pdf.plot=FALSE) {
  # Get the two signal data sets:
  print(mirna)
  print(experiment)
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
  print(head(bg1))
  # Get the colors for the points:
  col_bg <- ConvertRColortoRGB("gray")
  row_use <- paste(mirna, site, sep="_")
  print(row_use)
  row_ind <- which(rownames(bg1) == row_use)
  print(row_ind)
  col_site <- kSiteColors[site]
  # Make the log2fold-change for each rep, and average them:
  l2fc1 <- log(signal1/bg1, 2)
  l2fc2 <- log(signal2/bg2, 2)
  bg_use <- grep(sprintf("%s_", mirna), rownames(bg1))
  bg_use <- setdiff(1:nrow(bg1), bg_use)
  bg_use <- bg_use[seq(1, length(bg_use), length.out=15)]
  # bg_use <- seq(1, nrow(bg1), by=10)
  print(bg_use)
  print(dim(signal1))
  l2fc1 <<- l2fc1
  l2fc2 <<- l2fc2
  # Make a Kd vector, with 1 as the default value for each:
  SubfunctionCall(FigureSaveFile)
  xmin <- -5
  xmax <- 3
  ymin <- -5
  ymax <- 3
  BlankPlot()
  AddLinearAxis(1, label.space=1, tick.space=0.2, label="log2(fold change); rep 1")
  AddLinearAxis(2, label.space=1, tick.space=0.2, label="log2(fold change); rep 2")
  Points(l2fc1[bg_use,], l2fc2[bg_use,], col=col_bg)
  print(mean(l2fc1[bg_use, ], na.rm=TRUE))



  l2fc1_site_all <- l2fc1[row_ind, ]
  l2fc2_site_all <- l2fc2[row_ind, ]
  l2fc1_site <- mean(l2fc1_site_all[is.finite(l2fc1_site_all)])
  l2fc2_site <- mean(l2fc2_site_all[is.finite(l2fc2_site_all)])


  l2fc1_bg_all <- c(l2fc1[bg_use, ])
  l2fc2_bg_all <- c(l2fc2[bg_use, ])
  l2fc1_bg <- mean(l2fc1_bg_all[is.finite(l2fc1_bg_all)])
  l2fc2_bg <- mean(l2fc2_bg_all[is.finite(l2fc2_bg_all)])



  segments(l2fc1_bg,
           ymin,
           x1=l2fc1_bg,
           y1=ymax,
           col=col_bg)
  segments(xmin,
           l2fc2_bg,
           x1=xmax,
           y1=l2fc2_bg, col=col_bg)


  segments(0,
           ymin,
           x1=0,
           y1=ymax,
           col="darkgray", lty=2)
  segments(xmin,
           0,
           x1=xmax,
           y1=0, col="darkgray", lty=2)

  x <- l2fc1[row_ind,]
  y <- l2fc2[row_ind, ]

  Points(x, y, col=col_site)
  segments(l2fc1_site,
           ymin,
           x1=l2fc1_site,
           y1=ymax,
           col=col_site)
  segments(xmin,
           l2fc2_site,
           x1=xmax,
           y1=l2fc2_site, col=col_site)

  # Add the miRNA label to the plot.
  xy <- GetPlotFractionalCoords(0.95, 0.95)
  text(xy[1], xy[2], label=mirna, adj=1)
  # Add the site label to the plot.
  xy <- GetPlotFractionalCoords(0.95, 0.90)
  text(xy[1], xy[2], label=site, adj=1)
  # Add the experiment label to the plot.
  xy <- GetPlotFractionalCoords(0.05, 0.95)
  text(xy[1], xy[2], label=exp_label, adj=0)

  # Add pairwise correlation between the x and y values.
  xy <- GetPlotFractionalCoords(0.95, 0.05)
  inds_use <- (!is.na(x) & !is.na(y) & is.finite(x) & is.finite(y))
  AddCorrelationToPlot(x[inds_use], y[inds_use], xy[1], xy[2],
                       rsquared=FALSE, adj=1)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

CheckSiteDistributionHistogram <- function(mirna, site, xpos=xpos, ypos=ypos,
                                       combined=TRUE, tpm=TRUE, buffer=FALSE,
                                       pdf.plot=FALSE) {
  # Get the two signal data sets:
  signal1 <- ReporterCounts(mirna, "duplex", rep=1, tpm=tpm)
  signal2 <- ReporterCounts(mirna, "duplex", rep=2, tpm=tpm)
  mirna_site_inds <- grep(paste0(mirna, "_", site), rownames(signal1))
  site_counts_1 <- signal1[mirna_site_inds, ]
  site_counts_2 <- signal2[mirna_site_inds, ]
  ecdf1 <- ecdf(site_counts_1)
  ecdf2 <- ecdf(site_counts_2)

  totals <- c(sum(site_counts_1), sum(site_counts_1))
  SubfunctionCall(FigureSaveFile, width=10)
  par(mfrow=c(1, 2))
  xmin <- 1
  xmax <- max(c(signal1[mirna_site_inds, ], signal2[mirna_site_inds, ]))
  ymin <- 0
  ymax <- 1
  BlankPlot(log='x')
  x_ecdf <- seq(xmin, xmax)
  AddLogAxis(1, label="Counts")
  AddLinearAxis(2, tick.space=0.05, label.space=0.2, label="Cumulative frequency")
  lines(x_ecdf, ecdf1(x_ecdf), lwd=1, col="blue")
  lines(x_ecdf, ecdf2(x_ecdf), lwd=1, col="red")


  mirna_removed <- setdiff(kMirnas, mirna)
  for (mirna_bg in mirna_removed) {
    signal1 <- ReporterCounts(mirna_bg, "duplex", 1, tpm=tpm)
    signal2 <- ReporterCounts(mirna_bg, "duplex", 2, tpm=tpm)
    mirna_site_inds <- grep(paste0(mirna, "_", site), rownames(signal1))
    site_counts_1 <- signal1[mirna_site_inds, ]
    site_counts_2 <- signal2[mirna_site_inds, ]
    ecdf1 <- ecdf(site_counts_1)
    ecdf2 <- ecdf(site_counts_2)
    totals <- c(totals, sum(site_counts_1), sum(site_counts_1))
    lines(x_ecdf, ecdf1(x_ecdf), lwd=1, col="blue", lty=2)
    lines(x_ecdf, ecdf2(x_ecdf), lwd=1, col="red", lty=2)

  }
  xy <- GetPlotFractionalCoords(0.9, 0.95, log='x')
  text(xy[1], xy[2], label=mirna)
  xy <- GetPlotFractionalCoords(0.9, 0.90, log='x')
  text(xy[1], xy[2], label=site)

  xmin <- 0
  xmax <- 12
  ymin <- 1
  ymax <- 10^ceiling(log10(max(totals)))
  BlankPlot()
  kMirnas[length(kMirnas)] <- "miR-7"
  mirna_labs <- c(mirna, mirna_removed)
  mirna_cols <- kMirnaColors[mirna_labs]
  mirna_cols_2 <- ConvertRColortoRGB(mirna_cols, alpha=0.8)
  cols <- c(rbind(mirna_cols, mirna_cols_2))
  AddLinearAxis(1, alt_lab=mirna_labs, alt_lab_pos=(seq(length(kMirnas)) - 1)*2 + 2, label="",
                angled=TRUE, noline=TRUE)

  AddLinearAxis(2, tick.space=ymax/50, label.space=ymax/10,
                label="Total counts")
  x_l <- seq(xmin, xmax - 1)
  y_b <- 0
  x_r <- x_l + 1
  y_t <- totals
  rect(x_l, y_b, x_r, y_t, col=cols, border=NA)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

PlotCountHistogram <- function(mirna, pdf.plot=FALSE) {
  # Load the count data.
  counts_full <- read.table("ReporterScreen/rep2_counts.txt")
  counts_full_1 <- read.table("ReporterScreen/rep1_counts.txt")
  colnames(counts_full) <- gsub("\\.", "-", colnames(counts_full))
  # Remove the variants for which there are no counts in all seven samples:
  counts <- counts_full[which(rowSums(counts_full) != 0), ]
  num_zero_counts_all <-   length(which(rowSums(counts_full + counts_full_1) == 0))
  # Define the min and max of the count matrix, in order to specify the breaks
  # in the histogram function, and the mean, to put in the plot.
  c_min <- min(counts[, mirna])
  c_max <- max(counts[, mirna])
  c_mean <- mean(counts[, mirna])
  num_zeros <- length(which(counts[, mirna] == 0))
  num_means <- length(which(counts[, mirna] == round(c_mean, 0)))
  # Simulate a poisson distribution with the same mean.
  sample_counts <- rpois(nrow(counts), c_mean)
  # Set up the plot.
  pdf.plot_original <- pdf.plot
  if (class(pdf.plot) == "character") pdf.plot <- paste0(pdf.plot, "_left")
  SubfunctionCall(FigureSaveFile)
  # par(mfrow=c(1, 2))
  hist_counts <- hist(counts[, mirna], breaks=seq(c_min , c_max + 1),
                      right=FALSE, plot=FALSE)
  plot(hist_counts, main=NA, xlab=NA, ylab=NA)
  segments(0, num_zeros, 0.3*c_max, num_zeros, lty=2, xpd=NA)
  text(0.3*c_max, num_zeros, label=sprintf("%s zero-count variants (%1.2f%%)",
                                           num_zeros,
                                           100*num_zeros/nrow(counts)),
       xpd=NA, adj=0)
  title(xlab="Counts", line=1.5)
  title(ylab="Frequency", line=2)
  print(hist_counts$counts[1:10])
  segments(c_mean, 0, c_mean, 0.8*num_zeros, lty=2)
  segments(c_mean, 0.8*num_zeros, 0.3*c_max, lty=2)
  text(c_mean, 0.85*num_zeros, label=sprintf("Mean: %s", round(c_mean, 0)),
       xpd=NA, adj=0)

  text(0.3*c_max, 0.8*num_zeros, label=sprintf("%s mean-count variants",
                                               num_means),
       xpd=NA, adj=0)

  y_max <- max(hist_counts$counts)
  text(0.95*c_max, y_max*0.05, label=mirna, adj=1)

  hist_sim <- hist(sample_counts, breaks=seq(c_min - 1, c_max) + 0.5, right=FALSE,
                   plot=FALSE)
  if (class(pdf.plot) == "character") {
    dev.off()
  }  
  if (class(pdf.plot) == "character") pdf.plot <- paste0(pdf.plot_original,
                                                         "_right")
  SubfunctionCall(FigureSaveFile)
  plot(hist_sim, main=NA, xlab=NA, ylab=NA)
  title(xlab="Counts", line=1.5)
  title(ylab="Frequency", line=2)
  y_max <- max(hist_sim$counts)

  text(0.5*c_max, y_max*0.95,
       label="Simulated Poisson distribution", adj=0.5)
  text(0.8*c_max, y_max*0.05, label=mirna)
  if (class(pdf.plot) == "character") {
    dev.off()
  }  
}

PlotLuciferaseExperiment <- function(date, swap=TRUE, height=7, width=5,
                                     pdf.plot=FALSE) {
  path <- file.path("ReporterScreen", "LuciferaseExperiments",
                    "All_luciferase_data.txt")
  file <- as.data.frame(read.table(path, header=TRUE, row.names=1,
                                   stringsAsFactors=FALSE, sep="\t",
                                   comment.char=""),
                        stringsAsFactors=FALSE)
  file <- file[which(file$Date == date), ]
  file <- file[, 3:ncol(file)]
  # Make Renilla / Firefly Ratio column:
  file$Rel <- file$Renilla/file$Firefly
  if (date == "190227" & swap) {
    for (row_ind in 1:nrow(file)) {
      if (file$Site[row_ind] == "miR-1") {
        file$Sites[row_ind] <- "miR-124"
      } else {
        file$Sites[row_ind] <- "miR-1"
      }
    }
  }
  if (date == "190302") {
    file$Rel[6] <- NA
  }
  if (date == "190415") {
    file$Rel[6] <- NA
    file$Rel[11] <- NA
  }
  rep_columns <- grep("Rep", colnames(file))
  unique_col_inds <- sapply(1:(min(rep_columns) - 1), function(col_ind) {
  	length(unique(file[, col_ind]))
  	})
  var_inds <- which(unique_col_inds != 1)
  print(var_inds)
  print(colnames(file)[var_inds])
  if (date == "190415") {
    vec_pUC_temp <- file[, var_inds[3]]
    print(vec_pUC_temp)
    vec_pUC_temp_new <- vec_pUC_temp
    vec_pUC_temp_new[which(vec_pUC_temp == 0.6)] <- 1
    file[, var_inds[3]] <- vec_pUC_temp_new 
  }
  if (date == "190302") {
    var_inds <- var_inds[-3]
  }
  if (date == "190317") {
    average_df <- matrix(c(sapply(unique(file[, var_inds[2]]), function(ind_2) {
        sapply(unique(file[, var_inds[1]]), function(ind_1) {
          bool_1 <- file[, var_inds[1]] == ind_1
          bool_2 <- file[, var_inds[2]] == ind_2
          ind_all <- which(bool_1 & bool_2)
          print(ind_all)
          vals <- file$Rel[ind_all]
          row_mean <- mean(vals, na.rm=TRUE)
          row_sd <- sd(vals, na.rm=TRUE)/sqrt(sum((!is.na(vals))))
          out <- c(ind_1, ind_2, row_mean, row_sd)
          out
        })
      })), byrow=TRUE, nrow=8, ncol=4)
  } else {
    average_df <- matrix(c(sapply(unique(file[, var_inds[2]]), function(ind_2) {
      sapply(unique(file[, var_inds[3]]), function(ind_3) {
        sapply(unique(file[, var_inds[1]]), function(ind_1) {
          message("ind_1")
          message(ind_1)
          message("ind_2")
          message(ind_2)
          message("ind_3")
          message(ind_3)
          bool_1 <- file[, var_inds[1]] == ind_1
          bool_2 <- file[, var_inds[2]] == ind_2
          bool_3 <- file[, var_inds[3]] == ind_3
          ind_all <- which(bool_1 & bool_2 & bool_3)
          vals <- file$Rel[ind_all]
          row_mean <- mean(vals, na.rm=TRUE)
          row_sd <- sd(vals, na.rm=TRUE)/sqrt(sum((!is.na(vals))))
          out <- c(ind_1, ind_2, ind_3, row_mean, row_sd)
          out
          })
        })
      })), byrow=TRUE, nrow=8, ncol=5)  
    }


  average_df <- data.frame(average_df, stringsAsFactors=FALSE)
  colnames(average_df) <- c(colnames(file)[var_inds], "Mean", "SD")
  average_df$Mean <- as.numeric(average_df$Mean)
  average_df$SD <- as.numeric(average_df$SD)
  print(average_df)
  if (date == "190415") {
    average_df[["pUC..ug."]][7:8] <- c("0.6", "0.6")
    print(average_df[["pUC..ug."]][7:8])
  }
  print(average_df)
  norm_df <- data.frame(t(sapply(1:nrow(average_df), function(row_ind) {
    if (date == "190224") {
      mirna <- average_df$miRNA[row_ind]
      utr <- average_df$UTR[row_ind]
      norm_ind <- which(average_df$miRNA == mirna &
                       average_df$UTR == utr &
                       average_df$Sites != mirna)
    } else if (date == "190227") {
      norm_ind <- ceiling(row_ind/2)*2 - 1
    } else {
      norm_ind <- ceiling(row_ind/2)*2
    }
    mean_row <- average_df$Mean[row_ind]
    sd_row <- average_df$SD[row_ind]
    mean_ctrl <- average_df$Mean[norm_ind]
    sd_ctrl <- average_df$SD[norm_ind]
    mean_norm <- mean_row/mean_ctrl
    if (mean_norm == 1) {
      sd_norm <- 0
    } else {
      sd_norm <- mean_norm*sqrt((sd_row/mean_row)^2 + (sd_ctrl/mean_ctrl)^2)
    }
    out <- average_df[row_ind, ]
    if (date == "190317") {
      out[3] <- mean_norm
      out[4] <- sd_norm

    } else {
      out[4] <- mean_norm
      out[5] <- sd_norm          
    }
    out
  })))
  for (ind_col in 1:ncol(norm_df)) {
    norm_df[, ind_col] <- unlist(norm_df[, ind_col])
  }
  if (date != "190317") {
    norm_df <- norm_df[order(norm_df[, 1]), ]  
  }
  print(norm_df)
  if (date == "190224") {
    norm_df$miRNA_UTR <- paste(norm_df$miRNA, norm_df$UTR, sep="_")
    norm_df$UTR_site <- paste(norm_df$UTR, ", ", norm_df$Sites, " sites", sep="")    
    list_temp <- list(norm_df$UTR_site, norm_df$miRNA)
  } else if (date == "190317") {
    norm_df$Site <- paste0(norm_df$Site, " sites")
    norm_df$div <- norm_df[, 2]
    list_temp <- list(norm_df$Site, norm_df$div)
  } else {
  	norm_df$Site <- paste0(norm_df$Site, " sites")
    norm_df$div <- paste(norm_df[, 2], "\n", norm_df[, 3], sep=" ")
    list_temp <- list(norm_df$Site, norm_df$div)
  }
  print(list_temp)
  norm_df$Mean <- as.numeric(norm_df$Mean)
  tabbedMeans <- tapply(norm_df$Mean,
                        list_temp,
                        function(x) c(x = x))
  tabbedSDs <- tapply(norm_df$SD,
                        list_temp,
                        function(x) c(x = x))

  if (date == "190415") {
    tabbedMeans <- tabbedMeans[, c(1, 3, 2, 4)]
    tabbedSDs <- tabbedSDs[, c(1, 3, 2, 4)]
  } else if (date == "190317") {
    tabbedMeans <- tabbedMeans[c(1, 2), c(2, 1, 4, 3)]
    tabbedSDs <- tabbedSDs[c(1, 2), c(2, 1, 4, 3)]        
  } else if (date == "190224") {
    tabbedMeans <- tabbedMeans[c(2, 1, 4, 3), ]
    tabbedSDs <- tabbedSDs[c(2, 1, 4, 3), ]    
  } else if (date %in% c("190227")) {
    tabbedMeans <- tabbedMeans[, c(2, 4, 1, 3)]
    tabbedSDs <- tabbedSDs[, c(2, 4, 1, 3)]
  } else if (date %in% c("190228", "190303")) {
    tabbedMeans <- tabbedMeans[, c(3, 4, 1, 2)]
    tabbedSDs <- tabbedSDs[, c(3, 4, 1, 2)]
  } else if (date == "190302") {
    tabbedMeans <- tabbedMeans[, c(1, 3, 2, 4)]
    tabbedSDs <- tabbedSDs[, c(1, 3, 2, 4)]
  }
  if (date != "190224" && date != "190302" && date != "190415") {
    tabbedMeans <- tabbedMeans[c(2, 1), ]
    tabbedSDs <- tabbedSDs[c(2, 1), ]
  }

  if (date == "190224") {
    legend_title <- "3' UTR:"
    names_arg <- c("", "")
  } else {
    legend_title <- "3' UTR sites:"
    names_arg <- c("", "", "", "")
  }

  SubfunctionCall(FigureSaveFile)
  par(mar=c(10, 6, 2, 2))

  ymin <- 0
  ymax <- 1.8
  xmin <- 0
  xmax <- 12
  barCenters <- barplot(height=tabbedMeans, beside=TRUE, border=NA, names.arg=names_arg,
                        ylim=c(ymin, ymax), xlim=c(xmin, xmax), legend.text=TRUE, axes=FALSE,
                        args.legend=list(x="topleft", title=NULL, border=NA, box.lty=0))
  arrows(barCenters, tabbedMeans - 1.96*tabbedSDs,
         barCenters, tabbedMeans + 1.96*tabbedSDs, lwd=par()$cex, length=0.025,
         angle=90, code=3)

  text_labels <- c("miRNA:", "Duplex (pmol):", "RNAiMAX (uL):", "Renilla (ng):",
                "Firefly (ng):", "TRL (ng):", "pUC (ug):", "L2K (uL):", "Cells (1e6):")
  y_top <- GetPlotFractionalCoords(0, -0.05)[2]
  y_bottom <- GetPlotFractionalCoords(0, -0.5)[2]
  y_vals <- seq(y_top, y_bottom, length.out=length(text_labels))
  y_val_dif <- y_vals[1] - y_vals[2]
  text(0, y_vals,
       labels=text_labels, xpd=NA, adj=1)
  AddLinearAxis(2, label.space=0.2, tick.space=0.2,
                label="Relative Luciferase activity")

  cols_order <- c(1, 5, 6, 7, 8, 9, 10, 11, 12)
  if (date == "190224") {
    full_width <- c(1, 10)
  } else {
    full_width <- c(1, 12)
  }
  if (date == "190415") {
    file[, var_inds[3]] <- vec_pUC_temp
  }
  for (ind in seq(1,length(cols_order))) {
    print(colnames(file)[cols_order[ind]])
    all_conditions <- file[, cols_order[ind]]
    print(all_conditions)
    all_conditions[is.na(all_conditions)] <- "NA"
    conditions <- unique(all_conditions)
    print(conditions)
    if (length(conditions) == 1) {
      condition <- unique(conditions)
      if (is.na(condition)) {
        condition <- "NA"
      }
      x_ind <- mean(full_width)
      text(x_ind, y_vals[ind], labels=as.character(conditions), xpd=NA)
      segments(x0=full_width[1], y0=y_vals[ind] + y_val_dif/2, x1=full_width[2], xpd=NA)
    } else {
      box_r <- c(floor(full_width[2]/2), full_width[2])
      box_l <- c(full_width[1], floor(full_width[2]/2 + 1))
      mids <- (box_l + box_r)/2
      if (date == "190415" & conditions == c("0", "1", "0.6")) {
        box_l <- c(box_l[1], box_l[2], ceiling(mids[2]))
        box_r <- box_l + c(5, 2, 2)
        mids <- (box_l + box_r)/2
      } else if ((all_conditions[7] != all_conditions[1]) ||
          (date == "190227" & all_conditions[13] != all_conditions[1])) {
        box_l <- c(box_l[1], ceiling(mids[1]), box_l[2], ceiling(mids[2]))
        box_r <- box_l + 2
        mids <- (box_l + box_r)/2
        conditions <- rep(conditions, 2)
        if (date == "190228") {
          conditions <- rev(conditions)
        }
      }
        text(mids, y_vals[ind], labels=conditions, xpd=NA)
        segments(x0=box_l, y0=y_vals[ind] + y_val_dif/2, x1=box_r, xpd=NA)
    }
  }

  	trans_text <- unique(file$Format)
    xy <- GetPlotFractionalCoords(0.05, 1.05)
  text(xy[1], xy[2], labels=trans_text, xpd=NA, adj=0)
	if (class(pdf.plot) == "character") {
    dev.off()
  }  

}


PlotFlankKdsVsRepression(sitelist="topcanonical")

break




## FIGURES FOR RBNS EQUILIBRIUM PAPER ##########################################
MakeFigure1 <- function(uniq=FALSE) {
  message("Making Fig. 1")
  PlotEquilSiteWithInput("miR-1", 7, sitelist="canonical", combined=FALSE,
                         buffer=TRUE, pdf.plot="1.C")
  message("Done 1.C")
  PlotSiteEnrichments("miR-1", sitelist="canonical", combined=FALSE,
                      buffer=TRUE, write_kds=TRUE, pdf.plot="1.D")
  message("Done 1.D")
  PlotSiteEnrichments("miR-1", combined=FALSE, remove_sites=FALSE, buffer=TRUE,
                      pdf.plot="1.E")
  message("Done 1.E")
  PlotSiteKds("miR-1", combined=FALSE, remove_sites=FALSE,
              buffer=TRUE, pdf.plot="1.F")
  message("Done 1.F")
  PlotSiteOccupancy("miR-1", combined=FALSE, buffer=TRUE, pdf.plot="1.G")
  message("Done Fig. 1")
}

MakeFigure2 <- function() {
  message("Making Fig. 2")
  PlotSiteKds("let-7a", adjusted_height=TRUE, pdf.plot="2.Ai")
  PlotSiteOccupancy("let-7a", adjusted_height=TRUE, pdf.plot="2.Aii")
  PlotSiteKds("miR-155", adjusted_height=TRUE, pdf.plot="2.Bi")
  PlotSiteOccupancy("miR-155", adjusted_height=TRUE, pdf.plot="2.Bii")
  PlotSiteKds("miR-124", collapse_AA=TRUE, adjusted_height=TRUE,
              pdf.plot="2.Ci")
  PlotSiteKds("miR-124", collapse_AA=TRUE, adjusted_height=TRUE,
              pdf.plot="2.Ci")
  PlotSiteOccupancy("miR-124", collapse_AA=TRUE, compcorrect=FALSE,
                    adjusted_height=TRUE, pdf.plot="2.Cii")

  message("Done Fig. 2")
}

MakeFigure3 <- function() {
  message("Making Fig. 3")
  PlotPositionalKds(pdf.plot="3.A")
  Plot8merVs7merCombined(pdf.plot="3.B")
  PlotSiteKdsVsSPS(pdf.plot="3.C")
  PlotSiteKdsVsRepression("miR-1", combined=FALSE, buffer=TRUE, pdf.plot="3.D")
  PlotSiteKdsVsRepression("let-7a", pdf.plot="3.E")
  PlotSiteKdsVsRepression("miR-155", pdf.plot="3.F")
  PlotSiteKdsVsRepression("miR-124", pdf.plot="3.G")
  PlotSiteKdsVsRepression("lsy-6", pdf.plot="3.H")
  PlotSiteKdsVsRepression("miR-7-23nt", experiment="equilibrium2_nb",
                          pdf.plot="3.I")
  message("Done Fig. 3")
}

MakeFigure4 <- function() {
  message("Making Fig. 4")
  PlotSiteFlankEnrichments("miR-1", "8mer", combined=TRUE, combined_site=FALSE,
                           buffer=TRUE, pdf.plot="4.A")
  PlotSiteFlankKds("miR-1", combined=TRUE, combined_site=FALSE, buffer=TRUE,
                   pdf.plot="4.B")
  PlotFlankLinModel(pdf.plot="4.C_left")
  PlotFlankLinModelCoefficients(pdf.plot="4.C_right")
  PlotStructureVsFlankingKds("miR-1", "8mer", combined=TRUE, buffer=TRUE,
                             pdf.plot="4.D")
  message("Done Fig. 4")
}

MakeSupplementaryFigure1 <- function() {
  message("Making fig. S1")
  PlotAgoPrepPurity(pdf.plot="S1.A", no_marker=FALSE, no_adapter=FALSE)
  PlotMiR1_KmersCorrelation(pdf.plot="S1.B", kmer_len=9)
  PlotEnrichmentsAgainstKds("miR-1", sitelist="canonical", buffer=TRUE,
                            combined=FALSE, pdf.plot="S1.C")
  PlotPositionalEnrichment("miR-1", buffer=TRUE, pdf.plot="S1.Di") # checked
  PlotPositionalEnrichment("let-7a", pdf.plot="S1.Dii") # checked
  PlotPositionalEnrichment("miR-155", pdf.plot="S1.Diii") # checked
  PlotPositionalEnrichment("miR-124", pdf.plot="S1.Div") # checked
  PlotPositionalEnrichment("lsy-6", pdf.plot="S1.Dv") # checked
  PlotPositionalEnrichment("miR-7-23nt", exp="equilibrium2_nb", pdf.plot="S1.Dvi")
  PlotPositionalEnrichment("miR-155", sites=k3PSites, pdf.plot="S1.Ei") # checked
  PlotPositionalEnrichment("miR-124", sites=k3PSites, pdf.plot="S1.Eii") # checked
  PlotPositionalEnrichment("lsy-6", sites=k3PSites, pdf.plot="S1.Eiii") # checked  break
  PlotWorstSiteKdCrossValScatter("miR-1", sitelist="canonical", combined=FALSE,
                                 buffer=TRUE, pdf.plot="S1.F") # checked
  PlotSiteKdCrossValMatrix("miR-1", sitelist="canonical", combined=FALSE,
                           buffer=TRUE, pdf.plot="S1.G") # checked
  PlotSalomonComparison(pdf.plot="S1.H") # checked
  message("Done fig. S1")
}

MakeSupplementaryFigure2 <- function() {
  message("Making fig. S2.")
  PlotCompetitorOligoSiteSchematic(pdf.plot="S2.A")
  PlotSiteKds("lsy-6", adjusted_height=TRUE, pdf.plot="S2.Bi")
  PlotSiteOccupancy("lsy-6", sitelist="resubmissionfinal",  adjusted_height=TRUE,
                    pdf.plot="S2.Bii")
  PlotSiteKds("miR-7-23nt", exp="equilibrium2_nb", adjusted_height=TRUE,
              pdf.plot="S2.Ci")
  PlotSiteOccupancy("miR-7-23nt", exp="equilibrium2_nb", adjusted_height=TRUE,
                    pdf.plot="S2.Cii")
  message("Done fig. S2.")
}

MakeSupplementaryFigure3 <- function() {
  message("Making fig. S3.")
  PlotBaekKds("miR-1", combined=FALSE, buffer=TRUE, pdf.plot="S3.A")
  PlotBaekKds("let-7a", pdf.plot="S3.B")
  PlotBaekKds("miR-155", pdf.plot="S3.C")
  PlotBaekKds("miR-124", pdf.plot="S3.D")
  PlotBaekKds("lsy-6", pdf.plot="S3.E")
  PlotBaekKds("miR-7-23nt", experiment="equilibrium2_nb", pdf.plot="S3.F")
  message("Done fig. S3.")
}

MakeSupplementaryFigure4 <- function() {
  message("Making fig. S4.")
  PlotBulgeKds("miR-1", combined=FALSE, buffer=TRUE, pdf.plot="S4.B")
  PlotDelKds("miR-1", , combined=FALSE, buffer=TRUE, pdf.plot="S4.C")
  PlotBulgeKds("let-7a", pdf.plot="S4.D")
  PlotDelKds("let-7a", pdf.plot="S4.E")
  PlotBulgeKds("miR-155", pdf.plot="S4.F")
  PlotDelKds("miR-155", pdf.plot="S4.G")
  PlotBulgeKds("miR-124", pdf.plot="S4.H")
  PlotDelKds("miR-124", pdf.plot="S4.I")
  PlotBulgeKds("lsy-6", pdf.plot="S4.J")
  PlotDelKds("lsy-6", pdf.plot="S4.K")
  PlotBulgeKds("miR-7-23nt", experiment="equilibrium2_nb", pdf.plot="S4.L")
  PlotDelKds("miR-7-23nt", experiment="equilibrium2_nb",  pdf.plot="S4.M")
  message("Done fig. S4.")
}


MakeSupplementaryFigure5 <- function() {
  message("Making fig. S5.")
  PlotReporterAssaySchematicLetters(pdf.plot="S5.A")
  PlotReporterAssayKdVsL2fc("miR-1", pdf.plot="S5.Bi")
  PlotReporterAssayKdVsL2fc("let-7a", pdf.plot="S5.Bii")
  PlotReporterAssayKdVsL2fc("miR-155", pdf.plot="S5.Biii")
  PlotReporterAssayKdVsL2fc("miR-124", pdf.plot="S5.Biv")
  PlotReporterAssayKdVsL2fc("lsy-6", pdf.plot="S5.Bv")
  PlotReporterAssayKdVsL2fc("miR-7", pdf.plot="S5.Bvi")
  message("Done fig. S5.")
 }

MakeSupplementaryFigure6 <- function() {
  message("Making fig. S6.")
  PlotSiteFlankKds("let-7a", adjusted_height=TRUE, width=7.3,
                   pdf.plot="S6.A")
  PlotSiteFlankKds("miR-155", adjusted_height=TRUE, width=7.3,
                   pdf.plot="S6.B")
  PlotSiteFlankKds("miR-124", adjusted_height=TRUE, width=7.3,
                   pdf.plot="S6.C")
  PlotSiteFlankKds("lsy-6", adjusted_height=TRUE, width=7.3,
                   pdf.plot="S6.D")
  PlotSiteFlankKds("miR-7-23nt", experiment="equilibrium2_nb", width=7.3,
                   adjusted_height=TRUE, pdf.plot="S6.E")
  PlotFlankKdsVsRepression(sitelist="topcanonical", pdf.plot="S6.F")
  PlotAllSiteInputDistribution("miR-1", "8mer", combined=TRUE, buffer=TRUE,
                               pdf.plot="S6.G")
  PlotFlankPlCorrelationWindowsNew("miR-1", "8mer",
                                   buffer=TRUE, pdf.plot="S6.H")
  PlotFlanksSamplePl("miR-1", "8mer", "0.4", buffer=TRUE, matchdist=TRUE,
                     pdf.plot="S6.I")
  message("Done fig. S6.")
 }

MakeRefereeResponseFigure <- function() {
  message("Making referee response fig.")
  PlotFlowthrough("let-7a", pdf.plot="R.A")
  PlotFlowthrough("let-7a", pdf.plot="R.B", log_y=FALSE)
  PlotFlowthrough("let-7a", pdf.plot="R.C",
                  enrich_against_enrich=TRUE)
  message("Done referee response fig.")
 }



 
# MakeFigure1() # checked
# MakeFigure2() # checked
# MakeFigure3() # checked
# MakeFigure4() # checked
# MakeSupplementaryFigure1() # checked
# MakeSupplementaryFigure2()
# MakeSupplementaryFigure3()
# MakeSupplementaryFigure4()
# MakeSupplementaryFigure5()
# MakeSupplementaryFigure6()


# MakeRefereeResponseFigure()

