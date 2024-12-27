rm()
source("general/general.R")
source("general/ModelingFunctions.R")
source("general/PlotFunctions.R")

AGO_mir_label <- "[AGO2-miRNA] (pM)"
occupancy_mirna_seq_col <- "#958872"
occupancy_schematic_bg_col <- "#989898"
occupancy_schematic_bg_cex <- 0.45


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
      path <- paste0("McGearyLinEtAl_2019/Figure_", figure, "/", figurename, "_raw.pdf")
      print(path)
    }
    pdf(file=path, height=2.2/5*height, width=2.2/5*width, useDingbats=FALSE)
    par(kPDFParameters)
  } else {
    dev.new(xpos=xpos, ypos=ypos, height=height, width=width)
    par(kPlotParameters)    
  }  
}



# source("general/GenericFigures.R")

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
  print(23)
  SubfunctionCall(FigureSaveFile)
  print(24)
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
  # Note: The "remove_sites" boolean is whether or not to remove sites that have
  # no obvious name (e.g., CACACACA, etc.)
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

  # Get the model curves.
  model <- SubfunctionCall(EquilSingleSiteModelFreq, pars=pars.model,
                           A.dil=A.dil.model)
  model.R <- EquilEnrichments(model, l)


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
  print("here")
  cols <- c(cols, None="black")
  if (flowthrough) {
  	segments(A.pM.data, ymin, A.pM.data, ymax, lty=2, col="gray")
  }
  # Add the points and lines:
  sapply(rownames(data), function(site) {
  	if (!flowthrough) {
      print(A.pM.data)
      print(data.R[site, ])
	    Points(A.pM.data, data.R[site, ], col=cols[site], line=datalines)
		}
		if (!(datalines) & modellines) {
      print(A.pM.model)
      print(model.R[site, ])

	    lines(A.pM.model, model.R[site, ], col=cols[site], xpd=NA)    		
		}
  })
  print(262)
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
  print(kMirnas)
  print(mirna.sites)
  names(mirna.sites) <- c(kMirnas[1:6], "miR-7-22nt", "miR-7-24nt", "miR-7-25nt") # Had to modify this, kMirnas list is current 9 miRNAs long.
  xy <- GetPlotFractionalCoords(fx=0.65, fy=0, log='x', inv='x')
  print(mirna.sites[[mirna]])
  print(mirna.sites[[mirna]])
  print(legend.colors)
  print(legend.colors[mirna.sites[[mirna]]])
  print(legend.names[mirna.sites[[mirna]]])
  if (!(added.text)) {
    Legend(xy, legend=legend.names[mirna.sites[[mirna]]],
           col=legend.colors[mirna.sites[[mirna]]],
           yjust=0)    
  }
  if (class(pdf.plot) == "character") {
    dev.off()
  }
  print("removed sites:")
  print(removed_sites)
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
  # file_out <- sprintf("aesthetic_tables_for_kathy/%s_site_pastel_colors.txt",
  #                     mirna.str)
  # write.table(file_out, x=df.global, quote=FALSE, row.names=FALSE, sep="\t",
  #             col.names=TRUE)
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

  # kMirnas <- kMirnas
  tick <- 0
  max_kds <- 0
  none_kd_offset_start <- 0.8
  none_kd_offset <- none_kd_offset_start
  df.global <- data.frame(mirnas=c(), color=c())
  sapply(kMirnas[1:6], function(mirna) {
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
  # file_out <- sprintf("aesthetic_tables_for_kathy/all_mirna_colors.txt")
  # write.table(file_out, x=df.global, quote=FALSE, row.names=FALSE, col.names=TRUE,
  #             sep="\t")
  # Add text to plot:
  par(xpd=NA)
  print("here")
  xy_alt <- GetPlotFractionalCoords(fx=-0.01, fy=0.4, log='x', inv='x')
  text(x=xy_alt[1], y=max_kds - seq(max_kds) + 1, labels=max_kd_labels, adj=0,
       col="black")
  # Add legend to plot:
  xy <- GetPlotFractionalCoords(fx=0.75, fy=0.4, log='x', inv='x')
  kMirnas[length(kMirnas)] <- "miR-7"
  Legend(xy, legend=kMirnas[1:6], col=kMirnaColors[1:6])
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

# 3B____________________________________________________________________________
Plot8merVs7merCombined <- function(experiment="equilibrium", n_constant=5,
                                   sitelist="resubmissionfinal", buffer=FALSE,
                                   singleonly=TRUE, dG.table=3, height=4,
                                   width=5, pdf.plot=FALSE) {
  # dG.df <- read.table(file=paste0("canonical_sites_mfe_", dG.table, ".txt"))
  # colnames(dG.df) <- gsub("(^.*)\\.(.*$)", colnames(dG.df), replace="\\1-\\2")
  R <- 1.987e-3 # in kcal K-1 mol-1
  T <- 310.15 # in K
  total_x <- c(1)
  total_y <- c(1)
  print(1157)
  SubfunctionCall(FigureSaveFile)
  xmin <- 0
  xmax <- length(kMirnas[1:6])*4 + 7
  ymin <- -2.7
  ymax <- 0
  BlankPlot(inv="y")
  ymin <- -2.5
  mirna_labs <- kMirnas[1:6]
  mirna_labs[length(mirna_labs)] <- "miR-7"
  AddLinearAxis(1, alt_lab=mirna_labs, alt_lab_pos=(seq(length(kMirnas[1:6])) - 1)*4 + 2, label="",
                angled=TRUE, noline=TRUE)

  AddLinearAxis(2, tick.space=0.5, label.space=1,
                label=expression(Delta*Delta*italic(G)~"(kcal/mol)"))

  x_start <- 0
  print(1173)
  segments(x0=xmin, y0=-R*T*log(c(2, 10, 50)), x1=length(kMirnas[1:6])*4 + 0.5,
           lwd=0.25, xpd=NA)
  text(rep(length(kMirnas[1:6])*4 + 0.75, 3), -R*T*log(c(2, 10, 50)),
       labels=c("2-fold greater\nbinding affinity",
                "10-fold greater\nbinding affinity",
                "50-fold greater\nbinding affinity"),
       cex=0.8, adj=c(0, 0.5), xpd=NA)
  kds.all <- sapply(kMirnas[1:6], function(mirna) {
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
    print(1204)
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
  print(1215)
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
  dG.df <- read.table("general/miRNA_canonical_site_energies.txt")
  colnames(dG.df) <- gsub("(^.*)\\.(.*$)", colnames(dG.df), replace="\\1-\\2")
  print(dG.df)
  kMirnas <- kMirnas[1:6]
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
  Kd.matrix <- sapply(kMirnas[1:6], function(mirna) {
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
  dG.matrix <- sapply(kMirnas[1:6], function(mirna) {
    if (mirna == "miR-7-23nt") {
      mirna <- "miR-7"
    }
    print(dG.df[, mirna])
    as.numeric(dG.df[kSeedSites, mirna])
  })
  print(dG.matrix)
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
                                    new=FALSE, old=FALSE, cat_colors=FALSE,
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
  print(CombinedSiteAndFlankColors)
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




## FIGURES FOR RBNS EQUILIBRIUM PAPER ##########################################
MakeFigure1 <- function(uniq=FALSE) {
  ## To make this figure need run these commands in the command line:
  ## make PreprocessAllEquilibrium
  ## make AssignSites mirna=miR-1 exp=equilibrium n_constant=5 buffer=1 sitelist=canonical
  ## python SolveForKds/MakeSiteCountTable.py miR-1 equilibrium 5 canonical -buffer
  message("Making Fig. 1")
  PlotEquilSiteWithInput("miR-1", 7, sitelist="canonical", combined=FALSE,
                         buffer=TRUE, pdf.plot="1.C")
  message("Done 1.C")
  ## python SolveForKds/MakeMultiSiteCountTable.py miR-1 equilibrium 5 canonical -buffer3p
  ## Rscript SolveForKds/FitSiteKds.R miR-1 equilibrium 5 canonical -buffer3p -nocombI -single
  PlotSiteEnrichments("miR-1", sitelist="canonical", combined=FALSE,
                      buffer=TRUE, write_kds=TRUE, pdf.plot="1.D")
  message("Done 1.D")
  ## make AssignSites mirna=miR-1 exp=equilibrium n_constant=5 buffer=1 sitelist=resubmissionfinal
  ## python SolveForKds/MakeSiteCountTable.py miR-1 equilibrium 5 resubmissionfinal -buffer
  ## python SolveForKds/MakeMultiSiteCountTable.py miR-1 equilibrium 5 resubmissionfinal -buffer
  ## Rscript SolveForKds/FitSiteKds.R miR-1 equilibrium 5 resubmissionfinal -buffer -nocombI -single
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
  ## make AssignSites mirna=let-7a exp=equilibrium n_constant=5 sitelist=resubmissionfinal
  ## python SolveForKds/MakeSiteCountTable.py let-7a equilibrium 5 resubmissionfinal 
  ## python SolveForKds/MakeMultiSiteCountTable.py let-7a equilibrium 5 resubmissionfinal 
  ## Rscript SolveForKds/FitSiteKds.R let-7a equilibrium 5 resubmissionfinal -single
  ## make AssignSites mirna=lsy-6 exp=equilibrium n_constant=5 sitelist=resubmissionfinal
  ## python SolveForKds/MakeSiteCountTable.py lsy-6 equilibrium 5 resubmissionfinal 
  ## python SolveForKds/MakeMultiSiteCountTable.py lsy-6 equilibrium 5 resubmissionfinal 
  ## Rscript SolveForKds/FitSiteKds.R lsy-6 equilibrium 5 resubmissionfinal -single

  PlotSiteKds("let-7a", adjusted_height=TRUE, pdf.plot="2.Ai")
  PlotSiteOccupancy("let-7a", adjusted_height=TRUE, pdf.plot="2.Aii")
  ## make AssignSites mirna=miR-155 exp=equilibrium n_constant=5 sitelist=resubmissionfinal
  ## python SolveForKds/MakeSiteCountTable.py miR-155 equilibrium 5 resubmissionfinal 
  ## python SolveForKds/MakeMultiSiteCountTable.py miR-155 equilibrium 5 resubmissionfinal 
  ## Rscript SolveForKds/FitSiteKds.R miR-155 equilibrium 5 resubmissionfinal -single

  PlotSiteKds("miR-155", adjusted_height=TRUE, pdf.plot="2.Bi")
  PlotSiteOccupancy("miR-155", adjusted_height=TRUE, pdf.plot="2.Bii")
  ## make AssignSites mirna=miR-124 exp=equilibrium n_constant=5 sitelist=resubmissionfinal
  ## python SolveForKds/MakeSiteCountTable.py miR-124 equilibrium 5 resubmissionfinal 
  ## python SolveForKds/MakeMultiSiteCountTable.py miR-124 equilibrium 5 resubmissionfinal 
  ## Rscript SolveForKds/FitSiteKds.R miR-124 equilibrium 5 resubmissionfinal -single
  PlotSiteKds("miR-124", collapse_AA=TRUE, adjusted_height=TRUE,
              pdf.plot="2.Ci")
  PlotSiteOccupancy("miR-124", collapse_AA=TRUE, compcorrect=FALSE,
                    adjusted_height=TRUE, pdf.plot="2.Cii")

  message("Done Fig. 2")
}

MakeFigure3 <- function() {
  message("Making Fig. 3")
  ## make n_constant=5 sitelist=centered11 AssignSitesAllEquilibrium

  ## python SolveForKds/MakeSiteCountTable.py miR-1 equilibrium 5 centered11 -buffer
  ## python SolveForKds/MakeMultiSiteCountTable.py miR-1 equilibrium 5 centered11 -buffer
  ## Rscript SolveForKds/FitSiteKds.R miR-1 equilibrium 5 centered11 -buffer -nocombI -single

  ## python SolveForKds/MakeSiteCountTable.py let-7a equilibrium 5 centered11
  ## python SolveForKds/MakeMultiSiteCountTable.py let-7a equilibrium 5 centered11
  ## Rscript SolveForKds/FitSiteKds.R let-7a equilibrium 5 centered11 -single

  ## python SolveForKds/MakeSiteCountTable.py miR-155 equilibrium 5 centered11
  ## python SolveForKds/MakeMultiSiteCountTable.py miR-155 equilibrium 5 centered11
  ## Rscript SolveForKds/FitSiteKds.R miR-155 equilibrium 5 centered11 -single

  ## python SolveForKds/MakeSiteCountTable.py miR-124 equilibrium 5 centered11
  ## python SolveForKds/MakeMultiSiteCountTable.py miR-124 equilibrium 5 centered11
  ## Rscript SolveForKds/FitSiteKds.R miR-124 equilibrium 5 centered11 -single

  ## python SolveForKds/MakeSiteCountTable.py lsy-6 equilibrium 5 centered11
  ## python SolveForKds/MakeMultiSiteCountTable.py lsy-6 equilibrium 5 centered11
  ## Rscript SolveForKds/FitSiteKds.R lsy-6 equilibrium 5 centered11 -single

  ## make mirna=miR-7-23nt exp=equilibrium n_constant=5 sitelist=resubmissionfinal AssignSites
  ## python SolveForKds/MakeSiteCountTable.py miR-7-23nt equilibrium2_nb 5 centered11
  ## python SolveForKds/MakeMultiSiteCountTable.py miR-7-23nt equilibrium2_nb 5 centered11
  ## Rscript SolveForKds/FitSiteKds.R miR-7-23nt equilibrium2_nb 5 centered11 -single



  PlotPositionalKds(pdf.plot="3.A")

  ## python SolveForKds/MakeSiteCountTable.py miR-7-23nt equilibrium2_nb 5 resubmissionfinal
  ## python SolveForKds/MakeMultiSiteCountTable.py miR-7-23nt equilibrium2_nb 5 resubmissionfinal
  ## Rscript SolveForKds/FitSiteKds.R miR-7-23nt equilibrium2_nb 5 resubmissionfinal -single


  Plot8merVs7merCombined(pdf.plot="3.B")

  ## python SolveForKds/MakeEnergyTable.py
  ## Rscript SolveForKds/FitSiteKds.R miR-1 equilibrium 5 resubmissionfinal -buffer -nocombI
  ## Rscript SolveForKds/FitSiteKds.R let-7a equilibrium 5 resubmissionfinal
  ## Rscript SolveForKds/FitSiteKds.R miR-155 equilibrium 5 resubmissionfinal
  ## Rscript SolveForKds/FitSiteKds.R miR-124 equilibrium 5 resubmissionfinal
  ## Rscript SolveForKds/FitSiteKds.R lsy-6 equilibrium 5 resubmissionfinal
  ## Rscript SolveForKds/FitSiteKds.R miR-7-23nt equilibrium2_nb 5 resubmissionfinal
  PlotSiteKdsVsSPS(pdf.plot="3.C")
  # Currently unavailable because repression file not in repo.
  # PlotSiteKdsVsRepression("miR-1", combined=FALSE, buffer=TRUE, pdf.plot="3.D")
  # PlotSiteKdsVsRepression("let-7a", pdf.plot="3.E")
  # PlotSiteKdsVsRepression("miR-155", pdf.plot="3.F")
  # PlotSiteKdsVsRepression("miR-124", pdf.plot="3.G")
  # PlotSiteKdsVsRepression("lsy-6", pdf.plot="3.H")
  # PlotSiteKdsVsRepression("miR-7-23nt", experiment="equilibrium2_nb",
  #                         pdf.plot="3.I")
  message("Done Fig. 3")
}

MakeFigure4 <- function() {
  message("Making Fig. 4")
  ## make AssignFlanks mirna=miR-1 exp=equilibrium n_constant=5 buffer=1 sitelist=resubmissionfinal
  ## python SolveForKds/MakeFlankCountTable.py miR-1 equilibrium 5 resubmissionfinal -buffer
  #### Rscript SolveForKds/FitSiteKds.R miR-1 equilibrium 5 resubmissionfinal -buffer
  ## Rscript SolveForKds/FitFlankKds.R miR-1 equilibrium 5 resubmissionfinal 8mer -buffer
  PlotSiteFlankEnrichments("miR-1", "8mer", combined=TRUE, combined_site=FALSE,
                           buffer=TRUE, pdf.plot="4.A")

  ## make mirna=miR-1 exp=equilibrium n_constant=5 sitelist=resubmissionfinal buffer=1 FitFlankKds
  PlotSiteFlankKds("miR-1", combined=TRUE, combined_site=FALSE, buffer=TRUE,
                   pdf.plot="4.B")
  ## make AssignFlanks mirna=let-7a exp=equilibrium n_constant=5 sitelist=resubmissionfinal
  ## make AssignFlanks mirna=miR-155 exp=equilibrium n_constant=5 sitelist=resubmissionfinal
  ## make AssignFlanks mirna=miR-124 exp=equilibrium n_constant=5 sitelist=resubmissionfinal
  ## make AssignFlanks mirna=lsy-6 exp=equilibrium n_constant=5 sitelist=resubmissionfinal
  ## make AssignFlanks mirna=miR-7-23nt exp=equilibrium2_nb n_constant=5 sitelist=resubmissionfinal

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



 
# MakeFigure1()
# MakeFigure2()
# MakeFigure3()
MakeFigure4()
# MakeSupplementaryFigure1()
# MakeSupplementaryFigure2()
# MakeSupplementaryFigure3()
# MakeSupplementaryFigure4()
# MakeSupplementaryFigure5()
# MakeSupplementaryFigure6()



