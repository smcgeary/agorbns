PlotLet7aProgrammedLibraryReps <- function(
  average_sites=FALSE, height=5, width=5, xpos=20, ypos=20, pdf.plot=FALSE
) {
  kds_1 <- read.table("ThreePrimeTargetingPaper/temp_let-7a_kds.txt", sep="\t",
                       row.names=1, header=TRUE)
  kds_2 <- read.table("ThreePrimeTargetingPaper/temp_let-7a-21nt_kds.txt", sep="\t",
                         row.names=1, header=TRUE)
  if (average_sites) {
    kds_1 <- CollapseSeedSitesTable(kds_1)
    kds_2 <- CollapseSeedSitesTable(kds_2)
    alpha_sites <- 1    
  } else {
    alpha_sites <- 0.5
  }
  #TODO fix this, such that it doesn't need to be hardcoded into the data.
  # Realistically, I should extend to 5 nucleotides outside.
  # That makes the most sense.
  rownames(kds_1) <- gsub("^6mer\\|28\\|", x=rownames(kds_1), perl=TRUE,
                          replacement="7mer-m8\\|28\\|")
  rownames(kds_2) <- gsub("^6mer\\|28\\|", x=rownames(kds_2), perl=TRUE,
                          replacement="7mer-m8\\|28\\|")

  rownames(kds_1) <- gsub("^7mer-A1\\|27\\|", x=rownames(kds_1), perl=TRUE,
                          replacement="8mer\\|27\\|")
  rownames(kds_2) <- gsub("^7mer-A1\\|27\\|", x=rownames(kds_2), perl=TRUE,
                          replacement="8mer\\|27\\|")

  sites_shared <- intersect(rownames(kds_1), rownames(kds_2))
  sites_shared_single <- grep("&", sites_shared, invert=TRUE, value=TRUE)

  kds_shared_1 <- kds_1[sites_shared, 1]
  kds_shared_2 <- kds_2[sites_shared, 1]

  kds_single_1 <- kds_1[sites_shared_single, 1]
  kds_single_2 <- kds_2[sites_shared_single, 1]

  kds_shared <- t(rbind(kds_shared_1, kds_shared_2))
  kds_single <- 10^t(rbind(kds_single_1, kds_single_2))
  rownames(kds_shared) <- sites_shared
  rownames(kds_single) <- sites_shared_single
  SubfunctionCall(FigureSaveFile2)
  xmin=0.001
  xmax=3
  ymin=0.001
  ymax=3
  BlankPlot(log="xy")
  inds_mm <- sites_shared_single[1:18]
  inds_8 <- grep("^8mer\\|.*\\|.*_Kd", sites_shared_single, value=TRUE, perl=TRUE)
  inds_7m8 <- grep("^7mer-m8\\|.*\\|.*_Kd", sites_shared_single, value=TRUE, perl=TRUE)
  inds_7a1 <- grep("^7mer-A1\\|.*\\|.*_Kd", sites_shared_single, value=TRUE, perl=TRUE)
  inds_6 <- grep("^6mer\\|.*\\|.*_Kd", sites_shared_single, value=TRUE, perl=TRUE)
  inds_6m8 <- grep("^6mer-m8\\|.*\\|.*_Kd", sites_shared_single, value=TRUE, perl=TRUE)
  inds_6a1 <- grep("^6mer-A1\\|.*\\|.*_Kd", sites_shared_single, value=TRUE, perl=TRUE)



  inds_rest <- setdiff(sites_shared_single, c(inds_mm, inds_8, inds_7m8, inds_7a1,
                                       inds_6, inds_6m8, inds_6a1))
  segments(xmin, ymin, x1=xmax, y1=ymax, col="gray")
  Points(kds_single[inds_rest, 1], kds_single[inds_rest, 2],
         col=rgb(0, 0, 0, alpha=0.1))
  Points(kds_single[inds_8, 1], kds_single[inds_8, 2],
         col=ConvertRColortoRGB(kSiteColors["8mer"], alpha=alpha_sites))
  Points(kds_single[inds_7m8, 1], kds_single[inds_7m8, 2],
         col=ConvertRColortoRGB(kSiteColors["7mer-m8"], alpha=alpha_sites))
  Points(kds_single[inds_7a1, 1], kds_single[inds_7a1, 2],
         col=ConvertRColortoRGB(kSiteColors["7mer-A1"], alpha=alpha_sites))
  Points(kds_single[inds_6, 1], kds_single[inds_6, 2],
         col=ConvertRColortoRGB(kSiteColors["6mer"], alpha=alpha_sites))
  Points(kds_single[inds_6a1, 1], kds_single[inds_6a1, 2],
         col=ConvertRColortoRGB(kSiteColors["6mer-A1"], alpha=alpha_sites))
  Points(kds_single[inds_6m8, 1], kds_single[inds_6m8, 2],
         col=ConvertRColortoRGB(kSiteColors["6mer-m8"], alpha=alpha_sites))
  Points(kds_single[1:18, 1], kds_single[1:18, 2],
         col=ConvertRColortoRGB("forestgreen", alpha=alpha_sites))

  xy <- GetPlotFractionalCoords(0.1, 0.95, log="xy')")
  AddCorrelationToPlot(x=log(kds_single[, 1]), y=log(kds_single[, 2]),
                       xpos=xy[1], ypos=xy[2], rsquared=TRUE)

  AddLogAxis(1, "Relative Kd; replicate 1")
  AddLogAxis(2, "Relative Kd; replicate 2")
  # Finish the plot.
  if (class(pdf.plot) == "character") {
    dev.off()
  }

}

PlotScanningThreePrimeSite <- function(
  mirna="let-7a", seed_mmsite= "8mer-mmC7", height=5, width=5, xpos=20, ypos=20,
  pdf.plot=FALSE
) {
  kd_path <- sprintf("ThreePrimeTargetingPaper/temp_%s_kds.txt", mirna)
  kds <- read.table(kd_path, sep="\t", row.names=1, header=TRUE)
  kds_average <- CollapseSeedSitesTable(kds)


  rownames(kds) <- gsub("^6mer\\|28\\|", x=rownames(kds), perl=TRUE,
                          replacement="7mer-m8\\|28\\|")
  rownames(kds) <- gsub("^7mer-A1\\|27\\|", x=rownames(kds), perl=TRUE,
                          replacement="8mer\\|27\\|")
  # ROUGH CODE SECTION: Try plotting just the 5mer-sites, see if it checks out.
  sites_5mer <- GetSeedKdMatrix(kds, "5mer-m11.15")
  print(sites_5mer)
  SubfunctionCall(FigureSaveFile2)
  xmin <- min(as.integer(colnames(sites_5mer)))
  xmax <- 25
  ymin <- 0.001
  ymax <- 1
  BlankPlot(log="y")

  AddLinearAxis(1, tick.space=1, label.space=5,
                label="Position within random library")
  AddLogAxis(2, label="Relative Kd")
  mm_ind <- sprintf("%s_Kd", seed_mmsite)
  lines(colnames(sites_5mer), 10^sites_5mer[mm_ind, ],
        col="red")
  cols <- c(`6mer-m11.16`="green", `7mer-m11.17`="forestgreen",
            `8mer-m11.18`="blue", `9mer-m11.19`="purple")
  for (site_len in names(cols)) {
    sites_kds <- GetSeedKdMatrix(kds, site_len)
    lines(colnames(sites_kds), 10^sites_kds[mm_ind, ], col=cols[site_len])
  }
  print(kds[1:18, , drop=FALSE])
  mismatch_kd <- kds[sprintf("NA|NA|%s_Kd", seed_mmsite), 1]
  
  segments(xmin, 10^mismatch_kd, x1=xmax, y1=10^mismatch_kd, lty=2)
  kds_collapsed_seed <- CollapseSeedSitesTable(kds)
  s8mer_kd <- kds_collapsed_seed["8mer|NA|NA_Kd", 1]
  s6mer_kd <- kds_collapsed_seed["6mer|NA|NA_Kd", 1]
  segments(xmin, 10^s8mer_kd, x1=xmax, y1=10^s8mer_kd, lty=2,
           col=kSiteColors["8mer"])
  segments(xmin, 10^s6mer_kd, x1=xmax, y1=10^s6mer_kd, lty=2,
           col=kSiteColors["6mer"])
  # Finish the plot.
  if (class(pdf.plot) == "character") {
    dev.off()
  }

}



PlotOffSetAlignmentMatrix <- function(
  mirna, experiment, offset, n_constant=5, sitelist="programmed_collapsed",
  kdrel=TRUE, outto11mers=FALSE, collapsemm=FALSE,
  deldelG=FALSE, deldelGleft=FALSE, key=FALSE, xlabels=TRUE, mirna_label=FALSE,
  height=2.8, width=2.5, xpos=20, ypos=20, pdf.plot=FALSE
) {
  # Load the collapsed data (for the mismatch-only).
  message(sprintf("offset: %s", offset))
  mm8mer_sites <- GetAll8merMmSites(mirna)
  canon_sites <- c("8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-A1", "6mer-m8")
  # Get the Geometric average of the compensatory and supplemental pairing.
  if (sitelist == "programmed_suppcomp") {
    kds_thrp <- SubfunctionCall(EquilPars)
    if (collapsemm) {
      kds_comp <- GeoMean(kds_thrp["Comp_Kd", 2])
      kds_supp <- GeoMean(kds_thrp["Supp_Kd", 2])
    } else {
      kds_comp <- GeoMean(kds_thrp[sprintf("%s_Kd", mm8mer_sites), 2])
      kds_supp <- GeoMean(kds_thrp[sprintf("%s_Kd", canon_sites), 2])      
    }
  } else if (sitelist == "programmed_collapsed") {
    kds_collapsed <- SubfunctionCall(EquilPars)
    kds_supp <- kds_collapsed[sprintf("%s_Kd", canon_sites), 2]
    kds_comp <- kds_collapsed[sprintf("%s_Kd", mm8mer_sites), 2]
  }
  if (kdrel) {
    kd_base <- kds_comp
  } else {
    kd_base <- 1
  }
  print(kds_comp)
  print(kds_supp)
  len_mir <- nchar(kMirnaSeqs[mirna])
  # Define the possible 5-prime starting nucleotides and possible 3-prime
  # starting nucleotides, for overall matrix.
  nucs_5p <- 9:(len_mir - 5 + 1)
  nucs_3p <- 13:len_mir
  # Define the output matrix.
  out_matrix <- matrix(NA, nrow=length(nucs_3p), ncol=length(nucs_5p),
                       dimnames=list(nucs_3p, nucs_5p))
  if (outto11mers) {
    len_max <- 11
  } else {
    len_max <- 9
  }
  for (kmer_len in 5:len_max) {
    print(kmer_len)
    starts_i <- 9:(len_mir - kmer_len + 1)
    stops_i <- starts_i + kmer_len - 1
    sites <- sprintf("%smer-m%s.%s", kmer_len, starts_i, stops_i)
    print(sites)
    if (sitelist == "programmed_suppcomp") {
      sites_use <- sprintf("%s|%s|Comp_Kd", sites, starts_i + offset)
      print(sites_use)
      kds <- kds_thrp[sites_use, 2]/kd_base
      for (i in 1:length(starts_i)) {
        out_matrix[as.character(stops_i[i]), as.character(starts_i[i])] <- kds[i]
      }
      print(out_matrix)
    } else {
      vals <- t(sapply(sites, function(site) {
        kds <- SubfunctionCall(GetPositionalProgKds)
        num_pos <- 25 - kmer_len + 1
        mir_pos <- c("Seed", as.character(9:(9 + num_pos - 1)))
        kds_out <- matrix(0, nrow=18, ncol=1 + 25 - kmer_len + 1,
                          dimnames=list(mm8mer_sites, mir_pos))
        for (mm_site in mm8mer_sites) {
          for (pos in mir_pos) {
            if (pos == "Seed") {
              site_name <- sprintf("%s_Kd", mm_site)
              val <- kds_collapsed[site_name, 2]
            } else {
              site_name <- sprintf("%s|%s|%s_Kd", site, pos, mm_site)
              val <- kds[site_name, 2]
            }
            if (length(val) != 0) {
              kds_out[mm_site, as.character(pos)] <- val          
            }
          }
        }
        if (kdrel) {
          kds_out <- (t(t(kds_out/kds_out[, 1])))
        }
        kds_out <- kds_out[, 2:ncol(kds_out)]
        return(exp(colMeans(log(kds_out), na.rm=TRUE)))
      }))
      for (ind in 1:length(sites)) {
        # message("in site loop.")
        site <- sites[ind]
        start_stop <- unlist(strsplit(site, split="-m"))[2]
        start_stop <- unlist(strsplit(start_stop, split="\\."))
        start_out <- start_stop[1]
        stop_out <- start_stop[2]
        val <- vals[site, ][ind + offset]
        if (length(val) != 0) {
          out_matrix[stop_out, start_out] <- val      
        }
      }
    }

  }
  # print(out_matrix)
  mar1 <- 2.2
  mar2 <- 2.5
  mar3 <- 1.8
  mar4 <- 0.5
  # This is the conversion factor in order to have the correct height of the
  # plot for the specified width, that makes each box of equal height and width.

  R_mat <- -log10(out_matrix)
  R_mat <- t(R_mat)

  if (deldelGleft) {
    R_mat1 <- R_mat[-1, ]
    R_mat2 <- R_mat[-nrow(R_mat), ]
    R_mat <- t((R_mat2 - R_mat1)[, -1])

    R_mat_del <<- R_mat
  } else if (deldelG) {
    R_mat1 <- R_mat[, -ncol(R_mat)]
    R_mat2 <- R_mat[, -1]
    R_mat <- (R_mat2 - R_mat1)[-nrow(R_mat), ]
    R_mat_del <<- R_mat
  }
  mir_length <- nchar(kMirnaSeqs[mirna])
  SubfunctionCall(FigureSaveFile2)
  par(mar=c(mar1, mar2, mar3, mar4))
  xmin <- 0
  xmax <- nrow(R_mat)
  ymin <- 0
  ymax <- ncol(R_mat) 
  BlankPlot()
  ymin <- 0
  test_rows <- rowMeans(R_mat, na.rm=TRUE)
  test_rows <<- test_rows
  test_row_check <- sapply(1:(length(test_rows) - 1), function(ind) {
    mean(test_rows[ind:(ind+1)], na.rm=TRUE)  
  })
  xlefts <- rep(seq(0, ncol(R_mat) - 1), nrow(R_mat))
  xright <- xlefts + 1
  ybottom <- rep(seq(ymax - 1, ymin), each=ncol(R_mat))
  ytop <- ybottom + 1

  # Make the x-axis
  AddLinearAxis(2, 1, 1, label="5'-paired nt",
                alt_lab=rev(rownames(R_mat)),
                alt_lab_pos=ymin:(ymax - 1) + 0.5,
                alt_tick_pos=TRUE,
                line=0.8)
  AddLinearAxis(1, 1, 1, label="3'-paired nt",
                label_pos_ticks=TRUE,
                alt_lab=colnames(R_mat),
                alt_lab_y_dist=0.01,
                alt_lab_pos=xmin:(xmax - 1) + 0.5,
                alt_tick_pos=TRUE)
  # Add the label for the seed if not plotting relative to seed kds.
  if (!kdrel) {
    text(-1.5, ymin - 0.5, labels="Seed", srt=90, adj=c(1, 0.5), xpd=NA)
  }
  if (kdrel) {
    if (deldelGleft || deldelG) {
      col.inds <- round((R_mat)/log10(20)*100)
      x_lab_pos <- -0.75
    } else {
      col.inds <- round((R_mat)/log10(100)*100)
      x_lab_pos <- -0.75
    }
  } else {
    col.inds <- round((R_mat - log10(1/3))/3*100)
    x_lab_pos <- -2.5
  }
  # Load the seaborn cubeHelix color palette, and conver the indeces from the
  # and make them bounded between zero and 1
  if (kdrel) {
    if (deldelGleft) {
      start_col <- 0.4
      r_col <- -0.4
    } else if (deldelG) {
      start_col <- 0.4
      r_col <- -0.4
    } else {
      start_col <- 0
      r_col <- 0.4
    }
  } else {
    start_col <- 0.5
    r_col <- -0.75
  }
  # Make distribution of 100 colors and a 101th gray color for NA indeces.
  color.dist <- c(rev(cubeHelix(100, start=start_col, r=r_col, hue=0.8)), "gray")
  # Make the color index scale.
  col.inds <- sapply(t(col.inds), function(col.ind) {
    min(max(1, col.ind), 100)
  })
  col.inds[which(is.na(col.inds))] <- 101
  # Make the color rectangle.
  rect(xlefts, ybottom, xright, ytop, col=color.dist[col.inds], lwd=0.4, 
       xpd=NA, border="white")
  # Make the key.
  if (key) {
    y_div <- (ymax-0.5)/100
    kpl <- xmax + 1 # Left-hand position of the key
    kw <- 1.5                # Width of the key
    rect(xleft=kpl, ybottom=seq(100)*y_div - y_div + 0.5, xright=kpl + kw,
         ytop=seq(100)*y_div - y_div + 0.5 + y_div,
         col=color.dist[1:100], xpd=NA, border=NA)
    # Generate the axis for the legend and the label
    bottom <- ymin + y_div/2
    top <- ymax - y_div/2
    if (kdrel) {
      labels <- c(0.3, 1, 3, 10, 30, 100)
    } else {
      # if (site %in% c("8mer", "7mer-m8", "7mer-A1", "6mer")) {
      #   labels <- c(3, 1, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001, 0.0003)
      # } else {
        labels <- c(3, 1, 0.3, 0.1, 0.03, 0.01, 0.003)
      # }
    }
    pos_labels <- log(labels)
    centered_labels <- pos_labels - pos_labels[1]
    norm_labels <- centered_labels/(centered_labels[length(centered_labels)])
    # labels_span <- max(log(labels)) - min(log(labels))
    height_span <- top - bottom
    pos_labels <- norm_labels*height_span
    axis(4, at=pos_labels, labels=labels, lwd=par()$lwd, pos=kpl + kw, xpd=NA)
    if (kdrel) {
      text(x=kpl + kw + 3, y=ymax/2,
           labels=bquote(italic(K)[D]*.(" fold change")),
           srt=270, xpd=NA)    
    } else {
      text(x=kpl + kw + 3.5, y=ymax/2,
           labels=bquote(.("Relative")~italic(K)[D]),
           srt=270, xpd=NA)    
    }
  }
  # # Add the label indicating how many nucleotides of pairing.
  mtext(text=sprintf("%s nt offset", offset), side=3, at=xmax, adj=1,
        cex=par("cex"))
  # If the `mirna_label` conditional is true, add the label saying which miRNA
  # is being looked at.
  if (mirna_label) {
    if (mirna == "let-7a-21nt") {
      mirna_txt <- "let-7a"
    } else if (mirna == "let-7a_minus1") {
      mirna_txt <- "let-7a(-1)"
    } else if (mirna == "let-7a_plus1") {
      mirna_txt <- "let-7a(+1)"
    } else if (mirna == "let-7a_miR-155") {
      mirna_txt <- "let-7a-miR-155"
    } else if (mirna == "miR-155_let-7a") {
      mirna_txt <- "miR-155-let-7a"
    } else {
      mirna_txt <- mirna
    }
    mtext(text=mirna_txt, side=3, line=0.7, at=xmin, adj=0, cex=par("cex"))
  }

  # Finish the plot.
  if (class(pdf.plot) == "character") {
    dev.off()
  }
  return(R_mat)
}

PlotOffsetCoefficients <- function(
  n_constant=5, offsetmin=0, offsetmax=17, sitelist="programmed_suppcomp",
  let7_series=FALSE, chimera_series=FALSE, chimera_let7_series=FALSE,
  chimera_miR155_series=FALSE, outto11mers=TRUE, collapsemm=FALSE,
  supp_base=FALSE, xpos=20, ypos=20, height=4, width=4, pdf.plot=NULL
) {
  if (let7_series) {
    mirna_tables <- matrix(c("let-7a_minus1", "equil_c_nb",
                             "let-7a-21nt", "equil_c2_nb",
                             "let-7a_plus1", "equil_c_nb"), byrow=TRUE, nrow=3)

  } else if (chimera_series) {
    mirna_tables <- matrix(c("let-7a-21nt", "equil_c2_nb",
                             "miR-155", "equil_sc_nb",
                             "miR-155_let-7a", "equil_c_nb",
                             "let-7a_miR-155", "equil_c_nb"), byrow=TRUE,
                           nrow=4)
  } else if (chimera_let7_series) {
    mirna_tables <- matrix(c("let-7a-21nt", "equil_c2_nb",
                             "miR-155_let-7a", "equil_c_nb"), byrow=TRUE,
                           nrow=2)
  } else if (chimera_miR155_series) {
    mirna_tables <- matrix(c("miR-155", "equil_sc_nb",
                             "let-7a_miR-155", "equil_c_nb"), byrow=TRUE,
                           nrow=2)
  } else {
    mirna_tables <- matrix(c("let-7a-21nt", "equil_c2_nb",
                             "miR-1", "equil_c_nb",
                             "miR-155", "equil_sc_nb"), byrow=TRUE, nrow=3)    
  }
  offsets_df <- apply(mirna_tables, 1, function(row) {
    model <- FitOffsetAndPairingModel(
      mirna=row[1], exp=row[2], n_constant=n_constant, offsetmin=offsetmin,
      offsetmax=offsetmax, sitelist=sitelist, outto11mers=outto11mers,
      collapsemm=collapsemm, supp_base=supp_base
    )
    offsets <- c(0, model[["offsets"]])
    names(offsets)[1] <- sprintf("offset%s", offsetmin)
    offsets
  })
  print(offsets_df)
  # colnames(offsets_df) <- mirna_tables[, 1]
  SubfunctionCall(FigureSaveFile2)
  xmin <- offsetmin
  xmax <- offsetmax
  ymin <- 0.3
  ymax <- 10
  BlankPlot(, log="y")

  if (let7_series) {
    colnames(offsets_df) <- c("let-7a(-1)", "let-7a", "let-7a(+1)")
    # xmax <- 21 - len_k + 1.2
    cols <- c("dodgerblue", "#F7931D", "darkorchid3")
  } else if (chimera_series) {
    colnames(offsets_df) <- c("let-7a", "miR-155", "miR-155-let-7a", "let-7a-miR-155")
    # xmax <- 23 - len_k + 1.2
    cols <- c("#F7931D", "darkorchid3", "dodgerblue", "deeppink")
  } else if (chimera_let7_series) {
    colnames(offsets_df) <- c("let-7a", "miR-155-let-7a")
    # xmax <- 23 - len_k + 1.2
    cols <- c("dodgerblue", "deeppink")
  } else if (chimera_miR155_series) {
    colnames(offsets_df) <- c("miR-155", "let-7a-miR-155")
    # xmax <- 23 - len_k + 1.2
    cols <- c("darkorchid3", "deeppink")
  } else {
    colnames(offsets_df) <- c("let-7a", "miR-1", "miR-155")  
    # xmax <- 23 - len_k + 1.2
    cols <- c("#F7931D", "dodgerblue", "darkorchid3")
  }

  names(cols) <- colnames(offsets_df)  

  sapply(colnames(offsets_df), function(mirna) {
    lines(xmin:xmax, 10^offsets_df[, mirna], col=cols[mirna])  
    points(xmin:xmax, 10^offsets_df[, mirna], col=cols[mirna], pch=20)  
  })
  # Add the axes.
  AddLinearAxis(1, 1, 1, label="Offset (nt)")
  AddLogAxis(2, 1, 1, label="Kd fold change")

  xy <- GetPlotFractionalCoords(0.5, 0.5)
  Legend(xy, legend=colnames(offsets_df), col=cols)

  # Finish the plot.
  if (class(pdf.plot) == "character") {
    dev.off()
  }

}

PlotLet7aProgrammedLibraryAgainstRandom <- function(
  average_sites=FALSE, height=5, width=5, xpos=20, ypos=20, pdf.plot=FALSE
) {
  kds_programmed <- read.table("ThreePrimeTargetingPaper/temp_let-7a_kds.txt",
                               sep="\t", row.names=1, header=TRUE)
  # kds_2 <- read.table("ThreePrimeTargetingPaper/temp_let-7a-21nt_kds.txt", sep="\t",
  #                        row.names=1, header=TRUE)
  # if (average_sites) {
  kds_random <- EquilPars("let-7a")
  kds_programmed_collapsed <- CollapseSeedSitesTable(kds_programmed)
    # kds_2 <- CollapseSeedSitesTable(kds_2)
  alpha_sites <- 1
  print(kds_programmed_collapsed[1:20, , drop=FALSE])
  print(kds_random)


  row_names_random <- paste0(kSeedSites, "_Kd")
  row_names_programmed <- paste0(kSeedSites, "|NA|NA_Kd")

  # print(kds_programmed_collapsed[row_names_programmed, ])
  # print(kds_random[row_names_random, ])


  # break    
  # # } else {
  # #   alpha_sites <- 0.5
  # # }
  # #TODO fix this, such that it doesn't need to be hardcoded into the data.
  # # Realistically, I should extend to 5 nucleotides outside.
  # # That makes the most sense.
  # rownames(kds) <- gsub("^6mer\\|28\\|", x=rownames(kds), perl=TRUE,
  #                         replacement="7mer-m8\\|28\\|")
  # # rownames(kds_2) <- gsub("^6mer\\|28\\|", x=rownames(kds_2), perl=TRUE,
  # #                         replacement="7mer-m8\\|28\\|")

  # rownames(kds_1) <- gsub("^7mer-A1\\|27\\|", x=rownames(kds_1), perl=TRUE,
  #                         replacement="8mer\\|27\\|")
  # # rownames(kds_2) <- gsub("^7mer-A1\\|27\\|", x=rownames(kds_2), perl=TRUE,
  # #                         replacement="8mer\\|27\\|")

  # sites_shared <- intersect(rownames(kds_1), rownames(kds_2))
  # sites_shared_single <- grep("&", sites_shared, invert=TRUE, value=TRUE)

  # kds_shared_1 <- kds_1[sites_shared, 1]
  # kds_shared_2 <- kds_2[sites_shared, 1]

  # kds_single_1 <- kds_1[sites_shared_single, 1]
  # kds_single_2 <- kds_2[sites_shared_single, 1]

  # kds_shared <- t(rbind(kds_shared_1, kds_shared_2))
  # kds_single <- 10^t(rbind(kds_single_1, kds_single_2))
  # rownames(kds_shared) <- sites_shared
  # rownames(kds_single) <- sites_shared_single
  SubfunctionCall(FigureSaveFile2)
  xmin=0.0001
  xmax=1
  ymin=0.0001
  ymax=1
  BlankPlot(log="xy")

  AddLogAxis(1, label="Relative Kd; Programmed Library")
  AddLogAxis(2, label="Relative Kd; Random Library")
  points(10^kds_programmed_collapsed[row_names_programmed, ],
       kds_random[row_names_random, 2], col=kSiteColors[kSeedSites])

  xy <- GetPlotFractionalCoords(0.9, 0.1, log='xy')
  AddCorrelationToPlot(xy[1], xy[2],
                       x=kds_programmed_collapsed[row_names_programmed, 1],
                       y=log(kds_random[row_names_random, 2]), rsquared=TRUE,
                       adj=c(1, 0))
  # break
  # inds_mm <- sites_shared_single[1:18]
  # inds_8 <- grep("^8mer\\|.*\\|.*_Kd", sites_shared_single, value=TRUE, perl=TRUE)
  # inds_7m8 <- grep("^7mer-m8\\|.*\\|.*_Kd", sites_shared_single, value=TRUE, perl=TRUE)
  # inds_7a1 <- grep("^7mer-A1\\|.*\\|.*_Kd", sites_shared_single, value=TRUE, perl=TRUE)
  # inds_6 <- grep("^6mer\\|.*\\|.*_Kd", sites_shared_single, value=TRUE, perl=TRUE)
  # inds_6m8 <- grep("^6mer-m8\\|.*\\|.*_Kd", sites_shared_single, value=TRUE, perl=TRUE)
  # inds_6a1 <- grep("^6mer-A1\\|.*\\|.*_Kd", sites_shared_single, value=TRUE, perl=TRUE)



  # inds_rest <- setdiff(sites_shared_single, c(inds_mm, inds_8, inds_7m8, inds_7a1,
  #                                      inds_6, inds_6m8, inds_6a1))
  # segments(xmin, ymin, x1=xmax, y1=ymax, col="gray")
  # Points(kds_single[inds_rest, 1], kds_single[inds_rest, 2],
  #        col=rgb(0, 0, 0, alpha=0.1))
  # Points(kds_single[inds_8, 1], kds_single[inds_8, 2],
  #        col=ConvertRColortoRGB(kSiteColors["8mer"], alpha=alpha_sites))
  # Points(kds_single[inds_7m8, 1], kds_single[inds_7m8, 2],
  #        col=ConvertRColortoRGB(kSiteColors["7mer-m8"], alpha=alpha_sites))
  # Points(kds_single[inds_7a1, 1], kds_single[inds_7a1, 2],
  #        col=ConvertRColortoRGB(kSiteColors["7mer-A1"], alpha=alpha_sites))
  # Points(kds_single[inds_6, 1], kds_single[inds_6, 2],
  #        col=ConvertRColortoRGB(kSiteColors["6mer"], alpha=alpha_sites))
  # Points(kds_single[inds_6a1, 1], kds_single[inds_6a1, 2],
  #        col=ConvertRColortoRGB(kSiteColors["6mer-A1"], alpha=alpha_sites))
  # Points(kds_single[inds_6m8, 1], kds_single[inds_6m8, 2],
  #        col=ConvertRColortoRGB(kSiteColors["6mer-m8"], alpha=alpha_sites))
  # Points(kds_single[1:18, 1], kds_single[1:18, 2],
  #        col=ConvertRColortoRGB("forestgreen", alpha=alpha_sites))

  # xy <- GetPlotFractionalCoords(0.1, 0.95, log="xy')")
  # AddCorrelationToPlot(x=log(kds_single[, 1]), y=log(kds_single[, 2]),
  #                      xpos=xy[1], ypos=xy[2], rsquared=TRUE)

  # AddLogAxis(1, "Relative Kd; replicate 1")
  # AddLogAxis(2, "Relative Kd; replicate 2")
  # Finish the plot.
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

PlotScanningThreePrimeSiteOld <- function(
  mirna="let-7a-21nt", experiment="equil_c2_nb", seed_mmsite= "8mer-mmC7",
  n_constant=5, sitelist="programmed", equilibrium_nb=FALSE, height=5, width=5,
  xpos=20, ypos=20, pdf.plot=FALSE
) {
  thrp_site_list <- c("5mer-m11.15", "6mer-m11.16", "7mer-m11.17",
                         "8mer-m11.18", "9mer-m11.19")
  kds <- SubfunctionCall(EquilPars, sitelist="programmed_collapsed")
  SubfunctionCall(FigureSaveFile2)
  xmin <- 8
  xmax <- 25
  ymin <- 0.001
  ymax <- 1
  BlankPlot(log="y")

  AddLinearAxis(1, tick.space=1, label.space=5,
                label="Position with respect to miRNA")
  AddLogAxis(2, label="Relative Kd")
  mm_ind <- sprintf("%s_Kd", seed_mmsite)
  # lines(colnames(sites_5mer), 10^sites_5mer[mm_ind, ],
        # col="red")
  cols <- c(`5mer-m11.15`="red", `6mer-m11.16`="green",
            `7mer-m11.17`="forestgreen", `8mer-m11.18`="blue",
            `9mer-m11.19`="purple")
  offsets <- c(-0.1, -0.05, 0, 0.05, 0.1)
  names(offsets) <- names(cols)
  sapply(rev(thrp_site_list), function(site) {
    site_kds <- SubfunctionCall(GetPositionalProgKds)
    inds <- grep("8mer-mmC7", rownames(site_kds))
    df_use <- site_kds[inds, ]
    pos <- gsub("^(.*)\\|(.*)\\|(.*_Kd$)", replacement="\\2", rownames(df_use))
    rownames(df_use) <- pos
    x <- 9:min(max(as.integer(pos)), xmax)
    segments(x0=x + offsets[site], y0=df_use[as.character(x), 3],
             y1=df_use[as.character(x), 5],
             col=ConvertRColortoRGB(cols[site], alpha=0.5), lwd=1.5, xpd=NA)
    lines(x, df_use[as.character(x), 2], col=cols[site], lwd=1, xpd=NA)
    points(x, df_use[as.character(x), 2], col=cols[site], pch=19, xpd=NA)
  })
  mismatch_kd <- kds[sprintf("%s_Kd", seed_mmsite), 1]
  segments(xmin, mismatch_kd, x1=xmax, y1=mismatch_kd, lty=2)
  # kds_collapsed_seed <- CollapseSeedSitesTable(kds)
  s8mer_kd <- GeoMean(kds[grep("^8mer\\|8mer-mm[ACGT][2-7]_Kd$",
                       rownames(kds), perl=TRUE), 2])
  s6mer_kd <- GeoMean(kds[grep("^6mer\\|8mer-mm[ACGT][2-7]_Kd$",
                       rownames(kds), perl=TRUE), 2])
  mm_kd <- GeoMean(kds[grep("^8mer-mm[ACGT][2-7]_Kd$",
                            rownames(kds), perl=TRUE), 2])
  print(mm_kd)
  print(s8mer_kd)
  print(s6mer_kd)
  s8mer_kd_norm <- s8mer_kd/mm_kd
  s6mer_kd_norm <- s6mer_kd/mm_kd
  if (mirna %in% c("let-7a-21nt", "let-7a_plus1", "let-7a_minus1",
                   "let-7a_miR-155")) {
    if (equilibrium_nb) {
      mirna_rand <- "let-7a-21nt"
      experiment_rand <- "equilibrium_nb"
    } else {
      mirna_rand <- "let-7a"
      experiment_rand <- "equilibrium"
    }
  } else if (mirna == "miR-155_let-7a") {
    mirna_rand <- "miR-155"
    experiment_rand <- "equilibrium"
  } else {
    mirna_rand <- mirna
    experiment_rand <- "equilibrium"
  }
  if (mirna == "miR-1") {
    buffer_rand <- TRUE
    combined_rand <- FALSE
  } else {
    buffer_rand <- FALSE
    combined_rand <- TRUE
  }


  # kds_rand <- SubfunctionCall(EquilPars, mirna=mirna_rand,
  #                             experiment=experiment_rand, n_constant=5,
  #                             sitelist="programmed", buffer=buffer_rand,
  #                             combined=combined_rand)
  # kds_rand <<- kds_rand
  # s8mer_kd_rand <- kds_rand["8mer_Kd", 2]
  # s6mer_kd_rand <- kds_rand["6mer_Kd", 2]

  segments(xmin, s8mer_kd, x1=xmax, y1=s8mer_kd, lty=2,
           col=kSiteColors["8mer"])
  segments(xmin, s6mer_kd, x1=xmax, y1=s6mer_kd, lty=2,
           col=kSiteColors["6mer"])


  segments(xmin, s8mer_kd_norm, x1=xmax, y1=s8mer_kd_norm, lty=3,
           col=kSiteColors["8mer"])
  segments(xmin, s6mer_kd_norm, x1=xmax, y1=s6mer_kd_norm, lty=3,
           col=kSiteColors["6mer"])
  xy <- GetPlotFractionalCoords(0.1, 1, log="y")
  Legend(xy, legend=names(cols), col=cols, ncol=2)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

PlotSimNucleotideContributionMatrix <- function(
  mirna, experiment, offset_sim, n_constant=0, offsetmin=-4, offsetmax=17,
  sitelist="programmed_suppcomp", outto11mers=TRUE, downto2mers=FALSE,
  downto3mers=FALSE, downto4mers=FALSE, end3prand=FALSE, collapsemm=FALSE,
  deldelG=FALSE, deldelGleft=FALSE, residual=TRUE, globalmodel=FALSE,
  key=FALSE, xlabels=TRUE, mirna_label=FALSE, extralabel=FALSE, height=2.8,
  width=2.5, xpos=20, ypos=20, pdf.plot=FALSE
) {
  # Load the data matrix.
  if (globalmodel) {
    model <- SubfunctionCall(FitOffsetAndPairingModel)
    model <<- model    
  }
  R_mat <- t(model[["pairing"]])
  R_mat_model <<- R_mat
  if (offset_sim != offsetmin) {
    offset_name <- sprintf("offset%s", offset_sim)
    R_mat <- R_mat + model[["offsets"]][offset_name]
  }
  R_mat_data <- t(SubfunctionCall(MakeNucleotideContributionMatrix, offset=offset_sim))
  R_mat_data <<- R_mat_data
  print(dim(R_mat_data))
  print(dim(R_mat))
  # R_mat_data <- R_mat_data*0 + 1
  R_mat <- R_mat*(R_mat_data*0 + 1)
  if (residual) {
    R_mat <- R_mat_data - R_mat
    print(min(R_mat, na.rm=TRUE))
    print(max(R_mat, na.rm=TRUE))
  }
  R_mat <<- R_mat

  mar1 <- 2.2
  mar2 <- 2.5
  mar3 <- 1.8
  mar4 <- 0.5
  # This is the conversion factor in order to have the correct height of the
  # plot for the specified width, that makes each box of equal height and width.
  if (deldelGleft) {
    R_mat1 <- R_mat[-1, ]
    R_mat2 <- R_mat[-nrow(R_mat), ]
    R_mat <- t((R_mat2 - R_mat1)[, -1])

  } else if (deldelG) {
    R_mat1 <- R_mat[, -ncol(R_mat)]
    R_mat2 <- R_mat[, -1]
    R_mat <- (R_mat2 - R_mat1)[-nrow(R_mat), ]
    R_mat_del <<- R_mat
  }
  mir_length <- nchar(kMirnaSeqs[mirna])
  SubfunctionCall(FigureSaveFile2)
  par(mar=c(mar1, mar2, mar3, mar4))
  xmin <- 0
  xmax <- nrow(R_mat)
  ymin <- 0
  ymax <- ncol(R_mat) 
  BlankPlot()
  ymin <- 0
  test_rows <- rowMeans(R_mat, na.rm=TRUE)
  test_rows <<- test_rows
  test_row_check <- sapply(1:(length(test_rows) - 1), function(ind) {
    mean(test_rows[ind:(ind+1)], na.rm=TRUE)  
  })
  xlefts <- rep(seq(0, ncol(R_mat) - 1), nrow(R_mat))
  xright <- xlefts + 1
  ybottom <- rep(seq(ymax - 1, ymin), each=ncol(R_mat))
  ytop <- ybottom + 1

  # Make the x-axis
  AddLinearAxis(2, 1, 1, label="5'-paired nt",
                alt_lab=rev(rownames(R_mat)),
                alt_lab_pos=ymin:(ymax - 1) + 0.5,
                alt_tick_pos=TRUE,
                line=0.8)
  AddLinearAxis(1, 1, 1, label="3'-paired nt",
                label_pos_ticks=TRUE,
                alt_lab=colnames(R_mat),
                alt_lab_y_dist=0.01,
                alt_lab_pos=xmin:(xmax - 1) + 0.5,
                alt_tick_pos=TRUE)
  # Add the label for the seed if not plotting relative to seed kds.
  # if (kdrel) {
    if (deldelGleft || deldelG) {
      col.inds <- round((R_mat)/log10(20)*100)
      x_lab_pos <- -0.75
    } else if (residual) {
      col.inds <- round((R_mat)/log10(20)*50 + 50)
      print(min(col.inds, na.rm=TRUE))
      print(max(col.inds, na.rm=TRUE))
    } else {
      col.inds <- round((R_mat)/log10(100)*100)
      x_lab_pos <- -0.75
    }
  # } else {
  #   col.inds <- round((R_mat - log10(1/3))/3*100)
  #   x_lab_pos <- -2.5
  # }
  # Load the seaborn cubeHelix color palette, and conver the indeces from the
  # and make them bounded between zero and 1
  # if (kdrel) {
    if (deldelGleft) {
      start_col <- 0.4
      r_col <- -0.4
    } else if (deldelG) {
      start_col <- 0.4
      r_col <- -0.4
    } else {
      start_col <- 0
      r_col <- 0.4
    }
  # } else {
  #   start_col <- 0.5
  #   r_col <- -0.75
  # }
  # Make distribution of 100 colors and a 101th gray color for NA indeces.
  if (residual) {
    brbg <- RColorBrewer::brewer.pal(11, "BrBG")
    color.dist <- c(colorRampPalette(c(brbg[2], brbg[6]))(50), 
                    colorRampPalette(c(brbg[6], brbg[10]))(51)[-1], "gray")
    # color.dist <- c(rev(RColorBrewer::brewer.pal(100, "BrBG")), "gray")
  } else {
    color.dist <- c(rev(cubeHelix(100, start=start_col, r=r_col, hue=0.8)), "gray")
    # Make the color index scale.
    color.dist[100] <- "red"    
  }
  col.inds <- sapply(t(col.inds), function(col.ind) {
    min(max(1, col.ind), 100)
  })
  col.inds[which(is.na(col.inds))] <- 101

  # Make the color rectangle.
  rect(xlefts, ybottom, xright, ytop, col=color.dist[col.inds], lwd=0.4, 
       xpd=NA, border="white")
  # Make the key.
  if (key) {
    y_div <- (ymax-0.5)/100
    kpl <- xmax + 1 # Left-hand position of the key
    kw <- 1.5                # Width of the key
    rect(xleft=kpl, ybottom=seq(100)*y_div - y_div + 0.5, xright=kpl + kw,
         ytop=seq(100)*y_div - y_div + 0.5 + y_div,
         col=color.dist[1:100], xpd=NA, border=NA)
    # Generate the axis for the legend and the label
    bottom <- ymin + y_div/2
    top <- ymax - y_div/2
    # if (kdrel) {
    labels <- c(0.3, 1, 3, 10, 30, 100)
    # } else {
      # if (site %in% c("8mer", "7mer-m8", "7mer-A1", "6mer")) {
      #   labels <- c(3, 1, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001, 0.0003)
      # } else {
        # labels <- c(3, 1, 0.3, 0.1, 0.03, 0.01, 0.003)
      # }
    # }
    pos_labels <- log(labels)
    centered_labels <- pos_labels - pos_labels[1]
    norm_labels <- centered_labels/(centered_labels[length(centered_labels)])
    # labels_span <- max(log(labels)) - min(log(labels))
    height_span <- top - bottom
    pos_labels <- norm_labels*height_span
    axis(4, at=pos_labels, labels=labels, lwd=par()$lwd, pos=kpl + kw, xpd=NA)
    # if (kdrel) {
      text(x=kpl + kw + 3, y=ymax/2,
           labels=bquote(italic(K)[D]*.(" fold change")),
           srt=270, xpd=NA)    
    # } else {
    #   text(x=kpl + kw + 3.5, y=ymax/2,
    #        labels=bquote(.("Relative")~italic(K)[D]),
    #        srt=270, xpd=NA)    
    # }
  }
  # # Add the label indicating how many nucleotides of pairing.
  mtext(text=sprintf("%s nt offset", offset_sim), side=3, at=xmax, adj=1,
        cex=par("cex"))
  # If the `mirna_label` conditional is true, add the label saying which miRNA
  # is being looked at.
  if (mirna_label) {
    if (mirna == "let-7a-21nt") {
      mirna_txt <- "let-7a"
    } else if (mirna == "let-7a_minus1") {
      mirna_txt <- "let-7a(-1)"
    } else if (mirna == "let-7a_plus1") {
      mirna_txt <- "let-7a(+1)"
    } else if (mirna == "let-7a_miR-155") {
      mirna_txt <- "let-7a-miR-155"
    } else if (mirna == "miR-155_let-7a") {
      mirna_txt <- "miR-155-let-7a"
    } else {
      mirna_txt <- mirna
    }
    mtext(text=mirna_txt, side=3, line=0.7, at=xmin, adj=0, cex=par("cex"))
  }
  if (extralabel) {
    str.end3prand <- "r"
    if (collapsemm) {
      str.collapsemm <- "sum"
    } else {
      str.collapsemm <- "geo"
    }
    mtext(text=sprintf("%s_%s", str.end3prand, str.collapsemm), side=3,
          line=0.7, at=xmax, adj=1, cex=par("cex"))    
  }
  # Finish the plot.
  if (class(pdf.plot) == "character") {
    dev.off()
  }
  return(R_mat)
}

PlotSimNucleotideContributionMatrix <- function(
  mirna, experiment, offset_sim, n_constant=0, offsetmin=-4, offsetmax=17,
  sitelist="programmed_suppcomp", collapsemm=FALSE,
  deldelG=FALSE, deldelGleft=FALSE, residual=FALSE, globalmodel=FALSE,
  key=FALSE, xlabels=TRUE, mirna_label=FALSE, extralabel=FALSE, height=2.8,
  width=2.5, xpos=20, ypos=20, pdf.plot=FALSE
) {
  # Load the data matrix.
  if (globalmodel) {
    model <- SubfunctionCall(FitOffsetAndPairingNonlinModel, mm=FALSE)
    model <<- model    
  }
  R_mat <- t(model[["pairing"]])
  offset_name <- as.character(offset_sim)
  R_mat <- R_mat*(model[["offsets"]][offset_name]) + model[["base"]]
  R_mat_model <<- R_mat

  R_mat_data <- t(SubfunctionCall(MakeNucleotideContributionMatrix, offset=offset_sim))
  print("past R_mat_data")
  R_mat_data <<- R_mat_data
  print(R_mat)
  R_mat_data <- R_mat_data*0 + 1
  R_mat <- R_mat*(R_mat_data*0 + 1)
  if (residual) {
    R_mat <- R_mat_data - R_mat
    print(min(R_mat, na.rm=TRUE))
    print(max(R_mat, na.rm=TRUE))
  }
  R_mat <<- R_mat

  mar1 <- 2.2
  mar2 <- 0.5
  mar3 <- 1.8
  mar4 <- 2.5
  # This is the conversion factor in order to have the correct height of the
  # plot for the specified width, that makes each box of equal height and width.
  if (deldelGleft) {
    R_mat1 <- R_mat[-1, ]
    R_mat2 <- R_mat[-nrow(R_mat), ]
    R_mat <- t((R_mat2 - R_mat1)[, -1])

  } else if (deldelG) {
    R_mat1 <- R_mat[, -ncol(R_mat)]
    R_mat2 <- R_mat[, -1]
    R_mat <- (R_mat2 - R_mat1)[-nrow(R_mat), ]
    R_mat_del <<- R_mat
  }
  print("here")
  mir_length <- nchar(kMirnaSeqs[mirna])
  SubfunctionCall(FigureSaveFile2)
  par(mar=c(mar1, mar2, mar3, mar4))
  xmin <- 0
  xmax <- nrow(R_mat)
  ymin <- 0
  ymax <- ncol(R_mat)
  print(R_mat)
  BlankPlot()
  ymin <- 0
  test_rows <- rowMeans(R_mat, na.rm=TRUE)
  test_rows <<- test_rows
  test_row_check <- sapply(1:(length(test_rows) - 1), function(ind) {
    mean(test_rows[ind:(ind+1)], na.rm=TRUE)  
  })
  xlefts <- rep(seq(0, ncol(R_mat) - 1), nrow(R_mat))
  xright <- xlefts + 1
  ybottom <- rep(rev(seq(ymax - 1, ymin)), each=ncol(R_mat))
  ytop <- ybottom + 1

  pos_5p <- as.integer(rep(rownames(R_mat), each=ncol(R_mat)))
  pos_3p <- as.integer(rep(colnames(R_mat), nrow(R_mat)))

  impossible_cols <- which(pos_5p > pos_3p)
  # Make the x-axis
  AddLinearAxis(4, 1, 1, label="5'-paired nt",
                alt_lab=rev(rownames(R_mat)),
                alt_lab_pos=ymin:(ymax - 1) + 0.5,
                alt_tick_pos=TRUE,
                line=0.8)
  AddLinearAxis(1, 1, 1, label="3'-paired nt",
                label_pos_ticks=TRUE,
                alt_lab=colnames(R_mat),
                alt_lab_y_dist=0.01,
                alt_lab_pos=xmin:(xmax - 1) + 0.5,
                alt_tick_pos=TRUE)
  # Add the label for the seed if not plotting relative to seed kds.
  # if (kdrel) {
    if (deldelGleft || deldelG) {
      col.inds <- round((R_mat)/log10(20)*100)
      x_lab_pos <- -0.75
    } else if (residual) {
      col.inds <- round((R_mat)/log10(20)*50 + 50)
      print(min(col.inds, na.rm=TRUE))
      print(max(col.inds, na.rm=TRUE))
    } else {
      col.inds <- round((R_mat)/log10(100)*100)
      x_lab_pos <- -0.75
    }
  # Load the seaborn cubeHelix color palette, and conver the indeces from the
  # and make them bounded between zero and 1
  # if (kdrel) {
    if (deldelGleft) {
      start_col <- 0.4
      r_col <- -0.4
    } else if (deldelG) {
      start_col <- 0.4
      r_col <- -0.4
    } else {
      start_col <- 0
      r_col <- 0.4
    }
  # } else {
  #   start_col <- 0.5
  #   r_col <- -0.75
  # }
  # Make distribution of 100 colors and a 101th gray color for NA indeces.
  if (residual) {
    brbg <- RColorBrewer::brewer.pal(11, "BrBG")
    color.dist <- c(colorRampPalette(c(brbg[2], brbg[6]))(50), 
                    colorRampPalette(c(brbg[6], brbg[10]))(51)[-1], "gray90")
  } else {
    color.dist <- c(rev(cubeHelix(100, start=start_col, r=r_col, hue=0.8)), "gray90", "white")
    # Make the color index scale.
    color.dist[100] <- "red"    
  }
  col.inds <- sapply(t(col.inds), function(col.ind) {
    min(max(1, col.ind), 100)
  })
  col.inds[which(is.na(col.inds))] <- 101
  col.inds[impossible_cols] <- 102

  # Make the color rectangle.
  rect(xlefts, ybottom, xright, ytop, col=color.dist[col.inds], lwd=0.4, 
       xpd=NA, border="white")
  # Make the key.
  if (key) {
    y_div <- (ymax-0.5)/100
    kpl <- xmax + 1 # Left-hand position of the key
    kw <- 1.5                # Width of the key
    rect(xleft=kpl, ybottom=seq(100)*y_div - y_div + 0.5, xright=kpl + kw,
         ytop=seq(100)*y_div - y_div + 0.5 + y_div,
         col=color.dist[1:100], xpd=NA, border=NA)
    # Generate the axis for the legend and the label
    bottom <- ymin + y_div/2
    top <- ymax - y_div/2
    # if (kdrel) {
    labels <- c(0.3, 1, 3, 10, 30, 100)
    # } else {
      # if (site %in% c("8mer", "7mer-m8", "7mer-A1", "6mer")) {
      #   labels <- c(3, 1, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001, 0.0003)
      # } else {
        # labels <- c(3, 1, 0.3, 0.1, 0.03, 0.01, 0.003)
      # }
    # }
    pos_labels <- log(labels)
    centered_labels <- pos_labels - pos_labels[1]
    norm_labels <- centered_labels/(centered_labels[length(centered_labels)])
    # labels_span <- max(log(labels)) - min(log(labels))
    height_span <- top - bottom
    pos_labels <- norm_labels*height_span
    axis(4, at=pos_labels, labels=labels, lwd=par()$lwd, pos=kpl + kw, xpd=NA)
    # if (kdrel) {
      text(x=kpl + kw + 3, y=ymax/2,
           labels=bquote(italic(K)[D]*.(" fold change")),
           srt=270, xpd=NA)    
    # } else {
    #   text(x=kpl + kw + 3.5, y=ymax/2,
    #        labels=bquote(.("Relative")~italic(K)[D]),
    #        srt=270, xpd=NA)    
    # }
  }
  # # Add the label indicating how many nucleotides of pairing.
  mtext(text=sprintf("%s nt offset", offset_sim), side=3, at=xmax, adj=1,
        cex=par("cex"))
  # If the `mirna_label` conditional is true, add the label saying which miRNA
  # is being looked at.
  if (mirna_label) {
    if (mirna == "let-7a-21nt") {
      mirna_txt <- "let-7a"
    } else if (mirna == "let-7a_minus1") {
      mirna_txt <- "let-7a(-1)"
    } else if (mirna == "let-7a_plus1") {
      mirna_txt <- "let-7a(+1)"
    } else if (mirna == "let-7a_miR-155") {
      mirna_txt <- "let-7a-miR-155"
    } else if (mirna == "miR-155_let-7a") {
      mirna_txt <- "miR-155-let-7a"
    } else {
      mirna_txt <- mirna
    }
    mtext(text=mirna_txt, side=3, line=0.7, at=xmin, adj=0, cex=par("cex"))
  }
  if (extralabel) {
    str.end3prand <- "r"
    if (collapsemm) {
      str.collapsemm <- "sum"
    } else {
      str.collapsemm <- "geo"
    }
    mtext(text=sprintf("%s_%s", str.end3prand, str.collapsemm), side=3,
          line=0.7, at=xmax, adj=1, cex=par("cex"))    
  }
  # Finish the plot.
  if (class(pdf.plot) == "character") {
    dev.off()
  }
  return(R_mat)
}

PlotOffsetAndPairingModelMatrix <- function(
  mirna, experiment, n_constant=5, offsetmin=0, offsetmax=17,
  sitelist="programmed_suppcomp", collapsemm=FALSE,
  deldelG=FALSE, deldelGleft=FALSE, key=FALSE, xlabels=TRUE, mirna_label=TRUE,
  height=2.8, width=2.5, xpos=20, ypos=20, pdf.plot=FALSE
) {
  # Load the data matrix.
  model <- SubfunctionCall(FitOffsetAndPairingModel)
  R_mat <- t(model[["pairing"]])

  mar1 <- 2.2
  mar2 <- 2.5
  mar3 <- 1.8
  mar4 <- 0.5
  # This is the conversion factor in order to have the correct height of the
  # plot for the specified width, that makes each box of equal height and width.
  if (deldelGleft) {
    R_mat1 <- R_mat[-1, ]
    R_mat2 <- R_mat[-nrow(R_mat), ]
    R_mat <- t((R_mat2 - R_mat1)[, -1])

  } else if (deldelG) {
    R_mat1 <- R_mat[, -ncol(R_mat)]
    R_mat2 <- R_mat[, -1]
    R_mat <- (R_mat2 - R_mat1)[-nrow(R_mat), ]
    R_mat_del <<- R_mat
  }
  mir_length <- nchar(kMirnaSeqs[mirna])
  SubfunctionCall(FigureSaveFile2)
  par(mar=c(mar1, mar2, mar3, mar4))
  xmin <- 0
  xmax <- nrow(R_mat)
  ymin <- 0
  ymax <- ncol(R_mat) 
  BlankPlot()
  ymin <- 0
  test_rows <- rowMeans(R_mat, na.rm=TRUE)
  test_rows <<- test_rows
  test_row_check <- sapply(1:(length(test_rows) - 1), function(ind) {
    mean(test_rows[ind:(ind+1)], na.rm=TRUE)  
  })
  xlefts <- rep(seq(0, ncol(R_mat) - 1), nrow(R_mat))
  xright <- xlefts + 1
  ybottom <- rep(seq(ymax - 1, ymin), each=ncol(R_mat))
  ytop <- ybottom + 1

  # Make the x-axis
  AddLinearAxis(2, 1, 1, label="5'-paired nt",
                alt_lab=rev(rownames(R_mat)),
                alt_lab_pos=ymin:(ymax - 1) + 0.5,
                alt_tick_pos=TRUE,
                line=0.8)
  AddLinearAxis(1, 1, 1, label="3'-paired nt",
                label_pos_ticks=TRUE,
                alt_lab=colnames(R_mat),
                alt_lab_y_dist=0.01,
                alt_lab_pos=xmin:(xmax - 1) + 0.5,
                alt_tick_pos=TRUE)
  # Add the label for the seed if not plotting relative to seed kds.
  # if (kdrel) {
    message(sprintf("This is the max Kd fold change for %s: %s", mirna, 10^max(R_mat, na.rm=TRUE)))
    message(sprintf("This is the min Kd fold change for %s: %s", mirna, 10^min(R_mat, na.rm=TRUE)))
    if (deldelGleft || deldelG) {
      col.inds <- round((R_mat)/log10(20)*100)
      x_lab_pos <- -0.75
    } else {
      col.inds <- round((R_mat)/log10(100)*100)
      x_lab_pos <- -0.75
    }
  # } else {
  #   col.inds <- round((R_mat - log10(1/3))/3*100)
  #   x_lab_pos <- -2.5
  # }
  # Load the seaborn cubeHelix color palette, and conver the indeces from the
  # and make them bounded between zero and 1
  # if (kdrel) {
    if (deldelGleft) {
      start_col <- 0.4
      r_col <- -0.4
    } else if (deldelG) {
      start_col <- 0.4
      r_col <- -0.4
    } else {
      start_col <- 0
      r_col <- 0.4
    }
  # } else {
  #   start_col <- 0.5
  #   r_col <- -0.75
  # }
  # Make distribution of 100 colors and a 101th gray color for NA indeces.
  color.dist <- c(rev(cubeHelix(100, start=start_col, r=r_col, hue=0.8)), "white")
  # Make the color index scale.
  col.inds <- sapply(t(col.inds), function(col.ind) {
    min(max(1, col.ind), 100)
  })
  col.inds[which(is.na(col.inds))] <- 101
  # Make the color rectangle.
  rect(xlefts, ybottom, xright, ytop, col=color.dist[col.inds], lwd=0.4, 
       xpd=NA, border="white")
  # Make the key.
  if (key) {
    y_div <- (ymax-0.5)/100
    kpl <- xmax + 1 # Left-hand position of the key
    kw <- 1.5                # Width of the key
    rect(xleft=kpl, ybottom=seq(100)*y_div - y_div + 0.5, xright=kpl + kw,
         ytop=seq(100)*y_div - y_div + 0.5 + y_div,
         col=color.dist[1:100], xpd=NA, border=NA)
    # Generate the axis for the legend and the label
    bottom <- ymin + y_div/2
    top <- ymax - y_div/2
    # if (kdrel) {
    labels <- c(0.3, 1, 3, 10, 30, 100)
    # } else {
      # if (site %in% c("8mer", "7mer-m8", "7mer-A1", "6mer")) {
      #   labels <- c(3, 1, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001, 0.0003)
      # } else {
        # labels <- c(3, 1, 0.3, 0.1, 0.03, 0.01, 0.003)
      # }
    # }
    pos_labels <- log(labels)
    centered_labels <- pos_labels - pos_labels[1]
    norm_labels <- centered_labels/(centered_labels[length(centered_labels)])
    # labels_span <- max(log(labels)) - min(log(labels))
    height_span <- top - bottom
    pos_labels <- norm_labels*height_span
    axis(4, at=pos_labels, labels=labels, lwd=par()$lwd, pos=kpl + kw, xpd=NA)
    # if (kdrel) {
      text(x=kpl + kw + 3, y=ymax/2,
           labels=bquote(italic(K)[D]*.(" fold change")),
           srt=270, xpd=NA)    
    # } else {
    #   text(x=kpl + kw + 3.5, y=ymax/2,
    #        labels=bquote(.("Relative")~italic(K)[D]),
    #        srt=270, xpd=NA)    
    # }
  }
  # # Add the label indicating how many nucleotides of pairing.
  # mtext(text=sprintf("Model estimates", offset_sim), side=3, at=xmax, adj=1,
  #       cex=par("cex"))
  # If the `mirna_label` conditional is true, add the label saying which miRNA
  # is being looked at.
  if (mirna_label) {
    if (mirna == "let-7a-21nt") {
      mirna_txt <- "let-7a"
    } else if (mirna == "let-7a_minus1") {
      mirna_txt <- "let-7a(-1)"
    } else if (mirna == "let-7a_plus1") {
      mirna_txt <- "let-7a(+1)"
    } else if (mirna == "let-7a_miR-155") {
      mirna_txt <- "let-7a-miR-155"
    } else if (mirna == "miR-155_let-7a") {
      mirna_txt <- "miR-155-let-7a"
    } else {
      mirna_txt <- mirna
    }
    mtext(text=mirna_txt, side=3, line=0.7, at=xmin, adj=0, cex=par("cex"))
  }

  # Finish the plot.
  if (class(pdf.plot) == "character") {
    dev.off()
  }
  return(R_mat)
}

PlotSingleOffsetCoefficients <- function(
  mirna, experiment, n_constant=5, offsetmin=-4, offsetmax=16,
  sitelist="programmed_suppcomp", outto11mers=TRUE, downto2mers=FALSE,
  downto3mers=FALSE, downto4mers=TRUE, end3prand=TRUE, collapsemm=FALSE,
  supp_base=FALSE, xpos=20, ypos=20, height=4, width=4, pdf.plot=NULL
) {
  model <- SubfunctionCall(FitOffsetAndPairingModel)
  offsets <- c(0, model[["offsets"]])
  names(offsets)[1] <- sprintf("offset%s", offsetmin)
  # colnames(offsets_df) <- mirna_tables[, 1]
  SubfunctionCall(FigureSaveFile2)
  xmin <- offsetmin
  xmax <- offsetmax
  ymin <- 0.3
  ymax <- 10
  BlankPlot(, log="y")



  # if (kdrel) {
    start_col <- 0
    r_col <- 0.4
  # } else {
  #   start_col <- 0.5
  #   r_col <- -0.75
  # }
  # Make distribution of 100 colors and a 101th gray color for NA indeces.
  # if (kdrel) {
    col.inds <- round((offsets - log10(0.3))/log10(100/0.3)*100)  
    # col.inds <- round((log10(data))/log10(100)*100)  
    # x_lab_pos <- -1  
  # } else {
  #   if (site %in% c("8mer", "7mer-m8", "7mer-A1", "6mer")) {
  #     col.inds <- round((log10(data) - log10(1/3))/4*100)
  #   } else {
  #     col.inds <- round((log10(data) - log10(1/3))/3*100)    
  #   }
  # }

  # Make distribution of 100 colors and a 101th gray color for NA indeces.
  color.dist <- c(rev(cubeHelix(100, start=start_col, r=r_col, hue=0.8)), "gray")
  # Make the color index scale.
  col.inds <- sapply(col.inds, function(col.ind) {
    min(max(1, col.ind), 100)
  })





  lines(xmin:xmax, 10^offsets, col="gray")  
  points(xmin:xmax, 10^offsets, col=color.dist[col.inds], pch=19)  
  # Add the axes.
  AddLinearAxis(1, 1, 1, label="Offset (nt)")
  AddLogAxis(2, 1, 1, label="Kd fold change")

  # xy <- GetPlotFractionalCoords(0.5, 0.5)
  # Legend(xy, legend=colnames(offsets_df), col=cols)

  # Finish the plot.
  if (class(pdf.plot) == "character") {
    dev.off()
  }

}

PlotScanningThreePrimeSite <- function(
  mirna="let-7a-21nt", experiment="equil_c2_nb", n_constant=0,
  sitelist="programmed_suppcomp", starting_nuc=10, ending_nuc=NULL,
  collapsemm=TRUE, equilibrium_nb=FALSE,
  height=5, width=5, xpos=20, ypos=20, pdf.plot=FALSE
) {
  kmers <- 4:11
  if (class(ending_nuc) == "numeric") {
    thrp_site_list <- sprintf("%smer-m%s.%s", kmers, ending_nuc - kmers + 1,
                              rep(ending_nuc, length(kmers)))
  } else {
    thrp_site_list <- sprintf("%smer-m%s.%s", kmers,
                              rep(starting_nuc, length(kmers)),
                              starting_nuc + kmers - 1)    
  }
  mir_length = nchar(kMirnaSeqs[mirna])
  mir_ends <- sapply(thrp_site_list, function(site) {
    as.integer(gsub("^.*\\.(.*)$", replacement="\\1", site, perl=TRUE))
  })
  thrp_site_list <- thrp_site_list[which(mir_ends <= mir_length)]
  kds <- SubfunctionCall(EquilPars, sitelist="programmed_suppcomp")
  kds <<- kds
  SubfunctionCall(FigureSaveFile2)
  xmin <- 8
  xmax <- 25
  ymin <- 0.0003
  ymax <- 10
  BlankPlot(log="y")

  AddLinearAxis(1, tick.space=1, label.space=5,
                label="Position with respect to miRNA")
  AddLogAxis(2, label="Relative Kd")
  # mm_ind <- sprintf("_Kd", seed_mmsite)
  if (sitelist == "programmed_suppcomp" & collapsemm) {
    print(kds["Comp_Kd", 2])
    segments(x0=xmin, kds["Comp_Kd", 2], x1=xmax, lty=2, lwd=0.5)
  } else {
    print(kds[grep("^8mer-mm[ACTG][2-7]_Kd$",
                                       rownames(kds),
                                       perl=TRUE), 2])
    segments(x0=xmin, GeoMean(kds[grep("^8mer-mm[ACTG][2-7]_Kd$",
                                       rownames(kds),
                                       perl=TRUE), 2]), x1=xmax, lty=2, lwd=0.5)    
  }
  cols <- c(`5mer-m11.15`="red", `6mer-m11.16`="green",
            `7mer-m11.17`="forestgreen", `8mer-m11.18`="blue",
            `9mer-m11.19`="purple")
  cols <- rev(viridis::plasma(10))[3:10]
  cols <- rev(viridis::plasma(13))[4:11]
  
  names(cols) <- thrp_site_list
  print(cols)
  offsets <- seq(-0.1, 0.1, length.out=length(kmers))
  names(offsets) <- names(cols)
  sapply(rev(thrp_site_list), function(site) {
    print(site)
    site_grep <- gsub("\\.", replacement="\\\\\\.", site)
    site_grep <- gsub("\\|", replacement="\\\\\\|", site_grep)
    print(site_grep)
    # site_kds <- SubfunctionCall(GetPositionalProgKds)
    inds <- grep(sprintf("^%s\\|.*\\|Comp_Kd", site_grep), rownames(kds), perl=TRUE)
    print(inds)
    df_use <- kds[inds, ]
    print(df_use)
    pos <- gsub("^(.*)\\|(.*)\\|(Comp_Kd$)", replacement="\\2", rownames(df_use))
    print(pos)
    rownames(df_use) <- pos
    x <- 9:min(max(as.integer(pos)), xmax)
    print(df_use)
    # y_use <- which(as.integer(pos) <= 9 & as.integer(pos) <=)
    print(cols[site])
    segments(x0=x + offsets[site], y0=df_use[as.character(x), 3],
             y1=df_use[as.character(x), 5],
             col=ConvertRColortoRGB(cols[site], alpha=1), lwd=1, xpd=NA)
    lines(x, df_use[as.character(x), 2], col=cols[site], lwd=1, xpd=NA)
    points(x, df_use[as.character(x), 2], col=cols[site], pch=20, xpd=NA)
  })

  if (sitelist == "programmed_suppcomp" & collapsemm) {
    s8mer_kd <- GeoMean(kds[grep("^8mer\\|.*\\|Comp_Kd$",
                                 rownames(kds), perl=TRUE), 2])
    s6mer_kd <- GeoMean(kds[grep("^6mer\\|.*\\|Comp_Kd$",
                                 rownames(kds), perl=TRUE), 2])
    mm_kd <- GeoMean(kds[grep("^Comp_Kd$",
                              rownames(kds), perl=TRUE), 2])
  } else {
    s8mer_kd <- GeoMean(kds[grep("^8mer\\|.*\\|Comp_Kd$",
                                 rownames(kds), perl=TRUE), 2])
    s6mer_kd <- GeoMean(kds[grep("^6mer\\|.*\\|Comp_Kd$",
                                 rownames(kds), perl=TRUE), 2])
    mm_kd <- GeoMean(kds[grep("^8mer-mm[ACGT][2-7]_Kd$",
                              rownames(kds), perl=TRUE), 2])    
  }
  s8mer_kd_norm <- s8mer_kd*mm_kd/(s8mer_kd + mm_kd)
  s6mer_kd_norm <- s6mer_kd*mm_kd/(s6mer_kd + mm_kd)
  if (mirna %in% c("let-7a-21nt", "let-7a_plus1", "let-7a_minus1",
                   "let-7a_miR-155")) {
    if (equilibrium_nb) {
      mirna_rand <- "let-7a-21nt"
      experiment_rand <- "equilibrium_nb"
    } else {
      mirna_rand <- "let-7a"
      experiment_rand <- "equilibrium"
    }
  } else if (mirna == "miR-155_let-7a") {
    mirna_rand <- "miR-155"
    experiment_rand <- "equilibrium"
  } else {
    mirna_rand <- mirna
    experiment_rand <- "equilibrium"
  }
  if (mirna == "miR-1") {
    buffer_rand <- TRUE
    combined_rand <- FALSE
  } else {
    buffer_rand <- FALSE
    combined_rand <- TRUE
  }
  kds_rand <- EquilPars(mirna_rand, experiment_rand, n_constant=5,
                              sitelist="programmed", buffer=buffer_rand,
                              combined=combined_rand)
  s8mer_kd_rand <- kds_rand["8mer_Kd", 2]
  s6mer_kd_rand <- kds_rand["6mer_Kd", 2]

  segments(xmin, s8mer_kd, x1=xmax, y1=s8mer_kd,
           col=kSiteColors["8mer"])
  segments(xmin, s6mer_kd, x1=xmax, y1=s6mer_kd,
           col=kSiteColors["6mer"])

  # segments(xmin, s8mer_kd_norm, x1=xmax, y1=s8mer_kd_norm, lty=2,
  #          col=kSiteColors["8mer"], lwd=3)
  # segments(xmin, s6mer_kd_norm, x1=xmax, y1=s6mer_kd_norm, lty=2,
  #          col=kSiteColors["6mer"], lwd=3)

  # segments(xmin, s8mer_kd_rand, x1=xmax, y1=s8mer_kd_rand, lty=3,
  #          col=kSiteColors["8mer"], lwd=3)
  # segments(xmin, s6mer_kd_rand, x1=xmax, y1=s6mer_kd_rand, lty=3,
  #          col=kSiteColors["6mer"], lwd=3)

  xy <- GetPlotFractionalCoords(0, 1.05, log="y")
  Legend(xy, legend=names(cols), col=cols, ncol=2)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

PlotRegisterAndPositionKdFoldChange <- function(
  mirna, experiment, len_k=7, n_constant=5, sitelist="programmed", kdrel=FALSE,
  key=FALSE, xlabels=TRUE, mirna_label=FALSE, width=2.85, xpos=20,
  ypos=20, pdf.plot=FALSE
) {
  # Load the collapsed data (for the mismatch-only).
  mm_8mer_sites <- GetAll8merMmSites(mirna)
  kds_collapsed <- SubfunctionCall(EquilPars, sitelist="programmed_collapsed")
  print(dim(kds_collapsed))
  len_mir <- nchar(kMirnaSeqs[mirna])
  num_sites <- len_mir - 8 - len_k + 1
  print(num_sites)
  starts <- 9:(9 + num_sites - 1)
  stops <- starts + len_k - 1
  sites <- sprintf("%smer-m%s.%s", len_k, starts, stops)
    vals <- t(sapply(sites, function(site) {
      kds <- SubfunctionCall(GetPositionalProgKds)
      num_pos <- 25 - len_k + 1
      mir_pos <- c("Seed", as.character(9:(9 + num_pos - 1)))

      kds_out <- matrix(0, nrow=18, ncol=1 + 25 - len_k + 1,
                        dimnames=list(mm_8mer_sites, mir_pos))
      for (mm_site in mm_8mer_sites) {
        for (pos in mir_pos) {
          if (pos == "Seed") {
            site_name <- sprintf("%s_Kd", mm_site)

            val <- kds_collapsed[site_name, 2]
          } else {
            site_name <- sprintf("%s|%s|%s_Kd", site, pos, mm_site)
            val <- kds[site_name, 2]
          }
          kds_out[mm_site, as.character(pos)] <- val
        }
      }
      if (kdrel) {
        kds_out <- (t(t(kds_out/kds_out[, 1])))
      }
      kds_out <- kds_out[, 2:ncol(kds_out)]
      return(exp(colMeans(log(kds_out), na.rm=TRUE)))
  }))
  mar1 <- 2
  mar2 <- 3
  mar3 <- 2
  if (key) {
    mar4 <- 6
  } else {
    mar4 <- 0.5
  }
  # This is the conversion factor in order to have the correct height of the
  # plot for the specified width, that makes each box of equal height and width.
  mir_length <- nchar(kMirnaSeqs[mirna])
  height <- (mar1 + mar3 + 1.5*(mir_length - 8 - len_k + 1)*(5*width - mar2 - mar4)/17)/5  
  SubfunctionCall(FigureSaveFile2)
  par(mar=c(mar1, mar2, mar3, mar4))
  xmin <- 0
  xmax <- 16      
  ymin <- 0
  ymax <- mir_length - 8 - len_k + 1
  BlankPlot()
  ymin <- 0
  R_mat <- -log10(vals)
  R_mat <- R_mat[, 1:17]

  R_mat <<- R_mat
  test_rows <- rowMeans(R_mat, na.rm=TRUE)
  test_rows <<- test_rows
  test_row_check <- sapply(1:(length(test_rows) - 1), function(ind) {
    mean(test_rows[ind:(ind+1)], na.rm=TRUE)  
  })
  # print(test_row_check)
  # print(which(test_row_check == max(test_row_check)))
  # return()
  xlefts <- rep(seq(0, ncol(R_mat) - 1), ymax) - 0.5
  xright <- xlefts + 1
  ybottom <- rep(seq(ymax - 1, ymin), each=ncol(R_mat))
  ytop <- ybottom + 1

  # sapply(1:5, function(y) {
  #   segments(x0=xmin, y0=ymin - y, x1=xmax, col="red", xpd=NA)
  #   segments(x0=xmin, y0=ymax + y, x1=xmax, col="red", xpd=NA)
  #   segments(x0=xmin - y, y0=ymin, y1=ymax, col="red", xpd=NA)
  #   segments(x0=xmax + y, y0=ymin, y1=ymax, col="red", xpd=NA)
  # })
  # Convert the R values within into normalized values

  x_lib_labels <- c(9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29)
  x_lib_pos    <- c(0,  2,  4,  6,  8, 10, 12, 14, 16, 18, 20)
  # if (!kdrel) {
  #   # x_lib_labels <- c("Seed", x_lib_labels)
  #   x_lib_pos <- x_lib_pos + 2
  # }
  # Make the x-axis
  if (xlabels) {
    AddLinearAxis(1, 1, 1, label="Position", blank_lab=TRUE, pos=ymin)

  }
  AddLinearAxis(1, 1, 1, label=NA, label_pos_ticks=TRUE, alt_lab=x_lib_labels,
                alt_lab_pos=x_lib_pos, line=1.5, alt_tick_pos=TRUE, pos=ymin)
  # Add the label for the seed if not plotting relative to seed kds.
  if (!kdrel) {
    text(-1.5, ymin - 0.5, labels="Seed", srt=90, adj=c(1, 0.5), xpd=NA)
  }
  if (kdrel) {
    # col.inds <- round((R_mat - log10(0.3))/log10(100/0.3)*100)  
    col.inds <- round((R_mat)/log10(100)*100)  
    x_lab_pos <- -0.75  
  } else {
    col.inds <- round((R_mat - log10(1/3))/3*100)    
    x_lab_pos <- -2.5
  }
  col.inds <<- col.inds
  text(x=x_lab_pos, y=(ymax - 1):ymin + 0.5, labels=gsub("\\.", replacement="-",
                                        gsub("^.*mer-m",
                                             replacement="", rownames(R_mat),
                                             perl=TRUE)),
       adj=c(1, 0.5), xpd=NA)
  # Load the seaborn cubeHelix color palette, and conver the indeces from the
  # and make them bounded between zero and 1
  if (kdrel) {
    start_col <- 0
    r_col <- 0.4
  } else {
    start_col <- 0.5
    r_col <- -0.75
  }
  # Make distribution of 100 colors and a 101th gray color for NA indeces.
  color.dist <- c(rev(cubeHelix(100, start=start_col, r=r_col, hue=0.8)), "gray")
  # Make the color index scale.
  col.inds <- sapply(t(col.inds), function(col.ind) {
    min(max(1, col.ind), 100)
  })
  col.inds[which(is.na(col.inds))] <- 101
  # Make the color rectangle.
  rect(xlefts, ybottom, xright, ytop, col=color.dist[col.inds], lwd=0.20, 
       xpd=NA, border="white")
  # Make the key.
  if (key) {
    y_div <- (ymax-0.5)/100
    kpl <- xmax + 1 # Left-hand position of the key
    kw <- 1.5                # Width of the key
    rect(xleft=kpl, ybottom=seq(100)*y_div - y_div + 0.5, xright=kpl + kw,
         ytop=seq(100)*y_div - y_div + 0.5 + y_div,
         col=color.dist[1:100], xpd=NA, border=NA)
    # Generate the axis for the legend and the label
    bottom <- ymin + y_div/2
    top <- ymax - y_div/2
    if (kdrel) {
      labels <- c(0.3, 1, 3, 10, 30, 100)
    } else {
      # if (site %in% c("8mer", "7mer-m8", "7mer-A1", "6mer")) {
      #   labels <- c(3, 1, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001, 0.0003)
      # } else {
        labels <- c(3, 1, 0.3, 0.1, 0.03, 0.01, 0.003)
      # }
    }
    pos_labels <- log(labels)
    centered_labels <- pos_labels - pos_labels[1]
    norm_labels <- centered_labels/(centered_labels[length(centered_labels)])
    # labels_span <- max(log(labels)) - min(log(labels))
    height_span <- top - bottom
    pos_labels <- norm_labels*height_span
    axis(4, at=pos_labels, labels=labels, lwd=par()$lwd, pos=kpl + kw, xpd=NA)
    if (kdrel) {
      text(x=kpl + kw + 3, y=ymax/2,
           labels=bquote(italic(K)[D]*.(" fold change")),
           srt=270, xpd=NA)    
    } else {
      text(x=kpl + kw + 3.5, y=ymax/2,
           labels=bquote(.("Relative")~italic(K)[D]),
           srt=270, xpd=NA)    
    }
  }
  # Add the label indicating how many nucleotides of pairing.
  mtext(text=sprintf("%s nt of pairing", len_k), side=3, at=xmax, adj=1,
        cex=par("cex"))
  # If the `mirna_label` conditional is true, add the label saying which miRNA
  # is being looked at.
  if (mirna_label) {
    if (mirna == "let-7a-21nt") {
      mirna_txt <- "let-7a"
    } else if (mirna == "let-7a_minus1") {
      mirna_txt <- "let-7a(-1)"
    } else if (mirna == "let-7a_plus1") {
      mirna_txt <- "let-7a(+1)"
    } else {
      mirna_txt <- mirna
    }
    mtext(text=mirna_txt, side=3, line=1, at=xmin, adj=0, cex=par("cex"))
  }
  # Finish the plot.
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

PlotKdMatrixScatter <- function(
  mirna1, experiment1, mirna2, experiment2, site, n_constant=5,
  sitelist="programmed", kdrel=FALSE, height=5, width=5, xpos=20, ypos=20,
  pdf.plot=FALSE
) {
  # Load the collapsed data (for the mismatch-only).
  kds_1 <- SubfunctionCall(GetPositionalProgKds, mirna=mirna1,
                               experiment=experiment1)
  kds_2 <- SubfunctionCall(GetPositionalProgKds, mirna=mirna2,
                               experiment=experiment2)

  SubfunctionCall(FigureSaveFile2)
  xmin <- 1e-3
  xmax <- 10      
  ymin <- 1e-3
  ymax <- 10
  BlankPlot(log='xy')
  # Make the axes
  segments(xmin, ymin, xmax, ymax, col="gray", lty=2)
  AddLogAxis(1, label="Relative Kd; replicate 1")
  AddLogAxis(2, label="Relative Kd; replicate 2")
  x <- kds_1[, 2]
  y <- kds_2[rownames(kds_1), 2]

  cols <- ConvertRColortoRGB(c("red", "orange", "green",
                               "forestgreen", "blue", "purple"), alpha=0.5)
  col_inds <- as.integer(gsub("^.*mm[ACGT]([2-7])_Kd$", replacement="\\1", rownames(kds_1),
                   perl=TRUE)) - 1
  print(head(col_inds))
  Points(x, y, col=cols[col_inds])
  kds_1 <<- kds_1
  kds_2 <<- kds_2
  xy <- GetPlotFractionalCoords(0.05, 0.95, log='xy')
  AddCorrelationToPlot(xy[1], xy[2],
                       x=log(x), y=log(y), rsquared=TRUE, adj=c(0, 1))
  xy <- GetPlotFractionalCoords(0.65, 0.4, log='xy')
  Legend(xy, legend=c("8mer-mm2", "8mer-mm3", "8mer-mm4",
                      "8mer-mm5", "8mer-mm6", "8mer-mm7"), col=cols)
  # Finish the plot.
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

PlotBestRegisterKd <- function(
  len_k=7, average_win=3, n_constant=5, sitelist="programmed",
  let7_series=FALSE, chimera_series=FALSE, chimera_let7_series=FALSE,
  chimera_miR155_series=FALSE, kdrel=FALSE, key=TRUE, xlabels=TRUE,
  height=3, width=4, xpos=20, ypos=20, pdf.plot=FALSE
) {
  # Load the collapsed data (for the mismatch-only).
  if (let7_series) {
    mirna_tables <- matrix(c("let-7a_minus1", "equil_c_nb",
                             "let-7a-21nt", "equil_c2_nb",
                             "let-7a_plus1", "equil_c_nb"), byrow=TRUE, nrow=3)

  } else if (chimera_series) {
    mirna_tables <- matrix(c("let-7a-21nt", "equil_c2_nb",
                             "miR-155", "equil_sc_nb",
                             "miR-155_let-7a", "equil_c_nb",
                             "let-7a_miR-155", "equil_c_nb"), byrow=TRUE,
                           nrow=4)
  } else if (chimera_let7_series) {
    mirna_tables <- matrix(c("let-7a-21nt", "equil_c2_nb",
                             "miR-155_let-7a", "equil_c_nb"), byrow=TRUE,
                           nrow=2)
  } else if (chimera_miR155_series) {
    mirna_tables <- matrix(c("miR-155", "equil_sc_nb",
                             "let-7a_miR-155", "equil_c_nb"), byrow=TRUE,
                           nrow=2)
  } else {
    mirna_tables <- matrix(c("let-7a-21nt", "equil_c2_nb",
                             "miR-1", "equil_c_nb",
                             "miR-155", "equil_sc_nb"), byrow=TRUE, nrow=3)    
  }
  out <- apply(mirna_tables, 1, function(row) {
    mirna <- row[1]
    experiment <- row[2]
    mm_8mer_sites <- GetAll8merMmSites(mirna)

    kds_collapsed <- SubfunctionCall(EquilPars, sitelist="programmed_collapsed")
    print(dim(kds_collapsed))
    len_mir <- nchar(kMirnaSeqs[mirna])
    num_sites <- len_mir - 8 - len_k + 1
    # print(num_sites)
    starts <- 9:(9 + num_sites - 1)
    # print(starts)
    stops <- starts + len_k - 1
    # print(stops)
    sites <- sprintf("%smer-m%s.%s", len_k, starts, stops)
    print(sites)
    vals <- sapply(sites, function(site) {
      kds <- SubfunctionCall(GetPositionalProgKds)

      num_pos <- 25 - len_k + 1
      mir_pos <- c("Seed", as.character(9:(9 + num_pos - 1)))

      kds_out <- matrix(0, nrow=18, ncol=1 + 25 - len_k + 1,
                        dimnames=list(mm_8mer_sites, mir_pos))
      for (mm_site in mm_8mer_sites) {
        for (pos in mir_pos) {
          if (pos == "Seed") {
            site_name <- sprintf("%s_Kd", mm_site)

            val <- kds_collapsed[site_name, 2]
          } else {
            site_name <- sprintf("%s|%s|%s_Kd", site, pos, mm_site)
            val <- kds[site_name, 2]
          }
          kds_out[mm_site, as.character(pos)] <- val
        }
      }
      if (kdrel) {
        kds_out <- (t(t(kds_out/kds_out[, 1])))^(-1)
      }
      kds_out <- kds_out[, 2:ncol(kds_out)]
      kds_average <- exp(colMeans(log(kds_out), na.rm=TRUE))
      mean_val <- sapply(1:(length(kds_average) - average_win + 1), function(inds) {
        exp(mean(log(kds_average[inds:(inds + average_win - 1)]), na.rm=TRUE))
      })
      if (kdrel) {
        return(max(mean_val, na.rm=TRUE))
      } else {
        return(min(mean_val, na.rm=TRUE))      
      }
    })
    vals
  })
  if (class(out) == "matrix") {
    out <- split(out, col(out))
  }
  print(out)
  if (let7_series) {
    names(out) <- c("let-7a(-1)", "let-7a", "let-7a(+1)")
    xmax <- 21 - len_k + 1.2
    cols <- c("dodgerblue", "#F7931D", "darkorchid3")
  } else if (chimera_series) {
    names(out) <- c("let-7a", "miR-155", "miR-155-let-7a", "let-7a-miR-155")
    xmax <- 23 - len_k + 1.2
    cols <- c("#F7931D", "darkorchid3", "dodgerblue", "deeppink")
  } else if (chimera_let7_series) {
    names(out) <- c("let-7a", "miR-155-let-7a")
    xmax <- 23 - len_k + 1.2
    cols <- c("dodgerblue", "deeppink")
  } else if (chimera_miR155_series) {
    names(out) <- c("miR-155", "let-7a-miR-155")
    xmax <- 23 - len_k + 1.2
    cols <- c("darkorchid3", "deeppink")
  } else {
    names(out) <- c("let-7a", "miR-1", "miR-155")  
    xmax <- 23 - len_k + 1.2
    cols <- c("#F7931D", "dodgerblue", "darkorchid3")
  }
  names(cols) <- names(out)
  SubfunctionCall(FigureSaveFile2)
  xmin <- 9
  if (kdrel) {
    ymin <- 0.5
    if (chimera_series || chimera_let7_series || chimera_miR155_series) {
      ymax <- 20
    } else {
      ymax <- 100    
    }
    inv <- ""
    xlab <- "Kd fold change"
  } else {
    ymin <- 1e-2
    ymax <- 3
    inv <- "y"
    xlab <- "Relative Kd"
  }
  par(mar=c(3, 3, 2, 1))
  BlankPlot(log="y", inv=inv)
  xmax <- xmax - 0.2
  AddLogAxis(2, label=xlab)
  AddLinearAxis(1, 1, 1, "5'-most paired miRNA nucleotide")
  # xmax <- 17.2
  print(out)
  segments(x0=xmin, y0=1, x1=xmax, lty=2, col="gray")
  sapply(names(out), function(mirna) {
    y_vals <- out[[mirna]]
    lines(1:length(y_vals) + 8, y_vals, col=cols[mirna])
    points(1:length(y_vals) + 8, y_vals, col=cols[mirna], pch=19)
  })
  if (chimera_series || chimera_let7_series || chimera_miR155_series) {
    xy <- GetPlotFractionalCoords(0.4, 1.2, log="y", inv=inv)
  } else {
    xy <- GetPlotFractionalCoords(0.6, 1.1, log="y", inv=inv)
  }
  Legend(xy, legend=names(out), lwd=1, col=cols)
  xy <- GetPlotFractionalCoords(1, 0.025, log="y", inv=inv)
  text(xy[1], xy[2], labels=sprintf("%s nt of pairing", len_k), adj=c(1, 0),
       xpd=NA)
  if (chimera_let7_series || chimera_miR155_series) {
    xy <- GetPlotFractionalCoords(0.025, 0.025, log="y", inv=inv)
    AddCorrelationToPlot(log(out[[1]]), log(out[[2]]), xpos=xy[1], ypos=xy[2],
                         rsquared=TRUE, adj=c(0, 0))
  }
  # Finish the plot.
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

PlotKdMatrixCollapsed <- function(
  mirna, experiment, site, n_constant=5, lambda=0.01,
  ordermm=FALSE, fulloptim=FALSE, kdrel=FALSE, key=TRUE, xlabels=TRUE,
  outto11mers=FALSE, downto2mers=FALSE, downto3mers=FALSE, downto4mers=FALSE,
  end3prand=FALSE, new=FALSE, height=3, width=3, xpos=20, ypos=20,
  pdf.plot=FALSE
) {
  # Load the collapsed data (for the mismatch-only).
  kds_collapsed <- SubfunctionCall(EquilPars, sitelist="programmed_collapsed")  
  print("################################################################")
  # kds_site <<- kds_site
  # print(kds_site["7mer-m11.17|14|8mer-mmC7_Kd", , drop=FALSE])
  # return()
  mm_8mer_sites <- GetAll8merMmSites(mirna)
  len_site <- as.integer(substring(site, 1, 1))
  num_sites <- 25 - len_site + 1
  mir_pos <- c("Seed", "Collapsed")

  kds_out <- matrix(0, nrow=18, ncol=2,
                    dimnames=list(mm_8mer_sites, mir_pos))
  for (mm_site in mm_8mer_sites) {
    for (pos in mir_pos) {
      if (pos == "Seed") {
        site_name <- sprintf("%s_Kd", mm_site)

        val <- kds_collapsed[site_name, 2]
      } else {
        site_name <- sprintf("%s|%s_Kd", site, mm_site)
        val <- kds_collapsed[site_name, 2]
      }
      kds_out[mm_site, as.character(pos)] <- val
    }
  }
  if (ordermm) {
    kds_out <- kds_out[order(kds_out[, 1]), ]  
  }
  # Make the plot and define the limits.
  SubfunctionCall(FigureSaveFile2)
  if (key) {
    par(mar=c(3, 4, 2, 6))  
  } else {
    par(mar=c(3, 4, 2, 2))  
  }
  xmin <- -2
  xmax <- 2      
  ymin <- 0
  ymax <- 17.5
  BlankPlot()
  # Nucleotide positions to be labeled on the x-axis, relative to the miRNA.
  x_lib_labels <- c(9, 11, 13, 15, 17, 19, 21, 23, 25, 27)
  x_lib_pos <- c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18)
  # if (!kdrel) {
  #   # x_lib_labels <- c("Seed", x_lib_labels)
  #   x_lib_pos <- x_lib_pos + 2
  # }
  # Make the x-axis
  if (xlabels) {
    AddLinearAxis(1, 1, 1, label="Position", blank_lab=TRUE)

  }
  AddLinearAxis(1, 1, 1, label=NA, label_pos_ticks=TRUE,
                alt_lab=x_lib_labels, alt_lab_pos=x_lib_pos, line=1.5, alt_tick_pos=TRUE)
  # Add the label for the seed if not plotting relative to seed kds.
  if (!kdrel) {
    text(-1.5, -0.5, labels="Seed", srt=90, adj=c(1, 0.5), xpd=NA)
  }
  # Make the R matrix map from 0–4 being 10–0.001
  R_mat <- -log10(kds_out)
  if (kdrel) {
    R_mat <- R_mat - R_mat[, 1]
    R_mat <- R_mat[, 2:ncol(R_mat), drop=FALSE]
  }
  # Define the left, right, bottom, and top positions for the squares in the
  # heatmap.
  if (kdrel) {
    xlefts <- rep(seq(0, ncol(R_mat) - 1), ymax) - 0.5
  } else {
    xlefts <- rep(c(-1.5, seq(0, ncol(R_mat) - 2)), ymax) - 0.5
  }
  xright <- xlefts + 1
  ybottom <- rep(seq(17, 0), each=ncol(R_mat))
  ytop <- ybottom + 1
  # Convert the R values within into normalized values

  if (kdrel) {
    col.inds <- round((R_mat - log10(0.3))/log10(100/0.3)*100)  
    x_lab_pos <- -1  
  } else {
    if (site %in% c("8mer", "7mer-m8", "7mer-A1", "6mer")) {
      col.inds <- round((R_mat - log10(1/3))/4*100)
    } else {
      col.inds <- round((R_mat - log10(1/3))/3*100)    
    }
    x_lab_pos <- -2.5
  }
  col.inds <<- col.inds
  if (height == 3) {
    text_cex <- 0.75
  } else {
    text_cex <- 1
  }
  text(x=x_lab_pos, y=17:0 + 0.5, labels=gsub("T", replacement="U",
                                        gsub("8mer-mm",
                                             replacement="", rownames(R_mat))),
       adj=c(1, 0.5), xpd=NA, cex=text_cex)
  # Load the seaborn cubeHelix color palette, and conver the indeces from the
  # and make them bounded between zero and 1
  if (kdrel) {
    start_col <- 0
    r_col <- 0.4
  } else {
    start_col <- 0.5
    r_col <- -0.75
  }
  # Make distribution of 100 colors and a 101th gray color for NA indeces.
  color.dist <- c(rev(cubeHelix(100, start=start_col, r=r_col, hue=0.8)), "gray")
  # Make the color index scale.
  col.inds <- sapply(t(col.inds), function(col.ind) {
    min(max(1, col.ind), 100)
  })
  col.inds[which(is.na(col.inds))] <- 101
  # Make the color rectangle.
  rect(xlefts, ybottom, xright, ytop, col=color.dist[col.inds], lwd=0.20,
       xpd=NA, border="white")
  # Make the key.
  if (key) {
    y_div <- ymax/100
    kpl <- 17 # Left-hand position of the key
    kw <- 1.5                # Width of the key
    rect(xleft=kpl, ybottom=seq(100)*y_div - y_div + 0.5, xright=kpl + kw,
         ytop=seq(100)*y_div - y_div + 0.5 + y_div,
         col=color.dist[1:100], xpd=NA, border=NA)
    # Generate the axis for the legend and the label
    bottom <- ymin + y_div/2
    top <- ymax + 0.5 - y_div/2
    if (kdrel) {
      labels <- c(0.3, 1, 3, 10, 30, 100)
    } else {
      if (site %in% c("8mer", "7mer-m8", "7mer-A1", "6mer")) {
        labels <- c(3, 1, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001, 0.0003)
      } else {
        labels <- c(3, 1, 0.3, 0.1, 0.03, 0.01, 0.003)
      }
    }
    pos_labels <- log(labels)
    centered_labels <- pos_labels - pos_labels[1]
    norm_labels <- centered_labels/(centered_labels[length(centered_labels)])
    # labels_span <- max(log(labels)) - min(log(labels))
    height_span <- top - bottom
    pos_labels <- norm_labels*height_span
    axis(4, at=pos_labels, labels=labels, lwd=par()$lwd, pos=kpl + kw, xpd=NA)
    bquote(italic(r)^2~.(" with flanking dinucleotide")~italic(K)[D]*.("s"))
    if (kdrel) {
      text(x=kpl + kw + 3, y=10,
           labels=bquote(italic(K)[D]*.(" fold change")),
           srt=270, xpd=NA)    
    } else {
      text(x=kpl + kw + 3.5, y=10,
           labels=bquote(.("Relative")~italic(K)[D]*.("s")),
           srt=270, xpd=NA)    
    }
  }
  text(xmin, ymax + 1, labels=site, xpd=NA, adj=c(0, 0))
  # Finish the plot.
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

PlotKdMatrix <- function(
  mirna, experiment, site, n_constant=5, lambda=0.01, sitelist="programmed",
  ordermm=FALSE, fulloptim=FALSE, kdrel=FALSE, key=TRUE, xlabels=TRUE, suppcomp=FALSE,
  outto11mers=FALSE, downto2mers=FALSE, downto3mers=FALSE, downto4mers=FALSE,
  end3prand=FALSE, new=FALSE, height=3, width=3, xpos=20, ypos=20,
  pdf.plot=FALSE
) {
  print(ordermm)
  print(suppcomp)
  # Load the collapsed data (for the mismatch-only).
  if (suppcomp) {
    kds_collapsed <- SubfunctionCall(EquilPars, sitelist="programmed_suppcomp")  
  } else {
    kds_collapsed <- SubfunctionCall(EquilPars, sitelist="programmed_collapsed")  
  }
  if (fulloptim) {
    kds_site <- SubfunctionCall(EquilPars, sitelist="programmed")
  } else {
    kds_site <- SubfunctionCall(GetPositionalProgKds)  
  }
  print("################################################################")
  mm_8mer_sites <- GetAll8merMmSites(mirna)
  len_site <- as.integer(substring(site, 1, 1))
  num_sites <- 25 - len_site + 1
  mir_pos <- c("Seed", as.character(9:(9 + num_sites - 1)))

  kds_out <- matrix(0, nrow=18, ncol=1 + 25 - len_site + 1,
                    dimnames=list(mm_8mer_sites, mir_pos))
  for (mm_site in mm_8mer_sites) {
    for (pos in mir_pos) {
      if (pos == "Seed") {
        if (fulloptim) {
          site_name <- sprintf("NA|NA|%s_Kd", mm_site)
          val <- kds_site[site_name, 2]
        } else {
          site_name <- sprintf("%s_Kd", mm_site)
          val <- kds_collapsed[site_name, 2]
        }
      } else {
        site_name <- sprintf("%s|%s|%s_Kd", site, pos, mm_site)
        val <- kds_site[site_name, 2]
      }
      kds_out[mm_site, as.character(pos)] <- val
    }
  }
  if (ordermm) {
    kds_out <- kds_out[order(kds_out[, 1]), ]  
  }
  kds_out2 <<- kds_out
  # Make the plot and define the limits.
  SubfunctionCall(FigureSaveFile2)
  if (key) {
    par(mar=c(3, 4, 2, 6))  
  } else {
    par(mar=c(3, 4, 2, 2))  
  }
  xmin <- 0
  xmax <- 16      
  ymin <- 0
  ymax <- 17.5
  BlankPlot()
  # Nucleotide positions to be labeled on the x-axis, relative to the miRNA.
  x_lib_labels <- c(9, 11, 13, 15, 17, 19, 21, 23, 25, 27)
  x_lib_pos <- c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18)
  # if (!kdrel) {
  #   # x_lib_labels <- c("Seed", x_lib_labels)
  #   x_lib_pos <- x_lib_pos + 2
  # }
  # Make the x-axis
  if (xlabels) {
    AddLinearAxis(1, 1, 1, label="Position", blank_lab=TRUE)

  }
  AddLinearAxis(1, 1, 1, label=NA, label_pos_ticks=TRUE,
                alt_lab=x_lib_labels, alt_lab_pos=x_lib_pos, line=1.5, alt_tick_pos=TRUE)
  # Add the label for the seed if not plotting relative to seed kds.
  if (!kdrel) {
    text(-1.5, -0.5, labels="Seed", srt=90, adj=c(1, 0.5), xpd=NA)
  }
  # Make the R matrix map from 0–4 being 10–0.001
  R_mat <- -log10(kds_out[, 1:18])
  if (kdrel) {
    R_mat <- R_mat - R_mat[, 1]
    R_mat <- R_mat[, 2:ncol(R_mat)]
  }
  # Define the left, right, bottom, and top positions for the squares in the
  # heatmap.
  if (kdrel) {
    xlefts <- rep(seq(0, ncol(R_mat) - 1), ymax) - 0.5
  } else {
    xlefts <- rep(c(-1.5, seq(0, ncol(R_mat) - 2)), ymax) - 0.5
  }
  xright <- xlefts + 1
  ybottom <- rep(seq(17, 0), each=ncol(R_mat))
  ytop <- ybottom + 1
  # Convert the R values within into normalized values

  if (kdrel) {
    col.inds <- round((R_mat - log10(0.3))/log10(100/0.3)*100)  
    x_lab_pos <- -1  
  } else {
    if (site %in% c("8mer", "7mer-m8", "7mer-A1", "6mer")) {
      col.inds <- round((R_mat - log10(1/3))/4*100)
    } else {
      col.inds <- round((R_mat - log10(1/3))/3*100)    
    }
    x_lab_pos <- -2.5
  }
  col.inds <<- col.inds
  if (height == 3) {
    text_cex <- 0.75
  } else {
    text_cex <- 1
  }
  text(x=x_lab_pos, y=17:0 + 0.5, labels=gsub("T", replacement="U",
                                        gsub("8mer-mm",
                                             replacement="", rownames(R_mat))),
       adj=c(1, 0.5), xpd=NA, cex=text_cex)
  # Load the seaborn cubeHelix color palette, and conver the indeces from the
  # and make them bounded between zero and 1
  if (kdrel) {
    start_col <- 0
    r_col <- 0.4
  } else {
    start_col <- 0.5
    r_col <- -0.75
  }
  # Make distribution of 100 colors and a 101th gray color for NA indeces.
  color.dist <- c(rev(cubeHelix(100, start=start_col, r=r_col, hue=0.8)), "gray90")
  # Make the color index scale.
  col.inds <- sapply(t(col.inds), function(col.ind) {
    min(max(1, col.ind), 100)
  })
  col.inds[which(is.na(col.inds))] <- 101
  # Make the color rectangle.
  rect(xlefts, ybottom, xright, ytop, col=color.dist[col.inds], lwd=0.20,
       xpd=NA, border="white")
  # Make the key.
  if (key) {
    y_div <- ymax/100
    kpl <- 17 # Left-hand position of the key
    kw <- 1.5                # Width of the key
    rect(xleft=kpl, ybottom=seq(100)*y_div - y_div + 0.5, xright=kpl + kw,
         ytop=seq(100)*y_div - y_div + 0.5 + y_div,
         col=color.dist[1:100], xpd=NA, border=NA)
    # Generate the axis for the legend and the label
    bottom <- ymin + y_div/2
    top <- ymax + 0.5 - y_div/2
    if (kdrel) {
      labels <- c(0.3, 1, 3, 10, 30, 100)
    } else {
      if (site %in% c("8mer", "7mer-m8", "7mer-A1", "6mer")) {
        labels <- c(3, 1, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001, 0.0003)
      } else {
        labels <- c(3, 1, 0.3, 0.1, 0.03, 0.01, 0.003)
      }
    }
    pos_labels <- log(labels)
    centered_labels <- pos_labels - pos_labels[1]
    norm_labels <- centered_labels/(centered_labels[length(centered_labels)])
    # labels_span <- max(log(labels)) - min(log(labels))
    height_span <- top - bottom
    pos_labels <- norm_labels*height_span
    axis(4, at=pos_labels, labels=labels, lwd=par()$lwd, pos=kpl + kw, xpd=NA)
    bquote(italic(r)^2~.(" with flanking dinucleotide")~italic(K)[D]*.("s"))
    if (kdrel) {
      text(x=kpl + kw + 3, y=10,
           labels=bquote(italic(K)[D]*.(" fold change")),
           srt=270, xpd=NA)    
    } else {
      text(x=kpl + kw + 3.5, y=10,
           labels=bquote(.("Relative")~italic(K)[D]*.("s")),
           srt=270, xpd=NA)    
    }
  }
  text(xmin, ymax + 1, labels=site, xpd=NA, adj=c(0, 0))
  # Finish the plot.
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

PlotLet7aProgrammedLibraryReps <- function(
  n_constant=0, lambda=0.01, fulloptim=FALSE, rand=FALSE, doublerand=FALSE,
  combined_rand=TRUE, allmm=FALSE, suppcomp=TRUE, collapsemm=FALSE, alpha=0.5,
  prog_mm=FALSE, kmer_min=4, kmer_max=11, plotseed=TRUE, height=5, width=5,
  xpos=20, ypos=20, pdf.plot=FALSE

) {
  # Get the two sets of Kds
  if (suppcomp) {
    sitelist <- "programmed_suppcomp"
    grep_base_string <- "Comp"
    grep_search_string <- "(^.*)\\|(.*)\\|(%s_Kd)"
    bottom_label <- "Counts combined over\nseed-mismatch type"
    n_constant_1 <- 0
  } else {
    sitelist <- "programmed"
    grep_base_string <- "8mer-mm[ACTG][2-7]"
    grep_search_string <- "(^.*)\\|(.*)\\|(%s_Kd)"
    bottom_label <- "Counts"
    n_constant_1 <- 0
  }
  if (doublerand) {
    experiment_1 <- "equilibrium_nb"
    sitelist_1 <- "bipartite_random"
    if (suppcomp) {
      sitelist_1 <- sprintf("%s_suppcomp", sitelist_1)
    }
    n_constant_2 <- 5
  } else {
    experiment_1 <- "equil_c2_nb"
    sitelist_1 <- sitelist
    n_constant_2 <- n_constant    
  }
  kds_collapsed_1 <- SubfunctionCall(EquilPars, mirna="let-7a-21nt",
                                     experiment=experiment_1,
                                     n_constant=n_constant_1,
                                     sitelist=sitelist_1)
  if (rand) {
    mirna_2 <- "let-7a"
    experiment_2 <- "equilibrium"
    sitelist_2 <- "bipartite_random"

    grep_search_string <- "(^.*)\\|(.*)\\|(%s_Kd)|(^.*)(_Kd)"
    grep_replace_string <- "\\1\\4"

    if (suppcomp) {
      sitelist_2 <- sprintf("%s_suppcomp", sitelist_2)
    }
    n_constant_2 <- 5

  } else {
    mirna_2 <- "let-7a-21nt"
    experiment_2 <- "equil_c_nb"
    sitelist_2 <- sitelist
    n_constant_2 <- n_constant
    grep_replace_string <- "\\1"
    # combined <- TRUE
  }
  kds_collapsed_2 <- SubfunctionCall(EquilPars, mirna=mirna_2,
                                     experiment=experiment_2,
                                     n_constant=n_constant_2,
                                     sitelist=sitelist_2,
                                     combined=combined_rand)
  # Get the names of the sites to be plotted; specifically, all instances of a
  # seed site or a thrP site in the random region.
  if (rand) {
    # Need to 1.) Remove the programmed versions of the seed sites from the
    # programmed libraries. 2.) Need to log-average the random-region seed sites
    # in the programmed region. 3.) Relabel these as the
    print(head(rownames(kds_collapsed_1)))
    # Converts the "NA|NA|8mer_Kd" rownames to "8mer_Kd".
    if (!suppcomp) {
      rownames(kds_collapsed_1) <- gsub("^NA\\|NA\\|(.*)$", replacement="\\1",
                                        rownames(kds_collapsed_1), perl=TRUE)
    }
    seed_sites_in_programmed <- grep("\\|", rownames(kds_collapsed_1),
                                     invert=TRUE, perl=TRUE, value=TRUE)      
    seed_sites_in_programmed <- seed_sites_in_programmed[1:(length(seed_sites_in_programmed) - 2)]
    if (suppcomp) {
      regex_target_base <- "%s\\|.*\\|Comp_Kd"
    } else {
      regex_target_base <- "%s\\|.*\\|8mer-mm[ACGT][2-7]_Kd"
    }
    if (prog_mm) {
      mm_lim <- 24
    } else {
      mm_lim <- 24
    }
    kd_mm <- GeoMean(kds_collapsed_1[7:24, 2])

    for (site in seed_sites_in_programmed[1:mm_lim]) {
      site_trim <- unlist(strsplit(site, split="_"))[1]
      print(site_trim)
      regex_target <- sprintf(regex_target_base, site_trim)
      # This is the rows being averaged to get the single seed mismatch
      # Kd value.
      rows_ave <- grep(regex_target, rownames(kds_collapsed_1), perl=TRUE)
      if (site == "8mer-mmG4_Kd") {
        print(head(rownames(kds_collapsed_1)[rows_ave]))
        print(tail(rownames(kds_collapsed_1)[rows_ave]))
      }
      check_sites <<- rownames(kds_collapsed_1)[rows_ave]
      # This is supposed to be a single row.
      row_replace <- grep(sprintf("^%s$", site), rownames(kds_collapsed_1), perl=TRUE)
      replace_sites <<- rownames(kds_collapsed_1)[row_replace]
      kd_obs <- GeoMean(kds_collapsed_1[rows_ave, 2])
      print(kd_obs)
      print(kd_mm)
      print(kd_obs*kd_mm/(kd_mm - kd_obs))
      # if (substring(site, 1, 7) == "8mer-mm" & prog_mm) {
      #   base_kd <- kds_collapsed_1[site, 2]*mm_kd/(kds_collapsed_1[site, 2] + mm_kd)
      # }
      kds_collapsed_1[row_replace, 2] <- kd_obs*kd_mm/(kd_mm - kd_obs)
    }
    # kds_collapsed_1 <- kds_collapsed_1[grep("\\|", rownames(kds_collapsed_1),
    #                                         invert=TRUE, perl=TRUE), ]



  }
  sites <- intersect(rownames(kds_collapsed_1), rownames(kds_collapsed_2))
  sites <<- sites
  if (rand) {
    sites <- c(sites[1:24],
               grep("&", grep(sprintf("%s_Kd", grep_base_string), sites,
                              value=TRUE), invert=TRUE, value=TRUE))

  } else {
    sites <- grep("&", grep(sprintf("%s_Kd", grep_base_string), sites,
                            value=TRUE), invert=TRUE, value=TRUE)
  }
  # Get the positioning of the 3p site with respect to the programmed site.
  # sites_nums <- gsub("^.*\\|(.*)\\|.*$", replacement="\\1", sites)
  # Take only those sites that have a positioning at at least position 9.
  # sites <- sites[which(as.integer(sites_nums) >= 9)]

  site_names <- gsub(sprintf(grep_search_string, grep_base_string),
                     replacement=grep_replace_string, sites, perl=TRUE)
  site_names_unique <- unique(site_names)
  site_names_unique <- grep("Kd", site_names_unique, invert=TRUE, value=TRUE)

  seed_site_names <- site_names[which(site_names %in% site_names_unique[1:24])]
  thrp_site_names <- site_names[which(!(site_names %in% site_names_unique[1:24]))]
  thrp_site_lens <- as.integer(gsub("^(.*)mer.*", replacement="\\1", thrp_site_names))

  thrp_sites_use <- thrp_site_names[which(thrp_site_lens >= kmer_min & thrp_site_lens <= kmer_max)]
  sites_unique_use <- unique(c(seed_site_names, thrp_sites_use))
  sites_unique_use <- grep("Kd", sites_unique_use, invert=TRUE, value=TRUE)
  thrp_site_lens <- as.integer(gsub("^(.*)mer.*", replacement="\\1", thrp_sites_use))
  if (fulloptim) {
    kds_1 <- SubfunctionCall(EquilPars, mirna="let-7a-21nt",
                                       experiment="equil_c2_nb",
                                       sitelist="programmed")
    kds_2 <- SubfunctionCall(EquilPars, mirna="let-7a-21nt",
                                       experiment="equil_c_nb",
                                       sitelist="programmed")
    sites_use <- intersect(rownames(kds_1), rownames(kds_2))
    sites_use <- grep("&", sites_use, invert=TRUE, value=TRUE)
    sites_use <- grep("8mer-mm[ACTG][2-7]_Kd$", sites_use, perl=TRUE, value=TRUE)
    # Get the names of the sites within the random region.
    site_names <- gsub("(^.*)\\|(.*)\\|(Comp_Kd)", replacement="\\1",
                       sites_use, perl=TRUE)
    kds_1 <- kds_1[sites_use, ]
    kds_2 <- kds_2[sites_use, ]
  } else if (allmm) {
    kds_1 <- do.call("rbind", lapply(
      sites_unique_use, GetPositionalProgKds, mirna="let-7a-21nt",
      experiment="equil_c2_nb", n_constant=n_constant, lambda=lambda,
      suppcomp=suppcomp
    ))
    kds_2 <- do.call("rbind", lapply(
      sites_unique_use, GetPositionalProgKds, mirna="let-7a-21nt",
      experiment="equil_c_nb", n_constant=n_constant, lambda=lambda,
      suppcomp=suppcomp
    ))
    # Get the intersection of the two matrices of kds.
    sites_use <- intersect(rownames(kds_1), rownames(kds_2))
    # Get the names of the sites within the random region.
    site_names <- gsub("(^.*)\\|(.*)\\|(Comp_Kd)", replacement="\\1",
                       sites_use, perl=TRUE)
    # Subset the two kd matrices such that each only has shared rows.
    kds_1 <- kds_1[sites_use, ]
    kds_2 <- kds_2[sites_use, ]

  } else {
    sites_use <- sites[which(site_names %in% sites_unique_use)]
    kds_1 <- kds_collapsed_1[sites_use, ]
    kds_2 <- kds_collapsed_2[sites_use, ]
    if (rand) {
      # Swap kds_1 and kds_2 in order to have the random kds on the x-axis.
      kds_temp <- kds_1
      kds_1 <- kds_2
      kds_2 <- kds_temp
    }
  }
  sites_base <- sapply(sites_use, function(site_use) {
    unlist(strsplit(site_use, split="\\|"))[1]
  })
  if (rand) {
    sites_base <- sapply(sites_base, function(site_base) {
      unlist(strsplit(site_base, split="_Kd"))[1]
    })
  }
  inds_seed <- which(sites_base %in% c(kSeedSites, GetAll8merMmSites("let-7a")))
  inds_non_seed <- setdiff(1:length(sites_base), inds_seed)
  lens_colors <- sapply(sites_use[inds_non_seed], function(site_use) {
    unlist(strsplit(site_use, split="mer"))[1]  
  })
  print(lens_colors)
  cols_thrp <- ConvertRColortoRGB(kThrPLengthCols, alpha=alpha)
  names(cols_thrp) <- 4:11
  cols <- rep("black", length(sites_use))
  cols[inds_seed] <- kSiteColors[sites_base[inds_seed]]
  cols[inds_non_seed] <- cols_thrp[as.character(lens_colors)]
  print(cols)
  site_cols <- cols
  # Reorder the sites to make seed sites first and non-seed sites second.
  sites_use <- sites_use[c(inds_seed, inds_non_seed)]

  SubfunctionCall(FigureSaveFile2)
  if (rand) {
    xmin <- 0.0001
  } else {
    xmin <- 0.001
  }
  ymin <- xmin
  xmax <- 1
  ymax <- xmax
  BlankPlot(log="xy")
  if (rand & suppcomp) {
    CostFunction <- function(pars, x, y) {
      m <- pars[1]
      b <- pars[2]
      y_p <- log(m*x + b)
      res <- (log(y) - y_p)^2
      return(sum(res))
    }
    pars_init <- c(1, 1)
    x_s <- 1/kds_1[sites_use[inds_seed], 2]
    y_s <- 1/kds_2[sites_use[inds_seed], 2]
    x_ns <- 1/kds_1[sites_use[inds_non_seed], 2]
    y_ns <- 1/kds_2[sites_use[inds_non_seed], 2]

    x_sns <- 1/kds_1[sites_use, 2]
    y_sns <- 1/kds_2[sites_use, 2]



    print(x_s)
    print(y_s)
    # dev.new()
    # plot(x_s, y_s)
    opt <- optim(pars_init, CostFunction, x=x_s, y=y_s)
    pars <- opt$par
    m_s <- pars[1]
    b_s <- pars[2]

    opt <- optim(pars_init, CostFunction, x=x_ns, y=y_ns)
    pars <- opt$par
    m_ns <- pars[1]
    b_ns <- pars[2]

    opt <- optim(pars_init, CostFunction, x=x_sns, y=y_sns)
    pars <- opt$par
    m_sns <- pars[1]
    b_sns <- pars[2]

  }

  segments(xmin, ymin, x1=xmax, y1=ymax, col="gray")
  if (rand & suppcomp) {
    x_lines <- 10^seq(-4, 0, length.out=100)
    y_s_lines <- 1/(m_s/x_lines + b_s)
    y_ns_lines <- 1/(m_ns/x_lines + b_ns)
    y_sns_lines <- 1/(m_sns/x_lines + b_sns)

    lines(x_lines, y_s_lines, col="purple")
    lines(x_lines, y_ns_lines, col="forestgreen")
    lines(x_lines, y_sns_lines, col="black")
  }
  # Add the points to the plot.

  # points(kds_1[sites_use[inds_seed], 2], kds_2[sites_use[inds_seed], 2],
  #        col=site_cols[inds_seed], pch=20, xpd=NA)

  kds_1_global <<- kds_1[sites_use, ]
  kds_2_global <<- kds_2[sites_use, ]

  points(kds_1[sites_use[(length(inds_seed) + 1):length(sites_use)], 2],
         kds_2[sites_use[(length(inds_seed) + 1):length(sites_use)], 2],
         col=site_cols[inds_non_seed], pch=20, xpd=NA)

  if (plotseed) {
    points(kds_1[sites_use[1:length(inds_seed)], 2],
           kds_2[sites_use[1:length(inds_seed)], 2],
           col=site_cols[inds_seed], pch=20, xpd=NA)    
  }


  # if (ident & class(pdf.plot) != "character") {
  #   identify(kds_1[sites_use, 2],
  #            kds_2[sites_use, 2], labels=sites_use[c(inds_non_seed, inds_seed)])    
  # }
  xy <- GetPlotFractionalCoords(0.1, 0.95, log="xy")
  text(xy[1], xy[2], labels="let-7a", adj=0, xpd=NA)
  xy <- GetPlotFractionalCoords(0.1, 0.9, log="xy")
  if (plotseed) {
    sites_cor <- sites_use
  } else {
    sites_cor <- sites_use[inds_non_seed]
  }
  AddCorrelationToPlot(x=log(kds_1[sites_cor, 2]), y=log(kds_2[sites_cor, 2]),
                       xpos=xy[1], ypos=xy[2], rsquared=TRUE)
  xy <- GetPlotFractionalCoords(0.1, 0.85, log="xy")
  text(xy[1], xy[2], labels=sprintf("n = %s", length(sites_cor)), adj=0, xpd=NA)

  xy <- GetPlotFractionalCoords(1, 0.05, log="xy")

  text(xy[1], xy[2], labels=bottom_label, adj=c(1, 0), xpd=NA)

  if (rand) {
    x_string <- "Relative Kd; random-sequence library"
    y_string <- "Relative Kd; programmed library"
  } else {
    x_string <- "Relative Kd; replicate 1"
    y_string <- "Relative Kd; replicate 2"
  }

  xy <- GetPlotFractionalCoords(0.05, 0.1, log="xy")
  text(xy[1], xy[2], sprintf("%s-%s nt of pairing", kmer_min, kmer_max),
       adj=c(0, 0))
  # Text saying if the programmed
  if (rand) {
    xy <- GetPlotFractionalCoords(0.05, 0.05, log="xy")
    # if (prog_mm) {
    #   text(xy[1], xy[2], sprintf("Programmed 8mer-mm Kds"), adj=c(0, 0))
    # } else {
    #   text(xy[1], xy[2], sprintf("Random 8mer-mm Kds"), adj=c(0, 0))    
    # }    
  }

  AddLogAxis(1, x_string)
  AddLogAxis(2, y_string)
  # Finish the plot.
  if (class(pdf.plot) == "character") {
    dev.off()
  }

}

# PlotCanonicalSiteECDF <- function(
#   mirna="let-7a-21nt", experiment="equil_c2_nb", n_constant=0,
#   sitelist="programmed_suppcomp", collapsemm=FALSE, kd_fc=FALSE, height=5,
#   width=5, xpos=20, ypos=20, pdf.plot=FALSE
# ) {
#   kds <- SubfunctionCall(EquilPars, sitelist="programmed_suppcomp")
#   SubfunctionCall(FigureSaveFile2)
#   xmin <- 1e-4
#   xmax <- 1
#   ymin <- 0
#   ymax <- 1
#   BlankPlot(log="x")
#   can_sites <- c("8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-m8", "6mer-A1")
#   cols <- kSiteColors[can_sites]
#   names(cols) <- can_sites
#   if (kd_fc) {
#     kd_ref <- kds["Comp_Kd", 2]
#   } else {
#     kd_ref <- 1
#     grep_str <- "^8mer-mm[ACTG][2-7]_Kd$"
#     inds <- grep(grep_str, rownames(kds), perl=TRUE, value=TRUE)
#     df_use <- kds[inds, ]/kd_ref
#     x <- 10^seq(-3, 1, length.out=100)
#     y <- ecdf(df_use[, 2])(x)
#     lines(x, y, col="black", lwd=1, xpd=NA)
#   }
#   # Apply function over the canonical site names, that plots each of the lines. 
#   sapply(can_sites, function(can_site) {
#     grep_str <- sprintf("^%s\\|.*\\|Comp_Kd", can_site)
#     inds <- grep(grep_str, rownames(kds), perl=TRUE, value=TRUE)
#     df_use <- kds[inds, ]
#     x <- 10^seq(-3, 1, length.out=100)
#     y <- ecdf(df_use[, 2])(x)
#     lines(x, y, col=cols[can_site], lwd=1, xpd=NA)
#   })
#   # Add the x- and y-axes
#   AddLogAxis(1, label="Relative Kd")
#   AddLinearAxis(2, tick.space=0.05, label.space=0.2,
#                 label="CDF (%)", percent=TRUE)
#   # Add the legend to the plot.
#   xy <- GetPlotFractionalCoords(0, 1.05, log="x")
#   Legend(xy, legend=c("Seed mismatch only", can_sites),
#          col=c("black", cols), ncol=1)
#   # Finish the plot.
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }



# MismatchBarPlot <- function(
#   mirna, experiment, site, n_constant=3, kdrel=TRUE, offset_plot="0",
#   key=FALSE, average=FALSE, height=5, width=5, mirna_label=TRUE, xpos=20,
#   ypos=20, pdf.plot=FALSE
# ) {
#   # Load the collapsed data (for the mismatch-only).
#   kds_site <- SubfunctionCall(EquilPars, sitelist="programmed")
#   print(head(kds_site))
#   mm_8mer_sites <- GetAll8merMmSites(mirna)
#   ref_sites <- sprintf("NA|NA|%s_Kd", mm_8mer_sites)
#   ref_kd <- GeoMean(kds_site[ref_sites, 2])
#   print(ref_kd)
#   # Get the starting position of the site. 
#   start_pos <- as.integer(gsub("^.*mer-m(.*)\\..*$", replacement="\\1", site,
#                                perl=TRUE))
#   # Pre-allocate the matrix
#   out_matrix <- matrix(NA, nrow=18, ncol=length(-6:20), dimnames=list(mm_8mer_sites, -6:20))
#   out <- sapply(mm_8mer_sites, function(mm_site) {
#       # get the indeces of the rows that have the site in the function argument
#       # and the mismatch site of the function.
#       inds <- grep(sprintf("^%s\\|.*\\|%s_Kd$", site, mm_site),
#                    rownames(kds_site), perl=TRUE)
#       # Get the names of the sites.
#       names <- rownames(kds_site)[inds]
#       print(names)
#       offsets <- as.integer(gsub(sprintf("^%s\\|(.*)\\|%s_Kd$", site, mm_site),
#                                  replacement="\\1", names, perl=TRUE)) - start_pos
#       # Assign the kds to the appropriate row and column within the output matrix.
#       out_matrix[mm_site, as.character(offsets)] <<- kds_site[inds, 2]
#   })
#   out_matrix  <- out_matrix/c(kds_site[ref_sites, 2])
#   out_matrix <- out_matrix[, which(colnames(out_matrix) %in% as.character(-4:16))]
#   out_matrix <<- out_matrix
#   cols <- rainbow(length(-4:16), start=0.6, end=1)
#   names(cols) <- as.character(-4:16)
#   R_mat <- log10(1/out_matrix)
#   SubfunctionCall(FigureSaveFile2)

#   # par(mar=c(mar1, mar2, mar3, mar4))
#   xmin <- 0
#   xmax <- ncol(R_mat)
#   ymin <- 0
#   ymax <- nrow(R_mat) 
#   BlankPlot()
#   print(R_mat)
#   # Assign the positions of the corners of the boxes.
#   xlefts <- rep(seq(0, ncol(R_mat) - 1), nrow(R_mat))
#   xright <- xlefts + 1
#   ybottom <- rep(rev(seq(ymax - 1, ymin)), each=ncol(R_mat))
#   ytop <- ybottom + 1
#   # Make the x-axis
#   AddLinearAxis(1, 1, 1, label="5'-paired nt",
#                 alt_lab=colnames(R_mat),
#                 alt_lab_pos=xmin:(xmax - 1) + 0.5,
#                 alt_tick_pos=TRUE)
#   AddLinearAxis(2, 1, 1, label="3'-paired nt",
#                 label_pos_ticks=TRUE,
#                 alt_lab=gsub("^8mer-mm(.*)$", replacement="\\1",
#                              rownames(R_mat), perl=TRUE),
#                 alt_lab_y_dist=0.01,
#                 alt_lab_pos=ymin:(ymax - 1) + 0.5,
#                 alt_tick_pos=TRUE)
#   # Add the label for the seed if not plotting relative to seed kds.
#   col.inds <- round(t((t(R_mat))/apply(R_mat, 2, max, na.rm=TRUE))*99 + 1)
#   x_lab_pos <- -0.75
#   pos_5p <- as.integer(rep(rownames(R_mat), each=ncol(R_mat)))
#   pos_3p <- as.integer(rep(colnames(R_mat), nrow(R_mat)))
#   impossible_cols <- which(pos_5p > pos_3p)

#   start_col <- 0
#   r_col <- 0.4
#   # Make distribution of 100 colors and a 101th gray color for NA indeces.
#   color.dist <- c(rev(cubeHelix(100, start=start_col, r=r_col, hue=0.8)),
#                   "gray90", "white")
#   # Make the color index scale.
#   col.inds <- sapply(t(col.inds), function(col.ind) {
#     min(max(1, col.ind), 100)
#   })
#   col.inds[which(is.na(col.inds))] <- 101
#   col.inds[impossible_cols] <- 102
#   # Make the color rectangle.
#   rect(xlefts, ybottom, xright, ytop, col=color.dist[col.inds], lwd=0.4, 
#        xpd=NA, border="white")
#   xy <- GetPlotFractionalCoords(fx=0.05, fy=1.05)
#   text(xy[1], xy[2], labels=site, adj=c(0, 0), xpd=NA)
#   # Finish the plot.
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }



PlotEmpiricalMismatchCoefficients <- function(
  mirna, experiment, n_constant=3, offsetmin=-4, offsetmax=16, len_min=4,
  len_max=11, pos_3p_max=23, pos_3p_min=9, offset_pick=NA, sitelist="programmed",
  corrected_kds=TRUE, supp_base=FALSE, add_comp_sites=FALSE,
  collapsemm_addcomp=FALSE, kd_fc=TRUE, xpos=20, ypos=20, height=3.5,
  width=4, pdf.plot=FALSE, rand_base=FALSE, percentile=0.25
) {
  R_mat <<- SubfunctionCall(GetEmpiricalBestSeedMismatchCoefficients)
  SubfunctionCall(FigureSaveFile2)
  xmin <- 0
  xmax <- ncol(R_mat) + 2
  ymin <- 0
  ymax <- nrow(R_mat)
  if (supp_base) {
    par(mar=c(2, 4, 0.5, 1)) 
  } else {
    par(mar=c(2, 2, 0.5, 1)) 
  }
  BlankPlot()
  xlefts <- rep(seq(0, ncol(R_mat) - 1), each=nrow(R_mat))
  xright <- xlefts + 1
  ybottom <- rep(rev(seq(ymax - 1, ymin)), ncol(R_mat))
  ytop <- ybottom + 1

  # Make the x-axis
  text(x=-1, y=ymin:(ymax - 1) + 0.5,
       labels=gsub("(8mer-mm)(.*)", rownames(R_mat),
                   replacement="\\2",
                   perl=TRUE),
       xpd=NA, adj=c(1, 0.5))
  # # Add the label for the seed if not plotting relative to seed kds.
  col_min <- -1
  col_max <- 1    

  col.inds <- round((R_mat - col_min)/(col_max - col_min)*99 + 1)
  col.inds[which(col.inds <= 0)] <- 1
  col.inds[which(col.inds > 100)] <- 100
  print(col.inds) 
  # Make distribution of 100 colors and a 101th gray color for NA indeces.
  color.dist <- rev(colorspace::diverging_hcl(100, h=c(180, 50), c=80,
                                              l=c(20, 95), power=c(0.7, 1.3)))


  rect(xlefts, ybottom, xright, ytop, col=color.dist[col.inds], lwd=0.4, 
       xpd=NA, border="white")
  rect(xmin, ymin, round(percentile*xmax), ymax, col=NULL, border="black", xpd=NA)

  R_mat_ave <<- rowMeans(R_mat[, 1:(round(percentile*ncol(R_mat)))], na.rm=TRUE)
  xlefts <- ncol(R_mat) + 1
  xright <- xlefts + 1
  ybottom <- rev(seq(ymax - 1, ymin))
  ytop <- ybottom + 1
  col.inds <- round((R_mat_ave - col_min)/(col_max - col_min)*99 + 1)

  rect(xlefts, ybottom, xright, ytop, col=color.dist[col.inds], lwd=0.4,
       xpd=NA, border="black")
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}


PlotAllRandomLibrarySeedCoefficients <- function(
  experiment="equilibrium", n_constant=3, sitelist="bipartite_random",
  corrected_kds=FALSE, combined=TRUE, buffer=FALSE, supp_base=TRUE, key=TRUE,
  len_min=4, len_max=11, pos_3p_min=9, pos_3p_max=23, offsetmin=-4,
  offsetmax=16, kd_fc=TRUE, add_comp_sites=TRUE, collapsemm_addcomp=FALSE,
  intercept=FALSE, empirical=FALSE, percentile=0.25, height=4, width=4, xpos=20,
  ypos=20, pdf.plot=FALSE
) {
  # Load the collapsed data (for the mismatch-only).
  intercepts_global <<- list(`let-7a`=c(), `miR-1`=c(), `miR-155`=c(),
                            `miR-124`=c(), `lsy-6`=c(), `miR-7-23nt`=c())
  data_all <- sapply(kMirnas, function(mirna) {
    if (mirna == "miR-1") {

      combined <- FALSE
      buffer <- TRUE
    } else if (mirna == "miR-7-23nt") {
      experiment <- "equilibrium2_nb"
      combined <- FALSE
    }
    if (empirical) {
      out_mat <- SubfunctionCall(GetEmpiricalBestSeedMismatchCoefficients)
      data <- rowMeans(out_mat[, 1:round(percentile*ncol(out_mat))])
    } else {
      model <- SubfunctionCall(FitPairingOffsetAndMismatchModel)
      if (intercept) {
        intercepts_global[[mirna]] <<- model$b
      }
      data <- model$mm    
    }
  })
  data_all <<- data_all
  # Make the plot and define the limits.
  SubfunctionCall(FigureSaveFile2)
  par(mar=c(3, 3, 2, 0.5))
  xmin <- 0
  xmax <- dim(data_all)[1]*(dim(data_all)[2] + 1)
  if (empirical) {
    ymin <- -1
    ymax <- 1
  } else {
    ymin <- 0
    ymax <- 2
  }
  BlankPlot()
  # col.inds <- round((data - ymin)/(ymax - ymin)*99 + 1)
  # col.inds[which(col.inds <= 0)] <- 1  
  # # Make distribution of 100 colors and a 101th gray color for NA indeces.
  # color.dist <- rev(colorspace::diverging_hcl(100, h=c(180, 50), c=80,
  #                                             l=c(20, 95), power=c(0.7, 1.3)))
  # if (supp_base) {

  x_shifts <<- rep(0:(dim(data_all)[1] - 1), each=dim(data_all)[2])
  x_base <<- 1:(dim(data_all)[1]*dim(data_all)[2])
  # xleft <- 1:(xmax) - 1 + rep(0:(dim(data_all)[1] - 1), each=dim(data_all)[2])
  xleft <<- x_base + x_shifts

  xright <- xleft + 1
  # } else {
  #   xleft <- 0:17 + 0.05 + rep(c(0, 1, 2, 3, 4, 5), each=3)
  #   xright <- xleft + 1
  # }
  # Get the left- and right-hand positions of the rectangles.
  # Make the barplot rectangles ################################################
  rect(xleft=xleft, ybottom=mean(c(ymin, ymax)), xright=xright, ytop=c(t(data_all)),
       col=rep(kMirnaColors, nrow(data_all)), border=NA)
  AddLinearAxis(2, 0.1, 0.2, label=expression(Delta*Delta*italic(G)~"scaling-factor"))

  # # Get the each of the nucleotide mismatches.
  # if (!supp_base) {
  #   nucs <- gsub("T", "U", gsub("8mer-mm(.)(.)", replacement="\\1", names(data),
  #                               perl=TRUE))
  #   pos <- unique(gsub("8mer-mm(.)(.)", replacement="\\2", names(data), perl=TRUE))

  #   # Get miRNA nucleotides 2-7.
  #   mirna_nucs <- unlist(strsplit(kMirnaSeqs[mirna], split=""))[2:7]
  #   pos_cols <- c()
  #   for (i in 0:5) {
  #     nuc_ind <- i + 1
  #     mirna_nuc <- mirna_nucs[nuc_ind]
  #     mirna_nucs_start <- i*3 + 1
  #     mirna_nucs_stop <- mirna_nucs_start + 2
  #     cols_i <- rep("black", 3)
  #     if (mirna_nuc == "G") {
  #       ind_wobble <- which(nucs[mirna_nucs_start:mirna_nucs_stop] == "U")
  #       cols_i[ind_wobble] <- "blue"
  #     } else if (mirna_nuc == "U") {
  #       ind_wobble <- which(nucs[mirna_nucs_start:mirna_nucs_stop] == "G")
  #       cols_i[ind_wobble] <- "red"
  #     }
  #     pos_cols <- c(pos_cols, cols_i)
  #   }
  #   xy <- GetPlotFractionalCoords(0.05, -0.05)
  #   text(x=(xleft + xright)/2, y=xy[2], labels=nucs, col=pos_cols, xpd=NA)
  #   xy <- GetPlotFractionalCoords(0.05, -0.1)
  #   text(x=0:5*4 + 1.5, y=xy[2], labels=pos, xpd=NA)
  # } else {
  #   xy <- GetPlotFractionalCoords(0.05, 0)
  #   text(x=0:5*4 + 1.5, y=xy[2], labels=names(data), srt=45, xpd=NA, adj=c(1, 1))
  # }

  # if (mirna_label) {
  #   if (mirna == "let-7a-21nt") {
  #     mirna_txt <- "let-7a"
  #   } else if (mirna == "let-7a_minus1") {
  #     mirna_txt <- "let-7a(-1)"
  #   } else if (mirna == "let-7a_plus1") {
  #     mirna_txt <- "let-7a(+1)"
  #   } else if (mirna == "let-7a_miR-155") {
  #     mirna_txt <- "let-7a-miR-155"
  #   } else if (mirna == "miR-155_let-7a") {
  #     mirna_txt <- "miR-155-let-7a"
  #   } else {
  #     mirna_txt <- mirna
  #   }
  #   # x <- GetPlotFractionalCoords(0.0125, 0.5)[1]
  #   # mtext(text=mirna_txt, side=3, line=0, at=x, adj=c(0, 0), cex=par("cex"),
  #   #       xpd=NA)
  #   xy <- GetPlotFractionalCoords(0.025, 1.1)
  #   text(xy[1], xy[2], labels=sprintf(mirna_txt), adj=c(0, 1), xpd=NA)  

  # }
  ## Add label explaining that these are the model mismatch coefficients.
  # x <- GetPlotFractionalCoords(0.0125, 0.5)[1]
  xy <- GetPlotFractionalCoords(0.025, 1.025)
  if (empirical) {
    text(xy[1], xy[2], labels="Empirical variation", adj=c(0, 1), xpd=NA)  
  } else {
    text(xy[1], xy[2], labels="Mismatch coefficients", adj=c(0, 1), xpd=NA)  
  }

  # mtext(text="Mismatch coefficients", side=3, line=-1, at=x, adj=c(0, 0),
  #       cex=par("cex"), xpd=NA)
  # Finish the plot.
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

PlotMismatchCoefficientsAgainstKds <- function(
  mirna, experiment, n_constant=3, sitelist="programmed", corrected_kds=TRUE,
  globalmodel=TRUE, combined=TRUE, buffer=FALSE, average_all=FALSE,
  mirna_label=TRUE, len_min=4, len_max=11, pos_3p_min=9, pos_3p_max=23,
  offsetmin=-4, offsetmax=16, supp_base=FALSE, empirical=FALSE, kd_fc=TRUE,
  percentile=0.25, height=3.5, width=3.5, xpos=20, ypos=20, pdf.plot=FALSE
) {
  # Load the collapsed data (for the mismatch-only).
  if (experiment == "all_prog") {
    if (mirna == "let-7a") {
      mirnas <- c("let-7a-21nt", "let-7a-21nt", "let-7a-21nt",
                     "let-7a_plus1", "let-7a_minus1", "let-7a_miR-155")
      exps <- c("equil_c_nb", "equil_c2_nb", "equil_sc_nb",
                    "equil_c_nb", "equil_c_nb", "equil_c_nb")
    } else if (mirna == "miR-1") {
      mirnas <- c("miR-1", "miR-1")
      exps <- c("equil_c_nb", "equil_sc_nb")
    } else if (mirna == "miR-155") {
      mirnas <- c("miR-155", "miR-155_let-7a")
      exps <- c("equil_sc_nb", "equil_c_nb")
    }
    model_coefs_all <- apply(cbind(mirnas, exps), 1, function(row) {
      SubfunctionCall(FitPairingOffsetAndMismatchModel, mirna=row[1],
                      experiment=row[2])$mm
    })
    model_coefs_all <<- model_coefs_all
    mm_coefs <- rowMeans(model_coefs_all)
    names(mm_coefs) <- rownames(model_coefs_all)
    kds_all <- do.call("cbind", apply(cbind(mirnas, exps), 1, function(row) {
      if (corrected_kds) {
        kds <- SubfunctionCall(
          ApplyKdCorrection, mirna=row[1], experiment=row[2],
          rand_n_constant=n_constant, prog_n_constant=n_constant,
          prog_sitelist=sitelist, rand_sitelist="bipartite_random"
        )
      } else {
        kds <- SubfunctionCall(EquilPars, mirna=row[1], experiment=row[2]) 
      }
      kds[sprintf("%s_Kd", names(mm_coefs)), 2, drop=FALSE]
    }))
    kds_all <<- kds_all
    kds_programmed <- exp(rowMeans(log(kds_all)))
  } else {
    if (empirical) {
      out_mat <- SubfunctionCall(GetEmpiricalBestSeedMismatchCoefficients)
      out_mat_global <<- out_mat
      data <- rowMeans(out_mat[, 1:round(percentile*ncol(out_mat))])
      mm_coefs <- data - mean(data) + 1
      mm_coefs_check <<- mm_coefs

    } else {
      if (globalmodel) {
        model <<- SubfunctionCall(FitPairingOffsetAndMismatchModel)
      }    
      mm_coefs <- model$mm
    }
    mm_coefs <<- mm_coefs
    if (corrected_kds) {
      kds <- SubfunctionCall(
        ApplyKdCorrection, rand_n_constant=n_constant,
        prog_n_constant=n_constant, prog_sitelist=sitelist,
        rand_sitelist="bipartite_random"
      )
    } else {
      if (mirna == "miR-1" & experiment == "equilibrium") {
        buffer <- TRUE
        combined <- FALSE
      } else if (mirna == "miR-7-23nt" & experiment == "equilibrium2_nb") {
        buffer <- FALSE
        combined <- FALSE
      } else {
        buffer <- FALSE
        combined <- TRUE
      }
      kds <- SubfunctionCall(EquilPars)      
    }
    print(head(kds, n=24))
    kds_programmed <- kds[sprintf("%s_Kd", names(mm_coefs)), 2]
  }
  # Make the plot and define the limits.
  SubfunctionCall(FigureSaveFile2)
  par(mar=c(3, 3, 2, 2))
  if (supp_base) {
    xmin <- 1e-4
  } else if (corrected_kds | (experiment == "equilibrium") | (experiment == "equilibrium2_nb")) {
    xmin <-1e-2
  } else {
    xmin <- 1e-1
  }
  xmax <- 1      
  if ((mirna == "miR-7-23nt" | empirical) & !supp_base) {
    ymin <- 0
    ymax <- 2
  } else {
    ymin <- 0.4
    ymax <- 1.6
  }
  BlankPlot(log="x")

  col.inds <- floor((mm_coefs - ymin)/(ymax - ymin)*99 + 1)
  print(col.inds) 
  col.inds[which(col.inds <= 0)] <- 1
  # Make distribution of 100 colors and a 101th gray color for NA indeces.
  color.dist <- rev(colorspace::diverging_hcl(100, h=c(180, 50), c=80,
                                              l=c(20, 95), power=c(0.7, 1.3)))
  # Get the left- and right-hand positions of the rectangles.
  # Make the barplot rectangles ################################################
  kds_programmed <<- kds_programmed
  mm_coefs <<- mm_coefs
  Points(kds_programmed, mm_coefs, col=color.dist[col.inds])
  xy <- GetPlotFractionalCoords(0.9, 0.1, log="x")
  print(xy)
  AddCorrelationToPlot(xy[1], xy[2], x=log(kds_programmed), y=mm_coefs,
                       rsquared=TRUE, adj=1)
  AddLogAxis(1, label="Relative Kd")
  AddLinearAxis(2, 0.1, 0.2,
    label=expression(Delta*Delta*italic(G)~"scaling-factor"))
  # text(x=(xleft + xright)/2, y=0.6, labels=nucs, col=pos_cols, xpd=NA)
  # text(x=0:5*4 + 1.5, y=0.52, labels=pos, xpd=NA)

  if (mirna_label) {
    if (mirna == "let-7a-21nt") {
      mirna_txt <- "let-7a"
    } else if (mirna == "let-7a_minus1") {
      mirna_txt <- "let-7a(-1)"
    } else if (mirna == "let-7a_plus1") {
      mirna_txt <- "let-7a(+1)"
    } else if (mirna == "let-7a_miR-155") {
      mirna_txt <- "let-7a-miR-155"
    } else if (mirna == "miR-155_let-7a") {
      mirna_txt <- "miR-155-let-7a"
    } else {
      mirna_txt <- mirna
    }
    # x <- GetPlotFractionalCoords(0.0125, 0.5)[1]
    # mtext(text=mirna_txt, side=3, line=0, at=x, adj=c(0, 0), cex=par("cex"),
    #       xpd=NA)
    xy <- GetPlotFractionalCoords(0.025, 1.15)
    text(xy[1], xy[2], labels=sprintf(mirna_txt), adj=c(0, 1), xpd=NA)  

  }
  ## Add label explaining that these are the model mismatch coefficients.
  # x <- GetPlotFractionalCoords(0.0125, 0.5)[1]
  xy <- GetPlotFractionalCoords(0.025, 1.075)
  text(xy[1], xy[2], labels="Mismatch coefficients", adj=c(0, 1), xpd=NA)  

  # mtext(text="Mismatch coefficients", side=3, line=-1, at=x, adj=c(0, 0),
  #       cex=par("cex"), xpd=NA)
  # Finish the plot.
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}


MakeFigure7_trash <- function(collapsemm=TRUE) {

  PlotSimNucleotideContributionMatrixNew("let-7a-21nt", "equil_c2_nb", -4, globalmodel=TRUE, pdf.plot="8.Ai", mirna_label=TRUE, extralabel=TRUE)
  PlotSimNucleotideContributionMatrixNew("let-7a-21nt", "equil_c2_nb", -3, pdf.plot="8.Aii")
  PlotSimNucleotideContributionMatrixNew("let-7a-21nt", "equil_c2_nb", -2, pdf.plot="8.Aiii")
  PlotSimNucleotideContributionMatrixNew("let-7a-21nt", "equil_c2_nb", -1, pdf.plot="8.Aiv")
  PlotSimNucleotideContributionMatrixNew("let-7a-21nt", "equil_c2_nb",  0, pdf.plot="8.Av")
  PlotSimNucleotideContributionMatrixNew("let-7a-21nt", "equil_c2_nb",  1, pdf.plot="8.Avi")

  PlotSimNucleotideContributionMatrixNew("let-7a-21nt", "equil_c2_nb",  2, pdf.plot="8.Avii")
  PlotSimNucleotideContributionMatrixNew("let-7a-21nt", "equil_c2_nb",  3, pdf.plot="8.Aviii")
  PlotSimNucleotideContributionMatrixNew("let-7a-21nt", "equil_c2_nb",  4, pdf.plot="8.Aix")
  PlotSimNucleotideContributionMatrixNew("let-7a-21nt", "equil_c2_nb",  5, pdf.plot="8.Ax")
  PlotSimNucleotideContributionMatrixNew("let-7a-21nt", "equil_c2_nb",  6, pdf.plot="8.Axi")
  PlotSimNucleotideContributionMatrixNew("let-7a-21nt", "equil_c2_nb",  7, pdf.plot="8.Axii")

  PlotSimNucleotideContributionMatrixNew("let-7a-21nt", "equil_c2_nb",  8, pdf.plot="8.Axiii")
  PlotSimNucleotideContributionMatrixNew("let-7a-21nt", "equil_c2_nb",  9, pdf.plot="8.Axiv")
  PlotSimNucleotideContributionMatrixNew("let-7a-21nt", "equil_c2_nb", 10, pdf.plot="8.Axv")
  PlotSimNucleotideContributionMatrixNew("let-7a-21nt", "equil_c2_nb", 11, pdf.plot="8.Axvi")
  PlotSimNucleotideContributionMatrixNew("let-7a-21nt", "equil_c2_nb", 12, pdf.plot="8.Axvii")
  PlotSimNucleotideContributionMatrixNew("let-7a-21nt", "equil_c2_nb", 13, pdf.plot="8.Axviii")
}

MakeFigure8 <- function(collapsemm=TRUE) {
  PlotNucleotideContributionMatrix("miR-1", "equil_c_nb", -4, collapsemm=collapsemm, pdf.plot="7.Bi", mirna_label=TRUE, extralabel=TRUE)
  PlotNucleotideContributionMatrix("miR-1", "equil_c_nb", -3, collapsemm=collapsemm, pdf.plot="7.Bii")
  PlotNucleotideContributionMatrix("miR-1", "equil_c_nb", -2, collapsemm=collapsemm, pdf.plot="7.Biii")
  PlotNucleotideContributionMatrix("miR-1", "equil_c_nb", -1, collapsemm=collapsemm, pdf.plot="7.Biv")
  PlotNucleotideContributionMatrix("miR-1", "equil_c_nb",  0, collapsemm=collapsemm, pdf.plot="7.Bv")
  PlotNucleotideContributionMatrix("miR-1", "equil_c_nb",  1, collapsemm=collapsemm, pdf.plot="7.Bvi")

  PlotNucleotideContributionMatrix("miR-1", "equil_c_nb",  2, collapsemm=collapsemm, pdf.plot="7.Bvii")
  PlotNucleotideContributionMatrix("miR-1", "equil_c_nb",  3, collapsemm=collapsemm, pdf.plot="7.Bviii")
  PlotNucleotideContributionMatrix("miR-1", "equil_c_nb",  4, collapsemm=collapsemm, pdf.plot="7.Bix")
  PlotNucleotideContributionMatrix("miR-1", "equil_c_nb",  5, collapsemm=collapsemm, pdf.plot="7.Bx")
  PlotNucleotideContributionMatrix("miR-1", "equil_c_nb",  6, collapsemm=collapsemm, pdf.plot="7.Bxi")
  PlotNucleotideContributionMatrix("miR-1", "equil_c_nb",  7, collapsemm=collapsemm, pdf.plot="7.Bxii")

  PlotNucleotideContributionMatrix("miR-1", "equil_c_nb",  8, collapsemm=collapsemm, pdf.plot="7.Bxiii")
  PlotNucleotideContributionMatrix("miR-1", "equil_c_nb",  9, collapsemm=collapsemm, pdf.plot="7.Bxiv")
  PlotNucleotideContributionMatrix("miR-1", "equil_c_nb", 10, collapsemm=collapsemm, pdf.plot="7.Bxv")
  PlotNucleotideContributionMatrix("miR-1", "equil_c_nb", 11, collapsemm=collapsemm, pdf.plot="7.Bxvi")
  PlotNucleotideContributionMatrix("miR-1", "equil_c_nb", 12, collapsemm=collapsemm, pdf.plot="7.Bxvii")
  PlotNucleotideContributionMatrix("miR-1", "equil_c_nb", 13, collapsemm=collapsemm, pdf.plot="7.Bxviii")

  PlotSimNucleotideContributionMatrixNew("miR-1", "equil_c_nb", -4, globalmodel=TRUE, collapsemm=collapsemm, pdf.plot="8.Bi", mirna_label=TRUE, extralabel=TRUE)
  PlotSimNucleotideContributionMatrixNew("miR-1", "equil_c_nb", -3, collapsemm=collapsemm, pdf.plot="8.Bii")
  PlotSimNucleotideContributionMatrixNew("miR-1", "equil_c_nb", -2, collapsemm=collapsemm, pdf.plot="8.Biii")
  PlotSimNucleotideContributionMatrixNew("miR-1", "equil_c_nb", -1, collapsemm=collapsemm, pdf.plot="8.Biv")
  PlotSimNucleotideContributionMatrixNew("miR-1", "equil_c_nb",  0, collapsemm=collapsemm, pdf.plot="8.Bv")
  PlotSimNucleotideContributionMatrixNew("miR-1", "equil_c_nb",  1, collapsemm=collapsemm, pdf.plot="8.Bvi")

  PlotSimNucleotideContributionMatrixNew("miR-1", "equil_c_nb",  2, collapsemm=collapsemm, pdf.plot="8.Bvii")
  PlotSimNucleotideContributionMatrixNew("miR-1", "equil_c_nb",  3, collapsemm=collapsemm, pdf.plot="8.Bviii")
  PlotSimNucleotideContributionMatrixNew("miR-1", "equil_c_nb",  4, collapsemm=collapsemm, pdf.plot="8.Bix")
  PlotSimNucleotideContributionMatrixNew("miR-1", "equil_c_nb",  5, collapsemm=collapsemm, pdf.plot="8.Bx")
  PlotSimNucleotideContributionMatrixNew("miR-1", "equil_c_nb",  6, collapsemm=collapsemm, pdf.plot="8.Bxi")
  PlotSimNucleotideContributionMatrixNew("miR-1", "equil_c_nb",  7, collapsemm=collapsemm, pdf.plot="8.Bxii")

  PlotSimNucleotideContributionMatrixNew("miR-1", "equil_c_nb",  8, collapsemm=collapsemm, pdf.plot="8.Bxiii")
  PlotSimNucleotideContributionMatrixNew("miR-1", "equil_c_nb",  9, collapsemm=collapsemm, pdf.plot="8.Bxiv")
  PlotSimNucleotideContributionMatrixNew("miR-1", "equil_c_nb", 10, collapsemm=collapsemm, pdf.plot="8.Bxv")
  PlotSimNucleotideContributionMatrixNew("miR-1", "equil_c_nb", 11, collapsemm=collapsemm, pdf.plot="8.Bxvi")
  PlotSimNucleotideContributionMatrixNew("miR-1", "equil_c_nb", 12, collapsemm=collapsemm, pdf.plot="8.Bxvii")
  PlotSimNucleotideContributionMatrixNew("miR-1", "equil_c_nb", 13, collapsemm=collapsemm, pdf.plot="8.Bxviii")

}

MakeFigure9 <- function(collapsemm=TRUE) {
  PlotNucleotideContributionMatrix("miR-155", "equil_sc_nb", -4, collapsemm=collapsemm, pdf.plot="7.Ci", mirna_label=TRUE, extralabel=TRUE)
  PlotNucleotideContributionMatrix("miR-155", "equil_sc_nb", -3, collapsemm=collapsemm, pdf.plot="7.Cii")
  PlotNucleotideContributionMatrix("miR-155", "equil_sc_nb", -2, collapsemm=collapsemm, pdf.plot="7.Ciii")
  PlotNucleotideContributionMatrix("miR-155", "equil_sc_nb", -1, collapsemm=collapsemm, pdf.plot="7.Civ")
  PlotNucleotideContributionMatrix("miR-155", "equil_sc_nb",  0, collapsemm=collapsemm, pdf.plot="7.Cv")
  PlotNucleotideContributionMatrix("miR-155", "equil_sc_nb",  1, collapsemm=collapsemm, pdf.plot="7.Cvi")

  PlotNucleotideContributionMatrix("miR-155", "equil_sc_nb",  2, collapsemm=collapsemm, pdf.plot="7.Cvii")
  PlotNucleotideContributionMatrix("miR-155", "equil_sc_nb",  3, collapsemm=collapsemm, pdf.plot="7.Cviii")
  PlotNucleotideContributionMatrix("miR-155", "equil_sc_nb",  4, collapsemm=collapsemm, pdf.plot="7.Cix")
  PlotNucleotideContributionMatrix("miR-155", "equil_sc_nb",  5, collapsemm=collapsemm, pdf.plot="7.Cx")
  PlotNucleotideContributionMatrix("miR-155", "equil_sc_nb",  6, collapsemm=collapsemm, pdf.plot="7.Cxi")
  PlotNucleotideContributionMatrix("miR-155", "equil_sc_nb",  7, collapsemm=collapsemm, pdf.plot="7.Cxii")

  PlotNucleotideContributionMatrix("miR-155", "equil_sc_nb",  8, collapsemm=collapsemm, pdf.plot="7.Cxiii")
  PlotNucleotideContributionMatrix("miR-155", "equil_sc_nb",  9, collapsemm=collapsemm, pdf.plot="7.Cxiv")
  PlotNucleotideContributionMatrix("miR-155", "equil_sc_nb", 10, collapsemm=collapsemm, pdf.plot="7.Cxv")
  PlotNucleotideContributionMatrix("miR-155", "equil_sc_nb", 11, collapsemm=collapsemm, pdf.plot="7.Cxvi")
  PlotNucleotideContributionMatrix("miR-155", "equil_sc_nb", 12, collapsemm=collapsemm, pdf.plot="7.Cxvii")
  PlotNucleotideContributionMatrix("miR-155", "equil_sc_nb", 13, collapsemm=collapsemm, pdf.plot="7.Cxviii")

  PlotSimNucleotideContributionMatrixNew("miR-155", "equil_sc_nb", -4, globalmodel=TRUE, collapsemm=collapsemm, pdf.plot="8.Ci", mirna_label=TRUE, extralabel=TRUE)
  PlotSimNucleotideContributionMatrixNew("miR-155", "equil_sc_nb", -3, collapsemm=collapsemm, pdf.plot="8.Cii")
  PlotSimNucleotideContributionMatrixNew("miR-155", "equil_sc_nb", -2, collapsemm=collapsemm, pdf.plot="8.Ciii")
  PlotSimNucleotideContributionMatrixNew("miR-155", "equil_sc_nb", -1, collapsemm=collapsemm, pdf.plot="8.Civ")
  PlotSimNucleotideContributionMatrixNew("miR-155", "equil_sc_nb",  0, collapsemm=collapsemm, pdf.plot="8.Cv")
  PlotSimNucleotideContributionMatrixNew("miR-155", "equil_sc_nb",  1, collapsemm=collapsemm, pdf.plot="8.Cvi")

  PlotSimNucleotideContributionMatrixNew("miR-155", "equil_sc_nb",  2, collapsemm=collapsemm, pdf.plot="8.Cvii")
  PlotSimNucleotideContributionMatrixNew("miR-155", "equil_sc_nb",  3, collapsemm=collapsemm, pdf.plot="8.Cviii")
  PlotSimNucleotideContributionMatrixNew("miR-155", "equil_sc_nb",  4, collapsemm=collapsemm, pdf.plot="8.Cix")
  PlotSimNucleotideContributionMatrixNew("miR-155", "equil_sc_nb",  5, collapsemm=collapsemm, pdf.plot="8.Cx")
  PlotSimNucleotideContributionMatrixNew("miR-155", "equil_sc_nb",  6, collapsemm=collapsemm, pdf.plot="8.Cxi")
  PlotSimNucleotideContributionMatrixNew("miR-155", "equil_sc_nb",  7, collapsemm=collapsemm, pdf.plot="8.Cxii")

  PlotSimNucleotideContributionMatrixNew("miR-155", "equil_sc_nb",  8, collapsemm=collapsemm, pdf.plot="8.Cxiii")
  PlotSimNucleotideContributionMatrixNew("miR-155", "equil_sc_nb",  9, collapsemm=collapsemm, pdf.plot="8.Cxiv")
  PlotSimNucleotideContributionMatrixNew("miR-155", "equil_sc_nb", 10, collapsemm=collapsemm, pdf.plot="8.Cxv")
  PlotSimNucleotideContributionMatrixNew("miR-155", "equil_sc_nb", 11, collapsemm=collapsemm, pdf.plot="8.Cxvi")
  PlotSimNucleotideContributionMatrixNew("miR-155", "equil_sc_nb", 12, collapsemm=collapsemm, pdf.plot="8.Cxvii")
  PlotSimNucleotideContributionMatrixNew("miR-155", "equil_sc_nb", 13, collapsemm=collapsemm, pdf.plot="8.Cxviii")


}

MakeFigure6 <- function() {
  # PlotRegisterAndPositionKdFoldChange("let-7a-21nt", "equil_c2_nb", kdrel=TRUE,
  #                                     len_k=5, mirna_label=TRUE,
  #                                     pdf.plot="6.Bi")
  # PlotRegisterAndPositionKdFoldChange("let-7a-21nt", "equil_c2_nb", kdrel=TRUE,
  #                                     len_k=6, pdf.plot="6.Bii")
  # PlotRegisterAndPositionKdFoldChange("let-7a-21nt", "equil_c2_nb", kdrel=TRUE,
  #                                     len_k=7, pdf.plot="6.Biii")
  # PlotRegisterAndPositionKdFoldChange("let-7a-21nt", "equil_c2_nb", kdrel=TRUE,
  #                                     len_k=8, pdf.plot="6.Biv")
  # PlotRegisterAndPositionKdFoldChange("let-7a-21nt", "equil_c2_nb", kdrel=TRUE,
  #                                     len_k=9, pdf.plot="6.Bv")


  # PlotRegisterAndPositionKdFoldChange("miR-155_let-7a", "equil_c_nb", kdrel=TRUE,
  #                                     len_k=5, mirna_label=TRUE,
  #                                     pdf.plot="6.Bvi")
  # PlotRegisterAndPositionKdFoldChange("miR-155_let-7a", "equil_c_nb", kdrel=TRUE,
  #                                     len_k=6, pdf.plot="6.Bvii")
  # PlotRegisterAndPositionKdFoldChange("miR-155_let-7a", "equil_c_nb", kdrel=TRUE,
  #                                     len_k=7, pdf.plot="6.Bviii")
  # PlotRegisterAndPositionKdFoldChange("miR-155_let-7a", "equil_c_nb", kdrel=TRUE,
  #                                     len_k=8, pdf.plot="6.Bix")
  # PlotRegisterAndPositionKdFoldChange("miR-155_let-7a", "equil_c_nb", kdrel=TRUE,
  #                                     len_k=9, pdf.plot="6.Bx")

  # PlotRegisterAndPositionKdFoldChange("miR-155", "equil_sc_nb", kdrel=TRUE,
  #                                     len_k=5, mirna_label=TRUE, 
  #                                     pdf.plot="6.Bxi")
  # PlotRegisterAndPositionKdFoldChange("miR-155", "equil_sc_nb", kdrel=TRUE,
  #                                     len_k=6, pdf.plot="6.Bxii")
  # PlotRegisterAndPositionKdFoldChange("miR-155", "equil_sc_nb", kdrel=TRUE,
  #                                     len_k=7, pdf.plot="6.Bxiii")
  # PlotRegisterAndPositionKdFoldChange("miR-155", "equil_sc_nb", kdrel=TRUE,
  #                                     len_k=8, pdf.plot="6.Bxiv")
  # PlotRegisterAndPositionKdFoldChange("miR-155", "equil_sc_nb", kdrel=TRUE,
  #                                     len_k=9, pdf.plot="6.Bxv")

  # PlotRegisterAndPositionKdFoldChange("let-7a_miR-155", "equil_c_nb", kdrel=TRUE,
  #                                     len_k=5, mirna_label=TRUE,
  #                                     pdf.plot="6.Bxvi")
  # PlotRegisterAndPositionKdFoldChange("let-7a_miR-155", "equil_c_nb", kdrel=TRUE,
  #                                     len_k=6, pdf.plot="6.Bxvii")
  # PlotRegisterAndPositionKdFoldChange("let-7a_miR-155", "equil_c_nb", kdrel=TRUE,
  #                                     len_k=7, pdf.plot="6.Bxviii")
  # PlotRegisterAndPositionKdFoldChange("let-7a_miR-155", "equil_c_nb", kdrel=TRUE,
  #                                     len_k=8, pdf.plot="6.Bxix")
  # PlotRegisterAndPositionKdFoldChange("let-7a_miR-155", "equil_c_nb", kdrel=TRUE,
  #                                     len_k=9, pdf.plot="6.Bxx")

  PlotBestRegisterKd(average_win=5, chimera_series=TRUE, kdrel=TRUE,
                     pdf.plot="6.C")
  PlotBestRegisterKd(average_win=5, chimera_let7_series=TRUE, kdrel=TRUE,
                     pdf.plot="6.Ci")
  PlotBestRegisterKd(average_win=5, chimera_miR155_series=TRUE, kdrel=TRUE,
                     pdf.plot="6.Cii")

  PlotMismatchesByPositionstKd("let-7a-21nt", "equil_c2_nb", len_k=5,
                               kdrel=TRUE, pdf.plot="6.Di")
  PlotMismatchesByPositionstKd("let-7a-21nt", "equil_c2_nb", len_k=6,
                               kdrel=TRUE, pdf.plot="6.Dii")
  PlotMismatchesByPositionstKd("let-7a-21nt", "equil_c2_nb", len_k=7,
                               kdrel=TRUE, pdf.plot="6.Diii")
  PlotMismatchesByPositionstKd("let-7a-21nt", "equil_c2_nb", len_k=8,
                               kdrel=TRUE, pdf.plot="6.Div")
  PlotMismatchesByPositionstKd("let-7a-21nt", "equil_c2_nb", len_k=9,
                               kdrel=TRUE, pdf.plot="6.Dv")

  PlotMismatchesByPositionstKd("miR-155_let-7a", "equil_c_nb", len_k=5,
                               kdrel=TRUE, pdf.plot="6.Dvi")
  PlotMismatchesByPositionstKd("miR-155_let-7a", "equil_c_nb", len_k=6,
                               kdrel=TRUE, pdf.plot="6.Dvii")
  PlotMismatchesByPositionstKd("miR-155_let-7a", "equil_c_nb", len_k=7,
                               kdrel=TRUE, pdf.plot="6.Dviii")
  PlotMismatchesByPositionstKd("miR-155_let-7a", "equil_c_nb", len_k=8,
                               kdrel=TRUE, pdf.plot="6.Dix")
  PlotMismatchesByPositionstKd("miR-155_let-7a", "equil_c_nb", len_k=9,
                               kdrel=TRUE, pdf.plot="6.Dx")


  PlotMismatchesByPositionstKd("miR-155", "equil_sc_nb", len_k=5,
                               kdrel=TRUE, pdf.plot="6.Dxi")
  PlotMismatchesByPositionstKd("miR-155", "equil_sc_nb", len_k=6,
                               kdrel=TRUE, pdf.plot="6.Dxii")
  PlotMismatchesByPositionstKd("miR-155", "equil_sc_nb", len_k=7,
                               kdrel=TRUE, pdf.plot="6.Dxiii")
  PlotMismatchesByPositionstKd("miR-155", "equil_sc_nb", len_k=8,
                               kdrel=TRUE, pdf.plot="6.Dxiv")
  PlotMismatchesByPositionstKd("miR-155", "equil_sc_nb", len_k=9,
                               kdrel=TRUE, pdf.plot="6.Dxv")


  PlotMismatchesByPositionstKd("let-7a_miR-155", "equil_c_nb", len_k=5,
                               kdrel=TRUE, pdf.plot="6.Dxvi")
  PlotMismatchesByPositionstKd("let-7a_miR-155", "equil_c_nb", len_k=6,
                               kdrel=TRUE, pdf.plot="6.Dxvii")
  PlotMismatchesByPositionstKd("let-7a_miR-155", "equil_c_nb", len_k=7,
                               kdrel=TRUE, pdf.plot="6.Dxviii")
  PlotMismatchesByPositionstKd("let-7a_miR-155", "equil_c_nb", len_k=8,
                               kdrel=TRUE, pdf.plot="6.Dxix")
  PlotMismatchesByPositionstKd("let-7a_miR-155", "equil_c_nb", len_k=9,
                               kdrel=TRUE, pdf.plot="6.Dxx")
}

MakeFigure6 <- function() {
  message("Making Fig. 5")

  PlotNucleotideContributionMatrix("let-7a_minus1", "equil_c_nb", 3,
                                   pdf.plot="5.Bi")
  PlotNucleotideContributionMatrix("let-7a_minus1", "equil_c_nb", 3,
                                   model_values=TRUE, pdf.plot="5.Bii")
  PlotOffsetAndPairingModelMatrix("let-7a_minus1", "equil_c_nb",
                                  globalmodel=FALSE, pdf.plot="5.Biii")
  PlotSingleOffsetCoefficients("let-7a_minus1", "equil_c_nb", globalmodel=FALSE,
                                  pdf.plot="5.Biv")

  PlotNucleotideContributionMatrix("let-7a-21nt", "equil_c2_nb", 4,
                                   pdf.plot="5.Ci")
  PlotNucleotideContributionMatrix("let-7a-21nt", "equil_c2_nb", 4,
                                   model_values=TRUE, pdf.plot="5.Cii")
  PlotOffsetAndPairingModelMatrix("let-7a-21nt", "equil_c2_nb",
                                  globalmodel=FALSE, pdf.plot="5.Ciii")
  PlotSingleOffsetCoefficients("let-7a-21nt", "equil_c2_nb", globalmodel=FALSE,
                               pdf.plot="5.Civ")

  PlotNucleotideContributionMatrix("let-7a_plus1", "equil_c_nb", 1,
                                   pdf.plot="5.Di")
  PlotNucleotideContributionMatrix("let-7a_plus1", "equil_c_nb", 1,
                                   model_values=TRUE, pdf.plot="5.Dii")
  PlotOffsetAndPairingModelMatrix("let-7a_plus1", "equil_c_nb",
                                  globalmodel=FALSE, pdf.plot="5.Diii")
  PlotSingleOffsetCoefficients("let-7a_plus1", "equil_c_nb", globalmodel=FALSE,
                               pdf.plot="5.Div")
  message("Done Fig. 5")
}



# END OF SPECIFIC FIGURE PLANEL FUNCTIONS. #####################################
MakeFigureN1 <- function() {
  PlotPairwiseInputFrequencies("let-7a-21nt", "equil_c_nb",
                               "let-7a-21nt", "equil_s_nb",
                               pdf.plot="N1.A")
  PlotPairwiseInputFrequencies("let-7a-21nt", "equil_c_nb",
                               "let-7a-21nt", "equil_c2_nb",
                               pdf.plot="N1.B")
  PlotPairwiseInputFrequencies("let-7a-21nt", "equil_c_nb",
                               "let-7a_plus1", "equil_c_nb",
                               pdf.plot="N1.C")
  PlotPairwiseInputFrequencies("let-7a-21nt", "equil_c_nb",
                               "let-7a_minus1", "equil_c_nb",
                               pdf.plot="N1.D")
  PlotPairwiseInputFrequencies("let-7a-21nt", "equil_c_nb",
                               "let-7a-21nt", "equil_sc_nb",
                               pdf.plot="N1.E")
  PlotPairwiseInputFrequencies("let-7a-21nt", "equil_c_nb",
                               "let-7a_miR-155", "equil_c_nb",
                               pdf.plot="N1.F")
  PlotPairwiseInputFrequencies("miR-1", "equil_c_nb",
                               "miR-1", "equil_sc_nb",
                               pdf.plot="N1.G")
  PlotPairwiseInputFrequencies("miR-155", "equil_sc_nb",
                               "miR-155_let-7a", "equil_c_nb",
                               pdf.plot="N1.H")
}

MakeFigureN2 <- function() {
  PlotProgAndRandKds("let-7a-21nt", "equil_s_nb", pdf.plot="N2.A")
  PlotProgAndRandKds("let-7a-21nt", "equil_c_nb", pdf.plot="N2.B")
  PlotProgAndRandKds("let-7a-21nt", "equil_c2_nb", pdf.plot="N2.C")

  PlotProgAndRandKds("let-7a-21nt", "equil_sc_nb", pdf.plot="N2.D")
  PlotProgAndRandKds("let-7a-21nt", "equil_sc_nb", equilibrium_nb=TRUE,
                     pdf.plot="N2.E")
  PlotProgAndRandKds("let-7a_minus1", "equil_c_nb", pdf.plot="N2.F")
  
  PlotProgAndRandKds("let-7a_plus1", "equil_c_nb", pdf.plot="N2.G")
  PlotProgAndRandKds("let-7a_miR-155", "equil_c_nb", pdf.plot="N2.H")
  PlotProgAndRandKds("miR-1", "equil_sc_nb", pdf.plot="N2.I")

  PlotProgAndRandKds("miR-1", "equil_c_nb", pdf.plot="N2.J")
  PlotProgAndRandKds("miR-155", "equil_sc_nb", pdf.plot="N2.K")
  PlotProgAndRandKds("miR-155_let-7a", "equil_c_nb", pdf.plot="N2.L")

}

MakeFigureN3 <- function() {
  PlotProgAndRandKds("let-7a-21nt", "equil_s_nb", prog_mm=TRUE,
                     pdf.plot="N3.A")
  PlotProgAndRandKds("let-7a-21nt", "equil_c_nb", prog_mm=TRUE,
                     pdf.plot="N3.B")
  PlotProgAndRandKds("let-7a-21nt", "equil_c2_nb", prog_mm=TRUE,
                     pdf.plot="N3.C")

  PlotProgAndRandKds("let-7a-21nt", "equil_sc_nb", prog_mm=TRUE,
                     pdf.plot="N3.D")
  PlotProgAndRandKds("let-7a-21nt", "equil_sc_nb", equilibrium_nb=TRUE,
                     prog_mm=TRUE, pdf.plot="N3.E")
  PlotProgAndRandKds("let-7a_minus1", "equil_c_nb", prog_mm=TRUE,
                     pdf.plot="N3.F")
  
  PlotProgAndRandKds("let-7a_plus1", "equil_c_nb", prog_mm=TRUE,
                     pdf.plot="N3.G")
  PlotProgAndRandKds("let-7a_miR-155", "equil_c_nb", prog_mm=TRUE,
                     pdf.plot="N3.H")
  PlotProgAndRandKds("miR-1", "equil_sc_nb", prog_mm=TRUE, pdf.plot="N3.I")

  PlotProgAndRandKds("miR-1", "equil_c_nb", prog_mm=TRUE, pdf.plot="N3.J")
  PlotProgAndRandKds("miR-155", "equil_sc_nb", prog_mm=TRUE, pdf.plot="N3.K")
  PlotProgAndRandKds("miR-155_let-7a", "equil_c_nb", prog_mm=TRUE,
                     pdf.plot="N3.L")

}

MakeFigureN4 <- function() {
  PlotOneThreePrimeSite( 4, "let-7a-21nt", "equil_c2_nb", pdf.plot="N4.A")
  PlotOneThreePrimeSite( 5, pdf.plot="N4.B")
  PlotOneThreePrimeSite( 6, pdf.plot="N4.C")
  PlotOneThreePrimeSite( 7, pdf.plot="N4.D")
  PlotOneThreePrimeSite( 8, pdf.plot="N4.E")
  PlotOneThreePrimeSite( 9, pdf.plot="N4.F")
  PlotOneThreePrimeSite(10, pdf.plot="N4.G")
  PlotOneThreePrimeSite(11, pdf.plot="N4.H")
}

MakeFigureN5 <- function() {
  PlotOneThreePrimeSite( 4, mirna="miR-1", exp="equil_c_nb", pdf.plot="N5.A")
  PlotOneThreePrimeSite( 5, mirna="miR-1", exp="equil_c_nb", pdf.plot="N5.B")
  PlotOneThreePrimeSite( 6, mirna="miR-1", exp="equil_c_nb", pdf.plot="N5.C")
  PlotOneThreePrimeSite( 7, mirna="miR-1", exp="equil_c_nb", pdf.plot="N5.D")
  PlotOneThreePrimeSite( 8, mirna="miR-1", exp="equil_c_nb", pdf.plot="N5.E")
  PlotOneThreePrimeSite( 9, mirna="miR-1", exp="equil_c_nb", pdf.plot="N5.F")
  PlotOneThreePrimeSite(10, mirna="miR-1", exp="equil_c_nb", pdf.plot="N5.G")
  PlotOneThreePrimeSite(11, mirna="miR-1", exp="equil_c_nb", pdf.plot="N5.H")
}

MakeFigureN6 <- function() {
  PlotOneThreePrimeSite( 4, mirna="miR-155", exp="equil_sc_nb", pdf.plot="N6.A")
  PlotOneThreePrimeSite( 5, mirna="miR-155", exp="equil_sc_nb", pdf.plot="N6.B")
  PlotOneThreePrimeSite( 6, mirna="miR-155", exp="equil_sc_nb", pdf.plot="N6.C")
  PlotOneThreePrimeSite( 7, mirna="miR-155", exp="equil_sc_nb", pdf.plot="N6.D")
  PlotOneThreePrimeSite( 8, mirna="miR-155", exp="equil_sc_nb", pdf.plot="N6.E")
  PlotOneThreePrimeSite( 9, mirna="miR-155", exp="equil_sc_nb", pdf.plot="N6.F")
  PlotOneThreePrimeSite(10, mirna="miR-155", exp="equil_sc_nb", pdf.plot="N6.G")
  PlotOneThreePrimeSite(11, mirna="miR-155", exp="equil_sc_nb", pdf.plot="N6.H")
}

MakeFigureN7 <- function() {
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=4, suppcomp=FALSE, pdf.plot="N7.A")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=5, suppcomp=FALSE, pdf.plot="N7.B")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=6, suppcomp=FALSE, pdf.plot="N7.C")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=7, suppcomp=FALSE, pdf.plot="N7.D")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=8, suppcomp=FALSE, pdf.plot="N7.E")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=9, suppcomp=FALSE, pdf.plot="N7.F")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=10, suppcomp=FALSE, pdf.plot="N7.G")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=11, suppcomp=FALSE, pdf.plot="N7.H")
}

MakeFigureN8 <- function() {
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=4, suppcomp=FALSE, prog_mm=TRUE, pdf.plot="N8.A")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=5, suppcomp=FALSE, prog_mm=TRUE, pdf.plot="N8.B")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=6, suppcomp=FALSE, prog_mm=TRUE, pdf.plot="N8.C")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=7, suppcomp=FALSE, prog_mm=TRUE, pdf.plot="N8.D")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=8, suppcomp=FALSE, prog_mm=TRUE, pdf.plot="N8.E")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=9, suppcomp=FALSE, prog_mm=TRUE, pdf.plot="N8.F")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=10, suppcomp=FALSE, prog_mm=TRUE, pdf.plot="N8.G")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=11, suppcomp=FALSE, prog_mm=TRUE, pdf.plot="N8.H")
}

MakeFigureN9 <- function() {
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=4, pdf.plot="N9.A")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=5, pdf.plot="N9.B")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=6, pdf.plot="N9.C")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=7, pdf.plot="N9.D")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=8, pdf.plot="N9.E")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=9, pdf.plot="N9.F")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=10, pdf.plot="N9.G")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=11, pdf.plot="N9.H")
}

MakeFigureN10 <- function() {
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=4, prog_mm=TRUE, pdf.plot="N10.A")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=5, prog_mm=TRUE, pdf.plot="N10.B")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=6, prog_mm=TRUE, pdf.plot="N10.C")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=7, prog_mm=TRUE, pdf.plot="N10.D")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=8, prog_mm=TRUE, pdf.plot="N10.E")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=9, prog_mm=TRUE, pdf.plot="N10.F")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=10, prog_mm=TRUE, pdf.plot="N10.G")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, rand=TRUE, kmer_max=11, prog_mm=TRUE, pdf.plot="N10.H")
}

MakeFigureN11 <- function() {
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, suppcomp=FALSE, kmer_max=4, plotseed=FALSE, pdf.plot="N11.A")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, suppcomp=FALSE, kmer_max=5, plotseed=FALSE, pdf.plot="N11.B")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, suppcomp=FALSE, kmer_max=6, plotseed=FALSE, pdf.plot="N11.C")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, suppcomp=FALSE, kmer_max=7, plotseed=FALSE, pdf.plot="N11.D")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, suppcomp=FALSE, kmer_max=8, plotseed=FALSE, pdf.plot="N11.E")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, suppcomp=FALSE, kmer_max=9, plotseed=FALSE, pdf.plot="N11.F")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, suppcomp=FALSE, kmer_max=10, plotseed=FALSE, pdf.plot="N11.G")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, suppcomp=FALSE, kmer_max=11, plotseed=FALSE, pdf.plot="N11.H")
}
MakeFigureN12 <- function() {
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, kmer_max=4, plotseed=FALSE, pdf.plot="N12.A")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, kmer_max=5, plotseed=FALSE, pdf.plot="N12.B")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, kmer_max=6, plotseed=FALSE, pdf.plot="N12.C")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, kmer_max=7, plotseed=FALSE, pdf.plot="N12.D")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, kmer_max=8, plotseed=FALSE, pdf.plot="N12.E")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, kmer_max=9, plotseed=FALSE, pdf.plot="N12.F")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, kmer_max=10, plotseed=FALSE, pdf.plot="N12.G")
  PlotLet7aProgrammedLibraryReps(allmm=FALSE, kmer_max=11, plotseed=FALSE, pdf.plot="N12.H")  
}

MakeFigureN13 <- function() {
  PlotPairingRangeMatrix("let-7a-21nt", "equil_c2_nb", 0, loop=TRUE, pdf.plot="N13.A")
  PlotPairingRangeMatrix("let-7a-21nt", "equil_c2_nb", 2, loop=TRUE, pdf.plot="N13.B")
  PlotPairingRangeMatrix("let-7a-21nt", "equil_c2_nb", 4, loop=TRUE, pdf.plot="N13.C")
  PlotPairingRangeMatrix("let-7a-21nt", "equil_c2_nb", 6, loop=TRUE, pdf.plot="N13.D")
  PlotPairingRangeMatrix("let-7a-21nt", "equil_c2_nb", 8, loop=TRUE, pdf.plot="N13.E")
  PlotPairingRangeMatrix("let-7a-21nt", "equil_c2_nb", 10, loop=TRUE, pdf.plot="N13.F")
  PlotPairingRangeMatrix("let-7a-21nt", "equil_c2_nb", 12, loop=TRUE, pdf.plot="N13.G")
  PlotPairingRangeMatrix("let-7a-21nt", "equil_c2_nb", 14, loop=TRUE, pdf.plot="N13.H")
  PlotPairingRangeMatrix("let-7a-21nt", "equil_c2_nb", 16, loop=TRUE, pdf.plot="N13.I")
  PlotPairingRangeMatrix("let-7a-21nt", "equil_c2_nb", 18, loop=TRUE, pdf.plot="N13.J")
  PlotPairingRangeMatrix("let-7a-21nt", "equil_c2_nb", 20, loop=TRUE, pdf.plot="N13.K")
}

MakeFigureN14 <- function() {

  MismatchBarPlotModel("miR-155", "equil_sc_nb", pdf.plot="N14.A")
  PlotPairingModelCoefficients("miR-155", "equil_sc_nb", mm=TRUE,
                               globalmodel=FALSE, pdf.plot="N14.B")
  PlotOffsetCoefficients("miR-155", "equil_sc_nb", mm=TRUE, globalmodel=FALSE,
                         pdf.plot="N14.C")

  MismatchBarPlotModel("miR-155", "equil_sc_nb", pos_3p_max=18, len_max=7,
                       pdf.plot="N14.D")
  PlotPairingModelCoefficients("miR-155", "equil_sc_nb", mm=TRUE,
                               globalmodel=FALSE, pos_3p_max=18, len_max=7,
                               pdf.plot="N14.E")
  PlotOffsetCoefficients("miR-155", "equil_sc_nb", mm=TRUE, globalmodel=FALSE,
                         pos_3p_max=18, len_max=7, pdf.plot="N14.F")

  MismatchBarPlotModel("miR-155", "equil_sc_nb", pos_3p_min=21, len_min=8,
                       pdf.plot="N14.G")
  PlotPairingModelCoefficients("miR-155", "equil_sc_nb", globalmodel=FALSE,
                               pos_3p_min=21, len_min=8, pdf.plot="N14.H")
  PlotOffsetCoefficients("miR-155", "equil_sc_nb", mm=TRUE, globalmodel=FALSE,
                         pos_3p_min=21, len_min=8, pdf.plot="N14.I")
}

MakeFigureN15 <- function() {
  MismatchBarPlot("miR-155", "equil_sc_nb", "4mer-m13.16", pdf.plot="N15.A")
  MismatchBarPlot("miR-155", "equil_sc_nb", "5mer-m13.17", pdf.plot="N15.B")
  MismatchBarPlot("miR-155", "equil_sc_nb", "6mer-m13.18", pdf.plot="N15.C")
  MismatchBarPlot("miR-155", "equil_sc_nb", "7mer-m15.21", pdf.plot="N15.D")
  MismatchBarPlot("miR-155", "equil_sc_nb", "8mer-m15.22", pdf.plot="N15.E")
  MismatchBarPlot("miR-155", "equil_sc_nb", "9mer-m13.21", pdf.plot="N15.F")
  MismatchBarPlot("miR-155", "equil_sc_nb", "10mer-m13.22", pdf.plot="N15.G")
  MismatchBarPlot("miR-155", "equil_sc_nb", "11mer-m13.23", pdf.plot="N15.H")
}


# MakeFigureN4()
# MakeFigureN5()
# MakeFigureN6()
# MakeFigureN7()
# MakeFigureN8()
# MakeFigureN9()
# MakeFigureN10()
# MakeFigureN11()
# MakeFigureN12()
# MakeFigureN14()
# MakeFigureN15()



# PlotEmpiricalMismatchCoefficients(
#   "let-7a", "equilibrium", sitelist="bipartite_random", corrected_kds=FALSE,
#   supp_base=TRUE, add_comp_sites=TRUE, collapsemm_addcomp=TRUE, len_max=5,
#   offsetmin=0, offsetmax=10)



# MakeSupplementalFigure7 <- function() {
#   PlotPairingModelCoefficients(
#     "let-7a", "equilibrium", sitelist="bipartite_random_suppcomp", supp_base=TRUE,
#     corrected_kds=FALSE, pdf.plot="S7.Ai"
#   )
#   PlotOffsetCoefficients(
#     "let-7a", "equilibrium", sitelist="bipartite_random_suppcomp", supp_base=TRUE,
#     globalmodel=FALSE, corrected_kds=FALSE, pdf.plot="S7.Aii"
#   )
#   PlotPairingRangeMatrix(
#     "let-7a", "equilibrium", 3, sitelist="bipartite_random_suppcomp", supp_base=TRUE,
#     model_values=TRUE, globalmodel=FALSE, corrected_kds=FALSE,
#     pdf.plot="S7.Aiii"
#   )
#   PlotPairingRangeMatrix(
#     "let-7a", "equilibrium", 3, sitelist="bipartite_random_suppcomp", supp_base=TRUE,
#     globalmodel=FALSE, corrected_kds=FALSE, pdf.plot="S7.Aiv"
#   )

#   PlotPairingModelCoefficients(
#     "miR-1", "equilibrium", sitelist="bipartite_random_suppcomp", supp_base=TRUE,
#     corrected_kds=FALSE, pdf.plot="S7.Bi"
#   )
#   PlotOffsetCoefficients(
#     "miR-1", "equilibrium", sitelist="bipartite_random_suppcomp", supp_base=TRUE,
#     globalmodel=FALSE, corrected_kds=FALSE, pdf.plot="S7.Bii"
#   )
#   PlotPairingRangeMatrix(
#     "miR-1", "equilibrium", -1, sitelist="bipartite_random_suppcomp", supp_base=TRUE,
#     model_values=TRUE, globalmodel=FALSE, corrected_kds=FALSE,
#     pdf.plot="S7.Biii"
#   )
#   PlotPairingRangeMatrix(
#     "miR-1", "equilibrium", -1, sitelist="bipartite_random_suppcomp", supp_base=TRUE,
#     globalmodel=FALSE, corrected_kds=FALSE, pdf.plot="S7.Biv"
#   )

#   PlotPairingModelCoefficients(
#     "miR-155", "equilibrium", sitelist="bipartite_random_suppcomp", supp_base=TRUE,
#     corrected_kds=FALSE, pdf.plot="S7.Ci"
#   )
#   PlotOffsetCoefficients(
#     "miR-155", "equilibrium", sitelist="bipartite_random_suppcomp", supp_base=TRUE,
#     globalmodel=FALSE, corrected_kds=FALSE, pdf.plot="S7.Cii"
#   )
#   PlotPairingRangeMatrix(
#     "miR-155", "equilibrium", -1, sitelist="bipartite_random_suppcomp", supp_base=TRUE,
#     model_values=TRUE, globalmodel=FALSE, corrected_kds=FALSE,
#     pdf.plot="S7.Ciii"
#   )
#   PlotPairingRangeMatrix(
#     "miR-155", "equilibrium", -1, sitelist="bipartite_random_suppcomp", supp_base=TRUE,
#     globalmodel=FALSE, corrected_kds=FALSE, pdf.plot="S7.Civ"
#   )

#   PlotPairingModelCoefficients(
#     "miR-124", "equilibrium", sitelist="bipartite_random_suppcomp", supp_base=TRUE,
#     corrected_kds=FALSE, pdf.plot="S7.Di"
#   )
#   PlotOffsetCoefficients(
#     "miR-124", "equilibrium", sitelist="bipartite_random_suppcomp", supp_base=TRUE,
#     globalmodel=FALSE, corrected_kds=FALSE, pdf.plot="S7.Dii"
#   )
#   PlotPairingRangeMatrix(
#     "miR-124", "equilibrium", 3, sitelist="bipartite_random_suppcomp", supp_base=TRUE,
#     model_values=TRUE, globalmodel=FALSE, corrected_kds=FALSE,
#     pdf.plot="S7.Diii"
#   )
#   PlotPairingRangeMatrix(
#     "miR-124", "equilibrium", 3, sitelist="bipartite_random_suppcomp", supp_base=TRUE,
#     globalmodel=FALSE, corrected_kds=FALSE, pdf.plot="S7.Div"
#   )

#   PlotPairingModelCoefficients(
#     "lsy-6", "equilibrium", sitelist="bipartite_random_suppcomp", supp_base=TRUE,
#     corrected_kds=FALSE, pdf.plot="S7.Ei"
#   )
#   PlotOffsetCoefficients(
#     "lsy-6", "equilibrium", sitelist="bipartite_random_suppcomp", supp_base=TRUE,
#     globalmodel=FALSE, corrected_kds=FALSE, pdf.plot="S7.Eii"
#   )
#   PlotPairingRangeMatrix(
#     "lsy-6", "equilibrium", -1, sitelist="bipartite_random_suppcomp", supp_base=TRUE,
#     model_values=TRUE, globalmodel=FALSE, corrected_kds=FALSE,
#     pdf.plot="S7.Eiii"
#   )
#   PlotPairingRangeMatrix(
#     "lsy-6", "equilibrium", -1, sitelist="bipartite_random_suppcomp", supp_base=TRUE,
#     globalmodel=FALSE, corrected_kds=FALSE, pdf.plot="S7.Eiv"
#   )

#   PlotPairingModelCoefficients(
#     "miR-7-23nt", "equilibrium2_nb", sitelist="bipartite_random_suppcomp", supp_base=TRUE,
#     corrected_kds=FALSE, pdf.plot="S7.Fi"
#   )
#   PlotOffsetCoefficients(
#     "miR-7-23nt", "equilibrium2_nb", sitelist="bipartite_random_suppcomp", supp_base=TRUE,
#     globalmodel=FALSE, corrected_kds=FALSE, pdf.plot="S7.Fii"
#   )
#   PlotPairingRangeMatrix(
#     "miR-7-23nt", "equilibrium2_nb", 3, sitelist="bipartite_random_suppcomp", supp_base=TRUE,
#     model_values=TRUE, globalmodel=FALSE, corrected_kds=FALSE,
#     pdf.plot="S7.Fiii"
#   )
#   PlotPairingRangeMatrix(
#     "miR-7-23nt", "equilibrium2_nb", 3, sitelist="bipartite_random_suppcomp", supp_base=TRUE,
#     globalmodel=FALSE, corrected_kds=FALSE, pdf.plot="S7.Fiv"
#   )
#   message("Done fig. S7.")
# }

# MakeSupplementalFigure8 <- function() {
#   PlotPairingModelCoefficients(
#     "let-7a", "equilibrium", sitelist="bipartite_random_suppcomp", offset_base=TRUE,
#     corrected_kds=FALSE, pdf.plot="S8.Ai"
#   )
#   PlotOffsetCoefficients(
#     "let-7a", "equilibrium", sitelist="bipartite_random_suppcomp", offset_base=TRUE,
#     globalmodel=FALSE, corrected_kds=FALSE, pdf.plot="S8.Aii"
#   )
#   PlotPairingRangeMatrix(
#     "let-7a", "equilibrium", 3, sitelist="bipartite_random_suppcomp", offset_base=TRUE,
#     model_values=TRUE, globalmodel=FALSE, corrected_kds=FALSE,
#     pdf.plot="S8.Aiii"
#   )
#   PlotPairingRangeMatrix(
#     "let-7a", "equilibrium", 3, sitelist="bipartite_random_suppcomp", offset_base=TRUE,
#     globalmodel=FALSE, corrected_kds=FALSE, pdf.plot="S8.Aiv"
#   )

#   PlotPairingModelCoefficients(
#     "miR-1", "equilibrium", sitelist="bipartite_random_suppcomp", offset_base=TRUE,
#     corrected_kds=FALSE, pdf.plot="S8.Bi"
#   )
#   PlotOffsetCoefficients(
#     "miR-1", "equilibrium", sitelist="bipartite_random_suppcomp", offset_base=TRUE,
#     globalmodel=FALSE, corrected_kds=FALSE, pdf.plot="S8.Bii"
#   )
#   PlotPairingRangeMatrix(
#     "miR-1", "equilibrium", -1, sitelist="bipartite_random_suppcomp", offset_base=TRUE,
#     model_values=TRUE, globalmodel=FALSE, corrected_kds=FALSE,
#     pdf.plot="S8.Biii"
#   )
#   PlotPairingRangeMatrix(
#     "miR-1", "equilibrium", -1, sitelist="bipartite_random_suppcomp", offset_base=TRUE,
#     globalmodel=FALSE, corrected_kds=FALSE, pdf.plot="S8.Biv"
#   )

#   PlotPairingModelCoefficients(
#     "miR-155", "equilibrium", sitelist="bipartite_random_suppcomp", offset_base=TRUE,
#     corrected_kds=FALSE, pdf.plot="S8.Ci"
#   )
#   PlotOffsetCoefficients(
#     "miR-155", "equilibrium", sitelist="bipartite_random_suppcomp", offset_base=TRUE,
#     globalmodel=FALSE, corrected_kds=FALSE, pdf.plot="S8.Cii"
#   )
#   PlotPairingRangeMatrix(
#     "miR-155", "equilibrium", -1, sitelist="bipartite_random_suppcomp", offset_base=TRUE,
#     model_values=TRUE, globalmodel=FALSE, corrected_kds=FALSE,
#     pdf.plot="S8.Ciii"
#   )
#   PlotPairingRangeMatrix(
#     "miR-155", "equilibrium", -1, sitelist="bipartite_random_suppcomp", offset_base=TRUE,
#     globalmodel=FALSE, corrected_kds=FALSE, pdf.plot="S8.Civ"
#   )

#   PlotPairingModelCoefficients(
#     "miR-124", "equilibrium", sitelist="bipartite_random_suppcomp", offset_base=TRUE,
#     corrected_kds=FALSE, pdf.plot="S8.Di"
#   )
#   PlotOffsetCoefficients(
#     "miR-124", "equilibrium", sitelist="bipartite_random_suppcomp", offset_base=TRUE,
#     globalmodel=FALSE, corrected_kds=FALSE, pdf.plot="S8.Dii"
#   )
#   PlotPairingRangeMatrix(
#     "miR-124", "equilibrium", 3, sitelist="bipartite_random_suppcomp", offset_base=TRUE,
#     model_values=TRUE, globalmodel=FALSE, corrected_kds=FALSE,
#     pdf.plot="S8.Diii"
#   )
#   PlotPairingRangeMatrix(
#     "miR-124", "equilibrium", 3, sitelist="bipartite_random_suppcomp", offset_base=TRUE,
#     globalmodel=FALSE, corrected_kds=FALSE, pdf.plot="S8.Div"
#   )

#   PlotPairingModelCoefficients(
#     "lsy-6", "equilibrium", sitelist="bipartite_random_suppcomp", offset_base=TRUE,
#     corrected_kds=FALSE, pdf.plot="S8.Ei"
#   )
#   PlotOffsetCoefficients(
#     "lsy-6", "equilibrium", sitelist="bipartite_random_suppcomp", offset_base=TRUE,
#     globalmodel=FALSE, corrected_kds=FALSE, pdf.plot="S8.Eii"
#   )
#   PlotPairingRangeMatrix(
#     "lsy-6", "equilibrium", -1, sitelist="bipartite_random_suppcomp", offset_base=TRUE,
#     model_values=TRUE, globalmodel=FALSE, corrected_kds=FALSE,
#     pdf.plot="S8.Eiii"
#   )
#   PlotPairingRangeMatrix(
#     "lsy-6", "equilibrium", -1, sitelist="bipartite_random_suppcomp", offset_base=TRUE,
#     globalmodel=FALSE, corrected_kds=FALSE, pdf.plot="S8.Eiv"
#   )

#   PlotPairingModelCoefficients(
#     "miR-7-23nt", "equilibrium2_nb", sitelist="bipartite_random_suppcomp", offset_base=TRUE,
#     corrected_kds=FALSE, pdf.plot="S8.Fi"
#   )
#   PlotOffsetCoefficients(
#     "miR-7-23nt", "equilibrium2_nb", sitelist="bipartite_random_suppcomp", offset_base=TRUE,
#     globalmodel=FALSE, corrected_kds=FALSE, pdf.plot="S8.Fii"
#   )
#   PlotPairingRangeMatrix(
#     "miR-7-23nt", "equilibrium2_nb", 3, sitelist="bipartite_random_suppcomp", offset_base=TRUE,
#     model_values=TRUE, globalmodel=FALSE, corrected_kds=FALSE,
#     pdf.plot="S8.Fiii"
#   )
#   PlotPairingRangeMatrix(
#     "miR-7-23nt", "equilibrium2_nb", 3, sitelist="bipartite_random_suppcomp", offset_base=TRUE,
#     globalmodel=FALSE, corrected_kds=FALSE, pdf.plot="S8.Fiv"
#   )
#   message("Done fig. S8.")
# }




# MakeSupplementalFigure4 <- function() {
#   message("Making Fig. S4")

#   PlotPairingRangeMatrix("let-7a-21nt", "equil_c2_nb", -4, mirna_label=TRUE,
#                          model_values=TRUE, label_offset=TRUE, pdf.plot="S4.Ai")
#   PlotPairingRangeMatrix("let-7a-21nt", "equil_c2_nb", -2, model_values=TRUE,
#                          label_offset=TRUE, globalmodel=FALSE,
#                          pdf.plot="S4.Aii")
#   PlotPairingRangeMatrix("let-7a-21nt", "equil_c2_nb",  0, model_values=TRUE,
#                          label_offset=TRUE, globalmodel=FALSE,
#                          pdf.plot="S4.Aiii")
#   PlotPairingRangeMatrix("let-7a-21nt", "equil_c2_nb",  2, model_values=TRUE,
#                          label_offset=TRUE, globalmodel=FALSE,
#                          pdf.plot="S4.Aiv")
#   PlotPairingRangeMatrix("let-7a-21nt", "equil_c2_nb",  4, model_values=TRUE,
#                          label_offset=TRUE, globalmodel=FALSE, pdf.plot="S4.Av")
#   PlotPairingRangeMatrix("let-7a-21nt", "equil_c2_nb",  6, model_values=TRUE,
#                          label_offset=TRUE, globalmodel=FALSE,
#                          pdf.plot="S4.Avi")
#   PlotPairingRangeMatrix("let-7a-21nt", "equil_c2_nb",  8, model_values=TRUE,
#                          label_offset=TRUE, globalmodel=FALSE,
#                          pdf.plot="S4.Avii")
#   PlotPairingRangeMatrix("let-7a-21nt", "equil_c2_nb", 10, model_values=TRUE,
#                          label_offset=TRUE, globalmodel=FALSE,
#                          pdf.plot="S4.Aviii")
#   PlotPairingRangeMatrix("let-7a-21nt", "equil_c2_nb", 12, model_values=TRUE,
#                          label_offset=TRUE, globalmodel=FALSE,
#                          pdf.plot="S4.Aix")
#   PlotPairingRangeMatrix("let-7a-21nt", "equil_c2_nb", 14, model_values=TRUE,
#                          label_offset=TRUE, globalmodel=FALSE, pdf.plot="S4.Ax")
#   PlotPairingRangeMatrix("let-7a-21nt", "equil_c2_nb", 16, key=TRUE,
#                          model_values=TRUE, label_offset=TRUE,
#                          globalmodel=FALSE, pdf.plot="S4.Axi")

#   PlotPairingRangeMatrix("miR-1", "equil_c_nb", -4, mirna_label=TRUE,
#                          model_values=TRUE, label_offset=TRUE, pdf.plot="S4.Bi")
#   PlotPairingRangeMatrix("miR-1", "equil_c_nb", -2, model_values=TRUE,
#                          label_offset=TRUE, globalmodel=FALSE,
#                          pdf.plot="S4.Bii")
#   PlotPairingRangeMatrix("miR-1", "equil_c_nb",  0, model_values=TRUE,
#                          label_offset=TRUE, globalmodel=FALSE,
#                          pdf.plot="S4.Biii")
#   PlotPairingRangeMatrix("miR-1", "equil_c_nb",  2, model_values=TRUE,
#                          label_offset=TRUE, globalmodel=FALSE,
#                          pdf.plot="S4.Biv")
#   PlotPairingRangeMatrix("miR-1", "equil_c_nb",  4, model_values=TRUE,
#                          label_offset=TRUE, globalmodel=FALSE, pdf.plot="S4.Bv")
#   PlotPairingRangeMatrix("miR-1", "equil_c_nb",  6, model_values=TRUE,
#                          label_offset=TRUE, globalmodel=FALSE,
#                          pdf.plot="S4.Bvi")
#   PlotPairingRangeMatrix("miR-1", "equil_c_nb",  8, model_values=TRUE,
#                          label_offset=TRUE, globalmodel=FALSE,
#                          pdf.plot="S4.Bvii")
#   PlotPairingRangeMatrix("miR-1", "equil_c_nb", 10, model_values=TRUE,
#                          label_offset=TRUE, globalmodel=FALSE,
#                          pdf.plot="S4.Bviii")
#   PlotPairingRangeMatrix("miR-1", "equil_c_nb", 12, model_values=TRUE,
#                          label_offset=TRUE, globalmodel=FALSE,
#                          pdf.plot="S4.Bix")
#   PlotPairingRangeMatrix("miR-1", "equil_c_nb", 14, model_values=TRUE,
#                          label_offset=TRUE, globalmodel=FALSE, pdf.plot="S4.Bx")
#   PlotPairingRangeMatrix("miR-1", "equil_c_nb", 16, key=TRUE, model_values=TRUE,
#                          label_offset=TRUE, globalmodel=FALSE,
#                          pdf.plot="S4.Bxi")

#   PlotPairingRangeMatrix("miR-155", "equil_sc_nb", -4, mirna_label=TRUE,
#                          model_values=TRUE, label_offset=TRUE,
#                          pdf.plot="S4.Ci")
#   PlotPairingRangeMatrix("miR-155", "equil_sc_nb", -2, model_values=TRUE,
#                          label_offset=TRUE, globalmodel=FALSE,
#                          pdf.plot="S4.Cii")
#   PlotPairingRangeMatrix("miR-155", "equil_sc_nb",  0, model_values=TRUE,
#                          label_offset=TRUE, globalmodel=FALSE,
#                          pdf.plot="S4.Ciii")
#   PlotPairingRangeMatrix("miR-155", "equil_sc_nb",  2, model_values=TRUE,
#                          label_offset=TRUE, globalmodel=FALSE,
#                          pdf.plot="S4.Civ")
#   PlotPairingRangeMatrix("miR-155", "equil_sc_nb",  4, model_values=TRUE,
#                          label_offset=TRUE, globalmodel=FALSE,
#                          pdf.plot="S4.Cv")
#   PlotPairingRangeMatrix("miR-155", "equil_sc_nb",  6, model_values=TRUE,
#                          label_offset=TRUE, globalmodel=FALSE,
#                          pdf.plot="S4.Cvi")
#   PlotPairingRangeMatrix("miR-155", "equil_sc_nb",  8, model_values=TRUE,
#                          label_offset=TRUE, globalmodel=FALSE,
#                          pdf.plot="S4.Cvii")
#   PlotPairingRangeMatrix("miR-155", "equil_sc_nb", 10, model_values=TRUE,
#                          label_offset=TRUE, globalmodel=FALSE,
#                          pdf.plot="S4.Cviii")
#   PlotPairingRangeMatrix("miR-155", "equil_sc_nb", 12, model_values=TRUE,
#                          label_offset=TRUE, globalmodel=FALSE,
#                          pdf.plot="S4.Cix")
#   PlotPairingRangeMatrix("miR-155", "equil_sc_nb", 14, model_values=TRUE,
#                          label_offset=TRUE, globalmodel=FALSE,
#                          pdf.plot="S4.Cx")
#   PlotPairingRangeMatrix("miR-155", "equil_sc_nb", 16, key=TRUE,
#                          model_values=TRUE, label_offset=TRUE,
#                          globalmodel=FALSE, pdf.plot="S4.Cxi")

#   PlotPairingAndOffsetModelAgainstData("let-7a-21nt", "equil_c2_nb",
#                                        pdf.plot="S4.Di")
#   PlotPairingAndOffsetModelAgainstData("miR-1", "equil_c_nb", pdf.plot="S4.Dii")
#   PlotPairingAndOffsetModelAgainstData("miR-155", "equil_sc_nb",
#                                        pdf.plot="S4.Diii")

#   message("Done fig. S4.")
# }


MakeSupplementalFigure5 <- function() {
  message("Making Fig. 5")
  PlotKdFoldChangeAgainstDeltaG("let-7a-21nt", "equil_c2_nb", lengthnorm=FALSE,
                                pdf.plot="S5.Ai")
  PlotKdFoldChangeAgainstDeltaG("miR-1", "equil_c_nb", lengthnorm=FALSE,
                                pdf.plot="S5.Aii")
  PlotKdFoldChangeAgainstDeltaG("miR-155", "equil_sc_nb", lengthnorm=FALSE,
                                pdf.plot="S5.Aiii")

  PlotKdFoldChangeAgainstDeltaG("let-7a-21nt", "equil_c2_nb", model_values=FALSE,
                                offset=4, pdf.plot="S5.Bi")
  PlotKdFoldChangeAgainstDeltaG("miR-1", "equil_c_nb", model_values=FALSE, 
                                offset=1, pdf.plot="S5.Bii")
  PlotKdFoldChangeAgainstDeltaG("miR-155", "equil_sc_nb", model_values=FALSE,
                                offset=1, pdf.plot="S5.Biii")

  PlotKdFoldChangeAgainstDeltaG("let-7a-21nt", "equil_c2_nb", model_values=FALSE,
                                offset=4, lengthnorm=FALSE, pdf.plot="S5.Ci")
  PlotKdFoldChangeAgainstDeltaG("miR-1", "equil_c_nb", model_values=FALSE, 
                                offset=1, lengthnorm=FALSE, pdf.plot="S5.Cii")
  PlotKdFoldChangeAgainstDeltaG("miR-155", "equil_sc_nb", model_values=FALSE,
                                offset=1, lengthnorm=FALSE, pdf.plot="S5.Ciii")

  PlotPairingOffsetAndMismatchModelAgainstData("let-7a-21nt", "equil_c2_nb",
                                               pdf.plot="S5.Di")
  PlotPairingOffsetAndMismatchModelAgainstData("miR-1", "equil_c_nb",
                                               pdf.plot="S5.Dii")
  PlotPairingOffsetAndMismatchModelAgainstData("miR-155", "equil_sc_nb",
                                               pdf.plot="S5.Diii")

  PlotPairingAndOffsetModelAgainstData("let-7a-21nt", "equil_c2_nb",
                                       pdf.plot="S5.Ei")
  PlotPairingAndOffsetModelAgainstData("miR-1", "equil_c_nb", pdf.plot="S5.Eii")
  PlotPairingAndOffsetModelAgainstData("miR-155", "equil_sc_nb",
                                       pdf.plot="S5.Eiii")



  message("Done fig. S5.")
}



MakeSupplementalFigure25 <- function(
) {
  # PlotMismatchCoefficients("let-7a-21nt", "equil_c2_nb", corrected_kds=TRUE, intercept=FALSE, exponential=FALSE, pdf.plot="S25.A")

  # PlotMismatchCoefficients("let-7a-21nt", "equil_c2_nb", corrected_kds=TRUE, intercept=TRUE, exponential=FALSE, pdf.plot="S25.B")

  # PlotMismatchCoefficients("let-7a-21nt", "equil_c2_nb", corrected_kds=TRUE, intercept=FALSE, exponential=TRUE, pdf.plot="S25.C")

  # PlotMismatchCoefficients("let-7a-21nt", "equil_c2_nb", corrected_kds=TRUE, intercept=TRUE, exponential=TRUE, pdf.plot="S25.D")

  # PlotMismatchCoefficients("let-7a-21nt", "equil_c2_nb", empirical=TRUE, corrected_kds=TRUE, pdf.plot="S25.E")

  # PlotAllRandomLibrarySeedCoefficientsAgainstKds(
  #   len_max=5, offsetmin=0, offsetmax=10, pos_3p_max=17, add_comp_sites=TRUE,
  #   intercept=FALSE, exponential=FALSE, collapsemm_addcomp=TRUE,
  #   pdf.plot="S25.F"
  # )
  # PlotAllRandomLibrarySeedCoefficientsAgainstKds(
  #   len_max=5, offsetmin=0, offsetmax=10, pos_3p_max=17, add_comp_sites=TRUE,
  #   intercept=TRUE, exponential=FALSE, collapsemm_addcomp=TRUE,
  #   pdf.plot="S25.G"
  # )

  PlotAllRandomLibrarySeedCoefficientsAgainstKds(
    len_max=5, offsetmin=0, offsetmax=10, add_comp_sites=TRUE,
    intercept=FALSE, exponential=TRUE, collapsemm_addcomp=TRUE,
    pdf.plot="S25.H"
  )

  # PlotAllRandomLibrarySeedCoefficientsAgainstKds(
  #   len_max=5, offsetmin=0, offsetmax=10, pos_3p_max=17, add_comp_sites=TRUE,
  #   intercept=TRUE, exponential=TRUE, collapsemm_addcomp=TRUE,
  #   pdf.plot="S25.I"
  # )

  PlotAllRandomLibrarySeedCoefficientsAgainstKds(add_comp_sites=TRUE,
    collapsemm_addcomp=TRUE, empirical=TRUE, pdf.plot="S25.J"
  )

  PlotAllRandomLibrarySeedCoefficientsAgainstKds(
    len_max=5, offsetmin=0, offsetmax=10, add_comp_sites=TRUE,
    collapsemm_addcomp=TRUE, empirical=TRUE, pdf.plot="S25.K"
  )


}



  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 9, 19, model_values=model_values, pdf.plot="S12.Ai")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 10, 20, model_values=model_values, pdf.plot="S12.Aii")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 11, 21, model_values=model_values, pdf.plot="S12.Aiii")

  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 9,  18, model_values=model_values, pdf.plot="S12.Bi")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 10, 19, model_values=model_values, pdf.plot="S12.Bii")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 11, 20, model_values=model_values, pdf.plot="S12.Biii")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 12, 21, model_values=model_values, pdf.plot="S12.Biv")

  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 9,  17, model_values=model_values, pdf.plot="S12.Ci")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 10, 18, model_values=model_values, pdf.plot="S12.Cii")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 11, 19, model_values=model_values, pdf.plot="S12.Ciii")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 12, 20, model_values=model_values, pdf.plot="S12.Civ")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 13, 21, model_values=model_values, pdf.plot="S12.Cv")

  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 9,  16, model_values=model_values, pdf.plot="S12.Di")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 10, 17, model_values=model_values, pdf.plot="S12.Dii")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 11, 18, model_values=model_values, pdf.plot="S12.Diii")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 12, 19, model_values=model_values, pdf.plot="S12.Div")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 13, 20, model_values=model_values, pdf.plot="S12.Dv")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 14, 21, model_values=model_values, pdf.plot="S12.Dvi")

  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 9,  15, model_values=model_values, pdf.plot="S12.Ei")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 10, 16, model_values=model_values, pdf.plot="S12.Eii")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 11, 17, model_values=model_values, pdf.plot="S12.Eiii")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 12, 18, model_values=model_values, pdf.plot="S12.Eiv")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 13, 19, model_values=model_values, pdf.plot="S12.Ev")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 14, 20, model_values=model_values, pdf.plot="S12.Evi")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 15, 21, model_values=model_values, pdf.plot="S12.Evii")

  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 9,  14, model_values=model_values, pdf.plot="S12.Fi")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 10, 15, model_values=model_values, pdf.plot="S12.Fii")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 11, 16, model_values=model_values, pdf.plot="S12.Fiii")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 12, 17, model_values=model_values, pdf.plot="S12.Fiv")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 13, 18, model_values=model_values, pdf.plot="S12.Fv")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 14, 19, model_values=model_values, pdf.plot="S12.Fvi")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 15, 20, model_values=model_values, pdf.plot="S12.Fvii")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 16, 21, model_values=model_values, pdf.plot="S12.Fviii")

  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 9,  13, model_values=model_values, pdf.plot="S12.Gi")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 10, 14, model_values=model_values, pdf.plot="S12.Gii")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 11, 15, model_values=model_values, pdf.plot="S12.Giii")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 12, 16, model_values=model_values, pdf.plot="S12.Giv")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 13, 17, model_values=model_values, pdf.plot="S12.Gv")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 14, 18, model_values=model_values, pdf.plot="S12.Gvi")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 15, 19, model_values=model_values, pdf.plot="S12.Gvii")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 16, 20, model_values=model_values, pdf.plot="S12.Gviii")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 17, 21, model_values=model_values, pdf.plot="S12.Gix")

  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 9,  12, model_values=model_values, pdf.plot="S12.Hi")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 10, 13, model_values=model_values, pdf.plot="S12.Hii")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 11, 14, model_values=model_values, pdf.plot="S12.Hiii")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 12, 15, model_values=model_values, pdf.plot="S12.Hiv")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 13, 16, model_values=model_values, pdf.plot="S12.Hv")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 14, 17, model_values=model_values, pdf.plot="S12.Hvi")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 15, 18, model_values=model_values, pdf.plot="S12.Hvii")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 16, 19, model_values=model_values, pdf.plot="S12.Hviii")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 17, 20, model_values=model_values, pdf.plot="S12.Hix")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 18, 21, model_values=model_values, pdf.plot="S12.Hx")


# MakeSupplementalFigure13 <- function(model_values=TRUE) {
#   message("Making fig. S13")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 9, 19, model_values=model_values, pdf.plot="S13.Ai")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 10, 20, model_values=model_values, pdf.plot="S13.Aii")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 11, 21, model_values=model_values, pdf.plot="S13.Aiii")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 12, 22, model_values=model_values, pdf.plot="S13.Aiv")

#   PlotSiteMismatches("miR-1", "equil_c_nb", 9,  18, model_values=model_values, pdf.plot="S13.Bi")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 10, 19, model_values=model_values, pdf.plot="S13.Bii")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 11, 20, model_values=model_values, pdf.plot="S13.Biii")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 12, 21, model_values=model_values, pdf.plot="S13.Biv")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 13, 22, model_values=model_values, pdf.plot="S13.Bv")

#   PlotSiteMismatches("miR-1", "equil_c_nb", 9,  17, model_values=model_values, pdf.plot="S13.Ci")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 10, 18, model_values=model_values, pdf.plot="S13.Cii")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 11, 19, model_values=model_values, pdf.plot="S13.Ciii")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 12, 20, model_values=model_values, pdf.plot="S13.Civ")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 13, 21, model_values=model_values, pdf.plot="S13.Cv")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 14, 22, model_values=model_values, pdf.plot="S13.Cvi")

#   PlotSiteMismatches("miR-1", "equil_c_nb", 9,  16, model_values=model_values, pdf.plot="S13.Di")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 10, 17, model_values=model_values, pdf.plot="S13.Dii")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 11, 18, model_values=model_values, pdf.plot="S13.Diii")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 12, 19, model_values=model_values, pdf.plot="S13.Div")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 13, 20, model_values=model_values, pdf.plot="S13.Dv")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 14, 21, model_values=model_values, pdf.plot="S13.Dvi")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 15, 22, model_values=model_values, pdf.plot="S13.Dvii")

#   PlotSiteMismatches("miR-1", "equil_c_nb", 9,  15, model_values=model_values, pdf.plot="S13.Ei")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 10, 16, model_values=model_values, pdf.plot="S13.Eii")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 11, 17, model_values=model_values, pdf.plot="S13.Eiii")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 12, 18, model_values=model_values, pdf.plot="S13.Eiv")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 13, 19, model_values=model_values, pdf.plot="S13.Ev")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 14, 20, model_values=model_values, pdf.plot="S13.Evi")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 15, 21, model_values=model_values, pdf.plot="S13.Evii")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 16, 22, model_values=model_values, pdf.plot="S13.Eviii")

#   PlotSiteMismatches("miR-1", "equil_c_nb", 9,  14, model_values=model_values, pdf.plot="S13.Fi")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 10, 15, model_values=model_values, pdf.plot="S13.Fii")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 11, 16, model_values=model_values, pdf.plot="S13.Fiii")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 12, 17, model_values=model_values, pdf.plot="S13.Fiv")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 13, 18, model_values=model_values, pdf.plot="S13.Fv")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 14, 19, model_values=model_values, pdf.plot="S13.Fvi")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 15, 20, model_values=model_values, pdf.plot="S13.Fvii")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 16, 21, model_values=model_values, pdf.plot="S13.Fviii")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 17, 22, model_values=model_values, pdf.plot="S13.Fix")

#   PlotSiteMismatches("miR-1", "equil_c_nb", 9,  13, model_values=model_values, pdf.plot="S13.Gi")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 10, 14, model_values=model_values, pdf.plot="S13.Gii")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 11, 15, model_values=model_values, pdf.plot="S13.Giii")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 12, 16, model_values=model_values, pdf.plot="S13.Giv")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 13, 17, model_values=model_values, pdf.plot="S13.Gv")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 14, 18, model_values=model_values, pdf.plot="S13.Gvi")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 15, 19, model_values=model_values, pdf.plot="S13.Gvii")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 16, 20, model_values=model_values, pdf.plot="S13.Gviii")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 17, 21, model_values=model_values, pdf.plot="S13.Gix")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 18, 22, model_values=model_values, pdf.plot="S13.Gx")

#   PlotSiteMismatches("miR-1", "equil_c_nb", 9,  12, model_values=model_values, pdf.plot="S13.Hi")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 10, 13, model_values=model_values, pdf.plot="S13.Hii")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 11, 14, model_values=model_values, pdf.plot="S13.Hiii")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 12, 15, model_values=model_values, pdf.plot="S13.Hiv")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 13, 16, model_values=model_values, pdf.plot="S13.Hv")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 14, 17, model_values=model_values, pdf.plot="S13.Hvi")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 15, 18, model_values=model_values, pdf.plot="S13.Hvii")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 16, 19, model_values=model_values, pdf.plot="S13.Hviii")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 17, 20, model_values=model_values, pdf.plot="S13.Hix")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 18, 21, model_values=model_values, pdf.plot="S13.Hx")
#   PlotSiteMismatches("miR-1", "equil_c_nb", 19, 22, model_values=model_values, pdf.plot="S13.Hxi")

#   message("Done fig. S13.")
# }


# MakeSupplementalFigure14 <- function(model_values=TRUE) {
#   message("Making fig. S14")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 9, 19, model_values=model_values, pdf.plot="S14.Ai")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 10, 20, model_values=model_values, pdf.plot="S14.Aii")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 11, 21, model_values=model_values, pdf.plot="S14.Aiii")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 13, 22, model_values=model_values, pdf.plot="S14.Aiv")

#   PlotSiteMismatches("miR-155", "equil_sc_nb", 9,  18, model_values=model_values, pdf.plot="S14.Bi")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 10, 19, model_values=model_values, pdf.plot="S14.Bii")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 11, 20, model_values=model_values, pdf.plot="S14.Biii")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 12, 21, model_values=model_values, pdf.plot="S14.Biv")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 13, 22, model_values=model_values, pdf.plot="S14.Bv")

#   PlotSiteMismatches("miR-155", "equil_sc_nb", 9,  17, model_values=model_values, pdf.plot="S14.Ci")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 10, 18, model_values=model_values, pdf.plot="S14.Cii")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 11, 19, model_values=model_values, pdf.plot="S14.Ciii")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 12, 20, model_values=model_values, pdf.plot="S14.Civ")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 13, 21, model_values=model_values, pdf.plot="S14.Cv")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 14, 22, model_values=model_values, pdf.plot="S14.Cvi")

#   PlotSiteMismatches("miR-155", "equil_sc_nb", 9,  16, model_values=model_values, pdf.plot="S14.Di")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 10, 17, model_values=model_values, pdf.plot="S14.Dii")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 11, 18, model_values=model_values, pdf.plot="S14.Diii")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 12, 19, model_values=model_values, pdf.plot="S14.Div")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 13, 20, model_values=model_values, pdf.plot="S14.Dv")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 14, 21, model_values=model_values, pdf.plot="S14.Dvi")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 15, 22, model_values=model_values, pdf.plot="S14.Dvii")

#   PlotSiteMismatches("miR-155", "equil_sc_nb", 9,  15, model_values=model_values, pdf.plot="S14.Ei")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 10, 16, model_values=model_values, pdf.plot="S14.Eii")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 11, 17, model_values=model_values, pdf.plot="S14.Eiii")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 12, 18, model_values=model_values, pdf.plot="S14.Eiv")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 13, 19, model_values=model_values, pdf.plot="S14.Ev")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 14, 20, model_values=model_values, pdf.plot="S14.Evi")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 15, 21, model_values=model_values, pdf.plot="S14.Evii")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 16, 22, model_values=model_values, pdf.plot="S14.Eviii")

#   PlotSiteMismatches("miR-155", "equil_sc_nb", 9,  14, model_values=model_values, pdf.plot="S14.Fi")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 10, 15, model_values=model_values, pdf.plot="S14.Fii")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 11, 16, model_values=model_values, pdf.plot="S14.Fiii")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 12, 17, model_values=model_values, pdf.plot="S14.Fiv")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 13, 18, model_values=model_values, pdf.plot="S14.Fv")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 14, 19, model_values=model_values, pdf.plot="S14.Fvi")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 15, 20, model_values=model_values, pdf.plot="S14.Fvii")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 16, 21, model_values=model_values, pdf.plot="S14.Fviii")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 17, 22, model_values=model_values, pdf.plot="S14.Fix")

#   PlotSiteMismatches("miR-155", "equil_sc_nb", 9,  13, model_values=model_values, pdf.plot="S14.Gi")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 10, 14, model_values=model_values, pdf.plot="S14.Gii")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 11, 15, model_values=model_values, pdf.plot="S14.Giii")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 12, 16, model_values=model_values, pdf.plot="S14.Giv")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 13, 17, model_values=model_values, pdf.plot="S14.Gv")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 14, 18, model_values=model_values, pdf.plot="S14.Gvi")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 15, 19, model_values=model_values, pdf.plot="S14.Gvii")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 16, 20, model_values=model_values, pdf.plot="S14.Gviii")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 17, 21, model_values=model_values, pdf.plot="S14.Gix")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 18, 22, model_values=model_values, pdf.plot="S14.Gx")

#   PlotSiteMismatches("miR-155", "equil_sc_nb", 9,  12, model_values=model_values, pdf.plot="S14.Hi")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 10, 13, model_values=model_values, pdf.plot="S14.Hii")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 11, 14, model_values=model_values, pdf.plot="S14.Hiii")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 12, 15, model_values=model_values, pdf.plot="S14.Hiv")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 13, 16, model_values=model_values, pdf.plot="S14.Hv")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 14, 17, model_values=model_values, pdf.plot="S14.Hvi")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 15, 18, model_values=model_values, pdf.plot="S14.Hvii")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 16, 19, model_values=model_values, pdf.plot="S14.Hviii")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 17, 20, model_values=model_values, pdf.plot="S14.Hix")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 18, 21, model_values=model_values, pdf.plot="S14.Hx")
#   PlotSiteMismatches("miR-155", "equil_sc_nb", 19, 22, model_values=model_values, pdf.plot="S14.Hxi")

#   message("Done fig. S14.")
# }




# MakeSupplementalFigure15 <- function() {
#   message("Making fig. S15")
#   CompareSiteMismatches("let-7a-21nt", "equil_c2_nb",  9, 17, 10, 18, model_values=FALSE, pdf.plot="S15.Ai")
#   CompareSiteMismatches("let-7a-21nt", "equil_c2_nb", 10, 18, 11, 19, model_values=FALSE, pdf.plot="S15.Aii")
#   CompareSiteMismatches("let-7a-21nt", "equil_c2_nb", 11, 19, 12, 20, model_values=FALSE, pdf.plot="S15.Aiii")
#   # CompareSiteMismatches("let-7a-21nt", "equil_c2_nb", 12, 20, 13, 21, model_values=FALSE, pdf.plot="S15.Aiv")
#   # CompareSiteMismatches("let-7a-21nt", "equil_c2_nb", 13, 21, 14, 20, model_values=FALSE, pdf.plot="S15.Av")
#   # CompareSiteMismatches("let-7a-21nt", "equil_c2_nb", 14, 20, 15, 21, model_values=FALSE, pdf.plot="S15.Avi")

#   CompareSiteMismatches("let-7a-21nt", "equil_c2_nb",  9, 17, 10, 18, model_values=TRUE, pdf.plot="S15.Bi")
#   CompareSiteMismatches("let-7a-21nt", "equil_c2_nb", 10, 18, 11, 19, model_values=TRUE, pdf.plot="S15.Bii")
#   CompareSiteMismatches("let-7a-21nt", "equil_c2_nb", 11, 19, 12, 20, model_values=TRUE, pdf.plot="S15.Biii")
#   # CompareSiteMismatches("let-7a-21nt", "equil_c2_nb", 12, 18, 13, 19, model_values=TRUE, pdf.plot="S15.Biv")
#   # CompareSiteMismatches("let-7a-21nt", "equil_c2_nb", 13, 19, 14, 20, model_values=TRUE, pdf.plot="S15.Bv")
#   # CompareSiteMismatches("let-7a-21nt", "equil_c2_nb", 14, 20, 15, 21, model_values=TRUE, pdf.plot="S15.Bvi")


#   CompareSiteMismatches("miR-1", "equil_c_nb",  9, 17, 10, 18, model_values=FALSE, pdf.plot="S15.Ci")
#   CompareSiteMismatches("miR-1", "equil_c_nb", 10, 18, 11, 19, model_values=FALSE, pdf.plot="S15.Cii")
#   CompareSiteMismatches("miR-1", "equil_c_nb", 11, 19, 12, 20, model_values=FALSE, pdf.plot="S15.Ciii")
#   # CompareSiteMismatches("miR-1", "equil_c_nb", 12, 18, 13, 19, model_values=FALSE, pdf.plot="S15.Civ")
#   # CompareSiteMismatches("miR-1", "equil_c_nb", 13, 19, 14, 20, model_values=FALSE, pdf.plot="S15.Cv")
#   # CompareSiteMismatches("miR-1", "equil_c_nb", 14, 20, 15, 21, model_values=FALSE, pdf.plot="S15.Cvi")
#   # CompareSiteMismatches("miR-1", "equil_c_nb", 15, 21, 16, 22, model_values=FALSE, pdf.plot="S15.Cvii")

#   CompareSiteMismatches("miR-1", "equil_c_nb",  9, 17, 10, 18, model_values=TRUE, pdf.plot="S15.Di")
#   CompareSiteMismatches("miR-1", "equil_c_nb", 10, 18, 11, 19, model_values=TRUE, pdf.plot="S15.Dii")
#   CompareSiteMismatches("miR-1", "equil_c_nb", 11, 19, 12, 20, model_values=TRUE, pdf.plot="S15.Diii")
#   # CompareSiteMismatches("miR-1", "equil_c_nb", 12, 18, 13, 19, model_values=TRUE, pdf.plot="S15.Div")
#   # CompareSiteMismatches("miR-1", "equil_c_nb", 13, 19, 14, 20, model_values=TRUE, pdf.plot="S15.Dv")
#   # CompareSiteMismatches("miR-1", "equil_c_nb", 14, 20, 15, 21, model_values=TRUE, pdf.plot="S15.Dvi")
#   # CompareSiteMismatches("miR-1", "equil_c_nb", 15, 21, 16, 22, model_values=TRUE, pdf.plot="S15.Dvii")


#   CompareSiteMismatches("miR-155", "equil_sc_nb",  9, 17, 10, 18, model_values=FALSE, pdf.plot="S15.Ei")
#   CompareSiteMismatches("miR-155", "equil_sc_nb", 10, 18, 11, 19, model_values=FALSE, pdf.plot="S15.Eii")
#   CompareSiteMismatches("miR-155", "equil_sc_nb", 11, 19, 12, 20, model_values=FALSE, pdf.plot="S15.Eiii")
#   # CompareSiteMismatches("miR-155", "equil_sc_nb", 12, 18, 13, 19, model_values=FALSE, pdf.plot="S15.Eiv")
#   # CompareSiteMismatches("miR-155", "equil_sc_nb", 13, 19, 14, 20, model_values=FALSE, pdf.plot="S15.Ev")
#   # CompareSiteMismatches("miR-155", "equil_sc_nb", 14, 20, 15, 21, model_values=FALSE, pdf.plot="S15.Evi")
#   # CompareSiteMismatches("miR-155", "equil_sc_nb", 15, 21, 16, 22, model_values=FALSE, pdf.plot="S15.Evii")
#   # CompareSiteMismatches("miR-155", "equil_sc_nb", 16, 22, 17, 23, model_values=FALSE, pdf.plot="S15.Eviii")

#   CompareSiteMismatches("miR-155", "equil_sc_nb",  9, 17, 10, 18, model_values=TRUE, pdf.plot="S15.Fi")
#   CompareSiteMismatches("miR-155", "equil_sc_nb", 10, 18, 11, 19, model_values=TRUE, pdf.plot="S15.Fii")
#   CompareSiteMismatches("miR-155", "equil_sc_nb", 11, 19, 12, 20, model_values=TRUE, pdf.plot="S15.Fiii")
#   # CompareSiteMismatches("miR-155", "equil_sc_nb", 12, 18, 13, 19, model_values=TRUE, pdf.plot="S15.Fiv")
#   # CompareSiteMismatches("miR-155", "equil_sc_nb", 13, 19, 14, 20, model_values=TRUE, pdf.plot="S15.Fv")
#   # CompareSiteMismatches("miR-155", "equil_sc_nb", 14, 20, 15, 21, model_values=TRUE, pdf.plot="S15.Fvi")
#   # CompareSiteMismatches("miR-155", "equil_sc_nb", 15, 21, 16, 22, model_values=TRUE, pdf.plot="S15.Fvii")
#   # CompareSiteMismatches("miR-155", "equil_sc_nb", 16, 22, 17, 23, model_values=TRUE, pdf.plot="S15.Fviii")

#   message("Done fig. S15.")
# }

# MakeSupplementalFigure14 <- function() {
#   message("Making fig. S14")

#   PlotMismatchKdsAgainstDeltaG("let-7a-21nt", "equil_c2_nb", 10, 20, model_values=FALSE, just_complement=TRUE, pdf.plot="S14.Aii")
#   PlotMismatchKdsAgainstDeltaG("let-7a-21nt", "equil_c2_nb", 10, 20, model_values=TRUE, pdf.plot="S14.Aiii")
#   PlotMismatchKdsAgainstDeltaG("let-7a-21nt", "equil_c2_nb", 10, 20, model_values=TRUE, just_complement=TRUE, pdf.plot="S14.Aiv")

#   PlotMismatchKdsAgainstDeltaG("let-7a-21nt", "equil_c2_nb", 10, 19, model_values=FALSE, pdf.plot="S14.Bi")
#   PlotMismatchKdsAgainstDeltaG("let-7a-21nt", "equil_c2_nb", 10, 19, model_values=FALSE, just_complement=TRUE, pdf.plot="S14.Bii")
#   PlotMismatchKdsAgainstDeltaG("let-7a-21nt", "equil_c2_nb", 10, 19, model_values=TRUE, pdf.plot="S14.Biii")
#   PlotMismatchKdsAgainstDeltaG("let-7a-21nt", "equil_c2_nb", 10, 19, model_values=TRUE, just_complement=TRUE, pdf.plot="S14.Biv")

#   PlotMismatchKdsAgainstDeltaG("let-7a-21nt", "equil_c2_nb", 11, 19, model_values=FALSE, pdf.plot="S14.Ci")
#   PlotMismatchKdsAgainstDeltaG("let-7a-21nt", "equil_c2_nb", 11, 19, model_values=FALSE, just_complement=TRUE, pdf.plot="S14.Cii")
#   PlotMismatchKdsAgainstDeltaG("let-7a-21nt", "equil_c2_nb", 11, 19, model_values=TRUE, pdf.plot="S14.Ciii")
#   PlotMismatchKdsAgainstDeltaG("let-7a-21nt", "equil_c2_nb", 11, 19, model_values=TRUE, just_complement=TRUE, pdf.plot="S14.Civ")

#   PlotMismatchKdsAgainstDeltaG("let-7a-21nt", "equil_c2_nb", 12, 20, model_values=FALSE, pdf.plot="S14.Di")
#   PlotMismatchKdsAgainstDeltaG("let-7a-21nt", "equil_c2_nb", 12, 20, model_values=FALSE, just_complement=TRUE, pdf.plot="S14.Dii")
#   PlotMismatchKdsAgainstDeltaG("let-7a-21nt", "equil_c2_nb", 12, 20, model_values=TRUE, pdf.plot="S14.Diii")
#   PlotMismatchKdsAgainstDeltaG("let-7a-21nt", "equil_c2_nb", 12, 20, model_values=TRUE, just_complement=TRUE, pdf.plot="S14.Div")

#   PlotMismatchKdsAgainstDeltaG("let-7a-21nt", "equil_c2_nb", 13, 21, model_values=FALSE, pdf.plot="S14.Ei")
#   PlotMismatchKdsAgainstDeltaG("let-7a-21nt", "equil_c2_nb", 13, 21, model_values=FALSE, just_complement=TRUE, pdf.plot="S14.Eii")
#   PlotMismatchKdsAgainstDeltaG("let-7a-21nt", "equil_c2_nb", 13, 21, model_values=TRUE, pdf.plot="S14.Eiii")
#   PlotMismatchKdsAgainstDeltaG("let-7a-21nt", "equil_c2_nb", 13, 21, model_values=TRUE, just_complement=TRUE, pdf.plot="S14.Eiv")


#   message("Done fig. S14.")
# }

# MakeSupplementalFigure15 <- function() {
#   message("Making fig. S15")

#   PlotMismatchKdsAgainstDeltaG("miR-1", "equil_c_nb", 11, 21, model_values=FALSE, pdf.plot="S15.Ai")
#   PlotMismatchKdsAgainstDeltaG("miR-1", "equil_c_nb", 11, 21, model_values=FALSE, just_complement=TRUE, pdf.plot="S15.Aii")
#   PlotMismatchKdsAgainstDeltaG("miR-1", "equil_c_nb", 11, 21, model_values=TRUE, pdf.plot="S15.Aiii")
#   PlotMismatchKdsAgainstDeltaG("miR-1", "equil_c_nb", 11, 21, model_values=TRUE, just_complement=TRUE, pdf.plot="S15.Aiv")

#   PlotMismatchKdsAgainstDeltaG("miR-1", "equil_c_nb", 11, 20, model_values=FALSE, pdf.plot="S15.Bi")
#   PlotMismatchKdsAgainstDeltaG("miR-1", "equil_c_nb", 11, 20, model_values=FALSE, just_complement=TRUE, pdf.plot="S15.Bii")
#   PlotMismatchKdsAgainstDeltaG("miR-1", "equil_c_nb", 11, 20, model_values=TRUE, pdf.plot="S15.Biii")
#   PlotMismatchKdsAgainstDeltaG("miR-1", "equil_c_nb", 11, 20, model_values=TRUE, just_complement=TRUE, pdf.plot="S15.Biv")

#   PlotMismatchKdsAgainstDeltaG("miR-1", "equil_c_nb", 12, 20, model_values=FALSE, pdf.plot="S15.Ci")
#   PlotMismatchKdsAgainstDeltaG("miR-1", "equil_c_nb", 12, 20, model_values=FALSE, just_complement=TRUE, pdf.plot="S15.Cii")
#   PlotMismatchKdsAgainstDeltaG("miR-1", "equil_c_nb", 12, 20, model_values=TRUE, pdf.plot="S15.Ciii")
#   PlotMismatchKdsAgainstDeltaG("miR-1", "equil_c_nb", 12, 20, model_values=TRUE, just_complement=TRUE, pdf.plot="S15.Civ")

#   PlotMismatchKdsAgainstDeltaG("miR-1", "equil_c_nb", 13, 21, model_values=FALSE, pdf.plot="S15.Di")
#   PlotMismatchKdsAgainstDeltaG("miR-1", "equil_c_nb", 13, 21, model_values=FALSE, just_complement=TRUE, pdf.plot="S15.Dii")
#   PlotMismatchKdsAgainstDeltaG("miR-1", "equil_c_nb", 13, 21, model_values=TRUE, pdf.plot="S15.Diii")
#   PlotMismatchKdsAgainstDeltaG("miR-1", "equil_c_nb", 13, 21, model_values=TRUE, just_complement=TRUE, pdf.plot="S15.Div")

#   PlotMismatchKdsAgainstDeltaG("miR-1", "equil_c_nb", 14, 22, model_values=FALSE, pdf.plot="S15.Ei")
#   PlotMismatchKdsAgainstDeltaG("miR-1", "equil_c_nb", 14, 22, model_values=FALSE, just_complement=TRUE, pdf.plot="S15.Eii")
#   PlotMismatchKdsAgainstDeltaG("miR-1", "equil_c_nb", 14, 22, model_values=TRUE, pdf.plot="S15.Eiii")
#   PlotMismatchKdsAgainstDeltaG("miR-1", "equil_c_nb", 14, 22, model_values=TRUE, just_complement=TRUE, pdf.plot="S15.Eiv")

#   message("Done fig. S15.")
# }


# MakeSupplementalFigure16 <- function() {
#   message("Making fig. S16")

#   PlotMismatchKdsAgainstDeltaG("miR-155", "equil_sc_nb", 12, 22, model_values=FALSE, pdf.plot="S16.Ai")
#   PlotMismatchKdsAgainstDeltaG("miR-155", "equil_sc_nb", 12, 22, model_values=FALSE, just_complement=TRUE, pdf.plot="S16.Aii")
#   PlotMismatchKdsAgainstDeltaG("miR-155", "equil_sc_nb", 12, 22, model_values=TRUE, pdf.plot="S16.Aiii")
#   PlotMismatchKdsAgainstDeltaG("miR-155", "equil_sc_nb", 12, 22, model_values=TRUE, just_complement=TRUE, pdf.plot="S16.Aiv")

#   PlotMismatchKdsAgainstDeltaG("miR-155", "equil_sc_nb", 12, 21, model_values=FALSE, pdf.plot="S16.Bi")
#   PlotMismatchKdsAgainstDeltaG("miR-155", "equil_sc_nb", 12, 21, model_values=FALSE, just_complement=TRUE, pdf.plot="S16.Bii")
#   PlotMismatchKdsAgainstDeltaG("miR-155", "equil_sc_nb", 12, 21, model_values=TRUE, pdf.plot="S16.Biii")
#   PlotMismatchKdsAgainstDeltaG("miR-155", "equil_sc_nb", 12, 21, model_values=TRUE, just_complement=TRUE, pdf.plot="S16.Biv")

#   PlotMismatchKdsAgainstDeltaG("miR-155", "equil_sc_nb", 13, 21, model_values=FALSE, pdf.plot="S16.Ci")
#   PlotMismatchKdsAgainstDeltaG("miR-155", "equil_sc_nb", 13, 21, model_values=FALSE, just_complement=TRUE, pdf.plot="S16.Cii")
#   PlotMismatchKdsAgainstDeltaG("miR-155", "equil_sc_nb", 13, 21, model_values=TRUE, pdf.plot="S16.Ciii")
#   PlotMismatchKdsAgainstDeltaG("miR-155", "equil_sc_nb", 13, 21, model_values=TRUE, just_complement=TRUE, pdf.plot="S16.Civ")

#   PlotMismatchKdsAgainstDeltaG("miR-155", "equil_sc_nb", 14, 22, model_values=FALSE, pdf.plot="S16.Di")
#   PlotMismatchKdsAgainstDeltaG("miR-155", "equil_sc_nb", 14, 22, model_values=FALSE, just_complement=TRUE, pdf.plot="S16.Dii")
#   PlotMismatchKdsAgainstDeltaG("miR-155", "equil_sc_nb", 14, 22, model_values=TRUE, pdf.plot="S16.Diii")
#   PlotMismatchKdsAgainstDeltaG("miR-155", "equil_sc_nb", 14, 22, model_values=TRUE, just_complement=TRUE, pdf.plot="S16.Div")

#   PlotMismatchKdsAgainstDeltaG("miR-155", "equil_sc_nb", 15, 23, model_values=FALSE, pdf.plot="S16.Ei")
#   PlotMismatchKdsAgainstDeltaG("miR-155", "equil_sc_nb", 15, 23, model_values=FALSE, just_complement=TRUE, pdf.plot="S16.Eii")
#   PlotMismatchKdsAgainstDeltaG("miR-155", "equil_sc_nb", 15, 23, model_values=TRUE, pdf.plot="S16.Eiii")
#   PlotMismatchKdsAgainstDeltaG("miR-155", "equil_sc_nb", 15, 23, model_values=TRUE, just_complement=TRUE, pdf.plot="S16.Eiv")

#   message("Done fig. S16.")
# }







MakeSupplementalFigure17 <- function(plot_pred=FALSE, plot_emp=FALSE, ratio=FALSE, frac_dG=TRUE, frac_ddG=FALSE) {
  message("Making fig. S17")
  PlotThrPMismatchDeltaGResidual("let-7a-21nt", "equil_c2_nb", 11, plot_pred=plot_pred, plot_emp=plot_emp, ratio=ratio, frac_dG=frac_dG, frac_ddG=frac_ddG, pdf.plot="S17.Ai")
  PlotThrPMismatchDeltaGResidual("miR-1", "equil_c_nb", 11, plot_pred=plot_pred, plot_emp=plot_emp, ratio=ratio, frac_dG=frac_dG, frac_ddG=frac_ddG, pdf.plot="S17.Aii")
  PlotThrPMismatchDeltaGResidual("miR-155", "equil_sc_nb", 11, plot_pred=plot_pred, plot_emp=plot_emp, ratio=ratio, frac_dG=frac_dG, frac_ddG=frac_ddG, pdf.plot="S17.Aiii")

  PlotThrPMismatchDeltaGResidual("let-7a-21nt", "equil_c2_nb", 11, model_values=TRUE,plot_pred=plot_pred, plot_emp=plot_emp, ratio=ratio, frac_dG=frac_dG, frac_ddG=frac_ddG, pdf.plot="S17.Bi")
  PlotThrPMismatchDeltaGResidual("miR-1", "equil_c_nb", 11, model_values=TRUE, plot_pred=plot_pred, plot_emp=plot_emp, ratio=ratio, frac_dG=frac_dG, frac_ddG=frac_ddG, pdf.plot="S17.Bii")
  PlotThrPMismatchDeltaGResidual("miR-155", "equil_sc_nb", 11, model_values=TRUE, plot_pred=plot_pred, plot_emp=plot_emp, ratio=ratio, frac_dG=frac_dG, frac_ddG=frac_ddG, pdf.plot="S17.Biii")


  PlotThrPMismatchDeltaGResidual("let-7a-21nt", "equil_c2_nb", 8, plot_pred=plot_pred, plot_emp=plot_emp, ratio=ratio, frac_dG=frac_dG, frac_ddG=frac_ddG, pdf.plot="S17.Ci")
  PlotThrPMismatchDeltaGResidual("miR-1", "equil_c_nb", 8, plot_pred=plot_pred, plot_emp=plot_emp, ratio=ratio, frac_dG=frac_dG, frac_ddG=frac_ddG, pdf.plot="S17.Cii")
  PlotThrPMismatchDeltaGResidual("miR-155", "equil_sc_nb", 8, plot_pred=plot_pred, plot_emp=plot_emp, ratio=ratio, frac_dG=frac_dG, frac_ddG=frac_ddG, pdf.plot="S17.Ciii")


  PlotThrPMismatchDeltaGResidual("let-7a-21nt", "equil_c2_nb", 8, model_values=TRUE,plot_pred=plot_pred, plot_emp=plot_emp, ratio=ratio, frac_dG=frac_dG, frac_ddG=frac_ddG, pdf.plot="S17.Di")
  PlotThrPMismatchDeltaGResidual("miR-1", "equil_c_nb", 8, model_values=TRUE, plot_pred=plot_pred, plot_emp=plot_emp, ratio=ratio, frac_dG=frac_dG, frac_ddG=frac_ddG, pdf.plot="S17.Dii")
  PlotThrPMismatchDeltaGResidual("miR-155", "equil_sc_nb", 8, model_values=TRUE, plot_pred=plot_pred, plot_emp=plot_emp, ratio=ratio, frac_dG=frac_dG, frac_ddG=frac_ddG, pdf.plot="S17.Diii")


  PlotThrPMismatchDeltaGResidual("let-7a", "equilibrium", 8, plot_pred=plot_pred, plot_emp=plot_emp, ratio=ratio, frac_dG=frac_dG, frac_ddG=frac_ddG, pdf.plot="S17.Ei")
  PlotThrPMismatchDeltaGResidual("miR-1", "equilibrium", 8, plot_pred=plot_pred, plot_emp=plot_emp, ratio=ratio, frac_dG=frac_dG, frac_ddG=frac_ddG, pdf.plot="S17.Eii")
  PlotThrPMismatchDeltaGResidual("miR-155", "equilibrium", 8, plot_pred=plot_pred, plot_emp=plot_emp, ratio=ratio, frac_dG=frac_dG, frac_ddG=frac_ddG, pdf.plot="S17.Eiii")

  PlotThrPMismatchDeltaGResidual("let-7a", "equilibrium", 8, model_values=TRUE,plot_pred=plot_pred, plot_emp=plot_emp, ratio=ratio, frac_dG=frac_dG, frac_ddG=frac_ddG, pdf.plot="S17.Fi")
  PlotThrPMismatchDeltaGResidual("miR-1", "equilibrium", 8, model_values=TRUE, plot_pred=plot_pred, plot_emp=plot_emp, ratio=ratio, frac_dG=frac_dG, frac_ddG=frac_ddG, pdf.plot="S17.Fii")
  PlotThrPMismatchDeltaGResidual("miR-155", "equilibrium", 8, model_values=TRUE, plot_pred=plot_pred, plot_emp=plot_emp, ratio=ratio, frac_dG=frac_dG, frac_ddG=frac_ddG, pdf.plot="S17.Fiii")


  PlotThrPMismatchDeltaGResidual("miR-124", "equilibrium", 8, plot_pred=plot_pred, plot_emp=plot_emp, ratio=ratio, frac_dG=frac_dG, frac_ddG=frac_ddG, pdf.plot="S17.Gi")
  PlotThrPMismatchDeltaGResidual("lsy-6", "equilibrium", 8, plot_pred=plot_pred, plot_emp=plot_emp, ratio=ratio, frac_dG=frac_dG, frac_ddG=frac_ddG, pdf.plot="S17.Gii")
  PlotThrPMismatchDeltaGResidual("miR-7-23nt", "equilibrium2_nb", 8, plot_pred=plot_pred, plot_emp=plot_emp, ratio=ratio, frac_dG=frac_dG, frac_ddG=frac_ddG, pdf.plot="S17.Giii")

  PlotThrPMismatchDeltaGResidual("miR-124", "equilibrium", 8, model_values=TRUE,plot_pred=plot_pred, plot_emp=plot_emp, ratio=ratio, frac_dG=frac_dG, frac_ddG=frac_ddG, pdf.plot="S17.Hi")
  PlotThrPMismatchDeltaGResidual("lsy-6", "equilibrium", 8, model_values=TRUE, plot_pred=plot_pred, plot_emp=plot_emp, ratio=ratio, frac_dG=frac_dG, frac_ddG=frac_ddG, pdf.plot="S17.Hii")
  PlotThrPMismatchDeltaGResidual("miR-7-23nt", "equilibrium2_nb", 8, model_values=TRUE, plot_pred=plot_pred, plot_emp=plot_emp, ratio=ratio, frac_dG=frac_dG, frac_ddG=frac_ddG, pdf.plot="S17.Hiii")



  message("Done fig. S17.")

}



MakeSupplementalFigure18 <- function() {
  message("Making fig. S18")
  # PlotAverageMismatches("let-7a-21nt", "equil_c2_nb", 11, pdf.plot="S18.Ai")
  # PlotAverageMismatches("miR-1", "equil_c_nb", 11, pdf.plot="S18.Aii")
  # PlotAverageMismatches("miR-155", "equil_sc_nb", 11, pdf.plot="S18.Aiii")

  # PlotAverageMismatches("let-7a-21nt", "equil_c2_nb", 11, model_values=TRUE, pdf.plot="S18.Bi")
  # PlotAverageMismatches("miR-1", "equil_c_nb", 11, model_values=TRUE, pdf.plot="S18.Bii")
  # PlotAverageMismatches("miR-155", "equil_sc_nb", 11, model_values=TRUE, pdf.plot="S18.Biii")


  # PlotAverageMismatches("let-7a-21nt", "equil_c2_nb", 8, pdf.plot="S18.Ci")
  # PlotAverageMismatches("miR-1", "equil_c_nb", 8, pdf.plot="S18.Cii")
  # PlotAverageMismatches("miR-155", "equil_sc_nb", 8, pdf.plot="S18.Ciii")


  # PlotAverageMismatches("let-7a-21nt", "equil_c2_nb", 8, model_values=TRUE, pdf.plot="S18.Di")
  # PlotAverageMismatches("miR-1", "equil_c_nb", 8, model_values=TRUE, pdf.plot="S18.Dii")
  # PlotAverageMismatches("miR-155", "equil_sc_nb", 8, model_values=TRUE, pdf.plot="S18.Diii")


  PlotAverageMismatches("let-7a", "equilibrium", 8, pdf.plot="S18.Ei")
  PlotAverageMismatches("miR-1", "equilibrium", 8, pdf.plot="S18.Eii")
  PlotAverageMismatches("miR-155", "equilibrium", 8, pdf.plot="S18.Eiii")

  # PlotAverageMismatches("let-7a", "equilibrium", 8, model_values=TRUE, pdf.plot="S18.Fi")
  # PlotAverageMismatches("miR-1", "equilibrium", 8, model_values=TRUE, pdf.plot="S18.Fii")
  # PlotAverageMismatches("miR-155", "equilibrium", 8, model_values=TRUE, pdf.plot="S18.Fiii")


  PlotAverageMismatches("miR-124", "equilibrium", 8, pdf.plot="S18.Gi")
  PlotAverageMismatches("lsy-6", "equilibrium", 8, pdf.plot="S18.Gii")
  PlotAverageMismatches("miR-7-23nt", "equilibrium2_nb", 8, pdf.plot="S18.Giii")

  # PlotAverageMismatches("miR-124", "equilibrium", 8, model_values=TRUE, pdf.plot="S18.Hi")
  # PlotAverageMismatches("lsy-6", "equilibrium", 8, model_values=TRUE, pdf.plot="S18.Hii")
  # PlotAverageMismatches("miR-7-23nt", "equilibrium2_nb", 8, model_values=TRUE, pdf.plot="S18.Hiii")



  message("Done fig. S18.")

}





MakeSupplementalFigure20 <- function() {
  message("Making Fig. 20")

  # PlotTerminalMismatchAgainstBulge("let-7a-21nt", "equil_c2_nb", 10, pdf.plot="S20.Ai")
  # PlotTerminalMismatchAgainstBulge("let-7a-21nt", "equil_c2_nb", 10, threep=FALSE, pdf.plot="S20.Aii")
  # PlotTerminalMismatchAgainstBulge("let-7a-21nt", "equil_c2_nb", 10, bothp=TRUE, pdf.plot="S20.Aiii")

  # PlotTerminalMismatchAgainstBulge("miR-1", "equil_c_nb", 10, pdf.plot="S20.Bi")
  # PlotTerminalMismatchAgainstBulge("miR-1", "equil_c_nb", 10, threep=FALSE, pdf.plot="S20.Bii")
  # PlotTerminalMismatchAgainstBulge("miR-1", "equil_c_nb", 10, bothp=TRUE, pdf.plot="S20.Biii")

  # PlotTerminalMismatchAgainstBulge("miR-155", "equil_sc_nb", 10, pdf.plot="S20.Ci")
  # PlotTerminalMismatchAgainstBulge("miR-155", "equil_sc_nb", 10, threep=FALSE, pdf.plot="S20.Cii")
  # PlotTerminalMismatchAgainstBulge("miR-155", "equil_sc_nb", 10, bothp=TRUE, pdf.plot="S20.Ciii")

  PlotAllPositionalMismatches(8, rand_data=TRUE, pdf.plot="S20.E")



  message("Done Fig. S20")
}

MakeSupplementalFigure5 <- function(corrected_kds=TRUE) {
  message("Making fig. S5")
  mirnas <- c("miR-1", "miR-155")
  experiments <- c("equil_c_nb", "equil_sc_nb")
  letter <- "A"
  for (i_m in 1:2) {
    mirna <- mirnas[i_m]
    experiment <- experiments[i_m]
    for (i in 1:11) {
      offset <- i*2 - 6
      if (i == 1) mirna_label <- TRUE
      else        mirna_label <- FALSE
      if (i == 11) key <- TRUE
      else         key <- FALSE
      SubfunctionCall(PlotPairingMatrix,
                      pdf.plot=sprintf("S5.%s%s", letter, kRomanNumerals[i]))
    }
    letter <- kNextLetter[letter]
  }
  message("Done fig. S5")
}


MakeSupplementalFigure2 <- function() {
  message("Making Fig. S2")
  PlotBestThreePrimeSite(experiment="equil_c_nb", nbomitc=TRUE, height=4,
                         pdf.plot="S2.A")
  PlotPositionalCanonicalSites(pdf.plot="S2.B")
  PlotPositionalCanonicalSites(pdf.plot="S2.C", mirna="miR-1",
                               experiment="equil_c_nb")
  PlotPositionalCanonicalSites(pdf.plot="S2.D", mirna="miR-155",
                               experiment="equil_sc_nb")
  message("Done Fig. S2")
}


PlotAllPositionalMismatchesOld <- function(
  len, n_constant=3, model_values=FALSE, rand_data=FALSE, 
  new=TRUE, kd_fc=TRUE, win_average=3, corrected_kds=TRUE, titles=FALSE,
  frac_dG=FALSE, plot_kd_fc=FALSE, xpos=20, ypos=20, height=3, width=4.5, pdf.plot=FALSE
) {
  # Load the collapsed data (for the mismatch-only).
  df_use <<- SubfunctionCall(MakeAllLibraryMismatches)
  print(df_use)
  df_use$nuc <- gsub("U", replacement="T", as.character(df_use$nuc))
  SubfunctionCall(FigureSaveFile2)
  par(mar=c(3, 3, 1, 1))
  xmin <- 9
  if (rand_data) {
    xmax <- (21 - 9 + 1)*7 + 9 + 2
  } else {
    xmax <- (21 - 9 + 1)*4 + 9 + 2
  }
  ymin <- 0
  if (frac_dG) {
    ymax <- 1
  } else {
    ymax <- 2.5  
  }
  BlankPlot()
  x_l_global <- 9
  x_line_global <- c()
  y_line_global <- c()
  for (pos_i in 9:21) {
    inds_use <- which(df_use$pos == pos_i)
    print(df_use[inds_use, ])
    if (rand_data) {
      x_l_i <- 0:5
    } else {
      x_l_i <- 0:2
    }
    x_r_i <- x_l_i + 1
    vals <- df_use$dG_obs[inds_use]
    cols <- kNucleotideColors[df_use$nuc[inds_use]]
    rect(xleft=x_l_i + x_l_global, ybottom=0, xright=x_r_i + x_l_global,
         ytop=vals, col=cols, border=NA, xpd=NA)

    x_line_global <- c(x_line_global, x_l_global + (x_l_i[1] + x_r_i[length(x_r_i)])/2)
    y_line_global <- c(y_line_global, mean(vals))
    x_l_global <- x_l_global + length(x_l_i) + 1

  }
  lines(x_line_global, y_line_global, lty=2)
  AddLinearAxis(2, 1, 0.5, label=expression("Observed"~Delta*Delta*italic(G)))
  if (rand_data) {
    text(x=(0:12)*7 + 10.5, y=-0.2*ymax/2.5, labels=as.character(9:21), xpd=NA)
  } else {
    text(x=(0:12)*4 + 10.5, y=-0.2*ymax/2.5, labels=as.character(9:21), xpd=NA)
  }
  text(x=(xmin + xmax)/2, y=-0.50*ymax/2.5, labels="miRNA nucleotide", xpd=NA)
  ## Add the length of the sites. ##############################################
  xy <- GetPlotFractionalCoords(1, 1.1)
  text(x=xy[1], y=xy[2], labels=sprintf("%s-nt sites", len), xpd=NA,
       adj=c(1, 1))
  ## Add in the legend containing the nucleotides. #############################
  xy <- GetPlotFractionalCoords(0, 1.177)
  Legend(xy, col=kNucleotideColors, legend=ConvertTtoU(kNucs),
         legend_pch_use=15, ncol=4, x.intersp=0.5)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}



MakeSupplementalFigure14 <- function() {

  MakeLoopPSAM(13, 6, 0, titles2=TRUE, pdf.plot="S15.Ai")
  MakeLoopPSAM(13, 6, 1, pdf.plot="S15.Aii")
  MakeLoopPSAM(13, 6, 2, pdf.plot="S15.Aiii")
  MakeLoopPSAM(13, 6, 3, pdf.plot="S15.Aiv")
  MakeLoopPSAM(13, 6, 4, pdf.plot="S15.Av")
  MakeLoopPSAM(13, 6, 5, titles=TRUE, pdf.plot="S15.Avi")

  MakeLoopPSAM(13, 7, 0, titles3=TRUE, pdf.plot="S15.Avii")
  MakeLoopPSAM(13, 7, 1, pdf.plot="S15.Aviii")
  MakeLoopPSAM(13, 7, 2, pdf.plot="S15.Aix")
  MakeLoopPSAM(13, 7, 3, pdf.plot="S15.Ax")
  MakeLoopPSAM(13, 7, 4, pdf.plot="S15.Axi")
  MakeLoopPSAM(13, 7, 5, titles=TRUE, pdf.plot="S15.Axii")

  MakeLoopPSAM(13, 8, 0, titles2=TRUE, pdf.plot="S15.Axiii")
  MakeLoopPSAM(13, 8, 1, pdf.plot="S15.Axiv")
  MakeLoopPSAM(13, 8, 2, pdf.plot="S15.Axv")
  MakeLoopPSAM(13, 8, 3, pdf.plot="S15.Axvi")
  MakeLoopPSAM(13, 8, 4, pdf.plot="S15.Axvii")
  MakeLoopPSAM(13, 8, 5, titles=TRUE, pdf.plot="S15.Axviii")


  MakeLoopPSAM(12, 6, 0, titles2=TRUE, pdf.plot="S15.Bi")
  MakeLoopPSAM(12, 6, 1, pdf.plot="S15.Bii")
  MakeLoopPSAM(12, 6, 2, pdf.plot="S15.Biii")
  MakeLoopPSAM(12, 6, 3, pdf.plot="S15.Biv")
  MakeLoopPSAM(12, 6, 4, pdf.plot="S15.Bv")
  MakeLoopPSAM(12, 6, 5, titles=TRUE, pdf.plot="S15.Bvi")

  MakeLoopPSAM(12, 7, 0, titles3=TRUE, pdf.plot="S15.Bvii")
  MakeLoopPSAM(12, 7, 1, pdf.plot="S15.Bviii")
  MakeLoopPSAM(12, 7, 2, pdf.plot="S15.Bix")
  MakeLoopPSAM(12, 7, 3, pdf.plot="S15.Bx")
  MakeLoopPSAM(12, 7, 4, pdf.plot="S15.Bxi")
  MakeLoopPSAM(12, 7, 5, titles=TRUE, pdf.plot="S15.Bxii")

  MakeLoopPSAM(12, 8, 0, titles2=TRUE, pdf.plot="S15.Bxiii")
  MakeLoopPSAM(12, 8, 1, pdf.plot="S15.Bxiv")
  MakeLoopPSAM(12, 8, 2, pdf.plot="S15.Bxv")
  MakeLoopPSAM(12, 8, 3, pdf.plot="S15.Bxvi")
  MakeLoopPSAM(12, 8, 4, pdf.plot="S15.Bxvii")
  MakeLoopPSAM(12, 8, 5, titles=TRUE, pdf.plot="S15.Bxviii")


  MakeLoopPSAM(11, 6, 0, titles2=TRUE, pdf.plot="S15.Ci")
  MakeLoopPSAM(11, 6, 1, pdf.plot="S15.Cii")
  MakeLoopPSAM(11, 6, 2, pdf.plot="S15.Ciii")
  MakeLoopPSAM(11, 6, 3, pdf.plot="S15.Civ")
  MakeLoopPSAM(11, 6, 4, pdf.plot="S15.Cv")
  MakeLoopPSAM(11, 6, 5, titles=TRUE, pdf.plot="S15.Cvi")

  MakeLoopPSAM(11, 7, 0, titles3=TRUE, pdf.plot="S15.Cvii")
  MakeLoopPSAM(11, 7, 1, pdf.plot="S15.Cviii")
  MakeLoopPSAM(11, 7, 2, pdf.plot="S15.Cix")
  MakeLoopPSAM(11, 7, 3, pdf.plot="S15.Cx")
  MakeLoopPSAM(11, 7, 4, pdf.plot="S15.Cxi")
  MakeLoopPSAM(11, 7, 5, titles=TRUE, pdf.plot="S15.Cxii")

  MakeLoopPSAM(11, 8, 0, titles2=TRUE, pdf.plot="S15.Cxiii")
  MakeLoopPSAM(11, 8, 1, pdf.plot="S15.Cxiv")
  MakeLoopPSAM(11, 8, 2, pdf.plot="S15.Cxv")
  MakeLoopPSAM(11, 8, 3, pdf.plot="S15.Cxvi")
  MakeLoopPSAM(11, 8, 4, pdf.plot="S15.Cxvii")
  MakeLoopPSAM(11, 8, 5, titles=TRUE, pdf.plot="S15.Cxviii")

  message("Done Fig. S14")
}


MakeSupplementalFigure15 <- function() {

  MakeLoopPSAM(13, 6, 0, dinuc=TRUE, titles2=TRUE, pdf.plot="S16.Ai")
  MakeLoopPSAM(13, 6, 1, dinuc=TRUE, pdf.plot="S16.Aii")
  MakeLoopPSAM(13, 6, 2, dinuc=TRUE, pdf.plot="S16.Aiii")
  MakeLoopPSAM(13, 6, 3, dinuc=TRUE, pdf.plot="S16.Aiv")
  MakeLoopPSAM(13, 6, 4, dinuc=TRUE, pdf.plot="S16.Av")
  MakeLoopPSAM(13, 6, 5, dinuc=TRUE, titles=TRUE, pdf.plot="S16.Avi")

  MakeLoopPSAM(13, 6, 0, dinuc_sim=TRUE, titles3=TRUE, pdf.plot="S16.Avii")
  MakeLoopPSAM(13, 6, 1, dinuc_sim=TRUE, pdf.plot="S16.Aviii")
  MakeLoopPSAM(13, 6, 2, dinuc_sim=TRUE, pdf.plot="S16.Aix")
  MakeLoopPSAM(13, 6, 3, dinuc_sim=TRUE, pdf.plot="S16.Ax")
  MakeLoopPSAM(13, 6, 4, dinuc_sim=TRUE, pdf.plot="S16.Axi")
  MakeLoopPSAM(13, 6, 5, dinuc_sim=TRUE, titles=TRUE, pdf.plot="S16.Axii")



  MakeLoopPSAM(13, 7, 0, dinuc=TRUE, titles2=TRUE, pdf.plot="S16.Axiii")
  MakeLoopPSAM(13, 7, 1, dinuc=TRUE, pdf.plot="S16.Axiv")
  MakeLoopPSAM(13, 7, 2, dinuc=TRUE, pdf.plot="S16.Axv")
  MakeLoopPSAM(13, 7, 3, dinuc=TRUE, pdf.plot="S16.Axvi")
  MakeLoopPSAM(13, 7, 4, dinuc=TRUE, pdf.plot="S16.Axvii")
  MakeLoopPSAM(13, 7, 5, dinuc=TRUE, titles=TRUE, pdf.plot="S16.Axviii")

  MakeLoopPSAM(13, 7, 0, dinuc_sim=TRUE, titles2=TRUE, pdf.plot="S16.Bi")
  MakeLoopPSAM(13, 7, 1, dinuc_sim=TRUE, pdf.plot="S16.Bii")
  MakeLoopPSAM(13, 7, 2, dinuc_sim=TRUE, pdf.plot="S16.Biii")
  MakeLoopPSAM(13, 7, 3, dinuc_sim=TRUE, pdf.plot="S16.Biv")
  MakeLoopPSAM(13, 7, 4, dinuc_sim=TRUE, pdf.plot="S16.Bv")
  MakeLoopPSAM(13, 7, 5, dinuc_sim=TRUE, titles=TRUE, pdf.plot="S16.Bvi")



  MakeLoopPSAM(13, 8, 0, dinuc=TRUE, titles3=TRUE, pdf.plot="S16.Bvii")
  MakeLoopPSAM(13, 8, 1, dinuc=TRUE, pdf.plot="S16.Bviii")
  MakeLoopPSAM(13, 8, 2, dinuc=TRUE, pdf.plot="S16.Bix")
  MakeLoopPSAM(13, 8, 3, dinuc=TRUE, pdf.plot="S16.Bx")
  MakeLoopPSAM(13, 8, 4, dinuc=TRUE, pdf.plot="S16.Bxi")
  MakeLoopPSAM(13, 8, 5, dinuc=TRUE, titles=TRUE, pdf.plot="S16.Bxii")

  MakeLoopPSAM(13, 8, 0, dinuc_sim=TRUE, titles2=TRUE, pdf.plot="S16.Bxiii")
  MakeLoopPSAM(13, 8, 1, dinuc_sim=TRUE, pdf.plot="S16.Bxiv")
  MakeLoopPSAM(13, 8, 2, dinuc_sim=TRUE, pdf.plot="S16.Bxv")
  MakeLoopPSAM(13, 8, 3, dinuc_sim=TRUE, pdf.plot="S16.Bxvi")
  MakeLoopPSAM(13, 8, 4, dinuc_sim=TRUE, pdf.plot="S16.Bxvii")
  MakeLoopPSAM(13, 8, 5, dinuc_sim=TRUE, titles=TRUE, pdf.plot="S16.Bxviii")

  message("Done Fig. S15")

}





MakeSupplementalFigure16 <- function() {

  MakeLoopPSAM(12, 6, 0, dinuc=TRUE, titles2=TRUE, pdf.plot="S17.Ai")
  MakeLoopPSAM(12, 6, 1, dinuc=TRUE, pdf.plot="S17.Aii")
  MakeLoopPSAM(12, 6, 2, dinuc=TRUE, pdf.plot="S17.Aiii")
  MakeLoopPSAM(12, 6, 3, dinuc=TRUE, pdf.plot="S17.Aiv")
  MakeLoopPSAM(12, 6, 4, dinuc=TRUE, pdf.plot="S17.Av")
  MakeLoopPSAM(12, 6, 5, dinuc=TRUE, titles=TRUE, pdf.plot="S17.Avi")

  MakeLoopPSAM(12, 6, 0, dinuc_sim=TRUE, titles3=TRUE, pdf.plot="S17.Avii")
  MakeLoopPSAM(12, 6, 1, dinuc_sim=TRUE, pdf.plot="S17.Aviii")
  MakeLoopPSAM(12, 6, 2, dinuc_sim=TRUE, pdf.plot="S17.Aix")
  MakeLoopPSAM(12, 6, 3, dinuc_sim=TRUE, pdf.plot="S17.Ax")
  MakeLoopPSAM(12, 6, 4, dinuc_sim=TRUE, pdf.plot="S17.Axi")
  MakeLoopPSAM(12, 6, 5, dinuc_sim=TRUE, titles=TRUE, pdf.plot="S17.Axii")



  MakeLoopPSAM(12, 7, 0, dinuc=TRUE, titles2=TRUE, pdf.plot="S17.Axiii")
  MakeLoopPSAM(12, 7, 1, dinuc=TRUE, pdf.plot="S17.Axiv")
  MakeLoopPSAM(12, 7, 2, dinuc=TRUE, pdf.plot="S17.Axv")
  MakeLoopPSAM(12, 7, 3, dinuc=TRUE, pdf.plot="S17.Axvi")
  MakeLoopPSAM(12, 7, 4, dinuc=TRUE, pdf.plot="S17.Axvii")
  MakeLoopPSAM(12, 7, 5, dinuc=TRUE, titles=TRUE, pdf.plot="S17.Axviii")

  MakeLoopPSAM(12, 7, 0, dinuc_sim=TRUE, titles2=TRUE, pdf.plot="S17.Bi")
  MakeLoopPSAM(12, 7, 1, dinuc_sim=TRUE, pdf.plot="S17.Bii")
  MakeLoopPSAM(12, 7, 2, dinuc_sim=TRUE, pdf.plot="S17.Biii")
  MakeLoopPSAM(12, 7, 3, dinuc_sim=TRUE, pdf.plot="S17.Biv")
  MakeLoopPSAM(12, 7, 4, dinuc_sim=TRUE, pdf.plot="S17.Bv")
  MakeLoopPSAM(12, 7, 5, dinuc_sim=TRUE, titles=TRUE, pdf.plot="S17.Bvi")



  MakeLoopPSAM(12, 8, 0, dinuc=TRUE, titles3=TRUE, pdf.plot="S17.Bvii")
  MakeLoopPSAM(12, 8, 1, dinuc=TRUE, pdf.plot="S17.Bviii")
  MakeLoopPSAM(12, 8, 2, dinuc=TRUE, pdf.plot="S17.Bix")
  MakeLoopPSAM(12, 8, 3, dinuc=TRUE, pdf.plot="S17.Bx")
  MakeLoopPSAM(12, 8, 4, dinuc=TRUE, pdf.plot="S17.Bxi")
  MakeLoopPSAM(12, 8, 5, dinuc=TRUE, titles=TRUE, pdf.plot="S17.Bxii")

  MakeLoopPSAM(12, 8, 0, dinuc_sim=TRUE, titles2=TRUE, pdf.plot="S17.Bxiii")
  MakeLoopPSAM(12, 8, 1, dinuc_sim=TRUE, pdf.plot="S17.Bxiv")
  MakeLoopPSAM(12, 8, 2, dinuc_sim=TRUE, pdf.plot="S17.Bxv")
  MakeLoopPSAM(12, 8, 3, dinuc_sim=TRUE, pdf.plot="S17.Bxvi")
  MakeLoopPSAM(12, 8, 4, dinuc_sim=TRUE, pdf.plot="S17.Bxvii")
  MakeLoopPSAM(12, 8, 5, dinuc_sim=TRUE, titles=TRUE, pdf.plot="S17.Bxviii")

  message("Done Fig. S16")

}



MakeSupplementalFigure17 <- function() {

  MakeLoopPSAM(11, 6, 0, dinuc=TRUE, titles2=TRUE, pdf.plot="S18.Ai")
  MakeLoopPSAM(11, 6, 1, dinuc=TRUE, pdf.plot="S18.Aii")
  MakeLoopPSAM(11, 6, 2, dinuc=TRUE, pdf.plot="S18.Aiii")
  MakeLoopPSAM(11, 6, 3, dinuc=TRUE, pdf.plot="S18.Aiv")
  MakeLoopPSAM(11, 6, 4, dinuc=TRUE, pdf.plot="S18.Av")
  MakeLoopPSAM(11, 6, 5, dinuc=TRUE, titles=TRUE, pdf.plot="S18.Avi")

  MakeLoopPSAM(11, 6, 0, dinuc_sim=TRUE, titles3=TRUE, pdf.plot="S18.Avii")
  MakeLoopPSAM(11, 6, 1, dinuc_sim=TRUE, pdf.plot="S18.Aviii")
  MakeLoopPSAM(11, 6, 2, dinuc_sim=TRUE, pdf.plot="S18.Aix")
  MakeLoopPSAM(11, 6, 3, dinuc_sim=TRUE, pdf.plot="S18.Ax")
  MakeLoopPSAM(11, 6, 4, dinuc_sim=TRUE, pdf.plot="S18.Axi")
  MakeLoopPSAM(11, 6, 5, dinuc_sim=TRUE, titles=TRUE, pdf.plot="S18.Axii")



  MakeLoopPSAM(11, 7, 0, dinuc=TRUE, titles2=TRUE, pdf.plot="S18.Axiii")
  MakeLoopPSAM(11, 7, 1, dinuc=TRUE, pdf.plot="S18.Axiv")
  MakeLoopPSAM(11, 7, 2, dinuc=TRUE, pdf.plot="S18.Axv")
  MakeLoopPSAM(11, 7, 3, dinuc=TRUE, pdf.plot="S18.Axvi")
  MakeLoopPSAM(11, 7, 4, dinuc=TRUE, pdf.plot="S18.Axvii")
  MakeLoopPSAM(11, 7, 5, dinuc=TRUE, titles=TRUE, pdf.plot="S18.Axviii")

  MakeLoopPSAM(11, 7, 0, dinuc_sim=TRUE, titles2=TRUE, pdf.plot="S18.Bi")
  MakeLoopPSAM(11, 7, 1, dinuc_sim=TRUE, pdf.plot="S18.Bii")
  MakeLoopPSAM(11, 7, 2, dinuc_sim=TRUE, pdf.plot="S18.Biii")
  MakeLoopPSAM(11, 7, 3, dinuc_sim=TRUE, pdf.plot="S18.Biv")
  MakeLoopPSAM(11, 7, 4, dinuc_sim=TRUE, pdf.plot="S18.Bv")
  MakeLoopPSAM(11, 7, 5, dinuc_sim=TRUE, titles=TRUE, pdf.plot="S18.Bvi")


  MakeLoopPSAM(11, 8, 0, dinuc_sim=TRUE, titles3=TRUE, pdf.plot="S18.Bvii")
  MakeLoopPSAM(11, 8, 1, dinuc_sim=TRUE, pdf.plot="S18.Bviii")
  MakeLoopPSAM(11, 8, 2, dinuc_sim=TRUE, pdf.plot="S18.Bix")
  MakeLoopPSAM(11, 8, 3, dinuc_sim=TRUE, pdf.plot="S18.Bx")
  MakeLoopPSAM(11, 8, 4, dinuc_sim=TRUE, pdf.plot="S18.Bxi")
  MakeLoopPSAM(11, 8, 5, dinuc_sim=TRUE, titles=TRUE, pdf.plot="S18.Bxii")

  MakeLoopPSAM(11, 8, 0, dinuc_sim=TRUE, titles2=TRUE, pdf.plot="S18.Bxiii")
  MakeLoopPSAM(11, 8, 1, dinuc_sim=TRUE, pdf.plot="S18.Bxiv")
  MakeLoopPSAM(11, 8, 2, dinuc_sim=TRUE, pdf.plot="S18.Bxv")
  MakeLoopPSAM(11, 8, 3, dinuc_sim=TRUE, pdf.plot="S18.Bxvi")
  MakeLoopPSAM(11, 8, 4, dinuc_sim=TRUE, pdf.plot="S18.Bxvii")
  MakeLoopPSAM(11, 8, 5, dinuc_sim=TRUE, titles=TRUE, pdf.plot="S18.Bxviii")

  message("Done Fig. S17")
}






PlotKdsAgainstThreePScore <- function(
  mirna, experiment, plotscore=TRUE, label_mirna=TRUE,
  xpos=20, ypos=20, height=3, width=3, pdf.plot=FALSE
) {
  if (experiment %in% c("equilibrium", "equilibrium2_nb")) {
    len_lim <- c(4, 8)
    sitelist <- "randthrp_suppcomp"
    corrected_kds <- FALSE
  }
  model_use <- SubfunctionCall(FitPairingAndOffsetModelSingle)

  val_3p_score <- apply(model_use$data, 1, function(row) {
    CalculateThreePScore(
      c(as.integer(row[3]), as.integer(row[4])), as.integer(row[7])
    )
  })
  SubfunctionCall(FigureSaveFile2)
  par(mar=c(3, 3, 1, 1))
  ymin <- 0.2
  ymax <- 800
  if (plotscore) {
    xmin <- -6
    xmax <- 0.5*(nchar(kMirnaSeqs[mirna]) - 8 - 4) + 4
    log_str <- "y"
  } else {
    log_str <- "xy"
    xmin <- ymin
    xmax <- ymax
  }
  BlankPlot(log=log_str)
  if (plotscore) {
    AddLinearAxis(1, 1, 5, label="Grimson 3' Score")
    x <- val_3p_score
  } else {
    AddLogAxis(1, label="Model-estimated Kd fold change")
    x <- 10^model_use$values
  }
  y <- 10^model_use$data$logkd
  AddLogAxis(2, label="Measured Kd fold change")
  lens <- as.character(as.integer(as.character(model_use$data$pos3p))
                       - as.integer(as.character(model_use$data$pos5p)) + 1)

  print(lens)
  segments(xmin, ymin, x1=xmax, y1=ymax, lty=2)
  # Add the points.
  points(x, y, col=kThrPLengthCols[lens], lwd=0, xpd=NA)
  # pos_5p <- as.character(model_use$data$pos5p)
  # pos_3p <- as.character(model_use$data$pos3p)
  # offset_str <- as.character(model_use$data$offset)
  # str_all <- sprintf("%s_%s_%s", pos_5p, pos_3p, offset_str)
  # identify(x, y, labels=str_all)

  # Potentially add the miRNA label.
  if (label_mirna) {
    if (mirna == "let-7a-21nt") {
      if (experiment == "equil_c_nb") {
        mirna_txt <- "let-7a rep."
      } else {
        mirna_txt <- "let-7a"
      }
    } else if (mirna == "let-7a_minus1") {
      mirna_txt <- "let-7a(-1)"
    } else if (mirna == "let-7a_plus1") {
      mirna_txt <- "let-7a(+1)"
    } else {
      mirna_txt <- mirna
    }
    xy <- GetPlotFractionalCoords(0.05, 1, log=log_str)
    text(x=xy[1], y=xy[2], labels=mirna_txt, adj=0, xpd=NA)
  }

  xy <- GetPlotFractionalCoords(0.05, 0.9, log=log_str)
  if (plotscore) x_cor <- x
  else           x_cor <- log(x)
  y_cor <- log(y)
  AddCorrelationToPlot(x=x_cor, y=y_cor, xpos=xy[1], ypos=xy[2], rsquared=TRUE)


  # Finish the plot.
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}


PlotReporterSiteECDF <- function(
  mirna="let-7a", experiment="twist_reporter_assay_3p_tp", exp_type="series",
  location="left", bulge_nucs="both", rep="both", height=3.4, width=4, xpos=20,
  ypos=20, pdf.plot=FALSE
) {
    if (rep == "both") {
      log2fc_df <- SubfunctionCall(GetThreePrimeReporterFoldChanges, rep=1)
      log2fc_df_2 <- SubfunctionCall(GetThreePrimeReporterFoldChanges, rep=2)
      log2fc_df[, "log_fc"] <- (log2fc_df[, "log_fc"] + log2fc_df_2[, "log_fc"])/2
    } else {
      log2fc_df <- SubfunctionCall(GetThreePrimeReporterFoldChanges)
    }
    log2fc_df_global <<- log2fc_df
    SubfunctionCall(FigureSaveFile2)
    if (location == "dual") xmin <- -3
    else                    xmin <- -3
    xmax <- 0.5
    ymin <- 0
    ymax <- 1
    BlankPlot()

    log2fc_df[, "Seed"] <- gsub("^No_site_[0-9]{1,2}$", "None", log2fc_df[, "Seed"])
    log2fc_df_global <<- log2fc_df
    base_sites <- c("None", "8mer-w6", "6mer", "7mer-A1", "7mer-m8", "8mer")
    for (base_site_i in base_sites) {
      if (base_site_i == "None") {
        location_use <- ""
      } else {
        location_use <- location
      }
      inds <- which((log2fc_df$Seed == base_site_i) &
                    (log2fc_df$ThreePrime == "") &
                    (log2fc_df$Location == location_use))
      print(base_site_i)
      print(length(inds))
      log2fc_df_use <- log2fc_df[inds, ]
      if (nrow(log2fc_df_use) != 0) {
        x <- seq(xmin, xmax, length.out=100)
        y <- ecdf(log2fc_df_use[, "log_fc"])(x)
        lines(x, y, col=kSiteColors[base_site_i], lwd=1, xpd=NA)
      }
    }


    thrp_sites <- c("4mer-m13.16", "9mer-m11.19", "9mer-m13.21")
    if (bulge_nucs == "A") {
      bulges <- c("none", "A", "AAAA")
      bulge_len <- FALSE
    } else if (bulge_nucs == "T") {
      bulges <- c("none", "T", "TTTT")
      bulge_len <- FALSE
    } else if (bulge_nucs == "both") {
      bulges <- c("none", "N", "NNNN")
      bulge_len <- TRUE
    }
    cols <- c("yellowgreen", "olivedrab", "forestgreen")
    lwds <- c(1, 1, 1, 1, 1)
    ltys <- c(3, 2, 1)
    names(lwds) <- bulges
    names(ltys) <- bulges
    names(cols) <- thrp_sites
    for (thrp_site in thrp_sites) {
      for (bulge in bulges) {
        if (bulge == "none") {
          bulge_use <- ""
        } else {
          bulge_use <- bulge
        }
        if (bulge_len) {
          print(nchar(bulge_use))
          inds <- which((log2fc_df[, "Seed"] == "8mer-w6") &
                        (log2fc_df[, "ThreePrime"] == thrp_site) &
                        (log2fc_df[, "Location"] == location) &
                        (nchar(log2fc_df[, "Bulge"]) == nchar(bulge_use)))
        } else {
          inds <- which((log2fc_df[, "Seed"] == "8mer-w6") &
                        (log2fc_df[, "ThreePrime"] == thrp_site) &
                        (log2fc_df[, "Location"] == location) &
                        (log2fc_df[, "Bulge"] == bulge_use))
        }
        log2fc_df_use <- log2fc_df[inds, ]
        if (nrow(log2fc_df_use) != 0) {
          x <- seq(xmin, xmax, length.out=100)
          y <- ecdf(log2fc_df_use[, "log_fc"])(x)
          lines(x, y, col=cols[thrp_site], lwd=lwds[bulge], lty=ltys[bulge], xpd=NA)
        }        
      }
    }

    inds <- which(log2fc_df[, "Seed"] == "lin-41")
    log2fc_df_use <- log2fc_df[inds, ]
    if (nrow(log2fc_df_use) != 0) {
      x <- seq(xmin, xmax, length.out=100)
      y <- ecdf(log2fc_df_use[, "log_fc"])(x)
      lines(x, y, col="goldenrod", lwd=1, xpd=NA)
    }



    # ############################# Add lin-41_line ############################
    ind_lin_41 = which(log2fc_df[, "Seed"] == "lin-41" && log2fc_df[, "Context"] == "lin-41_context")
    segments(x0=log2fc_df[ind_lin_41, "log_fc"], y0=0, y1=1, col="goldenrod", lty=2)
    # ############################## Add axes ##################################
    AddLinearAxis(1, tick.space=0.5, label.space=1,
                  label="Fold change (log2)")
    AddLinearAxis(2, tick.space=0.05, label.space=0.2,
                  label="Cumulative fraction (%)", percent=TRUE)
    # ############################# Add labels ##############################
    exp_type_strings <- c("Serial trans.", "Co-trans.")
    names(exp_type_strings) <- c("series", "parallel")
    location_strings <- c("Upstream site", "Downstream site", "Dual site")
    names(location_strings) <- c("left", "right", "dual")
    bulge_nuc_strings <- c("A/U bulges", "A bulges", "U bulges")
    names(bulge_nuc_strings) <- c("both", "A", "T")
    xy <- GetPlotFractionalCoords(0.025, 1)
    text(xy[1], xy[2], labels=mirna, adj=c(0, 1), xpd=NA)
    xy <- GetPlotFractionalCoords(0.025, 0.90)
    text(xy[1], xy[2], labels=location_strings[location], adj=c(0, 1), xpd=NA)
    xy <- GetPlotFractionalCoords(0.025, 0.80)
    text(xy[1], xy[2], labels=exp_type_strings[exp_type], adj=c(0, 1), xpd=NA)
    xy <- GetPlotFractionalCoords(0.025, 0.70)
    text(xy[1], xy[2], labels=bulge_nuc_strings[bulge_nucs], adj=c(0, 1), xpd=NA)
    ############################# Finish the plot ##############################
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

PlotReporterSiteECDFLegend <- function(xpos=20, ypos=20, height=4.2, width=2.8,
                                       pdf.plot=FALSE) {
  SubfunctionCall(FigureSaveFile2)
  xmin <- 0
  xmax <- 1
  ymin <- 0
  ymax <- 1
  par(mar=c(0, 0, 0, 0))
  BlankPlot()
  xy <- GetPlotFractionalCoords(0, 1)

  base_sites <- c("Nosite", "8mer-w6", "6mer", "7mer-A1", "7mer-m8", "8mer")
  thrp_sites <- c("4mer-m13.16", "9mer-m11.19", "9mer-m13.21")
  cols <- c("yellowgreen", "olivedrab", "forestgreen")

  Legend(xy, legend=c("None", rev(base_sites[-1])), lwd=1,
         legend_pch=NA, col=c("black", kSiteColors[rev(base_sites[-1])]),
         ncol=1, seg.len=1.2)

  xy <- GetPlotFractionalCoords(0, 0.65)
  Legend(xy,
         legend=c(sprintf("%s; %s-nt offset", rep(thrp_sites, each=3),
                          rep(c(" 0", "+1", "+4"), times=3)), "lin-41 sites",
                  "lin-41 3'-UTR context"),
         lwd=1,
         legend_pch=NA, col=c(rep(cols, each=3), "goldenrod", "goldenrod"),
         lty=c(rep(c(3, 2, 1), times=3), 1, 2),
         ncol=1, seg.len=1.2)


  # xy <- GetPlotFractionalCoords(0, 0.4)
  # Legend(xy, legend=c(" 0", "+1", "+4"), lty=c(3, 2, 1), lwd=1, legend_pch=NA,
  #        col="gray", ncol=1)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}



MakeSupplementalFigure18 <- function(bulge_nucs="both") {
  letter_use <- "A"
  if (bulge_nucs == "both") {
    extra_str <- ""
  } else {
    extra_str <- bulge_nucs
  }
  for (exp_type in c("series", "parallel")) {
    for (location in c("left", "right", "dual")) {
      i <- 1
      for (mirna in c("let-7a", "let-7a-21nt", "miR-1")) {
        print(exp_type)
        print(location)
        print(mirna)
        print(letter_use)
        print(kRomanNumerals[i])

        PlotReporterSiteECDF(mirna=mirna, exp_type=exp_type, location=location,
                             bulge_nucs=bulge_nucs,
                             pdf.plot=sprintf("S19.%s%s%s",
                                              letter_use,
                                              kRomanNumerals[i],
                                              extra_str))
        i <- i + 1
      }
      letter_use <- kNextLetter[letter_use]
    }
  }
  PlotReporterSiteECDFLegend(pdf.plot="S18.Aiv")
  message("Done Fig. S19")
}


