source("general/ThrPFunctions.R")
library(wCorr)
source("general/GenericFigures.R")
source("Repression/repression_workspace.R")
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


slopes_lm <- rep(NA, 7)

names(slopes_lm) <- c(kMirnas[1:6], "flanks")


kThrPLengthCols <- rev(viridis::plasma(10))[3:10]
names(kThrPLengthCols) <- 4:11
kThrPPositionCols <- rev(viridis::inferno(17))[3:17]
names(kThrPPositionCols) <- 9:23

kThrpReporterColors <- c("#006400", "darkolivegreen2", "#3EB000")
kThrpReporterColors <- c("#6A8F24", "#7BFC01", "#3CB371")
names(kThrpReporterColors) <- c("4mer-m13.16", "9mer-m11.19", "9mer-m13.21")


graphics.off()

kMaxValueColor <- "coral"

## Global variable for whether or not the offset axis should be inverted.
kOffsetInv <- "x"

kNucleotideColorsThrP <- c(`A`="blue", `C`="purple", `G`="red", `U`="forestgreen", `T`="forestgreen")

                                  # blue         # magenta      # orange       # teal
kNucleotideColorsColorBlind <- c(`A`="#0077BB", `C`="#EE3377", `G`="#EE7733", `U`="#009988", `T`="#009988")

kNucleotideColorsColorBlind <- c(`A`="#785EF0", `C`="#DC267F", `G`="#FE6100", `U`="#FFB000", `T`="#FFB000")

kNucleotideColorsColorBlind <- c(`A`="#0077BB", `C`="#33BBEE", `G`="#009988", `U`="#EE7733", `T`="#EE7733")

kNucleotideColorsColorBlind <- c(`A`="#66CCEE", `C`="#EE6677", `G`="#CCBB44", `U`="#228833", `T`="#228833")

kNucleotideColorsColorBlind <- c(`A`="#7A6CFF", `C`="#AE4B57", `G`="#F6921E", `U`="#00A14B", `T`="#00A14B")


################################################################################
# FIGURE 1
################################################################################

# 1C____________________________________________________________________________
PlotPositionalEnrichmentForProgramedLibrary <- function(
  mirna, experiment, condition, n_constant, kmer_len, seedex=TRUE, ps=0,
  toppos=5, start=1, stop=18, height=4.2, width=5.5, xpos=20, ypos=20,
  pdf.plot=FALSE, outputdata=FALSE
) {
  # Load the two data tables.
  pc_I <- SubfunctionCall(GetPositionalKmersProgrammedLib, condition="I")[, start:stop]
  pc_A <- SubfunctionCall(GetPositionalKmersProgrammedLib)[, start:stop]
  pc_I_global <<- pc_I
  pc_A_global <<- pc_A
  # STUFF TO HELP WITH WRITING THE METHODS #
  message("normalized pseudo count")
  message(ps)
  message("number of kmer positions averaged:")
  message(toppos)
  message("range of positions considered:")
  message(sprintf("%s to %s", start, stop))
  # Give the psuedo counts
  ps_I <- sum(pc_I)/nrow(pc_I)*ps
  ps_A <- sum(pc_A)/nrow(pc_A)*ps
  # Add the pseudo counts and normalize.
  norm_I <- (pc_I + ps_I)/sum(pc_I + ps_I)
  norm_A <- (pc_A + ps_A)/sum(pc_A + ps_A)
  # Make the enrichment table.
  R_A <- norm_A / norm_I
  R_mat <- R_A
  # Get the sum of the top five enrichment positions for each kmer:
  top_n_enrichments <- apply(R_A, 1, function(row) {
    row <- row[start:stop]
    row <- row[!(is.na(row))]
    top <- sort(row, decreasing=TRUE)
    return(sum(top[1:min(toppos, length(top))]))
  })
  # Sort the list of enrichments to get the top 20 rows.
  kmers_sorted <- sort(top_n_enrichments, decreasing=TRUE)
  kmers_sorted <<- kmers_sorted
  R_mat <- R_A[names(kmers_sorted)[1:20], 1:18]
  R_mat <<- R_mat
  str_over_5 <- length(which(R_mat[1, ] >= 5))
  message("The number of positions for which the top k-mer")
  message(sprintf("is greater than 5-fold enriched: %s", str_over_5))
  print(which(R_mat[1, ] >= 5))
  if (outputdata) {
    return(R_mat)
  } else {
    # Make the plot and define the limits.
    SubfunctionCall(FigureSaveFile2)
    par(mar=c(2.5, 8, 2, 4))
    xmin <- 1
    xmax <- 18
    ymin <- 0.5
    ymax <- 20
    BlankPlot()
    # Nucleotide positions to be labeled on the x-axis, relative to the miRNA.
    x_lib_labels <- c(1, 5, 10, 15, 20, 25)
    x_lib_pos <- 27 - x_lib_labels
    # This corrects for wanting the labels to transition from 1 to -1, with no 0-
    # numbered position.
    x_lib_pos[which(x_lib_pos > 27)] <- x_lib_pos[which(x_lib_pos > 27)] - 1
    # Make the x-axis
    AddLinearAxis(1, 1, 1, label=expression(italic(k)*"-mer position"), label_pos_ticks=TRUE,
                  alt_lab=x_lib_labels, alt_lab_pos=x_lib_pos)
    # Make the plot that will be used in the figure.
    # Define the left, right, bottom, and top positions for the squares in the
    # heatmap.
    xlefts <- (
      rep(seq(1, ncol(R_mat)), ymax) - 0.5
    )
    xright <- xlefts + 1
    ybottom <- rep(seq(ymax, 1), each=ncol(R_mat)) - 0.5
    ytop <- ybottom + 1
    colors <- ColorViridisPalette(c(t(R_mat)), steps=100, min=0, max=20, palettemax=0.8)
    rect(xlefts, ybottom, xright, ytop, col=colors, lwd=0.4,
         xpd=NA, border="white")
    # PRINT THE KMERS TO THE LEFT WITH ALIGNED SEQUENCES
    # Kmers to be printed to the left of the rows
    text(-4.5, 21, labels=expression("Ranked"~italic(k)*"-mers:"), adj=c(0.5, 0),
         xpd=NA)
    kmers = names(kmers_sorted)[1:20]
    # Makea matrix of the longest substring between each pair-wise kmer
    consensus_seq <- sapply(kmers, function(kmer1) {
      sapply(kmers, function(kmer2) {
        LargestCommonSubstring(kmer1, kmer2)
      })
    })
    # Make a list of only the unique substrings between each pairwise kmer.
    uniq_sequences <- unique(c(consensus_seq))
    # Identify the kmer that is the shortest of all of these,
    # as this is the minimal element connecting all of them.
    min_seq <- uniq_sequences[which(
      nchar(uniq_sequences) == min(nchar(uniq_sequences))
    )]
    # Get the position of this kmer within all of the kmers.
    print(kmers)
    pos_min_seq <- sapply(kmers, function(kmer) {
      out <- gregexpr(min_seq, kmer)[[1]]
      out
    })
    # Vector of the possible X-positions in the figure.
    xpos = seq(-9, -0.5,
               length.out=nchar(kmers[1]) + max(pos_min_seq) - min(pos_min_seq))
    # Figure out the coloring of the kmer nucleotides, based on complementarity to
    # the miRNA.
    mirna_seq <- kMirnaSeqs[mirna]
    rc_mirna_seq <- RevComplement(mirna_seq)
    # Get the sequence of complementarity to the miRNA.
    complementary_seqs <- sapply(kmers, LargestCommonSubstring, seq2=rc_mirna_seq)
    # Make a matrix that has the range of complementary sequence of each kmer to
    # the miRNA.
    start_stop_mat <- apply(rbind(complementary_seqs, kmers), 2, function(row_i) {
      start <- gregexpr(row_i[1], row_i[2])[[1]]
      stop <- start + nchar(row_i[1]) - 1
      return(c(start, stop))  
    })
    # Transpose the start_stop_matrix just to be able to go over rows.

    start_stop_mat <- t(start_stop_mat)
    sapply(seq(1, length(kmers)), function(r_ind) {
      kmer = ConvertTtoU(kmers[r_ind])
      y_pos = 21 - r_ind
      kmer_list = unlist(strsplit(kmer, split=""))
      start = max(pos_min_seq) - pos_min_seq[r_ind] + 1
      stop = start + nchar(kmer) - 1
      xpos_i = xpos[start : stop]
      cols_i <- rep("black", length(kmer_list))
      start_stop_row <- start_stop_mat[r_ind, ]
      red_i <- setdiff(seq(length(kmer_list)),
                       seq(start_stop_row[1], start_stop_row[2]))
      cols_i[red_i] <- "#EC008B"
      text(x=xpos_i, y=y_pos, labels=kmer_list, col=cols_i, adj=c(0.5, 0.5),
           xpd=NA)
    })
    # Make the key.
    y_div <- ymax/100
    kpl <- ncol(R_mat) + 1 # Left-hand position of the key
    kw <- 1 # Width of the key
    colors <- ColorViridisPalette(seq(0, 20, length.out=100), 100, 0, 20,
                                  palettemax=0.8)
    rect(xleft=kpl, ybottom=seq(100)*y_div - y_div + 0.5, xright=kpl + kw,
         ytop=seq(100)*y_div - y_div + 0.5 + y_div,
         col=colors, xpd=NA, border=NA)
    # Generate the axis for the legend and the label
    axis(4, at=seq(ymin + y_div/2, ymax + 0.5 - y_div/2, length.out=5),
         labels=seq(0, 20, by=5), lwd=par()$lwd, pos=kpl + kw, xpd=NA)
    text(x= kpl + kw + 2.5, y=10, labels="Enrichment", srt=270, xpd=NA)
    # Finish the plot.
    if (class(pdf.plot) == "character") {
      dev.off()
    }
  }
}

# ################################################################################
# # FIGURE S1
# ################################################################################

# # S1A-B_________________________________________________________________________
# PlotPairwiseThrPKds <- function(
#   mirna_1, experiment_1, n_constant_1, sitelist_1, combined_1, buffer_1, 
#   mirna_2, experiment_2, n_constant_2, sitelist_2, combined_2, buffer_2,
#   corrected_kds=FALSE, sumseed=FALSE, alpha=0.3, span=10, len_lim=c(4, 11),
#   thrp_only=FALSE, height=5, width=5, xpos=20, ypos=20, pdf.plot=FALSE
# ) {
#   # Get the two sets of Kds
#   if (grepl("suppcomp", sitelist_1)) {
#     grep_base_string <- "Comp"
#     grep_search_string <- "(^.*)\\|(.*)\\|(%s_Kd)"
#     bottom_label <- "Counts combined over\nseed-mismatch type"
#   } else {
#     grep_base_string <- "8mer-mm[ACTG][2-7]"
#     grep_search_string <- "(^.*)\\|(.*)\\|(%s_Kd)"
#     bottom_label <- "Counts separated by\nseed-mismatch type"
#   }
#   if (experiment_1 == "equilibrium") {
#     grep_search_string <- "(^.*)\\|(.*)\\|(%s_Kd)|(^.*)(_Kd)"
#     grep_replace_string <- "\\1\\4"
#   } else {
#     grep_replace_string <- "\\1"
#   }
#   # Get the two sets of Kds:
#   kds_collapsed_1 <- SubfunctionCall(
#     EquilPars, mirna=mirna_1, experiment=experiment_1, n_constant=n_constant_1,
#     sitelist=sitelist_1, combined=combined_1, buffer=buffer_1
#   )
#   # Allows for the kds to be plotted either corrected or not corrected.
#   if (corrected_kds) {
#     kds_collapsed_2 <- SubfunctionCall(
#       ApplyKdCorrection, mirna=mirna_2, experiment=experiment_2,
#       prog_n_constant=n_constant_2, rand_n_constant=n_constant_1,
#       prog_sitelist=sitelist_2, rand_sitelist=sitelist_1,
#       combined_rand=combined_1, buffer_rand=buffer_1
#     )
#   } else {
#     kds_collapsed_2 <- SubfunctionCall(
#       EquilPars, mirna=mirna_2, experiment=experiment_2,
#       n_constant=n_constant_2, sitelist=sitelist_2, combined=combined_2
#     )
#   }
#   # This swaps the seed site kds for the averages of those values over the
#   # the varous values of it in the random region.
#   if (grepl("randthrp", sitelist_1) & (!corrected_kds)) {
#     kds_collapsed_2 <- SubfunctionCall(SwapProgrammedKds, kds=kds_collapsed_2,
#                                        suppcomp=grepl("suppcomp", sitelist_2))
#   }
#   kds_collapsed_1_global <<- kds_collapsed_2

#   # Plot only the intersection of the two kd sets.
#   if (experiment_1 == "equilibrium") {
#     sites_shared <- intersect(rownames(kds_collapsed_2),
#                               rownames(kds_collapsed_1))
#   } else {
#     sites_shared <- intersect(rownames(kds_collapsed_1),
#                               rownames(kds_collapsed_2))      
#   }
#   # This conditional is necessarily because I don't currently look for
#   # `8mer|NN|Comp_Kd` etc when performing the site counting in the random
#   # library. I rather just look for seed sites, and then average the
#   # `8mer|NN|Comp_Kd` in the programmed library. Because "grep_base_string"
#   # looks for the form `^site|NN|8mer-mmN#` or `site|NN|Comp`, it would drop the
#   # base sites in the random without the sites[1:24] part.
#   if (grepl("randthrp", sitelist_1)) {
#     sites <- c(sprintf("%s_Kd", c(kSeedSites, GetAll8merMmSites(mirna_1))),
#                grep("&", grep(sprintf("%s_Kd", grep_base_string), sites_shared,
#                               value=TRUE), invert=TRUE, value=TRUE))
#   } else {
#     sites <- grep("&", grep(sprintf("%s_Kd", grep_base_string), sites_shared,
#                             value=TRUE), invert=TRUE, value=TRUE)
#   }
#   # Get the names of the sites that isn't the seed/programmed base site.
#   site_names <- gsub(sprintf(grep_search_string, grep_base_string),
#                      replacement=grep_replace_string, sites, perl=TRUE)
#   # Get the unique set of these sites
#   site_names_unique <- unique(site_names)
#   # sites_unique_use <- site_names_unique

#   site_names_unique <- grep("Kd", site_names_unique, invert=TRUE, value=TRUE)


#   seed_site_names <- site_names[which(site_names %in% site_names_unique[1:24])]
#   thrp_site_names <- site_names[which(!(site_names %in% site_names_unique[1:24]))]
#   thrp_site_lens <- as.integer(gsub("^(.*)mer.*", replacement="\\1", thrp_site_names))

#   thrp_sites_use <- thrp_site_names[which(thrp_site_lens >= len_lim[1] & thrp_site_lens <= len_lim[2])]
#   print(thrp_sites_use)
#   if (thrp_only) {
#     sites_unique_use <- thrp_sites_use
#   } else {
#     sites_unique_use <- unique(c(thrp_sites_use, seed_site_names))
#   }

#   sites_unique_use <- grep("Kd", sites_unique_use, invert=TRUE, value=TRUE)
#   thrp_site_lens <- as.integer(gsub("^(.*)mer.*", replacement="\\1", thrp_sites_use))

#   sites_use <- sites[which(site_names %in% sites_unique_use)]
#   # message("head sites_use")
#   # print(head(sites_use))
#   # message("tail sites_use")
#   # print(tail(sites_use))
#   # message("head rownames(kds_collapsed_1")
#   # print(head(rownames(kds_collapsed_1)))

#   sXc_programmed <- SubfunctionCall(
#     SitesXCounts, mirna=mirna_1, experiment=experiment_1,
#     n_constant=n_constant_1, sitelist=sitelist_1, buffer=buffer_1
#   )
#   sites_use_counts <- gsub("^(.*)(_Kd)$", sites_use, replacement="\\1")
#   # Remove Kd values which have fewer than 5 counts in the input library.
#   sites_use <- sites_use[which(sXc_programmed[sites_use_counts, 1] > 5)]

#   kds_1 <- kds_collapsed_1[sites_use, ]
#   kds_2 <- kds_collapsed_2[sites_use, ]

#   # Gets the left hand portion of the bipartite sites, and also gets rid of the
#   # `_Kd` suffix, in the case of plotting against the random libraries.
#   sites_base <- sapply(sites_use, function(site_use) {
#     unlist(strsplit(unlist(strsplit(site_use, split="\\|"))[1], split="_Kd"))[1]
#   })
#   # Assign the colors to each row.
#   inds_seed <- which(sites_base %in% c(kSeedSites, GetAll8merMmSites(mirna_1)))
#   inds_non_seed <- setdiff(1:length(sites_base), inds_seed)

#   lens_colors <- sapply(sites_use[inds_non_seed], function(site_use) {
#     unlist(strsplit(site_use, split="mer"))[1]  
#   })
#   cols_thrp <- ConvertRColortoRGB(kThrPLengthCols, alpha=alpha)
#   names(cols_thrp) <- 4:11
#   cols <- rep("black", length(sites_use))
#   cols[inds_seed] <- kSiteColors[sites_base[inds_seed]]
#   cols[inds_non_seed] <- cols_thrp[as.character(lens_colors)]
#   site_cols <- cols
#   # Reorder the sites to make seed sites first and non-seed sites second.
#   if (sitelist_1 == "progthrp") {
#     sites_use <- sites_use[c(inds_seed, inds_non_seed)]
#   }
#   #
#   SubfunctionCall(FigureSaveFile2)
#   if (!(grepl("suppcomp", sitelist_1)) & !(grepl("suppcomp", sitelist_2))) {
#     xmin <- 0.0005
#     xmax <- 5
#   } else {
#     if (grepl("randthrp", sitelist_1)) {
#       xmin <- 0.0001
#     } else {
#       xmin <- 0.001
#     }
#     xmax <- 1  
#   }
#   ymin <- xmin
#   ymax <- xmax
#   BlankPlot(log="xy")
#   segments(xmin, ymin, x1=xmax, y1=ymax, col="gray")
#   if (grepl("randthrp", sitelist_1) & grepl("suppcomp", sitelist_2)) {
#     ################# Plot the loess curve line ########################
#     l_kds_2 <- log(kds_2[1:24, 2])
#     l_kds_1 <- log(kds_1[1:24, 2])
#     l_kds_dif <- l_kds_1 - l_kds_2
#     loess_object <- loess(l_kds_dif ~ l_kds_2, span=span,
#                           control=loess.control(surface="direct"))
#     y_line <- seq(min(l_kds_2) - 0.5, max(l_kds_2) + 0.5, length.out=50)
#     x_line <- y_line + predict(loess_object, data.frame(l_kds_2=y_line))
#     lines(exp(x_line), exp(y_line), col="purple")   
#     ########### Plot the loess curve line for the three-prime sites ############
#     # l_kds_2_thrp <- log(kds_2[25:nrow(kds_2), 2])
#     # l_kds_1_thrp <- log(kds_1[25:nrow(kds_1), 2])
#     # l_kds_dif_thrp <- l_kds_1_thrp - l_kds_2_thrp
#     # loess_object_global <<- loess_object
#     # loess_object_thrp <- loess(l_kds_dif_thrp ~ l_kds_2_thrp, span=span,
#     #                            control=loess.control(surface="direct"))
#     # y_line_thrp <- seq(min(l_kds_2_thrp) - 0.5, max(l_kds_2_thrp) + 0.5, length.out=50)
#     # x_line_thrp <- y_line_thrp + predict(loess_object_thrp, data.frame(l_kds_2_thrp=y_line_thrp))

#   }
#   ################ Plot the points ########################
#   points(kds_1[sites_use, 2],
#          kds_2[sites_use, 2],
#          col=site_cols, pch=20, xpd=NA)

#   ################ Add the miRNA label ########################
#   xy <- GetPlotFractionalCoords(0.1, 0.95, log="xy")
#   if (mirna_1 == "let-7a-21nt") {
#     mirna_1 <- "let-7a"
#   }
#   if (!corrected_kds){
#     text(xy[1], xy[2], labels=mirna_1, adj=0, xpd=NA)
#   }
#   ################ Add the r^2 and n ########################
#   xy <- GetPlotFractionalCoords(0.1, 0.9, log="xy")
#   AddCorrelationToPlot(x=log(kds_1[, 2]), y=log(kds_2[, 2]),
#                        xpos=xy[1], ypos=xy[2], rsquared=TRUE)
#   xy <- GetPlotFractionalCoords(0.1, 0.85, log="xy")
#   text(xy[1], xy[2], labels=sprintf("n = %s", nrow(kds_1)), adj=0, xpd=NA)
#   ################ Add the text labels ########################
#   xy <- GetPlotFractionalCoords(1, 0.05, log="xy")
#   text(xy[1], xy[2], labels=bottom_label, adj=c(1, 0), xpd=NA)
#   if (experiment_1 == "equilibrium") {
#     xy <- GetPlotFractionalCoords(1, 0.20, log="xy")
#     if (corrected_kds) {
#       text_label <- "Corrected values"
#     } else {
#       text_label <- "Raw values"
#     }
#     text(xy[1], xy[2], labels=text_label, adj=c(1, 0), xpd=NA)   
#   }
#   xy <- GetPlotFractionalCoords(1, 0.05, log="xy")
#   text(xy[1], xy[2], labels=bottom_label, adj=c(1, 0), xpd=NA)
#   ######################## Add the axes ########################
#   if (grepl("randthrp", sitelist_1)) {
#     x_string <- "Relative Kd; random-sequence library"
#     y_string <- "Relative Kd; programmed library"
#   } else {
#     x_string <- "Relative Kd; replicate 1"
#     y_string <- "Relative Kd; replicate 2"
#   }
#   # Add the axes:
#   AddLogAxis(1, x_string)
#   AddLogAxis(2, y_string)
#   # Finish the plot.
#   if (class(pdf.plot) == "character") {
#     print(pdf.plot)
#     print(corrected_kds)
#     dev.off()
#   }
#   print("done function")
# }

# ################################################################################
# # FIGURE 2
# ################################################################################

# # 2A, 3A&C______________________________________________________________________
# PlotThreePrimeSiteECDF <- function(
#   mirna="let-7a-21nt", experiment="equil_c2_nb", n_constant=3, sumseed=FALSE,
#   corrected_kds=TRUE, label_mirna=FALSE, height=5, width=5, nbomitc=FALSE,
#   xpos=20, ypos=20, pdf.plot=FALSE, outputdata=FALSE
# ) {
#   if (corrected_kds) {
#     log_xmin <- -4
#     kds <- SubfunctionCall(ApplyKdCorrection, prog_n_constant=n_constant,
#                            rand_n_constant=n_constant,
#                            prog_sitelist="progthrp_suppcomp")
#   } else {
#     log_xmin <- -3
#     kds <- SubfunctionCall(EquilPars, sitelist="progthrp_suppcomp")
#   }
#   if (outputdata) {
#     return(kds)
#   } else {
#     SubfunctionCall(FigureSaveFile2)
#     xmin <- 10^log_xmin
#     xmax <- 1
#     ymin <- 0
#     ymax <- 1
#     BlankPlot(log="x")
#     # cols <- kThrPLengthCols
#     # Plot the black EDF line and set up the geometric mean reference point for
#     # the statements being made.
#     grep_str <- "^8mer-mm[ACTG][2-7]_Kd$"
#     inds <- grep(grep_str, rownames(kds), perl=TRUE, value=TRUE)
#     df_use <- kds[inds, ]
#     ref_kd <- GeoMean(df_use[, 2])
#     x <- 10^seq(-3, 0, length.out=100)
#     y <- ecdf(df_use[, 2])(x)
#     lines(x, y, col="black", lwd=1, xpd=NA)
#     min_ecdf <- Inf
#     max_ecdf <- -1*Inf
#     min_site <- ""
#     max_site <- ""
#     total_plotted <- 0
#     names_use <<- c()
#     sapply(c(4:11), function(kmer) {
#       #### Get the relevant sites ###########
#       target <- sprintf("^%smer-m[0-9]{1,2}\\.[0-9]{1,2}\\|.*\\|Comp_Kd", kmer)
#       inds <- grep(target, rownames(kds), perl=TRUE, value=TRUE)
#       df_use <- kds[inds, ]
#       names_use <<- c(names_use, inds)
#       total_plotted <<- total_plotted + length(inds)
#       #### Plot the ECDF line ############################
#       x <- 10^seq(log_xmin, 0, length.out=100)
#       y <- ecdf(df_use[, 2])(x)
#       lines(x, y, col=kThrPLengthCols[as.character(kmer)], lwd=1, xpd=NA)
#       ###### Update the min value of the ECDF and the min site #####
#       min_x <- min(df_use[, 2])
#       max_x <- max(df_use[, 2])
#       if (min_x < min_ecdf) {
#         min_ecdf <<- min_x
#         min_site <<- rownames(df_use)[which(df_use[, 2] == min_x)]
#       }
#       if (max_x > max_ecdf) {
#         max_ecdf <<- max_x
#         max_site <<- rownames(df_use)[which(df_use[, 2] == max_x)]
#       }
#       ### Print the number of sites that are 2- and 10- fold over background. ##
#       n_2fold <- length(which(df_use[, 2] >= ref_kd/2))
#       n_10fold <- length(which(df_use[, 2] <= ref_kd/10))
#       n_sites <- length(inds)
#       message(sprintf("%smers: %s percent over 2-fold of %s", kmer,
#                       n_2fold/n_sites, n_sites))
#       message(sprintf("%smers: %s percent over 10-fold of %s", kmer,
#                       n_10fold/n_sites, n_sites))
#       message(sprintf("%smer median Kd fold-change value: %s", kmer, ref_kd/median(df_use[, 2])))
#     })
#     # Print messages related to the min, max, range and reference Kd for this
#     # experiment.
#     message(sprintf("Max is %s, for the %s.", max_ecdf, max_site))
#     message(sprintf("Min is %s, for the %s.", min_ecdf, min_site))
#     message(sprintf("Range is %s.", max_ecdf/min_ecdf))
#     message(sprintf("Ref Kd is %s.", ref_kd))
#     message(sprintf("Total number of sites potted: %s", total_plotted))
#     ############################## Add axes ############################
#     AddLogAxis(1, label="Relative Kd")
#     AddLinearAxis(2, tick.space=0.05, label.space=0.2,
#                   label="Cumulative fraction (%)", percent=TRUE)
#     ################# Add legend #############################################
#     xy <- GetPlotFractionalCoords(0.025, 1, log="x")
#     text(xy[1], xy[2], labels="Length of 3'\npairing (bp)", adj=c(0, 1),
#          xpd=NA)
#     xy <- GetPlotFractionalCoords(0.025, 0.925, log="x")
#     Legend(xy, legend=c("<4", names(kThrPLengthCols)),
#            col=c("black", kThrPLengthCols), ncol=1)
#     ############################# Add miRNA label ##############################
#     if (label_mirna) {
#       if (mirna == "let-7a-21nt") {
#         if (experiment == "equil_c_nb") {
#           mirna_txt <- "let-7a rep."
#         } else {
#           mirna_txt <- "let-7a"
#         }
#       } else if (mirna == "let-7a_minus1") {
#         mirna_txt <- "let-7a(-1)"
#       } else if (mirna == "let-7a_plus1") {
#         mirna_txt <- "let-7a(+1)"
#       } else {
#         mirna_txt <- mirna
#       }
#       mtext(text=mirna_txt, side=3, line=0.5, at=xmax, adj=1, cex=par("cex"))
#     }
#     ############################# Finish the plot ##############################
#     if (class(pdf.plot) == "character") {
#       dev.off()
#     }
#   }
# }


# PlotCombinedEnrichmentAndKd <- function(site,
#   mirna="let-7a-21nt", experiment="equil_c2_nb", n_constant_r=0,
#   n_constant_kd=3, sitelist="progthrp_suppcomp", nbomitc=FALSE,
#   corrected_kds=TRUE, len_lim=c(4, 11), offset_lim=c(-4, 16), collapsemm=FALSE,
#   equilibrium_nb=FALSE, justlegend=FALSE, condition=40, start=1, stop=18,
#   seedex=TRUE, ps=0, plot_enrich=FALSE, height=2, width=5, xpos=20, ypos=20,
#   pdf.plot=FALSE 
# ) {
#   ###### Get the enrichment values to plot alongside the Kd values. ############
#   # Load the two data tables.
#   kmer_len <- as.numeric(strsplit(site, split="mer")[[1]][1])
#   print(kmer_len)
#   n_constant <- n_constant_r
#   pc_I <- SubfunctionCall(GetPositionalKmersProgrammedLib,
#                           condition="I")[, start:stop]
#   pc_A <- SubfunctionCall(GetPositionalKmersProgrammedLib)[, start:stop]
#   pc_I_global <- pc_I
#   pc_A_global <- pc_A
#   # STUFF TO HELP WITH WRITING THE METHODS #
#   # Give the psuedo counts
#   ps_I <- sum(pc_I)/nrow(pc_I)*ps
#   ps_A <- sum(pc_A)/nrow(pc_A)*ps
#   # Add the pseudo counts and normalize.
#   norm_I <- (pc_I + ps_I)/sum(pc_I + ps_I)
#   norm_A <- (pc_A + ps_A)/sum(pc_A + ps_A)
#   # Make the enrichment table.
#   R_A <- norm_A / norm_I
#   R_mat <- R_A
#   print(head(R_mat))
#   ## GET THE KDs #############
#   n_constant <- n_constant_kd
#   if (corrected_kds) {
#     kds <- SubfunctionCall(ApplyKdCorrection, prog_n_constant=n_constant,
#                            rand_n_constant=n_constant, prog_sitelist=sitelist)
#   } else { 
#     if (mirna == "miR-1" & experiment == "equilibrium") {
#       combined <- FALSE
#       buffer <- TRUE
#     } else if (mirna == "miR-7-23nt" & experiment == "equilibrium2_nb") {
#       combined <- FALSE
#       buffer <- FALSE
#     } else {
#       combined <- TRUE
#       buffer <- FALSE
#     }
#     kds <- SubfunctionCall(EquilPars, sitelist=sitelist)
#   }

#   target <- sprintf("^%s\\|.*\\|Comp_Kd", site)
#   inds <- grep(target, rownames(kds), perl=TRUE, value=TRUE)
#   df_use <- kds[inds, ]
#   sites_inds <- grep(sprintf("^%s\\|.*\\|Comp_Kd", site), rownames(df_use), perl=TRUE)
#   site_pos <- as.integer(gsub("^.*mer-m(.*)\\..*$", replacement="\\1", site))
#   pos <- as.integer(gsub("^(.*)\\|(.*)\\|(Comp_Kd$)", replacement="\\2", rownames(df_use)))
#   offsets <- pos - site_pos
#   inds_pos <- which(offsets >= offset_lim[1] & offsets <= offset_lim[2])
#   offsets <- offsets[inds_pos]
#   df_use_pos <- df_use[inds_pos, ]
#   # df_use_pos <- df_use_pos[order(offsets), ]

#   site_sequence <- GetSiteSeq(mirna, site)
#   print(site_sequence)
#   R_mat_use <- R_mat[site_sequence, ]

#   print(R_mat_use)
#   names(R_mat_use) <- -1*(1:length(R_mat_use)) + 25 + n_constant_r - kmer_len + 1 + 9
#   print(R_mat_use)

#   rownames(df_use_pos) <- gsub("^(.*)\\|(.*)\\|(.*)$", replace="\\2", rownames(df_use_pos), perl=TRUE)
#   if (plot_enrich) height <- 1.75
#   else             height <- 2.25
#   SubfunctionCall(FigureSaveFile2)
#   xmin <- -4 + site_pos
#   xmax <- 16 + site_pos
#   if (plot_enrich) {
#     par(mar=c(0.5, 3, 0.5, 4))
#     ymin <- 1
#     ymax <- 100
#   } else {
#     par(mar=c(3, 3, 0.5, 4))
#     ymin <- 1e-3
#     ymax <- 1e-1
#   }

#   BlankPlot(log="y", inv="x")
#   x <- 9:26
#   if (plot_enrich) {
#     y_lab <- "Enrichment"    
#     colors <- ColorViridisPalette(R_mat_use[as.character(9:26)], steps=100,
#                                   min=0, max=20, palettemax=0.8)
#     lines(9:26, R_mat_use[as.character(9:26)], col="gray")
#     points(9:26, R_mat_use[as.character(9:26)], pch=19,
#            cex=2*par("cex"), col=colors)
#     segments(x0=xmin, y0=ymin, x1=xmax, lty=2, xpd=NA)
#     # Add the top line of the label, that gives the kmer motif itself. #########
#     site_lab <- ConvertTtoU(site_sequence)
#     xy <- GetPlotFractionalCoords(0.025, 0.95, log="y", inv="x")
#     text(xy[1], xy[2], labels=site_lab, adj=0)
#     # Add the second label saying "Piars to ....."
#     xy <- GetPlotFractionalCoords(0.025, 0.82, log="y", inv="x")
#     site_pos <- gsub("^.*mer-m(.*)*$", replacement="\\1", site)
#     site_pos <- gsub("\\.", replacement="-", site_pos)
#     text(xy[1], xy[2], labels=sprintf("Pairs to nt %s", site_pos), adj=0)

#   } else {
#     cols <- kThrPLengthCols[as.character(kmer_len)]
#     segments(x0=x, y0=df_use_pos[as.character(9:26), 3],
#              y1=df_use_pos[as.character(9:26), 5], col=cols, lwd=1, xpd=NA)
#     lines(x, df_use_pos[as.character(9:26), 2], col=cols, lwd=1, xpd=NA)
#     points(x, df_use_pos[as.character(9:26), 2], col=cols, pch=20, xpd=NA)      

#     # lines(9:26, df_use_pos[as.character(9:26), 2])
#     y_lab <- "Relative Kd"
#     AddLinearAxis(1, 1, 4, label="")
#     # ymin <- ymin/sqrt(10)
#     offset_labs <- c(-4, 0, 4, 8, 12, 16)
#     offset_x <- offset_labs + site_pos
#     text(x=offset_x, y=ymin/4, labels=offset_labs, adj=0.5, xpd=NA)
#     par(lheight=.8)
#     mtext(side=1, text="k-mer\nposition", cex=par("cex"), at=xmin - 1, adj=0, padj=0, line=0.25)

#     mtext(side=1, text="Offset", at=xmin - 1, cex=par("cex"), adj=0, padj=0, line=1.5)
#     segments(x0=xmin, y0=ymax, x1=xmax, lty=2, xpd=NA)
#   }
#   AddLogAxis(2, label=y_lab)




#   # kds_use_global <<- kds_use
#   # R_mat_use_global <<- R_mat_use
#   # print(site_pos)

#   # print(kds_use)
#   # print(R_mat_use)
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }


# # 2B, 3B&D______________________________________________________________________
# PlotBestThreePrimeSite <- function(
#   mirna="let-7a-21nt", experiment="equil_c2_nb", n_constant=3,
#   sitelist="progthrp_suppcomp", nbomitc=FALSE, corrected_kds=TRUE,
#   len_lim=c(4, 11), collapsemm=FALSE, equilibrium_nb=FALSE, justlegend=FALSE,
#   alt_height=FALSE, height=5, width=10, xpos=20, ypos=20, pdf.plot=FALSE
# ) {
#   ################### Load the Kds ########################################
#   # The text_bound is the fold-difference in comparison to the geometric mean of
#   # the 18 programmed kds, first was 1.5, but trying some other values.
#   text_bound <- 1.5

#   # These are calculations to allow multiple lengths to be shown in the
#   # justlegend=TRUE plots shown in Figure 7, and Supplemental Figure 12&13.
#   ymin_f <- 0.001/6
#   ymax_f <- 0.7/6
#   # fc_f <- 1.8
#   # fc_f <- 1.4
#   num_sites <- len_lim[2] - len_lim[1] + 1
#   y_mir_site <- (ymax_f*ymin_f^9)^(1/10)

#   # y_first_site <- fc_f*y_mir_site
#   # y_top_site <- fc_f^8*y_mir_site
#   # ymax_adjusted <- fc_f^(num_sites - 8)*ymax_f
#   # ymax_adjusted <- 0.015
#   ymax_adjusted <- 0.001*1.5*1.4^(num_sites + 2)
#   # print(sprintf("y_top_site:%s", y_top_site))
#   print(sprintf("ymax_adjusted:%s", ymax_adjusted))
#   # Find the fraction of the 10 units used up by the difference from the
#   # minimum to the first site, and the top site to the maximum:
#   # user_bottom <- log10(y_first_site/ymin_f)/log10(ymax_f/ymin_f)*10
#   # user_top <- log10(ymax_f/y_top_site)/log10(ymax_f/ymin_f)*10
#   # Adjust the height for a variable number of sites:
#   # height_adjusted <- ((10 - user_bottom - user_top)*num_sites/8 + user_bottom + user_top + 1)*2.5/11

#   height_adjusted <- 2* 1/11 * (1 + 10*((log(1.5) + (num_sites + 2)*log(1.4))/(log(1.5) + 7*log(1.4))))
#   print(sprintf("height_adjusted: %s", height_adjusted))
#   if (corrected_kds) {
#     kds <- SubfunctionCall(ApplyKdCorrection, prog_n_constant=n_constant,
#                            rand_n_constant=n_constant, prog_sitelist=sitelist)
#   } else { 
#     if (mirna == "miR-1" & experiment == "equilibrium") {
#       combined <- FALSE
#       buffer <- TRUE
#     } else if (mirna == "miR-7-23nt" & experiment == "equilibrium2_nb") {
#       combined <- FALSE
#       buffer <- FALSE
#     } else {
#       combined <- TRUE
#       buffer <- FALSE
#     }
#     kds <- SubfunctionCall(EquilPars, sitelist=sitelist)
#   }
#   kds_global <<- kds
#   break
#   names_cols_use <- as.character(len_lim[1]:len_lim[2])
#   cols <- kThrPLengthCols[names_cols_use]
#   if (justlegend) {
#     width <- 4
#     height <- 2.5 
#     height <- height_adjusted
#     SubfunctionCall(FigureSaveFile2)
#     par(mar=c(1, 0, 0, 0))  
#     xmin <- 21.5
#     xmax <- 41
#     ymin <- 0.001
#     ymax <- ymax_adjusted
#     # message(sprintf("ymin is %s", ymin))
#     # message(sprintf("ymax is %s", ymax))
#     BlankPlot(log="y")
#     # segments(x0=xmin, x1=xmax, y0=ymin, xpd=NA, lwd=1, col="black")
#     # segments(x0=xmin, x1=xmax, y0=ymin*1.5, xpd=NA, lwd=1, col="red")
#     # segments(x0=xmin, x1=xmax, y0=ymin*1.5*1.4, xpd=NA, lwd=1, col="purple")
#     # segments(x0=xmin, x1=xmax, y0=ymin*1.5*1.4*1.4, xpd=NA, lwd=1, col="purple")
#     # segments(x0=xmin, x1=xmax, y0=ymin*1.5*1.4*1.4*1.4, xpd=NA, lwd=1, col="purple")
#     # segments(x0=xmin, x1=xmax, y0=ymin*1.5*1.4*1.4*1.4*1.4, xpd=NA, lwd=1, col="purple")
#     # segments(x0=xmin, x1=xmax, y0=ymin*1.5*1.4*1.4*1.4*1.4*1.4, xpd=NA, lwd=1, col="purple")
#     # # segments(x0=xmin, x1=xmax, y0=y_top_site, xpd=NA, lwd=1, col="blue")
#     # # segments(x0=xmin, x1=xmax, y0=ymax_adjusted, xpd=NA, lwd=4, col="green")

#   } else {
#     if (alt_height) {
#       height <- 3.7
#     }
#     #################### Open the plot window ##############################
#     SubfunctionCall(FigureSaveFile2)
#     if (alt_height || height == 3.7*9/9.75) {
#       par(mar=c(3, 3.58, 0.7, 1.42))
#     }
#     xmin <- -4
#     xmax <- 42
#     if (corrected_kds) {
#       ymin <- 0.0001
#     } else {
#       ymin <- 0.0003
#     }
#     ymax <- 1
#     BlankPlot(log="y")
#     xmax <- 16 # This allows the axes to not extend all the way across the plot.

#     if (sitelist == "progthrp_suppcomp" & collapsemm) {
#       ref_kd <- kds["Comp_Kd", 2]
#     } else {
#       mm_inds <- grep("^8mer-mm[ACTG][2-7]_Kd$", rownames(kds),
#                                  perl=TRUE)
#       ref_kd <- GeoMean(kds[mm_inds, 2])
#     }
#     segments(x0=xmin, y0=ref_kd, x1=xmax + 0.5, lty=2, lwd=0.5)
#     # segments(x0=xmin, y0=ref_kd/text_bound, x1=xmax, lty=3, lwd=0.5)
#     x_skew <- seq(-0.1, 0.1, length.out=len_lim[2] - len_lim[1] + 1)
#     names(x_skew) <- names(cols)
#   }
#   print("made it here")
#   if (!justlegend) {
#     s8mer_kd <- GeoMean(kds[grep("^8mer\\|.*\\|Comp_Kd$",
#                                  rownames(kds), perl=TRUE), 2])
#     s6mer_kd <- GeoMean(kds[grep("^6mer\\|.*\\|Comp_Kd$",
#                                  rownames(kds), perl=TRUE), 2])
#     if (sitelist == "programmed_suppcomp" & collapsemm) {
#       mm_kd    <- GeoMean(kds[grep("^Comp_Kd$",
#                                    rownames(kds), perl=TRUE), 2])
#     } else {
#       mm_kd    <- GeoMean(kds[grep("^8mer-mm[ACGT][2-7]_Kd$",
#                                    rownames(kds), perl=TRUE), 2])    
#     }
#     # s8mer_kd_norm <- s8mer_kd*mm_kd/(s8mer_kd + mm_kd)
#     # s6mer_kd_norm <- s6mer_kd*mm_kd/(s6mer_kd + mm_kd)
#     if (mirna %in% c("let-7a-21nt", "let-7a_plus1", "let-7a_minus1",
#                      "let-7a_miR-155")) {
#       if (equilibrium_nb) {
#         mirna_rand <- "let-7a-21nt"
#         experiment_rand <- "equilibrium_nb"
#       } else {
#         mirna_rand <- "let-7a"
#         experiment_rand <- "equilibrium"
#       }
#     } else if (mirna == "miR-155_let-7a") {
#       mirna_rand <- "miR-155"
#       experiment_rand <- "equilibrium"
#     } else {
#       mirna_rand <- mirna
#       experiment_rand <- "equilibrium"
#     }
#     if (mirna == "miR-1") {
#       buffer_rand <- TRUE
#       combined_rand <- FALSE
#     } else {
#       buffer_rand <- FALSE
#       combined_rand <- TRUE
#     }
#     segments(xmin, s8mer_kd, x1=xmax + 0.5, col=kSiteColors["8mer"])
#     segments(xmin, s6mer_kd, x1=xmax + 0.5, col=kSiteColors["6mer"])  
#     text(xmax + 1, s8mer_kd, labels="8mer", adj=c(0, 0.5), xpd=NA)
#     text(xmax + 1, s6mer_kd, labels="6mer", adj=c(0, 0.5), xpd=NA)
#     text(xmax + 1, ref_kd, labels="Seed m.m.", adj=c(0, 0.5), xpd=NA)
#   }
#   # Change the xmax value to allow for the complex legend to be plotted on the
#   # right-hand-side of the plot.
#   legend_names <- c()
#   # Variable names corresponding to the total number of sites plotted, and the
#   # the number of sites with binding affinity better than 1.5 fold better than
#   # the mismatch site on its own.
#   total_sites_plotted <- 0
#   total_sites_better <- 0
#   worse_sites <- c()
#   sapply(len_lim[1]:len_lim[2], function(kmer) {
#     ##################### Get all sites of length k #####################
#     target <- sprintf("^%smer-m[0-9]{1,2}\\.[0-9]{1,2}\\|.*\\|Comp_Kd", kmer)
#     inds <- grep(target, rownames(kds), perl=TRUE, value=TRUE)
#     if (kmer == 4) {
#       if (mirna == "miR-1") {
#         inds <- grep("^4mer-m15.*|^4mer-m19.*", inds, value=TRUE, invert=TRUE,
#                       perl=TRUE)
#       } else if (mirna == "miR-7-23nt") {
#         inds <- grep("^4mer-m17.*|^4mer-m20.*", inds, value=TRUE, invert=TRUE,
#                       perl=TRUE)
#       }
#     }
#     df_use <- kds[inds, ]
#     best_site <- rownames(df_use)[which.min(df_use[, 2])]
#     site_pos <- gsub("^.*mer-m(.*)\\..*$", replacement="\\1", best_site)
#     best_site <- gsub("^(.*)\\|.*\\|Comp_Kd", replacement="\\1", best_site)
#     legend_names <<- c(legend_names, best_site)
#     best_site_grep <- gsub("\\.", replacement="\\\\\\.", best_site)
#     best_inds <- grep(sprintf("^%s\\|.*\\|Comp_Kd", best_site_grep), rownames(df_use), perl=TRUE)
#     df_use <- df_use[best_inds, ]
#     pos <- as.integer(gsub("^(.*)\\|(.*)\\|(Comp_Kd$)", replacement="\\2", rownames(df_use)))
#     offsets <- pos - as.integer(site_pos)
#     inds_pos <- which(offsets >= xmin & offsets <= xmax)
#     offsets <- offsets[inds_pos]
#     df_use_pos <- df_use[inds_pos, ]
#     df_use_pos <- df_use_pos[order(offsets), ]
#     offsets_original <- sort(offsets)
#     offsets <- 20 - sort(offsets) - 8
#     # Places the asterisk by the cleavage site
#     if (mirna == "let-7a-21nt" & !justlegend) {
#       if ((experiment == "equil_c2_nb" & kmer == 10) |
#           (experiment == "equil_c_nb"  & kmer == 11)) {
#         ind_use <- which(offsets_original == 0)
#         y_txt <- df_use_pos[ind_use, 5]*1.1
#         text(x=20 - 8, y=y_txt, labels="*", adj=c(0.5, 0.5), cex=par("cex")*5)
#       }
#     }
#     if (mirna == "let-7a-21nt" & !justlegend) {
#       if (kmer %in% c(4, 10)) {
#         value_use <- min(df_use_pos[, 2])
#         value_ind <- which(df_use_pos[, 2] == value_use)
#         Arrows(x0=offsets[value_ind], y0=0.4, y1=0.2, x1=offsets[value_ind],
#                col=ConvertRColortoRGB(cols[as.character(kmer)], alpha=1),
#                arr.type="triangle", arr.adj=0, arr.length=0.035, arr.width=0.08)
#         adj_x_dict <- c(0.65, 0.45)
#         names(adj_x_dict) <- c("10", "4")
#         text(offsets[value_ind], 0.5,
#              labels=sprintf("+%s nt", 12 - offsets[value_ind]),
#              col=ConvertRColortoRGB(cols[as.character(kmer)], alpha=1),
#              adj=c(adj_x_dict[as.character(kmer)], 0))
#              # adj=c(0.5, 0))
#       }
#     }
#     if (kmer == 11 & !justlegend) {
#       lkds <- log(ref_kd/df_use[inds_pos, 2])
#       # message(sprintf("Mean benefit of kd at this length: %s", exp(mean(lkds))))
#       message(sprintf("Max benefit of log-kd at this length: %s", max(lkds)))
#       message(sprintf("Standard deviation of kd at this length: %s", sd(lkds)))
#       message(sprintf("(Standard deviation of logkd)/max(logkd) at this length: %1.1f%%", 100*sd(lkds)/max(lkds)))
#     }
#     if (!justlegend) {
#       total_sites_plotted <<- total_sites_plotted + length(inds_pos)
#       total_sites_better  <<- total_sites_better  + length(which(df_use[inds_pos, 2] < (ref_kd/text_bound)))
#       worse_sites <<- c(worse_sites, rownames(df_use)[inds_pos][which(df_use[inds_pos, 2] >= (ref_kd/text_bound))])
#       segments(x0=offsets + x_skew[as.character(kmer)], y0=df_use_pos[, 3],
#                y1=df_use_pos[, 5],
#                col=ConvertRColortoRGB(cols[as.character(kmer)], alpha=0.5), lwd=1, xpd=NA)
#       lines(offsets, df_use_pos[, 2], col=cols[as.character(kmer)], lwd=1, xpd=NA)
#       points(offsets, df_use_pos[, 2], col=cols[as.character(kmer)], pch=20)      
#     }
#   })
#   if (!justlegend) {
#     message(sprintf("Number of sites better than %s fold above background kd: %s",
#                     text_bound, total_sites_better/total_sites_plotted))
#     print(worse_sites)
#     ################# Add axes ##############################
#     alt_lab_pos <- c(-4, 0, 4, 8, 12, 16)
#     AddLinearAxis(1, tick.space=1, label.space=5, alt_lab=rev(alt_lab_pos),
#                   alt_lab_pos=alt_lab_pos, label="Offset (nt)", adj_pos=0.2)
#     AddLogAxis(2, label="Relative Kd")

#   }
#   ##### Make the legend-schematic on the right hand side of the plot. ##########
#   mirnalist <- unlist(strsplit(kMirnaSeqs[mirna], split=""))
#   if (width == 7.4) {
#     pos_left <- 18.5
#     pos_right <- pos_left + (40 - pos_left)*length(mirnalist)/21
#   } else {
#     pos_left <- 22.5
#     pos_right <- pos_left + (37 - pos_left)*length(mirnalist)/21
#   }
#   mirnalist_rev <- c(rev(mirnalist), "-5'")
#   mir_pos_x <- seq(pos_left, pos_right, length.out=length(mirnalist))
#   del_mir_pos_x <- mir_pos_x[length(mirnalist)] - mir_pos_x[length(mirnalist) - 1]
#   mir_pos_x <- c(mir_pos_x, mir_pos_x[length(mirnalist)] + 1.5*del_mir_pos_x)
#   if (height == 3.7*9/9.75) {
#     foldchange_y <- 1.92
#     num_dip <- 1.74
#     y_mir_site <- ymin*3.1

#   } else if (!justlegend) {
#     foldchange_y <- 1.8
#     num_dip <- 1.7
#   } else {
#     # These get used in Figure 7 and Figure7-figure supplements.
#     y_mir_site <- ymin*1.5

#     # foldchange_y <- 1.8
#     foldchange_y <- 1.4

#     # num_dip <- 1.6
#     num_dip <- 1.3
#   }
#   if (corrected_kds & !justlegend) {
#     foldchange_y <- foldchange_y*log(1/ymin)/log(1/0.0002)
#   }
#   y_base <- y_mir_site
#   if (justlegend) {
#     y_dip <- y_mir_site/1.15
#   } else if (height == 3.7*9/9.75) {
#     y_dip <- y_mir_site/1.2
#   } else {
#     y_dip <- y_mir_site/1.25
#   }
#   mir_pos_y <- c(rep(y_base, length(mirnalist) - 1), rep(y_dip, 2))
#   cols_mirna <- rep(c("black", "black", "#F7931D", "black", "#ED1C24", "black"),
#               times=c(1, length(mirnalist) - 17, 4, 5, 6, 2))
#   text(x=mir_pos_x, y=mir_pos_y, labels=mirnalist_rev, xpd=NA,
#        adj=c(0.5, 0.5), col=cols_mirna)
#   mir_nums <- c(length(mirnalist), 16, 13, 8, 7, 6, 5, 4, 3, 2, 1)
#   mir_nums_x <- mir_pos_x[length(mirnalist) + 1 - mir_nums]
#   mir_nums_y <- (mir_pos_y[length(mirnalist) + 1 - mir_nums])/num_dip
#   text(x=mir_nums_x, y=mir_nums_y, labels=mir_nums, xpd=NA,
#        adj=c(0.5, 0.5), col=c("black", "#F7931D", "#F7931D", "black",
#                               rep("#ED1C24", 6), "black"))
#   dot_pos <- setdiff(mir_pos_x[-length(mir_pos_x)], mir_nums_x)
#   points(x=dot_pos, y=rep(mir_nums_y[1], length(dot_pos)), pch=20, cex=0.5,
#          col=rep(c("black", "#F7931D", "black"),
#                  times=c(length(mirnalist) - 17, 2, 4)),
#          xpd=NA)
#   y_val <- mir_pos_y[1]
#   print(legend_names)
#   for (name_i in rev(legend_names)) {
#     y_val <- y_val*foldchange_y
#     print(foldchange_y)
#     name_split <- unlist(strsplit(name_i, split="mer-m"))
#     len_i <- name_split[1]
#     name_split <- name_split[2]
#     name_split_2 <- unlist(strsplit(name_split, split="\\|"))[1]
#     start_stop <- as.integer(unlist(strsplit(name_split_2, split="\\.")))
#     start_positions <- mir_pos_x[length(mirnalist) + 1 - start_stop[1]:start_stop[2]]
#     segments(x0=start_positions[1], y0=y_val,
#              x1=start_positions[length(start_positions)], lwd=1.5,
#              col=cols[len_i])
#     points(x=c(start_positions[1], start_positions[length(start_positions)]),
#            y=rep(y_val, 2), pch=19, col=cols[len_i], cex=1.4)
#     dot_pos <- setdiff(mir_pos_x[1:(length(mirnalist) - 8)], start_positions)
#     points(x=dot_pos, y=rep(y_val, length(dot_pos)), col="gray", pch=20,
#            cex=0.3)
#     text(x=mir_pos_x[length(mirnalist) - 6], y=y_val,
#          labels=sprintf("%s-%s", start_stop[1], start_stop[2]), adj=c(0, 0.5))
#   }
#   if (!justlegend) {
#     text(x=mean(mir_pos_x), y=y_val*foldchange_y^2,
#          labels="Range of 3'\npairing")
#   }
#   # Add text indicating which miRNA is being plotted.
#   if (mirna == "let-7a-21nt") {
#     if (experiment == "equil_c_nb") {
#       mirna_txt <- "let-7a rep."
#     } else {
#       mirna_txt <- "let-7a"
#     }
#   } else if (mirna == "let-7a_minus1") {
#     mirna_txt <- "let-7a(-1)"
#   } else if (mirna == "let-7a_plus1") {
#     mirna_txt <- "let-7a(+1)"
#   } else if (mirna == "miR-7-23nt") {
#     mirna_txt <- "miR-7"
#   } else {
#     mirna_txt <- mirna
#   }
#   if (justlegend) {
#     str_line <- -1.5
#     str_at <- xmin + (xmax - xmin)*0.95
#     mtext(text=mirna_txt, side=3, line=str_line, at=str_at, adj=1,
#           cex=par("cex"))
#   } else {
#     text(x=mir_nums_x[length(mir_nums_x)] + strwidth("U-5'")*0.5, y=y_base/(num_dip^2.5), labels=mirna_txt,
#          adj=c(1, 1), xpd=NA)
#     xy <- GetPlotFractionalCoords(0, 1, log="y")
#     print(xy)
#     text(x=xy[1] + 0.5, y=xy[2], labels=mirna_txt,
#          adj=c(0, 1), xpd=NA)
#   }



#   line2user <- function(line, side) {
#     lh <- par('cin')[2] * par('cex') * par('lheight')
#     x_off <- diff(grconvertX(c(0, lh), 'inches', 'npc'))
#     y_off <- diff(grconvertY(c(0, lh), 'inches', 'npc'))
#     switch(side,
#            `1` = grconvertY(-line * y_off, 'npc', 'user'),
#            `2` = grconvertX(-line * x_off, 'npc', 'user'),
#            `3` = grconvertY(1 + line * y_off, 'npc', 'user'),
#            `4` = grconvertX(1 + line * x_off, 'npc', 'user'),
#            stop("Side must be 1, 2, 3, or 4", call.=FALSE))
#   }
#   # abline(h=line2user(0:-10, 1), lty=3, xpd=TRUE)
#   # segments(x0=xmin, y0=ymin, x1=xmax, xpd=NA)

#   # segments(x0=xmin, y0=y_first_site, x1=xmax, xpd=NA)

#   # segments(x0=xmin, y0=y_top_site, x1=xmax, xpd=NA)

#   # segments(x0=xmin, y0=ymax_adjusted, x1=xmax, lwd=4, xpd=NA)



#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }

# # 2C____________________________________________________________________________
# PlotPairingMatrix <- function(
#   mirna, experiment, offset, n_constant=3, sitelist="progthrp_suppcomp",
#   model_values=FALSE, makeglobalmodel=TRUE, offset_lim=c(-4, 16),
#   len_lim=c(4, 11), pos3p_lim=c(9, 12), corrected_kds=FALSE, sumseed=FALSE,
#   F_method=FALSE, exponential=FALSE, intercept=FALSE, fixed_offset=0, log_plus_one=FALSE,
#   supp_base=FALSE, offset_base=FALSE, site_base=NULL, residual=FALSE, key=FALSE,
#   xlabels=TRUE, suppress_label=FALSE, label_offset=FALSE, mirna_label=FALSE,
#   extralabel=FALSE, loop=FALSE, alt_top=FALSE, max_value_check=NULL,
#   height=2.8, width=2.47, xpos=20, ypos=20, pdf.plot=FALSE
# ) {
#   # Load the data matrix.
#   print("in pairing range matrix function")
#   if (model_values) {
#     if (makeglobalmodel) {
#       if (sitelist == "progthrp") {
#         model <<- SubfunctionCall(FitPairingOffsetAndMismatchModelWithError)
#       } else {
#         model <<- SubfunctionCall(FitPairingAndOffsetModelWithError)
#       }
#     }
#     offset_name <- as.character(offset)
#     if (log_plus_one) {
#       R_mat <- log10(exp(t(model$pairing$MLE) + model$offsets[offset_name, 1]) + 1)
#     } else {
#       R_mat <- t(model$pairing$MLE)
#       R_mat <- R_mat*model$offsets[offset_name, 1]
#       R_mat_use_temp_global <<- R_mat
#     }
#     # This just makes anything with data NA, even if there is a model value for
#     # it.
#     R_mat_data <- t(SubfunctionCall(MakePairingMatrix))
#     R_mat_data_use <- R_mat_data
#     R_mat_data_use_global <<- R_mat_data_use
#     R_mat_data <- R_mat_data*0 + 1
#     R_mat <- R_mat*(R_mat_data*0 + 1)
#     R_mat_global <<- R_mat
#     if (residual) {
#       R_mat <- R_mat_data - R_mat
#     }
#   } else {
#     R_mat <- t(SubfunctionCall(MakePairingMatrix))
#     R_mat_global_data <<- R_mat
#   }
#   mar1 <- 2.2
#   mar2 <- 0.5
#   mar3 <- 1.8
#   if(key) {
#     width <- 4
#     mar4 <- 9.23
#   } else {
#     mar4 <- 2.5
#   }
#   R_mat <<- R_mat
#   # This is the conversion factor in order to have the correct height of the
#   # plot for the specified width, that makes each box of equal height and width.
#   mir_length <- nchar(kMirnaSeqs[mirna])
#   SubfunctionCall(FigureSaveFile2)
#   par(mar=c(mar1, mar2, mar3, mar4))
#   xmin <- 0
#   xmax <- nrow(R_mat)
#   ymin <- 0
#   ymax <- ncol(R_mat) 
#   BlankPlot()
#   # Assign the positions of the corners of the boxes.
#   xlefts <- rep(seq(0, ncol(R_mat) - 1), nrow(R_mat))
#   xright <- xlefts + 1
#   ybottom <- rep(rev(seq(ymax - 1, ymin)), each=ncol(R_mat))
#   ytop <- ybottom + 1
#   # Make the x-axis
#   y_axis_pos <- seq(ymin, ymax - 1) + 0.5
#   y_axis_labs <- rownames(R_mat)
#   for (i in 1:length(y_axis_labs)) {
#     if (i %% 2 == 0) {
#       y_axis_labs[i] <- ""
#     }
#   }
#   AddLinearAxis(1, 1, 1, label="3'-paired nt",
#                 label_pos_ticks=TRUE,
#                 alt_lab=colnames(R_mat),
#                 alt_lab_y_dist=0.01,
#                 alt_lab_pos=xmin:(xmax - 1) + 0.5,
#                 alt_tick_pos=TRUE)
#   AddLinearAxis(4, 1, 2, label="5'-paired nt",
#                 alt_lab=y_axis_labs,
#                 alt_lab_pos=y_axis_pos,
#                 alt_tick_pos=TRUE,
#                 line=0.8)
#   # Add the label for the seed if not plotting relative to seed kds.
#   col.inds <- floor((R_mat)/log10(700)*99 + 1)
#   x_lab_pos <- -0.75
#   pos_5p <- as.integer(rep(rownames(R_mat), each=ncol(R_mat)))
#   pos_3p <- as.integer(rep(colnames(R_mat), nrow(R_mat)))
#   impossible_cols <- which(pos_5p > pos_3p)

#   start_col <- 0
#   r_col <- 0.4
#   # Make distribution of 100 colors and a 101th gray color for NA indeces.
#   color.dist <- c(rev(cubeHelix(100, start=start_col, r=r_col, hue=0.8)),
#                   "gray90", "white")
#   color.dist[100] <- kMaxValueColor
#   # Make the color index scale.
#   # col.inds <- sapply(t(col.inds), function(col.ind) {
#   #   min(max(1, col.ind), 100)
#   # })
#   col.inds <- t(col.inds)
#   col.inds[which(col.inds > 100)] <- 100
#   col.inds[which(col.inds < 1)] <- 1
#   col.inds[which(is.na(col.inds))] <- 101
#   col.inds[impossible_cols] <- 102
#   # Make the color rectangle.
#   rect(xlefts, ybottom, xright, ytop, col=color.dist[col.inds], lwd=0.4, 
#        xpd=NA, border="white")
#   # Make the key.
#   if (key) {
#     mir_scale <- (mir_length - 11)/(21 - 11)
#     if (alt_top) {
#       max_value_use <- max(max_value_check, max(10^R_mat, na.rm=TRUE))
#       if (max_value_use > 700) {
#         labels <- c(1, 3, 10, 30, 100, 300, 700, 100*round(max_value_use/100))
#       } else {
#         labels <- c(1, 3, 10, 30, 100, 300, 700)
#       }
#     } else {
#       labels <- c(1, 3, 10, 30, 100, 300, 700)
#     }
#     # if (alt_top) {
#       print(labels)
#       num_boxes <- ceiling(log(labels[length(labels)])/log(700)*100)
#     # } else {
#       # num_boxes <- 100
#     # }
#     print(num_boxes)
#     color.dist <- c(rev(cubeHelix(100, start=start_col, r=r_col, hue=0.8)),
#                     rep(kMaxValueColor, num_boxes - 100))
#     y_div <- (ymax - ymin)/num_boxes
#     kpl <- xmax + 5*mir_scale # Left-hand position of the key
#     kw <- mir_scale                # Width of the key
#     rect(xleft=kpl,
#          ybottom=seq(num_boxes)*y_div - y_div,
#          xright=kpl + kw,
#          ytop=seq(num_boxes)*y_div,
#          col=color.dist, xpd=NA, border=NA)
#     # Generate the axis for the legend and the label
#     bottom <- ymin + y_div/2
#     top <- ymax - y_div/2
#     # if (kdrel) {
#     pos_labels <- log(labels)
#     centered_labels <- pos_labels - pos_labels[1]
#     norm_labels <- centered_labels/(centered_labels[length(centered_labels)])
#     height_span <- top - bottom
#     pos_labels <- norm_labels*height_span + y_div/2
#     labels <- as.character(labels)
#     labels[1] <- expression(""<=1)
#     axis(4, at=pos_labels, labels=labels, lwd=par()$lwd, pos=kpl + kw, xpd=NA)
#     text(x=kpl + kw + 3*mir_scale, y=ymax/2,
#          labels=bquote(italic(K)[D]*.(" fold change")),
#          srt=270, xpd=NA)
#   }  
#   # # Add the label indicating how many nucleotides of pairing.
#   if (model_values) {
#     if (label_offset) {
#       mtext(text=sprintf("%s-nt offset; model", offset), side=3, at=xmax, adj=1,
#             cex=par("cex"))    
#     } else {
#       mtext(text="Model-estimated values", side=3,
#             at=xmin, adj=0, cex=par("cex"))      
#     }
#   } else if (!(suppress_label)) {
#     if (loop) {
#       unit <- "loop"
#     } else {
#       unit <- "offset"
#     }
#     mtext(text=sprintf("%s-nt %s", offset, unit), side=3,
#           at=xmax, adj=1, cex=par("cex"))    
#   }
#   # If the `mirna_label` conditional is true, add the label saying which miRNA
#   # is being looked at.
#   if (mirna_label) {
#     if (mirna == "let-7a-21nt") {
#       mirna_txt <- "let-7a"
#     } else if (mirna == "let-7a_minus1") {
#       mirna_txt <- "let-7a(-1)"
#     } else if (mirna == "let-7a_plus1") {
#       mirna_txt <- "let-7a(+1)"
#     } else if (mirna == "let-7a_miR-155") {
#       mirna_txt <- "let-7a-miR-155"
#     } else if (mirna == "miR-155_let-7a") {
#       mirna_txt <- "miR-155-let-7a"
#     } else {
#       mirna_txt <- mirna
#     }
#     x <- GetPlotFractionalCoords(0.125, 0.5)[1]
#     mtext(text=mirna_txt, side=3, line=0.8, at=x, adj=0, cex=par("cex"))
#   }

#   if (extralabel) {
#     str.end3prand <- "r"
#     if (collapsemm) {
#       str.collapsemm <- "sum"
#     } else {
#       str.collapsemm <- "geo"
#     }
#     mtext(text=sprintf("%s_%s", str.end3prand, str.collapsemm), side=3,
#           line=0.7, at=xmax, adj=1, cex=par("cex"))    
#   }

#   # Finish the plot.
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
#   return(R_mat)
# }

# ################################################################################
# # FIGURE S2
# ################################################################################

# # S2C___________________________________________________________________________
# PlotPositionalCanonicalSites <- function(
#   mirna="let-7a-21nt", experiment="equil_c2_nb", n_constant=3,
#   sitelist="progthrp_suppcomp", sumseed=FALSE, corrected_kds=TRUE, kd_fc=FALSE,
#   height=4, width=4, xpos=20, ypos=20, pdf.plot=FALSE
# ) {
#   ########################### Load the Kds #####################################
#   if (corrected_kds) {
#     kds <- SubfunctionCall(ApplyKdCorrection, prog_sitelist=sitelist,
#                            prog_n_constant=n_constant,
#                            rand_n_constant=n_constant)
#   } else {
#     kds <- SubfunctionCall(EquilPars, sitelist="progthrp_suppcomp")
#   }
#   ######################### Set up the plot window #############################
#   SubfunctionCall(FigureSaveFile2)
#   if (corrected_kds) {
#     ymin <- 1e-4
#   } else {
#     ymin <- 1e-3
#   }
#   xmin <- 5
#   xmax <- 25
#   ymax <- 1
#   BlankPlot(log="y", inv="x")
#   # cols <- kSiteColors[kSeedSites]

#   x_skew <- seq(-0.1, 0.1, length.out=length(kSeedSites))
#   names(x_skew) <- kSeedSites
#   ################### Put reference line ###############################
#   if (sitelist == "progthrp_suppcomp" & sumseed) {
#     segments(x0=xmin, kds["Comp_Kd", 2], x1=xmax, lty=2, lwd=0.5)
#   } else {
#     segments(x0=xmin, GeoMean(kds[grep("^8mer-mm[ACTG][2-7]_Kd$",
#                                        rownames(kds),
#                                        perl=TRUE), 2]), x1=xmax, lty=2, lwd=0.5)    
#   }
#   ################### Add the points and lines for each site ###################
#   sapply(kSeedSites, function(can_site) {
#     grep_str <- sprintf("^%s\\|.*\\|Comp_Kd", can_site)
#     inds <- grep(grep_str, rownames(kds), perl=TRUE, value=TRUE)
#     df_use <- kds[inds, ]
#     pos <- as.integer(gsub("^(.*)\\|(.*)\\|(Comp_Kd$)", replacement="\\2", rownames(df_use)))
#     inds <- which(pos >= xmin & pos <= xmax)
#     x <- pos[inds]
#     col_use <- kSiteColors[can_site]
#     segments(x0=x + x_skew[can_site], y0=df_use[inds, 3], y1=df_use[inds, 5],
#              col=col_use, lwd=1, xpd=NA)
#     lines(x, df_use[inds, 2], col=col_use, lwd=1, xpd=NA)
#     points(x, df_use[inds, 2], col=col_use, pch=20, xpd=NA)
#   })
#   ######################### Add the axes #############################
#   AddLinearAxis(1, 1, 5, label="Position")
#   AddLogAxis(2, label="Relative Kd")
#   xy <- GetPlotFractionalCoords(0.075, 1.175, log="y", inv="x")
#   ###################### Add the legend ################################
#   Legend(xy, legend=kSeedSites, col=kSiteColors[kSeedSites], ncol=2)
#   ######################### Add the miRNA label ################################
#   if (mirna == "let-7a-21nt") {
#     if (experiment == "equil_c_nb") {
#       mirna_txt <- "let-7a rep."
#     } else {
#       mirna_txt <- "let-7a"
#     }
#   } else if (mirna == "let-7a_minus1") {
#     mirna_txt <- "let-7a(-1)"
#   } else if (mirna == "let-7a_plus1") {
#     mirna_txt <- "let-7a(+1)"
#   } else {
#     mirna_txt <- mirna
#   }
#   mtext(text=mirna_txt, side=1, line=-1.5, at=xmin, adj=1, cex=par("cex"))
#   ########################## Finish the plot ###################################
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }

# # S2D-K, S3A&B__________________________________________________________________
# PlotOneThreePrimeSite <- function(
#   kmer, mirna, experiment, n_constant=3, sitelist="progthrp_suppcomp",
#   sumseed=FALSE, corrected_kds=TRUE, offset_lim=c(-4, 16), equilibrium_nb=FALSE,
#   mirna_label=FALSE, height=3.7, width=3.7, xpos=20, ypos=20, pdf.plot=FALSE
# ) {
#   if (corrected_kds) {
#     kds <- SubfunctionCall(ApplyKdCorrection, prog_n_constant=n_constant,
#                            rand_n_constant=n_constant,
#                            prog_sitelist="progthrp_suppcomp")
#     if (mirna == "miR-155" & kmer >= 8) ymin <- 1e-5
#     else                                ymin <- 1e-4

#   } else {
#     if (mirna == "miR-1" & experiment == "equilibrium") {
#       combined <- FALSE
#       buffer <- TRUE
#     } else if (mirna == "miR-7-23nt" & experiment == "equilibrium2_nb") {
#       combined <- FALSE
#       buffer <- FALSE
#     } else {
#       combined <- TRUE
#       buffer <- FALSE
#     }
#     kds <- SubfunctionCall(EquilPars)
#     ymin <- 3e-4
#   }
#   SubfunctionCall(FigureSaveFile2)
#   xmin <- offset_lim[1]
#   xmax <- offset_lim[2]
#   ymax <- 1
#   par(mar=c(3, 3, 1, 1))
#   BlankPlot(log="y", inv=kOffsetInv)
#   if (sitelist == "progthrp_suppcomp" & sumseed) {
#     print(kds["Comp_Kd", 2])
#     segments(x0=xmin, kds["Comp_Kd", 2], x1=xmax, lty=2, lwd=0.5)
#   } else {
#     print(kds[grep("^8mer-mm[ACTG][2-7]_Kd$",
#                                        rownames(kds),
#                                        perl=TRUE), 2])
#     segments(x0=xmin, GeoMean(kds[grep("^8mer-mm[ACTG][2-7]_Kd$",
#                                        rownames(kds),
#                                        perl=TRUE), 2]), x1=xmax, lty=2, lwd=0.5)    
#   }
#   grep_str <- sprintf("^%smer-m[0-9]{1,2}\\.[0-9]{1,2}\\|.*\\|Comp_Kd", kmer)
#   inds <- grep(grep_str, rownames(kds), perl=TRUE, value=TRUE)
#   kds <- kds[inds, ]
#   # Get all the thrp site names.
#   sites_thrp <- gsub("^(.*)\\|.*\\|.*$", replacement="\\1", rownames(kds))
#   unique_sites_thrp <- rev(unique(sites_thrp))
#   pos_1 <- as.integer(gsub("^.*mer-m(.*)\\..*$", replacement="\\1",
#                            unique_sites_thrp))
#   unique_sites_thrp <- unique_sites_thrp[order(pos_1)]
#   len_sites_thrp <- length(unique_sites_thrp)
#   x_skew <- seq(-0.1, 0.1, length.out=len_sites_thrp)
#   cols <- GetColors(n=12, scheme="smooth rainbow",
#                     stops=c(0.05, 1))[1:len_sites_thrp]
#   names(x_skew) <- unique_sites_thrp
#   names(cols) <- unique_sites_thrp
#   # Operate over each site.
#   # legend_names <- unique_sites_thrp
#   print(unique_sites_thrp)
#   if (mirna == "miR-1") {
#     sites_remove <- c("4mer-m15.18", "4mer-m19.22")
#   } else if (mirna == "miR-7-23nt") {
#     sites_remove <- c("4mer-m17.20", "4mer-m20.23")
#   } else {
#     sites_remove <- c()
#   }
#   unique_sites_thrp <- setdiff(unique_sites_thrp, sites_remove)
#   cols <- cols[unique_sites_thrp]
#   sapply(unique_sites_thrp, function(site) {
#     message("___________")
#     message(site)

#     grep_str <- sprintf("^%smer-m[0-9]{1,2}\\.[0-9]{1,2}\\|.*\\|Comp_Kd", kmer)
#     inds <- grep(site, rownames(kds), value=TRUE)
#     df_use <- kds[inds, ]
#     # Get the site position 5-prime most position
#     site_pos <- gsub("^.*mer-m(.*)\\..*$", replacement="\\1", site)
#     # Get the starting point of that site at each position in order to calculate
#     # the offset.
#     pos <- as.integer(gsub("^(.*)\\|(.*)\\|(Comp_Kd$)", replacement="\\2",
#                       rownames(df_use)))
#     offsets <- pos - as.integer(site_pos)
#     inds_pos <- which(offsets >= xmin & offsets <= xmax)
#     offsets <- offsets[inds_pos]
#     df_use_pos <- df_use[inds_pos, ]
#     # This make sure that the offsets are in order, because they are not in the
#     # random libraries, due to how the scripting is different.
#     df_use_pos <- df_use_pos[order(offsets), ]
#     offsets <- sort(offsets)
#     # Make the error bars, lines, and points for the figure.
#     # NOTE: The error bars have been given a 50% transparency.
#     segments(x0=offsets + x_skew[site], y0=df_use_pos[, 3], y1=df_use_pos[, 5],
#              col=ConvertRColortoRGB(cols[site], alpha=0.5), lwd=1, xpd=NA)
#     lines(offsets, df_use_pos[, 2], col=cols[site], lwd=1)
#     points(offsets, df_use_pos[, 2], col=cols[site], pch=20)
#     if (mirna == "let-7a-21nt") {
#       if (kmer >= 8 & as.integer(site_pos) == 9) {
#         ind_use <- which(offsets == 0)
#         y_txt <- df_use_pos[ind_use, 5]*1.1
#         text(x=0, y=y_txt, labels="*", adj=c(0.5, 0.5), cex=par("cex")*5)
#       }
#       if (kmer >= 6 & as.integer(site_pos) %in% c(11, 12)) {
#         x_use <- c(0.5, 3)
#         y0_use <- c(0.45, 0.25)
#         labels <- c("+1-0 nt", "+3 nt")

#         names(x_use) <- c("12", "11")
#         names(y0_use) <- c("12", "11")
#         names(labels) <- c("12", "11")
#         Arrows(x0=x_use[as.character(site_pos)], y0=y0_use[as.character(site_pos)], y1=0.2,
#                x1=x_use[as.character(site_pos)], col=cols[site],
#                arr.type="triangle", arr.adj=0, arr.length=0.035, arr.width=0.08)
#         adj_x_dict <- c(0.45, 0.65)
#         names(adj_x_dict) <- names(x_use)
#         text(x_use[as.character(site_pos)],
#              y0_use[as.character(site_pos)]*(6/5),
#              labels=labels[as.character(site_pos)],
#              col=cols[site],
#              adj=c(adj_x_dict[as.character(site_pos)], 0))
#       }
#     }
     
#   })

#   alt_lab_pos <- c(-4, 0, 4, 8, 12, 16)
#   AddLinearAxis(1, tick.space=1, label.space=5, alt_lab_pos=alt_lab_pos,
#                 label="Offset (nt)")
#   AddLogAxis(2, label="Relative Kd")

#   # Add title to the plot.
#   xy <- GetPlotFractionalCoords(0.025, 1, log="y", inv=kOffsetInv)
#   text(xy[1], xy[2], sprintf("%s bp of 3' pairing", kmer), adj=c(0, 0), xpd=NA)
#   xy <- GetPlotFractionalCoords(0.5, 0, log="y")
#   legend_names <- gsub(sprintf("^%smer-m(.*)\\.(.*)$", kmer), replacement="\\1-\\2",
#                        unique_sites_thrp, perl=TRUE)
#   Legend(xy, legend=legend_names, col=cols, ncol=3, xjust=0.5, yjust=0)
#   # Add the mirna label to the top of the plot.
#   if (mirna_label) {
#     if (mirna == "let-7a-21nt") {
#       if (experiment == "equil_c_nb") {
#         mirna_txt <- "let-7a rep."
#       } else {
#         mirna_txt <- "let-7a"
#       }
#     } else if (mirna == "let-7a_minus1") {
#       mirna_txt <- "let-7a(-1)"
#     } else if (mirna == "let-7a_plus1") {
#       mirna_txt <- "let-7a(+1)"
#     } else if (mirna == "miR-7-23nt") {
#       mirna_txt <- "miR-7"
#     } else {
#       mirna_txt <- mirna
#     }
#     xy <- GetPlotFractionalCoords(1, 1, log="y", inv=kOffsetInv)
#     text(xy[1], xy[2], mirna_txt, adj=c(1, 0), xpd=NA)
#     # mtext(text=mirna_txt, side=3, line=0, at=xmin, adj=1, cex=par("cex"))
#   }
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }

# ################################################################################
# # FIGURE 3
# ################################################################################
# PlotMeanReporterRepression <- function(mirna="let-7a",
#   experiment="twist_reporter_assay_3p_2_tp", exp_type="parallel",
#   dual_site=FALSE, all_lin41_UTR=FALSE, height=3.5, xpos=20,
#   ypos=20, pdf.plot=FALSE
# ) {
#   if (mirna == "let-7a") {
#     log2fc_df_1 <- SubfunctionCall(
#       GetThreePrimeReporterFoldChanges, mirna="let-7a", rep=1, aggregate=TRUE
#     )
#     log2fc_df_2 <- SubfunctionCall(
#       GetThreePrimeReporterFoldChanges, mirna="let-7a", rep=2, aggregate=TRUE
#     )
#     log2fc_df_3 <- SubfunctionCall(
#       GetThreePrimeReporterFoldChanges, mirna="let-7a-21nt", rep=1,
#       aggregate=TRUE
#     )
#     log2fc_df_4 <- SubfunctionCall(
#       GetThreePrimeReporterFoldChanges, mirna="let-7a-21nt", rep=2,
#       aggregate=TRUE
#     )
#     log2fc_df <- log2fc_df_1[, 1:4]
#     # Add each of the replicates to the dataframe.
#     log2fc_df$log2fc_1 <- log2fc_df_1$log_fc
#     log2fc_df$log2fc_2 <- log2fc_df_2$log_fc
#     log2fc_df$log2fc_3 <- log2fc_df_3$log_fc
#     log2fc_df$log2fc_4 <- log2fc_df_4$log_fc
#   } else if (mirna == "miR-1") {
#     log2fc_df_1 <- SubfunctionCall(
#       GetThreePrimeReporterFoldChanges, mirna="miR-1", rep=1, aggregate=TRUE
#     )
#     log2fc_df_2 <- SubfunctionCall(
#       GetThreePrimeReporterFoldChanges, mirna="miR-1", rep=2, aggregate=TRUE
#     )
#     log2fc_df <- log2fc_df_1[, 1:4]
#     # Add each of the replicates to the dataframe.
#     log2fc_df$log2fc_1 <- log2fc_df_1$log_fc
#     log2fc_df$log2fc_2 <- log2fc_df_2$log_fc   
#   }
#   if (all_lin41_UTR & dual_site)  width <- 11.46
#   else if (all_lin41_UTR)         width <- 11
#   else if (dual_site)             width <- 9.16
#   else                            width <- 8.8 
#   # Boolean in which a left-hand shift is given if both the total repression and
#   # the lin-41 UTR-specific repression is to be plotted.
#   if (all_lin41_UTR)  point_shift <- 0.2
#   else                point_shift <- 0

#   if (dual_site) x_denom <- 22
#   else                x_denom <- 21
#   # Boolean about where the right-hand side of lines ends.
#   if (all_lin41_UTR) line_r <- 0
#   else               line_r <- 0.25/x_denom

#   if (all_lin41_UTR) {
#     output_p_values <- data.frame(site=c(), p_vals=c())
#   }
#   SubfunctionCall(FigureSaveFile2)
#   xmin <- 0
#   xmax <- 1
#   ymin <- -3
#   ymax <- 0
#   mar_use <- c(1, 5, 6.5, 1)
#   par(mar=mar_use)
#   BlankPlot()
#   # Constants defining the vertical heights of the information above the points.
#   y_const <- abs(ymin)*0.6
#   y_seed <- 0.74*y_const
#   y_line_1 <- 0.6*y_const
#   y_thrp <- 0.40*y_const
#   y_line_2 <- 0.32*y_const
#   y_bulge <- 0.18*y_const
#   segments(x0=xmin, y0=0, x1=xmax, col="gray")
#   if (dual_site)  config_str <- "dual"
#   else            config_str <- "single"
#   x_start <- 0.5/x_denom
#   base_sites <- c("None", "8mer", "7mer-m8", "7mer-A1", "6mer", "8mer-w6")
#   if (dual_site) {
#     log2fc_df <- log2fc_df[which(log2fc_df$Seed == "None" | log2fc_df$Seed == "None_lin41_UTR" | log2fc_df$Location == "dual"), ]
#   } else {
#     log2fc_df <- log2fc_df[which(log2fc_df$Seed == "None" | log2fc_df$Location != "dual"), ]
#   }
#   print(log2fc_df)
#   for (base_site in base_sites) {
#     inds <- which(log2fc_df$Seed == base_site & log2fc_df$ThreePrime == "")
#     log2fc_df_use <- log2fc_df[inds, 5:ncol(log2fc_df)]
#     # Calculate background subtraction. ________________________________________
#     if (base_site == "None")  bg_subtract <- mean(unlist(log2fc_df_use))
#     x <- rep(x_start, length(unlist(log2fc_df_use)))
#     y <- unlist(log2fc_df_use) - bg_subtract
#     # Plot the points associated with the site. ________________________________
#     points(x - point_shift/x_denom, y, col=kSiteColors[base_site], lwd=1, pch=1, xpd=NA)
#     # Plot the horizontal line giving the mean repression. _____________________
#     segments(x0=x[1] - 0.25/x_denom - point_shift/x_denom, x1=x[1] + line_r, y0=mean(y), xpd=NA)
#     if (all_lin41_UTR) {
#       y_old <- y
#       mean_old <- mean(y)
#       inds <- which(log2fc_df$Seed == paste0(base_site, "_lin41_UTR") & log2fc_df$ThreePrime == "")
#       log2fc_df_use <- log2fc_df[inds, 5:ncol(log2fc_df)]
#       x <- rep(x_start, length(unlist(log2fc_df_use)))
#       y <- unlist(log2fc_df_use) - bg_subtract
#       y_difs <- y - y_old
#       se_difference <- sd(y_difs)/sqrt(length(y_difs))
#       t_statistic <- mean(y_difs)/se_difference
#       p_val <- 2*pt(t_statistic, df=length(y) - 1)
#       # This versin of it does it treating the reps as unpaired, which is more
#       # similar to how the Tukey tests were done.
#       p_val <- t.test(y, y_old, alternative="two.sided", var.equal=FALSE)$p.value
#       p_value_row <- data.frame(site=sprintf("%s_%s", base_site, config_str),
#                                 p_vals=p_val)
#       output_p_values <- rbind(output_p_values, p_value_row)
#       colnames(output_p_values) <- c("site", "p_vals")
#       points(x + point_shift/x_denom, y, col=kSiteColors[base_site], lwd=1, pch=1, xpd=NA)
#       segments(x0=x[1], x1=x[1] + 0.25/x_denom + point_shift/x_denom, y0=mean(y), xpd=NA)
#       segments(x0=x[1], x1=x[1], y0=mean_old, y1=mean(y), xpd=NA)
#     }

#     # Label the information associated with this site. _________________________
#     if (base_site != "8mer-w6") {
#     	if (base_site == "None")	base_site <- "No site"
#       text(x_start, y_seed, label=base_site, adj=0, srt=45, xpd=NA)
#     }
#     text(x=x_start, y=y_bulge, label="-", xpd=NA)
#     text(x=x_start, y=y_thrp, label="-", xpd=NA, adj=c(0.5, 0))
#     # Update the x-axis value. _________________________________________________
#     x_start <- x_start + 1/x_denom
#   }
#   ########################## Plot the 3-prime sites ############################
#   thrp_sites <- c("4mer-m13.16", "9mer-m11.19", "9mer-m13.21")
#   bulges <- c("none", "A", "T", "AAAA", "TTTT")
#   # Plot horizontal line above plot.
#   segments(x0=x_start - 1.4/x_denom, y0=y_line_1,
#            x1=x_start + 14.4/x_denom, xpd=NA)
#   # Label the base site. _______________________________________________________
#   text(x=x_start + 6.5/x_denom, y=y_seed, label="8mer-w6", xpd=NA)
#   for (thrp_site in thrp_sites) {
#     # Label the pairing range. _________________________________________________
#     text(x_start + 2/x_denom, y_thrp,
#          label=gsub("^.*mer-m(.*)\\.(.*)$", replacement="\\1-\\2",
#                     thrp_site, perl=TRUE),
#          srt=0, xpd=NA, adj=c(0.5, 0))
#     # Plot the line segment above the pairing range. ___________________________
#     segments(x0=x_start - 0.4/x_denom, y0=y_line_2,
#              x1=x_start + 4.4/x_denom, xpd=NA)
#     for (bulge in bulges) {
#       if (bulge == "none")  bulge_use <- ""
#       else                  bulge_use <- bulge
#       # Subset the log2fc matrix. ______________________________________________
#       inds <- which((log2fc_df[, "Seed"] == "8mer-w6") &
#                     (log2fc_df[, "ThreePrime"] == thrp_site) &
#                     (log2fc_df[, "Bulge"] == bulge_use))
#       log2fc_df_use <- log2fc_df[inds, 5:ncol(log2fc_df)]
#       x <- rep(x_start, length(unlist(log2fc_df_use)))
#       y <- unlist(log2fc_df_use) - bg_subtract
#       # Plot the points. _______________________________________________________
#       points(x - point_shift/x_denom, y, col=kThrpReporterColors[thrp_site], lwd=1, pch=1, xpd=NA)
#       # Plot the horizontal lne associated with the average. ___________________
#       segments(x0=x[1] - 0.25/x_denom - point_shift/x_denom, y0=mean(y),
#                x1=x[1] + line_r, xpd=NA)
#       # Optionally add the lin-41 UTR repression.
#       if (all_lin41_UTR) {
#         y_old <- y
#         mean_old <- mean(y)
#         inds <- which((log2fc_df[, "Seed"] == "8mer-w6_lin41_UTR") &
#                       (log2fc_df[, "ThreePrime"] == thrp_site) &
#                       (log2fc_df[, "Bulge"] == bulge_use))
#         log2fc_df_use <- log2fc_df[inds, 5:ncol(log2fc_df)]
#         x <- rep(x_start, length(unlist(log2fc_df_use)))
#         y <- unlist(log2fc_df_use) - bg_subtract

#         y_difs <- y - y_old
#         se_difference <- sd(y_difs)/sqrt(length(y_difs))
#         t_statistic <- mean(y_difs)/se_difference
#         p_val <- 2*pt(t_statistic, df=length(y) - 1)
#         # This versin of it does it treating the reps as unpaired, which is more
#         # similar to how the Tukey tests were done.
#         p_val <- t.test(y, y_old, alternative="two.sided", var.equal=FALSE)$p.value
#         p_value_row <- data.frame(site=sprintf("%s_%s_%s", thrp_site, bulge_use, config_str), p_vals=p_val)
#         output_p_values <- rbind(output_p_values, p_value_row)
#         points(x + point_shift/x_denom, y, col=kThrpReporterColors[thrp_site], lwd=1, pch=1, xpd=NA)
#         segments(x0=x[1], x1=x[1] + 0.25/x_denom + point_shift/x_denom, y0=mean(y), xpd=NA)
#         segments(x0=x[1], x1=x[1], y0=mean_old, y1=mean(y), xpd=NA)
#       }

#       # Format the bulge string. _______________________________________________
#       if (bulge_use == "") {
#         bulge_use <- "0"
#       } else {
#         bulge_use <- gsub("T", "U", bulge_use)
#         len_bulge <- nchar(bulge_use)
#         bulge_nuc <- strsplit(bulge_use, split="")[[1]]
#         bulge_use <- bquote(
#           .("+")*.(len_bulge)[.(bulge_nuc)]
#         )
#       }
#       # Add the bulge label. ___________________________________________________
#       text(x_start, y_bulge, label=bulge_use, adj=0.5, xpd=NA)
#       # Add the double- or single-, or no asterisk for significance. ___________
#       site_label <- sprintf("%s-%s", thrp_site, bulge)
#       # Update the x-axis position. ____________________________________________
#       x_start <- x_start + 1/x_denom
#     }
#   }
#   ############################## Add lin-41 sites ##############################
#   if (dual_site) {
#     ################### Plot the lin-41 sites in all contexts. #################
#     inds <- which(log2fc_df[, "Seed"] == "lin-41")
#     log2fc_df_use <- log2fc_df[inds, 5:ncol(log2fc_df)]
#     x <- rep(x_start, length(unlist(log2fc_df_use)))
#     y <- unlist(log2fc_df_use) - bg_subtract
#     # Plot the points. _________________________________________________________
#     points(x - point_shift/x_denom, y, col="goldenrod", lwd=1, pch=1, xpd=NA)
#     # Plot the horizontal lne associated with the average. ___________________
#     segments(x0=x[1] - 0.25/x_denom - point_shift/x_denom, y0=mean(y),
#              x1=x[1] + line_r, xpd=NA)
#     text(x_start, y_bulge, label="lin-41 sites", srt=90, adj=0, xpd=NA)

#     if (all_lin41_UTR) {
#       y_old <- y
#       mean_old <- mean(y)
#       inds <- which(log2fc_df[, "Seed"] == "lin-41_lin41_UTR")
#       log2fc_df_use <- log2fc_df[inds, 5:ncol(log2fc_df)]
#       x <- rep(x_start, length(unlist(log2fc_df_use)))
#       y <- unlist(log2fc_df_use) - bg_subtract
#       y_difs <- y - y_old
#       se_difference <- sd(y_difs)/sqrt(length(y_difs))
#       t_statistic <- mean(y_difs)/se_difference
#       p_val <- 2*pt(t_statistic, df=length(y) - 1)
#       # This versin of it does it treating the reps as unpaired, which is more
#       # similar to how the Tukey tests were done.
#       p_val <- t.test(y, y_old, alternative="two.sided", var.equal=FALSE)$p.value
#       p_value_row <- data.frame(site="lin-41", p_vals=p_val)
#       output_p_values <- rbind(output_p_values, p_value_row)
#       points(x + point_shift/x_denom, y, col="goldenrod", lwd=1, pch=1, xpd=NA)
#       segments(x0=x[1], x1=x[1] + 0.25/x_denom + point_shift/x_denom, y0=mean(y), xpd=NA)
#       segments(x0=x[1], x1=x[1], y0=mean_old, y1=mean(y), xpd=NA)
#     }
#   }
#   ######################### Add the site category labels. ######################
#   text(0, y_seed, label="Seed match:", adj=1, xpd=NA)
#   text(0, y_thrp, label="3' pairing:", adj=1, xpd=NA)
#   text(0, y_bulge, label="Offset (nt):", adj=1, xpd=NA)
#   # # ############################## Add axes ##################################
#   AddLinearAxis(2, tick.space=0.5, label.space=1,
#                 label="Fold change (log2)")
#   if (dual_site)  label <- "Dual sites"
#   else            label <- "Single sites"
#   y_label <- 1.15*y_const
#   if (all_lin41_UTR)  x_use <- -1.75
#   else                x_use <- -2.25
#   text(x=x_use/x_denom, y=y_label, label=label, xpd=NA, adj=c(0, 0))
#   ############################# Finish the plot ################################
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
#   if (all_lin41_UTR) {
#     output_p_values_global <<- output_p_values
#   }
# }




# PlotReporterEfficacyAgainstProgKds <- function(
#   mirna="let-7a-21nt", experiment="equil_c2_nb",
#   experiment_rep="twist_reporter_assay_3p_2_tp", n_constant=3,
#   sitelist="progthrp_suppcomp", nbomitc=FALSE, corrected_kds=TRUE,
#   collapsemm=FALSE, equilibrium_nb=FALSE, exp_type="parallel",
#   height=4.5, width=4.5, xpos=20, ypos=20, pdf.plot=FALSE
# ) {
#   kds <- SubfunctionCall(ApplyKdCorrection, prog_n_constant=n_constant,
#                          rand_n_constant=n_constant, prog_sitelist=sitelist,
#                          new2=TRUE)
#   log2fc_df_1 <- SubfunctionCall(
#     GetThreePrimeReporterFoldChanges, experiment=experiment_rep, mirna="let-7a",
#     rep=1, aggregate=TRUE
#   )
#   log2fc_df_2 <- SubfunctionCall(
#     GetThreePrimeReporterFoldChanges, experiment=experiment_rep, mirna="let-7a",
#     rep=2, aggregate=TRUE
#   )
#   log2fc_df_3 <- SubfunctionCall(
#     GetThreePrimeReporterFoldChanges, experiment=experiment_rep,
#     mirna="let-7a-21nt", rep=1, aggregate=TRUE
#   )
#   log2fc_df_4 <- SubfunctionCall(
#     GetThreePrimeReporterFoldChanges, experiment=experiment_rep,
#     mirna="let-7a-21nt", rep=2, aggregate=TRUE
#   )
#   # Subset the information from the first matrix (could have used any, they
#   # are redundant).
#   log2fc_df <- log2fc_df_1[, 1:4]
#   # Add each of the replicates to the dataframe.
#   log2fc_df$log2fc_1 <- log2fc_df_1$log_fc
#   log2fc_df$log2fc_2 <- log2fc_df_2$log_fc
#   log2fc_df$log2fc_3 <- log2fc_df_3$log_fc
#   log2fc_df$log2fc_4 <- log2fc_df_4$log_fc
#   fold_change_values <- c()
#   colors <- c()
#   pchs <- c()

#   # Iterate over each of the site types to get the average log2fc to plot
#   # the Kd values.
#   seed_sites <- c("None", "8mer", "7mer-m8", "7mer-A1", "6mer", "8mer-w6")
#   label_index = 1
#   sites_track <- c()
#   bulges_track <- c()
#   for (seed_site in seed_sites) {
#     inds <- which((log2fc_df$Seed == seed_site) &
#                   (log2fc_df$ThreePrime == "") &
#                   log2fc_df$Location != "dual")
#     log2fc_df_use <- log2fc_df[inds, 5:8]
#     if (seed_site == "None") {
#       bg_subtract <- mean(unlist(log2fc_df_use))
#     }
#     # Correct the selected fold change values for the background.
#     fc_use <- mean(unlist(log2fc_df_use)) - bg_subtract
#     fold_change_values <- c(fold_change_values, fc_use)
#     if (seed_site == "8mer-w6") seed_site_use <- "8mer-mmG6"
#     else                        seed_site_use <- seed_site
#     names(fold_change_values)[label_index] <- sprintf("%s_Kd", seed_site_use)
#     colors <- c(colors, kSiteColors[seed_site])
#     pchs <- c(pchs, 19)
#     label_index <- label_index + 1
#     sites_track <- c(sites_track, seed_site)
#     bulges_track <- c(bulges_track, "-")
#   }

#   # Define the relevant parameters for the 3-prime sites
#   thrp_sites <- c("4mer-m13.16", "9mer-m11.19", "9mer-m13.21")
#   bulges <- c("None", "A", "T", "AAAA", "TTTT")
#   thrp_pchs <- c(4, 1, 0, 5, 2)
#   names(thrp_pchs) <- bulges
#   # Predefine list of labels that will need to be modified for the Kd correction
#   # for the 9mer-m11.19|AAAA| and 9mer-m13.21|TTTT| sites.
#   labels_shift <- c()
#   for (thrp_site in thrp_sites) {
#     for (bulge in bulges) {
#       if (bulge == "None")  bulge_use <- ""
#       else                  bulge_use <- bulge
#       # First calculate which 
#       # Determine the row indeces of the sites.
#       inds <- which((log2fc_df$Seed == "8mer-w6") &
#                     (log2fc_df$ThreePrime == thrp_site) &
#                     (log2fc_df$Bulge == bulge_use) &
#                     (log2fc_df$Location != "dual"))
#       # Subset the log2fc matrix to just have the fold-change values for the
#       # appropriate rows.
#       log2fc_df_use <- log2fc_df[inds, c("log2fc_1", "log2fc_2",
#                                          "log2fc_3", "log2fc_4")]
#       # Calculate the average log2 fold-change, using the same bg_subtract value
#       # as determined when iterating over the seed sequences.
#       fc_use <- mean(unlist(log2fc_df_use)) - bg_subtract
#       # Append the fold-change value to the growing list of fold-change values.
#       fold_change_values <- c(fold_change_values, fc_use)
#       mir_position <- gsub("^.*mer-m(.*)\\..*$", replacement="\\1", thrp_site)
#       tar_position <- as.integer(mir_position) + nchar(bulge_use)
#       if (bulge_use == "") {
#         names(fold_change_values)[label_index] <- sprintf("%s|%s|Comp_Kd", thrp_site, tar_position)
#       } else if ((thrp_site %in% c("9mer-m11.19", "9mer-m13.21")) &
#           (bulge %in% c("AAAA", "TTTT"))) {
#         names(fold_change_values)[label_index] <- sprintf("%s|%s|WWWW|Comp_Kd", thrp_site, tar_position)
#         labels_shift <- c(labels_shift, label_index)
#         names(labels_shift)[length(labels_shift)] <- sprintf("%s|%s", thrp_site, bulge_use)
#       } else {
#         names(fold_change_values)[label_index] <- sprintf("%s|%s|%s|Comp_Kd", thrp_site, tar_position, bulge_use)        
#       }
#       colors <- c(colors, kThrpReporterColors[thrp_site])
#       pchs <- c(pchs, thrp_pchs[bulge])
#       label_index <- label_index + 1
#       sites_track <- c(sites_track, thrp_site)
#       bulges_track <- c(bulges_track, bulge)
#     }
#   }
#   kds <- rbind(kds, None_Kd=c(1, 1, 1, 1, 1))
#   kds_use <- kds[names(fold_change_values), ]
#   kds_use_pre <- kds_use
#   ############################################################################## 
#   # Portion of figure related to rescaling the three-prime sites based on their
#   # loop sequences.
#   sXc_rescale <- SitesXCounts("let-7a-21nt", experiment="equil_c2_nb", n_constant=3,
#                       sitelist="progthrp_suppcomp", new2=TRUE, include_zeros=TRUE)
#   print(dim(sXc_rescale))
#   print(colSums(sXc_rescale))
#   sXc_rescale <- as.data.frame(t(t(sXc_rescale)/colSums(sXc_rescale)))
#   print(dim(sXc_rescale))
#   sXc_rescale_global <<- sXc_rescale
#   # Average the reps (this is only appropriate for the second let-7a-21nt
#   # experiment, since it has two 12.65% reps).
#   sXc_rescale[, c("12.65")] <- rowMeans(sXc_rescale[, c("12.65", "12.65_2")], na.rm=TRUE)
#   # Subsample the counts matrix to get the relevant sequences for the calculation
#   # for the 9mer-m11.19.
#   sXc_rescale <- sXc_rescale[, c("0.4", "1.265", "4", "12.65", "40")]
#   s_11_19_A <- sXc_rescale[grep("9mer-m11\\.19\\|15\\|AAAA", perl=TRUE, rownames(sXc_rescale), value=TRUE), ]
#   s_11_19_U <- sXc_rescale[grep("9mer-m11\\.19\\|15\\|TTTT", perl=TRUE, rownames(sXc_rescale), value=TRUE), ]
#   s_11_19_AU_ratio <- sum(s_11_19_A)/sum(s_11_19_U)

#   # Should really perform the grep every time.
#   AAAA_I <- as.integer(system("cut -c -25 /lab/solexa_bartel/mcgeary/AgoRBNS/let-7a-21nt/equil_c2_nb/reads/I.txt | grep AAAA | wc -l", intern=TRUE))
#   UUUU_I <- as.integer(system("cut -c -25 /lab/solexa_bartel/mcgeary/AgoRBNS/let-7a-21nt/equil_c2_nb/reads/I.txt | grep TTTT | wc -l", intern=TRUE))
#   print(AAAA_I)
#   print(UUUU_I)
#   AAAA_I <- 3870743
#   UUUU_I <- 3256850
#   print(AAAA_I)
#   print(UUUU_I)
#   s_11_19_AU_ratio <- s_11_19_AU_ratio/(AAAA_I/UUUU_I)
#   for (label_shift_name in names(labels_shift)) {
#     print(label_shift_name)
#     if (label_shift_name == "9mer-m11.19|AAAA") {
#       kds_use[labels_shift[label_shift_name], ] <- kds_use[labels_shift[label_shift_name], ]/sqrt(s_11_19_AU_ratio)
#     } else if (label_shift_name == "9mer-m11.19|TTTT") {
#       kds_use[labels_shift[label_shift_name], ] <- kds_use[labels_shift[label_shift_name], ]*sqrt(s_11_19_AU_ratio)
#     }
#   }

#   SubfunctionCall(FigureSaveFile2)
#   xmin <- 3e-4
#   xmax <- 3
#   ymin <- -1.5
#   ymax <- 0.25
#   par(mar=c(3, 3, 1, 1))
#   BlankPlot(log="x", inv="x")
#   x <- kds_use$Mean
#   print(x)
#   print(labels_shift)
#   print(cbind(sites_track, bulges_track, x, kds_use_pre$Mean))
#   y <- fold_change_values

#   ############ Fit trend line to biochemical model #############################
#   RepressionModel <- function(kds, a, b) {
#     1/(1 + b*(a/(a + kds)))
#   }
#   LossFunction <- function(pars) {
#     a <- exp(pars[1])
#     b <- exp(pars[2])
#     rep_predicted <- log(RepressionModel(x, a, b), 2)
#     return(sum((y - rep_predicted)^2))
#   }
#   # Optimize the model for `a` and `b`.
#   opt <- optim(c(0, 0), LossFunction, method="L-BFGS-B", lower=c(-10, -10), 
#                upper=c(10, 0))
#   # Extract the `a` and `b` parameters.
#   a <- exp(opt$par[1])
#   b <- exp(opt$par[2])
#   a_global <<- a
#   b_global <<- b
#   # Get equally-spaced lines between the min and max of the x-axis and apply the
#   # Model function with optimized parameters to get the trend line.
#   x_lines <- 10^seq(log10(xmin), log10(xmax), length.out=100)
#   out_model <- log(RepressionModel(x_lines, a, b), 2)
#   lines(x_lines, out_model, lwd=2, col="gray")
#   ###################### Plot the points #######################################
#   Points(x, y, col=colors, pch=pchs, pt.lwd=1)
#   ############################ Add axes. #######################################
#   AddLogAxis(1, label="Relative Kd")
#   AddLinearAxis(2, 0.25, 0.5, label="Fold change (log2)")
#   ################ Add r-squred and corresponding text labels. #################
#   x_text <- GetPlotFractionalCoords(0.35, 0.15, log="x", inv="x")[1]
#   # xy <- GetPlotFractionalCoords(0.05, 0.15, log="x", inv="x")
#   # AddCorrelationToPlot(x=log(x), y=y, xpos=xy[1], ypos=xy[2], rsquared=TRUE)
#   # text(x=x_text, y=xy[2], label="Linear; all points", adj=0)

#   # xy <- GetPlotFractionalCoords(0.05, 0.10, log="x", inv="x")
#   # AddCorrelationToPlot(x=log(x[-1]), y=y[-1], xpos=xy[1], ypos=xy[2], rsquared=TRUE)
#   # text(x=x_text, y=xy[2], label="Linear; Excluding None", adj=0)

#   xy <- GetPlotFractionalCoords(0.8, 0.05, log="x", inv="x")
#   AddCorrelationToPlot(x=y, y=log(RepressionModel(x, a, b), 2), xpos=xy[1], ypos=xy[2], rsquared=TRUE)
#   # text(x=x_text, y=xy[2], label="Model fit", adj=0)

#   xy <- GetPlotFractionalCoords(0.05, 1.05, log="x", inv="x")
#   Legend(xy, legend=seed_sites, col=kSiteColors[seed_sites], y.intersp=0.75,
#          ncol=3)

#   # Add the 3-prime legend.
#   xy <- GetPlotFractionalCoords(0.05, 0.25, log="x", inv="x")
#   text(xy[1], xy[2], label="Seed m.m. with", adj=c(0, 0))
#   xy <- GetPlotFractionalCoords(0.05, 0.20, log="x", inv="x")
#   text(xy[1], xy[2], label="3' pairing:", adj=0)
#   xy <- GetPlotFractionalCoords(0.0, 0.20, log="x", inv="x")
#   Legend(xy, legend=gsub("^.*mer-m(.*)\\.(.*)$", replacement="\\1-\\2",
#                          thrp_sites, perl=TRUE),
#          col=kThrpReporterColors, y.intersp=0.75)
#   xy <- GetPlotFractionalCoords(0.35, 0.20, log="x", inv="x")
#   text(xy[1], xy[2], label="Offset:", adj=0)
#   xy <- GetPlotFractionalCoords(0.30, 0.20, log="x", inv="x")
#   # Reformatting of the bulge strings.
#   bulges_legend <- gsub("T", "U", bulges)
#   bulges_legend <- gsub("None", "-", bulges_legend)
#   bulges_legend <- sapply(bulges_legend, function(bulge) {
#     if (bulge == "-") {
#       return("  0")
#     } else {
#       nuc <- strsplit(bulge, split="")[[1]][1]
#       len_off <- nchar(bulge)
#       bulge_use <- bquote(.(nuc)[.(len_off)])
#       bulge_use <- bquote(.("+")*.(len_off)[.(nuc)])
#       return(as.expression(bulge_use))
#     }
#   })
#   legend(xy[1], xy[2], legend=bulges_legend, pch=thrp_pchs, bty="n", pt.lwd=1,
#          y.intersp=0.75, ncol=2)
#   # Finish the plot.
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }

# }


# # Functions found in FIGURE 2___________________________________________________

# ################################################################################
# # FIGURE S3
# ################################################################################

# # Functions found in FIGURE S2__________________________________________________

# ################################################################################
# # FIGURE 4
# ################################################################################

# # 4A-C,_________________________________________________________________________
# PlotPairingCoefficients <- function(
#   mirna, experiment, n_constant=3, sitelist="progthrp_suppcomp",
#   offset_lim=c(-4, 16), len_lim=c(4, 11), pos_lim=c(9, 23), corrected_kds=TRUE,
#   makeglobalmodel=TRUE, cutoff=TRUE, supp_base=FALSE, fixed_offset=0,
#   log_plus_one=FALSE, F_method=FALSE, exponential=FALSE, intercept=FALSE, LowerCI=FALSE,
#   UpperCI=FALSE, offset_base=FALSE, site_base=NULL, mm=FALSE, delG=FALSE,
#   deldelG=FALSE, deldelGleft=FALSE, key=FALSE, xlabels=TRUE, mirna_label=TRUE,
#   decomp_pairing=FALSE, fit_max_len=5, keep_data=TRUE,
#   height=2.8, width=2.5, xpos=20, ypos=20, pdf.plot=FALSE
# ) {
#   print("PlotPairingCoefficients")
#   # Thermodynamic constants for conversion to delta G.
#   R <- 1.987e-3 # in kcal K-1 mol-1
#   T <- 310.15 # in K
#   # Global model is a parameter to speed up fitting the model.
#   if (makeglobalmodel) {
#     if (mm) {
#       model <<- SubfunctionCall(FitPairingOffsetAndMismatchModel)  
#     } else {
#       model <<- SubfunctionCall(FitPairingAndOffsetModelWithError)  
#     }
#   }
#   if (LowerCI) {
#     R_mat <- t(model$pairing$LowerCI)
#   } else if (UpperCI) {
#     R_mat <- t(model$pairing$UpperCI)
#   } else if (decomp_pairing) {
#     R_mat_orig <- t(model$pairing$MLE)
#     print(R_mat_orig)
#     R_mat <- t(
#       DecomposePairingCoefficients(model$pairing$MLE, fit_max_len=fit_max_len,
#                                    keep_data=keep_data)
#     )
#   } else {
#     R_mat <- t(model$pairing$MLE)
#   }
#   mar1 <- 2.2
#   mar2 <- 0.5
#   mar3 <- 1.8
#   mar4 <- 2.5
#   R_mat <<- R_mat
#   # This is the conversion factor in order to have the correct height of the
#   # plot for the specified width, that makes each box of equal height and width.
#   if (deldelGleft) {
#     R_mat1 <- R_mat[-1, ]
#     R_mat2 <- R_mat[-nrow(R_mat), ]
#     R_mat <- t((R_mat2 - R_mat1)[, -1])

#   } else if (deldelG) {
#     R_mat1 <- R_mat[, -ncol(R_mat)]
#     R_mat2 <- R_mat[, -1]
#     R_mat <- (R_mat2 - R_mat1)[-nrow(R_mat), ]
#     R_mat_del <<- R_mat
#   }
#   mir_length <- nchar(kMirnaSeqs[mirna])
#   SubfunctionCall(FigureSaveFile2)
#   par(mar=c(mar1, mar2, mar3, mar4))
#   xmin <- 0
#   xmax <- nrow(R_mat)
#   ymin <- 0
#   ymax <- ncol(R_mat) 
#   BlankPlot()
#   ymin <- 0

#   xlefts <- rep(seq(0, ncol(R_mat) - 1), nrow(R_mat))
#   xright <- xlefts + 1
#   ybottom <- rep(rev(seq(ymax - 1, ymin)), each=ncol(R_mat))
#   ytop <- ybottom + 1
#   # Make every other label on the y-axis blank.
#   y_axis_pos <- seq(ymin, ymax - 1) + 0.5
#   y_axis_labs <- rownames(R_mat)
#   for (i in 1:length(y_axis_labs)) {
#     if (i %% 2 == 0) {
#       y_axis_labs[i] <- ""
#     }
#   }
#   # Make the x-axis
#   AddLinearAxis(1, 1, 1, label="3'-paired nt",
#                 label_pos_ticks=TRUE,
#                 alt_lab=colnames(R_mat),
#                 alt_lab_y_dist=0.01,
#                 alt_lab_pos=xmin:(xmax - 1) + 0.5,
#                 alt_tick_pos=TRUE)
#   # Make the y-axis
#   AddLinearAxis(4, 1, 2, label="5'-paired nt",
#                 alt_lab=y_axis_labs,
#                 alt_lab_pos=y_axis_pos,
#                 alt_tick_pos=TRUE,
#                 line=0.8)

#   # Add the label for the seed if not plotting relative to seed kds.
#     message(sprintf("This is the min Kd fold change for %s: %s", mirna, 10^min(R_mat, na.rm=TRUE)))
#     if (deldelGleft || deldelG) {
#       col.inds <- floor((R_mat)/log10(20)*100)
#       x_lab_pos <- -0.75
#     } else if (delG) {
#       # R_mat <<- R_mat
#       R_mat_adjusted <- R_mat/log10(exp(1))*R*T
#       R_mat_adjusted <<- R_mat_adjusted
#       col.inds <- floor(R_mat_adjusted/4*99 + 1)
#       col.inds <<- col.inds
#       x_lab_pos <- -0.75
#     } else {
#       col.inds <- floor((R_mat)/log10(700)*99 + 1)
#       x_lab_pos <- -0.75
#     }
#   ################### Assign colors to each square in the matrix ###############
#   # Make distribution of 100 colors and a 101th gray color for NA indeces.
#   start_col <- 0.5
#   r_col <- -0.75
#   color.dist <- c(rev(cubeHelix(100, start=start_col, r=r_col, hue=0.8)),
#                   "gray90", "white")
#   color.dist[100] <- kMaxValueColor
#   # Make the color index scale.
#   col.inds <- t(col.inds)
#   # Assign values of "NA" to gray, but assign impossible positions to white.
#   pos_5p <- as.integer(rep(rownames(R_mat), each=ncol(R_mat)))
#   pos_3p <- as.integer(rep(colnames(R_mat), nrow(R_mat)))
#   if (decomp_pairing) {
#     impossible_cols <- which((pos_3p - pos_5p + 1 < 1 )| (pos_3p - pos_5p + 1 > 11))
#   } else {
#     impossible_cols <- which((pos_3p - pos_5p + 1 < 4 )| (pos_3p - pos_5p + 1 > 11))
#   }
#   # impossible_cols <- c()
#   col.inds[which(col.inds > 100)] <- 100
#   col.inds[which(col.inds < 1)] <- 1
#   col.inds[which(is.na(col.inds))] <- 101
#   col.inds[impossible_cols] <- 102
#   #################### Add the rectangles ############################
#   rect(xlefts, ybottom, xright, ytop, col=color.dist[col.inds], lwd=0.4, 
#        xpd=NA, border="white")
#   ####################### Add the key ##########################################
#   if (key) {
#     mir_scale <- (mir_length - 11)/(21 - 11)
#     num_boxes <- 100
#     y_div <- (ymax - ymin)/num_boxes
#     kpl <- xmax + 5*mir_scale # Left-hand position of the key
#     kw <- mir_scale                # Width of the key
#     rect(xleft=kpl,
#          ybottom=seq(100)*y_div - y_div,
#          xright=kpl + kw,
#          ytop=seq(100)*y_div,
#          col=color.dist[1:100], xpd=NA, border=NA)
#     # Generate the axis for the legend and the label
#     bottom <- ymin + y_div/2
#     top <- ymax - y_div/2
#     if (delG) {
#       labels <- c(0, 1, 2, 3, 4)    
#       pos_labels <- labels
#     } else {
#       labels <- c(1, 3, 10, 30, 100)
#       pos_labels <- log(labels)
#     }
#     centered_labels <- pos_labels - pos_labels[1]
#     norm_labels <- centered_labels/(centered_labels[length(centered_labels)])
#     height_span <- top - bottom
#     pos_labels <- norm_labels*height_span

#     axis(4, at=pos_labels, labels=labels, lwd=par()$lwd, pos=kpl + kw, xpd=NA)
#     text(x=kpl + kw + 3, y=ymax/2,
#          labels=bquote(italic(K)[D]*.(" fold change")),
#          srt=270, xpd=NA)    
#   }
#   mtext(text="Pairing coefficients", side=3, at=xmax, adj=1, cex=par("cex"))
#   ##################### Add the miRNA label #####################
#   if (mirna_label) {
#     if (mirna == "let-7a-21nt") {
#       mirna_txt <- "let-7a"
#     } else if (mirna == "let-7a_minus1") {
#       mirna_txt <- "let-7a(-1)"
#     } else if (mirna == "let-7a_plus1") {
#       mirna_txt <- "let-7a(+1)"
#     } else if (mirna == "let-7a_miR-155") {
#       mirna_txt <- "let-7a-miR-155"
#     } else if (mirna == "miR-155_let-7a") {
#       mirna_txt <- "miR-155-let-7a"
#     } else if (mirna == "miR-7-23nt") {
#       mirna_txt <- "miR-7"
#     } else {
#       mirna_txt <- mirna
#     }
#     x <- GetPlotFractionalCoords(0.125, 0.5)[1]
#     mtext(text=mirna_txt, side=3, line=0.8, at=x, adj=0, cex=par("cex"))
#   }
#   ############### Finish the plot #############################
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
#   return(R_mat)
# }

# # 4A-C__________________________________________________________________________
# PlotOffsetCoefficients <- function(
#   mirna, experiment, n_constant=3, offset_lim=c(-4, 16), makeglobalmodel=TRUE,
#   sitelist="progthrp_suppcomp", corrected_kds=TRUE,  F_method=FALSE,
#   exponential=FALSE, intercept=FALSE, cutoff=TRUE, sumseed=FALSE,
#   supp_base=FALSE, offset_base=FALSE, site_base=NULL, mm=FALSE,
#   len_lim=c(4, 11), pos_lim=c(9, 23), weighted_ave=FALSE, log_plus_one=FALSE,
#   xpos=20, ypos=20, height=2.8, width=3, pdf.plot=NULL, weights_new=FALSE
# ) {
#   if (makeglobalmodel) {
#     if (mm) {
#       model <<- SubfunctionCall(FitPairingOffsetAndMismatchModel)
#     } else {
#       model <<- SubfunctionCall(FitPairingAndOffsetModelWithError)
#     }
#   }
#   print("got here")
#   offsets <- model$offsets
#   offsets <<- offsets
#   SubfunctionCall(FigureSaveFile2)
#   xmin <- -4
#   xmax <- 16
#   if (log_plus_one) {
#     if (experiment %in% c("equilibrium", "equilibrium2_nb")) {
#       ymin <- -4
#     } else {
#       ymin <- -2.5
#     }
#     ymax <- 0.5
#   } else {
#     ymin <- 0
#     ymax <- 1.05
#   }
#   mar1 <- 2.2
#   mar2 <- 2.8
#   mar3 <- 1.8
#   mar4 <- 0.5
#   par(mar=c(mar1, mar2, mar3, mar4))  
#   BlankPlot(inv="x")
#   ymax <- 1
#   start_col <- 0.5
#   r_col <- -0.75
#   if (log_plus_one) {
#     col.inds <- floor((offsets[, 1] + ymin)/(-ymin)*99) + 1  
#   } else {
#     col.inds <- floor(offsets[, 1]*99) + 1  
#   }
#   # Make distribution of 100 colors and a 101th gray color for NA indeces.
#   # color.dist <- c(rev(cubeHelix(100, start=start_col, r=r_col, hue=0.8)), "gray")
#   color.dist <- rep(cubeHelix(100, start=0.25*pi, r=r_col, hue=0.8)[50], 101)
#   # Make the color index scale.
#   col.inds <- sapply(col.inds, function(col.ind) {
#     min(max(1, col.ind), 100)
#   })
#   if (weighted_ave) {
#     offsets_use <- offsets[which(as.integer(rownames(offsets)) <= 16
#                                  & as.integer(rownames(offsets)) >= -4)]
#     weighted_average <- sum(as.integer(names(offsets_use))*offsets_use)/sum(offsets_use)
#     print(weighted_average)
#     segments(x0=weighted_average, y0=0, y1=1, lty=3)

#   }
#   max_offset <<- as.integer(rownames(offsets)[which(offsets[, 1] == max(offsets[, 1]))])
#   lines(as.integer(rownames(offsets)), offsets[, 1], col="gray", xpd=NA)  
#   lines(as.integer(rownames(offsets)), offsets[, 2], col="gray", lty=2, xpd=NA)  
#   lines(as.integer(rownames(offsets)), offsets[, 3], col="gray", lty=2, xpd=NA)  
#   points(as.integer(rownames(offsets)), offsets[, 1], col=color.dist[col.inds],
#          pch=19, xpd=NA)  
#   # Add the axes.
#   AddLinearAxis(1, 1, 4, label="Offset (nt)")
#   AddLinearAxis(2, 0.1, 0.2, label="Fraction maximum")
#   # Add the label indicating the type of plot that it is.
#   mtext(side=3, text="Offset coefficients", at=xmin, adj=1, cex=par()$cex)
#   # Finish the plot.
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }

# # 4A____________________________________________________________________________
# PlotPairingMatrixKey <- function(
#   model=TRUE, delG=FALSE, alt_top=FALSE, short_margins=FALSE, horizontal=FALSE,
#   height=2.8, max_value=700, width=1, xpos=20, ypos=20, pdf.plot=NULL
# ) {
#   if (horizontal) {
#     height <- 0.95
#     width <- 3.5
#   }
#   SubfunctionCall(FigureSaveFile2)
#   if (short_margins) {
#     mar1 <- 0.5
#     mar3 <- 0.5
#   } else {
#     mar1 <- 2.2
#     mar3 <- 1.8
#   }
#   mar2 <- 0.5
#   mar4 <- 2.5
#   if (horizontal) {
#     par(mar=c(0.5, 1, 2.5, 1))
#   } else {
#     par(mar=c(mar1, mar2, mar3, mar4))
#   }
#   xmin <- 0
#   xmax <- 1
#   ymin <- 0
#   ymax <- 1
#   BlankPlot()
#   if (model) {
#     if (delG) {
#       labels <- 0:-4
#       pos_labels <- labels
#       axis_label <- expression(Delta*Delta*italic(G)~"(kcal/mol)")
#       axis_space <- 1.2
#     } else {
#       if (alt_top) {
#         labels <- c(1, 3, 10, 30, 100, 300, 700, 100*round(max_value/100))
#       } else {
#         labels <- c(1, 3, 10, 30, 100, 300, 700)
#       }
#       pos_labels <- log(labels)
#       axis_label <- bquote(italic(K)[D]*.(" fold change"))
#       axis_space <- 1.6
#     }
#     start_col <- 0.5
#     r_col <- -0.75
#   } else {
#     if (alt_top) {
#       labels <- c(1, 3, 10, 30, 100, 300, 700, 100*round(max_value/100))
#     } else {
#       labels <- c(1, 3, 10, 30, 100, 300, 700)
#     }
#     pos_labels <- log(labels)
#     axis_label <- bquote(italic(K)[D]*.(" fold change"))
#     axis_space <- 1.6
#     start_col <- 0
#     r_col <- 0.4
#   }
#   if (alt_top) {
#     num_boxes <- ceiling(log(labels[length(labels)])/log(700)*100)
#   } else {
#     num_boxes <- 100
#   }
#   color.dist <- c(rev(cubeHelix(100, start=start_col, r=r_col, hue=0.8)),
#                   rep(kMaxValueColor, num_boxes - 100))
#   y_div <- (ymax - ymin)/num_boxes
#   kpl <- xmin # Left-hand position of the key
#   kw <- 0.6
#   if (horizontal) {
#     rect(ybottom=kpl, xleft=seq(num_boxes)*y_div - y_div, ytop=kpl + kw,
#          xright=seq(num_boxes)*y_div - y_div + y_div,
#          col=color.dist, xpd=NA, border=NA)
#   } else {
#     rect(xleft=kpl, ybottom=seq(num_boxes)*y_div - y_div, xright=kpl + kw,
#          ytop=seq(num_boxes)*y_div - y_div + y_div,
#          col=color.dist, xpd=NA, border=NA)
#   }               # Width of the key
#   # Generate the axis for the legend and the label
#   bottom <- ymin + y_div/2
#   top <- ymax - y_div/2
#   centered_labels <- pos_labels - pos_labels[1]
#   norm_labels <- centered_labels/(centered_labels[length(centered_labels)])
#   height_span <- top - bottom
#   pos_labels <- norm_labels*height_span + y_div/2
#   if (labels[1] == 1) {
#     labels <- as.character(labels)
#     labels[1] <- expression(""<=1)
#   }
#   if (horizontal) {
#     axis(3, at=pos_labels, labels=labels, lwd=par()$lwd, pos=kpl + kw, xpd=NA)
#     text(y=kpl + kw + axis_space, x=ymax/2,
#          labels=axis_label,
#          xpd=NA)
#   } else {
#     axis(4, at=pos_labels, labels=labels, lwd=par()$lwd, pos=kpl + kw, xpd=NA)
#     text(x=kpl + kw + axis_space, y=ymax/2,
#          labels=axis_label,
#          srt=270, xpd=NA)
#   }
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }  

# # 4D ___________________________________________________________________________
# PlotDeltaGMatrix <- function(
#   mirna, just_complement=TRUE, wobble=TRUE, deldelG=FALSE, deldelGleft=FALSE, key=FALSE, xlabels=TRUE,
#   mirna_label=TRUE, height=2.8, width=2.5, xpos=20, ypos=20, pdf.plot=FALSE
# ) {
#   # Load the data matrix.
#   R_mat <- t(SubfunctionCall(GetDeltaGPredictedMatrix))
#   DeltaG_mat_global <<- R_mat
#   mar1 <- 2.2
#   mar2 <- 0.5
#   mar3 <- 1.8
#   if(key) {
#     width <- 4
#     mar4 <- 9
#   } else {
#     mar4 <- 2.5
#   }
#   # This is the conversion factor in order to have the correct height of the
#   # plot for the specified width, that makes each box of equal height and width.
#   mir_length <- nchar(kMirnaSeqs[mirna])
#   SubfunctionCall(FigureSaveFile2)
#   par(mar=c(mar1, mar2, mar3, mar4))
#   xmin <- 0
#   xmax <- nrow(R_mat)
#   ymin <- 0
#   ymax <- ncol(R_mat) 
#   BlankPlot()
#   xlefts <- rep(seq(0, ncol(R_mat) - 1), nrow(R_mat))
#   xright <- xlefts + 1
#   ybottom <- rep(rev(seq(ymax - 1, ymin)), each=ncol(R_mat))
#   ytop <- ybottom + 1

#   # Make the x-axis
#   # Make the x-axis
#   y_axis_pos <- seq(ymin, ymax - 1) + 0.5
#   y_axis_labs <- rownames(R_mat)
#   for (i in 1:length(y_axis_labs)) {
#     if (i %% 2 == 0) {
#       y_axis_labs[i] <- ""
#     }
#   }
#   AddLinearAxis(1, 1, 1, label="3'-paired nt",
#                 label_pos_ticks=TRUE,
#                 alt_lab=colnames(R_mat),
#                 alt_lab_y_dist=0.01,
#                 alt_lab_pos=xmin:(xmax - 1) + 0.5,
#                 alt_tick_pos=TRUE)
#   AddLinearAxis(4, 1, 2, label="5'-paired nt",
#                 alt_lab=y_axis_labs,
#                 alt_lab_pos=y_axis_pos,
#                 alt_tick_pos=TRUE,
#                 line=0.8)

#   # AddLinearAxis(4, 1, 1, label="5'-paired nt",
#   #               alt_lab=rownames(R_mat),
#   #               alt_lab_pos=ymin:(ymax - 1) + 0.5,
#   #               alt_tick_pos=TRUE,
#   #               line=0.8)
#   # AddLinearAxis(1, 1, 1, label="3'-paired nt",
#   #               label_pos_ticks=TRUE,
#   #               alt_lab=colnames(R_mat),
#   #               alt_lab_y_dist=0.01,
#   #               alt_lab_pos=xmin:(xmax - 1) + 0.5,
#   #               alt_tick_pos=TRUE)
#   # Add the label for the seed if not plotting relative to seed kds.
#   col.inds <- floor(R_mat/-22*99 + 1)
#   # x_lab_pos <- -0.75
#   # Load the seaborn cubeHelix color palette, and conver the indeces from the
#   # and make them bounded between zero and 1
#   start_col <- 0.5
#   r_col <- -0.75
#   # Make distribution of 100 colors and a 101th gray color for NA indeces.
#   color.dist <- c(rev(cubeHelix(100, start=start_col, r=r_col, hue=0.8)),
#                   "white")
#   ## Don't need to change the max value to kMaxValueColor because none of the
#   ## values of the heatmap actually have this value, so it's better to not have
#   ## the key have this value at the top.
#   col.inds <- t(col.inds)
#   col.inds[which(col.inds > 100)] <- 100
#   col.inds[which(is.na(col.inds))] <- 101
#   col.inds[which(col.inds < 1)] <- 1
#   # Make the color rectangle.
#   rect(xlefts, ybottom, xright, ytop, col=color.dist[col.inds], lwd=0.4, 
#        xpd=NA, border="white")
#   # Make the key.
#   if (key) {
#     mir_scale <- (mir_length - 11)/(21 - 11)
#     num_boxes <- 100
#     y_div <- (ymax - ymin)/num_boxes
#     kpl <- xmax + 5*mir_scale # Left-hand position of the key
#     kw <- mir_scale                # Width of the key
#     rect(xleft=kpl,
#          ybottom=seq(num_boxes)*y_div - y_div,
#          xright=kpl + kw,
#          ytop=seq(num_boxes)*y_div,
#          col=color.dist[1:100], xpd=NA, border=NA)
#     # Generate the axis for the legend and the label
#     bottom <- ymin + y_div/2
#     top <- ymax - y_div/2
#     labels <- c(0, -5, -10, -15, -20, -22)
#     pos_labels <- labels
#     centered_labels <- pos_labels - pos_labels[1]
#     norm_labels <- centered_labels/(centered_labels[length(centered_labels)])
#     height_span <- top - bottom
#     pos_labels <- norm_labels*height_span + y_div/2
#     axis(4, at=pos_labels, labels=labels, lwd=par()$lwd, pos=kpl + kw, xpd=NA)
#     text(x=kpl + kw + 2.7*mir_scale, y=ymax/2,
#          labels=expression(Delta*italic(G)~"(kcal/mol)"), srt=270, xpd=NA)
#   }  
#   # # Add the label indicating how many nucleotides of pairing.
#   mtext(text="Predicted stability", side=3, at=xmax, adj=1,
#         cex=par("cex"))
#   # If the `mirna_label` conditional is true, add the label saying which miRNA
#   # is being looked at.
#   if (mirna_label) {
#     if (mirna == "let-7a-21nt") {
#       mirna_txt <- "let-7a"
#     } else if (mirna == "let-7a_minus1") {
#       mirna_txt <- "let-7a(-1)"
#     } else if (mirna == "let-7a_plus1") {
#       mirna_txt <- "let-7a(+1)"
#     } else if (mirna == "let-7a_miR-155") {
#       mirna_txt <- "let-7a-miR-155"
#     } else if (mirna == "miR-155_let-7a") {
#       mirna_txt <- "miR-155-let-7a"
#     } else {
#       mirna_txt <- mirna
#     }
#     x <- GetPlotFractionalCoords(0.125, 0.5)[1]
#     mtext(text=mirna_txt, side=3, line=0.8, at=x, adj=0, cex=par("cex"))
#   }

#   # Finish the plot.
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
#   return(R_mat)
# }

# # 4E ___________________________________________________________________________
# PlotKdFoldChangeAgainstDeltaG <- function(
#   mirna, experiment, n_constant=3, offset_lim=c(-4, 16), len_lim=c(4, 11),
#   pos_lim=c(9, 23), pos_min=23, pos_max=9, pos_inv=FALSE, corrected_kds=TRUE,
#   model_values=TRUE, log_plus_one=FALSE, lambda_p=0.1, makeglobalmodel=TRUE,
#   offset=0, delG=FALSE, deldelG=FALSE, deldelGleft=FALSE, lengthnorm=TRUE,
#   just_complement=TRUE, wobble=TRUE, xlabels=TRUE, mirna_label=TRUE,
#   height=3, width=3, xpos=20, ypos=20, pdf.plot=FALSE
# ) {
#   R <- 1.987e-3 # in kcal K-1 mol-1
#   T <- 310.15 # in K  

#   if (model_values) {
#     if (makeglobalmodel) {
#       model <<- SubfunctionCall(FitPairingAndOffsetModelWithError)
#     }    
#     pairing_mat <- model$pairing$MLE
#   } else {
#     R_mat <- SubfunctionCall(MakePairingMatrix)
#     pairing_mat <- R_mat  
#   }
#   dG_mat <- SubfunctionCall(GetDeltaGPredictedMatrix)
#   if (lengthnorm) {
#     pairing_mat <- NormalizeByLength(pairing_mat)
#     dG_mat <- NormalizeByLength(dG_mat)    
#   }

#   if (delG) {
#     pairing_mat <- -pairing_mat/log10(exp(1))*R*T
#     logstr <- ""
#     ylab <- expression("Model"~Delta*Delta*italic(G)~"(kcal/mol)")
#     inv <- "xy"
#   } else {
#     if (log_plus_one & model_values) {
#       pairing_mat <- exp(pairing_mat)
#     } else {
#       pairing_mat <- 10^pairing_mat
#     }
#     logstr <- "y"
#     if (model_values) {
#       ylab <- expression("Model"~italic(K)[D]~"fold change")
#       ylab <- "Pairing coefficients"
#     } else {
#       ylab <- expression("Measured"~italic(K)[D]~"fold change")
#     }
#     inv <- "x"
#   }

#   if (pos_inv) {
#     if (pos_min != 23) {
#       rows_use <- 1:nrow(pairing_mat)
#       cols_use <- which(as.numeric(colnames(pairing_mat)) > pos_min)
#     } else if (pos_max != 9) {
#       rows_use <- which(as.numeric(rownames(pairing_mat)) < pos_max)
#       cols_use <- 1:ncol(pairing_mat)
#     } else {
#       rows_use <- 1:nrow(pairing_mat)
#       cols_use <- 1:ncol(pairing_mat)
#     }
#   } else {
#     if (pos_min != 23) {
#       rows_use <- 1:nrow(pairing_mat)
#       cols_use <- which(as.numeric(colnames(pairing_mat)) <= pos_min)
#     } else if (pos_max != 9) {
#       rows_use <- which(as.numeric(rownames(pairing_mat)) >= pos_max)
#       cols_use <- 1:ncol(pairing_mat)
#     } else {
#       rows_use <- 1:nrow(pairing_mat)
#       cols_use <- 1:ncol(pairing_mat)
#     }
#   }
#   pairing_mat <- pairing_mat[rows_use, cols_use, drop=FALSE]
#   dG_mat <- dG_mat[rows_use, cols_use, drop=FALSE]
#   pairing_mat_global <<- pairing_mat
#   dG_mat_global <<- dG_mat
#   ################################ Create the plot #############################
#   SubfunctionCall(FigureSaveFile2)
#   if (lengthnorm) {
#     xmin <- -6
#     xmax <- 6
#     x_lab_space <- 2
#     if (delG) {
#       ymin <- -2
#       ymax <- 2
#     } else {
#       ymin <- 0.1
#       ymax <- 10
#     }
#   } else {
#     xmin <- -25
#     xmax <- 0
#     x_lab_space <- 5
#     if (delG) {
#       ymin <- -4
#       ymax <- 0
#     } else {
#       ymin <- 1
#       ymax <- 1000
#     }
#   }
#   xlab <- expression("Predicted"~Delta*italic(G)~"(kcal/mol)")
#   par(mar=c(3, 3, 1, 1))
#   BlankPlot(inv=inv, log=logstr)
#   # AddLinearAxis(1, 1, 2, label=xlab, alt_lab_pos=c(-6, -4, -2, 0, 2, 4, 6))
#   AddLinearAxis(1, 1, x_lab_space, label=xlab)
#   ################ Add the trend line and confidence interval ##################
#   x.line <- seq(xmin, xmax, length.out=100)
#   if (delG) {
#     AddLinearAxis(2, 1, 0.5, label=ylab)
#     segments(max(xmin, ymin), max(xmin, ymin), min(xmax, ymax), min(xmax, ymax),
#              lty=2)
#     linmodel <- lm(as.numeric(pairing_mat) ~ as.numeric(dG_mat),
#                    na.action=na.omit)
#     predict <- predict(linmodel, data.frame(dG_mat=x.line),
#                        interval="confidence")
#     predict[which(predict[, 2] < ymin), 2] <- ymin
#     predict[which(predict[, 3] < ymin), 3] <- ymin
#     predict[which(predict[, 2] > ymax), 2] <- ymax
#     predict[which(predict[, 3] > ymax), 3] <- ymax
#     segments(x0=-6, y0=1, x1=6)
#     polygon(c(x.line, rev(x.line)), c(predict[,2], rev(predict[, 3])),
#             col=ConvertRColortoRGB("gray", alpha=0.3), border=NA)

#   } else {
#     AddLogAxis(2, label=ylab)
#     linmodel <- lm(as.numeric(log10(pairing_mat)) ~ as.numeric(dG_mat),
#                    na.action=na.omit)
#     predict <- predict(linmodel, data.frame(dG_mat=x.line),
#                        interval="confidence")
#     # This is required to make the intervals bounded once the clipping masks are
#     # removed in Illustrator.
#     predict[which(10^predict[, 2] < ymin), 2] <- log10(ymin)
#     predict[which(10^predict[, 3] < ymin), 3] <- log10(ymin)
#     predict[which(10^predict[, 2] > ymax), 2] <- log10(ymax)
#     predict[which(10^predict[, 3] > ymax), 3] <- log10(ymax)
#     polygon(c(x.line, rev(x.line)), c(10^predict[, 2], rev(10^predict[, 3])),
#             col=ConvertRColortoRGB("gray", alpha=0.3), border=NA)
#     if (!lengthnorm) {
#       ymax_use <- 700
#     } else {
#       ymax_use <- ymax
#     }
#     xmax_convert <- -log(ymax_use)*R*T
#     xmin_convert <- -log(ymin)*R*T
#     segments(xmin_convert, ymin, xmax_convert, ymax_use, lty=2)
#   }
#   ########################### Make the points ##################################
#   pos3p <- rep(rownames(pairing_mat), times=ncol(pairing_mat))
#   pos5p <- rep(colnames(pairing_mat), each=nrow(pairing_mat))
#   len_sites <- as.character(as.integer(pos3p) - as.integer(pos5p) + 1)
#   points(dG_mat, pairing_mat, col=kThrPLengthCols[len_sites])
#   ##################### Add correlation to the plot ############################
#   xy <- GetPlotFractionalCoords(1, 1, inv=inv)
#   if (!delG) {
#     pairing_mat <- log(pairing_mat)
#   }
#   AddCorrelationToPlot(xpos=xy[1], ypos=xy[2], x=c(dG_mat),
#                        y=c(pairing_mat), rsquared=TRUE, adj=1)
#   # Add miRNA label ############################################################
#   if (mirna_label) {
#     if (mirna == "let-7a-21nt") {
#       mirna_txt <- "let-7a"
#     } else if (mirna == "let-7a_minus1") {
#       mirna_txt <- "let-7a(-1)"
#     } else if (mirna == "let-7a_plus1") {
#       mirna_txt <- "let-7a(+1)"
#     } else if (mirna == "let-7a_miR-155") {
#       mirna_txt <- "let-7a-miR-155"
#     } else if (mirna == "miR-155_let-7a") {
#       mirna_txt <- "miR-155-let-7a"
#     } else {
#       mirna_txt <- mirna
#     }
#     xy <- GetPlotFractionalCoords(0.025, 1, inv=inv)
#     text(x=xy[1], y=xy[2], labels=mirna_txt, adj=0, xpd=NA)
#   }
#   # Finish the plot. ###########################################################
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }

# # 4F____________________________________________________________________________
# PlotMismatchCoefficients <- function(
#   mirna, experiment, n_constant=3, sitelist="progthrp", corrected_kds=TRUE,
#   combined=TRUE, buffer=FALSE, makeglobalmodel=TRUE, fullpairs=FALSE,
#   len_lim=c(4, 11), pos_lim=c(9, 23), offset_lim=c(-4, 16), supp_base=FALSE,
#   key=TRUE, mirna_label=TRUE, kd_fc=TRUE, exponential=FALSE, intercept=FALSE,
#   log_plus_one=FALSE, lambda_p=0.1, empirical=FALSE, percentile=0.25, ylab=TRUE,
#   height=3.5, width=3.5, xpos=20, ypos=20, pdf.plot=FALSE
# ) {
#   # Load the collapsed data (for the mismatch-only).
#   if (empirical) {
#     out_mat <<- SubfunctionCall(GetEmpiricalBestSeedMismatchCoefficients)
#     n_cols <- round(percentile*ncol(out_mat))
#     data <- matrix(rowMeans(out_mat[, 1:n_cols]), ncol=1)
#     se_means <- apply(out_mat[, 1: n_cols], 1, sd)/sqrt(n_cols)
#     lower_CI <- data - se_means
#     upper_CI <- data + se_means
#     data <- cbind(data, lower_CI, upper_CI)
#     rownames(data) <- rownames(out_mat)
#     empirical_mismatch_global <<- data
#   } else {
#     if (makeglobalmodel) {
#       print("in makeglobalmodel conditional")
#       model <<- SubfunctionCall(FitPairingOffsetAndMismatchModelWithError)  
#     }
#     mismatch_model_global <<- model
#     # model <- mismatch_model_global
#     data <- model$mm    
#   }
#   # Make the plot and define the limits.
#   # if (!ylab) width <- 3.2
#   # if (!ylab) width <- 3.15
#   if (!ylab) width <- 3.175
#   SubfunctionCall(FigureSaveFile2)
#   if (ylab) par(mar=c(3, 3, 2, 2))
#   else      par(mar=c(3, 1.6, 2, 2))
#   xmin <- 0
#   xmax <- 23
#   if (empirical | log_plus_one) {
#     ymin <- -1
#     ymax <- 1
#   } else if (mirna == "miR-7-23nt") {
#     ymin <- 0
#     ymax <- 2
#   } else {
#     ymin <- 0.5
#     ymax <- 1.5    
#   }
#   BlankPlot()
#   col.inds <- floor((data[, 1] - ymin)/(ymax - ymin)*99 + 1)
#   col.inds[which(col.inds <= 0)] <- 1  
#   # Make distribution of 100 colors and a 101th gray color for NA indeces.
#   color.dist <- rev(colorspace::diverging_hcl(100, h=c(180, 50), c=80,
#                                               l=c(20, 95), power=c(0.7, 1.3)))
#   if (supp_base) {
#     xleft <- 0:6*3 + c(0, 1, 2, 3, 4, 5, 6)
#     xright <- xleft + 2
#   } else {
#     xleft <- 0:17 + 0.05 + rep(c(0, 1, 2, 3, 4, 5), each=3)
#     xright <- xleft + 1
#   }
#   # Get the left- and right-hand positions of the rectangles.
#   # Make the barplot rectangles ################################################
#   rect(xleft=xleft, ybottom=mean(c(ymin, ymax)), xright=xright, ytop=data[, 1],
#        col=color.dist[col.inds], border=NA)
#   # if (!empirical) {
#     arrows(x0=(xleft + xright)/2, y0=data[, 2], x1=(xleft + xright)/2,
#            y1=data[, 3], length=0.02, angle=90, code=3, xpd=NA)
#   # }

#   # Get the each of the nucleotide mismatches.
#   if (!supp_base) {
#     nucs <- gsub("T", "U", gsub("8mer-mm(.)(.)", replacement="\\1", rownames(data),
#                                 perl=TRUE))
#     pos <- unique(gsub("8mer-mm(.)(.)", replacement="\\2", rownames(data), perl=TRUE))

#     # Get miRNA nucleotides 2-7.
#     mirna_nucs <- unlist(strsplit(kMirnaSeqs[mirna], split=""))[2:7]
#     pos_cols <- c()
#     for (i in 0:5) {
#       nuc_ind <- i + 1
#       mirna_nuc <- mirna_nucs[nuc_ind]
#       mirna_nucs_start <- i*3 + 1
#       mirna_nucs_stop <- mirna_nucs_start + 2
#       cols_i <- rep("black", 3)
#       if (mirna_nuc == "G") {
#         ind_wobble <- which(nucs[mirna_nucs_start:mirna_nucs_stop] == "U")
#         cols_i[ind_wobble] <- "blue"
#       } else if (mirna_nuc == "U") {
#         ind_wobble <- which(nucs[mirna_nucs_start:mirna_nucs_stop] == "G")
#         cols_i[ind_wobble] <- "red"
#       }
#       pos_cols <- c(pos_cols, cols_i)
#     }
#     # Add a negative, zero, and positive x-shift to the nucleotide letters to
#     # make them more spread out. 
#     # Get the y coordinate for the nucleotide labels, then print the labels.
#     # The frac_scale variable deals with the different size mismatch plots in
#     # Figure 5 in comparison to Figures 4 and 6.
#     if (height == 4) {
#       frac_scale <- 0.85
#     } else {
#       frac_scale <- 1
#     }
#     x_shifts <- rep(c(-0.2, 0, 0.2), 6)*frac_scale
#     xy <- GetPlotFractionalCoords(0.05, -0.04*frac_scale)
#     text(x=(xleft + xright)/2 + x_shifts, y=xy[2], labels=nucs, col=pos_cols, xpd=NA)
#     # Get the y coordinate for the seed positions, then print the positions.
#     xy <- GetPlotFractionalCoords(0.05, -0.11*frac_scale)
#     text(x=0:5*4 + 1.5, y=xy[2], labels=pos, xpd=NA)
#     # Add mismatch label to the plot:
#     mtext(text="Mismatch identity & position", side=1, line=1.4, cex=par("cex"), xpd=NA)
#   } else {
#     xy <- GetPlotFractionalCoords(0.05, 0)
#     if (supp_base) {
#       text(x=0:6*4 + 1.5, y=xy[2], labels=rownames(data), srt=45, xpd=NA, adj=c(1, 1))
#     } else {
#       text(x=0:5*4 + 1.5, y=xy[2], labels=rownames(data), srt=45, xpd=NA, adj=c(1, 1))
#     }
#   }

#   if (mirna_label) {
#     if (mirna == "let-7a-21nt") {
#       mirna_txt <- "let-7a"
#     } else if (mirna == "let-7a_minus1") {
#       mirna_txt <- "let-7a(-1)"
#     } else if (mirna == "let-7a_plus1") {
#       mirna_txt <- "let-7a(+1)"
#     } else if (mirna == "let-7a_miR-155") {
#       mirna_txt <- "let-7a-miR-155"
#     } else if (mirna == "miR-155_let-7a") {
#       mirna_txt <- "miR-155-let-7a"
#     } else {
#       mirna_txt <- mirna
#     }
#     # x <- GetPlotFractionalCoords(0.0125, 0.5)[1]
#     # mtext(text=mirna_txt, side=3, line=0, at=x, adj=c(0, 0), cex=par("cex"),
#     #       xpd=NA)
#     xy <- GetPlotFractionalCoords(0.025, 1.1)
#     text(xy[1], xy[2], labels=sprintf(mirna_txt), adj=c(0, 1), xpd=NA)  

#   }
#   ## Add label explaining that these are the model mismatch coefficients.
#   xy <- GetPlotFractionalCoords(0.025, 1.025)
#   if (empirical) {
#     AddLinearAxis(2, 0.1, 0.2, label=expression(Delta~italic(K)[D]~"fold change ("*log[10]*")"))
#     text(xy[1], xy[2], labels="Empirical variation", adj=c(0, 1), xpd=NA)
#     xy <- GetPlotFractionalCoords(0.025, 0.95)
#     text(xy[1], xy[2], labels=sprintf("Top %sth percentile", round(percentile*100, 0)), adj=c(0, 1), xpd=NA)
#   } else {
#     if (ylab) y_label <- "Mismatch coefficients"
#     else      y_label <- ""
#     AddLinearAxis(2, 0.1, 0.2, label=y_label)
#     # text(xy[1], xy[2], labels="Mismatch coefficients", adj=c(0, 1), xpd=NA)  
#     # xy <- GetPlotFractionalCoords(0.025, 0.95)
#     # if (intercept) {
#     #   label_i <- "With intercept term"
#     # } else {
#     #   label_i <- "No intercept term"
#     # }
#     # text(xy[1], xy[2], labels=label_i, adj=c(0, 1), xpd=NA)  
#     # xy <- GetPlotFractionalCoords(0.025, 0.875)
#     # if (exponential) {
#     #   label_i <- "Exponential coefs"
#     # } else {
#     #   label_i <- "Linear coefs"
#     # }
#     # text(xy[1], xy[2], labels=label_i, adj=c(0, 1), xpd=NA)  
#   }
#   # mtext(text="Mismatch coefficients", side=3, line=-1, at=x, adj=c(0, 0),
#   #       cex=par("cex"), xpd=NA)
#   # Finish the plot.
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }

# # 4G____________________________________________________________________________
# PlotAllProgrammedbraryMismatchCoefficientsAgainstKds <- function(
#   n_constant=3, sitelist="progthrp", corrected_kds=TRUE, supp_base=FALSE,
#   sumseed=FALSE, key=TRUE, len_lim=c(4, 11), pos_lim=c(9, 23),
#   offset_lim=c(-4, 16), kd_fc=TRUE, intercept=FALSE, exponential=FALSE,
#   log_plus_one=FALSE, empirical=FALSE, single_error=FALSE, percentile=0.25,
#   height=3.5, width=3.5, xpos=20, ypos=20, pdf.plot=FALSE
# ) {
#   print("Inside PlotAllProgrammedLibrarySeedCoefficientsAgainstKds")
#   # Load the collapsed data (for the mismatch-only).
#   mirnas_use <- c("let-7a-21nt", "miR-1", "miR-155")
#   exps <- c("equil_c2_nb", "equil_c_nb", "equil_sc_nb")
#   names(exps) <- mirnas_use
#   data_all <<- sapply(mirnas_use, function(mirna) {
#     experiment <- exps[mirna]
#     if (empirical) {
#       out_mat <- SubfunctionCall(GetEmpiricalBestSeedMismatchCoefficients)
#       print(out_mat)
#       print("Just made out_mat")
#       data <- rowMeans(out_mat[, 1:round(percentile*ncol(out_mat))])
#     } else if (single_error) {
#       model <- SubfunctionCall(FitPairingOffsetAndMismatchModelSingle)
#       data <- model$coefs[grep("8mer", rownames(model$coefs)), 1]
#     } else {
#       model <- SubfunctionCall(FitPairingOffsetAndMismatchModelWithError)
#       data <- model$mm[, 1]
#     }
#   })

#   kds_all <<- sapply(mirnas_use, function(mirna) {
#     experiment <- exps[mirna]
#     kds <- SubfunctionCall(EquilPars)
#     kds_use <- kds[sprintf("%s_Kd", GetAll8merMmSites(mirna)), 2]
#     names(kds_use) <- GetAll8merMmSites(mirna)
#     kds_use
#   })

#   print(data_all)
#   print(kds_all)
#   # Make the plot and define the limits.
#   SubfunctionCall(FigureSaveFile2)
#   par(mar=c(3, 3, 2, 2))
#   xmin <- 1e-1
#   xmax <- 3
#   if (empirical) {
#     ymin <- -1
#     ymax <- 1
#   } else if (log_plus_one) {
#     ymin <- -2
#     ymax <- 2
#   } else {
#     ymin <- 0
#     ymax <- 2
#   }
#   BlankPlot(log="x")

#   x.line <- seq(log10(xmin), log10(xmax), length.out=20)

#   log10kds <- log10(c(kds_all))
#   linmodel <- lm(c(data_all) ~ log10kds)
#   b <- linmodel$coefficients[1]
#   m <- linmodel$coefficients[2]
#   y.line <- m*x.line + b
#   coefs <- summary(linmodel)$coefficients[2, 1:2]

#   predict_y <- predict(linmodel, data.frame(log10kds=x.line), interval = 'confidence')
#   polygon(c(10^x.line, rev(10^x.line)), c(predict_y[,2], rev(predict_y[, 3])),
#           col=ConvertRColortoRGB("gray", alpha=0.3), border=NA)
#   lines(10^x.line, y.line, lty=2)

#   Points(c(kds_all), c(data_all), col=rep(kMirnaColors[c(2, 1, 3)], each=nrow(data_all)))
#   AddLogAxis(1, label="Relative Kd")
#   xy <- GetPlotFractionalCoords(0.025, 1.025, log="x")
#   if (empirical) {
#     AddLinearAxis(2, 0.1, 0.2, label=expression("Average difference in"~italic(K)[D]~"("*log[10]*")"))
#   } else {
#     AddLinearAxis(2, 0.1, 0.2, label="Mismatch coefficients")
#     # text(xy[1], xy[2], labels="Mismatch coefficients", adj=c(0, 0.95), xpd=NA)  
#     # xy <- GetPlotFractionalCoords(0.025, 0.95, log="x")
#     # if (intercept) {
#     #   label_i <- "With intercept term"
#     # } else {
#     #   label_i <- "No intercept term"
#     # }
#     # text(xy[1], xy[2], labels=label_i, adj=c(0, 1), xpd=NA)  
#     # xy <- GetPlotFractionalCoords(0.025, 0.875, log="x")
#     # if (exponential) {
#     #   label_i <- "Exponential coefs"
#     # } else {
#     #   label_i <- "Linear coefs"
#     # }
#     # text(xy[1], xy[2], labels=label_i, adj=c(0, 1), xpd=NA)  
#   }
#   xy <- GetPlotFractionalCoords(0, 0, log="x")
#   legend_vals <- c("let-7a", "miR-1", "miR-155", "miR-124", "lsy-6", "miR-7")
#   Legend(
#     xy, legend=legend_vals, col=kMirnaColors[legend_vals], yjust=0, ncol=2, text.width=1.5)
#   # mtext(text="Mismatch coefficients", side=3, line=-1, at=x, adj=c(0, 0),
#   #       cex=par("cex"), xpd=NA)
#   # Finish the plot.
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }

# # 4H____________________________________________________________________________
# PlotAllRandomLibrarySeedCoefficientsAgainstKds <- function(
#   experiment="equilibrium", n_constant=3, sitelist="randthrp_comp",
#   corrected_kds=FALSE, combined=TRUE, buffer=FALSE, supp_base=TRUE,
#   sumseed=FALSE, key=TRUE, len_lim=c(4, 11), pos_lim=c(9, 23),
#   offset_lim=c(-4, 16), pos_5p_lim=FALSE, kd_fc=TRUE, intercept=FALSE,
#   exponential=FALSE, empirical=FALSE, percentile=0.25, height=3.5, width=3.5,
#   xpos=20, ypos=20, pdf.plot=FALSE
# ) {
#   print("Inside PlotAllRandomLibrarySeedCoefficientsAgainstKds")
#   # Load the collapsed data (for the mismatch-only).
#   data_all <- sapply(kMirnas[1:6], function(mirna) {
#     if (mirna == "miR-1") {

#       combined <- FALSE
#       buffer <- TRUE
#     } else if (mirna == "miR-7-23nt") {
#       experiment <- "equilibrium2_nb"
#       combined <- FALSE
#     }
#     if (empirical) {
#       out_mat <- SubfunctionCall(GetEmpiricalBestSeedMismatchCoefficients)
#       print(out_mat)
#       out_mat <<- out_mat
#       print("Just made out_mat")
#       data <- rowMeans(out_mat[, 1:round(percentile*ncol(out_mat))])
#     } else {
#       model <- SubfunctionCall(FitPairingOffsetAndMismatchModel)
#       print(model$offsets)
#       data <- model$mm    
#     }
#   })
#   kds_all <- sapply(kMirnas[1:6], function(mirna) {
#     if (mirna == "miR-1") {
#       combined <- FALSE
#       buffer <- TRUE
#     } else if (mirna == "miR-7-23nt") {
#       experiment <- "equilibrium2_nb"
#       combined <- FALSE
#     }
#       kds <- SubfunctionCall(EquilPars)
#       kds_use <- kds[sprintf("%s_Kd", c(kSeedSites, "Comp")), 2]
#       names(kds_use) <- c(kSeedSites, "Comp")
#     kds_use
#   })
#   data_all <- data_all[-nrow(data_all), ]
#   kds_all <- kds_all[-nrow(kds_all), ]
#   data_all <<- data_all
#   kds_all <<- kds_all
#   # Make the plot and define the limits.
#   SubfunctionCall(FigureSaveFile2)
#   par(mar=c(3, 3, 2, 2))
#   xmin <- 1e-4
#   xmax <- 1
#   if (empirical) {
#     ymin <- -1
#     ymax <- 1
#   } else {
#     ymin <- 0
#     ymax <- 2
#   }
#   BlankPlot(log="x")

#   x.line <- seq(log10(xmin), log10(xmax), length.out=20)

#   log10kds <- log10(c(kds_all))
#   linmodel <- lm(c(data_all) ~ log10kds)
#   b <- linmodel$coefficients[1]
#   m <- linmodel$coefficients[2]
#   y.line <- m*x.line + b
#   coefs <- summary(linmodel)$coefficients[2, 1:2]

#   predict_y <- predict(linmodel, data.frame(log10kds=x.line), interval = 'confidence')
#   polygon(c(10^x.line, rev(10^x.line)), c(predict_y[,2], rev(predict_y[, 3])), col=ConvertRColortoRGB("gray", alpha=0.3), border=NA)
#   lines(10^x.line, y.line, lty=2)

#   Points(c(kds_all), c(data_all), col=rep(kMirnaColors, each=nrow(data_all)))
#   AddLogAxis(1, label="Relative Kd")
#   xy <- GetPlotFractionalCoords(0.025, 1.025, log="x")
#   if (empirical) {
#     AddLinearAxis(2, 0.1, 0.2, label=expression("Average difference in"~italic(K)[D]~"("*log[10]*")"))
#   } else {
#     AddLinearAxis(2, 0.1, 0.2, label=expression(italic(K)[D]~"scaling factor ("*log[10]*")"))
#     text(xy[1], xy[2], labels="Mismatch coefficients", adj=c(0, 0.95), xpd=NA)  
#     xy <- GetPlotFractionalCoords(0.025, 0.95, log="x")
#     if (intercept) {
#       label_i <- "With intercept term"
#     } else {
#       label_i <- "No intercept term"
#     }
#     text(xy[1], xy[2], labels=label_i, adj=c(0, 1), xpd=NA)  
#     xy <- GetPlotFractionalCoords(0.025, 0.875, log="x")
#     if (exponential) {
#       label_i <- "Exponential coefs"
#     } else {
#       label_i <- "Linear coefs"
#     }
#     text(xy[1], xy[2], labels=label_i, adj=c(0, 1), xpd=NA)  
#   }
#   xy <- GetPlotFractionalCoords(0, 0, log="x")
#   legend_vals <- c("let-7a", "miR-1", "miR-155", "miR-124", "lsy-6", "miR-7")
#   Legend(
#     xy, legend=legend_vals,
#     col=kMirnaColors[legend_vals], yjust=0, ncol=2, text.width=1.5)
#   # mtext(text="Mismatch coefficients", side=3, line=-1, at=x, adj=c(0, 0),
#   #       cex=par("cex"), xpd=NA)
#   # Finish the plot.
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }

# ################################################################################
# # FIGURE 5
# ################################################################################

# # 5D-G__________________________________________________________________________
# PlotChimeraModelCoefficientsScatter <- function(
#   mirna, experiment, n_constant=3, offset_lim=c(-4, 16), len_lim=c(4, 11),
#   pos_lim=c(9, 23), makeglobalmodel=TRUE, sitelist="progthrp_suppcomp",
#   corrected_kds=TRUE, sumseed=FALSE, supp_base=FALSE, delG=TRUE,
#   pairing_coefs=TRUE, seed_compare=FALSE, reg_compare=1, log_plus_one=FALSE,
#   lengthnorm=TRUE, mirna_label=FALSE, legend=FALSE, height=4, width=4, xpos=20,
#   ypos=20, pdf.plot=FALSE
# ) {
#   if (mirna == "let-7a-21nt") {
#     if (seed_compare) {
#       mirna_chim <- "let-7a_miR-155"
#       mirna_chim_str <- "let-7a-miR-155"
#     } else {
#       mirna_chim <- "miR-155_let-7a"
#       mirna_chim_str <- "miR-155-let-7a"
#     }
#   } else if (mirna == "miR-155") {
#     if (seed_compare) {
#       mirna_chim <- "miR-155_let-7a"
#       mirna_chim_str <- "miR-155-let-7a"
#     } else {
#       mirna_chim <- "let-7a_miR-155"
#       mirna_chim_str <- "let-7a-miR-155"
#     }
#   }
#   exp_chim <- "equil_c_nb"
#   if (makeglobalmodel) {
#     model_base <<- SubfunctionCall(FitPairingAndOffsetModelWithError)
#     model_chim <<- SubfunctionCall(FitPairingAndOffsetModelWithError,
#                                    mirna=mirna_chim, experiment=exp_chim)
#   }
#   # Make the color index scale.  
#   SubfunctionCall(FigureSaveFile2)
#   if (pairing_coefs) {
#     if (lengthnorm) {
#       xmin <- 0.1
#       xmax <- 10
#     } else {
#       xmin <- 1
#       xmax <- 400      
#     }
#     logstr <- "xy"
#   } else {
#     if (log_plus_one) {
#       xmin <- -2.5
#       xmax <- 0.5
#     } else {
#       xmin <- 0
#       xmax <- 1
#     }
#     logstr <- ""
#   }
#   ymin <- xmin
#   ymax <- xmax    
#   BlankPlot(log=logstr)
#   if (pairing_coefs) {
#     pairing_mat <- model_base$pairing$MLE
#     pairing_mat_lCI <- model_base$pairing$LowerCI
#     pairing_mat_uCI <- model_base$pairing$UpperCI

#     pairing_mat_chim <- model_chim$pairing$MLE
#     pairing_mat_chim_lCI <- model_chim$pairing$LowerCI
#     pairing_mat_chim_uCI <- model_chim$pairing$UpperCI

#     pairing_mat <<- pairing_mat
#     pairing_mat_chim <<- pairing_mat_chim
#     n_base <- nrow(pairing_mat)
#     n_chim <- nrow(pairing_mat_chim)
#     if (lengthnorm) {
#       # Normalize the rows
#       pairing_mat_norm <- NormalizeByLength(pairing_mat)
#       pairing_mat_chim_norm <- NormalizeByLength(pairing_mat_chim)
#       # Adjust the error values
#       pairing_mat_lCI <- pairing_mat_lCI - pairing_mat + pairing_mat_norm
#       pairing_mat_uCI <- pairing_mat_uCI - pairing_mat + pairing_mat_norm

#       pairing_mat_chim_lCI <- pairing_mat_chim_lCI - pairing_mat_chim + pairing_mat_chim_norm
#       pairing_mat_chim_uCI <- pairing_mat_chim_uCI - pairing_mat_chim + pairing_mat_chim_norm

#       pairing_mat <- pairing_mat_norm
#       pairing_mat_chim <- pairing_mat_chim_norm

#     }
#     if (n_base < n_chim) {
#       pairing_mat_chim <- pairing_mat_chim[(1:n_base + reg_compare - 1), (1:n_base + reg_compare - 1)]
#       pairing_mat_chim_lCI <- pairing_mat_chim_lCI[(1:n_base + reg_compare - 1), (1:n_base + reg_compare - 1)]
#       pairing_mat_chim_uCI <- pairing_mat_chim_uCI[(1:n_base + reg_compare - 1), (1:n_base + reg_compare - 1)]
#     } else if (n_base > n_chim) {
#       pairing_mat <- pairing_mat[(1:n_chim + reg_compare - 1), (1:n_chim + reg_compare - 1)]
#       pairing_mat_lCI <- pairing_mat_lCI[(1:n_chim + reg_compare - 1), (1:n_chim + reg_compare - 1)]
#       pairing_mat_uCI <- pairing_mat_uCI[(1:n_chim + reg_compare - 1), (1:n_chim + reg_compare - 1)]
#     }
#     # Portion of figure that calculates the linear model and its associated
#     # error.
#     x.line <- seq(log10(xmin), log10(xmax), length=100)

#     x_lm <- c(pairing_mat)
#     y_lm <- c(pairing_mat_chim)
#     linmodel <- lm(y_lm ~ x_lm)

#     predict <- predict(linmodel, data.frame(x_lm=x.line), interval = 'confidence')
#     predict[which(10^predict[, 2] < ymin), 2] <- log10(ymin)
#     predict[which(10^predict[, 3] > ymax), 3] <- log10(ymax)
#     polygon(c(10^x.line, 10^rev(x.line)), c(10^predict[, 2], rev(10^predict[, 3])), col=ConvertRColortoRGB("gray90", alpha=0.3), border=NA)

#     ## TODO get rid of this.
#     if (log_plus_one) {
#       pairing_mat <- pairing_mat/log(10)
#       pairing_mat_chim <- pairing_mat_chim/log(10)

#       pairing_mat_lCI <- pairing_mat_lCI/log(10)
#       pairing_mat_uCI <- pairing_mat_uCI/log(10)

#       pairing_mat_chim_lCI <- pairing_mat_chim_lCI/log(10)
#       pairing_mat_chim_uCI <- pairing_mat_chim_uCI/log(10)
#     }

#     x <- c(10^pairing_mat)
#     y <- c(10^pairing_mat_chim)
#     x_l <- c(10^pairing_mat_lCI)
#     x_r <- c(10^pairing_mat_uCI)
#     y_l <- c(10^pairing_mat_chim_lCI)
#     y_u <- c(10^pairing_mat_chim_uCI)



#     stops <- rep(rownames(pairing_mat), times=ncol(pairing_mat))
#     starts <- rep(colnames(pairing_mat), each=nrow(pairing_mat))
#     len_sites <- as.integer(stops) - as.integer(starts) + 1
#     cols <- kThrPLengthCols
#     cols_use <- cols[as.character(len_sites)]

#   } else {

#     x.line <- seq(xmin, xmax, length=100)

#     x <- model_base$offsets[, 1]
#     x_l <- model_base$offsets[, 2]
#     x_r <- model_base$offsets[, 3] 

#     y <- model_chim$offsets[, 1]
#     y_l <- model_chim$offsets[, 2]
#     y_u <- model_chim$offsets[, 3]

#     linmodel <- lm(y ~ x)

#     names(x) <- rownames(model_base$offsets)

#     predict <- predict(linmodel, data.frame(x=x.line), interval = 'confidence')
#     predict[which(predict[, 2] < ymin), 2] <- ymin
#     predict[which(predict[, 3] > ymax), 3] <- ymax

#     polygon(c(x.line, rev(x.line)), c(predict[,2], rev(predict[, 3])), col=ConvertRColortoRGB("gray90", alpha=0.3), border=NA)
#     r_col <- -0.75
#     color.dist <- cubeHelix(length(x) + 8, start=0.25*pi, r=r_col, hue=0.8)[5:(length(x) + 4)]
#     cols_use <- color.dist
#     names(cols_use) <- names(x)
#   }

#   print(linmodel)
#   linmodel_global <<- linmodel
#   if (pairing_coefs) {
#     AddLogAxis(1, label="Native miRNA")
#     AddLogAxis(2, label="Chimeric miRNA")
#     label_str <- "Normalized\npairing coefficients"
#   } else {
#     AddLinearAxis(1, 0.1, 0.2, label="Native miRNA")
#     AddLinearAxis(2, 0.1, 0.2, label="Chimeric miRNA")
#     label_str <- "Offset coefficients"
#   }
#   segments(x0=x_l, x1=x_r, y0=y, lwd=0.25, xpd=NA)
#   segments(x0=x, y0=y_l, y1=y_u, lwd=0.25, xpd=NA)

#   points(x, y, col=cols_use, xpd=NA)
#   # Label the miRNA in the plot.
#   if (mirna_label) {
#     if (mirna == "let-7a-21nt") {
#         mirna_txt <- "let-7a"
#     } else {
#       mirna_txt <- mirna
#     }
#     xy <- GetPlotFractionalCoords(-0.1, 1.14, log=logstr)
#     text(xy[1], xy[2], labels=sprintf("%s versus %s", mirna_chim_str, mirna_txt),
#          adj=c(0, 1), xpd=NA)  
#   }
#   ############## Label the coefficients ########################################
#   xy <- GetPlotFractionalCoords(0.025, 0.025, log=logstr)
#   text(xy[1], xy[2], labels=label_str, adj=c(0, 0), xpd=NA)
#   # Add the r-squared to the plot
#   if (pairing_coefs) {
#     x <- log(x)
#     y <- log(y)
#   }
#   ######################### Add Pearson r ######################################
#   xy <- GetPlotFractionalCoords(1, 0.025, log=logstr)
#   AddCorrelationToPlot(x=x, y=y, xpos=xy[1], ypos=xy[2], rsquared=FALSE,
#                        adj=c(1, 0))
#   ######################### Add legend #########################################
#   if (legend) {
#     if (pairing_coefs) {
#       legend_title <- "Pairing length (bp)"
#       legend_labels <- names(kThrPLengthCols)
#       legend_cols <- kThrPLengthCols
#       ncol <- 3
#     } else {
#       legend_title <- "Offset (nt)"
#       legend_labels <- names(cols_use)
#       legend_cols <- cols_use
#       ncol <- 4
#     }
#     print("making legend")
#     xy_title <- GetPlotFractionalCoords(fx=0.025, fy=1.06, log=logstr)
#     xy_legend <- GetPlotFractionalCoords(fx=0, fy=1.035, log=logstr)
#     text(xy_title[1], xy_title[2], labels=legend_title, adj=c(0, 1), xpd=NA)
#     Legend(
#         xy_legend, legend=legend_labels, col=legend_cols, x.intersp=0.45,
#         text.width=0.07, ncol=ncol)
#   }
#   ######################### Finish the plot ####################################
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }

# PlotChimeraMismatchModelCoefficientsScatter <- function(
#   mirna, experiment, n_constant=3, offset_lim=c(-4, 16), len_lim=c(4, 11), 
#   pos_lim=c(9, 23), makeglobalmodel=TRUE, sitelist="progthrp",
#   corrected_kds=TRUE, sumseed=FALSE, supp_base=FALSE, mirna_label=TRUE,
#   exponential=FALSE, intercept=FALSE, log_plus_one=FALSE, legend=FALSE,
#   height=4, width=4, xpos=20, ypos=20, pdf.plot=FALSE
# ) {
#   if (mirna == "let-7a-21nt") {
#     mirna_chim <- "let-7a_miR-155"
#     mirna_chim_str <- "let-7a-miR-155"
#   } else if (mirna == "miR-155") {
#     mirna_chim <- "miR-155_let-7a"
#     mirna_chim_str <- "miR-155-let-7a"
#   }
#   exp_chim <- "equil_c_nb"
#   if (makeglobalmodel) {
#     model_base <<- SubfunctionCall(FitPairingOffsetAndMismatchModelSingle)
#     model_chim <<- SubfunctionCall(FitPairingOffsetAndMismatchModelSingle,
#                                    mirna=mirna_chim, experiment=exp_chim)
#   }
#   # Make the ordered pairs:
#   x <<- model_base$coefs[grep("mer", rownames(model_base$coefs)), 1]
#   y <<- model_chim$coefs[grep("mer", rownames(model_chim$coefs)), 1]
#   SubfunctionCall(FigureSaveFile2)
#   if (log_plus_one) {
#     xmin <- -1
#     xmax <- 1
#   } else {
#     xmin <- 0.6
#     xmax <- 1.4
#   }
#   ymin <- xmin
#   ymax <- xmax    
#   BlankPlot()
#   AddLinearAxis(1, 0.1, 0.2, label="Native miRNA")
#   AddLinearAxis(2, 0.1, 0.2, label="Chimeric miRNA")
#   label_str <- "Mismatch coefficients"
#   # Pick the colors for the points #############################################
#   r_col <- -0.50
#   cols <- cubeHelix(length(x) + 8, start=0.25*pi, r=r_col,
#                     hue=0.8)[5:(length(x) + 4)]
#   print(x)
#   print(y)
#   names(cols) <- sapply(names(y), function(name_i) {
#     unlist(strsplit(name_i, split="mm"))[2]
#   })



#   nucs <- substr(names(cols), start=1, stop=1)
#   pos <- substr(names(cols), start=2, stop=2)
#   pch_use <- c(`A`=1, `C`=2, `G`=3, `T`=4)
#   col_use <- c("#D73027", "#FC8D59", "#FEE090", "#E0F3F8", "#91BFDB", "#4575B4")

#   # col_use <- c("#ffffb2",
#   # "#fed976",
#   # "#feb24c",
#   # "#fd8d3c",
#   # "#f03b20",
#   # "#bd0026")

#   # col_use <- c(`2`="#7BBCE7",
#   #               `3`="#88A5DD",
#   #               `4`="#9B8AC4",
#   #               `5`="#9A709E",
#   #               `6`="#805770",
#   #               `7`="#46353A")

#   col_use <- c(`2`="#1965B0", #10
#                `3`="#7BAFDE", #12
#                `4`="#4EB265", #15
#                `5`="#CAE0AB", #17
#                # `6`="#F7F056", #18
#                # `6`="#EE8026", #23
#                `6`="#F6C141", #23
#                `7`="#DC050C") #26

#   col_use <- c(`2`="#0077BB", #blue
#                `3`="#33BBEE", #cyan
#                `4`="#009988", #teal
#                `5`="#EE7733", #orange
#                `6`="#CC3311", #red
#                `7`="#EE3377") #magenta

#   col_use <- c(`2`="#4477AA", #blue
#                `3`="#66CCEE", #cyan
#                `4`="#228833", #green
#                `5`="#CCBB44", #yellow
#                `6`="#EE6677", #red
#                `7`="#AA3377") #purple


#   print(col_use[pos])
#   print(pch_use[nucs])
#   # Plot the points ############################################################
#   segments(x0=xmin, y0=1, x1=xmax, lty=2, col="gray90")
#   segments(x0=1, y0=ymin, y1=ymax, lty=2, col="gray90")
#   points(x, y, col=col_use[pos], pch=pch_use[nucs], xpd=NA)
#   # Label the miRNA in the plot. ###############################################
#   if (mirna_label) {
#     if (mirna == "let-7a-21nt") {
#       mirna_txt <- "let-7a"
#     } else {
#       mirna_txt <- mirna
#     }
#     xy <- GetPlotFractionalCoords(0.025, 1.1)
#     text(xy[1], xy[2], labels=sprintf("%s versus %s", mirna_chim_str,
#                                       mirna_txt),
#          adj=c(0, 1), xpd=NA)  
#   }
#   # Label if these are the pairing coefficients or the offset coefficients. ####
#   xy <- GetPlotFractionalCoords(0.025, 0.025)
#   text(xy[1], xy[2], labels=label_str, adj=c(0, 0), xpd=NA)
#   # Add the r-squared to the plot
#   xy <- GetPlotFractionalCoords(1, 0.025)
#   AddCorrelationToPlot(x=x, y=y, xpos=xy[1], ypos=xy[2], rsquared=TRUE,
#                        adj=c(1, 0))

#   if (legend) {

#     xy_title <- GetPlotFractionalCoords(0.025, 1)
#     text(xy_title[1], xy_title[2], labels="M.m. pos.", adj=c(0, 1), xpd=NA)
#     xy_legend <- GetPlotFractionalCoords(0, 0.975)
#     # Legend(xy_legend, legend=unique(pos), col=col_use[unique(pos)])
#     legend(x=xy_legend[1], y=xy_legend[2], unique(pos), col=col_use[unique(pos)],
#            lwd=1, bty="n", xpd=NA, ncol=2, seg.len=0.6, x.intersp=0.5,
#            y.intersp=0.8)

#     xy_title <- GetPlotFractionalCoords(0.025, 0.7)
#     text(xy_title[1], xy_title[2], labels="M.m. ident.", adj=c(0, 1), xpd=NA)
#     xy_legend <- GetPlotFractionalCoords(0.025, 0.675)
#     # Legend(xy_legend, legend=names(pch_use), pch=pch_use)
#     legend(x=xy_legend[1], y=xy_legend[2], ConvertTtoU(names(pch_use)),
#            pch=pch_use, bty="n", xpd=NA, ncol=2, x.intersp=0.5, y.intersp=0.8)
#   }

#   # Finish the plot.
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }




# # S5D___________________________________________________________________________
# PlotPairingAndOffsetModelAgainstData <- function(
#   mirna, experiment, n_constant=3, sitelist="progthrp_suppcomp",
#   offset_lim=c(-4, 16), len_lim=c(4, 11), pos_lim=c(9, 23), loop=FALSE,
#   intercept=FALSE, log_plus_one=FALSE, supp_base=FALSE, offset_base=FALSE,
#   site_base=NULL, alpha=1, makeglobalmodel=TRUE, mirna_label=TRUE, height=3.5,
#   width=3.5, xpos=20, ypos=20, pdf.plot=FALSE
# ) {
#   if (makeglobalmodel) {
#     model <- SubfunctionCall(FitPairingAndOffsetModelWithError)  
#     model <<- model  
#   }
#   data <- model$data
#   model_sim <- model$values
#   num_global <<- length(model_sim)
#   SubfunctionCall(FigureSaveFile2)
#   xmin <- 0.5
#   xmax <- 1000
#   ymin <- xmin
#   ymax <- xmax
#   BlankPlot(log="xy")
#   AddLogAxis(1, label="Model-estimated Kd fold change")
#   AddLogAxis(2, label="Measured Kd fold change")
#   lens <- as.character(as.integer(as.character(data$pos3p))
#                        - as.integer(as.character(data$pos5p)) + 1)
#   segments(xmin, ymin, x1=xmax, y1=ymax, lty=2)
#   points(10^model_sim, 10^data$logkd, col=kThrPLengthCols[lens], lwd=0, xpd=NA)
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
#   xy <- GetPlotFractionalCoords(0.05, 0.95, log="xy")
#   text(xy[1], xy[2], labels=mirna_txt, adj=0)
#   xy <- GetPlotFractionalCoords(0.05, 0.875, log="xy")
#   AddCorrelationToPlot(xy[1], xy[2], x=model_sim, y=data$logkd, rsquared=TRUE,
#                        adj=0)


#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }

# # S10G__________________________________________________________________________
# PlotProgrammedAndRandomPairingCoefficients <- function(
#   mirna_1, experiment_1, mirna_2, experiment_2, n_constant_1=3,
#   sitelist_1="progthrp_suppcomp", n_constant_2=3,
#   sitelist_2="randthrp_suppcomp", offset_lim=c(-4, 16), len_lim=c(4, 11), 
#   pos_lim=c(9, 23), loop=FALSE, supp_base=FALSE,
#   offset_base=FALSE, site_base=NULL, alpha=1, exponential=FALSE,
#   intercept=FALSE, additive=FALSE, log_plus_one=FALSE, mirna_label=TRUE,
#   height=4, width=4, xpos=20, ypos=20, pdf.plot=FALSE
# ) {
#   model_1 <<- SubfunctionCall(
#     FitPairingAndOffsetModelWithError, mirna=mirna_1,
#     experiment=experiment_1, n_constant=n_constant_1, sitelist=sitelist_1,
#     corrected_kds=TRUE
#   )$pairing
#   model_2 <<- SubfunctionCall(
#     FitPairingAndOffsetModelWithError, mirna=mirna_2,
#     experiment=experiment_2, n_constant=n_constant_2, sitelist=sitelist_2,
#     corrected_kds=FALSE
#   )$pairing

#   # Generate the len values
#   pos5p <- rep(colnames(model_1$MLE), each=nrow(model_1$MLE))
#   pos3p <- rep(rownames(model_1$MLE), ncol(model_1$MLE))
#   lens <- as.character(as.integer(as.character(pos3p))
#                        - as.integer(as.character(pos5p)) + 1)
#   # Make sure the second model is the same size as the first model (only
#   # important for let-7a, because the random library is one nucleotide longer).
#   model_2$MLE <- model_2$MLE[rownames(model_1$MLE), colnames(model_1$MLE)]
#   model_2$LowerCI <- model_2$LowerCI[rownames(model_1$MLE),
#                                      colnames(model_1$MLE)]
#   model_2$UpperCI <- model_2$UpperCI[rownames(model_1$MLE),
#                                      colnames(model_1$MLE)]

#   data <- data.frame(`lens`=lens, `kds_prog`=c(model_1$MLE),
#                      `kds_prog_lCI`=c(model_1$LowerCI),
#                      `kds_prog_uCI`=c(model_1$UpperCI),
#                      `kds_rand`=c(model_2$MLE),
#                      `kds_rand_lCI`=c(model_2$LowerCI),
#                      `kds_rand_uCI`=c(model_2$UpperCI))

#   inds <- which(!is.na(data$kds_prog) & !(is.na(data$kds_rand)))
#   data <- data[inds, ]
#   ######################### Make the plot window ###############################
#   SubfunctionCall(FigureSaveFile2)
#   if(log_plus_one) {
#     xmin <- 0.3
#     xmax <- 1e2
#   } else {
#     xmin <- 1
#     xmax <- 1e3
#   }
#   ymin <- xmin
#   ymax <- xmax
#   par(mar=c(3, 3, 1, 1))
#   BlankPlot(log="xy")
#   ############################# Make the axes ##################################
#   AddLogAxis(1, label="Programmed library")
#   AddLogAxis(2, label="Random library")
#   ############################## Addd x = y line ###############################
#   segments(xmin, ymin, x1=xmax, y1=ymax, lty=2)
#   ####################### Add error bars and points ############################
#   if (log_plus_one) {
#     segments(x0=exp(data$kds_prog), y0=exp(data$kds_rand_lCI),
#              y1=exp(data$kds_rand_uCI), lwd=0.5, xpd=NA)
#     segments(x0=exp(data$kds_prog_lCI), y0=exp(data$kds_rand),
#              x1=exp(data$kds_prog_uCI), lwd=0.5, xpd=NA)
#     points(exp(data$kds_prog), exp(data$kds_rand),
#            col=kThrPLengthCols[as.character(data$lens)], lwd=0, xpd=NA)
#     } else {
#       segments(x0=10^data$kds_prog, y0=10^data$kds_rand_lCI,
#                y1=10^data$kds_rand_uCI, lwd=0.5, xpd=NA)
#       segments(x0=10^data$kds_prog_lCI, y0=10^data$kds_rand,
#                x1=10^data$kds_prog_uCI, lwd=0.5, xpd=NA)
#       points(10^data$kds_prog, 10^data$kds_rand,
#              col=kThrPLengthCols[as.character(data$lens)], lwd=0, xpd=NA)

#     }
#   ############################ Add miRNA label #################################
#   if (mirna_1 == "let-7a-21nt") {
#     mirna_txt <- "let-7a"
#   } else if (mirna_1 == "let-7a_minus1") {
#     mirna_txt <- "let-7a(-1)"
#   } else if (mirna_1 == "let-7a_plus1") {
#     mirna_txt <- "let-7a(+1)"
#   } else if (mirna_1 == "let-7a_miR-155") {
#     mirna_txt <- "let-7a-miR-155"
#   } else if (mirna_1 == "miR-155_let-7a") {
#     mirna_txt <- "miR-155-let-7a"
#   } else {
#     mirna_txt <- mirna_1
#   }
#   xy <- GetPlotFractionalCoords(0.05, 1.05, log="xy")
#   text(xy[1], xy[2], labels=mirna_txt, adj=c(0, 1), xpd=NA)
#   ########################### Add pairing coefficients label ###################
#   xy <- GetPlotFractionalCoords(0.05, 0.95, log="xy")
#   text(xy[1], xy[2], labels="Pairing coefficients", adj=c(0, 1))
#   ####################### Add weighted r squared ###############################
#   xy <- GetPlotFractionalCoords(1, 0.025, log="xy")
#   AddCorrelationToPlot(
#     xpos=xy[1], ypos=xy[2], x=data$kds_prog, y=data$kds_rand,
#     weights=c(1/(data$kds_rand - data$kds_rand_lCI)), rsquared=TRUE, adj=c(1, 0)
#   )
#   ############################ Finish the plot #################################
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }

# # S10H__________________________________________________________________________
# PlotProgrammedAndRandomOffsetCoefficients <- function(
#   mirna_1, experiment_1, mirna_2, experiment_2, n_constant_1=3,
#   sitelist_1="progthrp_suppcomp", n_constant_2=3,
#   sitelist_2="randthrp_suppcomp", offset_lim=c(-4, 16), len_lim=c(4, 11),
#   pos_lim=c(9, 23), loop=FALSE, log_plus_one=FALSE, supp_base=FALSE,
#   offset_base=FALSE, site_base=NULL, alpha=1, globalmodel=TRUE,
#   mirna_label=TRUE, height=4, width=4, xpos=20, ypos=20, pdf.plot=FALSE
# ) {
#   offsets_p <<- SubfunctionCall(
#     FitPairingAndOffsetModelWithError, mirna=mirna_1,
#     experiment=experiment_1, n_constant=n_constant_1, sitelist=sitelist_1,
#     corrected_kds=TRUE
#   )$offsets
#   offsets_r <<- SubfunctionCall(
#     FitPairingAndOffsetModelWithError, mirna=mirna_2,
#     experiment=experiment_2, n_constant=n_constant_2, sitelist=sitelist_2,
#     corrected_kds=FALSE
#   )$offsets

#   SubfunctionCall(FigureSaveFile2)
#   if (log_plus_one) {
#     xmin <- -4
#     xmax <- 0.5
#     label.space <- 0.5
#     tick.space <- 0.5
#   } else {
#     xmin <- 0
#     xmax <- 1.2
#     tick.space <- 0.1
#     label.space <- 0.2
#   }
#   ymin <- xmin
#   ymax <- xmax
#   par(mar=c(3, 3, 1, 1))
#   BlankPlot()
#   # xmax <- 1
#   # ymax <- 1
#   AddLinearAxis(1, label.space=label.space, tick.space=tick.space, label="Programmed library")
#   AddLinearAxis(2, label.space=label.space, tick.space=tick.space, label="Random library")
#   segments(xmin, ymin, x1=xmax, y1=ymax, lty=2)
#   segments(x0=offsets_p[, 1], y0=offsets_r[, 2], y1=offsets_r[, 3], lwd=0.5,
#            xpd=NA)
#   segments(x0=offsets_p[, 2], y0=offsets_r[, 1], x1=offsets_p[, 3], lwd=0.5,
#            xpd=NA)
#   x <- offsets_p[, 1]
#   y <- offsets_r[, 1]
#   ######################### Set up the colors for the points ###################
#   r_col <- -0.75
#   color.dist <- cubeHelix(length(x) + 8,
#                           start=0.25*pi, r=r_col, hue=0.8)[5:(length(x) + 4)]

#   points(offsets_p[, 1], offsets_r[, 1], col=color.dist, lwd=0, xpd=NA)

#   if (mirna_1 == "let-7a-21nt") {
#     mirna_txt <- "let-7a"
#   } else if (mirna_1 == "let-7a_minus1") {
#     mirna_txt <- "let-7a(-1)"
#   } else if (mirna_1 == "let-7a_plus1") {
#     mirna_txt <- "let-7a(+1)"
#   } else if (mirna_1 == "let-7a_miR-155") {
#     mirna_txt <- "let-7a-miR-155"
#   } else if (mirna_1 == "miR-155_let-7a") {
#     mirna_txt <- "miR-155-let-7a"
#   } else {
#     mirna_txt <- mirna_1
#   }
#   xy <- GetPlotFractionalCoords(0.05, 1.05)
#   text(xy[1], xy[2], labels=mirna_txt, adj=c(0, 1), xpd=NA)
#   xy <- GetPlotFractionalCoords(0.05, 0.95)
#   text(xy[1], xy[2], labels="Offset coefficients", adj=c(0, 1), xpd=NA)
#   # xy <- GetPlotFractionalCoords(0.05, 0.875)
#   xy <- GetPlotFractionalCoords(1, 0.025)

#   # AddCorrelationToPlot(xy[1], xy[2], x=x, y=y, rsquared=TRUE,
#   #                      adj=0)
#   # xy <- GetPlotFractionalCoords(0.05, 0.775)
#   AddCorrelationToPlot(xy[1], xy[2], x=x, y=y,
#                        weights=c(1/(y - offsets_r[, 2])), rsquared=TRUE,
#                        adj=c(1, 0))
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }

# # S7____________________________________________________________________________
# PlotSiteBySeedMismatchesAndOffset <- function(
#   mirna, experiment, len, pos5p, n_constant=3, sitelist="progthrp",
#   corrected_kds=TRUE, kd_fc=TRUE, model_values=FALSE, log_plus_one=FALSE,
#   modelnew=FALSE, makeglobalmodel=TRUE, supp_base=FALSE, rand_base=FALSE,
#   show_base=FALSE, offset_lim=c(-4, 16), sumseed=FALSE, show_site_names=TRUE,
#   show_offsets=TRUE, print_mirna=TRUE, highlight_column=FALSE, percentile=0.25,
#   xpos=20, ypos=20, height=2, width=1.51, pdf.plot=FALSE
# ) {
#   if (model_values) {
#     if (modelnew) {
#       model <- SubfunctionCall(FitPairingAndOffsetByMismatchModel)
#       pairing <- model$pairing
#       offsetXmm <- model$offsetXmm
#       pos_5p <- as.character(start_pos)
#       pos_3p <- as.character(start_pos + site_len - 1)
#       dG_p <- model$pairing[pos_3p, pos_5p]
#       mat_out <- dG_p*offsetXmm
#     } else {
#       if (makeglobalmodel) {
#         model <<- SubfunctionCall(FitPairingOffsetAndMismatchModelSingle)
#       }
#       row_name <- sprintf("%s|%s", pos5p, pos5p + len - 1)
#       if (row_name %in% rownames(model$coefs)) {
#         pairing <- model$coefs[sprintf("%s|%s", pos5p, pos5p + len - 1), 1]
#       } else {
#         pairing <- NA
#       }
#       mm <- model$coefs[grep("mer", rownames(model$coefs)), 1]
#       names(mm) <- grep("mer", rownames(model$coefs), value=TRUE)
#       print(pairing)
#       offsets <- model$coefs[grep("(mer|\\||b)", rownames(model$coefs), perl=TRUE, invert=TRUE), 1]
#       names(offsets) <- grep("(mer|\\||b)", rownames(model$coefs), perl=TRUE, invert=TRUE, value=TRUE)
#       # dG_p <- model$pairing[pos_3p, pos_5p]
#       if (log_plus_one) {
#         mm_mat <- mm %o% (offsets* 0 + 1)
#         offset_mat <- (mm*0 + 1) %o% offsets
#         mat_out <- log10(exp(pairing + mm_mat + offset_mat) + 1)
#       } else {
#         mat_out <- pairing*(mm %o% offsets)
#       }
#     }
#     # print("1913")
#     # mat_out <- cbind(mat_out[, 1, drop=FALSE]*0, mat_out)
#     # colnames(mat_out)[1] <- "None"
#     # print(mat_out)
#     R_mat <<- mat_out
#     R_mat_data <- SubfunctionCall(MakeMismatchByOffsetMatrix)
#     R_mat_data <<- R_mat_data
#     R_mat <- R_mat * (R_mat_data*0 + 1)
#   } else {
#     print("here")
#     R_mat <<- SubfunctionCall(MakeMismatchByOffsetMatrix)
#     # R_mat <<- log10(1/mat_out)
#   }
#   if (supp_base) {
#     if (show_offsets) height <- 1.89
#     else              height <- 1.4
#     if (show_site_names)  width <- 2.25
#     else                  width <- 1.4
#   } else {
#     if (!show_offsets)  height <- 1.52
#     if (!show_site_names) width <- 0.97
#   }
#   if (!supp_base) {
#     R_mat <- R_mat[nrow(R_mat):1, ]
#     colnames_keep <- as.character(seq(-4, 16, by=2))
#   } else {
#     colnames_keep <- c(as.character(offset_lim[1]:offset_lim[2]))
#   }
#   R_mat <- R_mat[, colnames_keep]
#   if (!show_base & colnames(R_mat)[1] == "None") {
#     R_mat <- R_mat[, -1]
#   }
#   SubfunctionCall(FigureSaveFile2)
#   xmin <- 0
#   xmax <- ncol(R_mat)
#   ymin <- 0
#   ymax <- nrow(R_mat)
#   if (show_site_names) {
#     if (supp_base) mar2 <- 4
#     else           mar2 <- 2.60
#   } else {
#     mar2 <- 0.25
#   }
#   if (show_offsets) {
#     mar1 <- 2.5
#   } else {
#     mar1 <- 0.42
#   }
#   par(mar=c(mar1, mar2, 1, 0.125))
#   BlankPlot()
#   xlefts <- rep(seq(0, ncol(R_mat) - 1), each=nrow(R_mat))
#   xright <- xlefts + 1
#   ybottom <- rep(seq(ymin, ymax - 1), ncol(R_mat))
#   ytop <- ybottom + 1
#   coords <<- cbind(xlefts, ybottom)
#   # Make the x-axis
#   AddLinearAxis(1, 1, 1, label=NA,
#                 alt_lab=NA,
#                 alt_lab_pos=xmin:(xmax - 1) + 0.5,
#                 alt_tick_pos=TRUE,
#                 blank_lab=TRUE,
#                 line=1.4)

#   if (show_offsets) {
#     label <- "Offsets"
#     if (supp_base) alt_lab <- c(0, 3, 6, 9)
#     else           alt_lab <- c(0, 10)
#     AddLinearAxis(1, 1, 1, label="Offsets",
#                   alt_lab=alt_lab,
#                   # alt_lab_pos=xmin:(xmax - 1) + 0.5,
#                   alt_lab_pos=(xmin:(xmax - 1) + 0.5)[which(as.integer(colnames(R_mat)) %in% alt_lab)],
#                   noline=TRUE,
#                   line=1.4)
#   }
#   if (show_site_names) {
#     if (supp_base) {
#       x_pos <- -0.5
#       text(x=x_pos,
#            y=ymin:(ymax - 1) + 0.5,
#            labels=gsub("(8mer-mm)(.*)", rownames(R_mat), replacement="\\2",
#                        perl=TRUE),
#            xpd=NA,
#            adj=c(1, 0.5))
#     } else {
#       x_pos <- -6.2 + rep(c(0, 1, 2)*1.8, 6)
#       nucs <- gsub("(8mer-mm)([ACTG])([2-7])", rownames(R_mat), replacement="\\2",
#                        perl=TRUE)
#       nucs <- gsub("T", nucs, replacement="U")
#       positions <- 7:2

#       text(x=x_pos,
#            y=ymin:(ymax - 1) + 0.5,
#            labels=nucs,
#            xpd=NA,
#            adj=c(0.5, 0.5))
#       text(x=-0.4,
#            y=(0:5)*3 + 1.5,
#            labels=positions,
#            xpd=NA,
#            adj=c(1, 0.5))
#     }
#   }
#   # # Add the label for the seed if not plotting relative to seed kds.
#   col.inds <- floor((R_mat + 25/33)/(3 + 25/33)*124 + 1)
#   col.inds <- floor((R_mat)/log10(700)*99 + 1)

#   x_lab_pos <- -0.75
#   # Load the seaborn cubeHelix color palette, and conver the indeces from the
#   # and make them bounded between zero and 1
#   start_col <- 0.5
#   r_col <- -0.75
#   # Make distribution of 100 colors and a 101th gray color for NA indeces.
#   color.dist <- c(rev(cubeHelix(100, start=start_col, r=r_col, hue=0.8)),
#                   "gray90")
#   color.dist[100] <- kMaxValueColor

#   col.inds[which(col.inds > 100)] <- 100
#   col.inds[which(col.inds < 1)] <- 1
#   col.inds[which(is.na(col.inds))] <- 101
#   # Make the color rectangle.
#   rect(xlefts, ybottom, xright, ytop, col=color.dist[col.inds], lwd=0.4, 
#        xpd=NA, border="white")
#   if (highlight_column) {
#     # Get the the name of the column that contributed to the Main text
#     # calculation of the random-library seed-coefficient value.
#     cols_best <- grep(sprintf("l%sp%so", len, pos5p), colnames(table_use_global)[1:round(percentile*ncol(table_use_global))],
#                       value=TRUE)
#     # If this is one of the top percentile columns:
#     if (length(cols_best) != 0) {
#       # Get the name of the offset
#       offset <- unlist(strsplit(cols_best, split="o"))[2]
#       x_left <- which(colnames(R_mat) == offset) - 1
#       x_right <- x_left + 1
#       # Make a horizontal line above that column in the figure.
#       segments(x0=x_left, y0=ymax + 0.1, x1=x_right, lwd=1, col="black",
#                lend=1, xpd=NA)
#     }
#   }
#   mtext(sprintf("%s-%s", pos5p, pos5p + len - 1), side=3, at=xmax,
#         adj=1, cex=par("cex"))
#   if (print_mirna) {
#     if (mirna == "miR-7-23nt")        mirna_text <- "miR-7"
#     else if (mirna == "let-7a-21nt")  mirna_text <- "let-7a"
#     else                              mirna_text <- mirna
#     if (supp_base) at_use <- xmin
#     else           at_use <- -4
#     mtext(mirna_text, side=3, at=at_use, adj=0, cex=par("cex"))
#     xy <- GetPlotFractionalCoords(fx=0, fy=1, log="y")
#   }
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
#   return(R_mat)
# } 



# # S6D___________________________________________________________________________
# PlotPairingOffsetAndMismatchModelAgainstData <- function(
#   mirna, experiment, n_constant=3, modelnew=FALSE, offset_lim=c(-4, 16),
#   len_lim=c(4, 11), pos_lim=c(9, 23), makeglobalmodel=TRUE, log_plus_one=FALSE,
#   mirna_label=TRUE, height=3, width=3, xpos=20, ypos=20, pdf.plot=FALSE
# ) {
#   if (modelnew) {
#     model <- SubfunctionCall(FitPairingAndOffsetByMismatchModel)
#     pairing <- model$pairing
#     offsetXmm <- model$offsetXmm
#     print(offsetXmm)
#     data <- model$data
#     model_sim <- apply(data, 1, function(row) {
#       pairing[as.character(as.integer(row[4])), as.character(as.integer(row[3]))]*offsetXmm[row[7], row[6]]
#     })
#   } else {
#     if (makeglobalmodel) {
#       model <<- SubfunctionCall(FitPairingOffsetAndMismatchModelSingle)
#       print("done with model")
#     }
#     data <- model$data
#     model_sim <- model$values
#   }
#   SubfunctionCall(FigureSaveFile2)
#   xmin <- 0.2
#   xmax <- 800
#   ymin <- xmin
#   ymax <- xmax
#   par(mar=c(3, 3, 1, 1))
#   BlankPlot(log="xy")
#   AddLogAxis(1, label="Model-estimated Kd fold change")
#   AddLogAxis(2, label="Measured Kd fold change")
#   lens <- as.character(as.integer(as.character(data$pos3p))
#                        - as.integer(as.character(data$pos5p)) + 1)
#   segments(xmin, ymin, x1=xmax, y1=ymax, lty=2)
#   cols <- sapply(kThrPLengthCols, function(color) {
#     rgb_i <- c(col2rgb(color)/255)
#     rgb(red=rgb_i[1], green=rgb_i[2], blue=rgb_i[3], alpha=0.2)
#   })
#   names(cols) <- names(kThrPLengthCols)
#   points(10^model_sim, 10^data$logkd, col=cols[lens], lwd=0)
#   num_global <<- length(model_sim)
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
#   xy <- GetPlotFractionalCoords(0.05, 0.95, log="xy")
#   text(xy[1], xy[2], labels=mirna_txt, adj=0)
#   xy <- GetPlotFractionalCoords(0.05, 0.875, log="xy")
#   AddCorrelationToPlot(xy[1], xy[2], x=model_sim, y=data$logkd, rsquared=TRUE,
#                        adj=0)
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }


# PlotPairingCorrelationForShiftedLet7 <- function(mirna,
#   offset_lim=c(-4, 16), len_lim=c(4, 11), pos_lim=c(9, 23), corrected_kds=TRUE,
#   log_plus_one=FALSE, lengthnorm=FALSE, height=3.5, width=3.5, xpos=20, ypos=20, pdf.plot=FALSE
# ) {
#   # pairing_l7 <<- SubfunctionCall(FitPairingAndOffsetModelWithError,
#   #                           mirna="let-7a-21nt",
#   #                           experiment="equil_c2_nb")$pairing
#   # pairing_l7_shift <<- SubfunctionCall(FitPairingAndOffsetModelWithError,
#   #                           experiment="equil_c_nb")$pairing
#   if (mirna == "let-7a_minus1") {
#     pairing_l7_use <<- pairing_l7$MLE[c(-1, -2), c(-1, -2)]
#     pairing_l7_shift_use <<- pairing_l7_shift$MLE[c(-1, -nrow(pairing_l7_shift$MLE)), c(-1, -ncol(pairing_l7_shift$MLE))]
#   } else if (mirna == "let-7a_plus1") {
#     pairing_l7_use <- pairing_l7$MLE[c(-1, -2, -nrow(pairing_l7_shift$MLE)), c(-1, -2, -ncol(pairing_l7_shift$MLE))]
#     pairing_l7_shift_use <- pairing_l7_shift$MLE[c(-1, -2, -3), c(-1, -2, -3)]

#   }
#   print(pairing_l7_use)
#   print(pairing_l7_shift_use)
#   if (log_plus_one) {
#     pairing_l7_use <- pairing_l7_use/log(10)
#     pairing_l7_shift_use <- pairing_l7_shift_use/log(10)
#   }
#   if (lengthnorm) {
#     pairing_l7_use <- NormalizeByLength(pairing_l7_use)
#     pairing_l7_shift_use <- NormalizeByLength(pairing_l7_shift_use)
#   }

#   SubfunctionCall(FigureSaveFile2)
#   if (lengthnorm) {
#     xmin <- 0.1
#     xmax <- 10
#   } else {
#     xmin <- 1
#     xmax <- 400      
#   }
#   ymin <- xmin
#   ymax <- xmax
#   BlankPlot(log="xy")
#   stops <- rep(rownames(pairing_l7_use), times=ncol(pairing_l7_use))
#   starts <- rep(colnames(pairing_l7_use), each=nrow(pairing_l7_use))
#   print(head(stops))
#   print(head(starts))
#   len_sites <- as.integer(stops) - as.integer(starts) + 1
#   cols <- kThrPLengthCols
#   cols_use <- cols[as.character(len_sites)]
#   cols_use_global <<- cols_use

#   AddLogAxis(1, label="Kd fold change; native let-7a")
#   AddLogAxis(2, label="Kd fold change; shifted let-7a")
#   x.line <- seq(log10(xmin), log10(xmax), length=20)

#   if (mirna == "let-7a_minus1") {
#     x_lm_G <- c(pairing_l7_use[, c(1, 2)])
#     y_lm_G <- c(pairing_l7_shift_use[, c(1, 2)])
#     # print(pairing_l7_use)
#     # print(pairing_l7_shift_use)
#     x_lm_nonG <- c(pairing_l7_use[, c(-1, -2)])
#     y_lm_nonG <- c(pairing_l7_shift_use[, c(-1, -2)])
#     linmodel_G <- lm(y_lm_G ~ x_lm_G - 1)
#     linmodel_nonG <- lm(y_lm_nonG ~ x_lm_nonG - 1)
#     predict_G <- predict(linmodel_G, data.frame(x_lm_G=x.line), interval = 'confidence')
#     predict_nonG <- predict(linmodel_nonG, data.frame(x_lm_nonG=x.line), interval = 'confidence')
#     polygon(c(10^x.line, 10^rev(x.line)), c(10^predict_G[, 2], rev(10^predict_G[, 3])), col=ConvertRColortoRGB("red", alpha=0.4), border=NA)
#     polygon(c(10^x.line, 10^rev(x.line)), c(10^predict_nonG[, 2], rev(10^predict_nonG[, 3])), col=ConvertRColortoRGB("green", alpha=0.4), border=NA)

#   }

#   x_lm <- c(pairing_l7_use)
#   y_lm <- c(pairing_l7_shift_use)
#   linmodel <- lm(y_lm ~ x_lm - 1)
#   print(x_lm)
#   print(y_lm)
#   predict_all <- predict(linmodel, data.frame(x_lm=x.line), interval = 'confidence')
#   print(predict_all)
#   segments(xmin, ymin, xmax, ymax, lty=2)
#   # polygon(c(10^x.line, 10^rev(x.line)), c(10^predict_all[, 2], rev(10^predict_all[, 3])), col=ConvertRColortoRGB("gray90", alpha=0.4), border=NA)
#   polygon(c(10^x.line, 10^rev(x.line)), c(10^predict_all[, 2], rev(10^predict_all[, 3])), col=ConvertRColortoRGB("gray20", alpha=0.4), border=NA)
#   print(cols_use)
#   # points(c(10^pairing_l7_use), c(10^pairing_l7_shift_use), col=cols_use)
#   points(c(10^pairing_l7_use[, c(-1, -2)]), c(10^pairing_l7_shift_use[, c(-1, -2)]), col="forestgreen")
#   points(c(10^pairing_l7_use[, c(1, 2)]), c(10^pairing_l7_shift_use[, c(1, 2)]), col="red")

#   if (mirna == "let-7a_minus1") {
#     mirna_txt <- "let-7a(-1)"
#   } else if (mirna == "let-7a_plus1") {
#     mirna_txt <- "let-7a(+1)"
#   }
#   xy <- GetPlotFractionalCoords(0.025, 1.1, log="xy")
#   text(xy[1], xy[2], labels=mirna_txt, adj=c(0, 1), xpd=NA)

#   xy <- GetPlotFractionalCoords(1, 0.025, log="xy")
#   AddCorrelationToPlot(c(pairing_l7_use), c(pairing_l7_shift_use), xpos=xy[1],
#                        ypos=xy[2], rsquared=TRUE, adj=c(1, 0))


#   # xy <- GetPlotFractionalCoords(0.025, 1.1)
#   # Legend(xy, legend=c("let-7a(-1)", "let-7a(+1)"), col=c("blue", "red"))
#   # Finish the plot.
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }




# # 6E____________________________________________________________________________
# PlotOffsetCrossCorrelation <- function(
#   offset_lim=c(-4, 16), len_lim=c(4, 11), pos_lim=c(9, 23), corrected_kds=TRUE,
#   log_plus_one=FALSE, height=4, width=4, xpos=20, ypos=20, pdf.plot=FALSE
# ) {
#   off_l7 <- SubfunctionCall(FitPairingAndOffsetModelWithError,
#                             mirna="let-7a-21nt",
#                             experiment="equil_c2_nb")$offsets[, 1]
#   off_l7_m1 <- SubfunctionCall(FitPairingAndOffsetModelWithError,
#                             mirna="let-7a_minus1",
#                             experiment="equil_c_nb")$offsets[, 1]
#   off_l7_p1 <- SubfunctionCall(FitPairingAndOffsetModelWithError,
#                             mirna="let-7a_plus1",
#                             experiment="equil_c_nb")$offsets[, 1]

#   p1_lag <- ccf(off_l7_p1, off_l7, plot=FALSE, lag.max=5)
#   # p1_lag_alt <- ccf(off_l7, off_l7_p1, plot=FALSE, lag.max=5)

#   m1_lag <- ccf(off_l7_m1, off_l7, plot=FALSE, lag.max=5)
#   SubfunctionCall(FigureSaveFile2)
#   par(mar=c(3, 3, 2, 1))
#   xmin <- -5
#   xmax <- 5
#   ymin <- -0.4
#   ymax <- 1
#   BlankPlot()
#   AddLinearAxis(1, 1, 2, label="Offset difference (nt)", alt_lab_pos=c(-4, -2, 0, 2, 4))
#   AddLinearAxis(2, 0.1, 0.2, label="Cross-correlation")
#   points(-5:5, m1_lag$acf, xlim=c(-5, 5), ylim=c(-0.5, 1), type="o", col="blue")
#   points(-5:5, p1_lag$acf, col="red", type="o")
#   xy <- GetPlotFractionalCoords(0.025, 1.2)
#   Legend(xy, legend=c("let-7a(-1)", "let-7a(+1)"), col=c("blue", "red"))
#   # Finish the plot.
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }



# ########### FUNCTIONS THAT MAKE PLOTS NOT SPECIFICALLY INTENDED ################
# ########################### TO BE IN THE PAPER #################################

# # Figure N1 ####################################################################
# PlotPairwiseInputFrequencies <- function(
#   mirna1, experiment1, mirna2, experiment2, ident=FALSE, xpos=20, ypos=20,
#   height=5, width=5, pdf.plot=FALSE
# ) {
#   # First, load each of the kmer files.
#   counts1 <- SubfunctionCall(GetProgrammedKmerFrequencies, mirna=mirna1,
#                              experiment=experiment1, condition="I")
#   counts2 <- SubfunctionCall(GetProgrammedKmerFrequencies, mirna=mirna2,
#                              experiment=experiment2, condition="I")
#   # Determine the names of each of the possible 8mer sites (1), 7mer-m8 sites 
#   # (3), 7mer-A1 sites (3), 6mer sites (9), and 8mer-mm sites (18), for both
#   # miRNAs.
#   sites1 <- c(GetAll8mersCanonicalSite(mirna1), GetAll8merMmSites(mirna1))
#   sites2 <- c(GetAll8mersCanonicalSite(mirna2), GetAll8merMmSites(mirna2))
#   kmers1 <- sapply(sites1, GetSiteSeq, mirna=mirna1)
#   kmers2 <- sapply(sites2, GetSiteSeq, mirna=mirna2)
#   x <- counts1[c(kmers1, setdiff(rownames(counts1), kmers1)), 1]
#   y <- counts2[c(kmers2, setdiff(rownames(counts2), kmers2)), 1]
#   cols_seed <- rep(
#     kSiteColors[c("8mer", "7mer-m8", "7mer-A1", "6mer")],
#     times=c(1, 3, 3, 9)
#   )
#   cols_seed_mm_and_none <- rep(
#     c("green", "orange", "black"),
#     times=c(18, nrow(counts1) - length(kmers1) - 1, 1)
#   )
#   cols_seed_mm_and_none <- sapply(cols_seed_mm_and_none, ConvertRColortoRGB, alpha=0.5)
#   # cols_seed <- sapply(cols_seed_1mm, ConvertRColortoRGB, alpha=0.8)
  
#   cols <- c(cols_seed, cols_seed_mm_and_none)


#   SubfunctionCall(FigureSaveFile2)
#   xmin <- 1e-9
#   xmax <- 1
#   ymin <- xmin
#   ymax <- xmax
#   BlankPlot(log="xy")
#   # Add the axes.
#   AddLogAxis(1, label=sprintf("%s; %s", mirna1, experiment1))
#   AddLogAxis(2, label=sprintf("%s; %s", mirna2, experiment2))
#   # Label the points
#   segments(xmin, ymin, x1=xmax, y1=ymax, lty=2)
#   Points(x, y, col=cols)
#   xy <- GetPlotFractionalCoords(0, 1.1, log='xy')
#   Legend(xy,
#          legend=c("8mer", "7mer-m8", "7mer-A1", "6mer",
#                   "8mer-1mm", "8mer-2mm", "None"),
#          col=unique(cols))
#   # Finish the plot.
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }

# # Figure N2 & N3 ###############################################################
# PlotProgAndRandKds <- function(mirna, experiment, equilibrium_nb=FALSE,
#                                xpos=20, ypos=20, height=5, width=5,
#                                prog_mm=FALSE,  pdf.plot=FALSE) {
#   if (mirna %in% c("let-7a-21nt", "let-7a_plus1", "let-7a_minus1",
#                    "let-7a_miR-155")) {
#     if (equilibrium_nb) {
#       mirna_rand <- "let-7a-21nt"
#       experiment_rand <- "equilibrium_nb"
#     } else {
#       mirna_rand <- "let-7a"
#       experiment_rand <- "equilibrium"
#     }
#   } else if (mirna == "miR-155_let-7a") {
#     mirna_rand <- "miR-155"
#     experiment_rand <- "equilibrium"
#   } else {
#     mirna_rand <- mirna
#     experiment_rand <- "equilibrium"
#   }
#   if (mirna == "miR-1") {
#     buffer_rand <- TRUE
#     combined_rand <- FALSE
#   } else {
#     buffer_rand <- FALSE
#     combined_rand <- TRUE
#   }
#   # Load the two Kd files.
#   kds_rand <- SubfunctionCall(EquilPars, mirna=mirna_rand,
#                               experiment=experiment_rand, n_constant=n_constant,
#                               sitelist="programmed", buffer=buffer_rand,
#                               combined=combined_rand)
#   kds_prog <- SubfunctionCall(EquilPars, n_constant=3,
#                               sitelist="programmed_collapsed")
#   # Get the names of the sites in `kds_rand`, as these are the ones to look for
#   # in the `programmed_collapsed` library.
#   sites_base <- gsub("_Kd", replacement="",
#                      rownames(kds_rand)[1:(nrow(kds_rand) - 3)])
#   # Calculate, for each of the sites represented in the random library, the
#   # log average of its value over the 18 mismatch sites in the
#   # `programmed_collapsed` library.
#   kds_prog_ave <- sapply(sites_base, function(site) {
#     # Get the indeces of all 18 rows with the site.
#     if (prog_mm) {
#       row_inds <- grep(sprintf("^%s_Kd", site),
#                        rownames(kds_prog), perl=TRUE)
#       GeoMean(kds_prog[row_inds, 2])      

#     } else {
#       row_inds <- grep(sprintf("%s\\|8mer-mm[ACTG][2-7]_Kd", site),
#                        rownames(kds_prog), perl=TRUE)
#       GeoMean(kds_prog[row_inds, 2])      
#     }
#   })
#   SubfunctionCall(FigureSaveFile2)
#   xmin <- 1e-4
#   xmax <- 10
#   ymin <- xmin
#   ymax <- xmax
#   BlankPlot(log="xy")
#   # Add the axes.
#   AddLogAxis(1, label="Relative Kd; Random Library")
#   AddLogAxis(2, label="Relative Kd; Programmed Library")

#   # Label the points
#   segments(xmin, ymin, x1=xmax, y1=ymax, lty=2)
#   cols <- c(kSiteColors[sites_base[1:6]], rep("green", 18))
#   cols <- ConvertRColortoRGB(cols, alpha=rep(c(1, 0.5), times=c(6, 18)))
#   cols <<- cols
#   print(kds_rand[1:length(kds_prog_ave), 2])
#   Points(y=kds_prog_ave, x=kds_rand[1:length(kds_prog_ave), 2],
#          col=cols)
#   xy <- GetPlotFractionalCoords(1, 0, log='xy')
#   Legend(xy,
#          legend=c(sites_base[1:6], "8mer-1mm"),
#          col=cols[1:7], xjust=1, yjust=0)

#   xy <- GetPlotFractionalCoords(0.05, 1, log='xy')
#   text(x=xy[1], y=xy[2], labels=mirna, adj=c(0, 1))
#   xy <- GetPlotFractionalCoords(0.05, 0.95, log='xy')
#   text(x=xy[1], y=xy[2], labels=experiment, adj=c(0, 1))
#   if (equilibrium_nb) {
#     xy <- GetPlotFractionalCoords(0.05, 0.90, log='xy')
#     text(x=xy[1], y=xy[2], labels="NB random", adj=c(0, 1))
#   }
#   # Finish the plot.
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }


# PlotKmerCountsTheoryVersusData <- function(
#   mirna, experiment, condition, n_constant=0, sitelist="programmed_suppcomp",
#   outto11mers=TRUE, downto2mers=TRUE, downto3mers=FALSE, downto4mers=FALSE,
#   end3prand=TRUE, collapsemm=FALSE, xpos=20,
#   ypos=20, height=5, width=5, pdf.plot=FALSE
# ) {
#   probs <- sapply(2:11, function(k) {
#     prob_k <- ProbKmerInRandom(k, 25)
#     prob_k_plus_1 <- ProbKmerInRandom(k + 1, 25)
#     prob_3p_k <- 1 - (1 - prob_k)^length(GetMirna3pKmers(mirna, k))
#     prob_not_3p_k_plus_1 <- (1 - prob_k_plus_1)^length(GetMirna3pKmers(mirna, k + 1))
#     return(prob_3p_k * prob_not_3p_k_plus_1)
#    })
#   sXc <- SubfunctionCall(SitesXCounts)
#   sXc_norm <- t(t(sXc)/colSums(sXc))

#   average_kmers <- sapply(2:11, function(k) {
#     grep_string <- sprintf("^%smer-m[0-9]{1,2}\\.[0-9]{1,2}\\|[0-9]{1,2}\\|Comp", k)
#     kmers <- sXc_norm[grep(grep_string, rownames(sXc_norm), perl=TRUE), condition]
#     nonzero_kmers <- kmers[which(kmers != 0)]
#     return(GeoMean(nonzero_kmers))
#   })
#   SubfunctionCall(FigureSaveFile2)
#   xmin <- 1e-7
#   xmax <- 1
#   ymin <- xmin
#   ymax <- xmax
#   BlankPlot(, log="xy")
#   AddLogAxis(1, label="Theoretical library fraction")
#   AddLogAxis(2, label="Observed library fraction")
#   cols <- rainbow(n=10) 
#   points(probs, average_kmers, pch=19, col=cols)

#   xy <- GetPlotFractionalCoords(0.95, 0.05, log="xy")
#   text(xy[1], xy[2], labels=sprintf("%s\n%s", mirna, experiment), adj=c(1, 0))

#   segments(xmin, ymin, x1=xmax, y1=ymax, lty=2, lwd=0.5)
#   xy <- GetPlotFractionalCoords(0.05, 0.95, log="xy")
#   Legend(xy, legend=sprintf("%smer", 2:11), col=cols)

#   # Finish the plot.
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }

# PlotGeoMeanKmerKdsVersusEnrichments <- function(
#   mirna, experiment, n_constant=0, sitelist="programmed_suppcomp",
#   downto2mers=TRUE, outto11mers=TRUE, end3prand=TRUE, collapsemm=FALSE,
#   xpos=20, ypos=20, height=5, width=5, pdf.plot=NULL
# ) {
#   sXc <- SubfunctionCall(SitesXCounts)
#   sXc_norm <- t(t(sXc)/colSums(sXc))
#   kds <- SubfunctionCall(EquilPars)

#   fracts_I <- sapply(2:11, function(k) {
#     grep_string <- sprintf(
#       "^%smer-m[0-9]{1,2}\\.[0-9]{1,2}\\|[0-9]{1,2}\\|Comp", k
#     )
#     kmers <- sXc_norm[grep(grep_string, rownames(sXc_norm), perl=TRUE), 1]
#     nonzero_kmers <- kmers[which(kmers != 0)]
#     return(GeoMean(nonzero_kmers))
#   })

#   kds_geoMean <- sapply(2:11, function(k) {
#     grep_string <- sprintf(
#       "^%smer-m[0-9]{1,2}\\.[0-9]{1,2}\\|[0-9]{1,2}\\|Comp_Kd", k
#     )
#     return(GeoMean(kds[grep(grep_string, rownames(kds), perl=TRUE), 2]))
#   })

#   fracts_A <- sapply(3:8, function(column) {
#     sapply(2:11, function(k) {
#     grep_string <- sprintf(
#       "^%smer-m[0-9]{1,2}\\.[0-9]{1,2}\\|[0-9]{1,2}\\|Comp", k
#     )
#     kmers <- sXc_norm[grep(grep_string, rownames(sXc_norm), perl=TRUE), column]
#     nonzero_kmers <- kmers[which(kmers != 0)]
#     return(GeoMean(nonzero_kmers))
#     })
#   })
#   fracts_A <<- fracts_A
#   fracts_I <<- fracts_I
#   Rs <- fracts_A/fracts_I
#   Rs_norm <- t(t(Rs)/Rs[1, ])
#   Rs_norm <<- Rs_norm
#   SubfunctionCall(FigureSaveFile2)
#   xmin <- 0.5
#   xmax <- 100
#   ymin <- 0.1
#   ymax <- 100
#   BlankPlot(log="xy")
#   AddLogAxis(1, label="Kd fold change")
#   AddLogAxis(2, label="Geometric mean enrichment")

#   print(Rs_norm)
#   kds_geoMean <<- kds_geoMean
#   colors <- c("red", "orangered", "forestgreen", "blue", "purple", "gray")
#   for (i in 1:ncol(Rs_norm)) {
#     print(kds_geoMean[1]/kds_geoMean)
#     print(Rs_norm[, i])
#     lines(kds_geoMean[1]/kds_geoMean, Rs[, i], col=colors[i], pch=1, lwd=1)
#     points(kds_geoMean[1]/kds_geoMean, Rs[, i], col=colors[i], pch=19)
#   }
#   xy <- GetPlotFractionalCoords(0.95, 0.05, log="xy")
#   text(xy[1], xy[2], labels=sprintf("%s\n%s", mirna, experiment), adj=c(1, 0))
#   xy <- GetPlotFractionalCoords(0.05, 0.95, log="xy")
#   Legend(xy, legend=c("40", "12.6", "4", "1.26", "0.4", "0"), col=colors)
#   # Finish the plot.
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }




# ################################################################################
# # FIGURE 7
# ################################################################################

# # 7A-C__________________________________________________________________________
# PlotSiteMismatches <- function(
#   mirna, experiment, start_mm, stop_mm, n_constant=3, bulge=FALSE,
#   model_values=FALSE, mm_and_bulge=TRUE, offset_lim=c(-2, 6), new=TRUE,
#   kd_fc=TRUE, win_average=3, corrected_kds=TRUE, titles=FALSE, deltaG=FALSE,
#   xpos=20, ypos=20, height=1.75, pdf.plot=FALSE
# ) {
#   if (model_values) {
#     parameters <- SubfunctionCall(GetThreePModelMismatchCoefficients)
#     kdfc_mmdb <- parameters[[2]]
#     if (deltaG) {
#       dG_base <- GetThreePSiteDeltaG(mirna, start_mm, stop_mm)
#       dG_mm <- sapply(names(kdfc_mmdb)[-1], GetThreePSiteDeltaG, mirna=mirna,
#                   start=start_mm, stop=stop_mm)
#       dG_all <- c(dG_base, dG_mm)
#       names(dG_all) <- names(kdfc_mmdb)
#       kdfc_mmdb <- dG_all
#     }
#   } else {
#     kdfc_mmdb <<- log10(SubfunctionCall(GetThreePrimeMmBDKds, best_average=TRUE))
#   }
#   # make the new matrix that will be used to visualize the mismatches/bulges,
#   # where the rows are the four possible nucleotides, and the columns are the
#   # positions of the mismatch or bulge.
#   if (bulge) {
#     ncol_mat <- stop_mm - start_mm
#     colnames_mat <- as.character((start_mm + 1):stop_mm)
#   } else {
#     ncol_mat <- stop_mm - start_mm + 1
#     colnames_mat <- as.character(start_mm:stop_mm)
#   }
#   if (bulge) {
#     nrow_mat <- 5
#     rownames_mat <- c("A", "C", "G", "T", "d")
#   } else {
#     nrow_mat <- 4
#     rownames_mat <- c("A", "C", "G", "T")
#   }
#   # Pre-allocate the output values of the matrix.
#   logkdfc_mat <- matrix(NaN, nrow=nrow_mat, ncol=ncol_mat,
#                         dimnames=list(rownames_mat, colnames_mat))
#   for (mm_i in names(kdfc_mmdb)[-1]) {
#     nuc <- substr(mm_i, start=1, stop=1)
#     if (!(nuc %in% c("A", "C", "G", "T"))) {
#       nuc <- "d"
#       par_pos <- 1
#     } else {
#       par_pos <- 2
#     }
#     if (substr(mm_i, par_pos, par_pos) == "(") {
#       splits <- unlist(strsplit(unlist(strsplit(mm_i, split="\\("))[2],
#                        split="\\)"))[1]
#       range <- unlist(strsplit(splits, split="\\."))
#       pos_list <- as.character(seq(as.integer(range[1]), as.integer(range[2])))
#       for (pos in pos_list) {
#         logkdfc_mat[nuc, pos] <- kdfc_mmdb[mm_i]
#       }
#     } else {
#       pos <- substr(mm_i, start=par_pos, stop=nchar(mm_i))
#       logkdfc_mat[nuc, pos] <- kdfc_mmdb[mm_i]
#     }
#   }
#   if (!bulge) {
#     for (i in colnames(logkdfc_mat)) {
#       mir_nuc <- RevComplement(substr(kMirnaSeqs[mirna], start=i, stop=i))
#       logkdfc_mat[mir_nuc, i] <- kdfc_mmdb["wt"]
#     }
#   }
#   rownames(logkdfc_mat)[4] <- "U"
#   if (bulge) {
#     rownames(logkdfc_mat)[5] <- "Del."
#     logkdfc_mat <- logkdfc_mat[c("Del.", "A", "C", "G", "U"), ]
#   }
#   logkdfc_mat_global <<- logkdfc_mat
#   if (mirna %in% c("let-7a-21nt", "miR-155")) base_width <- 12
#   else                                        base_width <- 11    
#   if (bulge)  base_width <- base_width - 1
#   if (titles) mar1 <- 2.2
#   else        mar1 <- 1.7
#   mar2 <- 1.7
#   mar3 <- 0.5
#   mar4 <- 0.2
#   height <- 1.1 + (mar1 + mar3)/5
#   width <- base_width/5*0.8 + (mar2 + mar4)/5
#   SubfunctionCall(FigureSaveFile2)
#   par(mar=c(mar1, mar2, mar3, mar4))
#   xmin <- 0
#   if (bulge) {
#     xmax_lists <- c(`let-7a-21nt`=11, `miR-1`=10, `miR-155`=10, `let-7a`=11, `miR-124`=13, `lsy-6`=13, `miR-7-23nt`=13)
#   } else {
#     xmax_lists <- c(`let-7a-21nt`=12, `miR-1`=11, `miR-155`=11, `let-7a`=12, `miR-124`=14, `lsy-6`=14, `miR-7-23nt`=14)
#     # xmax_lists <- c(`let-7a-21nt`=12, `miR-1`=11, `miR-155`=14, `let-7a`=12, `miR-124`=14, `lsy-6`=14, `miR-7-23nt`=14)
#   }
#   offset_lists <- c(`let-7a-21nt`=20, `miR-1`=21, `miR-155`=23, `let-7a`=20, `miR-124`=22, `lsy-6`=22, `miR-7-23nt`=23)
#   xmax <- xmax_lists[mirna]
#   col_offset <- offset_lists[mirna] - max(as.integer(colnames(logkdfc_mat)))
#   ymin <- 0
#   ymax <- 5 
#   BlankPlot()
#   ymin <- 0
#   if (!bulge) {
#     ymax <- 4
#   }
#   # Define the positions of the rectangles in the matrix.
#   xlefts <- rev(rep(seq(0, ncol(logkdfc_mat) - 1), nrow(logkdfc_mat))) + col_offset
#   if (bulge) {
#     xlefts <- xlefts + rep(c(-0.5, 0), times=c(ncol(logkdfc_mat), ncol(logkdfc_mat)*(nrow(logkdfc_mat) - 1)))
#     x_axis_label_pos <- xmin:(xmin + ncol(logkdfc_mat) - 1) + col_offset    
#   } else {
#     x_axis_label_pos <- xmin:(xmin + ncol(logkdfc_mat) - 1) + 0.5 + col_offset    
#   }
#   xright <- xlefts + 1
#   ybottom <- rep(rev(seq(ymax - 1, ymin)), each=ncol(logkdfc_mat))
#   ytop <- ybottom + 1
#   # Make the x-axis
#   if (titles) {
#     if (bulge)  x_label <- "Bulge/del. position"
#     else        x_label <- "Position"
#   } else {
#     x_label <- ""
#   }
#   # labels_inds <- seq(1, length(x_axis_label_pos), by=3)
#   labels_inds <- seq(1, length(x_axis_label_pos))
#   print(labels_inds)
#   AddLinearAxis(1, 1, 3, label=x_label,
#                 label_pos_ticks=TRUE,
#                 alt_lab=rev(colnames(logkdfc_mat))[labels_inds],
#                 alt_lab_y_dist=0.01,
#                 alt_lab_pos=x_axis_label_pos[labels_inds],
#                 alt_tick_pos=TRUE,
#                 alt_tick_lab_space=TRUE)
#   # Make the y-axis ############################################################
#   # This part was added in allow the period to be a bit to the right of the
#   # nucleotide letters. It shrinks par(mgp)[2] from 0.4 to 0.2, which brings the
#   # axis labels closer to the tick marks by 50%. This allows the period to be
#   # slightly below the axis lines. The gsub adds a space in the string for each
#   # nucleotide letter, so that each of them will also still look justified with
#   # the `Del.` string at the letter `l` rather than the period.
#   par_mgp <- par("mgp")
#   par_mgp[2] <- 0.2
#   par(mgp=par_mgp)
#   yaxis_labs <- gsub("(A|C|G|U)", rownames(logkdfc_mat), replacement="\\1 ", perl=TRUE)
#   AddLinearAxis(2, 1, 1, label="",
#                 alt_lab=yaxis_labs,
#                 alt_lab_pos=ymin:(ymax - 1) + 0.5,
#                 alt_tick_pos=TRUE)
#   if (deltaG) {
#     col.inds <- floor(t(logkdfc_mat)/-22*99 + 1)
#   } else {
#     col.inds <- floor(t(logkdfc_mat)/log10(700)*99 + 1)
#   }
#   start_col <- 0
#   r_col <- 0.4
#   # Part where the rectangles are modified for the bulge/deletion plots, to
#   # extend some of them in the case of the ambiguous-position bulges. (i.e.,
#   # the a T(17.18) bulge).
#   if (bulge) {
#     # Make the list of indeces to exclude.
#     inds_exclude <- c()
#     # Get the positions and nucleotides.
#     positions <- rep(rownames(col.inds), ncol(col.inds))
#     nucs <- rep(colnames(col.inds), each=nrow(col.inds))
#     full_matrix <- do.call(
#       "cbind", 
#       list(xlefts, xright, ybottom, ytop, c(col.inds), positions, nucs)
#     )
#     colnames(full_matrix) <- c("xleft", "xright", "ybottom", "ytop", "col.inds",
#                                "pos", "nuc")
#     # Remove the row associated with the deletion at the 3'-most position,
#     # because it is indistinguishable from a site shorter by nucleotide..
#     inds_del_exclude <- which(full_matrix[, "nuc"] == "Del." &
#                               full_matrix[, "pos"] == stop_mm)
#     full_matrix <- full_matrix[-inds_del_exclude, ]
#     # Iterate over the sites from the original matrix that are of the type
#     # A(13.14).
#     for (name in grep("\\(", names(kdfc_mmdb), value=TRUE, perl=TRUE)) {
#       # Get the nucleotide idenetity, replace T with U.
#       nuc <- gsub("^(.*)\\(.*\\..*\\)$", replacement="\\1", name, perl=TRUE)
#       if (nuc == "") nuc <- "Del."
#       else if (nuc == "T") nuc <- "U"
#       pos_l <- gsub("^.*\\((.*)\\..*\\)$", replacement="\\1", name, perl=TRUE)
#       pos_r <- gsub("^.*\\(.*\\.(.*)\\)$", replacement="\\1", name, perl=TRUE)
#       inds_use <- sapply(seq(as.integer(pos_l), as.integer(pos_r)), function(pos) {
#         which(full_matrix[, "nuc"] == nuc & full_matrix[, "pos"] == pos)
#       })
#       full_matrix[inds_use[1], "xleft"] <- full_matrix[inds_use[length(inds_use)], "xleft"]
#       inds_exclude <- c(inds_exclude, inds_use[-1])
#     }
#     prev_nan <- FALSE
#     for (ind in seq(nrow(full_matrix))) {
#       if (full_matrix[ind, "col.inds"] == "NaN") {
#         xleft_i <- full_matrix[ind, "xleft"]
#         ybottom_i <- full_matrix[ind, "ybottom"]
#         if (prev_nan) {
#           if (ybottom_i == ybottom_save & as.numeric(xleft_i) == xleft_save - 1) {
#             full_matrix[ind_save, "xleft"] <- xleft_i
#             inds_exclude <- c(inds_exclude, ind)
#           } else {
#             ind_save <- ind
#           }
#         } else {
#           prev_nan <- TRUE
#           ind_save <- ind
#         }
#         xleft_save <- as.numeric(xleft_i)
#         ybottom_save <- ybottom_i
#       } else {
#         prev_nan <- FALSE
#       }
#     }
#     # Remove the rows that are now redundant with the widened rectangles.
#     full_matrix <- full_matrix[-inds_exclude, ]
#     # Split the matrix back into the constituent information necessary to plot
#     # the cells of the heatmap.
#     xlefts <- as.numeric(full_matrix[, "xleft"])
#     xright <- as.numeric(full_matrix[, "xright"])
#     ybottom <- as.numeric(full_matrix[, "ybottom"])
#     ytop <- as.numeric(full_matrix[, "ytop"])
#     col.inds <- as.integer(full_matrix[, "col.inds"])
#   }
#   # Make distribution of 100 colors and a 101th gray color for NA indeces.
#   color.dist <- c(rev(cubeHelix(100, start=start_col, r=r_col, hue=0.8)),
#                   "gray90")
#   color.dist[100] <- kMaxValueColor
#   # Make the color index scale.
#   col.inds[which(col.inds > 100)] <- 100
#   # col.inds[which(col.inds > 200)] <- 102
#   col.inds[which(col.inds < 1)] <- 1
#   col.inds[which(is.na(col.inds))] <- 101
#   # Make the color rectangle.
#   rect(xlefts, ybottom, xright, ytop, col=color.dist[col.inds], lwd=0.4, 
#        border="white", xpd=NA)
#   # Add colors around the perfect Watson-Crick pairs.
#   if (!bulge) {
#     xlefts <- rev(seq(0, ncol(logkdfc_mat) - 1)) + col_offset
#     xright <- xlefts + 1
#     ybottom <- sapply(colnames(logkdfc_mat), function(col_name) {
#       nuc <- RevComplement(
#         substr(kMirnaSeqs[mirna], start=col_name, stop=col_name),
#         RNA=TRUE
#       )
#       ybottom <- which(rownames(logkdfc_mat) == nuc) - 1
#     })
#     ytop <- ybottom + 1
#     rect(xlefts, ybottom, xright, ytop, col=NA, lwd=0.6, xpd=NA,
#          border="skyblue1")
#     # Add the range of nucleotides for the site.
#     xy <- GetPlotFractionalCoords(0, 1.05)
#     text(xy[1], xy[2], labels=sprintf("%s-%s", start_mm, stop_mm), adj=c(0, 0),
#          xpd=NA)
#   }
#   if (model_values & titles & !bulge) {
#     xy <- GetPlotFractionalCoords(0.9, 1.05)
#     text(xy[1], xy[2], labels="Model", adj=c(1, 0), xpd=NA)
#   }
#   # Finish the plot.
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }

# # 7A-C__________________________________________________________________________
# PlotSiteMismatchesBar <- function(
#   mirna, experiment, start_mm, stop_mm, n_constant=3, bulge=FALSE,
#   model_values=FALSE, mm_and_bulge=TRUE, offset_lim=c(-4, 16), new=TRUE,
#   kd_fc=TRUE, win_average=3, corrected_kds=TRUE, legend=TRUE,
#   xpos=20, ypos=20, height=2.5, width=3.5, pdf.plot=FALSE
# ) {
#   # Get the matrix with the dG values and the error.
#   dG_df <- SubfunctionCall(GetDeltaDeltaGResidual, obs_only=TRUE, error=TRUE)
#   # Replace the wobble symbols in the data frame
#   dG_df[, "mm_type"] <- gsub(pattern="w", replacement="", x=as.character(dG_df[, "mm_type"]))
#   # Define the plot region.
#   SubfunctionCall(FigureSaveFile2)
#   xmin <- 0
#   xmax <- (stop_mm - start_mm + 1)*4
#   ymin <- 0
#   ymax <- 3
#   par(mar=c(2, 2.8, 0.6, 0.1))
#   BlankPlot()
#   # Assign the x-axis variables 
#   x_l_global <- 0
#   x_l_i <- 0:2
#   x_r_i <- x_l_i + 1
#   for (pos_i in start_mm:stop_mm) {
#     inds_use <- which(dG_df$pos == pos_i)
#     vals <- dG_df$ddG_obs[inds_use]
#     errors <- dG_df$ddG_obs_errors[inds_use]
#     cols <- ConvertRColortoRGB(kNucleotideColorsColorBlind[dG_df$mm_type[inds_use]], alpha=0.7)
#     if (length(vals) != 0) {
#       # Add the bars of the bar-plot, followed by the error bars.
#       rect(xleft=x_l_i + x_l_global, ybottom=0, xright=x_r_i + x_l_global,
#            ytop=vals, col=cols, border=NA, xpd=NA)
#       arrows(x0=x_l_global + (x_l_i + x_r_i)/2, y0=vals - errors,
#              y1=vals + errors, length=0.01, angle=90, lwd=0.5, code=3, xpd=NA)
#     }
#     # Update the procession across the axis.
#     x_l_global <- x_l_global + 4
#   }
#   # Add the y-axis.
#   AddLinearAxis(2, 1, 0.5, label=expression(Delta*Delta*italic(G)))
#   # Add x-axis miRNA positions.
#   text(x=(0:(stop_mm - start_mm))*4 + 1.5, y=-0.25,
#        labels=as.character(start_mm:stop_mm), xpd=NA)
#   text(x=(xmin + xmax)/2, y=-0.6, labels="miRNA position", xpd=NA)
#   if (legend) {
#     ## Add in the legend containing the nucleotides. ###########################
#     xy <- GetPlotFractionalCoords(0.7, 1.1)
#     Legend(xy, col=ConvertRColortoRGB(kNucleotideColorsColorBlind[kNucs],
#                                       alpha=0.7),
#            legend=kNucs, legend_pch_use=15, ncol=2, x.intersp=0.5,
#            text.width=1.5)
#   }
#   # Finish the plot.
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }



# # 7D____________________________________________________________________________
# PlotTerminalMismatchAgainstBulgeDifferenceDist <- function(
#   len, library_type="programmed", n_constant=3, model_values=FALSE,
#   modelweights=TRUE, new=TRUE, kd_fc=TRUE, win_average=3, corrected_kds=TRUE,
#   titles=FALSE, threep=TRUE, bothp=FALSE, limited=FALSE, xpos=20, ypos=20,
#   height=2.5, width=4, pdf.plot=FALSE
# ) {
#   # Load the collapsed data (for the mismatch-only).
#   if (library_type == "programmed") {
#     ratio.df <<- SubfunctionCall(GetAllProgrammedBulgeAndMismatchDifferences)
#   }
#   ratio.df$mirna <- as.character(ratio.df$mirna)
#   ratio.df$end <- as.character(ratio.df$end)
#   ratio.df$mirna[which(ratio.df$mirna == "let-7a-21nt")] <- "let-7a"
#   ratio.df <- cbind(ratio.df, rank=as.character(paste0(ratio.df$mirna, "_", ratio.df$end)),
#                     cols=rep("blue", nrow=ratio.df))
#   ratio.df$cols <- as.character(ratio.df$cols)
#   #### Order the factors such that they correspond to the appropriate order
#   if (library_type == "programmed") {
#     ratio.df$rank <- factor(
#       ratio.df$rank ,
#       levels=c("let-7a_5p", "let-7a_3p", "miR-1_5p", "miR-1_3p", "miR-155_5p",
#                "miR-155_3p")
#     )
#   }
#   if (limited) {
#     inds_keep <- which(ratio.df$mirna == "miR-1" & ratio.df$end == "5p")
#     ratio.df <- ratio.df[inds_keep, ]
#   }
#   ######## Convert the Kd values to be log10 Kd because the bar plot can't be #3
#   ##### done on a log axis #####################################################
#   ratio.df$ratio <- log10(ratio.df$ratio)
#   ###################### Define the axis limits ################################
#   xmin <- log10(0.1)
#   xmax <- log10(4)
#   ymin <- 0.5
#   ymax <- 6.5
#   SubfunctionCall(FigureSaveFile2)

#   par(mar=c(3, 1, 1, 1))
#   xpd=NA
#   BlankPlot(inv="y")
#   # Define the background box showing a Kd fold change of 1. ###################
#   segments(x0=0, y0=ymin, y1=ymax, lty=1, lwd=0.5, col="gray")
#   ################## Add the box plot ##########################################
#   boxplot(ratio ~ rank,
#           data       = ratio.df,
#           axes       = FALSE,
#           xaxt       = "n",
#           col        = "white",
#           horizontal = TRUE,
#           range      = 0,
#           outline    = FALSE,
#           xlim       = c(ymin, ymax),
#           ylim       = c(xmin, xmax),
#           ann        = FALSE,
#           add        = TRUE)
#   AddLogAxis(1, label="Kd fold change attributable", boxplot=TRUE)
#   mtext(side=1, text=expression("to terminal bulges"), line=2.2, cex=par("cex"))

#   par(xpd=NA)
#   ########################### Add the beeswarm points. #########################
#   beeswarm(ratio ~ rank,
#            data        = ratio.df,
#            add         = TRUE,
#            method      = "swarm",
#            corral      = "random",
#            corralWidth = 0.5,
#            pch         = 1,
#            lwd         =1.2,
#            cex         = 0.8,
#            horizontal  = TRUE,  
#            pwcol       = cols,
#            axes        = FALSE,
#            xpd=NA)
#   ############ Plot the Wilcoxon paired test for the ends separateely. #########
#   for (rank_i in unique(ratio.df$rank)) {
#     df_sub <- subset(ratio.df, rank == rank_i, select = c("mm","bu"))
#     test_out <- wilcox.test(log(df_sub$bu), log(df_sub$mm), paired=TRUE,
#                             alternative="greater")
#     print(sprintf("%s: %9.8f", rank_i, test_out$p.value))
#   }
#   ############ Plot the Wilcoxon paired test when combining bot termini for ####
#   ############ each miRNA. #####################################################
#   for (mirna_i in unique(ratio.df$mirna)) {
#     df_sub <- subset(ratio.df, mirna == mirna_i, select = c("mm","bu", "ratio"))
#     test_out <- wilcox.test(log(df_sub$bu), log(df_sub$mm), paired=TRUE,
#       alternative="greater")

#     test_out_alt <- wilcox.test(df_sub$ratio)

#     print(sprintf("%s: %9.8f", mirna_i, test_out$p.value))
#     print(sprintf("%s: %9.8f", mirna_i, test_out_alt$p.value))
#     print(test_out$p.value)
#   }
#   # Add the miRNA names:  
#   xy <- GetPlotFractionalCoords(-0.05, 0.9)
#   text(x=xy[1], y= (0:2)*2 + 1.5, labels=unique(ratio.df$mirna), xpd=NA,
#        adj=c(0, 0.5))
#   # Add the 5′ and 3′ termini labels.
#   xy <- GetPlotFractionalCoords(0.20, 0.9)
#   text(x=xy[1], y=1:6, labels=unique(ratio.df$end), xpd=NA, adj=c(0, 0.5))
#   xy <- GetPlotFractionalCoords(0.18, 0.9)
#   segments(x0=xy[1], y0=c(0.5, 2.5, 4.5) + 0.2, y1=c(2.5, 4.5, 6.5) - 0.2)
  
#   xy <- GetPlotFractionalCoords(1, 1.1, inv="y")
#   # if (library_type == "programmed") {
#   #   library_label <- "Programmed libraries"
#   # } else {
#   #   library_label <- "Random libraries"
#   # }
#   text(xy[1], xy[2], labels=sprintf("%s-nt sites", len), adj=c(1, 1))
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }

# # 7E____________________________________________________________________________
# PlotAllAverageMismatchesAndDeltaGs <- function(
#   library_type="programmed", len, n_constant=3, bulge=FALSE, model_values=FALSE,
#   new=TRUE, kd_fc=TRUE, win_average=3, corrected_kds=TRUE, titles=FALSE,
#   ratio=FALSE, frac_dG=FALSE, frac_ddG=FALSE, xpos=20, ypos=20, height=3.5,
#   width=3.5, pdf.plot=FALSE
# ) {
#   # Load the collapsed data (for the mismatch-only).
#   df_obs_pred <<- SubfunctionCall(MakeAllMirnaDeltaDeltaGMismatchDf)
#   print(df_obs_pred)
#   colors <- c(`A`="#0077BB", `C`="#EE3377", `G`="#EE7733", `U`="#009988")
#   colors <- kNucleotideColorsColorBlind
#   pch_use <- c(`A`=1, `C`=2, `G`=3, `U`=4)
#   mir_nuc <- substr(df_obs_pred$Group.1, start=2, stop=2)
#   tar_nuc <- substr(df_obs_pred$Group.1, start=4, stop=4)
#   # Make the plot and define the limits.
#   SubfunctionCall(FigureSaveFile2)
#   par(mar=c(3, 3, 1, 0.8))
#   if (frac_dG) {
#     xmin <- 0
#     xmax <- 1
#   } else {
#     xmin <- -1
#     xmax <- 5
#   }
#   ymin <- xmin
#   ymax <- xmax
#   BlankPlot()
#   if (frac_dG) {
#     # Add x-axis._______________________________________________________________
#     AddLinearAxis(1, 0.1, 0.2, label="Predicted fractional", percent=TRUE)
#     mtext(side=1, text=expression("reduction in"~Delta*italic(G)~"(%)"),
#           line=2.2, cex=par("cex"))
#     # Add y-axis._______________________________________________________________
#     AddLinearAxis(2, 0.1, 0.2, label="Observed fractional", percent=TRUE, line=1.45)
#     mtext(side=2, text=expression("reduction in"~Delta*italic(G)~"(%)"), line=1.15, cex=par("cex"), las=0)
#   } else {
#     # Add x-axis._______________________________________________________________
#     AddLinearAxis(1, 1, 1, label=expression("Predicted"~Delta*Delta*italic(G)))
#     # Add y-axis._______________________________________________________________
#     AddLinearAxis(2, 1, 1, label=expression("Observed"~Delta*Delta*italic(G)))
#     segments(x0=xmin, y0=0, x1=xmax, lty=2)
#     if (frac_dG) seg_y_max <- ymax
#     else         seg_y_max <- 2.5
#     segments(x0=0, y0=ymin, y1=seg_y_max, lty=2)
#   }
#   if (frac_dG) {
#     x <- 1 - df_obs_pred$ddG_pred
#     y <- 1 - df_obs_pred$ddG_obs
#   } else {
#     x <- df_obs_pred$ddG_pred
#     y <- df_obs_pred$ddG_obs
#   }
#   points(x, y, col=colors[mir_nuc], pch=pch_use[tar_nuc], cex=1.5, lwd=1.25,
#          xpd=NA)

#   x_global <<- df_obs_pred$ddG_pred
#   y_global <<- df_obs_pred$ddG_obs

#   if (!frac_dG) {
#     # Make the first legend.____________________________________________________
#     xy <- GetPlotFractionalCoords(0.06, 1)
#     text(xy[1], xy[2], labels="miRNA", adj=c(0, 0.5), xpd=NA)
#     legend_text_1 <- unique(mir_nuc)
#     xy <- GetPlotFractionalCoords(-0.02, 1)
#     legend(x=xy[1], y=xy[2], kNucs, col=colors[kNucs], lwd=2,
#            bty="n", xpd=NA, ncol=2, seg.len=0.6, x.intersp=0.3, y.intersp=0.8,
#            text.width=0.3)
#     # Make the second legend.___________________________________________________
#     xy <- GetPlotFractionalCoords(0.06, 0.77)
#     text(xy[1], xy[2], labels="Target", adj=c(0, 0.5), xpd=NA)
#     xy <- GetPlotFractionalCoords(0, 0.77)
#     legend(x=xy[1], y=xy[2], kNucs, pch=pch_use[kNucs], bty="n", 
#            xpd=NA, ncol=2, x.intersp=0.6, y.intersp=0.8, pt.lwd=1.25, pt.cex=3*par()$cex, text.width=0.4)
#   }
#   # Print the frac_dG (presumably used for troubleshooting).____________________
#   if (frac_dG) {
#     print(df_obs_pred$ddG_pred/df_obs_pred$ddG_obs)
#   }
#   # Add label of the length of the 3-prime sites used in the plot.______________
#   xy <- GetPlotFractionalCoords(1, 1)
#   text(xy[1], xy[2], labels=sprintf("%s-nt sites", len), adj=c(1, 0.5), xpd=NA)  
#   # Finish the plot.____________________________________________________________
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }

# # 7F____________________________________________________________________________
# PlotAllPositionalMismatches <- function(
#   mirna, len, n_constant=3, model_values=FALSE, rand_data=FALSE, bulge=FALSE,
#   legend=FALSE, new=TRUE, kd_fc=TRUE, win_average=3, offset_lim=c(-4, 16),
#   alt_y_max=FALSE, corrected_kds=TRUE, titles=FALSE, x_axis=FALSE,
#   frac_dG=FALSE, plot_kd_fc=FALSE, xpos=20, ypos=20, height=1, width=4.5,
#   pdf.plot=FALSE
# ) {
#   # Load the collapsed data (for the mismatch-only).
#   mirnas <- c("let-7a-21nt", "miR-1", "miR-155", "let-7a_plus1",
#            "let-7a_minus1", "let-7a_miR-155", "miR-155_let-7a")
#   experiments <- c("equil_c2_nb", "equil_c_nb", "equil_sc_nb", "equil_c_nb",
#                 "equil_c_nb", "equil_c_nb", "equil_c_nb")
#   names(experiments) <- mirnas
#   df_use <<- SubfunctionCall(MakePositionalMismatchDf,
#                              experiment=experiments[mirna])
#   if (x_axis) {
#     height <- 1.35 # max
#     height <- 1.3375 # min
#     height <- 1.34
#   }
#   df_use$nuc <- gsub("U", replacement="T", as.character(df_use$nuc))
#   SubfunctionCall(FigureSaveFile2)
#   if (x_axis) {
#     # par(mar=c(2, 4, 0.6, 0.1))
#     par(mar=c(2, 3, 0.6, 1))
#   } else {
#     par(mar=c(0.5, 3, 0.6, 1))
#   }
#   xmin <- 0
#   xmax <- 15.5*4 # (14 = 23 - 9 + 1)
#   ymin <- 0
#   if (class(alt_y_max) %in% c("integer", "numeric")) {
#     ymax <- alt_y_max
#   } else if (mirna == "miR-1") {
#     ymax <- 1
#   } else {
#     ymax <- 2
#   }
#   BlankPlot()
#   # Convert the strings for the bulge dataframe into strings compatible with the
#   # mismatch based formatting.
#   if (bulge) {
#     df_use$pos <- gsub("^(.*)(A|C|G|U)d$", replacement="\\1",
#                        as.character(df_use$pos))
#     df_use$pos <- gsub("^(.*)(A|C|G|U)t$", replacement="\\1\\.5",
#                        as.character(df_use$pos))
#     df_use$nuc <- ConvertUtoT(gsub("^.*((A|C|G|U))(t|d)$", replacement="\\1",
#                               as.character(df_use$pos_nuc)))
#     inds_change <- grep("\\.5", df_use$pos, perl=TRUE)
#     df_use$pos[inds_change] <- as.character(as.numeric(df_use$pos[inds_change]) - 1)
#     df_use$nuc[inds_change] <- "Bulge"
#   }
#   # Define the left-hand sides of rectangles.
#   x_lefts <- (as.numeric(as.character(df_use$pos)) - 9)*4
#   if (bulge) {
#     x_lefts <- x_lefts + 0.5
#     x_rights <- x_lefts + 2   
#     x_lefts <- x_lefts + 0.2
#     x_rights <- x_rights - 0.2 
#   } else {
#     x_rights <- x_lefts + 3
#   }
#   inds <- order(as.numeric(as.character(df_use$pos)))
#   cols <- ConvertRColortoRGB(kNucleotideColorsColorBlind, alpha=0.7)[df_use$nuc]
#   cols[which(is.na(cols))] <- "#CCCCCC"
#   vals <- df_use$dG_obs
#   rect(xleft=x_lefts, ybottom=0, xright=x_rights,
#        ytop=vals, col=cols, border=NA, xpd=NA)

#   # title(ylab="Observed", line=2.5)
#   AddLinearAxis(2, 1, 0.5, label=expression(Delta*Delta*italic(G)))

#   if (x_axis) {
#     text(x=(0:14)*4 + 1.5, y=-0.30*ymax/2.5, labels=as.character(9:23), xpd=NA)
#     text(x=(xmin + xmax)/2, y=-0.9*ymax/2.5, labels="miRNA position", xpd=NA)
#   }
#   ## Add the miRNA name.
#   if (mirna == "let-7a-21nt") mirna_txt <- "let-7a"
#   else                        mirna_txt <- mirna
#   xy <- GetPlotFractionalCoords(0.025, 1.12)
#   text(x=xy[1], y=xy[2], labels=mirna_txt, xpd=NA, adj=c(0, 1))

#   ## Add the length of the sites. ##############################################
#   if (legend) {
#     if (bulge) {
#       xy <- GetPlotFractionalCoords(0.95, 1.12)
#       labels_use <- "Del."
#       legend_cols <- c(ConvertRColortoRGB(kNucleotideColorsColorBlind[kNucs[1:2]], alpha=0.7),
#                        "#CCCCCC",
#                        ConvertRColortoRGB(kNucleotideColorsColorBlind[kNucs[3:4]], alpha=0.7))
#       legend_names <- c(kNucs[1:2], "Bulge", kNucs[3:4])
#     } else {
#       xy <- GetPlotFractionalCoords(1, 1.12)
#       labels_use <- sprintf("%s-nt sites", len)
#       labels_use <- "Mismatch"
#       legend_cols <- ConvertRColortoRGB(kNucleotideColorsColorBlind[kNucs],
#                                         alpha=0.7)
#       legend_names <- kNucs
#     }
#     text(x=xy[1], y=xy[2], labels=labels_use, xpd=NA, adj=c(1, 1))
#     ## Add in the legend containing the nucleotides. #############################
#     xy <- GetPlotFractionalCoords(0.79, 1.05)
#     Legend(xy, col=legend_cols, legend=legend_names, legend_pch_use=15, ncol=2,
#            x.intersp=0.5, text.width=1.5)
#   }
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }

# # 7G____________________________________________________________________________
# PlotAllAverageSeedMismatchesAndDeltaGs <- function(n_constant=3, A1=FALSE,
#   m8_sites=FALSE, frac_dG=FALSE, xpos=20, ypos=20, height=3.5, width=3.5,
#   pdf.plot=FALSE
# ) {
#   # Load the collapsed data (for the mismatch-only).
#   df_obs_pred <<- SubfunctionCall(GetAllRandomSeedSiteInternalDeltaG)
#   df_obs_pred$Group.1 <- ConvertTtoU(as.character(df_obs_pred$Group.1))
#   colors <- kNucleotideColorsColorBlind
#   pch_use <- c(`A`=1, `C`=2, `G`=3, `U`=4)
#   mir_nuc <- substr(df_obs_pred$Group.1, start=2, stop=2)
#   tar_nuc <- substr(df_obs_pred$Group.1, start=4, stop=4)
#   # Make the plot and define the limits.
#   SubfunctionCall(FigureSaveFile2)
#   par(mar=c(3, 3, 1, 0.8))
#   if (frac_dG) {
#     xmin <- 0
#     xmax <- 1
#   } else {
#     xmin <- -1
#     xmax <- 8
#   }
#   ymin <- xmin
#   ymax <- xmax
#   BlankPlot()
#   if (frac_dG) {
#     AddLinearAxis(1, 0.1, 0.2, label="Predicted fractional", percent=TRUE)
#     mtext(side=1, text=expression("reduction in"~Delta*italic(G)~"(%)"),
#           line=2.2, cex=par("cex"))
#     AddLinearAxis(2, 0.1, 0.2, label="Observed fractional", percent=TRUE,
#                   line=1.45)
#     mtext(side=2, text=expression("reduction in"~Delta*italic(G)~"(%)"),
#           line=1.15, cex=par("cex"), las=0)

#   } else {
#     AddLinearAxis(1, 1, 1, label=expression("Predicted"~Delta*Delta*italic(G)))
#     AddLinearAxis(1, 1, 1, label="",
#                   blank_lab=FALSE, noline=TRUE, alt_lab_pos=c(0, 2, 4, 6, 8))

#     AddLinearAxis(2, 1, 1, label=expression("Observed"~Delta*Delta*italic(G)))
#     segments(x0=xmin, y0=0, x1=xmax, lty=2)
#     segments(x0=0, y0=ymin, y1=ymax, lty=2)
#   }
#   if (frac_dG) {
#     x <- df_obs_pred$frac_ddG_pred
#     y <- df_obs_pred$frac_ddG_obs
#   } else {
#     x <- df_obs_pred$ddG_pred
#     y <- df_obs_pred$ddG_obs
#   }
#   points(x, y, col=colors[mir_nuc], pch=pch_use[tar_nuc], cex=1.5, lwd=1.25,
#          xpd=NA)

#   # Assign global values.
#   x_global <<- x
#   y_global <<- y
#   # Print ratio (presumably while troubleshooting).
#   # if (frac_dG) {
#   #   print(df_obs_pred$ddG_pred/df_obs_pred$ddG_obs)
#   # }
#   # Add the "Seed sites" label to the plot.
#   xy <- GetPlotFractionalCoords(1, 1)
#   text(xy[1], xy[2], labels="Seed sites", adj=c(1, 0.5), xpd=NA)  
#   # Finish the plot.
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }


# CompareSiteMismatches <- function(
#   mirna, experiment, start_mm_1, stop_mm_1, start_mm_2, stop_mm_2, n_constant=3,
#   bulge=FALSE, model_values=FALSE, mm_and_bulge=TRUE, offset_lim=c(0, 10),
#   new=TRUE, kd_fc=TRUE, win_average=3, corrected_kds=TRUE, xpos=20, ypos=20,
#   titles=FALSE, ident=FALSE, mirna_label=TRUE, height=3, width=3, pdf.plot=FALSE
# ) {
#   if (model_values) {
#     pairing_1 <- SubfunctionCall(GetThreePModelMismatchCoefficients,
#                                  start_mm=start_mm_1, stop_mm=stop_mm_1)$pairing
#     pairing_2 <- SubfunctionCall(GetThreePModelMismatchCoefficients,
#                                  start_mm=start_mm_2, stop_mm=stop_mm_2)$pairing
#   } else {
#     pairing_1 <- log10(SubfunctionCall(GetThreePrimeMmBDKds,
#                                        start_mm=start_mm_1, stop_mm=stop_mm_1,
#                                        best_average=TRUE))
#     pairing_2 <- log10(SubfunctionCall(GetThreePrimeMmBDKds,
#                                        start_mm=start_mm_2, stop_mm=stop_mm_2,
#                                        best_average=TRUE))
#   }
#   del_pairing_1 <- pairing_1["wt"] - pairing_1[2:length(pairing_1)]
#   del_pairing_2 <- pairing_2["wt"] - pairing_2[2:length(pairing_2)]
#   inds <- intersect(names(del_pairing_1), names(del_pairing_2))
#   SubfunctionCall(FigureSaveFile2)
#   # par(mar=c(mar1, mar2, mar3, mar4))
#   xmin <- -0.5
#   xmax <- 3
#   ymin <- xmin
#   ymax <- xmax
#   par(mar=c(3, 3, 1, 1))
#   BlankPlot()
#   # Define the positions of the rectangles in the matrix.
#   pos <- as.integer(gsub("^([ACGT])(.*)$", inds, replacement="\\2"))
#   inds_1 <- which(pos == start_mm_1 | pos == stop_mm_1)
#   inds_2 <- which(pos == start_mm_2 | pos == stop_mm_2)
#   inds_1_2 <- intersect(inds_1, inds_2)

#   cols <- kThrPPositionCols[as.character(pos)]
#   cols[inds_1] <- "blue"
#   cols[inds_2] <- "red"
#   cols[inds_1_2] <- "purple"
#   cols[which(is.na(pos))] <- "black"

#   pch_use <- rep(1, length(cols))
#   pch_use[inds_1] <- 19
#   pch_use[inds_2] <- 19
#   pch_use[which(is.na(pos))] <- 19

#   x <- del_pairing_1[inds]
#   y <- del_pairing_2[inds]
#   segments(x0=0, y0=ymin, y1=0.85*(ymax - ymin) + ymin, lty=3)
#   segments(x0=xmin, y0=0, x1=xmax, lty=3)
#   points(x, y, col=cols, pch=pch_use)
#   site_name_1 <- sprintf("%smer-m%s.%s", stop_mm_1 - start_mm_1 + 1, start_mm_1,
#                          stop_mm_1)
#   site_name_2 <- sprintf("%smer-m%s.%s", stop_mm_2 - start_mm_2 + 1, start_mm_2,
#                          stop_mm_2)

#   AddLinearAxis(1, 1, 0.5, label=site_name_1)
#   AddLinearAxis(2, 1, 0.5, label=site_name_2)
#   xy <- GetPlotFractionalCoords(0.95, 0.03)
#   AddCorrelationToPlot(x=x, y=y, xpos=xy[1], ypos=xy[2], rsquared=TRUE,
#                        adj=c(1, 0))
#   if (ident)  identify(x, y, labels=inds)

#   if (mirna_label) {
#     if (mirna == "let-7a-21nt") {
#       mirna_txt <- "let-7a"
#     } else if (mirna == "let-7a_minus1") {
#       mirna_txt <- "let-7a(-1)"
#     } else if (mirna == "let-7a_plus1") {
#       mirna_txt <- "let-7a(+1)"
#     } else if (mirna == "let-7a_miR-155") {
#       mirna_txt <- "let-7a-miR-155"
#     } else if (mirna == "miR-155_let-7a") {
#       mirna_txt <- "miR-155-let-7a"
#     } else {
#       mirna_txt <- mirna
#     }
#     xy <- GetPlotFractionalCoords(0.05, 1)
#     text(xy[1], xy[2], label=mirna_txt, adj=c(0, 0), xpd=NA)
#   }
#   xy <- GetPlotFractionalCoords(1, 1)
#   if (model_values) analysis_label <- sprintf("Model %s-%s", offset_lim[1], offset_lim[2])
#   else              analysis_label <- sprintf("Best %s-nt", win_average)
#   text(xy[1], xy[2], label=analysis_label, adj=c(1, 0), xpd=NA)

#   xy <- GetPlotFractionalCoords(0, 1.07)
#   legend_cols_use <- unique(cols)
#   labels_legend <- unique(pos)
#   pch_legend <- rep(1, length(labels_legend))
#   inds_1 <- which(legend_cols_use %in% c("blue", "red", "purple"))
#   pch_legend[inds_1] <- 19
#   legend(xy[1], xy[2], unique(pos), col=unique(cols), text.width=0.27,
#          ncol=ceiling(length(unique(pos))/2), pch=pch_legend, bty="n", xpd=NA,
#          y.intersp=0.75, x.intersp=0.5, xjust=0)

#   # Finish the plot.
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }

# # S12
# PlotMismatchKdsAgainstDeltaG <- function(
#   mirna, experiment, start_mm, stop_mm, n_constant=3, bulge=FALSE,
#   model_values=FALSE, new=TRUE, kd_fc=TRUE, win_average=3, corrected_kds=TRUE,
#   mirna_residual=FALSE, mirna_res_len=11, just_complement=TRUE, titles=FALSE,
#   xpos=20, ypos=20, height=3.5, width=3.5, pdf.plot=FALSE
# ) {
#   if (model_values) {
#     parameters <- SubfunctionCall(GetThreePModelMismatchCoefficients)
#     kdfc_mmdb <- parameters[[2]]
#   } else {
#     kdfc_mmdb <- log10(SubfunctionCall(GetThreePrimeMmBDKds, best_average=TRUE))
#   }
#   # # kdfc_mmdb_global <<- kdfc_mmdb
#   # kdfc_mmdb <- kdfc_mmdb_global
#   R <- 1.987e-3 # in kcal K-1 mol-1
#   T <- 310.15 # in K
#   dG_base <- GetThreePSiteDeltaG(mirna, start_mm, stop_mm,
#                                  just_complement=just_complement)
#   dG_mm <- sapply(names(kdfc_mmdb)[-1], GetThreePSiteDeltaG, mirna=mirna,
#                   start=start_mm, stop=stop_mm, just_complement=just_complement)
#   dG_all <- c(dG_base, dG_mm)
#   dG_all <- dG_all
#   if (class(mirna_residual) == "character") {
#     exps <- c(`let-7a-21nt`="equil_c2_nb", `miR-1`="equil_c_nb",
#               `miR-155`="equil_sc_nb")
#     residual_matrix <- SubfunctionCall(
#       MakeAllDeltaDeltaGMismatchDf, mirna=mirna_residual,
#       experiment=exps[mirna_residual], len=mirna_res_len
#     )
#   }
#   SubfunctionCall(FigureSaveFile2)
#   if (mirna == "miR-155") {
#     xmin <- -26
#     ymax <- 1e4 
#   } else {
#     xmin <- -18
#     ymax <- 1e3
#   }
#   xmax <- 0
#   ymin <- 1
#   BlankPlot(inv="x", log="y")
#   ############## Get the nucleotide and position of each mismatch #############
#   target_nucs <- sapply(names(dG_all)[-1], substr, start=1, stop=1)
#   pos <- sapply(names(dG_all)[-1], function(name_i) {
#     substr(name_i, start=2, stop=nchar(name_i))
#   })
#   pch_use <- c(`A`=3, `C`=4, `G`=0, `T`=2)
#   miR_seq_list <- unlist(strsplit(substr(kMirnaSeqs[mirna], start_mm, stop_mm),
#                          split=""))
#   names(miR_seq_list) <- start_mm:stop_mm
#   miR_nucs <- miR_seq_list[pos]
#   pairing_str <- sprintf("m%st%s", miR_nucs, gsub("T", target_nucs, replacement="U"))


#   if (class(mirna_residual) == "character") {
#     corrections <- residual_matrix[, 2]
#     names(corrections) <- residual_matrix[, 1]
#     corrections_global <<- corrections
#     dG_all_pre_global <<- dG_all
#     dG_all[-1] <- dG_all[-1] + corrections[pairing_str]
#     dG_all_post_global <<- dG_all
#   }


#   wobbleGU_inds <- which(miR_nucs == "U" & target_nucs == "G")
#   wobbleUG_inds <- which(miR_nucs == "G" & target_nucs == "T")

#   pch_all <- pch_use[target_nucs]
#   pch_all[wobbleGU_inds] <- 15
#   pch_all[wobbleUG_inds] <- 17
#   ############################## Add the points ################################
#   points(dG_all, 10^kdfc_mmdb, pch=c(19, pch_all),
#          col=c("black", kThrPPositionCols[pos]))
#   ############################### Add the axes #################################
#   AddLinearAxis(1, 1, 5, label=expression(Delta*Delta*italic(G)~"(kcal/mol)"))
#   AddLogAxis(2, label="Kd fold-change")
#   ################# Plot the trendline for the predicted ∆∆G ###################
#   xmax_convert <- -log(ymax)*R*T
#   xmin_convert <- -log(ymin)*R*T
#   segments(xmin_convert, ymin, xmax_convert, ymax, lty=2)
#   xmax_convert <- -log(ymax)*R*T + dG_all[1] + kdfc_mmdb["wt"]*log(10)*R*T
#   xmin_convert <- -log(ymin)*R*T + dG_all[1] + kdfc_mmdb["wt"]*log(10)*R*T
#   segments(xmin_convert, ymin, xmax_convert, ymax, lty=4)
#   xmax_convert <- -log(ymax)*R*T + dG_all[1]
#   xmin_convert <- -log(ymin)*R*T + dG_all[1]
#   segments(xmin_convert, ymin, xmax_convert, ymax, lty=3)

#   ########################### Add r-squared to plot ############################
#   xy <- GetPlotFractionalCoords(0.05, 0.05, log="y", inv="x")
#   AddCorrelationToPlot(dG_all[-1], kdfc_mmdb[-1], xy[1], xy[2], rsquared=TRUE)


#   ######################### Add the nucleotide legend ##########################
#   xy <- GetPlotFractionalCoords(0.85, 0.65, log="y", inv="x")
#   legend(xy[1], xy[2], legend=c("A", "C", "G", "U", "wG", "wU"), col="black",
#          pch=c(3, 4, 1, 2, 15, 17), bty="n", ncol=1, text.width=1, xpd=NA)
#   ###################### Add the pairing position legend #######################
#   xy <- GetPlotFractionalCoords(0, 1.25, log="y", inv="x")
#   Legend(xy, legend=c("wt", unique(pos)), ncol=6, text.width=1.25,
#          col=c("black", kThrPPositionCols[unique(pos)]))
#   ########################### Add the miRNA label ##############################
#   xy <- GetPlotFractionalCoords(0, 1, log="y", inv="x")
#   if (mirna == "let-7a-21nt") {
#     mirna_text <- "let-7a"
#   } else {
#     mirna_text <- mirna
#   }
#   text(xy[1], xy[2], label=mirna_text, adj=c(0, 1))
#   # Finish the plot. ###########################################################
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }


# PlotAverageMismatches <- function(
#   mirna, experiment, len, n_constant=3, bulge=FALSE, model_values=FALSE,
#   new=TRUE, kd_fc=TRUE, win_average=3, corrected_kds=TRUE, titles=FALSE,
#   ratio=FALSE, frac_dG=FALSE, frac_ddG=FALSE,
#   xpos=20, ypos=20, height=2.5, width=2.5, pdf.plot=FALSE
# ) {
#   # Load the collapsed data (for the mismatch-only).
#   df_mm_pred <- SubfunctionCall(MakeAllDeltaDeltaGMismatchDf, plot_pred=TRUE)
#   df_mm_emp <- SubfunctionCall(MakeAllDeltaDeltaGMismatchDf, plot_emp=TRUE)
#   print(df_mm_pred)
#   print(df_mm_emp)

#   colors <- c(`A`="blue", `C`="purple", `G`="red", `U`="forestgreen")
#   pch_use <- c(`A`=1, `C`=2, `G`=3, `U`=4)
#   mir_nuc <- substr(df_mm_pred$Group.1, start=2, stop=2)
#   tar_nuc <- substr(df_mm_pred$Group.1, start=4, stop=4)
#   # Make the plot and define the limits.
#   SubfunctionCall(FigureSaveFile2)
#   par(mar=c(3, 3, 1, 1))
#   if (frac_dG | frac_ddG) {
#     xmin <- 0
#     ymin <- 0
#     xmax <- 1
#     ymax <- 1
#   } else {
#     xmin <- -0.5
#     ymin <- -0.5
#     xmax <- 5
#     ymax <- 5
#   }
#   BlankPlot()
#   if (frac_dG) {
#     AddLinearAxis(1, 0.2, 0.1, label=expression("Pred. %"~Delta*italic(G)))
#     AddLinearAxis(2, 0.2, 0.1, label=expression("Obs. %"~Delta*italic(G)))
#   } else if (frac_ddG) {
#     AddLinearAxis(1, 0.2, 0.1, label=expression("Pred. %"~Delta*Delta*italic(G)))
#     AddLinearAxis(2, 0.2, 0.1, label=expression("Obs. %"~Delta*Delta*italic(G)))
#   } else {
#     AddLinearAxis(1, 1, 0.5, label=expression("Predicted"~Delta*Delta*italic(G)))
#     AddLinearAxis(2, 1, 0.5, label=expression("Observed"~Delta*Delta*italic(G)))
#     segments(x0=xmin, y0=0, x1=xmax, lty=2)
#     segments(x0=0, y0=ymin, y1=4.5, lty=2)
#   }
#   points(df_mm_pred$x, df_mm_emp$x, col=colors[mir_nuc], pch=pch_use[tar_nuc],
#          xpd=NA)
#   legend_text <- sprintf("%s:%s", mir_nuc, tar_nuc)
#   if (frac_dG) {
#     xy <- GetPlotFractionalCoords(0, 0.97)
#     legend(x=xy[1], y=xy[2], legend_text, col=colors[mir_nuc],
#            pch=pch_use[tar_nuc], bty="n", xpd=NA, cex=0.55, ncol=1, x.intersp=0.5)
#   } else if (frac_ddG) {
#     xy <- GetPlotFractionalCoords(0.9, 0.95)
#     legend(x=xy[1], y=xy[2], legend_text, col=colors[mir_nuc],
#            pch=pch_use[tar_nuc], bty="n", xpd=NA, cex=0.55, ncol=1, x.intersp=0.5)
#   } else {
#     xy <- GetPlotFractionalCoords(0.1, 0.9)
#     legend(x=xy[1], y=xy[2], legend_text, col=colors[mir_nuc],
#          pch=pch_use[tar_nuc], bty="n", xpd=NA, cex=0.6, ncol=3, x.intersp=0.6)
#   }

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
#   } else if (mirna == "miR-7-23nt") {
#     mirna_txt <- "miR-7"
#   } else {
#     mirna_txt <- mirna
#   }
#   xy <- GetPlotFractionalCoords(0.025, 1.1)
#   text(xy[1], xy[2], labels=sprintf(mirna_txt), adj=c(0, 1), xpd=NA, cex=0.8)  

#   xy <- GetPlotFractionalCoords(0.025, 1)
#   text(xy[1], xy[2], labels=sprintf("%s-nt sites", len), adj=c(0, 1), xpd=NA, cex=0.8)  

#   xy <- GetPlotFractionalCoords(1, 1.1)
#   if (experiment %in% c("equilibrium", "equilibrium2_nb")) {
#     experiment_label <- "Rand. lib."
#   } else {
#     experiment_label <- "Prog. lib."
#   }
#   text(xy[1], xy[2], labels=experiment_label, adj=c(1, 1), xpd=NA, cex=0.8)  

#   xy <- GetPlotFractionalCoords(1, 1)
#   if (model_values) analysis_label <- "Model"
#   else              analysis_label <- sprintf("Best %s", win_average)
#   text(xy[1], xy[2], labels=analysis_label, adj=c(1, 1), xpd=NA, cex=0.8)  




#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }


# PlotThrPMismatchDeltaGResidual <- function(
#   mirna, experiment, len, n_constant=3, bulge=FALSE, model_values=FALSE,
#   new=TRUE, kd_fc=TRUE, win_average=3, corrected_kds=TRUE, titles=FALSE,
#   plot_pred=FALSE, plot_emp=FALSE, ratio=FALSE, frac_dG=FALSE, frac_ddG=FALSE,
#   xpos=20, ypos=20, height=3, width=2.5, pdf.plot=FALSE
# ) {
#   # Load the collapsed data (for the mismatch-only).
#   df_mm <- SubfunctionCall(MakeAllDeltaDeltaGMismatchDf)
#   # Make the plot and define the limits.
#   SubfunctionCall(FigureSaveFile2)
#   par(mar=c(3, 3.5, 2, 0))
#   xmin <- 0
#   xmax <- nrow(df_mm) + nrow(df_mm)/3 - 1
#   if (ratio) {
#     if (frac_ddG) {
#       ymin <- 0
#       ymax <- 20
#     } else if (frac_dG) {
#       ymin <- 0
#       ymax <- 2
#     } else {
#       ymin <- 0
#       ymax <- 1
#     }
#   } else if (frac_ddG | frac_dG) {
#     if (plot_pred | plot_emp) {
#       ymin <- -1
#       ymax <- 6    
#     } else {
#       ymin <- -1
#       ymax <- 1    
#     }
#   } else {
#     if (plot_pred | plot_emp) {
#       ymin <- -1
#       ymax <- 6    
#     } else {
#       ymin <- -5
#       ymax <- 1.5    
#     }
#   }
#   BlankPlot()
#   col.inds <- round((df_mm[, 2] + max(abs(ymin), abs(ymax)))/(2*max(abs(ymin), abs(ymax)))*99 + 1)
#   print(col.inds)
#   col.inds[which(col.inds > 100)] <- 101  
#   col.inds[which(col.inds <= 0)] <- 102  
  
#   # Make distribution of 100 colors and a 101th gray color for NA indeces.
#   color.dist <- c(rev(colorspace::diverging_hcl(100, h=c(180, 50), c=80,
#                                               l=c(20, 95), power=c(0.7, 1.3))), "red", "blue")
#   xleft <- 0:(nrow(df_mm) - 1) + rep(0:(nrow(df_mm)/3 - 1), each=3)
#   xright <- xleft + 1
#   print(xleft)
#   print(xright)
#   # Get the left- and right-hand positions of the rectangles.
#   print(mean(c(ymin, ymax)))
#   print(df_mm[, 2])

#   inds_max <- which(df_mm[, 2] > ymax)
#   df_mm[inds_max, 2] <- ymax
#   inds_min <- which(df_mm[, 2] < ymin)
#   df_mm[inds_min, 2] <- ymin

#   rect(xleft=xleft, ybottom=0, xright=xright, ytop=df_mm[, 2],
#        col=color.dist[col.inds], border=NA)
#   # arrows(x0=(xleft + xright)/2, y0=data[, 2], x1=(xleft + xright)/2,
#   #        y1=data[, 3], length=0.02, angle=90, code=3, xpd=NA)
#   if (ratio) {
#     AddLinearAxis(2, 0.1, 0.2, label=expression(Delta*Delta*italic(G)~"ratio"))
#   } else {
#     AddLinearAxis(2, 0.5, 1, label=expression(Delta*Delta*italic(G)~"residual"))
#   }

#   # # Get the each of the nucleotide mismatches.
#   tar_nucs <- gsub("^m(.*)t(.*)$", replacement="\\2", df_mm$Group.1, perl=TRUE)
#   mir_nucs <- gsub("^m(.*)t(.*)$", replacement="\\1", df_mm$Group.1, perl=TRUE)

#   text(x=-0.1, y = ymin - (ymax - ymin)*0.1, labels="Target:", xpd=NA, adj=c(1, 0.5), xpd=NA)
#   text(x=-0.1, y = ymin - (ymax - ymin)*0.2, labels="miRNA:", xpd=NA, adj=c(1, 0.5), xpd=NA)

#   text(x=xleft + 0.5, y = ymin - (ymax - ymin)*0.1, labels=tar_nucs, xpd=NA)
#   text(x=xleft + 0.5, y = ymin - (ymax - ymin)*0.2, labels=mir_nucs, xpd=NA)

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
#   } else if (mirna == "miR-7-23nt") {
#     mirna_txt <- "miR-7"
#   } else {
#     mirna_txt <- mirna
#   }
#   xy <- GetPlotFractionalCoords(0.025, 1.2)
#   text(xy[1], xy[2], labels=sprintf(mirna_txt), adj=c(0, 1), xpd=NA)  

#   xy <- GetPlotFractionalCoords(0.025, 1.1)
#   text(xy[1], xy[2], labels=sprintf("%s-nt sites", len), adj=c(0, 1), xpd=NA)  

#   if (plot_pred) {
#     xy <- GetPlotFractionalCoords(1, 1.2)
#     text(xy[1], xy[2], labels=expression(Delta*Delta*italic(G)), adj=c(1, 1), xpd=NA)  
#   } else {
#     xy <- GetPlotFractionalCoords(1, 1.2)
#     if (experiment %in% c("equilibrium", "equilibrium2_nb")) {
#       experiment_label <- "Rand. lib."
#     } else {
#       experiment_label <- "Prog. lib."
#     }
#     text(xy[1], xy[2], labels=experiment_label, adj=c(1, 1), xpd=NA)  

#     xy <- GetPlotFractionalCoords(1, 1.1)
#     if (model_values) analysis_label <- "Model"
#     else              analysis_label <- sprintf("Best %s", win_average)
#     text(xy[1], xy[2], labels=analysis_label, adj=c(1, 1), xpd=NA)  
#   }





#   # ## Add label explaining that these are the model mismatch coefficients.
#   # xy <- GetPlotFractionalCoords(0.025, 1.025)
#   # if (empirical) {
#   #   AddLinearAxis(2, 0.1, 0.2, label=expression("Average difference in"~Delta*Delta*italic(G)))
#   #   text(xy[1], xy[2], labels="Empirical variation", adj=c(0, 1), xpd=NA)
#   #   xy <- GetPlotFractionalCoords(0.025, 0.95)
#   #   text(xy[1], xy[2], labels=sprintf("Top %sth percentile", round(percentile*100, 0)), adj=c(0, 1), xpd=NA)

#   # } else {
#   #   AddLinearAxis(2, 0.1, 0.2, label=expression(Delta*Delta*italic(G)~"scaling-factor"))
#   #   text(xy[1], xy[2], labels="Mismatch coefficients", adj=c(0, 1), xpd=NA)  
#   # }
#   # Finish the plot.
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }






# PlotTerminalMismatchAgainstBulge <- function(
#   mirna, experiment, len, n_constant=3, model_values=FALSE,
#   new=TRUE, kd_fc=TRUE, win_average=3, corrected_kds=TRUE, titles=FALSE,
#   threep=TRUE, bothp=FALSE, xpos=20, ypos=20, height=4, width=4, pdf.plot=FALSE
# ) {
#   # Load the collapsed data (for the mismatch-only).
#   if (bothp) {
#     mm_and_bu_5p <- SubfunctionCall(GetBulgeAndMismatchValues, threep=TRUE)
#     mm_and_bu_3p <- SubfunctionCall(GetBulgeAndMismatchValues, threep=FALSE)
#     mm <- c(mm_and_bu_5p$mm, mm_and_bu_3p$mm)
#     bu <- c(mm_and_bu_5p$bu, mm_and_bu_3p$bu)
#   } else {
#     mm_and_bu <- SubfunctionCall(GetBulgeAndMismatchValues)
#     mm <- mm_and_bu$mm
#     bu <- mm_and_bu$bu
#   }
#   print(mm)
#   print(bu)
#   # Make the plot and define the limits.
#   SubfunctionCall(FigureSaveFile2)
#   par(mar=c(3, 3, 1, 1))
#   xmin <- 0.01
#   xmax <- 3
#   ymin <- 0
#   ymax <- 1


#   BlankPlot(log="x")
#   x <- 10^seq(log10(xmin), log10(xmax), length.out=100)
#   y <- ecdf(mm)(x)
#   y2 <- ecdf(bu)(x)
#   lines(x, y, lwd=1, xpd=NA)
#   lines(x, y2, col="blue", lwd=1, xpd=NA)

#   AddLogAxis(1, label="Kd fold change")
#   AddLinearAxis(2, 0.2, 0.1, label="ECDF")
#   xy <- GetPlotFractionalCoords(0.05, 1, log="x")
#   text(xy[1], xy[2], labels=mirna, adj=c(0, 1), xpd=NA)
#   xy <- GetPlotFractionalCoords(0.05, 0.9, log="x")
#   if (bothp) {
#     analysis_label <- "Both 5' and 3' termini"
#   } else if (threep) {
#     analysis_label <- "3' termini only"
#   } else {
#     analysis_label <- "5' termini only "
#   }
#   text(xy[1], xy[2], labels=analysis_label, adj=c(0, 1), xpd=NA)
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }


# PlotBeckerEtAlMismatchData <- function(
#   mirna, xpos=20, ypos=20, height=5, width=5, pdf.plot=FALSE
# ) {
#   data <- LoadBeckerEtAlData(mirna, with_names=TRUE)
#   sites_8mer <- GetAll8merMmSites(mirna)
#   print(sites_8mer)
#   kds.df <- do.call("rbind", lapply(sites_8mer, function(site_8mer) {
#     base_inds <- grep(sprintf("^%s$", site_8mer), data[, 1], perl=TRUE)
#     data[base_inds, ]
#   }))
#   kds.df[["Kd_median (pM)"]] <- log10(as.numeric(kds.df[["Kd_median (pM)"]]))
#   kds.df <- cbind(kds.df, cols=rep("violet", nrow(kds.df)))
#   colnames(kds.df)[3] <- "Kd"
#   colnames(kds.df)[1] <- "site_name"
#   kds.df$site_name <- factor(kds.df$site_name, levels=intersect(sites_8mer, unique(kds.df$site_name)))
#   print(dim(kds.df))
#   kds.df <<- kds.df
#   xmin <- log10(10000)
#   xmax <- log10(1)
#   ymin <- 0.5
#   ymax <- 16.5
#   SubfunctionCall(FigureSaveFile2)

#   par(mar=c(3, 1, 1, 1))
#   xpd=NA
#   BlankPlot(inv="y")
#   # Define the background box showing a Kd fold change of 1. ###################
#   # segments(x0=0, y0=ymin, y1=ymax, lty=1, lwd=0.5, col="gray")
#   segments(x0=xmin, y0=0:16, x1=xmax, y1=0:16, lty=1, lwd=0.5, col="gray")

#   ################## Add the box plot ##########################################
#   boxplot(Kd ~ site_name,
#           data       = kds.df,
#           axes       = FALSE,
#           xaxt       = "n",
#           col        = "white",
#           horizontal = TRUE,
#           range      = 0,
#           outline    = FALSE,
#           xlim       = c(ymin, ymax),
#           ylim       = c(xmin, xmax),
#           ann        = FALSE,
#           add        = TRUE)
#   AddLogAxis(1, label="Kd fold change of terminal wobbles", boxplot=TRUE)
#   par(xpd=NA)
#   ########################### Add the beeswarm points. #########################
#   beeswarm(Kd ~ site_name,
#            data        = kds.df,
#            add         = TRUE,
#            method      = "swarm",
#            corral      = "random",
#            corralWidth = 0.5,
#            pch         = 1,
#            lwd         =1.2,
#            cex         = 0.8,
#            horizontal  = TRUE,  
#            pwcol       = cols,
#            axes        = FALSE,
#            xpd=NA)


#   xy <- GetPlotFractionalCoords(1, 0.9)
#   text(x=xy[1], y= 1:16, labels=unique(kds.df$site_name), xpd=NA,
#        adj=c(1, 0.5))


#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }

# }

# PlotBeckerMismatchAgainstProgrammedData <- function(
#   mirna, use_rand=FALSE, xpos=20, ypos=20, height=3.5, width=3.5, pdf.plot=FALSE
# ) {
#   if (use_rand) {
#     kds_prog <- EquilPars(mirna="let-7a", experiment="equilibrium", n_constant=3, sitelist="randthrp")
#   } else {
#     kds_prog <- ApplyKdCorrection(mirna="let-7a-21nt", experiment="equil_c2_nb", prog_sitelist="progthrp")
#   }
#   kds_prog <<- kds_prog

#   data <- LoadBeckerEtAlData(mirna, with_names=TRUE)
#   all_sites_becker <- as.character(unique(data[, 2]))
#   all_sites_becker_compare <<- sprintf("%s_Kd", all_sites_becker)
#   sites_use <<- intersect(rownames(kds_prog), all_sites_becker_compare)
#   print(length(sites_use))

#   print(sites_use)
#   sites_use <- gsub("_Kd", replacement="", sites_use)
#   print(head(all_sites_becker_compare))
#   sites_8mer <- c("8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-m8", "6mer-A1", GetAll8merMmSites(mirna))
#   sites_8mer <- sites_use
#   print(sites_8mer)
#   kds.df <- do.call("rbind", lapply(sites_8mer, function(site_8mer) {
#     print(site_8mer)
#     base_inds <- grep(sprintf("^%s$", site_8mer), data[, 1], perl=TRUE)
#     data[base_inds, ]
#   }))
#   kds.df[["Kd_median (pM)"]] <- log10(as.numeric(kds.df[["Kd_median (pM)"]]))
#   # kds.df[["Kd_median (pM)"]] <- as.numeric(kds.df[["Kd_median (pM)"]])

#   kds.df <- cbind(kds.df, cols=rep("violet", nrow(kds.df)))
#   colnames(kds.df)[3] <- "Kd"
#   colnames(kds.df)[1] <- "site_name"
#   kds.df$site_name <- factor(kds.df$site_name, levels=sites_use)
#   average_kd <- aggregate(kds.df$Kd, list(kds.df$site_name), mean, na.rm=TRUE)
#   average_kd <<- average_kd
#   # average_kd[, 2] <- -average_kd[, 2]

#   kds_prog_use <- kds_prog[sprintf("%s_Kd", as.character(average_kd[, 1])), 2, drop=FALSE]
#   kds_prog_use <<- kds_prog_use
#   xmin <- 1e-3
#   xmax <- 1e3
#   ymin <- 1e-4
#   ymax <- 1
#   SubfunctionCall(FigureSaveFile2)

#   # par(mar=c(3, 3, 1, 1))
#   BlankPlot(log="xy", inv="xy")
#   points(10^average_kd[, 2]/1000, kds_prog_use[, 1], xpd=NA)
#   AddLogAxis(1, label="Becker et al.")
#   AddLogAxis(2, label="Our study")

#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }

# }


# PlotBeckerPairingMatrix <- function(
#   mirna, offset_, seed_site_, get_global_data=FALSE, offset_lim=c(-4, 16),
#   len_lim=c(4, 11), no_mm=TRUE, pos3p_lim=c(9, 12), kd_fc=TRUE,
#   remove_multi=TRUE, remove_CCC=FALSE, F_method=FALSE, exponential=FALSE,
#   intercept=FALSE, key=FALSE, xlabels=TRUE, suppress_label=FALSE,
#   label_offset=FALSE, mirna_label=FALSE, extralabel=FALSE, sitelabel=FALSE,
#   height=2.8, width=2.5, xpos=20, ypos=20, pdf.plot=FALSE
# ) {
#   # Load the data matrix.
#   print("in pairing range matrix function")
#   R_mat <- t(SubfunctionCall(MakeBeckerPairingMatrix))
#   mar1 <- 2.2
#   mar2 <- 0.5
#   mar3 <- 1.8
#   if(key) {
#     width <- 4
#     mar4 <- 9
#   } else {
#     mar4 <- 2.5
#   }
#   # This is the conversion factor in order to have the correct height of the
#   # plot for the specified width, that makes each box of equal height and width.
#   if (mirna == "let-7a")  mirna_use <- "let-7a-21nt"
#   else                    mirna_use <- mirna
#   mir_length <- nchar(kMirnaSeqs[mirna_use])
#   SubfunctionCall(FigureSaveFile2)
#   par(mar=c(mar1, mar2, mar3, mar4))
#   xmin <- 0
#   xmax <- nrow(R_mat)
#   ymin <- 0
#   ymax <- ncol(R_mat) 
#   BlankPlot()
#   # Assign the positions of the corners of the boxes.
#   xlefts <- rep(seq(0, ncol(R_mat) - 1), nrow(R_mat))
#   xright <- xlefts + 1
#   ybottom <- rep(rev(seq(ymax - 1, ymin)), each=ncol(R_mat))
#   ytop <- ybottom + 1

#   # Make every other label on the y-axis blank.
#   y_axis_pos <- seq(ymin, ymax - 1) + 0.5
#   y_axis_labs <- rownames(R_mat)
#   for (i in 1:length(y_axis_labs)) {
#     if (i %% 2 == 0) {
#       y_axis_labs[i] <- ""
#     }
#   }

#   # Make the x-axis.
#   AddLinearAxis(1, 1, 3, label="3'-paired nt",
#                 label_pos_ticks=TRUE,
#                 alt_lab=colnames(R_mat),
#                 alt_lab_y_dist=0.01,
#                 alt_lab_pos=xmin:(xmax - 1) + 0.5,
#                 alt_tick_pos=TRUE,
#                 alt_tick_lab_space=TRUE)
#   # Make the y-axis.
#   AddLinearAxis(4, 1, 2, label="5'-paired nt",
#                 alt_lab=y_axis_labs,
#                 alt_lab_pos=y_axis_pos,
#                 alt_tick_pos=TRUE,
#                 line=0.8)

#   # Add the label for the seed if not plotting relative to seed kds.
#   col.inds <- floor((R_mat)/log10(700)*99 + 1)
#   x_lab_pos <- -0.75
#   pos_5p <- as.integer(rep(rownames(R_mat), each=ncol(R_mat)))
#   pos_3p <- as.integer(rep(colnames(R_mat), nrow(R_mat)))
#   impossible_cols <- which(pos_5p > pos_3p)

#   start_col <- 0
#   r_col <- 0.4
#   # Make distribution of 100 colors and a 101th gray color for NA indeces.
#   color.dist <- c(rev(cubeHelix(100, start=start_col, r=r_col, hue=0.8)),
#                   "gray90", "white")
#   # color.dist[100] <- "red"
#   # Make the color index scale.
#   col.inds <- sapply(t(col.inds), function(col.ind) {
#     min(max(1, col.ind), 100)
#   })
#   col.inds[which(is.na(col.inds))] <- 101
#   col.inds[impossible_cols] <- 102
#   # Make the color rectangle.
#   rect(xlefts, ybottom, xright, ytop, col=color.dist[col.inds], lwd=0.4, 
#        xpd=NA, border="white")
#   # Make the key.
#   if (key) {
#     mir_scale <- (mir_length - 11)/(21 - 11)
#     num_boxes <- 100
#     y_div <- (ymax - ymin)/num_boxes
#     kpl <- xmax + 5*mir_scale # Left-hand position of the key
#     kw <- mir_scale                # Width of the key
#     rect(xleft=kpl,
#          ybottom=seq(num_boxes)*y_div - y_div,
#          xright=kpl + kw,
#          ytop=seq(num_boxes)*y_div,
#          col=color.dist[1:100], xpd=NA, border=NA)
#     # Generate the axis for the legend and the label
#     bottom <- ymin + y_div/2
#     top <- ymax - y_div/2
#     # if (kdrel) {
#     labels <- c(1, 3, 10, 30, 100, 300, 700)
#     pos_labels <- log(labels)
#     centered_labels <- pos_labels - pos_labels[1]
#     norm_labels <- centered_labels/(centered_labels[length(centered_labels)])
#     height_span <- top - bottom
#     pos_labels <- norm_labels*height_span + y_div/2
#     axis(4, at=pos_labels, labels=labels, lwd=par()$lwd, pos=kpl + kw, xpd=NA)
#     text(x=kpl + kw + 3*mir_scale, y=ymax/2,
#          labels=bquote(italic(K)[D]*.(" fold change")),
#          srt=270, xpd=NA)
#   }  
#   # # Add the label indicating how many nucleotides of pairing.
#   # if (model_values) {
#   #   if (label_offset) {
#   #     mtext(text=sprintf("%s-nt offset; model", offset), side=3, at=xmax, adj=1,
#   #           cex=par("cex"))    
#   #   } else {
#   #     mtext(text="Model-estimated values", side=3,
#   #           at=xmin, adj=0, cex=par("cex"))      
#   #   }
#   # } else if (!(suppress_label)) {
#   if (!(suppress_label)) {
#     mtext(text=sprintf("%s-nt offset", offset_), side=3,
#           at=xmax, adj=1, cex=par("cex"))    
#   }
#   # If the `mirna_label` conditional is true, add the label saying which miRNA
#   # is being looked at.
#   if (mirna_label) {
#     if (mirna == "let-7a-21nt") {
#       mirna_txt <- "let-7a"
#     } else if (mirna == "let-7a_minus1") {
#       mirna_txt <- "let-7a(-1)"
#     } else if (mirna == "let-7a_plus1") {
#       mirna_txt <- "let-7a(+1)"
#     } else if (mirna == "let-7a_miR-155") {
#       mirna_txt <- "let-7a-miR-155"
#     } else if (mirna == "miR-155_let-7a") {
#       mirna_txt <- "miR-155-let-7a"
#     } else {
#       mirna_txt <- mirna
#     }
#     x <- GetPlotFractionalCoords(0.125, 0.5)[1]
#     mtext(text=mirna_txt, side=3, line=0.8, at=x, adj=0, cex=par("cex"))
#   }
#   # Boolean label giving the name of the site, in the top right-hand corner,
#   # next to the miRNA.
#   if (sitelabel) {
#     x <- GetPlotFractionalCoords(1.2, 0.5)[1]
#     seed_site_ <- gsub("mm", replacement="x", x=seed_site_)
#     mtext(text=seed_site_, side=3, line=0.8, at=x, adj=1, cex=par("cex"))
#   }
#   if (extralabel) {
#     str.end3prand <- "r"
#     if (collapsemm) {
#       str.collapsemm <- "sum"
#     } else {
#       str.collapsemm <- "geo"
#     }
#     mtext(text=sprintf("%s_%s", str.end3prand, str.collapsemm), side=3,
#           line=0.7, at=xmax, adj=1, cex=par("cex"))    
#   }

#   # Finish the plot.
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
#   return(R_mat)
# }

# PlotBeckerOffsets <- function(
#   mirna, site_type="compensatory", model_values=FALSE, make_global_df=TRUE,
#   remove_multi=TRUE, remove_CCC=FALSE, offset_lim=c(-4, 16), len_lim=c(4, 11),
#   pos3p_lim=c(9, 12), kd_fc=TRUE, F_method=FALSE, exponential=FALSE,
#   intercept=FALSE, key=FALSE, xlabels=TRUE, suppress_label=FALSE,
#   label_offset=FALSE, mirna_label=FALSE, site_label=TRUE, extralabel=FALSE,
#   sitelabel=FALSE, height=2.8, width=2.5, xpos=20, ypos=20, pdf.plot=FALSE
# ) {
#   # Load the data matrix.

#   df.out <- SubfunctionCall(MakeBeckerThreePDataFrame)
#   print(head(df.out))
#   if (site_type == "compensatory") {
#     num_offsets <- sapply(unique(df.out$offset), function(offset_i) {
#       nrow(subset(df.out, offset == offset_i & mm %in% GetAll8merMmSites(mirna)))
#     })
#   } else {
#     num_offsets <- sapply(unique(df.out$offset), function(offset_i) {
#       nrow(subset(df.out, offset == offset_i & mm == site_type))
#     })
#   }

#   names(num_offsets) <- unique(df.out$offset)
#   num_offsets <- num_offsets[order(as.integer(names(num_offsets)))]
#   print(num_offsets)
#   df_use <- data.frame(`offset`=as.integer(names(num_offsets)),
#                        `n`=as.integer(num_offsets))
#   print(df_use)
#   SubfunctionCall(FigureSaveFile2)
#   par(mar=c(3, 3, 1, 1))
#   xmin <- 0
#   xmax <- nrow(df_use)
#   ymin <- 0
#   if (site_type == "compensatory") {
#     ymax <- 600
#     ytick <- 50
#     ylabel <- 100
#   } else {
#     ymax <- 100
#     ytick <- 10
#     ylabel <- 20
#   }
#   BlankPlot()
#   # Assign the positions of the corners of the boxes.
#   xlefts <- xmin:(xmax - 1)
#   xright <- xlefts + 1
#   ybottom <- rep(0, nrow(df_use))
#   ytop <- df_use$n
#   # Make the x-axis
#   AddLinearAxis(1, 1, 1, label="5'-paired nt",
#                 alt_lab=df_use$offset,
#                 alt_lab_pos=xmin:(xmax - 1) + 0.5,
#                 alt_tick_pos=TRUE)
#   AddLinearAxis(2, ytick, ylabel, label="#")

#   rect(xlefts, ybottom, xright, ytop, col="black", lwd=0.4, 
#        xpd=NA, border="white")
 
#   if (mirna_label) {
#     if (mirna == "let-7a-21nt") {
#       mirna_txt <- "let-7a"
#     } else if (mirna == "let-7a_minus1") {
#       mirna_txt <- "let-7a(-1)"
#     } else if (mirna == "let-7a_plus1") {
#       mirna_txt <- "let-7a(+1)"
#     } else if (mirna == "let-7a_miR-155") {
#       mirna_txt <- "let-7a-miR-155"
#     } else if (mirna == "miR-155_let-7a") {
#       mirna_txt <- "miR-155-let-7a"
#     } else {
#       mirna_txt <- mirna
#     }
#     x <- GetPlotFractionalCoords(0.125, 0.5)[1]
#     mtext(text=mirna_txt, side=3, line=0.8, at=x, adj=0, cex=par("cex"))
#   }
#   if (site_label) {
#     if (site_type == "compensatory") {
#       site_txt <- "All 8mer-mm sites"
#     } else {
#       site_txt <- site_type
#     }
#     x <- GetPlotFractionalCoords(0.8, 0.5)[1]
#     mtext(text=site_txt, side=3, line=0, at=x, adj=1, cex=par("cex"), xpd=NA)
#   }
#   # Boolean label giving the name of the site, in the top right-hand corner,
#   # next to the miRNA.

#   # Finish the plot.
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
#   return(R_mat)
# }



# MakeLoopPSAM <- function(pos, len, offset, conditions="all", n_constant=3,
#                          sitelist="progthrp", dinuc=FALSE, dinuc_sim=FALSE,
#                          trinuc=FALSE, trinuc_sim=FALSE, titles=FALSE,
#                          titles2=FALSE, titles3=FALSE, height=2.8, width=2.5,
#                          xpos=20, ypos=20, pdf.plot=FALSE) {
#   # Define the extension string.
#   R_mat <- log2(SubfunctionCall(MakeLoopNucleotideEnrichmentMatrix))
#   rownames(R_mat) <- gsub("T", "U", rownames(R_mat))
#   base_width <- pos - 8 - 1 + offset
#   if (trinuc | trinuc_sim)    base_width <- base_width - 2
#   else if (dinuc | dinuc_sim) base_width <- base_width - 1
#   # if (bulge)  base_width <- base_width - 1
#   if (titles) mar1 <- 2.2
#   else        mar1 <- 1.7
#   if (trinuc | trinuc_sim)  mar2 <- 2.5
#   else                      mar2 <- 1.7
#   if (titles3 | titles2) mar3 <- 2.7
#   else                   mar3 <- 1
#   mar4 <- 0.2
#   height <- 2.35^(2.6*(trinuc | trinuc_sim) + ((dinuc & !trinuc_sim) | dinuc_sim)) + (mar1 + mar3)/5
#   width <- base_width/5 + (mar2 + mar4)/5
#   print(base_width)
#   print(height)
#   print(width)
#   SubfunctionCall(FigureSaveFile2)
#   par(mar=c(mar1, mar2, mar3, mar4))
#   xmin <- 0
#   xmax <- base_width
#   ymin <- 0
#   if (trinuc | trinuc_sim)    ymax <- 64
#   else if (dinuc | dinuc_sim) ymax <- 16
#   else                        ymax <- 4
#   BlankPlot()
#   # if (dinuc | dinuc_sim) ymax <- 16
#   # else                   ymax <- 4
#   # Define the positions of the rectangles in the matrix.
#   xlefts <- rep(seq(0, ncol(R_mat) - 1), nrow(R_mat))
#   x_axis_label_pos <- xmin:(xmin + ncol(R_mat) - 1) + 0.5    
#   xright <- xlefts + 1
#   ybottom <- rep(rev(seq(ymax - 1, ymin)), each=ncol(R_mat))
#   ytop <- ybottom + 1
#   # Make the x-axis
#   if (titles) x_label <- "Position"
#   else        x_label <- ""
#   labels_inds <- seq(1, length(x_axis_label_pos))
#   AddLinearAxis(1, 1, 3, label=x_label,
#                 label_pos_ticks=TRUE,
#                 alt_lab=colnames(R_mat),
#                 alt_lab_y_dist=0.01,
#                 alt_lab_pos=x_axis_label_pos[labels_inds],
#                 alt_tick_pos=TRUE,
#                 alt_tick_lab_space=TRUE)
#   print('made it past axis')
#   # Make the y-axis ############################################################
#   par_mgp <- par("mgp")
#   par_mgp[2] <- 0.2
#   par(mgp=par_mgp)
#   yaxis_labs <- paste0(rownames(R_mat), " ")
#   AddLinearAxis(2, 1, 1, label="",
#                 alt_lab=yaxis_labs,
#                 alt_lab_pos=ymin:(ymax - 1) + 0.5,
#                 alt_tick_pos=TRUE)
#   if (trinuc | trinuc_sim)      col.inds <- floor((t(R_mat) + 3.5)/7*99 + 1)
#   else if (dinuc | dinuc_sim)   col.inds <- floor((t(R_mat) + 2  )/4*99 + 1)
#   else                          col.inds <- floor((t(R_mat) + 1.5)/3*99 + 1)    
#   # Make distribution of 100 colors and a 101th gray color for NA indeces.
#   color.dist <- c(rev(colorspace::diverging_hcl(100, h=c(180, 50), c=80,
#                                               l=c(20, 95), power=c(0.7, 1.3))),
#                   "gray90", kMaxValueColor)
#   color.dist[100] <- "forestgreen"
#   # Make the color index scale.
#   col.inds[which(col.inds > 100)] <- 100
#   col.inds[which(col.inds < 1)] <- 102
#   col.inds[which(is.na(col.inds))] <- 101
#   # Make the color rectangle.
#   rect(xlefts, ybottom, xright, ytop, col=color.dist[col.inds], lwd=0.4, 
#        border="white", xpd=NA)
#   # Add the titles for the length, position, and offset of pairing shown.
#   if (titles3) {
#     mtext(sprintf("pos %s", pos), 3, line=1.8, adj=1, cex=par("cex"), xpd=NA)
#   }
#   if (titles2 | titles3) {
#     mtext(sprintf("%s bp", len), 3, line=0.9, adj=1, cex=par("cex"), xpd=NA)
#   }
#   mtext(sprintf("%s-nt offset", offset), 3, line=0, adj=1, cex=par("cex"),
#         xpd=NA)
#   # Finish the plot.
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
#   return(R_mat)
# }





# PlotModelRsquareds <- function(
#   cell_line, mirnas, height=2.9, width=3.5, xpos=20, ypos=20, pdf.plot=FALSE
# ) {
#   if (mirnas == "sixteen") {
#     conditions_matrix <- matrix(c("biochem", FALSE,
#                                   "biochem", TRUE,
#                                   "biochemplus", FALSE,
#                                   "biochemplus", TRUE),
#                                 byrow=TRUE,
#                                 ncol=2)
#     out <- apply(conditions_matrix, 1, function(row) {
#       GetModelPredictionCor(cell_line=cell_line, model=row[1], kds="predicted",
#                             mirnas=mirnas, passenger=row[2])
#     })

#   } else {
#     conditions_matrix <- matrix(c("biochem", "measured",
#                                   "biochem", "predicted",
#                                   "biochemplus", "measured",
#                                   "biochemplus", "predicted"),
#                                 byrow=TRUE,
#                                 ncol=2)
#     out <- apply(conditions_matrix, 1, function(row) {
#       GetModelPredictionCor(cell_line=cell_line, model=row[1], kds=row[2],
#                             mirnas=mirnas)
#     })
#   }
#   print(conditions_matrix)
#   print(out)
#   SubfunctionCall(FigureSaveFile2)
#   par(mar=c(3, 3, 0.5, 0.5))
#   xmin <- 0
#   xmax <- 4
#   ymin <- 0
#   ymax <- 0.5
#   BlankPlot()
#   xleft <- 0:3 + 0.25
#   xright <- xleft + 0.5
#   AddLinearAxis(2, 0.05, 0.1, label="r2")
#   # Get the left- and right-hand positions of the rectangles.
#   # Make the barplot rectangles ################################################
#   rect(xleft=xleft, ybottom=0, xright=xright, ytop=out, col="black", border=NA)
#   if (mirnas == "sixteen") {
#     text(x = -0.1, y=-0.05, labels=c("Star:"), adj=c(1, 0.5), xpd=NA)
#     text(x = 0:3 + 0.5, y=-0.05, labels=c("-", "+", "-", "+"), xpd=NA)
#   } else {
#     text(x = -0.1, y=-0.05, labels=c("Kds:"), adj=c(1, 0.5), xpd=NA)
#     text(x = 0:3 + 0.5, y=-0.05, labels=c("RBNS", "CNN", "RBNS", "CNN"), xpd=NA)
#   }

#   text(x = -0.1, y=-0.11, labels=c("Model:"), adj=c(1, 0.5), xpd=NA)
#   segments(x0=c(0, 2) + 0.125, y0=-0.08, x1=c(2, 4) - 0.125, xpd=NA)
#   text(x = c(1, 3), y=-0.11, labels=c("biochem", "biochem+"), xpd=NA)
#   xy <- GetPlotFractionalCoords(fx=0.05, fy=0.95)
#   text(x=xy[1], y=xy[2], labels=cell_line, xpd=NA, adj=c(0, 1))
#   xy <- GetPlotFractionalCoords(fx=0.95, fy=0.95)
#   text(x=xy[1], y=xy[2], labels=mirnas, xpd=NA, adj=c(1, 1))


#   # Finish the plot.
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }


# PlotThreePrimeScoreFeatures <- function(
#   cell_line, model, kds, mirnas, passenger=FALSE, model_coefs=FALSE, two_bmodes=NULL,
#   min_pairing=NULL, offset_opt=NULL, supp_l=NULL, supp_r=NULL, offset_tol=NULL,
#   w_pairing=NULL, w_offset=NULL, w_supp=NULL, plot_labels=FALSE, height=2.95,
#   width=2.95, xpos=20, ypos=20, pdf.plot=FALSE
# ) {
#   if (model_coefs) {
#     x_use <- c(1, 2, 3)
#     x_axis_lab <- "Model coefficients"
#     tick_space <- 1
#     lab_space <- 1
#     xmin <- 0
#     xmax <- 4
#     x_default <- 1
#     outs <- c(
#       GetModelPredictionCor(
#         cell_line, "biochemplus", kds, mirnas, passenger=passenger
#       ),
#       GetModelPredictionCor(
#         cell_line, "biochemplus", kds, mirnas, passenger=passenger,
#         model_coefs="programmed"
#       ),
#       GetModelPredictionCor(
#         cell_line, "biochemplus", kds, mirnas, passenger=passenger,
#         model_coefs="random"
#       )
#       )


#   }
#   if (!is.null(two_bmodes)) {
#     x_use <- c(1, 2)
#     x_axis_lab <- "Two binding modes"
#     tick_space <- 1
#     lab_space <- 1
#     xmin <- 0
#     xmax <- 3
#     x_default <- 1
#     outs <- c(
#       GetModelPredictionCor(
#         cell_line, "biochemplus", kds, mirnas, passenger=passenger
#       ),
#       GetModelPredictionCor(
#         cell_line, "biochemplus", kds, mirnas, passenger=passenger,
#         two_bmodes=TRUE
#       )
#     )
#   } else if (!is.null(min_pairing)) {
#     x_use <- min_pairing
#     x_axis_lab <- "Minimum pairing length"
#     xmin <- 0
#     xmax <- 8
#     tick_space <- 1
#     lab_space <- 1
#     x_default <- 2
#     outs <- sapply(min_pairing, function(min_pairing_i) {
#       # if (min_pairing_i != 2) {
#         GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#                               passenger=passenger, min_pairing=min_pairing_i)
#       # } else {
#       #   GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#       #                         passenger=passenger)

#       # }
#     })
#   } else if (!is.null(offset_opt)) {
#     x_use <- offset_opt
#     x_axis_lab <- "Optimal offset"
#     xmin <- 0
#     xmax <- 8
#     tick_space <- 1
#     lab_space <- 1
#     x_default <- 0
#     outs <- sapply(offset_opt, function(offset_opt_i) {
#       if (offset_opt_i != 0) {
#         GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#                               passenger=passenger, offset_opt=offset_opt_i)
#       } else {
#         GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#                               passenger=passenger)

#       }
#     })
#   } else if (!is.null(supp_l) & is.null(supp_r)) {
#     x_use <- supp_l
#     x_axis_lab <- "5' end of supp. region"
#     xmin <- 9
#     xmax <- 13
#     tick_space <- 1
#     lab_space <- 1
#     x_default <- 13
#     outs <- sapply(supp_l, function(supp_l_i) {
#       if (supp_l_i != 13) {
#         GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#                               passenger=passenger, supp_l=supp_l_i)
#       } else {
#         GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#                               passenger=passenger)

#       }
#     })
#   } else if (!is.null(supp_r) & is.null(supp_l)) {
#     x_use <- supp_r
#     x_axis_lab <- "3' end of supp. region"
#     xmin <- 16
#     xmax <- 22
#     tick_space <- 1
#     lab_space <- 1
#     x_default <- 16
#     outs <- sapply(supp_r, function(supp_r_i) {
#       if (supp_r_i != 16) {
#         GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#                               passenger=passenger, supp_r=supp_r_i)
#       } else {
#         GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#                               passenger=passenger)

#       }
#     })
#   } else if (!is.null(supp_r) & !is.null(supp_l)) {
#     x_use <- supp_l
#     x_axis_lab <- "Entire supp. region"
#     xmin <- 9
#     xmax <- 16
#     tick_space <- 1
#     lab_space <- 1
#     supp_lr <- cbind(supp_l, supp_r)
#     x_default <- 13
#     outs <- apply(supp_lr, 1, function(supp_lr_i) {
#       if (supp_lr_i[1] != 13) {
#         GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#                               passenger=passenger, supp_l=supp_lr_i[1],
#                               supp_r=supp_lr_i[2])
#       } else {
#         GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#                               passenger=passenger)

#       }
#     })
#   } else if (!is.null(offset_tol)) {
#     x_use <- offset_tol
#     x_axis_lab <- "Offset tolerance"
#     xmin <- 0
#     xmax <- 10
#     tick_space <- 1
#     lab_space <- 1
#     x_default <- 2
#     outs <- sapply(offset_tol, function(offset_tol_i) {
#       if (offset_tol_i != 2) {
#         GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#                               passenger=passenger, offset_tol=offset_tol_i)
#       } else {
#         GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#                               passenger=passenger)

#       }
#     })
#   } else if (!is.null(w_pairing)) {
#     x_use <- w_pairing
#     x_axis_lab <- "Weight of pairing"
#     xmin <- 0
#     xmax <- 1
#     tick_space <- 0.1
#     lab_space <- 0.1
#     x_default <- 0.5
#     outs <- sapply(w_pairing, function(w_pairing_i) {
#       if (w_pairing_i != 0.5) {
#         GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#                               passenger=passenger, w_pairing=w_pairing_i)
#       } else {
#         GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#                               passenger=passenger)

#       }
#     })
#   } else if (!is.null(w_offset)) {
#     x_use <- w_offset
#     x_axis_lab <- "Weight of offset"
#     xmin <- 0
#     xmax <- 1
#     tick_space <- 0.1
#     lab_space <- 0.1
#     x_default <- 0.5
#     outs <- sapply(w_offset, function(w_offset_i) {
#       if (w_offset_i != 0.5) {
#         GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#                               passenger=passenger, w_offset=w_offset_i)
#       } else {
#         GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#                               passenger=passenger)

#       }
#     })
#   } else if (!is.null(w_supp)) {
#     x_use <- w_supp
#     x_axis_lab <- "Weight of supp. region"
#     xmin <- 0
#     xmax <- 1
#     tick_space <- 0.1
#     lab_space <- 0.1
#     x_default <- 0.5
#     outs <- sapply(w_supp, function(w_supp_i) {
#       if (w_supp_i != 0.5) {
#         GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#                               passenger=passenger, w_supp=w_supp_i)
#       } else {
#         GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#                               passenger=passenger)

#       }
#     })
#   }
#   outs <- outs*100
#   print(outs)
#   SubfunctionCall(FigureSaveFile2)
#   par(mar=c(3, 3, 0.5, 0.5))
#   if (mirnas == "five") {
#     if (kds == "measured") {
#       ymin <- 36
#       ymax <- 37
#     } else {
#       ymin <- 39
#       ymax <- 40  
#     }
#   } else if (mirnas == "six") {
#     ymin <- 33.7
#     ymax <- 34.7
#   } else if (mirnas == "sixteen") {
#     ymin <- 32
#     ymax <- 33
#   }
#   BlankPlot()
#   AddLinearAxis(1, tick_space, lab_space, label=x_axis_lab)
#   AddLinearAxis(2, .05, .1, label="r2")
#   cols <- rep("black", length(x_use))
#   cols[which(x_use == x_default)] <- "red"
#   lines(x_use, outs)
#   points(x_use, outs, pch=19, col=cols)
#   if (plot_labels) {
#     xy <- GetPlotFractionalCoords(0.95, 1)
#     text(xy[1], xy[2], labels=sprintf("%s miRNAs", mirnas), adj=c(1, 1), xpd=NA)
#     xy <- GetPlotFractionalCoords(0.95, 0.90)
#     if (kds == "measured") kd_label <- "RBNS Kds"
#     if (kds == "predicted") kd_label <- "CNN Kds"
#     text(xy[1], xy[2], labels=kd_label, adj=c(1, 1), xpd=NA)
#   }
#   # Finish the plot.
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }




# PlotBiochemicalModelFeatureWeights <- function(xpos=20, ypos=20, height=3,
#   width=3.5, pdf.plot=FALSE
# ) {
#   model1 <- GetModelFit("HeLa", "biochemplus", "measured", "five")
#   model2 <- GetModelFit("HeLa", "biochemplus", "predicted", "five")
#   model3 <- GetModelFit("HeLa", "biochemplus", "predicted", "sixteen")
#   out <- do.call("rbind", lapply(list(model1, model2, model3), function(model_i) {
#     print(model_i)
#     return(model_i[c("logSA_diff", "Threep_canon", "PCT")])
#   }))
#   print(out)
#   out <- c(t(out))
#   SubfunctionCall(FigureSaveFile2)
#   par(mar=c(3, 3.2, 0.5, 0.5))
#   xmin <- 0
#   xmax <- 11
#   ymin <- 0
#   ymax <- 2
#   BlankPlot()
#   xleft <- c(0, 1, 2, 4, 5, 6, 8, 9, 10)
#   xright <- xleft + 1
#   AddLinearAxis(2, 0.1, 0.1, label="Feature weight")
#   # Get the left- and right-hand positions of the rectangles.
#   # Make the barplot rectangles ################################################
#   rect(xleft=xleft, ybottom=0, xright=xright, ytop=out, col=c("red", "purple", "blue"), border=NA)
#     text(x=0.2, y=-0.15, labels=c("Kds:"), adj=c(1, 0.5), xpd=NA)
#     text(x = c(1.5, 5.5, 9.5), y=-0.15, labels=c("RBNS", "CNN", "CNN"), xpd=NA)

#   text(x=0.2, y=-0.3, labels=c("miRNAs:"), adj=c(1, 0.5), xpd=NA)
#   text(x = c(1.5, 5.5, 9.5), y=-0.3, labels=c("five", "five", "sixteen"), xpd=NA)
#   xy <- GetPlotFractionalCoords(fx=0.5, fy=0.5)

#   text(x=11, y=c(2, 1.8, 1.6), labels=c("SA", "Threep", "PCT"), col=c("red", "purple", "blue"), adj=c(1, 1), xpd=NA)
#   # Finish the plot.
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }

# }




# PlotTargetScanThreePrimeScoreFeatures <- function(
#   cell_line, mirnas, train_mirnas, alt=FALSE, rescaled=FALSE, bounded=FALSE, two_bmodes=NULL,
#   min_pairing=NULL, offset_opt=NULL, supp_l=NULL, supp_r=NULL, offset_tol=NULL,
#   w_pairing=NULL, w_offset=NULL, w_supp=NULL, plot_labels=FALSE, height=2.95,
#   width=2.95, xpos=20, ypos=20, pdf.plot=FALSE
# ) {
#   if (!is.null(two_bmodes)) {
#     x_use <- c(1, 2)
#     x_axis_lab <- "Two binding modes"
#     tick_space <- 1
#     lab_space <- 1
#     xmin <- 0
#     xmax <- 3
#     x_default <- 1
#     outs <- c(
#       GetTargetScanR2(
#         cell_line, mirnas, train_mirnas, alt=alt, rescaled=rescaled, bounded=bounded
#       ),
#       GetTargetScanR2(
#         cell_line, mirnas, train_mirnas, alt=alt, rescaled=rescaled, bounded=bounded, two_bmodes=TRUE
#       )
#     )
#   } else if (!is.null(min_pairing)) {
#     x_use <- min_pairing
#     x_axis_lab <- "Minimum pairing length"
#     xmin <- 0
#     xmax <- 8
#     tick_space <- 1
#     lab_space <- 1
#     x_default <- 2
#     outs <- sapply(min_pairing, function(min_pairing_i) {
#       # if (min_pairing_i != 2) {
#         # GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#         #                       passenger=passenger, min_pairing=min_pairing_i)


#         GetTargetScanR2(
#           cell_line, mirnas, train_mirnas, alt=alt, rescaled=rescaled, bounded=bounded,
#           min_pairing=min_pairing_i
#         )

#       # } else {
#       #   GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#       #                         passenger=passenger)

#       # }
#     })
#   } else if (!is.null(offset_opt)) {
#     x_use <- offset_opt
#     x_axis_lab <- "Optimal offset"
#     xmin <- 0
#     xmax <- 8
#     tick_space <- 1
#     lab_space <- 1
#     x_default <- 0
#     outs <- sapply(offset_opt, function(offset_opt_i) {
#       if (offset_opt_i != 0) {
#         # GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#         #                       passenger=passenger, offset_opt=offset_opt_i)

#         GetTargetScanR2(
#           cell_line, mirnas, train_mirnas, alt=alt, rescaled=rescaled, bounded=bounded,
#           offset_opt=offset_opt_i
#         )

#       } else {
#         # GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#         #                       passenger=passenger)

#         GetTargetScanR2(
#           cell_line, mirnas, train_mirnas, alt=alt, rescaled=rescaled, bounded=bounded,
#         )



#       }
#     })
#   } else if (!is.null(supp_l) & is.null(supp_r)) {
#     x_use <- supp_l
#     x_axis_lab <- "5' end of supp. region"
#     xmin <- 9
#     xmax <- 13
#     tick_space <- 1
#     lab_space <- 1
#     x_default <- 13
#     outs <- sapply(supp_l, function(supp_l_i) {
#       if (supp_l_i != 13) {
#         # GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#         #                       passenger=passenger, supp_l=supp_l_i)
#         GetTargetScanR2(
#           cell_line, mirnas, train_mirnas, alt=alt, rescaled=rescaled, bounded=bounded,
#           supp_l=supp_l_i
#         )
#       } else {
#         # GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#         #                       passenger=passenger)
#         GetTargetScanR2(
#           cell_line, mirnas, train_mirnas, alt=alt, rescaled=rescaled, bounded=bounded
#         )

#       }
#     })
#   } else if (!is.null(supp_r) & is.null(supp_l)) {
#     x_use <- supp_r
#     x_axis_lab <- "3' end of supp. region"
#     xmin <- 16
#     xmax <- 22
#     tick_space <- 1
#     lab_space <- 1
#     x_default <- 16
#     outs <- sapply(supp_r, function(supp_r_i) {
#       if (supp_r_i != 16) {
#         # GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#         #                       passenger=passenger, supp_r=supp_r_i)
#         GetTargetScanR2(
#           cell_line, mirnas, train_mirnas, alt=alt, rescaled=rescaled, bounded=bounded,
#           supp_r=supp_r_i
#         )
#       } else {
#         # GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#         #                       passenger=passenger)
#         GetTargetScanR2(
#           cell_line, mirnas, train_mirnas, alt=alt, rescaled=rescaled, bounded=bounded
#         )

#       }
#     })
#   } else if (!is.null(supp_r) & !is.null(supp_l)) {
#     x_use <- supp_l
#     x_axis_lab <- "Entire supp. region"
#     xmin <- 9
#     xmax <- 16
#     tick_space <- 1
#     lab_space <- 1
#     supp_lr <- cbind(supp_l, supp_r)
#     x_default <- 13
#     outs <- apply(supp_lr, 1, function(supp_lr_i) {
#       if (supp_lr_i[1] != 13) {
#         # GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#         #                       passenger=passenger, supp_l=supp_lr_i[1],
#         #                       supp_r=supp_lr_i[2])
#         GetTargetScanR2(
#           cell_line, mirnas, train_mirnas, alt=alt, rescaled=rescaled, bounded=bounded,
#           supp_l=supp_lr_i[1], supp_r=supp_lr_i[2]
#         )
#       } else {
#         # GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#         #                       passenger=passenger)
#         GetTargetScanR2(
#           cell_line, mirnas, train_mirnas, alt=alt, rescaled=rescaled, bounded=bounded
#         )

#       }
#     })
#   } else if (!is.null(offset_tol)) {
#     x_use <- offset_tol
#     x_axis_lab <- "Offset tolerance"
#     xmin <- 0
#     xmax <- 10
#     tick_space <- 1
#     lab_space <- 1
#     x_default <- 2
#     outs <- sapply(offset_tol, function(offset_tol_i) {
#       if (offset_tol_i != 2) {
#         # GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#         #                       passenger=passenger, offset_tol=offset_tol_i)
#         GetTargetScanR2(
#           cell_line, mirnas, train_mirnas, alt=alt, rescaled=rescaled, bounded=bounded,
#           offset_tol=offset_tol_i
#         )
#       } else {
#         # GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#         #                       passenger=passenger)
#         GetTargetScanR2(
#           cell_line, mirnas, train_mirnas, alt=alt, rescaled=rescaled, bounded=bounded
#         )

#       }
#     })
#   } else if (!is.null(w_pairing)) {
#     x_use <- w_pairing
#     x_axis_lab <- "Weight of pairing"
#     xmin <- 0
#     xmax <- 1
#     tick_space <- 0.1
#     lab_space <- 0.1
#     x_default <- 0.5
#     outs <- sapply(w_pairing, function(w_pairing_i) {
#       if (w_pairing_i != 0.5) {
#         # GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#         #                       passenger=passenger, w_pairing=w_pairing_i)
#         GetTargetScanR2(
#           cell_line, mirnas, train_mirnas, alt=alt, rescaled=rescaled, bounded=bounded,
#           w_pairing=w_pairing_i
#         )
#       } else {
#         # GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#         #                       passenger=passenger)
#         GetTargetScanR2(
#           cell_line, mirnas, train_mirnas, alt=alt, rescaled=rescaled, bounded=bounded
#         )

#       }
#     })
#   } else if (!is.null(w_offset)) {
#     x_use <- w_offset
#     x_axis_lab <- "Weight of offset"
#     xmin <- 0
#     xmax <- 1
#     tick_space <- 0.1
#     lab_space <- 0.1
#     x_default <- 0.5
#     outs <- sapply(w_offset, function(w_offset_i) {
#       if (w_offset_i != 0.5) {
#         # GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#         #                       passenger=passenger, w_offset=w_offset_i)
#         GetTargetScanR2(
#           cell_line, mirnas, train_mirnas, alt=alt, rescaled=rescaled, bounded=bounded,
#           w_offset=w_offset_i
#         )
#       } else {
#         # GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#         #                       passenger=passenger)
#         GetTargetScanR2(
#           cell_line, mirnas, train_mirnas, alt=alt, rescaled=rescaled, bounded=bounded
#         )

#       }
#     })
#   } else if (!is.null(w_supp)) {
#     x_use <- w_supp
#     x_axis_lab <- "Weight of supp. region"
#     xmin <- 0
#     xmax <- 1
#     tick_space <- 0.1
#     lab_space <- 0.1
#     x_default <- 0.5
#     outs <- sapply(w_supp, function(w_supp_i) {
#       if (w_supp_i != 0.5) {
#         # GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#         #                       passenger=passenger, w_supp=w_supp_i)
#         GetTargetScanR2(
#           cell_line, mirnas, train_mirnas, alt=alt, rescaled=rescaled, bounded=bounded,
#           w_supp=w_supp_i
#         )
#       } else {
#         # GetModelPredictionCor(cell_line, "biochemplus", kds, mirnas,
#         #                       passenger=passenger)
#         GetTargetScanR2(
#           cell_line, mirnas, train_mirnas, alt=alt, rescaled=rescaled, bounded=bounded
#         )

#       }
#     })
#   }
#   outs <- outs*100
#   print(outs)
#   SubfunctionCall(FigureSaveFile2)
#   par(mar=c(3, 3, 0.5, 0.5))
#   if (mirnas == "five") {
#     ymin <- 27
#     ymax <- 28
#   } else if (mirnas == "six") {
#     ymin <- 25
#     ymax <- 26
#   } else if (mirnas == "sixteen") {
#     ymin <- 21
#     ymax <- 22
#   }
#   BlankPlot()
#   AddLinearAxis(1, tick_space, lab_space, label=x_axis_lab)
#   AddLinearAxis(2, .05, .1, label="r2")
#   cols <- rep("black", length(x_use))
#   cols[which(x_use == x_default)] <- "red"
#   lines(x_use, outs)
#   points(x_use, outs, col=cols, pch=19)
#   if (plot_labels) {
#     xy <- GetPlotFractionalCoords(0.95, 1)
#     text(xy[1], xy[2], labels=sprintf("%s miRNAs", mirnas), adj=c(1, 1), xpd=NA)
#     xy <- GetPlotFractionalCoords(0.95, 0.90)
#     # if (kds == "measured") kd_label <- "RBNS Kds"
#     # if (kds == "predicted") kd_label <- "CNN Kds"
#     # text(xy[1], xy[2], labels=kd_label, adj=c(1, 1), xpd=NA)
#   }
#   # Finish the plot.
#   if (class(pdf.plot) == "character") {
#     dev.off()
#   }
# }


## MAIN TEXT FIGURES ###########################################################
MakeFigure1 <- function() {
  message("Making Fig. 1.")
  # A.__________________________________________________________________________
  # Illustrator schematic.
  # B.__________________________________________________________________________
  # Illustrator schematic.
  # C.__________________________________________________________________________
  # Illustrator schematic.
  # D.__________________________________________________________________________
  PlotPositionalEnrichmentForProgramedLibrary(
    "let-7a-21nt", "equil_c2_nb", "40", 0, 8, stop=31, ps=0, pdf.plot="1.D"
  )
  message("Done Fig. 1.")
}


MakeFigure2 <- function(corrected_kds=TRUE) {
  message("Making Fig. 2")
  # A.__________________________________________________________________________
  PlotCombinedEnrichmentAndKd("8mer-m11.18", plot_enrich=TRUE, pdf.plot="2.Ai")
  PlotCombinedEnrichmentAndKd("8mer-m11.18", pdf.plot="2.Aii")
  # B.__________________________________________________________________________
  PlotBestThreePrimeSite(height=4, pdf.plot="2.B")
  # C.__________________________________________________________________________
  PlotOneThreePrimeSite(4, "let-7a-21nt", "equil_c2_nb", pdf.plot="2.Ci")
  PlotOneThreePrimeSite(5, "let-7a-21nt", "equil_c2_nb", pdf.plot="2.Cii")
  PlotOneThreePrimeSite(6, "let-7a-21nt", "equil_c2_nb", pdf.plot="2.Ciii")
  PlotOneThreePrimeSite(7, "let-7a-21nt", "equil_c2_nb", pdf.plot="2.Civ")
  PlotOneThreePrimeSite(8, "let-7a-21nt", "equil_c2_nb", pdf.plot="2.Cv")
  PlotOneThreePrimeSite(9, "let-7a-21nt", "equil_c2_nb", pdf.plot="2.Cvi")
  # D.__________________________________________________________________________
  PlotPairingMatrix("let-7a-21nt", "equil_c2_nb", -2,
                    corrected_kds=corrected_kds, pdf.plot="2.Di")
  PlotPairingMatrix("let-7a-21nt", "equil_c2_nb", 0,
                    corrected_kds=corrected_kds, pdf.plot="2.Dii")
  PlotPairingMatrix("let-7a-21nt", "equil_c2_nb",  2,
                    corrected_kds=corrected_kds, pdf.plot="2.Diii")
  PlotPairingMatrix("let-7a-21nt", "equil_c2_nb",  4,
                    corrected_kds=corrected_kds, pdf.plot="2.Div")
  R_mat4 <- R_mat
  message(sprintf("Kd fold change of 11mer-m10.20 a offset 4: %s",
                  10^R_mat["10", "20"]))
  message(sprintf("If that is the max this value should be the same: %s",
                  10^max(R_mat, na.rm=TRUE)))
  PlotPairingMatrix("let-7a-21nt", "equil_c2_nb",  6,
                    corrected_kds=corrected_kds, pdf.plot="2.Dv")
  PlotPairingMatrix("let-7a-21nt", "equil_c2_nb",  8,
                    corrected_kds=corrected_kds, pdf.plot="2.Dvi")
  PlotPairingMatrix("let-7a-21nt", "equil_c2_nb",  10,
                    corrected_kds=corrected_kds, pdf.plot="2.Dvii")
  PlotPairingMatrix("let-7a-21nt", "equil_c2_nb", 12, key=TRUE,
                    corrected_kds=corrected_kds, pdf.plot="2.Dviii")
  message("Done Fig. 2")
}

MakeFigure3 <- function(experiment="twist_reporter_assay_3p_2_tp",
                        bulge_nucs="both") {
  message("Making Fig. 3")
  # A.__________________________________________________________________________
  # Illustrator schematic.
  # B.__________________________________________________________________________
  # PlotMeanReporterRepression(experiment=experiment, pdf.plot="3.Bi")
  # PlotMeanReporterRepression(experiment=experiment, dual_site=TRUE, pdf.plot="3.Bii")
  # PrintTukeyPValues()
  # PrintTukeyPValues(all_lin41_UTR=TRUE)

  # # C.__________________________________________________________________________  
  PlotReporterEfficacyAgainstProgKds(experiment_rep=experiment, pdf.plot="3.C")
  # message("Done Fig. 3")
}


MakeFigure4 <- function(corrected_kds=TRUE) {
  message("Making Fig. 4.")
  # A___________________________________________________________________________
  # PlotBestThreePrimeSite(mirna="miR-1", exp="equil_c_nb", alt_height=TRUE,
  #                        pdf.plot="4.A")
  # B___________________________________________________________________________
  PlotBestThreePrimeSite(mirna="miR-155", exp="equil_sc_nb", alt_height=TRUE,
                         pdf.plot="4.B")
  # C___________________________________________________________________________
  PlotPairingMatrix("miR-1", "equil_c_nb", -2, mirna_label=TRUE,
                         corrected_kds=corrected_kds, pdf.plot="4.Ci")
  data_max_value <<- max(10^R_mat, na.rm=TRUE)

  PlotPairingMatrix("miR-1", "equil_c_nb", 0,
                         corrected_kds=corrected_kds, pdf.plot="4.Cii")
  data_max_value <<- max(data_max_value, max(10^R_mat, na.rm=TRUE))
  PlotPairingMatrix("miR-1", "equil_c_nb",  2,
                         corrected_kds=corrected_kds, pdf.plot="4.Ciii")
  data_max_value <<- max(data_max_value, max(10^R_mat, na.rm=TRUE))
  PlotPairingMatrix("miR-1", "equil_c_nb",  4, corrected_kds=corrected_kds,
                    pdf.plot="4.Civ")
  R_mat4 <- R_mat
  message(sprintf("Kd fold change of 11mer-m10.20 a offset 4: %s",
                  10^R_mat["10", "20"]))
  message(sprintf("If that is the max this value should be the same: %s",
                  10^max(R_mat, na.rm=TRUE)))
  PlotPairingMatrix("miR-1", "equil_c_nb",  6, corrected_kds=corrected_kds,
                    pdf.plot="4.Cv")
  data_max_value <<- max(data_max_value, max(10^R_mat, na.rm=TRUE))
  PlotPairingMatrix("miR-1", "equil_c_nb",  8, corrected_kds=corrected_kds,
                    pdf.plot="4.Cvi")
  data_max_value <<- max(data_max_value, max(10^R_mat, na.rm=TRUE))
  PlotPairingMatrix("miR-1", "equil_c_nb",  10, corrected_kds=corrected_kds,
                    pdf.plot="4.Cvii")
  data_max_value <<- max(data_max_value, max(10^R_mat, na.rm=TRUE))
  PlotPairingMatrix("miR-1", "equil_c_nb", 12, corrected_kds=corrected_kds,
                    key=TRUE, alt_top=TRUE, max_value_check=data_max_value,
                    pdf.plot="4.Cviii")
  # D___________________________________________________________________________
  PlotPairingMatrix("miR-155", "equil_sc_nb", -2, mirna_label=TRUE,
                         corrected_kds=corrected_kds, pdf.plot="4.Di")
  print(data_max_value)
  data_max_value <<- max(10^R_mat, na.rm=TRUE)

  PlotPairingMatrix("miR-155", "equil_sc_nb", 0, corrected_kds=corrected_kds,
                    pdf.plot="4.Dii")
  data_max_value <<- max(data_max_value, max(10^R_mat, na.rm=TRUE))
  print(data_max_value)
  PlotPairingMatrix("miR-155", "equil_sc_nb",  2, corrected_kds=corrected_kds,
                    pdf.plot="4.Diii")
  data_max_value <<- max(data_max_value, max(10^R_mat, na.rm=TRUE))
  print(data_max_value)

  PlotPairingMatrix("miR-155", "equil_sc_nb",  4, corrected_kds=corrected_kds,
                    pdf.plot="4.Div")
  data_max_value <<- max(data_max_value, max(10^R_mat, na.rm=TRUE))
  print(data_max_value)
  R_mat4 <- R_mat
  message(sprintf("Kd fold change of 11mer-m10.20 a offset 4: %s",
                  10^R_mat["10", "20"]))
  message(sprintf("If that is the max this value should be the same: %s",
                  10^max(R_mat, na.rm=TRUE)))
  PlotPairingMatrix("miR-155", "equil_sc_nb",  6, corrected_kds=corrected_kds,
                    pdf.plot="4.Dv")
  print(data_max_value)
  data_max_value <<- max(data_max_value, max(10^R_mat, na.rm=TRUE))
  PlotPairingMatrix("miR-155", "equil_sc_nb",  8, corrected_kds=corrected_kds,
                    pdf.plot="4.Dvi")
  print(data_max_value)
  data_max_value <<- max(data_max_value, max(10^R_mat, na.rm=TRUE))
  PlotPairingMatrix("miR-155", "equil_sc_nb",  10, corrected_kds=corrected_kds,
                    pdf.plot="4.Dvii")
  print(data_max_value)
  data_max_value <<- max(data_max_value, max(10^R_mat, na.rm=TRUE))
  PlotPairingMatrix("miR-155", "equil_sc_nb", 12, corrected_kds=corrected_kds,
                    key=TRUE, alt_top=TRUE, max_value_check=data_max_value,
                    pdf.plot="4.Dviii")
  print(data_max_value)

  message("Done Fig. 4.")
}


MakeFigure5 <- function(corrected_kds=TRUE) {
  # message("Making Fig. 5.")
  # # A___________________________________________________________________________
  # PlotPairingCoefficients("let-7a-21nt", "equil_c2_nb", pdf.plot="5.Ai")
  # model_l7 <<- model
  # PlotOffsetCoefficients("let-7a-21nt", "equil_c2_nb", makeglobalmodel=FALSE,
  #                        pdf.plot="5.Aii")
  # PlotPairingMatrix("let-7a-21nt", "equil_c2_nb", 4, model_values=TRUE,
  #                   makeglobalmodel=FALSE, corrected_kds=corrected_kds,
  #                   pdf.plot="5.Aiii")
  # # R_let7_model <<- R_mat
  # PlotPairingMatrix("let-7a-21nt", "equil_c2_nb", 4,
  #                   corrected_kds=corrected_kds, pdf.plot="5.Aiv")
  # R_let7_data <<- R_mat
  # PlotPairingMatrixKey(pdf.plot="5.Av")
  # PlotPairingMatrixKey(model=FALSE, pdf.plot="5.Avi")
  # # B___________________________________________________________________________
  # PlotPairingCoefficients("miR-1", "equil_c_nb", pdf.plot="5.Bi")
  # model_m1 <<- model
  # PlotOffsetCoefficients("miR-1", "equil_c_nb", pdf.plot="5.Bii")
  # PlotPairingMatrix("miR-1", "equil_c_nb", 1, model_values=TRUE,
  #                   makeglobalmodel=FALSE, corrected_kds=corrected_kds,
  #                   pdf.plot="5.Biii")
  # # R_miR1_model <<- R_mat
  # PlotPairingMatrix("miR-1", "equil_c_nb", 1, corrected_kds=corrected_kds,
  #                   pdf.plot="5.Biv")
  # R_miR1_data <<- R_mat
  # # C___________________________________________________________________________
  # PlotPairingCoefficients("miR-155", "equil_sc_nb", pdf.plot="5.Ci")
  # model_m155 <<- model
  # PlotOffsetCoefficients("miR-155", "equil_sc_nb", pdf.plot="5.Cii")
  # PlotPairingMatrix("miR-155", "equil_sc_nb", 1, model_values=TRUE,
  #                   makeglobalmodel=FALSE, corrected_kds=corrected_kds,
  #                   pdf.plot="5.Ciii")
  # # R_miR155_model <<- R_mat
  # PlotPairingMatrix("miR-155", "equil_sc_nb", 1, corrected_kds=corrected_kds,
  #                   pdf.plot="5.Civ")
  # R_miR155_data <<- R_mat
  # D___________________________________________________________________________
  PlotDeltaGMatrix("let-7a-21nt", pdf.plot="5.Di")
  # DeltaG_mat_let7 <<- DeltaG_mat_global
  PlotDeltaGMatrix("miR-1", pdf.plot="5.Dii")
  # DeltaG_mat_miR1 <<- DeltaG_mat_global
  PlotDeltaGMatrix("miR-155", key=TRUE, pdf.plot="5.Diii")
  # DeltaG_mat_miR155 <<- DeltaG_mat_global
  # E___________________________________________________________________________
  # PlotKdFoldChangeAgainstDeltaG("let-7a-21nt", "equil_c2_nb", pdf.plot="5.Ei")
  # PlotKdFoldChangeAgainstDeltaG("miR-1", "equil_c_nb", pdf.plot="5.Eii")
  # PlotKdFoldChangeAgainstDeltaG("miR-155", "equil_sc_nb", pdf.plot="5.Eiii")
  # # F___________________________________________________________________________
  # PlotMismatchCoefficients("let-7a-21nt", "equil_c2_nb", pdf.plot="5.Fi")
  # mismatches_let7 <<- mismatch_model_global
  # # mismatches_let7 <<- empirical_mismatch_global
  # PlotMismatchCoefficients("miR-1", "equil_c_nb", pdf.plot="5.Fii")
  # mismatches_miR1 <<- mismatch_model_global
  # # mismatches_miR1 <<- empirical_mismatch_global
  # PlotMismatchCoefficients("miR-155", "equil_sc_nb", pdf.plot="5.Fiii")
  # mismatches_miR155 <<- mismatch_model_global
  # # mismatches_miR155 <<- empirical_mismatch_global
  # # Within-figure analysis looking at the average contribution of each position
  # # of three prime pairing.
  # all_mismatches <- do.call("rbind", list(mismatches_let7$mm,
  #                                         mismatches_miR1$mm,
  #                                         mismatches_miR155$mm))
  # all_mismatches_global <<- all_mismatches
  # pos <- substr(rownames(all_mismatches), start=9, stop=9)
  # df_new <<- data.frame(all_mismatches, `pos`=as.numeric(pos))

  # average_pos_effect <<- aggregate(df_new$MLE, list(df_new$pos), mean,
  #                                 na.rm=TRUE)
  # print(average_pos_effect)
  # sd_pos_effect <<- aggregate(df_new$MLE, list(df_new$pos), sd, na.rm=TRUE)
  # # G___________________________________________________________________________
  # PlotAllProgrammedbraryMismatchCoefficientsAgainstKds(single_error=TRUE,
  #   pdf.plot="5.G"
  # )
  # # H___________________________________________________________________________
  # PlotAllRandomLibrarySeedCoefficientsAgainstKds(
  #   len_lim=c(4, 5), offset_lim=c(0, 10), pos_lim=c(9, 18), pos_5p_lim=TRUE,
  #   empirical=TRUE, sumseed=TRUE, pdf.plot="5.H"
  # )
  message("Done Fig. 5.")
}


# No Figure 6 function, since it is entirely a schematic and lifted r values.


MakeFigure7 <- function(model_values=FALSE, mm_and_bulge=TRUE) {
  message("Making Fig. 7.")
  # A__________________________________________________________________________
  # PlotBestThreePrimeSite("let-7a-21nt", len_lim=c(8, 11), justlegend=TRUE,
  #                        pdf.plot="7.Ai")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 11, 18,
  #                     pdf.plot="7.Aii")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb",  9, 17,
  #                     pdf.plot="7.Aiii")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb",  9, 18,
  #                     pdf.plot="7.Aiv")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 10, 20, titles=TRUE,
  #                     pdf.plot="7.Av")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 11, 18, bulge=TRUE,
  #                     pdf.plot="7.Avi")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb",  9, 17, bulge=TRUE,
  #                     pdf.plot="7.Avii")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb",  9, 18, bulge=TRUE,
  #                     pdf.plot="7.Aviii")
  # PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 10, 20, bulge=TRUE,
  #                     titles=TRUE, pdf.plot="7.Aix")
  # # B___________________________________________________________________________
  # PlotBestThreePrimeSite("miR-1", experiment="equil_c_nb", len_lim=c(8, 11),
  #                        justlegend=TRUE, pdf.plot="7.Bi")
  # PlotSiteMismatches("miR-1", "equil_c_nb", 11, 18, 
  #                    pdf.plot="7.Bii")
  # PlotSiteMismatches("miR-1", "equil_c_nb", 11, 19, 
  #                    pdf.plot="7.Biii")
  # PlotSiteMismatches("miR-1", "equil_c_nb", 11, 20, 
  #                    pdf.plot="7.Biv")
  # PlotSiteMismatches("miR-1", "equil_c_nb", 11, 21, titles=TRUE,
  #                     pdf.plot="7.Bv")
  # PlotSiteMismatches("miR-1", "equil_c_nb", 11, 18, bulge=TRUE,
  #                     pdf.plot="7.Bvi")
  # PlotSiteMismatches("miR-1", "equil_c_nb", 11, 19, bulge=TRUE,
  #                     pdf.plot="7.Bvii")
  # PlotSiteMismatches("miR-1", "equil_c_nb", 11, 20, bulge=TRUE,
  #                     pdf.plot="7.Bviii")
  # PlotSiteMismatches("miR-1", "equil_c_nb", 11, 21, bulge=TRUE, titles=TRUE,
  #                     pdf.plot="7.Bix")
  # # C___________________________________________________________________________
  # PlotBestThreePrimeSite("miR-155", experiment="equil_sc_nb", len_lim=c(8, 11),
  #                        justlegend=TRUE, pdf.plot="7.Ci")
  # PlotSiteMismatches("miR-155", "equil_sc_nb", 15, 22,
  #                     pdf.plot="7.Cii")
  # PlotSiteMismatches("miR-155", "equil_sc_nb", 13, 21,
  #                     pdf.plot="7.Ciii")
  # PlotSiteMismatches("miR-155", "equil_sc_nb", 13, 22,
  #                     pdf.plot="7.Civ")
  # PlotSiteMismatches("miR-155", "equil_sc_nb", 13, 23, titles=TRUE,
  #                     pdf.plot="7.Cv")
  # PlotSiteMismatches("miR-155", "equil_sc_nb", 15, 22, bulge=TRUE,
  #                     pdf.plot="7.Cvi")
  # PlotSiteMismatches("miR-155", "equil_sc_nb", 13, 21, bulge=TRUE,
  #                     pdf.plot="7.Cvii")
  # PlotSiteMismatches("miR-155", "equil_sc_nb", 13, 22, bulge=TRUE,
  #                     pdf.plot="7.Cviii")
  # PlotSiteMismatches("miR-155", "equil_sc_nb", 13, 23, bulge=TRUE, titles=TRUE,
  #                     pdf.plot="7.Cix")
  # PlotPairingMatrixKey(model=FALSE, height=1.9, short_margins=TRUE, pdf.plot="7.Cx")

  # # D___________________________________________________________________________
  # PlotAllPositionalMismatches("let-7a-21nt", 10, 
  #                             legend=TRUE, pdf.plot="7.Di")
  # PlotAllPositionalMismatches("miR-1", 10, 
  #                             pdf.plot="7.Dii")
  # PlotAllPositionalMismatches("miR-155", 10, 
  #                             x_axis=TRUE, pdf.plot="7.Diii")
  # # E___________________________________________________________________________
  # PlotAllPositionalMismatches("let-7a-21nt", 10, bulge=TRUE,
  #                              legend=TRUE,
  #                             pdf.plot="7.Ei")
  # PlotAllPositionalMismatches("miR-1", 10, bulge=TRUE,
  #                              pdf.plot="7.Eii")
  # PlotAllPositionalMismatches("miR-155", 10, bulge=TRUE,
  #                              x_axis=TRUE,
  #                             pdf.plot="7.Eiii")
  # # F___________________________________________________________________________
  # PlotTerminalMismatchAgainstBulgeDifferenceDist(10, 
  #                                                pdf.plot="7.F")
  # # G___________________________________________________________________________
  PlotAllAverageMismatchesAndDeltaGs("programmed", 10,
                                     
                                     pdf.plot="7.Gi")
  x_global_abs <- x_global
  y_global_abs <- y_global
  PlotAllAverageMismatchesAndDeltaGs("programmed", 10, frac_dG=TRUE,                                     
                                     pdf.plot="7.Gii")

  PlotAllAverageSeedMismatchesAndDeltaGs(pdf.plot="7.Hi")
  PlotAllAverageSeedMismatchesAndDeltaGs(frac_dG=TRUE, pdf.plot="7.Hii")
  x_global_frac <- x_global
  y_global_frac <- y_global
  print(x_global_frac)
  print(y_global_frac)
  print(cor(x_global_abs, y_global_abs)^2)
  print(cor(x_global_frac, y_global_frac)^2)
  message("Done Fig. 7.")
}


## SUPPLEMENTAL FIGURES ########################################################
MakeFigure2s1 <- function() {
  ## TODO, This function doesn't work for B when sitelist_1 == "randthrp_comp",
  ## i.e., it doesn't display any of the three-prime sites. It's also true that
  ## this figure function currently can't handle performing the correction with
  ## the sumseed flag set to TRUE.
  message("Making Fig. 2s1")
  # A___________________________________________________________________________
  PlotPairwiseThrPKds(
    "let-7a-21nt", "equil_c2_nb", 3, "progthrp_suppcomp", TRUE, FALSE,
    "let-7a-21nt", "equil_c_nb", 3, "progthrp_suppcomp", TRUE, FALSE,
    thrp_only=FALSE, pdf.plot="2s1.Ai")
  PlotPairwiseThrPKds(
    "let-7a-21nt", "equil_c2_nb", 3, "progthrp", TRUE, FALSE,
    "let-7a-21nt", "equil_c_nb", 3, "progthrp", TRUE, FALSE,
    thrp_only=FALSE, pdf.plot="2s1.Aii")
  # B___________________________________________________________________________
  PlotPairwiseThrPKds(
    "let-7a", "equilibrium", 3, "randthrp_suppcomp", TRUE, FALSE,
    "let-7a-21nt", "equil_c2_nb", 3, "progthrp_suppcomp", TRUE, FALSE,
    pdf.plot="2s1.Bi")
  PlotPairwiseThrPKds(
    "let-7a", "equilibrium", 3, "randthrp_suppcomp", TRUE, FALSE,
    "let-7a-21nt", "equil_c2_nb", 3, "progthrp_suppcomp", TRUE, FALSE,
    corrected_kds=TRUE, pdf.plot="2s1.Bii")
  PlotPairwiseThrPKds(
    "miR-1", "equilibrium", 3, "randthrp_suppcomp", FALSE, TRUE,
    "miR-1", "equil_c_nb", 3, "progthrp_suppcomp", TRUE, FALSE,
    pdf.plot="2s1.Biii")
  PlotPairwiseThrPKds(
    "miR-1", "equilibrium", 3, "randthrp_suppcomp", FALSE, TRUE,
    "miR-1", "equil_c_nb", 3, "progthrp_suppcomp", TRUE, FALSE,
    corrected_kds=TRUE, pdf.plot="2s1.Biv")
  PlotPairwiseThrPKds(
    "miR-155", "equilibrium", 3, "randthrp_suppcomp", TRUE, FALSE,
    "miR-155", "equil_sc_nb", 3, "progthrp_suppcomp", TRUE, FALSE,
    pdf.plot="2s1.Bv")
  PlotPairwiseThrPKds(
    "miR-155", "equilibrium", 3, "randthrp_suppcomp", TRUE, FALSE,
    "miR-155", "equil_sc_nb", 3, "progthrp_suppcomp", TRUE, FALSE,
    corrected_kds=TRUE, pdf.plot="2s1.Bvi")
  message("Done Fig. 2s1")
}

MakeFigure2s2 <- function() {
  message("Making Fig. 2s2")
  # A___________________________________________________________________________
  PlotBeckerPairingMatrix("let-7a", 0, "8mer", mirna_label=TRUE, sitelabel=TRUE,
                          get_global_data=TRUE, pdf.plot="2s2.Ai")
  PlotBeckerPairingMatrix("let-7a", 1, "8mer", pdf.plot="2s2.Aii")
  PlotBeckerPairingMatrix("let-7a", 2, "8mer", pdf.plot="2s2.Aiii")
  PlotBeckerPairingMatrix("let-7a", 3, "8mer", pdf.plot="2s2.Aiv")
  PlotBeckerPairingMatrix("let-7a", 4, "8mer", key=TRUE, pdf.plot="2s2.Av")
  # B___________________________________________________________________________
  PlotBeckerPairingMatrix("let-7a", 0, "7mer-m8", mirna_label=TRUE,
                          sitelabel=TRUE, get_global_data=TRUE,
                          pdf.plot="2s2.Bi")
  PlotBeckerPairingMatrix("let-7a", 1, "7mer-m8", pdf.plot="2s2.Bii")
  PlotBeckerPairingMatrix("let-7a", 2, "7mer-m8", pdf.plot="2s2.Biii")
  PlotBeckerPairingMatrix("let-7a", 3, "7mer-m8", pdf.plot="2s2.Biv")
  PlotBeckerPairingMatrix("let-7a", 4, "7mer-m8", pdf.plot="2s2.Bv")
  # C___________________________________________________________________________
  PlotBeckerPairingMatrix("let-7a", 0, "7mer-A1", mirna_label=TRUE,
                          sitelabel=TRUE, get_global_data=TRUE,
                          pdf.plot="2s2.Ci")
  PlotBeckerPairingMatrix("let-7a", 1, "7mer-A1", pdf.plot="2s2.Cii")
  PlotBeckerPairingMatrix("let-7a", 2, "7mer-A1", pdf.plot="2s2.Ciii")
  PlotBeckerPairingMatrix("let-7a", 3, "7mer-A1", pdf.plot="2s2.Civ")
  PlotBeckerPairingMatrix("let-7a", 4, "7mer-A1", pdf.plot="2s2.Cv")
  # D___________________________________________________________________________
  PlotBeckerPairingMatrix("let-7a", 0, "6mer", mirna_label=TRUE, sitelabel=TRUE,
                          get_global_data=TRUE, pdf.plot="2s2.Di")
  PlotBeckerPairingMatrix("let-7a", 1, "6mer", pdf.plot="2s2.Dii")
  PlotBeckerPairingMatrix("let-7a", 2, "6mer", pdf.plot="2s2.Diii")
  PlotBeckerPairingMatrix("let-7a", 3, "6mer", pdf.plot="2s2.Div")
  PlotBeckerPairingMatrix("let-7a", 4, "6mer", pdf.plot="2s2.Dv")

  # E___________________________________________________________________________
  # Determine which 8mer-mm sites have the max and min of the number of Kds with
  # offsets between 0 and +4.
  becker_df  <- MakeBeckerThreePDataFrame("let-7a")
  counts <- sapply(GetAll8merMmSites("let-7a"), function(site) {
    nrow(subset(becker_df, mm == site & offset >= 0 & offset <= 4))
  })
  max_site <- names(counts)[which.max(counts)]
  min_site <- names(counts)[which.min(counts)]
  PlotBeckerPairingMatrix("let-7a", 0, max_site, mirna_label=TRUE,
                          sitelabel=TRUE, get_global_data=TRUE,
                          pdf.plot="2s2.Ei")
  PlotBeckerPairingMatrix("let-7a", 1, max_site, pdf.plot="2s2.Eii")
  PlotBeckerPairingMatrix("let-7a", 2, max_site, pdf.plot="2s2.Eiii")
  PlotBeckerPairingMatrix("let-7a", 3, max_site, pdf.plot="2s2.Eiv")
  PlotBeckerPairingMatrix("let-7a", 4, max_site, pdf.plot="2s2.Ev")
  # F___________________________________________________________________________
  PlotBeckerPairingMatrix("let-7a", 0, min_site, mirna_label=TRUE,
                          sitelabel=TRUE, get_global_data=TRUE,
                          pdf.plot="2s2.Fi")
  PlotBeckerPairingMatrix("let-7a", 1, min_site, pdf.plot="2s2.Fii")
  PlotBeckerPairingMatrix("let-7a", 2, min_site, pdf.plot="2s2.Fiii")
  PlotBeckerPairingMatrix("let-7a", 3, min_site, pdf.plot="2s2.Fiv")
  PlotBeckerPairingMatrix("let-7a", 4, min_site, pdf.plot="2s2.Fv")

  message("Done Fig. 2s2.")
}

MakeFigure3s1 <- function() {
  message("Making Fig. 3s1.")
  # Calculate the p-values ahead of time (due to the need for all of them in 
  # order to perform the Benjamini-Hochberg method).
  PlotMeanReporterRepression(mirna="miR-1", pdf.plot="3s1.A")
  PlotMeanReporterRepression(mirna="miR-1", dual_site=TRUE,
                             pdf.plot="3s1.B")
  message("Done Fig. 3s1.")
}

MakeFigure3s2 <- function(experiment="twist_reporter_assay_3p_2_tp") {
  message("Making Fig. 3s2.")
  # Calculate the p-values ahead of time (due to the need for all of them in 
  # order to perform the Benjamini-Hochberg method).
  PlotMeanReporterRepression(experiment=experiment, all_lin41_UTR=TRUE,
                             pdf.plot="3s2.Ai")
  single_site_p_values <- output_p_values_global
  PlotMeanReporterRepression(experiment=experiment, dual_site=TRUE, 
                             all_lin41_UTR=TRUE, pdf.plot="3s2.Aii")
  dual_site_p_values <- output_p_values_global
  print(single_site_p_values)
  print(dual_site_p_values)
  full_df <- rbind(single_site_p_values, dual_site_p_values[-1, ])
  full_df <- cbind(full_df, p.adjust(full_df$p_vals, method="fdr"))
  colnames(full_df)[3] <- "p.adjust"
  full_df <- cbind(full_df, full_df$p.adjust < 0.05)
  print(full_df)
  message("Done Fig. 3s2.")
}


MakeFigure4s1 <- function() {
  message("Making Fig. 4s1.")
  # A___________________________________________________________________________
  PlotOneThreePrimeSite( 4, "miR-1", "equil_c_nb", pdf.plot="4s1.Ai")
  PlotOneThreePrimeSite( 5, "miR-1", "equil_c_nb", pdf.plot="4s1.Aii")
  PlotOneThreePrimeSite( 6, "miR-1", "equil_c_nb", mirna_label=TRUE,
                        pdf.plot="4s1.Aiii")
  PlotOneThreePrimeSite( 7, "miR-1", "equil_c_nb", pdf.plot="4s1.Aiv")
  PlotOneThreePrimeSite( 8, "miR-1", "equil_c_nb", pdf.plot="4s1.Av")
  PlotOneThreePrimeSite( 9, "miR-1", "equil_c_nb", pdf.plot="4s1.Avi")
  # B___________________________________________________________________________
  PlotOneThreePrimeSite( 4, "miR-155", "equil_sc_nb", pdf.plot="4s1.Bi")
  PlotOneThreePrimeSite( 5, "miR-155", "equil_sc_nb", pdf.plot="4s1.Bii")
  PlotOneThreePrimeSite( 6, "miR-155", "equil_sc_nb", mirna_label=TRUE,
                        pdf.plot="4s1.Biii")
  PlotOneThreePrimeSite( 7, "miR-155", "equil_sc_nb", pdf.plot="4s1.Biv")
  PlotOneThreePrimeSite( 8, "miR-155", "equil_sc_nb", pdf.plot="4s1.Bv")
  PlotOneThreePrimeSite( 9, "miR-155", "equil_sc_nb", pdf.plot="4s1.Bvi")

  message("Done Fig. 4s1.")
}


MakeFigure5s1 <- function(
  corrected_kds=TRUE, additive=FALSE, intercept=FALSE
) {
  message("Making fig. 5s1.")
  # A-C_________________________________________________________________________
  model_values <- TRUE
  label_offset <- TRUE
  mirnas <- c("let-7a-21nt", "miR-1", "miR-155")
  experiments <- c("equil_c2_nb", "equil_c_nb", "equil_sc_nb")
  letter <- "A"
  for (i_m in 1:3) {
    mirna <- mirnas[i_m]
    experiment <- experiments[i_m]
    for (i in 1:8) {
      offset <- i*2 - 4
      if (i == 1) mirna_label <- TRUE
      else        mirna_label <- FALSE
      if (i == 1) makeglobalmodel <- TRUE
      else        makeglobalmodel <- FALSE
      if (i == 8) key <- TRUE
      else         key <- FALSE
      SubfunctionCall(PlotPairingMatrix,
                      pdf.plot=sprintf("5s1.%s%s", letter, kRomanNumerals[i]))
    }
    letter <- kNextLetter[letter]
  }
  # D___________________________________________________________________________
  PlotPairingAndOffsetModelAgainstData("let-7a-21nt", "equil_c2_nb",
                                       intercept=intercept,
                                       log_plus_one=log_plus_one,
                                       pdf.plot="5s1.Di")
  num_global_l7 <- num_global
  PlotPairingAndOffsetModelAgainstData("miR-1", "equil_c_nb",
                                       intercept=intercept,
                                       log_plus_one=log_plus_one,
                                       pdf.plot="5s1.Dii")
  num_global_m1 <- num_global
  PlotPairingAndOffsetModelAgainstData("miR-155", "equil_sc_nb",
                                       intercept=intercept,
                                       log_plus_one=log_plus_one,
                                       pdf.plot="5s1.Diii")
  num_global_m155 <- num_global
  print(num_global_l7)
  print(num_global_m1)
  print(num_global_m155)
  message("Done fig. 5s1.")
}

MakeFigure5s2 <- function() {
  message("Making fig. 5s2.")
  # A___________________________________________________________________________
  PlotKdFoldChangeAgainstDeltaG("let-7a-21nt", "equil_c2_nb", lengthnorm=FALSE,
                                pdf.plot="5s2.Ai")
  mat_global_l7 <<- pairing_mat_global
  dG_global_l7 <<- dG_mat_global
  PlotKdFoldChangeAgainstDeltaG("miR-1", "equil_c_nb", lengthnorm=FALSE,
                                pdf.plot="5s2.Aii")
  mat_global_m1 <<- pairing_mat_global
  dG_global_m1 <<- dG_mat_global
  PlotKdFoldChangeAgainstDeltaG("miR-155", "equil_sc_nb", lengthnorm=FALSE,
                                pdf.plot="5s2.Aiii")
  mat_global_m155 <<- pairing_mat_global
  dG_global_m155 <<- dG_mat_global
  # B___________________________________________________________________________
  PlotKdFoldChangeAgainstDeltaG("let-7a-21nt", "equil_c2_nb", lengthnorm=FALSE,
                                pos_min=11, pdf.plot="5s2.Bi")
  PlotKdFoldChangeAgainstDeltaG("miR-1", "equil_c_nb", lengthnorm=FALSE,
                                pos_min=12, pdf.plot="5s2.Bii")
  PlotKdFoldChangeAgainstDeltaG("miR-155", "equil_sc_nb", lengthnorm=FALSE,
                                pos_max=20, pdf.plot="5s2.Biii")
  # C___________________________________________________________________________
  PlotKdFoldChangeAgainstDeltaG("let-7a-21nt", "equil_c2_nb", lengthnorm=FALSE,
                                pos_min=11, pos_inv=TRUE, pdf.plot="5s2.Ci")
  PlotKdFoldChangeAgainstDeltaG("miR-1", "equil_c_nb", lengthnorm=FALSE,
                                pos_min=12, pos_inv=TRUE, pdf.plot="5s2.Cii")
  PlotKdFoldChangeAgainstDeltaG("miR-155", "equil_sc_nb", lengthnorm=FALSE,
                                pos_max=20, pos_inv=TRUE, pdf.plot="5s2.Ciii")
  # D___________________________________________________________________________
  PlotKdFoldChangeAgainstDeltaG("let-7a-21nt", "equil_c2_nb",
                                model_values=FALSE, offset=4, pdf.plot="5s2.Di")
  PlotKdFoldChangeAgainstDeltaG("miR-1", "equil_c_nb", model_values=FALSE, 
                                offset=1, pdf.plot="5s2.Dii")
  PlotKdFoldChangeAgainstDeltaG("miR-155", "equil_sc_nb", model_values=FALSE,
                                offset=1, pdf.plot="5s2.Diii")
  # E___________________________________________________________________________
  PlotKdFoldChangeAgainstDeltaG("let-7a-21nt", "equil_c2_nb",
                                model_values=FALSE, offset=4, lengthnorm=FALSE,
                                pdf.plot="5s2.Ei")
  PlotKdFoldChangeAgainstDeltaG("miR-1", "equil_c_nb", model_values=FALSE, 
                                offset=1, lengthnorm=FALSE, pdf.plot="5s2.Eii")
  PlotKdFoldChangeAgainstDeltaG("miR-155", "equil_sc_nb", model_values=FALSE,
                                offset=1, lengthnorm=FALSE, pdf.plot="5s2.Eiii")
  # F___________________________________________________________________________
  PlotPairingOffsetAndMismatchModelAgainstData("let-7a-21nt", "equil_c2_nb",
                                               pdf.plot="5s2.Fi")
  num_global_l7 <- num_global
  PlotPairingOffsetAndMismatchModelAgainstData("miR-1", "equil_c_nb",
                                               pdf.plot="5s2.Fii")
  num_global_m1 <- num_global
  PlotPairingOffsetAndMismatchModelAgainstData("miR-155", "equil_sc_nb",
                                               pdf.plot="5s2.Fiii")
  num_global_m155 <- num_global
  print(num_global_l7)
  print(num_global_m1)
  print(num_global_m155)

  message("Done fig. 5s2.")
}

MakeFigure5s3 <- function(
  exponential=FALSE, intercept=FALSE, additive=FALSE, len_lim=c(4, 8),
  offset_lim=c(-4, 16), pos_lim=c(9, 23)
) {
  message("Making fig. 5s3.")
  # A-F_________________________________________________________________________
  sitelist <- "randthrp_suppcomp"
  mm <- FALSE
  corrected_kds <- FALSE
  # globalmodel <- FALSE
  letter <- "A"
  mirnas <- c(kMirnas[2], kMirnas[1], kMirnas[3:6])
  for (mirna in mirnas) {
    if (mirna == "miR-7-23nt") {
      experiment <- "equilibrium2_nb"
    } else {
      experiment <- "equilibrium"
    }
    SubfunctionCall(PlotPairingCoefficients, makeglobalmodel=TRUE,
                    pdf.plot=sprintf("5s3.%si", letter))
    SubfunctionCall(PlotOffsetCoefficients, makeglobalmodel=FALSE,
                    pdf.plot=sprintf("5s3.%sii", letter))
    # Makes the keys for the first row (against let-7a).
    if (mirna == "let-7a") {
      PlotPairingMatrixKey(pdf.plot="5s3.Aiii")
    }
    if (mirna == "lsy-6") {
      model_lsy6 <<- model
    } else if (mirna == "miR-124") {
      model_miR124 <<- model
    } else if (mirna == "miR-155") {
      model_miR155 <<- model
    }
    letter <- kNextLetter[letter]
  }
  # G___________________________________________________________________________
  SubfunctionCall(PlotProgrammedAndRandomPairingCoefficients,
    mirna_1="let-7a-21nt", experiment_1="equil_c2_nb",
    mirna_2="let-7a", experiment_2="equilibrium", pdf.plot="5s3.Gi"
  )
  SubfunctionCall(PlotProgrammedAndRandomPairingCoefficients,
    mirna_1="miR-1", experiment_1="equil_c_nb",
    mirna_2="miR-1", experiment_2="equilibrium", pdf.plot="5s3.Gii"
  )
  SubfunctionCall(PlotProgrammedAndRandomPairingCoefficients,
    mirna_1="miR-155", experiment_1="equil_sc_nb",
    mirna_2="miR-155", experiment_2="equilibrium", pdf.plot="5s3.Giii"
  )
  # H___________________________________________________________________________
  SubfunctionCall(PlotProgrammedAndRandomOffsetCoefficients,
    mirna_1="let-7a-21nt", experiment_1="equil_c2_nb",
    mirna_2="let-7a", experiment_2="equilibrium", pdf.plot="5s3.Hi"
  )
  SubfunctionCall(PlotProgrammedAndRandomOffsetCoefficients,
    mirna_1="miR-1", experiment_1="equil_c_nb",
    mirna_2="miR-1", experiment_2="equilibrium", pdf.plot="5s3.Hii"
  )
  SubfunctionCall(PlotProgrammedAndRandomOffsetCoefficients,
    mirna_1="miR-155", experiment_1="equil_sc_nb",
    mirna_2="miR-155", experiment_2="equilibrium", pdf.plot="5s3.Hiii"
  )

  message("Done fig. 5s3.")
}


MakeFigure5s4 <- function(height=3.7*9/9.75) {
  message("Making fig. 5s4.")
  mirna <- "let-7a"
  experiment <- "equilibrium"
  # A___________________________________________________________________________
  PlotBestThreePrimeSite(mirna, experiment, sitelist="randthrp_suppcomp",
                         corrected_kds=FALSE, len_lim=c(4, 8), height=height,
                         width=7.4, pdf.plot="5s4.A")
  # B___________________________________________________________________________
  PlotOneThreePrimeSite(mirna=mirna, experiment=experiment,
                        sitelist="randthrp_suppcomp", corrected_kds=FALSE,
                        kmer=4, height=height, pdf.plot="5s4.Bi")
  PlotOneThreePrimeSite(mirna=mirna, experiment=experiment,
                        sitelist="randthrp_suppcomp", corrected_kds=FALSE,
                        kmer=5, height=height, pdf.plot="5s4.Bii")
  mirna <- "miR-124"
  experiment <- "equilibrium"
  # C___________________________________________________________________________
  PlotBestThreePrimeSite(mirna, experiment, sitelist="randthrp_suppcomp",
                         corrected_kds=FALSE, len_lim=c(4, 8), height=height,
                         width=7.4, pdf.plot="5s4.C")
  # D___________________________________________________________________________
  PlotOneThreePrimeSite(mirna=mirna, experiment=experiment,
                        sitelist="randthrp_suppcomp", corrected_kds=FALSE,
                        kmer=4, height=height, pdf.plot="5s4.Di")
  PlotOneThreePrimeSite(mirna=mirna, experiment=experiment,
                        sitelist="randthrp_suppcomp", corrected_kds=FALSE,
                        kmer=5, height=height, pdf.plot="5s4.Dii")

  mirna <- "lsy-6"
  experiment <- "equilibrium"
  # E___________________________________________________________________________
  PlotBestThreePrimeSite(mirna, experiment, sitelist="randthrp_suppcomp",
                         corrected_kds=FALSE, len_lim=c(4, 8), height=height,
                         width=7.4, pdf.plot="5s4.E")
  # F___________________________________________________________________________
  PlotOneThreePrimeSite(mirna=mirna, experiment=experiment,
                        sitelist="randthrp_suppcomp", corrected_kds=FALSE,
                        kmer=4, height=height, pdf.plot="5s4.Fi")
  PlotOneThreePrimeSite(mirna=mirna, experiment=experiment,
                        sitelist="randthrp_suppcomp", corrected_kds=FALSE,
                        kmer=5, height=height, pdf.plot="5s4.Fii")

  mirna <- "miR-7-23nt"
  experiment <- "equilibrium2_nb"
  # G___________________________________________________________________________
  PlotBestThreePrimeSite(mirna, experiment, sitelist="randthrp_suppcomp",
                         corrected_kds=FALSE, len_lim=c(4, 8), height=height,
                         width=7.4, pdf.plot="5s4.G")
  # H___________________________________________________________________________
  PlotOneThreePrimeSite(mirna=mirna, experiment=experiment,
                        sitelist="randthrp_suppcomp", corrected_kds=FALSE,
                        kmer=4, height=height, pdf.plot="5s4.Hi")
  PlotOneThreePrimeSite(mirna=mirna, experiment=experiment,
                        sitelist="randthrp_suppcomp", corrected_kds=FALSE,
                        kmer=5, height=height, pdf.plot="5s4.Hii")

  message("Making fig. 5s4.")
}

MakeFigure5s5 <- function(
  len_lim=c(4, 10), offset_lim=c(-4, 16), pos_lim=c(9, 23)
) {
  message("Making fig. 5s5")
  mirna <- "let-7a-21nt"
  experiment <- "equil_c2_nb"
  len_mir <- nchar(kMirnaSeqs[mirna])
  sitelist <- "progthrp"
  corrected_kds <- TRUE
  kd_fc <- TRUE
  supp_base <- FALSE
  letter_use <- "A"
  model_values <- FALSE
  # Outer loop over the lengths.
  print_mirna <- TRUE
  for (len in seq(len_lim[1], len_lim[2])) {
    max_pos <- len_mir - 8 - len + 1
    if (len == len_lim[2])  show_offsets <- TRUE
    else                    show_offsets <- FALSE
    # Inner loop over the 5-prime positions of pairing within the miRNA.
    for (i in 1:max_pos) {
      if (i == 1) {
        show_site_names <- TRUE
      } else {
        show_site_names <- FALSE
      }
      SubfunctionCall(PlotSiteBySeedMismatchesAndOffset, pos5p=i + 8,
                      pdf.plot=sprintf("5s5.%s%s", letter_use,
                                       kRomanNumerals[i]))
      print_mirna <- FALSE
    }
    letter_use <- kNextLetter[letter_use]
  }
  model_values <- TRUE
  makeglobalmodel <- TRUE
  for (len in seq(len_lim[2], len_lim[1])) {
    max_pos <- len_mir - 8 - len + 1
    if (len == 4) show_offsets <- TRUE
    else          show_offsets <- FALSE
    # Inner loop over the 5-prime positions of pairing within the miRNA.
    for (i in max_pos:1) {
      if (i == max_pos) {
        show_site_names <- TRUE
      } else {
        show_site_names <- FALSE
      }
      SubfunctionCall(PlotSiteBySeedMismatchesAndOffset, pos5p=i + 8,
                      pdf.plot=sprintf("5s5.%s%s", letter_use,
                                       kRomanNumerals[max_pos - i + 1]))
      makeglobalmodel <- FALSE
    }
  letter_use <- kNextLetter[letter_use]
  }
  PlotPairingMatrixKey(horizontal=TRUE, pdf.plot="5s5.O")
  message("Done fig. 5s5.")  
}

MakeFigure5s6 <- function(
  len_lim=c(4, 10), offset_lim=c(-4, 16), pos_lim=c(9, 23)
) {
  message("Making fig. 5s6.")
  mirna <- "miR-1"
  experiment <- "equil_c_nb"
  len_mir <- nchar(kMirnaSeqs[mirna])
  sitelist <- "progthrp"
  corrected_kds <- TRUE
  kd_fc <- TRUE
  supp_base <- FALSE
  letter_use <- "A"
  print_mirna <- TRUE
  model_values <- FALSE

  # Outer loop over the lengths.
  for (len in seq(len_lim[1], len_lim[2])) {
    max_pos <- len_mir - 8 - len + 1
    if (len == len_lim[2])  show_offsets <- TRUE
    else                    show_offsets <- FALSE
    # Inner loop over the 5-prime positions of pairing within the miRNA.
    for (i in 1:max_pos) {
      if (i == 1) show_site_names <- TRUE
      else        show_site_names <- FALSE
      if (i == max_pos) show_offsets <- FALSE
      SubfunctionCall(PlotSiteBySeedMismatchesAndOffset, pos5p=i + 8,
                      pdf.plot=sprintf("5s6.%s%s", letter_use,
                                       kRomanNumerals[i]))
      print_mirna <- FALSE
    }
    letter_use <- kNextLetter[letter_use]
  }
  model_values <- TRUE
  makeglobalmodel <- TRUE
  for (len in seq(len_lim[2], len_lim[1])) {
    max_pos <- len_mir - 8 - len + 1
    if (len == 4) show_offsets <- TRUE
    else          show_offsets <- FALSE
    # Inner loop over the 5-prime positions of pairing within the miRNA.
    for (i in max_pos:1) {
      if (i == max_pos) show_site_names <- TRUE
      else              show_site_names <- FALSE
      SubfunctionCall(PlotSiteBySeedMismatchesAndOffset, pos5p=i + 8,
                      pdf.plot=sprintf("5s6.%s%s", letter_use,
                                       kRomanNumerals[max_pos - i + 1]))
      makeglobalmodel <- FALSE
    }
  letter_use <- kNextLetter[letter_use]
  }
  PlotPairingMatrixKey(pdf.plot="5s6.O")
  message("Done fig. 5s6.")  
}

MakeFigure5s7 <- function(
  len_lim=c(4, 10), offset_lim=c(-4, 16), pos_lim=c(9, 23)
) {
  message("Making fig. 5s7.")
  mirna <- "miR-155"
  experiment <- "equil_sc_nb"
  len_mir <- nchar(kMirnaSeqs[mirna])
  sitelist <- "progthrp"
  corrected_kds <- TRUE
  kd_fc <- TRUE
  supp_base <- FALSE
  letter_use <- "A"
  print_mirna <- TRUE
  model_values <- FALSE
  # Outer loop over the lengths.
  for (len in seq(len_lim[1], len_lim[2])) {
    max_pos <- len_mir - 8 - len + 1
    if (len == len_lim[2])  show_offsets <- TRUE
    else                    show_offsets <- FALSE
    # Inner loop over the 5-prime positions of pairing within the miRNA.
    for (i in 1:max_pos) {
      if (i == 1) show_site_names <- TRUE
      else        show_site_names <- FALSE
      if (i == max_pos) show_offsets <- FALSE
      SubfunctionCall(PlotSiteBySeedMismatchesAndOffset, pos5p=i + 8,
                      pdf.plot=sprintf("5s7.%s%s", letter_use,
                                       kRomanNumerals[i]))
      print_mirna <- FALSE
    }
    letter_use <- kNextLetter[letter_use]
  }
  model_values <- TRUE
  makeglobalmodel <- TRUE
  for (len in 10:4) {
    max_pos <- len_mir - 8 - len + 1
    if (len == 4) show_offsets <- TRUE
    else          show_offsets <- FALSE
    # Inner loop over the 5-prime positions of pairing within the miRNA.
    for (i in max_pos:1) {
      if (i == max_pos) show_site_names <- TRUE
      else              show_site_names <- FALSE
      SubfunctionCall(PlotSiteBySeedMismatchesAndOffset, pos5p=i + 8,
                      pdf.plot=sprintf("5s7.%s%s", letter_use,
                                       kRomanNumerals[max_pos - i + 1]))
      makeglobalmodel <- FALSE
    }
  letter_use <- kNextLetter[letter_use]
  }
  PlotPairingMatrixKey(pdf.plot="5s7.O")
  message("Done fig. 5s7.")  
}


MakeFigure5s8 <- function() {
  # Make the random library coefficients
  message("Making fig. 5s8.")
  experiment <- "equilibrium"
  sitelist <- "randthrp_comp"
  corrected_kds <- FALSE
  kd_fc <- TRUE
  supp_base <- TRUE
  model <- FALSE
  offset_lim <- c(0, 10)
  len_lim <- c(4, 5)
  sumseed <- TRUE
  highlight_column <- TRUE
  model_values <- FALSE
  letter_use <- "A"
  mirnas <- c("let-7a", "miR-1", "miR-155", "miR-124", "lsy-6", "miR-7-23nt")
  # Outer loop over the miRNAS.
  for (mirna in mirnas) {
    i_use <- 1
    if (mirna == "miR-1") {
      experiment <- "equilibrium"
      combined <- FALSE
      buffer <- TRUE
    } else if (mirna == "miR-7-23nt") {
      experiment <- "equilibrium2_nb"
      combined <- FALSE
      buffer <- FALSE
    } else {
      experiment <- "equilibrium"
      combined <- TRUE
      buffer <- FALSE
    }
    table_use_global <<- SubfunctionCall(
      GetEmpiricalBestSeedMismatchCoefficients, offset_lim=c(0, 10),
      len_lim=c(4, 5), pos_lim=c(9, 18), pos_5p_lim=TRUE
    )
    # Middle loop over the two lengths.
    for (len in 4:5) {
      if (len == 4) {
        show_offsets <- FALSE
      } else {
        show_offsets <- TRUE
      }
      # Inner loop over the 5-prime positions of pairing within the miRNA.
      for (i in 1:10) {
        if (i == 1) {
          print_mirna <- TRUE
          show_site_names <- TRUE
        } else {
          print_mirna <- FALSE
          show_site_names <- FALSE
        }
        SubfunctionCall(PlotSiteBySeedMismatchesAndOffset, pos5p=i + 8,
                        pdf.plot=sprintf("5s8.%s%s", letter_use,
                                         kRomanNumerals[i_use]))
        i_use <- i_use + 1
      }
    }
    letter_use <- kNextLetter[letter_use]
  }
  PlotPairingMatrixKey(horizontal=TRUE, pdf.plot="5s8.G")
  message("Done fig. 5s8.")
}

MakeFigure6s1 <- function(corrected_kds=TRUE) {
  message("Making fig. 6s1.")
  # A___________________________________________________________________________
  # Illustrator schematic.
  # B___________________________________________________________________________
  PlotPairingCoefficients("let-7a_miR-155", "equil_c_nb", pdf.plot="6s1.Bi")
  # model_coefficients_used_global <<- 10^model$pairing$MLE
  model_max_value <<- max(10^model$pairing$MLE, na.rm=TRUE)
  # model_save <- model$pairing$MLE  
  PlotOffsetCoefficients("let-7a_miR-155", "equil_c_nb", makeglobalmodel=FALSE,
                         pdf.plot="6s1.Bii")

  # Note: shouldn't need to calculate the max from this panel because it should
  # be identical to the pairing coefficients (because the fit is happening at
  # the optimum offset coefficient, which is definitionally equal to 1).
  PlotPairingMatrix("let-7a_miR-155", "equil_c_nb", 3, model_values=TRUE,
                    makeglobalmodel=FALSE, corrected_kds=corrected_kds,
                    pdf.plot="6s1.Biii")
  model_used_global <<- 10^R_mat
  PlotPairingMatrix("let-7a_miR-155", "equil_c_nb", 3,
                    corrected_kds=corrected_kds, pdf.plot="6s1.Biv")
  # data_used_global <<- 10^R_mat
  # Get the maximum value of the first three plots.
  data_max_value <<- max(10^R_mat, na.rm=TRUE)
  max_value <- max(model_max_value, data_max_value)
  PlotPairingMatrixKey(alt_top=TRUE, max_value=max_value, pdf.plot="6s1.Bv")
  PlotPairingMatrixKey(model=FALSE, alt_top=TRUE, max_value=max_value,
                       pdf.plot="6s1.Bvi")
  # C___________________________________________________________________________
  PlotPairingCoefficients("miR-155_let-7a", "equil_c_nb", pdf.plot="6s1.Ci")
  model_m155l7 <<- model
  PlotOffsetCoefficients("miR-155_let-7a", "equil_c_nb", makeglobalmodel=FALSE,
                         pdf.plot="6s1.Cii")
  PlotPairingMatrix("miR-155_let-7a", "equil_c_nb", 4, model_values=TRUE,
                    makeglobalmodel=FALSE, corrected_kds=corrected_kds,
                    pdf.plot="6s1.Ciii")
  PlotPairingMatrix("miR-155_let-7a", "equil_c_nb", 4,
                    corrected_kds=corrected_kds, pdf.plot="6s1.Civ")
  # D___________________________________________________________________________
  PlotChimeraModelCoefficientsScatter(
    "miR-155", "equil_sc_nb", mirna_label=TRUE,
    legend=TRUE, pdf.plot="6s1.Di"
  )
  # linmodel_global1 <<- linmodel_global
  PlotChimeraModelCoefficientsScatter(
    "miR-155", "equil_sc_nb", makeglobalmodel=FALSE, pairing_coefs=FALSE,
    legend=TRUE, pdf.plot="6s1.Dii"
  )
  # linmodel_global2 <<- linmodel_global
  # E___________________________________________________________________________
  PlotChimeraModelCoefficientsScatter(
    "let-7a-21nt", "equil_c2_nb", mirna_label=TRUE, pdf.plot="6s1.Ei"
  )
  # linmodel_global3 <<- linmodel_global
  PlotChimeraModelCoefficientsScatter(
    "let-7a-21nt", "equil_c2_nb", makeglobalmodel=FALSE, pairing_coefs=FALSE,
    pdf.plot="6s1.Eii"
  )
  # linmodel_global4 <<- linmodel_global
  # F___________________________________________________________________________
  PlotChimeraModelCoefficientsScatter(
    "miR-155", "equil_sc_nb", seed_compare=TRUE, mirna_label=TRUE,
    pdf.plot="6s1.Fi"
  )
  PlotChimeraModelCoefficientsScatter(
    "miR-155", "equil_sc_nb", makeglobalmodel=FALSE, seed_compare=TRUE,
    pairing_coefs=FALSE, pdf.plot="6s1.Fii"
  )
  # G___________________________________________________________________________
  PlotChimeraModelCoefficientsScatter(
    "let-7a-21nt", "equil_c2_nb", seed_compare=TRUE, mirna_label=TRUE,
    pdf.plot="6s1.Gi"
  )
  PlotChimeraModelCoefficientsScatter(
    "let-7a-21nt", "equil_c2_nb", makeglobalmodel=FALSE, seed_compare=TRUE,
    pairing_coefs=FALSE, pdf.plot="6s1.Gii"
  )
  # H___________________________________________________________________________
  PlotMismatchCoefficients("let-7a_miR-155", "equil_c_nb", height=4, width=4,
                           pdf.plot="6s1.Hi")
  PlotMismatchCoefficients("miR-155_let-7a", "equil_c_nb", height=4, width=4,
                           pdf.plot="6s1.Hii")

  # I___________________________________________________________________________
  PlotChimeraMismatchModelCoefficientsScatter(
    "let-7a-21nt", "equil_c2_nb", pdf.plot="6s1.Ii", legend=TRUE
  )
  # # x_1 <<- x
  # # y_1 <<- y
  # # ind <- which.max((x_1 - y_1)^2)
  # # cor(x_1, y_1)^2
  # # cor(x_1[-ind], y_1[-ind])^2
  PlotChimeraMismatchModelCoefficientsScatter(
    "miR-155", "equil_sc_nb", pdf.plot="6s1.Iii"
  )
  # x_2 <<- x
  # y_2 <<- y
  # ind <- which.max((x_2 - y_2)^2)
  # print(cor(x_2, y_2)^2)
  # print(cor(x_2[-ind], y_2[-ind])^2)
  message("Done fig. 6s1.")
}

MakeFigure6s2 <- function(corrected_kds=TRUE) {
  message("Making fig. 6s2.")
  # A___________________________________________________________________________
  # Illustrator schematic.
  # B___________________________________________________________________________
  PlotPairingCoefficients("let-7a_minus1", "equil_c_nb", pdf.plot="6s2.Bi")
  PlotOffsetCoefficients("let-7a_minus1", "equil_c_nb", makeglobalmodel=FALSE,
                         pdf.plot="6s2.Bii")
  PlotPairingMatrix("let-7a_minus1", "equil_c_nb", 3,
                    corrected_kds=corrected_kds, model_values=TRUE,
                    makeglobalmodel=FALSE, pdf.plot="6s2.Biii")
  PlotPairingMatrix("let-7a_minus1", "equil_c_nb", 3,
                    corrected_kds=corrected_kds, pdf.plot="6s2.Biv")
  PlotPairingMatrixKey(height=1.8, short_margins=TRUE, pdf.plot="6s2.Bv")
  PlotPairingMatrixKey(height=1.8, model=FALSE, short_margins=TRUE,
                       pdf.plot="6s2.Bvi")
  # C___________________________________________________________________________
  PlotPairingCoefficients("let-7a-21nt", "equil_c2_nb", pdf.plot="6s2.Ci")
  PlotOffsetCoefficients("let-7a-21nt", "equil_c2_nb", makeglobalmodel=TRUE,
                         pdf.plot="6s2.Cii")
  PlotPairingMatrix("let-7a-21nt", "equil_c2_nb", 4, model_values=TRUE,
                    makeglobalmodel=FALSE, corrected_kds=corrected_kds,
                    pdf.plot="6s2.Ciii")
  PlotPairingMatrix("let-7a-21nt", "equil_c2_nb", 4,
                    corrected_kds=corrected_kds, pdf.plot="6s2.Civ")
  # D___________________________________________________________________________
  PlotPairingCoefficients("let-7a_plus1", "equil_c_nb", pdf.plot="6s2.Di")
  PlotOffsetCoefficients("let-7a_plus1", "equil_c_nb", makeglobalmodel=TRUE,
                         pdf.plot="6s2.Dii")
  PlotPairingMatrix("let-7a_plus1", "equil_c_nb", 2, model_values=TRUE,
                    makeglobalmodel=FALSE, corrected_kds=corrected_kds,
                    pdf.plot="6s2.Diii")
  PlotPairingMatrix("let-7a_plus1", "equil_c_nb", 2,
                    corrected_kds=corrected_kds, pdf.plot="6s2.Div")
  # E___________________________________________________________________________
  PlotOffsetCrossCorrelation(height=3.5, width=3, pdf.plot="6s2.E")
  # F___________________________________________________________________________
  PlotMismatchCoefficients("let-7a_minus1", "equil_c_nb",
                           corrected_kds=corrected_kds, pdf.plot="6s2.Fi")
  mm_l7m1 <<- mismatch_model_global
  PlotMismatchCoefficients("let-7a-21nt", "equil_c2_nb", ylab=FALSE,
                           pdf.plot="6s2.Fii")
  mm_l7 <<- mismatch_model_global
  PlotMismatchCoefficients("let-7a_plus1", "equil_c_nb", ylab=FALSE,
                           pdf.plot="6s2.Fiii")
  mm_l7p1 <<- mismatch_model_global

  message("Done fig. 6s2.")
}


MakeFigure7s1 <- function(model_values=FALSE, mm_and_bulge=TRUE) {
  message("Making fig. 7s1.")
  # A___________________________________________________________________________
  PlotBestThreePrimeSite("let-7a-21nt", len_lim=c(4, 7), justlegend=TRUE,
                         pdf.plot="7s1.Ai")
  PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 11, 14, pdf.plot="7s1.Aii")
  PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 11, 15, pdf.plot="7s1.Aiii")
  PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 11, 16, pdf.plot="7s1.Aiv")
  PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 11, 17, titles=TRUE,
                     pdf.plot="7s1.Av")
  PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 11, 17, bulge=TRUE,
                     titles=TRUE, pdf.plot="7s1.Aix")
  PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 11, 16, bulge=TRUE,
                     pdf.plot="7s1.Aviii")
  PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 11, 15, bulge=TRUE,
                     pdf.plot="7s1.Avii")
  PlotSiteMismatches("let-7a-21nt", "equil_c2_nb", 11, 14, bulge=TRUE,
                     pdf.plot="7s1.Avi")
  # B___________________________________________________________________________
  PlotBestThreePrimeSite("miR-1", experiment="equil_c_nb", len_lim=c(4, 7),
                         justlegend=TRUE, pdf.plot="7s1.Bi")
  PlotSiteMismatches("miR-1", "equil_c_nb", 12, 18, titles=TRUE,
                     pdf.plot="7s1.Bv")
  PlotSiteMismatches("miR-1", "equil_c_nb", 12, 17, pdf.plot="7s1.Biv")
  PlotSiteMismatches("miR-1", "equil_c_nb", 12, 16, pdf.plot="7s1.Biii")
  PlotSiteMismatches("miR-1", "equil_c_nb", 12, 15, pdf.plot="7s1.Bii")
  PlotSiteMismatches("miR-1", "equil_c_nb", 12, 18, bulge=TRUE, titles=TRUE,
                     pdf.plot="7s1.Bix")
  PlotSiteMismatches("miR-1", "equil_c_nb", 12, 17, bulge=TRUE,
                     pdf.plot="7s1.Bviii")
  PlotSiteMismatches("miR-1", "equil_c_nb", 12, 16, bulge=TRUE,
                     pdf.plot="7s1.Bvii")
  PlotSiteMismatches("miR-1", "equil_c_nb", 12, 15, bulge=TRUE,
                     pdf.plot="7s1.Bvi")
  # C________________________________________________________________________
  PlotBestThreePrimeSite("miR-155", experiment="equil_sc_nb", len_lim=c(4, 7),
                         justlegend=TRUE, pdf.plot="7s1.Ci")
  PlotSiteMismatches("miR-155", "equil_sc_nb", 15, 21, titles=TRUE,
                     pdf.plot="7s1.Cv")
  PlotSiteMismatches("miR-155", "equil_sc_nb", 13, 18, pdf.plot="7s1.Civ")
  PlotSiteMismatches("miR-155", "equil_sc_nb", 13, 17, pdf.plot="7s1.Ciii")
  PlotSiteMismatches("miR-155", "equil_sc_nb", 13, 16, pdf.plot="7s1.Cii")
  PlotSiteMismatches("miR-155", "equil_sc_nb", 15, 21, bulge=TRUE, titles=TRUE,
                     pdf.plot="7s1.Cix")
  PlotSiteMismatches("miR-155", "equil_sc_nb", 13, 18, bulge=TRUE,
                     pdf.plot="7s1.Cviii")
  PlotSiteMismatches("miR-155", "equil_sc_nb", 13, 17, bulge=TRUE,
                     pdf.plot="7s1.Cvii")
  PlotSiteMismatches("miR-155", "equil_sc_nb", 13, 16, bulge=TRUE,
                     pdf.plot="7s1.Cvi")
  PlotPairingMatrixKey(model=FALSE, height=1.9, short_margins=TRUE,
                       pdf.plot="7s1.Cx")
  # D___________________________________________________________________________
  PlotSiteMismatchesBar("let-7a-21nt", "equil_c2_nb", 11, 19,
                        offset_lim=c(1, 1), win_average=1, pdf.plot="7s1.D")

  message("Done fig. 7s1.")
}

MakeFigure7s2 <- function(model_values=FALSE, mm_and_bulge=TRUE) {
  message("Making fig. 7s2.")
  # A___________________________________________________________________________
  PlotBestThreePrimeSite("let-7a", experiment="equilibrium", len_lim=c(4, 8),
                         sitelist="randthrp_suppcomp", corrected_kds=FALSE,
                         justlegend=TRUE, pdf.plot="7s2.Ai")
  PlotSiteMismatches("let-7a", "equilibrium", corrected_kds=FALSE, 10, 17,
                     titles=TRUE, pdf.plot="7s2.Avi")
  PlotSiteMismatches("let-7a", "equilibrium", corrected_kds=FALSE, 11, 17,
                     pdf.plot="7s2.Av")
  PlotSiteMismatches("let-7a", "equilibrium", corrected_kds=FALSE, 11, 16,
                     pdf.plot="7s2.Aiv")
  PlotSiteMismatches("let-7a", "equilibrium", corrected_kds=FALSE, 11, 15,
                     pdf.plot="7s2.Aiii")
  PlotSiteMismatches("let-7a", "equilibrium", corrected_kds=FALSE, 12, 15,
                     pdf.plot="7s2.Aii")
  PlotSiteMismatches("let-7a", "equilibrium", corrected_kds=FALSE, 10, 17,
                     bulge=TRUE, titles=TRUE, pdf.plot="7s2.Axi")
  PlotSiteMismatches("let-7a", "equilibrium", corrected_kds=FALSE, 11, 17,
                     bulge=TRUE, pdf.plot="7s2.Ax")
  PlotSiteMismatches("let-7a", "equilibrium", corrected_kds=FALSE, 11, 16,
                     bulge=TRUE, pdf.plot="7s2.Aix")
  PlotSiteMismatches("let-7a", "equilibrium", corrected_kds=FALSE, 11, 15,
                     bulge=TRUE, pdf.plot="7s2.Aviii")
  PlotSiteMismatches("let-7a", "equilibrium", corrected_kds=FALSE, 12, 15,
                     bulge=TRUE, pdf.plot="7s2.Avii")
  # B___________________________________________________________________________
  PlotBestThreePrimeSite("miR-1", experiment="equilibrium", len_lim=c(4, 8),
                         sitelist="randthrp_suppcomp", corrected_kds=FALSE,
                         justlegend=TRUE, pdf.plot="7s2.Bi")
  PlotSiteMismatches("miR-1", "equilibrium", corrected_kds=FALSE, 12, 19,
                     titles=TRUE, pdf.plot="7s2.Bvi")
  PlotSiteMismatches("miR-1", "equilibrium", corrected_kds=FALSE, 12, 18,
                     pdf.plot="7s2.Bv")
  PlotSiteMismatches("miR-1", "equilibrium", corrected_kds=FALSE, 13, 18,
                     pdf.plot="7s2.Biv")
  PlotSiteMismatches("miR-1", "equilibrium", corrected_kds=FALSE, 12, 16,
                     pdf.plot="7s2.Biii")
  PlotSiteMismatches("miR-1", "equilibrium", corrected_kds=FALSE, 12, 15,
                     pdf.plot="7s2.Bii")
  PlotSiteMismatches("miR-1", "equilibrium", corrected_kds=FALSE, 12, 19,
                     bulge=TRUE, titles=TRUE, pdf.plot="7s2.Bxi")
  PlotSiteMismatches("miR-1", "equilibrium", corrected_kds=FALSE, 12, 18,
                     bulge=TRUE, pdf.plot="7s2.Bx")
  PlotSiteMismatches("miR-1", "equilibrium", corrected_kds=FALSE, 13, 18,
                     bulge=TRUE, pdf.plot="7s2.Bix")
  PlotSiteMismatches("miR-1", "equilibrium", corrected_kds=FALSE, 12, 16,
                     bulge=TRUE, pdf.plot="7s2.Bviii")
  PlotSiteMismatches("miR-1", "equilibrium", corrected_kds=FALSE, 12, 15,
                     bulge=TRUE, pdf.plot="7s2.Bvii")
  # C___________________________________________________________________________
  PlotBestThreePrimeSite("miR-155", experiment="equilibrium", len_lim=c(4, 8),
                         sitelist="randthrp_suppcomp", corrected_kds=FALSE,
                         justlegend=TRUE, pdf.plot="7s2.Ci")
  PlotSiteMismatches("miR-155", "equilibrium", corrected_kds=FALSE, 15, 22,
                     titles=TRUE, pdf.plot="7s2.Cvi")
  PlotSiteMismatches("miR-155", "equilibrium", corrected_kds=FALSE, 15, 21,
                     pdf.plot="7s2.Cv")
  PlotSiteMismatches("miR-155", "equilibrium", corrected_kds=FALSE, 13, 18,
                     pdf.plot="7s2.Civ")
  PlotSiteMismatches("miR-155", "equilibrium", corrected_kds=FALSE, 13, 17,
                     pdf.plot="7s2.Ciii")
  PlotSiteMismatches("miR-155", "equilibrium", corrected_kds=FALSE, 13, 16,
                     pdf.plot="7s2.Cii")
  PlotSiteMismatches("miR-155", "equilibrium", corrected_kds=FALSE, 15, 22,
                     bulge=TRUE, titles=TRUE, pdf.plot="7s2.Cxi")
  PlotSiteMismatches("miR-155", "equilibrium", corrected_kds=FALSE, 15, 21,
                     bulge=TRUE, pdf.plot="7s2.Cx")
  PlotSiteMismatches("miR-155", "equilibrium", corrected_kds=FALSE, 13, 18,
                     bulge=TRUE, pdf.plot="7s2.Cix")
  PlotSiteMismatches("miR-155", "equilibrium", corrected_kds=FALSE, 13, 17,
                     bulge=TRUE, pdf.plot="7s2.Cviii")
  PlotSiteMismatches("miR-155", "equilibrium", corrected_kds=FALSE, 13, 16,
                     bulge=TRUE, pdf.plot="7s2.Cvii")
  PlotPairingMatrixKey(model=FALSE, height=1.9, short_margins=TRUE,
                       pdf.plot="7s2.Cxii")

  # D___________________________________________________________________________
  PlotBestThreePrimeSite("miR-124", experiment="equilibrium", len_lim=c(4, 8),
                         sitelist="randthrp_suppcomp", corrected_kds=FALSE,
                         justlegend=TRUE, pdf.plot="7s2.Di")
  PlotSiteMismatches("miR-124", "equilibrium", corrected_kds=FALSE, 11, 18,
                     titles=TRUE, pdf.plot="7s2.Dvi")  
  PlotSiteMismatches("miR-124", "equilibrium", corrected_kds=FALSE, 11, 17,
                     pdf.plot="7s2.Dv")
  PlotSiteMismatches("miR-124", "equilibrium", corrected_kds=FALSE, 11, 16,
                     pdf.plot="7s2.Div")
  PlotSiteMismatches("miR-124", "equilibrium", corrected_kds=FALSE, 11, 15,
                     pdf.plot="7s2.Diii")
  PlotSiteMismatches("miR-124", "equilibrium", corrected_kds=FALSE, 11, 14,
                     pdf.plot="7s2.Dii")
  PlotSiteMismatches("miR-124", "equilibrium", corrected_kds=FALSE, 11, 18,
                     titles=TRUE, bulge=TRUE, pdf.plot="7s2.Dxi")
  PlotSiteMismatches("miR-124", "equilibrium", corrected_kds=FALSE, 11, 17,
                     bulge=TRUE, pdf.plot="7s2.Dx")
  PlotSiteMismatches("miR-124", "equilibrium", corrected_kds=FALSE, 11, 16,
                     bulge=TRUE, pdf.plot="7s2.Dix")
  PlotSiteMismatches("miR-124", "equilibrium", corrected_kds=FALSE, 11, 15,
                     bulge=TRUE, pdf.plot="7s2.Dviii")
  PlotSiteMismatches("miR-124", "equilibrium", corrected_kds=FALSE, 11, 14,
                     bulge=TRUE, pdf.plot="7s2.Dvii")
  # E___________________________________________________________________________
  PlotBestThreePrimeSite("lsy-6", experiment="equilibrium", len_lim=c(4, 8),
                         sitelist="randthrp_suppcomp", corrected_kds=FALSE,
                         justlegend=TRUE, pdf.plot="7s2.Ei")
  PlotSiteMismatches("lsy-6", "equilibrium", corrected_kds=FALSE, 9, 16,
                     titles=TRUE, pdf.plot="7s2.Evi")
  PlotSiteMismatches("lsy-6", "equilibrium", corrected_kds=FALSE, 10, 16,
                     pdf.plot="7s2.Ev")
  PlotSiteMismatches("lsy-6", "equilibrium", corrected_kds=FALSE, 12, 17,
                     pdf.plot="7s2.Eiv")
  PlotSiteMismatches("lsy-6", "equilibrium", corrected_kds=FALSE, 12, 16,
                     pdf.plot="7s2.Eiii")
  PlotSiteMismatches("lsy-6", "equilibrium", corrected_kds=FALSE, 12, 15,
                     pdf.plot="7s2.Eii")
  PlotSiteMismatches("lsy-6", "equilibrium", corrected_kds=FALSE, 9, 16,
                     titles=TRUE, bulge=TRUE, pdf.plot="7s2.Exi")
  PlotSiteMismatches("lsy-6", "equilibrium", corrected_kds=FALSE, 10, 16,
                     bulge=TRUE, pdf.plot="7s2.Ex")
  PlotSiteMismatches("lsy-6", "equilibrium", corrected_kds=FALSE, 12, 17,
                     bulge=TRUE, pdf.plot="7s2.Eix")
  PlotSiteMismatches("lsy-6", "equilibrium", corrected_kds=FALSE, 12, 16,
                     bulge=TRUE, pdf.plot="7s2.Eviii")
  PlotSiteMismatches("lsy-6", "equilibrium", corrected_kds=FALSE, 12, 15,
                     bulge=TRUE, pdf.plot="7s2.Evii")
  # F___________________________________________________________________________
  PlotBestThreePrimeSite("miR-7-23nt", experiment="equilibrium2_nb",
                         len_lim=c(4, 8), sitelist="randthrp_suppcomp",
                         corrected_kds=FALSE, justlegend=TRUE,
                         pdf.plot="7s2.Fi")
  PlotSiteMismatches("miR-7-23nt", "equilibrium2_nb", corrected_kds=FALSE, 10,
                     17, titles=TRUE, pdf.plot="7s2.Fvi")
  PlotSiteMismatches("miR-7-23nt", "equilibrium2_nb", corrected_kds=FALSE, 11,
                     17, pdf.plot="7s2.Fv")
  PlotSiteMismatches("miR-7-23nt", "equilibrium2_nb", corrected_kds=FALSE, 13,
                     18, pdf.plot="7s2.Fiv")
  PlotSiteMismatches("miR-7-23nt", "equilibrium2_nb", corrected_kds=FALSE, 11,
                     15, pdf.plot="7s2.Fiii")
  PlotSiteMismatches("miR-7-23nt", "equilibrium2_nb", corrected_kds=FALSE, 13,
                     16, pdf.plot="7s2.Fii")
  PlotSiteMismatches("miR-7-23nt", "equilibrium2_nb", corrected_kds=FALSE, 10,
                     17, titles=TRUE, bulge=TRUE, pdf.plot="7s2.Fxi")
  PlotSiteMismatches("miR-7-23nt", "equilibrium2_nb", corrected_kds=FALSE, 11,
                     17, bulge=TRUE, pdf.plot="7s2.Fx")
  PlotSiteMismatches("miR-7-23nt", "equilibrium2_nb", corrected_kds=FALSE, 13,
                     18, bulge=TRUE, pdf.plot="7s2.Fix")
  PlotSiteMismatches("miR-7-23nt", "equilibrium2_nb", corrected_kds=FALSE, 11,
                     15, bulge=TRUE, pdf.plot="7s2.Fviii")
  PlotSiteMismatches("miR-7-23nt", "equilibrium2_nb", corrected_kds=FALSE, 13,
                     16, bulge=TRUE, pdf.plot="7s2.Fvii")

  message("Done fig. 7s2.")
}










############## END FIGURES FOR THE PAPER #######################################
MakeFigure1()
# MakeFigure2()
# MakeFigure3()
# MakeFigure4()
# MakeFigure5()
# MakeFigure7()

# MakeFigure2s1()                 #1
# MakeFigure2s2()                 #2
# MakeFigure3s1()                 #3
# MakeFigure3s2()                 #3
# MakeFigure4s1()                 #4
# MakeFigure5s1()                 #5
# MakeFigure5s2()                 #6
# MakeFigure5s3()                 #7
# MakeFigure5s4()                 #8
# MakeFigure5s5()                 #9
# MakeFigure5s6()                 #10
# MakeFigure5s7()                 #11
# MakeFigure5s8()                 #12
# MakeFigure6s1()                 #13
# MakeFigure6s2()                 #14
# MakeFigure7s1()                 #15
# MakeFigure7s2()                 #16





