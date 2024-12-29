#23456789!123456789@123456789#123456789$123456789%123456789^123456789&123456789*


# # Assignment of consistent colors for all figures:
# kEquilCondColors <- c("black", "black", "gray60", "brown4", "darkgoldenrod1",
#                       "red", "orange", "forestgreen", "blue", "violet")
# names(kEquilCondColors) <- c("I", "I_combined", "0", "0.04", "0.126", "0.4", "1.26", "4", "12.6",
#                              "40")

# kBaekColors <- c("purple1", "firebrick", "blue", "cyan", "purple3", "purple2",
#                  "lightblue", "darkslategrey", "darkslategray4",
#                  "darkslategray3", "darkslategray2", "black")

# names(kBaekColors) <- c("8mer", "7mer-m8", "7mer-A1", "6mer", "7mer-m3.9",
#                         "6mer-m8", "6mer-A1", "CDNST 1", "CDNST 2", "CDNST 3",
#                         "CDNST 4", "None")

kMirnaColors <- c("deepskyblue2", "black", "red", "forestgreen",
                  "purple", "orange", "black", "orange")

names(kMirnaColors) <- kMirnas

kNucleotideColors <- c("blue", "green", "purple", "red")
names(kNucleotideColors) <- c("A", "T", "C", "G")
# kRNucleotideColors <- kNucleotideColors
# names(kRNucleotideColors) <- c("A", "U", "C", "G")

kSiteColors <- read.table(
  "general/site_info.txt",
  row.names=1, header=FALSE, sep="\t", stringsAsFactors=FALSE)[, 1, drop = FALSE]
names.temp <- rownames(kSiteColors)
kSiteColors <- kSiteColors[,1]
names(kSiteColors) <- names.temp
kSiteColors2 <- kSiteColors
names(kSiteColors) <- paste0(names(kSiteColors), "_Kd")
kSiteColors <- c(kSiteColors, kSiteColors2)
kSiteColorsbg <- rep(kSiteColors["bg"], length(kMirnas))
kSiteColorsAGO <- rep(kSiteColors["AGO"], length(kMirnas))
names(kSiteColorsbg) <- paste0("bg_", kMirnas)
names(kSiteColorsAGO) <- paste0("AGO_", kMirnas)

kSiteColors <- c(kSiteColors, kSiteColorsbg, kSiteColorsAGO)


# kSiteCategoryColors <- read.table(
#   "/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/site_info.txt",
#   row.names=1, header=FALSE, sep="\t", stringsAsFactors=FALSE)[, 6, drop = FALSE]
# names.temp <- rownames(kSiteCategoryColors)
# kSiteCategoryColors <- kSiteCategoryColors[,1]
# names(kSiteCategoryColors) <- names.temp

kSiteCategoryColors <- c("purple2",
                         "cyan",
                         "royalblue", 
                         "green3",
                         "gray",
                         "black",
                         "violet")
names(kSiteCategoryColors) <- c(unique(kSiteCategories), "Noncanonical")

# kPlotParameters <- list(
#   cex.main  = 1.5,
#   lwd       = 2,
#   pch       = 18,
#   cex.lab   = 1.5,
#   cex.axis  = 1.5,
#   ann       = FALSE,
#   font      = 1,
#   las       = 1,
#   mar       = c(5, 5, 4, 2) + 0.1,
#   font.main = 1,
#   bty       = "n",
#   mgp       = c(2.2, 1, 0))

# # Assignment of the plotting parameters for all figures:
# kPlotParameters <- list(
#   cex.main  = 1,
#   lwd       = 1.5,
#   lheight   = 1,
#   pch       = 16,
#   cex.lab   = 1,
#   lab       = c(5, 5, 2),
#   cex.axis  = 1,
#   ann       = FALSE,
#   font      = 1,
#   las       = 1,
#   mar       = c(3, 3, 2, 2),
#   font.main = 1,
#   tcl       = -0.2,
#   bty       = "n",
#   mgp       = c(0.5, 0.3, 0),
#   xaxs      = "i",
#   yaxs      = "i")


kPDFParameters <- list(
  cex       = 0.5,
  # cex.axis  = 0.4,
  # cex.lab   = 0.4,
  lwd       = 0.7,
  lheight   = 1, # The line height multiplier.
  pch       = 16,
  lab       = c(2, 2, 4/5), # doesn't matter
  ann       = FALSE,
  font      = 1,
  las       = 1, # Makes axis labels always horizontal.
  mar       = c(3, 3, 2, 2),
  font.main = 1,
  tcl       = -0.2,
  bty       = "n",
  mgp       = c(0.7, 0.4, 0), # Sets the axis label location relative to inner plot window
  xaxs      = "i",
  yaxs      = "i")

# # kCairoPDFParameters <- list(
# #   cex       = 0.5,
# #   # cex.axis  = 0.4,
# #   # cex.lab   = 0.4,
# #   lwd       = 0.55,
# #   ps        = 11,
# #   lheight   = 1, # The line height multiplier.
# #   pch       = 20,
# #   lab       = c(2, 2, 4/5), # doesn't matter
# #   ann       = FALSE,
# #   font      = 1,
# #   las       = 1, # Makes axis labels always horizontal.
# #   mar       = c(3, 3, 2, 2),
# #   font.main = 1,
# #   tcl       = -0.2,
# #   bty       = "n",
# #   mgp       = c(1, 0.35, 0), # Sets the axis label location relative to inner plot window
# #   xaxs      = "i",
# #   yaxs      = "i")

error_height_final <- 0.05 * 1.2
line_dash_length <- 2
pt_cex_final <- 1.35
legend_pch <- 19




#23456789!123456789@123456789#123456789$123456789%123456789^123456789&123456789*
GetDynLims <- function(side, inv=FALSE) {
  lim.matrix <- matrix(c("xmin", "xmax", "ymin", "ymax", "xmin", "xmax", "ymin", "ymax"), nrow=2, ncol=4)
  lims <- sapply(lim.matrix[, side], dynGet)
  if (inv) lims <- rev(lims)
  return(lims)
}
#23456789!123456789@123456789#123456789$123456789%123456789^123456789&123456789*

ExtendLims <- function(lims, extension=0.0, log=FALSE){
  pmin <- lims[1]
  pmax <- lims[2]
  if (log) return(c(pmin, pmin*(pmax/pmin)^(1 + extension)))
  else     return(c(pmin, pmin + (pmax - pmin)*(1 + extension)))
}


#23456789!123456789@123456789#123456789$123456789%123456789^123456789&123456789*
BlankPlot <- function(log="", inv="", adjusted=FALSE) {
  inv.l <- sapply(c("x", "y"), grepl, x=inv)
  xlim <- ExtendLims(GetDynLims(1, inv.l[1]), log=grepl("x", log))
  ylim <- ExtendLims(GetDynLims(2, inv.l[2]), log=grepl("y", log))
  if (adjusted) {
    xmin <- xlim[1]
    xmax <- xlim[2]
    xlim[2] <- BoxAdjustedXMax(log=grepl("x", log), inv=grepl("x", inv))
  }
  plot(1, type="n", log=log, xlim=xlim, ylim=ylim, axes=FALSE, ann=FALSE)
}
#23456789!123456789@123456789#123456789$123456789%123456789^123456789&123456789*
GetXYPlotDimRatio <- function() {
  y.plot.dist <- dev.size()[2]*5 - sum(par()$mar[c(1, 3)])*par()$cex
  x.plot.dist <- dev.size()[1]*5 - sum(par()$mar[c(2, 4)])*par()$cex
  return(x.plot.dist/y.plot.dist)
}
#23456789!123456789@123456789#123456789$123456789%123456789^123456789&123456789*
BoxAdjustedXMax <- function(log=FALSE, inv=FALSE){
  lims <- sapply(c("xmin", "xmax"), dynGet)
  xmin <- lims[1]
  xmax <- lims[2]
  xy.ratio <- GetXYPlotDimRatio()
  if (log) return(xmin*(xmax/xmin)^xy.ratio)
  else     return(xmin + (xmax - xmin)*xy.ratio)
}

GetPlotFractionalCoords <- function(fx, fy, log="", inv="") {
  # Places text within a plot based on the fractional position of the text
  # with respect to the axes. This should allow streamlining of plotting
  # figures.
  #
  # Args:
  # fx:      The fractional placement of the text along the x axis.
  # fy:      The fractional placement of the text along the y axis.
  # log:   A string to be one of "", "x", "y", or "xy", specifying which
  #          axes are log-transformed.
  inv.l <- sapply(c("x", "y"), grepl, x=inv)
  xlim <- GetDynLims(1, inv.l[1])
  ylim <- GetDynLims(2, inv.l[2])
  if (length(grep("x", log)) == 1) xlim <- log10(xlim)
  if (length(grep("y", log)) == 1) ylim <- log10(ylim)
  x.range <- xlim[2] - xlim[1]
  y.range <- ylim[2] - ylim[1]
  xpos <- xlim[1] + fx*x.range
  ypos <- ylim[1] + fy*y.range
  if (length(grep("x", log)) == 1) xpos <- 10^xpos
  if (length(grep("y", log)) == 1) ypos <- 10^ypos
  return(c(xpos, ypos))
}

XLabelAdj <- function(label, log=FALSE, ...) {
  dots <- list(...)
  names_dots <- names(dots)
  dots <- unlist(dots)
  dim.xy <- dev.size() # Get plot dimensions
  mar <- par()$mar*par()$cex # Get margin dimensions
  x.dist <- dim.xy[1]*5 - sum(mar[c(2, 4)]) # Get plot x dimension.
  y.dist <- dim.xy[2]*5 - sum(mar[c(1, 3)]) # Get plot y-dimension.
  label.dist <- strwidth(label, units="figure")*dev.size()[1]*5
  xlim1 <- GetDynLims(1)
  if (length(dots) > 0) {
    for (i in 1:length(dots)) {
      assign(names_dots[i], dots[i])
    }
  }
  xlim2 <- GetDynLims(1)
  if (xlim1[2] != xlim2[2]) {
    if (log) {
      y.dist <- y.dist*log(xlim2[2]/xlim2[1])/log(xlim1[2]/xlim1[1])
    }
  }
  xy <- SubfunctionCall(GetPlotFractionalCoords, fx=1, fy=1)
  return((0.5*y.dist - label.dist/2)/(x.dist - label.dist))
}

YLabelAdj <- function(label, log=FALSE, ...) {
  dots <- list(...)
  names_dots <- names(dots)
  dots <- unlist(dots)
  dim.xy <- dev.size() # Get plot dimensions
  mar <- par()$mar*par()$cex # Get margin dimensions
  y.dist <- dim.xy[2]*5 - sum(mar[c(1, 3)]) # Get plot y-dimension.
  label.dist <- strwidth(label, units="figure")*dev.size()[1]*5
  ylim1 <- GetDynLims(2)
  if (length(dots) > 0) {
    for (i in 1:length(dots)) {
      assign(names_dots[i], dots[i])
    }
  }
  ylim2 <- GetDynLims(2)
  if (ylim1[2] != ylim2[2]) {
    if (log) {
      y.dist2 <- y.dist*log(ylim2[2]/ylim2[1])/log(ylim1[2]/ylim1[1])
    } else {
      y.dist2 <- y.dist*(ylim2[2] - ylim2[1])/(ylim1[2] - ylim1[1])
    }
  } else {
    y.dist2 <- y.dist
  }
  return((0.5*y.dist2 - label.dist/2)/(y.dist - label.dist))
}


FormatKdLabel <- function(label){
  if (grepl("Kd", label)) {
    pre.string <- gsub("(^.*)Kd(.*$)", label, replace="\\1", perl=TRUE)
    post.string <- gsub("(^.*)Kd(.*$)", label, replace="\\2", perl=TRUE)
    if (grepl("log2", pre.string)) {
      mid.string <- gsub("(^.*)log2(.*$)", pre.string, replace="\\2", perl=TRUE)
      pre.string <- gsub("(^.*)log2(.*$)", pre.string, replace="\\1", perl=TRUE)
      return(bquote(.(pre.string)*log[2]*.(mid.string)*italic(K)[D]*.(post.string)))
    } else if (grepl("log2", post.string)) {
      mid.string <- gsub("(^.*)log2(.*$)", post.string, replace="\\1", perl=TRUE)
      post.string <- gsub("(^.*)log2(.*$)", post.string, replace="\\2", perl=TRUE)
      return(bquote(.(pre.string)*italic(K)[D]*.(mid.string)*log[2]*.(post.string)))
    } else if (grepl("log10", pre.string)) {
      mid.string <- gsub("(^.*)log10(.*$)", pre.string, replace="\\2", perl=TRUE)
      pre.string <- gsub("(^.*)log10(.*$)", pre.string, replace="\\1", perl=TRUE)
      return(bquote(.(pre.string)*log[2]*.(mid.string)*italic(K)[D]*.(post.string)))
    } else if (grepl("log10", post.string)) {
      mid.string <- gsub("(^.*)log10(.*$)", post.string, replace="\\1", perl=TRUE)
      post.string <- gsub("(^.*)log110(.*$)", post.string, replace="\\2", perl=TRUE)
      return(bquote(.(pre.string)*italic(K)[D]*.(mid.string)*log[2]*.(post.string)))
    } else {
      return(bquote(.(pre.string)*italic(K)[D]*.(post.string)))
    }
  } else if (grepl("log10", label)) {
    pre.string <- gsub("(^.*)log10(.*$)", label, replace="\\1", perl=TRUE)
    post.string <- gsub("(^.*)log10(.*$)", label, replace="\\2", perl=TRUE)
    return(bquote(.(pre.string)*log[2]*.(post.string)))
  } else {
    pre.string <- gsub("(^.*)log2(.*$)", label, replace="\\1", perl=TRUE)
    post.string <- gsub("(^.*)log2(.*$)", label, replace="\\2", perl=TRUE)
    return(bquote(.(pre.string)*log[2]*.(post.string)))
  }
}
#23456789!123456789@123456789#123456789$123456789%123456789^123456789&123456789*
AddLogAxis <- function(side, label, maglabel=1, line=1.2,
                       percent=FALSE, adj=FALSE, boxplot=FALSE, blank_lab=FALSE, ...) {
  if (sum(sapply(c("Kd", "log2", "log10"), grepl, x=label))) {
    label <- FormatKdLabel(label)
  }
  if (adj) adj <- XLabelAdj(label, log=TRUE, ...)
  else     adj <- 0.5
  dots <- list(...)
  names_dots <- names(dots)
  dots <- unlist(dots)
  if (length(dots) > 0) {
    for (i in 1:length(dots)) {
      assign(names_dots[i], dots[i])
    }
  }
  ticks.all <- sapply(seq(-10, 6), function(x) seq(9)*10^(x - 1))
  lims <- GetDynLims(side)
  if (boxplot) {
    lims <- 10^lims
  }
  ticks.pos <- c(lims[1], ticks.all[which(ticks.all > lims[1] &
                                          ticks.all < lims[2])], lims[2])
  label.pos <- 10^seq(ceiling(log10(min(ticks.pos))),
                      floor(log10(max(ticks.pos))), by=maglabel)
  if (boxplot) {
    ticks.pos <- log10(ticks.pos)
    label.pos <- log10(label.pos)
  }
  if (blank_lab) {
    labels <- rep("", length(label.pos))
  } else if (percent) {
    labels <- label.pos*100
  } else if (boxplot) {
    labels <- sapply(label.pos, function(name) {
    #   if (name < 0) {
    #     name <- paste0("\u2212", substr(name, 2, nchar(name)))
    #   }
      eval(substitute(expression(10^x), list(x=name)))
    })
  } else {
    labels <- sapply(label.pos, function(name) {
      str_name <- log10(name)
      # if (str_name < 0) {
      #   str_name <- paste0("\u2212", substr(str_name, 2, nchar(str_name)))
      # }
      eval(substitute(expression(10^x), list(x=str_name)))
    })
  }
  if (side==2) {
    tick.adj <- 0
    tick.line <- 1.5
  } else {
    tick.adj <- 0.5
    tick.line <- NA
  }
  # Make line.
  axis(side, at=c(min(ticks.pos), max(ticks.pos)), labels=FALSE, lwd=par()$lwd,
       lwd.ticks=0)
  # Make ticks.
  axis(side, at=label.pos, labels=labels, lwd=0, lwd.ticks=par()$lwd,
       tcl=par()$tcl*1.5)
  # Make minor ticks.
  axis(side, at=setdiff(ticks.pos, label.pos), labels=FALSE, lwd=0,
       lwd.ticks=par()$lwd)
  # condition for if the box requires adjustment:
  l_2 <- 0
  if (class(label) %in% c("call", "expression")) l_2 <- -0.16
  if (side==1) title(xlab=label, line=line - l_2, adj=adj)
  else if (side==4) title(ylab=label, line=line - 19.5, adj=adj, xpd=NA)
  else         title(ylab=label, line=line + l_2 + 0.70, adj=adj)
}


#23456789!123456789@123456789#123456789$123456789%123456789^123456789&123456789*
AddLinearAxis <- function(side, tick.space, label.space, label, line=1.2,
                          percent=FALSE, adj=FALSE, adj_pos=NULL, alt_lab=NULL,
                          alt_lab_pos=NULL, angled=FALSE, label_pos_ticks=TRUE,
                          alt_lab_y_dist=0.05, blank_lab=FALSE,
                          alt_tick_pos=FALSE, noline=FALSE,
                          alt_tick_lab_space=FALSE, ...) {
  # print(label.space)
  if (!(side %in% c(1, 2, 4))) {
    print("Error: 'side' argument must be 1 or 2.")
  } else {
    
    if (adj & side==1) adj <- XLabelAdj(label, ...)
    else if (adj & side==2) adj <- YLabelAdj(label, ...)
    else     adj <- 0.5
    # message(sprintf("adj is %s", adj))
    if (adj < 0) {
      adj <- 0
    }
    # message(sprintf("adj is %s", adj))
    dots <- list(...)
    names_dots <- names(dots)
    dots <- unlist(dots)
    if (length(dots) > 0) {
      for (i in 1:length(dots)) {
        assign(names_dots[i], dots[i])
      }
    }
    if (class(adj_pos) == "numeric") adj <- adj_pos
    lims <- GetDynLims(side)
    ticks.pos <- unique(c(lims[1], seq(lims[1], lims[2], by=tick.space),
                   lims[2]))
    if (length(alt_lab_pos) != 0) label.pos <- alt_lab_pos
    else                      label.pos <- seq(lims[1], lims[2], by=label.space)
    if (length(alt_lab) != 0) labels <- alt_lab
    else if (blank_lab)       labels <- rep("", length(label.pos))
    else if (length(alt_lab_pos) !=0)   labels <- alt_lab_pos    
    else if (percent)         labels <- round(label.pos*100, digits=1)
    else                      labels <- round(label.pos, digits=1)
    if (angled) labels <- rep("", length(label.pos))
    # Where the labels are written.
    if (alt_tick_lab_space) {
        ind_labs_use <- seq(1, length(label.pos), by=label.space)
        axis(side, at=label.pos[ind_labs_use], labels=labels[ind_labs_use], lwd=0, gap.axis=0)
    } else {
      axis(side, at=label.pos, labels=labels, lwd=0, gap.axis=0)
    }
    if (!noline) {
      if (alt_tick_pos) {
        axis(side, at=c(min(label.pos), max(label.pos)), lwd.ticks=0, labels=FALSE, lwd=par()$lwd)      
      } else {
        axis(side, at=c(min(ticks.pos), max(ticks.pos)), lwd.ticks=0, labels=FALSE, lwd=par()$lwd)      
      }
      if (label_pos_ticks) {
        if (alt_tick_lab_space) {
          axis(side, at=label.pos[ind_labs_use], tcl=par()$tcl*1.5, labels=FALSE, lwd=0, lwd.ticks=par()$lwd)      
          axis(side, at=label.pos[-ind_labs_use], tcl=par()$tcl, labels=FALSE, lwd=0, lwd.ticks=par()$lwd)      
        } else {
          axis(side, at=label.pos, tcl=par()$tcl*1.5, labels=FALSE, lwd=0, lwd.ticks=par()$lwd)      
        }
      }
      if (!(alt_tick_pos)) {
        axis(side, at=setdiff(ticks.pos, label.pos), labels=FALSE, lwd=0, lwd.ticks=par()$lwd)      
      }
    }
    if (sum(sapply(c("Kd", "log2", "log10"), grepl, x=label))) label <- FormatKdLabel(label)
    # condition for if the box requires adjustment:
    l_2 <- 0
    if (class(label) %in% c("call", "expression")) l_2 <- -0.16
    if (side==4) mtext(text=label, side=4, line=line + l_2 + 0.7, adj=adj, cex=par()$cex, las=0)
    else if (side==1) title(xlab=label, line=line - l_2, adj=adj)
    else         title(ylab=label, line=line + l_2 + 0.70, adj=adj)
    if (angled) {
      xy <- text(x=label.pos - 0.5, alt_lab_y_dist, labels=alt_lab, srt=40,
                     xpd=NA, adj=c(1, 1))
    }
  }
}
# #23456789!123456789@123456789#123456789$123456789%123456789^123456789&123456789*
# GetPlotFractionalCoords <- function(fx, fy, log="", inv="") {
#   # Places text within a plot based on the fractional position of the text
#   # with respect to the axes. This should allow streamlining of plotting
#   # figures.
#   #
#   # Args:
#   # fx:      The fractional placement of the text along the x axis.
#   # fy:      The fractional placement of the text along the y axis.
#   # log:   A string to be one of "", "x", "y", or "xy", specifying which
#   #          axes are log-transformed.
#   inv.l <- sapply(c("x", "y"), grepl, x=inv)
#   xlim <- GetDynLims(1, inv.l[1])
#   ylim <- GetDynLims(2, inv.l[2])
#   if (length(grep("x", log)) == 1) xlim <- log10(xlim)
#   if (length(grep("y", log)) == 1) ylim <- log10(ylim)
#   x.range <- xlim[2] - xlim[1]
#   y.range <- ylim[2] - ylim[1]
#   xpos <- xlim[1] + fx*x.range
#   ypos <- ylim[1] + fy*y.range
#   if (length(grep("x", log)) == 1) xpos <- 10^xpos
#   if (length(grep("y", log)) == 1) ypos <- 10^ypos
#   return(c(xpos, ypos))
# }
# #23456789!123456789@123456789#123456789$123456789%123456789^123456789&123456789*
AddCorrelationToPlot <- function(x, y, xpos, ypos, rsquared=FALSE, adj=0,
                                 col="black", round_dec=2, weights=NULL,
                                 method="pearson") {
  if (length(weights) != 0) {
    r <- weightedCorr(x, y, weights=weights, method="pearson")  
  } else {
    r <- cor(x, y, use="pairwise.complete.obs", method=method)  
  }
  if (rsquared) {
    r <- format(round(r^2, round_dec), nsmall=round_dec)
    text(x=xpos, y=ypos, bquote(italic(r)^2 == .(r)), col=col, adj=adj, xpd=NA)
  } else {
    r <- format(round(r, round_dec), nsmall=round_dec)
    if (method == "spearman") {
      text(x=xpos, y=ypos, bquote(italic(r)[S] == .(r)), col=col, adj=adj, xpd=NA)      
    } else {
      text(x=xpos, y=ypos, bquote(italic(r) == .(r)), col=col, adj=adj, xpd=NA)
    }
  }
}

# AddPValueToPlot <- function(p, xpos, ypos, adj=0,
#                                  col="black", round_dec=2, weights=NULL) {
# p <- format(round(p, round_dec), nsmall=round_dec)  
# text(x=xpos, y=ypos, bquote(italic(p) == .(p)), col=col, adj=adj, xpd=NA)
# }


ConvertRColortoRGB <- function(color.vector, alpha=1) {
  # Get the rgb values for each color:
  rgb.vals <- rbind(col2rgb(color.vector)/255, alpha)

  # Convert each of the columns in the rgb.vals matrix into its RGB code:
  return(apply(rgb.vals, 2, function(column) {
    r <- column[1]
    g <- column[2]
    b <- column[3]
    return(rgb(r, g, b, alpha=column[4]))
  }))
}

# GetColorFunction <- function(flanks, alpha=1) {
#   cols.counts <- c(-2, 1, 1, -1, 0, 0)
#   names(cols.counts) <- c("G", "A", "T", "C", ".", "|")
#   flank.colors <- rev(topo.colors(15, alpha=alpha)[1:13])
#   sapply(flanks, function(flank) {
#     if (flank %in% c(names(kSiteColors), "11mer-m3.13", "12mer-m3.14",
#                     "11mer-m4.14", "12mer-m4.15")) {
#       return("gray")
#     } else {
#       counts <- 9 + sum(cols.counts[unlist(strsplit(flank, NULL))])
#       return(flank.colors[counts])
#     }
#   })
# }

CombinedSiteAndFlankColors <- function(sXc) {
  colors <- rep("gray", nrow(sXc))
  flanks.inds <- grep("\\|", rownames(sXc), perl=TRUE)
  print(rownames(sXc)[flanks.inds])
  flanks     <- gsub("(^.*)\\|(.*$)", rownames(sXc)[flanks.inds], replace="\\2",
                     perl=TRUE)
  colors[flanks.inds] <- GetColorFunction(flanks)
  colors
}


GetColorFunction <- function(flanks, alpha=1) {
  cols.counts <- c(0, 1, 1, 0, 0, 0)
  names(cols.counts) <- c("G", "A", "T", "C", ".", "|")
  flank.colors <- rev(topo.colors(200, alpha=alpha)[c(1, 20,60,70,120)])
  sapply(flanks, function(flank) {
    if (flank %in% c(names(kSiteColors), "11mer-m3.13", "12mer-m3.14",
                    "11mer-m4.14", "12mer-m4.15")) {
      return("gray")
    } else {
      counts <- 1 + sum(cols.counts[unlist(strsplit(flank, NULL))])
      return(flank.colors[counts])
    }
  })
}

Points <- function(x, y, x_error=NULL, y_error=NULL, ln.lwd=par()$lwd, pt.lwd=0, cex=pt_cex_final,
                   line=FALSE, ...) {
  if (!is.null(x_error)) {
    inds_error <- which(x_error[[1]] - x_error[[2]] != 0)
    arrows(x0=x_error[[1]][inds_error], y0=y[inds_error],
                            x1=x_error[[2]][inds_error], y1=y[inds_error],
           length=error_height_final*par()$cex, angle=90, code=3, xpd=NA)
  }
  if (!is.null(y_error)) {
    inds_error <- which(y_error[[1]] - y_error[[2]] != 0)
    arrows(x0=x[inds_error], y0=y_error[[1]][inds_error],
                            x1=x[inds_error], y1=y_error[[2]][inds_error],
           length=error_height_final*par()$cex, angle=90, code=3, xpd=NA)
  }
  if (line) {
    lines(x, y, lwd=ln.lwd, xpd=NA, ...)
  }
    points(x, y, lwd=pt.lwd, cex=cex, xpd=NA, ...)  
}

Legend <- function(xy, legend, legend_pch_use=legend_pch, ...) {

    legend(x=xy[1], y=xy[2], legend=legend, pt.cex=pt_cex_final, pt.lwd=0,
           bty="n", pch=legend_pch_use, xpd=NA, ...)
  }

# ColorViridisPalette <- function(values,
#   steps, min, max, log=FALSE, nacolor="gray90", mincolor=NULL, maxcolor=NULL,
#   start=0.0, r=0.4, hue=0.8, palettemax=1
# ) {
#   color_ramp <- rev(cubeHelix(round(steps*(1/palettemax)), start=start, r=r, hue=hue))[1:steps]
#   if (class(nacolor) == "character") {
#     color_ramp <- c(color_ramp, nacolor)
#   }
#   if (class(mincolor) == "character") {
#     color_ramp[1] <- mincolor
#   }
#   if (class(maxcolor) == "character") {
#     color_ramp[steps] <- maxcolor
#   }
#   if (log) {
#     values <- log10(values)
#     min <- log10(min)
#     max <- log10(max)
#   }
#   # Transform the data such that the max and min values are 0 and 1
#   values <- (values - min)/(max - min)
#   # Further transform the data such that the max and min values are `1` and
#   # `steps`, respectively.
#   col_inds <- round(values*(steps - 1) + 1)
#   col_inds <- sapply(col_inds, function(col_ind) {
#     min(max(1, col_ind), steps)
#   })
#   col_inds[which(is.na(col_inds))] <- steps + 1
#   return(color_ramp[col_inds])
# }

