options(warn=0)

library(beeswarm)

mirnas_all <- c("miR-1", "let-7a", "miR-155", "miR-124", "lsy-6")

centered_sites <- c("11mer-m3.13", "12mer-m3.14",
                    "11mer-m4.14", "12mer-m4.15")

canonical_sites <- c("8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-A1", "6mer-m8")

kmer_list_names <- c("8mers", "9mers", "10mers", "11mers", "12mers")

exp_conditions <- c("I", "0", "0.4", "1.26", "4", "12.6", "40")

exp_cond_colors <- c("black", "gray60", "red", "orange", "forestgreen", "blue", "violet")
names(exp_cond_colors) <- exp_conditions

baek_colors <- c("purple1","firebrick","blue","cyan","purple3","purple2","lightblue","darkslategrey","darkslategray4","darkslategray3","darkslategray2","black")
names(baek_colors) <- c("8mer", "7mer-m8", "7mer-A1", "6mer", "7mer-m3.9", "6mer-m8", "6mer-A1",
                        "CDNST 1", "CDNST 2","CDNST 3", "CDNST 4","None")
mirna_colors <- c("deepskyblue2", "black", "red", "forestgreen", "purple")
names(mirna_colors) <- mirnas_all
repression.df <- data.frame(read.table("RepressionData/all_flanking_kds_and_repression.txt",header=TRUE,na.strings="",sep="\t", stringsAsFactors=FALSE))

mirna_sequences <- c("UUAAUGCUAAUCGUGAUAGGGGU")
names(mirna_sequences) <- "miR-155"

RC_vector <- c("A", "C", "G", "T")
names(RC_vector) <- c("U", "G", "C", "A")

GiveTargetSequence <- function(mirna, start, stop) {
  sequence <- mirna_sequences[mirna]
  sequence_window <- unlist(strsplit(sequence, split = ""))[start:stop]
  sequence_converted <- RC_vector[rev(sequence_window)]
  sequence_new <- paste0(sequence_converted, collapse = "")
  return(sequence_new)
}
kNucleotideColors <- c("blue", "green", "purple", "red")
names(kNucleotideColors) <- c("A", "T", "C", "G")
kSiteColors <- read.table(
  "/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/site_info.txt",
  row.names=1, header=FALSE, sep="\t", stringsAsFactors=FALSE)[, 1, drop = FALSE]

kSiteCategoryColors <- read.table(
  "/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/site_info.txt",
  row.names=1, header=FALSE, sep="\t", stringsAsFactors=FALSE)[, 6, drop = FALSE]


stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
                       "/SolveForKds/k_c_stockago.txt"), row.names=1,
                       header=TRUE, sep="\t")


GetSingleFlankPosition <- function(flanks, position) {
  if (position > 2) {
    position <- position + 1
  }
  output <- sapply(flanks, function(flank_name) {
    unlist(strsplit(flank_name, split=""))[position]
    })
return(output)

}

kPlotParameters <- list(
  cex.main  = 1.5,
  lwd       = 2,
  pch       = 20,
  cex.lab   = 1.5,
  cex.axis  = 1.5,
  ann       = FALSE,
  font      = 1,
  las       = 1,
  mar       = c(5, 5, 4, 2) + 0.1,
  font.main = 1,
  bty       = "n",
  mgp       = c(2.2, 1, 0))

kPlotParameters <- list(
  cex.main  = 1,
  lwd       = 1.5,
  lheight   = 1,
  pch       = 19,
  cex.lab   = 1,
  lab       = c(5, 5, 2),
  cex.axis  = 1,
  ann       = FALSE,
  font      = 1,
  las       = 1,
  mar       = c(3, 3, 2, 1),
  font.main = 1,
  tcl       = -0.2,
  bty       = "n",
  mgp       = c(0.5, 0.3, 0))




Norm <- function(vector) {
  return(vector/ sum(vector))
}

Cumul <- function(vector) {
  norm <- Norm(vector)
  tot <- 0
  out <- sapply(norm, function(x){
    tot <<- tot + x
    return(tot)

    })
  return(out)
}

Logistic <- function(vector, max) {
  return(max / (1 + exp(-vector)))
}

Logit <- function(vector, max) {
  print(vector)
  return(-log(max / vector - 1))
}


## FUNCTIONS USED IN EQUILIBRIUM MODELING (4 total)
# 1.________
GetOccupancy <- function(c.freeago, kds) {
  # Generates occupancy matrix for the entire input matrix which has rows
  # that are all possible flnaking nucleotide combinations, and columns that
  # are the probabilities of being unpaired in the window.
  #
  # Args:
  # c.freeago: The concentration of free AGO in the binding reaction
  # kds: a list of all kds corresponding to the individual site types including
  #   a non-specific binding kd.
  # IDs: The matrix which gives the site types in column 1 and the flanking
  #    nucleotide identity in column 2.
  # freq.flankXp: The matrix which splits the input into all site and flanking
  #   identities, as well as along the possible average windows for being
  #   unpaired.
  # p.exponent: The exponent used for the weight of geometric average
  #   probability of being unpaired.
  #
  # Returns:
  # A matrix giving the coresponding occuapncy (between 0 and 1) for each
  #   position in the matrix, where 
  return(c.freeago * (c.freeago + kds)^(-1))
}



GetContamResidual <- function(frees, kds, kd.contam, c.tots, c.ago, c.contam){
  c.freeago <- 10^frees[1]
  c.freecontam <- 10^frees[2]
  occs <- GetOccupanciesContaminant(c.freeago, c.freecontam, kds, kd.contam)
  residual_ago <- (c.ago - c.freeago - sum(occs[[1]]*c.tots))^2
  residual_contam <- (c.contam - c.freecontam - sum(occs[[2]]*c.tots))^2
  # print(residual_ago)
  # print(residual_contam)
  # print(c.freeago)
  # print(c.freecontam)
  return(sum(residual_ago, residual_contam))
}

GetFreeAgoAndContam <- function(kds,kd.contam,c.tots,c.ago,c.contam) {
  # print(c.ago)
  # print(c.contam)
  solution <- try(optim(log10(c(c.ago,c.contam)),
    GetContamResidual,
    kds=kds,
    kd.contam=kd.contam,
    c.tots=c.tots,
    c.ago=c.ago,
    c.contam=c.contam))
  
  return(10^solution$par)

}


GetContamResidual <- function(c.freecontam, kd.contam, c.tots, c.ago, c.freeago, c.contam){
  return(c.freecontam^2 + (sum(c.tots) - c.ago - c.contam + c.freeago + kd.contam)*c.freecontam - c.contam*kd.contam)
}

GetContamResidual2 <- function(c.freecontam, kd.contam, c.tots, c.ago, c.freeago, c.contam){
  return((c.freecontam^2 + (sum(c.tots) - c.ago - c.contam + c.freeago + kd.contam)*c.freecontam - c.contam*kd.contam)^2)
}


GetFreeContam <- function(kds,kd.contam,c.tots,c.ago,c.freeago,c.contam) {
  solution <- NaN
  try(solution <- uniroot(GetContamResidual,
    c(0,c.contam),
    kd.contam=kd.contam,
    c.tots=c.tots,
    c.ago=c.ago,
    c.contam=c.contam,
    c.freeago=c.freeago,
    tol = 0.00001 * .Machine$double.eps^0.25)$root,silent=TRUE)
  if (is.na(solution)) {
    solution <- optimize(GetContamResidual2,
    c(0,c.contam),
    kd.contam=kd.contam,
    c.tots=c.tots,
    c.ago=c.ago,
    c.contam=c.contam,
    c.freeago=c.freeago,
    tol = 0.00001 * .Machine$double.eps^0.25)$minimum
  }
  return(solution)

}

FreeAgoAndContamCombinedRoot <- function(c.freeago,kds,kd.contam, c.tots,c.ago,c.contam){
    c.freecontam <- GetFreeContam(kds,kd.contam,c.tots,c.ago,c.freeago,c.contam)

    occs_I <- GetOccupancy(c.freeago,kds*(c.freecontam + kd.contam)*(kd.contam)^(-1))
    root <- c.ago - c.freeago - sum(occs_I*c.tots)
    return(root)
}

GetFreeAgoAndContam <- function(kds, kd.contam, c.tots, c.ago, c.contam) {
  c.freeago <- uniroot(FreeAgoAndContamCombinedRoot,
                      c(0,c.ago),
                      kds=kds,
                      kd.contam = kd.contam,
                      c.tots = c.tots,
                      c.ago = c.ago,
                      c.contam = c.contam,
                      tol = 0.00001 * .Machine$double.eps^0.25)$root
  c.freecontam <- GetFreeContam(kds,kd.contam,c.tots,c.ago,c.freeago,c.contam)
  return(c(c.freeago,c.freecontam))
}

GetOccupanciesContaminant <- function(c.freeago, c.freecontam, kds, kd.contam) {
  return(list(c.freeago * (kds * (1 + c.freecontam / kd.contam) + c.freeago)^(-1),
              c.freecontam * (kd.contam * (1 + c.freeago / kds) + c.freecontam)^(-1)))
}


# 2._______
GetFreeRoot <- function(c.freeago, kds, c.tots, c.ago) {

  occs_I <- GetOccupancy(c.freeago, kds)

  root <- c.ago - c.freeago - sum(occs_I*c.tots)
  return(root)
}

# 3.___________
GetFreeResidual <- function(c.freeago, kds, c.tots, c.ago) {
  occs_I <- GetOccupancy(c.freeago, kds)
  residual <- (c.ago - c.freeago - sum(occs_I*c.tots))^2
  return(residual)
}


# 4._______
GetBoundRNA <- function(kds, c.tots, c.ago) {
  if (c.ago > 0) {
    c.free <- NaN
    try(c.free <- uniroot(GetFreeRoot,
                          c(0, c.ago),
                          kds = kds,
                          c.tots = c.tots,
                          c.ago = c.ago,
                          tol = 0.0001 * .Machine$double.eps^0.25)$root)
    if (is.na(c.free)) {
      c.free <- optimize(GetFreeResidual,
                         c(0, c.ago),
                         kds = kds,
                         c.tots = c.tots,
                         c.ago = c.ago,
                         tol = 0.00001*.Machine$double.eps^0.25)$minimum
    } 
    c.bound <- GetOccupancy(c.free, kds)*c.tots
    return(c.bound)
  } else {
    return(rep(0, length(c.tots)))
  }
}

GetFreeAgo <- function(kds, c.tots, c.ago) {
  if (c.ago > 0) {
    c.free <- NaN
    try(c.free <- uniroot(GetFreeRoot,
                          c(0, c.ago),
                          kds = kds,
                          c.tots = c.tots,
                          c.ago = c.ago,
                          tol = 0.00001 * .Machine$double.eps^0.25)$root)
    if (is.na(c.free)) {
      c.free <- optimize(GetFreeResidual,
                         c(0, c.ago),
                         kds = kds,
                         c.tots = c.tots,
                         c.ago = c.ago,
                         tol = 0.00001*.Machine$double.eps^0.25)$minimum
    } 
    return(c.free)
  } else {
    return(0)
  }
}


## Functions that load data. (__ total)
# 1.___________
GetSitesXCounts <- function(mirna, exp, n_constant, sitelist,
                            mirna.start=NULL, mirna.stop=NULL) {
  if (sitelist %in% kmer_list_names) {
   sites_file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
                             mirna, "/", exp,
                             "/site_count_tables/all_sites_",
                             n_constant, "_", sitelist, "_", mirna.start, "-",
                             mirna.stop, ".txt")
    sitesXcounts <- read.table(sites_file_name, stringsAsFactors=FALSE)
    colnames.temp <- colnames(sitesXcounts)[3 : ncol(sitesXcounts)]
    colnames.new <- sapply(colnames.temp, function(name) {
        return(unlist(strsplit(name, split="A", fixed=TRUE))[2])
      })
    colnames(sitesXcounts)[3 : ncol(sitesXcounts)] <- colnames.new
  } else {
    sites_file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
                         mirna, "/",exp,"/site_count_tables/all_sites_",
                         n_constant,"_", sitelist, ".txt")
    sitesXcounts <- read.table(sites_file_name, stringsAsFactors=FALSE)
    colnames.temp <- colnames(sitesXcounts)[4 : ncol(sitesXcounts)]
    colnames.new <- sapply(colnames.temp, function(name) {
        return(unlist(strsplit(name, split="A", fixed=TRUE))[2])
      })
    colnames(sitesXcounts)[4 : ncol(sitesXcounts)] <- colnames.new
  }
  return(sitesXcounts)
}

GetSitesXCountsUnique <- function(mirna, experiment, start, stop, sitelist,
                            mirna.start=NULL, mirna.stop=NULL) {
  if (sitelist %in% kmer_list_names) {
   sites_file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
                             mirna, "/", experiment,
                             "/full_site_count_tables_unique/all_sites_",
                             start, "-", stop, "_", sitelist, "_", mirna.start, "-",
                             mirna.stop, ".txt")
    sitesXcounts <- read.table(sites_file_name, stringsAsFactors=FALSE)
    colnames.temp <- colnames(sitesXcounts)[2 : ncol(sitesXcounts)]
    colnames.new <- sapply(colnames.temp, function(name) {
        return(unlist(strsplit(name, split="A", fixed=TRUE))[2])
      })
    colnames(sitesXcounts)[2 : ncol(sitesXcounts)] <- colnames.new
  } else {
    sites_file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
                         mirna, "/",experiment,"/full_site_count_tables_unique/all_sites_",
                         start, "-", stop,"_", sitelist, ".txt")
    sitesXcounts <- read.table(sites_file_name, stringsAsFactors=FALSE)
    colnames.temp <- colnames(sitesXcounts)[3 : ncol(sitesXcounts)]
    colnames.new <- sapply(colnames.temp, function(name) {
        return(unlist(strsplit(name, split="A", fixed=TRUE))[2])
      })
    colnames(sitesXcounts)[3 : ncol(sitesXcounts)] <- colnames.new
  }
  return(sitesXcounts)
}



GetSitesXCounts.Kinetics <- function(mirna, experiment, start, stop, sitelist) {
  sites_file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
                       mirna, "/",experiment,"/full_site_count_tables/all_sites_",
                       start, "-", stop, "_", sitelist, "_pulse.txt")
  sitesXcounts.pulse <- read.table(sites_file_name)
  sites_file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
                       mirna, "/",experiment,"/full_site_count_tables/all_sites_",
                       start, "-", stop, "_", sitelist, "_chase.txt")
  sitesXcounts.chase <- read.table(sites_file_name)

  return(list(sitesXcounts.pulse,sitesXcounts.chase))
}



# 3.___________
GetSiteFlanksXCounts <- function(mirna, experiment, n_constant,
                                 sitelist, site) {
  sites_file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
                       mirna, "/",experiment,"/site_count_tables/", site,
                       "_flanking_", n_constant, "_", sitelist, ".txt")
  print(sites_file_name)
  sitesXcounts <- read.table(sites_file_name, stringsAsFactors=FALSE)
  print(dim(sitesXcounts))
  colnames.temp <- colnames(sitesXcounts)[3 : ncol(sitesXcounts)]
  colnames.new <- sapply(colnames.temp, function(name) {
    return(unlist(strsplit(name, split="A", fixed=TRUE))[2])
  })
  colnames(sitesXcounts)[3 : ncol(sitesXcounts)] <- colnames.new
  return(sitesXcounts)
}

GetSiteFlanksXCountsKinetics <- function(mirna, experiment, site, start, stop,
                                 sitelist) {
  sites_file_name_pulse <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
                       mirna, "/",experiment,"/full_site_count_tables/", site,
                       "_flanking_", start, "-", stop,"_", sitelist, "_pulse.txt")
  sitesXcounts.pulse <- read.table(sites_file_name_pulse, stringsAsFactors=FALSE)
  sites_file_name_chase <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
                       mirna, "/",experiment,"/full_site_count_tables/", site,
                       "_flanking_", start, "-", stop,"_", sitelist, "_chase.txt")
  sitesXcounts.chase <- read.table(sites_file_name_chase, stringsAsFactors=FALSE)

  return(list(sitesXcounts.pulse, sitesXcounts.chase))
}

GetSiteKmersXCountsKinetics <- function(mirna, experiment, site, start, stop,
                                 sitelist, k) {
  sites_file_name_pulse <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
                       mirna, "/",experiment,"/full_site_count_tables/", site,
                       "_sitekmers_", start, "-", stop,"_", sitelist, "_k", k, "_pulse.txt")
  sitesXcounts.pulse <- read.table(sites_file_name_pulse, stringsAsFactors=FALSE)
  sites_file_name_chase <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
                       mirna, "/",experiment,"/full_site_count_tables/", site,
                       "_sitekmers_", start, "-", stop,"_", sitelist, "_k", k, "_chase.txt")
  sitesXcounts.chase <- read.table(sites_file_name_chase, stringsAsFactors=FALSE)

  return(list(sitesXcounts.pulse, sitesXcounts.chase))
}


# 4.
LinearModelInputFlanks <- function(flank_matrix,col) {
  flanks <- matrix(sapply(c(1,2,4,5), function(index) {
  sapply(rownames(flank_matrix), function(flank_name) {
    unlist(strsplit(flank_name, split=""))[index]
    })
  }), nrow=nrow(flank_matrix), ncol=4, byrow=FALSE)
  input_flanks <- data.frame(I = as.numeric(flank_matrix[,col]),
                             f5p.o = flanks[, 1],
                             f5p.i = flanks[, 2],
                             f3p.i = flanks[, 3],
                             f3p.o = flanks[, 4],
                             stringsAsFactors = FALSE)
  return(input_flanks)

}

GetColorFunction <- function(name,alpha=1) {
  if (name %in% rownames(kSiteColors) | name %in% c("11mer-m3.13", "12mer-m3.14",
                    "11mer-m4.14", "12mer-m4.15")) {
    return("gray")
  } else {
    i <- c(sapply(unlist(strsplit(name,".",fixed=TRUE)), function(x) {
      unlist(strsplit(x,"",fixed=TRUE))
      }))
    cols_counts <- list(c(231,41,138)/255,c(50,200,159)/255,c(50,200,159)/255,c(231,41,138)/255)
    cols_counts <- c(-1,0,1,1)

    names(cols_counts) <- c("G","A","T","C")
    counts <- c(0,0,0)
    counts <- 5
    # names(counts) <- c("red","green","blue")
    for (nucleotide in i) {
      counts <- counts + cols_counts[nucleotide]
    }
    # red = counts[1]
    # green = counts[2]
    # blue = counts[3]
    return(rev(topo.colors(15)[1:9])[counts])
    # return(rev(terrain.colors(9, alpha = 0.8)[1:9])[counts])

  }
}

ConvertTtoU <- function(site) {
  sites.split <- strsplit(site, "")
  out <- unlist(lapply(sites.split, function(site.split) {
    ind.T <- grep("T",site.split)
    site.split[ind.T] <- "U"
    return(paste0(site.split, collapse = ""))
    }))
  return(out)
}

MakeIterationPlot <- function(out, type, extension="") {
  setEPS()
  postscript(file=paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
                         "/figures/kds/", mirna, "/", type, "/iterations/",
                         site, "_", k.c.stockago, extension, ".eps")
            )
  par(kPlotParameters)
  x = seq(dim(out)[1])
  probs = out[ ,"-logp"]
  out.print <- out[ , seq(dim(out)[2] - 1)]
  ys <- 10^c(floor(min(out.print)), ceiling(max(out.print)))

  probs.scaled <- probs / max(probs) * (log10(ys[2]) - log10(ys[1])) +
                  log10(ys[1])

  out.print <- 10^out.print

  plot(x , 10^probs.scaled, log='y', axes=FALSE, type="l", ylim=ys,
       lwd=2, ann=FALSE,
       col="black")
  title(main=mirna, line=-1, adj=0.1)
  title(main=site, col.main=kSiteColors[site,], line=-2.5, adj=0.1)

  title(xlab = "Iteration")
  title(ylab = "Parameter values (nM)")
  axis(1, at=seq(1, max(x), by=max(1, floor(max(x) / 20))), pos=ys[1], lwd=2,
       labels=FALSE, tck=-0.01)
  axis(1, at=seq(1, max(x), by=max(1, floor(max(x) / 20))*5), pos=ys[1], lwd=2)
  axis(2, at=10^seq(log10(ys[1]), log10(ys[2])),
       labels=sapply(seq(log10(ys[1]), log10(ys[2])), function(name) {
         eval(substitute(expression(10^x), list(x=name)))
       }),
       pos=1, lwd=2, hadj=0.8)
  axis(2, at=c(sapply(seq(log10(ys[1]), log10(ys[2])-1), function(x) seq(9)*10^x),ys[2]), labels=FALSE,
       pos=1, lwd=2)
  sapply(colnames(out.print), function(name) {
          lines(x, out.print[, name], lwd=2, col=GetColorFunction(name))
        }
        )
  lines(x, 10^probs.scaled, type="l", col="black")
  dev.off()
}

PlotFlankKdRepression <- function(mirna,cutoff=FALSE,merge=FALSE,noncanon=TRUE){
    par(kPlotParameters)
  data<- GetSitesXCounts(mirna, "equilibrium", 5, 5, sitelist="paper")[,1,drop=FALSE]
  kds <- GetKds(mirna,"equilibrium", 5, 5, "paper",nosite=TRUE)
  data_seqs <- data[,1]
  names(data_seqs) <- rownames(data)

  print(data_seqs)

  sites_all <- unlist(unique(subset(repression.df,mir==mirna,select=site_type)))
  sites_all <- sites_all[sites_all!="nosite"]
  site_seqs <- sapply(sites_all, function(site) {
    data_seqs[site]
  })
  names(site_seqs) <- sites_all
  print('sites_all')
  print(sites_all)

  print(site_seqs)
  print(data_seqs)
  site_seqs_noncanonical <- site_seqs
  site_seqs_noncanonical <- site_seqs_noncanonical[grep(data_seqs["6mer"], site_seqs_noncanonical, invert=TRUE)]
  site_seqs_noncanonical <- site_seqs_noncanonical[grep(data_seqs["6mer-A1"], site_seqs_noncanonical, invert=TRUE)]
  site_seqs_noncanonical <- site_seqs_noncanonical[grep(data_seqs["6mer-m8"], site_seqs_noncanonical, invert=TRUE)]
  site_seqs_noncanonical <- site_seqs_noncanonical[names(site_seqs_noncanonical)!="8mer-bG(6.7)"]
  site_seqs_noncanonical <- site_seqs_noncanonical[names(site_seqs_noncanonical)!="7mer-m8bG(6.7)"]

  site_seqs_canonical <- setdiff(sites_all,names(site_seqs_noncanonical))
  print("canonical")
  print(site_seqs_canonical)
  print("noncanonical")
  print(site_seqs_noncanonical)
  if (merge == TRUE){
    repression.df$site_type[which(repression.df$site_type %in% names(site_seqs_noncanonical))] <- "Noncanonical"
      print(unique(repression.df$site_type))
    sites_all <- c(site_seqs_canonical, "Noncanonical")
  }

  out <- sapply(sites_all, function(site) {
    print(site)
    if (cutoff != FALSE & site=="Noncanonical") {
      reduced_frame <- subset(repression.df, mir==mirna & site_type==site & log_kd<=cutoff ,select=c(log_fc,log_kd))
    } else {
      reduced_frame <- subset(repression.df, 
        mir==mirna & site_type==site,
        select=c(log_fc,log_kd)
      )      
    }
    print(dim(reduced_frame))
    mean_values <- colMeans(reduced_frame)
    sd_values <- apply(reduced_frame,2,sd)/sqrt(nrow(reduced_frame))
    return(c(mean_values,sd_values))
  })
  par(mfrow=c(1,1))
  colnames(out) <- sites_all
  print(out)
  if (noncanon == FALSE) {
    out <- out[,c("8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-A1", "6mer-m8")]
    sites_all <- colnames(out)
  }


  print(out)
      xmin <- 1
    xmax <- 2000
    ymin <- -1
    ymax <- 0.1
    nosite_rep <- mean(subset(repression.df,mir==mirna & site_type=="nosite",select=c(log_fc))[,1])
  plot(c(1,1),c(1,1),xlim=c(xmin, xmax), ylim = c(ymin,ymax),log='x',col="white",cex=2,pch=19, lwd=2,ann=FALSE,axes=FALSE)
  arrows(1/(2^out[2,]), out[1,]+out[3,] - nosite_rep,1/(2^out[2,]),out[1,]-out[3,]-nosite_rep,length=0.05, lwd=1.5,angle=90, code=3)
  arrows(1/(2^(out[2,]+out[4,])), out[1,]-nosite_rep,1/(2^(out[2,]-out[4,])),out[1,]-nosite_rep,length=0.05, lwd=1.5,angle=90, code=3)
  points(1/(2^out[2,]),out[1,]-nosite_rep,col="black", bg=kSiteColors[sites_all,],pch=21,cex=1.5)


    xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
    ys <- seq(ymin,ymax,by=0.05)

    xs <- xs[xs >= xmin & xs <= xmax]
    ys <- ys[ys >= ymin & ys <= ymax]

    # xmin <- min(xs)
    # xmax <- max(xs)
    # ymin <- min(ys)
    # ymax <- max(ys)
    xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
    yl <- seq(ymin,ymax, by= 0.1)





    axis(1, at=xl,
         labels=sapply(xl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=ymin, lwd=0, cex.axis=1.7, padj=0.2)
    axis(1, at=xs, labels=FALSE,
         pos=ymin, lwd=2)
    # Label the axis at each order of magnitude.

    axis(2, at=yl,
         labels=round(yl,2),
         pos=xmin, las=2, lwd=0, cex.axis=1.7,hadj=0.8)
    axis(2, at=ys, labels=FALSE,
         pos=xmin, lwd=2)


  remove <- unique(which(is.na(out[1,])), which(is.na(out[1,])))
  print(remove)
  print(length(remove))
  if (length(remove) > 0) {
    out <- out[,-remove]  
  }
  text(x=800, y=0.05, mirna,cex=1.5)

  text(x=800, y=-0.02, eval(substitute(expression(italic(r) == x), 
            list(x = round(cor(out[2,],out[1,]),3)))),cex=1.5)
  if (merge == TRUE) {
  text(x=500, y=-0.09, eval(substitute(expression(log[2](italic(K)[D][italic(noncanon)]) <= x),
            list(x = cutoff))), cex=1.5)
  }

  # Second plot
  # kds <- GetKds(mirna,"equilibrium", 5, 5, "paper",nosite=TRUE)
  # sites_all <- unlist(unique(subset(repression.df,mir==mirna,select=site_type)))
  # sites_all <- sites_all[sites_all!="nosite"]
  # out <- sapply(sites_all, function(site) {
  #   reduced_frame <- subset(repression.df, mir==mirna & site_type==site,select=c(log_fc))
  #   mean_values <- colMeans(reduced_frame)
  #   sd_values <- apply(reduced_frame,2,sd)/sqrt(nrow(reduced_frame))
  #   return(c(mean_values,log(kds[site],base=2),sd_values))
  # })
  # title(main=mirna,font.main=1,cex=2)
  # plot(out[2,],out[1,],xlim=c(-10,1),ylim=c(-1,0.25),col=kSiteColors[sites_all,],cex=2,pch=seq(1,6),ann=FALSE,axes=FALSE)
  # arrows(out[2,], out[1,]+out[3,],out[2,],out[1,]-out[3,],length=0.05, angle=90, code=3, lwd=1.5)
  # points(out[2,],out[1,],col=kSiteColors[sites_all,],pch=seq(1,6),cex=2)
  # axis(1,at=seq(-10,1),pos = -1,lwd = 2)
  # axis(2, at = seq(-1,0.25,by=0.25),pos = -10,lwd=2)
  out <<- out
  sites_all <- sites_all[order(out[2,])]
  legend(x=1.5,y=-0.5, legend=sites_all, bty="n", pch=19,col=kSiteColors[sites_all,],cex=1.1, ncol=1)
  sites_all <<- sites_all
  title(xlab=expression(italic(K)[D]))
  title(ylab=expression(log[2](paste("fold change"))))
  # text(x=30,y=0.1,round(cor(out[2,],out[1,]),3),cex=1.5)


}

PlotKdRepression <- function(mirna){
    par(kPlotParameters)
  kds <- GetKds(mirna,"equilibrium", 5, 5, "paper",nosite=TRUE)
  sites_all <- unlist(unique(subset(repression.df,mir==mirna,select=site_type)))
  sites_all <- sites_all[sites_all!="nosite"]
  print(sites_all)
  out <- sapply(sites_all, function(site) {
    print(site)
    reduced_frame <- subset(repression.df, mir==mirna & site_type==site,select=c(log_fc))
    print(reduced_frame)
    mean_values <- colMeans(reduced_frame)
    sd_values <- apply(reduced_frame,2,sd)/sqrt(nrow(reduced_frame))
    return(c(mean_values,log(kds[site],base=2),sd_values))
  })
  colnames(out) <- sites_all
  print(out)
  plot(out[2,],out[1,],xlim=c(-10,1),ylim=c(-1,0.25),col=kSiteColors[sites_all,],cex=2,pch=seq(1,5),ann=FALSE,axes=FALSE)
  arrows(out[2,], out[1,]+out[3,],out[2,],out[1,]-out[3,],length=0.05, angle=90, code=3, lwd=1.5)
  points(out[2,],out[1,],col=kSiteColors[sites_all,],pch=seq(1,5),cex=2)
  axis(1,at=seq(-10,1),pos = -1,lwd = 2)
  axis(2, at = seq(-1,0.25,by=0.25),pos = -10,lwd=2)
  legend(x=, legend=sites_all,pch=seq(1,5),col=kSiteColors[sites_all,])

}


WriteIterationFile <- function(out,extension="") {
  out.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                     "/equilibrium/kds/", site,"_flanking_",
                     k.c.stockago, extension,".txt")
  write.table(file=out.file, out, sep="\t", quote=FALSE, row.names=FALSE,
                col.names=TRUE)
}

WriteFinalParameterFile <- function(out,extension="") {
  out.final <- out[dim(out)[1], ]
  names(out.final) <- colnames(out)
  out.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                           "/equilibrium/kds/final_", site,
                           "_flanking_", k.c.stockago,extension, ".txt")
  write.table(file=out.file, out.final, sep="\t", quote=FALSE, row.names=TRUE,
              col.names=FALSE)
}

## Helpful Fnctions

CompareSiteInputsN <- function(mirna, experiment, n_constant_1, n_constant_2, sitelist,column) {
  column_1 <- GetSitesXCounts(mirna, experiment, n_constant_1, sitelist)[,column,drop=FALSE]
  column_2 <- GetSitesXCounts(mirna, experiment, n_constant_2, sitelist)[,column,drop=FALSE]
  out <- cbind(column_1, column_2)
  plot(out[,1], out[,2], log='xy')

  return(out)
}

CompareSiteFlanksInputsN <- function(mirna, experiment, n_constant_1, n_constant_2, sitelist, site, column) {
  column_1 <- GetSiteFlanksXCounts(mirna, experiment, n_constant_1, sitelist, site)[,column,drop=FALSE]
  column_2 <- GetSiteFlanksXCounts(mirna, experiment, n_constant_2, sitelist, site)[,column,drop=FALSE]
  out <- cbind(column_1, column_2)
  plot(out[,1], out[,2], log='xy')
  return(out)
}










MakeSiteIterationPlot <- function(out,method,mirna, num.kds, num.bgs, colors = FALSE) {
  setEPS()
  postscript(file=paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/figures/kds/",
    mirna,"/", method, ".eps"))
  par(kPlotParameters)
  x = seq(1,dim(out)[1],length = max(ncol(out),1000))
  out <- out[x,]
  probs = out[ ,"-logp"]
  out.print <- out[ , seq(dim(out)[2] - 1)]
  out.print.kds <- Logistic(out[,1:num.kds],max = 1)
  out.print.bgs <- 10^out.print[,(num.kds + 1) : (ncol(out.print))]
  out.print <- out.print <- cbind(out.print.kds,out.print.bgs)
  ys <- 10^c(max(floor(log10(min(out.print))), -5), ceiling(log10(max(out.print))))
  probs.scaled <- probs / max(probs) * (log10(ys[2]) - log10(ys[1])) + log10(ys[1])


  plot(x    = x,
       y    = 10^probs.scaled,
       log  = "y",
       axes = FALSE,
       type = "l",
       col = "white",
       ylim = ys)
  axis(side = 1,
       at   = seq(1, max(x), by=max(1, floor(max(x) / 20))),
       pos  = ys[1], labels = FALSE,
       tck  = -0.01)
  axis(1, at=seq(1, max(x), by=max(1, floor(max(x) / 20))*5), pos=ys[1], lwd=2)
  axis(2, at=10^seq(log10(ys[1]), log10(ys[2])),
       labels=sapply(seq(log10(ys[1]), log10(ys[2])), function(name) {
         eval(substitute(expression(10^x), list(x=name)))
       }),
       pos=1, lwd=2)
  title(main = mirna, font=1)
  axis(2, at=c(sapply(seq(log10(ys[1]), log10(ys[2])-1), function(x) seq(9)*10^x),ys[2]), labels=FALSE,
       pos=1)
  if (colors != FALSE) {

  cols <- c(colors,kSiteColors[colnames(out.print)[(length(colors)+1):length(colnames(out.print))],])
  names(cols) <- colnames(out.print)
  cols["AGO"] <- "grey"

  sapply(colnames(out.print),function(name) {
    lines(x, out.print[, name], lwd=1, col=cols[name])
    })
  } else  {
    sapply(colnames(out.print),function(name) {
    lines(x, out.print[, name], lwd=1, col=kSiteColors[name,])
    })
  }
  lines(x, 10^probs.scaled, type="l", lwd=2, col="black")
  dev.off()
}


# Made FOR PAPER NOW
GetSiteKds <- function(mirna, experiment, n_constant, sitelist) {
    params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                 experiment, "/kds_PAPER/", n_constant, "_", 
                 sitelist, "_PAPER.txt")

    params <- data.frame(read.table(params.file, header=TRUE, row.names=1,
                         stringsAsFactors=FALSE))
    return(params)
}

GetFlankKds <- function(mirna, experiment, n_constant, sitelist, site) {
    params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                 experiment, "/kds_PAPER/", n_constant, "_", 
                 sitelist, "_", site, "_PAPER.txt")

    params <- data.frame(read.table(params.file, header=TRUE, row.names=1,
                         stringsAsFactors=FALSE))
    return(params)
}



PlotSiteKds <- function(mirna, experiment, n_constant, sitelist, plotlist,
                        xpos = 20, ypos = 20) {
  # Get kds for all site-types of the mirna.
  dev.new(xpos = xpos, ypos = ypos, height = 5, width = 8)

  kds <- GetSiteKds(mirna, experiment, n_constant, sitelist)
  if (length(plotlist) == 0) {
    site_list <- rownames(data)
  } else if (class(plotlist) == "character") {
    site_list <- c(read.table(file = paste0("/lab/bartel1_ata/mcgeary/",
                                            "computation/AgoRBNS/",
                                            "AssignSiteTypes/sites.", mirna,
                                            "_", plotlist, ".txt"),
                              stringsAsFactors=FALSE)[,1], "None")
  } else {
    site_list <- c(rownames(data)[order(kds)][1:site_list], "None")
  }
  kds <- kds[-which(rownames(kds) %in% c("bg", "AGO")),]
  kds <- kds[order(kds$Full),]

  par(kPlotParameters)
  xs <- kds
  ys <- nrow(kds) - seq(nrow(kds)) + 1
  plot(1, type ="n",
          axes    = FALSE,
          log = 'x',
          ylim       = c(0, 22),
          xlim       = rev(c(0.00003, 3)))
 arrows(kds$Upper_CI, nrow(kds) - seq(nrow(kds)) + 1,
        kds$Lower_CI, nrow(kds) - seq(nrow(kds)) + 1, length=0.05, angle=90, code=3)

  title(main = mirna,
        line = -2,
        adj  = 0.1)
  title(xlab = expression(italic(K)[D]))
 points(kds$Mean,nrow(kds) - seq(nrow(kds)) + 1,
          col = "black",
          bg = kSiteCategoryColors[rownames(kds),],
          pch = 21,
          lwd=1,
          cex = 1.2)
  ymin=0.0001
  ymax=3
  ys <- c(sapply(seq(-5,6), function(x) seq(9)*10^(x-1)))
  ys <- ys[which(ys>=ymin & ys <= ymax)]
  # ys <- ys[ys >= ymin & ys <= ymax]
  yl <- 10^seq(ceiling(log10(min(ys))), floor(log10(max(ys))))
  axis(1, at=yl,
       labels=sapply(yl, function(name) {
         eval(substitute(expression(10^x), list(x=log10(name))))
       }),
       pos=0, lwd=0)
  axis(1, at=ys, labels=FALSE,
       pos=0, lwd=1)
# points(rep(0.0006, nrow(kds)),nrow(kds) - seq(nrow(kds)) + 1,col=kSiteCategoryColors[rownames(kds),],cex=1.5,pch=19)


  text((kds$Full)/1.6,nrow(kds) - seq(nrow(kds)) + 1, labels=ConvertTtoU(rownames(kds)), adj=0, cex = 0.8, col= "black")


  if (mirna == "miR-124") {
      legend(x = 10^-2.1, y = 8,legend = c("7-8-nt canonical site",
                                 "6-nt canonical site",
                                 "Enhanced 6mer-containing sites",
                                 "Noncanonical sites",
                                 "3' sites"),
         bty="n",
         col="black",
         pt.bg = c("purple2", "deepskyblue2", "cyan", "violet", "green3"),
         pch = 21,
         pt.cex = 1.2,
         pt.lwd=1)

      } else if (mirna == "miR-155") {
  legend(x = 10^-2.1, y = 8,legend = c("7-8-nt canonical site",
                                 "6-nt canonical site",
                                 "Noncanonical sites",
                                 "3' sites",
                                 "??"),
         bty="n",
         col="black",
         pt.bg = c("purple2", "deepskyblue2", "violet", "green3","gray"),
         pch = 21,
         pt.cex = 1.2,
         pt.lwd=1)
} else if (mirna == "miR-1") {
    legend(x = 10^-2.1, y = 7,legend = c("7-8-nt canonical site",
                                 "6-nt canonical site",
                                 "Noncanonical sites",
                                 "??"),
         bty="n",
         col="black",
         pt.bg = c("purple2", "deepskyblue2", "violet", "gray"),
         pch = 21,
         pt.cex = 1.2,
         pt.lwd=1)

} else if (mirna == "let-7a") {
    legend(x = 10^-2.1, y = 6.5,legend = c("7-8-nt canonical site",
                                 "6-nt canonical site",
                                 "Noncanonical sites",
                                 "??"),
         bty="n",
         col="black",
         pt.bg = c("purple2", "deepskyblue2", "violet", "gray"),
         pch = 21,
         pt.cex = 1.2,
         pt.lwd=1)

} else {
  legend(x = 10^-2.1, y = 9.5,legend = c("7-8-nt canonical site",
                                 "6-nt canonical site",
                                 "Enhanced 6mer-containing sites",
                                 "Noncanonical sites",
                                 "3' sites",
                                 "??"),
         bty="n",
         col="black",
         pt.bg = c("purple2", "deepskyblue2", "cyan", "violet", "green3","gray"),
         pch = 21,
         pt.cex = 1.2,
         pt.lwd=1)

  }
}

PlotSiteEnrichments <- function(mirna, experiment, n_constant, sitelist,
                                plotlist, xpos = 20, ypos = 20,
                                bgoff = FALSE) {
  params <- GetSiteKds(mirna, experiment, n_constant, sitelist)
  
  sitesXcounts <- GetSitesXCounts(mirna, experiment, n_constant, sitelist)
  sitesXcounts <- sitesXcounts[,-1]

  # Remove any rows with no reads whoatsoever.
  sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
  # Remove any rows for which there are no input reads:
  sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]
  kds <- params[1:(nrow(params)-2),]$Mean
  names(kds) <- rownames(sitesXcounts)

  bgs <- rep(params["bg",]$Mean, 5)
  k.c.stockago <- params["AGO",]$Mean
  k.c.stockago.plot <- stockago[mirna, "equilibrium"]
  c.I.tots <- Norm(sitesXcounts[,2])*100
  names(c.I.tots) <- rownames(sitesXcounts)
  data <- sitesXcounts[,3:7]
  names(c.I.tots) <- rownames(data)
  x_model <- seq(min(as.numeric(colnames(data))) * 0.7,
                 max(as.numeric(colnames(data))) / 0.7,
                 length=100)

  c.totals <- matrix(c.I.tots, nrow=nrow(data), ncol = 100, byrow=FALSE)
  colnames(c.totals) <- colnames(x_model)
  rownames(c.totals) <- rownames(data)

  c.agos <- sapply(x_model, function(x) {
      as.numeric(x) * k.c.stockago / 100
    }
  )

  c.bounds <- as.matrix(
    sapply(c.agos, function(x) {
        return(GetBoundRNA(kds, c.I.tots, x))
      }
    )
  )
  c.frees <- c.totals - c.bounds
  c.bgs <- t(t(c.frees) * bgs / colSums(c.frees))

  if (bgoff == TRUE) {
    c.all <- c.bounds
  } else {
    c.all <- c.bounds + c.bgs
  }

  c.final <- data.frame(t(t(c.all) / colSums(c.all)))
  rownames(c.final) <- rownames(data)
  x <- c(40,12.65,4,1.265,0.4)*k.c.stockago.plot/100*1000
  y <- c(1,1,1,1,1)
  sites.norm <- Norm(c.I.tots)
  data.norm <- t(t(data)/colSums(data))

  data.R <- data.norm/(sites.norm)
  model.R <- c.final/(sites.norm)
  xmin <- min(x)*0.3
  xmax <- max(x)*3
  ymin <- 0.2
  ymax <- 300
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
  dev.new(xpos = xpos, ypos = ypos, height = 5, width = 6)
  par(kPlotParameters)

  plot(x, y,log='xy', xlim=c(xmin, xmax*(xmax/xmin)^0.4), ylim=c(ymin, ymax), type="l",
     col="white", axes=FALSE, ann=FALSE)        
  # Generate tickmarks for axis.

  axis(1, at=xl,
       labels=sapply(xl, function(name) {
         eval(substitute(expression(10^x), list(x=log10(name))))
       }),
       pos=ymin, lwd=0)
  axis(1, at=xs, labels=FALSE,
       pos=ymin)
  # Label the axis at each order of magnitude.

  axis(2, at=yl,
       labels=sapply(yl, function(name) {
         eval(substitute(expression(10^x), list(x=log10(name))))
       }),
       pos=xmin, las=2, lwd=0)
  axis(2, at=ys, labels=FALSE,
       pos=xmin)

  title(main = mirna, font.main=1, line=-2, adj=0.1)
  title(xlab = "[AGO2-miRNA] (pM)", adj=0.3)
  title(ylab = "Enrichment")

  plotlist <- c(read.table(file = paste0("/lab/bartel1_ata/mcgeary/",
                                          "computation/AgoRBNS/",
                                          "AssignSiteTypes/sites.", mirna,
                                            "_", plotlist, ".txt"),
                              stringsAsFactors=FALSE)[,1], "None")

  legend.names <- rownames(data)[order(kds)]
  legend.names <- legend.names[which(legend.names %in% plotlist)]
  ordered_list <- legend.names


  legend(x=xmax, y=ymax, legend=ConvertTtoU(legend.names), pch=19,
         col=kSiteColors[ordered_list, ], bty="n")
  for (name in plotlist) {
    type = "p"
    points(x, data.R[name, ], col=kSiteColors[name, ], type=type, pch=19, cex=1.2)
    lines(x_model*k.c.stockago.plot/100*1000, model.R[name, ], col=kSiteColors[name, ],lwd=2)      
  }
}

PlotSiteOccupancy <- function(mirna, experiment, n_constant, sitelist,
                                plotlist, xpos = 20, ypos = 20,
                                bgoff = FALSE) {
  params <- GetSiteKds(mirna, experiment, n_constant, sitelist)
  
  sitesXcounts <- GetSitesXCounts(mirna, experiment, n_constant, sitelist)
  sitesXcounts <- sitesXcounts[,-1]

  # Remove any rows with no reads whoatsoever.
  sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
  # Remove any rows for which there are no input reads:
  sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]
  kds <- params[1:(nrow(params)-2),]$Mean
  names(kds) <- rownames(sitesXcounts)

  bgs <- rep(params["bg",]$Mean, 5)
  k.c.stockago <- params["AGO",]$Mean
  k.c.stockago.plot <- stockago[mirna, "equilibrium"]

  c.I.tots <- Norm(sitesXcounts[,2])*100
  names(c.I.tots) <- rownames(sitesXcounts)
  data <- sitesXcounts[,3:7]
  names(c.I.tots) <- rownames(data)
  x_model <- seq(min(as.numeric(colnames(data))) * 0.7,
                 max(as.numeric(colnames(data))) / 0.7,
                 length=100)
  x_points <- c(40,12.65,4,1.265,0.4)

  c.totals <- matrix(c.I.tots, nrow=nrow(data), ncol = 100, byrow=FALSE)
  colnames(c.totals) <- colnames(x_model)
  rownames(c.totals) <- rownames(data)

  c.agos_model <- sapply(x_model, function(x) {
      as.numeric(x) * k.c.stockago / 100
    }
  )

  c.agos_points <- sapply(x_points, function(x) {
      as.numeric(x) * k.c.stockago / 100
    }
  )
  c.bounds_model <- as.matrix(
    sapply(c.agos_model, function(x) {
        return(GetBoundRNA(kds, c.I.tots, x))
      }
    )
  )

  c.bounds_points <- as.matrix(
    sapply(c.agos_points, function(x) {
        return(GetBoundRNA(kds, c.I.tots, x))
      }
    )
  )

  c.bounds_model_norm <- apply(c.bounds_model, 2, Norm)
  c.bounds_points_norm <- apply(c.bounds_points, 2, Norm)
  rownames(c.bounds_model_norm) <- rownames(sitesXcounts)
  rownames(c.bounds_points_norm) <- rownames(sitesXcounts)
  xmin <- min(x_points*k.c.stockago.plot/100*1000)*0.3
  xmax <- max(x_points*k.c.stockago.plot/100*1000)*3
  ymin <- 0
  ymax <- 0.7
  yextension <- (ymax/ymin)
  xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
  xs <- xs[xs >= xmin & xs <= xmax]
  xmin <- min(xs)
  xmax <- max(xs)
  xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
  ys <- seq(0,7)/10
  yl <- ys
  dev.new(xpos = xpos, ypos = ypos, height = 5, width = 6)
  par(kPlotParameters)

  plot(1,
       type = "n",
       log='x',
       xlim=c(xmin, xmax*(xmax/xmin)^0.4),
       ylim=c(ymin, ymax),
       axes=FALSE,
       ann=FALSE)        
  # Generate tickmarks for axis.

  axis(1, at=xl,
       labels=sapply(xl, function(name) {
         eval(substitute(expression(10^x), list(x=log10(name))))
       }),
       pos=ymin, lwd=0)
  axis(1, at=xs, labels=FALSE,
       pos=ymin)
  # Label the axis at each order of magnitude.

  axis(2, at=yl,
       labels=yl,
       pos=xmin, las=2, lwd=0)
  axis(2, at=ys, labels=FALSE,
       pos=xmin)

  title(main = mirna, font.main=1, line=-2, adj=0.1)
  title(xlab = "[AGO2-miRNA] (pM)", adj=0.3)
  title(ylab = "Enrichment")

  plotlist <- c(read.table(file = paste0("/lab/bartel1_ata/mcgeary/",
                                          "computation/AgoRBNS/",
                                          "AssignSiteTypes/sites.", mirna,
                                            "_", plotlist, ".txt"),
                              stringsAsFactors=FALSE)[,1], "None")

  legend.names <- rownames(data)[order(kds)]
  legend.names <- legend.names[which(legend.names %in% plotlist)]
  ordered_list <- legend.names


  legend(x=xmax, y=ymax, legend=ConvertTtoU(legend.names), pch=19,
         col=kSiteColors[ordered_list, ], bty="n")
  for (name in plotlist) {
    type = "p"
    points(x_points*k.c.stockago.plot/100*1000, c.bounds_points_norm[name,], col=kSiteColors[name, ], type=type, pch=19,cex=1.2)
    lines(x_model*k.c.stockago.plot/100*1000, c.bounds_model_norm[name,], col=kSiteColors[name, ], lwd = 2)      
  }
}



PlotSiteScatterWithInput <- function(mirna, experiment, n_constant, sitelist, column,
                                     xpos = 20, ypos = 20) {
  
  sitesXcounts <- GetSitesXCounts(mirna, experiment, n_constant, sitelist)
  sitesXcounts <- sitesXcounts[,-1]

  # Remove any rows with no reads whoatsoever.
  sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
  # Remove any rows for which there are no input reads:
  sitesXcounts <- sitesXcounts[sitesXcounts[,2]>0,]

  x <- Norm(sitesXcounts[,2])*100
  y <- Norm(sitesXcounts[,column])*100

  xymin <- 0.05
  xymax <- 100
  yextension <- (xymax/xymin)
  xys <- c(sapply(seq(floor(log10(xymin)), ceiling(log10(xymax))), function(x) seq(10)*10^x))
  xys <- xys[xys >= xymin & xys <= xymax]
  xymin <- min(xys)
  xymax <- max(xys)
  xyl <- 10^seq(ceiling(log10(min(xys))), floor(log10(max(xys))))

  dev.new(xpos = xpos, ypos = ypos, height = 5, width = 6)
  par(kPlotParameters)
  # Make dummy plot:
  plot(1, type = "n",log='xy',
       xlim=c(xymin, xymax*(xymax/xymin)^0.4),
       ylim=c(xymin, xymax),
      axes=FALSE, ann=FALSE)
  # Make x=y line:
  segments(xymin, xymin, xymax, xymax, lty = 2)
  # Make the lines connecting the points to the x = y line:
  segments(x,x,x,y,lt = 2, col = kSiteColors[rownames(sitesXcounts), ])
  # Make axes:
  axis(1, at=xyl,
       labels=xyl,
       pos=xymin, lwd=0)
  axis(1, at=xys, labels=FALSE,
       pos=xymin)
  axis(2, at=xyl,
       labels=xyl,
       pos=xymin, las=2, lwd=0)
  axis(2, at=xys, labels=FALSE,
       pos=xymin)
  # Add the points to the plot:
  points(x, y,
     col=kSiteColors[rownames(sitesXcounts),])        

  ago.percent <- as.numeric(colnames(sitesXcounts)[column])/100
  title(main = paste0(ago.percent*stockago[mirna,"equilibrium"]*1000,' pM AGO2-', mirna), line=-2.5, adj=0.1)
  title(xlab = "Input library (%)", adj=0.4)
  title(ylab = "AGO-bound library(%)")


  legend(x=50, y=2, legend=ConvertTtoU(rownames(sitesXcounts)), pch=19,
         col=kSiteColors[rownames(sitesXcounts), ], bty="n", y.intersp=0.9)
  }




PlotSiteFlankEnrichments <- function(mirna, experiment, n_constant, sitelist, plotlist, site, bgoff = FALSE, xpos = 20, ypos = 20) {
  dev.new(xpos = xpos, ypos = ypos, height = 5, width = 6)
  par(kPlotParameters)
  params <- GetSiteKds(mirna, experiment, n_constant, sitelist)
  flank.kds <- GetFlankKds(mirna, experiment, n_constant, sitelist, site)
  sitesXcounts <- GetSitesXCounts(mirna, experiment, n_constant, sitelist)
  sitesXcounts <- sitesXcounts[,-1]

  s.c <- as.numeric(sitesXcounts[site, ])

  sfXc <- GetSiteFlanksXCounts(mirna, experiment, n_constant, sitelist, site)

  colnames(sfXc)[1] <- "I"
  sfXc <- sfXc[rowSums(sfXc[, 2:6]) > 0,]
  sfXc <- sfXc[which(sfXc[,1] > 0),,drop=FALSE]

  sfXc <- round(t((t(sfXc) / colSums(sfXc)) * s.c))
  sfXc[is.na(sfXc)] <- 0

  # Remove any rows with no reads whoatsoever.
  sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
  # Remove any rows for which there are no input reads:
  sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]
  kds <- params[1:(nrow(params)-2),]$Mean
  names(kds) <- rownames(sitesXcounts)
  num.sf <- nrow(sfXc)

  kds.s <- params$Mean[1:num.kds]
  names(kds.s) <- rownames(params)[1:num.kds]
  kd.site <- kds.s[names(kds.s)==site]
  # Omit the site kd for which the flanking sites are being fit.
  kds.s <- kds.s[names(kds.s) != site]

  sitesXcounts <- rbind(sitesXcounts, sfXc)

  sitesXcounts <- sitesXcounts[rownames(sitesXcounts) != site, ]


  bgs <- rep(params["bg",]$Mean, 5)
  k.c.stockago <- params["AGO",]$Mean

  c.I.tots <- Norm(sitesXcounts[,2])*100
  names(c.I.tots) <- rownames(sitesXcounts)
  data <- sitesXcounts[,3:7]
  names(c.I.tots) <- rownames(data)
  x_model <- seq(min(as.numeric(colnames(data))) * 0.7, max(as.numeric(colnames(data))) / 0.7,
                 length=100)

  c.totals <- matrix(c.I.tots, nrow=nrow(data), ncol = 100, byrow=FALSE)
  colnames(c.totals) <- colnames(x_model)
  rownames(c.totals) <- rownames(data)

  c.agos <- sapply(x_model, function(x) {
      as.numeric(x) * k.c.stockago / 100
    }
  )
  kds <- c(kds.s,flank.kds$Mean)
  c.bounds <- as.matrix(
    sapply(c.agos, function(x) {
        return(GetBoundRNA(kds, c.I.tots, x))
      }
    )
  )
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
  x <- c(40,12.65,4,1.265,0.4)*k.c.stockago/100*1000/2
  print(x)
  y <- c(1,1,1,1,1)
  print(x)
  sites.norm <- Norm(c.I.tots)
  data.norm <- t(t(data)/colSums(data))

  data.R <- data.norm/(sites.norm)
  model.R <- c.final/(sites.norm)

  xmin <- min(x)*0.3
  xmax <- max(x)*3
  ymin <- 0.2
  ymax <- 300
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

  plot(x, y,log='xy', xlim=c(xmin, xmax*(xmax/xmin)^0.4), ylim=c(ymin, ymax), type="l",
     col="white", axes=FALSE, ann=FALSE)        
  # Generate tickmarks for axis.

  axis(1, at=xl,
       labels=sapply(xl, function(name) {
         eval(substitute(expression(10^x), list(x=log10(name))))
       }),
       pos=ymin, lwd=0, )
  axis(1, at=xs, labels=FALSE,
       pos=ymin)
  # Label the axis at each order of magnitude.

  axis(2, at=yl,
       labels=sapply(yl, function(name) {
         eval(substitute(expression(10^x), list(x=log10(name))))
       }),
       pos=xmin, las=2, lwd=0, )
  axis(2, at=ys, labels=FALSE,
       pos=xmin)

  title(main = mirna, font.main=1, line = -2, adj=0.1)
  title(xlab = "[AGO2-miRNA] (pM)", adj=0.3)
  title(ylab = "Enrichment")

  plotlist <- c(read.table(file = paste0("/lab/bartel1_ata/mcgeary/",
                                          "computation/AgoRBNS/",
                                          "AssignSiteTypes/sites.", mirna,
                                            "_", plotlist, ".txt"),
                              stringsAsFactors=FALSE)[,1], "None")

  legend.names <- rownames(data)[order(kds)]
  legend.names <- legend.names[which(legend.names %in% plotlist)]
  ordered_list <- legend.names


  site_colors <- rep("gray", length(plotlist)-1)
  flank_colors <- sapply(rownames(flank.kds),GetColorFunction)
  print(flank_colors)

  colors_all <- c(site_colors, flank_colors)
  names(colors_all) <- rownames(data)
  for (name in rownames(data)) {
    type = "p"
    points(x, data.R[name, ], col=colors_all[name], type=type, pch=19,
           cex=1.2, lwd=3)
    lines(x_model*k.c.stockago/100*1000/2, model.R[name, ], col=colors_all[name], lwd=2)      
  }
}


PlotBaekKds <- function(mirna, experiment, n_constant, xpos = 20, ypos =20) {
  # Make plot window.
  dev.new(xpos = xpos, ypos = ypos, height = 5, width = 8)
  # Get the kds.
  kds <- GetSiteKds(mirna, experiment, n_constant, "baek")
  # Get the list of sites (Redundant with kd names?)
  site_list <- c(
    read.table(
      file = paste0("/lab/bartel1_ata/mcgeary/",
                    "computation/AgoRBNS/",
                    "AssignSiteTypes/sites.", mirna,
                    "_", "baek", ".txt"),
      stringsAsFactors=FALSE)[,1], "None")
  # Remove Ago and bg parameters from kd list.
  kds <- kds[-which(rownames(kds) %in% c("bg", "AGO")),]
  # Make a data frame where each kd is associated with a sitetype
  # category.
  baek_categories <- rownames(kds)
  baek_categories[which(baek_categories=="5mer-m2.6")] <- "CDNST 1"
  baek_categories[grep("7mer-A1mm",baek_categories)] <- "CDNST 2"
  baek_categories[grep("8mer-m2.9", baek_categories)] <- "CDNST 3"
  baek_categories[grep("8mer-mm", baek_categories)] <- "CDNST 4"
  kds <- data.frame(Mean = kds$Mean,
                    Sitetype = baek_categories,
                    stringsAsFactors = FALSE)
  # Make dataframe with means for each category.
  sitetype.means.df = aggregate(
    kds$Mean, list(kds$Sitetype),function(x) {10^mean(log10(x))})

  # Get the ordering vector (for colors).
  site.order = order(sitetype.means.df[,2])
  # Converts the order to the the per-site positions.

  # I.E.: If the vector is c(7, 6, 5, 1, etc., then the
  # new vector will have a 1 at position 7, a 2, at position 6,
  # a 3 at position 5, a 4 at position 1, etc.)
  site.rank = sapply(
    seq(
      length(site.order)), function(ind) {
        return(which(site.order == ind))
    }
  )
  par(kPlotParameters)
  # Initial plot.
  plot(1, type ="n",
          axes    = FALSE,
          ann = FALSE,
          log = 'x',
          ylim       = c(0, 22),
          xlim       = rev(c(0.00003, 3)))
  # Makes the mirna title.
  title(main = mirna,
        line = -2,
        adj  = 0.1)
  # Makes the Kd axis title.
  title(xlab = expression(italic(K)[D]))
  
  # AXES:
  ymin=0.0001
  ymax=3

  # Makes the axes ticks.
  # ys refers to the short ticks, of which there are more
  # yl refers to the long ticks, of which there are fewer.
  ys <- c(sapply(seq(-5,6), function(x) seq(9)*10^(x-1)))
  ys <- ys[which(ys>=ymin & ys <= ymax)]
  yl <- 10^seq(ceiling(log10(min(ys))), floor(log10(max(ys))))
  axis(1, at=yl,
       labels=sapply(yl, function(name) {
         eval(substitute(expression(10^x), list(x=log10(name))))
       }),
       pos=0, lwd=0)
  axis(1, at=ys, labels=FALSE,
       pos=0)

  # Plot the beeswarm.
  beeswarm(Mean ~ Sitetype,
          data=kds,
          pch = 21,
          col = "black",
          bg = baek_colors[levels(factor(kds$Sitetype))],
          method = "swarm",
          pt.lwd= 1,
          cex = 1.2,
          add = TRUE,
          corral = "random",
          corralWidth = 0.5,
          # ylim = c(2,0.00001),
          at = length(site.rank) -site.rank +1,
          horizontal = TRUE)

  # Plots the names of the site types to the right of each category.
  labelpositions.x <- aggregate(kds$Mean, list(kds$Sitetype), min)[,2]/1.6
  text(x      = labelpositions.x,
       y      = length(site.rank)-site.rank +1, 
       labels = sitetype.means.df[,1],
       adj    = 0,
       cex = 0.9,
       col    = "black")
}
SortKdsFile <- function(mirna, sitelist) {
    # This prints out a new site list that is ordered as per the mean
    # Kd value from the original sitelist. This allows the Kds to be fit a
    # second time to make sure that the values are robust to the ordering.
    kds <- GetSiteKds(mirna, "equilibrium", 5, sitelist)
    kds <- kds[-which(rownames(kds) %in% c("bg", "AGO", "None")),]

    write.table(file=paste0("AssignSiteTypes/sites.", mirna, "_", sitelist,
                            ",ordered.txt"),
                x=rownames(kds)[order(kds$Mean)], col.names = FALSE, row.names= FALSE,quote=FALSE)
    print(kds)
    print(kds[order(kds$Mean),])
}


PlotPositionalKds <- function(experiment, n_constant, sitelist, xpos = 20, ypos = 20) {
  # Get kds for all site-types of the mirna.

  positional_sites <- c("8mer", "6mer", "6mer-m8",
                        "11mer-m3.13",
                        "11mer-m4.14",
                        "11mer-m5.15",
                        "11mer-m6.16",
                        "11mer-m7.17",
                        "11mer-m8.18",
                        "11mer-m9.19",
                        "11mer-m10.20",
                        "11mer-m11.21",
                        "11mer-m12.22",
                        "11mer-m13.23",
                        "None")
  dev.new(xpos = xpos, ypos = ypos, height = 5, width = 8)
  kds <- GetSiteKds("miR-1", experiment, n_constant, sitelist)    
  kds <- kds[positional_sites,]
  print(kds)
  par(kPlotParameters)

  # Make plot with miR-1 data:
  ind_p <- c(1,2,3,nrow(kds))
  y_p <- (nrow(kds) - seq(nrow(kds)) + 1)[ind_p]
  ind_l <- seq(nrow(kds))[-ind_p]
  y_l <- (nrow(kds) - seq(nrow(kds)) + 1)[ind_l]
  plot(1, type = "n",
        col = "white",
        axes    = FALSE,
        log = 'x',
        ylim       = c(0, 22),
        xlim       = rev(c(0.00003, 5)))

  title(xlab = expression(italic(K)[D]))
  ymin=0.0001
  ymax=5
  ys <- c(sapply(seq(-5,6), function(x) seq(9)*10^(x-1)))
  ys <- ys[which(ys>=ymin & ys <= ymax)]
  yl <- 10^seq(ceiling(log10(min(ys))), floor(log10(max(ys))))
  axis(1, at=yl,
       labels=sapply(yl, function(name) {
         eval(substitute(expression(10^x), list(x=log10(name))))
       }),
       pos=0)
  axis(1, at=ys, labels=FALSE,
       pos=0)

  for (mirna in c("miR-1","let-7a", "miR-124", "lsy-6", "miR-155")) {
    kds <- GetSiteKds(mirna, experiment, n_constant, sitelist)    
  kds <- kds[positional_sites,]
  points(kds$Mean[ind_p], y_p,
        col = mirna_colors[mirna],
        pch=19)
  lines(kds$Mean[ind_l], y_l,
        col = mirna_colors[mirna],
        lwd = 2,
        type = "o")
  lines(kds$Lower_CI[ind_l], y_l,
        col = mirna_colors[mirna],
        lwd = 1,
        lty=2,
        cex = 1)
  lines(kds$Upper_CI[ind_l], y_l,
        col = mirna_colors[mirna],
        lwd = 1,
        lty=2,
        cex = 1)
  }
  text(6.5,
       nrow(kds) - seq(nrow(kds)) + 1,
       labels=rownames(kds),
       cex=0.9,
       adj=0,
       col= "black")
  legend(x = 10^-3,
         y = 8,
         legend = mirnas_all,
         col = mirna_colors,
         bty="n",
         pch = 19)
}


PlotSiteFlankKds <- function(mirna, experiment, n_constant, sitelist,
                                plotlist, xpos = 20, ypos = 20) {

  if (length(plotlist) == 0) {
    site_list <- rownames(data)
  } else if (class(plotlist) == "character") {
    site_list <- c(read.table(file = paste0("/lab/bartel1_ata/mcgeary/",
                                            "computation/AgoRBNS/",
                                            "AssignSiteTypes/sites.", mirna,
                                            "_", plotlist, ".txt"),
                              stringsAsFactors=FALSE)[,1], "None")
  } else {
    site_list <- c(rownames(data)[order(kds)][1:site_list], "None")
  }
  num.sites <- length(site_list)
  # Get kds for all site-types of the mirna.
  params <- GetSiteKds(mirna, experiment, n_constant, sitelist)
  site.kds <- params$Mean[1 : (nrow(params) - 2)]
  names(site.kds) <- rownames(params)[1 : (nrow(params) - 2)]
  # # Use the 8mer flanks to just get the flank strings.
  flanks.temp <- GetFlankKds(mirna, experiment, n_constant, sitelist, "8mer")
  # Get the flank string names
  flanks <- rownames(flanks.temp)

  # Pre-allocate the matrix with the flanking kds.
  flank.kds <- matrix(NaN,nrow=length(flanks),ncol=length(site_list)-1)
  # Name the rows and columns.
  rownames(flank.kds) <- flanks
  colnames(flank.kds) <- site_list[-length(site_list)]
  # Remove the teporary flank data structure.
  rm(flanks.temp)

  # Fill in the flanking kd matrix for each sitetype and flanking dinucleotide
  # combination that exists.
  for (site in site_list[1:(length(site_list)-1)]) {
    flank.kds.site <- GetFlankKds(mirna, experiment, n_constant, sitelist, site)
    for (flank in rownames(flank.kds.site)) {
      flank.kds[flank,site] <- flank.kds.site[flank,]$Mean
    }
  }
  # Removes site types for which there are literally zero flank Kds.
  site.kds <- site.kds[which(colSums(is.na(flank.kds))!=256)]
  site.kds <- site.kds[site_list[-length(site_list)]]
  flank.kds.trim <- flank.kds[,site_list[-length(site_list)]]
  flank.kds.trim <- flank.kds.trim[,order(site.kds)]
  kds.flanks <- c(flank.kds.trim)
  data.sites <- rep(colnames(flank.kds.trim), each=nrow(flank.kds.trim))
  data.ranks <- rep(num.sites-seq(num.sites-1), each=nrow(flank.kds.trim))
  data.colors <- rep(sapply(rownames(flank.kds), GetColorFunction, alpha=1),
                     ncol(flank.kds.trim))
  flanks.df <- data.frame(kds=log10(kds.flanks), rank=as.numeric(data.ranks), sites=data.sites,
                          cols=data.colors,
                           stringsAsFactors=FALSE)

  print(unique(flanks.df$site))
  flanks.df <<- flanks.df
  dev.new(xpos = xpos, ypos = ypos, height = 5, width = 8)

  par(kPlotParameters)
  boxplot(kds ~ rank,
          data       = flanks.df,
          axes       = FALSE,
          horizontal = TRUE,
          outline    = FALSE,
          xlim       = c(0, 22),
          ylim       = log10(rev(c(0.00003, 10))))
  print(c(0,num.sites+5))
  title(main = mirna,
        line = -2,
        adj  = 0.1)
  title(xlab = expression(italic(K)[D]))
    beeswarm(kds ~ rank,
             data       = flanks.df,
             add        = TRUE,
             method     = "swarm",
             corral     = "random",
             corralWidth = 0.5,
             pch        = 1,
             cex = 0.8,
             horizontal = TRUE,  
             pwcol = cols)




  ymin=0.0001
  ymax=10
  ys <- c(sapply(seq(floor(min(flanks.df$kds, na.rm=TRUE)),
                     ceiling(max(flanks.df$kds, na.rm=TRUE))), function(x) seq(10)*10^x))
  ys <- ys[ys >= ymin & ys <= ymax]
  ymin <- min(ys)
  ymax <- max(ys)
  yl <- seq(ceiling(log10(min(ys))), floor(log10(max(ys))))
  axis(1, at=yl,
       labels=sapply(yl, function(name) {
         eval(substitute(expression(10^x), list(x=name)))
       }),
       pos=0, lwd=0)
  axis(1, at=log10(ys), labels=FALSE,
       pos=0)
  # points(rep(-3.4, num.sites-1),num.sites - seq(num.sites-1),col= kSiteColors[names(site.kds),][order(site.kds)],cex=1.5,pch=19)

  flank_mins <- apply(flank.kds.trim,2,function(col) {min(col[!is.na(col)])})
    text(log10(flank_mins) - 0.2,
         num.sites - seq(num.sites-1),
         labels=ConvertTtoU(unique(flanks.df$sites)),
         adj=0, col= "black", cex = 0.8)

  arrows(log10(max(flank.kds[,"8mer"])), 16.8,
         log10(min(flank.kds[,"8mer"])), 16.8,
         length=0.05, lty=1, angle=90, code=3)
  text(-2.1,
       17.3, "Range of 8mer binding affinity across dinucleotide context", cex = 0.6)

  arrows(log10(site.kds["6mer"]), 18,
        log10(site.kds["8mer"]), 18,
        length=0.05, lty=1, angle=90, code=3)
  text(-2.1,
       18.5, "Difference in average binding betwen 6mer and 8mer site", cex = 0.6)

}

PlotSiteFlankKdsScatter <- function(mirna, experiment, n_constant_1, n_constant_2, sitelist,
                                site, xpos = 20, ypos = 20) {


  flanks.kds.1 <- GetFlankKds(mirna, experiment, n_constant_1, sitelist, "8mer")
  flanks.kds.2 <- GetFlankKds(mirna, experiment, n_constant_2, sitelist, "8mer")

  # Get the flank string names


  dev.new(xpos = xpos, ypos = ypos, height = 5, width = 5)
    par(kPlotParameters)

  plot(flanks.kds.1$Mean,
    flanks.kds.2$Mean,
    xlim = c(0.0001, 5),
    ylim = c(0.0001, 5),
      axes = FALSE,
      col=sapply(rownames(flanks.kds.1), GetColorFunction, alpha=1),
      ann= FALSE,
      log ='xy')
  title(main = mirna,
        line = -2,
        adj  = 0.1)
  title(xlab = expression(italic(K)[D]))


  ymin=0.0001
  ymax=5


  ys <- c(sapply(seq(-5,6), function(x) seq(9)*10^(x-1)))
  ys <- ys[which(ys>=ymin & ys <= ymax)]
  yl <- 10^seq(ceiling(log10(min(ys))), floor(log10(max(ys))))
  print(yl)
  print(ys)
  axis(1, at=yl,
       labels=sapply(yl, function(name) {
         eval(substitute(expression(10^x), list(x=log10(name))))
       }),
       pos=ymin)
  axis(1, at=ys, labels=FALSE,
       pos=ymin)

  # points(rep(-3.4, num.sites-1),num.sites - seq(num.sites-1),col= kSiteColors[names(site.kds),][order(site.kds)],cex=1.5,pch=19)


}

GetCanonicalSiteFlankingKdMatrix <- function(experiment, n_constant, sitelist) {

  flanks.temp <- GetFlankKds("miR-1", experiment, n_constant, sitelist, "8mer")
  flank.strings <- rownames(flanks.temp)

  data.matrix <- matrix(NA, ncol = 258, nrow =30)
  colnames(data.matrix) <- c("miRNA", "site", flank.strings)
  row <- 1
  for (mirna in mirnas_all) {
    for (site in canonical_sites) {
      flanks <- GetFlankKds(mirna, experiment, n_constant, sitelist, site)
      data.matrix[row, 1] <- mirna
      data.matrix[row, 2] <- site
      row_mean <- mean(log10(flanks$Mean))
      for (flank in rownames(flanks)) {
        data.matrix[row,which(colnames(data.matrix)==flank)] <- log10(as.numeric(flanks[flank,]$Mean)) - row_mean
      }
      row <- row + 1
    }
  }
  return(data.matrix)
}

PlotCanonicalSiteFlankingKdHeatmap <- function(experiment, n_constant, sitelist, xpos = 20, ypos = 20) {
  data.matrix <- GetCanonicalSiteFlankingKdMatrix(experiment, n_constant, sitelist)
  data.flanks <- data.frame(data.matrix[,3:258],stringsAsFactors=FALSE)
  rownames(data.flanks) <- apply(data.matrix,1,function(row) {
    paste0(c(row[1],row[2]), collapse=" | ")
    })
  dev.new(xpos = xpos, ypos = ypos, width = 14, height = 6)
  heatmap.2(data.matrix(data.flanks),
    col=rainbow(100,s=0.8,v=1,start=0.45,end=0.9),
    trace = "none",
    scale = "row",
    cexRow = 0.9,
    margins = c(3.5, 8),
    na.color="gray",
    lhei=c(1,10),
    lwid = c(1, 9),
    key=FALSE)

}

GenerateLMForFlankingKds <- function(experiment, n_constant, sitelist) {
  data.matrix <- GetCanonicalSiteFlankingKdMatrix(experiment, n_constant, sitelist)
  print(data.matrix[1:5, 1:5])
  print(c(data.matrix[3:5,3:5]))
}

GetPairingFlankData <- function(mirna, experiment, condition, n_constant,
                                sitelist, site, mir_start, mir_stop,
                                absolute=FALSE, noconstant=TRUE) {
  if (absolute == TRUE) {
    folder = "/structural_analysis_PAPER_realfinal/"
  } else {
    folder = "/structural_analysis_PAPER_final/"
  }
  file = paste0(condition, "_", n_constant, "_", mir_start, "-", mir_stop,
                collapse="")
  if (noconstant == TRUE) {
    file = paste0(file, "_noconstant", collapse="")
  }
  data <- read.table(paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                            experiment, folder, site, "/", file, ".txt"),
                     header=TRUE)
  data <- data.frame(data)
  return(data)
}


MakeFlankScatterPlotSimple <- function(mirna, experiment, condition, n_constant, sitelist, site, mir_start, mir_stop) {
  data <- data.frame(GetPairingFlankData(mirna, experiment, condition, n_constant, sitelist, site, mir_start, mir_stop))

  print(data[1:10,])
  print(data$Flank)
  attach(data)
  par(kPlotParameters)
  plot(data$plfold^(1/(mir_stop - mir_start + 1)),
       data$accessibility,
       col=sapply(unlist(lapply(data$Flank,as.character)), GetColorFunction, alpha=0.1))
  detach(data)
}

MakeFlankScatterPlot <- function(mirna, experiment, condition, n_constant, sitelist, site, mir_start, mir_stop) {
  data <- data.frame(GetPairingFlankData(mirna, experiment, condition, n_constant, sitelist, site, mir_start, mir_stop))

  print(data[1:10,])
  print(data$Flank)
  attach(data)
  par(kPlotParameters)
  plot(by(data,Flank,function(x){mean(x$plfold^(1/(mir_stop - mir_start + 1)))}),
       by(data,Flank,function(x){mean(x$accessibility)}),
       col=sapply(unique(unlist(lapply(data$Flank,as.character))), GetColorFunction, alpha=1),
       log='xy',
       xlim = c(0.03,1),
       ylim = c(0.03,1))
  detach(data)
}

MakeFlankECDFPlot <- function(mirna, experiment, n_constant, sitelist, site,
                              mir_start, mir_stop, absolute=TRUE,
                              noconstant=TRUE, xpos = 20, ypos = 20) {
  dev.new(xpos = xpos, ypos = ypos, height = 5, width = 10)
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
  xs = seq(0, 1, by = 0.2)
  x_ecdf <- seq(0, 1, by = 0.02)

  par(kPlotParameters)
  par(mfrow=c(1,2))
  par(lwd = 1.5)

  legend_func <- function() {
      legend.names <- c("Input library",
                    "0% AG02-miR-1",
                    "0.4% AG02-miR-1",
                    "1.26% AG02-miR-1",
                    "4% AGO2-miR-1",
                    "12.6% AGO2-miR-1",
                    "40% AGO2-miR-1")
  legend(x=0.02,y=1, legend=legend.names, lwd = 1, col=exp_cond_colors,
         bty="n", ncol = 1)

  }
  win = (mir_stop - mir_start + 1)
  # win = 1
  data <- data.frame(GetPairingFlankData(mirna, experiment, "I", n_constant,
                                         sitelist, site, mir_start, mir_stop,
                                         absolute=absolute,
                                         noconstant=noconstant))
  data.intermediate <- data$plfold^(1/win)
  ecdf.I <- ecdf(data.intermediate)

  plot(x = x_ecdf, y = ecdf.I(x_ecdf),
       type = "l",
       xlim = c(0,1),
       axes = FALSE,
       ann = FALSE)

  sapply(c("0", "0.4", "1.26", "4", "12.6", "40"), function(condition) {
    data <- data.frame(GetPairingFlankData(mirna, experiment, condition,
                                           n_constant, sitelist, site,
                                           mir_start, mir_stop,
                                           absolute=absolute,
                                           noconstant=noconstant))
    data.intermediate <- data$plfold^(1/win)
    ecdf.d <- ecdf(data.intermediate)
    lines(x_ecdf, ecdf.d(x_ecdf), col = exp_cond_colors[condition])
  })

    axis(1, at=xs,
         pos=ymin)
    # Label the axis at each order of magnitude.
    axis(2, at=xs,
         pos=ymin)

  title(xlab = "Per-nucleotide accessibility across from miRNA nt 1–15")
  title(ylab = "ECDF")

  # text(x=0.15, y=0.95, mirna)

  legend_func()

  data <- data.frame(GetPairingFlankData(mirna, experiment, "I", n_constant,
                                         sitelist, site, mir_start, mir_stop,
                                         absolute=absolute,
                                         noconstant=noconstant))
  ecdf.I <- ecdf(data$accessibility^((mir_stop - mir_start + 1) /win))
  plot(x = x_ecdf, y = ecdf.I(x_ecdf),
       type = "l",
       xlim = c(0,1),
       axes = FALSE,
       ann = FALSE)

  sapply(c("0", "0.4", "1.26", "4", "12.6", "40"), function(condition) {
    data   <- data.frame(GetPairingFlankData(mirna, experiment, condition,
                                             n_constant, sitelist, site,
                                             mir_start, mir_stop,
                                             absolute=absolute,
                                             noconstant=noconstant))
    ecdf.d <- ecdf(data$accessibility^((mir_stop - mir_start + 1) /win))
    lines(x_ecdf, ecdf.d(x_ecdf), col = exp_cond_colors[condition])
  })

    axis(1, at=xs,
         pos=ymin)
    # Label the axis at each order of magnitude.
    axis(2, at=xs,
         pos=ymin)

  title(xlab = "Average nucleotide accessibility across from miRNA nt 1–15")
  title(ylab = "ECDF")

  legend_func()

}

MakeFlankECDFPlot_AU <- function(mirna, experiment, n_constant, sitelist, site,
                                 mir_start, mir_stop, absolute=TRUE,
                                 noconstant=TRUE, AU_score="AU_cs", xpos = 20,
                                 ypos = 20) {
  dev.new(xpos = xpos, ypos = ypos, height = 5, width = 10)
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
  xs = seq(0, 1, by = 0.2)
  x_ecdf <- seq(0, 1, by = 0.01)

  par(kPlotParameters)
  par(mfrow=c(1,2))
  par(lwd = 1.5)

  legend_func <- function() {
      legend.names <- c("Input library",
                    "0% AG02-miR-1",
                    "0.4% AG02-miR-1",
                    "1.26% AG02-miR-1",
                    "4% AGO2-miR-1",
                    "12.6% AGO2-miR-1",
                    "40% AGO2-miR-1")
  legend(x=0.02,y=1, legend=legend.names, lwd = 1, col=exp_cond_colors,
         bty="n", ncol = 1)

  }
  win = (mir_stop - mir_start + 1)
  # win = 1
  data <- data.frame(GetPairingFlankData(mirna, experiment, "I", n_constant,
                                         sitelist, site, mir_start, mir_stop))
  data.intermediate <- data[[AU_score]]
  ecdf.I <- ecdf(data.intermediate)

  plot(x = x_ecdf, y = ecdf.I(x_ecdf),
       type = "l",
       xlim = c(0,1),
       axes = FALSE,
       ann = FALSE)

  sapply(c("0", "0.4", "1.26", "4", "12.6", "40"), function(condition) {
    data <- data.frame(GetPairingFlankData(mirna, experiment, condition,
                                           n_constant, sitelist, site,
                                           mir_start, mir_stop))
    data.intermediate <- data[[AU_score]]
    ecdf.d <- ecdf(data.intermediate)
    lines(x_ecdf, ecdf.d(x_ecdf), col = exp_cond_colors[condition])
  })

    axis(1, at=xs,
         pos=ymin)
    # Label the axis at each order of magnitude.
    axis(2, at=xs,
         pos=ymin)

  title(xlab = "Per-nucleotide accessibility across from miRNA nt 1–15")
  title(ylab = "ECDF")

  # text(x=0.15, y=0.95, mirna)

  legend_func()

}


MakeKdVsAccessibilityPlot <- function(mirna, experiment, n_constant,
                                      sitelist, site, mir_start, mir_stop,
                                      xpos = 20, ypos = 20, absolute=TRUE,
                                      noconstant=TRUE) {
  dev.new(xpos = xpos, ypos = ypos, width = 12, height = 5)
  par(kPlotParameters)
  par(mfrow = c(2, 6))
  # Window size for normalizing the pl_fold with the accessibility:
  win = 1/(mir_stop - mir_start + 1)

  # Get the flanking kds:
  flank.kds <- GetFlankKds(mirna, experiment, n_constant, sitelist, site)
  kds <- flank.kds$Mean
  names(kds) <- sapply(rownames(flank.kds), function(name) {
    paste0(unlist(strsplit(name, split = ""))[-3], collapse = "")
    })

  plot_subfunction <- function(condition,factor,ylimit, logparam) {
        data <- data.frame(GetPairingFlankData(mirna, experiment, condition,
                           n_constant, sitelist, site, mir_start, mir_stop,
                           absolute=absolute, noconstant=noconstant))
    print(data[1 : 10, ])
    attach(data)

    p_access <- c(by(data, Flank, function(x) {
      print(mean(log(x[[factor]])))
        if (factor == "plfold"){
          10^(mean(log(x[[factor]])))
          } else if (factor == "accessibility") {
          10^(mean(log(x[[factor]])))
          } else {
          mean(x[[factor]])
          }
      }
    ))

    detach(data)

    plot(1/kds,
         p_access[names(kds)],
         col=sapply(names(kds), GetColorFunction, alpha=0.1),
         log = logparam,
         cex = 0.8,
         xlim = c(5,3000),
         ylim = ylimit)
    text(x=20, y = 0.95, condition)
    if (logparam == "xy") {
      cor_text <- round(cor(-log(kds),log(p_access[names(kds)]),
                                use = "pairwise.complete.obs")^2,2)      
    } else {
      cor_text <- round(cor(-log(kds),p_access[names(kds)],
                                use = "pairwise.complete.obs")^2,2)      

    }
    text(x = 500, y = 0.95, bquote(italic(r)^2 ==  .(cor_text)))

  }

  sapply(c("I_combined", "0.4", "1.26", "4","12.6", "40"), plot_subfunction, factor = "plfold", ylimit = 10^c(-0.9,0),logparam='xy')
  sapply(c("I_combined", "0.4", "1.26", "4","12.6", "40"), plot_subfunction, factor = "accessibility", ylimit = c(0.001, 1),logparam='xy')

}

MakeAccessibilityvsInputPlot <- function(mirna, experiment, condition, n_constant,
                                      sitelist, site, mir_start, mir_stop,pl=TRUE) {
  dev.new(xpos = 20, ypos = 220, height = 5, width = 6)
  par(kPlotParameters)
  # Window size for normalizing the pl_fold with the accessibility:
  win = 1/(mir_stop - mir_start + 1)

  # Get the flanking kds:
  flank.kds <- GetFlankKds(mirna, experiment, n_constant, sitelist, site)
  kds <- flank.kds$Mean
  names(kds) <- sapply(rownames(flank.kds), function(name) {
    paste0(unlist(strsplit(name, split = ""))[-3], collapse = "")
    })
    subfunction <- function(condition_temp) {
    data <- data.frame(GetPairingFlankData(mirna, experiment, condition_temp,
                       n_constant, sitelist, site, mir_start, mir_stop))
    attach(data)

    p_access<- c(by(data, Flank, function(x) {
          10^(mean(log(x$plfold^win)))
      }
    ))
    detach(data)

    return(p_access)
  }
    p_I <- subfunction("I")
    p_condition <- subfunction(condition)
    plot(p_I,
         p_condition[names(p_I)],
         col=sapply(names(p_I), GetColorFunction, alpha=0.1),
         log = 'xy',
         xlim = c(0.1, 1),
         ylim = c(0.1, 1))
  

 

}

PlotPlfoldByAUScoreBin <- function(mirna, experiment, condition, n_constant,
                                 sitelist, site, mir_start, mir_stop, xpos=20,
                                 ypos=20, num.bins=11,AU_score="AU_win",p_score="plfold",
                                 noconstant=FALSE,absolute=TRUE) {
  dev.new(xpos=xpos, ypos=ypos, width=5, height=5)
  par(kPlotParameters)


  plot(1, type = "n",
         xlim = c(0, 1),
         ylim = c(0, 1))

  bin.limits <- seq(0,1, length = num.bins)
  print(bin.limits)
  b.l <- length(bin.limits)
  x.values <- (bin.limits[2 : b.l] + bin.limits[1 : (b.l - 1)]) / 2
  print(x.values)
  # Window size for normalizing the pl_fold with the accessibility:
  win = 1/(mir_stop - mir_start + 1)
  subfunction <- function(condition, factor, cond.color = "black") {
    data <- data.frame(GetPairingFlankData(mirna, experiment, condition,
                       n_constant, sitelist, site, mir_start, mir_stop,
                       noconstant=noconstant, absolute=absolute), stringsAsFactors=FALSE)

    data$AU_bins <- as.integer(data[[factor]]*num.bins)
    print(unique(data$AU_bins))
    attach(data)
      y.values <- c(by(data, AU_bins, function(x) {
        mean(x[[p_score]])
        }))
      y.error <- c(by(data, AU_bins, function(x) {
        sd(x[[p_score]])/sqrt(length(x[[p_score]])-1)
      }))
    detach(data)

      agg <- aggregate(. ~ AU_bins, data, function(x) {
      c(mean = mean(x),
        se = sd(x),
        len = length(x))
      })

    print(agg)

    x.values <- agg$AU_bins/num.bins
    y.values <- agg[[p_score]][,1]
    y.error <- agg[[p_score]][,2]/sqrt(agg[[p_score]][,3]-1)
    print(x.values)
    points(x.values,
         y.values,
         type = "o",
         col=cond.color)
    arrows(x.values, y.values-y.error, x.values, y.values + y.error,
           col=cond.color, length=0.07, angle=90, code=3)
  }



    text(x=20, y = 0.95, condition)

    subfunction("I", AU_score)
    subfunction("4", AU_score,cond.color="red")

}

PlotAUScoreByPlfold <- function(mirna, experiment, condition, n_constant,
                                 sitelist, site, mir_start, mir_stop, noconstant=FALSE,
                                 absolute=TRUE, xpos=20,
                                 ypos=20, num.bins=11,AU_score="AU_win", p_score="plfold") {
  dev.new(xpos=xpos, ypos=ypos, width=5, height=5)
  par(kPlotParameters)


  plot(1, type = "n",
         xlim = c(0.0, 1),
         ylim = c(0, 1))

  bin.limits <- seq(0.2,1, length = num.bins)
  print(bin.limits)
  b.l <- length(bin.limits)
  x.values <- (bin.limits[2 : b.l] + bin.limits[1 : (b.l - 1)]) / 2
  print(x.values)
  # Window size for normalizing the pl_fold with the accessibility:
  win = 1/(mir_stop - mir_start + 1)
  subfunction <- function(condition, factor, cond.color = "black") {
    data <- data.frame(GetPairingFlankData(mirna, experiment, condition,
                       n_constant, sitelist, site, mir_start, mir_stop,
                       noconstant=noconstant, absolute=absolute))

    print(data[1:10,])
    print("done")
    x.values = as.integer(data[[p_score]]*num.bins)
    data$plfold_bins = as.integer(data[[p_score]]*num.bins)


    print(data[1:10,])
    agg <- aggregate(. ~ plfold_bins, data, function(x) {
      c(mean = mean(x),
        se = sd(x),
        len = length(x))
      })

    attach(data)
      y.values <- c(by(data, plfold_bins, function(x) {
        mean(x[[AU_score]])
        }))
      y.error <- c(by(data, plfold_bins, function(x) {
        sd(x[[AU_score]])/sqrt(length(x[[AU_score]])-1)
      }))
    detach(data)
    print(agg)
    # print(y.error)

    x <- agg$plfold_bins/num.bins
    y <- agg[[AU_score]][,1]
    y.e <- agg[[AU_score]][,2]/sqrt(agg[[AU_score]][,3]-1)
    print("length:")
    print(agg$AU_cs[,3])

    print(cbind(y.values, y))
    print(x)
    print(y)
    # print(y_error)
    points(x, y,
         type = "o",
         col=cond.color)
    arrows(x, y - y.e, x, y + y.e,
           col=cond.color, length=0.07, angle=90, code=3)

    text(x=0.1,y=0.9, p_score)
    text(x=0.1,y=0.85, AU_score)
  }



    text(x=20, y = 0.95, condition)

    subfunction("I", AU_score)
    subfunction(condition, AU_score,cond.color="red")

}


# graphics.off()
# MakeFlankECDFPlot("miR-1", "equilibrium", -3, "paper", "8mer", 1, 15, xpos = 20, ypos =20, noconstant=FALSE)
# MakeFlankECDFPlot_AU("miR-1", "equilibrium", 5, "paper", "8mer", 1, 15, xpos = 20, ypos =20)
# MakeFlankECDFPlot_AU("miR-1", "equilibrium", 5, "paper", "8mer", 1, 15, xpos = 20, ypos =20, AU_score="AU_win")
# MakeFlankECDFPlot_AU("miR-1", "equilibrium", 5, "paper", "8mer", 1, 15, xpos = 20, ypos =20, AU_score = "AU_read")
# MakeKdVsAccessibilityPlot("miR-1", "equilibrium", 5, "paper", "8mer", 1, 15, xpos = 20, ypos = 220, absolute=TRUE,noconstant=FALSE)
graphics.off()
# MakeKdVsAccessibilityPlot("miR-1", "equilibrium", -3, "paper", "8mer", 1, 15, xpos = 20, ypos = 220, absolute=TRUE,noconstant=FALSE)
# PlotPlfoldByAUScoreBin("miR-1", "equilibrium", 0.4, 5, "paper", "8mer", 1, 15, xpos = 420, ypos = 20, num.bins=12)
# PlotPlfoldByAUScoreBin("miR-1", "equilibrium", 4, -3, "paper", "8mer", 1, 15, xpos = 20, ypos = 420, num.bins=11, p_score="plfold", AU_score ="AU_win_wo")
# PlotPlfoldByAUScoreBin("miR-1", "equilibrium", 4, -3, "paper", "8mer", 1, 15, xpos = 420, ypos = 420, num.bins=11, p_score="accessibility",AU_score ="AU_win_wo")

# PlotPlfoldByAUScoreBin("miR-1", "equilibrium", 4, -3, "paper", "8mer", 1, 15, xpos = 20, ypos = 420, num.bins=11, p_score="plfold", AU_score ="AU_read")
# PlotPlfoldByAUScoreBin("miR-1", "equilibrium", 4, -3, "paper", "8mer", 1, 15, xpos = 420, ypos = 420, num.bins=11, p_score="accessibility",AU_score ="AU_read")


# PlotPlfoldByAUScoreBin("miR-1", "equilibrium", 4, -3, "paper", "8mer", 1, 15, xpos = 20, ypos = 420, num.bins=11, p_score="plfold", AU_score ="AU_read_wo")
PlotPlfoldByAUScoreBin("miR-1", "equilibrium", 4, -3, "paper", "8mer", 1, 15, xpos = 20, ypos = 20, num.bins=16, p_score="accessibility",AU_score ="AU_win")
PlotPlfoldByAUScoreBin("miR-1", "equilibrium", 4, -3, "paper", "8mer", 1, 15, xpos = 420, ypos = 20, num.bins=38, p_score="accessibility",AU_score ="AU_read")

PlotAUScoreByPlfold("miR-1", "equilibrium", 4, -3, "paper", "8mer", 1, 15, xpos = 20, ypos = 420, num.bins=11, AU_score = "AU_read", p_score="accessibility")
PlotAUScoreByPlfold("miR-1", "equilibrium", 4, -3, "paper", "8mer", 1, 15, xpos = 420, ypos = 420, num.bins=11, AU_score ="AU_win", p_score="accessibility")
# PlotAUScoreByPlfold("miR-1", "equilibrium", 4, -3, "paper", "8mer", 1, 15, xpos = 820, ypos = 20, num.bins=11, AU_score ="AU_read_wo",p_score="accessibility")
# PlotAUScoreByPlfold("miR-1", "equilibrium", 4, -3, "paper", "8mer", 1, 15, xpos = 1220, ypos = 20, num.bins=11, AU_score ="AU_win_wo",p_score="accessibility")

# PlotAUScoreByPlfold("miR-1", "equilibrium", 4, -3, "paper", "8mer", 1, 15, xpos = 20, ypos = 620, num.bins=11, AU_score = "AU_read", p_score="plfold")
# PlotAUScoreByPlfold("miR-1", "equilibrium", 4, -3, "paper", "8mer", 1, 15, xpos = 420, ypos = 620, num.bins=11, AU_score ="AU_win", p_score="plfold")
# PlotAUScoreByPlfold("miR-1", "equilibrium", 4, -3, "paper", "8mer", 1, 15, xpos = 820, ypos = 620, num.bins=11, AU_score ="AU_read_wo",p_score="plfold")
# PlotAUScoreByPlfold("miR-1", "equilibrium", 4, -3, "paper", "8mer", 1, 15, xpos = 1220, ypos = 620, num.bins=11, AU_score ="AU_win_wo",p_score="plfold")


# PlotAUScoreByPlfold("miR-1", "equilibrium", 4, -3, "paper", "8mer", 1, 15, xpos = 20, ypos = 520, num.bins=11, absolute=TRUE, noconstant=FALSE,p_score="accessibility")
# PlotAUScoreByPlfold("miR-1", "equilibrium", 0.4, 5, "paper", "8mer", 1, 15, xpos = 320, ypos = 520, num.bins=10, AU_score ="AU_read", p_score="accessibility")
# PlotAUScoreByPlfold("miR-1", "equilibrium", 0.4, 5, "paper", "8mer", 1, 15, xpos = 620, ypos = 520, num.bins=10, AU_score ="AU_cs", p_score="accessibility")



# MakeAccessibilityvsInputPlot("miR-1", "equilibrium", 40, 5, "paper", "8mer", 1, 15)
# MakeAccessibilityvsInputPlot("miR-1", "equilibrium", 4, 5, "paper", "8mer", 1, 15)
# MakeAccessibilityvsInputPlot("miR-1", "equilibrium", 1.26, 5, "paper", "8mer", 1, 15)

# FOR PAPER
PlotpairwiseKds <- function(mirna, experiment, start, stop, sitelist, combined.input=TRUE) {
  kds.1 <- GetKds(mirna, experiment, start, stop, sitelist, combined.input = combined.input)
  kds.2 <- GetKds(mirna, experiment, start, stop, sitelist, log.residual = TRUE, combined.input=combined.input)
  par(kPlotParameters)
    xmin <- 10^floor(min(log10(kds.1), log10(kds.2)))
    xmax <- 10^ceiling(max(log10(as.numeric(unlist(c(kds.1, kds.2))))))
    xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
    xs <- xs[xs >= xmin & xs <= xmax]
    xmin <- min(xs)
    xmax <- max(xs)
    xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
    plot(kds.1,
      kds.2,
      col = kSiteColors[c(names(kds.1), "bg", "Ago"),],
      pch = c(rep(20, length(kds.1) - 2), 1, 1),
      log = 'xy',
      xlim = c(xmin,xmax*(xmax/xmin)^0.55), 
      ylim = c(xmin, xmax),
      axes = FALSE,
      ann = FALSE)
    segments(xmin, xmin, xmax, xmax, lty = 2)
    axis(1, at=xl,
         labels=sapply(xl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=xmin, lwd=0, cex.axis=1.5, padj=0.2)
    axis(1, at=xs, labels=FALSE,
         pos=xmin, lwd=2)
    # Label the axis at each order of magnitude.

    axis(2, at=xl,
         labels=sapply(xl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
    axis(2, at=xs, labels=FALSE,
         pos=xmin, lwd=2)

        title(main = mirna, font.main=1, cex.main=1.5, line=-2, adj=0.1)
    title(xlab = paste0("Multinomial optimized parameters (nM)"), cex.lab=1.5, line=2, adj=0.3)
    title(ylab = paste0("Log-normal optimized parameters (nM)"), cex.lab=1.5, line=2)
    legend.names <- c(rownames(sitesXcounts), "bg", "Ago")
    # centered_sites <- c("11mer-m3.13", "12mer-m3.14",
    #                 "11mer-m4.14", "12mer-m4.15")
    # centered_index <- sapply(centered_sites, function(site){
    #                        which(legend.names == site)
    #                      }
    #                      )
    # legend.names[centered_index] <- "Centered"
    # legend.names <- unique(legend.names)
    legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=0.9, bty="n", ncol = 1)

}

PlotpairwiseKdsInput <- function(mirna, experiment, start, stop, sitelist, log.residual=FALSE) {
  kds.1 <- GetKds(mirna, experiment, start, stop, sitelist, log.residual=log.residual)
  kds.2 <- GetKds(mirna, experiment, start, stop, sitelist, combined.input = FALSE, log.residual=log.residual)
  par(kPlotParameters)
    xmin <- 10^floor(min(log10(kds.1), log10(kds.2)))
    xmax <- 10^ceiling(max(log10(as.numeric(unlist(c(kds.1, kds.2))))))
    xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
    xs <- xs[xs >= xmin & xs <= xmax]
    xmin <- min(xs)
    xmax <- max(xs)
    xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))


    plot(kds.1,
      kds.2,
      col = kSiteColors[c(names(kds.1), "bg", "Ago"),],
      pch = c(rep(20, length(kds.1) - 2), 1, 1),
      log = 'xy',
      xlim = c(xmin,xmax*(xmax/xmin)^0.55), 
      ylim = c(xmin, xmax),
      axes = FALSE,
      ann = FALSE)
    segments(xmin, xmin, xmax, xmax, lty = 2)
    axis(1, at=xl,
         labels=sapply(xl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=xmin, lwd=0, cex.axis=1.5, padj=0.2)
    axis(1, at=xs, labels=FALSE,
         pos=xmin, lwd=2)
    # Label the axis at each order of magnitude.

    axis(2, at=xl,
         labels=sapply(xl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
    axis(2, at=xs, labels=FALSE,
         pos=xmin, lwd=2)

        title(main = mirna, font.main=1, cex.main=1.5, line=-2, adj=0.1)
    title(xlab = paste0("Combined input parameters (nM)"), cex.lab=1.5, line=2, adj=0.3)
    title(ylab = paste0("Single input parameters (nM)"), cex.lab=1.5, line=2)
    legend.names <- c(rownames(sitesXcounts), "bg", "Ago")
    # centered_sites <- c("11mer-m3.13", "12mer-m3.14",
    #                 "11mer-m4.14", "12mer-m4.15")
    # centered_index <- sapply(centered_sites, function(site){
    #                        which(legend.names == site)
    #                      }
    #                      )
    # legend.names[centered_index] <- "Centered"
    # legend.names <- unique(legend.names)
    legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=0.9, bty="n", ncol = 1)

}

PlotpairwiseFlankKds <- function(mirna, experiment, start, stop, site,
                                 sitelist, combined.input=TRUE) {
  kds.1 <- GetFlankKds(mirna, experiment, start, stop, site, sitelist,
                       combined.input=combined.input)
  kds.2 <- GetFlankKds(mirna, experiment, start, stop, site, sitelist,
                       combined.input=combined.input, log.residual=TRUE)
  par(kPlotParameters)
    xmin <- 10^floor(min(log10(kds.1), log10(kds.2)))
    xmax <- 10^ceiling(max(log10(as.numeric(unlist(c(kds.1, kds.2))))))
    xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
    xs <- xs[xs >= xmin & xs <= xmax]
    xmin <- min(xs)
    xmax <- max(xs)
    xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))


    plot(kds.1,
      kds.2,
      col = sapply(names(kds.1), GetColorFunction),
      pch = 20,
      log = 'xy',
      xlim = c(xmin,xmax*(xmax/xmin)^0.55), 
      ylim = c(xmin, xmax),
      axes = FALSE,
      ann = FALSE)
    segments(xmin, xmin, xmax, xmax, lty = 2)
    axis(1, at=xl,
         labels=sapply(xl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=xmin, lwd=0, cex.axis=1.5, padj=0.2)
    axis(1, at=xs, labels=FALSE,
         pos=xmin, lwd=2)
    # Label the axis at each order of magnitude.

    axis(2, at=xl,
         labels=sapply(xl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
    axis(2, at=xs, labels=FALSE,
         pos=xmin, lwd=2)

        title(main = mirna, font.main=1, cex.main=1.5, line=-2, adj=0.1)
    title(xlab = paste0("Multinomial optimized parameters (nM)"), cex.lab=1.5, line=2, adj=0.3)
    title(ylab = paste0("Log-normal optimized parameters (nM)"), cex.lab=1.5, line=2)
    legend.names <- c(rownames(sitesXcounts), "bg", "Ago")
    # centered_sites <- c("11mer-m3.13", "12mer-m3.14",
    #                 "11mer-m4.14", "12mer-m4.15")
    # centered_index <- sapply(centered_sites, function(site){
    #                        which(legend.names == site)
    #                      }
    #                      )
    # legend.names[centered_index] <- "Centered"
    # legend.names <- unique(legend.names)
    legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=0.9, bty="n", ncol = 1)

}


PlotpairwiseFlankKdsInput <- function(mirna, experiment, start, stop, site,
                                      sitelist, log.residual=FALSE) {
  print("622")
  kds.1 <- GetFlankKds(mirna, experiment, start, stop, site, sitelist, log.residual=log.residual)
  print("624")
  kds.2 <- GetFlankKds(mirna, experiment, start, stop, site, sitelist, combined.input = FALSE, log.residual=log.residual)
  print(kds.1)
  print(kds.2)
  par(kPlotParameters)
    xmin <- 10^floor(min(log10(kds.1), log10(kds.2)))
    xmax <- 10^ceiling(max(log10(as.numeric(unlist(c(kds.1, kds.2))))))
    xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
    xs <- xs[xs >= xmin & xs <= xmax]
    xmin <- min(xs)
    xmax <- max(xs)
    xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))


    plot(kds.1,
      kds.2,
      col = sapply(names(kds.1), GetColorFunction),
      pch = 20,
      log = 'xy',
      xlim = c(xmin,xmax*(xmax/xmin)^0.55), 
      ylim = c(xmin, xmax),
      axes = FALSE,
      ann = FALSE)
    segments(xmin, xmin, xmax, xmax, lty = 2)
    axis(1, at=xl,
         labels=sapply(xl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=xmin, lwd=0, cex.axis=1.5, padj=0.2)
    axis(1, at=xs, labels=FALSE,
         pos=xmin, lwd=2)
    # Label the axis at each order of magnitude.

    axis(2, at=xl,
         labels=sapply(xl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
    axis(2, at=xs, labels=FALSE,
         pos=xmin, lwd=2)

        title(main = mirna, font.main=1, cex.main=1.5, line=-2, adj=0.1)
    title(xlab = paste0("Combined input parameters (nM)"), cex.lab=1.5, line=2, adj=0.3)
    title(ylab = paste0("Single input parameters (nM)"), cex.lab=1.5, line=2)
    legend.names <- c(rownames(sitesXcounts), "bg", "Ago")
    # centered_sites <- c("11mer-m3.13", "12mer-m3.14",
    #                 "11mer-m4.14", "12mer-m4.15")
    # centered_index <- sapply(centered_sites, function(site){
    #                        which(legend.names == site)
    #                      }
    #                      )
    # legend.names[centered_index] <- "Centered"
    # legend.names <- unique(legend.names)
    legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=0.9, bty="n", ncol = 1)

}





PlotpairwiseFlankKdsSites <- function(mirna, experiment, start, stop, site1, site2, sitelist, combined.input=TRUE, log.residual=FALSE,position=FALSE,nosite=TRUE) {
  kds.1 <- GetFlankKds(mirna, experiment, site1, start, stop, sitelist, combined.input = combined.input, log.residual=log.residual)
  kds.2 <- GetFlankKds(mirna, experiment, site2, start, stop, sitelist, combined.input = combined.input, log.residual=log.residual)
 par(kPlotParameters)
      inds <- intersect(names(kds.1), names(kds.2))
      kds.1 <- kds.1[inds]
      kds.2 <- kds.2[inds]
    xmin <- 10^(min(log10(kds.1)) - 0.1)
    xmax <- 10^(max(log10(kds.1)) + 0.1)
    ymin <- 10^(min(log10(kds.2)) - 0.1)
    ymax <- 10^(max(log10(kds.2)) + 0.1)


    xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
    ys <- c(sapply(seq(floor(log10(ymin)), ceiling(log10(ymax))), function(x) seq(10)*10^x))

    xs <- xs[xs >= xmin & xs <= xmax]
    ys <- ys[ys >= ymin & ys <= ymax]

    # xmin <- min(xs)
    # xmax <- max(xs)
    # ymin <- min(ys)
    # ymax <- max(ys)
    xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
    yl <- 10^seq(ceiling(log10(min(ys))), floor(log10(max(ys))))
    if (position != FALSE) {
      cex1 <- 1
      pos1 <- position
    } else {
      cex1 <- 3.5
      pos1 <- 1
    }
    plot(kds.1,
      kds.2,
      col = kNucleotideColors[GetSingleFlankPosition(names(kds.1), pos1)],
      cex = cex1,
      pch = 19,
      log = 'xy',
      xlim = c(xmin,xmax*(xmax/xmin)^0.55), 
      ylim = c(ymin, ymax),
      axes = FALSE,
      ann = FALSE)
    if (position == FALSE){
      points(kds.1,  kds.2, pch = 19, cex = 2.5,
             col = kNucleotideColors[GetSingleFlankPosition(names(kds.1), 2)])
      points(kds.1,  kds.2, pch = 19, cex = 1.5,
             col = kNucleotideColors[GetSingleFlankPosition(names(kds.1), 3)])
      points(kds.1,  kds.2, pch = 19, cex = 0.5,
             col = kNucleotideColors[GetSingleFlankPosition(names(kds.1), 4)])
    }
    plot_min <- max(xmin, ymin)
    plot_max <- min(xmax, ymax)
    segments(plot_min, plot_min, plot_max, plot_max, lty = 2)
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

        title(mirna, font.main = 1, cex.main = 1.5, line=-2, adj=0.1)
        cor_text <- round(
                      cor(log(kds.1[inds]), log(kds.2[inds])),
                      digits = 3
                    )

        title(main = bquote(italic(r) ==  .(cor_text)), font.main=1, cex.main=1.5, line=-3.5, adj=0.1)

    title(xlab = paste0(site1," parameters (nM)"), cex.lab=1.5, line=2, adj=0.3)
    title(ylab = paste0(site2, " parameters (nM)"), cex.lab=1.5, line=2)
    legend.names <- c(rownames(sitesXcounts), "bg", "Ago")
    # centered_sites <- c("11mer-m3.13", "12mer-m3.14",
    #                 "11mer-m4.14", "12mer-m4.15")
    # centered_index <- sapply(centered_sites, function(site){
    #                        which(legend.names == site)
    #                      }
    #                      )
    # legend.names[centered_index] <- "Centered"
    # legend.names <- unique(legend.names)
    legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=0.9, bty="n", ncol = 1)

}



PlotpairwiseFlankKdsmiRNAs <- function(mirna1, mirna2, experiment, start, stop, site, sitelist, combined.input=TRUE, log.residual=FALSE) {
  kds.1 <- GetFlankKds(mirna1, experiment, start, stop, site, sitelist, combined.input = combined.input, log.residual=log.residual)
  kds.2 <- GetFlankKds(mirna2, experiment, start, stop, site, sitelist, combined.input = combined.input, log.residual=log.residual)
  par(kPlotParameters)
      inds <- intersect(names(kds.1), names(kds.2))
      kds.1 <- kds.1[inds]
      kds.2 <- kds.2[inds]
    xmin <- 10^(min(log10(kds.1)) - 0.1)
    xmax <- 10^(max(log10(kds.1)) + 0.1)
    ymin <- 10^(min(log10(kds.2)) - 0.1)
    ymax <- 10^(max(log10(kds.2)) + 0.1)


    xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
    ys <- c(sapply(seq(floor(log10(ymin)), ceiling(log10(ymax))), function(x) seq(10)*10^x))

    xs <- xs[xs >= xmin & xs <= xmax]
    ys <- ys[ys >= ymin & ys <= ymax]

    # xmin <- min(xs)
    # xmax <- max(xs)
    # ymin <- min(ys)
    # ymax <- max(ys)
    xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
    yl <- 10^seq(ceiling(log10(min(ys))), floor(log10(max(ys))))
    plot(kds.1,
      kds.2,
      col = kNucleotideColors[GetSingleFlankPosition(names(kds.1), 1)],
      cex = 3.5,
      pch = 19,
      log = 'xy',
      xlim = c(xmin,xmax*(xmax/xmin)^0.55), 
      ylim = c(ymin, ymax),
      axes = FALSE,
      ann = FALSE)
      points(kds.1,  kds.2, pch = 19, cex = 2.5, col = kNucleotideColors[GetSingleFlankPosition(names(kds.1), 2)])
      points(kds.1,  kds.2, pch = 19, cex = 1.5, col = kNucleotideColors[GetSingleFlankPosition(names(kds.1), 3)])
      points(kds.1,  kds.2, pch = 19, cex = 0.5, col = kNucleotideColors[GetSingleFlankPosition(names(kds.1), 4)])
    plot_min <- max(xmin, ymin)
    plot_max <- min(xmax, ymax)
    segments(plot_min, plot_min, plot_max, plot_max, lty = 2)
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

        title(site, font.main = 1, cex.main = 1.5, line=-2, adj=0.1)
        cor_text <- round(
                      cor(log(kds.1[inds]), log(kds.2[inds])),
                      digits = 3
                    )

        title(main = bquote(italic(r) ==  .(cor_text)), font.main=1, cex.main=1.5, line=-3.5, adj=0.1)

    title(xlab = paste0(mirna1," parameters (nM)"), cex.lab=1.5, line=2, adj=0.3)
    title(ylab = paste0(mirna2, " parameters (nM)"), cex.lab=1.5, line=2)
    legend.names <- c(rownames(sitesXcounts), "bg", "Ago")
    # centered_sites <- c("11mer-m3.13", "12mer-m3.14",
    #                 "11mer-m4.14", "12mer-m4.15")
    # centered_index <- sapply(centered_sites, function(site){
    #                        which(legend.names == site)
    #                      }
    #                      )
    # legend.names[centered_index] <- "Centered"
    # legend.names <- unique(legend.names)
    legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=0.9, bty="n", ncol = 1)

}



GetAllMirnaFlankKds <- function(mirna, experiment, start, stop, sitelist,
                          log.residual=FALSE, combined.input=TRUE, 
                          site_list=NULL, nosite=TRUE) {
  # Get Kds for the miRNA.
  kds <- GetKds(mirna, experiment, start, stop, sitelist,
                log.residual=log.residual, combined.input=combined.input,
                nosite=nosite)
  # Get the 8mer flanks, so as to populate the flank names as row names.
  kds.8mer <- GetFlankKds(mirna, experiment, "8mer", start, stop, sitelist,
                          log.residual=log.residual,
                          combined.input=combined.input, nosite=nosite)
  print(kds.8mer)
  flanks <- names(kds.8mer)
  # Pre-allocate the matrix with "NA".
  flanks.all <- matrix(NA, nrow=length(kds.8mer), ncol=length(kds) - 3)

  rownames(flanks.all) <- names(kds.8mer)
  colnames(flanks.all) <- names(kds)[1:(length(kds) - 3)]

  # Iterate over all the kds names. Exclude the last three since this is
  # the background parameter, the Ago concentration vector, and the
  # loglikelihood or sum of squares for the log-transformed fit.
  sapply(names(kds)[1:(length(kds) - 3)], function(site) {
    kds_flanks <- try(GetFlankKds(mirna, experiment, site, start, stop, sitelist, log.residual=log.residual, combined.input=combined.input))
    if (length(kds_flanks) > 0) {
    flanks.all[names(kds_flanks), site] <<- as.numeric(kds_flanks)
  }
    })

  flanks.new <- matrix(mapply(flanks.all, FUN=as.numeric),
                       ncol=ncol(flanks.all), nrow=nrow(flanks.all))

  rownames(flanks.new) <- rownames(flanks.all)
  colnames(flanks.new) <- colnames(flanks.all)

  return(flanks.new)
}



MakeFigure1 <- function() {
  graphics.off()
  PlotSiteScatterWithInput("miR-1", "equilibrium", 5, "canonical", 7, xpos = 20, ypos = 20)
  dev.copy2eps(file = "2017_Paper/Figure1B_raw_v6.eps")
  PlotSiteEnrichments("miR-1", "equilibrium", 5, "canonical", "canonical", xpos = 620, ypos = 20)
  dev.copy2eps(file = "2017_Paper/Figure1C_raw_v6.eps")
  PlotSiteOccupancy("miR-1", "equilibrium", 5, "canonical", "canonical", xpos = 20, ypos = 420)
  dev.copy2eps(file = "2017_Paper/Figure1D_raw_v6.eps")
  PlotSiteEnrichments("miR-1", "equilibrium", 5, "paper", "paper", xpos = 620, ypos = 420)
  dev.copy2eps(file = "2017_Paper/Figure1E_raw_v6.eps")
  PlotSiteKds("miR-1", "equilibrium", 5, "paper", "paper", xpos = 20, ypos = 820)
  dev.copy2eps(file = "2017_Paper/Figure1F_raw_v6.eps")
}

MakeFigure2 <- function() {
  graphics.off()
  PlotSiteKds("let-7a", "equilibrium", 5, "paper", "paper", xpos = 20, ypos = 20)
  dev.copy2eps(file = "2017_Paper/Figure2A_raw_v6.eps")
  PlotSiteKds("miR-155", "equilibrium", 5, "paper", "paper", xpos = 820, ypos = 20)
  dev.copy2eps(file = "2017_Paper/Figure2B_raw_v6.eps")
  PlotSiteKds("miR-124", "equilibrium", 5, "paper", "paper", xpos = 20, ypos = 620)
  dev.copy2eps(file = "2017_Paper/Figure2C_raw_v6.eps")
  PlotSiteKds("lsy-6", "equilibrium", 5, "paper", "paper", xpos = 820, ypos = 620)
  dev.copy2eps(file = "2017_Paper/Figure2D_raw_v6.eps")
}

MakeSupplementalFigure2 <- function() {
  graphics.off()
  PlotBaekKds("miR-1", "equilibrium", 5, xpos = 20, ypos = 20)
  dev.copy2eps(file = "2017_Paper/FigureS2A_raw_v6.eps")
  PlotBaekKds("let-7a", "equilibrium", 5, xpos = 820, ypos = 20)
  dev.copy2eps(file = "2017_Paper/FigureS2B_raw_v6.eps")
  PlotBaekKds("miR-155", "equilibrium", 5, xpos = 20, ypos = 220)
  dev.copy2eps(file = "2017_Paper/FigureS2C_raw_v6.eps")
  PlotBaekKds("miR-124", "equilibrium", 5, xpos = 820, ypos = 220)
  dev.copy2eps(file = "2017_Paper/FigureS2D_raw_v6.eps")
  PlotBaekKds("lsy-6", "equilibrium", 5, xpos = 20, ypos = 420)
  dev.copy2eps(file = "2017_Paper/FigureS2E_raw_v6.eps")
  PlotPositionalKds("equilibrium", 5, "centered11", xpos = 820, 420)
  dev.copy2eps(file = "2017_Paper/FigureS2F_raw_v6.eps")
}

MakeFigure3 <- function() {
  graphics.off()
  # PlotSiteFlankEnrichments("miR-1", "equilibrium", 5, "paper", "paper","8mer", xpos = 20, ypos = 20)
  # dev.copy2eps(file = "2017_Paper/Figure3A_raw_v6.eps")
  PlotSiteFlankKds("miR-1", "equilibrium", 5, "paper", "paper", xpos = 820, ypos = 20)
  # dev.copy2eps(file = "2017_Paper/Figure3B_raw_v6.eps")
  PlotCanonicalSiteFlankingKdHeatmap("equilibrium", 5, "paper", xpos = 20, ypos = 420)
  print("hi")
  dev.copy2eps(file = "2017_Paper/Figure3C_raw_v6.eps")

}



# graphics.off()
# # MakeFigure1()
# # MakeFigure2()

# # MakeSupplementalFigure2()

# MakeFigure3()




GetAllFlanks <- function(start, stop, sitelist, log.residual=FALSE, combined.input=TRUE) {
  output <- GetAllMirnaFlankKds("let-7a", "equilibrium", start, stop, sitelist, log.residual=log.residual, combined.input=combined.input)
  colnames(output) <-sapply(colnames(output), function(name) {paste0("let-7a_", name)})
  for (mirna in c("miR-1", "miR-155", "miR-124", "lsy-6")) {
    output.temp <- GetAllMirnaFlankKds(mirna, "equilibrium", start, stop, sitelist, log.residual=log.residual, combined.input=combined.input)
    colnames(output.temp) <-sapply(colnames(output.temp), function(name) {paste0(mirna, "_", name)})

    output <- cbind(output, output.temp)
  }
  output <- log10(output)
  output <- t(t(output) - colMeans(output, na.rm =TRUE))
  return(output)

}

MakeKdandSequenceTable <- function(mirna, experiment, start, stop, sitelist,
                                   log.residual=FALSE, combined.input=TRUE,
                                   nosite=TRUE){
  kds <- GetKds(mirna, experiment, start, stop, sitelist,
                log.residual=log.residual, combined.input=combined.input,
                nosite=nosite)
  print(kds)
  data <- GetSitesXCountsCombined(mirna, experiment, start, stop, sitelist)
  print(data[,1:2])
  out <- data.frame(seq = data[,1], kd = as.numeric(kds[rownames(data)]))
  rownames(out) <- rownames(data)
  return(out)
}

MakeSiteKdBeeswarms <- function(mirna, experiment, start, stop, sitelist,
                                site_list=NULL,flankdata=NULL,
                                colorByPoints=FALSE,
                                sitelist.print=FALSE, log.residual=FALSE,
                                combined.input=TRUE, nosite=TRUE) {

  dev.new(xpos = 20, ypos = 20, height = 8, width = 12)
  if (length(site_list) == 0) {
    site_list <- rownames(data)
  } else if (class(site_list) == "character") {
    site_list <- c(read.table(file = paste0("/lab/bartel1_ata/mcgeary/",
                                            "computation/AgoRBNS/",
                                            "AssignSiteTypes/sites.", mirna,
                                            "_", site_list, ".txt"),
                              stringsAsFactors=FALSE)[,1], "None")
  } else {
    site_list <- c(rownames(data)[order(kds)][1:site_list], "None")
  }
  num.sites <- length(site_list)
  # Get kds for all site-types of the mirna.
  site.kds <- GetKds(mirna, experiment, start, stop, sitelist,
                     log.residual=log.residual, combined.input=combined.input,
                     nosite=nosite)
  print(site.kds)
  site.kds <- site.kds[1 : (length(site.kds) - 2)]  
  if (nosite == FALSE) {
    site.kds <- c(site.kds, 1)
    names(site.kds)[length(site.kds)] <- "None"
  }
  # # Get flanking kds for the mirna.
  if (length(flankdata) == 0){
      flank.kds <- GetAllMirnaFlankKds(mirna, experiment, start, stop, sitelist,
                                   log.residual=log.residual,
                                   combined.input=combined.input,nosite=nosite)
  } else {
    flank.kds <- flankdata
  }
  site.kds <- site.kds[which(colSums(is.na(flank.kds))!=256)]
  site.kds <- site.kds[site_list[-length(site_list)]]
  site.kds <<- site.kds
  site_list <<- site_list
  flank.kds.trim <- flank.kds[,site_list[-length(site_list)]]
  flank.kds.trim <- flank.kds.trim[,order(site.kds)]
  flank.kds.trim <<- flank.kds.trim
  kds.flanks <- c(flank.kds.trim)
  data.sites <- rep(colnames(flank.kds.trim), each=nrow(flank.kds.trim))
  data.ranks <- rep(num.sites-seq(num.sites-1), each=nrow(flank.kds.trim))
  data.colors <- rep(sapply(rownames(flank.kds), GetColorFunction, alpha=0.7),
                     ncol(flank.kds.trim))
  flanks.df <- data.frame(kds=log10(kds.flanks), rank=data.ranks, sites=data.sites,
                          cols=data.colors,
                           stringsAsFactors=FALSE)
  print(unique(flanks.df$rank))
  print(unique(flanks.df$site))

  print(unique(flanks.df$rank))
  print(unique(flanks.df$site))
  flanks.df <<- flanks.df
  par(kPlotParameters)
  boxplot(kds ~ rank,
          data       = flanks.df,
          axes       = FALSE,
          horizontal = TRUE,
          outline    = FALSE,
          xlim       = c(0, num.sites+5),
          ylim       = rev(c(-4,1)))
  print(c(0,num.sites+5))
  title(main = mirna,
        line = -2,
        adj  = 0.1)
  title(xlab = expression(K[D]))
  if (colorByPoints == TRUE){
    beeswarm(kds ~ rank,
             data       = flanks.df,
             add        = TRUE,
             method     = "swarm",
             corral     = "random",
             pch        = 1,
             horizontal = TRUE,  
             pwcol = cols)
  } else {
    beeswarm(kds ~ rank,
         data       = flanks.df,
         add        = TRUE,
         method     = "swarm",
         corral     = "random",
         pch        = 20,
         horizontal = TRUE,  
         col        = kSiteColors[rev(unique(flanks.df$sites)), ])
  }
    # points(log10(site.kds[order(site.kds)]),num.sites - seq(num.sites-1),col=kSiteColors[names(site.kds),][order(site.kds)],cex=3)

  ymin=0.0001
  ymax=10
  ys <- c(sapply(seq(floor(min(flanks.df$kds, na.rm=TRUE)),
                     ceiling(max(flanks.df$kds, na.rm=TRUE))), function(x) seq(10)*10^x))
  ys <- ys[ys >= ymin & ys <= ymax]
  ymin <- min(ys)
  ymax <- max(ys)
  yl <- seq(ceiling(log10(min(ys))), floor(log10(max(ys))))
  axis(1, at=yl,
       labels=sapply(yl, function(name) {
         eval(substitute(expression(10^x), list(x=log10(name))))
       }),
       pos=0, lwd=0, cex.axis=1.5, padj=0.2)
  axis(1, at=log10(ys), labels=FALSE,
       pos=0, lwd=2)
  points(rep(-3.4, num.sites-1),num.sites - seq(num.sites-1),col= kSiteColors[names(site.kds),][order(site.kds)],cex=1.5,pch=19)

  for (ind in seq(num.sites)) {
    site = unique(flanks.df$sites)[ind]
    text(-3.5,num.sites - ind, labels=site, adj=0, cex=1.2, col= "black")
  }
}

MakeSiteBarPlotsOrig <- function(mirna, experiment, start, stop, sitelist, site_list, 
                                colorByPoints=FALSE,
                                sitelist.print=FALSE, log.residual=FALSE,
                                combined.input=TRUE) {
  # Get kds for all site-types of the mirna.
  dev.new(xpos = 20, ypos = 20, height = 8, width = 12)
  site.kds <- GetKdsErrorOrig(mirna, experiment, start, stop, sitelist)
  print(site.kds)
      site_list <- c(read.table(file = paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/AssignSiteTypes/sites.", mirna, "_", site_list,".txt"),
  stringsAsFactors=FALSE)[,1], "None")
  print(site.kds[1:10,])
  print(site_list)
  print(rownames(site.kds))
  print(site.kds)
  site.kds <- site.kds[which(rownames(site.kds) %in% site_list),]
  site.kds <- site.kds / c(site.kds["None",])
  flanks.df <- data.frame(kds=site.kds[,2],upper=site.kds[,2],lower=site.kds[,3], site = rownames(site.kds),stringsAsFactors=FALSE)
  flanks.df <- flanks.df[order(flanks.df$kds),]
  print(flanks.df)
  num.sites <- nrow(site.kds)
  print(dim(flanks.df))
  par(kPlotParameters)
  xs <- site.kds
  ys <- nrow(site.kds) - seq(nrow(site.kds)) + 1
  xs <<- xs
  ys <<- ys
  print(length(flanks.df$kds))
  print(kSiteColors[flanks.df$site,])
  print(length(nrow(site.kds) - seq(nrow(site.kds)) + 1))
  plot(flanks.df$kds,nrow(site.kds) - seq(nrow(site.kds)) + 1,
          col = "white",
          axes    = FALSE,
          log = 'x',
          pch = 20,
          cex = 2,
          ylim       = c(0, num.sites+5),
          xlim       = rev(c(0.00001, 1.1)))
 arrows(flanks.df$upper, nrow(site.kds) - seq(nrow(site.kds)) + 1,
        flanks.df$lower, nrow(site.kds) - seq(nrow(site.kds)) + 1, length=0.07, angle=90, code=3)

  title(main = mirna,
        line = -2,
        adj  = 0.1)
  title(xlab = expression(K[D]))
 points(flanks.df$kds,nrow(site.kds) - seq(nrow(site.kds)) + 1,
          bg = kSiteColors[flanks.df$site,],
          col = "black",
          pch = 21,
          cex = 1.5)
  ymin=0.0001
  ymax=1
  ys <- c(sapply(seq(log10(ymin),
                     log10(ymax)), function(x) seq(10)*10^x))
  # ys <- ys[ys >= ymin & ys <= ymax]
  ymin <- min(ys)
  ymax <- max(ys)
  yl <- 10^seq(ceiling(log10(min(ys))), floor(log10(max(ys))))
  axis(1, at=yl,
       labels=sapply(yl, function(name) {
         eval(substitute(expression(10^x), list(x=log10(name))))
       }),
       pos=0, lwd=0, cex.axis=1.5, padj=0.2)
  axis(1, at=ys, labels=FALSE,
       pos=0, lwd=2)
print(unique(flanks.df$rank))
print(unique(flanks.df$sites))
print(flanks.df[1:10,])
points(rep(0.00006, nrow(site.kds)),nrow(site.kds) - seq(nrow(site.kds)) + 1,col="black",bg=kSiteColors[flanks.df$site,],cex=1.5,pch=21)
sapply(seq(num.sites), function(ind) {
  site = flanks.df$site[ind]

  text(0.00005,num.sites - ind + 1, labels=site, adj=0, cex=1.2, col= "black")
  })
}






PlotFlankLinearModel <- function(flanks, site_fit, site_y,  int2=FALSE, int3=FALSE) {

  flanks <- log10(flanks)
  flanks <- t(t(flanks) - colMeans(flanks, na.rm=TRUE))
  fit <- LinearModelInputFlanks(flanks, which(colnames(flanks) == site_fit))
  if (int3 == TRUE) {
    fit2 <- lm(I ~ f5p.i*f5p.o*f3p.i*f3p.o, data=fit)      
  } else if (int2 == TRUE) {
    fit2 <- lm(I ~ f5p.i*f5p.o + f3p.i*f3p.o, data=fit)      
  } else {
    fit2 <- lm(I ~ f5p.i + f5p.o + f3p.i + f3p.o, data=fit)  
  }
  plot(  fit2$fitted.values,  flanks[,site_y], xlab = site_fit, pch = 19, cex = 2.0, ylab = site_y, col = kNucleotideColors[fit$f5p.o])
  points(fit2$fitted.values,  flanks[,site_y],                  pch = 19, cex = 1.5, col = kNucleotideColors[fit$f5p.i])
  points(fit2$fitted.values,  flanks[,site_y],                  pch = 19, cex = 1.0, col = kNucleotideColors[fit$f3p.i])
  points(fit2$fitted.values,  flanks[,site_y],                  pch = 19, cex = 0.5, col = kNucleotideColors[fit$f3p.o])
  cor_text <- round(
                cor(fit2$fitted.values, flanks[, site_y], use="complete"),
                digits = 3
              )

  title(main = bquote(italic(r) ==  .(cor_text)), font.main=1, cex.main=1.5, line=-3.5, adj=0.1)

}

PlotFlankLinearModelAverage <- function(flanks, site_y, int2=FALSE, int3=FALSE) {

  flanks <- log10(flanks)
  flanks <- t(t(flanks) - colMeans(flanks, na.rm=TRUE))
  flanks <- cbind(rowMeans(flanks, na.rm=TRUE), flanks)
  fit <- LinearModelInputFlanks(flanks, 1)
  if (int3 == TRUE) {
    fit2 <- lm(I ~ f5p.i*f5p.o*f3p.i*f3p.o, data=fit)      
  } else if (int2 == TRUE) {
    fit2 <- lm(I ~ f5p.i*f5p.o + f3p.i*f3p.o, data=fit)      
  } else {
    fit2 <- lm(I ~ f5p.i + f5p.o + f3p.i + f3p.o, data=fit)  
  }
  # print(fit)
  # print(fit2)
  # # print(flanks[, site_y])
  plot(  fit2$fitted.values,  flanks[,site_y], xlab = "average", pch = 19, cex = 2.0, ylab = site_y, col = kNucleotideColors[fit$f5p.o])
  points(fit2$fitted.values,  flanks[,site_y],                  pch = 19, cex = 1.5, col = kNucleotideColors[fit$f5p.i])
  points(fit2$fitted.values,  flanks[,site_y],                  pch = 19, cex = 1.0, col = kNucleotideColors[fit$f3p.i])
  points(fit2$fitted.values,  flanks[,site_y],                  pch = 19, cex = 0.5, col = kNucleotideColors[fit$f3p.o])
  cor_text <- round(
                cor(fit2$fitted.values, flanks[, site_y], use="complete"),
                digits = 3
              )

  title(main = bquote(italic(r) ==  .(cor_text)), font.main=1, cex.main=1.5, line=-3.5, adj=0.1)

}



# CompareWobbles()



MakeWithWithoutProteinScatter <- function(mirna, experiment, start, stop, sitelist) {
    sitesXcounts <- GetSitesXCounts(mirna, experiment, start, stop, sitelist)
    colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
      colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
        return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
      }
    )
      params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/", experiment, "/kds_with_structure/",
                   start, "-", stop, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg_free_protein.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.1 <- Logistic(params[1:nrow(sitesXcounts)], 10)
      k.c.stockago.1 <- 10^params["AGO"]
      bgs.1 <- 10^params["bg"]

      params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/", experiment, "/kds_with_structure/",
                   start, "-", stop, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.2 <- Logistic(params[1:nrow(sitesXcounts)], 10)
      stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
                       "/SolveForKds/k_c_stockago.txt"), row.names=1,
                       header=TRUE, sep="\t")
      k.c.stockago.2 <- stockago[mirna,experiment]

      bgs.2 <- 10^params["bg"]

    print(kds.1)
    print(kds.2)
    print(k.c.stockago.1)
    print(k.c.stockago.2)
    xmin <- 10^floor(min(log10(kds.1), log10(kds.2)))
    xmax <- 10^ceiling(max(log10(as.numeric(unlist(c(kds.1, kds.2, bgs.2, bgs.1, k.c.stockago.2, k.c.stockago.1))))))
    print(xmin)
    print(xmax)
    xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
    xs <- xs[xs >= xmin & xs <= xmax]
    xmin <- min(xs)
    xmax <- max(xs)
    xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))


    plot(c(kds.1, bgs.1, k.c.stockago.1),
      c(kds.2, bgs.2, k.c.stockago.2),
      col = kSiteColors[c(rownames(sitesXcounts), "bg", "Ago"),],
      pch = 19,
      log = 'xy',
      xlim = c(xmin,xmax*(xmax/xmin)^0.55), 
      ylim = c(xmin, xmax),
      axes = FALSE,
      ann = FALSE)
        segments(xmin, xmin, xmax, xmax, lty = 2)

    axis(1, at=xl,
         labels=sapply(xl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=xmin, lwd=0, cex.axis=1.5, padj=0.2)
    axis(1, at=xs, labels=FALSE,
         pos=xmin, lwd=2)
    # Label the axis at each order of magnitude.

    axis(2, at=xl,
         labels=sapply(xl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
    axis(2, at=xs, labels=FALSE,
         pos=xmin, lwd=2)

        title(main = mirna, font.main=1, cex.main=1.5, line=-2, adj=0.1)
    title(xlab = "Fixed protein parameters (nM)", cex.lab=1.5, line=2, adj=0.3)
    title(ylab = "Fit protein parameters (nM)", cex.lab=1.5, line=2)
    legend.names <- c(rownames(sitesXcounts), "bg", "Ago")
    # centered_sites <- c("11mer-m3.13", "12mer-m3.14",
    #                 "11mer-m4.14", "12mer-m4.15")
    # centered_index <- sapply(centered_sites, function(site){
    #                        which(legend.names == site)
    #                      }
    #                      )
    # legend.names[centered_index] <- "Centered"
    # legend.names <- unique(legend.names)
    legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=0.9, bty="n", ncol = 1)

    
}


MakeScatterAcrossWindows <- function(mirna, experiment, lim1, lim2, sitelist, method) {
    sitesXcounts <- GetSitesXCounts(mirna, experiment, lim1, lim1, sitelist)
    colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
      colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
        return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
      }
    )
    if (method == "free_protein") {

      params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/", experiment, "/kds_with_structure/",
                   lim1, "-", lim1, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg_free_protein.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.1 <- Logistic(params[1:nrow(sitesXcounts)], 10)
      k.c.stockago.1 <- 10^params["AGO"]
      bgs.1 <- 10^params["bg"]

      params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/", experiment, "/kds_with_structure/",
                   lim2, "-", lim2, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg_free_protein.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.2 <- Logistic(params[1:nrow(sitesXcounts)], 10)
      k.c.stockago.2 <- 10^params["AGO"]
      bgs.2 <- 10^params["bg"]

      } else {
      params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/", experiment, "/kds_with_structure/",
                   lim1, "-", lim1, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.1 <- Logistic(params[1:nrow(sitesXcounts)], 10)
      stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
                       "/SolveForKds/k_c_stockago.txt"), row.names=1,
                       header=TRUE, sep="\t")
      k.c.stockago.1 <- stockago[mirna,experiment]

      bgs.1 <- 10^params["bg"]


      params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/", experiment, "/kds_with_structure/",
                   lim2, "-", lim2, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.2 <- Logistic(params[1:nrow(sitesXcounts)], 10)
      stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
                       "/SolveForKds/k_c_stockago.txt"), row.names=1,
                       header=TRUE, sep="\t")
      k.c.stockago.2 <- stockago[mirna,experiment]

      bgs.2 <- 10^params["bg"]
}
    print(kds.1)
    print(kds.2)
    xmin <- 10^floor(min(log10(kds.1), log10(kds.2)))
    xmax <- 10^ceiling(max(log10(as.numeric(unlist(c(kds.1, kds.2, bgs.2, bgs.1, k.c.stockago.2, k.c.stockago.1))))))
    print(xmin)
    print(xmax)
    xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
    xs <- xs[xs >= xmin & xs <= xmax]
    xmin <- min(xs)
    xmax <- max(xs)
    xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))


    plot(c(kds.1, bgs.1, k.c.stockago.1),
      c(kds.2, bgs.2, k.c.stockago.2),
      col = kSiteColors[c(rownames(sitesXcounts), "bg", "Ago"),],
      pch = 19,
      log = 'xy',
      xlim = c(xmin,xmax*(xmax/xmin)^0.55), 
      ylim = c(xmin, xmax),
      axes = FALSE,
      ann = FALSE)
    segments(xmin, xmin, xmax, xmax, lty = 2)
    axis(1, at=xl,
         labels=sapply(xl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=xmin, lwd=0, cex.axis=1.5, padj=0.2)
    axis(1, at=xs, labels=FALSE,
         pos=xmin, lwd=2)
    # Label the axis at each order of magnitude.

    axis(2, at=xl,
         labels=sapply(xl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
    axis(2, at=xs, labels=FALSE,
         pos=xmin, lwd=2)

        title(main = mirna, font.main=1, cex.main=1.5, line=-2, adj=0.1)
    title(xlab = paste0(lim1, " parameters (nM)"), cex.lab=1.5, line=2, adj=0.3)
    title(ylab = paste0(lim2, " parameters (nM)"), cex.lab=1.5, line=2)
    legend.names <- c(rownames(sitesXcounts), "bg", "Ago")
    # centered_sites <- c("11mer-m3.13", "12mer-m3.14",
    #                 "11mer-m4.14", "12mer-m4.15")
    # centered_index <- sapply(centered_sites, function(site){
    #                        which(legend.names == site)
    #                      }
    #                      )
    # legend.names[centered_index] <- "Centered"
    # legend.names <- unique(legend.names)
    legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=0.9, bty="n", ncol = 1)

    
}

MakeScatterAcrossOrder <- function(mirna, experiment, lim, sitelist1, sitelist2, method) {
    sitesXcounts <- GetSitesXCounts(mirna, experiment, lim, lim, sitelist1)
    colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
      colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
        return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
      }
    )
    sitesXcounts1 <- sitesXcounts
        sitesXcounts <- GetSitesXCounts(mirna, experiment, lim, lim, sitelist2)
    colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
      colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
        return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
      }
    )
    sitesXcounts2 <- sitesXcounts
    print(sitesXcounts2)
    if (method == "free_protein") {

      params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/", experiment, "/kds_with_structure/",
                   lim, "-", lim, "_", 1, "-", 15, "_", sitelist1, "_Basic_singlebg_free_protein.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.1 <- Logistic(params[1:nrow(sitesXcounts1)], 10)
      k.c.stockago.1 <- 10^params["AGO"]
      bgs.1 <- 10^params["bg"]

      params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/", experiment, "/kds_with_structure/",
                   lim, "-", lim, "_", 1, "-", 15, "_", sitelist2, "_Basic_singlebg_free_protein.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.2 <- Logistic(params[1:nrow(sitesXcounts2)], 10)
      k.c.stockago.2 <- 10^params["AGO"]
      bgs.2 <- 10^params["bg"]

      } else {
      params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/", experiment, "/kds_with_structure/",
                   lim, "-", lim, "_", 1, "-", 15, "_", sitelist1, "_Basic_singlebg.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.1 <- Logistic(params[1:nrow(sitesXcounts1)], 10)
      stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
                       "/SolveForKds/k_c_stockago.txt"), row.names=1,
                       header=TRUE, sep="\t")
      k.c.stockago.1 <- stockago[mirna,experiment]

      bgs.1 <- 10^params["bg"]


      params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/", experiment, "/kds_with_structure/",
                   lim, "-", lim, "_", 1, "-", 15, "_", sitelist2, "_Basic_singlebg.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.2 <- Logistic(params[1:nrow(sitesXcounts2)], 10)
      stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
                       "/SolveForKds/k_c_stockago.txt"), row.names=1,
                       header=TRUE, sep="\t")
      k.c.stockago.2 <- stockago[mirna,experiment]

      bgs.2 <- 10^params["bg"]
}
    names(kds.1) <- rownames(sitesXcounts1)
    names(kds.2) <- rownames(sitesXcounts2)
    kds.2 <- kds.2[names(kds.1)]
    print(kds.1)
    print(kds.2)
    xmin <- 10^floor(min(log10(kds.1), log10(kds.2)))
    xmax <- 10^ceiling(max(log10(as.numeric(unlist(c(kds.1, kds.2, bgs.2, bgs.1, k.c.stockago.2, k.c.stockago.1))))))
    print(xmin)
    print(xmax)
    print(bgs.1)
    print(bgs.2)
    print(k.c.stockago.1)
    print(k.c.stockago.2)
    xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
    xs <- xs[xs >= xmin & xs <= xmax]
    xmin <- min(xs)
    xmax <- max(xs)
    xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))


    plot(c(kds.1, bgs.1, k.c.stockago.1),
      c(kds.2, bgs.2, k.c.stockago.2),
      col = kSiteColors[c(rownames(sitesXcounts1), "bg", "Ago"),],
      pch = 19,
      log = 'xy',
      xlim = c(xmin,xmax*(xmax/xmin)^0.55), 
      ylim = c(xmin, xmax),
      axes = FALSE,
      ann = FALSE)
    segments(xmin, xmin, xmax, xmax, lty = 2)
    axis(1, at=xl,
         labels=sapply(xl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=xmin, lwd=0, cex.axis=1.5, padj=0.2)
    axis(1, at=xs, labels=FALSE,
         pos=xmin, lwd=2)
    # Label the axis at each order of magnitude.

    axis(2, at=xl,
         labels=sapply(xl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
    axis(2, at=xs, labels=FALSE,
         pos=xmin, lwd=2)

        title(main = mirna, font.main=1, cex.main=1.5, line=-2, adj=0.1)
    title(xlab = paste0(sitelist1, " parameters (nM)"), cex.lab=1.5, line=2, adj=0.3)
    title(ylab = paste0(sitelist2, " parameters (nM)"), cex.lab=1.5, line=2)
    legend.names <- c(rownames(sitesXcounts1), "bg", "Ago")
    # centered_sites <- c("11mer-m3.13", "12mer-m3.14",
    #                 "11mer-m4.14", "12mer-m4.15")
    # centered_index <- sapply(centered_sites, function(site){
    #                        which(legend.names == site)
    #                      }
    #                      )
    # legend.names[centered_index] <- "Centered"
    # legend.names <- unique(legend.names)
    legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=0.9, bty="n", ncol = 1)

    
}






MakeSiteIterationPlotKinetics <- function(out,method,mirna, num.kds, num.bgs, colors = FALSE) {
  setEPS()
  postscript(file=paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/figures/koffs/",
    mirna,"/", method, ".eps"))
  par(kPlotParameters)
  x = seq(1,dim(out)[1],length = max(ncol(out),1000))
  out <- out[x,]
  # print(out)
  probs = out[ ,"-logp"]
  out.print <- out[ , seq(dim(out)[2] - 1)]
  out.print.koffs <- Logistic(out[,1:num.kds],max = 10)
  out.print.percentages <- Logistic(out[,(num.kds+1):(2*num.kds)],max = 1)
  out.print.setpoint <- 10^out[,2*num.kds + 1]
  out.print.bgs <- 10^out.print[,(2 * num.kds + 2) : (2 * num.kds + num.bgs + 1)]

  out.print <- cbind(out.print.koffs,out.print.percentages,out.print.setpoint, out.print.bgs)
  ys <- 10^c(max(floor(log10(min(out.print))), -6), ceiling(log10(max(out.print))))
  probs.scaled <- probs / max(probs) * (log10(ys[2]) - log10(ys[1])) + log10(ys[1])


  plot(x    = x,
       y    = 10^probs.scaled,
       log  = "y",
       axes = FALSE,
       ann = FALSE, 
       type = "l",
       col = "white",
       ylim = ys)
  axis(side = 1,
       at   = seq(1, max(x), by=max(1, floor(max(x) / 20))),
       pos  = ys[1],
       tck  = -0.01,
       labels = FALSE)
  axis(1, at=seq(1, max(x), by=max(1, floor(max(x) / 20))*5), pos=ys[1], lwd=2)
  axis(2, at=10^seq(log10(ys[1]), log10(ys[2])),
       labels=sapply(seq(log10(ys[1]), log10(ys[2])), function(name) {
         eval(substitute(expression(10^x), list(x=name)))
       }),
       pos=1, lwd=2)
  title(main = mirna, font=1)
  axis(2, at=c(sapply(seq(log10(ys[1]), log10(ys[2])-1), function(x) seq(9)*10^x),ys[2]), labels=FALSE,
       pos=1)
  if (colors != FALSE) {

  cols <- c(colors,kSiteColors[colnames(out.print)[(length(colors)+1):length(colnames(out.print))],])
  names(cols) <- colnames(out.print)
  cols["AGO"] <- "grey"

  sapply(colnames(out.print),function(name) {
    lines(x, out.print[, name], lwd=1, col=cols[name])
    })
  } else  {
    sapply(colnames(out.print),function(name) {
    lines(x, out.print[, name], lwd=1, col=kSiteColors[name,])
    })
  }
  lines(x, 10^probs.scaled, type="l", lwd=2, col="black")
  dev.off()
}

centered_sites <- c("11mer-m3.13", "12mer-m3.14",
                    "11mer-m4.14", "12mer-m4.15")


kSiteColors <- read.table(
  "/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/site_info.txt",
  row.names=1, header=FALSE, sep="\t", stringsAsFactors=FALSE)[, 1, drop = FALSE]





Norm <- function(vector) {
  return(vector/ sum(vector))
}

Cumul <- function(vector) {
  norm <- Norm(vector)
  tot <- 0
  out <- sapply(norm, function(x){
    tot <<- tot + x
    return(tot)

    })
  return(out)
}

Logistic <- function(vector, max) {
  return(max / (1 + exp(-vector)))
}

Logit <- function(vector, max) {
  return(-log(max / vector - 1))
}


MakeIterationPlot <- function(out,type,extension = "") {
  setEPS()
  postscript(file=paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
                         "/figures/kds/", mirna, "/", type,
                         "/iterations/", site, "_", k.c.stockago, extension,".eps"))
  par(kPlotParameters)
  x = seq(dim(out)[1])
  probs = out[ ,"-logp"]
  out.print <- out[ , seq(dim(out)[2] - 1)]
  ys <- 10^c(floor(min(out.print)), ceiling(max(out.print)))
  probs.scaled <- probs / max(probs) * (log10(ys[2]) - log10(ys[1])) + log10(ys[1])
  out.print <- 10^out.print
  plot(x , 10^probs.scaled, log='y', axes=FALSE, type="l", ylim=ys,
       lwd=2, ann=FALSE,
       col="black")
  title(main = mirna, line=-1, adj=0.1)
  title(main = site, col.main=kSiteColors[site,], line=-2.5, adj=0.1)

  title(xlab = "Iteration")
  title(ylab = "Parameter values (nM)")
  axis(1, at=seq(1, max(x), by=max(1, floor(max(x) / 20))), pos=ys[1], lwd=2,
       labels=FALSE, tck=-0.01)
  axis(1, at=seq(1, max(x), by=max(1, floor(max(x) / 20))*5), pos=ys[1], lwd=2)
  axis(2, at=10^seq(log10(ys[1]), log10(ys[2])),
       labels=sapply(seq(log10(ys[1]), log10(ys[2])), function(name) {
         eval(substitute(expression(10^x), list(x=name)))
       }),
       pos=1, lwd=2, hadj=0.8)
  axis(2, at=c(sapply(seq(log10(ys[1]), log10(ys[2])-1), function(x) seq(9)*10^x),ys[2]), labels=FALSE,
       pos=1, lwd=2)
  sapply(colnames(out.print), function(name) {
          lines(x, out.print[, name], lwd=2, col=GetColorFunction(name))
        }
        )
  lines(x, 10^probs.scaled, type="l", col="black")
  dev.off()
}



WriteIterationFile <- function(out,extension="") {
  out.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                     "/equilibrium/kds/", site,"_flanking_",
                     k.c.stockago, extension,".txt")
  write.table(file=out.file, out, sep="\t", quote=FALSE, row.names=FALSE,
                col.names=TRUE)
}

WriteFinalParameterFile <- function(out,extension="") {
  out.final <- out[dim(out)[1], ]
  names(out.final) <- colnames(out)
  out.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                           "/equilibrium/kds/final_", site,
                           "_flanking_", k.c.stockago,extension, ".txt")
  write.table(file=out.file, out.final, sep="\t", quote=FALSE, row.names=TRUE,
              col.names=FALSE)
}












# FIXED FOR PAPER
PlotSiteKdOptimization <- function(out, specific_string, mirna, num.kds,
                                   num.bgs, colors = FALSE) {
  # Assign file name to the figure.
  setEPS()
  filename <- paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/",
                         "2017_Paper/Kd_fits/", mirna, "_", specific_string, ".eps")
  postscript(file=filename)
  # Load standard plot settings.
  par(kPlotParameters)
  # Define x as either the list of all rows, or alternatively 1000 equally
  # spaced rows.
  x <- seq(1, nrow(out), length=min(nrow(out), 1000))
  out <- out[x, ]
  # Split output matrix into loglikelihood and parameters:
  probs <- out[, "-logp"]
  pars <- out[, seq(ncol(out) - 1)]
  kds <- Logistic(pars[, 1:num.kds], max = 10)
  bgs <- 10^pars[,(num.kds + 1) : (num.kds + num.bgs)]
  AGO <- 10^pars[,(num.kds + num.bgs + 1) : ncol(pars)]
  out.print <- cbind(kds, bgs, AGO)
  ys <- c(10^-10, 10^5)
  probs.scaled <- probs / probs[1]

  plot(x    = x,
       y    = 10^probs.scaled,
       log  = "y",
       axes = FALSE,
       type = "l",
       col = "white",
       ylim = ys)
  axis(side = 1,
       at   = seq(1, max(x), by=max(1, floor(max(x) / 20))),
       pos  = ys[1], labels = FALSE,
       tck  = -0.01)
  axis(1, at=seq(1, max(x), by=max(1, floor(max(x) / 20))*5), pos=ys[1], lwd=2)
  axis(2, at=10^seq(log10(ys[1]), log10(ys[2])),
       labels=sapply(seq(log10(ys[1]), log10(ys[2])), function(name) {
         eval(substitute(expression(10^x), list(x=name)))
       }),
       pos=1, lwd=2)
  title(main = mirna, font=1)
  axis(2, at=c(sapply(seq(log10(ys[1]), log10(ys[2])-1), function(x) seq(9)*10^x),ys[2]), labels=FALSE,
       pos=1)
  if (colors != FALSE) {

  cols <- colors
  names(cols) <- colnames(out.print)
  cols["AGO"] <- "grey"

  sapply(colnames(out.print),function(name) {
    lines(x, out.print[, name], lwd=1, col=cols[name])
    })
  } else  {
    sapply(colnames(out.print),function(name) {
    lines(x, out.print[, name], lwd=1, col=kSiteColors[name,])
    })
  }
  lines(x, ys[1]*(ys[2]/ys[1])^probs.scaled, type="l", lwd=2, col="black")
  dev.off()
}


PlotSiteFlankKdOptimization <- function(out, specific_string, mirna, num.kds,
                                   num.bgs, colors) {
  # Assign file name to the figure.
  setEPS()
  postscript(file=paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/",
                         "2017_Paper/Kd_fits/", mirna, "_", specific_string, ".eps"))
  # Load standard plot settings.
  par(kPlotParameters)
  # Define x as either the list of all rows, or alternatively 1000 equally
  # spaced rows.
  x <- seq(1, nrow(out), length=min(nrow(out), 1000))
  out <- out[x, ]

  # Split output matrix into loglikelihood and parameters:
  probs <- out[, "-logp"]
  pars <- out[, seq(ncol(out) - 1)]
  out.print <- Logistic(pars, max = 10)

  ys <- 10^c(max(floor(log10(min(out.print))), -10), ceiling(log10(max(out.print))))
  probs.scaled <- probs / max(probs) * (log10(ys[2]) - log10(ys[1])) + log10(ys[1])

  plot(x    = x,
       y    = 10^probs.scaled,
       log  = "y",
       axes = FALSE,
       type = "l",
       col = "white",
       ylim = ys)
  axis(side = 1,
       at   = seq(1, max(x), by=max(1, floor(max(x) / 20))),
       pos  = ys[1], labels = FALSE,
       tck  = -0.01)
  axis(1, at=seq(1, max(x), by=max(1, floor(max(x) / 20))*5), pos=ys[1], lwd=2)
  axis(2, at=10^seq(log10(ys[1]), log10(ys[2])),
       labels=sapply(seq(log10(ys[1]), log10(ys[2])), function(name) {
         eval(substitute(expression(10^x), list(x=name)))
       }),
       pos=1, lwd=2)
  title(main = mirna, font=1)
  axis(2, at=c(sapply(seq(log10(ys[1]), log10(ys[2])-1), function(x) seq(9)*10^x),ys[2]), labels=FALSE,
       pos=1)

  cols <- colors
  print(dim(out.print))
  print(length(cols))
  names(cols) <- colnames(out.print)

  sapply(colnames(out.print),function(name) {
    lines(x, out.print[, name], lwd=1, col=cols[name])
    })
  lines(x, 10^probs.scaled, type="l", lwd=2, col="black")
  dev.off()
}







MakeEquilibriumEnrichmentPlots <- function(mirna, experiment, start, stop, sitelist, site_list=NULL, combined.input=TRUE, log.residual=FALSE, bgoff=FALSE, model=TRUE, connected.points=FALSE,nosite=FALSE) {
  dev.new(xpos = 20, ypos = 20, height = 6.9, width = 8.157)
  params <- GetKds(mirna, experiment, start, stop, sitelist, combined.input=combined.input, log.residual=log.residual, scaled=FALSE,nosite=nosite)
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
  if (nosite == FALSE){
    kds <- c(unlist(Logistic(params[1:(nrow(sitesXcounts)-1)], 10)),1)
    } else{
      kds <- Logistic(params[1:nrow(sitesXcounts)], 10)
    }
  names(kds) <- rownames(sitesXcounts)

  bgs <- rep(10^params["bg"], 5)
  k.c.stockago <- 10^params["AGO"]

  c.I.tots <- Norm(sitesXcounts[,1])*100
  names(c.I.tots) <- rownames(sitesXcounts)
  data <- sitesXcounts[,2:6]
  names(c.I.tots) <- rownames(data)
  x_model <- seq(min(as.numeric(colnames(data))) * 0.7, max(as.numeric(colnames(data))) / 0.7,
                 length=100)

  c.totals <- matrix(c.I.tots, nrow=nrow(data), ncol = 100, byrow=FALSE)
  colnames(c.totals) <- colnames(x_model)
  rownames(c.totals) <- rownames(data)

  c.agos <- sapply(x_model, function(x) {
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
  x <- c(40,12.65,4,1.265,0.4)*k.c.stockago/100*1000/2
  print(x)
  y <- c(1,1,1,1,1)
  print(x)
  sites.norm <- Norm(c.I.tots)
  data.norm <- t(t(data)/colSums(data))

  data.R <- data.norm/(sites.norm)
  model.R <- c.final/(sites.norm)

  xmin <- min(x)*0.3
  xmax <- max(x)*3
  ymin <- 0.2
  ymax <- 300
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

  plot(x, y,log='xy', xlim=c(xmin, xmax*(xmax/xmin)^0.5), ylim=c(ymin, ymax), type="l",
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

  if (length(site_list) == 0) {
    site_list <- rownames(data)
  } else if (class(site_list) == "character") {
    site_list <- c(read.table(file = paste0("/lab/bartel1_ata/mcgeary/",
                                            "computation/AgoRBNS/",
                                            "AssignSiteTypes/sites.", mirna,
                                            "_", site_list, ".txt"),
                              stringsAsFactors=FALSE)[,1], "None")
  } else {
    site_list <- c(rownames(data)[order(kds)][1:site_list], "None")
  }

  legend.names <- rownames(data)[order(kds)]
  print(order(kds))
  legend.names <- legend.names[which(legend.names %in% site_list)]
  ordered_list <- legend.names


  centered_sites <- c("11mer-m3.13", "12mer-m3.14", "11mer-m4.14",
                      "12mer-m4.15", "13mer-m3.15", "11mer-m5.15")
  centered_rename <- c("Centered-m3.13", "Centered-m3.14", "Centered-m4.14",
                       "Centered.m4-15", "Centered-m3.15", "Centered-m5.15")
  names(centered_rename) <- centered_sites

  legend.names <- sapply(legend.names, function(site) {
                           if (site %in% centered_sites) {
                             return(centered_rename[site])
                           } else {
                             return(site)
                          }
                         }
                         )

  legend(x=xmax, y=ymax, legend=legend.names, pch=19,
         col=kSiteColors[ordered_list, ], cex=1, bty="n", ncol=1)
  for (name in site_list) {
    if (connected.points == TRUE) {
      type = "o"
    } else {
      type = "p"
    }
    points(x, data.R[name, ], col=kSiteColors[name, ], type=type, pch=19,
           cex=1.5, lwd=3)
    if (model == TRUE) {
    lines(x_model*k.c.stockago/100*1000/2, model.R[name, ], col=kSiteColors[name, ], lwd=2)      
    }
  }
}

MakeOccupancyPlots <- function(mirna, experiment, start, stop, sitelist, site_list=NULL, combined.input=TRUE, log.residual=FALSE, bgoff=FALSE, model=TRUE, connected.points=FALSE,nosite=FALSE) {
  params <- GetKds(mirna, experiment, start, stop, sitelist, combined.input=combined.input, log.residual=log.residual, scaled=FALSE, nosite=nosite)
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

  if (nosite == FALSE){
    kds <- c(unlist(Logistic(params[1:(nrow(sitesXcounts)-1)], 10)),1)
    } else{
      kds <- Logistic(params[1:nrow(sitesXcounts)], 10)
    }
  names(kds) <- rownames(sitesXcounts)

  bgs <- rep(10^params["bg"], 5)
  k.c.stockago <- 10^params["AGO"]
  print(k.c.stockago)
  c.I.tots <- Norm(sitesXcounts[,1])*100
  names(c.I.tots) <- rownames(sitesXcounts)
  data <- sitesXcounts[,2:6]
  names(c.I.tots) <- rownames(data)
  x <- c(40,12.65,4,1.265,0.4)*k.c.stockago/100
  print(x)
  x_model <- seq(min(as.numeric(colnames(data))) * 0.7, max(as.numeric(colnames(data))) / 0.7,
                 length=100)

  c.totals <- matrix(c.I.tots, nrow=nrow(data), ncol = 100, byrow=FALSE)
  colnames(c.totals) <- colnames(x_model)
  rownames(c.totals) <- rownames(data)

  c.agos <- sapply(x_model, function(x) {
      as.numeric(x) * k.c.stockago / 100
    }
  )
  c.bounds <- as.matrix(
    sapply(c.agos, function(x) {
        return(GetBoundRNA(kds, c.I.tots, x))
      }
    )
  )
  c.bounds.points <- as.matrix(sapply(x,function(i) {
    return(GetBoundRNA(kds,c.I.tots,i))
    }
    )
  )
  c.frees <- c.totals - c.bounds


  c.bgs <- t(t(c.frees) * bgs / colSums(c.frees))

  if (bgoff == TRUE) {
    c.all <- c.bounds
  } else {
    c.all <- c.bounds + c.bgs
  }

  c.final <- data.frame(t(t(c.all) / colSums(c.all)))
  rownames(c.final) <- rownames(data)
  print(x)
  y <- c(1,1,1,1,1)
  print(x)
  sites.norm <- Norm(c.I.tots)
  data.norm <- t(t(data)/colSums(data))

  data.R <- data.norm/(sites.norm)
  model.ago.occupancy <- t(t(c.bounds) / colSums(c.bounds))
  model.ago.occupancy.points <- t(t(c.bounds.points) / colSums(c.bounds.points))
  xmin <- min(x)*0.3
  xmax <- max(x)*3
  xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
  xs <- xs[xs >= xmin & xs <= xmax]
  xmin <- min(xs)
  xmax <- max(xs)
  xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
  ymin <- 0
  ymax <- 0.7
  ys <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7)
  yl <- ys

  plot(x, y,log='x', xlim=c(xmin, xmax*(xmax/xmin)^0.55), ylim=c(ymin, ymax), type="l",
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
       labels=yl,
       pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
  axis(2, at=ys, labels=FALSE,
       pos=xmin, lwd=2)

  title(main = mirna, font.main=1, cex.main=1.5, line=-2, adj=0.1)
  title(xlab = "[AGO2-miRNA] (pM)", cex.lab=1.5, line=2, adj=0.3)
  title(ylab = "Enrichment", cex.lab=1.5, line=2)
  print(site_list)
  if (length(site_list) == 0) {
    site_list <- rownames(data)
  } else if (class(site_list) == "character") {
    site_list <- c(read.table(file = paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/AssignSiteTypes/sites.", mirna, "_", site_list,".txt"),
  stringsAsFactors=FALSE)[,1], "None")
  } else {
    site_list <- c(rownames(data)[order(kds)][1:site_list], "None")
  }
  print(site_list)
  legend.names <- rownames(data)[order(kds)]

  legend.names <- legend.names[which(legend.names %in% site_list)]
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
  legend(x=xmax,y=ymax,legend=legend.names, pch=19, col=kSiteColors[ordered_list, ], cex=1.1, bty="n", ncol = 1)
  for (name in site_list) {
    if (connected.points == TRUE) {
      type = "o"
    } else {
      type = "p"
    }
    if (model == TRUE) {
    lines(x_model*k.c.stockago/100*1000/2, model.ago.occupancy[name, ], col=kSiteColors[name, ], lwd=2)      
    points(x,model.ago.occupancy.points[name,],col=kSiteColors[name, ], pch=19,cex=1.5)
    }
  }
}








MakeFlankingEquilibriumEnrichmentPlots <- function(mirna, experiment, start,
                                                   stop, site, sitelist, 
                                                   site_list=NULL,
                                                   combined.input=TRUE,
                                                   log.residual=FALSE,
                                                   bgoff=FALSE, model=TRUE,
                                                   connected.points=FALSE,
                                                   nosite=TRUE) {
  params <- GetKds(mirna, experiment, start, stop, sitelist,
                   combined.input=combined.input, log.residual=log.residual,
                   scaled=FALSE, nosite=nosite)
  params.flanks <- GetFlankKds(mirna, experiment, site, start, stop, sitelist,
                               combined.input=combined.input,
                               log.residual=log.residual,scaled=FALSE,nosite=nosite)

  if (combined.input == TRUE) {
    sitesXcounts <- GetSitesXCountsCombined(mirna, experiment, start, stop,
                                            sitelist)
  } else {
    sitesXcounts <- GetSitesXCounts(mirna, experiment, start, stop, sitelist)
  }
  sitesXcounts <- sitesXcounts[,-1]
  print(params)
  print(params.flanks)
  # Remove any rows with no reads whoatsoever.
  sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
  # Remove any rows for which there are no input reads:
  sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]
  if (nosite==TRUE){
  kds.s <- unlist(Logistic(params[1:nrow(sitesXcounts)], 10))
  } else {
  kds.s <- c(unlist(Logistic(params[1:(nrow(sitesXcounts)-1)], 10)), 1)

  }
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
  print(kds)
  bgs <- rep(10^params["bg"], 5)
  k.c.stockago <- 10^params["AGO"]
  print(bgs)
  print(k.c.stockago)
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
  x_model <- seq(min(as.numeric(colnames(data))) * 0.7, max(as.numeric(colnames(data))) / 0.7,
                 length=100)

  c.totals <- matrix(c.I.tots, nrow=nrow(data), ncol = 100, byrow=FALSE)
  colnames(c.totals) <- colnames(x_model)
  rownames(c.totals) <- rownames(data)

  c.agos <- sapply(x_model, function(x) {
      as.numeric(x) * k.c.stockago / 100
    }
  )

  c.bounds <- as.matrix(
    sapply(c.agos, function(x) {
        return(GetBoundRNA(kds, c.I.tots, x))
      }
    )
  )
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
  data.R <<- data.R
  model.R <<- model.R
  xmin <- min(x)*0.3
  xmax <- max(x)*3
  ymin <- 0.2
  ymax <- 300
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
  stringsAsFactors=FALSE)[,1], "None")
  } else {
    site_list_real <- c(rownames(sitesXcounts.sites)[order(kds.s)][1:site_list], "None")
  }
  site_list_real <- site_list_real[site_list_real!=site]
  print(site_list_real)
  legend.names <- rownames(sfXc)

  ordered_list <- legend.names
  print(site_list_real)
  legend.names <- unique(legend.names)
  print(legend.names)
  legend(x=xmax,y=ymax,legend=legend.names, text.font=2, text.col= sapply(legend.names,GetColorFunction), cex=0.9, bty="n", ncol = 8)
  for (name in c(site_list_real,rownames(sfXc))) {
    print(name)
    if (connected.points == TRUE) {
      type = "o"
    } else {
      type = "p"
    }
    points(x, data.R[name, ], col=colors_all[name], type = type, pch=19, cex=1.5, lwd=3)
    if (model == TRUE) {
    lines(x_model*k.c.stockago/100*1000, model.R[name, ], col=colors_all[name], lwd=2)      
    }
  }
}

PlotAllCenteredKds <- function(start, stop, sitelist, combined.input=TRUE,
                               log.residual=FALSE, nosite=TRUE) {
  # This pulls the last two digits out of the sitelist
  # ("13" from "centered13", etc.)
  split.sitelist <- unlist(strsplit(sitelist, split=""))
  split.sitelist.last <- split.sitelist[(nchar(sitelist)-1):nchar(sitelist)]
  split.sitelist.last <<- split.sitelist.last
  length_sitelist <- paste(split.sitelist.last,collapse="")

  kds_centered <- matrix(unlist(sapply(mirnas_all, function(mirna) {
    kds <- unlist(GetKds(mirna,"equilibrium", start, stop, sitelist,
                  combined.input=combined.input, log.residual=log.residual,
                  nosite=nosite))
    print(kds)
    kds <- c(kds["8mer"],kds["6mer"], kds["6mer-m8"],kds[grep("\\.",names(kds)),drop=FALSE],kds["None"])
    kds <- kds[grep("23", names(kds), invert=TRUE),drop=FALSE]
    return(kds)
  })),nrow=length(mirnas_all),byrow=TRUE)
  rownames(kds_centered) <- mirnas_all
  temp_kds <- GetKds("let-7a", "equilibrium", start, stop, sitelist,
                    combined.input=combined.input, log.residual=log.residual,
                  nosite=nosite)
  temp_kd_names <- names(temp_kds)[grep("\\.", names(temp_kds))]
  colnames(kds_centered) <- c("8mer","6mer", "6mer-m8", temp_kd_names, "None")
  print(kds_centered)
  ylims = c(0.5*min(kds_centered), 2*max(kds_centered))
  plot(1:ncol(kds_centered), kds_centered[1,],
       type="o",log='y',
       ylim=ylims,
       axes=FALSE, 
       ann=FALSE,
       pch=19,lwd=2)
  lines(1:ncol(kds_centered),pch=19,lwd=2,kds_centered[2,],type="o",col="blue")
  lines(1:ncol(kds_centered),pch=19,lwd=2,kds_centered[3,],type="o",col="red")
  lines(1:ncol(kds_centered),pch=19,lwd=2,kds_centered[4,],type="o",col="forestgreen")
  lines(1:ncol(kds_centered),pch=19,lwd=2,kds_centered[5,],type="o",col="purple")
  print(colnames(kds_centered))
  axis(1,pos=ylims[1],at=c(1:ncol(kds_centered)),labels=colnames(kds_centered),las=2)
  axis(2,pos=0.5,at=10^seq(ceiling(log10(min(ylims))-1), floor(log10(max(ylims)))+1))
  legend("topleft",legend = mirnas_all, col=c("black", "blue", "red", "forestgreen", "purple"), lwd=2, pch=19)

}



MakeInputScatter <- function(mirna, experiment, n_constant, ombined.input=TRUE,
                             log.residual=FALSE, bgoff=FALSE,
                             vertical.lines=FALSE) {
  params <- GetSiteKds(mirna, experiment, start, stop, sitelist,
                   combined.input=combined.input, log.residual=log.residual,
                   scaled=FALSE)
  if (combined.input == TRUE) {
    sitesXcounts <- GetSitesXCountsCombined(mirna, experiment, start, stop, sitelist)
  } else {
    sitesXcounts <- GetSitesXCounts(mirna, experiment, start, stop, sitelist)
  }
  print(sitesXcounts)
  sitesXcounts <- sitesXcounts[,-1]
  # Remove any rows with no reads whoatsoever.
  sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
  # Remove any rows for which there are no input reads:
  sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]

  kds <- c(unlist(Logistic(params[1:(nrow(sitesXcounts)-1)], 10)),1)
  names(kds) <- rownames(sitesXcounts)

  if (length(site_list) == 0) {
    site_list <- rownames(sitesXcounts)
  } else if (class(site_list) == "character") {
    site_list <- c(read.table(file = paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/AssignSiteTypes/sites.", mirna, "_", site_list,".txt"),
  stringsAsFactors=FALSE)[,1], "None")
  } else {
    site_list <- c(rownames(data)[order(kds)][1:site_list], "None")
  }
  print(site_list)
  x <- Norm(sitesXcounts[site_list, 1])
  y <- Norm(sitesXcounts[site_list, column.plot])
  xmin <- 0.5 * min(x, y)
  xmax <- 2 * max(x, y)
  ymin <- 0.5 * min(x, y)
  ymax <- 2 * max(x, y)
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

  plot(x, y,log='xy', xlim=c(xmin, xmax*(xmax/xmin)^0.3), ylim=c(ymin, ymax),pch=19, cex = 1.5, ann=FALSE, axes=FALSE,
     col=kSiteColors[site_list, ])        
  # Generate tickmarks for axis.
  axis(1, at=xl,
       labels=xl*100,
       pos=ymin, lwd=0, cex.axis=1.5, padj=0.2)
  axis(1, at=xs, labels=FALSE,
       pos=ymin, lwd=2)
  # Label the axis at each order of magnitude.

  axis(2, at=yl,
       labels=yl*100,
       pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
  axis(2, at=ys, labels=FALSE,
       pos=xmin, lwd=2)
  segments(xmin, xmin, xmax, xmax, lwd = 2, lty = 2)
  if (vertical.lines == TRUE) {
    sapply(1:length(x), function(index) {
      segments(x[index], x[index], x[index], y[index], col=kSiteColors[site_list,][index])
      })

  }
  ago.percent <- as.numeric(colnames(sitesXcounts)[column.plot])/100
  title(main = paste0(ago.percent*stockago[mirna,"equilibrium"]*100,' pM AGO2-', mirna), font.main=1, cex.main=1.5, line=-2.5, adj=0.1)

  title(xlab = "Input library (%)", cex.lab=1.5, line=2, adj=0.3)
  title(ylab = "AGO-bound library (%)", cex.lab=1.5, line=2)
  print(site_list)
  print(site_list)
  legend.names <- rownames(sitesXcounts)[order(kds)]

  legend.names <- legend.names[which(legend.names %in% site_list)]
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
  legend(x=2 * xmax,y=ymax,legend=legend.names, pch=19, col=kSiteColors[ordered_list, ], cex=1.2, bty="n", ncol = 1)
}


MakeWithWithoutProteinScatter <- function(mirna, experiment, start, stop, sitelist) {
    sitesXcounts <- GetSitesXCounts(mirna, experiment, start, stop, sitelist)
    colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
      colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
        return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
      }
    )
      params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/", experiment, "/kds_with_structure/",
                   start, "-", stop, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg_free_protein.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.1 <- Logistic(params[1:nrow(sitesXcounts)], 10)
      k.c.stockago.1 <- 10^params["AGO"]
      bgs.1 <- 10^params["bg"]

      params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/", experiment, "/kds_with_structure/",
                   start, "-", stop, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.2 <- Logistic(params[1:nrow(sitesXcounts)], 10)
      stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
                       "/SolveForKds/k_c_stockago.txt"), row.names=1,
                       header=TRUE, sep="\t")
      k.c.stockago.2 <- stockago[mirna,experiment]

      bgs.2 <- 10^params["bg"]

    xmin <- 10^floor(min(log10(kds.1), log10(kds.2)))
    xmax <- 10^ceiling(max(log10(as.numeric(unlist(c(kds.1, kds.2, bgs.2, bgs.1, k.c.stockago.2, k.c.stockago.1))))))
    print(xmin)
    print(xmax)
    xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
    xs <- xs[xs >= xmin & xs <= xmax]
    xmin <- min(xs)
    xmax <- max(xs)
    xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))


    plot(c(kds.1, bgs.1, k.c.stockago.1),
      c(kds.2, bgs.2, k.c.stockago.2),
      col = kSiteColors[c(rownames(sitesXcounts), "bg", "Ago"),],
      pch = 19,
      log = 'xy',
      xlim = c(xmin,xmax*(xmax/xmin)^0.55), 
      ylim = c(xmin, xmax),
      axes = FALSE,
      ann = FALSE)
        segments(xmin, xmin, xmax, xmax, lty = 2)

    axis(1, at=xl,
         labels=sapply(xl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=xmin, lwd=0, cex.axis=1.5, padj=0.2)
    axis(1, at=xs, labels=FALSE,
         pos=xmin, lwd=2)
    # Label the axis at each order of magnitude.

    axis(2, at=xl,
         labels=sapply(xl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
    axis(2, at=xs, labels=FALSE,
         pos=xmin, lwd=2)

        title(main = mirna, font.main=1, cex.main=1.5, line=-2, adj=0.1)
    title(xlab = "Fixed protein parameters (nM)", cex.lab=1.5, line=2, adj=0.3)
    title(ylab = "Fit protein parameters (nM)", cex.lab=1.5, line=2)
    legend.names <- c(rownames(sitesXcounts), "bg", "Ago")
    # centered_sites <- c("11mer-m3.13", "12mer-m3.14",
    #                 "11mer-m4.14", "12mer-m4.15")
    # centered_index <- sapply(centered_sites, function(site){
    #                        which(legend.names == site)
    #                      }
    #                      )
    # legend.names[centered_index] <- "Centered"
    # legend.names <- unique(legend.names)
    legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=1.2, bty="n", ncol = 1)

    
}

MakeSMNBScatter <- function(start, stop, sitelist) {
    sitesXcounts <- GetSitesXCounts("let-7", "equilibrium", start, stop, sitelist)
    colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
      colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
        return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
      }
    )
      params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/", experiment, "/kds_with_structure/",
                   start, "-", stop, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg_free_protein.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.1 <- Logistic(params[1:nrow(sitesXcounts)], 10)
      k.c.stockago.1 <- 10^params["AGO"]
      bgs.1 <- 10^params["bg"]

      params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/", experiment, "/kds_with_structure/",
                   start, "-", stop, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.2 <- Logistic(params[1:nrow(sitesXcounts)], 10)
      stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
                       "/SolveForKds/k_c_stockago.txt"), row.names=1,
                       header=TRUE, sep="\t")
      k.c.stockago.2 <- stockago[mirna,experiment]

      bgs.2 <- 10^params["bg"]

    xmin <- 10^floor(min(log10(kds.1), log10(kds.2)))
    xmax <- 10^ceiling(max(log10(as.numeric(unlist(c(kds.1, kds.2, bgs.2, bgs.1, k.c.stockago.2, k.c.stockago.1))))))
    print(xmin)
    print(xmax)
    xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
    xs <- xs[xs >= xmin & xs <= xmax]
    xmin <- min(xs)
    xmax <- max(xs)
    xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))


    plot(c(kds.1, bgs.1, k.c.stockago.1),
      c(kds.2, bgs.2, k.c.stockago.2),
      col = kSiteColors[c(rownames(sitesXcounts), "bg", "Ago"),],
      pch = 19,
      log = 'xy',
      xlim = c(xmin,xmax*(xmax/xmin)^0.55), 
      ylim = c(xmin, xmax),
      axes = FALSE,
      ann = FALSE)
        segments(xmin, xmin, xmax, xmax, lty = 2)

    axis(1, at=xl,
         labels=sapply(xl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=xmin, lwd=0, cex.axis=1.5, padj=0.2)
    axis(1, at=xs, labels=FALSE,
         pos=xmin, lwd=2)
    # Label the axis at each order of magnitude.

    axis(2, at=xl,
         labels=sapply(xl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
    axis(2, at=xs, labels=FALSE,
         pos=xmin, lwd=2)

        title(main = mirna, font.main=1, cex.main=1.5, line=-2, adj=0.1)
    title(xlab = "Fixed protein parameters (nM)", cex.lab=1.5, line=2, adj=0.3)
    title(ylab = "Fit protein parameters (nM)", cex.lab=1.5, line=2)
    legend.names <- c(rownames(sitesXcounts), "bg", "Ago")
    # centered_sites <- c("11mer-m3.13", "12mer-m3.14",
    #                 "11mer-m4.14", "12mer-m4.15")
    # centered_index <- sapply(centered_sites, function(site){
    #                        which(legend.names == site)
    #                      }
    #                      )
    # legend.names[centered_index] <- "Centered"
    # legend.names <- unique(legend.names)
    legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=1.2, bty="n", ncol = 1)

    
}

MakeScatterAcrossWindows <- function(mirna, experiment, lim1, lim2, sitelist, method) {
    sitesXcounts <- GetSitesXCounts(mirna, experiment, lim1, lim1, sitelist)
    colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
      colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
        return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
      }
    )
    if (method == "free_protein") {

      params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/", experiment, "/kds_with_structure/",
                   lim1, "-", lim1, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg_free_protein.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.1 <- Logistic(params[1:nrow(sitesXcounts)], 10)
      k.c.stockago.1 <- 10^params["AGO"]
      bgs.1 <- 10^params["bg"]

      params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/", experiment, "/kds_with_structure/",
                   lim2, "-", lim2, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg_free_protein.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.2 <- Logistic(params[1:nrow(sitesXcounts)], 10)
      k.c.stockago.2 <- 10^params["AGO"]
      bgs.2 <- 10^params["bg"]

      } else {
      params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/", experiment, "/kds_with_structure/",
                   lim1, "-", lim1, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.1 <- Logistic(params[1:nrow(sitesXcounts)], 10)
      stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
                       "/SolveForKds/k_c_stockago.txt"), row.names=1,
                       header=TRUE, sep="\t")
      k.c.stockago.1 <- stockago[mirna,experiment]

      bgs.1 <- 10^params["bg"]


      params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/", experiment, "/kds_with_structure/",
                   lim2, "-", lim2, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.2 <- Logistic(params[1:nrow(sitesXcounts)], 10)
      stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
                       "/SolveForKds/k_c_stockago.txt"), row.names=1,
                       header=TRUE, sep="\t")
      k.c.stockago.2 <- stockago[mirna,experiment]

      bgs.2 <- 10^params["bg"]
}
    print(kds.1)
    print(kds.2)
    xmin <- 10^floor(min(log10(kds.1), log10(kds.2)))
    xmax <- 10^ceiling(max(log10(as.numeric(unlist(c(kds.1, kds.2, bgs.2, bgs.1, k.c.stockago.2, k.c.stockago.1))))))
    print(xmin)
    print(xmax)
    xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
    xs <- xs[xs >= xmin & xs <= xmax]
    xmin <- min(xs)
    xmax <- max(xs)
    xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))


    plot(c(kds.1, bgs.1, k.c.stockago.1),
      c(kds.2, bgs.2, k.c.stockago.2),
      col = kSiteColors[c(rownames(sitesXcounts), "bg", "Ago"),],
      pch = 19,
      log = 'xy',
      xlim = c(xmin,xmax*(xmax/xmin)^0.55), 
      ylim = c(xmin, xmax),
      axes = FALSE,
      ann = FALSE)
    segments(xmin, xmin, xmax, xmax, lty = 2)
    axis(1, at=xl,
         labels=sapply(xl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=xmin, lwd=0, cex.axis=1.5, padj=0.2)
    axis(1, at=xs, labels=FALSE,
         pos=xmin, lwd=2)
    # Label the axis at each order of magnitude.

    axis(2, at=xl,
         labels=sapply(xl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
    axis(2, at=xs, labels=FALSE,
         pos=xmin, lwd=2)

        title(main = mirna, font.main=1, cex.main=1.5, line=-2, adj=0.1)
    title(xlab = paste0(lim1, " parameters (nM)"), cex.lab=1.5, line=2, adj=0.3)
    title(ylab = paste0(lim2, " parameters (nM)"), cex.lab=1.5, line=2)
    legend.names <- c(rownames(sitesXcounts), "bg", "Ago")
    # centered_sites <- c("11mer-m3.13", "12mer-m3.14",
    #                 "11mer-m4.14", "12mer-m4.15")
    # centered_index <- sapply(centered_sites, function(site){
    #                        which(legend.names == site)
    #                      }
    #                      )
    # legend.names[centered_index] <- "Centered"
    # legend.names <- unique(legend.names)
    legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=0.9, bty="n", ncol = 1)

    
}

MakeScatterAcrossOrder <- function(mirna, experiment, lim, sitelist1, sitelist2, method) {
    sitesXcounts <- GetSitesXCounts(mirna, experiment, lim, lim, sitelist1)
    colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
      colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
        return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
      }
    )
    sitesXcounts1 <- sitesXcounts
        sitesXcounts <- GetSitesXCounts(mirna, experiment, lim, lim, sitelist2)
    colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
      colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
        return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
      }
    )
    sitesXcounts2 <- sitesXcounts
    print(sitesXcounts2)
    if (method == "free_protein") {

      params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/", experiment, "/kds_with_structure/",
                   lim, "-", lim, "_", 1, "-", 15, "_", sitelist1, "_Basic_singlebg_free_protein.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.1 <- Logistic(params[1:nrow(sitesXcounts1)], 10)
      k.c.stockago.1 <- 10^params["AGO"]
      bgs.1 <- 10^params["bg"]

      params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/", experiment, "/kds_with_structure/",
                   lim, "-", lim, "_", 1, "-", 15, "_", sitelist2, "_Basic_singlebg_free_protein.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.2 <- Logistic(params[1:nrow(sitesXcounts2)], 10)
      k.c.stockago.2 <- 10^params["AGO"]
      bgs.2 <- 10^params["bg"]

      } else {
      params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/", experiment, "/kds_with_structure/",
                   lim, "-", lim, "_", 1, "-", 15, "_", sitelist1, "_Basic_singlebg.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.1 <- Logistic(params[1:nrow(sitesXcounts1)], 10)
      stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
                       "/SolveForKds/k_c_stockago.txt"), row.names=1,
                       header=TRUE, sep="\t")
      k.c.stockago.1 <- stockago[mirna,experiment]

      bgs.1 <- 10^params["bg"]


      params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                   "/", experiment, "/kds_with_structure/",
                   lim, "-", lim, "_", 1, "-", 15, "_", sitelist2, "_Basic_singlebg.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.2 <- Logistic(params[1:nrow(sitesXcounts2)], 10)
      stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
                       "/SolveForKds/k_c_stockago.txt"), row.names=1,
                       header=TRUE, sep="\t")
      k.c.stockago.2 <- stockago[mirna,experiment]

      bgs.2 <- 10^params["bg"]
}
    names(kds.1) <- rownames(sitesXcounts1)
    names(kds.2) <- rownames(sitesXcounts2)
    kds.2 <- kds.2[names(kds.1)]
    print(kds.1)
    print(kds.2)
    xmin <- 10^floor(min(log10(kds.1), log10(kds.2)))
    xmax <- 10^ceiling(max(log10(as.numeric(unlist(c(kds.1, kds.2, bgs.2, bgs.1, k.c.stockago.2, k.c.stockago.1))))))
    print(xmin)
    print(xmax)
    print(bgs.1)
    print(bgs.2)
    print(k.c.stockago.1)
    print(k.c.stockago.2)
    xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
    xs <- xs[xs >= xmin & xs <= xmax]
    xmin <- min(xs)
    xmax <- max(xs)
    xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))


    plot(c(kds.1, bgs.1, k.c.stockago.1),
      c(kds.2, bgs.2, k.c.stockago.2),
      col = kSiteColors[c(rownames(sitesXcounts1), "bg", "Ago"),],
      pch = 19,
      log = 'xy',
      xlim = c(xmin,xmax*(xmax/xmin)^0.55), 
      ylim = c(xmin, xmax),
      axes = FALSE,
      ann = FALSE)
    segments(xmin, xmin, xmax, xmax, lty = 2)
    axis(1, at=xl,
         labels=sapply(xl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=xmin, lwd=0, cex.axis=1.5, padj=0.2)
    axis(1, at=xs, labels=FALSE,
         pos=xmin, lwd=2)
    # Label the axis at each order of magnitude.

    axis(2, at=xl,
         labels=sapply(xl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=xmin, las=2, lwd=0, cex.axis=1.5,hadj=0.8)
    axis(2, at=xs, labels=FALSE,
         pos=xmin, lwd=2)

        title(main = mirna, font.main=1, cex.main=1.5, line=-2, adj=0.1)
    title(xlab = paste0(sitelist1, " parameters (nM)"), cex.lab=1.5, line=2, adj=0.3)
    title(ylab = paste0(sitelist2, " parameters (nM)"), cex.lab=1.5, line=2)
    legend.names <- c(rownames(sitesXcounts1), "bg", "Ago")
    # centered_sites <- c("11mer-m3.13", "12mer-m3.14",
    #                 "11mer-m4.14", "12mer-m4.15")
    # centered_index <- sapply(centered_sites, function(site){
    #                        which(legend.names == site)
    #                      }
    #                      )
    # legend.names[centered_index] <- "Centered"
    # legend.names <- unique(legend.names)
    legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=0.9, bty="n", ncol = 1)

    
}






MakeSiteIterationPlotKineticsSingle <- function(out,method,mirna, num.kds, num.bgs, colors = FALSE) {
  setEPS()
  postscript(file=paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/figures/koffs/",
    mirna,"/", method, ".eps"))
  par(kPlotParameters)
  x = seq(1,dim(out)[1],length = max(ncol(out),1000))
  out <- out[x,]
  # print(out)
  probs = out[ ,"-logp"]
  out.print <- out[ , seq(dim(out)[2] - 1)]
  out.print.koffs <- Logistic(out[,1:num.kds],max = 10)
  out.print.setpoint <- 10^out[,num.kds + 1]
  out.print.bgs <- 10^out.print[,(num.kds + 2) : (num.kds + num.bgs + 1)]

  out.print <- cbind(out.print.koffs,out.print.setpoint, out.print.bgs)
  ys <- 10^c(max(floor(log10(min(out.print))), -6), ceiling(log10(max(out.print))))
  probs.scaled <- probs / max(probs) * (log10(ys[2]) - log10(ys[1])) + log10(ys[1])


  plot(x    = x,
       y    = 10^probs.scaled,
       log  = "y",
       axes = FALSE,
       ann = FALSE, 
       type = "l",
       col = "white",
       ylim = ys)
  axis(side = 1,
       at   = seq(1, max(x), by=max(1, floor(max(x) / 20))),
       pos  = ys[1],
       tck  = -0.01,
       labels = FALSE)
  axis(1, at=seq(1, max(x), by=max(1, floor(max(x) / 20))*5), pos=ys[1], lwd=2)
  axis(2, at=10^seq(log10(ys[1]), log10(ys[2])),
       labels=sapply(seq(log10(ys[1]), log10(ys[2])), function(name) {
         eval(substitute(expression(10^x), list(x=name)))
       }),
       pos=1, lwd=2)
  title(main = mirna, font=1)
  axis(2, at=c(sapply(seq(log10(ys[1]), log10(ys[2])-1), function(x) seq(9)*10^x),ys[2]), labels=FALSE,
       pos=1)
  if (colors != FALSE) {

  cols <- c(colors,kSiteColors[colnames(out.print)[(length(colors)+1):length(colnames(out.print))],])
  names(cols) <- colnames(out.print)
  cols["AGO"] <- "grey"

  sapply(colnames(out.print),function(name) {
    lines(x, out.print[, name], lwd=1, col=cols[name])
    })
  } else  {
    sapply(colnames(out.print),function(name) {
    lines(x, out.print[, name], lwd=1, col=kSiteColors[name,])
    })
  }
  lines(x, 10^probs.scaled, type="l", lwd=2, col="black")
  dev.off()
}


MakeSiteIterationPlotKineticsSingleNoSetPoint <- function(out,method,mirna, num.kds, num.bgs, colors = FALSE) {
  setEPS()
  postscript(file=paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/figures/koffs/",
    mirna,"/", method, ".eps"))
  par(kPlotParameters)
  x = seq(1,dim(out)[1],length = max(ncol(out),1000))
  out <- out[x,]
  # print(out)
  probs = out[ ,"-logp"]
  out.print <- out[ , seq(dim(out)[2] - 1)]
  out.print.koffs <- Logistic(out[,1:num.kds],max = 10)
  out.print.bgs <- 10^out.print[,(num.kds + 1) : (num.kds + num.bgs)]

  out.print <- cbind(out.print.koffs,out.print.bgs)
  ys <- 10^c(max(floor(log10(min(out.print))), -5), 5)
  probs.scaled <- probs*0.9 / probs[1] * (log10(ys[2]) - log10(ys[1])) + log10(ys[1])

  colnames <- c(rep(sapply(colnames(out.print)[1:num.kds], function(name) {
    trim = unlist(strsplit(name,split = "_"))
    return(trim[1])
    }), 3), colnames(out.print)[(3*num.kds + 1):length(colnames(out.print))])
  names(colnames) <- colnames(out.print)

  plot(x    = x,
       y    = 10^probs.scaled,
       log  = "y",
       axes = FALSE,
       ann = FALSE, 
       type = "l",
       col = "white",
       ylim = ys)
  axis(side = 1,
       at   = seq(1, max(x), by=max(1, floor(max(x) / 20))),
       pos  = ys[1],
       tck  = -0.01,
       labels = FALSE)
  axis(1, at=seq(1, max(x), by=max(1, floor(max(x) / 20))*5), pos=ys[1], lwd=2)
  axis(2, at=10^seq(log10(ys[1]), log10(ys[2])),
       labels=sapply(seq(log10(ys[1]), log10(ys[2])), function(name) {
         eval(substitute(expression(10^x), list(x=name)))
       }),
       pos=1, lwd=2)
  title(main = mirna, font=1)
  axis(2, at=c(sapply(seq(log10(ys[1]), log10(ys[2])-1), function(x) seq(9)*10^x),ys[2]), labels=FALSE,
       pos=1)
  if (colors != FALSE) {

  cols <- c(colors,kSiteColors[colnames(out.print)[(length(colors)+1):length(colnames(out.print))],])
  names(cols) <- colnames(out.print)
  cols["AGO"] <- "grey"

  sapply(colnames(out.print),function(name) {
    lines(x, out.print[, name], lwd=1, col="blue")
    })
  } else  {
    sapply(colnames(out.print),function(name) {
    lines(x, out.print[, name], lwd=1, col=kSiteColors[colnames[name],])
    })
  }
  lines(x, 10^probs.scaled, type="l", lwd=2, col="black")
  dev.off()
  dev.set(3)
}

MakeSiteIterationPlotKineticsDoubleNoSetPoint <- function(out,method,mirna, num.sites, inds_double, num.bgs, colors = FALSE) {
  setEPS()
  postscript(file=paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/figures/koffs/",
    mirna,"/", method, ".eps"))  
  par(kPlotParameters)
  x = seq(1,dim(out)[1],length = min(nrow(out),1000))
  out <- out[x,]
  probs = out[ ,"-logp"]
  out.print <- out[ , seq(ncol(out) - 1)]
  out.print.koffs1 <- 10^out[,1 : num.sites]
  out.print.koffs2 <- 10^out[,(num.sites + 1) : (2 * num.sites)]
  out.print.fwd <- 10^out[,(2 * num.sites + 1): (2 * num.sites + length(inds_double))]
  out.print.rev <- 10^out[,(2 * num.sites + length(inds_double) + 1) : (2 * (num.sites + length(inds_double)))]

  out.print.bgs <- 10^out.print[,(2 * (num.sites + length(inds_double)) + 1) : (2 * (num.sites + length(inds_double)) + num.bgs)]

  out.print <- cbind(out.print.koffs1, out.print.koffs2, out.print.fwd, out.print.rev,out.print.bgs)
  print(dim(out.print))
  ys <- 10^c(max(floor(log10(min(out.print))), -20), ceiling(log10(max(out.print))))
  probs.scaled <- probs/(ys[2])
  colnames <- c(rep(sapply(colnames(out.print)[1:num.sites], function(name) {
    trim = unlist(strsplit(name,split = "_"))
    return(trim[1])
    }), 2),
    rep(sapply(colnames(out.print)[1:num.sites], function(name) {
    trim = unlist(strsplit(name,split = "_"))
    return(trim[1])
    })[inds.double], 2),
    colnames(out.print)[(2 * (num.sites + length(inds_double)) + 1):length(colnames(out.print))])
  print(length(colnames))
  names(colnames) <- colnames(out.print)
  out.print <<- out.print
  plot(x    = x,
       y    = 10^probs.scaled,
       log  = "y",
       axes = FALSE,
       ann = FALSE, 
       type = "l",
       col = "white",
       ylim = ys)
  axis(side = 1,
       at   = seq(1, max(x), by=max(1, floor(max(x) / 20))),
       pos  = ys[1],
       tck  = -0.01,
       labels = FALSE)
  axis(1, at=seq(1, max(x), by=max(1, floor(max(x) / 20))*5), pos=ys[1], lwd=2)
  axis(2, at=10^seq(log10(ys[1]), log10(ys[2])),
       labels=sapply(seq(log10(ys[1]), log10(ys[2])), function(name) {
         eval(substitute(expression(10^x), list(x=name)))
       }),
       pos=1, lwd=2)
  title(main = mirna, font=1)
  axis(2, at=c(sapply(seq(log10(ys[1]), log10(ys[2])-1), function(x) seq(9)*10^x),ys[2]), labels=FALSE,
       pos=1)
  print("2821")
  if (colors != FALSE) {

  cols <- c(colors,kSiteColors[colnames(out.print)[(length(colors)+1):length(colnames(out.print))],])
  names(cols) <- colnames(out.print)
  cols["AGO"] <- "grey"

  sapply(colnames(out.print),function(name) {
    lines(x, out.print[, name], lwd=1, col=cols[name])
    })
  } else  {
    lwds <- c(rep(1, 2*num.sites), rep(1, 2*length(inds.double)), rep(1, num.bgs))
    ltys <- c(rep(c(1, 2),each=num.sites),rep(c(3, 3),each=length(inds.double)),rep(1, num.bgs))
    names(lwds) <- colnames(out.print)
    names(ltys) <- names(lwds)
    lwds <<- lwds
    ltys <<- ltys
    sapply(colnames(out.print),function(name) {

    lines(x, out.print[, name], lwd=lwds[name],lty = ltys[name] , col=kSiteColors[colnames[name],])
    })
  }
  print("2842")
  lines(x, probs.scaled, type="l", lwd=2, col="black")
  dev.off()
  dev.set(2)
}
