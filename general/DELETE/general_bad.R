options(warn=0)

library(beeswarm)
library(gplots)
library(wrswoR)
mirnas_all <- c("miR-1", "let-7a", "miR-155", "miR-124", "lsy-6")

centered_sites <- c("11mer-m3.13", "12mer-m3.14",
                    "11mer-m4.14", "12mer-m4.15")

canonical_sites <- c("8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-A1", "6mer-m8")

kmer_list_names <- c("8mers", "9mers", "10mers", "11mers", "12mers")

exp_conditions <- c("I", "I_combined", "0", "0.4", "1.26", "4", "12.6", "40")

exp_cond_colors <- c("black", "black", "gray60", "red", "orange", "forestgreen", "blue", "violet")
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

GetFlankbgKds <- function(mirna, experiment, n_constant, sitelist, site) {
    params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                 experiment, "/kds_PAPER/", n_constant, "_", 
                 sitelist, "_", site, "_controlplfold_PAPER.txt")

    params <- data.frame(read.table(params.file, header=TRUE, row.names=1,
                         stringsAsFactors=FALSE))
    return(params)
}


GetFlankbyPlFoldKds <- function(mirna, experiment, n_constant, sitelist, site) {
    params.file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                 experiment, "/kds_PAPER/", n_constant, "_", 
                 sitelist, "_", site, "_withprob_PAPER.txt")

    params <- data.frame(read.table(params.file, header=TRUE, row.names=1,
                         stringsAsFactors=FALSE))
    return(params)
}

PlotFlankPlFoldKds <- function(mirna, experiment, n_constant, sitelist, site) {
  kds <- GetFlankKds(mirna, experiment, n_constant, sitelist, site)
  kds.plot <- kds$Mean
  kds.sub <- GetFlankbyPlFoldKds(mirna, experiment, n_constant, sitelist, site)
  sub.divisions <- unique(sapply(rownames(kds.sub), function(name) {
    temp_name <- paste0("_",unlist(strsplit(name, split = "_"))[2], collapse = "")
    }))
  print(sub.divisions)
  par(kPlotParameters)
  average_vector <- matrix(NaN,nrow=length(sub.divisions), ncol=length(kds.plot))
  plot(kds.plot, rnorm(length(kds.plot),mean=0,sd=0.1),
       log='x', xlim = c(0.0001, 10), ylim = c(0, length(sub.divisions)+5),
       col=sapply(rownames(kds), GetColorFunction))
  starting_sd <- sd(log10(kds.plot))^2
  text(10,0,sd(log10(kds.plot))^2/starting_sd)
  sapply(1:length(sub.divisions), function(ind) {
    print(ind)
    print(sub.divisions)
    print(sub.divisions[ind])
    print(paste0(sub.divisions[ind],"$",collapse=""))
    sub.division <- paste0(sub.divisions[ind],"$",collapse="")
    print(sub.division)
    kds.plot <- kds.sub[grep(sub.division,rownames(kds.sub),perl=TRUE),]$Median
    print(length(kds.plot))
    print(length(log10(kds.plot)-mean(log10(kds.plot))))
    print(length(average_vector[ind,]))
    average_vector[ind,] <<- log10(kds.plot)-mean(log10(kds.plot))+mean(log10(kds$Mean))
    kd.plot_alt <- 10^(log10(kds.plot)-mean(log10(kds.plot))+mean(log10(kds$Mean)))
    points(kds.plot, ind+rnorm(length(kds.plot),mean=0,sd=0.1), col=sapply(rownames(kds), GetColorFunction))
     text(10,ind,round(sd(log10(kds.plot))^2/starting_sd,2))
  })
  final <- 10^colMeans(average_vector)
  points(10^colMeans(average_vector),
         length(sub.divisions)+2+rnorm(length(kds.plot), mean=0, sd=0.1),
         col=sapply(rownames(kds), GetColorFunction))
       text(10,length(sub.divisions)+2,round(sd(log10(final))^2/starting_sd,2))

}

PlotFlankControlKds <- function(mirna, experiment, n_constant, sitelist, site) {
  kds <- GetFlankKds(mirna, experiment, n_constant, sitelist, site)
  kds.plot <- kds$Mean
  kds.sub <- GetFlankbgKds(mirna, experiment, n_constant, sitelist, site)$Mean

  par(kPlotParameters)
  # par(mfrow=c(2,1))
  # plot(kds.plot, 1+rnorm(length(kds.plot),mean=0,sd=0.1),
  #      log='x', xlim = c(0.0001, 10), ylim = c(0, 5),
  #      col=sapply(rownames(kds), GetColorFunction))

  # points(kds.sub, 2 + rnorm(length(kds.plot),mean=0,sd=0.1),
  #        col=sapply(rownames(kds), GetColorFunction))
      xmin = 0.0001
      xmax = 10
        plot( kds.sub, kds.plot,
       log='xy', xlim = c(0.0001, 10), ylim = c(0.0001, 10),
       col=sapply(rownames(kds), GetColorFunction))
      linear_model <- lm(log10(kds.plot) ~ log10(kds.sub))
      m <- linear_model$coefficients[2]
      b <- linear_model$coefficients[1]

    x_line <- 10^seq(log10(xmin), log10(xmax),length = 20)
    y_line <- 10^(m*log10(x_line) + b)

    lines(x_line, y_line, lty = 2,lwd = 0.5)
    text(1e-1, 1e-3, bquote(log(y) == log(x)*.(round(m,2)) + .(round(b,2))))
    text(1e-1, 2e-3, round(cor(log(kds.sub), log(kds.plot))^2,2))
    # text(x = 0.001, y = 0.95, bquote(italic(r)^2 ==  .(cor_text)))


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


GetModel <- function(mirna, experiment, n_constant, sitelist) {
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

  c.I.tots <- Norm(sitesXcounts[,2])*100
  names(c.I.tots) <- rownames(sitesXcounts)
  data <- sitesXcounts[,3:7]
  names(c.I.tots) <- rownames(data)
  x_points <- c(40,12.65,4,1.265,0.4)

  c.totals <- matrix(c.I.tots, nrow=nrow(data), ncol = length(x_points), byrow=FALSE)
  colnames(c.totals) <- colnames(x_points)
  rownames(c.totals) <- rownames(data)

  c.agos <- sapply(x_points, function(x) {
      as.numeric(x) * k.c.stockago / 100
    }
  )
  names(c.agos) <- x_points

  c.bounds <- as.matrix(
    sapply(c.agos, function(x) {
        return(GetBoundRNA(kds, c.I.tots, x))
      }
    )
  )
  return(c.bounds / c.totals)

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
 