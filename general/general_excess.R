kAgoStock <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
                       "/SolveForKds/k_c_stockago.txt"), row.names=1,
                       header=TRUE, sep="\t")
# Repression files:
salomon.df <- data.frame(read.table("outside_data/salomon2015.txt", header=1,
                                    sep="\t", stringsAsFactors=FALSE))
# Energy files:
dG.df <- read.table("canonical_sites_mfe.txt")
dG.df.new <- read.table("canonical_sites_mfe_new.txt")
colnames(dG.df) <- gsub("(^.*)\\.(.*$)", colnames(dG.df), replace="\\1-\\2")
colnames(dG.df.new) <- gsub("(^.*)\\.(.*$)", colnames(dG.df.new), replace="\\1-\\2")
canon_2.7 <- c(2, 4, 6)

dG.df <- read.table("canonical_sites_mfe.txt")

UTR_h <- read.table("RepressionData/UTRs_hg19.txt",
                    header=FALSE,
                    row.names=NULL,
                    stringsAsFactors=FALSE,
                    sep="\t")
UTR_h <- structure(UTR_h[,2], names=UTR_h[,1])

UTR_m <- read.table("RepressionData/UTRs_mm10.txt",
                    header=FALSE,
                    row.names=NULL,
                    stringsAsFactors=FALSE,
                    sep="\t")
UTR_m <- structure(UTR_m[,2], names=UTR_m[,1])

UTRS <- list(human=UTR_h, mouse=UTR_m)
GetUTRS <- function(list, species="human") {
	UTRS[[species]][list] 
}


UTR_kathy <- read.table("RepressionData/Lin-Shi_transfection_data/utr3.txt")



FindUTRSite <- function(utrs, mirna, sites=c("7mer-A1", "7mer-m8")) {
	sevenA1 <- MirnaTargetSequence(mirna, 1, 7)
	sevenM8 <- MirnaTargetSequence(mirna, 2, 8)
	targets <- unique(unlist(sapply(c(sevenA1, sevenM8), grep, x=utrs)))
	return(targets)
}


GetLinShiRepressionData <- function(best=FALSE, old=FALSE, threePseq=TRUE) {
	if (best) str1 <- "best_"
	else str1 <- "all_"
	if (old) old <- "_old"
	else old <- ""
	if (threePseq) threePseq <- "_3pseq"
	else threePseq <- ""
	path <- paste0("RepressionData/Lin-Shi_transfection_data/", str1,
	               "flanking_kds_and_repression", threePseq, old, ".txt")
	fread(path, data.table=FALSE)
}

LinModelSeedPairing <- function(dG.table=3) {
  dG.df <- read.table(file=paste0("canonical_sites_mfe_", dG.table, ".txt"))
  colnames(dG.df) <- gsub("(^.*)\\.(.*$)", colnames(dG.df), replace="\\1-\\2")
  R <- 1.987e-3 # in kcal K-1 mol-1
  T <- 310.15 # in K

  kds <- GetSiteKds(mirna)
  names_kds <- rownames(kds)
  output <- kds["7mer-A1_Kd",]$Mean/kds["7mer-m8_Kd", ]$Mean 
  return(output)
}


GetRisslandRepData <- function(line="HEK293", qn=FALSE) {
	if (qn) {
		qn <- ".qn"
		skip=1
	} else {
		qn <- ""
		skip=0
	}
	path1 <- "RepressionData/Rissland-Nam_transfection_data/GSE52530_"
	path2 <- ".expData"
	suffix <- ".txt"
	fullpath <- paste0(path1, line, path2, qn, suffix)
	table <- data.frame(fread(fullpath, header=FALSE, sep="\t",
	                             stringsAsFactors=FALSE, showProgress=FALSE))
	out <- data.frame(gene=table[, 1], len=table[, 2],
	           ctr.1=as.numeric(gsub("^(.*),(.*)$", table[, 3], replace="\\1", perl=TRUE)),
	           ctr.2=as.numeric(gsub("^(.*),(.*)$", table[, 3], replace="\\2", perl=TRUE)),
	           m124.1=as.numeric(gsub("^(.*),(.*)$", table[, 4], replace="\\1", perl=TRUE)),
	           m124.2=as.numeric(gsub("^(.*),(.*)$", table[, 4], replace="\\2", perl=TRUE)),
	           m155.1=as.numeric(gsub("^(.*),(.*)$", table[, 5], replace="\\1", perl=TRUE)),
	           m155.2=as.numeric(gsub("^(.*),(.*)$", table[, 5], replace="\\2", perl=TRUE)),
	           stringsAsFactors=FALSE)
	colnames(out) <- c("gene", "length", "control.1", "control.2", "miR-124_1",
	                     "miR-124_2", "miR-155_1", "miR-155_2")
	out
}


OverlappingGenes <- function(df.1, df.2) {
	intersect(df.1$gene, df.2$gene)
}

GetORFoldChange <- function(line, mirna, qn=FALSE, cutoff=3) {
	if (mirna %in% c("miR-155", "miR-124")) {
	args <- as.list(match.call())
	rep.df <- SubfunctionCall(GetRisslandRepData, args)
	print(dim(rep.df))
	cutoff.inds <- which(rep.df$control.1*rep.df$length/1000 > cutoff & rep.df$control.2*rep.df$length/1000 > cutoff)
	print(length(cutoff.inds))
	control <- rowMeans(rep.df[,grep("control", colnames(rep.df))])
	signal <- rowMeans(rep.df[,grep(mirna, colnames(rep.df))])
	genes <- rep.df$gene
	utrs <- GetUTRS(rep.df$gene)
	target.inds <- FindUTRSite(utrs, "miR-155")
	print(length(target.inds))
	final.inds <- union(cutoff.inds, target.inds)
	out <- signal/control
	names(out) <- rep.df$gene
	return(out[final.inds])
	return(out)
	} else {
		print("not a miRNA studied by Olivia")
		break
	}
}

PlotTwoFoldChanges <- function(line1, line2, mirna, qn=FALSE, cutoff=0,
                               xpos=20, ypos=20) {
	width <- 5
	height <- 5
	args <- as.list(match.call())
	fc.1 <- SubfunctionCall(GetORFoldChange, args, line=line1)
	fc.2 <- SubfunctionCall(GetORFoldChange, args, line=line2)
	inds.both <- intersect(names(fc.1), names(fc.2))
	fc.1 <- fc.1[inds.both]
	fc.2 <- fc.2[inds.both]
	print(head(fc.1))
	print(head(fc.2))
	xmin <- -4
	xmax <- 2
	ymin <- -4
	ymax <- 2
	dev.new(xpos=xpos, ypos=ypos, height=height, width=width)
	BlankPlot()
	abline(h=0)
	abline(v=0)
	AddLinearAxis(1, tick.space=1, label.space=1, label=line1)
	AddLinearAxis(2, tick.space=1, label.space=1, label=line2)
	points(log2(fc.1), log2(fc.2), col=ConvertRColortoRGB("black", alpha=0.2), pch=1)
}

GetRisslandRepressionData <- function(best=FALSE, old=FALSE, threePseq=TRUE) {
	if (best) str1 <- "best_"
	else str1 <- "all_"
	if (old) old <- "_old"
	else old <- ""
	if (threePseq) threePseq <- "_3pseq"
	else threePseq <- ""
	path <- paste0("RepressionData/Lin-Shi_transfection_data/", str1,
	               "flanking_kds_and_repression", threePseq, old, ".txt")
	fread(path, data.table=FALSE)
}


GetAnalysisPath <- function(mirna, experiment, condition, analysis_type, ext="",
                            suffix="txt") {
  subdirectory <- paste(mirna, experiment, analysis_type, sep="/")
  directory <- paste0(kSolexaDir, subdirectory)
  full_path <- paste0(directory, "/", condition, ext, ".", suffix)
  return(full_path)
}



MakeKmerList <- function(n) {
  if (n == 1)
    sort(kDNucs)
  else
    c(sapply(MakeKmerList(n - 1), function(x) {
      paste0(x, sort(kDNucs))
    }))
}
print(MakeKmerList(1))
print(MakeKmerList(2))
print(MakeKmerList(3))

GetInputKmers <- function(mirna, experiment, condition, n_constant, kmer_length, sites_removed) {
  ext <- sprintf("_%s_k%s_%ssites", n_constant, kmer_length, sites_removed)
  path <- SubfunctionCall(GetAnalysisPath, analysis_type = "kmers")
  print(path)
  out <- fread(path, sep="\t", fill=TRUE,header=FALSE,
                       stringsAsFactors=FALSE, showProgress=FALSE, data.table=FALSE)
  structure(out[, 2], names = out[,1])
}

KmerMatrix <- function(mirna, experiment, n_constant, kmer_length, sites_removed) {
  kmers <- MakeKmerList(kmer_length)
  conditions <- c("I", "40", "12.6", "4", "1.26", "0.4", "0")
  output <- matrix(NaN, nrow=length(kmers), ncol=length(conditions),
                   dimnames=list(kmers, conditions))
  sapply(conditions, function(condition) {
    SubfunctionCall(GetInputKmers)
  })
}

PlotEnrichmentMatrix <- function(mirna, experiment, n_constant, n_kmer, n_sites, kmer) {
  matrix <- KmerMatrix(mirna, experiment, n_constant, n_kmer, n_sites)
  data.R <- EquilEnrichments(matrix[, -1], matrix[, 1])
  color <- rep(ConvertRColortoRGB("gray", alpha=0.5), nrow(matrix))
  ind.kmer <- grep(kmer, rownames(matrix))
  dev.new(xpos = 20, ypos=20, height=8, width=12)
  par(kPlotParameters)
  color[ind.kmer] <- "blue"
  par(mfrow=c(2, 3))
  xmin=1e-8
  xmax=1e-3
  ymin=5e-3
  ymax=200
  sapply(colnames(data.R), function(condition) {
    BlankPlot(log='xy')
    AddLogAxis(1, label="R")
    AddLogAxis(2, label="input freq")

    points(Norm(matrix[,1]), data.R[, condition], col=color)
    points(Norm(matrix[,1])[ind.kmer], data.R[ind.kmer, condition], col="blue")
    title(main=condition)  
  })
}

# PlotEnrichmentMatrix("lsy-6", "equilibrium", 5, 9, 5, "ATGACAAAA")

#23456789!123456789@123456789#123456789$123456789%123456789^123456789&123456789*
RemoveStringPrefix <- function(x, prefix) {
  gsub(paste0("^", prefix, "(.*)$"), x, replace="\\1", perl=TRUE)
}

#23456789!123456789@123456789#123456789$123456789%123456789^123456789&123456789*

#23456789!123456789@123456789#123456789$123456789%123456789^123456789&123456789*
SeqComplement <- function(x, RNA=FALSE) {
  # For a single string or a vector of strings, returns the complementary RNA
  # sequence (not reversed!).
  complement.base <- c(U="A", T="A", G="C", C="G")
  if (RNA) complement <- c(complement.base, A="U")
  else complement <- c(complement.base, A="T")
  return(sapply(lapply(strsplit(x, NULL), function(y) complement[y]), paste,
                collapse=""))
}

RevComplement <- function(x, RNA=FALSE) {
  StrRev(SeqComplement(x, RNA=RNA))
}

MirnaTargetSequence <- function(mirna, start, stop, RNA=FALSE) {
  # Gives the reverse-compement sequence to the miRNA positions, using R
  # indexing (not pythonic).
  mirna.seq <- substr(kMirnaSeqs[mirna], start, stop)
  return(RevComplement(mirna.seq, RNA=RNA))
}

GetSingleFlankPosition <- function(flanks, position) {
  if (position > 2 & substr(c(flanks)[1], 3, 3) %in% c(".", "|")) {
    position <- position + 1
  }
  return(sapply(flanks, substr, start=position, stop=position))
}

# Basic mathematical functions:
Norm <- function(vector) {
  # Returns the normalized entries of the vector.
  vector/sum(vector)
}

MatNorm <- function(matrix, dim=2) {
  # Returns the normalized entries of the vector.
  apply(matrix, dim, Norm)
}

GeoMean <- function(vector) {
  exp(mean(log(vector), na.rm=TRUE))
}

Cumul <- function(vector) {
  # Returns the cumulative sum of a vector across its length.
  len <- length(vector)
  if (len == 1) return(vector)
  else return(c(Cumul(vector[1:(len - 1)]), sum(vector)))
}

Logistic <- function(vector, max) {
  return(max/(1 + exp(-vector)))
}

#23456789!123456789@123456789#123456789$123456789%123456789^123456789&123456789*
Logit <- function(vector, max) {
  return(-log(max/vector - 1))
}

# 1.________
GetAverageOfReplicates <- function(time, times, data) {
  return(rowMeans(data[, which(times == time), drop=FALSE]))
}

ZScore <- function(x) {
  (x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
}

ZScoreMatrix <- function(matrix, dim) {
  t(apply(matrix, dim, ZScore))
}
####################### I/O ####################################################
## Functions that load data. (__ total) # 1.___________
SitesXCounts <- function(mirna, experiment="equilibrium", n_constant=5,
                         sitelist="paper", mirna.start=NULL,
                         multisite=FALSE) {
  if (multisite) site.dir1 <- "/multisite_count_tables"
  else           site.dir1 <- "/site_count_tables"
  site.dir <- paste0(site.dir1, "/")
  file.prefix <- paste0(kSolexaDir, mirna, "/", experiment, site.dir, "all_sites_",
                        n_constant, "_", sitelist)
  if (sitelist %in% kKmerSiteLists) {
    file.suffix <- paste0("_", mirna.start, "-", mirna.start + 3, ".txt")
  } else {
    file.suffix <- ".txt"
  }
  file <- paste0(file.prefix, file.suffix)
  sXc <- read.table(file, sep="\t", header=TRUE, row.names=1,
                    stringsAsFactors=FALSE)
  sXc <- sXc[, grep("Seq", colnames(sXc), inv=TRUE)]
  colnames(sXc) <- gsub("^A", colnames(sXc), replace="", perl=TRUE)
  # Check the two conditions, that the site exists in the input library:
  sXc <- sXc[rowSums(sXc, na.rm=TRUE) > 0,]
  if (!(multisite)) sXc <- sXc[sXc[, 2] > 0,]
  return(sXc)
}

SitesXCounts12mersOld <- function(mirna, experiment="equilibrium", n_left=5,
                            n_right=5, sitelist="12mers", mirna.start=2) {
  extension <- sprintf("_%s-%s_%s_%s-%s", n_left, n_right, sitelist,
                       mirna.start, as.integer(mirna.start) + 3)
  condition <- "all_sites"
  path <- SubfunctionCall(GetAnalysisPath,
                          analysis_type="full_site_count_tables",
                          ext=extension)
  print(path)
  out <- data.frame(fread(path), row.names=1, stringsAsFactors=FALSE)
  colnames(out) <- gsub("^A(.*)$", colnames(out), replace="\\1", perl=TRUE)
  out
}

# Made FOR PAPER NOW
GetSiteKds <- function(mirna, experiment="equilibrium", n_constant=5,
                       sitelist="paper", mirna.start=FALSE,
                       combined.input=TRUE) {
  print(mirna.start)
  if (mirna.start) {
    str.mir.start <- sprintf("_%s-%s", mirna.start, as.integer(mirna.start) + 3)
  } else {
    str.mir.start <- ""
  }
  if (combined.input) {
  str.combined <- "_nocombInput"
  } else {
    str.combined <- ""
  }
  params.file <- paste0(kSolexaDir, mirna, "/",
               experiment, "/kds_PAPER/", n_constant, "_", 
               sitelist, str.mir.start, str.combined,"_PAPER.txt")
  params <- data.frame(read.table(params.file, header=TRUE, row.names=1,
                       sep= "\t", stringsAsFactors=FALSE))
  return(params)
}

KdVector <- function(mirna, experiment="equilibrium", n_constant=5,
                       sitelist="paper", mirna.start=FALSE,
                       combined.input=TRUE) {
  pars <- SubfunctionCall(GetSiteKds)
  kds <- pars[grep("_Kd", rownames(pars)), ]$Mean
  names(kds) <- gsub("^(.*)_Kd$", rownames(pars), replace="\\1",
                         perl=TRUE)[1 : length(kds)]
  kds
}


KdMatrix <- function(mirna, experiment="equilibrium", n_constant=5,
                       sitelist="paper", mirna.start=FALSE) {
  pars <- SubfunctionCall(GetSiteKds)
  pars <- pars[grep("_Kd", rownames(pars)), ]
  rownames(pars) <- gsub("^(.*)_Kd$", rownames(pars), replace="\\1",
                         perl=TRUE)
  pars
}

ConvertTtoU <- function(site) {
  return(gsub("T", c(site), replace="U"))
}


Get12merKdsOld <- function(mirna, experiment="equilibrium", n_left=5,
                            n_right=5, sitelist="12mers", mirna.start=2) {
  condition <- sprintf("%s-%s_%s_%s-%s_singlebg_multinomial_PAPER_last",
                       n_left, n_right, sitelist, mirna.start,
                       as.integer(mirna.start) + 3)
  path <- SubfunctionCall(GetAnalysisPath, analysis_type="kds_PAPER")
  print(path)
  t(fread(path, header=TRUE, sep="\t", stringsAsFactors))
}

Get12merKdsFull <- function(mirna, experiment="equilibrium", n_constant=5,
                            sitelist="12mers", mirna.start=2, combined.input=TRUE) {
  if (combined.input) {
    combined <- FALSE
    str.combined <- "_nocombInput"
  } else {
    combined <- TRUE
    str.combined <- ""
  }
  condition <- sprintf("%s_12mers_%s-%s%s_PAPER_full", n_constant, mirna.start, 
                       mirna.start + 3, str.combined)
  path <- SubfunctionCall(GetAnalysisPath, analysis_type="kds_PAPER")
  print(path)
  10^t(fread(path, header=TRUE, sep="\t", stringsAsFactors))
}

Get12merKdsMean <- function(mirna, experiment="equilibrium", n_constant=5,
                            sitelist="12mers", mirna.start=2, combined.input=TRUE) {
  if (combined.input) {
    combined <- FALSE
    str.combined <- "_nocombInput"
  } else {
    combined <- TRUE
    str.combined <- ""
  }
  condition <- sprintf("%s_12mers_%s-%s%s_PAPER_mean", n_constant, mirna.start, 
                       mirna.start + 3, str.combined)
  path <- SubfunctionCall(GetAnalysisPath, analysis_type="kds_PAPER")
  print(path)
  10^t(fread(path, header=TRUE, sep="\t", stringsAsFactors))
}

GetSiteSeqForKathy <- function(mirna, experiment="equilibrium", n_constant=5,
                               sitelist="paper") {
  condition <- sprintf("site_and_flank_kds_%s_%s",
                       n_constant, sitelist)
  path <- SubfunctionCall(GetAnalysisPath, analysis_type="kds_PAPER")
  print(path)
  t(fread(path, header=TRUE, sep="\t", stringsAsFactors))
}


Compare12merKds <- function(mirna, experiment="equilibrium", n_constant=5,
                            sitelist="12mers") {
  n_left = n_constant
  n_right = n_constant
  dev.new(xpos=20, ypos=20, height=10, width=15)
  par(mfrow=c(4, 5))
  sapply(seq(2, 5), function(mirna.start) {
    kd12mers_old <- SubfunctionCall(Get12merKdsOld)
    kd12mers_old <- kd12mers_old[1:(nrow(kd12mers_old) - 3), ]
    kd12mers_new <- SubfunctionCall(KdVector)
    kd12mers_full <- SubfunctionCall(Get12merKdsFull)
    kd12mers_full_nocombI <- SubfunctionCall(Get12merKdsFull, nocombI=TRUE)
    kd12mers_mean <- SubfunctionCall(Get12merKdsMean)
    kd12mers_mean_nocombI <- SubfunctionCall(Get12merKdsMean, nocombI=TRUE)
    shared <- intersect(names(kd12mers_old), names(kd12mers_new))
    x <- Logistic(kd12mers_old[shared], 10)
    y <- kd12mers_new[shared]
    y2 <- kd12mers_full[shared, ]
    y3 <- kd12mers_full_nocombI[shared, ]
    y4 <- kd12mers_mean[shared, ]
    y5 <- kd12mers_mean_nocombI[shared, ]
    xmin <- 1e-4
    xmax <- 200
    ymin <- xmin
    ymax <- xmax
    par(kPlotParameters)
    lapply(list(y2, y3, y4, y5), function(y) {
      BlankPlot(log='xy')
      AddLogAxis(1, label="Old 12mer Kds")
      AddLogAxis(2, label="New 12mer Kds")
      points(x, y, col=ConvertRColortoRGB("black", alpha=0.1))
      abline(0, 1)
    })
    plot(density(log(x)), ylim=c(0, 0.6))
    lines(density(log(y2)), col="red")
    lines(density(log(y3)), col="blue")
    lines(density(log(y4)), col="green")
    lines(density(log(y5)), col="purple")
    legend("topleft", bty="n",
           legend=c("old", "new_combined", "new_single", "new_mean_combined",
                    "new_mean_single"), lwd=1, col=c("black", "red", "blue",
                    "green", "purple"))
  })
}


Compare12merKdsmiR7 <- function(experiment="equilibrium2_nb", n_constant=5,
                                sitelist="12mers", type=FALSE, mirna.start=2,
                                combined.input=TRUE) {
  n_left = n_constant
  n_right = n_constant
  dev.new(xpos=20, ypos=20, height=8, width=10)
  par(mfrow=c(3, 3))
  sapply(c("miR-7-23nt", "miR-7-24nt", "miR-7-25nt"), function(mirna) {
    print(mirna)
    if (type=="full") {
      kd12mers_x <- SubfunctionCall(Get12merKdsFull, mirna=mirna)
    } else if (type=="new") {
      kd12mers_x <- SubfunctionCall(Get12merKdsNew, mirna=mirna)
    } else {
      kd12mers_x <- SubfunctionCall(KdVector, mirna=mirna)      
    }
    sapply(seq(3, 5), function(mirna.start) {
      if (type=="full") {
        kd12mers_y <- SubfunctionCall(Get12merKdsFull, mirna=mirna,
                                      mirna.start=mirna.start)
        namefunc <- rownames
      } else if (type=="new") {
        kd12mers_y <- SubfunctionCall(Get12merKdsNew, mirna=mirna,
                                      mirna.start=mirna.start)
        namefunc <- rownames
      } else {
        kd12mers_y <- SubfunctionCall(KdVector, mirna=mirna,
                                      mirna.start=mirna.start)      
        namefunc <- names
      }
      print(length(kd12mers_x))
      print(length(kd12mers_y))
      print(head(kd12mers_x))
      print(head(kd12mers_y))
      shared <- intersect(namefunc(kd12mers_x), namefunc(kd12mers_y))
      print(length(shared))
      if (type == "full" | type =="new"){
        x <- kd12mers_x[shared, ]
        y <- kd12mers_y[shared, ]        
      } else {
        x <- kd12mers_x[shared]
        y <- kd12mers_y[shared]        
      }
      # y4 <- kd12mers_mean_nocombI[, 1]
      xmin <- 1e-6
      xmax <- 2000
      ymin <- xmin
      ymax <- xmax
      par(kPlotParameters)
        BlankPlot(log='xy')
        AddLogAxis(1, label=paste0("nt ", as.integer(mirna.start) - 1, "-",
                                   as.integer(mirna.start) + 2, "12mer Kds"))
        AddLogAxis(2, label=paste0("nt ", mirna.start, "-",
                                   as.integer(mirna.start) + 3, "12mer Kds"))
        points(x, y, col=ConvertRColortoRGB("black", alpha=0.1))
        abline(0, 1)  
      xy <- GetPlotFractionalCoords(fx=0.1, fy=0.9, log='xy')
      text(xy[1], xy[2], paste0(mirna, "\n", type), adj=0)
      kd12mers_x <<- kd12mers_y
    })
  })
}

Compare12merKdsmiR72 <- function(experiment="equilibrium2_nb", n_constant=5,
                                sitelist="12mers", type=FALSE, mirna.start=2) {
  n_left = n_constant
  n_right = n_constant
  dev.new(xpos=20, ypos=20, height=8, width=10)
  par(mfrow=c(3, 3))
  sapply(seq(2, 5), function(mirna.start) {
    mirna <- "miR-7-23nt"
    print(mirna)
    if (type=="full") {
      kd12mers_x <- SubfunctionCall(Get12merKdsFull, mirna=mirna,
                                    mirna.start=mirna.start)
    } else if (type=="new") {
      kd12mers_x <- SubfunctionCall(Get12merKdsNew, mirna=mirna,
                                    mirna.start=mirna.start)
    } else {
      kd12mers_x <- SubfunctionCall(KdVector, mirna=mirna,
                                    mirna.start=mirna.start)      
    }
    sapply(c("miR-7-24nt", "miR-7-25nt"), function(mirna) {
      if (type=="full") {
        kd12mers_y <- SubfunctionCall(Get12merKdsFull, mirna=mirna,
                                      mirna.start=mirna.start)
        namefunc <- rownames
      } else if (type=="new") {
        kd12mers_y <- SubfunctionCall(Get12merKdsNew, mirna=mirna,
                                      mirna.start=mirna.start)
        namefunc <- rownames
      } else {
        kd12mers_y <- SubfunctionCall(KdVector, mirna=mirna,
                                      mirna.start=mirna.start)      
        namefunc <- names
      }
      print(length(kd12mers_x))
      print(length(kd12mers_y))
      print(head(kd12mers_x))
      print(head(kd12mers_y))
      shared <- intersect(namefunc(kd12mers_x), namefunc(kd12mers_y))
      print(length(shared))
      if (type == "full" | type =="new"){
        x <- kd12mers_x[shared, ]
        y <- kd12mers_y[shared, ]        
      } else {
        x <- kd12mers_x[shared]
        y <- kd12mers_y[shared]        
      }
      # y4 <- kd12mers_mean_nocombI[, 1]
      xmin <- 1e-4
      xmax <- 200
      ymin <- xmin
      ymax <- xmax
      par(kPlotParameters)
        BlankPlot(log='xy')
        AddLogAxis(1, label=paste0("nt ", as.integer(mirna.start) - 1, "-",
                                   as.integer(mirna.start) + 2, "12mer Kds"))
        AddLogAxis(2, label=paste0("nt ", mirna.start, "-",
                                   as.integer(mirna.start) + 3, "12mer Kds"))
        points(x, y, col=ConvertRColortoRGB("black", alpha=0.1))
        abline(0, 1)  
      xy <- GetPlotFractionalCoords(fx=0.1, fy=0.9, log='xy')
      text(xy[1], xy[2], paste0(mirna, "\n", type), adj=0)
      kd12mers_x <<- kd12mers_y
    })
  })
}


Compare12merKdsDatamiR7 <- function(experiment="equilibrium2_nb", n_constant=5,
                                sitelist="12mers", mirna.start=2, cond="I") {
  n_left = n_constant
  n_right = n_constant
  dev.new(xpos=20, ypos=20, height=8, width=10)
  par(mfrow=c(3, 3))
  sapply(c("miR-7-23nt", "miR-7-24nt", "miR-7-25nt"), function(mirna) {
    print(mirna)
    data_x <- SubfunctionCall(SitesXCounts, mirna=mirna,
                              mirna.start=mirna.start)
    data_x <<- data_x
    sapply(seq(3, 5), function(mirna.start) {
      data_y <- SubfunctionCall(SitesXCounts, mirna=mirna,
                                mirna.start=mirna.start)
      print(head(data_x))
      print(head(data_y))
      shared <- intersect(rownames(data_x), rownames(data_y))
      print(length(shared))
      x <- data_x[shared, cond]
      y <- data_y[shared, cond]        
      xmin <- 1
      xmax <- 1e7
      ymin <- xmin
      ymax <- xmax
      par(kPlotParameters)
        BlankPlot(log='xy')
        AddLogAxis(1, label=paste0("nt ", as.integer(mirna.start) - 1, "-",
                                   as.integer(mirna.start) + 2, "12mer counts"))
        AddLogAxis(2, label=paste0("nt ", mirna.start, "-",
                                   as.integer(mirna.start) + 3, "12mer counts"))
        points(x, y, col=ConvertRColortoRGB("black", alpha=0.1))
        abline(0, 1)  
      xy <- GetPlotFractionalCoords(fx=0.1, fy=0.9, log='xy')
      text(xy[1], xy[2], paste0(mirna, "\n", cond), adj=0)
      data_x <<- data_y
    })
  })
}

  # xy <- GetPlotFractionalCoords(0.9, 0.1, log='xy')
  # AddCorrelationToPlot(x=log(x), y=log(y), xpos=xy[1], ypos=xy[2], adj=1)
  # xmin <- 1
  # xmax <- 1e6
  # ymin <- xmin
  # ymax <- xmax
  # sapply(1:6, function(i) {
  #   BlankPlot(log='xy')
  # AddLogAxis(1, label="Old 12mer counts")
  # AddLogAxis(2, label="New 12mer counts")
  # x <- sXc[, i + 1]
  # y <- sXc_old[, i]
  # points(x, y, col=ConvertRColortoRGB("black", alpha=0.1))
  # xy <- GetPlotFractionalCoords(0.9, 0.1, log='xy')
  # AddCorrelationToPlot(x=log(x[(x > 0 & y > 0)]), y=log(y[(x > 0 & y > 0)]), xpos=xy[1], ypos=xy[2], adj=1)
  #   })
  # BlankPlot(log='xy')
  # AddLogAxis(1, label="Old 12mer counts")
  # AddLogAxis(2, label="New 12mer counts")
  # x <- sXc[, i]
  # y <- sXc_old[, i]
  # x <<- x
  # y <<- y
  # points(x, y, col=ConvertRColortoRGB("black", alpha=0.1))
  # xy <- GetPlotFractionalCoords(0.9, 0.1, log='xy')
  # AddCorrelationToPlot(x=log(x[(x > 0 & y > 0)]), y=log(y[(x > 0 & y > 0)]), xpos=xy[1], ypos=xy[2], adj=1)


AllMirnaHist12merKds <- function(experiment="equilibrium", n_constant=5,
                            sitelist="12mers") {
  n_left = n_constant
  n_right = n_constant
  mirna.start <- 2
  dev.new(xpos=20, ypos=20, height=8, width=12)
  par(mfcol=c(3, 4))
  sapply(kMirnas, function(mirna) {
    kd12mers_old <- SubfunctionCall(Get12merKdsOld)
    kd12mers_old <- kd12mers_old[1:(nrow(kd12mers_old) - 3), ]
    # kd12mers_new <- SubfunctionCall(KdVector)
    kd12mers_full <- SubfunctionCall(Get12merKdsFull)
    kd12mers_full <- SubfunctionCall(Get12merKdsFull)

    sXc <- SubfunctionCall(SitesXCounts)
    print(dim(sXc))
    sXc_old <- SubfunctionCall(SitesXCounts12mersOld)
    # print(head(sXc_old))
    # print(head(sXc))
    # print(head(kd12mers_old))
    print(head(kd12mers_new))
    print(head(kd12mers_full))

    shared <- intersect(names(kd12mers_old), names(kd12mers_new))
    print(length(shared))
    print(length(names(kd12mers_old)))
    print(length(names(kd12mers_new)))


    x <- Logistic(kd12mers_old[shared], 10)
    y <- kd12mers_new[shared]
    y2 <- kd12mers_full[shared, ]
    print(head(y2))
    xmin <- 1e-4
    xmax <- 10
    ymin <- xmin
    ymax <- xmax
    par(kPlotParameters)
    BlankPlot(log='xy')
    AddLogAxis(1, label="Old 12mer Kds")
    AddLogAxis(2, label="New 12mer Kds")
    points(x, y, col=ConvertRColortoRGB("black", alpha=0.1))
    BlankPlot(log='xy')
    AddLogAxis(1, label="Old 12mer Kds")
    AddLogAxis(2, label="New 12mer Kds")
    points(x, y2, col=ConvertRColortoRGB("black", alpha=0.1))
    plot(density(log(x)))
    plot(density(log(y2)), add=TRUE)

})
  # xy <- GetPlotFractionalCoords(0.9, 0.1, log='xy')
  # AddCorrelationToPlot(x=log(x), y=log(y), xpos=xy[1], ypos=xy[2], adj=1)
  # xmin <- 1
  # xmax <- 1e6
  # ymin <- xmin
  # ymax <- xmax
  # sapply(1:6, function(i) {
  #   BlankPlot(log='xy')
  # AddLogAxis(1, label="Old 12mer counts")
  # AddLogAxis(2, label="New 12mer counts")
  # x <- sXc[, i + 1]
  # y <- sXc_old[, i]
  # points(x, y, col=ConvertRColortoRGB("black", alpha=0.1))
  # xy <- GetPlotFractionalCoords(0.9, 0.1, log='xy')
  # AddCorrelationToPlot(x=log(x[(x > 0 & y > 0)]), y=log(y[(x > 0 & y > 0)]), xpos=xy[1], ypos=xy[2], adj=1)
  #   })
  # BlankPlot(log='xy')
  # AddLogAxis(1, label="Old 12mer counts")
  # AddLogAxis(2, label="New 12mer counts")
  # x <- sXc[, i]
  # y <- sXc_old[, i]
  # x <<- x
  # y <<- y
  # points(x, y, col=ConvertRColortoRGB("black", alpha=0.1))
  # xy <- GetPlotFractionalCoords(0.9, 0.1, log='xy')
  # AddCorrelationToPlot(x=log(x[(x > 0 & y > 0)]), y=log(y[(x > 0 & y > 0)]), xpos=xy[1], ypos=xy[2], adj=1)
}



SiteFlanksXCounts <- function(mirna, site, experiment="equilibrium", n_constant=5,
                         sitelist="paper") {
  sites_file_name <- paste0(kSolexaDir,
                       mirna, "/",experiment,"/site_count_tables/", site,
                       "_flanking_", n_constant, "_", sitelist, ".txt")
  sitesXcounts <- read.table(sites_file_name, stringsAsFactors=FALSE)
  colnames.temp <- colnames(sitesXcounts)[3 : ncol(sitesXcounts)]
  colnames.new <- sapply(colnames.temp, function(name) {
    return(unlist(strsplit(name, split="A", fixed=TRUE))[2])
  })
  colnames(sitesXcounts)[3 : ncol(sitesXcounts)] <- colnames.new
  return(sitesXcounts)
}


SitesXCountsUnique <- function(mirna, experiment, start, stop, sitelist,
                            mirna.start=NULL, mirna.stop=NULL) {
  if (sitelist %in% kKmerSiteLists) {
   sites_file_name <- paste0(kSolexaDir,
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
    sites_file_name <- paste0(kSolexaDir,
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

SitesXCountsKinetics <- function(mirna, experiment, n_constant, sitelist, n_lag=0) {
  if (n_lag > 0) {
    sites_file_name <- paste0(kSolexaDir, mirna, "/", experiment,
                              "/site_count_tables/all_sites_", n_constant, "_",
                              sitelist, "_ex", n_lag, "_pulse.txt")    
  } else {
    sites_file_name <- paste0(kSolexaDir, mirna, "/",experiment, 
                              "/site_count_tables/all_sites_", n_constant, "_",
                              sitelist, "_pulse.txt")

  }
  print(sites_file_name)
  sitesXcounts.pulse <- read.table(sites_file_name)
  # sitesXcounts.pulse <- sitesXcounts.pulse[,-4]
  # colnames(sitesXcounts.pulse)[3] <- "I_combined"
  colnames(sitesXcounts.pulse) <- sapply(colnames(sitesXcounts.pulse),
                                         RemoveStringPrefix, prefix="X")
  colnames(sitesXcounts.pulse)[(
  ncol(sitesXcounts.pulse)-2):ncol(sitesXcounts.pulse)] <- c(14400,"Equil-","0-")

  if (n_lag > 0) {
    sites_file_name <- paste0(kSolexaDir, mirna, "/", experiment,
                              "/site_count_tables/all_sites_", n_constant, "_",
                              sitelist, "_ex", n_lag, "_chase.txt")    
  } else {
    sites_file_name <- paste0(kSolexaDir, mirna, "/",experiment, 
                              "/site_count_tables/all_sites_", n_constant, "_",
                              sitelist, "_chase.txt")

  }
  sitesXcounts.chase <- read.table(sites_file_name)
  # sitesXcounts.chase <- sitesXcounts.chase[,-3]
  # colnames(sitesXcounts.chase)[3] <- "I_combined"
  colnames(sitesXcounts.chase) <- sapply(colnames(sitesXcounts.chase),
                                         RemoveStringPrefix, prefix="X")
  colnames(sitesXcounts.chase)[(
  ncol(sitesXcounts.chase)-2):ncol(sitesXcounts.chase)] <- c(14400,"Equil-","0-")

  return(list(sitesXcounts.pulse,sitesXcounts.chase))
}

SiteFlanksXCountsKinetics <- function(mirna, experiment, n_constant, sitelist, site) {
  sites_file_name_pulse <- paste0(kSolexaDir,
                       mirna, "/",experiment,"/site_count_tables/", site,
                       "_flanking_pulse_", n_constant,"_", sitelist, ".txt")
  sitesXcounts.pulse <- read.table(sites_file_name_pulse, stringsAsFactors=FALSE)
  sitesXcounts.pulse <- sitesXcounts.pulse[,-3]

  colnames(sitesXcounts.pulse)[2] <- "I_combined"
  colnames(sitesXcounts.pulse) <- sapply(colnames(sitesXcounts.pulse),
                                         RemoveStringPrefix, prefix="X")
  colnames(sitesXcounts.pulse)[(
  ncol(sitesXcounts.pulse)-2):ncol(sitesXcounts.pulse)] <- c(14400,"Equil-","0-")

  sites_file_name_chase <- paste0(kSolexaDir,
                       mirna, "/",experiment,"/site_count_tables/", site,
                       "_flanking_chase_", n_constant, "_", sitelist, ".txt")
  sitesXcounts.chase <- read.table(sites_file_name_chase, stringsAsFactors=FALSE)
  sitesXcounts.chase <- sitesXcounts.chase[,-2]
  colnames(sitesXcounts.chase)[2] <- "I_combined"
  colnames(sitesXcounts.chase) <- sapply(colnames(sitesXcounts.chase),
                                         RemoveStringPrefix, prefix="X")
  colnames(sitesXcounts.chase)[(
  ncol(sitesXcounts.chase)-2):ncol(sitesXcounts.chase)] <- c(14400,"Equil-","0-")

  return(list(sitesXcounts.pulse, sitesXcounts.chase))
}

GetSiteKmersXCountsKinetics <- function(mirna, experiment, site, start, stop,
                                 sitelist, k) {
  sites_file_name_pulse <- paste0(kSolexaDir,
                       mirna, "/",experiment,"/full_site_count_tables/", site,
                       "_sitekmers_", start, "-", stop,"_", sitelist, "_k", k, "_pulse.txt")
  sitesXcounts.pulse <- read.table(sites_file_name_pulse, stringsAsFactors=FALSE)
  sites_file_name_chase <- paste0(kSolexaDir,
                       mirna, "/",experiment,"/full_site_count_tables/", site,
                       "_sitekmers_", start, "-", stop,"_", sitelist, "_k", k, "_chase.txt")
  sitesXcounts.chase <- read.table(sites_file_name_chase, stringsAsFactors=FALSE)

  return(list(sitesXcounts.pulse, sitesXcounts.chase))
}




ConvertTtoU <- function(site) {
  return(gsub("T", c(site), replace="U"))
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




WriteIterationFile <- function(out,extension="") {
  out.file <- paste0(kSolexaDir, mirna,
                     "/equilibrium/kds/", site,"_flanking_",
                     k.c.stockago, extension,".txt")
  write.table(file=out.file, out, sep="\t", quote=FALSE, row.names=FALSE,
                col.names=TRUE)
}

WriteFinalParameterFile <- function(out,extension="") {
  out.final <- out[dim(out)[1], ]
  names(out.final) <- colnames(out)
  out.file <- paste0(kSolexaDir, mirna,
                           "/equilibrium/kds/final_", site,
                           "_flanking_", k.c.stockago,extension, ".txt")
  write.table(file=out.file, out.final, sep="\t", quote=FALSE, row.names=TRUE,
              col.names=FALSE)
}

## Helpful Fnctions





GetSiteKoffs <- function(mirna, experiment, n_constant, sitelist, costfunc, subset=FALSE, dil=FALSE) {
    if (dil !=FALSE) {
      dil.ext <- "_plusdil"
    } else {
      dil.ext <- ""
    }
    if (subset != FALSE) {
      subset.ext <- paste0("_", subset, "-sites")
    } else {
      subset.ext <- ""
    }
    params.file <- paste0(kSolexaDir, mirna, "/",
                 experiment, "/koffs_PAPER/", n_constant, "_", 
                 sitelist, "_", costfunc, subset.ext, "_round2", dil.ext, ".txt")
    print(params.file)
    params <- data.frame(read.table(params.file, header=TRUE, row.names=NULL,
                         stringsAsFactors=FALSE))
    colnames(params) <- sapply(colnames(params),
                                         RemoveStringPrefix, prefix="X")
    params.final <- as.numeric(params[nrow(params), -ncol(params)])
    names(params.final) <- colnames(params)[1:length(params.final)]
    # print(params.final)
    # kds <- params.final[grep("_Kd", names(params.final))]
    # koffs <- params.final[grep("_koff", names(params.final))]
    # kinetics_params <- cbind(kds, koffs)
    # contam_params <- params.final["contam"]
    return(params.final)
}

GetFlankKoffs <- function(mirna, experiment, n_constant, sitelist, site, costfunc, dil=FALSE) {
    if (dil !=FALSE) {
      dil.ext <- "_plusdil"
    } else {
      dil.ext <- ""
    }
    params.file <- paste0(kSolexaDir, mirna, "/",
                 experiment, "/koffs_PAPER/", n_constant, "_", 
                 sitelist, "_", costfunc, "_", site, "_round2", dil.ext, ".txt")
    params <- data.frame(read.table(params.file, header=TRUE, row.names=NULL,
                         stringsAsFactors=FALSE))
    colnames(params) <- sapply(colnames(params),
                                         RemoveStringPrefix, prefix="X")
    params.final <- as.numeric(params[nrow(params), -ncol(params)])
    names(params.final) <- colnames(params)[1:length(params.final)]
    return(params.final)
}

CompareFlankingKoffs <- function(mirna, experiment, n_constant, sitelist, site1, site2, costfunc, dil=FALSE) {
  flanks1 <- GetFlankKoffs(mirna, experiment, n_constant, sitelist, site1, costfunc, dil=FALSE)
  flanks2 <- GetFlankKoffs(mirna, experiment, n_constant, sitelist, site2, costfunc, dil=FALSE)
}

PlotKoffsVsKds <- function(mirna, n_constant, sitelist, ident=FALSE, costfunc="multinom", subset=FALSE, dil=FALSE, alpha=0.5) {
    dev.new(xpos = 20, ypos = 20, height = 8, width = 12)
    par(kPlotParameters)
    par(mfrow = c(2, 3))
    par(mgp = c(0.5, 0.6, 0.0))
    tempkds <- GetSiteKds(mirna, "equilibrium", n_constant, sitelist)
    tempkoffs <- GetSiteKoffs(mirna, "kinetics", n_constant, sitelist, costfunc, subset=subset, dil = dil)
    print(tempkoffs)
    if (subset != FALSE) {
      tempkds <- tempkds[c(1:subset, nrow(tempkds)-2), ]
    } else {
      tempkds <- tempkds[1:(nrow(tempkds)-2), ]
    }
    print(tempkoffs)
      kinetics.kds <- 10^tempkoffs[1:nrow(tempkds)]    
      kinetics.koff <- 10^tempkoffs[(nrow(tempkds) + 1):(2*nrow(tempkds))]
      print(kinetics.kds)
    tempkds <<- tempkds
    tempkoffs <<- tempkoffs
    kinetics.koff <<- kinetics.koff
    kinetics.kds <<- kinetics.kds
      # break

    if (sitelist %in% c("biological12mers", "biological12mersnew")){
      colors <- rep("black", length(kinetics.kds))
      ind.6mers <- grep(MirnaTargetSequence("miR-1", 2, 7), rownames(tempkds), fixed=TRUE)
      colors[ind.6mers] <- ConvertRColortoRGB(kSiteColors["6mer"], alpha = alpha)
      ind.7mersA1 <- grep(MirnaTargetSequence("miR-1", 1, 7), rownames(tempkds), fixed=TRUE)
      colors[ind.7mersA1] <- ConvertRColortoRGB(kSiteColors["7mer-A1", ], alpha = alpha)
      ind.7mersm8 <- grep(MirnaTargetSequence("miR-1", 2, 8), rownames(tempkds), fixed=TRUE)
      colors[ind.7mersm8] <- ConvertRColortoRGB(kSiteColors["7mer-m8"], alpha = alpha)
      ind.8mers <- grep(MirnaTargetSequence("miR-1", 1, 8), rownames(tempkds), fixed=TRUE)
      colors[ind.8mers] <- ConvertRColortoRGB(kSiteColors["8mer"], alpha = alpha)
    } else {
      colors <- kSiteColors[rownames(tempkds)]
    }

    kdl    <- seq(-4, 3)
    kds    <- sapply(kdl, function(i) {10^i * seq(9)})
    kdl <- 10^kdl
    kdmin <- kdl[1]
    kdmax <- kdl[length(kdl)]
    kdlims <- c(kdmin, kdmax)
    koffl    <- seq(-3, 2)
    koffs    <- sapply(koffl, function(i) {10^i * seq(9)})
    koffl    <- 10^koffl
    koffmin <- koffl[1]
    koffmax <- koffl[length(koffl)]
    kofflims <- c(koffmin, koffmax)

    plot(tempkds$Mean,
      kinetics.kds,
      log = 'xy',
      col=colors,
      xlim = kdlims,
      ylim = kdlims,
      ann=FALSE,
      axes=FALSE)
    mtext(costfunc, side = 3, line = 1)
    abline(0, 1, lty = 2, lwd = 0.5)
    # axis(1, at = ticks, labels = FALSE, pos = min(xlims))
    # axis(2, at = ticks, labels = FALSE, pos = min(xlims))
    # axis(1, at = 10^seq(-3, 3), lwd = 0, pos = min(xlims))
    # axis(2, at = 10^seq(-3, 3), lwd = 0, pos = min(xlims))
    axis(1, at=kdl,
         labels=sapply(kdl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=kdmin, lwd=0, line = 3, cex.axis = 1.5)
    axis(1, at=kds, labels=FALSE,
         pos=kdmin)
    # Label the axis at each order of magnitude.

    axis(2, at=kdl,
         labels=sapply(kdl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=kdmin, las=2, line = 1, lwd=0, cex.axis = 1.5)
    axis(2, at=kds, labels=FALSE,
         pos=kdmin)


    title(xlab = "Equilibrium KD (nM)", line= 1, cex.lab = 1.5)
    title(ylab = "Kinetics KD (nM)", line = 1.5, cex.lab = 1.5)



    plot(tempkds$Mean,
      kinetics.koff,
      log = 'xy',
      col=colors,
      xlim = kdlims,
      ylim = kofflims,
      ann = FALSE,
      axes = FALSE)
    xstart <- kdmin
    range <- kdmax/kdmin
    yoffset <- (kinetics.koff["None_koff"])/(tempkds$Mean[nrow(tempkds)])
    segments(xstart, xstart*yoffset, xstart*range, xstart*yoffset*range, lty = 2, lwd = 0.5)

    axis(1, at=kdl,
         labels=sapply(kdl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=koffmin, lwd=0, line = 3, cex.axis = 1.5)
    axis(1, at=kds, labels=FALSE,
         pos=koffmin)
    # Label the axis at each order of magnitude.

    axis(2, at=koffl,
         labels=sapply(koffl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=kdmin, las=2, line = 1, lwd=0, cex.axis = 1.5)
    axis(2, at=koffs, labels=FALSE,
         pos=kdmin)


    title(xlab = "Equilibrium KD (nM)", line= 1, cex.lab = 1.5)
    title(ylab = "Kinetics koff (min^-1)", line = 1.5, cex.lab = 1.5)





    if (ident == TRUE) {
      identify(tempkds$Mean[1:(nrow(tempkds)-2)],
      kinetics.koff,
      labels = rownames(tempkds)[1:(nrow(tempkds)-2)])
    }
    title(main = mirna)
    plot(kinetics.kds,
      kinetics.koff,
      log = 'xy',
      col=colors,
      xlim = kdlims,
      ylim = kofflims,
      ann = FALSE,
      axes = FALSE)
    xstart <- kdmin
    range <- kdmax/kdmin
    print(kinetics.koff)
    yoffset <- (kinetics.koff["None_koff"])/(kinetics.kds["None_Kd"])
    segments(xstart, xstart*yoffset, xstart*range, xstart*yoffset*range, lty = 2, lwd = 0.5)

    axis(1, at=kdl,
         labels=sapply(kdl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=koffmin, lwd=0, line = 3, cex.axis = 1.5)
    axis(1, at=kds, labels=FALSE,
         pos=koffmin)
    # Label the axis at each order of magnitude.

    axis(2, at=koffl,
         labels=sapply(koffl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=kdmin, las=2, line = 1, lwd=0, cex.axis = 1.5)
    axis(2, at=koffs, labels=FALSE,
         pos=kdmin)


    title(xlab = "Kinetics KD (nM)", line= 1, cex.lab = 1.5)
    title(ylab = "Kinetics koff (min^-1)", line = 1.5, cex.lab = 1.5)






    if (ident == TRUE) {
      identify(tempkds$Mean[1:(nrow(tempkds)-2)],
      kinetics.koff,
      labels = rownames(tempkds)[1:(nrow(tempkds)-2)])
    }
    costfunc <- "multinom"
    tempkoffs <- GetSiteKoffs(mirna, "kinetics", n_constant, sitelist, costfunc, subset=subset, dil = dil)
    print(tempkoffs)
      kinetics.kds <- 10^tempkoffs[1:nrow(tempkds)]
      print(kinetics.kds)
      break
      kinetics.koff <- 10^tempkoffs[(nrow(tempkds) + 1):(2*nrow(tempkds))]
    tempkds <<- tempkds
    tempkoffs <<- tempkoffs
    kinetics.koff <<- kinetics.koff
    kinetics.kds <<- kinetics.kds
    if (sitelist %in% c("biological12mers", "biological12mersnew")){
      colors <- rep("black", length(kinetics.kds))
      ind.6mers <- grep(MirnaTargetSequence("miR-1", 2, 7), rownames(tempkds), fixed=TRUE)
      colors[ind.6mers] <- ConvertRColortoRGB(kSiteColors["6mer" ], alpha = alpha)
      ind.7mersA1 <- grep(MirnaTargetSequence("miR-1", 1, 7), rownames(tempkds), fixed=TRUE)
      colors[ind.7mersA1] <- ConvertRColortoRGB(kSiteColors["7mer-A1" ], alpha = alpha)
      ind.7mersm8 <- grep(MirnaTargetSequence("miR-1", 2, 8), rownames(tempkds), fixed=TRUE)
      colors[ind.7mersm8] <- ConvertRColortoRGB(kSiteColors["7mer-m8" ], alpha = alpha)
      ind.8mers <- grep(MirnaTargetSequence("miR-1", 1, 8), rownames(tempkds), fixed=TRUE)
      colors[ind.8mers] <- ConvertRColortoRGB(kSiteColors["8mer" ], alpha = alpha)
    } else {
      colors <- kSiteColors[rownames(tempkds)]
    }

    plot(tempkds$Mean,
      kinetics.kds,
      log = 'xy',
      col=colors,
      xlim = kdlims,
      ylim = kdlims,
      ann=FALSE,
      axes=FALSE)
    mtext(costfunc, side = 3, line = 1)
    abline(0, 1, lty = 2, lwd = 0.5)
    # axis(1, at = ticks, labels = FALSE, pos = min(xlims))
    # axis(2, at = ticks, labels = FALSE, pos = min(xlims))
    # axis(1, at = 10^seq(-3, 3), lwd = 0, pos = min(xlims))
    # axis(2, at = 10^seq(-3, 3), lwd = 0, pos = min(xlims))
    axis(1, at=kdl,
         labels=sapply(kdl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=kdmin, lwd=0, line = 3, cex.axis = 1.5)
    axis(1, at=kds, labels=FALSE,
         pos=kdmin)
    # Label the axis at each order of magnitude.

    axis(2, at=kdl,
         labels=sapply(kdl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=kdmin, las=2, line = 1, lwd=0, cex.axis = 1.5)
    axis(2, at=kds, labels=FALSE,
         pos=kdmin)


    title(xlab = "Equilibrium KD (nM)", line= 1, cex.lab = 1.5)
    title(ylab = "Kinetics KD (nM)", line = 1.5, cex.lab = 1.5)



    plot(tempkds$Mean,
      kinetics.koff,
      log = 'xy',
      col=colors,
      xlim = kdlims,
      ylim = kofflims,
      ann = FALSE,
      axes = FALSE)
    xstart <- kdmin
    range <- kdmax/kdmin
    yoffset <- (kinetics.koff["None_koff"])/(tempkds$Mean[nrow(tempkds)])
    segments(xstart, xstart*yoffset, xstart*range, xstart*yoffset*range, lty = 2, lwd = 0.5)

    axis(1, at=kdl,
         labels=sapply(kdl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=koffmin, lwd=0, line = 3, cex.axis = 1.5)
    axis(1, at=kds, labels=FALSE,
         pos=koffmin)
    # Label the axis at each order of magnitude.

    axis(2, at=koffl,
         labels=sapply(koffl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=kdmin, las=2, line = 1, lwd=0, cex.axis = 1.5)
    axis(2, at=koffs, labels=FALSE,
         pos=kdmin)


    title(xlab = "Equilibrium KD (nM)", line= 1, cex.lab = 1.5)
    title(ylab = "Kinetics koff (min^-1)", line = 1.5, cex.lab = 1.5)





    if (ident == TRUE) {
      identify(tempkds$Mean[1:(nrow(tempkds)-2)],
      kinetics.koff,
      labels = rownames(tempkds)[1:(nrow(tempkds)-2)])
    }
    title(main = mirna)
    plot(kinetics.kds,
      kinetics.koff,
      log = 'xy',
      col=colors,
      xlim = kdlims,
      ylim = kofflims,
      ann = FALSE,
      axes = FALSE)
    xstart <- kdmin
    range <- kdmax/kdmin
    print(kinetics.koff)
    yoffset <- (kinetics.koff["None_koff"])/(kinetics.kds["None_Kd"])
    segments(xstart, xstart*yoffset, xstart*range, xstart*yoffset*range, lty = 2, lwd = 0.5)

    axis(1, at=kdl,
         labels=sapply(kdl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=koffmin, lwd=0, line = 3, cex.axis = 1.5)
    axis(1, at=kds, labels=FALSE,
         pos=koffmin)
    # Label the axis at each order of magnitude.

    axis(2, at=koffl,
         labels=sapply(koffl, function(name) {
           eval(substitute(expression(10^x), list(x=log10(name))))
         }),
         pos=kdmin, las=2, line = 1, lwd=0, cex.axis = 1.5)
    axis(2, at=koffs, labels=FALSE,
         pos=kdmin)


    title(xlab = "Kinetics KD (nM)", line= 1, cex.lab = 1.5)
    title(ylab = "Kinetics koff (min^-1)", line = 1.5, cex.lab = 1.5)






    if (ident == TRUE) {
      identify(tempkds$Mean[1:(nrow(tempkds)-2)],
      kinetics.koff,
      labels = rownames(tempkds)[1:(nrow(tempkds)-2)])
    }

}


# 2F____________________________________________________________________________
PlotRelKdsVsDDG <- function(experiment, n_constant, sitelist, xpos=20, ypos=20,
                            new=FALSE, height=5, width=5, pdf.plot=FALSE,
                            dG.table=3, square=FALSE){
  dG.df <- read.table(file=paste0("canonical_sites_mfe_", dG.table, ".txt"))
  colnames(dG.df) <- gsub("(^.*)\\.(.*$)", colnames(dG.df), replace="\\1-\\2")   
  R <- 1.987e-3 # in kcal K-1 mol-1
  T <- 310.15 # in K
  total_x <- c()
  total_y <- c()
  if (class(pdf.plot) == "character") {
    pdf(file=paste0("2017_Paper/", pdf.plot, ".pdf"), height=height*2/5, width=width*2/5)
    par(kPDFParameters)
  } else {
    dev.new(xpos=xpos, ypos=ypos, height=height, width=width)
    par(kPlotParameters)    
  }
  xmin <- 0.005
  xmax <- 2
  ymin <- -4
  ymax <- 0.5
  if (square) {
    adjusted = TRUE
  } else {
    adjusted = FALSE
  }
  BlankPlot(log="x", adjusted=adjusted, inv='xy')
  AddLogAxis(1, label="Kd relative to 6mer")
  AddLinearAxis(2, tick.space=0.5, label=expression(Delta*Delta*"G; predicted"))
  x.line <- exp(seq(-8, 6, length=20))
  y.line <- R*T*log(x.line)
  lines(x.line, y.line, col="gray60", lwd=2)
  for (mirna in kMirnas) {
    color = kMirnaColors[mirna]
    kds   <- GetSiteKds(mirna, "equilibrium", 5, "paper")[kCanonicalSites,]$Mean
    kds.0 <- GetSiteKds(mirna, "equilibrium", 5, "paper")["6mer",]$Mean
    dG    <- dG.df[kCanonicalSites, mirna]
    dG.0  <- dG.df["6mer", mirna]
    x.data <- kds/kds.0
    y.data  <- dG - dG.0
    total_x <- c(total_x, x.data)
    total_y <- c(total_y, y.data)
    lmodel <- lm(y.data ~ log(x.data))
    m <- lmodel$coefficients[2]
    b <- lmodel$coefficients[1]
    y.line <- m*log(x.line) + b
    lines(x.line, y.line, col=color, lty=2)
    points(x.data, y.data, col=c(rep(color, 3), "gray"), cex=1.5, lwd=1.2, pch = c(1, 2, 0, 20)) 
  }
  text(3, -3.5, round(cor(total_x, total_y)^2, 3))
  leg.xy1 <- GetPlotFractionalCoords(fx=0, fy=1, log='x', inv='xy')
  legend(x=leg.xy1[1], y=leg.xy1[2], seg.len=1, legend=kMirnas, col=kMirnaColors[kMirnas], bty="n", lty=1)
  leg.xy2 <- GetPlotFractionalCoords(fx=0.3, fy=1, log='x', inv='xy')
  legend(x=leg.xy2[1], y=leg.xy2[2], legend = kCanonicalSites, bty="n",pch=c(1, 2, 0, 20))
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

#23456789!123456789@123456789#123456789$123456789%123456789^123456789&123456789*
# Figure 2H
Plot8merVs7merEnergyCombined <- function(experiment, n_constant, sitelist, xpos=20,
                                   ypos=20, height=5, width=5, pdf.plot=FALSE,
                                   dG.table=3){
  dG.df <- read.table(file=paste0("canonical_sites_mfe_", dG.table, ".txt"))
  colnames(dG.df) <- gsub("(^.*)\\.(.*$)", colnames(dG.df), replace="\\1-\\2")   
  R <- 1.987e-3 # in kcal K-1 mol-1
  T <- 310.15 # in K
  total_x <- c(0)
  total_y <- c(0)
  if (class(pdf.plot) == "character") {
    pdf(file=paste0("2017_Paper/", pdf.plot, ".pdf"), height=height*2/5, width=width*2/5)
    par(kPDFParameters)
  } else {
    dev.new(xpos=xpos, ypos=ypos, height=height, width=width)
    par(kPlotParameters)    
  }
  xmin <- -4
  xmax <- 0.5
  ymin <- -4
  ymax <- 0.5
  BlankPlot(inv='xy')
  AddLinearAxis(1, tick.space=0.5, label=expression("7mer-m8"~Delta*Delta*"G + 7mer-A1 "~Delta*Delta*"G (relative to 6mer)"))
  AddLinearAxis(2, tick.space=0.5, label=expression("8mer"~Delta*Delta*"G (relative to 6mer)"))

  abline(0, 1, lty=2, col="gray")
  for (mirna in kMirnas) {
    dG    <- dG.df[kSeedSites, mirna]
    names(dG) <- rownames(dG.df[kSeedSites, ])
    ddG    <-  dG - dG["6mer"]
    ddG.8   <- ddG["8mer"]
    ddG.7m8 <- ddG["7mer-m8"]
    ddG.7a1 <- ddG["7mer-A1"]
    total_x <- c(total_x, ddG.7m8 + ddG.7a1)
    total_y <- c(total_y, ddG.8)
    arrows(0, ddG.8, ddG.7a1, ddG.8, length=0.05*par()$cex, angle=90,
         code=3, col=kMirnaColors[mirna], lty=1)
    arrows(ddG.7a1, ddG.8, ddG.7a1 + ddG.7m8, ddG.8, length=0.05*par()$cex, angle=90,
         code=3, lty=2, col=kMirnaColors[mirna])
  }
  points(total_x, total_y, col=c("gray60",kMirnaColors[kMirnas]), pch=19)    
  cor.xy <- GetPlotFractionalCoords(fx=0.9, fy=0.9, inv='xy')
  AddCorrelationToPlot(total_x, total_y, xpos=cor.xy[1], ypos=cor.xy[2], adj=1, rsquared=TRUE)
  leg.xy <- GetPlotFractionalCoords(fx=0.7, fy=0.35, inv='xy')
  legend(x=leg.xy[1], y=leg.xy[2], legend = kMirnas, col=kMirnaColors[kMirnas], bty="n",pch=19)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}


GetFlankKds <- function(mirna, site, experiment="equilibrium", n_constant=5,
                        sitelist="paper") {
  extension <- paste0("_", paste(sitelist, site, "PAPER", sep="_"))
  path <- GetAnalysisPath(mirna, experiment, n_constant, "kds_PAPER",
                        ext=extension)
  data <- fread(path, fill=TRUE,header=TRUE,
                       stringsAsFactors=FALSE, showProgress=FALSE)
  colnames(data) <- c("", colnames(data)[-ncol(data)])
  out <- data.frame(data, row.names=1)
  out
}


GetFlankbgKds <- function(mirna, experiment, n_constant, sitelist, site) {
    params.file <- paste0(kSolexaDir, mirna, "/",
                 experiment, "/kds_PAPER/", n_constant, "_", 
                 sitelist, "_", site, "_controlplfold_PAPER.txt")
    params <- data.frame(read.table(params.file, header=TRUE, row.names=1,
                         stringsAsFactors=FALSE))
    return(params)
}

GetFlankbyPlFoldKds <- function(mirna, experiment, n_constant, sitelist, site) {
    params.file <- paste0(kSolexaDir, mirna, "/",
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

SitesAndSingleFlanksXCounts <- function(mirna, site, experiment="equilibrium",
                                        n_constant=5, sitelist="paper") {
  print(mirna)
  sXc <- SubfunctionCall(SitesXCounts)
  print(sXc)
  fXc <- SubfunctionCall(SiteFlanksXCounts)
  site.ind <- which(rownames(sXc) == site)
  left.flanks <- substr(rownames(fXc), 1, 2)
  right.flanks <- substr(rownames(fXc), 4, 5)
  rownames(fXc) <- paste0(site, "|", rownames(fXc))
  rbind(sXc[0:(site.ind - 1), ],
        fXc,
        sXc[(site.ind + 1): nrow(sXc), ])
}

SitesAndSingleFlanksKds <- function(mirna, experiment, n_constant,
                                              sitelist, site) {
  pars.sites <- SubfunctionCall(GetSiteKds)
  pars.flanks <- SubfunctionCall(GetFlankKds)
  rownames(pars.flanks) <- paste0(site, "|", rownames(pars.flanks))
  site.ind <- which(rownames(pars.sites) == paste0(site, "_Kd"))
  rbind(pars.sites[0:(site.ind - 1), ],
        pars.flanks,
        pars.sites[(site.ind + 1): nrow(pars.sites), ])
}

#KEEP, for figure 3A
CombindSiteAndFlankColors <- function(sXc) {
  colors <- rep("gray", nrow(sXc))
  flanks.inds <- grep("\\|", rownames(sXc), perl=TRUE)
  flanks     <- gsub("(^.*)\\|(.*$)", rownames(sXc)[flanks.inds], replace="\\2",
                     perl=TRUE)
  colors[flanks.inds] <- GetColorFunction(flanks)
  colors
}

GetPlFoldDataNew <- function(mirna, experiment, condition, n_constant, sitelist,
                            site, win) {
  extension <- paste0("_", n_constant, "_", sitelist, "_", site,  "_win", win-1)
  path <- GetAnalysisPath(mirna, experiment, condition, "plfold_2018_PAPER",
                          ext=extension)
  print(path)
  data <- fread(file=path, header=TRUE, stringsAsFactors=FALSE, sep="\t")
  return(data)
}

MakeFlankAveragesNewSingle <- function(mirna, site, win, experiment="equilibrium",
                                    condition="I_combined", n_constant=5,
                                    sitelist="paper") {
  pl_data <- SubfunctionCall(GetPlFoldDataNew)
  flanks <- pl_data[, 2]
  pl_data <- t(t(pl_data[,-c(1, 2)]^(1/win)))
  out <- aggregate(pl_data,by=flanks,FUN=GeoMean)
  flanks.new <- out[,1]
  out <- out[,-1]
  rownames(out) <- flanks.new
  out <- out[order(flanks.new),]
  return(out)
}

MakeFlankAveragesNew <- function(mirna, site, experiment="equilibrium",
                                    condition="I_combined", n_constant=5,
                                    sitelist="paper") {
  pl_data_matrix <- sapply(seq(2), function(win_i) {
    SubfunctionCall(MakeFlankAveragesNewSingle, win=win_i)
  })
  return(pl_data_matrix)
}

GetCorrelationMatrixNew <- function(mirna, site, experiment="equilibrium",
                                    condition="I_combined", n_constant=5,
                                    sitelist="paper") {

  kds <- SubfunctionCall(GetFlankKds)
  kd_names <- gsub("^(.*)_Kd", rownames(kds), replace="\\1")
  kds <- kds$Mean
  names(kds) <- kd_names
  cormatrix <- sapply(seq(30), function(win_i) {
    data <- SubfunctionCall(MakeFlankAveragesNewSingle, win=win_i)
    out <- apply(log(data), 2, cor,
           y=-log(kds))
    return(out)
  })
  t(cormatrix)
}


MakeAveragePlFoldData <- function(pl_data, log=FALSE, per.nt = FALSE) {
  flanks <- pl_data[, 1]
  pl_data <- pl_data[,-1]
  if (per.nt == TRUE) {
    wins = as.numeric(gsub("^.*w(.*)$", colnames(pl_data), replace="\\1", perl=TRUE))
    pl_data = t(t(pl_data)^(1/wins))
  }
  if (log == TRUE) {
    pl_data <- log10(pl_data)
  }
  pl_data_summary <- colMeans(pl_data, na.rm=TRUE)
  pos.strings <- unique(gsub("(^.*)\\|(.*$)", colnames(pl_data), replace="\\1", perl=TRUE))
  win.strings <- unique(gsub("(^.*)\\|(.*$)", colnames(pl_data), replace="\\2", perl=TRUE))
  pl.data.matrix <- matrix(pl_data_summary, nrow=length(win.strings), ncol=length(pos.strings), byrow=FALSE, dimnames=list(win.strings, pos.strings))
  return(pl.data.matrix)
}



PlotContourPlFold <- function(pl_average_A, pl_average_I, xpos=20, ypos=20, pdf.plot=FALSE) {
  height <- 5
  width <- 10
  if (class(pdf.plot) == "character") {
    pdf(file=paste0("2017_Paper/", pdf.plot, ".pdf"), height=height*2/5, width=width*2/5)
    par(kPDFParameters)
  } else {
    dev.new(xpos=xpos, ypos=ypos, height=height, width=width)
    par(kPlotParameters)    
  }
  pl_average <- pl_average_A - pl_average_I
  xmin <- 0
  xmax <- ncol(pl_average) + 1
  ymin <- 0
  ymax <- 1
  BlankPlot()

  axis(1, at=seq(ncol(pl_average)), labels=colnames(pl_average), xlab="Position")
  AddLinearAxis(2, label="log(target site accessibility")


  sapply(seq(20), function(row) {
    print(pl_average[row,])
    lines(1:ncol(pl_average), pl_average[row, ], col = rainbow(20)[row])
    })
}



CheckAllPlFoldWindows <- function(mirna, experiment, condition, n_constant,
                                  sitelist, site, win_start, win_stop, win,
                                  logdata=FALSE,
                                  xpos=20, ypos=20, pdf.plot=FALSE) {
  height <- 5
  width <- 5
  Data.A <- MakeAveragePlFoldData(GetPlFoldData(mirna, experiment, condition, n_constant, sitelist,
                           site, win_start, win_stop, win, logdata=logdata))
  Data.I <- MakeAveragePlFoldData(GetPlFoldData(mirna, experiment, "I", n_constant, sitelist,
                           site, win_start, win_stop, win, logdata=logdata))

  # print(length(win.strings))
  # print(length(pos.strings))
  # print(pl.data.matrix[1:3, 1:3])
  # print(pl_data_summary[1:9])
  # print(pl_data_summary[31:39])
  dev.new(xpos = xpos, ypos=ypos, height=height, width=width)
  par(kPlotParameters)
  image(t(Data.A - Data.I))
}

SubsetPLFold <- function(pl_data_raw, start, stop, win, flank) {
  grep.symbol <- paste0("\\.w", win, "$")
  col_greps <- grep(grep.symbol, colnames(pl_data_raw), perl=TRUE, value=TRUE)
  row_greps <- grep("AA.AA", pl_data_raw[,1])
  return(pl_data_raw[row_greps,col_greps])
}


# Figure 3B
GetFullMirnaSiteFlankKds <- function(mirna, site, experiment="equilibrium",
                                     n_constant=5, sitelist="paper") {
  out <- setNames(rep(NaN, length(kFlanks)), kFlanks)
  f_kds <- SubfunctionCall(GetFlankKds)
  flanks <- rownames(f_kds)
  out[flanks] <- f_kds[flanks, ]$Mean
  return(out)
}

# Useful; keep:
PlotFlankKdScatter <-function(mirna_1, site_1, mirna_2=mirna_1, site_2=site_1,
                              experiment_1="equilibrium",
                              experiment_2=experiment_1, n_constant_1=5,
                              n_constant_2=nconstant_1, sitelist_1="paper",
                              sitelist_2=sitelist_1, identif=FALSE) {
  f_kds_1 <- SubfunctionCall(GetFullMirnaSiteFlankKds, mirna=mirna_1,
                              site=site_1, experiment=experiment_1,
                              sitelist=sitelist_1, n_constant=n_constant_1)
  f_kds_2 <- SubfunctionCall(GetFullMirnaSiteFlankKds, mirna=mirna_2,
                              site=site_2, experiment=experiment_2,
                              sitelist=sitelist_2, n_constant=n_constant_2)
  dev.new(xpos=20, ypos=20, height=5, width=5)
  par(kPlotParameters)
  xmin <- 1e-4
  ymin <- xmin
  xmax <- 10
  ymax <- xmax
  BlankPlot(log='xy', inv='xy')
  AddLogAxis(1, label=paste0(site_1, " flank Kd"))
  AddLogAxis(2, label=paste0(site_2, " flank Kd"))
  points(x=f_kds_1, y=f_kds_2, col=GetColorFunction(kFlanks))
  if (identif) {
    identify(x=f_kds_1, y=f_kds_2, labels=kFlanks)
  }
}

# Function important for 3C:
GetCanonicalSiteFlankingKdMatrix <- function(experiment, n_constant, sitelist) {
  mirsXsites <- data.frame(miRNA=rep(kMirnas, each=length(kSeedSites)),
                      site=rep(kSeedSites, length(kMirnas)))
  flank.data <- data.frame(t(apply(mirsXsites, 1, function(row) {
    GetFullMirnaSiteFlankKds(row[1], experiment, n_constant, sitelist, row[2])
  })))
  data.matrix <- cbind(mirsXsites, flank.data)
  colnames(data.matrix) <- c("miRNA", "site", kFlanks)
  return(data.matrix)
}

RelativeToNoneKds <- function(mirna) {
  kds <- GetSiteKds(mirna, sitelist="baek")
  names_kds <- rownames(kds)
  output <- kds["None_Kd",]$Mean/kds$Mean 
  names(output) <- gsub("^(.*)_Kd$", names_kds, replace="\\1", perl=TRUE)
  return(output)
}

Ratioof7merKds <- function(mirna) {
  kds <- GetSiteKds(mirna)
  names_kds <- rownames(kds)
  output <- kds["7mer-A1_Kd",]$Mean/kds["7mer-m8_Kd", ]$Mean 
  return(output)
}

GetSPS <- function(mirna, sitetype, dG.table=3) {
  dG.df <- read.table(file=paste0("canonical_sites_mfe_", dG.table, ".txt"))
  colnames(dG.df) <- gsub("(^.*)\\.(.*$)", colnames(dG.df), replace="\\1-\\2")
  R <- 1.987e-3 # in kcal K-1 mol-1
  T <- 310.15 # in K
  as.numeric(dG.df[sitetype, mirna])

}

GetSingleKd <- function(mirna, sitetype, experiment="equilibrium", n_constant=5,
                        sitelist="paper") {
  kds <- SubfunctionCall(GetSiteKds)
  out <- kds[paste0(sitetype, "_Kd"), ]$Mean/kds[paste0("None_Kd"), ]$Mean
}

TestSPSRelationship <- function(sitetype, experiment="equilibrium", n_constant=5,
                                sitelist="paper", dG.table=3) {
  tempSPS <- sapply(kMirnas, function(mirna) {
    SubfunctionCall(GetSPS)
  })
  tempKDS <- log(sapply(kMirnas, function(mirna) {
    SubfunctionCall(GetSingleKd)
  }))
  print(cor.test(tempSPS, tempKDS))
}

GetFlankRange <- function(mirna, site, experiment="equilibrium", n_constant=5,
                          sitelist="paper") {
  kds <- SubfunctionCall(GetFlankKds)
  flank_names <- rownames(kds)
  kds <- kds$Mean
  names(kds) <- flank_names
  print(kds[which(kds == min(kds))])
  print(kds[which(kds == max(kds))])
}

GetFlankSpan <- function(mirna, site, experiment="equilibrium", n_constant=5,
                          sitelist="paper") {
  kds <- SubfunctionCall(GetFlankKds)
  flank_names <- rownames(kds)
  kds <- kds$Mean
  max(kds)/min(kds)
}

GetFlankSD <- function(mirna, site, experiment="equilibrium", n_constant=5,
                          sitelist="paper") {
  kds <- SubfunctionCall(GetFlankKds)
  flank_names <- rownames(kds)
  kds <- kds$Mean
  exp(sd(log(kds)))
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


PlotPositionalKdsMiR7 <- function(n_constant, sitelist, xpos=20, ypos=20,
                                  pdf.plot=FALSE) {
  # Get kds for all site-types of the mirna.

  dev.new(xpos = xpos, ypos = ypos, height = 5, width = 8)
  kds <- GetSiteKds("miR-7-22nt", "equilibrium_nb", n_constant, sitelist)    
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
  kPositionalSites <- rownames(kds)[1:(nrow(kds) - 2)]
  kMirnaColors <- topo.colors(8)[2:5]
  names(kMirnaColors) <- c("miR-7-22nt","miR-7-23nt", "miR-7-24nt", "miR-7-25nt")
  kMirnas <- names(kMirnaColors)
  for (mirna in c("miR-7-22nt","miR-7-23nt", "miR-7-24nt", "miR-7-25nt")) {
    kds <- GetSiteKds(mirna, "equilibrium_nb", n_constant, sitelist)    
  kds <- kds[kPositionalSites,]
  print(kds)
  points(kds$Mean[ind_p], y_p,
        col = kMirnaColors[mirna],
        pch=19)
  lines(kds$Mean[ind_l], y_l,
        col = kMirnaColors[mirna],
        lwd = 2,
        type = "o")
  lines(kds$Lower_CI[ind_l], y_l,
        col = kMirnaColors[mirna],
        lwd = 1,
        lty=2,
        cex = 1)
  lines(kds$Upper_CI[ind_l], y_l,
        col = kMirnaColors[mirna],
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
         legend = kMirnas,
         col = kMirnaColors,
         bty="n",
         pch = 19)
}

ExtractSingleFlanks <- function(flank.names, side="5p") {
  if (side == "5p") {
    replace.string <- "\\1"
  } else if (side == "3p") {
    replace.string <- "\\2"
  }
  return(gsub("^(.{2,2}).*(.{2,2})$", flank.names, replace=replace.string,
         perl=TRUE))
}

GetMirnaFlankDataFrame <- function(mirna, experiment="equilibrium",
                                   n_constant=5, sitelist="paper") {
  rbindlist(lapply(kSeedSites, function(site) {
      # Get each site and miRNA Kd table:
      kds <- SubfunctionCall(GetFlankKds)
      flank.5p <- ExtractSingleFlanks(rownames(kds))
      flank.3p <- ExtractSingleFlanks(rownames(kds), side="3p")
      data.frame(logkd=log(kds$Mean),
                 mirna=rep(mirna, nrow(kds)),
                 site=rep(site, nrow(kds)),
                 f5p1=sapply(flank.5p, substr, start=1, stop=1),
                 f5p2=sapply(flank.5p, substr, start=2, stop=2),
                 f3p1=sapply(flank.3p, substr, start=1, stop=1),
                 f3p2=sapply(flank.3p, substr, start=2, stop=2))
  }))
}
# All part of linear leave one out model in figure 3:
################################################################################
GetSingleMirnaModel <- function(mirna, experiment="equilibrium",
                                   n_constant=5, sitelist="paper") {
  lm(logkd ~ site, data=SubfunctionCall(GetMirnaFlankDataFrame))
}

GetFlankLinearModel <- function(experiment="equilibrium", n_constant=5,
                                       sitelist="paper", leaveout=FALSE) {
  # Output matrix:
  if (leaveout != FALSE) {
    kMirnas <- setdiff(kMirnas, leaveout)
  }
  kd.data <- rbindlist(lapply(kMirnas, function(mirna) {
    rbindlist(lapply(kSeedSites, function(site) {
      # Get each site and miRNA Kd table:
      kds <- SubfunctionCall(GetFlankKds)
      flank.5p <- ExtractSingleFlanks(rownames(kds))
      flank.3p <- ExtractSingleFlanks(rownames(kds), side="3p")
      data.frame(logkd=log(kds$Mean),
                 mirna=rep(mirna, nrow(kds)),
                 site=rep(site, nrow(kds)),
                 f5p1=substr(flank.5p, 1, 1),
                 f5p2=substr(flank.5p, 2, 2),
                 f3p1=substr(flank.3p, 1, 1),
                 f3p2=substr(flank.3p, 2, 2))
    }))
  }))
  out <- lm(logkd ~ mirna*site + f5p1*f5p2 + f3p1*f3p2, data=kd.data)
  # print(out)
  out
}

ParseFlankModelCoefs <- function(model) {
  coefs <- c(f5p1A=0, f5p2A=0, f3p1A=0, f3p2A=0, model$coefficients)
  lin <- sapply(c("f5p1", "f5p2", "f3p1", "f3p2"), function(string) {
    grep.key <- paste0("^", string, ".{1}$")
    out <- coefs[grep(grep.key, names(coefs), perl=TRUE)]
    names(out) <- sapply(names(out), substr, start=5, stop=5)
    out
  })
  ints <- lapply(c(5, 3), function(flank){
    sapply(c("1C", "1G", "1T"), function(string) {
    grep.key <- paste0("^f", flank, "p", string, ":")
    out <- coefs[grep(grep.key, names(coefs), perl=TRUE)]
    len <- nchar(names(out))
    names(out) <- sapply(names(out), substr, start=len - 1, stop=len)
    out
    })  
  })
  list(lin=lin, int.5p=ints[[1]], int.3p=ints[[2]])
}

GetFlankLMCoefs <- function(experiment="equilibrium", n_constant=5,
                            sitelist="paper") {
  model <- SubfunctionCall(GetFlankLinearModel)
  SubfunctionCall(ParseFlankModelCoefs)
}

LeaveOneOutFlankModel <- function(mirna, experiment="equilibrium",
                                   n_constant=5, sitelist="paper") {
  flank.model <- SubfunctionCall(GetFlankLinearModel, leaveout=mirna)
  site.mir.model <- SubfunctionCall(GetSingleMirnaModel)
  # Get the site coefficients:
  site.coeffs <- c(site8mer=0, site.mir.model$coefficients)
  names(site.coeffs) <- gsub("^site(.*)$", names(site.coeffs),
                             replace="\\1", perl=TRUE)
  # Get the flank coefficients:
  flank.coeffs <- SubfunctionCall(ParseFlankModelCoefs, model=flank.model)
  flank.lin <- flank.coeffs[[1]]
  flank.5int <- flank.coeffs[[2]]
  flank.3int <- flank.coeffs[[3]]
  data.flanks <- SubfunctionCall(GetMirnaFlankDataFrame)
  test.model <- apply(data.flanks, 1, function(row) {
    site <- paste0("site", row["site"])
    f5p1 <- row["f5p1"]
    f5p2 <- row["f5p2"]
    f3p1 <- row["f3p1"]
    f3p2 <- row["f3p2"]
    site.val <- site.coeffs[row["site"]]
    f5p1.val <- flank.lin[f5p1, "f5p1"]
    f5p2.val <- flank.lin[f5p2, "f5p2"]
    f3p1.val <- flank.lin[f3p1, "f3p1"]
    f3p2.val <- flank.lin[f3p2, "f3p2"]
    if (f5p1 != "A" & f5p2 != "A") {
      int5.val <- flank.5int[paste0(2, f5p2), paste0(1, f5p1)]    
    } else {
      int5.val <- 0
    }
    if (f3p1 != "A" & f3p2 != "A") {
      int3.val <- flank.3int[paste0(2, f3p2), paste0(1, f3p1)]    
    } else {
      int3.val <- 0
    }
    sum(site.val, f5p1.val, f5p2.val, f3p1.val, f3p2.val, int5.val, int3.val)
  })
  test.model - mean(test.model, na.rm=TRUE) + mean(data.flanks$logkd, na.rm=TRUE)
}


PlotFlankCoefLinModel <- function(experiment="equilibrium", n_constant=5,
                                        sitelist="paper", xpos=20, ypos=20,
                                        pdf.plot=FALSE) {
  height <- 5
  width <- 5

  args <- as.list(match.call())
  flank.kds.all <- SubfunctionCall(GetFlankLinearModel, args)
  # Data analysis for the 1st plot:
  # First linear model; splitting up all flanks by 5p and 3p sequence.
  lm.split <- lm(logkd ~ mirna*site + f5p + f3p, data=flank.kds.all)

  # Data analysis for 2nd plot:
  # Get the coefficients for each of the 16 flanking dinucleotides to both
  # the 5p and 3p ends, to be able to plot them against one another to assess
  # their relative magnitude:
  model.coeffs <- c(lm.split$coefficients, f5pAA=0, f3pAA=0)
  f5p.coeffs <- model.coeffs[grep("f5p", names(model.coeffs))]
  f3p.coeffs <- model.coeffs[grep("f3p", names(model.coeffs))]
  f5p.coeffs <- f5p.coeffs[order(names(f5p.coeffs))]
  f3p.coeffs <- f3p.coeffs[order(names(f3p.coeffs))]
  names(f5p.coeffs) <- gsub("^f5p(.{2,2})$", names(f5p.coeffs), replace="\\1",
                            perl=TRUE)
  names(f3p.coeffs) <- gsub("^f3p(.{2,2})$", names(f3p.coeffs), replace="\\1",
                            perl=TRUE)
  names(f3p.coeffs) <- StrRev(names(f3p.coeffs))
  inds.intersect <- intersect(names(f5p.coeffs), names(f3p.coeffs))

  f5p.df <- data.frame(kd=f5p.coeffs,
                              f1=substr(names(f5p.coeffs), 1, 1),
                              f2=substr(names(f5p.coeffs), 2, 2))
  f5p.lm<- lm(kd ~ f1 + f2, data=f5p.df)
  f3p.df <- data.frame(kd=f3p.coeffs,
                              f1=substr(names(f3p.coeffs), 1, 1),
                              f2=substr(names(f3p.coeffs), 2, 2))
  f3p.lm <- lm(kd ~ f1 + f2, data=f3p.df)
  colors.nucs <- c(A="blue", C="purple", G="red", T="green")
  # 3rd plot:
  xmin <- -1
  xmax <- 3
  ymin <- xmin
  ymax <- xmax

  if (class(pdf.plot) == "character") {
    pdf(file=paste0("2017_Paper/", pdf.plot, ".pdf"), height=height*2/5, width=width*2/5)
    par(kPDFParameters)
  } else {
    dev.new(xpos=xpos, ypos=ypos, height=height, width=width)
    par(kPlotParameters)    
  }


  plot(x=f5p.lm$fitted.values, y=f5p.df$kd, pch=1, col="red", ann=FALSE,
       axes=FALSE, xlim=c(xmin, xmax), ylim=c(ymin, ymax))
  AddLinearAxis(1, tick.space=0.2, label.space=1, label="5′ dinucleotide model")
  AddLinearAxis(2, tick.space=0.2, label.space=1, label="5′ dinucleotide data")
  corr.coords <- GetPlotFractionalCoords(fx=0.9, fy=0.1)
  AddCorrelationToPlot(x=f5p.lm$fitted.values, y=f5p.df$kd, xpos=corr.coords[1],
                       ypos=corr.coords[2], rsquared=TRUE, col="red")
  points(x=f3p.lm$fitted.values, y=f3p.df$kd, pch=1, col="blue", lwd=2.5)
  corr.coords <- GetPlotFractionalCoords(fx=0.9, fy=0.05)
  AddCorrelationToPlot(x=f3p.lm$fitted.values, y=f3p.df$kd, xpos=corr.coords[1],
                       ypos=corr.coords[2], rsquared=TRUE, col="blue")
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

PlotCanonicalSiteFlankingKdHeatmap <- function(experiment="equilibrium", n_constant=5,
                             sitelist="paper", xpos=20, ypos=20, pdf.plot=FALSE) {
  width <- 18
  height <- 8
  data.matrix <- GetCanonicalSiteFlankingKdMatrix(experiment, n_constant, sitelist)
  z.mat <- ZScoreMatrix(as.matrix(log(data.matrix[,-c(1, 2)])), 1)
  rownames(z.mat) <- paste0(data.matrix$miRNA, "|", data.matrix$site)
  if (class(pdf.plot) == "character") {
    pdf(file=paste0("2017_Paper/", pdf.plot, ".pdf"), height=height*2/5, width=width*2/5)
    par(kPDFParameters)
  } else {
    dev.new(xpos=xpos, ypos=ypos, height=height, width=width)
    par(kPlotParameters)    
  }
  par(mar = c(8, 3, 2, 4.5))
  h.clust <- hclust(dist(t(z.mat)))
  dd <- as.dendrogram(h.clust)
  dd.reorder <- reorder(dd, wts=colMeans(z.mat, na.rm=TRUE))
  h.clust <- as.hclust(dd.reorder)
  h.clust <<- h.clust
  colors=rainbow( 100, s=0.8, v=1, start=0.40, end=1)
  zlims <- c(-max(abs(z.mat), na.rm=TRUE), max(abs(z.mat), na.rm=TRUE))
  # Adds the blank rows to separate the heat map.
  for (i in seq(4, 1, by=-1)) {
    z.mat <- rbind(z.mat[1:(i*6), ], rep(NaN, ncol(z.mat)), z.mat[(i*6 + 1):nrow(z.mat), ])
  }
  # Inverts the order for plotting the image device.
  y.order <- seq(nrow(z.mat), 1, by=-1)
  # Plot the image:
  image(t(z.mat[y.order, h.clust$order]), col=colors, zlim=zlims, ann=FALSE,
        axes=FALSE, useRaster=FALSE)
  kFlanks.printed <- ConvertTtoU(colnames(z.mat)[h.clust$order])
  colors <- GetColorFunction(colnames(z.mat)[h.clust$order])
  # print the flank axes
  for(i in seq(4)){
    mtext(sapply(kFlanks.printed[seq(i, 256, by=4)], substr,start=1, stop=2),
          side=1, at=seq(i, 256,by=4)*1/(256-1) + 1/(1-256),line=i*1.6-1.4,cex=0.37, col=colors[seq(i, 256, by=4)])
    mtext(sapply(kFlanks.printed[seq(i, 256, by=4)], substr,start=4, stop=5),
          side=1, at=seq(i, 256,by=4)*1/(256-1) + 1/(1-256),line=i*1.6-1.4+0.7,cex=0.37, col=colors[seq(i, 256, by=4)])
  }
  conv1 <- seq(5)*28/4+3.5-28/4
  conv2 <- (conv1 - 1)/33
  axis(2, at=conv2, labels=rev(kMirnas), las=0, lwd=0)
  y.notch <- (seq(1, 34)[-c(seq(4)*7)] - 1)/33
  axis(4, at=y.notch, labels=rev(rep(kSeedSites, 5)), las=2, lwd=0)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}


  PlotFlankLinModelNtPosCoefs <- function(experiment="equilibrium", n_constant=5,
                                        sitelist="paper", xpos=20, ypos=20,
                                        pdf.plot=FALSE) {
  width <- 5
  height <- 5
  # 4th plot:
  args <- as.list(match.call())
  flank.kds.all <- SubfunctionCall(GetFlankLinearModel, args)
  # Data analysis for the 1st plot:
  # First linear model; splitting up all flanks by 5p and 3p sequence.
  lm.split <- lm(logkd ~ mirna*site + f5p + f3p, data=flank.kds.all)

  # Data analysis for 2nd plot:
  # Get the coefficients for each of the 16 flanking dinucleotides to both
  # the 5p and 3p ends, to be able to plot them against one another to assess
  # their relative magnitude:
  model.coeffs <- c(lm.split$coefficients, f5pAA=0, f3pAA=0)
  f5p.coeffs <- model.coeffs[grep("f5p", names(model.coeffs))]
  f3p.coeffs <- model.coeffs[grep("f3p", names(model.coeffs))]
  f5p.coeffs <- f5p.coeffs[order(names(f5p.coeffs))]
  f3p.coeffs <- f3p.coeffs[order(names(f3p.coeffs))]
  names(f5p.coeffs) <- gsub("^f5p(.{2,2})$", names(f5p.coeffs), replace="\\1",
                            perl=TRUE)
  names(f3p.coeffs) <- gsub("^f3p(.{2,2})$", names(f3p.coeffs), replace="\\1",
                            perl=TRUE)
  names(f3p.coeffs) <- StrRev(names(f3p.coeffs))
  inds.intersect <- intersect(names(f5p.coeffs), names(f3p.coeffs))

  f5p.df <- data.frame(kd=f5p.coeffs,
                              f1=substr(names(f5p.coeffs), 1, 1),
                              f2=substr(names(f5p.coeffs), 2, 2))
  f5p.lm<- lm(kd ~ f1 + f2, data=f5p.df)
  f3p.df <- data.frame(kd=f3p.coeffs,
                              f1=substr(names(f3p.coeffs), 1, 1),
                              f2=substr(names(f3p.coeffs), 2, 2))
  f3p.lm <- lm(kd ~ f1 + f2, data=f3p.df)

  f5p.nuc.coeffs <- c(f5p.lm$coefficients, f1A=0, f2A=0)
  f3p.nuc.coeffs <- c(f3p.lm$coefficients, f1A=0, f2A=0)
  n1.coeffs <- f5p.nuc.coeffs[grep("f1", names(f5p.nuc.coeffs))]
  n2.coeffs <- f5p.nuc.coeffs[grep("f2", names(f5p.nuc.coeffs))]
  n3.coeffs <- f3p.nuc.coeffs[grep("f2", names(f3p.nuc.coeffs))]
  n4.coeffs <- f3p.nuc.coeffs[grep("f1", names(f3p.nuc.coeffs))]
  names(n1.coeffs) <- gsub("^.*(.$)", names(n1.coeffs), replace="\\1",
                           perl=TRUE)
  names(n2.coeffs) <- names(n1.coeffs)
  names(n3.coeffs) <- names(n1.coeffs)
  names(n4.coeffs) <- names(n1.coeffs)
  coeffs.all <- cbind(n1.coeffs, n2.coeffs, n3.coeffs, n4.coeffs)
  coeffs.all <- coeffs.all[c("A", "C", "G", "T"),]
  colnames(coeffs.all) <- c("f5p1", "f5p2", "f3p1", "f3p2")
  colors.nucs <- c(A="blue", C="purple", G="red", T="green")
  if (class(pdf.plot) == "character") {
    pdf(file=paste0("2017_Paper/", pdf.plot, ".pdf"), height=height*2/5, width=width*2/5)
    par(kPDFParameters)
  } else {
    dev.new(xpos=xpos, ypos=ypos, height=height, width=width)
    par(kPlotParameters)    
  }

  barplot(coeffs.all, col=colors.nucs[rownames(coeffs.all)], beside=TRUE)
  print(coeffs.all)
  legend("topright", bty="n", legend=names(colors.nucs), col=colors.nucs,lwd=2)
  if (class(pdf.plot) == "character") {
    dev.off()
  }
}

CompareSwappedNucleotideFlanks <- function (mirna, experiment, n_constant,
                                       sitelist, site, nt1, nt2, xpos=20, ypos=20) {
  # Extract the flanking dinucleotides:
  flank.kds.data <- GetFlankKds(mirna, experiment, n_constant, sitelist, site)
  flank.seqs <- rownames(flank.kds.data)
  N1 <- grep(paste0("^", nt1), flank.seqs)
  N2 <- grep(paste0("^.", nt1), flank.seqs)
  N3 <- grep(paste0(nt1, ".$"), flank.seqs)
  N4 <- grep(paste0(nt1, "$"), flank.seqs)
  N.inds <- c(N1, N2, N3, N4)
  M1 <- gsub(paste0("^", nt1), flank.seqs[N1], replace=nt2, perl=TRUE)
  M2 <- gsub(paste0("(^.)", nt1), flank.seqs[N2], replace=paste0("\\1", nt2), perl=TRUE)
  M3 <- gsub(paste0(nt1, "(.$)"), flank.seqs[N3], replace=paste0(nt2, "\\1"), perl=TRUE)
  M4 <- gsub(paste0(nt1, "$"), flank.seqs[N4], replace=nt2, perl=TRUE)
  M.inds <- sapply(c(M1, M2, M3, M4), function(flank) which(flank.seqs == flank))
  xmin <- 5e-4
  xmax <- 2e-1
  ymin <- xmin
  ymax <- xmax
  dev.new(xpos=xpos, ypos=ypos, height=5, width=5)
  par(kPlotParameters)
  # print(cbind(flank.seqs[N1], flank.seqs[M.inds[1:length(N1)]]))

  colors <- c(rep("blue", length(N1)), rep("green", length(N2)), rep("purple", length(N3)), rep("red", length(N4)))
  # colors <- topo.colors(length(N1))

  # len.temp <- 40
  # print(cbind(flank.seqs[N1], flank.kds.data$Mean[N1], flank.seqs[M.inds[1:length(N1)]], flank.kds.data$Mean[M.inds[1:length(N1)]])[1:len.temp,])
  plot(x=flank.kds.data$Mean[N.inds], y=flank.kds.data$Mean[M.inds], log='xy',
       axes=FALSE, ann=FALSE, xlim=c(xmax, xmin), ylim=c(ymax, ymin), col=colors)
  print(exp(mean(log(flank.kds.data$Mean[N.inds])) - mean(log(flank.kds.data$Mean[M.inds]))))
  AddLogXAxis(xmin=xmin, xmax=xmax, ymin=ymax, xlabel=nt1)
  AddLogYAxis(ymin=ymin, ymax=ymax, xmin=xmax, ylabel=nt2)
  abline(0, 1, lty = 2)
}


#USED IN FIGURES:
GetPairingFlankData <- function(mirna, site, condition,
                                experiment="equilibrium", n_constant=5,
                                sitelist="paper", mir.start=1, mir.stop=15,
                                absolute=TRUE, noconstant=FALSE) {
  print(site)
  if (absolute == TRUE) {
    folder = "/structural_analysis_PAPER_realfinal/"
  } else {
    folder = "/structural_analysis_PAPER_final/"
  }
  file = paste0(condition, "_", n_constant, "_", mir.start, "-", mir.stop,
                collapse="")
  if (noconstant == TRUE) {
    file = paste0(file, "_noconstant", collapse="")
  }
  data <- read.table(paste0(kSolexaDir, mirna, "/",
                            experiment, folder, site, "/", file, ".txt"),
                     header=TRUE)
  data <- data.frame(data)
  return(data)
}


MakeFlankScatterPlotSimple <- function(mirna, experiment, condition, n_constant, sitelist, site, mir.start, mir.stop) {
  data <- data.frame(GetPairingFlankData(mirna, experiment, condition, n_constant, sitelist, site, mir.start, mir.stop))

  print(data[1:10,])
  print(data$Flank)
  attach(data)
  par(kPlotParameters)
  plot(data$plfold^(1/(mir.stop - mir.start + 1)),
       data$accessibility,
       col=sapply(unlist(lapply(data$Flank,as.character)), GetColorFunction, alpha=0.1))
  detach(data)
}

MakeFlankScatterPlot <- function(mirna, experiment, condition, n_constant, sitelist, site, mir.start, mir.stop) {
  data <- data.frame(GetPairingFlankData(mirna, experiment, condition, n_constant, sitelist, site, mir.start, mir.stop))

  print(data[1:10,])
  print(data$Flank)
  attach(data)
  par(kPlotParameters)
  plot(by(data,Flank,function(x){mean(x$plfold^(1/(mir.stop - mir.start + 1)))}),
       by(data,Flank,function(x){mean(x$accessibility)}),
       col=sapply(unique(unlist(lapply(data$Flank,as.character))), GetColorFunction, alpha=1),
       log='xy',
       xlim = c(0.03,1),
       ylim = c(0.03,1))
  detach(data)
}

#23456789!123456789@123456789#123456789$123456789%123456789^123456789&123456789*
PlotStructECDF <- function(mirna, experiment, n_constant, sitelist, site,
                           mir.start, mir.stop, score="plfold", flank=FALSE,
                           absolute=TRUE, noconstant=FALSE, xpos=20, ypos=20,
                           xlims=c(0.2, 1), ylims=c(0, 1)) {
  # Define ranges and open plot window:
  xmin <- 1e-8
  xmax <- 1
  ymin <- 0
  ymax <- 1
  dev.new(xpos=xpos, ypos=ypos, height=5, width=5)
  par(kPlotParameters)
  par(lwd=1.5)
  # Put down empty plot:
  plot(1, type = "n", xlim=c(xmin, xmax), ylim=c(ymin, ymax), log="x",
       axes=FALSE, ann=FALSE)
  AddLogXAxis(xmin=xmin, xmax=xmax, ymin=ymin, maglabel=2,
              xlabel="Target site accessibility across from miRNA nt 1-15")
  AddLinearYAxis(ymin=ymin, ymax=ymax, xmin=xmin, tick.space=0.1,
                 label.space=0.2, ylab="Cumulative Fraction")
  conditions <- c("I_combined", "0", "0.4", "1.26", "4", "12.6", "40")
  # Assign the range for the ECDF distribution:
  ecdf.range <- exp(seq(log(xmin), log(xmax), by=0.02))
  sapply(conditions, function(condition) {
    data <- GetPairingFlankData(mirna, experiment, condition, n_constant,
                                sitelist, site, mir.start, mir.stop,
                                absolute=absolute, noconstant=noconstant)
    if (flank != FALSE) {
      data <- subset(data, Flank==flank)
    }
    win <- mir.stop - mir.start + 1
    data.intermediate <- data[[score]]^win
    print(condition)
    print(GeoMean(data.intermediate))
    ecdf.d <- ecdf(data.intermediate)
    lines(ecdf.range, ecdf.d(ecdf.range), col=kEquilCondColors[condition])
  })
  legend.names <- c("Input library",
                    "0% AG02-miR-1",
                    "0.4% AG02-miR-1",
                    "1.26% AG02-miR-1",
                    "4% AGO2-miR-1",
                    "12.6% AGO2-miR-1",
                    "40% AGO2-miR-1")
  names(legend.names) <- c("I_combined", "0", "0.4", "1.26", "4", "12.6", "40")
  legend.coords <- GetPlotFractionalCoords(xleft=xmin, xright=xmax,
                                           ybottom=ymin, ytop=ymax, fx=0.05,
                                           fy=0.95, log='x')
  legend(x=legend.coords[1], y=legend.coords[2], cex=0.9,
         legend=legend.names[conditions], lwd=1,
         col=kEquilCondColors[conditions],
         bty="n", ncol=1)
}
#23456789!123456789@123456789#123456789$123456789%123456789^123456789&123456789*

#23456789!123456789@123456789#123456789$123456789%123456789^123456789&123456789*
PlotAllFlankStructECDF <- function(mirna, experiment, condition, n_constant,
                                   sitelist, site, mir.start, mir.stop,
                                   score="plfold", absolute=TRUE,
                                   noconstant=FALSE, xpos=20, ypos=20,
                                   xlims=c(0.2, 1), ylims=c(0, 1)) {
  # Define ranges and open plot window:
  xmin <- 1e-8
  xmax <- 1
  ymin <- 0
  ymax <- 1
  dev.new(xpos=xpos, ypos=ypos, height=5, width=5)
  par(kPlotParameters)
  par(lwd=1.5)
  # Put down empty plot:
  plot(1, type = "n", xlim=c(xmin, xmax), ylim=c(ymin, ymax), log="x",
       axes=FALSE, ann=FALSE)
  AddLogXAxis(xmin=xmin, xmax=xmax, ymin=ymin, maglabel=2,
              xlabel="Target site accessibility across from miRNA nt 1-15")
  AddLinearYAxis(ymin=ymin, ymax=ymax, xmin=xmin, tick.space=0.1,
                 label.space=0.2, ylab="Cumulative Fraction")
  # Assign the range for the ECDF distribution:
  ecdf.range <- exp(seq(log(xmin), log(xmax), by=0.02))
  data.all <- GetPairingFlankData(mirna, experiment, condition, n_constant,
                                sitelist, site, mir.start, mir.stop,
                                absolute=absolute, noconstant=noconstant)
  flanks.all <- unique(data.all$Flank)
  print(flanks.all)
  for (flank in flanks.all) {
      data <- subset(data.all, Flank==flank)
    win <- mir.stop - mir.start + 1
    data.intermediate <- data[[score]]^win
    print(condition)
    print(GeoMean(data.intermediate))
    ecdf.d <- ecdf(data.intermediate)
    lines(ecdf.range, ecdf.d(ecdf.range), col=GetColorFunction(flank, alpha=0.4))
  }
  # legend.names <- c("Input library",
  #                   "0% AG02-miR-1",
  #                   "0.4% AG02-miR-1",
  #                   "1.26% AG02-miR-1",
  #                   "4% AGO2-miR-1",
  #                   "12.6% AGO2-miR-1",
  #                   "40% AGO2-miR-1")
  # names(legend.names) <- c("I_combined", "0", "0.4", "1.26", "4", "12.6", "40")
  # legend.coords <- GetPlotFractionalCoords(xleft=xmin, xright=xmax,
  #                                          ybottom=ymin, ytop=ymax, fx=0.05,
  #                                          fy=0.95, log='x')
  # legend(x=legend.coords[1], y=legend.coords[2], cex=0.9,
  #        legend=legend.names[conditions], lwd=1,
  #        col=kEquilCondColors[conditions],
  #        bty="n", ncol=1)
}
#23456789!123456789@123456789#123456789$123456789%123456789^123456789&123456789*






MakAccessibilityvsInputPlot <- function(mirna, experiment, condition, n_constant,
                                      sitelist, site, mir.start, mir.stop,pl=TRUE,
                                      noconstant=FALSE, absolute=TRUE) {
  dev.new(xpos = 20, ypos = 220, height = 5, width = 6)
  par(kPlotParameters)
  # Window size for normalizing the pl_fold with the accessibility:
  win = 1/(mir.stop - mir.start + 1)

  # Get the flanking kds:
  flank.kds <- GetFlankKds(mirna, experiment, n_constant, sitelist, site)
  kds <- flank.kds$Mean
  names(kds) <- sapply(rownames(flank.kds), function(name) {
    paste0(unlist(strsplit(name, split = ""))[-3], collapse = "")
    })
    subfunction <- function(condition_temp) {
    data <- data.frame(GetPairingFlankData(mirna, experiment, condition_temp,
                       n_constant, sitelist, site, mir.start, mir.stop,
                       noconstant=noconstant, absolute=absolute))
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

PlotPlByAUScoreBin <- function(mirna, experiment, condition, n_constant,
                                   sitelist, site, mir.start, mir.stop,
                                   xpos=20, ylimit = c(0, 1), xlimit = c(0, 1),
                                   num.bins=10, inside.range = c(0.3, 0.8),
                                   ypos=20, AU_score="AU_win",
                                   p_score="plfold", noconstant=FALSE,
                                   absolute=TRUE) 
{
  dev.new(xpos=xpos, ypos=ypos, width=5, height=5)
  par(kPlotParameters)

  xs <- seq(xlimit[1], xlimit[2], by = 0.2)
  ys <- seq(ylimit[1], ylimit[2], by = 0.1)
  plot(1, type = "n",
       xlim = xlimit,
       ylim = ylimit,
       ann = FALSE,
       axes = FALSE)

  bin.limits <- c(0, seq(inside.range[1], inside.range[2], length = num.bins), 1)
  print(bin.limits)
  b.l <- length(bin.limits)
  x <- (bin.limits[2 : b.l] + bin.limits[1 : (b.l - 1)]) / 2
  print(x)
  # Window size for normalizing the pl_fold with the accessibility:
  subfunction <- function(condition, factor, cond.color = "black") {
    data <- data.frame(GetPairingFlankData(mirna, experiment, condition,
                       n_constant, sitelist, site, mir.start, mir.stop,
                       noconstant=noconstant, absolute=absolute))

    data$plfold_bins = cut(data[[AU_score]], bin.limits,labels=x)

    agg <- aggregate(. ~ plfold_bins, data, function(x) {
      c(mean = mean(x),
        se = sd(x)/sqrt(length(x)-1))
      })

    x <- as.numeric(as.character(agg$plfold_bins))
    y <- agg[[p_score]][,1]
    y.e <- agg[[p_score]][,2]
    segments(x0=xlimit[1], y0=mean(data[[p_score]]), x1=xlimit[2], y1=mean(data[[p_score]]),col=cond.color,lty=2)

    points(x, y,
         type = "o",
         col=cond.color)
    arrows(x, y - y.e, x, y + y.e,
           col=cond.color, length=0.15, angle=90, code=3)

  }




    subfunction("I_combined", AU_score)
    subfunction(condition, AU_score,cond.color="forestgreen")
    axis(1, xs, pos = ylimit[1])
    axis(2, ys, pos = xlimit[1])
    title(xlab = "Mean AU context-score across from miRNA nt 1-15")
    title(ylab = "Per-nucleotide accessibility across from miRNA nt 1-15", line=1)
    legend(x=0.025, y = 1, bty = "n", legend = c("Input library", "4% AGO2-miR-1"), col = c("black","forestgreen"), lwd=1)

}

PlotAUScoreByPlBin <- function(mirna, experiment, condition, n_constant,
                                 sitelist, site, mir.start, mir.stop, 
                                 noconstant=FALSE, absolute=TRUE, xpos=20,
                                 ypos=20, xlimit=c(0,1), ylimit=c(0,1),num.bins=10,
                                 inside.range = c(0.3, 0.8),
                                 AU_score="AU_win", p_score="plfold")
{
  dev.new(xpos=xpos, ypos=ypos, width=5, height=5)
  par(kPlotParameters)

  xs <- seq(xlimit[1], xlimit[2], by = 0.2)

  ys <- seq(ylimit[1], ylimit[2], by = 0.05)
  plot(1, type = "n",
       xlim = xlimit,
       ylim = ylimit,
       ann = FALSE,
       axes = FALSE)

  bin.limits <- c(0, seq(inside.range[1], inside.range[2], length = num.bins), 1)
  print(bin.limits)
  b.l <- length(bin.limits)
  x <- (bin.limits[2 : b.l] + bin.limits[1 : (b.l - 1)]) / 2
  print(x)
  # Window size for normalizing the pl_fold with the accessibility:
  subfunction <- function(condition, factor, cond.color = "black") {
    data <- data.frame(GetPairingFlankData(mirna, experiment, condition,
                       n_constant, sitelist, site, mir.start, mir.stop,
                       noconstant=noconstant, absolute=absolute))

    data$plfold_bins = cut(data[[p_score]], bin.limits,labels=x)

    agg <- aggregate(. ~ plfold_bins, data, function(x) {
      c(mean = mean(x),
        se = sd(x)/sqrt(length(x)-1))
      })

    x <- as.numeric(as.character(agg$plfold_bins))
    y <- agg[[AU_score]][,1]
    y.e <- agg[[AU_score]][,2]
    segments(x0=xlimit[1], y0=mean(data[[AU_score]]), x1=xlimit[2], y1=mean(data[[AU_score]]),col=cond.color,lty=2)

    points(x, y,
         type = "o",
         col=cond.color)
    arrows(x, y - y.e, x, y + y.e,
           col=cond.color, length=0.15, angle=90, code=3)


  }




    subfunction("I_combined", AU_score)
    subfunction(condition, AU_score,cond.color="forestgreen")
    axis(1, xs, pos = ylimit[1])
    axis(2, ys, pos = xlimit[1])
    title(xlab = "Per-nucleotide accessibility across from miRNA nt 1-15")
    title(ylab = "Mean AU context-score across from miRNA nt 1-15", line=1)
    legend(x=0.025, y = 0.6, bty = "n", legend = c("Input library", "4% AGO2-miR-1"), col = c("black","forestgreen"), lwd=1)

}

MakeDoubleSubsetMatrix <- function(mirna, experiment, condition, n_constant,
                                 sitelist, site, mir.start, mir.stop, 
                                 noconstant=FALSE, absolute=TRUE, xpos=20,
                                 ypos=20, xlimit=c(0,1), ylimit=c(0,1),num.bins.AU=8,
                                 num.bins.pl=6,
                                 inside.range.AU = c(0.3, 0.7),
                                 inside.range.pl = c(0.3,0.7),
                                 AU_score="AU_cs", p_score="plfold")
{
  dev.new(xpos=xpos, ypos=ypos, width=10, height=10)

  par(kPlotParameters)
  par(mfrow=c(2,2))
  xs <- seq(xlimit[1], xlimit[2], by = 0.2)

  ys <- seq(ylimit[1], ylimit[2], by = 0.05)
  # plot(1, type = "n",
  #      xlim = xlimit,
  #      ylim = ylimit,
  #      ann = FALSE,
  #      axes = FALSE)

  bin.limits.AU <- c(0, seq(inside.range.AU[1], inside.range.AU[2], length = num.bins.AU), 1)
  bin.limits.pl <- c(0, seq(inside.range.pl[1], inside.range.pl[2], length = num.bins.pl), 1)

  print(bin.limits)
  b.l.pl <- length(bin.limits.pl)
  x.pl <- (bin.limits.pl[2 : b.l.pl] + bin.limits.pl[1 : (b.l.pl - 1)]) / 2

  b.l.AU <- length(bin.limits.AU)
  x.AU <- (bin.limits.AU[2 : b.l.AU] + bin.limits.AU[1 : (b.l.AU - 1)]) / 2


  # Window size for normalizing the pl_fold with the accessibility:
  subfunction <- function(condition, factor, cond.color = "black") {
    data <- data.frame(GetPairingFlankData(mirna, experiment, condition,
                       n_constant, sitelist, site, mir.start, mir.stop,
                       noconstant=noconstant, absolute=absolute))

    data$plfold_bins = cut(data[[p_score]], bin.limits.pl,labels=x.pl)
    data$AU_bins = cut(data[[AU_score]], bin.limits.AU,labels=x.AU)

    output = matrix(0,nrow=length(levels(data$plfold_bins)),
                      ncol=length(levels(data$AU_bins)))
    rownames(output) <- as.character(levels(data$plfold_bins))
    colnames(output) <- as.character(levels(data$AU_bins))
    for (pl_bin in levels(data$plfold_bins)) {
      for (AU_bin in levels(data$AU_bins)) {
        count <- nrow(subset(data, (plfold_bins == pl_bin & AU_bins == AU_bin)))+0.1
        output[as.character(pl_bin), as.character(AU_bin)] <- count
      }
    }
    output <- output/sum(output)
    # x <- as.numeric(as.character(agg$plfold_bins))
    # y <- agg[[AU_score]][,1]
    # y.e <- agg[[AU_score]][,2]
    # segments(x0=xlimit[1], y0=mean(data[[AU_score]]), x1=xlimit[2], y1=mean(data[[AU_score]]),col=cond.color,lty=2)

    # points(x, y,
    #      type = "o",
    #      col=cond.color)
    # arrows(x, y - y.e, x, y + y.e,
    #        col=cond.color, length=0.15, angle=90, code=3)

    return(output)
  }




    matrix_I <- subfunction("I_combined", AU_score)
    print(matrix_I)
    matrix_A <- subfunction(condition, AU_score,cond.color="forestgreen")
    print(matrix_A)
    final <- matrix_A / matrix_I
    rownames(final) <- round(as.numeric(rownames(matrix_A)), 2)
    colnames(final) <- round(as.numeric(colnames(matrix_A)), 2)
    # heatmap.2(final,Rowv=FALSE, Colv=FALSE,
    #       trace = "none",
    # cexRow = 0.9,
    # margins = c(3.5, 3),
    # na.color="gray",
    # lhei=c(1,10),
    # lwid = c(1, 9),
    # key=FALSE)
    plot(colnames(final), final[1,], col = "white",xlim = c(0, 1), ylim=c(0.01, 5),type="l")
        colramp = rainbow(nrow(final),start=0,end=0.5)

    for (row in 1:nrow(final)){
      lines(colnames(final), final[row,],col=colramp[row])
    }

    plot(rownames(final), final[,1],col="white", xlim = c(0, 1), ylim=c(0.01, 5),type="l")
    colramp = rainbow(ncol(final),start=0,end=0.5)
    for (column in 1:ncol(final)){
      lines(rownames(final), final[,column],col=colramp[column])
    }

    plot(colnames(final), final[1,], col = "white",xlim = c(0, 1),log='y', ylim=c(0.01, 5),type="l")
        colramp = rainbow(nrow(final),start=0,end=0.5)

    for (row in 1:nrow(final)){
      lines(colnames(final), final[row,],col=colramp[row])
    }

    plot(rownames(final), final[,1],col="white", xlim = c(0, 1), log='y',ylim=c(0.01, 5),type="l")
    colramp = rainbow(ncol(final),start=0,end=0.5)
    for (column in 1:ncol(final)){
      lines(rownames(final), final[,column],col=colramp[column])
    }
    dev.new()
    heatmap.2(final,Rowv=FALSE, Colv=FALSE,
          trace = "none",
    cexRow = 0.9,
    margins = c(3.5, 3),
    na.color="gray",
    lhei=c(1,10),
    lwid = c(1, 9),
    key=FALSE)

}

# Functions that work with structure-flank dataframes:
NumFlanks <- function(data) {
  flank_abundances <- aggregate(. ~ Flank, data, function(x) {
    length(x)
  })
  flank_num <- flank_abundances[, 2]
  names(flank_num) <- flank_abundances$Flank
  return(flank_num)
}

SampleByDinucleotideEnrichment <- function(dist.sample, dist.enriched, samplesize) {
  weights <- NumFlanks(dist.enriched)/NumFlanks(dist.sample)
  inds <- sample(1:nrow(dist.sample), replace=TRUE, size=samplesize,
                 prob=weights[dist.sample$Flank])
  return(inds)
}




SampleBySiteAccess <- function(dist.sample, samplesize,
                               prob.factor, exponent=1) {
  inds <- sample_int_rej(nrow(dist.sample), size=samplesize,
                         prob=dist.sample[[prob.factor]]^exponent)
  return(inds)
}


PlotFlanksSamplePl <- function(mirna, experiment, condition, n_constant,
                               sitelist, site, mir.start, mir.stop,
                               noconstant=FALSE, absolute=TRUE, xpos=20,
                               depth=20000, ypos=20, p_score="plfold", matchdist=FALSE) {
  # Window size for normalizing the pl_fold with the accessibility:
  # Window size for normalizing the pl_fold with the accessibility:
  data.I <- GetPairingFlankData(mirna, experiment, "I_combined", n_constant,
                                sitelist, site, mir.start, mir.stop,
                                noconstant=noconstant, absolute=absolute)

  data.A <- GetPairingFlankData(mirna, experiment, condition, n_constant,
                                sitelist, site, mir.start, mir.stop,
                                noconstant=noconstant, absolute=absolute)
  if (matchdist) {
    # Calculate the mean and sd of the target distribution:
    dist.A.mean <- mean(15*log(data.A[["plfold"]]))
    dist.A.sd <- sd(15*log(data.A[["plfold"]]))
    # print(dist.A.sd)
    # plot(ecdf(data.I[["plfold"]]^15), xlim=c(1e-8, 1), log='x')
    # plot(ecdf(data.A[["plfold"]]^15), col="blue", add=TRUE)

    SampleCostFunction <- function(pars) {
      tick <<- tick + 1
      par.exponent <- pars[1]
      par.size <- ceiling(Logistic(pars[2], nrow(data.I)))
      inds <- SampleBySiteAccess(data.I, par.size, "plfold", exponent=15*par.exponent)
      residual.mean <- (mean(15*log(data.I[["plfold"]][inds])) - dist.A.mean)^2
      residual.sd <- (sd(15*log(data.I[["plfold"]][inds])) - dist.A.sd)^2
      residual <- residual.mean + residual.sd
      # print(c(residual.mean, residual.sd))
      if (tick%%1000 == 0) {
        plot(ecdf(data.I[["plfold"]][inds]^15), col=ConvertRColortoRGB("red", alpha=0.1), add=TRUE)
        print(pars)
        print(residual)        
      }
      # print(c(Logistic(pars[1], 1), pars[2], 10^(pars[3])))
      return(residual)
    }
    tick <- 1
    exp.par <- optim(c(0, 0),SampleCostFunction)$par
    par.exponent <- exp.par[1]
    par.size <- ceiling(Logistic(exp.par[2], nrow(data.I)))
    print(pars)
  } else {
    par.exponent <- 1
    par.size <- 20000
  }
  inds <- SampleBySiteAccess(data.I, par.size, "plfold", exponent=15*par.exponent)
  while(min(NumFlanks(data.I[inds, ])) < 50) {
    print(min(NumFlanks(data.I[inds, ])))
    print(length(NumFlanks(data.I[inds, ])))
    inds <- c(inds, SampleBySiteAccess(data.I, par.size, "plfold", exponent=15*par.exponent))
  }
  dev.new(xpos=xpos, ypos=ypos , height=5, width=9)
  par(kPlotParameters)
  layout(matrix(c(1, 2), 1, 2, byrow=TRUE), widths=c(1, 1))
  xmin <- 1e-5
  xmax <- 1e-1
  ymin <- xmin
  ymax <- xmax
  plot(x=Norm(NumFlanks(data.I[inds, ])), y=Norm(NumFlanks(data.A)),
       col=sapply(names(NumFlanks(data.A)), GetColorFunction), log='xy',
       xlim=c(xmin, xmax), ylim=c(ymin, ymax), ann=FALSE, axes=FALSE)
  AddLogXAxis(xmin, xmax, ymin, "Sample miR-1 8mer flanking dinucleotide frequencies")
  AddLogYAxis(ymin, ymax, xmin, "Observed miR-1 8mer flanking dinucleotide frequencies")
  AddCorrelationToPlot(x=log(Norm(NumFlanks(data.A))), y=log(Norm(NumFlanks(data.I[inds, ]))),
                       xpos=1e-2, ypos=3e-2, rsquared=TRUE)
  abline(0, 1, lty = 2)
  # ECDF plot:
  xmin <- 1e-8
  xmax <- 1
  ymin <- 0
  ymax <- 1
  plot(1, type="n", xlim=c(xmin, xmax),log='x', ylim=c(ymin, ymax), axes=FALSE,
       ann=FALSE)
  # Make the ECDF for the input, 0.4 Ago, and the resampled:
  x_ecdf <- 10^seq(-8, 0, by=0.01)
  ecdf.I <- ecdf(data.I[["plfold"]]^15)
  ecdf.A <- ecdf(data.A[["plfold"]]^15)
  ecdf.I.sample <- ecdf(data.I[["plfold"]][inds]^15)
  lines(x_ecdf, ecdf.I(x_ecdf), lwd=1)
  lines(x_ecdf, ecdf.A(x_ecdf),
        col=kEquilCondColors[as.character(condition)], lwd=2)
  lines(x_ecdf, ecdf.I.sample(x_ecdf),
        col="blue", lty=2)
  AddLogXAxis(xmin, xmax, ymin,
              "8mer site accessibility across from miRNA nt 1-15", maglabel=2)
  AddLinearYAxis(ymin, ymax, 0.2, 0.2, xmin, "Cumulative fraction")
  }



PlotAllSiteKdsVsRepression <- function(experiment, n_constant, sitelist,
                                 cutoff=FALSE, bulk=FALSE,
                                 noncanon=FALSE, repression_df = repression.df, merge=FALSE, xpos=20, 
                                 ypos=20, pdf.plot=FALSE) {
  dev.new(xpos=xpos, ypos=ypos , height=5, width=5)

  par(kPlotParameters)
    xmin <- 1/2000
    xmax <- 1
    ymin <- -1.4
    ymax <- 0.4
  plot(1, type = "n",xlim=c(xmax, xmin), ylim = c(ymin,ymax),log='x',pch=19, lwd=2,ann=FALSE,axes=FALSE)

  for (mirna in c("let-7a", "miR-1", "miR-155", "miR-124", "lsy-6")) {
  data <- SitesXCounts(mirna, experiment, n_constant, sitelist)
  kds <- GetSiteKds(mirna, experiment, n_constant, sitelist)
  print(kds)
  data_seqs <- data[,1]
  names(data_seqs) <- rownames(data)

  print(data_seqs)
  site_type = "8mer"
  sites_all <- unlist(unique(subset(repression_df,mir==mirna,select=site_type)))
  print(sites_all)
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
    repression_df$site_type[which(repression_df$site_type %in% names(site_seqs_noncanonical))] <- "Noncanonical"
      print(unique(repression_df$site_type))
    sites_all <- c(site_seqs_canonical, "Noncanonical")
  }

  out <- sapply(sites_all, function(site) {
    print(site)
    if (cutoff != FALSE &
        (site %in% names(site_seqs_noncanonical) |
          site == "Noncanonical")) {
      reduced_frame <- subset(repression_df, 
        mir==mirna & site_type==site & log_kd<=cutoff,
        select=c(log_fc, log_kd))
      if (site == "9mer-m13.21") {
          print(reduced_frame)
      }
    } else {
      reduced_frame <- subset(repression_df, 
        mir == mirna & site_type == site,
        select=c(log_fc, log_kd)
      )      
    }
    print(dim(reduced_frame))
    mean_values <- colMeans(reduced_frame)
    sd_values <- apply(reduced_frame,2,sd)/sqrt(nrow(reduced_frame)-1)
    print(sd_values)
    if (bulk == TRUE) {
      mean_values[2] <- log(kds[site,]$Mean,base=2)
      sd_values[2] <- log(kds[site,]$Mean/(kds[site,]$Lower_CI),base=2)
    }
    return(c(mean_values,sd_values))
  })

  colnames(out) <- sites_all
  print(out)
  ind_remove <- !is.na(out[3,])
  print(ind_remove)
  out <- out[,ind_remove]
  print(out)
  if (noncanon == FALSE) {
    out <- out[,c("8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-A1", "6mer-m8")]
    sites_all <- colnames(out)
  }
  nosite_rep <- mean(subset(repression_df,
                            mir==mirna & site_type=="nosite",
                            select=c(log_fc))[, 1])
  if (mirna == "let-7a") {
    color <- "gray70"
    kMirnaColors[mirna] <- "gray70"
  } else {
    color <- kMirnaColors[mirna]    
  }
  arrows(2^out[2,],
         out[1,]+out[3,] - nosite_rep,
         2^out[2,],
         out[1,] - out[3,] - nosite_rep, length=0.05, lwd=1, angle=90,
         col=color, code=3)
  arrows(2^(out[2,] + out[4,]),
         out[1,] - nosite_rep,
         2^(out[2,] - out[4,]),
         out[1,] - nosite_rep, length=0.05, lwd=1, angle=90, col=color,
         code=3)
  points(2^out[2,],
         out[1,] - nosite_rep,
         col=color,
         pch=19)

  remove <- unique(which(is.na(out[1,])), which(is.na(out[1,])))
  if (length(remove) > 0) {
    out <- out[,-remove]  
  }
  sites_all <- sites_all[order(out[2,])]
  if (mirna == "let-7a") {
    all.sites <- out
  } else if (mirna == "miR-1") {
    exclude.l7.sites <- out

  } else {
    all.sites <- cbind(all.sites, out)
    exclude.l7.sites <- cbind(exclude.l7.sites, out)
  }
  }

  xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
  print(xs)
  ys <- seq(ymin,ymax,by=0.1)

  xs <- xs[xs >= xmin & xs <= xmax]
  ys <- ys[ys >= ymin & ys <= ymax]

  xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
  yl <- seq(ymin,ymax, by= 0.2)

  axis(1, at=xl,
       labels=sapply(xl, function(name) {
         eval(substitute(expression(10^x), list(x=log10(name))))
       }),
       pos=ymin, lwd=0)
  axis(1, at=xs, labels=FALSE,
       pos=ymin, lwd=2)
  # Label the axis at each order of magnitude.

  axis(2, at=yl,
       labels=round(yl,2),
       pos=xmax, las=2, lwd=0, )
  axis(2, at=ys, labels=FALSE,
       pos=xmax, lwd=2)

  legend(x=0.8,y=-0.8, legend=kMirnas, bty="n", pch=19,col=kMirnaColors[kMirnas], ncol=1)
  sites_all <<- sites_all

  text(x=1e-3, y=0, eval(substitute(expression(italic(r) == x), 
            list(x = round(cor(all.sites[2,],all.sites[1,]),3)))),
            col="gray50")
  text(x=1e-3, y=-0.07, eval(substitute(expression(italic(r) == x), 
            list(x = round(cor(exclude.l7.sites[2,],exclude.l7.sites[1,]),3)))))
  if (cutoff != FALSE) {
  text(x=1/500, y=0.1, eval(substitute(expression(log[2](italic(K)[D][italic(nc)]) <= x),
            list(x = cutoff))), cex = 0.9)
  }

  fit <- lm(all.sites[1, ] ~ log10(2^all.sites[2, ]))
  m <- fit$coefficients[2]
  b <- fit$coefficients[1]
  x_line <- 10^(seq(log10(xmin), log10(xmax),length = 20))
  y_line <- m*log10(x_line) + b
  lines(x_line, y_line, lty = 2,lwd = 1,col="gray50")
  
  fit <- lm(exclude.l7.sites[1, ] ~ log10(2^exclude.l7.sites[2, ]))
  m <- fit$coefficients[2]
  b <- fit$coefficients[1]
  x_line <- 10^(seq(log10(xmin), log10(xmax),length = 20))
  y_line <- m*log10(x_line) + b
  lines(x_line, y_line, lty = 2,lwd = 1)
  # lowess_fit <<- lowess(log_fc ~ log_kd,
  #                      data=subset(repression.df,
  #                                  site_type %in% kSeedSites,
  #                                  select = c(log_kd,log_fc)))
  # print(lowess_fit)
  # lines(2^lowess_fit$x,lowess_fit$y)
  # lowess_fit <<- lowess(log_fc ~ log_kd,
  #                      data=subset(repression.df,
  #                                  mir != "let-7a" & site_type %in% kSeedSites,
  #                                  select = c(log_kd,log_fc)))
  # print(lowess_fit)
  # lines(2^lowess_fit$x,lowess_fit$y)

  # lowess_fit <<- lowess(log_fc ~ log_kd,
  #                      data=subset(repression.df,
  #                                  mir != "let-7a" & site_type == "8mer",
  #                                  select = c(log_kd,log_fc)))
  # print(lowess_fit)
  # lines(2^lowess_fit$x,lowess_fit$y)


  if (merge == TRUE) {
  text(x=1/500, y=-0.09, eval(substitute(expression(log[2](italic(K)[D][italic(noncanon)]) <= x),
            list(x = cutoff))))
  }


  title(xlab=expression(italic(K)[D]))
  title(ylab=expression(log[2](paste("fold change"))))
  # text(x=30,y=0.1,round(cor(out[2,],out[1,]),3),cex=1.5)

}

CompareTwoSiteTypesForRepression <- function(n_constant, sitelist, sample_text,
                                 num_boxes=10,
                                 cutoff=FALSE, bulk=FALSE,
                                 repression_df = repression.df, merge=FALSE, xpos=20, 
                                 ypos=20) {
  dev.new(xpos=xpos, ypos=ypos , height=6, width=10)
  par(kPlotParameters)
  kSeedSites <- c("8mer", "7mer-m8", "7mer-A1", "6mer")
  par(mfrow=c(2, 3))
  for (mirna in kMirnas) {
    rep.df <- subset(repression_df, mir==mirna & site_type %in% kSeedSites)
    round_kd <- round(rep.df$log_kd)    
    rep.df$breaks <- round_kd
    print(breaks)
    breaks <<- breaks
    print(unique(breaks))
    rep.df$breaksbysite <- paste0(rep.df$site_type,"|", rep.df$breaks)
    print(unique(rep.df$breaksbysite))
    new_matrix <- aggregate(log_fc ~ site_type + breaks, data=rep.df, mean)
    print(new_matrix)
    output <- matrix(0,nrow=length(unique(rep.df$site_type)), ncol = length(unique(rep.df$breaks)))
    rownames(output) <- kSeedSites
    print(output)
    print("hi")
    print(sort(unique(breaks)))
    colnames(output) <- sort(unique(round_kd))
    print("hi")

    print(output)
    for (row in seq(nrow(new_matrix))) {
      print(row)
      print(output)
      print(new_matrix[row,])
      output[as.character(new_matrix[row,1]), as.character(new_matrix[row, 2])] <- new_matrix[row, 3]
    }
    print(output)
    barplot(output, beside=TRUE,col=kSiteColors[kSeedSites,],ylim = c(-2, 0.5))
    text(5, 0.2, mirna)
}

    text(10, -1.9, sample_text)
}
# CompareTwoSiteTypesForRepression(5, "paper", "new repression", num_boxes=check, xpos = 20, ypos = 20)
# dev.copy2pdf(file = "2017_Paper/Figure4G1_raw_v8.pdf")

# repression.df <- repression_old.df
# CompareTwoSiteTypesForRepression(5, "paper", "old repression", num_boxes=check, xpos = 20, ypos = 20)
# dev.copy2pdf(file = "2017_Paper/Figure4G2_raw_v8.pdf")


PlotPairwisePlot <- function(mirna, line1, line2) {
	data.line1 <- GetRisslandRepData(line=line1)
	data.line2 <- GetRisslandRepDa(line=line2)
	utrs <- grep()
}
########### TO DO FIX THIS


PlotSingleFlankKdsVsRepression <- function(mirna, experiment, n_constant,
                                           sitelist, site,
                                           cutoff=FALSE, bulk=FALSE,
                                           noncanon=FALSE, repression_df = repression.df, merge=FALSE, xpos=20, 
                                           ypos=20) {
  dev.new(xpos=xpos, ypos=ypos , height=5, width=5)
  par(kPlotParameters)
    xmin <- 1/10000
    xmax <- 10
    ymin <- -4
    ymax <- 2
  plot(1, type = "n", 
       xlim=c(xmax, xmin),
       ylim=c(ymin, ymax),
       log='x',
       pch=19,
       lwd=2,
       ann=FALSE,
       axes=FALSE)
  segments(xmin, 0, xmax, 0, lty = 2)
  # Get the flanking kds for the site:
  # kds <- GetFlankKds(mirna, experiment, n_constant, sitelist, site)
  # print(kds)
  # Subset the matrix for only the mirna and site type:
    if (mirna == "all"){
  rep.df <- data.frame(site_sequence = c(), log_fc = c(), log_kd = c(), site_type = c())
   for (mirna in kMirnas) {
           rep.df_temp <- subset(repression_df,
    mir == mirna & site_type == site, select=c(site_sequence, log_fc, log_kd, site_type))   
           rep.df_temp$log_fc <- rep.df_temp$log_fc - mean(unlist(subset(repression_df, mir==mirna & site == "nosite",select=c(log_fc))))
    rep.df <- rbind(rep.df, rep.df_temp)
   }
   } else {
      rep.df <- subset(repression_df,
    mir==mirna & site_type == site, select=c(site_sequence, log_fc, log_kd, site_type))

   }



  # Convert the sequence on the left into the flanking dinucleotides,
  # for the coloration of the points.
  ConvertSequenceToFlanks <- function(sequence){
    left <- substr(sequence, 1, 2)
    right <- substr(sequence, nchar(sequence) - 1, nchar(sequence))
    out <- paste0(left, right, collapse="")
    return(out)
  }
  rep.df$flanks <- sapply(rep.df$site_sequence, ConvertSequenceToFlanks)
  agg_fc <- aggregate(log_fc ~ flanks, rep.df, function(x) {
    c(mean = mean(c(x)),
      se = sd(c(x))/sqrt(length(c(x))-1))
    })
  agg_kd <- aggregate(log_kd ~ flanks, rep.df, function(x) {
    c(mean = mean(x))
    })
  agg_kd <<- agg_kd
  rep.df <<- rep.df
  agg_fc <<- agg_fc
  # Make the color vector:
  colors <- sapply(rep.df$flanks, GetColorFunction, alpha=0.9)
  colors_sum <- sapply(agg_fc$flanks, GetColorFunction)
  points(2^rep.df$log_kd, rep.df$log_fc, pch=19, col = colors)
    xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
  ys <- seq(ymin,ymax,by=0.1)

  xs <- xs[xs >= xmin & xs <= xmax]
  ys <- ys[ys >= ymin & ys <= ymax]

  xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
  yl <- seq(ymin,ymax, by= 0.5)

  axis(1, at=xl,
       labels=sapply(xl, function(name) {
         eval(substitute(expression(10^x), list(x=log10(name))))
       }),
       pos=ymin, lwd=0)
  axis(1, at=xs, labels=FALSE,
       pos=ymin, lwd=2)
  # Label the axis at each order of magnitude.

  axis(2, at=yl,
       labels=round(yl,2),
       pos=xmax, las=2, lwd=0, )
  axis(2, at=ys, labels=FALSE,
       pos=xmax, lwd=2)
  title(xlab=expression(italic(K)[D]))
  title(ylab=expression(log[2](paste("fold change"))))
  text(x = 10^-3, y = 1.5, mirna)
  text(x = 10^-3, y = 1.2, site)

  # plot(1, type = "n", 
  #      xlim=c(xmax, xmin),
  #      ylim=c(ymin, ymax),
  #      log='x',
  #      pch=19,
  #      lwd=2,
  #      ann=FALSE,
  #      axes=FALSE)
  # segments(xmin, 0, xmax, 0, lty = 2)

  # points(2^agg_kd[,2], agg_fc$log_fc[,1], col = colors_sum)


  # axis(1, at=xl,
  #      labels=sapply(xl, function(name) {
  #        eval(substitute(expression(10^x), list(x=log10(name))))
  #      }),
  #      pos=ymin, lwd=0)
  # axis(1, at=xs, labels=FALSE,
  #      pos=ymin, lwd=2)
  # # Label the axis at each order of magnitude.

  # axis(2, at=yl,
  #      labels=round(yl,2),
  #      pos=xmax, las=2, lwd=0, )
  # axis(2, at=ys, labels=FALSE,
  #      pos=xmax, lwd=2)


}


PlotAllCanonicalKdsVsRepression <- function(mirna, experiment, n_constant,
                                           sitelist, 
                                           cutoff=FALSE, bulk=FALSE,
                                           noncanon=FALSE, repression_df = repression.df, merge=FALSE, xpos=20, 
                                           ypos=20) {
  dev.new(xpos=xpos, ypos=ypos , height=5, width=5)
  par(kPlotParameters)
    xmin <- 1/10000
    xmax <- 10
    ymin <- -4
    ymax <- 2
  plot(1, type = "n", 
       xlim=c(xmax, xmin),
       ylim=c(ymin, ymax),
       log='x',
       pch=19,
       lwd=2,
       ann=FALSE,
       axes=FALSE)
  segments(xmin, 0, xmax, 0, lty = 2)
  # Get the flanking kds for the site:
  # Subset the matrix for only the mirna and site type:
  if (mirna == "all"){
   rep.df <- subset(repression_df,
    site_type %in% kSeedSites, select=c(site_sequence, log_fc, log_kd, site_type))
   
   } else {
      rep.df <- subset(repression_df,
    mir==mirna & site_type %in% kSeedSites, select=c(site_sequence, log_fc, log_kd, site_type))

   }
  # Convert the sequence on the left into the flanking dinucleotides,
  # for the coloration of the points.
  ConvertSequenceToFlanks <- function(sequence){
    left <- substr(sequence, 1, 2)
    right <- substr(sequence, nchar(sequence) - 1, nchar(sequence))
    out <- paste0(left, right, collapse="")
    return(out)
  }
  rep.df$flanks <- sapply(rep.df$site_sequence, ConvertSequenceToFlanks)
  agg_fc <- aggregate(log_fc ~ flanks, rep.df, function(x) {
    c(mean = mean(c(x)),
      se = sd(c(x))/sqrt(length(c(x))-1))
    })
  agg_kd <- aggregate(log_kd ~ flanks, rep.df, function(x) {
    c(mean = mean(x))
    })
  agg_kd <<- agg_kd
  rep.df <<- rep.df
  agg_fc <<- agg_fc
  # Make the color vector:
  colors.site <- sapply(kSiteColors[rep.df$site_type,],function(name){
    rgb_values <- c(col2rgb(name))/255
    red <- rgb_values[1]
    green <- rgb_values[2]
    blue <- rgb_values[3]

    return(rgb(red, green, blue ,alpha=0.2))})

    colors.flanks <- sapply(rep.df$flanks, GetColorFunction, alpha=0.1)

  colors_sum <- sapply(agg_fc$flanks, GetColorFunction)
  points(2^rep.df$log_kd, rep.df$log_fc, pch=19, col = colors.site)
    xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
  print(xs)
  ys <- seq(ymin,ymax,by=0.1)

  xs <- xs[xs >= xmin & xs <= xmax]
  ys <- ys[ys >= ymin & ys <= ymax]

  xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
  yl <- seq(ymin,ymax, by= 0.5)

  axis(1, at=xl,
       labels=sapply(xl, function(name) {
         eval(substitute(expression(10^x), list(x=log10(name))))
       }),
       pos=ymin, lwd=0)
  axis(1, at=xs, labels=FALSE,
       pos=ymin, lwd=2)
  # Label the axis at each order of magnitude.
  text(x=1, y=-3.5, eval(substitute(expression(italic(r) == x), 
            list(x = round(cor(rep.df$log_kd,rep.df$log_fc),3)))))

  axis(2, at=yl,
       labels=round(yl,2),
       pos=xmax, las=2, lwd=0, )
  axis(2, at=ys, labels=FALSE,
       pos=xmax, lwd=2)
  title(xlab=expression(italic(K)[D]))
  title(ylab=expression(log[2](paste("fold change"))))
  text(x = 10^-3, y = 1.5, mirna)
  for (site in kSeedSites) {
    lowess_object <- lowess(log_fc ~ log_kd, subset(rep.df, site_type==site),f=1)
    lines(2^lowess_object$x,lowess_object$y, col=kSiteColors[site,])
  }
  # plot(1, type = "n", 
  #      xlim=c(xmax, xmin),
  #      ylim=c(ymin, ymax),
  #      log='x',
  #      pch=19,
  #      lwd=2,
  #      ann=FALSE,
  #      axes=FALSE)
  # segments(xmin, 0, xmax, 0, lty = 2)


  # points(2^rep.df$log_kd, rep.df$log_fc, pch=1, col = colors.flanks)
  #   xs <- c(sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x))
  # print(xs)
  # ys <- seq(ymin,ymax,by=0.1)

  # xs <- xs[xs >= xmin & xs <= xmax]
  # ys <- ys[ys >= ymin & ys <= ymax]

  # xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
  # yl <- seq(ymin,ymax, by= 0.5)

  # axis(1, at=xl,
  #      labels=sapply(xl, function(name) {
  #        eval(substitute(expression(10^x), list(x=log10(name))))
  #      }),
  #      pos=ymin, lwd=0)
  # axis(1, at=xs, labels=FALSE,
  #      pos=ymin, lwd=2)
  # # Label the axis at each order of magnitude.

  # axis(2, at=yl,
  #      labels=round(yl,2),
  #      pos=xmax, las=2, lwd=0, )
  # axis(2, at=ys, labels=FALSE,
  #      pos=xmax, lwd=2)
  #   xmin <- 1/2000
  #   xmax <- 10
  #   ymin <- -1
  #   ymax <- 0.5

  # plot(1, type = "n", 
  #      xlim=c(xmax, xmin),
  #      ylim=c(ymin, ymax),
  #      log='x',
  #      pch=19,
  #      lwd=2,
  #      ann=FALSE,
  #      axes=FALSE)
  # segments(xmin, 0, xmax, 0, lty = 2)

  # points(2^agg_kd[,2], agg_fc$log_fc[,1], col = colors_sum)


  # axis(1, at=xl,
  #      labels=sapply(xl, function(name) {
  #        eval(substitute(expression(10^x), list(x=log10(name))))
  #      }),
  #      pos=ymin, lwd=0)
  # axis(1, at=xs, labels=FALSE,
  #      pos=ymin, lwd=2)
  # # Label the axis at each order of magnitude.

  # axis(2, at=yl,
  #      labels=round(yl,2),
  #      pos=xmax, las=2, lwd=0, )
  # axis(2, at=ys, labels=FALSE,
  #      pos=xmax, lwd=2)


}

PlotLinearModelKdVsRepressionCanonical <- function(
  experiment, n_constant, sitelist, cutoff=FALSE, bulk=FALSE, noncanon=FALSE,
  merge=FALSE, repression_df=repression.df, xpos=20, ypos=20) {
  # Get the flanking kds for the site:
  # Subset the matrix for only the mirna and site type:
   rep.df <- subset(repression_df,
                    site_type %in% kSeedSites,
                    select=c(mir, site_sequence, log_fc, log_kd, site_type))
  # Convert the sequence on the left into the flanking dinucleotides,
  # for the coloration of the points.
  # Define subfunction:
  ConvertSequenceToFlanks <- function(sequence) {
    left  <- substr(sequence,
                    1,
                    2)
    right <- substr(sequence,
                    nchar(sequence) - 1,
                    nchar(sequence))
    out   <- paste0(left,
                    right,
                    collapse="")
    return(out)
  }
  # Use function.
  rep.df$flank <- sapply(rep.df$site_sequence,
                         ConvertSequenceToFlanks)
  # Solve all four linear models:
  attach(rep.df)
  lm.kd  <- lm(log_kd ~ mir + site_type + flank)
  lm.kd2 <- lm(log_kd ~ mir * site_type + flank)
  lm.fc  <- lm(log_fc ~ mir + site_type + flank)
  lm.fc2 <- lm(log_fc ~ mir * site_type + flank)
  detach(rep.df)
  # Get colors for plot:
  colors.site <- sapply(kSiteColors[rep.df$site_type,],
                        ConvertRColortoRGB,
                        alpha = 0.1)
  colors.flank <- sapply(rep.df$flank,
                         GetColorFunction,
                         alpha=0.1)
  colors.mat <- cbind(colors.site, colors.flank)
  colors.site <<- colors.site
  colors.flank <<- colors.flank
  colors.mat <<- colors.mat
  colors_sum <- sapply(agg_fc$flank, GetColorFunction)
  # points(2^model.kds,2^rep.df$log_kd, pch=1, col = colors.site)

  dev.new(xpos   = xpos,
          ypos   = ypos,
          height = 8,
          width  = 12)
  par(kPlotParameters)
  par(mfrow=c(3, 4))
  xmin <- 1/10000
  xmax <- 10
  ymin <- -4.5
  ymax <- 1.5
  for (lm.ind in seq(2)) {
    for (color.ind in seq(2)) {
    x_data <- fitted(list(lm.kd, lm.kd2)[[lm.ind]])
    y_data <- rep.df$log_kd
    plot(2^x_data,
         2^y_data,
         xlim = c(xmax, xmin),
         ylim = c(xmax, xmin),
         log  = 'xy',
         pch  = 19,
         col  = colors.mat[,color.ind],
         lwd  = 2,
         ann  = FALSE,
         axes = FALSE)
    segments(xmin,
             0,
             xmax,
             0,
             lty = 2)
    xs <- c(sapply(seq(floor(log10(xmin)),
                       ceiling(log10(xmax))),
                   function(x) seq(10)*10^x))
    ys <- seq(ymin,
              ymax,
              by=0.1)
    xs <- xs[xs >= xmin & xs <= xmax]
    ys <- ys[ys >= ymin & ys <= ymax]

    xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
    yl <- seq(ymin,ymax, by= 0.5)

    axis(1,
         at     = xl,
         labels = sapply(xl,
                         function(name) {
           eval(substitute(expression(10^x),
                           list(x=log10(name))))
           }),
         pos    =xmax,
         lwd    = 0)
    axis(1,
         at     = xs,
         labels = FALSE,
         pos    = xmax,
         lwd    = 2)
    # Label the axis at each order of magnitude.
    text(x=1e-3,
         y=0.5,
         eval(substitute(expression(italic(r^2) == x), 
              list(x = round(cor(x_data, y_data)^2,3)))))
    axis(2,
         at     = xl,
         labels = sapply(xl,
                         function(name) {
           eval(substitute(expression(10^x),
                           list(x=log10(name))))
           }),
         pos    =xmax,
         lwd    = 0)
    axis(2,
         at     = xs,
         labels = FALSE,
         pos    = xmax,
         lwd    = 2)
  title(xlab=expression(italic(K)[D],predicted))
  title(ylab=expression(italic(K)[D],measured))

}}

  sapply(seq(2), function(lm.ind) {

  sapply(seq(2), function(color.ind) {
    x_data <- fitted(list(lm.fc, lm.fc2)[[lm.ind]])
    y_data <- rep.df$log_fc
    plot(0, type = "n",
         xlim = c(ymax, ymin),
         ylim = c(ymax, ymin),
         pch  = 19,
         ann  = FALSE,
         axes = FALSE)
    segments(ymin,
             0,
             ymax,
             0,
             lty = 2)
    segments(0,
             ymin,
             0,
             ymax,
             lty = 2)
    points(x_data,
           y_data,
           col = colors.mat[,color.ind],
           pch=19)

    ys <- seq(ymin,
              ymax,
              by=0.5)
    ys <- ys[ys >= ymin & ys <= ymax]

    yl <- seq(ymin, ymax, by= 1)

    axis(1,
         at     = yl,
         pos    =ymax,
         lwd    = 0)
    axis(1,
         at     = ys,
         labels = FALSE,
         pos    = ymax,
         lwd    = 2)
    # Label the axis at each order of magnitude.
    text(x=-3,
         y=0.5,
         eval(substitute(expression(italic(r^2) == x), 
              list(x = round(cor(x_data, y_data)^2,3)))))
    axis(2,
         at     = yl,
         pos    =ymax,
         lwd    = 0)
    axis(2,
         at     = ys,
         labels = FALSE,
         pos    = ymax,
         lwd    = 2)
  title(xlab=expression(log[2](paste("fold change"))))

  title(ylab=expression(log[2](paste("fold change"))))

})})


    GetPredictedValuesSiteMiR <- function(model,error=FALSE){
      if (error == TRUE) {
        extract <- "se.fit"
      } else {
        extract <- "fit"
      }
      sapply(kSeedSites, function(site) {
        sapply(kMirnas, function(mirna) {
          prediction <- predict(model,data.frame(mir=c(mirna), site_type=c(site),flank=c("AAAA")),se.fit=TRUE)
          # predict(lm.kd,data.frame(mir=c("miR-1"), site_type=c("8mer"),flank=c("AAAA")),se.fit=TRUE)
          return(c(prediction[[extract]]))
          })
      })
    }

    GetPredictedValuesFlanks <- function(model,error=FALSE){
      if (error == TRUE) {
        extract <- "se.fit"
      } else {
        extract <- "fit"
      }
      sapply(sort(unique(rep.df$flank)), function(flank) {
          prediction <- predict(model,data.frame(mir=c("let-7a"), site_type=c("6mer"),flank=c(flank)),se.fit=TRUE)
          # predict(lm.kd,data.frame(mir=c("miR-1"), site_type=c("8mer"),flank=c("AAAA")),se.fit=TRUE)
          return(c(prediction[[extract]]))
          })
    }


    y.ms.kd <<- GetPredictedValuesSiteMiR(lm.kd2)
    e.ms.kd <<- GetPredictedValuesSiteMiR(lm.kd2, error=TRUE)
    y.ms.fc <<- GetPredictedValuesSiteMiR(lm.fc2)
    e.ms.fc <<- GetPredictedValuesSiteMiR(lm.fc2, error=TRUE)

    y.f.kd <<- GetPredictedValuesFlanks(lm.kd2)
    e.f.kd <<- GetPredictedValuesFlanks(lm.kd2, error=TRUE)
    y.f.fc <<- GetPredictedValuesFlanks(lm.fc2)
    e.f.fc <<- GetPredictedValuesFlanks(lm.fc2, error=TRUE)




    color.flanks <- sapply(sort(unique(rep.df$flank)),GetColorFunction,alpha=0.5)
    color.site <- kSiteColors[rep(kSeedSites,each=5),]
    plot(1, type="n", xlim=c(-12, -2), ylim=c(-1.5, 1))
    arrows(y.ms.kd,
      y.ms.fc-e.ms.fc,
      y.ms.kd,
      y.ms.fc+e.ms.fc,
,length=0.07, angle=90, code=3,lwd=0.5)

    arrows(y.ms.kd - e.ms.kd,
      y.ms.fc,
      y.ms.kd + e.ms.kd,
      y.ms.fc,
,length=0.07, angle=90, code=3,lwd=0.5)


    kd_rep_mirsite <<- lm(c(y.ms.fc) ~ c(y.ms.kd), weights=1/c(e.ms.fc))
    kd_rep_flank <<- lm(c(y.f.fc) ~ c(y.f.kd), weights=1/c(e.f.fc))

    abline(kd_rep_flank,lty=2,col="gray")
    abline(kd_rep_mirsite,lty=2)

    points(y.ms.kd, y.ms.fc,col=color.site,pch=19)

    axis(1,
         at     = xl,
         labels = sapply(xl,
                         function(name) {
           eval(substitute(expression(10^x),
                           list(x=log10(name))))
           }),
         pos    =xmax,
         lwd    = 0)
    axis(1,
         at     = xs,
         labels = FALSE,
         pos    = xmax,
         lwd    = 2)

    axis(2,
         at     = yl,
         pos    =ymax,
         lwd    = 0)
    axis(2,
         at     = ys,
         labels = FALSE,
         pos    = ymax,
         lwd    = 2)
  title(xlab=expression(italic(K)[D],predicted))
  title(ylab=expression(log[2](paste("fold change"))))



    plot(1, type="n", xlim=c(-12, -2), ylim=c(-1.5, 1))

    arrows(y.f.kd,
      y.f.fc-e.f.fc,
      y.f.kd,
      y.f.fc+e.f.fc,
,length=0.07, angle=90, code=3,lwd=0.2)

    arrows(y.f.kd - e.f.kd,
      y.f.fc,
      y.f.kd + e.f.kd,
      y.f.fc,
,length=0.07, angle=90, code=3,lwd=0.2)

      abline(kd_rep_flank,lty=2)
      abline(kd_rep_mirsite,lty=2,col="gray")

    points(y.f.kd, y.f.fc,col=color.flanks,pch=19)
  title(xlab=expression(italic(K)[D],predicted))

  title(ylab=expression(log[2](paste("fold change"))))

    print(kd_rep_mirsite)
    print(kd_rep_flank)
}


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
    # sites.centered <- c("11mer-m3.13", "12mer-m3.14",
    #                 "11mer-m4.14", "12mer-m4.15")
    # centered_index <- sapply(sites.centered, function(site){
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
    # sites.centered <- c("11mer-m3.13", "12mer-m3.14",
    #                 "11mer-m4.14", "12mer-m4.15")
    # centered_index <- sapply(sites.centered, function(site){
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
    # sites.centered <- c("11mer-m3.13", "12mer-m3.14",
    #                 "11mer-m4.14", "12mer-m4.15")
    # centered_index <- sapply(sites.centered, function(site){
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
    # sites.centered <- c("11mer-m3.13", "12mer-m3.14",
    #                 "11mer-m4.14", "12mer-m4.15")
    # centered_index <- sapply(sites.centered, function(site){
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
    # sites.centered <- c("11mer-m3.13", "12mer-m3.14",
    #                 "11mer-m4.14", "12mer-m4.15")
    # centered_index <- sapply(sites.centered, function(site){
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
    # sites.centered <- c("11mer-m3.13", "12mer-m3.14",
    #                 "11mer-m4.14", "12mer-m4.15")
    # centered_index <- sapply(sites.centered, function(site){
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
K
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

GetKineticsData <- function(mirna, experiment, n_constant, sitelist, kDilRatio, kSubSet=FALSE) {
  sitesXcounts <- SitesXCounts(mirna, "equilibrium", n_constant, sitelist)
  print(sitesXcounts)
  sitesXcounts.kinetics <- SitesXCountsKinetics(mirna, experiment, n_constant,
                                              sitelist)
  print(sitesXcounts.kinetics)
  sitesXcounts.p <- sitesXcounts.kinetics[[1]]
  sitesXcounts.c <- sitesXcounts.kinetics[[2]]
  # 2. Load the Kds from the equilibrium fits
  params.e <- GetSiteKds(mirna, "equilibrium", n_constant, sitelist)
  # I never use the sequences, realistically
  seqs <- sitesXcounts[,1]
  # Subset the equilibrium and kinetic data matrices, such that the equilibrium
  # Matrix includes the combined input through to the zero protein, and the
  # kinetics matrices include the pulse through to the equilibrium sample.
  sitesXcounts <- sitesXcounts[, c(-1, -2)]
  indeces.to.use <- which(rowSums(sitesXcounts) > 0 & sitesXcounts[, 1] > 0)
  print(indeces.to.use)
  kData.e <- sitesXcounts[indeces.to.use, 1 : 7]
  kData.p <- sitesXcounts.p[indeces.to.use, 4 : (ncol(sitesXcounts.p) - 2)]
  kData.c <- sitesXcounts.c[indeces.to.use, 4 : (ncol(sitesXcounts.c) - 2)]
  kData.c[, 1] <- 0
  # Conditionally subset the data if flag is not FALSE:
  if (kSubSet != FALSE) {
    print(kSubSet)
    kData.e.sites <- kData.e[1:kSubSet, ]
    kData.e.none  <- colSums(kData.e[(kSubSet + 1):nrow(kData.e), ])
    kData.e       <- rbind(kData.e.sites, kData.e.none)
    rownames(kData.e)[kSubSet + 1] <- "None"
    kData.p.sites <- kData.p[1:kSubSet, ]
    kData.p.none  <- colSums(kData.p[(kSubSet + 1):nrow(kData.p), ])
    kData.p       <- rbind(kData.p.sites, kData.p.none)
    rownames(kData.p)[kSubSet + 1] <- "None"
    kData.c.sites <- kData.c[1:kSubSet, ]
    kData.c.none  <- colSums(kData.c[(kSubSet + 1):nrow(kData.c), ])
    kData.c       <- rbind(kData.c.sites, kData.c.none)
    rownames(kData.c)[kSubSet + 1] <- "None"
  }
  if (mirna == "let-7a") {
  colnames(kData.p)[2:5] <- c("2.1", "5.1", "2.2", "5.2")
  colnames(kData.c) <- colnames(kData.p)

  }
  times <- as.integer(floor(as.numeric(colnames(kData.p))))
  kData <- sapply(unique(times), GetAverageOfReplicates,
                   times=times, data=rbind(kData.p, kData.c))
  print("Kdata:")
  print(kData)
  print("done kData")
  rownames(kData) <- c(sapply(c(".p", ".c"), function(i) {
    sapply(rownames(kData.p), function(site) {
      paste0(site, i , collapse = "")
      })
    }))
  kTimes <- unique(times)/60
  colnames(kData) <- kTimes
  kInputPulseReads <- sitesXcounts.p[indeces.to.use, c("I_combined")]
kInputChaseReads <- sitesXcounts.c[indeces.to.use, c("I_combined")]
if (kSubSet != FALSE) {
  kInputPulseReads.sites <- kInputPulseReads[1:kSubSet]
  kInputPulseReads.none  <- sum(kInputPulseReads[(kSubSet + 1):length(kInputPulseReads)])
  kInputPulseReads       <- c(kInputPulseReads.sites, kInputPulseReads.none)
  kInputChaseReads.sites <- kInputChaseReads[1:kSubSet]
  kInputChaseReads.none  <- sum(kInputChaseReads[(kSubSet + 1):length(kInputChaseReads)])
  kInputChaseReads       <- c(kInputChaseReads.sites, kInputChaseReads.none)
}
# Fit the ratio of pulse to chase, in the native library,
# to correct for the actual difference in concentration from 1:1 in the
# experiment.
kPulseChaseRatio <- (sum(sitesXcounts.c[indeces.to.use,c("I")]) /
                     sum(sitesXcounts.p[indeces.to.use,c("I")]))
kInputPulseConc <- Norm(kInputPulseReads) * kLibraryConcInRxn
kInputChaseConc <- Norm(kInputChaseReads) * kLibraryConcInRxn * kPulseChaseRatio
names(kInputPulseConc) <- rownames(kData.p)
names(kInputChaseConc) <- rownames(kData.c)
# Constants throughout the optimization routine:
kAgoStockConc <- params.e["AGO","Mean"]
kIndsKDs        <<- 1:kNumSites
kIndsKoffs      <<- (kNumSites + 1):(2 * kNumSites)
kIndsContamKD   <<- 2 * kNumSites + 1
kIndsContamKoff <<- 2 * kNumSites + 2
kIndsAgo        <<- 2 * kNumSites + 3
kIndsContam     <<- 2 * kNumSites + 4
kIndsBgs        <<- (2 * kNumSites + 5):(2 * kNumSites + 16)
# Make the matrix to multiply through the adjoint equations:
kFDotXMatrix <- diag(x=1,nrow=4*kNumSites+2)
kFDotXMatrix[nrow(kFDotXMatrix) -1, 1:(2 * kNumSites)] <- 1
kFDotXMatrix[nrow(kFDotXMatrix), (2 * kNumSites + 1):(4 * kNumSites)] <- 1
# Define the total concentration in the INITIAL BINDING (kInputInitial), and in the
# chase (L).
kInputInitial <- c(kInputPulseConc, kInputChaseConc * 0)
kInput <- ((kInputInitial + c(kInputPulseConc * 0, kInputChaseConc) * kDilRatio)
           / (kDilRatio + 1))
kInputMatrix <- matrix(c(kInputInitial / (1 + kDilRatio), 
                         rep(kInput, kNumBgs - 1)),
                       nrow=length(kInput),
                       ncol=kNumBgs,
                       byrow=FALSE)
rownames(kInputMatrix) <- rownames(kData)
colnames(kInputMatrix) <- colnames(kData)
kInputTotals <- colSums(kInputMatrix)
  return(list(kData,kInputMatrix))
}

GetBgs <- function(pars, kNumSites) {
  kIndsBgs        <- (2 * kNumSites + 5):(2 * kNumSites + 16)
  return(10^pars[kIndsBgs])
}

MakeODEPars <- function(pars, kNumSites, l. = kInput) {
  pars <- 10^pars
  kIndsKDs        <- 1:kNumSites
  kIndsKoffs      <- (kNumSites + 1):(2 * kNumSites)
  kIndsContamKD   <- 2 * kNumSites + 1
  kIndsContamKoff <- 2 * kNumSites + 2
  # The site type on and off rates, the background terms for each column,
  # the and the total concentration of Ago (a) and the contaminant (b):
  print(kIndsKDs)
  print(kIndsKoffs)
  kds.a.   <- pars[kIndsKDs]
  koffs.a. <- pars[kIndsKoffs]
  kd.b     <- pars[kIndsContamKD]
  koff.b   <- pars[kIndsContamKoff]
  # Assignment of the Kds (being the ratio of the on and off rates)
  kons.a. <- koffs.a. / kds.a.
  kon.b   <- koff.b / kd.b
  # Removal of Kd parameter as it is no longer useful:
  # Calculate the vector of diluted, bound RNA:
  # Assign the parameter vector and initial conditions vector for the ODE solver:
  parms.k <- c(koffs.a., kons.a., l., koff.b, kon.b)
  return(parms.k)
}

MakeX0 <-function(pars, kNumSites, l0. = kInputInitial) {
  pars <- 10^pars
  kIndsKDs        <- 1:kNumSites
  kIndsKoffs      <- (kNumSites + 1):(2 * kNumSites)
  kIndsContamKD   <- 2 * kNumSites + 1
  kIndsContamKoff <- 2 * kNumSites + 2
  kIndsAgo        <- 2 * kNumSites + 3
  kIndsContam     <- 2 * kNumSites + 4

  # The site type on and off rates, the background terms for each column,
  # the and the total concentration of Ago (a) and the contaminant (b):
  kds.a.   <- pars[kIndsKDs]
  koffs.a. <- pars[kIndsKoffs]
  kd.b     <- pars[kIndsContamKD]
  koff.b   <- pars[kIndsContamKoff]
  A0       <- 0.4 * pars[kIndsAgo]
  B0       <- 0.4 * pars[kIndsContam]
  # Assignment of the Kds (being the ratio of the on and off rates)
  kons.a. <- koffs.a. / kds.a.
  kon.b   <- koff.b / kd.b
  # print(kons.a.)
  # print("kons.a.")
  # print(kon.b)
  # print("kon.b")
  # print(kds.a.)
  # print("kds.a.")
  # print(koffs.a.)
  # print("koffs.a.")
  # print(kd.b)
  # print("kd.b")
  # print(A0)
  # print("A0")
  # print(B0)
  # print("B0")
  # print(l0.)
  # print("l0.")
  # Get the free amount of Ago and contaminant:
  a0.b0. <- GetFreeAgoAndContam(kds.a., kd.b, l0., A0, B0)
  # print("a0.b0.")
  # print(a0.b0.)
  a0 <- a0.b0.[1]
  b0 <- a0.b0.[2]
  a0 <<- a0
  b0 <<- b0
  l0. <<- l0.
  # print("a0")
  # print(a0)
  # print("b0")
  # print(b0)

  # Calculate the initial Ago and contaminant occupancies:
  a.occs0.b.occs0 <- GetOccupanciesContaminant(a0, b0, kds.a., kd.b)
  xap0.xac0. <- l0. * a.occs0.b.occs0[[1]]
  xbp0.xbc0. <- l0. * a.occs0.b.occs0[[2]]
  # print("xap0.xac0.")
  # print(xap0.xac0.)
  # print("xbp0.xbc0.")
  # print(xbp0.xbc0.)
  # Removal of Kd parameter as it is no longer useful:
  # Calculate the vector of diluted, bound RNA:
  # Assign the parameter vector and initial conditions vector for the ODE solver:
  x0. <- c(xap0.xac0., xbp0.xbc0., a0.b0.) / (kDilRatio + 1)
  return(x0.)
}

MakeXTimeCourse <- function(x0, parms.k, kNumSites, kCScriptDir, kCScriptName, kTimesAll, kTimes, verbose.=FALSE) {
  # Load and run the ODE:
  dyn.load(kCScriptDir)
  # print(x0)
  # print(parms.k)
  x.t <- t(ode(y          = x0,
               times      = kTimesAll,
               func       = "ode_deriv",
               parms      = parms.k,
               method     = "lsodes",
               dllname    = kCScriptName,
               sparsetype = "sparseint",
               initfunc   = "ode_p_init",
               nout       = 1,
               verbose    = verbose.)[,2:(4*kNumSites+3)])
  dyn.unload(kCScriptDir)
  x <- x.t[, which(kTimesAll %in% kTimes)]
  print(x)
  print(kTimes)
  colnames(x) <- kTimes
  # Calculate total bound pulse and chase sites:
  return(x)
}


MakeModelPrediction <- function(x, bgs, data = kData, l = kInputMatrix) {
  # Combine the ago-bound and contaminant-bound site types:
  # Written out the numerator for maximum accuracy:
  # Form of equation is :        x(L - X) + B(l - x)
  #                          D * –––––––––––––––––––
  #                                (L - X)(X + B)
  kNumSites = nrow(kData)/2
  # print(dim(x))
  # print(dim(l))
  # print(length(bgs))
  # print(dim(data))
  x.ago <- x[1:(2 * kNumSites), ] 
  x.con <- x[(2 * kNumSites + 1):(4 * kNumSites), ]
  x     <- x.ago + x.con
  X     <- colSums(x)
  L     <- colSums(l)
  D     <- colSums(kData)
  B     <- bgs
  B <<- B
  L <<- L
  X <<- X
  # Transpose D-multiply x and l for row-multiplication:
  Dx.     <- D * t(x)
  Dl.     <- D * t(l)
  Dl. <<- Dl.
  Dx. <<- Dx.
  # print(Dx.)
  # print(Dl.)
  # print(B)
  pred.  <- (Dx.*L - Dx.*X + Dl.*B - Dx.*B) /
                (L*B + L*X -  B*X - X^2)
  pred   <- t(pred.)
  return(pred)
}


GetKineticsModelGeneral <- function(mirna, experiment, n_constant, sitelist, costfunction,dil=FALSE,pars=FALSE) {
  if (pars == FALSE) {
    pars <- GetSiteKoffs(mirna, experiment, n_constant, sitelist, costfunction, dil=dil)
  }
  if (dil == TRUE) {
    kDilRatio <- 10^pars[length(pars)-12]
    print(kDilRatio)
  } else {
    kDilRatio <- 11
  }
  kData_kInput <- GetKineticsData(mirna, experiment, n_constant, sitelist, kDilRatio=11)
  kData <- kData_kInput[[1]]
  kNumSites <- nrow(kData)/2
  print(kNumSites)
  # system(paste0("python SolveForOffRates/",
  #             "MakeCScriptSingleExponential.py ",
  #             kNumSites))
  # # Performa a system command to compile both of these scripts.
  # system(paste0("R CMD SHLIB SolveForOffRates/SingleODEContam_",
  #               kNumSites, ".c"))

  kInputMatrix <- kData_kInput[[2]]

# Assign names to their name and directory to be used within the optimization
# routine.
  kCScriptName <- paste0("SingleODEContam_", kNumSites)
  kCScriptDir <- paste0("SolveForOffRates/", kCScriptName, ".so")
  try(dyn.unload(kCScriptDir))
  print("hi")
  print(names(pars))
  temp_names <- names(pars)
  pars <- as.numeric(pars)

  names(pars) <- temp_names
  if (dil != FALSE) {
    kDilRatio <- 10^pars[2 * kNumSites + 5]
    pars <- pars[-(2*kNumSites + 5)]
  } else {
    kDilRatio <- 11
  }
  x0      <- MakeX0(pars, kNumSites, l0. = (kDilRatio+1) * kInputMatrix[,1])
  x0 <<- x0
  print("done x0")
  # print(x0)
  parms.k <- MakeODEPars(pars, kNumSites, l. = kInputMatrix[, 2])
  print("done parms.k")
  parms.k <<- parms.k
  # print(parms.k)
  x       <- MakeXTimeCourse(x0, parms.k, kNumSites, kCScriptDir, kCScriptName,kTimesAll, kTimes)
  print("done x")
  x <<- x
  kData <<- kData
  bgs     <- GetBgs(pars, kNumSites)
  bgs <<- bgs
  model     <- MakeModelPrediction(x, bgs, data = kData, l = kInputMatrix)
  print("done model")
  return(model)
}

PlotKineticsFit <- function(mirna, num_sites, sitelist, costfunc, dil=FALSE, col=FALSE) {
  kinetics.data <- GetKineticsData(mirna, "kinetics", num_sites, sitelist,kDilRatio=11)
  kData <- kinetics.data[[1]]
  kInputMatrix <- kinetics.data[[2]]
  data.e <- SitesXCounts(mirna, "equilibrium", 5, "paper")
  colors <- kSiteColors[rownames(data.e), ]
  model <- GetKineticsModelGeneral(mirna, "kinetics", 5, "paper", costfunc, dil=dil)
  model <<- model
  par(kPlotParameters)
  model.norm <- t(t(model + 1) / colSums(model + 1))
  data.norm <- t(t(kData + 1) / colSums(kData + 1))
  if (col != FALSE) {
    model.norm <- model.norm[, col, drop=FALSE]
    data.norm <- data.norm[, col, drop=FALSE]
  }
  plot(c(model.norm), c(data.norm),
       xlim = c(1e-8, 0.8),
       ylim = c(1e-8, 0.8),
       log  = 'xy',
       col  = colors,
       pch  = rep(c(20,1),
       each = kNumSites),
       cex = 1) 
    abline(0, 1, lty = 3)
}

PlotKineticsFitInidividual <- function(mirna, n_constant, sitelist, costfunc, log. = "", dil=FALSE) {
  site_list <- c(read.table(file = paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/AssignSiteTypes/sites.", mirna, "_", sitelist,".txt"),
    stringsAsFactors=FALSE)[,1], "None")
  print(site_list)
  kinetics.data <- GetKineticsData(mirna, "kinetics", n_constant, sitelist, kDilRatio=11)
  kData <- kinetics.data[[1]]
  kInputMatrix <- kinetics.data[[2]]
  data.e <- SitesXCounts(mirna, "equilibrium", 5, sitelist)
  model <- GetKineticsModelGeneral(mirna, "kinetics", 5, sitelist, costfunc, dil = dil)
  model <<- model
  koffs <- GetSiteKoffs(mirna, "kinetics", n_constant, sitelist, costfunc, dil=dil)
  print(koffs)
  koffs <- 10^koffs[(nrow(kData)/2+1):(nrow(kData))]
  names(koffs) <- rownames(data.e)
  print(koffs)

  times <- as.numeric(colnames(model))
  print(times)
  if (sitelist == "canonical") {
    dev.new(xpos = 20, ypos = 20, height = 7, width = 10)
    par(mfrow = c(2, 4))
  } else if (mirna == "miR-124") {

    dev.new(xpos = 20, ypos = 20, height = 8, width = 14)
    par(mfrow = c(4, 6))  
  } else {
    dev.new(xpos = 20, ypos = 20, height = 8, width = 12)
    par(mfrow = c(4, 5))  
  }
  par(kPlotParameters)

  kInput.Pulse <- Norm(kInputMatrix[1:(nrow(kData)/2),1])
  kInput.Chase <- Norm(kInputMatrix[(nrow(kData)/2 + 1):nrow(kData),2])

  kInput.Pulse <<- kInput.Pulse
  kInput.Chase <<- kInput.Chase
  kInput.Pulse <- 1
  kInput.Chase <- 1
  for (site in site_list) {
    print(site)
    color <- kSiteColors[site, ]
    times[1] <- 0.5/60
    site.inds <- grep(paste0(site,"."), rownames(model),fixed = TRUE)
    print(rownames(model))
    print(site.inds)
    model.p <- model[site.inds[1], ]/colSums(model)
    model.c <- model[site.inds[2], ]/colSums(model)
    data.p <- kData[site.inds[1], ]/colSums(kData)
    data.c <- kData[site.inds[2], ]/colSums(kData)
    # model.c[1] <- 1
    # data.c[1] <- 1

    print(model.p)
    print(model.c)
    print(data.p)
    print(data.c)

    ylim. = c(min(model.p, model.c[-1], data.p, data.c[-1]),
             max(model.p, model.c, data.p, data.c))
    # ylim. <- c(0.1, max(model.p, model.c, data.p, data.c)*2)
    print(ylim.)
    plot(times, data.p,
           xlim = c(0.3, 20000)/60,
           ylim = ylim.,
           log  = log.,
           col  = color,
           pch  = 20) 
    points(times, data.c, pch = 1,col=color)
    lines(times, model.p, col=color)
    lines(times, model.c, lty = 2, col=color)
    mtext(site, 3, line = 0.5, adj=0.4)
    print(grep(paste0(site),names(koffs),fixed=TRUE))
    dwell = 1/(koffs[site])
    segments(dwell,ylim.[1],dwell,ylim.[2], lty=2,col=color)
  }
    plot(1, type = "n", ann=FALSE, axes=FALSE)
  title(main = costfunc, cex = 2)

}
# graphics.off()
# PlotKoffsVsKds("miR-1", 5, "paper", dil=TRUE)
# dev.copy2pdf(file = "171208_miR-1_koffs_logres.pdf")
# graphics.off()
# PlotKoffsVsKds("let-7a", 5, "paper", dil=TRUE)
# dev.copy2pdf(file = "171208_let-7a_koffs_logres.pdf")
# graphics.off()
# PlotKoffsVsKds("miR-124", 5, "paper", dil=TRUE)
# dev.copy2pdf(file = "171208_miR-124_koffs_logres.pdf")
# graphics.off()
# PlotKoffsVsKds("lsy-6", 5, "paper", dil=TRUE)
# dev.copy2pdf(file = "171208_lsy-6_koffs_logres.pdf")
# graphics.off()
# PlotKineticsFitInidividual("miR-1", 5, "paper", "logres", log.="xy", dil=TRUE)
# dev.copy2pdf(file = "171208_miR-1_fits_logres.pdf")
# graphics.off()
# PlotKineticsFitInidividual("let-7a", 5, "paper", "logres", log.="xy", dil=TRUE)
# dev.copy2pdf(file = "171208_let-7a_fits_logres.pdf")
# graphics.off()
# PlotKineticsFitInidividual("miR-124", 5, "paper", "logres", log.="xy", dil=TRUE)
# dev.copy2pdf(file = "171208_miR-124_fits_logres.pdf")
# graphics.off()
# PlotKineticsFitInidividual("lsy-6", 5, "paper", "logres", log.="xy", dil=TRUE)
# dev.copy2pdf(file = "171208_lsy-6_fits_logres.pdf")
# graphics.off()

# PlotKineticsFitInidividual("miR-1", 5, "paper", "multinom", log.="xy", dil=TRUE)
# dev.copy2pdf(file = "171208_miR-1_fits_multinom.pdf")
# graphics.off()
# PlotKineticsFitInidividual("let-7a", 5, "paper", "multinom", log.="xy", dil=TRUE)
# dev.copy2pdf(file = "171208_let-7a_fits_multinom.pdf")
# graphics.off()
# PlotKineticsFitInidividual("miR-124", 5, "paper", "multinom", log.="xy", dil=TRUE)
# dev.copy2pdf(file = "171208_miR-124_fits_multinom.pdf")
# graphics.off()
# PlotKineticsFitInidividual("lsy-6", 5, "paper", "multinom", log.="xy", dil=TRUE)
# dev.copy2pdf(file = "171208_lsy-6_fits_multinom.pdf")
# graphics.off()



PlotKineticsFitInidividualPars <- function(mirna, n_constant, sitelist, costfunc, log. = "", dil=FALSE,pars=FALSE) {
  site_list <- c(read.table(file = paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/AssignSiteTypes/sites.", mirna, "_", sitelist,".txt"),
    stringsAsFactors=FALSE)[,1], "None")
  print(site_list)
  kinetics.data <- GetKineticsData(mirna, "kinetics", n_constant, sitelist,kDilRatio=11)
  kData <- kinetics.data[[1]]
  kInputMatrix <- kinetics.data[[2]]
  data.e <- SitesXCounts(mirna, "equilibrium", 5, sitelist)
  model <- GetKineticsModelGeneral(mirna, "kinetics", 5, sitelist, costfunc, dil = dil, pars=pars)
  print(model)
  print("done model")
  koffs <- GetSiteKoffs(mirna, "kinetics", n_constant, sitelist, costfunc, dil=dil)
  print(koffs)
  koffs <- 10^koffs[(nrow(kData)/2+1):(nrow(kData))]
  names(koffs) <- rownames(data.e)
  print(koffs)

  times <- colnames(model)
  if (sitelist == "canonical") {
    dev.new(xpos = 20, ypos = 20, height = 7, width = 10)
    par(mfrow = c(2, 4))
  } else {
    dev.new(xpos = 20, ypos = 20, height = 10, width = 17)
    par(mfrow = c(4, 6))  
  }
  par(kPlotParameters)

  for (site in site_list) {
    print(site)
    color <- kSiteColors[site, ]
    times[1] <- 0.5/60
    site.inds <- grep(paste0(site,"."), rownames(model),fixed = TRUE)
    print(rownames(model))
    print(site.inds)
    model.p <- model[site.inds[1], ]/colSums(model)
    model.c <- model[site.inds[2], ]/colSums(model)
    data.p <- kData[site.inds[1], ]/colSums(kData)
    data.c <- kData[site.inds[2], ]/colSums(kData)
    # print(model.p)
    # print(model.c)
    # print(data.p)
    # print(data.c)
    ylim. = c(min(model.p, model.c[-1], data.p, data.c[-1]),
             max(model.p, model.c, data.p, data.c))
    print(ylim.)

    plot(times, data.p,
           xlim = c(0.3, 20000)/60,
           ylim = ylim.,
           log  = log.,
           col  = color,
           pch  = 20) 
    points(times, data.c, pch = 1,col=color)
    lines(times, model.p, col=color)
    lines(times, model.c, lty = 2, col=color)
    mtext(site, 3, line = 2, adj=0.4)
    print(grep(paste0(site),names(koffs),fixed=TRUE))
    dwell = 1/(koffs[site])
    segments(dwell,ylim.[1],dwell,ylim.[2], lty=2,col=color)
  }
}


tSNE <- function(data_probs) {
  len <- nrow(data_probs)
  # plot(1, type = "n", xlim = c(1, 38), ylim = c(0, 1))
  output <- matrix(NaN,nrow=len,ncol = 38)
  for (row in seq(len)) {
    mir_index <- as.numeric(data_probs[row,1] + 1)
    range_initial <- c(mir_index - 15, mir_index + 15 + 7)
    seed_positions <- data_probs[row, (range_initial[1]+5):(range_initial[2]+5)]
    output[row,1:38] <- as.numeric(unlist(seed_positions))
  }
  return(output)
}

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
  data <- SitesXCountsCombined(mirna, experiment, start, stop, sitelist)
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






MakeWithWithoutProteinScatter <- function(mirna, experiment, start, stop, sitelist) {
    sitesXcounts <- SitesXCounts(mirna, experiment, start, stop, sitelist)
    colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
      colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
        return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
      }
    )
      params.file <- paste0(kSolexaDir, mirna,
                   "/", experiment, "/kds_with_structure/",
                   start, "-", stop, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg_free_protein.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.1 <- Logistic(params[1:nrow(sitesXcounts)], 10)
      k.c.stockago.1 <- 10^params["AGO"]
      bgs.1 <- 10^params["bg"]

      params.file <- paste0(kSolexaDir, mirna,
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
    # sites.centered <- c("11mer-m3.13", "12mer-m3.14",
    #                 "11mer-m4.14", "12mer-m4.15")
    # centered_index <- sapply(sites.centered, function(site){
    #                        which(legend.names == site)
    #                      }
    #                      )
    # legend.names[centered_index] <- "Centered"
    # legend.names <- unique(legend.names)
    legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=0.9, bty="n", ncol = 1)

    
}


MakeScatterAcrossWindows <- function(mirna, experiment, lim1, lim2, sitelist, method) {
    sitesXcounts <- SitesXCounts(mirna, experiment, lim1, lim1, sitelist)
    colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
      colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
        return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
      }
    )
    if (method == "free_protein") {

      params.file <- paste0(kSolexaDir, mirna,
                   "/", experiment, "/kds_with_structure/",
                   lim1, "-", lim1, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg_free_protein.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.1 <- Logistic(params[1:nrow(sitesXcounts)], 10)
      k.c.stockago.1 <- 10^params["AGO"]
      bgs.1 <- 10^params["bg"]

      params.file <- paste0(kSolexaDir, mirna,
                   "/", experiment, "/kds_with_structure/",
                   lim2, "-", lim2, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg_free_protein.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.2 <- Logistic(params[1:nrow(sitesXcounts)], 10)
      k.c.stockago.2 <- 10^params["AGO"]
      bgs.2 <- 10^params["bg"]

      } else {
      params.file <- paste0(kSolexaDir, mirna,
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


      params.file <- paste0(kSolexaDir, mirna,
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
    # sites.centered <- c("11mer-m3.13", "12mer-m3.14",
    #                 "11mer-m4.14", "12mer-m4.15")
    # centered_index <- sapply(sites.centered, function(site){
    #                        which(legend.names == site)
    #                      }
    #                      )
    # legend.names[centered_index] <- "Centered"
    # legend.names <- unique(legend.names)
    legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=0.9, bty="n", ncol = 1)

    
}

MakeScatterAcrossOrder <- function(mirna, experiment, lim, sitelist1, sitelist2, method) {
    sitesXcounts <- SitesXCounts(mirna, experiment, lim, lim, sitelist1)
    colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
      colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
        return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
      }
    )
    sitesXcounts1 <- sitesXcounts
        sitesXcounts <- SitesXCounts(mirna, experiment, lim, lim, sitelist2)
    colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
      colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
        return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
      }
    )
    sitesXcounts2 <- sitesXcounts
    print(sitesXcounts2)
    if (method == "free_protein") {

      params.file <- paste0(kSolexaDir, mirna,
                   "/", experiment, "/kds_with_structure/",
                   lim, "-", lim, "_", 1, "-", 15, "_", sitelist1, "_Basic_singlebg_free_protein.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.1 <- Logistic(params[1:nrow(sitesXcounts1)], 10)
      k.c.stockago.1 <- 10^params["AGO"]
      bgs.1 <- 10^params["bg"]

      params.file <- paste0(kSolexaDir, mirna,
                   "/", experiment, "/kds_with_structure/",
                   lim, "-", lim, "_", 1, "-", 15, "_", sitelist2, "_Basic_singlebg_free_protein.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.2 <- Logistic(params[1:nrow(sitesXcounts2)], 10)
      k.c.stockago.2 <- 10^params["AGO"]
      bgs.2 <- 10^params["bg"]

      } else {
      params.file <- paste0(kSolexaDir, mirna,
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


      params.file <- paste0(kSolexaDir, mirna,
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
    # sites.centered <- c("11mer-m3.13", "12mer-m3.14",
    #                 "11mer-m4.14", "12mer-m4.15")
    # centered_index <- sapply(sites.centered, function(site){
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

sites.centered <- c("11mer-m3.13", "12mer-m3.14",
                    "11mer-m4.14", "12mer-m4.15")


kSiteColors <- read.table(
  "/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/site_info.txt",
  row.names=1, header=FALSE, sep="\t", stringsAsFactors=FALSE)[, 1, drop = FALSE]
temp.names <- rownames(kSiteColors)
kSiteColors <- kSiteColors[, 1]
names(kSiteColors) <- temp.names
rm(temp.names)





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
  out.file <- paste0(kSolexaDir, mirna,
                     "/equilibrium/kds/", site,"_flanking_",
                     k.c.stockago, extension,".txt")
  write.table(file=out.file, out, sep="\t", quote=FALSE, row.names=FALSE,
                col.names=TRUE)
}

WriteFinalParameterFile <- function(out,extension="") {
  out.final <- out[dim(out)[1], ]
  names(out.final) <- colnames(out)
  out.file <- paste0(kSolexaDir, mirna,
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
    sitesXcounts <- SitesXCountsCombined(mirna, experiment, start, stop, sitelist)
  } else {
    sitesXcounts <- SitesXCounts(mirna, experiment, start, stop, sitelist)
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


  sites.centered <- c("11mer-m3.13", "12mer-m3.14", "11mer-m4.14",
                      "12mer-m4.15", "13mer-m3.15", "11mer-m5.15")
  centered_rename <- c("Centered-m3.13", "Centered-m3.14", "Centered-m4.14",
                       "Centered.m4-15", "Centered-m3.15", "Centered-m5.15")
  names(centered_rename) <- sites.centered

  legend.names <- sapply(legend.names, function(site) {
                           if (site %in% sites.centered) {
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
    sitesXcounts <- SitesXCountsCombined(mirna, experiment, start, stop, sitelist)
  } else {
    sitesXcounts <- SitesXCounts(mirna, experiment, start, stop, sitelist)
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
  sites.centered <- c("11mer-m3.13", "12mer-m3.14",
                  "11mer-m4.14", "12mer-m4.15", "13mer-m3.15", "11mer-m5.15")
  centered_rename <- c("Centered-m3.13", "Centered-m3.14", "Centered-m4.14", "Centered.m4-15", "Centered-m3.15", "Centered-m5.15")
  names(centered_rename) <- sites.centered
  print(kSiteColors[site_list, ])
  print(kSiteColors[ordered_list, ])

  legend.names <- sapply(legend.names, function(site){
                        if (site %in% sites.centered) {
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
    sitesXcounts <- SitesXCountsCombined(mirna, experiment, start, stop,
                                            sitelist)
  } else {
    sitesXcounts <- SitesXCounts(mirna, experiment, start, stop, sitelist)
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
  sfXc <- SiteFlanksXCounts(mirna, experiment, site, start, stop, sitelist)
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

  kds_centered <- matrix(unlist(sapply(kMirnas, function(mirna) {
    kds <- unlist(GetKds(mirna,"equilibrium", start, stop, sitelist,
                  combined.input=combined.input, log.residual=log.residual,
                  nosite=nosite))
    print(kds)
    kds <- c(kds["8mer"],kds["6mer"], kds["6mer-m8"],kds[grep("\\.",names(kds)),drop=FALSE],kds["None"])
    kds <- kds[grep("23", names(kds), invert=TRUE),drop=FALSE]
    return(kds)
  })),nrow=length(kMirnas),byrow=TRUE)
  rownames(kds_centered) <- kMirnas
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
  legend("topleft",legend = kMirnas, col=c("black", "blue", "red", "forestgreen", "purple"), lwd=2, pch=19)

}



MakeInputScatter <- function(mirna, experiment, n_constant, ombined.input=TRUE,
                             log.residual=FALSE, bgoff=FALSE,
                             vertical.lines=FALSE) {
  params <- GetSiteKds(mirna, experiment, start, stop, sitelist,
                   combined.input=combined.input, log.residual=log.residual,
                   scaled=FALSE)
  if (combined.input == TRUE) {
    sitesXcounts <- SitesXCountsCombined(mirna, experiment, start, stop, sitelist)
  } else {
    sitesXcounts <- SitesXCounts(mirna, experiment, start, stop, sitelist)
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
  sites.centered <- c("11mer-m3.13", "12mer-m3.14",
                  "11mer-m4.14", "12mer-m4.15", "13mer-m3.15", "11mer-m5.15")
  centered_rename <- c("Centered-m3.13", "Centered-m3.14", "Centered-m4.14", "Centered.m4-15", "Centered-m3.15", "Centered-m5.15")
  names(centered_rename) <- sites.centered
  print(kSiteColors[site_list, ])
  print(kSiteColors[ordered_list, ])

  legend.names <- sapply(legend.names, function(site){
                        if (site %in% sites.centered) {
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
    sitesXcounts <- SitesXCounts(mirna, experiment, start, stop, sitelist)
    colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
      colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
        return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
      }
    )
      params.file <- paste0(kSolexaDir, mirna,
                   "/", experiment, "/kds_with_structure/",
                   start, "-", stop, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg_free_protein.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.1 <- Logistic(params[1:nrow(sitesXcounts)], 10)
      k.c.stockago.1 <- 10^params["AGO"]
      bgs.1 <- 10^params["bg"]

      params.file <- paste0(kSolexaDir, mirna,
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
    # sites.centered <- c("11mer-m3.13", "12mer-m3.14",
    #                 "11mer-m4.14", "12mer-m4.15")
    # centered_index <- sapply(sites.centered, function(site){
    #                        which(legend.names == site)
    #                      }
    #                      )
    # legend.names[centered_index] <- "Centered"
    # legend.names <- unique(legend.names)
    legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=1.2, bty="n", ncol = 1)

    
}

MakeSMNBScatter <- function(start, stop, sitelist) {
    sitesXcounts <- SitesXCounts("let-7", "equilibrium", start, stop, sitelist)
    colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
      colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
        return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
      }
    )
      params.file <- paste0(kSolexaDir, mirna,
                   "/", experiment, "/kds_with_structure/",
                   start, "-", stop, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg_free_protein.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.1 <- Logistic(params[1:nrow(sitesXcounts)], 10)
      k.c.stockago.1 <- 10^params["AGO"]
      bgs.1 <- 10^params["bg"]

      params.file <- paste0(kSolexaDir, mirna,
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
    # sites.centered <- c("11mer-m3.13", "12mer-m3.14",
    #                 "11mer-m4.14", "12mer-m4.15")
    # centered_index <- sapply(sites.centered, function(site){
    #                        which(legend.names == site)
    #                      }
    #                      )
    # legend.names[centered_index] <- "Centered"
    # legend.names <- unique(legend.names)
    legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=1.2, bty="n", ncol = 1)

    
}

MakeScatterAcrossWindows <- function(mirna, experiment, lim1, lim2, sitelist, method) {
    sitesXcounts <- SitesXCounts(mirna, experiment, lim1, lim1, sitelist)
    colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
      colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
        return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
      }
    )
    if (method == "free_protein") {

      params.file <- paste0(kSolexaDir, mirna,
                   "/", experiment, "/kds_with_structure/",
                   lim1, "-", lim1, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg_free_protein.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.1 <- Logistic(params[1:nrow(sitesXcounts)], 10)
      k.c.stockago.1 <- 10^params["AGO"]
      bgs.1 <- 10^params["bg"]

      params.file <- paste0(kSolexaDir, mirna,
                   "/", experiment, "/kds_with_structure/",
                   lim2, "-", lim2, "_", 1, "-", 15, "_", sitelist, "_Basic_singlebg_free_protein.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.2 <- Logistic(params[1:nrow(sitesXcounts)], 10)
      k.c.stockago.2 <- 10^params["AGO"]
      bgs.2 <- 10^params["bg"]

      } else {
      params.file <- paste0(kSolexaDir, mirna,
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


      params.file <- paste0(kSolexaDir, mirna,
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
    # sites.centered <- c("11mer-m3.13", "12mer-m3.14",
    #                 "11mer-m4.14", "12mer-m4.15")
    # centered_index <- sapply(sites.centered, function(site){
    #                        which(legend.names == site)
    #                      }
    #                      )
    # legend.names[centered_index] <- "Centered"
    # legend.names <- unique(legend.names)
    legend(x=xmax,y=xmax,legend=legend.names, pch=19, col=kSiteColors[legend.names, ], cex=0.9, bty="n", ncol = 1)

    
}

MakeScatterAcrossOrder <- function(mirna, experiment, lim, sitelist1, sitelist2, method) {
    sitesXcounts <- SitesXCounts(mirna, experiment, lim, lim, sitelist1)
    colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
      colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
        return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
      }
    )
    sitesXcounts1 <- sitesXcounts
        sitesXcounts <- SitesXCounts(mirna, experiment, lim, lim, sitelist2)
    colnames(sitesXcounts)[3:dim(sitesXcounts)[2]] <- sapply(
      colnames(sitesXcounts)[3:dim(sitesXcounts)[2]], function(x) {
        return(unlist(strsplit(x, split="A", fixed=TRUE))[2])
      }
    )
    sitesXcounts2 <- sitesXcounts
    print(sitesXcounts2)
    if (method == "free_protein") {

      params.file <- paste0(kSolexaDir, mirna,
                   "/", experiment, "/kds_with_structure/",
                   lim, "-", lim, "_", 1, "-", 15, "_", sitelist1, "_Basic_singlebg_free_protein.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.1 <- Logistic(params[1:nrow(sitesXcounts1)], 10)
      k.c.stockago.1 <- 10^params["AGO"]
      bgs.1 <- 10^params["bg"]

      params.file <- paste0(kSolexaDir, mirna,
                   "/", experiment, "/kds_with_structure/",
                   lim, "-", lim, "_", 1, "-", 15, "_", sitelist2, "_Basic_singlebg_free_protein.txt")
      params <- read.table(params.file, header=TRUE, stringsAsFactors = FALSE)
      params <- params[nrow(params), - ncol(params)]
      kds.2 <- Logistic(params[1:nrow(sitesXcounts2)], 10)
      k.c.stockago.2 <- 10^params["AGO"]
      bgs.2 <- 10^params["bg"]

      } else {
      params.file <- paste0(kSolexaDir, mirna,
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


      params.file <- paste0(kSolexaDir, mirna,
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
    # sites.centered <- c("11mer-m3.13", "12mer-m3.14",
    #                 "11mer-m4.14", "12mer-m4.15")
    # centered_index <- sapply(sites.centered, function(site){
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


ModelFunctionDist <- function(pars, c.I.tots, data) {
  # print(c.totals[1,1])
  # print(c.I.tots[1])
  # Split up the parameters into the kd and background parameters.
  num.kds <- length(c.I.tots)
  kds  <- 10^pars[1 : num.kds]
  bgs  <- rep(10^pars[(num.kds + 1)], ncol(data))
  stock.ago <- 10^pars[num.kds + 2]
  alpha <- 10^pars[num.kds + 3]
  beta <- 10^pars[num.kds + 4]
  k <- 1

  # Get the bound counts in each exp:
  c.agos <- sapply(colnames(data), function(x) {
    as.numeric(x) / 100
  })
  c.totals <- matrix(
                rep(c.I.tots,ncol(data)), nrow=length(c.I.tots), ncol=ncol(data), byrow=FALSE)


  n_dist <- 100
  kds.dist <- c(sapply(1:length(kds), function(kdind) {
    p.dist <- rbeta(n_dist, shape1=alpha, shape2=beta)
    k.dist <- (kds[kdind])/(p.dist^k)
    return(k.dist)
    }))
  kds <<- kds
  kds.dist <<- kds.dist
  # kds.dist[(length(kds.dist)-n_dist+1):length(kds.dist)] <- kds[length(kds)]
  c.I.tots.dist <- c(sapply(c.I.tots, function(conc) {
    rep(conc, n_dist)/n_dist
  }))
  c.totals.dist <- matrix(
                rep(c.I.tots.dist,ncol(data)), nrow=length(c.I.tots.dist), ncol=ncol(data), byrow=FALSE)
    c.bounds.dist <- as.matrix(
    sapply(c.agos, function(percent) {
      return(GetBoundRNA(kds.dist, c.I.tots.dist, percent * stock.ago))
    }
  ))
  c.bounds.dist.total <- matrix(0, nrow=nrow(data), ncol = ncol(data))
    for (row in seq(length(c.I.tots))) {
    rowinds <- (1:n_dist)+(row-1)*n_dist
    c.bounds.dist.total[row,] <- colSums(c.bounds.dist[rowinds,])
  }
  c.frees.dist <- c.totals - c.bounds.dist.total
  c.bgs.dist <- t(t(c.frees.dist) * bgs / colSums(c.frees.dist))
  c.all.dist <- c.bounds.dist.total + c.bgs.dist
  c.final <- data.frame(t(t(c.all.dist) / colSums(c.all.dist) * colSums(data)))

  return(c.final)
}

plotBeta <- function(pars) {
  num.kds <- length(pars) - 5
      alpha <- 10^(pars[num.kds + 3])
    beta  <- 10^(pars[num.kds + 4])
    x.int <- seq(0, 1, length = 100)

    plot(x.int, dbeta(x.int, shape1 = alpha, shape2 = beta), type = "l")

}


ModelLikelihoodDist <- function(pars, data, input){
  pars.current <<- pars
  model <<- ModelFunctionDist(pars, input, data)
  model.norm <- t(t(model) / colSums(model))
  data.norm <- t(t(data) / colSums(data))
  loglikelihood <- -sum(data*log(model.norm))
  model.R <- model.norm/Norm(input)
  data.R <- data.norm/Norm(input)
  loglikelihood <- sum((model.R - data.R)^2)
  if (tick %% 100 == 0) {
    par(mfrow=c(1, 2))
    plot(0, type = "n", log = 'xy', xlim = c(.1,100), ylim=c(0.1, 300))
    colors.sites <- c(sapply(rownames(data.R)[1:256], GetColorFunction), kSiteColors[rownames(data.R)[257:nrow(data.R)], ])
    # print(colors.sites)
    # print(data.R[1:10, ])
    # print(model.R[1:10, ])
  sapply(c(seq(1, 256, by =16), 257:nrow(model)), function(row) {
    points(colnames(data),c(data.R[row,]),col=colors.sites[row], pch = 20)
    lines(colnames(data),c(model.R[row,]),col=colors.sites[row], pch = 20)
    num.kds <- nrow(model)
  })
    plotBeta(pars)
    # plot(c(model.R), c(data.R), log = 'xy',col=kSiteColors[rownames(data.R),])
    print(loglikelihood)
    print(pars[(length(pars)-3):length(pars)])

  }

# }
tick <<- tick + 1
  return(loglikelihood)
}



CompareDist <- function(pars=GetSiteKds("miR-1", "equilibrium", 5, "paper")) {
  sitesXcounts <- SitesXCounts("miR-1", "equilibrium", 5, "paper")
  print(sitesXcounts)
  # Separate site sequences from data file.
  if (sitelist %in% kKmerSiteLists) {
    seqs <- rownames(sitesXcounts)
  } else {
    seqs <- sitesXcounts[,1]
    sitesXcounts <- sitesXcounts[,-1]
  }

  # Remove any rows with no reads whoatsoever.
  sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
  # Remove any rows for which there are no input reads:
  sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]
  sitesXcounts <<- sitesXcounts
  sfXc <- SiteFlanksXCounts("miR-1", "equilibrium", 5, "paper", "8mer")
  pars.flanks <- GetFlankKds("miR-1", "equilibrium", 5, "paper", "8mer")
  sfXc <<- sfXc
  sitesXcounts <- sitesXcounts[-1,]
  sitesXcounts <- rbind(sfXc, sitesXcounts)
  data <- sitesXcounts[,3:7]
  data0 <- sitesXcounts[,8]
  names(data0) <- rownames(data)
  data <- data.matrix(data)
  data <<- data
  data0 <<- data0
  colnames(data) <- colnames(sitesXcounts)[3:7]
  data.norm <- t(t(data) / colSums(data))
  I.seq <- sitesXcounts[,2]
  c.I.tots <- Norm(I.seq)*100
  c.I.tots <<- c.I.tots
  c.totals <- matrix(
                rep(c.I.tots,ncol(data)), nrow=length(c.I.tots), ncol=ncol(data), byrow=FALSE)
  colnames(c.totals) <- colnames(data)
  pars.init <- log10(c(c(pars.flanks$Mean, (pars$Mean)[2:(nrow(pars)-2)])^(0.8), ((pars$Mean)[(nrow(pars)-1)])/10,pars$Mean[nrow(pars)], 1, 2))
  


  names(pars.init) <- c(rownames(data), "bg", "AGO", "alpha", "beta")
  print(pars.init)
  pars.init <<- pars.init
  opt  <- optim(pars.init, ModelLikelihoodDist, method="BFGS",data=data, input = c.I.tots, control=list(maxit = 100000, parscale = rep(10,length(pars.init))))
  par.round1 <<- opt$par
  opt2 <- optim(opt$par, ModelLikelihoodDist, data=data, input = c.I.tots, control=list(maxit = 100000))
  par.round2 <<- opt2$par
  opt3 <- optim(opt2$par, ModelLikelihoodDist, data=data, input = c.I.tots, control=list(maxit = 100000))
  par.round3 <<- opt3$par
  opt4 <- optim(opt3$par, ModelLikelihoodDist, data=data, input = c.I.tots, control=list(maxit = 100000))
  par.round4 <<- opt4$par
  opt5 <- optim(opt4$par, ModelLikelihoodDist, data=data, input = c.I.tots, control=list(maxit = 100000))
  par.round5 <<- opt5$par

  print(opt5$par)


}