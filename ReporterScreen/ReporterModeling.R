source("general/general.R")

source("general/GenericFigures.R")


graphics.off()



counts_v1_1 <- MakeReporterMatrix(rep=1)
counts_v1_2 <- MakeReporterMatrix(rep=2)


counts_v2 <- MakeReporterMatrixV2test()



break



SimulateTransformation <- function(colonies, variants=29992, reps=10) {
	outs <- sapply(1:reps, function(rep_i) {
		col_vec <- rep(0, variants)
		inds <- sample(1:variants, colonies, replace=TRUE)
		for (ind in inds) {
			col_vec[ind] <- col_vec[ind] + 1
		}
		length(which(col_vec == 0))		
	})
	c(mean(outs), sd(outs))
}

GetColonyHistogram <- function(colonies, variants=29992, reps=10, xpos=20,
                               ypos=20, height=5, width=5, pdf.plot=FALSE) {
	outs <- sapply(1:reps, function(rep_i) {
		col_vec <- rep(0, variants)
		inds <- sample(1:variants, colonies, replace=TRUE)
		for (ind in inds) {
			col_vec[ind] <- col_vec[ind] + 1
		}
		sapply(0:round(colonies/variants*10, 0), function(count_i) {
			length(which(col_vec == count_i))
		})
	})
	# c(mean(outs), sd(outs))
	outs

	nonzeros <- which(rowSums(outs) != 0)
	xmin <- 0
	xmax <- max(nonzeros)
	ymin <- 0
	ymax <- round(max(outs*1.2), 0)
	SubfunctionCall(FigureSaveFile)
	BlankPlot()
	par(kPlotParameters)
	AddLinearAxis(1, tick.space=1, label.space=5, label="Counts")
	AddLinearAxis(2, tick.space=round(ymax/10, 0),
	              label.space=round(ymax/10, 0), label="Number of Variants")
	apply(outs, 2, function(column) {
		print(column)
		print(xmin:xmax)
		lines(xmin:xmax, column[(xmin + 1):(xmax + 1)])
		})
	xy <- GetPlotFractionalCoords(0.05, 0.95)
	text(xy[1], xy[2], label=sprintf("%s colonies", colonies), adj=0)
	outs
}




# x <- c(seq(1, 9)*1e4, seq(1, 10)*1e5)

# y <- sapply(x, SimulateTransformation)

# plot(x, y[1, ], log='x')
# arrows(x, y[1, ] + y[2, ], x, y[1, ] - y[2, ],
#        length=error_height_final*par()$cex, lwd=par()$cex, angle=90, code=3)

# hist <- GetColonyHistogram(2e5, xpos=20, ypos=20)
# hist <- GetColonyHistogram(5e5, xpos=520, ypos=20)
hist1 <- GetColonyHistogram(1e6, xpos=1020, ypos=20)
hist2 <- GetColonyHistogram(2e6, xpos=20, ypos=520)
hist3 <- GetColonyHistogram(3e6, xpos=20, ypos=520)
# hist <- GetColonyHistogram(5e6, xpos=520, ypos=520)
# print(hist)




break




CountDinucsSeq <- function(sequence) {
	output <- rep(0, 16)
	names(output) <- GetKmerList(2)
	strlist <- unlist(strsplit(sequence, split=""))
	# print(sequence)
	# print(strlist)
	dinucs <- paste0(strlist[1:(nchar(sequence) - 1)], strlist[2:nchar(sequence)])
	sapply(dinucs, function(dinuc) {
		output[dinuc] <<- output[dinuc] + 1
	})
	output
}


CountDinucs <- function(row, mirna, site, proofread=TRUE) {
	# print(row)
	seq <- row[4]
	# message(seq)
	if (length(grep(",", row[5])) != 0) {
		sites <- unlist(strsplit(row[5], split=","))
		if (proofread) {
			l <- sapply(sites, function(i) {
				as.integer(unlist(strsplit(
					unlist(strsplit(i, split="\\|"))[3],
					split="\\."))[1])
				})
			r <- sapply(sites, function(i) {
				as.integer(unlist(strsplit(
					unlist(strsplit(i, split="\\|"))[3],
					split="\\."))[2])
			})
			l <- min(l)
			r <- max(r)
		} else {
			limits <- sapply(sites, function(i) {
					unlist(strsplit(i, split="\\|"))[3]
				})
			# mir_i <- sapply(sites, function(i) {
			# 		unlist(strsplit(i, split="\\|"))[1]
			# 	})
			# sites_i <- sapply(sites, function(i) {
			# 		unlist(strsplit(i, split="\\|"))[2]
			# 	})
			# ind_use <- which(mir_i == mirna)
			# site_use <- which(sites_i == site & mir_i == mirna)
			# if (length(site_use) > 1) {
			# 	limits_use <- limits[site_use]
			# 	# print(limits_use)
			# 	# print(site_use)
			# 	cents <- sapply(limits_use, function(limit_use) {
			# 		lims <- as.numeric(unlist(strsplit(limit_use, split="\\.")))
			# 		mean(lims)
			# 		})
			# 	print(cents)

			# 	ind_use <- which((cents - nchar(seq)/2)^2 == min((cents - nchar(seq)/2)^2))
			# 	print(ind_use)
			# 	site_use <- site_use[ind_use]

			# }
			cents <- sapply(limits, function(limit_use) {
				lims <- as.numeric(unlist(strsplit(limit_use, split="\\.")))
				mean(lims)
			})

			site_use <- which((cents - nchar(seq)/2)^2 == min((cents - nchar(seq)/2)^2))
			limit <- unlist(strsplit(limits[site_use], split="\\."))
			# print(limit)
			l <- limit[1]
			r <- limit[2]
		}
	} else {
		sites <- as.character(row[5])
		limits <- unlist(strsplit(sites, split="\\|"))[3]
		l <- unlist(strsplit(limits, split="\\."))[1]
		r <- unlist(strsplit(limits, split="\\."))[2]
	}
	seq_ll = substr(seq, 1, 17)
	# message(seq_ll)
	seq_l = substr(seq, 18, as.integer(l))
	# message(paste0(paste0(rep(" ", 17), collapse=""), seq_l))
	seq_site = substr(seq, as.integer(l) + 1, as.integer(r))
	# message(paste0(paste0(rep(" ", as.integer(l)), collapse=""), seq_site))
	seq_r = substr(seq, as.integer(r) + 1, nchar(seq) - 17)
	# message(paste0(paste0(rep(" ", as.integer(r)), collapse=""), seq_r))
	seq_rr = substr(seq, nchar(seq) - 16, nchar(seq))
	# message(paste0(paste0(rep(" ", nchar(seq) - 17), collapse=""), seq_rr))

	count_l <- CountDinucsSeq(seq_l)
	count_r <- CountDinucsSeq(seq_r)
	if (length(count_l) != length(count_r)) {
		print("error")
		print(seq_l)
		print(seq_r)
		break
	}
	Norm(count_l[names(count_l)] + count_r[names(count_r)])
}


GetmRNAUTRFreqs <- function() {

    utr_path <- file.path("/lab/bartel4_ata/kathyl/RNA_Seq/analysis", 
			              "data/no_baseline_analysis/merged.txt")
    data <- read.table(utr_path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
    utrs <- data$sequence
    utr_lengths <- data$utr_length
    utr_length_manual <- sapply(utrs, nchar)
    plot(utr_lengths, utr_length_manual)
    dinuc_freqs <- sapply(utrs, CountDinucsSeq)
    colnames(dinuc_freqs) <- seq(ncol(dinuc_freqs))
    dinuc_freqs
}

freqs_raw <- read.table("variants/mRNA_UTR_frequences.txt", header=TRUE, row.names=1)
freqs <- freqs_raw[sort(rownames(freqs_raw)), , drop=FALSE]



# freqs_alt <- GetmRNAUTRFreqs()

# break
print(freqs)
print(freqs_alt)
freqs_alt_norm <- t(t(freqs_alt) / colSums(freqs_alt))
print(colSums(freqs_alt_norm))
freqs_norm_mean <- rowMeans(freqs_alt_norm)
freqs_norm_sd <- apply(freqs_alt_norm, 1, sd)
# print(freqs_alt_norm)
freqs_alt_with_UTR_length <- Norm(rowSums(freqs_alt))
plot(freqs_norm_mean, freqs_alt_with_UTR_length)

arrows(freqs_norm_mean - freqs_norm_sd, freqs_alt_with_UTR_length,
       freqs_norm_mean + freqs_norm_sd, freqs_alt_with_UTR_length, code=0)
# arrows(unlist(freqs), dinucs_mean_nr - dinucs_sd_nr, unlist(freqs), dinucs_mean_nr + dinucs_sd_nr, code=0, col="red")
# points(unlist(freqs), dinucs_mean_nr, col="red")



PlotDinucsFreqs <- function(mirna, site) {
	freqs_raw <- read.table("variants/mRNA_UTR_frequences.txt", header=TRUE, row.names=1)
	freqs <- freqs_raw[sort(rownames(freqs_raw)), , drop=FALSE]
	path <- paste0("variants/oct12_final/", mirna, "/", site, "_library_variants_oct12_final.txt")
	path_nr <- paste0("variants/oct12_final/", mirna, "/", site, "_library_variants_oct12_final_noproofread.txt")

	variants <- read.table(path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
	variants_nr <- read.table(path_nr, sep="\t", header=TRUE, stringsAsFactors=FALSE)

	output_global <- rep(0, 16)
	names(output_global) <- GetKmerList(2)

	dinucs_all <- apply(variants, 1, CountDinucs, mirna=mirna, site=site)
	dinucs_all_nr <- apply(variants_nr, 1, CountDinucs, proofread=FALSE,
	                       mirna=mirna, site=site)

	dinucs_mean <- rowMeans(dinucs_all)
	dinucs_mean_nr <- rowMeans(dinucs_all_nr)

	dinucs_sd <- apply(dinucs_all, 1, sd)
	dinucs_sd_nr <- apply(dinucs_all_nr, 1, sd)
	dev.new(xpos=20, ypos=20, height=5, width=5)
	par(kPlotParameters)
	xmin <- 0
	xmax <- 0.15
	ymin <- 0
	ymax <- 0.15
	BlankPlot()
	AddLinearAxis(1, tick.space=0.01, label.space=0.05,
	              label="Variant dinucleotide frequencies (%)", percent=TRUE)
	AddLinearAxis(2, tick.space=0.01, label.space=0.05,
	              label="3'UTR dinucleotide frequencies (%)", percent=TRUE)
	abline(0, 1, lty=2)
	arrows(unlist(freqs), dinucs_mean - dinucs_sd, unlist(freqs), dinucs_mean + dinucs_sd,
	       length=error_height_final*par()$cex, lwd=par()$cex, angle=90, code=3)
	arrows(unlist(freqs), dinucs_mean_nr - dinucs_sd_nr, unlist(freqs), dinucs_mean_nr + dinucs_sd_nr,
	       length=error_height_final*par()$cex, lwd=par()$cex, angle=90, code=3, col="red")
	arrows(unlist(freqs) - freqs_norm_sd, dinucs_mean,
       unlist(freqs) + freqs_norm_sd, dinucs_mean,
       length=error_height_final*par()$cex, lwd=par()$cex, angle=90, code=3)
	points(unlist(freqs), dinucs_mean)

	points(unlist(freqs), dinucs_mean_nr, col="red")
	xy <- GetPlotFractionalCoords(0.05, 0.95)
	text(xy[1], xy[2], label=mirna, adj=0)
	xy <- GetPlotFractionalCoords(0.05, 0.90)
	text(xy[1], xy[2], label=site, adj=0)


	xy <- GetPlotFractionalCoords(0.4, 0.1)
	text(xy[1], xy[2], label="Selected variants", adj=c(0, 0))
	xy <- GetPlotFractionalCoords(0.95, 0.1)
	AddCorrelationToPlot(unlist(freqs), dinucs_mean, xy[1], xy[2], adj=c(1, 0),
	                     rsquared=TRUE)
	xy <- GetPlotFractionalCoords(0.4, 0.05)
	text(xy[1], xy[2], label="Control variants", adj=c(0, 0), col="red")
	xy <- GetPlotFractionalCoords(0.95, 0.05)
	AddCorrelationToPlot(unlist(freqs), dinucs_mean_nr, xy[1], xy[2],
	                     adj=c(1, 0), rsquared=TRUE, col="red")




}


# graphics.off()

# PlotDinucsFreqs("miR-1", "8mer")
# dev.copy2pdf(file="/lab/bartel1_ata/mcgeary/Presentations/seminars/181024\ Group\ Meeting/miR_1_8mer_dinucs.pdf")

GetAllDinucsFreqs <- function(mirna, sitelist="paperfinal", experiment="equilibrium") {
	if (mirna == "miR-1") {
		buffer <- TRUE
	} else {
		buffer <- FALSE
	}
	if (mirna == "miR-7") {
		mirna <- "miR-7-23nt"
		experiment <- "equilibrium2_nb"
	} else {
		experiment <- "equilibrium"
	}

	sites <- rownames(SubfunctionCall(SitesXCounts))
	if (mirna == "miR-7-23nt") {
		mirna <- "miR-7"
	}
	sites <- sites[-length(sites)]
	# print(sites)
	cors <- sapply(sites, function(site) {
		print(site)
		freqs_raw <- read.table("variants/mRNA_UTR_frequences.txt", header=TRUE, row.names=1)
		freqs <- freqs_raw[sort(rownames(freqs_raw)), , drop=FALSE]
		path <- paste0("variants/oct12_final/", mirna, "/", site, "_library_variants_oct12_final.txt")
		path_nr <- paste0("variants/oct12_final/", mirna, "/", site, "_library_variants_oct12_final_noproofread.txt")

		variants <- read.table(path, sep="\t", header=TRUE, stringsAsFactors=FALSE)
		variants_nr <- read.table(path_nr, sep="\t", header=TRUE, stringsAsFactors=FALSE)

		output_global <- rep(0, 16)
		names(output_global) <- GetKmerList(2)

		dinucs_all <- apply(variants, 1, CountDinucs, mirna=mirna, site=site)
		dinucs_all_nr <- apply(variants_nr, 1, CountDinucs, proofread=FALSE,
		                       mirna=mirna, site=site)

		dinucs_mean <- rowMeans(dinucs_all)
		dinucs_mean_nr <- rowMeans(dinucs_all_nr)

		dinucs_sd <- apply(dinucs_all, 1, sd)
		dinucs_sd_nr <- apply(dinucs_all_nr, 1, sd)

		# print(cor(dinucs_all, unlist(freqs))^2)
		return(c(cor(dinucs_mean, unlist(freqs))^2, cor(dinucs_mean_nr, unlist(freqs))^2))
	})
	print(cors)
	cors
}

# PlotDinucsFreqs("lsy-6", "8mer-bG7")
# dev.copy2pdf(file="/lab/bartel1_ata/mcgeary/Presentations/seminars/181024\ Group\ Meeting/lsy6_8merbG7_dinucs.pdf")

# PlotDinucsFreqs("miR-124", "CCGCCAC")
# dev.copy2pdf(file="/lab/bartel1_ata/mcgeary/Presentations/seminars/181024\ Group\ Meeting/mir124_CCGCCAC_dinucs.pdf")
# break


# cors_mir1 <- GetAllDinucsFreqs("miR-1")
# cors_let7a <- GetAllDinucsFreqs("let-7a")
# cors_mir155 <- GetAllDinucsFreqs("miR-155")
# cors_mir124 <- GetAllDinucsFreqs("miR-124")
# cors_lsy6 <- GetAllDinucsFreqs("lsy-6")
# cors_miR7 <- GetAllDinucsFreqs("miR-7")


cors_all <- cbind(cors_mir1, cors_let7a, cors_mir155, cors_mir124, cors_lsy6, cors_miR7)
# plot(ecdf(cors_mir1[1,]))

PlotAllCorrelations <- function() {
	dev.new(xpos=20, ypos=20, height=5, width=5)
	cors_all <- cors_all[, order(cors_all[1, ])]
	xmin <- 0
	xmax <- ncol(cors_all)
	ymin <- 0.90
	ymax <- 1
	par(kPlotParameters)
	BlankPlot()
	AddLinearAxis(1, tick.space=10, label.space=20, label="Sites")
	AddLinearAxis(2, tick.space=0.01, label.space=0.01, label=expression(italic(r)^2), percent=TRUE)
	lines(seq(ncol(cors_all)), cors_all[2, ], col="red")
	lines(seq(ncol(cors_all)), cors_all[1, ])
	xy <- GetPlotFractionalCoords(0.6, 0.10)
	text(xy[1], xy[2], labels="Selected variants", adj=0)
	xy <- GetPlotFractionalCoords(0.6, 0.05)
	text(xy[1], xy[2], labels="Control variants", col="red", adj=0)
	dev.copy2pdf(file="/lab/bartel1_ata/mcgeary/Presentations/seminars/181024\ Group\ Meeting/all_dinucs.pdf")

}

# PlotAllCorrelations()


PlotLikelihoods <- function() {
	dev.new(xpos=20, ypos=20, height=5, width=5)
	N_v <- 29992
	N <- 1e-15
	Na <- 6.022e23
	fracts <- 10^seq(-10, 1, length.out=100)

	xmin <- 1e-10
	xmax <- 1
	ymin <- 1e-10
	ymax <- 1
	par(kPlotParameters)
	BlankPlot(log='xy')
	AddLogAxis(1, label="Fraction of synthesis sampled")
	AddLogAxis(2, label="Probability of losing 0 variants")
	total_fracts <- Na*N*fracts
	pois_0 <- exp(-total_fracts)
	pois_greater_than_0 <- 1 - pois_0
	lose_none <- pois_greater_than_0^N_v
	lines(fracts, lose_none)
	dev.copy2pdf(file="/lab/bartel1_ata/mcgeary/Presentations/seminars/181024\ Group\ Meeting/library_dilution.pdf")

}

PlotError <- function() {
	dev.new(xpos=20, ypos=20, height=5, width=5)
	N_v <- 29992
	N <- 1e-15
	Na <- 6.022e23
	fracts <- 10^seq(-10, 1, length.out=100)

	xmin <- 1e-10
	xmax <- 1
	ymin <- 1e-5
	ymax <- 1
	par(kPlotParameters)
	BlankPlot(log='xy')
	AddLogAxis(1, label="Fraction of synthesis sampled")
	AddLogAxis(2, label="Error in variant fractional abundance (%)", percent=TRUE)
	total_fracts <- Na*N*fracts
	print(total_fracts)
	ind_min <- which((total_fracts - 1111.4)^2 == min((total_fracts - 1111.4)^2))
	print(ind_min)
	error <- 1/sqrt(total_fracts)
	# pois_greater_than_0 <- 1 - pois_0
	# lose_none <- pois_greater_than_0^N_v
	error_alt <- 4.1e-5/sqrt(fracts)
	lines(fracts, error)
	points(c(fracts[ind_min]), c(error[ind_min]), col="blue")
	# lines(fracts, error_alt, col="red")
	dev.copy2pdf(file="/lab/bartel1_ata/mcgeary/Presentations/seminars/181024\ Group\ Meeting/library_error.pdf")

}


PlotError()

