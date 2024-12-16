# graphics.off()
source("general/general.R")
source("general/GenericFigures.R")
library(mvmeta)
library(grid)
library(VennDiagram)
options(width=as.integer(200))
KConditions <- list(c("I", "0", "2,1", "2,2", "5,1", "5,2", "8,1", "8,2",
                      "15,1", "15,2", "30", "30", "90", "300", "900", "2700",
                      "7200,1", "7200,2","Equil", "0-", "Equil-"),
					c("I", "0", "2,1", "2,2", "5,1", "5,2", "8,1", "8,2",
					  "15,1", "15,2", "30,1", "30,2", "90", "300", "900",
					  "2700", "7200", "Equil", "0-", "Equil-"),
					c("I", "0", "2,1", "2,2", "5,1", "5,2", "8,1", "8,2",
					  "15,1", "15,2", "30,1", "30,2", "90", "300", "900",
					  "2700", "7200", "Equil", "0-", "Equil-"),
					c("I", "0", "2,1", "2,2", "5,1", "5,2", "8,1", "8,2",
					  "15,1", "15,2", "30,1", "30,2", "90", "300", "900",
					  "2700", "7200", "Equil", "0-", "Equil-"),
				    c("I", "0", "2", "5", "8", "15", "30", "90", "300", "900",
				      "2700", "7200","Equil", "0-", "Equil-"))
names(KConditions) <- c("let-7a", "miR-1", "miR-124", "lsy-6", "noreps")

ReadInput <- function(condition, mirna="miR-1", experiment="kinetics",
                      n_constant=5, len_k=6, n_ex=0, pulse=TRUE) {
	print(n_ex)
	analysis_type <- "kmers_cutoff_final"
	ext <- sprintf("_%s_k%s", n_constant, len_k)
	if (n_ex != 0) {
		ext <- sprintf("%s_ex%s", ext, n_ex)
	}
	print(ext)
	file_path <- SubfunctionCall(GetAnalysisPath)
	print(file_path)
	output <- read.table(file_path, header=TRUE, row.names=1)
	if (pulse) {
		out <- data.frame(output$counts_p)
	} else {
		out <- data.frame(output$counts_c)
	}
	rownames(out) <- rownames(output)
	out
}

RecursiveFlankLoop <- function(n) {
	if (n < 0) {
		print("Only can use integers >= 0")
	} else if (n == 0) {
		c("")
	} else {
		c(sapply(RecursiveFlankLoop(n - 1), paste0, y=c("A", "C", "G", "T")))
	}
}


GetKmerColors <- function(mirna, len_k, alpha=0.1) {
	kmer_cols <- rep(ConvertRColortoRGB("gray", alpha=alpha), 4^len_k)
	names(kmer_cols) <- RecursiveFlankLoop(len_k)
	for (site in rev(kSeedSites)) {
		message(site)
		seq <- SubfunctionCall(GetSiteSeq)
		message(seq)
		if (nchar(seq) <= len_k) {
			inds <- grep(seq, names(kmer_cols))
			kmer_cols[inds] <- kSiteColors[site]
		}
	}
	kmer_cols
}

MakePulseAndChaseKmerTables <- function(mirna, len_k, experiment="kinetics",
                                        analysis_type="kmers_cutoff_final",
                                        n_constant=5, n_ex=0, 
                                        combined_reps=TRUE) {
	print(n_ex)
	conditions <- KConditions[[mirna]]
	print(conditions)
	nrow <- 4^len_k
	ncol <- length(conditions)
	dimnames <- list(RecursiveFlankLoop(len_k), conditions)	
	pulse_counts <- matrix(NaN, nrow=nrow, ncol=ncol, dimnames=dimnames)
	chase_counts <- pulse_counts
	print(n_ex)
	for (condition in conditions) {
		print(n_ex)
		col_p <- SubfunctionCall(ReadInput)
		col_c <- SubfunctionCall(ReadInput, pulse=FALSE)
		colnames(col_p) <- condition
		colnames(col_c) <- condition
		pulse_counts[, condition] <- col_p[, 1]
		chase_counts[, condition] <- col_c[, 1]
		print(pulse_counts["CTTCCT", ])
	}
	if (combined_reps) {
		rep_cols <- grep(",", conditions, value=TRUE)
		reps_split <- matrix(unlist(strsplit(rep_cols, split=",")),
		                     ncol=2, byrow=TRUE)
		print(reps_split)
		conditions_noreps <- KConditions[["noreps"]]
		dimnames[[2]] <- conditions_noreps
		ncol <- length(conditions_noreps)
		pulse_counts_new <- matrix(0, nrow=nrow, ncol=ncol, dimnames=dimnames)
		chase_counts_new <- pulse_counts_new
		for (condition in conditions) {
			if (condition %in% rep_cols) {
				condition_norep <- unlist(strsplit(condition, ","))[1]
			} else {
				condition_norep <- condition
			}
			pulse_counts_new[, condition_norep] <- (
			    pulse_counts_new[, condition_norep] + pulse_counts[, condition]
		    )
			chase_counts_new[, condition_norep] <- (
			    chase_counts_new[, condition_norep] + chase_counts[, condition]
            )
		}
	pulse_counts <- pulse_counts_new
	chase_counts <- chase_counts_new
	}

	return(list(pulse_counts, chase_counts))
}

PlotKmerDists <- function(mirna, len_k, experiment="kinetics",
                          analysis_type="kmers_cutoff_final", n_constant=5,
                          n_ex=0, make.pdf=False) {
	kmer_cols <- SubfunctionCall(GetKmerColors)

	counts <- SubfunctionCall(MakePulseAndChaseKmerTables)

	pulse_counts <- counts[[1]]
	chase_counts <- counts[[2]]

	pulse_norm <- t(t(pulse_counts) / colSums(pulse_counts + chase_counts))
	chase_norm <- t(t(chase_counts) / colSums(pulse_counts + chase_counts))

	plot_name <- sprintf("%s_k%s_enrichments", mirna, len_k)
	if (make.pdf) {
		KineticsScrapFigureSaveFile(pdf.plot=plot_name, height=20, width=20)	
	} else {
		dev.new(xpos=20, ypos=20, height=20, width=20)
		par(kPlotParameters)
	}
	par(mfrow=c(4, 4))
	xmin <- 1e-2
	xmax <- 1e2
	sapply(colnames(pulse_norm)[-1], function(condition) {
		print(i)
		if (condition %in% c("0", "0-")) {
			ymin <- 1e-5
			ymax <- 0.1
		} else {
			ymin <- xmin
			ymax <- xmax
		}
		BlankPlot(log='xy')
		points(pulse_norm[, condition]/Norm(pulse_norm[, "I"]),
		     chase_norm[, condition]/Norm(chase_norm[, "I"]), col=kmer_cols)
		AddLogAxis(1, label="Pulse enrichment")
		AddLogAxis(2, label="Chase enrichment")
		title(main=condition)
		abline(0, 1, lty=2)
	})
	dev.off()
}

GetLaggingKmerDists <- function(mirna, len_k, experiment="kinetics",
                                analysis_type="kmers_cutoff_final",
                                n_constant=5, n_ex=0, condition="7200",
                                pdf.plot=FALSE, make_plots=FALSE) {
	kmer_cols <- SubfunctionCall(GetKmerColors, alpha=0.01)
	kmer_cols[1:length(kmer_cols)] <- rep(ConvertRColortoRGB("gray", alpha=0.5), length(kmer_cols))

	print(head(kmer_cols))
	inds_kmer <- grep("CTTCC", names(kmer_cols), value=TRUE)
	kmer_cols[inds_kmer] <- "blue"
	counts <- SubfunctionCall(MakePulseAndChaseKmerTables)

	pulse_counts <- counts[[1]]
	chase_counts <- counts[[2]]

	print(head(pulse_counts))
	pulse_norm <- t(t(pulse_counts) / colSums(pulse_counts))
	chase_norm <- t(t(chase_counts) / colSums(chase_counts))


	print(head(pulse_norm))
	# xmin <- 0.4
	# xmax <- 50
	# ymin <- xmin
	# ymax <- xmax
	# graphics.off()
	# dev.new(xpos=20, ypos=20, height=10, width=10)
	# par(kPlotParameters)
	# par(mfrow=c(2, 2))
	# BlankPlot(log='xy')
	# points(pulse_norm[, "7200"]/pulse_norm[, "I"],
	#      chase_norm[, "7200"]/chase_norm[, "I"], pch=1, col=kmer_cols)
	# points(pulse_norm[inds_kmer, "7200"]/pulse_norm[inds_kmer, "I"],
	#      chase_norm[inds_kmer, "7200"]/chase_norm[inds_kmer, "I"], pch=1, col="blue")


	# AddLogAxis(1, label="pulse 2h")
	# AddLogAxis(2, label="chase 2h")
	# xy <- GetPlotFractionalCoords(fx=0.05, fy=0.9, log='xy')
	# text(xy[1], xy[2], label=mirna, adj=0)

	# BlankPlot(log='xy')
	# points(pulse_norm[, "Equil"]/pulse_norm[, "I"],
	#      chase_norm[, "Equil"]/chase_norm[, "I"], pch=1, col=kmer_cols)
	# points(pulse_norm[inds_kmer, "Equil"]/pulse_norm[inds_kmer, "I"],
	#      chase_norm[inds_kmer, "Equil"]/chase_norm[inds_kmer, "I"], pch=1, col="blue")
	# AddLogAxis(1, label="pulse Equil")
	# AddLogAxis(2, label="chase Equil")

	# # identify(pulse_norm[, "Equil"]/pulse_norm[, "I"],
	# #      chase_norm[, "Equil"]/chase_norm[, "I"], labels=rownames(pulse_norm))




	# BlankPlot(log='xy')
	# points(pulse_norm[, "7200"]/pulse_norm[, "I"],
	#      pulse_norm[, "Equil"]/pulse_norm[, "I"], pch=1, col=kmer_cols)
	# points(pulse_norm[inds_kmer, "7200"]/pulse_norm[inds_kmer, "I"],
	#      pulse_norm[inds_kmer, "Equil"]/pulse_norm[inds_kmer, "I"], pch=1, col="blue")
	# AddLogAxis(1, label="pulse 2h")
	# AddLogAxis(2, label="pulse Equil")

	# BlankPlot(log='xy')
	# points(chase_norm[, "7200"]/chase_norm[, "I"],
	#      chase_norm[, "Equil"]/chase_norm[, "I"], pch=1, col=kmer_cols)
	# points(chase_norm[inds_kmer, "7200"]/chase_norm[inds_kmer, "I"],
	#      chase_norm[inds_kmer, "Equil"]/chase_norm[inds_kmer, "I"], pch=1, col="blue")
	# AddLogAxis(1, label="chase 2h")
	# AddLogAxis(2, label="chase Equil")


	# break
	pc_ratio_2h <- chase_norm[, condition]/Norm(chase_norm[, "I"])/(
	    pulse_norm[, condition]/Norm(pulse_norm[, "I"])
    )


	mean_pc_ratio_2h <- GeoMean(pc_ratio_2h)
    sd_pc_ratio_2h <- exp(sd(log(pc_ratio_2h), na.rm=TRUE))

    lower_CI <- mean_pc_ratio_2h/(sd_pc_ratio_2h^2)
    upper_CI <- mean_pc_ratio_2h*(sd_pc_ratio_2h^2)


	p_2he_ratio <- pulse_norm[, "Equil"]/pulse_norm[, condition]

	mean_p_2he_ratio <- GeoMean(p_2he_ratio)
	sd_p_2he_ratio <- exp(sd(log(p_2he_ratio), na.rm=TRUE))
	lower_CI_2 <- mean_p_2he_ratio/(sd_p_2he_ratio^2)
	upper_CI_2 <- mean_p_2he_ratio*(sd_p_2he_ratio^2)


	inds_ratio <- which(pc_ratio_2h <= lower_CI & p_2he_ratio <= lower_CI_2)


	kmer_cols[inds_ratio] <- "green"
	if (make_plots) {
		SubfunctionCall(KineticsScrapFigureSaveFile, height=8, width=8)
		par(mfrow=c(2, 2))
		xmin <- 1e-2
		xmax <- 1e2
		ymin <- xmin
		ymax <- xmax
		BlankPlot(log='xy')
		points(pulse_norm[, condition]/Norm(pulse_norm[, "I"]),
		     chase_norm[, condition]/Norm(chase_norm[, "I"]), col=kmer_cols)
		AddLogAxis(1, label="Pulse enrichment")
		AddLogAxis(2, label="Chase enrichment")

		xline <- exp(log(seq(xmin, xmax, length.out=100)))
		lines(xline, xline * median(pc_ratio_2h), lty=2, col="red")
		lines(xline, xline * GeoMean(pc_ratio_2h), lty=2, col="blue")
		lines(xline, xline * mean(pc_ratio_2h), lty=2, col="green")

		xmin <- 3e0
		xmax <- 3e1
		ymin <- 0
		ymax <- 1

		ecdf_1 <- ecdf(pc_ratio_2h)
		BlankPlot(log='x')
		lines(ecdf_1)
		AddLogAxis(1, label="Pulse/chase; 2 h")
		AddLinearAxis(2, label="CDF")
		segments(c(lower_CI, mean_pc_ratio_2h, upper_CI), c(0, 0, 0),
		         c(lower_CI, mean_pc_ratio_2h, upper_CI), c(1, 1, 1),
		         lty=c(2, 1, 2))
		title(main=condition)



		xmin <- 1e-2
		xmax <- 1e2
		ymin <- xmin
		ymax <- xmax
		BlankPlot(log='xy')
		points(pulse_norm[, condition]/Norm(pulse_norm[, "I"]),
		     pulse_norm[, "Equil"]/Norm(pulse_norm[, "I"]), col=kmer_cols)
		AddLogAxis(1, label="Pulse enrichment")
		AddLogAxis(2, label="Chase enrichment")
		title(main=condition)

		lines(xline, xline * median(p_2he_ratio), lty=2, col="red")
		lines(xline, xline * GeoMean(p_2he_ratio), lty=2, col="blue")
		lines(xline, xline * mean(p_2he_ratio), lty=2, col="green")

		xmin <- 3e-1
		xmax <- 3e0
		ymin <- 0
		ymax <- 1
		BlankPlot(log='x')

		ecdf_2 <- ecdf(p_2he_ratio)
		lines(ecdf_2)
		AddLogAxis(1, label="Pulse; 2 h / equilibrium")
		AddLinearAxis(2, label="CDF")
		segments(c(lower_CI_2, mean_p_2he_ratio, upper_CI_2), c(0, 0, 0),
		         c(lower_CI_2, mean_p_2he_ratio, upper_CI_2), c(1, 1, 1),
		         lty=c(2, 1, 2))
		if (class(pdf.plot) == "character") {
			dev.off()	
		}
	}

	out <- cbind(pc_ratio_2h, p_2he_ratio)
	out <- t(t(out) - colMeans(out, na.rm=TRUE))
	out <- t(t(out)/apply(out, 2, sd, na.rm=TRUE))
	out
}

PlotLaggingKmers <- function(len_k, n_ex=0) {
	k_l7 <- GetLaggingKmerDists("let-7a", len_k, n_ex=n_ex)
	k_m1 <- GetLaggingKmerDists("miR-1", len_k, n_ex=n_ex)
	k_m124 <- GetLaggingKmerDists("miR-124", len_k, n_ex=n_ex)
	k_l6 <- GetLaggingKmerDists("lsy-6", len_k, n_ex=n_ex)

	# print(k_l7)
	k_list <- list(k_l7, k_m1, k_m124, k_l6)
	# k_list <- list(k_m1, k_m124)
	names(k_list) <- c("let-7a", "miR-1", "miR-124", "lsy-6")
	# names(k_list) <- c("miR-1", "miR-124")
	# k_all <- unique(c(rownames(k_l7), rownames(k_m1), rownames(k_m124),
	# 				  rownames(k_l6)))



	# in_l7 <- (k_all %in% rownames(k_l7))
	# in_m1 <- (k_all %in% rownames(k_m1))
	# in_m124 <- (k_all %in% rownames(k_m124))
	# in_l6 <- (k_all %in% rownames(k_l6))


	# in_alls <- cbind(in_l7, in_m1, in_m124, in_l6)

	pc_ratio_all <- cbind(k_l7[, 1],
	                  k_m1[, 1],
	                  k_m124[, 1],
	                  k_l6[, 1, drop=FALSE])

	# pc_ratio_all <- cbind(k_m1[, 1],
	#                   	  k_m124[, 1])

	equil_2h_ratio_all <- cbind(k_l7[, 2],
	                            k_m1[, 2],
				                k_m124[, 2],
				                k_l6[, 2])

	colnames(pc_ratio_all) <- names(k_list)

	colnames(equil_2h_ratio_all) <- colnames(pc_ratio_all)



	pc_ratio_all <- cbind(pc_ratio_all, rowMeans(pc_ratio_all))
	equil_2h_ratio_all <- cbind(equil_2h_ratio_all, rowMeans(equil_2h_ratio_all))

	print(head(pc_ratio_all))
	print(head(equil_2h_ratio_all))
	

	dev.new(xpos=20, ypos=20, height=5, width=5)
	xmin <- -6
	xmax <- 4
	ymin <- -6
	ymax <- 4
	par(kPlotParameters)
	BlankPlot()
	print(dim(pc_ratio_all))
	segments(xmin, 0, xmax, 0)
	segments(0, ymin, 0, ymax)
	abline(0, 1, lty=2)
	points(pc_ratio_all[, 5], equil_2h_ratio_all[, 5], col=ConvertRColortoRGB("black", alpha=0.2))
	AddLinearAxis(1, tick.space=1, label.space=2, label="Pulse/Chase; 2 h")
	AddLinearAxis(2, tick.space=1, label.space=2, label="2 h/Equilibrium; Pulse reads")
	# identify(pc_ratio_all[, 5], equil_2h_ratio_all[, 5], labels=rownames(pc_ratio_all))
	dev.copy2pdf(
  		file=paste0("/lab/bartel1_ata/mcgeary/Presentations/seminars/181024\ Group\ Meeting/",
                  "Lagging_Kmer_distribution_unlabeled.pdf")
	)


	print(head(pc_ratio_all))
	out <- cbind(pc_ratio_all[, 5], equil_2h_ratio_all[, 5],
	             rowMeans(cbind(pc_ratio_all[, 5], equil_2h_ratio_all[, 5])))
	colnames(out) <- c("pc_2h_ratio", "p_2he_ratio", "mean")
	out <- out[order(out[, 3]), ]
	print(head(out))
	write.table(out, file="SolveForOffRates/LaggingKmers.txt", quote=FALSE,
	            sep="\t")
	out
}

# k_l7 <- GetLaggingKmerDists("let-7a", 6)
# k_m1 <- GetLaggingKmerDists("miR-1", 6)
# k_m124 <- GetLaggingKmerDists("miR-124", 6)
# k_l6 <- GetLaggingKmerDists("lsy-6", 6)

out <- PlotLaggingKmers(6)
m_out <- mean(out[,3])
sd <- sd(out[,3])

len_2sd <- length(which(out[, 3] <= m_out - 2*sd))
len_sd  <- length(which(out[, 3] <= m_out - sd))
len_mean <- length(which(out[, 3] <= m_out))

print(len_2sd)
print(len_sd)
print(len_mean)

break
ex0 <- ReadInput("miR-1", condition="7200", n_constant=5, experiment="kinetics",
                 len_k=6, pulse=TRUE)

ex1 <- ReadInput("miR-1", condition="7200", n_constant=5, experiment="kinetics",
                 len_k=6, pulse=TRUE, n_ex="CTTCCT")

ex2 <- ReadInput("miR-1", condition="7200", n_constant=5, experiment="kinetics",
                 len_k=6, pulse=TRUE, n_ex="CTTCCT,CCGCCA")


ex3 <- ReadInput("miR-1", condition="7200", n_constant=5, experiment="kinetics",
                 len_k=6, pulse=TRUE, n_ex="CTTCCT,CCGCCA,GCTTCC")

print(colSums(ex0))
print(colSums(ex1))
print(colSums(ex2))
print(colSums(ex3))



break

mirna <- "let-7a"
vec_t <- c("0", "2", "5", "8", "15", "30", "90", "300",
           "900", "2700", "7200", "Equil", "0-", "Equil-")

dups <- list('let-7a'=c("2", "5", "8", "15", "7200"),
             'miR-1' =c("2", "5", "8", "15", "30"))

PlotPositionalTimeCourse <- function(mirna, n_constant, len_k) {
	mat_I <- GetPositionalKmers(mirna, "kinetics", "I", 5, 6)
	mat_I_p <- (mat_I[[1]] + 1)/sum(mat_I[[1]] + mat_I[[2]] + 1 + 1)
	mat_I_c <- (mat_I[[2]] + 1)/sum(mat_I[[1]] + mat_I[[2]] + 1 + 1)
	graphics.off()
	dev.new(xpos=20, ypos=20, height=10, width=10)
	par(mfrow=c(4, 4))
	sapply(vec_t, function(condition) {
		print(condition)
		if (condition %in% dups[[mirna]]) {
			conditions <- sprintf("%s,%s", condition, c("1", "2"))
			print(conditions)
			mat_cond_1 <- GetPositionalKmers(mirna, "kinetics", conditions[1], 5, 6)
			mat_cond_2 <- GetPositionalKmers(mirna, "kinetics", conditions[2], 5, 6)
			mat_cond <- list(mat_cond_1[[1]] + mat_cond_2[[1]],
			                 mat_cond_1[[2]] + mat_cond_2[[2]])
		} else {
			print(condition)
			mat_cond <- GetPositionalKmers(mirna, "kinetics", condition, 5, 6)	
		}
		mat_cond_p <- mat_cond[[1]]/sum(mat_cond[[1]] + mat_cond[[2]])
		mat_cond_c <- mat_cond[[2]]/sum(mat_cond[[1]] + mat_cond[[2]])

		mat_R_cond_p <- mat_cond_p["CTTCCT",]/mat_I_p["CTTCCT",]
		mat_R_cond_c <- mat_cond_c["CTTCCT",]/mat_I_c["CTTCCT",]
		mat_R_cond_p <- mat_R_cond_p / mean(as.numeric(mat_R_cond_p))
		mat_R_cond_c <- mat_R_cond_c / mean(as.numeric(mat_R_cond_c))
		print(mat_R_cond_p)
		print(mat_R_cond_c)
		ymax <- ceiling(max(c(max(mat_R_cond_p), max(mat_R_cond_c))))
		plot(  1:ncol(mat_I_p), mat_R_cond_p, pch=20, type="o", ylim=c(0, 10))
		points(1:ncol(mat_I_p), mat_R_cond_c, pch=20, type="o", col="blue")
		title(main=condition)
	})
}

mirna <- "miR-1"
n_constant <- 5
len_k <- 6


kmers_full <- MakePulseAndChaseKmerTables(mirna, len_k=len_k,
                                          n_constant=n_constant, exp="kinetics")

kmers_full_ex1 <- MakePulseAndChaseKmerTables(mirna, len_k=len_k,
                                          n_constant=n_constant, exp="kinetics",
                                          n_ex="CTTCCT")



kmers_p <- kmers_full[[1]]
kmers_c <- kmers_full[[2]]

kmers_ex1_p <- kmers_full_ex1[[1]]
kmers_ex1_c <- kmers_full_ex1[[2]]


kmers_p_norm <- t(t(kmers_p) / (colSums(kmers_p) + colSums(kmers_c)))
kmers_c_norm <- t(t(kmers_c) / (colSums(kmers_p) + colSums(kmers_c)))

kmers_ex1_p_norm <- t(t(kmers_ex1_p) / (colSums(kmers_ex1_p) + colSums(kmers_ex1_c)))
kmers_ex1_c_norm <- t(t(kmers_ex1_c) / (colSums(kmers_ex1_p) + colSums(kmers_ex1_c)))

print(head(kmers_p_norm))

graphics.off()
dev.new(xpos=20, ypos=20, height=5, width=5)

plot(c(1, 2, 5, 8, 15, 30, 90, 300, 900, 2700, 7200, 10000),
     (kmers_p_norm["CTTCCT", 2:13]/(kmers_p_norm["CTTCCT", 1])), log='x',
     type="o")
points(c(1, 2, 5, 8, 15, 30, 90, 300, 900, 2700, 7200, 10000),
     (kmers_c_norm["CTTCCT", 2:13]/(kmers_c_norm["CTTCCT", 1])), type="o",
     col="red")


dev.new(xpos=20, ypos=520, height=5, width=5)
plot(c(1, 2, 5, 8, 15, 30, 90, 300, 900, 2700, 7200, 10000),
     (kmers_p_norm["CATTCC", 2:13]/(kmers_p_norm["CATTCC", 1])), log='x',
     type="o", ylim=c(1, 100))
points(c(1, 2, 5, 8, 15, 30, 90, 300, 900, 2700, 7200, 10000),
     (kmers_c_norm["CATTCC", 2:13]/(kmers_c_norm["CATTCC", 1])), type="o",
     col="red")
# points(c(1, 2, 5, 8, 15, 30, 90, 300, 900, 2700, 7200, 10000),
#      (kmers_ex1_p_norm["CATTCC", 2:13]/(kmers_ex1_p_norm["CATTCC", 1])),
#      type="o", lty=2)

points(c(1, 2, 5, 8, 15, 30, 90, 300, 900, 2700, 7200, 10000),
     (kmers_ex1_c_norm["CATTCC", 2:13]/(kmers_ex1_c_norm["CATTCC", 1])), type="o",
     col="red", lty=2)



dev.new(xpos=520, ypos=20, height=5, width=5)

plot(c(1, 2, 5, 8, 15, 30, 90, 300, 900, 2700, 7200, 10000),
     (colSums(kmers_p_norm[, 2:13])/sum(kmers_p_norm[, 1])), log='x',
     type="o")

points(c(1, 2, 5, 8, 15, 30, 90, 300, 900, 2700, 7200, 10000),
     (colSums(kmers_c_norm[, 2:13])/sum(kmers_c_norm[, 1])), type="o",
     col="red")

dev.new(xpos=520, ypos=520, height=5, width=5)


kmers_exclude <- c("CTTCCT")
inds_exclude <- sapply(kmers_exclude, function(kmer) {
	which(rownames(kmers_c_norm) == kmer)
})


plot(c(1, 2, 5, 8, 15, 30, 90, 300, 900, 2700, 7200, 10000),
     (colSums(kmers_p_norm[-inds_exclude, 2:13])/sum(kmers_p_norm[-inds_exclude, 1])), log='x',
     type="o")

points(c(1, 2, 5, 8, 15, 30, 90, 300, 900, 2700, 7200, 10000),
     (colSums(kmers_c_norm[-inds_exclude, 2:13])/sum(kmers_c_norm[-inds_exclude, 1])), type="o",
     col="red")

points(c(1, 2, 5, 8, 15, 30, 90, 300, 900, 2700, 7200, 10000),
     (colSums(kmers_ex1_p_norm[, 2:13])/sum(kmers_ex1_p_norm[, 1])), log='x',
     type="o", col="gray", lty=2, lwd=2)

points(c(1, 2, 5, 8, 15, 30, 90, 300, 900, 2700, 7200, 10000),
     (colSums(kmers_ex1_c_norm[, 2:13])/sum(kmers_ex1_c_norm[, 1])), type="o",
     col="purple", lty=2, lwd=2)



# plot(c(1, 2, 5, 8, 15, 30, 90, 300, 900, 2700, 7200, 10000),
#      (kmers_p_norm["CTTCCT", 2:13]/(kmers_p_norm["CTTCCT", 1])), log='xy',
#      type="o")

# points(c(1, 2, 5, 8, 15, 30, 90, 300, 900, 2700, 7200, 10000),
#      (kmers_c_norm["CTTCCT", 2:13]/(kmers_c_norm["CTTCCT", 1])), type="o")



break

kmers_2h <- GetPositionalKmers(mirna, "kinetics", condition="7200",
                                 n_constant=n_constant, kmer_length=len_k)

kmers_2h_p <- kmers_2h[[1]]
kmers_2h_c <- kmers_2h[[2]]

print(kmers_p["CTTCCT", ])

print(sum(kmers_2h_p["CTTCCT", ]))
print(kmers_c["CTTCCT", ])
print(sum(kmers_2h_c["CTTCCT", ]))


break

PlotPositionalTimeCourse("miR-1", 5, 6)


break
mat_I_e <- GetPositionalKmers("let-7a", "equilibrium", "I", 5, 6)

mat_I <- GetPositionalKmers("let-7a", "kinetics", "I", 5, 6)
mat_2h.1 <- GetPositionalKmers("let-7a", "kinetics", "7200,1", 5, 6)
mat_2h.2 <- GetPositionalKmers("let-7a", "kinetics", "7200,2", 5, 6)
mat_2h <- list(mat_2h.1[[1]] + mat_2h.2[[1]], mat_2h.1[[2]] + mat_2h.2[[2]])
mat_eq <- GetPositionalKmers("let-7a", "kinetics", "Equil", 5, 6)



mat_I_e <- Norm(mat_I_e + 1)
mat_I_p <- (mat_I[[1]] + 1)/sum(mat_I[[1]] + mat_I[[2]] + 1 + 1)
mat_I_c <- (mat_I[[2]] + 1)/sum(mat_I[[1]] + mat_I[[2]] + 1 + 1)
mat_2h_p <- mat_2h[[1]]/sum(mat_2h[[1]] + mat_2h[[2]])
mat_2h_c <- mat_2h[[2]]/sum(mat_2h[[1]] + mat_2h[[2]])
mat_eq_p <- mat_eq[[1]]/sum(mat_eq[[1]] + mat_eq[[2]])
mat_eq_c <- mat_eq[[2]]/sum(mat_eq[[1]] + mat_eq[[2]])




dev.new(xpos=20, ypos=20, height=5, width=10)
par(mfrow=c(1, 2))
plot(1:ncol(mat_I_p), mat_I_p["CTTCCT",], pch=20, type="o", log='y')
points(1:ncol(mat_I_p), mat_I_c["CTTCCT",], pch=20, type="o", col="blue")
points(1:ncol(mat_I_p), mat_I_e["CTTCCT",], pch=20, type="o", col="green")

mat_R_2h_p <- mat_2h_p["CTTCCT",]/mat_I_p["CTTCCT",]
mat_R_2h_c <- mat_2h_c["CTTCCT",]/mat_I_c["CTTCCT",]

mat_R_eq_p <- mat_eq_p["CTTCCT",]/mat_I_p["CTTCCT",]
mat_R_eq_c <- mat_eq_c["CTTCCT",]/mat_I_c["CTTCCT",]



plot(  1:ncol(mat_I_p), mat_R_2h_p["CTTCCT",], pch=20, type="o", log='y', ylim=c(0.1, 10))
points(1:ncol(mat_I_p), mat_R_2h_c["CTTCCT",], pch=20, type="o", col="blue")
points(1:ncol(mat_I_p), mat_R_eq_p["CTTCCT",], pch=20, type="o", col="red")
points(1:ncol(mat_I_p), mat_R_eq_c["CTTCCT",], pch=20, type="o", col="green")





break

dev.new(xpos=20, ypos=20, height=8, width=8)
par(kPlotParameters)
par(mfrow=c(3, 3))





sapply(names(k_list)[1:3], function(mirna1) {
	k_mir1 <- k_list[[mirna1]]
	new_names <- setdiff(names(k_list), mirna1)
	print(mirna1)
	sapply(new_names, function(mirna2) {
		print(mirna2)
		k_mir2 <- k_list[[mirna2]]
		xmin <- 0
		xmax <- 8
		ymin <- xmin
		ymax <- xmax

		rows_use <- intersect(rownames(k_mir1), rownames(k_mir2))
		print(rows_use)
		print(k_shared)
		SubfunctionCall(BlankPlot)
		colors <- rep("gray", length(rows_use))
		names(colors) <- rows_use
		colors[k_shared] <- "blue"
		points(k_mir1[rows_use, 1], k_mir2[rows_use, 1], col=colors)
		# identify(k_mir1[rows_use, 1], k_mir2[rows_use, 1], labels=rows_use)
		AddLinearAxis(1, tick.space=1, label.space=2, label=mirna1)
		AddLinearAxis(2, tick.space=1, label.space=2, label=mirna2)
		xy <- GetPlotFractionalCoords(0.1, 0.95)
		text(xy[1], xy[2], mirna)
	})
})


break

dev.new(xpos=20, ypos=20, height=8, width=8)
par(kPlotParameters)
par(mfrow=c(2, 2))


sapply(names(k_list), function(mirna) {
	k_mir <- k_list[[mirna]]
	xmin <- 0
	xmax <- 8
	ymin <- 0
	ymax <- 1
	SubfunctionCall(BlankPlot)
	colors <- rep("gray", nrow(k_mir))
	names(colors) <- rownames(k_mir)
	colors[k_shared] <- "blue"
	points(k_mir[, 1], k_mir[, 2], col=colors)
	AddLinearAxis(1, tick.space=1, label.space=2, label="pulse/chase ratio at 2h")
	AddLinearAxis(2, tick.space=0.1, label.space=0.2, label="pulse 2h/equil ratio")
	xy <- GetPlotFractionalCoords(0.1, 0.95)
	text(xy[1], xy[2], mirna)
})



break

# a <- vennCounts(in_alls)



venn_temp <- venn.diagram(list(rownames(k_l7), rownames(k_m1), rownames(k_m124),
                               rownames(k_l6)),
						  filename=NULL,
						  category=c("let-7a", "miR-1", "miR-124", "lsy-5"))
grid.draw(venn_temp)

break


dev.new(xpos=10, ypos=20, height=10, width=10)
par(mfrow=c(2, 2))

plot(pulse_norm[, 15]/Norm(pulse_norm[, 1]),
     chase_norm[, 15]/Norm(chase_norm[, 1]), log='xy',
     xlim=c(0.03, 100), ylim=c(0.03, 100), col=kmer_cols)
title(main=colnames(pulse_norm)[15])
abline(0, 1, lty=2)

# identify(pulse_norm[, 15]/Norm(pulse_norm[, 1]),
#      	 chase_norm[, 15]/Norm(chase_norm[, 1]), labels=rownames(pulse_norm))

# dev.new(xpos=610, ypos=20, height=7, width=7)
plot(pulse_norm[, 16]/Norm(pulse_norm[, 1]),
     chase_norm[, 16]/Norm(chase_norm[, 1]), log='xy',
     xlim=c(0.03, 100), ylim=c(0.03, 100), col=kmer_cols)
title(main=colnames(pulse_norm)[16])
abline(0, 1, lty=2)

# identify(pulse_norm[, 16]/Norm(pulse_norm[, 1]),
#      	 chase_norm[, 16]/Norm(chase_norm[, 1]), labels=rownames(pulse_norm))

plot(pulse_norm[, 15]/Norm(pulse_norm[, 1]),
     pulse_norm[, 16]/Norm(pulse_norm[, 1]), log='xy',
     xlim=c(0.03, 100), ylim=c(0.03, 100), col=kmer_cols)
abline(0, 1, lty=2)

plot(chase_norm[, 15]/Norm(chase_norm[, 1]),
     chase_norm[, 16]/Norm(chase_norm[, 1]), log='xy',
     xlim=c(0.03, 100), ylim=c(0.03, 100), col=kmer_cols)
abline(0, 1, lty=2)


dev.new(xpos=720, ypos=20, height=5, width=5)
plot(pulse_norm[, 15]/pulse_norm[, 16],
     pulse_norm[, 15]/Norm(pulse_norm[, 1])*20/
     (chase_norm[, 15]/Norm(chase_norm[, 1])),
      log='xy', xlim=c(0.3, 10), ylim=c(0.3, 10), col=kmer_cols)
abline(0, 1, lty=2)

deviate_1 <- pulse_norm[, 15]/pulse_norm[, 16]
deviate_2 <- (pulse_norm[, 15]/Norm(pulse_norm[, 1])*20/
			  (chase_norm[, 15]/Norm(chase_norm[, 1])))

deviate_both <- cbind(deviate_1, deviate_2)
rownames(deviate_both) <- rownames(chase_norm)

deviate_both_order1 <- deviate_both[order(-deviate_both[, 1]), ]
deviate_both_order2 <- deviate_both[order(-deviate_both[, 2]), ]

print(head(deviate_both_order1, 20))
print(head(deviate_both_order2, 20))
break



# dev.new(xpos=720, ypos=20, height=10, width=10)
# par(mfrow=c(4, 4))
# sapply(2:16, function(i) {
# 	plot(Norm(chase_norm[, 2])/Norm(chase_norm[, 1]),
# 	     Norm(chase_norm[, i])/Norm(chase_norm[, 1]), log='xy',
# 	     xlim=c(0.1, 100), ylim=c(0.1, 100), col=kmer_cols)
# 	title(main=colnames(chase_norm)[i])
# 	abline(0, 1, lty=2)

# })






