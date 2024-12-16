graphics.off()
source("general/general.R")
source("general/GenericFigures.R")





base_path <- file.path("lab", "bartel1_ata", "nbisaria", "AgoRBNS", "figures",
                       "kds_loopplots_longloop")

mirna_path <- c(`let-7a`="let-7a_mix2_rep", `miR-1`="miR-1_mix2",
                 `miR-155`="miR-155_mix2_3_3")

GetKdMatrix <- function(mirna, register, kmer) {

	full_path <- sprintf("/%s/%s/matdfchange%smer%s_reg%s.txt", base_path,
	                     mirna_path[mirna], kmer, mirna_path[mirna], register)
	print(full_path)
	dG_table <- read.table(full_path, sep=",", row.names=1, header=TRUE)
	colnames(dG_table) <- gsub("^(X)", colnames(dG_table), replace="")
	dG_table
}



GetDeltaGMatrix <- function(mirna, kmer_len, register, xpos=10, ypos=20, height=4,
                            width=4, pdf.plot=FALSE) {

	dG_tables_list <- lapply(c(kmer_len, kmer_len + 1), GetKdMatrix, mirna=mirna,
	                         register=register)
	ddG_table <- dG_tables_list[[2]] - dG_tables_list[[1]]
	SubfunctionCall(GenericFigureSaveFile)
	xmin <- 0
	xmax <- 16
	ymin <- -1
	ymax <- 5
	BlankPlot()
	AddLinearAxis(1, tick.space=1, label.space=4, label="Loop length")
	x_vals <- 1:ncol(ddG_table) - 1
	AddLinearAxis(2, tick.space=0.2, label.space=1, label = "log(2) 7mer/6mer Kd")
	lines(x_vals, colMeans(dG_tables_list[[1]]), col="blue")
	lines(x_vals, colMeans(dG_tables_list[[2]]), col="red")
	xy <- GetPlotFractionalCoords(0.9, 0.9)
	legend_labels <- sprintf("%s-nt kmer to position %s", c())
	legend(x=xy[1], y=xy[2], legend=c())
	lines(x_vals, colMeans(ddG_table), lty=2)
	abline(0, 0, lty=2, col="gray")
	xy <- GetPlotFractionalCoords(0.05, 0.95)
	text(xy[1], xy[2], mirna, adj=0)
	xy <- GetPlotFractionalCoords(0.05, 0.85)
	text(xy[1], xy[2], sprintf("Beginning at position %s", register), adj=0)

	if (class(pdf.plot) == "character") {
    dev.off()
	}
}

GetPositionalDeltaG <- function(mirna, kmer_len, position, orientation,
                                ave_window=5, xpos=10, ypos=20, height=4,
                                width=4, pdf.plot=FALSE) {
	# Convert position, orientation, and kmer_len to register_1 and register_2.
	if (orientation == "right") {
		register_1 <- position - kmer_len
		register_2 <- register_1
	} else if (orientation == "left") {
		register_1 <- position + 1
		register_2 <- position

	}
	# Make the two dG tables, and subtract them.
	dG_table_1 <- GetKdMatrix(mirna, register_1, kmer_len)
	dG_table_2 <- GetKdMatrix(mirna, register_2, kmer_len + 1)
	# Adjust the columns corresponding to loop length if the left-hand-position is
	# the position being queried.
	if (orientation == "left") {
		dG_table_1 <- dG_table_1[, -1]
		dG_table_2 <- dG_table_2[, -ncol(dG_table_2)]
	}
	# Subtract the two matrices.
	ddG_table <- dG_table_2 - dG_table_1
	# This is the average delta G of each loop length, across the two.
	dG_ll_means <- colMeans(dG_table_1 + dG_table_2)/2
	ddG_ll <- colMeans(ddG_table)
	# Define the starting positions for the averaging window:
	pos_l_vec <- seq(1, length(dG_ll_means) - ave_window + 1)
	# Define the values of the averaging window
	dG_ll_means_win <- sapply(pos_l_vec, function(pos) {
		mean(dG_ll_means[pos:(pos + ave_window - 1)])
	})
	# Calculate the actual matrix of values corresponding to the centered windows.
	dG_window <- (
    which(dG_ll_means_win == max(dG_ll_means_win)) + ceiling(ave_window/2) - 1
	)
	# SubfunctionCall(GenericFigureSaveFile)
	# xmin <- 0
	# xmax <- 16
	# ymin <- -1
	# ymax <- 5
	# BlankPlot()
	# AddLinearAxis(1, tick.space=1, label.space=4, label="Loop length")
	# x_vals <- 1:ncol(ddG_table) - 1
	# x_vals_ave_win <- seq(ceiling(ave_window/2),
	#                       ceiling(ave_window/2) + length(dG_ll_means_win) - 1) - 1
	# AddLinearAxis(2, tick.space=0.2, label.space=1, label = "log(2) 7mer/6mer Kd")
	# lines(x_vals, colMeans(dG_table_1), col="blue")
	# lines(x_vals, colMeans(dG_table_2), col="red")
	# xy <- GetPlotFractionalCoords(0.9, 0.9)
	# # legend_labels <- sprintf("%s-nt kmer to position %s", c())
	# # legend(x=xy[1], y=xy[2], legend=c())
	# lines(x_vals, ddG_ll, lty=2)
	# lines(x_vals, dG_ll_means, lty=3)
	# points(x_vals_ave_win, dG_ll_means_win, lty=3, pch=20)
	# abline(0, 0, lty=2, col="gray")
	# xy <- GetPlotFractionalCoords(0.05, 0.95)
	# text(xy[1], xy[2], mirna, adj=0)
	# xy <- GetPlotFractionalCoords(0.05, 0.85)
	# text(xy[1], xy[2],
	#      sprintf("Querying position %s\n from %s", position, orientation), adj=0)
	# if (class(pdf.plot) == "character") {
 #    dev.off()
	# }
	return(ddG_ll[dG_window])

}

GetAllPositionsOfLength <- function(mirna, kmer_len, orientation) {
	mirna_str <- kMirnaSeqs[mirna]
	print(mirna_str)
	print(nchar(mirna_str))
	print(kmer_len)
	max_registers <- nchar(mirna_str) - 8 - kmer_len + 1
	if (orientation == "right") {
		positions_all <- (9 + kmer_len):nchar(mirna_str)	
	} else if (orientation == "left") {
		positions_all <- 9:(nchar(mirna_str) - kmer_len)
	}
	print(positions_all)
	ddGs <- sapply(positions_all, function(position) {
		print(mirna)
		print(kmer_len)
		print(position)
		SubfunctionCall(GetPositionalDeltaG) 	
	})
	names(ddGs) <- positions_all
	print(ddGs)
}

GetAllPositionsAllLength <- function(mirna, orientation, xpos=20, ypos=20,
                                     pdf.plot=FALSE) {
	mirna_str <- kMirnaSeqs[mirna]
	print(mirna_str)
	positions <- 9:nchar(mirna_str)

	Matrix_all <- matrix(NaN, nrow=4, ncol=length(positions),
	                     dimnames=list(c(5:8), positions))
	print(Matrix_all)
	for (kmer_len in rownames(Matrix_all)) {
		kmer_len_int <- as.numeric(kmer_len)
		output <- SubfunctionCall(GetAllPositionsOfLength, kmer_len=kmer_len_int)
		print(output)
		Matrix_all[kmer_len, names(output)] <- output
	}
	Matrix_all

	SubfunctionCall(GenericFigureSaveFile)
	xmin <- 8
	xmax <- 24
	ymin <- -2
	ymax <- 3
	print(xmin)
	print(xmax)
	BlankPlot()
	AddLinearAxis(1, tick.space=1, label.space=4, label="Position")
	# x_vals <- 1:ncol(ddG_table) - 1
	# x_vals_ave_win <- seq(ceiling(ave_window/2),
	#                       ceiling(ave_window/2) + length(dG_ll_means_win) - 1) - 1
	AddLinearAxis(2, tick.space=0.2, label.space=1, label = "log(2) Fold-improvement with position")

	kRowColors <- c(`5`="red", `6`="orange", `7`="forestgreen", `8`="blue")
	segments(xmin, 0, xmax, 0, lty=line_dash_length)
	sapply(rownames(Matrix_all), function(row_name) {
		print(colnames(Matrix_all))
		print(row_name)
		lines(colnames(Matrix_all), Matrix_all[row_name, ],
		       col=kRowColors[row_name], type="o", lwd=2, pch=20)
	})
	# lines(x_vals, colMeans(dG_table_1), col="blue")
	# lines(x_vals, colMeans(dG_table_2), col="red")
	if (orientation == "right") {
		xy <- GetPlotFractionalCoords(0.05, 0.955)
		text.adj <- 0
	} else if (orientation == "left") {
		xy <- GetPlotFractionalCoords(0.75, 0.95)
		text.adj <- 0.5
	}
	text.xy <- GetPlotFractionalCoords(0, 0.955)[2]
	text(xy[1], text.xy, labels="Background k-mer context:", adj=text.adj)
	legend_labels <- sprintf("%s nt", 6:9)
	legend(x=xy[1], y=xy[2], legend=legend_labels, lwd=2, pch=20, col=kRowColors,
	       bty="n")

	# lines(x_vals, ddG_ll, lty=2)
	# lines(x_vals, dG_ll_means, lty=3)
	# points(x_vals_ave_win, dG_ll_means_win, lty=3, pch=20)
	# abline(0, 0, lty=2, col="gray")
	xy <- GetPlotFractionalCoords(0.05, 0.05)
	text(xy[1], xy[2], mirna, adj=0)
	# xy <- GetPlotFractionalCoords(0.05, 0.85)
	# text(xy[1], xy[2],
	#      sprintf("Querying position %s\n from %s", position, orientation), adj=0)
	if (class(pdf.plot) == "character") {
    dev.off()
	}


}


# GetPositionalDeltaG("miR-155", 5, 16, "right")
# GetPositionalDeltaG("miR-155", 5, 16, "left")

# GetAllPositionsOfLength("miR-155", 5, "right")
# GetAllPositionsOfLength("miR-155", 5, "left")

GetAllPositionsAllLength("miR-155", "left", pdf.plot="190717/miR-155_left")
GetAllPositionsAllLength("miR-155", "right", pdf.plot="190717/miR-155_right")
GetAllPositionsAllLength("miR-1", "left", pdf.plot="190717/miR-1_left")
GetAllPositionsAllLength("miR-1", "right", pdf.plot="190717/miR-1_right")
GetAllPositionsAllLength("let-7a", "left", pdf.plot="190717/let-7a_left")
GetAllPositionsAllLength("let-7a", "right", pdf.plot="190717/let-7a_right")



break
# GetDeltaGMatrix("let-7a",  9, pdf.plot="190625/190625 let-7a_3p_register9")
# GetDeltaGMatrix("let-7a", 10, pdf.plot="190625/190625 let-7a_3p_register10")
# GetDeltaGMatrix("let-7a", 11, pdf.plot="190625/190625 let-7a_3p_register11")
# GetDeltaGMatrix("let-7a", 12, pdf.plot="190625/190625 let-7a_3p_register12")
# GetDeltaGMatrix("let-7a", 13, pdf.plot="190625/190625 let-7a_3p_register13")
# GetDeltaGMatrix("let-7a", 14, pdf.plot="190625/190625 let-7a_3p_register14")
# GetDeltaGMatrix("let-7a", 15, pdf.plot="190625/190625 let-7a_3p_register15")
# GetDeltaGMatrix("let-7a", 16, pdf.plot="190625/190625 let-7a_3p_register16")


# GetDeltaGMatrix("miR-1",  9, pdf.plot="190625/190625 miR-1_3p_register9")
# GetDeltaGMatrix("miR-1", 10, pdf.plot="190625/190625 miR-1_3p_register10")
# GetDeltaGMatrix("miR-1", 11, pdf.plot="190625/190625 miR-1_3p_register11")
# GetDeltaGMatrix("miR-1", 12, pdf.plot="190625/190625 miR-1_3p_register12")
# GetDeltaGMatrix("miR-1", 13, pdf.plot="190625/190625 miR-1_3p_register13")
# GetDeltaGMatrix("miR-1", 14, pdf.plot="190625/190625 miR-1_3p_register14")
# GetDeltaGMatrix("miR-1", 15, pdf.plot="190625/190625 miR-1_3p_register15")
# GetDeltaGMatrix("miR-1", 16, pdf.plot="190625/190625 miR-1_3p_register16")

# GetDeltaGMatrix("miR-155",  9, pdf.plot="190625/190625 miR-155_3p_register9")
# GetDeltaGMatrix("miR-155", 10, pdf.plot="190625/190625 miR-155_3p_register10")
# GetDeltaGMatrix("miR-155", 11, pdf.plot="190625/190625 miR-155_3p_register11")
# GetDeltaGMatrix("miR-155", 12, pdf.plot="190625/190625 miR-155_3p_register12")
# GetDeltaGMatrix("miR-155", 13, pdf.plot="190625/190625 miR-155_3p_register13")
# GetDeltaGMatrix("miR-155", 14, pdf.plot="190625/190625 miR-155_3p_register14")
# GetDeltaGMatrix("miR-155", 15, pdf.plot="190625/190625 miR-155_3p_register15")
# GetDeltaGMatrix("miR-155", 16, pdf.plot="190625/190625 miR-155_3p_register16")
# GetDeltaGMatrix("miR-155", 17, pdf.plot="190625/190625 miR-155_3p_register17")

# break

GetFracPlot <- function(base_size=7, xpos=10, ypos=20, height=4,
                            width=4, pdf.plot=FALSE) {


	SubfunctionCall(GenericFigureSaveFile)
	xmin <- 9
	xmax <- 23 - base_size
	ymin <- 0
	if (base_size == 8) {
		ymax <- 2
	} else {
		ymax <- 2	
	}
	# ymax <- 10
	BlankPlot()
	AddLinearAxis(1, tick.space=1, label.space=1, label="Register")
	AddLinearAxis(2, tick.space=0.2, label.space=1,
	              label=sprintf("log(2) %smer/%smer Kd", base_size + 1,
	                            base_size))


	sapply(c("miR-1", "let-7a", "miR-155"), function(mirna) {

		if (mirna == "miR-155") {
			registers <- 9:xmax
		} else {
			registers <- 9:(xmax - 1)
		}
		del_vector <- sapply(registers, function(register) {
			full_paths <- sapply(base_size + c(0, 1), function(k) {
				sprintf("/%s/%s/matdfchange%smer%s_reg%s.txt", base_path,
				        mirna_path[mirna], k, mirna_path[mirna], register)
			})	
			dG_tables_list <- lapply(full_paths, function(path) {
				dG_table <- read.table(path, sep=",", row.names=1, header=TRUE)
				colnames(dG_table) <- gsub("^(X)", colnames(dG_table), replace="")
				if (base_size == 8) {
					dG_table <- dG_table[, 1:(ncol(dG_table) - 1)]
				}
				dG_table
			})

			print(register)
			ddG_table <- colMeans(dG_tables_list[[2]] - dG_tables_list[[1]])

			# ddG_table <- dG_tables_list[[2]] - dG_tables_list[[1]]
			ddG_average <- Norm(colMeans(dG_tables_list[[2]] + dG_tables_list[[1]]))
			print(ddG_average)
			print(sum(ddG_average))

			sum(ddG_table*ddG_average)
			# print(ddG_table)
			# dG_7mer <- colMeans(dG_tables_list[[2]])
			# dG_7mer <- 1

			# mean((ddG_table/dG_7mer)[3:8])
		})

		print(del_vector)

		lines(registers, del_vector, col=kMirnaColors[mirna], pch=20, type="o", xpd=NA)



	})
	legend("topleft", legend=c("miR-1", "let-7a", "miR-155"), pch=20,
	       col=kMirnaColors[c("miR-1", "let-7a", "miR-155")], bty="n")

	if (class(pdf.plot) == "character") {
    dev.off()
	}


}


# m155_r9_k7 <- GetKdMatrix("miR-155", 9, 7, pdf.plot="190625 ")
# m155_r9_k8 <- GetKdMatrix("miR-155", 9, 8, pdf.plot="190625 ")
# m155_r9_k9 <- GetKdMatrix("miR-155", 9, 9, pdf.plot="190625 ")

# break



GetFracPlot(base_size=5, pdf.plot="190625/base_5nt")
GetFracPlot(base_size=6, pdf.plot="190625/base_6nt")
GetFracPlot(pdf.plot="190625/base_7nt")
GetFracPlot(base_size=8, pdf.plot="190625/base_8nt")



