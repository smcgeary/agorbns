################################################################################
#Analyzestructure.py
################################################################################
# library(gplots)
source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")
# library(wrswoR)
get_total_structures <- function(mirna, site, condition,start,width) {
	file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
					    "/equilibrium/probability_unpaired_by_site/reduced/", site, "/", condition, "_-3--3_reduced_", start, "-",width,"_constant-False_mir_center-True.txt")
		out <- read.table(file_name,header=FALSE,skip=1,sep="\t",stringsAsFactors=FALSE)
	return(out)
}

get_geo_average_structures <- function(mirna, site, condition,start,width) {
	file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
					    "/equilibrium/probability_unpaired_geometric_mean_by_site/reduced/", site, "/", condition, "_-3--3_reduced_", start, "-",width,"_constant-False_mir_center-True.txt")
		out <- read.table(file_name,header=FALSE,skip=1,sep="\t",stringsAsFactors=FALSE)
	return(out)
}

get_site_flanks <- function(data) {
	apply(data, 1, function(row) {
		return(paste0(row[2:5],collapse=""))
		})
}
TotalFlanks <- function(flanks) {
		sapply(unique(flanks_A), function(flank) {
		length(which(flanks[,1]==flank))/length(flanks[,1])
		})}
 
# A_flanks <- get_geo_average_structures("let-7a", "8mer", "0.4", 1, 1)

total_flanks_A <- TotalFlanks(A_flanks)[sort(unique(flanks_A))]

# sample <- get_total_structures(mirna, site, "I", 9, 5) 
# x <-TotalFlanks(sample[sample_int_R(20*nrow(sample),
#  									   round(0.6*20*nrow(sample)),
#  									   prob = rep(sample[,2], 20)), ])[sort(unique(flanks_A))]
# y <- total_flanks_A[sort(unique(flanks_A))]
#  plot(x,y, log = 'xy')
#  inds <- which(x > 0 & y > 0)
#  title(main = cor(log(x[inds]),log(y[inds])))


GetCorrelationGeo <- function(condition, position, window_length, samplingdepth, foldexcess, k = 1, input_norm = FALSE) {
	sample <- get_geo_average_structures(mirna, site, condition, position, window_length)
	sample <- sample[!(is.na(sample[,2])),]
	I_initial <- TotalFlanks(sample)
	I_flanks_selection <- TotalFlanks(sample[sample_int_R(foldexcess*nrow(sample),
 									   round(samplingdepth*foldexcess*nrow(sample)),
 									   prob = rep((sample[,2])^k, foldexcess)), ])[sort(unique(flanks_A))]
	ind_use <- which(I_initial > 0 & I_flanks_selection > 0 & total_flanks_A > 0)
	# plot(I_flanks_selection[ind_use]/I_initial[ind_use], total_flanks_A[ind_use]/I_initial[ind_use], log = 'xy')
	if (input_norm == TRUE){
		out <- cor(log10(I_flanks_selection[ind_use]) - log10(I_initial[ind_use]),
				log10(total_flanks_A[ind_use]) - log10(I_initial[ind_use]))
	} else {
		out <- cor(log10(I_flanks_selection[ind_use]),
				log10(total_flanks_A[ind_use]))
	}
	print(out)
	return(out)
}



GetCorrelationTotal <- function(condition, position, window_length, samplingdepth, foldexcess, input_norm = FALSE) {
	sample <- get_total_structures(mirna, site, condition, position, window_length)
	sample <- sample[!(is.na(sample[,2])),]

	I_initial <- TotalFlanks(sample)
	I_flanks_selection <- TotalFlanks(sample[sample_int_R(foldexcess*nrow(sample),
 									   round(samplingdepth*foldexcess*nrow(sample)),
 									   prob = rep(sample[,2], foldexcess)), ])[sort(unique(flanks_A))]
	ind_use <- which(I_initial > 0 & I_flanks_selection > 0 & total_flanks_A > 0)
	# plot(I_flanks_selection[ind_use]/I_initial[ind_use], total_flanks_A[ind_use]/I_initial[ind_use], log = 'xy')

	if (input_norm == TRUE){
		out <- cor(log10(I_flanks_selection[ind_use]) - log10(I_initial[ind_use]),
				log10(total_flanks_A[ind_use]) - log10(I_initial[ind_use]))
	} else {
		out <- cor(log10(I_flanks_selection[ind_use]),
				log10(total_flanks_A[ind_use]))
	}

	print(out)
	return(out)
}
# matrix_final_addition <- matrix(0, nrow = 38, ncol = 5)
# sapply(seq(16,20), function(j) {
# 	print(j)
# 	sapply(seq(-17, 20), function(i) {
# 		print(i)
# 		attemp <- try(GetCorrelationTotal("I", i, j, 0.6, 5), silent = TRUE)
# 		if (class(attemp) == "numeric") {
# 		matrix_final_addition[i + 18, j-15] <<- try(GetCorrelationTotal("I", i, j, 0.6, 5), silent = TRUE)

# 		}
# 		})
# 	})

# matrix_final_addition_2 <- matrix(0, nrow = 38, ncol = 5)
# sapply(seq(16,20), function(j) {
# 	print(j)
# 	sapply(seq(-17, 20), function(i) {
# 		print(i)
# 		attemp <- try(GetCorrelationTotal("I", i, j, 0.6, 5, input_norm = TRUE), silent = TRUE)
# 		if (class(attemp) == "numeric") {
# 		matrix_final_addition_2[i + 18, j-15] <<- try(GetCorrelationTotal("I", i, j, 0.6, 5, input_norm = TRUE), silent = TRUE)
			
# 		}
# 		})
# 	})

# matrix_final_addition_3 <- matrix(0, nrow = 38, ncol = 5)
# sapply(seq(16,20), function(j) {
# 	print(j)
# 	sapply(seq(-17, 20), function(i) {
# 		print(i)
# 		attemp <- try(GetCorrelationGeo("I", i, j, 0.6, 5), silent = TRUE)
# 		if (class(attemp) == "numeric") {
# 		matrix_final_addition_3[i + 18, j-15] <<- try(GetCorrelationGeo("I", i, j, 0.6, 5), silent = TRUE)
			
# 		}
# 		})
# 	})

# matrix_final_addition_4 <- matrix(0, nrow = 38, ncol = 5)
# sapply(seq(16,20), function(j) {
# 	print(j)
# 	sapply(seq(-17, 20), function(i) {
# 		print(i)
# 		attemp <- try(GetCorrelationGeo("I", i, j, 0.6, 5, input_norm = TRUE), silent = TRUE)
# 		if (class(attemp) == "numeric") {
# 		matrix_final_addition_4[i + 18, j-15] <<- try(GetCorrelationGeo("I", i, j, 0.6, 5, input_norm = TRUE), silent = TRUE)
			
# 		}
# 		})
# 	})



# matrix_final_addition_2 <- sapply(seq(16,20), function(j) {
# 	print(j)
# 	sapply(seq(-17, 20), function(i) {
# 		return(GetCorrelationTotal("I", i, j, 0.6, 5, input_norm = TRUE))
# 		})
# 	})

# matrix_final_addition_3 <- sapply(seq(16,20), function(j) {
# 	print(j)
# 	sapply(seq(-17, 20), function(i) {
# 		return(GetCorrelationGeo("I", i, j, 0.6, 5))
# 		})
# 	})

# matrix_final_5 <- matrix(0, nrow = 38, ncol = 20)
# sapply(seq(1,20), function(j) {
# 	print(j)
# 	sapply(seq(-17, 20), function(i) {
# 		print(i)
# 		attemp <- try(GetCorrelationGeo("I", i, j, 0.6, 5, k = 3), silent = TRUE)
# 		if (class(attemp) == "numeric") {
# 		matrix_final_5[i + 18, j] <<- try(GetCorrelationGeo("I", i, j, 0.6, 5, k = 3), silent = TRUE)
			
# 		}
# 		})
# 	})

# matrix_final_6 <- matrix(0, nrow = 38, ncol = 20)
# sapply(seq(1,20), function(j) {
# 	print(j)
# 	sapply(seq(-17, 20), function(i) {
# 		print(i)
# 		attemp <- try(GetCorrelationGeo("I", i, j, 0.6, 5, k = 3, input_norm = TRUE), silent = TRUE)
# 		if (class(attemp) == "numeric") {
# 		matrix_final_6[i + 18, j] <<- try(GetCorrelationGeo("I", i, j, 0.6, 5, k = 3, input_norm = TRUE), silent = TRUE)
			
# 		}
# 		})
# 	})





# matrix_final <- sapply(seq(15), function(j) {
# 	sapply(seq(-17, 20), function(i) {
# 		return(GetCorrelationTotal("I", i, j, 0.6, 5))
# 		})
# 	})

# matrix_final_2 <- sapply(seq(15), function(j) {
# 	sapply(seq(-17, 20), function(i) {
# 		return(GetCorrelationTotal("I", i, j, 0.6, 5, input_norm = TRUE))
# 		})
# 	})

# matrix_final_3 <- sapply(seq(15), function(j) {
# 	sapply(seq(-17, 20), function(i) {
# 		return(GetCorrelationGeo("I", i, j, 0.6, 5))
# 		})
# 	})

# matrix_final_4 <- sapply(seq(15), function(j) {
# 	sapply(seq(-17, 20), function(i) {
# 		return(GetCorrelationGeo("I", i, j, 0.6, 5, input_norm = TRUE))
# 		})
# # 	})
# matrix_round <- cbind(matrix_final,matrix_final_addition)^2 - min(cbind(matrix_final,matrix_final_addition)^2)
# matrix_round_round <- matrix_round / max(matrix_round)

# dev.new(xpos = 20, ypos = 20, width = 10, height = 8)
# plot(seq(-20,20),c(seq(-20,19)*0,20),col="white")

# seeds <- c("NNDTACCTCBNN", "NNBCATTCCBNN", "NNBGCATTABNN", "NNHTGCCTTBNN",
# 				   "NNBTACAAABNN")

# names(seeds) <- c("let-7a", "miR-1", "miR-155", "miR-124", "lsy-6")

# seed_string <- strsplit(seeds[mirna], split = "")[[1]]

# sapply(-2:10, function(i) {
# 	       text(x = i,
# 	       	    y = 0.1,
# 	       	    labels = seed_string[i + 3])
# 	     }
# 	     )

# segments(-0.5, 0, -0.5, 20, lwd = 2, lty = 2)
# segments(7.5, 0, 7.5, 20, lwd = 2, lty = 2)
# segments(-0.5, 0, -0.5, 20, lwd = 2, lty = 2)
# segments(-2.5, 0, -2.5, 20, lwd = 1, lty = 2)
# segments(9.5, 0, 9.5, 20, lwd = 1, lty = 2)

# for (j in seq(20)) {
# 	x_tot <- c()
# 	y_tot <- c()
# 	for (i in seq(-17, 20)) {
# 		print(j)
# 		print(i)
# 		left <- 7 - i - j + 1
# 		right <- left+j + 1
# 		x_now <- mean(c(left, right))
# 		print(x_now)
# 		col_ <- "black"
# 		col_ <- try(rainbow(100,start = 0.05,end=0.7)[100-round(matrix_round_round[i + 18,j]*99)])
# 		 rect(x_now-0.5,j-0.5,x_now+0.5,j+0.5,col=col_,
# 					 border=NA)

# }
# }
# matrix_round_2 <- cbind(matrix_final_2,matrix_final_addition_2)^2 - min(cbind(matrix_final_2,matrix_final_addition_2)^2)
# matrix_round_round_2 <- matrix_round_2 / max(matrix_round_2)
# dev.new(xpos = 20, ypos = 20, width = 10, height = 8)
# plot(seq(-20,20),c(seq(-20,19)*0,20),col="white")

# seeds <- c("NNDTACCTCBNN", "NNBCATTCCBNN", "NNBGCATTABNN", "NNHTGCCTTBNN",
# 				   "NNBTACAAABNN")

# names(seeds) <- c("let-7a", "miR-1", "miR-155", "miR-124", "lsy-6")

# seed_string <- strsplit(seeds[mirna], split = "")[[1]]

# sapply(-2:10, function(i) {
# 	       text(x = i,
# 	       	    y = 0.1,
# 	       	    labels = seed_string[i + 3])
# 	     }
# 	     )

# segments(-0.5, 0, -0.5, 20, lwd = 2, lty = 2)
# segments(7.5, 0, 7.5, 20, lwd = 2, lty = 2)
# segments(-0.5, 0, -0.5, 20, lwd = 2, lty = 2)
# segments(-2.5, 0, -2.5, 20, lwd = 1, lty = 2)
# segments(9.5, 0, 9.5, 20, lwd = 1, lty = 2)

# for (j in seq(20)) {
# 	x_tot <- c()
# 	y_tot <- c()
# 	for (i in seq(-17, 20)) {
# 		print(j)
# 		print(i)
# 		left <- 7 - i - j + 1
# 		right <- left+j + 1
# 		x_now <- mean(c(left, right))
# 		print(x_now)
# 		col_ <- "black"
# 		col_ <- try(rainbow(100,start = 0.05,end=0.7)[100-round(matrix_round_round_2[i + 18, j]*99)])
# 		 rect(x_now-0.5,j-0.5,x_now+0.5,j+0.5,col=col_,
# 					 border=NA)

# }
# }


# matrix_round_3 <- cbind(matrix_final_3,matrix_final_addition_3)^2 - min(cbind(matrix_final_3,matrix_final_addition_3)^2)
# matrix_round_round_3 <- matrix_round_3 / max(matrix_round_3)


# dev.new(xpos = 20, ypos = 20, width = 10, height = 8)
# plot(seq(-20,20),c(seq(-20,19)*0,20),col="white")

# seeds <- c("NNDTACCTCBNN", "NNBCATTCCBNN", "NNBGCATTABNN", "NNHTGCCTTBNN",
# 				   "NNBTACAAABNN")

# names(seeds) <- c("let-7a", "miR-1", "miR-155", "miR-124", "lsy-6")

# seed_string <- strsplit(seeds[mirna], split = "")[[1]]

# sapply(-2:10, function(i) {
# 	       text(x = i,
# 	       	    y = 0.1,
# 	       	    labels = seed_string[i + 3])
# 	     }
# 	     )

# segments(-0.5, 0, -0.5, 20, lwd = 2, lty = 2)
# segments(7.5, 0, 7.5, 20, lwd = 2, lty = 2)
# segments(-0.5, 0, -0.5, 20, lwd = 2, lty = 2)
# segments(-2.5, 0, -2.5, 20, lwd = 1, lty = 2)
# segments(9.5, 0, 9.5, 20, lwd = 1, lty = 2)

# for (j in seq(20)) {
# 	x_tot <- c()
# 	y_tot <- c()
# 	for (i in seq(-17, 20)) {
# 		print(j)
# 		print(i)
# 		left <- 7 - i - j + 1
# 		right <- left+j + 1
# 		x_now <- mean(c(left, right))
# 		print(x_now)
# 		col_ <- "black"
# 		col_ <- try(rainbow(100,start = 0.05,end=0.7)[100-round(matrix_round_round_3[i + 18, j]*99)])
# 		 rect(x_now-0.5,j-0.5,x_now+0.5,j+0.5,col=col_,
# 					 border=NA)

# }
# }




# matrix_round_4 <- cbind(matrix_final_4,matrix_final_addition_4)^2 - min(cbind(matrix_final_4,matrix_final_addition_4)^2)
# matrix_round_round_4 <- matrix_round_4 / max(matrix_round_4)


# dev.new(xpos = 20, ypos = 20, width = 10, height = 8)
# plot(seq(-20,20),c(seq(-20,19)*0,20),col="white")

# seeds <- c("NNDTACCTCBNN", "NNBCATTCCBNN", "NNBGCATTABNN", "NNHTGCCTTBNN",
# 				   "NNBTACAAABNN")

# names(seeds) <- c("let-7a", "miR-1", "miR-155", "miR-124", "lsy-6")

# seed_string <- strsplit(seeds[mirna], split = "")[[1]]

# sapply(-2:10, function(i) {
# 	       text(x = i,
# 	       	    y = 0.1,
# 	       	    labels = seed_string[i + 3])
# 	     }
# 	     )

# segments(-0.5, 0, -0.5, 20, lwd = 2, lty = 2)
# segments(7.5, 0, 7.5, 20, lwd = 2, lty = 2)
# segments(-0.5, 0, -0.5, 20, lwd = 2, lty = 2)
# segments(-2.5, 0, -2.5, 20, lwd = 1, lty = 2)
# segments(9.5, 0, 9.5, 20, lwd = 1, lty = 2)

# for (j in seq(20)) {
# 	x_tot <- c()
# 	y_tot <- c()
# 	for (i in seq(-17, 20)) {
# 		print(j)
# 		print(i)
# 		left <- 7 - i - j + 1
# 		right <- left+j + 1
# 		x_now <- mean(c(left, right))
# 		print(x_now)
# 		col_ <- "black"
# 		col_ <- try(rainbow(100,start = 0.05,end=0.7)[100-round(matrix_round_round_4[i + 18, j]*99)])
# 		 rect(x_now-0.5,j-0.5,x_now+0.5,j+0.5,col=col_,
# 					 border=NA)

# }
# # }
matrix_round_5 <- matrix_final_5^2 - min(matrix_final_5^2)
matrix_round_round_5 <- matrix_round_5 / max(matrix_round_5)

dev.new(xpos = 20, ypos = 20, width = 10, height = 8)
plot(seq(-20,20),c(seq(-20,19)*0,20),col="white")

seeds <- c("NNDTACCTCBNN", "NNBCATTCCBNN", "NNBGCATTABNN", "NNHTGCCTTBNN",
				   "NNBTACAAABNN")

names(seeds) <- c("let-7a", "miR-1", "miR-155", "miR-124", "lsy-6")

seed_string <- strsplit(seeds[mirna], split = "")[[1]]

sapply(-2:10, function(i) {
	       text(x = i,
	       	    y = 0.1,
	       	    labels = seed_string[i + 3])
	     }
	     )

segments(-0.5, 0, -0.5, 20, lwd = 2, lty = 2)
segments(7.5, 0, 7.5, 20, lwd = 2, lty = 2)
segments(-0.5, 0, -0.5, 20, lwd = 2, lty = 2)
segments(-2.5, 0, -2.5, 20, lwd = 1, lty = 2)
segments(9.5, 0, 9.5, 20, lwd = 1, lty = 2)

for (j in seq(20)) {
	x_tot <- c()
	y_tot <- c()
	for (i in seq(-17, 20)) {
		print(j)
		print(i)
		left <- 7 - i - j + 1
		right <- left+j + 1
		x_now <- mean(c(left, right))
		print(x_now)
		col_ <- "black"
		col_ <- try(rainbow(100,start = 0.05,end=0.7)[100-round(matrix_round_round_5[i + 18, j]*99)])
		 rect(x_now-0.5,j-0.5,x_now+0.5,j+0.5,col=col_,
					 border=NA)

}
}



dev.new(xpos = 20, ypos = 20, width = 10, height = 8)
plot(seq(-20,20),c(seq(-20,19)*0,20),col="white")

seeds <- c("NNDTACCTCBNN", "NNBCATTCCBNN", "NNBGCATTABNN", "NNHTGCCTTBNN",
				   "NNBTACAAABNN")

names(seeds) <- c("let-7a", "miR-1", "miR-155", "miR-124", "lsy-6")

seed_string <- strsplit(seeds[mirna], split = "")[[1]]
matrix_round_6 <- matrix_final_6^2 - min(matrix_final_6^2)
matrix_round_round_6 <- matrix_round_6 / max(matrix_round_6)

sapply(-2:10, function(i) {
	       text(x = i,
	       	    y = 0.1,
	       	    labels = seed_string[i + 3])
	     }
	     )

segments(-0.5, 0, -0.5, 20, lwd = 2, lty = 2)
segments(7.5, 0, 7.5, 20, lwd = 2, lty = 2)
segments(-0.5, 0, -0.5, 20, lwd = 2, lty = 2)
segments(-2.5, 0, -2.5, 20, lwd = 1, lty = 2)
segments(9.5, 0, 9.5, 20, lwd = 1, lty = 2)

for (j in seq(20)) {
	x_tot <- c()
	y_tot <- c()
	for (i in seq(-17, 20)) {
		print(j)
		print(i)
		left <- 7 - i - j + 1
		right <- left+j + 1
		x_now <- mean(c(left, right))
		print(x_now)
		col_ <- "black"
		col_ <- try(rainbow(100,start = 0.05,end=0.7)[100-round(matrix_round_round_6[i + 18, j]*99)])
		 rect(x_now-0.5,j-0.5,x_now+0.5,j+0.5,col=col_,
					 border=NA)

}
}




break


# s_I_8mer <- get_structures(mirna, site, "I", 1, 1)
# s_0.4_8mer <- get_structures(mirna, site, "0.4")

# r_I_8mer <- get_reads(mirna,site,"I")
# r_0.4_8mer <- get_reads(mirna,site,"0.4")

# rlogitnormal <- function(...){
#   x <- exp(rnorm(...))
#   x / (1+x)
# }
# print("hi")


get_site_probs <- function(data,temp_start=0,temp_stop=7) {
	apply(data, 1, function(row) {
		start <- as.numeric(row[1])
		# print(as.numeric(row[(start+5+temp_start):(start+5+7+temp_stop)]))
		# print(1-as.numeric(row[(start+5+temp_start):(start+5+7+temp_stop)]))
		# print(10^sum(log10(1-as.numeric(row[(start+5+temp_start):(start+5+7+temp_stop)]))))
		probs <- as.numeric(row[6:length(row)])
		if (temp_start+start+1>=1&temp_stop+start+1<=length(probs)){
		return(10^sum(log10(1-probs[(start+1+temp_start):(start+1+temp_stop)])))
		} else {
			return(0)
		}

		})
}

# get_site_read_probs_position <- function(data,read_data,position) {
# 	data <- data[1:10,]
# 	read_data <- read_data[1:10,,drop=FALSE]

# 	pos <- as.numeric(data[,1]+position)
# 	return(sapply(seq(dim(data)[1]), function(row){
# 		c(data[row,pos[row]+5],substr(read_data[row,1],pos,pos))
# 		}))
# }

x <- seq(-9,15)

# probs_I_pos <- sapply(x,function(i){get_site_probs(s_I_8mer,i,i)})
# probs_A_pos <- sapply(x,function(i){get_site_probs(s_0.4_8mer,i,i)})
# probs_I_pos_full <- probs_I_pos
# probs_A_pos_full <- probs_A_pos
# dev.new(width=10,height=7)
# par(mfrow=c(4,5))
# par(mar=c(2,2,2,2))

# # for (i in seq(dim(probs_I_pos)[2])) {
# # 	I_dat <- probs_I_pos[,i]
# # 	A_dat <- probs_A_pos[,i]
# # 	plot(ecdf(1-I_dat),main=x[i])
# # 	segments(mean(1-I_dat), 0, mean(1 - I_dat), 1)
# # 	plot(ecdf(1-A_dat),col="red",add=TRUE)
# # 	segments(mean(1-A_dat), 0, mean(1 - A_dat), 1)
# # }

# # dev.new(width=10,height=7)
# # par(mar=c(2,2,2,2))
# # par(mfrow=c(4,5))
# # for (i in seq(dim(probs_I_pos)[2])) {
# # 	I_dat <- probs_I_pos[,i]
# # 	A_dat <- probs_A_pos[,i]

# # 	plot(ecdf(log10(1-I_dat)),main=x[i],xlim=c(-4,0))
# # 	segments(mean(log10(1-I_dat)), 0, mean(log10(1 - I_dat)), 1)
# # 	plot(ecdf(log10(1-A_dat)),col="red",add=TRUE)
# # 	segments(mean(log10(1-A_dat)), 0, mean(log10(1 - A_dat)), 1)
# # }
# # dev.new(width=10,height=7)
# # par(mfrow=c(4,5))
# # par(mar=c(2,2,2,2))

# # for (i in seq(dim(probs_I_pos)[2])) {
# # 		I_dat <- probs_I_pos[,i]
# # 	A_dat <- probs_A_pos[,i]

# # 	plot(ecdf(log10(I_dat)),main=x[i],xlim=c(-2,0))
# # 	segments(mean(log10(I_dat)), 0, mean(log10(I_dat)), 1)
# # 	plot(ecdf(log10(A_dat)),col="red",add=TRUE)
# # 	segments(mean(log10(A_dat)), 0, mean(log10(A_dat)), 1)
# # }
# dev.new(width=9,height=10)
pos <- 30:50
f1 <- "T"
f2 <- "T"
# ind_I <- which(sapply(get_site_flanks(s_I_8mer),function(flank){strsplit(flank,split="")[[1]][1]})==f1&sapply(get_site_flanks(s_I_8mer),function(flank){strsplit(flank,split="")[[1]][2]})==f2&as.integer(s_I_8mer[, 1]) %in% pos)
# ind_A <- which(sapply(get_site_flanks(s_0.4_8mer),function(flank){strsplit(flank,split="")[[1]][1]})==f1&sapply(get_site_flanks(s_0.4_8mer),function(flank){strsplit(flank,split="")[[1]][2]})==f2&as.integer(s_0.4_8mer[, 1]) %in% pos)
ind_I <- which(as.integer(s_I_8mer_use[, 1]) %in% pos)
ind_A <- which(as.integer(s_0.4_8mer_use[, 1]) %in% pos)


probs_I_pos_temp <- probs_I_pos_full[ind_I,]
probs_A_pos_temp <- probs_A_pos_full[ind_A,]

par(par_plots)
par(mfrow=c(4,3))
text_labels <- c("N","N","N","N","N","N","N","N","N","N","C","U","A","C","C","U","C","A","N","N","N","N","N","N","N")
y_i <- colMeans(1-probs_I_pos_temp)
y_a <- colMeans(1-probs_A_pos_temp)
plot(x,y_i,type="l",ylim=c(-0.1,1),axes=FALSE,ann=FALSE)
axis(1,at=x,labels=FALSE,cex=0.9,pos=0,lwd=2)
axis(2,at=seq(0,1,0.2),pos=-9,lwd=2)
sapply(x,function(i) { text(x=i,-0.07,labels=text_labels[i+10])})
title(main=paste0(f1,f2))
title(ylab="Paired probability",line=3)
lines(x,y_a,col="red")
legend("topleft",legend=c("Input","0.4% AGO"),col=c("black","red"),lwd=2,bty="n",cex=1.5)

plot(x,y_a/y_i,type="l",ylim=c(0.45,1),axes=FALSE,ann=FALSE)
axis(1,at=x,labels=FALSE,cex=0.9,pos=0.5,lwd=2)
axis(2,at=seq(0.5,1,0.1),pos=-9,lwd=2)
sapply(x,function(i) { text(x=i,0.47,labels=text_labels[i+10])})
title(ylab="AGO/Input")

plot(x,log(y_a/y_i,base=2),type="l",ylim=c(-1.1,0.2),axes=FALSE,ann=FALSE)
axis(1,at=x,labels=FALSE,cex=0.9,pos=-1,lwd=2)
axis(2,at=seq(-1,0.2,0.3),pos=-9,lwd=2)
sapply(x,function(i) { text(x=i,-1.08,labels=text_labels[i+10])})
title(ylab="log2(AGO/Input)")

plot(x,1-y_i,type="l",ylim=c(-0.1,1),axes=FALSE,ann=FALSE)
axis(1,at=x,labels=FALSE,cex=0.9,pos=0,lwd=2)
axis(2,at=seq(0,1,0.2),pos=-9,lwd=2)
sapply(x,function(i) { text(x=i,-0.07,labels=text_labels[i+10])})
title(ylab="Unpaired probability",line=2)
lines(x,1-y_a,col="red")
legend("topleft",legend=c("Input","0.4% AGO"),col=c("black","red"),lwd=2,bty="n",cex=1.5)

plot(x,(1-y_a)/(1-y_i),type="l",ylim=c(0.95,1.5),axes=FALSE,ann=FALSE)
axis(1,at=x,labels=FALSE,cex=0.9,pos=1,lwd=2)
axis(2,at=seq(1,1.5,0.1),pos=-9,lwd=2)
sapply(x,function(i) { text(x=i,0.97,labels=text_labels[i+10])})
title(ylab="AGO/Input")

plot(x,log((1-y_a)/(1-y_i),base=2),type="l",ylim=c(-0.05,0.5),axes=FALSE,ann=FALSE)
axis(1,at=x,labels=FALSE,cex=0.9,pos=0,lwd=2)
axis(2,at=seq(0,0.5,0.1),pos=-9,lwd=2)
sapply(x,function(i) { text(x=i,-0.02,labels=text_labels[i+10])})
title(ylab="log2(AGO/Input)")
y_i <- 10^colMeans(log10(1-probs_I_pos_temp))
y_a <- 10^colMeans(log10(1-probs_A_pos_temp))

plot(x,y_i,type="l",ylim=c(-0.1,1),axes=FALSE,ann=FALSE)
axis(1,at=x,labels=FALSE,cex=0.9,pos=0,lwd=2)
axis(2,at=seq(0,1,0.2),pos=-9,lwd=2)
sapply(x,function(i) { text(x=i,-0.07,labels=text_labels[i+10])})
title(ylab="Paired probability")
lines(x,y_a,col="red")
legend("topleft",legend=c("Input","0.4% AGO"),col=c("black","red"),lwd=2,bty="n",cex=1.5)

plot(x,y_a/y_i,type="l",ylim=c(0.45,1),axes=FALSE,ann=FALSE)
axis(1,at=x,labels=FALSE,cex=0.9,pos=0.5,lwd=2)
axis(2,at=seq(0.5,1,0.1),pos=-9,lwd=2)
sapply(x,function(i) { text(x=i,0.47,labels=text_labels[i+10])})
title(ylab="AGO/Input")

plot(x,log(y_a/y_i,base=2),type="l",ylim=c(-1.1,0.2),axes=FALSE,ann=FALSE)
axis(1,at=x,labels=FALSE,cex=0.9,pos=-1,lwd=2)
axis(2,at=seq(-1,0.2,0.3),pos=-9,lwd=2)
sapply(x,function(i) { text(x=i,-1.08,labels=text_labels[i+10])})
title(ylab="log2(AGO/Input)")

y_i <- 10^colMeans(log10(probs_I_pos_temp))
y_a <- 10^colMeans(log10(probs_A_pos_temp))

plot(x,y_i,type="l",ylim=c(-0.1,1),axes=FALSE,ann=FALSE)
axis(1,at=x,labels=FALSE,cex=0.9,pos=0,lwd=2)
axis(2,at=seq(0,1,0.2),pos=-9,lwd=2)
sapply(x,function(i) { text(x=i,-0.07,labels=text_labels[i+10])})
title(ylab="Unpaired probability")
lines(x,y_a,col="red")
legend("topleft",legend=c("Input","0.4% AGO"),col=c("black","red"),lwd=2,bty="n",cex=1.5)

plot(x,y_a/y_i,type="l",ylim=c(0.85,2.5),axes=FALSE,ann=FALSE)
axis(1,at=x,labels=FALSE,cex=0.9,pos=1,lwd=2)
axis(2,at=seq(1,2.5,0.3),pos=-9,lwd=2)
sapply(x,function(i) { text(x=i,0.97,labels=text_labels[i+10])})
title(ylab="AGO/Input")

plot(x,log(y_a/y_i,base=2),type="l",ylim=c(-0.1,1.5),axes=FALSE,ann=FALSE)
axis(1,at=x,labels=FALSE,cex=0.9,pos=0,lwd=2)
axis(2,at=seq(0,1.5,0.3),pos=-9,lwd=2)
sapply(x,function(i) { text(x=i,-0.08,labels=text_labels[i+10])})
title(ylab="log2(AGO/Input)")

get_read_positions <- function(data,read_data,position) {
	apply(data, 1, function(row) {
		start <- as.numeric(row[1])
		# print(as.numeric(row[(start+5+temp_start):(start+5+7+temp_stop)]))
		# print(1-as.numeric(row[(start+5+temp_start):(start+5+7+temp_stop)]))
		# print(10^sum(log10(1-as.numeric(row[(start+5+temp_start):(start+5+7+temp_stop)]))))

		as.numeric(row[start+5+position-1])
		})
}

