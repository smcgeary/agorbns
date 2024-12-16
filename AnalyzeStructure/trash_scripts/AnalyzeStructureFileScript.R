################################################################################
#Analyzestructure.py
################################################################################

source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")
library(wrswoR)

get_structures <- function(mirna, site, condition) {
	file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
                      "/equilibrium/structures_bp_prob/", site, "/", condition,
                      "_0-0.txt")
	out <- read.table(file_name, header = FALSE, skip = 1, sep = "\t",
		                stringsAsFactors = FALSE)
	return(out)
}

get_reads <- function(mirna, site, condition) {
	file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
					            "/equilibrium/reads_by_site/", site, "/", condition,
					            "_0-0.txt")
	out <- read.table(file_name, header = FALSE, skip = 1, sep = "\t",
		                stringsAsFactors = FALSE)
	return(out)
}

get_site_flanks <- function(data) {
	apply(data, 1, function(row) {
		return(paste0(row[2:5],collapse=""))
		})
}
# args = commandArgs(trailingOnly=TRUE)
# mirna = args[1]
# site = args[2]
# condition = args[3]
# percent = as.numeric(args[4])

s_I <- get_structures(mirna, site, "I")
s_A <- get_structures(mirna, site, condition)


flanks_I <- get_site_flanks(s_I)
flanks_A <- get_site_flanks(s_A)


get_site_probs <- function(data,temp_start=0,temp_stop=7) {
	apply(data, 1, function(row) {
		start <- as.numeric(row[1])
		probs <- as.numeric(row[6:length(row)])
		if ((temp_start + start + 1) >= 1 &
			  (temp_stop + start + 1) <= length(probs)){
		return(10^sum(log10(1-probs[(start+1+temp_start):(start+1+temp_stop)])))
		} else {
			return(0)
		}

		})
}



pos <- 26:55
ind_A <- which(as.integer(s_A[, 1]) %in% pos)
flanks_A <- get_site_flanks(s_A[ind_A,])
# setEPS(width = 10)

# postscript(file = paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
#                          "/figures/structures/structure_flank_correlations/",
#                          mirna, "/", site, "_", condition, "_",
#                          format(percent, nsmall = 2), ".eps"))
par(par_plots)

plot(seq(-20,20),c(seq(-20,19)*0,20),col="white")

seeds <- c("NNDTACCTCBNN", "NNBCATTCCBNN", "NNBGCATTABNN", "NNHTGCCTTBNN",
				   "NNBTACAAABNN")

names(seeds) <- c("let-7a", "miR-1", "miR-155", "miR-124", "lsy-6")

seed_string <- strsplit(seeds[mirna], split = "")[[1]]

sapply(-3:8, function(i) {
	       text(x = i,
	       	    y = 0.1,
	       	    labels = seed_string[i + 4])
	     }
	     )

segments(-0.5, 0, -0.5, 20, lwd = 2, lty = 2)
segments(7.5, 0, 7.5, 20, lwd = 2, lty = 2)
segments(-0.5, 0, -0.5, 20, lwd = 2, lty = 2)
segments(-2.5, 0, -2.5, 20, lwd = 1, lty = 2)
segments(9.5, 0, 9.5, 20, lwd = 1, lty = 2)

tot <- 10
ind_A_all <- rep(seq(1, dim(s_I)[1]), tot)

total_flanks_A <- sapply(unique(flanks_A), function(flank) {
	length(which(flanks_A == flank)) / length(ind_A)
	})

for (j in seq(20)) {
	x_tot <- c()
	y_tot <- c()
	for (i in seq(-17, 20)) {
		left <- i-floor((j-1)/2)
		right <- left+j-1
		x_now <- mean(c(left, right))
		x_tot <- c(x_tot,x_now)
		print(c(left,right,x_now))

		print("pre sample")
		ind_A_sample <- ind_A_all[sample_int_R(length(ind_A_all),round(length(ind_A_all)*percent),prob=get_site_probs(s_I[ind_A_all,],left,right))]

		s_A_sample <- s_I[ind_A_sample,]	
		print("post sample")
		ind_A_sample <- which(as.integer(s_A_sample[, 1]) %in% pos)
		flanks_A_sample <- get_site_flanks(s_A_sample[ind_A_sample,])
		total_flanks_A_sample <- sapply(unique(flanks_A), function(flank) {
			length(which(flanks_A_sample == flank))/length(ind_A_sample)
			})
		ind_keep <- which(total_flanks_A_sample > 0 & total_flanks_A > 0)
		y_new <- cor(log(total_flanks_A_sample[ind_keep]),
			           log(total_flanks_A[ind_keep])
			           )
		y_tot <- c(y_tot, y_new)
		print(y_new)
					 rect(x_now-0.5,j-0.5,x_now+0.5,j+0.5,col=rainbow(100,end=0.67)[100-round(y_new*99)],
					 border=NA)

	}
}

# dev.off()