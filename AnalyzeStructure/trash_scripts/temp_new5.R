# wI2l <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/dinucleotide_frequencies_by_motif/CT/I_left.txt")
# wI3l <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/dinucleotide_frequencies_by_motif/CTA/I_left.txt")
# wI4l <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/dinucleotide_frequencies_by_motif/CTAC/I_left.txt")
# wI5l <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/dinucleotide_frequencies_by_motif/CTACC/I_left.txt")
# wI6l <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/dinucleotide_frequencies_by_motif/CTACCT/I_left.txt")
# wI7l <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/dinucleotide_frequencies_by_motif/CTACCTC/I_left.txt")
# wI8l <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/dinucleotide_frequencies_by_motif/CTACCTCA/I_left.txt")

# wI2r <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/dinucleotide_frequencies_by_motif/CA/I_right.txt")
# wI3r <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/dinucleotide_frequencies_by_motif/TCA/I_right.txt")
# wI4r <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/dinucleotide_frequencies_by_motif/CTCA/I_right.txt")
# wI5r <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/dinucleotide_frequencies_by_motif/CCTCA/I_right.txt")
# wI6r <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/dinucleotide_frequencies_by_motif/ACCTCA/I_right.txt")
# wI7r <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/dinucleotide_frequencies_by_motif/TACCTCA/I_right.txt")
# wI8r <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/dinucleotide_frequencies_by_motif/CTACCTCA/I_right.txt")

# w02l <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/dinucleotide_frequencies_by_motif/CT/0_left.txt")
# w03l <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/dinucleotide_frequencies_by_motif/CTA/0_left.txt")
# w04l <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/dinucleotide_frequencies_by_motif/CTAC/0_left.txt")
# w05l <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/dinucleotide_frequencies_by_motif/CTACC/0_left.txt")
# w06l <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/dinucleotide_frequencies_by_motif/CTACCT/0_left.txt")
# w07l <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/dinucleotide_frequencies_by_motif/CTACCTC/0_left.txt")
# w08l <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/dinucleotide_frequencies_by_motif/CTACCTCA/0_left.txt")

# w02r <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/dinucleotide_frequencies_by_motif/CA/0_right.txt")
# w03r <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/dinucleotide_frequencies_by_motif/TCA/0_right.txt")
# w04r <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/dinucleotide_frequencies_by_motif/CTCA/0_right.txt")
# w05r <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/dinucleotide_frequencies_by_motif/CCTCA/0_right.txt")
# w06r <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/dinucleotide_frequencies_by_motif/ACCTCA/0_right.txt")
# w07r <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/dinucleotide_frequencies_by_motif/TACCTCA/0_right.txt")
get_site_probs_alt <- function(data,temp_start=0,temp_stop=7,constant=TRUE) {
	apply(data, 1, function(row) {
		start <- as.numeric(row[1])
		# print(as.numeric(row[(start+5+temp_start):(start+5+7+temp_stop)]))
		# print(1-as.numeric(row[(start+5+temp_start):(start+5+7+temp_stop)]))
		# print(10^sum(log10(1-as.numeric(row[(start+5+temp_start):(start+5+7+temp_stop)]))))
		probs <- as.numeric(row[6:length(row)])
		if (constant == FALSE) {
			full_left <- 27
			full_right <- 27+37
		} else {
			full_left <- 1
			full_right <- length(probs)
		}
		left <- max(start+1+temp_start,full_left)
		right <- min(start + 1 + temp_stop,full_right)
		return(10^mean(log10(1-probs[left:right])))
		})
}

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


get_site_probs_mean <- function(data,temp_start=0,temp_stop=7) {
	apply(data, 1, function(row) {
		start <- as.numeric(row[1])
		# print(as.numeric(row[(start+5+temp_start):(start+5+7+temp_stop)]))
		# print(1-as.numeric(row[(start+5+temp_start):(start+5+7+temp_stop)]))
		# print(10^sum(log10(1-as.numeric(row[(start+5+temp_start):(start+5+7+temp_stop)]))))
		probs <- as.numeric(row[6:length(row)])
		if (temp_start+start+1>=1&temp_stop+start+1<=length(probs)){
		return(10^mean(log10(1-probs[(start+1+temp_start):(start+1+temp_stop)])))
		} else {
			return(0)
		}

		})
}
# I_n0_7 <-  get_site_probs(s_I_8mer, 0,7)^(1/8)
# I_n1_7 <-  get_site_probs(s_I_8mer,-1,7)^(1/9)
# I_n2_7 <-  get_site_probs(s_I_8mer,-2,7)^(1/10)
# I_n3_7 <-  get_site_probs(s_I_8mer,-3,7)^(1/11)
# I_n4_7 <-  get_site_probs(s_I_8mer,-4,7)^(1/12)
# I_n5_7 <-  get_site_probs(s_I_8mer,-5,7)^(1/13)
# print("I 7")
# I_n0_8 <- get_site_probs(s_I_8mer, 0,8)^(1/9)
# I_n1_8 <- get_site_probs(s_I_8mer,-1,8)^(1/10)
# I_n2_8 <- get_site_probs(s_I_8mer,-2,8)^(1/11)
# I_n3_8 <- get_site_probs(s_I_8mer,-3,9)^(1/12)
# I_n4_8 <- get_site_probs(s_I_8mer,-4,9)^(1/13)
# I_n5_8 <- get_site_probs(s_I_8mer,-5,9)^(1/14)
# print("I 8")
# I_n0_9 <- get_site_probs(s_I_8mer, 0,9)^(1/10)
# I_n1_9 <- get_site_probs(s_I_8mer,-1,9)^(1/11)
# I_n2_9 <- get_site_probs(s_I_8mer,-2,9)^(1/12)
# I_n3_9 <- get_site_probs(s_I_8mer,-3,9)^(1/13)
# I_n4_9 <- get_site_probs(s_I_8mer,-4,9)^(1/14)
## NEW STUFF
# I_n5_9 <- get_site_probs(s_I_8mer,-5,9)
# I_n5_9_alt <- get_site_probs(s_I_8mer,-5,9,constant=FALSE)
# A1_n5_9 <- get_site_probs(s_0.4_8mer,-5,9)
# A1_n5_9_alt <- get_site_probs(s_0.4_8mer,-5,9,constant=FALSE)


# print("I 9")
# A1_n0_7 <-  get_site_probs(s_0.4_8mer, 0,7)^(1/8)
# A1_n1_7 <-  get_site_probs(s_0.4_8mer,-1,7)^(1/9)
# A1_n2_7 <-  get_site_probs(s_0.4_8mer,-2,7)^(1/10)
# A1_n3_7 <-  get_site_probs(s_0.4_8mer,-3,7)^(1/11)
# A1_n4_7 <-  get_site_probs(s_0.4_8mer,-4,7)^(1/12)
# A1_n5_7 <-  get_site_probs(s_0.4_8mer,-5,7)^(1/13)
# print("A1 7")
# A1_n0_8 <- get_site_probs(s_0.4_8mer, 0,8)^(1/9)
# A1_n1_8 <- get_site_probs(s_0.4_8mer,-1,8)^(1/10)
# A1_n2_8 <- get_site_probs(s_0.4_8mer,-2,8)^(1/11)
# A1_n3_8 <- get_site_probs(s_0.4_8mer,-3,9)^(1/12)
# A1_n4_8 <- get_site_probs(s_0.4_8mer,-4,9)^(1/13)
# A1_n5_8 <- get_site_probs(s_0.4_8mer,-5,9)^(1/14)
# print("A 8")
# A1_n0_9 <- get_site_probs(s_0.4_8mer, 0,9)^(1/10)
# A1_n1_9 <- get_site_probs(s_0.4_8mer,-1,9)^(1/11)
# A1_n2_9 <- get_site_probs(s_0.4_8mer,-2,9)^(1/12)
# A1_n3_9 <- get_site_probs(s_0.4_8mer,-3,9)^(1/13)
# A1_n4_9 <- get_site_probs(s_0.4_8mer,-4,9)^(1/14)
# A1_n5_9 <- get_site_probs(s_0.4_8mer,-5,9)^(1/15)
# print("A 9")

# A0_n0_7 <- get_site_probs(s_0_8mer, 0,7)^(1/8)
# A0_n1_7 <- get_site_probs(s_0_8mer,-1,7)^(1/9)
# A0_n2_7 <- get_site_probs(s_0_8mer,-2,7)^(1/10)
# A0_n3_7 <- get_site_probs(s_0_8mer,-3,7)^(1/11)
# A0_n4_7 <- get_site_probs(s_0_8mer,-4,7)^(1/12)
# A0_n5_7 <- get_site_probs(s_0_8mer,-5,7)^(1/13)
# A0_n6_7 <- get_site_probs(s_0_8mer,-6,7)^(1/14)
# A0_n7_7 <- get_site_probs(s_0_8mer,-7,7)^(1/15)
# A0_n8_7 <- get_site_probs(s_0_8mer,-8,7)^(1/16)
# A0_n0_8 <- get_site_probs(s_0_8mer, 0,8)^(1/9)
# A0_n1_8 <- get_site_probs(s_0_8mer,-1,8)^(1/10)
# A0_n2_8 <- get_site_probs(s_0_8mer,-2,8)^(1/11)
# A0_n3_8 <- get_site_probs(s_0_8mer,-3,8)^(1/12)
# A0_n4_8 <- get_site_probs(s_0_8mer,-4,8)^(1/13)
# A0_n5_8 <- get_site_probs(s_0_8mer,-5,8)^(1/14)
# A0_n6_8 <- get_site_probs(s_0_8mer,-6,8)^(1/15)
# A0_n7_8 <- get_site_probs(s_0_8mer,-7,8)^(1/16)
# A0_n8_8 <- get_site_probs(s_0_8mer,-8,8)^(1/17)
# A0_n0_9 <- get_site_probs(s_0_8mer, 0,9)^(1/10)
# A0_n1_9 <- get_site_probs(s_0_8mer,-1,9)^(1/11)
# A0_n2_9 <- get_site_probs(s_0_8mer,-2,9)^(1/12)
# A0_n3_9 <- get_site_probs(s_0_8mer,-3,9)^(1/13)
# A0_n4_9 <- get_site_probs(s_0_8mer,-4,9)^(1/14)
# A0_n5_9 <- get_site_probs(s_0_8mer,-5,9)^(1/15)
# A0_n6_9 <- get_site_probs(s_0_8mer,-6,9)^(1/16)
# A0_n7_9 <- get_site_probs(s_0_8mer,-7,9)^(1/17)
# A0_n8_9 <- get_site_probs(s_0_8mer,-8,9)^(1/18)



# I_n5_6 <-  get_site_probs(s_I_8mer,-5,6)^(1/12)
# A1_n5_6 <-  get_site_probs(s_0.4_8mer,-5,6)^(1/12)

# I_n6_7 <-  get_site_probs(s_I_8mer,-6,7)^(1/14)
# A1_n6_7 <-  get_site_probs(s_0.4_8mer,-6,7)^(1/14)

# I_n6_8 <-  get_site_probs(s_I_8mer,-6,8)^(1/15)
# A1_n6_8 <-  get_site_probs(s_0.4_8mer,-6,8)^(1/15)

# I_n6_9 <-  get_site_probs(s_I_8mer,-6,9)^(1/16)
# A1_n6_9 <-  get_site_probs(s_0.4_8mer,-6,9)^(1/16)

# I_n7_7 <-  get_site_probs(s_I_8mer,-7,7)^(1/15)
# A1_n7_7 <-  get_site_probs(s_0.4_8mer,-7,7)^(1/15)

# I_n7_8 <-  get_site_probs(s_I_8mer,-7,8)^(1/16)
# A1_n7_8 <-  get_site_probs(s_0.4_8mer,-7,8)^(1/16)

# I_n7_9 <-  get_site_probs(s_I_8mer,-7,9)^(1/17)
# A1_n7_9 <-  get_site_probs(s_0.4_8mer,-7,9)^(1/17)

# I_n8_7 <-  get_site_probs(s_I_8mer,-8,7)^(1/16)
# A1_n8_7 <-  get_site_probs(s_0.4_8mer,-8,7)^(1/16)

# I_n8_8 <-  get_site_probs(s_I_8mer,-8,8)^(1/17)
# A1_n8_8 <-  get_site_probs(s_0.4_8mer,-8,8)^(1/17)

# I_n8_9 <-  get_site_probs(s_I_8mer,-8,9)^(1/18)
# A1_n8_9 <-  get_site_probs(s_0.4_8mer,-8,9)^(1/18)
occ <- function(n,sat){
	return(1 / (1 + 1/(sat*(x^n))))
}

# total_flanks_I <- sapply(names(total_flanks_A),function(flank){
# 	left <- substr(flank,1,2)
# 	right <- substr(flank,3,4)
# 	print(left)
# 	print(right)
# 	return((rowSums(w_l)[left])*(rowSums(w_r)[right]))
# })
# total_flanks_I <- total_flanks_I / sum(total_flanks_I)
# names(total_flanks_I) <- names(total_flanks_A)
# med_0_7  <- median(I_n0_7) - median(A1_n0_7)
# med_n1_7 <- median(I_n1_7) - median(A1_n1_7)
# med_n2_7 <- median(I_n2_7) - median(A1_n2_7)
# med_n3_7 <- median(I_n3_7) - median(A1_n3_7)
# med_n4_7 <- median(I_n4_7) - median(A1_n4_7)
# med_n5_7 <- median(I_n5_7) - median(A1_n5_7)
# med_n6_7 <- median(I_n6_7) - median(A1_n6_7)
# med_0_8  <- median(I_n0_8) - median(A1_n0_8)
# med_n1_8 <- median(I_n1_8) - median(A1_n1_8)
# med_n2_8 <- median(I_n2_8) - median(A1_n2_8)
# med_n3_8 <- median(I_n3_8) - median(A1_n3_8)
# med_n4_8 <- median(I_n4_8) - median(A1_n4_8)
# med_n5_8 <- median(I_n5_8) - median(A1_n5_8)
# med_n6_8 <- median(I_n6_8) - median(A1_n6_8)
# med_0_9  <- median(I_n0_9) - median(A1_n0_9)
# med_n1_8 <- median(I_n1_8) - median(A1_n1_8)
# med_n1_9 <- median(I_n1_9) - median(A1_n1_9)
# med_n2_8 <- median(I_n2_8) - median(A1_n2_8)
# med_n2_9 <- median(I_n2_9) - median(A1_n2_9)
# med_n3_9 <- median(I_n3_9) - median(A1_n3_9)
# med_n4_9 <- median(I_n4_9) - median(A1_n4_9)
# med_n5_9 <- median(I_n5_9) - median(A1_n5_9)
# med_n6_9 <- median(I_n6_9) - median(A1_n6_9)


# mean_0_7  <- mean(I_n0_7) - mean(A1_n0_7)
# mean_n1_7 <- mean(I_n1_7) - mean(A1_n1_7)
# mean_n2_7 <- mean(I_n2_7) - mean(A1_n2_7)
# mean_n3_7 <- mean(I_n3_7) - mean(A1_n3_7)
# mean_n4_7 <- mean(I_n4_7) - mean(A1_n4_7)
# mean_n5_7 <- mean(I_n5_7) - mean(A1_n5_7)
# mean_n6_7 <- mean(I_n6_7) - mean(A1_n6_7)
# mean_0_8  <- mean(I_n0_8) - mean(A1_n0_8)
# mean_n1_8 <- mean(I_n1_8) - mean(A1_n1_8)
# mean_n2_8 <- mean(I_n2_8) - mean(A1_n2_8)
# mean_n3_8 <- mean(I_n3_8) - mean(A1_n3_8)
# mean_n4_8 <- mean(I_n4_8) - mean(A1_n4_8)
# mean_n5_8 <- mean(I_n5_8) - mean(A1_n5_8)
# mean_n6_8 <- mean(I_n6_8) - mean(A1_n6_8)
# mean_0_9  <- mean(I_n0_9) - mean(A1_n0_9)
# mean_n1_8 <- mean(I_n1_8) - mean(A1_n1_8)
# mean_n1_9 <- mean(I_n1_9) - mean(A1_n1_9)
# mean_n2_8 <- mean(I_n2_8) - mean(A1_n2_8)
# mean_n2_9 <- mean(I_n2_9) - mean(A1_n2_9)
# mean_n3_9 <- mean(I_n3_9) - mean(A1_n3_9)
# mean_n4_9 <- mean(I_n4_9) - mean(A1_n4_9)
# mean_n5_9 <- mean(I_n5_9) - mean(A1_n5_9)
# mean_n6_9 <- mean(I_n6_9) - mean(A1_n6_9)

# mean_n7_7 <- mean(I_n7_7) - mean(A1_n7_7)
# mean_n7_8 <- mean(I_n7_8) - mean(A1_n7_8)
# mean_n7_9 <- mean(I_n7_9) - mean(A1_n7_9)

# mean_n8_7 <- mean(I_n8_7) - mean(A1_n8_7)
# mean_n8_8 <- mean(I_n8_8) - mean(A1_n8_8)
# mean_n8_9 <- mean(I_n8_9) - mean(A1_n8_9)

# mean_n9_7 <- mean(I_n9_7) - mean(A1_n9_7)
# mean_n9_8 <- mean(I_n9_8) - mean(A1_n9_8)
# mean_n9_9 <- mean(I_n9_9) - mean(A1_n9_9)



# par(mfrow=c(3,9))

# plot(x,ecdf(I_n0_7)(x),type="l",lwd=2)
# text(mean(I_n0_7),y=1,labels=format(round(mean(I_0_7),3),nsmall=3))
# text((mean(I_n0_7)+mean(A1_n0_7))/2,0.9,labels=format(round(mean(I_n0_7)-mean(A1_n0_7),3),nsmall=3))
# lines(x,ecdf(A1_n0_7)(x),type="l",lwd=2,col="red");text(mean(A1_0_7),y=0.95,labels=format(round(mean(A1_n0_7),3),nsmall=3))
# #4
# plot(x,ecdf(I_n1_7)(x),type="l",lwd=2)
# text(mean(I_n1_7),y=1,labels=format(round(mean(I_n1_7),3),nsmall=3))
# text((mean(I_n1_7)+mean(A1_n1_7))/2,0.9,labels=format(round(mean(I_n1_7)-mean(A1_n1_7),3),nsmall=3))
# lines(x,ecdf(A1_n1_7)(x),type="l",lwd=2,col="red");text(mean(A1_n1_7),y=0.95,labels=format(round(mean(A1_n1_7),3),nsmall=3))
# #6
# plot(x,ecdf(I_n2_7)(x),type="l",lwd=2)
# text(mean(I_n2_7),y=1,labels=format(round(mean(I_n2_7),3),nsmall=3))
# text((mean(I_n2_7)+mean(A1_n2_7))/2,0.9,labels=format(round(mean(I_n2_7)-mean(A1_n2_7),3),nsmall=3))
# lines(x,ecdf(A1_n2_7)(x),type="l",lwd=2,col="red");text(mean(A1_n2_7),y=0.95,labels=format(round(mean(A1_n2_7),3),nsmall=3))
# #8
# plot(x,ecdf(I_n3_7)(x),type="l",lwd=2)
# text(mean(I_n3_7),y=1,labels=format(round(mean(I_n3_7),3),nsmall=3))
# text((mean(I_n3_7)+mean(A1_n3_7))/2,0.9,labels=format(round(mean(I_n3_7)-mean(A1_n3_7),3),nsmall=3))
# lines(x,ecdf(A1_n3_7)(x),type="l",lwd=2,col="red");text(mean(A1_n3_7),y=0.95,labels=format(round(mean(A1_n3_7),3),nsmall=3))
# #9
# plot(x,ecdf(I_n4_7)(x),type="l",lwd=2)
# text(mean(I_n4_7),y=1,labels=format(round(mean(I_n4_7),3),nsmall=3))
# text((mean(I_n4_7)+mean(A1_n4_7))/2,0.9,labels=format(round(mean(I_n4_7)-mean(A1_n4_7),3),nsmall=3))
# lines(x,ecdf(A1_n4_7)(x),type="l",lwd=2,col="red");text(mean(A1_n4_7),y=0.95,labels=format(round(mean(A1_n4_7),3),nsmall=3))
# #10=
# plot(x,ecdf(I_n5_7)(x),type="l",lwd=2)
# text(mean(I_n5_7),y=1,labels=format(round(mean(I_n5_7),3),nsmall=3))
# text((mean(I_n5_7)+mean(A1_n5_7))/2,0.9,labels=format(round(mean(I_n5_7)-mean(A1_n5_7),3),nsmall=3))
# lines(x,ecdf(A1_n5_7)(x),type="l",lwd=2,col="red");text(mean(A1_n5_7),y=0.95,labels=format(round(mean(A1_n5_7),3),nsmall=3))

# plot(x,ecdf(I_n6_7)(x),type="l",lwd=2)
# text(mean(I_n6_7),y=1,labels=format(round(mean(I_n6_7),3),nsmall=3))
# text((mean(I_n6_7)+mean(A1_n6_7))/2,0.9,labels=format(round(mean(I_n6_7)-mean(A1_n6_7),3),nsmall=3))
# lines(x,ecdf(A1_n6_7)(x),type="l",lwd=2,col="red");text(mean(A1_n6_7),y=0.95,labels=format(round(mean(A1_n6_7),3),nsmall=3))

# plot(x,ecdf(I_n7_7)(x),type="l",lwd=2)
# text(mean(I_n7_7),y=1,labels=format(round(mean(I_n7_7),3),nsmall=3))
# text((mean(I_n7_7)+mean(A1_n7_7))/2,0.9,labels=format(round(mean(I_n7_7)-mean(A1_n7_7),3),nsmall=3))
# lines(x,ecdf(A1_n7_7)(x),type="l",lwd=2,col="red");text(mean(A1_n7_7),y=0.95,labels=format(round(mean(A1_n7_7),3),nsmall=3))

# plot(x,ecdf(I_n8_7)(x),type="l",lwd=2)
# text(mean(I_n8_7),y=1,labels=format(round(mean(I_n8_7),3),nsmall=3))
# text((mean(I_n8_7)+mean(A1_n8_7))/2,0.9,labels=format(round(mean(I_n8_7)-mean(A1_n8_7),3),nsmall=3))
# lines(x,ecdf(A1_n8_7)(x),type="l",lwd=2,col="red");text(mean(A1_n8_7),y=0.95,labels=format(round(mean(A1_n8_7),3),nsmall=3))





# plot(x,ecdf(I_n0_8)(x),type="l",lwd=2)
# text(mean(I_n0_8),y=1,labels=format(round(mean(I_0_8),3),nsmall=3))
# text((mean(I_n0_8)+mean(A1_n0_8))/2,0.9,labels=format(round(mean(I_n0_8)-mean(A1_n0_8),3),nsmall=3))
# lines(x,ecdf(A1_n0_8)(x),type="l",lwd=2,col="red");text(mean(A1_0_8),y=0.95,labels=format(round(mean(A1_n0_8),3),nsmall=3))
# #4
# plot(x,ecdf(I_n1_8)(x),type="l",lwd=2)
# text(mean(I_n1_8),y=1,labels=format(round(mean(I_n1_8),3),nsmall=3))
# text((mean(I_n1_8)+mean(A1_n1_8))/2,0.9,labels=format(round(mean(I_n1_8)-mean(A1_n1_8),3),nsmall=3))
# lines(x,ecdf(A1_n1_8)(x),type="l",lwd=2,col="red");text(mean(A1_n1_8),y=0.95,labels=format(round(mean(A1_n1_8),3),nsmall=3))
# #6
# plot(x,ecdf(I_n2_8)(x),type="l",lwd=2)
# text(mean(I_n2_8),y=1,labels=format(round(mean(I_n2_8),3),nsmall=3))
# text((mean(I_n2_8)+mean(A1_n2_8))/2,0.9,labels=format(round(mean(I_n2_8)-mean(A1_n2_8),3),nsmall=3))
# lines(x,ecdf(A1_n2_8)(x),type="l",lwd=2,col="red");text(mean(A1_n2_8),y=0.95,labels=format(round(mean(A1_n2_8),3),nsmall=3))
# #8
# plot(x,ecdf(I_n3_8)(x),type="l",lwd=2)
# text(mean(I_n3_8),y=1,labels=format(round(mean(I_n3_8),3),nsmall=3))
# text((mean(I_n3_8)+mean(A1_n3_8))/2,0.9,labels=format(round(mean(I_n3_8)-mean(A1_n3_8),3),nsmall=3))
# lines(x,ecdf(A1_n3_8)(x),type="l",lwd=2,col="red");text(mean(A1_n3_8),y=0.95,labels=format(round(mean(A1_n3_8),3),nsmall=3))
# #9
# plot(x,ecdf(I_n4_8)(x),type="l",lwd=2)
# text(mean(I_n4_8),y=1,labels=format(round(mean(I_n4_8),3),nsmall=3))
# text((mean(I_n4_8)+mean(A1_n4_8))/2,0.9,labels=format(round(mean(I_n4_8)-mean(A1_n4_8),3),nsmall=3))
# lines(x,ecdf(A1_n4_8)(x),type="l",lwd=2,col="red");text(mean(A1_n4_8),y=0.95,labels=format(round(mean(A1_n4_8),3),nsmall=3))
# #10
# plot(x,ecdf(I_n5_8)(x),type="l",lwd=2)
# text(mean(I_n5_8),y=1,labels=format(round(mean(I_n5_8),3),nsmall=3))
# text((mean(I_n5_8)+mean(A1_n5_8))/2,0.9,labels=format(round(mean(I_n5_8)-mean(A1_n5_8),3),nsmall=3))
# lines(x,ecdf(A1_n5_8)(x),type="l",lwd=2,col="red");text(mean(A1_n5_8),y=0.95,labels=format(round(mean(A1_n5_8),3),nsmall=3))

# plot(x,ecdf(I_n6_8)(x),type="l",lwd=2)
# text(mean(I_n6_8),y=1,labels=format(round(mean(I_n6_8),3),nsmall=3))
# text((mean(I_n6_8)+mean(A1_n6_8))/2,0.9,labels=format(round(mean(I_n6_8)-mean(A1_n6_8),3),nsmall=3))
# lines(x,ecdf(A1_n6_8)(x),type="l",lwd=2,col="red");text(mean(A1_n6_8),y=0.95,labels=format(round(mean(A1_n6_8),3),nsmall=3))

# plot(x,ecdf(I_n6_8)(x),type="l",lwd=2)
# text(mean(I_n7_8),y=1,labels=format(round(mean(I_n7_8),3),nsmall=3))
# text((mean(I_n7_8)+mean(A1_n7_8))/2,0.9,labels=format(round(mean(I_n7_8)-mean(A1_n7_8),3),nsmall=3))
# lines(x,ecdf(A1_n7_8)(x),type="l",lwd=2,col="red");text(mean(A1_n7_8),y=0.95,labels=format(round(mean(A1_n7_8),3),nsmall=3))

# plot(x,ecdf(I_n8_8)(x),type="l",lwd=2)
# text(mean(I_n8_8),y=1,labels=format(round(mean(I_n8_8),3),nsmall=3))
# text((mean(I_n8_8)+mean(A1_n8_8))/2,0.9,labels=format(round(mean(I_n8_8)-mean(A1_n8_8),3),nsmall=3))
# lines(x,ecdf(A1_n8_8)(x),type="l",lwd=2,col="red");text(mean(A1_n8_8),y=0.95,labels=format(round(mean(A1_n8_8),3),nsmall=3))






# #3
# plot(x,ecdf(I_n0_9)(x),type="l",lwd=2)
# text(mean(I_n0_9),y=1,labels=format(round(mean(I_0_9),3),nsmall=3))
# text((mean(I_n0_9)+mean(A1_n0_9))/2,0.9,labels=format(round(mean(I_n0_9)-mean(A1_n0_9),3),nsmall=3))
# lines(x,ecdf(A1_n0_9)(x),type="l",lwd=2,col="red");text(mean(A1_0_9),y=0.95,labels=format(round(mean(A1_n0_9),3),nsmall=3))
# #4
# plot(x,ecdf(I_n1_9)(x),type="l",lwd=2)
# text(mean(I_n1_9),y=1,labels=format(round(mean(I_n1_9),3),nsmall=3))
# text((mean(I_n1_9)+mean(A1_n1_9))/2,0.9,labels=format(round(mean(I_n1_9)-mean(A1_n1_9),3),nsmall=3))
# lines(x,ecdf(A1_n1_9)(x),type="l",lwd=2,col="red");text(mean(A1_n1_9),y=0.95,labels=format(round(mean(A1_n1_9),3),nsmall=3))
# #6
# plot(x,ecdf(I_n2_9)(x),type="l",lwd=2)
# text(mean(I_n2_9),y=1,labels=format(round(mean(I_n2_9),3),nsmall=3))
# text((mean(I_n2_9)+mean(A1_n2_9))/2,0.9,labels=format(round(mean(I_n2_9)-mean(A1_n2_9),3),nsmall=3))
# lines(x,ecdf(A1_n2_9)(x),type="l",lwd=2,col="red");text(mean(A1_n2_9),y=0.95,labels=format(round(mean(A1_n2_9),3),nsmall=3))
# #8
# plot(x,ecdf(I_n3_9)(x),type="l",lwd=2)
# text(mean(I_n3_9),y=1,labels=format(round(mean(I_n3_9),3),nsmall=3))
# text((mean(I_n3_9)+mean(A1_n3_9))/2,0.9,labels=format(round(mean(I_n3_9)-mean(A1_n3_9),3),nsmall=3))
# lines(x,ecdf(A1_n3_9)(x),type="l",lwd=2,col="red");text(mean(A1_n3_9),y=0.95,labels=format(round(mean(A1_n3_9),3),nsmall=3))
# #9
# plot(x,ecdf(I_n4_9)(x),type="l",lwd=2)
# text(mean(I_n4_9),y=1,labels=format(round(mean(I_n4_9),3),nsmall=3))
# text((mean(I_n4_9)+mean(A1_n4_9))/2,0.9,labels=format(round(mean(I_n4_9)-mean(A1_n4_9),3),nsmall=3))
# lines(x,ecdf(A1_n4_9)(x),type="l",lwd=2,col="red");text(mean(A1_n4_9),y=0.95,labels=format(round(mean(A1_n4_9),3),nsmall=3))
# #10
# plot(x,ecdf(I_n5_9)(x),type="l",lwd=2)
# text(mean(I_n5_9),y=1,labels=format(round(mean(I_n5_9),3),nsmall=3))
# text((mean(I_n5_9)+mean(A1_n5_9))/2,0.9,labels=format(round(mean(I_n5_9)-mean(A1_n5_9),3),nsmall=3))
# lines(x,ecdf(A1_n5_9)(x),type="l",lwd=2,col="red");text(mean(A1_n5_9),y=0.95,labels=format(round(mean(A1_n5_9),3),nsmall=3))

# plot(x,ecdf(I_n6_9)(x),type="l",lwd=2)
# text(mean(I_n6_9),y=1,labels=format(round(mean(I_n6_9),3),nsmall=3))
# text((mean(I_n6_9)+mean(A1_n6_9))/2,0.9,labels=format(round(mean(I_n6_9)-mean(A1_n6_9),3),nsmall=3))
# lines(x,ecdf(A1_n6_9)(x),type="l",lwd=2,col="red");text(mean(A1_n6_9),y=0.95,labels=format(round(mean(A1_n6_9),3),nsmall=3))

# plot(x,ecdf(I_n7_9)(x),type="l",lwd=2)
# text(mean(I_n7_9),y=1,labels=format(round(mean(I_n7_9),3),nsmall=3))
# text((mean(I_n7_9)+mean(A1_n7_9))/2,0.9,labels=format(round(mean(I_n7_9)-mean(A1_n7_9),3),nsmall=3))
# lines(x,ecdf(A1_n7_9)(x),type="l",lwd=2,col="red");text(mean(A1_n7_9),y=0.95,labels=format(round(mean(A1_n7_9),3),nsmall=3))

# plot(x,ecdf(I_n8_9)(x),type="l",lwd=2)
# text(mean(I_n8_9),y=1,labels=format(round(mean(I_n8_9),3),nsmall=3))
# text((mean(I_n8_9)+mean(A1_n8_9))/2,0.9,labels=format(round(mean(I_n8_9)-mean(A1_n8_9),3),nsmall=3))
# lines(x,ecdf(A1_n8_9)(x),type="l",lwd=2,col="red");text(mean(A1_n8_9),y=0.95,labels=format(round(mean(A1_n8_9),3),nsmall=3))
col_2 <- rep(rainbow(1,start=0,end=1,alpha=0.6),256)
col_2[which(substr(names(total_flanks_A), 1, 2) == "CT" | substr(names(total_flanks_A), 1, 2) == "CC")] <- rainbow(1,start=0.6,end=1,alpha=0.9)
mean_exponent<- function(data,n){sapply(names(total_flanks_A),function(flank){mean((data^n)[which(flanks_temp==flank)])})}
total_flanks_I <- total_flanks_I_sub
get_flank_correlation_window <- function(data_I,data_A,flanks_data,n,sat){
	parameters <- rbind(sapply(names(total_flanks_A),function(flank){
	get_beta_params(data_I[which(flanks_I==flank)])}),total_flanks_I)
	input_beta_dist <- matrix(apply(parameters,2,function(col){col[3]*dbeta(x,col[1],col[2])}),nrow=dim(parameters)[2],byrow=TRUE)

	output_beta_dist <- matrix(apply(input_beta_dist,1,function(row){row*occ(n,sat)}),ncol=dim(input_beta_dist)[2],byrow=TRUE)
	par(mfrow=c(2,2))
	total_flanks <- sapply(names(total_flanks_A),function(flank){length(which(flanks_data==flank))})/length(flanks_data)
	inds <- which(substr(names(total_flanks_A), 1, 2) != "CT" & substr(names(total_flanks_A), 1, 2) != "CC")
	# inds <- seq(1,256)
	plot((rowSums(output_beta_dist)/sum(output_beta_dist))[inds],total_flanks[inds],log='xy',xlim=c(.00015,.02),ylim=c(.00015,.02),col=col_2[inds],lwd=3,cex=2)
	title(main=paste0("exponent: ",n),sub=sat)
	text(.0002,.015,labels=format(round(cor(log(total_flanks[inds]),log(rowSums(output_beta_dist)[inds])),4),nsmall=4))

	plot(x,ecdf(data_I)(x),type="l",lwd=2)
	lines(x,ecdf(data_A)(x),type="l",lwd=3,col="red")
	lines(x,ecdf(data_A[which(substr(flanks_data,1,2)!="CT" & substr(flanks_data,1,2)!="CC")])(x),type="l",lwd=2,col="purple")
	title(main=paste0("exponent: ",n),sub=sat)

	lines(x,ecdf_fun(colSums(output_beta_dist[inds,])),lwd=1,col="red")
	lines(x,ecdf_fun(colSums(output_beta_dist)),lwd=1,col="purple")

	# return(solution)
	I_means <- sapply(names(total_flanks_A),function(flank){mean(data_I[which(flanks_I==flank)])})
	A_means <- sapply(names(total_flanks_A),function(flank){mean(data_A[which(flanks_data==flank)])})

	A_predict <- apply(output_beta_dist,1,function(row){sum(x*row)/sum(row)})
	print(length(A_means))
	print(length(A_predict))
	plot(A_predict,I_means,xlim=c(0,1),ylim=c(0,1),col=col_2,lwd=3,cex=2)
	segments(0,0,1,1,lty=2)
	plot(A_predict,A_means,xlim=c(0,1),ylim=c(0,1),col=col_2,lwd=3,cex=2)
	segments(0,0,1,1,lty=2)


}
get_flank_correlation_window(I_n5_9,A1_n5_9,flanks_A,2.15,1)
total_flanks_I <- total_flanks_I_min


get_flank_correlation_window(I_n5_9_alt,A1_n5_9_alt,flanks_A,2.15,1)
