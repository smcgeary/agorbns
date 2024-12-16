library("Hmisc")
mirna <- "miR-1"
site <- "8mer"
dinucs <- apply(expand.grid(c("A", "C", "G", "T"),c("A","C","G","T")),
                 1, function(row) {
                 paste(rev(row),collapse = "")})

double_dinucs <- apply(expand.grid(dinucs,dinucs),
                 1, function(row) {
                 paste(rev(row),collapse = "")})

Get_structure_Probs <- function(condition, site, mirna) {
	out <- c()
	for (dinuc in double_dinucs) {
		name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",mirna,"/equilibrium/site_flank_unpaired_data/",site,"/",dinuc,"/",condition,"_-3--3_1-15.txt")
		testcon <- file(name,open="r")
		if (length(readLines(testcon)) > 0) {
			temp <- read.table(name)
			out <- c(out,temp[,1])
		}
		if (dinuc == "GTCG") {
			print(temp)
		}
		close(testcon)
	}
	return(out)
}
Get_structure_Probs_flank <- function(condition, site, mirna,flank) {
	out <- c()
		name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",mirna,"/equilibrium/site_flank_unpaired_data/",site,"/",flank,"/",condition,"_-3--3_1-15.txt")
		testcon <- file(name,open="r")
		if (length(readLines(testcon)) > 0) {
			out <- read.table(name)
		}
		close(testcon)
	return(out)
}

# Get_structure_Probs_average <- function(condition, site, mirna) {
# 	out <- rep(NA,length(double_dinucs))
# 	names(out) <- double_dinucs
# 	for (dinuc in double_dinucs) {
# 		name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/equilibrium/site_flank_unpaired_data/",site,"/",dinuc,"/",condition,"_-3--3_1-15.txt")
# 		testcon <- file(name,open="r")
# 		if (length(readLines(testcon)) > 0) {
# 			temp <- read.table(name)[,1]
# 			out[dinuc] <- mean(temp)
# 		}
# 		close(testcon)
# 	}
# 	return(out)
# }


# Get_flank_number <- function(condition, site, mirna) {
# 	out <- rep(NA,length(double_dinucs))
# 	names(out) <- double_dinucs
# 	for (dinuc in double_dinucs) {
# 		name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/equilibrium/site_flank_unpaired_data/",site,"/",dinuc,"/",condition,"_-3--3_1-15.txt")
# 		testcon <- file(name,open="r")
# 		out[dinuc] <- length(readLines(testcon))
# 		close(testcon)
# 	}
# 	return(out)
# }

# Get_AU_win <- function(condition, site, mirna) {
# 	out <- c()
# 	for (dinuc in double_dinucs) {
# 		name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/equilibrium/AU_win/",site,"/",dinuc,"/",condition,"_-3--3_1-15.txt")
# 		testcon <- file(name,open="r")
# 		if (length(readLines(testcon)) > 0) {
# 			temp <- read.table(name)[,1]
# 			out <- c(out,temp)
# 		}
# 		close(testcon)
# 	}
# 	return(out)
# }

# Get_AU_win_average <- function(condition, site, mirna) {
# 	out <- rep(NA,length(double_dinucs))
# 	names(out) <- double_dinucs
# 	for (dinuc in double_dinucs) {
# 		name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/equilibrium/AU_win/",site,"/",dinuc,"/",condition,"_-3--3_1-15.txt")
# 		testcon <- file(name,open="r")
# 		if (length(readLines(testcon)) > 0) {
# 			temp <- read.table(name)[,1]
# 			out[dinuc] <- mean(temp)
# 		}
# 		close(testcon)
# 	}
# 	return(out)
# }


# Get_AU_read <- function(condition, site, mirna) {
# 	out <- c()
# 	for (dinuc in double_dinucs) {
# 		name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/equilibrium/AU_read/",site,"/",dinuc,"/",condition,"_-3--3_1-15.txt")
# 		testcon <- file(name,open="r")
# 		if (length(readLines(testcon)) > 0) {
# 			temp <- read.table(name)[,1]
# 			out <- c(out,temp)
# 		}
# 		close(testcon)
# 	}
# 	return(out)
# }

# Get_AU_read_average <- function(condition, site, mirna) {
# 	out <- rep(NA,length(double_dinucs))
# 	names(out) <- double_dinucs
# 	for (dinuc in double_dinucs) {
# 		name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/equilibrium/AU_read/",site,"/",dinuc,"/",condition,"_-3--3_1-15.txt")
# 		testcon <- file(name,open="r")
# 		if (length(readLines(testcon)) > 0) {
# 			temp <- read.table(name)[,1]
# 			out[dinuc] <- mean(temp)
# 		}
# 		close(testcon)
# 	}
# 	return(out)
# }

# Get_AU_cs <- function(condition, site, mirna) {
# 	out <- c()
# 	for (dinuc in double_dinucs) {
# 		name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/equilibrium/AU_cs/",site,"/",dinuc,"/",condition,"_-3--3_1-15.txt")
# 		testcon <- file(name,open="r")
# 		if (length(readLines(testcon)) > 0) {
# 			temp <- read.table(name)[,1]
# 			out <- c(out,temp)
# 		}
# 		close(testcon)
# 	}
# 	return(out)
# }

flanks.I <- Get_flank_number("I", site, mirna)
probs.I <- Get_structure_Probs("I", site, mirna)
probs.I.byflank <- Get_structure_Probs_average("I", site, mirna)

AU_win.I <- Get_AU_win("I", site, mirna)
AU_read.I <- Get_AU_read("I", site, mirna)
AU_read.I.byflank <- Get_AU_read_average("I", site, mirna)

AU_cs.I <- Get_AU_cs("I", site, mirna)
flanks.A <- Get_flank_number("4", site, mirna)
probs.A <- Get_structure_Probs("4", site, mirna)
probs.A.byflank <- Get_structure_Probs_average("4", site, mirna)

AU_win.A <- Get_AU_win("4", site, mirna)
AU_read.A <- Get_AU_read("4", site, mirna)
AU_read.A.byflank <- Get_AU_read_average("4", site, mirna)
AU_cs.A <- Get_AU_cs("4", site, mirna)
ints.I <- sample(c(1:length(AU_read.I)), 100)



# dev.new(xpos = 20, ypos = 20, height = 10, width = 10)
# par(mfrow = c(2, 2))



# find_best_n <- function(n) {
# 	probs_sample_I <- sample(1:length(probs.I),
# 						 size = floor(length(probs.I)/10),
# 						 prob = probs.I^n)
# 	return((mean(probs.I[probs_sample_I]) - mean(probs.A))^2)
# }
# n <- optimize(find_best_n,c(-5, 5))$minimum

# probs_sample_I <- sample(1:length(probs.I),
# 						 size = floor(length(probs.I)/10),
# 						 prob = probs.I^n)


# ecdf_1 <- ecdf(probs.I)
# ecdf_2 <- ecdf(probs.A)
# ecdf_3 <- ecdf(probs.I[probs_sample_I])

# #PLOT 1 new
# x_range <- c(0:1000)/1000
# plot(x_range,
# 	ecdf_1(x_range),
# 	type = "l",
# 	xlim = c(0, 1),
# 	ylim = c(0, 1),
# 	lwd = 2,
# 	axes = FALSE,
# 	ann = FALSE)

# lines(x_range,
# 	ecdf_2(x_range),
# 	lwd = 2,
# 	col = "red")

# lines(x_range,
# 	ecdf_3(x_range),
# 	lwd = 2,
# 	col = "blue")


# axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), pos = 0, cex.axis = 1.2, lwd = 2)
# axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), pos = 0, cex.axis = 1.2, lwd = 2)
# title(ylab = "CDF", line = 1.5, cex.lab = 1.2)
# title(xlab = "Mean per-nucleotide probability of window across from\nmiRNA pos 1-15 being unpaired", line = 2.5, cex.lab = 1.2)

# legend("topleft",
# 	legend = c("Input", "0.4% AGO", "Input, matched\nfor unpaired probability"),
# 	lwd = 2,
# 	col = c("black", "red", "blue"),
# 	bty = "n",
# 	cex  = 1.2)



# ecdf_1 <- ecdf(AU_cs.I)
# ecdf_2 <- ecdf(AU_cs.A)
# ecdf_3 <- ecdf(AU_cs.I[probs_sample_I])

# effect <- (mean(AU_cs.I[probs_sample_I]) - mean(AU_cs.I)) / (mean(AU_cs.A) - mean(AU_cs.I))

# #PLOT 2 new
# x_range <- c(0:7000)/1000
# plot(x_range,
# 	ecdf_1(x_range),
# 	type = "l",
# 	xlim = c(2, 6),
# 	ylim = c(0, 1),
# 	lwd = 2,
# 	axes = FALSE,
# 	ann = FALSE)

# lines(x_range,
# 	ecdf_2(x_range),
# 	lwd = 2,
# 	col = "red")

# lines(x_range,
# 	ecdf_3(x_range),
# 	lwd = 2,
# 	col = "blue")


# axis(1, at = c(2, 3, 4, 5, 6), pos = 0, cex.axis = 1.2, lwd = 2)
# axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), pos = 2, cex.axis = 1.2, lwd = 2)
# title(ylab = "CDF", line = 1.5, cex.lab = 1.2)
# title(xlab = "AU context score", line = 1.5, cex.lab = 1.2)
# text(x = 5, y = 0.1, paste0("Effect: ",round(effect, 2)))




# find_best_n <- function(n) {
# 	sample_I <- sample(1:length(probs.I),
# 						 size = floor(length(probs.I)/10),
# 						 prob = AU_cs.I^n)
# 	return((mean(AU_cs.I[sample_I]) - mean(AU_cs.A))^2)
# }
# n <- optimize(find_best_n,c(1, 7))$minimum

# sample_I <- sample(1:length(probs.I),
# 						 size = floor(length(probs.I)/10),
# 						 prob = AU_cs.I^n)


# ecdf_1 <- ecdf(probs.I)
# ecdf_2 <- ecdf(probs.A)
# ecdf_3 <- ecdf(probs.I[sample_I])

# effect <- (mean(probs.I[sample_I]) - mean(probs.I)) / (mean(probs.A) - mean(probs.I))


# #PLOT 3 new
# x_range <- c(0:1000)/1000
# plot(x_range,
# 	ecdf_1(x_range),
# 	type = "l",
# 	xlim = c(0, 1),
# 	ylim = c(0, 1),
# 	lwd = 2,
# 	axes = FALSE,
# 	ann = FALSE)

# lines(x_range,
# 	ecdf_2(x_range),
# 	lwd = 2,
# 	col = "red")

# lines(x_range,
# 	ecdf_3(x_range),
# 	lwd = 2,
# 	col = "forestgreen")


# axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), pos = 0, cex.axis = 1.2, lwd = 2)
# axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), pos = 0, cex.axis = 1.2, lwd = 2)
# title(ylab = "CDF", line = 1.5, cex.lab = 1.2)
# title(xlab = "Mean per-nucleotide probability of window across from\nmiRNA pos 1-15 being unpaired", line = 2.5, cex.lab = 1.2)

# legend("topleft",
# 	legend = c("Input", "0.4% AGO", "Input, matched\nfor AU context score"),
# 	lwd = 2,
# 	col = c("black", "red", "forestgreen"),
# 	bty = "n",
# 	cex  = 1.2)

# text(x = 0.7, y = 0.1, paste0("Effect: ",round(effect,2)))


# ecdf_1 <- ecdf(AU_read.I)
# ecdf_2 <- ecdf(AU_read.A)
# ecdf_3 <- ecdf(AU_read.I[sample_I])



# #PLOT 4 new
# x_range <- c(0:7000)/1000
# plot(x_range,
# 	ecdf_1(x_range),
# 	type = "l",
# 	xlim = c(2, 6),
# 	ylim = c(0, 1),
# 	lwd = 2,
# 	axes = FALSE,
# 	ann = FALSE)

# lines(x_range,
# 	ecdf_2(x_range),
# 	lwd = 2,
# 	col = "red")

# lines(x_range,
# 	ecdf_3(x_range),
# 	lwd = 2,
# 	col = "forestgreen")



# axis(1, at = c(2, 3, 4, 5, 6), pos = 0, cex.axis = 1.2, lwd = 2)
# axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), pos = 2, cex.axis = 1.2, lwd = 2)
# title(ylab = "CDF", line = 1.5, cex.lab = 1.2)
# title(xlab = "AU context score", line = 1.5, cex.lab = 1.2)








IndexMatrix <- function(probs_new,AU_new) {
	probs <- as.integer(100*probs_new)
	AU <- as.integer(15*AU_new)
	out <- matrix(0,nrow = 100, ncol = 15)
	pairs <- cbind(probs,AU)
	for (i in 1:nrow(pairs)) {
		row <- pairs[i, 1]
		col <- pairs[i, 2]
		out[row,col] <- out[row,col] + 1
	}
	return(out)
}
frame_I <- data.frame(cbind(as.integer(15*AU_win.I), as.integer(37*AU_read.I), probs.I))
frame_A <- data.frame(cbind(as.integer(15*AU_win.A), as.integer(37*AU_read.A), probs.A))

frame_A_sort <- frame_A[order(frame_A[, 1]), ]

cutoff_inds <- seq(1, nrow(frame_A), length = 5)

cutoffs <- frame_A_sort[cutoff_inds,1]

colnames(frame_I) <- c("AU_win","AU_read", "prob")
colnames(frame_A) <- colnames(frame_I)

# bins <- sapply(frame_A_sort[,1], function(x) {
# 	length(which(cutoffs > = )
# 	})

ag_I <- aggregate(. ~ AU_win, frame_I, function(x) c( mean = mean(x), se = sd(x)/sqrt(length(x)-1)))
ag_A <- aggregate(. ~ AU_win, frame_A, function(x) c( mean = mean(x), se = sd(x)/sqrt(length(x)-1)))

par(mfrow = c(2,3))
##PLOT 1
plot(ag_I$AU_win/15, ag_I$prob[,1], xlim = c(0,1), ylim = c(0,1), type = "o", lwd = 2, axes = FALSE, ann = FALSE)
title(main = paste0(mirna, "\n", site, collapse=""), font.main = 1, cex.main = 2)
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), pos = 0, lwd = 2, cex.axis = 1.5)
title(xlab = "AU content across from miRNA pos 1-15", cex.lab = 1.5, line = 1.5)
axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), pos = 0, lwd = 2, cex.axis = 1.5)
title(ylab = "Mean per-nucleotide probability of window across from\nmiRNA pos 1-15 being unpaired", line = 1.5, cex.lab = 1.5)


with (
  data = ag_I
  , expr = errbar(AU_win/15, prob[,1], prob[,1]+prob[,2], prob[,1]-prob[,2], add=T, pch=1, lwd = 2, cap=.1)
)

lines(ag_A$AU_win/15, ag_A$prob[,1], type = "o", col = "red", lwd = 2)
par(fg = "red")
with (
  data = ag_A
  , expr = errbar(AU_win/15, prob[,1], prob[,1]+prob[,2], prob[,1]-prob[,2], add=T, pch=1, lwd = 2, cap=.1)
)
par(fg = "black")
par(fg = "black")

ag_I <- aggregate(. ~ AU_read, frame_I, function(x) c( mean = mean(x), se = sd(x)/sqrt(length(x)-1)))
ag_A <- aggregate(. ~ AU_read, frame_A, function(x) c( mean = mean(x), se = sd(x)/sqrt(length(x)-1)))


##PLOT 2
plot(ag_I$AU_read/37, ag_I$prob[,1], xlim = c(0,1), ylim = c(0,1), type = "o", lwd = 2, axes = FALSE, ann = FALSE)
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), pos = 0, lwd = 2, cex.axis = 1.5)
title(xlab = "AU content within entire randomized region", line = 1.5, cex.lab = 1.5)
axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), pos = 0, lwd = 2, cex.axis = 1.5)
title(ylab = "Mean per-nucleotide probability of window across from\nmiRNA pos 1-15 being unpaired", line = 1.5, cex.lab = 1.5)
with (
  data = ag_I
  , expr = errbar(AU_read/37, prob[,1], prob[,1]+prob[,2], prob[,1]-prob[,2], add=T, pch=1, lwd = 2, cap=.1)
)
lines(ag_A$AU_read/37, ag_A$prob[,1], type = "o", col = "red", lwd = 2)
par(fg = "red")
with (
  data = ag_A
  , expr = errbar(AU_read/37, prob[,1], prob[,1]+prob[,2], prob[,1]-prob[,2], add=T, pch=1, lwd = 2, cap=.1)
)
par(fg = "black")

ecdf_1 <- ecdf(frame_I$prob)
ecdf_2 <- ecdf(frame_A$prob)

#PLOT 3
x_range <- c(0:1000)/1000
plot(x_range,
	ecdf_1(x_range),
	type = "l",
	xlim = c(0, 1),
	ylim = c(0, 1),
	lwd = 2,
	axes = FALSE,
	ann = FALSE)

lines(x_range,
	ecdf_2(x_range),
	lwd = 2,
	col = "red")

legend("right",
	legend = c("Input", "0.4% AGO"),
	lwd = 2,
	col = c("black", "red"),
	bty = "n",
	cex  = 2.0)

axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), pos = 0, cex.axis = 1.5, lwd = 2)
axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), pos = 0, cex.axis = 1.5, lwd = 2)
title(ylab = "CDF", line = 1.5, cex.lab = 1.5)
title(xlab = "Mean per-nucleotide probability of window cross from\nmiRNA pos 1-15 being unpaired", line = 3, cex.lab = 1.5)


frame_I_new <- data.frame(cbind(AU_win.I, AU_read.I, round(probs.I, 1)))
frame_A_new <- data.frame(cbind(AU_win.A, AU_read.A, round(probs.A, 1)))

colnames(frame_I_new) <- c("AU_win","AU_read", "prob")
colnames(frame_A_new) <- colnames(frame_I)

ag_I <- aggregate(. ~ prob, frame_I_new, function(x) c( mean = mean(x), se = sd(x)/sqrt(length(x)-1), len = length(x)))
ag_A <- aggregate(. ~ prob, frame_A_new, function(x) c( mean = mean(x), se = sd(x)/sqrt(length(x)-1), len = length(x)))


#PLOT 4
plot(ag_I$prob, ag_I$AU_win[,1], xlim = c(0,1), ylim = c(0,1), type = "o", lwd = 2, axes = FALSE, ann = FALSE)
with (
  data = ag_I
  , expr = errbar(prob, AU_win[,1], AU_win[,1] + AU_win[,2], AU_win[, 1] - AU_win[, 2], add=T, pch=1, lwd = 2, cap=.1)
)
lines(ag_A$prob, ag_A$AU_win[,1], type = "o", col = "red", lwd = 2)
par(fg = "red")
with (
  data = ag_A
  , expr = errbar(prob, AU_win[, 1], AU_win[, 1]+AU_win[, 2], AU_win[,1]-AU_win[, 2], add=T, pch=1, lwd = 2, cap=.1)
)
par(fg = "black")
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), pos = 0, lwd = 2, cex.axis = 1.5)
title(xlab = "Mean per-nucleotide probability of window across from\nmiRNA pos 1-15 being unpaired", line = 3, cex.lab = 1.5)
axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), pos = 0, lwd = 2, cex.axis = 1.5)
title(ylab = "AU content across from miRNA pos 1-15", line = 1.5, cex.lab = 1.5)

#PLOT 5

plot(ag_I$prob, ag_I$AU_read[,1], xlim = c(0,1), ylim = c(0,1), type = "o", lwd = 2, axes = FALSE, ann = FALSE)
with (
  data = ag_I
  , expr = errbar(prob, AU_read[,1], AU_read[,1] + AU_read[,2], AU_read[, 1] - AU_read[, 2], add=T, pch=1, lwd = 2, cap=.1)
)
lines(ag_A$prob, ag_A$AU_read[,1], type = "o", col = "red", lwd = 2)
par(fg = "red")
with (
  data = ag_A
  , expr = errbar(prob, AU_read[, 1], AU_read[, 1]+AU_read[, 2], AU_read[,1]-AU_read[, 2], add=T, pch=1, lwd = 2, cap=.1)
)

par(fg = "black")
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), pos = 0, lwd = 2, cex.axis = 1.5)
title(xlab = "Mean per-nucleotide probability of window across from\nmiRNA pos 1-15 being unpaired", line = 3, cex.lab = 1.5)
axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), pos = 0, lwd = 2, cex.axis = 1.5)
title(ylab = "AU content within entire randomized sequence", cex.lab = 1.5)

#PLOT 6


ecdf_1 <- ecdf(frame_I_new$AU_win)
ecdf_2 <- ecdf(frame_A_new$AU_win)
plot(x_range,
	ecdf_1(x_range),
	xlim = c(0, 1),
	ylim = c(0, 1),
	lwd = 2,
	type = "l",
	axes = FALSE,
	ann = FALSE)

title(xlab = "AU content across from miRNA pos 1-15", line = 1.5, cex.lab = 1.5)
title(ylab = "CDF", line = 1.5, cex.lab = 1.5)
# plot(ecdf(frame_I_new$AU_read), add = TRUE)
lines(x_range,
	ecdf_2(x_range),
	lwd = 2,
	col = "red")
	# plot(ecdf(frame_A_new$AU_read),add = TRUE, col = "red")


axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), pos = 0, lwd = 2, cex.axis = 1.5)
axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), pos = 0, lwd = 2, cex.axis = 1.5)


# par(mfrow = c(1,2))
# plot(AU_win.I[ints.I], probs.I[ints.I],
# 	axes = FALSE,
# 	ann = FALSE)
# axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), pos = 0, lwd = 2, cex.axis = 1.5)
# title(xlab = "Mean per-nucleotide probability of window across from\nmiRNA pos 1-15 being unpaired", line = 2, cex.lab = 1.5)
# axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), pos = 0, lwd = 2, cex.axis = 1.5)
# title(ylab = "AU content within entire randomized sequence", cex.lab = 1.5)

# plot(AU_win[ints], probs[ints])
