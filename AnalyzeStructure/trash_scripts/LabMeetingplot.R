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
graphics.off()

# s_I <- get_structures(mirna, site, "I")
# s_A <- get_structures(mirna, site, condition)

# s_I_8 <- get_site_probs(s_I)
# s_I_15<- get_site_probs(s_I, temp_stop = 14)

# s_A_8 <- get_site_probs(s_A)
# s_A_15<- get_site_probs(s_A, temp_stop = 14)

# flanks_I <- get_site_flanks(s_I)
# flanks_A <- get_site_flanks(s_A)

x_range <- c(0:1000)/1000

# dev.new(xpos = 20, ypos = 20 )
# par(par_plots)
# x_range <- c(0:1000)/1000
# plot(x_range,
# 	ecdf(s_I_8[flanks_I=="AAAA"])(x_range),
# 	type = "l",
# 	xlim = c(0, 1),
# 	ylim = c(0, 1),
# 	lwd = 2,
# 	axes = FALSE,
# 	ann = FALSE)
# lines(x_range,
# 	ecdf(s_A_8[flanks_A=="AAAA"])(x_range),
# 	lwd = 2,
# 	col = "red")
# axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), pos = 0, cex.axis = 1.2, lwd = 2)
# axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), pos = 0, cex.axis = 1.2, lwd = 2)

# title(main = "8mer target site nt 1–8", font.main = 1, cex.main = 1.2, adj = 0.5, line = 0)
# title(xlab = "Probability of target site window being unpaired", line = 1.5, cex.lab = 1.2)
# title(ylab = "CDF", line = 1.5, cex.lab = 1.2)


# dev.new(xpos = 1000, ypos = 20 )
# par(par_plots)
# x_range <- c(0:1000)/1000
# plot(x_range,
# 	ecdf(s_I_15[flanks_I=="AAAA"])(x_range),
# 	type = "l",
# 	xlim = c(0, 1),
# 	ylim = c(0, 1),
# 	lwd = 2,
# 	axes = FALSE,
# 	ann = FALSE)
# lines(x_range,
# 	ecdf(s_A_15[flanks_A=="AAAA"])(x_range),
# 	lwd = 2,
# 	col = "red")
# axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), pos = 0, cex.axis = 1.2, lwd = 2)
# axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), pos = 0, cex.axis = 1.2, lwd = 2)

# title(main = "8mer target site nt 1–15", font.main = 1, cex.main = 1.2, adj = 0.5, line = 0)
# title(xlab = "Probability of target site window being unpaired", line = 1.5, cex.lab = 1.2)
# title(ylab = "CDF", line = 1.5, cex.lab = 1.2)

# dev.new(xpos = 20, ypos = 500 )
# par(par_plots)
# x_range <- c(0:1000)/1000
# plot(x_range,
# 	ecdf((s_I_8^(1/8))[flanks_I=="AAAA"])(x_range),
# 	type = "l",
# 	xlim = c(0, 1),
# 	ylim = c(0, 1),
# 	lwd = 2,
# 	axes = FALSE,
# 	ann = FALSE)
# lines(x_range,
# 	ecdf((s_A_8^(1/8))[flanks_A=="AAAA"])(x_range),
# 	lwd = 2,
# 	col = "red")
# axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), pos = 0, cex.axis = 1.2, lwd = 2)
# axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), pos = 0, cex.axis = 1.2, lwd = 2)

# title(main = "8mer target site nt 1–8\ngeometric mean", font.main = 1, cex.main = 1.2, adj = 0.5, line = 0)
# title(xlab = "Per-nucleotide average probability of target site window being unpaired", line = 1.5, cex.lab = 1.2)
# title(ylab = "CDF", line = 1.5, cex.lab = 1.2)


dev.new(xpos = 20, ypos = 20 )
par(par_plots)
x_range <- c(0:1000)/1000
plot(x_range,
	ecdf((s_I_15^(1/15))[flanks_I=="AAAA"])(x_range),
	type = "l",
	xlim = c(0, 1),
	ylim = c(0, 1),
	lwd = 2,
	axes = FALSE,
	ann = FALSE)
lines(x_range,
	ecdf((s_A_15^(1/15))[flanks_A=="AAAA"])(x_range),
	lwd = 2,
	col = "red")
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), pos = 0, cex.axis = 1.2, lwd = 2)
axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), pos = 0, cex.axis = 1.2, lwd = 2)

title(main = "AA-8mer-AA target site nt 1–15\ngeometric mean", font.main = 1, cex.main = 1.2, adj = 0.5, line = 0)
title(xlab = "Per-nucleotide average probability of target site window being unpaired", line = 1.5, cex.lab = 1.2)
title(ylab = "CDF", line = 1.5, cex.lab = 1.2)

dev.new(xpos = 20, ypos = 500 )
par(par_plots)
x_range <- c(0:1000)/1000
plot(x_range,
	ecdf((s_I_15^(1/15))[flanks_I=="TTTT"])(x_range),
	type = "l",
	xlim = c(0, 1),
	ylim = c(0, 1),
	lwd = 2,
	axes = FALSE,
	ann = FALSE)
lines(x_range,
	ecdf((s_A_15^(1/15))[flanks_A=="TTTT"])(x_range),
	lwd = 2,
	col = "red")
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), pos = 0, cex.axis = 1.2, lwd = 2)
axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), pos = 0, cex.axis = 1.2, lwd = 2)

title(main = "TT-8mer-TT target site nt 1–15\ngeometric mean", font.main = 1, cex.main = 1.2, adj = 0.5, line = 0)
title(xlab = "Per-nucleotide average probability of target site window being unpaired", line = 1.5, cex.lab = 1.2)
title(ylab = "CDF", line = 1.5, cex.lab = 1.2)

dev.new(xpos = 1000, ypos = 20 )
par(par_plots)
x_range <- c(0:1000)/1000
plot(x_range,
	ecdf((s_I_15^(1/15))[flanks_I=="CCCC"])(x_range),
	type = "l",
	xlim = c(0, 1),
	ylim = c(0, 1),
	lwd = 2,
	axes = FALSE,
	ann = FALSE)
lines(x_range,
	ecdf((s_A_15^(1/15))[flanks_A=="CCCC"])(x_range),
	lwd = 2,
	col = "red")
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), pos = 0, cex.axis = 1.2, lwd = 2)
axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), pos = 0, cex.axis = 1.2, lwd = 2)

title(main = "CC-8mer-CC target site nt 1–15\ngeometric mean", font.main = 1, cex.main = 1.2, adj = 0.5, line = 0)
title(xlab = "Per-nucleotide average probability of target site window being unpaired", line = 1.5, cex.lab = 1.2)
title(ylab = "CDF", line = 1.5, cex.lab = 1.2)


dev.new(xpos = 1000, ypos = 500 )
par(par_plots)
x_range <- c(0:1000)/1000
plot(x_range,
	ecdf((s_I_15^(1/15))[flanks_I=="GGGG"])(x_range),
	type = "l",
	xlim = c(0, 1),
	ylim = c(0, 1),
	lwd = 2,
	axes = FALSE,
	ann = FALSE)
lines(x_range,
	ecdf((s_A_15^(1/15))[flanks_A=="GGGG"])(x_range),
	lwd = 2,
	col = "red")
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), pos = 0, cex.axis = 1.2, lwd = 2)
axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), pos = 0, cex.axis = 1.2, lwd = 2)

title(main = "GG-8mer-GG target site nt 1–15\ngeometric mean", font.main = 1, cex.main = 1.2, adj = 0.5, line = 0)
title(xlab = "Per-nucleotide average probability of target site window being unpaired", line = 1.5, cex.lab = 1.2)
title(ylab = "CDF", line = 1.5, cex.lab = 1.2)




dev.set(2)
dev.copy2pdf(file = "170402 AA8merAA_15ntgeo.pdf")

dev.set(3)
dev.copy2pdf(file = "170402 TT8merTT_15ntgeo.pdf")

dev.set(4)
dev.copy2pdf(file = "170402 CC8merCC_15ntgeo.pdf")

dev.set(5)
dev.copy2pdf(file = "170402 GG8merGG_15ntgeo.pdf")





