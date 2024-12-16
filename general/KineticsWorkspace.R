source("general/general.R")

graphics.off()

MakeRepressionLine <- function(dwelltime, occupancy) {
	graphics.off()

	dev.new(xpos=1020, ypos=20, height=5, width=5)
	par(kPlotParameters)
	xmin <- 0
	xmax <- 4
	ymin <- 0
	ymax <- 30

	BlankPlot()
	AddLinearAxis(1, tick.space=1, label.space=1, label="Time (h)")
	AddLinearAxis(2, tick.space=1, label.space=5, label="Poly(A) tail length")
	xpoints <- seq(0, 4)
	ypoints <- seq(24, 0, by=-6) + 2
	ypoints2 <- c(seq(24, 0, by=-8) + 2, 0)

	xpoints_dense <- seq(0, 4, length.out=4*60)
	ypoints_dense <- rep(26, 100)
	ypoints_dense[1] <- 26

	slope_unbound <- -6/60
	slope_bound_app <- -8/60
	slope_bound_del <- slope_bound_app - slope_unbound
	print(slope_bound_del)
	slope_bound_enrich <- slope_bound_del/occupancy + slope_unbound
	print(slope_bound_enrich)
	print(slope_unbound)
	slope_inactive <- slope_unbound
	slope_active <- -6/60
	offtime <- dwelltime/occupancy - dwelltime
	print("offtime")
	print(offtime)
	bound <- TRUE
	total_time <- 4*60
	time <- 0
	transition_time <- dwelltime
	waiting_time <- offtime
	x_l <- 0
	for (i in seq(2, total_time)) {
		# print(i)
		# print(ypoints_dense[i])
		if (bound) {
			ypoints_dense[i] <- ypoints_dense[i - 1] + slope_bound_enrich
		} else {
			ypoints_dense[i] <- ypoints_dense[i - 1] + slope_unbound					
		}
		time <- time + 1
		if (time == transition_time) {
			time <- 0
			if (bound) {
				bound <- FALSE
				transition_time <- offtime
				rect(x_l, 0, i/60, 30, lwd=NA, col="gray90")
			} else {
				bound <- TRUE
				transition_time <- dwelltime
				x_l <- i/60
			}
		}
		if (i%%30 == 0) {
			slope_temp <- slope_active
			slope_active <- slope_inactive
			slope_inactive <- slope_temp
		}
	}

	points(xpoints, ypoints, pch=20, type="o", col="blue")
	points(xpoints, ypoints2, pch=20, type="o")
	lines(xpoints_dense, ypoints_dense, col="red")
	dev.copy2pdf(file=file.path("/lab/bartel1_ata/mcgeary/Presentations/seminars/",
                            "180604 Group Meeting/plot_", occupancy, "_",
                            dwelltime, ".pdf", fsep=""))


}


MakeRepressionLineDeadTime <- function(dwelltime, occupancy, deadtime) {
	graphics.off()

	dev.new(xpos=1020, ypos=20, height=5, width=5)
	par(kPlotParameters)
	xmin <- 0
	xmax <- 4
	ymin <- 0
	ymax <- 30

	BlankPlot()
	AddLinearAxis(1, tick.space=1, label.space=1, label="Time (h)")
	AddLinearAxis(2, tick.space=1, label.space=5, label="Poly(A) tail length")
	xpoints <- seq(0, 4)
	ypoints <- seq(24, 0, by=-6) + 2
	ypoints2 <- c(seq(24, 0, by=-8) + 2, 0)

	xpoints_dense <- seq(0, 4, length.out=4*60)
	ypoints_dense <- rep(26, 100)
	ypoints_dense[1] <- 26


	slope_unbound <- -6/60
	slope_bound_app <- -8/60
	slope_bound_del <- slope_bound_app - slope_unbound
	print(slope_bound_del)
	slope_bound_enrich <- slope_bound_del/occupancy + slope_unbound
	print(slope_bound_enrich)
	print(slope_unbound)
	slope_inactive <- slope_unbound
	slope_active <- -6/60
	offtime <- dwelltime/occupancy - dwelltime
	print("offtime")
	print(offtime)
	bound <- TRUE
	total_time <- 4*60
	time <- 0
	transition_time <- dwelltime
	waiting_time <- offtime
	x_l <- 0
	for (i in seq(2, total_time)) {
		# print(i)
		# print(ypoints_dense[i])
		if (bound & time >= deadtime) {
			ypoints_dense[i] <- ypoints_dense[i - 1] + slope_bound_enrich
		} else {
			ypoints_dense[i] <- ypoints_dense[i - 1] + slope_unbound					
		}
		time <- time + 1
		if (time == transition_time) {
			time <- 0
			if (bound) {
				bound <- FALSE
				transition_time <- offtime
				rect(x_l, 0, i/60, 30, lwd=NA, col="gray90")
			} else {
				bound <- TRUE
				dead <- TRUE
				transition_time <- dwelltime
				x_l <- i/60
			}
		}
		if (i%%30 == 0) {
			slope_temp <- slope_active
			slope_active <- slope_inactive
			slope_inactive <- slope_temp
		}
	}

	points(xpoints, ypoints, pch=20, type="o", col="blue")
	points(xpoints, ypoints2, pch=20, type="o")
	lines(xpoints_dense, ypoints_dense, col="red")
	dev.copy2pdf(file=file.path("/lab/bartel1_ata/mcgeary/Presentations/seminars/",
                            "180604 Group Meeting/plot_", occupancy, "_",
                            dwelltime, "_dead", deadtime,".pdf", fsep=""))


}

SimulateLog2FC <- function(dwelltime, kd, deadtime, a=0.05) {
	xpoints <- seq(0, 4)
	ypoints <- seq(24, 0, by=-6) + 2
	ypoints2 <- c(seq(24, 0, by=-8) + 2, 0)

	xpoints_dense <- seq(0, 4, length.out=4*60)
	ypoints_dense <- rep(26, 100)
	ypoints_dense[1] <- 26


	occupancy <- a / (a + kd)

	slope_unbound <- -6/60
	slope_bound <- -10/60


	active_bound_time <- max(c(dwelltime - deadtime)/dwelltime, 0)

	active_occupancy <- active_bound_time * occupancy

	apparent_rate <- active_occupancy*slope_bound + (1 - active_occupancy)* slope_unbound

	print(apparent_rate)

	time_to_lose_tail <- -24/apparent_rate/60

	print(time_to_lose_tail)

	log2fc <- log(time_to_lose_tail/4, 2)

	return(log2fc)
}

random_dwelltime <- exp(runif(100, min=log(0.1), max=log(30)))
random_kds <- exp(runif(100, min=log(0.001), max=log(3)))

categories_dwelltime <- log(random_dwelltime) - min(log(random_dwelltime))
categories_dwelltime <- ceiling(categories_dwelltime/max(categories_dwelltime)*99 + 1)

categories_kds <- log(random_kds) - min(log(random_kds))
categories_kds <- ceiling(categories_kds/max(categories_kds)*99 + 1)

temp <- SimulateLog2FC(10, 0.5, 2)

random_matrix <- cbind(random_dwelltime, random_kds)

PlotScatterWithKdAndKoff <- function(deadtime) {
	colors <- rainbow(100, start=0.4, end=0.9)

	fcs <- apply(random_matrix, 1, function(row) {
		SimulateLog2FC(row[1], row[2] , deadtime)
	})
	data_all <- cbind(random_matrix, fcs)

	dev.new(xpos=20, ypos=20, height=5, width=5)
	par(kPlotParameters)
	xmin <- 1e-3
	xmax <- 3
	ymin <- -1
	ymax <- 0.2
	BlankPlot(log="x", inv="x")
	AddLogAxis(1, label="Kd")
	AddLinearAxis(2, tick.space=0.1, label.space=0.2, label="Log2 fold-change")
	points(data_all[, 2], data_all[, 3], col=colors[categories_dwelltime], pch=19)
	dev.copy2pdf(file=file.path("/lab/bartel1_ata/mcgeary/Presentations/seminars/",
                            "180604 Group Meeting/KdVsRep_dead", deadtime,".pdf", fsep=""))


	dev.new(xpos=520, ypos=20, height=5, width=5)
	par(kPlotParameters)
	xmin <- 0.1
	xmax <- 30
	ymin <- -1
	ymax <- 0.2
	BlankPlot(log="x")
	AddLogAxis(1, label="Dwelltime")
	AddLinearAxis(2, tick.space=0.1, label.space=0.2, label="Log2 fold-change")
	points(data_all[, 1], data_all[, 3], col=rev(colors)[categories_kds], pch=19)
	dev.copy2pdf(file=file.path("/lab/bartel1_ata/mcgeary/Presentations/seminars/",
                            "180604 Group Meeting/DwelltimeVsRep_dead", deadtime,".pdf", fsep=""))



	fcs_categorical <- -1*ceiling(fcs*100) + 1
	colors_fc <- c("gray", rev(rainbow(99, start=0.4, end=0.1)))

	dev.new(xpos=1020, ypos=20, height=5, width=5)
	par(kPlotParameters)
	ymin <- 1e-3
	ymax <- 3
	xmin <- 0.1
	xmax <- 30







	BlankPlot(log="xy", inv="y")
	AddLogAxis(2, label="Relative Kd")
	AddLogAxis(1, label="Dwelltime (min)")
	points(data_all[,1], data_all[, 2], col=colors_fc[fcs_categorical], pch=19)
	dev.copy2pdf(file=file.path("/lab/bartel1_ata/mcgeary/Presentations/seminars/",
                            "180604 Group Meeting/DwelltimeVsKd_dead", deadtime,".pdf", fsep=""))
}

repression_matrix <- data.frame(fread(sprintf("/lab/solexa_bartel/mcgeary/AgoRBNS/%s/repression_hela_cs/lin_model_df/%s.txt",
                                     "miR-1", "paper"), sep="\t"), row.names=1)

data_cb <- read.table("HeLa_half_lives.txt", sep="\t", header=FALSE,
                      row.names=1, stringsAsFactors=FALSE)
data_cb_conversion <- read.table("Metadata_HS_Refseq_to_Ensembl.txt", sep="\t",
                                 stringsAsFactors=FALSE)
names_cb_refseq <- rownames(data_cb)
ensembl_ids <- data_cb_conversion[, 2]
names(ensembl_ids) <- data_cb_conversion[, 3]
halflife_cb <- unlist(data_cb)
names(halflife_cb) <- ensembl_ids[rownames(data_cb)]

halflife_cb <- halflife_cb[!is.na(names(halflife_cb))]

data_tani <- read.table("Tani_Supp_Tables_revised2.txt", skip=3, sep="\t",
                        header=TRUE, row.names=1, stringsAsFactors=FALSE)
halflife_tani <- sapply(data_tani[, 3], function(x) {
	as.numeric(substr(x, 1, nchar(x) - 1))
})

names(halflife_tani) <- sapply(rownames(data_tani), function(x) {
	names_all <- strsplit(x, split=",")[[1]]
	for (name in names_all) {
		ind <- which(rownames(repression_matrix) == name)
		ind2 <- which(names(halflife_cb) == name)
		if (length(ind) != 0) {
			return(rownames(repression_matrix)[ind])
		} else if (length(ind2) != 0) {
			return(names(halflife_cb)[ind2])
		}
	}
	return(paste0(names_all, collapse="/"))
})

names_all <- union(names(halflife_tani), names(halflife_cb))

halflife_all <- matrix(NaN, nrow=length(names_all), ncol=2)

rownames(halflife_all) <- names_all

halflife_all[names(halflife_cb), 1] <-  halflife_cb[names(halflife_cb)]
halflife_all[names(halflife_tani), 2] <-  halflife_tani[names(halflife_tani)]

dev.new(xpos=10, ypos=20, height=4, width=4)
par(kPlotParameters)
xmin <- 0.1
xmax <- 100
ymin <- 0.1
ymax <- 100
BlankPlot(log='xy')
AddLogAxis(1, label="Half-life (CB)")
AddLogAxis(2, label="Half-life (Tani et al)")
points(halflife_all[, 1], halflife_all[, 2], pch=1, col=rgb(0, 0, 0, alpha=0.2))
abline(0, 1, lty=2)
xy <- GetPlotFractionalCoords(fx=0.05, fy=0.95, log='xy')

AddCorrelationToPlot(log(halflife_all[, 1]), log(halflife_all[, 2]), xy[1], xy[2])

xy <- GetPlotFractionalCoords(fx=0.05, fy=0.90, log='xy')
AddCorrelationToPlot(log(halflife_all[, 1]), log(halflife_all[, 2]), 
                     xy[1], xy[2], method="spearman")



dev.copy2pdf(file=file.path("/lab/bartel1_ata/mcgeary/Presentations/seminars/",
                            "180604 Group Meeting/Correlation_ALL.pdf", fsep=""))


names_intersect <- intersect(rownames(halflife_all), rownames(repression_matrix))
halflife_use <- halflife_all[names_intersect, ]
print(dim(repression_matrix))
repression_matrix <- repression_matrix[names_intersect, ]
print(dim(repression_matrix))
print(dim(halflife_use))

dev.new(xpos=400, ypos=20, height=4, width=4)
par(kPlotParameters)
xmin <- 0.1
xmax <- 100
ymin <- 0.1
ymax <- 100
BlankPlot(log='xy')
AddLogAxis(1, label="Half-life (CB)")
AddLogAxis(2, label="Half-life (Tani et al)")
points(halflife_use[, 1], halflife_use[, 2], pch=1, col=rgb(0, 0, 0, alpha=0.2))
abline(0, 1, lty=2)
xy <- GetPlotFractionalCoords(fx=0.05, fy=0.95, log='xy')

AddCorrelationToPlot(log(halflife_use[, 1]), log(halflife_use[, 2]), xy[1], xy[2])

xy <- GetPlotFractionalCoords(fx=0.05, fy=0.90, log='xy')
AddCorrelationToPlot(log(halflife_use[, 1]), log(halflife_use[, 2]), 
                     xy[1], xy[2], method="spearman")

dev.copy2pdf(file=file.path("/lab/bartel1_ata/mcgeary/Presentations/seminars/",
                            "180604 Group Meeting/Correlation_with_FC.pdf", fsep=""))



# sites_8me <- which(rowSums(repression_matrix[, -ncol(repression_matrix)]) == 1 & repression_matrix$X8mer == 1 & !is.na(halflife_use))
# sites_7m8 <- which(rowSums(repression_matrix[, -ncol(repression_matrix)]) == 1 & repression_matrix$X7mer.m8 == 1 & !is.na(halflife_use))
# sites_7A1 <- which(rowSums(repression_matrix[, -ncol(repression_matrix)]) == 1 & repression_matrix$X7mer.A1 == 1 & !is.na(halflife_use))
# sites_6me <- which(rowSums(repression_matrix[, -ncol(repression_matrix)]) == 1 & repression_matrix$X6mer == 1 & !is.na(halflife_use))
# sites_6m8 <- which(rowSums(repression_matrix[, -ncol(repression_matrix)]) == 1 & repression_matrix$X6mer.m8 == 1 & !is.na(halflife_use))
# sites_6A1 <- which(rowSums(repression_matrix[, -ncol(repression_matrix)]) == 1 & repression_matrix$X6mer.A1 == 1 & !is.na(halflife_use))

sites_8me <- which(rowSums(repression_matrix[, -ncol(repression_matrix)]) == 1 & repression_matrix$X8mer == 1)
sites_7m8 <- which(rowSums(repression_matrix[, -ncol(repression_matrix)]) == 1 & repression_matrix$X7mer.m8 == 1)
sites_7A1 <- which(rowSums(repression_matrix[, -ncol(repression_matrix)]) == 1 & repression_matrix$X7mer.A1 == 1)
sites_6me <- which(rowSums(repression_matrix[, -ncol(repression_matrix)]) == 1 & repression_matrix$X6mer == 1)
sites_6m8 <- which(rowSums(repression_matrix[, -ncol(repression_matrix)]) == 1 & repression_matrix$X6mer.m8 == 1)
sites_6A1 <- which(rowSums(repression_matrix[, -ncol(repression_matrix)]) == 1 & repression_matrix$X6mer.A1 == 1)


# print(repression_matrix[sites_8me, ])
# print(halflife_use[sites_8me, ])



# break
colors <- c(rep(kSiteColors["8mer"], length(sites_8me)),
            rep(kSiteColors["7mer-m8"], length(sites_7m8)),
            rep(kSiteColors["7mer-A1"], length(sites_7A1)),
            rep(kSiteColors["6mer"], length(sites_6me)),
            rep(kSiteColors["6mer-m8"], length(sites_6m8)),
            rep(kSiteColors["6mer-A1"], length(sites_6A1)))


sites_all <- c(sites_8me,
			   sites_7m8,
			   sites_7A1,
			   sites_6me,
			   sites_6m8,
			   sites_6A1)

averages <- c(mean(repression_matrix$fc[sites_8me]),
			  mean(repression_matrix$fc[sites_7m8]),
			  mean(repression_matrix$fc[sites_7A1]),
			  mean(repression_matrix$fc[sites_6me]),
			  mean(repression_matrix$fc[sites_6m8]),
			  mean(repression_matrix$fc[sites_6A1]))

residuals <- c(repression_matrix$fc[sites_8me]-mean(repression_matrix$fc[sites_8me]),
			  repression_matrix$fc[sites_7m8]-mean(repression_matrix$fc[sites_7m8]),
			  repression_matrix$fc[sites_7A1]-mean(repression_matrix$fc[sites_7A1]),
			  repression_matrix$fc[sites_6me]-mean(repression_matrix$fc[sites_6me]),
			  repression_matrix$fc[sites_6m8]-mean(repression_matrix$fc[sites_6m8]),
			  repression_matrix$fc[sites_6A1]-mean(repression_matrix$fc[sites_6A1]))





print(averages)

fcs <- 2^averages

bg_half_life_cb <- exp(mean(log(halflife_use[, 1]/log(2)), na.rm=TRUE))
bg_half_life_tani <- exp(mean(log(halflife_use[, 2]/log(2)), na.rm=TRUE))

kds <- EquilPars("miR-1", sitelist="canonical", buffer=TRUE, combined=FALSE)

kd_names <- rownames(kds)

kds <- kds[1:6, 2]
names(kds) <- kd_names[1:6]

fcs <- GetRepressionLinearModel("miR-1", sitelist="canonical")[,1]

Model <- function(a, b, bg) {
	out <- log(1/(1/bg + a/(a + kds)*b)/(1/(1/bg)), 2)
	# print(out)
	# plot(fcs, out, col=kSiteColors[names(kds)], xlim=c(-1, 0.2), ylim=c(-1, 0.2))
	out
}

ResidualFunction <- function(pars, bg) {
	pars <- 10^pars
	y <- Model(pars[1], pars[2], bg)
	return(sum((y - fcs)^2))
}

matrix_output_cb <- matrix(NaN, nrow=100, ncol=3)
matrix_output_tani <- matrix(NaN, nrow=100, ncol=3)

matrix_input <- matrix(rnorm(100*4, mean=0, sd=2), nrow=100, ncol=4)

for (i in seq(100)) {
	solution_cb <- optim(matrix_input[i, 1:2], ResidualFunction,
	                  bg=bg_half_life_cb)
	solution_tani <- optim(matrix_input[i, 3:4], ResidualFunction,
	                    bg=bg_half_life_tani)

	matrix_output_cb[i, 1:2] <- solution_cb$par
	matrix_output_cb[i, 3] <- solution_cb$value
	matrix_output_tani[i, 1:2] <- solution_tani$par
	matrix_output_tani[i, 3] <- solution_tani$value

}


matrix_output_cb <- matrix_output_cb[which(matrix_output_cb[, 3] < 0.01), ]
matrix_output_tani <- matrix_output_tani[which(matrix_output_tani[, 3] < 0.01), ]



# print(pars_all$par)

pars_cb <- 10^colMeans(matrix_output_cb[, c(1, 2)])
pars_tani <- 10^colMeans(matrix_output_tani[, c(1, 2)])


a_cb <- pars_cb[1]
b_cb <- pars_cb[2]

a_tani <- pars_tani[1]
b_tani <- pars_tani[2]

occs_cb <- a_cb / (a_cb + kds)
occs_tani <- a_tani / (a_tani + kds)

kds_all <- rep.int(kds, c(length(sites_8me), length(sites_7m8),
                           length(sites_7A1), length(sites_6me),
                           length(sites_6m8), length(sites_6A1)))

inds_use <- 1:length(sites_all)

model_all_cb <- log(1/(1/halflife_use[, 1][sites_all][inds_use] + a_cb/(a_cb + kds_all[inds_use])*b_cb)/(1/(1/halflife_use[, 1][sites_all][inds_use])), 2)
model_all_tani <- log(1/(1/halflife_use[, 2][sites_all][inds_use] + a_tani/(a_tani + kds_all[inds_use])*b_tani)/(1/(1/halflife_use[, 2][sites_all][inds_use])), 2)



model_averages_cb <- log(1/(1/bg_half_life_cb + a_cb/(a_cb + kds)*b_cb)/(1/(1/bg_half_life_cb)), 2)
model_averages_tani <- log(1/(1/bg_half_life_tani + a_tani/(a_tani + kds)*b_tani)/(1/(1/bg_half_life_tani)), 2)


model_residual_all_cb <- model_all_cb - rep.int(fcs, c(length(sites_8me), length(sites_7m8),
                           length(sites_7A1), length(sites_6me),
                           length(sites_6m8), length(sites_6A1)))
model_residual_all_tani <- model_all_tani - rep.int(fcs, c(length(sites_8me), length(sites_7m8),
                           length(sites_7A1), length(sites_6me),
                           length(sites_6m8), length(sites_6A1)))

###### COFFEE #########################################

dev.new(xpos=10, ypos=20, height=4, width=4)
par(kPlotParameters)
xmin <- 0
xmax <- 1
ymin <- -1
ymax <- 0.2
BlankPlot()
AddLinearAxis(1, tick.space=0.1, label.space=0.2, label="Predicted repression (log2 fold change)")
AddLinearAxis(2, tick.space=0.1, label.space=0.2, label="Observed repression (log2 fold change)")
points(-model_averages_cb, fcs, col=kSiteColors[names(fcs)])
abline(0, -1, lty=2)
xy <- GetPlotFractionalCoords(fx=0.6, fy=0.95)
legend(xy[1], xy[2], legend=kSeedSites, col=kSiteColors[kSeedSites], pch=19, bty="n")

dev.copy2pdf(file=file.path("/lab/bartel1_ata/mcgeary/Presentations/seminars/",
                            "180604 Group Meeting/CB_1.pdf", fsep=""))


dev.new(xpos=400, ypos=20, height=4, width=4)
par(kPlotParameters)
xmin <- -3
xmax <- 1
ymin <- 0.1
ymax <- 100
BlankPlot(log='y', inv='x')
AddLinearAxis(1, tick.space=0.1, label.space=0.5, label="Predicted repression (log2 fold change)")
AddLogAxis(2, label="miRNA half-life (h)")
points(model_all_cb, halflife_use[, 1][sites_all][inds_use], col=colors, pch=19)
segments(averages, ymin, averages, ymax, col=kSiteColors[c("8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-m8", "6mer-A1")], lty=2)
xy <- GetPlotFractionalCoords(fx=0.7, fy=0.4, log='y', inv='x')
legend(xy[1], xy[2], legend=kSeedSites, col=kSiteColors[kSeedSites], pch=19, bty="n")

dev.copy2pdf(file=file.path("/lab/bartel1_ata/mcgeary/Presentations/seminars/",
                            "180604 Group Meeting/CB_2.pdf", fsep=""))

dev.new(xpos=790, ypos=20, height=4, width=4)
df_new <- data.frame(res=residuals, halflife=log10(halflife_use[sites_all, 1]))
df_model <- data.frame(res=model_residual_all_cb, halflife=log10(halflife_use[sites_all, 1]))
par(kPlotParameters)
xmin <- 0.1
xmax <- 100
ymin <- -1.5
ymax <- 1
BlankPlot(log='x', inv='y')
AddLogAxis(1, label="miRNA half-life (h)")
AddLinearAxis(2, tick.space=0.1, label.space=0.5, label="Residual repression (log2 fold change)")
abline(0, 0, lty=2, col="gray")
points(halflife_use[, 1][sites_all], residuals, col=colors, pch=1)
lowess_fit_8me <- lowess(res ~ halflife, data=df_new[1:length(sites_8me), ])
lowess_fit_7m8 <- lowess(res ~ halflife, data=df_new[(length(sites_8me)+1):(length(sites_8me) + length(sites_7m8)), ])
lowess_fit_7A1 <- lowess(res ~ halflife, data=df_new[(length(sites_8me) + length(sites_7m8) + 1):(length(sites_8me) + length(sites_7m8) + length(sites_7A1)), ])
lowess_fit_6me <- lowess(res ~ halflife, data=df_new[(length(sites_8me) + length(sites_7m8) + length(sites_7A1) + 1):(length(sites_8me) + length(sites_7m8) + length(sites_7A1) + length(sites_6me)), ])
lowess_fit_6m8 <- lowess(res ~ halflife, data=df_new[(length(sites_8me) + length(sites_7m8) + length(sites_7A1) + length(sites_6me) + 1):(length(sites_8me) + length(sites_7m8) + length(sites_7A1) + length(sites_6me) + length(sites_6m8)), ])
lowess_fit_6A1 <- lowess(res ~ halflife, data=df_new[(length(sites_8me) + length(sites_7m8) + length(sites_7A1) + length(sites_6me) + length(sites_6m8) + 1):length(sites_all), ])
lines(10^lowess_fit_8me$x, lowess_fit_8me$y, col=kSiteColors["8mer"], lwd=2)
lines(10^lowess_fit_7m8$x, lowess_fit_7m8$y, col=kSiteColors["7mer-m8"], lwd=2)
lines(10^lowess_fit_7A1$x, lowess_fit_7A1$y, col=kSiteColors["7mer-A1"], lwd=2)
lines(10^lowess_fit_6me$x, lowess_fit_6me$y, col=kSiteColors["6mer"], lwd=2)
lines(10^lowess_fit_6m8$x, lowess_fit_6m8$y, col=kSiteColors["6mer-m8"], lwd=2)
lines(10^lowess_fit_6A1$x, lowess_fit_6A1$y, col=kSiteColors["6mer-A1"], lwd=2)
xy <- GetPlotFractionalCoords(fx=0.7, fy=0.4, log='x', inv='y')
legend(xy[1], xy[2], legend=kSeedSites, col=kSiteColors[kSeedSites], pch=19, bty="n")

dev.copy2pdf(file=file.path("/lab/bartel1_ata/mcgeary/Presentations/seminars/",
                            "180604 Group Meeting/CB_3.pdf", fsep=""))


dev.new(xpos=1180, ypos=20, height=4, width=4)
par(kPlotParameters)
xmin <- -3
xmax <- 1
ymin <- 0.1
ymax <- 100
BlankPlot(log='y', inv='x')
AddLinearAxis(1, tick.space=0.1, label.space=0.5, label="Predicted repression (log2 fold change)")
AddLogAxis(2, label="miRNA half-life (h)")
points(repression_matrix$fc[sites_all][inds_use], halflife_use[, 1][sites_all][inds_use],  col=colors, pch=19)
segments(averages, ymin, averages, ymax, col=kSiteColors[c("8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-m8", "6mer-A1")], lty=2)
# segments(fcs, 2, fcs, 25, col=kSiteColors[c("8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-m8", "6mer-A1")], lty=2)
xy <- GetPlotFractionalCoords(fx=0.7, fy=0.4, log='y', inv='x')
legend(xy[1], xy[2], legend=kSeedSites, col=kSiteColors[kSeedSites], pch=19, bty="n")

dev.copy2pdf(file=file.path("/lab/bartel1_ata/mcgeary/Presentations/seminars/",
                            "180604 Group Meeting/CB_4.pdf", fsep=""))


dev.new(xpos=1570, ypos=20, height=4, width=4)
par(kPlotParameters)
xmin <- 0.1
xmax <- 100
ymin <- -1.5
ymax <- 1
BlankPlot(log='x', inv='y')
AddLogAxis(1, label="miRNA half-life (h)")
AddLinearAxis(2, tick.space=0.1, label.space=0.5, label="Residual repression (log2 fold change)")
abline(0, 0, lty=2, col="gray")
points(halflife_use[sites_all, 1], model_residual_all_cb, col=colors, pch=1)
lowess_fit2_8me <- lowess(res ~ halflife, data=df_model[1:length(sites_8me), ])
lowess_fit2_7m8 <- lowess(res ~ halflife, data=df_model[(length(sites_8me)+1):(length(sites_8me) + length(sites_7m8)), ])
lowess_fit2_7A1 <- lowess(res ~ halflife, data=df_model[(length(sites_8me) + length(sites_7m8) + 1):(length(sites_8me) + length(sites_7m8) + length(sites_7A1)), ])
lowess_fit2_6me <- lowess(res ~ halflife, data=df_model[(length(sites_8me) + length(sites_7m8) + length(sites_7A1) + 1):(length(sites_8me) + length(sites_7m8) + length(sites_7A1) + length(sites_6me)), ])
lowess_fit2_6m8 <- lowess(res ~ halflife, data=df_model[(length(sites_8me) + length(sites_7m8) + length(sites_7A1) + length(sites_6me) + 1):(length(sites_8me) + length(sites_7m8) + length(sites_7A1) + length(sites_6me) + length(sites_6m8)), ])
lowess_fit2_6A1 <- lowess(res ~ halflife, data=df_model[(length(sites_8me) + length(sites_7m8) + length(sites_7A1) + length(sites_6me) + length(sites_6m8) + 1):length(sites_all), ])
lines(10^lowess_fit2_8me$x, lowess_fit2_8me$y, col=kSiteColors["8mer"], lwd=2)
lines(10^lowess_fit2_7m8$x, lowess_fit2_7m8$y, col=kSiteColors["7mer-m8"], lwd=2)
lines(10^lowess_fit2_7A1$x, lowess_fit2_7A1$y, col=kSiteColors["7mer-A1"], lwd=2)
lines(10^lowess_fit2_6me$x, lowess_fit2_6me$y, col=kSiteColors["6mer"], lwd=2)
lines(10^lowess_fit2_6m8$x, lowess_fit2_6m8$y, col=kSiteColors["6mer-m8"], lwd=2)
lines(10^lowess_fit2_6A1$x, lowess_fit2_6A1$y, col=kSiteColors["6mer-A1"], lwd=2)
xy <- GetPlotFractionalCoords(fx=0.7, fy=0.4, log='x', inv='y')
legend(xy[1], xy[2], legend=kSeedSites, col=kSiteColors[kSeedSites], pch=19, bty="n")

dev.copy2pdf(file=file.path("/lab/bartel1_ata/mcgeary/Presentations/seminars/",
                            "180604 Group Meeting/CB_5.pdf", fsep=""))



########### TANI ET AL ################################################

dev.new(xpos=10, ypos=420, height=4, width=4)
par(kPlotParameters)
xmin <- 0
xmax <- 1
ymin <- -1
ymax <- 0.2
BlankPlot()
AddLinearAxis(1, tick.space=0.1, label.space=0.2, label="Predicted repression (log2 fold change)")
AddLinearAxis(2, tick.space=0.1, label.space=0.2, label="Observed repression (log2 fold change)")
points(-model_averages_2, fcs, col=kSiteColors[names(fcs)])
abline(0, -1, lty=2)
xy <- GetPlotFractionalCoords(fx=0.6, fy=0.95)
legend(xy[1], xy[2], legend=kSeedSites, col=kSiteColors[kSeedSites], pch=19, bty="n")

dev.copy2pdf(file=file.path("/lab/bartel1_ata/mcgeary/Presentations/seminars/",
                            "180604 Group Meeting/TANI_1.pdf", fsep=""))


dev.new(xpos=400, ypos=420, height=4, width=4)
par(kPlotParameters)
xmin <- -3
xmax <- 1
ymin <- 0.1
ymax <- 100
BlankPlot(log='y', inv='x')
AddLinearAxis(1, tick.space=0.1, label.space=0.5, label="Predicted repression (log2 fold change)")
AddLogAxis(2, label="miRNA half-life (h)")
points(model_all_tani, halflife_use[, 2][sites_all][inds_use], col=colors, pch=19)
segments(averages, ymin, averages, ymax, col=kSiteColors[c("8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-m8", "6mer-A1")], lty=2)
xy <- GetPlotFractionalCoords(fx=0.7, fy=0.4, log='y', inv='x')
legend(xy[1], xy[2], legend=kSeedSites, col=kSiteColors[kSeedSites], pch=19, bty="n")

dev.copy2pdf(file=file.path("/lab/bartel1_ata/mcgeary/Presentations/seminars/",
                            "180604 Group Meeting/TANI_2.pdf", fsep=""))

dev.new(xpos=790, ypos=420, height=4, width=4)
df_new <- data.frame(res=residuals, halflife=log10(halflife_use[sites_all, 2]))
df_model <- data.frame(res=model_residual_all_tani, halflife=log10(halflife_use[sites_all, 2]))


par(kPlotParameters)
xmin <- 0.1
xmax <- 100
ymin <- -1.5
ymax <- 1
BlankPlot(log='x', inv='y')
AddLogAxis(1, label="miRNA half-life (h)")
AddLinearAxis(2, tick.space=0.1, label.space=0.5, label="Residual repression (log2 fold change)")
abline(0, 0, lty=2, col="gray")
points(halflife_use[sites_all, 2], residuals, col=colors, pch=1)
lowess_fit_8me <- lowess(res ~ halflife, data=df_new[1:length(sites_8me), ])
lowess_fit_7m8 <- lowess(res ~ halflife, data=df_new[(length(sites_8me)+1):(length(sites_8me) + length(sites_7m8)), ])
lowess_fit_7A1 <- lowess(res ~ halflife, data=df_new[(length(sites_8me) + length(sites_7m8) + 1):(length(sites_8me) + length(sites_7m8) + length(sites_7A1)), ])
lowess_fit_6me <- lowess(res ~ halflife, data=df_new[(length(sites_8me) + length(sites_7m8) + length(sites_7A1) + 1):(length(sites_8me) + length(sites_7m8) + length(sites_7A1) + length(sites_6me)), ])
lowess_fit_6m8 <- lowess(res ~ halflife, data=df_new[(length(sites_8me) + length(sites_7m8) + length(sites_7A1) + length(sites_6me) + 1):(length(sites_8me) + length(sites_7m8) + length(sites_7A1) + length(sites_6me) + length(sites_6m8)), ])
lowess_fit_6A1 <- lowess(res ~ halflife, data=df_new[(length(sites_8me) + length(sites_7m8) + length(sites_7A1) + length(sites_6me) + length(sites_6m8) + 1):length(sites_all), ])
lines(10^lowess_fit_8me$x, lowess_fit_8me$y, col=kSiteColors["8mer"], lwd=2)
lines(10^lowess_fit_7m8$x, lowess_fit_7m8$y, col=kSiteColors["7mer-m8"], lwd=2)
lines(10^lowess_fit_7A1$x, lowess_fit_7A1$y, col=kSiteColors["7mer-A1"], lwd=2)
lines(10^lowess_fit_6me$x, lowess_fit_6me$y, col=kSiteColors["6mer"], lwd=2)
lines(10^lowess_fit_6m8$x, lowess_fit_6m8$y, col=kSiteColors["6mer-m8"], lwd=2)
lines(10^lowess_fit_6A1$x, lowess_fit_6A1$y, col=kSiteColors["6mer-A1"], lwd=2)
xy <- GetPlotFractionalCoords(fx=0.7, fy=0.4, log='x', inv='y')
legend(xy[1], xy[2], legend=kSeedSites, col=kSiteColors[kSeedSites], pch=19, bty="n")

dev.copy2pdf(file=file.path("/lab/bartel1_ata/mcgeary/Presentations/seminars/",
                            "180604 Group Meeting/TANI_3.pdf", fsep=""))


dev.new(xpos=1180, ypos=420, height=4, width=4)
par(kPlotParameters)
xmin <- -3
xmax <- 1
ymin <- 0.1
ymax <- 100
BlankPlot(log='y', inv='x')
AddLinearAxis(1, tick.space=0.1, label.space=0.5, label="Predicted repression (log2 fold change)")
AddLogAxis(2, label="miRNA half-life (h)")
points(repression_matrix$fc[sites_all][inds_use], halflife_use[, 2][sites_all][inds_use],  col=colors, pch=19)
segments(averages, ymin, averages, ymax, col=kSiteColors[c("8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-m8", "6mer-A1")], lty=2)
# segments(fcs, 2, fcs, 25, col=kSiteColors[c("8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-m8", "6mer-A1")], lty=2)
xy <- GetPlotFractionalCoords(fx=0.7, fy=0.4, log='y', inv='x')
legend(xy[1], xy[2], legend=kSeedSites, col=kSiteColors[kSeedSites], pch=19, bty="n")

dev.copy2pdf(file=file.path("/lab/bartel1_ata/mcgeary/Presentations/seminars/",
                            "180604 Group Meeting/TANI_4.pdf", fsep=""))


dev.new(xpos=1570, ypos=420, height=4, width=4)
par(kPlotParameters)
xmin <- 0.1
xmax <- 100
ymin <- -1.5
ymax <- 1
BlankPlot(log='x', inv='y')
AddLogAxis(1, label="miRNA half-life (h)")
AddLinearAxis(2, tick.space=0.1, label.space=0.5, label="Residual repression (log2 fold change)")
abline(0, 0, lty=2, col="gray")
points(halflife_use[, 2][sites_all], model_residual_all_tani, col=colors, pch=1)
lowess_fit2_8me <- lowess(res ~ halflife, data=df_model[1:length(sites_8me), ])
lowess_fit2_7m8 <- lowess(res ~ halflife, data=df_model[(length(sites_8me)+1):(length(sites_8me) + length(sites_7m8)), ])
lowess_fit2_7A1 <- lowess(res ~ halflife, data=df_model[(length(sites_8me) + length(sites_7m8) + 1):(length(sites_8me) + length(sites_7m8) + length(sites_7A1)), ])
lowess_fit2_6me <- lowess(res ~ halflife, data=df_model[(length(sites_8me) + length(sites_7m8) + length(sites_7A1) + 1):(length(sites_8me) + length(sites_7m8) + length(sites_7A1) + length(sites_6me)), ])
lowess_fit2_6m8 <- lowess(res ~ halflife, data=df_model[(length(sites_8me) + length(sites_7m8) + length(sites_7A1) + length(sites_6me) + 1):(length(sites_8me) + length(sites_7m8) + length(sites_7A1) + length(sites_6me) + length(sites_6m8)), ])
lowess_fit2_6A1 <- lowess(res ~ halflife, data=df_model[(length(sites_8me) + length(sites_7m8) + length(sites_7A1) + length(sites_6me) + length(sites_6m8) + 1):length(sites_all), ])
lines(10^lowess_fit2_8me$x, lowess_fit2_8me$y, col=kSiteColors["8mer"], lwd=2)
lines(10^lowess_fit2_7m8$x, lowess_fit2_7m8$y, col=kSiteColors["7mer-m8"], lwd=2)
lines(10^lowess_fit2_7A1$x, lowess_fit2_7A1$y, col=kSiteColors["7mer-A1"], lwd=2)
lines(10^lowess_fit2_6me$x, lowess_fit2_6me$y, col=kSiteColors["6mer"], lwd=2)
lines(10^lowess_fit2_6m8$x, lowess_fit2_6m8$y, col=kSiteColors["6mer-m8"], lwd=2)
lines(10^lowess_fit2_6A1$x, lowess_fit2_6A1$y, col=kSiteColors["6mer-A1"], lwd=2)
xy <- GetPlotFractionalCoords(fx=0.7, fy=0.4, log='x', inv='y')
legend(xy[1], xy[2], legend=kSeedSites, col=kSiteColors[kSeedSites], pch=19, bty="n")

dev.copy2pdf(file=file.path("/lab/bartel1_ata/mcgeary/Presentations/seminars/",
                            "180604 Group Meeting/TANI_5.pdf", fsep=""))







break
sXc <- SitesXCounts("miR-1", exp="kinetics")
print(sXc)

dim_p <- nrow(sXc)/2
pos_I   <- GetPositionalSites("miR-1", exp="kinetics", condition="I")/sum(sXc[1:dim_p, 1])
pos_0   <- GetPositionalSites("miR-1", exp="kinetics", condition="0")/sum(sXc[1:dim_p, 2])
pos_2.1 <- GetPositionalSites("miR-1", exp="kinetics", condition="2,1")/sum(sXc[1:dim_p, 3])
pos_2.2 <- GetPositionalSites("miR-1", exp="kinetics", condition="2,2")/sum(sXc[1:dim_p, 4])
pos_5.1 <- GetPositionalSites("miR-1", exp="kinetics", condition="5,1")/sum(sXc[1:dim_p, 5])
pos_5.2 <- GetPositionalSites("miR-1", exp="kinetics", condition="5,2")/sum(sXc[1:dim_p, 6])
pos_8.1 <- GetPositionalSites("miR-1", exp="kinetics", condition="8,1")/sum(sXc[1:dim_p, 7])
pos_8.2 <- GetPositionalSites("miR-1", exp="kinetics", condition="8,2")/sum(sXc[1:dim_p, 8])
pos_15.1 <- GetPositionalSites("miR-1", exp="kinetics", condition="15,1")/sum(sXc[1:dim_p, 9])
pos_15.2 <- GetPositionalSites("miR-1", exp="kinetics", condition="15,2")/sum(sXc[1:dim_p, 10])
pos_30.1 <- GetPositionalSites("miR-1", exp="kinetics", condition="30,1")/sum(sXc[1:dim_p, 11])
pos_30.2 <- GetPositionalSites("miR-1", exp="kinetics", condition="30,2")/sum(sXc[1:dim_p, 12])
pos_90 <- GetPositionalSites("miR-1", exp="kinetics", condition="90")/sum(sXc[1:dim_p, 13])
pos_300 <- GetPositionalSites("miR-1", exp="kinetics", condition="300")/sum(sXc[1:dim_p, 14])
pos_900 <- GetPositionalSites("miR-1", exp="kinetics", condition="900")/sum(sXc[1:dim_p, 15])
pos_2700 <- GetPositionalSites("miR-1", exp="kinetics", condition="2700")/sum(sXc[1:dim_p, 16])
pos_7200 <- GetPositionalSites("miR-1", exp="kinetics", condition="7200")/sum(sXc[1:dim_p, 17])
pos_Equil <- GetPositionalSites("miR-1", exp="kinetics", condition="Equil")/sum(sXc[1:dim_p, 18])

R_0 <-   pos_0  /pos_I
R_2.1 <- pos_2.1/pos_I
R_2.2 <- pos_2.2/pos_I


x <- 1:ncol(pos_I)
y_p <- R_0["8mer_p", ]
y_c <- R_0["8mer_c", ]

plot(1:ncol(pos_I), y_p, type="o", log="y", ylim=c(0.01, 100))
points(1:ncol(pos_I), pos_2.1["8mer_p", ]/(pos_I["8mer_p", ]), type="o")
points(1:ncol(pos_I), pos_2.2["8mer_p", ]/(pos_I["8mer_p", ]), type="o")
points(1:ncol(pos_I), pos_5.1["8mer_p", ]/(pos_I["8mer_p", ]), type="o")
points(1:ncol(pos_I), pos_5.2["8mer_p", ]/(pos_I["8mer_p", ]), type="o")
points(1:ncol(pos_I), pos_8.1["8mer_p", ]/(pos_I["8mer_p", ]), type="o")
points(1:ncol(pos_I), pos_8.2["8mer_p", ]/(pos_I["8mer_p", ]), type="o")
points(1:ncol(pos_I), pos_15.1["8mer_p",]/(pos_I["8mer_p", ]), type="o")
points(1:ncol(pos_I), pos_15.2["8mer_p",]/(pos_I["8mer_p", ]), type="o")
points(1:ncol(pos_I), pos_30.1["8mer_p",]/(pos_I["8mer_p", ]), type="o")
points(1:ncol(pos_I), pos_30.2["8mer_p",]/(pos_I["8mer_p", ]), type="o")
points(1:ncol(pos_I),  pos_90["8mer_p", ]/(pos_I["8mer_p", ]), type="o")
points(1:ncol(pos_I), pos_300["8mer_p", ]/(pos_I["8mer_p", ]), type="o")
points(1:ncol(pos_I), pos_900["8mer_p", ]/(pos_I["8mer_p", ]), type="o")
points(1:ncol(pos_I), pos_2700["8mer_p",]/(pos_I["8mer_p", ]), type="o")
points(1:ncol(pos_I), pos_7200["8mer_p",]/(pos_I["8mer_p", ]), type="o")
points(1:ncol(pos_I),pos_Equil["8mer_p",]/(pos_I["8mer_p", ]), type="o")

print(sXck)