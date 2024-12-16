################################################################################
#GenerateSiteTypeKds.py
################################################################################

# # Initial parameters and constants.
args = commandArgs(trailingOnly=TRUE)
mirna = args[1]
method = args[2]
source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")

site_cols <- read.table(
	"/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/site_colors.txt",
	row.names=1, header=FALSE, sep="\t", stringsAsFactors=FALSE)
print(site_cols)
stockago = read.table("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/SolveForKds/k_c_stockago.txt",
	row.names=1,header=FALSE,sep="\t")
if (length(args)==3) {
	k.c.stockago <- as.numeric(args[3])
} else {
	k.c.stockago = stockago[mirna,1]
}
k.c.lib = 100
print(k.c.stockago)
## Functions ###################################################################
## I/O functions
GetSitesXcounts <- function(mirna) {
	sites_file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
                       mirna,"/equilibrium/site_count_tables/all_sites.txt")
	sitesXcounts <- read.table(sites_file_name)
	return(sitesXcounts)
}

GetSiteKds <- function(mirna,k.c.stockago) {
	sites_file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
	     						 "/equilibrium/kds/final_", k.c.stockago, ".txt")
	site.kds <- read.table(sites_file_name,row.names=1,header=FALSE)
	return(site.kds)
}

# Modeling Functions

GetOccupancy <- function(c.freeago, kds) {
	return(c.freeago / (c.freeago + kds))
}
GetFreeResidual <- function(c.freeago, kds,c.sites,c.ago) {
	residual <- (c.ago - c.freeago - sum(GetOccupancy(c.freeago, kds)*c.sites))^2
	return(residual)
}
GetFreeAgo <- function(kds, c.sites, c.ago) {
	c.free <- optimize(GetFreeResidual,
	         					 c(0,c.ago),
	         					 kds=kds,
	       					   c.sites=c.sites,
	         					 c.ago=c.ago)$minimum
	return(c.free)
}
# Output Functions:
CheckMaxDifference <- function(out) {
	row.last <- dim(out)[1]
	row.secondtolast <- row.last-1
	diffs <- sapply(1:dim(out)[2],function(x) {
		abs((out[row.last,x]-out[row.secondtolast,x])/out[row.secondtolast,x])
	})
	return(max(diffs))
}




# MAIN #########################################################################

# Get data table:

sitesXcounts <- GetSitesXcounts(mirna)
colnames(sitesXcounts)[2:dim(sitesXcounts)[2]]=sapply(
	colnames(sitesXcounts)[2:dim(sitesXcounts)[2]], function(x) {
		return(unlist(strsplit(x,split="A",fixed=TRUE))[2])
	}
)
# Assign parameters for Kds and background, subracting 1 from rows due to "None"
# Kd assigned to 1 nM, and the fact that the I and A0 columns are not modeled
# in the experiment, so there's no background term assigned ot them.

# Assign the number of kd parameters and background parameters:
num.kds <- dim(sitesXcounts)[1]-1
num.bgs <- dim(sitesXcounts)[2]-2

# Assign the total site concentration in each experiment, and initialize the
# data column to be used.
# k.c.lib should be 100, for 100 nM.
c.sites <- sitesXcounts["I"]/sum(sitesXcounts["I"])*k.c.lib

# Remove the I and A0 columns from the data to be fit to the model. 
data <- sitesXcounts[, 2:6]
print(data)

# Initialize starting Kds, which are set to 1/the enrichment of each site type
# in the A12.6 experiment.
r_numerator <- data[, 2] / sum(data[, 2])
r_denomenator <- sitesXcounts["I"] / sum(sitesXcounts["I"])
kds.init <- (r_numerator / r_denomenator)[1 : num.kds, 1]
pars.init <- c(-log10(kds.init) - 1, -1, -1, -1, -1, -1)

site_colors = sample(colors(),dim(data)[1])
tick <- 0
# Define function of just kds and pars.

GetModelFrequencies <- function(pars) {
	# Split up the parameters into the kd and background parameters.
	kds <- 10 ^ c(pars[1 : num.kds], 0)

	bgs <- 10 ^ pars[(num.kds + 1):(num.kds + num.bgs)]



	# Solve for the free Ago concentration in each experiment.
	c.freeagos <- sapply(colnames(data), function(x) {
		c.ago <- as.numeric(x) / 100 * k.c.stockago
			return(GetFreeAgo(kds, c.sites, c.ago))
		}
	)
	
	# Initialize a matrix with the same total concentration of each site type
	# for each experiment.
	c.totals <- as.matrix(
		sapply(colnames(data), function(x) {
					   return(c.sites[[1]])
					 }
					 ))

	rownames(c.totals) <- rownames(data)
	
	# Use the free Ago concentrations to get the amount of each complex bound
	# to Ago.
	c.bounds <- as.matrix(
		sapply(c.freeagos, function(x) {
					   return((GetOccupancy(x, kds) * c.sites)[[1]])
					 }
					 ))
	rownames(c.bounds) <- rownames(data)

	# Get the amount of background binding by subtracting the bound from the
	# total sites in each experiment, normalizing. Must transpose to multiply
	# each column.
	c.frees <- c.totals - c.bounds
	c.bgs <- t(t(c.frees) * bgs / colSums(c.frees))
	c.all <- c.bounds + c.bgs
	c.final <- data.frame(t(t(c.all) / colSums(c.all)))
	colnames(c.final) <- colnames(data)

	prob <- sapply(seq(1, 5), function(x) {
				   			 dmultinom(data[, x],
									 				 prob=as.numeric(c.final[, x]), log=TRUE)
								 }
								 )


	if (tick%%100 == 0) {
		print(-sum(prob))
		setEPS(width=10)
		postscript(file=paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/figures/kds/",
	mirna,"/temp_",method,"_", k.c.stockago, ".eps"))

		x <- c(40,12.65,4,1.265,0.4)*k.c.stockago/100*1000
		y <- c(1,1,1,1,1)
		sites.norm <- c.sites / sum(c.sites)

		data.norm <- t(t(data)/colSums(data))

		data.R <- data.norm/(sites.norm[,1])

		model.R <- c.final/(sites.norm[,1])


		xmin <- floor(0.5*min(x))
		xmax <- ceiling(2*max(x))

		ymin <- 0.2
		ymax <- ceiling(max(data.R)*2)
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
		legend_names <- rownames(data)
		# centered_sites <- c("11mer-m3.13", "12mer-m3.14",
		# 								"11mer-m4.14", "12mer-m4.15")
		# centered_index <- sapply(centered_sites, function(site){
		# 										   which(legend_names == site)
		# 										 }
		# 										 )
		# legend_names[centered_index] <- "Centered"
		# legend_names <- unique(legend_names)
		legend(x=xmax,y=ymax,legend=legend_names, pch=19, col=site_cols[legend_names, ], cex=1.2, bty="n")
		for (name in rownames(data)) {
			points(x, data.R[name, ], col=site_cols[name, ], pch=19, lwd=3)
			lines(x, model.R[name, ], col=site_cols[name, ], lwd=2)
			
		}


		dev.off()

	}
	tick <<- tick + 1
	return(-sum(prob))
}
val <- GetModelFrequencies(pars.init)
print("NOW IT'S OPTIMIZING!!")
# Solve the first run of the function, and create the output matrix.
solution <- optim(pars.init, GetModelFrequencies,
									method="Nelder-Mead",
									control=c("maxit"<-20000))
out <- rbind(c(pars.init, val),
						 c(solution$par, solution$value))
colnames(out) <- c(rownames(data)[1:num.kds],colnames(data),"-logp")

# Get maximum difference in output.
converg <- CheckMaxDifference(out)
print(converg)

# Assign stopping criteria for continuing the optimizagin.
while (converg >= 0.000001) {
	solution <- optim(solution$par, GetModelFrequencies,
										method=method,
										control=c("maxit"<-50000))
	out <- rbind(out,c(solution$par,solution$value))
	converg <- CheckMaxDifference(out)
	print(converg)
	out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
	     						 "/equilibrium/kds/", method, "_", k.c.stockago, ".txt")
	write.table(file=out_file, out, sep="\t", quote=FALSE, row.names=FALSE,
							col.names=TRUE)

	# Assign output file For the final parameters and write to it.
	out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
		     						 "/equilibrium/kds/final_", method, "_", k.c.stockago, ".txt")
	out_final <- out[dim(out)[1],]
	names(out_final) <- colnames(out)

	write.table(file=out_file, out_final, sep="\t", quote=FALSE, row.names=TRUE,
							col.names=FALSE)

	MakeSiteIterationPlot(out,method)	
}

# Assign output file for the entire sequnce of the optimization and write
# to it.

