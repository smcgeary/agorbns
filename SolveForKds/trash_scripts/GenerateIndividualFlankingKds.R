################################################################################
#GenerateSiteTypeKds.py
################################################################################

# Initial parameters and constants.
args = commandArgs(trailingOnly=TRUE)
mirna = args[1]
site = args[2]
site_cols <- read.table(
	"/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/site_colors.txt",
	row.names=1, header=FALSE, sep="\t", stringsAsFactors=FALSE)

source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")
stockago = read.table("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/SolveForKds/k_c_stockago.txt",
	row.names=1,header=FALSE,sep="\t")

k.c.stockago = stockago[mirna,1]

k.c.lib = 100

## Functions ###################################################################
## I/O functions
GetSitesXcounts <- function(mirna) {
	sites_file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
                       mirna,"/equilibrium/site_count_tables/all_sites.txt")
	sitesXcounts <- read.table(sites_file_name)
	return(sitesXcounts)
}

GetSiteFlanksXcounts <- function(mirna,site) {
	sites_file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
                       mirna,"/equilibrium/site_count_tables/",
                       site,"_flanking.txt")
	siteflanksXcounts <- read.table(sites_file_name)
	return(siteflanksXcounts)
}


GetSiteKds <- function(mirna,k.c.stockago) {
	sites_file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
	     						 "/equilibrium/kds/final_Nelder-Mead_", k.c.stockago, ".txt")
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



tick <- 0

# MAIN #########################################################################

# Get siteXcounts, called sXc:
sXc <- GetSitesXcounts(mirna)
colnames(sXc)[2:dim(sXc)[2]]=sapply(
	colnames(sXc)[2:dim(sXc)[2]], function(x) {
		return(unlist(strsplit(x,split="A",fixed=TRUE))[2])
	}
)

# Get vector of single site counts s.c.
s.c <- as.numeric(sXc[site, ])

sfXc <- GetSiteFlanksXcounts(mirna,site)


sfXc <- round(t((t(sfXc) / colSums(sfXc)) * s.c))
colnames(sfXc) <- colnames(sXc)
# Assign parameters for Kds and background, subracting 1 from rows due to "None"
# Kd assigned to 1 nM, and the fact that the I and A0 columns are not modeled
# in the experiment, so there's no background term assigned ot them.
c.freeagos <- 0
pars <- GetSiteKds(mirna,k.c.stockago)


# Get the number of flanking sites
num.sf <- dim(sfXc)[1]
# Get the number of parameters (kds + bg + 1 log(prob))
num.pars <- dim(pars)[1]

bgs <- 10 ^ pars[(num.pars - 5): (num.pars - 1), ]

# Get the site kds
kds.s <- pars[1:(num.pars - 6), , drop=FALSE]
# Omit the site kd for which the flanking sites are being fit.
kds.s <- kds.s[rownames(kds.s) != site, ]

kds.s <- as.numeric(kds.s)

# Assign the total site concentration in each experiment, and initialize the
# data column to be used.
# k.c.lib should be 100, for 100 nM.
# print(kds.s)
sXc <- rbind(sfXc, sXc)
sXc <- sXc[rownames(sXc) != site, ]

c.sites <- sXc["I"]/sum(sXc["I"])*k.c.lib


# Remove the I and A0 columns from the data to be fit to the model. 
data <- sXc[, 2:6]

plot_range <- c(sample(1:num.sf,10),(num.sf+1):(dim(data)[1]))

# Initialize starting Kds, which are set to 1/the enrichment of each site type
# in the A12.6 experiment.
site_colors = sample(colors(),dim(data)[1])
pars.init <- as.numeric(rep(pars[site,][[1]],num.sf))
names(pars.init) <- rownames(sfXc)

# Define function of just kds and pars.
GetModelFrequencies <- function(par, index=FALSE, pars=FALSE) {

	# Split up the parameters into the kd and background parameters.
	if (index!=FALSE) {
		kds <- 10 ^ c(pars,kds.s, 0)
		kds[index] <- 10^par
	} else {
		kds <- 10 ^ c(par,kds.s, 0)		
	}
	names(kds) <- rownames(data)
	names(bgs) <- colnames(data)

	# Solve for the free Ago concentration in each experiment.
	c.freeagos <<- sapply(colnames(data), function(x) {
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
								 dmultinom(as.integer(data[, x]),
									 				 prob=as.numeric(c.final[, x]), log=TRUE)
								 }
								 )

	if (tick%%100 == 0) {

		setEPS(width=10)
		postscript(file=paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/figures/kds/",
	mirna,"/flanking_individual_fits/model_fits/",site,"_", k.c.stockago, ".eps"))
		par(par_plots)
		x <- c(40,12.65,4,1.265,0.4)*k.c.stockago/100*1000
		y <- c(1,1,1,1,1)
		sites.norm <- c.sites / sum(c.sites)

		data.norm <- t(t(data)/colSums(data))

		data.R <- data.norm/(sites.norm[,1])

		model.R <- c.final/(sites.norm[,1])


		xmin <- floor(0.5*min(x))
		xmax <- ceiling(2*max(x))

		ymin <- 0.2
		ymax <- ceiling(max(data.R[data.R!=Inf],na.rm=TRUE)*2)
		yextension <- (ymax/ymin)

		plot(x, y,log='xy', xlim=c(xmin, xmax*(xmax/xmin)^0.55), ylim=c(ymin, ymax), type="l",
			 col="white", axes=FALSE)				
		# Generate tickmarks for axis.
		xs <- c(xmin,sapply(seq(floor(log10(xmin)), ceiling(log10(xmax))), function(x) seq(10)*10^x),xmax)
		xs <- xs[xs >= xmin & xs <= xmax]
		xl <- 10^seq(ceiling(log10(min(xs))), floor(log10(max(xs))))
		axis(1, at=xl,
				 labels=sapply(xl, function(name) {
				   eval(substitute(expression(10^x), list(x=log10(name))))
				 }),
				 pos=ymin, lwd=0, padj=0.2)
		axis(1, at=xs, labels=FALSE,
				 pos=ymin, lwd=2)
		# Label the axis at each order of magnitude.
		ys <- c(ymin,sapply(seq(floor(log10(ymin)), ceiling(log10(ymax))), function(x) seq(10)*10^x),ymax)
		ys <- ys[ys >= ymin & ys <= ymax]
		yl <- 10^seq(ceiling(log10(min(ys))), floor(log10(max(ys))))
		axis(2, at=yl,
				 labels=sapply(yl, function(name) {
				   eval(substitute(expression(10^x), list(x=log10(name))))
				 }),
				 pos=xmin, las=2, lwd=0, hadj=0.8)
		axis(2, at=ys, labels=FALSE,
				 pos=xmin, lwd=2)

		title(main = mirna,line=-1, adj=0.1)
		title(main = site, line= -2.5, col.main=site_cols[site,], adj=0.1)

		title(xlab = "[AGO2-miRNA] (pM)", adj=0.3)
		title(ylab = "Enrichment")
		legend_names <- rownames(data)

		legend(x=xmax,y=ymax,legend=legend_names, pch=19, col=sapply(legend_names,GetColorFunction), cex=0.5, bty="n",ncol=4)
		for (name in rev(rownames(data))) {
			points(x, data.R[name, ], col=GetColorFunction(name), pch=19)
			lines(x, model.R[name, ], col=GetColorFunction(name))
			
		}
		dev.off()
	}
	tick <<- tick + 1
	return(-sum(prob))
}
pars.optim <- pars.init
val <- GetModelFrequencies(pars.init[1],1,pars.init)
solution <- optimize(GetModelFrequencies, c(-10,10), index= 1, pars=pars.optim)
pars.optim[1] <- solution$minimum
out <- rbind(c(pars.init, val),
							 c(pars.optim, solution$objective))	
colnames(out)[dim(out)[2]] <- rownames(pars)[dim(pars)[1]]
out.init <- out
converg <- 10


# Solve the first run of the function, and create the output matrix.
while (converg >= 0.0001) {
	for (i in 1:length(pars.init)) {
		solution <- optimize(GetModelFrequencies, c(-10,10), index= i, pars=pars.optim)
		pars.optim[i] <- solution$minimum

		out <- rbind(out,c(pars.optim,solution$objective))


		if (dim(out)[1]%%10 == 1) {
			inds <- seq(1,dim(out)[1],by=10)
			WriteIterationFile(out[inds,])
			WriteFinalParameterFile(out)
			MakeIterationPlot(out[inds,],"flanking_individual_fits")		
		}

	}
	out.init <- rbind(out.init,c(pars.optim,solution$objective))
	converg <- CheckMaxDifference(out.init)
	print(converg)
}
converg <- 1

while (converg >= .0000001) {
	solution <- optim(pars.optim, GetModelFrequencies,
										method="Nelder-Mead" , control=c("maxit"=1000000))

	pars.optim <- solution$par
	out <- rbind(out, c(pars.optim, solution$value))
	converg <- CheckMaxDifference(out)
	print(converg)
	WriteIterationFile(out)
	WriteFinalParameterFile(out)
	MakeIterationPlot(out,"flanking_individual_fits")
}


