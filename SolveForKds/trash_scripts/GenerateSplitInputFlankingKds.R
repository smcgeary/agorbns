################################################################################
#GenerateSiteTypeKds.py
################################################################################

# Initial parameters and constants.
args = commandArgs(trailingOnly=TRUE)
mirna = args[1]
site = args[2]

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

GetSiteLeftFlanksXcounts <- function(mirna,site) {
	sites_file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
                       mirna,"/equilibrium/site_count_tables/",
                       site,"_leftflanking.txt")
	siteflanksXcounts <- read.table(sites_file_name)
	return(siteflanksXcounts)
}

GetSiteRightFlanksXcounts <- function(mirna,site) {
	sites_file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/",
                       mirna,"/equilibrium/site_count_tables/",
                       site,"_rightflanking.txt")
	siteflanksXcounts <- read.table(sites_file_name)
	return(siteflanksXcounts)
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
slfXc <- GetSiteLeftFlanksXcounts(mirna,site)
srfXc <- GetSiteRightFlanksXcounts(mirna,site)


sfXc <- round(t((t(sfXc) / colSums(sfXc)) * s.c))
sfXc[,"I"] <- c(sapply(slfXc[,"I"], function(each) { each * srfXc[,"I"]}))/sum(srfXc[,"I"])
colnames(sfXc) <- colnames(sXc)
# Assign parameters for Kds and background, subracting 1 from rows due to "None"
# Kd assigned to 1 nM, and the fact that the I and A0 columns are not modeled
# in the experiment, so there's no background term assigned ot them.
c.freeagos <- 0
pars <- GetSiteKds(mirna,k.c.stockago)
print(pars)

# print(pars)
num.sf <- dim(sfXc)[1]

num.pars <- dim(pars)[1]
bgs <- 10 ^ pars[(num.pars - 5): (num.pars - 1), ]

kds.s <- pars[1:(num.pars - 6), , drop=FALSE]
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
print(plot_range)
# Initialize starting Kds, which are set to 1/the enrichment of each site type
# in the A12.6 experiment.
site_colors = sample(colors(),dim(data)[1])
pars.init <- as.numeric(rep(pars[site,][[1]],num.sf))
# Define function of just kds and pars.
GetModelFrequencies <- function(pars) {
	# Split up the parameters into the kd and background parameters.

	kds <- 10 ^ c(pars,kds.s, 0)
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
		print(-sum(prob))
		setEPS()
		postscript(file=paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/figures/kds/",
	mirna,"/temp_",site,"_altflanking_", k.c.stockago, ".eps"))

		x <- c(40,12.65,4,1.265,0.4)
		y <- c(1,1,1,1,1)

		plot(x,y,log='xy',ylim=c(0.3,200),xlim=c(0.2,50),type="l",col="white")				

		for (i in plot_range) {
			points(x,(data[i,])/colSums(data)/((c.sites[i,])/sum(c.sites)),col=site_colors[i])
			lines(x,(c.final[i,])/((c.sites[i,])/sum(c.sites)),col=site_colors[i])
		}
		dev.off()
	}
	tick <<- tick + 1
	return(-sum(prob))
}
val <- GetModelFrequencies(pars.init)

# Solve the first run of the function, and create the output matrix.
solution <- optim(pars.init, GetModelFrequencies,
									method="Nelder-Mead")
out <- rbind(c(pars.init, val),
						 c(solution$par, solution$value))

colnames(out) <- c(rownames(sfXc),"-logp")

# Get maximum difference in output.
converg <- CheckMaxDifference(out)
print(converg)

# Assign stopping criteria for continuing the optimization.
while (converg >= 0.0001) {
	solution <- optim(solution$par, GetModelFrequencies,
										method="Nelder-Mead")
	out <- rbind(out,c(solution$par,solution$value))
	converg <- CheckMaxDifference(out)
	# Assign output file for the entire sequnce of the optimization and write
	# to it.
	out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
		     						 "/equilibrium/kds/",site,"_altflanking_", k.c.stockago, ".txt")
	write.table(file=out_file, out, sep="\t", quote=FALSE, row.names=FALSE,
							col.names=TRUE)

	# Assign output file For the final parameters and write to it.
	out_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
		     						 "/equilibrium/kds/final_",site,"_altflanking_", k.c.stockago, ".txt")
	out_final <- out[dim(out)[1],]
	names(out_final) <- colnames(out)

	write.table(file=out_file, out_final, sep="\t", quote=FALSE, row.names=TRUE,
							col.names=FALSE)

	# Printing figures.
	# Print the convergence picture:
	setEPS()
	postscript(file=paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/figures/kds/",
		mirna,"/all_",site,"_altflanking_", k.c.stockago, ".eps"))
	x = c(1:(dim(out)[1]))
	probs = out[ ,dim(out)[2]]
	probs = probs / max(probs)

	out_print <- out[ , 1:(dim(out)[2]-1)]
	ys <- c(min(out_print),max(out_print))
	probs <- probs * (ys[2]-ys[1]) + ys[1]


	plot(x,probs,axes=FALSE,type="l",ylim=ys,
		main="",xlab="Iteration", ylab="Normalized parameter estimate",lwd=2,col="red")
	axis(1,at=x,pos=ys[1],lwd=2)
	axis(2,at=seq(ys[1],ys[2],length=5),pos=1,lwd=2)
	apply(out_print,2, function(y) {
		lines(x,y,lwd=1,col="blue")
		})
	lines(x,probs,type="l",lwd=2,col="red")
	dev.off()
}


