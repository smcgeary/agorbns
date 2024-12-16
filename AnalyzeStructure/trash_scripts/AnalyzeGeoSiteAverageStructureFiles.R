################################################################################
#AnalyzeSiteAverageStructureFiles.R
################################################################################

# # Initial parameters and constants.
library(colorspace)
site_cols <- read.table(
  "/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/site_colors.txt",
  row.names=1, header=FALSE, sep="\t", stringsAsFactors=FALSE)
site_cols <- read.table(
  "/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/site_colors.txt",
  row.names=1, header=FALSE, sep="\t", stringsAsFactors=FALSE)
# args = commandArgs(trailingOnly=TRUE)
# mirna = args[1]
# site = args[2]


source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")

getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

fix_names <- function(x){
				temp_ <- unlist(strsplit(x,split=""))
				if (temp_[1]=="X") {
					x <- paste0(temp_[-1],collapse="")
				}
				temp_2 <- unlist(strsplit(x,"mer\\."))
				if (length(temp_2)==2) {
					x <- paste0(temp_2,collapse="mer-")
				}
				return(x)
			}
color_palette <- rainbow_hcl(n = 100,
				 			c = 90,
				 			l = 60,
				 			start = 0,
				 			end = 4 / 6 * 360)
col_ind <- function(x) {
				 	if (is.na(x)) {
			 			return("white")
			 		} else {
			 			return(round(min(max((x)/-1.9*(150)+110,1),100)))
				 	}
				}
get_color <- function(x) {
	return(color_palette[col_ind(x)])
}

mag<-1.5
# # setPDF(width=20,height=15)
# pdf(file=paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/figures/structures/",
# 	"/geo_average_structure_enrichments/",mirna,"/",site,".pdf"),width=20,height=15)

# par(mfrow=c(3,4),mai = c(1, 0.1, 0.1, 0.1))
# par(par_plots)
conditions <- c("0","0.4","1.26","4","12.6","40")
conditions <- c("0.4")
# dev.new()
sapply(conditions,function(condition) {
	xlim = c(-19,1)
	# ylim = c(0.5,20.5)
	# plot(c(1,10),c(1,10),type="l",xlim=xlim,ylim=ylim,axes=FALSE,col="white")
	# sapply(seq(1,20),function(win_size) {
	# 	print(win_size)
	# 		window_I <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
	# 	                  "/equilibrium/structures_bp_prob/","I",
	# 	                  "_0-0_geo_average_window_",win_size,".txt")
	# 		total_I <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
	# 	                  "/equilibrium/structures_bp_prob/","I",
	# 	                  "_0-0_geo_totals_window_",win_size,".txt")
	# 		probs_I <- read.table(window_I,header=TRUE,sep="\t",stringsAsFactors=FALSE,row.names=1)
	# 		totals_I <- read.table(total_I,header=TRUE,sep="\t",stringsAsFactors=FALSE,row.names=1)


	# 		window_AGO <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
	# 	                  "/equilibrium/structures_bp_prob/",condition,
	# 	                  "_0-0_geo_average_window_",win_size,".txt")
	# 		total_AGO <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
	# 	                  "/equilibrium/structures_bp_prob/",condition,
	# 	                  "_0-0_geo_totals_window_",win_size,".txt")

	# 		probs_AGO <- read.table(window_AGO,header=TRUE,sep="\t",stringsAsFactors=FALSE,row.names=1)
	# 		totals_AGO <- read.table(total_AGO,header=TRUE,sep="\t",stringsAsFactors=FALSE,row.names=1)
	# 		colnames(probs_I) <- sapply(colnames(probs_I),fix_names)
	# 		print("colnames I")
	# 		print(colnames(probs_I))
	# 		colnames(probs_AGO) <- sapply(colnames(probs_AGO),fix_names)
	# 		colnames(totals_I) <- sapply(colnames(probs_I),fix_names)
	# 		colnames(totals_AGO) <- sapply(colnames(probs_AGO),fix_names)
	# 		print("colnames AGO")
	# 		print(colnames(probs_AGO))
	# 		print("probs AGO")
	# 		print(probs_AGO)
	# 		print("probs_I")
	# 		print(probs_I)
	# 		x_dat = as.numeric(rownames(probs_I))
	# 		inds = which(x_dat >= xlim[1]+1 & x_dat <= xlim[2])
	# 		print(x_dat[inds])
	# 		print(x_dat)
	# 		data = (probs_AGO[,site]-probs_I[,site])/log10(2)
	# 		print("data")
	# 		print(data[inds])
	# 		 rect(x_dat[inds]-0.5,rep(win_size-0.5,length(inds)),x_dat[inds]+0.5,rep(win_size+0.5,length(inds)),col=sapply(data[inds],get_color),
	# 			 border=NA)

	# 	})
	# 	axis(side=1,at=seq(xlim[1],xlim[2]),tck=0,labels=FALSE,lwd=2,pos=ylim[1])
	# 	axis(side=1,at=c(-21,-15,-12,-7,0),labels=FALSE,lwd=2,pos=ylim[1])
	# 	axis(side=1,at=-23:0,labels=FALSE,lwd=1,pos=ylim[1])

	# 	axis(side=1,at=c(-13.5,-3.5),labels=c("3'","seed"),tck=0,pos=ylim[1])

	# 	axis(side=2,at=seq(ylim[1],ylim[2]),labels=FALSE,lwd=2,pos=xlim[1])
	# 	axis(side=2,at=seq(ylim[1],ylim[2],by=5),lwd=2,pos=xlim[1])
		ylim=c(0,1.5)
		plot(c(1,10),c(1,10),type="l",xlim=xlim,ylim=ylim,axes=FALSE,col="white")
		sapply(seq(1),function(win_size) {
			print(win_size)
			window_I <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
		                  "/equilibrium/structures_bp_prob/","I",
		                  "_0-0_geo_average_window_",win_size,".txt")
			total_I <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
		                  "/equilibrium/structures_bp_prob/","I",
		                  "_0-0_geo_totals_window_",win_size,".txt")
			probs_I <- read.table(window_I,header=TRUE,sep="\t",stringsAsFactors=FALSE,row.names=1)
			totals_I <- read.table(total_I,header=TRUE,sep="\t",stringsAsFactors=FALSE,row.names=1)


			window_AGO <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
		                  "/equilibrium/structures_bp_prob/",condition,
		                  "_0-0_geo_average_window_",win_size,".txt")
			total_AGO <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
		                  "/equilibrium/structures_bp_prob/",condition,
		                  "_0-0_geo_totals_window_",win_size,".txt")

			probs_AGO <- read.table(window_AGO,header=TRUE,sep="\t",stringsAsFactors=FALSE,row.names=1)
			totals_AGO <- read.table(total_AGO,header=TRUE,sep="\t",stringsAsFactors=FALSE,row.names=1)
			colnames(probs_I) <- sapply(colnames(probs_I),fix_names)
			colnames(probs_AGO) <- sapply(colnames(probs_AGO),fix_names)
			colnames(totals_I) <- sapply(colnames(probs_I),fix_names)
			colnames(totals_AGO) <- sapply(colnames(probs_AGO),fix_names)
			x_dat = as.numeric(rownames(probs_I))
			data = (probs_AGO[,site]-probs_I[,site])/log10(2)
			lines(x_dat,data,lwd=2)

		})
		axis(side=1,at=seq(xlim[1],xlim[2]),tck=0,labels=FALSE,lwd=2,pos=ylim[1])
		axis(side=1,at=c(-21,-15,-12,-7,0),labels=FALSE,lwd=2,pos=ylim[1])
		axis(side=1,at=-23:0,labels=FALSE,lwd=1,pos=ylim[1])

		axis(side=1,at=c(-13.5,-3.5),labels=c("3'","seed"),tck=0,pos=ylim[1])

		axis(side=2,at=seq(ylim[1],ylim[2],length=5),lwd=2,tck=-0.025,pos=xlim[1])
		# axis(side=4,at=seq(ylim[1],ylim[2],length=10),labels <- sapply(seq(ylim[1],ylim[2],length=10), col_ind),lwd=2,tck=-0.025,pos=xlim[2])

		yrange = ylim[2] - ylim[1]
		y_step = yrange / 30
		# rect(rep(xlim[1], 30),
		# 	 seq(ylim[1], ylim[2], length = 30),
		# 	 rep(xlim[1], 30) + 2,
		# 	 seq(ylim[1], ylim[2], length = 30) + y_step,
		# 	 col=sapply(seq(ylim[1], ylim[2], length=30), get_color),
		# 	 border = NA)

})

		# dev.off()
