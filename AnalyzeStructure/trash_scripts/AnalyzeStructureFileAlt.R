################################################################################
#Analyzestructure.py
################################################################################
# library(gplots)
source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")
# library(wrswoR)
get_structures <- function(mirna, site, condition) {
	file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
					    "/equilibrium/structures_bp_prob/", site, "/", condition, "_0-0.txt")
		out <- read.table(file_name,header=FALSE,skip=1,sep="\t",stringsAsFactors=FALSE)
	return(out)
}

get_reads <- function(mirna, site, condition) {
	file_name <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna,
					    "/equilibrium/reads_by_site/", site, "/", condition, "_0-0.txt")
		out <- read.table(file_name,header=FALSE,skip=1,sep="\t",stringsAsFactors=FALSE)
	return(out)
}

get_site_flanks <- function(data) {
	apply(data, 1, function(row) {
		return(paste0(row[2:5],collapse=""))
		})
}


# s_I_8mer <- get_structures(mirna, site, "I")
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

