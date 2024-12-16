################################################################################
#Analyzestructure.py
################################################################################
# library(gplots)
source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")
library(wrswoR)
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


s_I <- get_structures(mirna, site, "I")
s_0.4 <- get_structures(mirna, site, "0.4")
s_1.26 <- get_structures(mirna, site, "1.26")
s_4 <- get_structures(mirna, site, "4")
s_12.6 <- get_structures(mirna, site, "12.6")
s_40 <- get_structures(mirna, site, "40")


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



get_read_positions <- function(data,read_data,position) {
	apply(data, 1, function(row) {
		start <- as.numeric(row[1])
		# print(as.numeric(row[(start+5+temp_start):(start+5+7+temp_stop)]))
		# print(1-as.numeric(row[(start+5+temp_start):(start+5+7+temp_stop)]))
		# print(10^sum(log10(1-as.numeric(row[(start+5+temp_start):(start+5+7+temp_stop)]))))

		as.numeric(row[start+5+position-1])
		})
}
pos <- 26:55
data <- s_0.4_8mer
Ind_A <- which(as.integer(data[, 1]) %in% pos)
flanks_A <- get_site_flanks(data[Ind_A,])
plot(seq(-20,20),c(seq(-20,19)*0,20),col="white")
# seed_string <- strsplit("NNATACAAAANN",split="")[[1]]
sapply(-2:9,function(i) { text(x=i,0.1,labels=seed_string[i+3])})
segments(-0.5,0,-0.5,20,lwd=2,lty=2)
segments(7.5,0,7.5,20,lwd=2,lty=2)
segments(-0.5,0,-0.5,20,lwd=2,lty=2)
segments(-2.5,0,-2.5,20,lwd=1,lty=2)
segments(9.5,0,9.5,20,lwd=1,lty=2)
# legend("topleft",legend=seq(20),pch=1,lwd=2,col=rainbow(20),bty="n")

tot <- 3
inds_all <- rep(seq(1,(dim(s_I_8mer)[1])),tot)

total_flanks_A <- sapply(unique(flanks_A), function(flank) {
		length(which(flanks_A==flank))/length(Ind_A)
		})

# fit_flank_parameters <- function(flank){
# 	optim(c(0,0),function(pars){
# 		data <- ecdf((I^(1/9))[which(flanks_temp==flank)])(x)
# 		model <- 
# 		})
# 	x,ecdf((I^(1/9))[which(flanks_temp==flank)])(x)
# }
final_mat <- matrix(0,ncol=38,nrow=20)
for (j in seq(20)){
	x_tot <- c()
	y_tot <- c()
for (i in seq(-17,20)) {
	left <- i-floor((j-1)/2)
	right <- left+j-1
	x_now <- mean(c(left, right))
	x_tot <- c(x_tot,x_now)
	print(c(left,right,x_now))

	inds <- inds_all[sample_int_R(length(inds_all),round(length(inds_all)*0.60),prob=get_site_probs(s_I_8mer[inds_all,],left,right))]

	s_0.4_8mer_sample <- s_I_8mer[inds,]	
	print("post sample")
	Ind_A_s <- which(as.integer(s_0.4_8mer_sample[, 1]) %in% pos)
	# means_A_s <- get_site_probs(s_0.4_8mer_sample[Ind_A_s, ],left,right)
	flanks_A_s <- get_site_flanks(s_0.4_8mer_sample[Ind_A_s,])
# 	print(flanks_A_s[1:10])
# 	print(flanks_A[1:10])
# 	print(unique(flanks_A_s)[1:10])
# 	# print(unique(flanks_A)[1:10])
	# flank_means_A_s <- sapply(unique(flanks_A), function(flank) {
	# 	mean(means_A_s[which(flanks_A_s==flank)])
	# 	})
	total_flanks_A_s <- sapply(unique(flanks_A), function(flank) {
		length(which(flanks_A_s==flank))/length(Ind_A_s)
		})
	# flank_means_A <- sapply(unique(flanks_A), function(flank) {
	# 	mean(means_A[which(flanks_A==flank)])
	# 	})

# 		# mean_flanks_predict <- sapply(unique(flanks_A_s), function(flank) {
# 		# mean(means_A_s[which(flanks_A[temp_inds]==flank)])/length(temp_inds)
# 		# })

	y_new <- cor(log(total_flanks_A_s)[which(total_flanks_A_s>0&total_flanks_A>0)],log(total_flanks_A)[which(total_flanks_A_s>0&total_flanks_A>0)])
	final_mat[j,i] <- y_new
	y_tot <- c(y_tot, y_new)
	print(y_new)
	# par(mfrow=c(2,3))
	# flank_cols <- sapply(unique(flanks_A_s), function(x) {GetColorFunction(paste0(c(unlist(strsplit(x,split=""))[1:4]),collapse=""))})
	# plot(flank_means_A_s,flank_means_A,pch=19,cex=2,col=flank_cols,log='xy')
	# title(main=left)

	# abline()
	# plot(total_flanks_A_s,total_flanks_A,pch=19,cex=2,col= flank_cols)
	# title(main=right)

	# plot(flank_means_A/flank_means_A_s,total_flanks_A/total_flanks_A_s,pch=19,cex=2,col= flank_cols,log='xy')
	# plot(total_flanks_A,total_flanks_A/total_flanks_A_s,pch=19,cex=2,col= flank_cols,log='xy')
	# plot(flank_means_A*total_flanks_A_s,total_flanks_A,pch=19,cex=2,col= flank_cols,log='xy')
	# title(main=cor(log(flank_means_A*total_flanks_A_s),log(total_flanks_A)))
	# flank_cols <- sapply(unique(flanks_A), function(x) {GetColorFunction(paste0(c(unlist(strsplit(x,split=""))[1:4]),collapse=""))})
	# dev.new(mfrow
	# legend_text <- c("A","C","G","T")
	# legendcols <- c("AAAA","CCCC","GGGG","TTTT")
	# # points(x_now,y_new,col=rainbow(20)[j],pch=19)
				 rect(x_now-0.5,j-0.5,x_now+0.5,j+0.5,col=rainbow(100)[100-round(y_new*99)],
				 border=NA)

	# title(main=cor(log(total_flanks_A_s),log(total_flanks_A)))
}
	# lines(x_tot,y_tot,col=rainbow(20)[j],lwd=2)

}

	# indsfinal <- inds_all[sample(inds_all,round(length(inds_all)*0.40),replace=FALSE,prob=get_site_probs(s_I_8mer[inds_all,],-2,-1)*get_site_probs(s_I_8mer[inds_all,],8,9))]
	# s_0.4_8mer_sample <- s_I_8mer[indsfinal,]	
	# print("post sample")
	# Ind_A_s <- which(as.integer(s_0.4_8mer_sample[, 1]) %in% pos)
	# flanks_A_s <- get_site_flanks(s_0.4_8mer_sample[Ind_A_s,])
	# total_flanks_A_s <- sapply(unique(flanks_A), function(flank) {
	# 	length(which(flanks_A_s==flank))/length(Ind_A_s)
	# 	})



# }
# graphics.off()
# print("done")
# s_0_8mer <- get_structures(mirna, site, "0")
# print("done")
# print("done")
# s_1.26_8mer <- get_structures(mirna, site, "1.26")
# print("done")
# s_4_8mer <- get_structures(mirna, site, "4")
# print("done")
# s_12.6_8mer <- get_structures(mirna, site, "12.6")
# print("done")
# s_40_8mer <- get_structures(mirna, site, "40")
# print("done")
# dev.new()
# dev.new()

# heatmap(as.matrix(s_I_8mer[, -1]),Colv=NA,Rowv=NA)
# dev.new()
# heatmap(as.matrix(s_0_8mer[, -1]),Colv=NA,Rowv=NA)
# dev.new()
# heatmap(as.matrix(s_0.4_8mer[, -1]),Colv=NA,Rowv=NA)
# dev.new()
# heatmap(as.matrix(s_1.26_8mer[, -1]),Colv=NA,Rowv=NA)
# dev.new()
# heatmap(as.matrix(s_4_8mer[, -1]),Colv=NA,Rowv=NA)
# dev.new()
# heatmap(as.matrix(s_12.6_8mer[, -1]),Colv=NA,Rowv=NA)
# dev.new()
# heatmap(as.matrix(s_40_8mer[, -1]),Colv=NA,Rowv=NA)

# heatmap.2(as.matrix(rbind(s_1.26_8mer[sample(which(s_1.26_8mer[,1]==40),1000),-1],
# 	                      s_I_8mer[which(s_I_8mer[,1]==40),-1])),
# 						  density.info="none",
# 						  trace="none",
# 						  margins=c(12,9),
# 						  dendrogram="row",
# 						  Colv="NA",
# 						  RowSideColors= c(rep("blue",1000),rep("red",length(which(s_I_8mer[,1]==40)))))
# I_ind <- grep("8mer:",read.table(
# 	"/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/full_sites/I_5-5.txt",
# 	header=FALSE,nrows = 400000,sep="\t",
# 	stringsAsFactors = FALSE)[, 1])
# I_pos <- grep("8mer:",read.table(
# 	"/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/full_sites/I_5-5.txt",
# 	header=FALSE,nrows = 400000,sep="\t",
# 	stringsAsFactors = FALSE)[, 1],
# 	value = TRUE)
# I_pos_real <- sapply(I_pos,function(x) {
# 	pos_full <- unlist(strsplit(x,"8mer:"))[2]
# 	pos_full <- unlist(strsplit(pos_full,","))[1]
# 	pos <- unlist(strsplit(pos_full,"-"))

# 	return(pos)
# 	})
# I_pos_final <- matrix(as.numeric(I_pos_real),ncol=2,byrow=TRUE)
# I_reads <- read.table(
# 		"/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/full_reads/I.txt",
# 	header=FALSE,nrows = 400000,sep="\t",
# 	stringsAsFactors = FALSE)[I_ind, 1]
# I_left_1 <- sapply(1:length(I_ind),function(x){
# 	substring(I_reads[x],22-2+I_pos_final[x,1],22-2+I_pos_final[x,1])
# 	})
# I_left_2 <- sapply(1:length(I_ind),function(x){
# 	substring(I_reads[x],22-1+I_pos_final[x,1],22-1+I_pos_final[x,1])
# 	})

# I_right_1 <- sapply(1:length(I_ind),function(x){
# 	substring(I_reads[x],22+1+I_pos_final[x,2],22+1+I_pos_final[x,2])
# 	})
# I_right_2 <- sapply(1:length(I_ind),function(x){
# 	substring(I_reads[x],22+2+I_pos_final[x,2],22+2+I_pos_final[x,2])
# 	})
# # flanks.df <- data.frame(pairing = I_8_site_means,
# # 						pos = I_pos_final[,1],
# # 						flank.1.5p = I_left_2,
# # 						flank.2.5p = I_left_1,
# # 						flank.1.3p = I_right_1,
# # 						flank.2.3p = I_right_2)
# # print("4")
# A_ind <- grep("8mer:",read.table(
# 	"/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/full_sites/4_5-5.txt",
# 	header=FALSE,nrows = 400000,sep="\t",
# 	stringsAsFactors = FALSE)[, 1])
# print("done")
# A_pos <- grep("8mer:",read.table(
# 	"/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/full_sites/4_5-5.txt",
# 	header=FALSE,nrows = 400000,sep="\t",
# 	stringsAsFactors = FALSE)[, 1],
# 	value = TRUE)
# A_pos_real <- sapply(A_pos,function(x) {
# 	pos_full <- unlist(strsplit(x,"8mer:"))[2]
# 	pos_full <- unlist(strsplit(pos_full,","))[1]
# 	pos <- unlist(strsplit(pos_full,"-"))
# 	return(pos)
# 	})
# A_pos_final <- matrix(as.numeric(A_pos_real),ncol=2,byrow=TRUE)

# I_file_name <- "/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/structures_bp_prob/I.txt"
# A_file_name <- "/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/structures_bp_prob/4.txt"
# I_struc <- read.table(
# 	I_file_name,
# 	header=FALSE,nrows = 400000)

# A_struc <- read.table(
# 	A_file_name,
# 	header=FALSE,nrows = 400000)

# # A_reads <- read.table(
# # 		"/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/full_reads/I.txt",
# # 	header=FALSE,nrows = 400000,sep="\t",
# # 	stringsAsFactors = FALSE)[A_ind, 1]
# # A_left_1 <- sapply(1:length(A_ind),function(x){
# # 	substring(A_reads[x],22-2+A_pos_final[x,1],22-2+A_pos_final[x,1])
# # 	})
# # A_left_2 <- sapply(1:length(A_ind),function(x){
# # 	substring(A_reads[x],22-1+A_pos_final[x,1],22-1+A_pos_final[x,1])
# # 	})

# # A_right_1 <- sapply(1:length(A_ind),function(x){
# # 	substring(A_reads[x],22+1+A_pos_final[x,2],22+1+A_pos_final[x,2])
# # 	})
# # A_right_2 <- sapply(1:length(A_ind),function(x){
# # 	substring(A_reads[x],22+2+A_pos_final[x,2],22+2+A_pos_final[x,2])
# # 	})
# # Aflanks.df <- data.frame(pairing = A_8_site_means,
# # 						pos = A_pos_final[,1],
# # 						flank.1.5p = A_left_2,
# # 						flank.2.5p = A_left_1,
# # 						flank.1.3p = A_right_1,
# # 						flank.2.3p = A_right_2)


# mir1_word  <- unlist(strsplit("GGGCAGAGUUCUACAGUCCGACGAUCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNUCGUAUGCCGUCUUCUGCUUG",split=""))
# other_word <- unlist(strsplit("GGGCAGAGUUCUACAGUCCGACGAUCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNUGUUCGUAUGCCGUCUUCUGCUUG",split=""))

# par(mfrow=c(2,2))
# plot(seq(1,87),
# 	A_struc[A_ind[1],],
# 	type = "p",
# 	col  = "white",
# 	axes = FALSE,
# 	ann  = FALSE)
# axis(side   = 1,
# 	 at     = seq(1, 87, 1),
# 	 labels = FALSE,
# 	 pos    = 0)
# axis(side   = 2,
# 	 at     = seq(0, 1, 0.1),
# 	 pos    = 0)
# sapply(1:100,function(x){
# 	if (x%%10==0){
# 		print(x)
# 	}
# 	points(seq(1,87),A_struc[A_ind[x],],col=rgb(0,0,0,alpha=0.1),pch=20)
# 	})

# plot(seq(-40,40),
# 	c(rep(0,80),1),
# 	type = "p",
# 	col  = "white",
# 	axes = FALSE,
# 	ann  = FALSE)
# axis(side   = 1,
# 	 at     = seq(1, 87, 1),
# 	 labels = unlist(strsplit(other_word, split = "")),
# 	 pos    = 0)
# axis(side   = 2,
# 	 at     = seq(0, 1, 0.1),
# 	 pos    = -40)
# rect(0,0,7,1,border=NA,col="gray60")
# sapply(1:100,function(x){
# 	if (x%%100==0){
# 		print(x)
# 	}
# 	points(seq(0,37)-A_pos_final[x,1],A_struc[A_ind[x],seq(26,26+37)],col=rgb(1,0,0,alpha=0.5),pch=20) 
# 	})

# plot(seq(1,87),
# 	A_struc[A_ind[1],],
# 	type = "p",
# 	col  = "white",
# 	axes = FALSE,
# 	ann  = FALSE)
# axis(side   = 1,
# 	 at     = seq(1, 87, 1),
# 	 labels = unlist(strsplit(other_word, split = "")),
# 	 pos    = 0)
# axis(side   = 2,
# 	 at     = seq(0, 1, 0.1),
# 	 pos    = 0)
# axis(side   = 1,
# 	 at     = seq(1, 87, 1),
# 	 labels = unlist(strsplit(other_word, split = "")),
# 	 pos    = 0)
# axis(side   = 2,
# 	 at     = seq(0, 1, 0.1),
# 	 pos    = 0)
# sapply(1:100,function(x){
# 	if (x%%100==0){
# 		print(x)
# 	}
# 	points(seq(1,87),I_struc[I_ind[x],],col=rgb(0,0,0,alpha=0.1),pch=20)
# 	})
# plot(seq(-40,40),
# 	c(rep(0,80),1),
# 	type = "p",
# 	col  = "white",
# 	axes = FALSE,
# 	ann  = FALSE)
# axis(side   = 1,
# 	 at     = seq(1, 87, 1),
# 	 labels = unlist(strsplit(other_word, split = "")),
# 	 pos    = 0)
# axis(side   = 2,
# 	 at     = seq(0, 1, 0.1),
# 	 pos    = -40)
# sapply(1:100,function(x){
# 	if (x%%100==0){
# 		print(x)
# 	}
# 	points(seq(0,37)-I_pos_final[x,1],I_struc[I_ind[x],seq(26,26+37)],col=rgb(1,0,0,alpha=0.5),pch=20) 
# 	})

# I_8_means <- rowMeans(I_struc[I_ind,])
# I_8_site_means <- sapply(1:length(I_ind),function(ind){
# 	mean(unlist(I_struc[I_ind[ind],(22+I_pos_final[ind,1]):(22+I_pos_final[ind,2])]))
# 	})

# A_8_means <- rowMeans(A_struc[A_ind,])
# A_8_site_means <- sapply(1:length(A_ind),function(ind){
# 	mean(unlist(A_struc[A_ind[ind],(22+A_pos_final[ind,1]):(22+A_pos_final[ind,2])]))
# 	})
# par(mfrow=c(1,1))
# com = sample(seq(length(A_8_means)), size=length(I_8_means))
# plot(I_8_means, I_8_site_means, xlim=c(0, 1), ylim=c(0, 1))
# plot(A_8_means[com], A_8_site_means[com], xlim=c(0, 1), ylim = c(0, 1))
