
# A0_f <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/4ntflank_5p_unpaired_probs/p/0_8mer.txt",row.names=1,skip=1,header=FALSE)
# A0.4_f <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/4ntflank_5p_unpaired_probs/p/0.4_8mer.txt",row.names=1,skip=1,header=FALSE)
# A4_f <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/4ntflank_5p_unpaired_probs/p/4_8mer.txt",row.names=1,skip=1,header=FALSE)
# I_f <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/4ntflank_5p_unpaired_probs/p/I_8mer.txt",row.names=1,skip=1,header=FALSE)
I_p <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/4ntflank_5p_unpaired_probs/p/I_8mer.txt",row.names=1,header=TRUE)
A0.4_p <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/4ntflank_5p_unpaired_probs/p/0.4_8mer.txt",row.names=1,skip=1,header=FALSE)
A0.4_lp <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/4ntflank_5p_unpaired_probs/logp/0.4_8mer.txt",row.names=1,skip=1,header=FALSE)
A0.4_t <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/4ntflank_5p_unpaired_probs/totals/0.4_8mer.txt",row.names=1,skip=1,header=FALSE)

I_lp <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/4ntflank_5p_unpaired_probs/logp/I_8mer.txt",row.names=1,skip=1,header=FALSE)
I_t <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/4ntflank_5p_unpaired_probs/totals/I_7mer-m8.txt",row.names=1,skip=1,header=FALSE)
I_t_3p <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/4ntflank_3p_unpaired_probs/totals/I_7mer-A1.txt",row.names=1,skip=1,header=FALSE)


# I_lp <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/4ntflank_3p_unpaired_probs/logp/I_8mer.txt",row.names=1,header=TRUE)
# I_p <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/4ntflank_3p_unpaired_probs/p/I_8mer.txt",row.names=1,header=TRUE)
# I_t <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/4ntflank_3p_unpaired_probs/totals/I_8mer.txt",row.names=1,header=TRUE)

# A0.4_lp <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/4ntflank_3p_unpaired_probs/logp/0.4_8mer.txt",row.names=1,header=TRUE)
# A0.4_p <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/4ntflank_3p_unpaired_probs/p/0.4_8mer.txt",row.names=1,header=TRUE)
# A0.4_t <- read.table("/lab/solexa_bartel/mcgeary/AgoRBNS/let-7a/equilibrium/4ntflank_3p_unpaired_probs/totals/0.4_8mer.txt",row.names=1,header=TRUE)
# flank_data_full <- data.frame(I_lp,log10(I_t/sum(I_t))-log10(A0.4_t/sum(A0.4_t)))
# flank_data_full_3p <- data.frame(I_lp_3p,log10(I_t_3p/sum(I_t_3p))-log10(A0.4_t_3p/sum(A0.4_t_3p)))

print(flank_data_full[1,])
print(flank_data_full_3p[1,])

# break
dev.new(height=7,width=15)
par(mfcol = c(3,6))
par(ann=FALSE)
par(lwd=1.5)
par(mai=c(0.5,0.5,0.2,0.2))
###################
y1 = colSums(I_p * matrix(rep(I_t[,1],12),ncol=12))/sum(I_t)
plot(1:12,y1,type="l",ylim=c(0,1),lty=2)
lines(1:8,y1[1:8],lwd=2)
sapply(1:12,function(i) { text(x=i,0,labels=seed_string[i])})
y2 = colSums(A0.4_p *  matrix(rep(A0.4_t[,1],12),ncol=12))/sum(A0.4_t)
lines(1:12,y2,lty=2,col="red")
lines(1:8,y2[1:8],col="red",lwd=2)


plot(1:12,y2/y1,type="l",lty=2,ylim=c(1,1.4))
sapply(1:12,function(i) { text(x=i,1,labels=seed_string[i])})
lines(1:8,(y2/y1)[1:8],lwd=2)

plot(1:12,y2 - y1, lty = 2,type="l",ylim=c(0,0.2))
sapply(1:12,function(i) { text(x=i,0,labels=seed_string[i])})
lines(1:8,(y2-y1)[1:8],lwd=2)


###################

y1 = 10^(colSums(I_lp * matrix(rep(I_t[,1],12),ncol=12))/sum(I_t))
plot(1:12,y1,type="l",ylim=c(0,1),lty=2)
lines(1:8,y1[1:8],lwd=2)
sapply(1:12,function(i) { text(x=i,0,labels=seed_string[i])})
y2 = 10^(colSums(A0.4_lp *  matrix(rep(A0.4_t[,1],12),ncol=12))/sum(A0.4_t))
lines(1:12,y2,lty=2,col="red")
lines(1:8,y2[1:8],lwd=2,col="red")


plot(1:12,y2/y1,type="l",lty=2,ylim=c(1,2.5))
sapply(1:12,function(i) { text(x=i,1,labels=seed_string[i])})
lines(1:8,(y2/y1)[1:8])

plot(1:12,y2 - y1, lty = 2,type="l",ylim=c(0,0.3))
sapply(1:12,function(i) { text(x=i,0,labels=seed_string[i])})
lines(1:8,(y2-y1)[1:8])


###################
y1 = colMeans(I_p)
plot(1:12,y1,type="l",ylim=c(0,1),lty=2)
lines(1:8,y1[1:8])
sapply(1:12,function(i) { text(x=i,0,labels=seed_string[i])})
y2 = colMeans(A0.4_p)
lines(1:12,y2,lty=2,col="red")
lines(1:8,y2[1:8],col="red")


plot(1:12,y2/y1,type="l",lty=2,ylim=c(1,1.4))
sapply(1:12,function(i) { text(x=i,1,labels=seed_string[i])})
lines(1:8,(y2/y1)[1:8])

plot(1:12,y2 - y1, lty = 2,type="l",ylim=c(0,0.2))
sapply(1:12,function(i) { text(x=i,0,labels=seed_string[i])})
lines(1:8,(y2-y1)[1:8])



###################
y1 = 10^colMeans(I_lp)
plot(1:12,y1,type="l",ylim=c(0,1),lty=2)
lines(1:8,y1[1:8])
sapply(1:12,function(i) { text(x=i,0,labels=seed_string[i])})
y2 = 10^colMeans(A0.4_lp)
lines(1:12,y2,lty=2,col="red")
lines(1:8,y2[1:8],col="red")


plot(1:12,y2/y1,type="l",lty=2,ylim=c(1,2.5))
sapply(1:12,function(i) { text(x=i,1,labels=seed_string[i])})
lines(1:8,(y2/y1)[1:8])

plot(1:12,y2 - y1, lty = 2,type="l",ylim=c(0,0.3))
sapply(1:12,function(i) { text(x=i,0,labels=seed_string[i])})
lines(1:8,(y2-y1)[1:8])


###################
y1 <- colSums(I_p * matrix(rep(I_t[,1],12),ncol=12))/sum(I_t)
plot(1:12,y1,type="l",ylim=c(0,1),lty=2)
lines(1:8,y1[1:8],type="l",ylim=c(0,1))

sapply(1:12,function(i) { text(x=i,0,labels=seed_string[i])})
y2 <- colSums(A0.4_p * matrix(rep(I_t[,1],12),ncol=12))/sum(I_t)
lines(1:12,y2,lty=2,col="red")
lines(1:8,y2[1:8],col="red")

y3 <- colSums(I_p * matrix(rep(A0.4_t[,1],12),ncol=12))/sum(A0.4_t)
plot(1:12,y3,type="l",ylim=c(0,1),lty=2)
sapply(1:12,function(i) { text(x=i,0,labels=seed_string[i])})
lines(1:8,y3[1:8],lwd=2)
y4 <- colSums(A0.4_p *  matrix(rep(A0.4_t[,1],12),ncol=12))/sum(A0.4_t)
lines(1:12,y4,lty=2,col="red")
lines(1:8,y4[1:8],lwd=2, col="red")

y = colSums(A0.4_p * matrix(rep(I_t[,1],12),ncol=12))/sum(I_t)-colSums(I_p * matrix(rep(I_t[,1],12),ncol=12))/sum(I_t)
plot(1:12,colSums(A0.4_p * matrix(rep(I_t[,1],12),ncol=12))/sum(I_t)-colSums(I_p * matrix(rep(I_t[,1],12),ncol=12))/sum(I_t),type="l",lty=2,ylim=c(0,0.2))
lines(1:8,y[1:8])
sapply(1:12,function(i) { text(x=i,0,labels=seed_string[i])})
y = colSums(A0.4_p *  matrix(rep(A0.4_t[,1],12),ncol=12))/sum(A0.4_t)-(colSums(I_p *  matrix(rep(A0.4_t[,1],12),ncol=12))/sum(A0.4_t))
lines(1:12,y,col="red")
lines(1:8,y[1:8],col="red")

###################
y1 <- colSums(I_p * matrix(rep(I_t[,1],12),ncol=12))/sum(I_t)
plot(1:12,y1,type="l",ylim=c(0,1),lty=2)
lines(1:8,y1[1:8],type="l",ylim=c(0,1))

sapply(1:12,function(i) { text(x=i,0,labels=seed_string[i])})
y2 <- colSums(A0.4_p * matrix(rep(I_t[,1],12),ncol=12))/sum(I_t)
lines(1:12,y2,lty=2,col="red")
lines(1:8,y2[1:8],col="red")

y3 <- colSums(I_p * matrix(rep(A0.4_t[,1],12),ncol=12))/sum(A0.4_t)
plot(1:12,y3,type="l",ylim=c(0,1),lty=2)
sapply(1:12,function(i) { text(x=i,0,labels=seed_string[i])})
lines(1:8,y3[1:8],lwd=2)
y4 <- colSums(A0.4_p *  matrix(rep(A0.4_t[,1],12),ncol=12))/sum(A0.4_t)
lines(1:12,y4,lty=2,col="red")
lines(1:8,y4[1:8],lwd=2, col="red")

y = colSums(A0.4_p * matrix(rep(I_t[,1],12),ncol=12))/sum(I_t)-colSums(I_p * matrix(rep(I_t[,1],12),ncol=12))/sum(I_t)
plot(1:12,colSums(A0.4_p * matrix(rep(I_t[,1],12),ncol=12))/sum(I_t)-colSums(I_p * matrix(rep(I_t[,1],12),ncol=12))/sum(I_t),type="l",lty=2,ylim=c(0,0.2))
lines(1:8,y[1:8])
sapply(1:12,function(i) { text(x=i,0,labels=seed_string[i])})
y = colSums(A0.4_p *  matrix(rep(A0.4_t[,1],12),ncol=12))/sum(A0.4_t)-(colSums(I_p *  matrix(rep(A0.4_t[,1],12),ncol=12))/sum(A0.4_t))
lines(1:12,y,col="red")
lines(1:8,y[1:8],col="red")
break
# dev.new(height=10,width=12)
par(mfrow=c(3,4))
sapply(1:10,function(i){
	rows <- unique(c(i:(i+2),4,6))
	x <- I_norm[,1]*10^(0.45*rowSums(log10(I_p[,rows])))
	plot(x,A0.4_norm[,1],col=rainbow(256,s=1,v=0.8),lwd=2,log='xy',main=c(cor(log10(x),log10(A0.4_norm[,1])),"\n",paste0(rows,collapse=" ")))
	})

break
data <- A4_f
# data <- data[rownames(A0_f),]
plot(ecdf((data-A0_f)[,1]),do.points=FALSE,lwd=2)
plot(ecdf((data-A0_f)[,2]),do.points=FALSE,add=TRUE,lwd=2,col="gray20")
plot(ecdf((data-A0_f)[,3]),do.points=FALSE,add=TRUE,lwd=2,col="gray45")
plot(ecdf((data-A0_f)[,4]),do.points=FALSE,add=TRUE,lwd=2,col="gray60")
plot(ecdf((data-A0_f)[,5]),do.points=FALSE,add=TRUE,col="red",lwd=2)
plot(ecdf((data-A0_f)[,6]),do.points=FALSE,add=TRUE,col="orangered",lwd=2)
plot(ecdf((data-A0_f)[,7]),do.points=FALSE,add=TRUE,col="green",lwd=2)
plot(ecdf((data-A0_f)[,8]),do.points=FALSE,add=TRUE,col="blue",lwd=2)
plot(ecdf((data-A0_f)[,9]),do.points=FALSE,add=TRUE,col="purple",lwd=2)
plot(ecdf((data-A0_f)[,10]),do.points=FALSE,add=TRUE,col="magenta",lwd=2)
# plot(ecdf((data-A0_f)[,11]),do.points=FALSE,add=TRUE,col="red",lwd=2)
# plot(ecdf((data-A0_f)[,12]),do.points=FALSE,add=TRUE,col="red",lwd=2)
