information_content <- function(row){
	row <- row/sum(row)
	return(-sum(sapply(row,function(p){p*log(p,base=length(row))}),na.rm=TRUE))
}


average_position <- function(row){
	print(row)
	row <- row/sum(row)
	print(row)
	print(seq(length(row)))
	return(sum(sapply(seq(length(row)),function(i){row[i]*i}),na.rm=TRUE))
}
plot(0:1,0:1,xlim=c(1,24),ylim=c(-1,20),col="white",axes=FALSE,ann=FALSE)
axis(side=1,pos=-1,at=c(seq(1,22)),labels=seq(22,1))
labels_mirna <- strsplit("UUGAUAUGUUGGAUGAUGGAGU",split="")[[1]]
cols = rep("black",23)
cols[7:10] <- "blue"
cols[15:22] <- "red"
sapply(1:22,function(i){
	text(x = i, y = 0, labels_mirna[i],cex=2.6,col=cols[i],family="mono")
	})
sapply(1:20,function(i){
	ind <- order(-ratio_total)[i]
	x <- ic_pos_A[ind]
	sapply(1:5,function(pos){
	text(x+2+pos,i,strsplit(names(ic_pos_A)[ind],split="")[[1]][pos],cex=2.6,family="mono")

	})
	text(1,i,round(ratio_total[ind],2),cex=1.5,pos=4)

	})