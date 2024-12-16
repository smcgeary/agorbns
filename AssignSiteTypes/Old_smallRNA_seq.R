## sRS_analyze_read_data
library(som)
graphics.off()

## FUNCTIONS for SCRIPT:

# Removes rows from expression matrix in which there is not a single read.
Filter_zeros <-function (mat) {
	return(mat[which(rowSums(mat)>0),])
	}

# Creates a log(10) transformation of the data, adding a pseudocount of + 1 to each read.	
Lognorm <- function (X) {log(X+1,10)}

# Generates a matrix in which each value is converted to a percentage of the sum of that column.
Col_frac <- function(mat) {
	# print("Col_frac")
	output<- as.matrix(apply(mat,1,function(row) {
		# print(row)
		# print(colSums(mat))
		row_out <- row/colSums(mat)
		# return(row_out)
	}),
	nrow=dim(mat)[1],ncol=dim(mat[2],byrow=TRUE))
	return(t(output))
}

Col_div <- function(mat,tot) {
	# print("Col_frac")
	output<- as.matrix(apply(mat,1,function(row) {
		# print(row)
		# print(colSums(mat))
		row_out <- row/tot
		# return(row_out)
	}),
	nrow=dim(mat)[1],ncol=dim(mat[2],byrow=TRUE))
	return(t(output))
}

# Row normalizes a matrix "mat"
Row_norm <- function(mat) {
	# print("Row_norm")
	output<- as.matrix(apply(mat,1,function(row) {
		row_out <- row/mean(row)
		return(row_out)
	}),
	nrow=dim(mat)[1],ncol=dim(mat[2],byrow=TRUE))
	return(t(output))
}

# Generates a matrix from another matrix, sorted numerically according to a chosen column
Sort <- function(mat,col) {
	# print("Sort")
	return(mat[order(as.vector(mat[,col]),decreasing = TRUE),])
}

# Generates a matrix in which each row is the sum, along each column, of all the values from the first row to that row. If the columns of the input matrix sum to 1, then this will produce a cumulative function from 0 to 1 along each column.
Cum_col <- function(mat,col) {
	print("Cum_col")
	out <- as.matrix(
		sapply(rownames(mat),function(row)  {
			ind <- which(rownames(mat)==row)		
			# print(ind)
			x <- mat[ind,col]
			# print("x")
			# print(x)
			# print(sum(mat[1:ind,col]))
			return(sum(mat[1:ind,col]))
			}
			),ncol=1, nrow=length(mat[,col]), byrow=FALSE)
	names(out) <- rownames(mat)
	return(out)
}

# Generates a matrix in which each row is the sum, along each column, of all the values from the first row to that row. If the columns of the input matrix sum to 1, then this will produce a cumulative function from 0 to 1 along each column.
Cum_mat <- function(mat) {
	print("Cum_mat")
	print(dim(mat))
	mat <- Col_frac(mat)
	out <- as.matrix(
		sapply(rownames(mat)[-1],function(row)  {
			ind <- which(rownames(mat)==row)		
			
			x <- as.matrix(mat[1:ind,],ncol=dim(mat)[2])

			return(colSums(x))
			}
			),ncol=dim(mat)[2], nrow=dim(mat)[1]-1, byrow=TRUE)
	print("out:")
	out <- t(out)
	out <- rbind(mat[1,],out)
	rownames(out) <- rownames(mat)
	# print(out[1:10])
	# print(out[1:10,])
	return(out)
}

# This function will remove the bottom segment of a row, depending on a percentage chosen as a "cutoff" from a cumulative distribution function.
Remove_frac <- function(mat,perc,col) {
	print("Remove_frac")
	orig <- mat
	sort <- Sort(mat,col)
	print("sorted:")
	print(sort[1:10,])
	mat_sorted <- Cum_mat(sort)
	# print("82 in remove frac")
	# print(mat_sorted[1:10,])
	ind <- max(which(as.matrix(mat_sorted)[,col]<=perc))
	# print(ind)
	# print(Sort(mat,col))
	mat_new <- Sort(mat,col)[1:ind,]
	mat_l <- dim(mat)[1]
	# print(mat_l)
	rem <- Sort(mat,col)[(ind+1):(dim(mat)[1]),]
	# print("rem")
	# print(rem)
	# print(colSums(rem))
	mat_new <- rbind(mat_new,colSums(rem))
	rownames(mat_new)[ind+1] <- "remaining miRNAs"
	return(mat_new)
}










D <- read.table("140603_WIGTC-HISEQA_C4W8TACXX_expression_matrix.txt", row.names=1, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors=FALSE)

D <- D[-1,]
D <- data.matrix(D)
print(head(D))
print("break41")
# Removes the REF column from the data, places the miR name as a row name. Then removes the "unknown" column, and divies up the matrix in to totals ("D_TOT"), size markers ("D_size"), quantification markers ("D_q"), capture / competitor oligos ("D_capt"), and mapped miRNAs ("D_miR")
D <- D[, order(as.integer((colnames(D))))]
D <- D[,-10]

D_T <- D[seq((nrow(D)-1),nrow(D),1),]
D <- D[-seq((nrow(D)-1),nrow(D),1),]

D_size <- D[1:2,]
D_q <- D[3:5,]
D_capt <- D[6:9,]
D_miR <- Filter_zeros(D[10:nrow(D),])

# Total assignment
T_size <- colSums(D_size)
T_q <- colSums(D_q)
T_capt <- colSums(D_capt)
T_miR <- colSums(D_miR)
T_all <- colSums(rbind(T_size,T_q,T_capt,T_miR))
T_lib <- colSums(rbind(T_q,T_capt,T_miR))




Pt_size <- T_size/T_all
Pt_q <- T_q/T_all
Pt_capt <- T_capt/T_all
Pt_miR <- T_miR/T_all
Pt_lib <- T_lib/T_all

Pl_q <- T_q/T_lib
Pl_capt <- T_capt/T_lib
Pl_miR <- T_miR/T_lib
Pl_all <- T_all/T_lib

miR_q <- matrix(T_miR/T_q,nrow=1)

break
rownames(miR_q) <- "miRNA/quant read ratio"
Percents_tot <- rbind(Pt_size,Pt_q,Pt_capt,Pt_miR,Pt_lib)
rownames(Percents_tot) <- c("% size markers","% quant standards","% capt/comp oligo","% human miRNAs","% total mapped")
Percents_lib <- rbind(Pl_q,Pl_capt,Pl_miR)
rownames(Percents_lib) <- c("% quant standards","% capt/comp oligo","% human miRNAs")

print("Percents_tot")
print(Percents_tot)
print("Percents_tot, normalized")
print(Row_norm(Percents_tot))


barplot(as.matrix(Percents_tot[c(2,4),]), main="small RNA-Seq Library percentages",
	xlab="Sample", col=c("blue","red"),
	legend = rownames(Percents_tot[c(2,4,5),]), beside=TRUE)
dev.new()	
barplot(as.matrix(Percents_lib[c(1,3),]), main="small RNA-Seq Library percentages",
	xlab="Sample", col=c("blue","red"),
	legend = rownames(Percents_lib[c(1,3),]), beside=TRUE)
dev.new()	
barplot(as.matrix(miR_q), main="small RNA-Seq Library percentages",
	xlab="Sample", col=c("black"),
	legend = rownames(as.matrix(miR_q)), beside=TRUE)
	
# Analysis of miRNAs:

miRs <- c('hsa-miR-1-5p','hsa-miR-1-3p','hsa-miR-155-5p','hsa-miR-155-3p')
ind <- which(!rownames(D_miR) %in% miRs)
D_q <- D_q[which(rownames(D_q)!='dme-miR-14-3p'),]
bg <- Remove_frac(D_miR[ind,],.5,1)
DATA_all <- rbind(D_q,D_miR[miRs,],bg)
DATA_miRs <- rbind(D_miR[miRs,],bg)

len_bg <- dim(bg)[1]
start <- 1/3
light <- 2.5*len_bg
cols_bg = sapply(seq(1,len_bg), function(x) {	
	return(rgb(light-x,light-x,light-x, maxColorValue = 3*len_bg))
}
)

cols_miRs <- c("#DF00DF","#FF00FF","#00FFFF","#00DFDF")
cols_q <- c("red","blue")
cols_all <- c(cols_q,cols_miRs,cols_bg)
cols_miRs <-c(cols_miRs,cols_bg)


dev.new(width = 8, height=6)
barplot(as.matrix(DATA_all), main="Reads",
	xlab="Sample", xlim = c(1,46), beside=FALSE, col = cols_all,
	legend.text = rownames(DATA_all),
	width = 3,
	space = 0.2,
	names.arg =  c("-","AGO1\n1","155","-","AGO2\n1","155","-","GFP\n1","155"))



dev.new(width = 8, height=6)
barplot(as.matrix(Col_frac(DATA_miRs)), main="% of mapped library",
	xlab="Sample", xlim = c(1,46), beside=FALSE, col = cols_miRs,
	legend.text = rownames(DATA_miRs),
	width = 3,
	space = 0.2,
	names.arg =  c("-","AGO1\n1","155","-","AGO2\n1","155","-","GFP\n1","155"))

dev.new(width = 8, height=6)
barplot(as.matrix(Col_div(DATA_all,as.matrix(T_q))), main="Relative abundance to spikes",
	xlab="Sample", xlim = c(1,46), beside=FALSE, col = cols_all,
	legend.text = rownames(DATA_all),
	width = 3,
	space = 0.2,
	names.arg =  c("-","AGO1\n1","155","-","AGO2\n1","155","-","GFP\n1","155"))

	
dev.new(width = 8, height=9)
barplot(as.matrix(Col_div(DATA_all,as.matrix(T_q))), main="Relative abundance to spikes",
	xlab="Sample", xlim = c(1,46), ylim = c(0,200), beside=FALSE, col = cols_all,
	legend.text = rownames(DATA_all),
	width = 3,
	space = 0.2,
	names.arg =  c("-","AGO1\n1","155","-","AGO2\n1","155","-","GFP\n1","155"))
	