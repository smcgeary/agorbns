## sRS_analyze_read_data


D <- read.table("140603_WIGTC-HISEQA_C4W8TACXX_expression_matrix.txt", header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors=FALSE)
print(colnames(D))

# Ordering table according to chromosome & position (without consideration of strandedness)
rownames(D)<- D$REF

D <- D[,-1]
D_tags <- D[1,]
D <- D[-1,]
D <- D[, order(as.integer((colnames(D))))]

print(head(D))

D_size = D[1:2,]
D_q = D[3:5,]
D_capt = D[6:9,]
D_miR = D[10:nrow(D)-6,]
D_tots = D[nrow(D)-6:nrow(D),]

filter_zeros <-function (mat) {
	keep <- which(rowSums(mat)>0)
	print(keep)
	output <- mat[keep,]
	output
	}

D_miR <- filter_zeros(D_miR)

# Function to column-normalize expression data, in reads-per-million.
Col_norm <- function(X) t(apply(X,1,function(x) x/colSums(X)))*10^6

A_tot = matrix(colSums(A_exp),nrow = nrow(A_exp), ncol = ncol(A_exp),byrow=TRUE)
A_norm = t(t(A_exp)/rowSums(t(A_exp)))*10^6

# miRNA info
A_inf = A[,1:6]
rownames(A_inf) = A$miRNA.Sample

# FARAZI + ENCODE REMAPPED	
B = read.table("/Volumes/bartel_lab1/mcgeary/computation/miR_coverage_vikram/expression_clustering/FARAZIreanalysis+ENCODE_output.gtf", header = TRUE, sep = "\t", comment.char = "")
# Raw table	
B = B[order(A$chr, A$start),]
B_exp = B[,7:ncol(B)]
rownames(B_exp) = B$miRNA.Sample
colnames(B_exp) = EXP_$V1
B_norm = t(t(B_exp)/rowSums(t(B_exp)))*10^6

B_n_enc = B_norm[,246:ncol(B_norm)]
B_n_far = B_norm[,1:245]

rownames(B_n_enc) = B$miRNA.Sample
rownames(B_n_far) = B$miRNA.Sample


B_n_noreps = matrix(0,nrow = nrow(B_n_enc),ncol = length(UNIQUE$V1))
NEWnames = as.vector(sapply(colnames(B_n_enc),function (x) strsplit(x,"Rep")[[1]][1]))
EXP_ = read.table("/Volumes/bartel_lab1/mcgeary/computation/miR_coverage_vikram/expression_clustering/FARAZI+ENCODE_exp_names.txt", header = FALSE, sep = "\t", comment.char = "")
UNIQUE = read.table("/Volumes/bartel_lab1/mcgeary/computation/miR_coverage_vikram/expression_clustering/ENCODE_exp_names_unique.txt", header = FALSE, sep = "\t", comment.char = "")

for (i in 1:length(UNIQUE$V1)) {

	temp = B_n_enc[,NEWnames==UNIQUE$V1[i]]

	if (class(temp)=="matrix"){
		temp = as.matrix(rowMeans(temp),ncol=1)
		}
		B_n_noreps[,i] = temp

		}
rownames(B_n_noreps) = B$miRNA.Sample
colnames(B_n_noreps) = UNIQUE$V1


B_inf = B[,1:6]

A_ave = rowMeans(A_norm)
B_ave = rowMeans(B_norm)
B_ave_far = rowMeans(B_n_far)
B_ave_enc = rowMeans(B_n_enc)


# Determines all unique chromosomal occurences of data
CHROM <- function(inf) {
chr<- unique(inf$chr., use.names = TRUE)
str <- unique(inf$strand, use.names = TRUE)
output <- na.omit(expand.grid(chr,str))
colnames(output) <- c("chr.","strand")
output
}

# Finds all unique miRNA pairs per chromosome / strand
SELECT <- function(X,Y) {
	strand <- rownames(na.omit(X[X$chr.==Y$chr&X$strand==Y$strand,]))
	print(strand)
	n <- length(strand)
	print(n)
	print((n*(n-1))/2)
	pairs <- expand.grid(strand,strand)
	pairs <- pairs[pairs[,1]!=pairs[,2],]
	print(dim(unique(t(apply(pairs,1,sort)))))		
	unique(t(apply(pairs,1,sort)))

}

# Finds all unique miRNA pairs connected by a physical chromosome
ALL <- function(X,Y) {
	output <- matrix(data=NA,nrow=0,ncol=2)
	print(nrow(Y))
	for (i in 1:nrow(Y)) {
		print(Y[i,])
		new <- SELECT(X,Y[i,])
		print('nrow')
		print(nrow(new))
		if (nrow(new)>0) {
		output <- rbind(output,new)
		}
	}
	output
}

PAIRS <- ALL(A_inf,CHROM)

PAIR_CORR <- function(D,X) {
	temp <- cor(t(D), method = "pearson")
	 temp[rownames(temp)==X[1],colnames(temp)==X[2]]
	}

ALL_CORR <- function(D,PAIRS) {
	D <- D
	apply(PAIRS,1,function(PAIRS) PAIR_CORR(D,PAIRS))
}

ALL_DIST <- function(PAIRS) {
	temp <- as.matrix(apply(PAIRS,1:2,function(X) (A_inf$start[rownames(A_inf)==X]+A_inf$stop[rownames(A_inf)==X])/2),ncol=2)
	abs(temp[,1]-temp[,2])
	
	}



# Plots showing correlation of Vikram's analysis with original farazi analysis

dev.new(width=7, height=7)

plot(A_ave,B_ave, type = 'p', log = "xy")
dev.new(width=7, height=7)
plot(A_ave,B_ave_enc, type = 'p',log = "xy")
dev.new(width=7, height=7)
plot(A_ave,B_ave_far, type = 'p', log = "xy")

Dom_MEAN = NULL
Dom_STD = NULL
tot_MEAN = NULL
tot_STD = NULL




filter <-function (exp,inf) {cut <- which(rowSums(exp)>cutoff/ncol(exp));
	
	
	inf = inf[cut,]
	exp = exp[cut,]

	
	Keep = matrix(1,nrow(inf),1)
rownames(Keep) = inf$miRNA.Sample

dup1 = sort(nrow(inf)-which(duplicated(rev(inf$start))))+1
dup2 = which(duplicated(inf$start))
DUP = cbind(dup1,dup2)

DUPdif = exp[DUP[,1],]-exp[DUP[,2],]
DUPsum = exp[DUP[,1],]+exp[DUP[,2],]
DOM = DUPdif/DUPsum
DOM_mean = rowMeans(DOM, na.rm = TRUE)
DOM_sd = apply(DOM,1,sd,na.rm = TRUE)
SUM_mean = rowMeans(DUPsum, na.rm = TRUE)
SUM_sd = apply(DUPsum,1,sd,na.rm = TRUE)
stats = cbind(DOM_mean,DOM_sd,SUM_mean,SUM_sd)
colnames(stats) = c('DOM mean', 'DOM SD', 'TOTAL mean', 'TOTAL SD')

DOM_round = round(1.5-DOM_mean/2)
USE = sort(c(c(1:nrow(inf))[-DUP],DUP[which(DOM_round==2)]+1,DUP[which(DOM_round==1)]))

output <- list(exp[USE,],inf[USE,])
	}

Bdom_far = filter(B_n_far,B_inf)
Bdom_enc = filter(B_n_enc,B_inf)
Bdom = filter(B_norm,B_inf)
Bdom_enc_noreps = filter(B_n_noreps,B_inf)






distpear <- function (X) {as.dist(1-cor(t(X), method = "pearson"))}

lognorm <- function (X) {normalize(log(X+1,10))}

EXP_ <- function (X) {X[1][[1]]}

INF_ <- function (X) {X[2][[1]]}

MAPS <- function (X,Y) {
	X <- EXP_(X)
	X <- lognorm(X)
	dev.new(1,1)
	heatmap(X,distfun = distpear,col = greenred(11), main = Y,breaks=c(min(X),seq(-3,3,length=10),max(X)))
}



MAPS(Bdom_enc,"ENCODE")
MAPS(Bdom,"ENCODE & Farazi")
MAPS(Bdom_far,"Farazi")
MAPS(Bdom_enc_noreps,"ENCODE, averaged replicates")

CHROM <- function(X,chr,str,title) {
	dev.new()
	coord <- INF_(X)
	X <- EXP_(X)
	X <- lognorm(X)
	i <- which(coord$chr.==chr & coord$strand==str)
	INF <- coord[i,]
	EXP <- X[i,]
	heatmap(EXP,distfun = distpear,col = greenred(11), main = title, seq(-3,3,length=10),max(X))
}

PRCOMP <- function(X,n) {
	X <- EXP_(X)
	X <- lognorm(X)
	pca <- prcomp(X)
	EXP <- pca$x[,1:n]
	# dev.new(3,1)
	pdf(file="a.pdf", height=20, width=10)
	heatmap(EXP,Colv = NA, scale = "none", distfun = distpear,col = greenred(11), main = "PCA clustering",breaks=c(min(X),seq(-3,3,length=10),max(X)))
	dev.off()
}

PRCOMP(Bdom,10)

