
library(ggplot2)
library(ggseqlogo)

GetKmerEnrichments <- function(data) {
	data_norm <- t(t(data) / colSums(data))
	print(colSums(data_norm))
	data_R <- data_norm / data_norm[, 1]
	data_R[, -1]
}

graphics.off()

rbp <- "FUBP3"

condition_I <- "I"

file_path_input <- paste0("/lab/solexa_bartel/mcgeary/RBNS_DominguezFreeseAlexis/",
                          rbp, "/equilibrium/kmers/", condition_I,
                          "_5-mers.txt")
kmer_data_input <- read.table(file_path_input, sep="\t", row.names=1)
kmer_data_full <- kmer_data_input

print(condition_I)
for (i in seq(5)) {
	condition <- as.character(i)
	print(condition)
	file_path_rbp <- paste0("/lab/solexa_bartel/mcgeary/RBNS_DominguezFreeseAlexis/",
	                        rbp, "/equilibrium/kmers/", condition,
	                        "_5-mers.txt")
	kmer_data_rbp <- read.table(file_path_rbp, sep="\t", row.names=1)
	kmer_data_full <- cbind(kmer_data_full, kmer_data_rbp)
}

colnames(kmer_data_full) <- c("I", "1", "2", "3", "4", "5")

data_R <- GetKmerEnrichments(kmer_data_full)

data_R_average <- rowMeans(data_R)
ind_max <- which(data_R_average == max(data_R_average))

test_data_path <- paste0("/lab/solexa_bartel/mcgeary/RBNS_DominguezFreeseAlexis/",
                         "raw_data/ENCFF826BAQ.tsv")

test_data <- read.table(test_data_path, row.names=1, skip=1)

test_data_ordered <- test_data[rownames(data_R), ]

data_R_ranked <- data_R[rownames(test_data), ]

path_enrichment_file <- file.path("/lab/solexa_bartel/mcgeary",
                                   "RBNS_DominguezFreeseAlexis/raw_data",
                                   "Logo_motifs.txt")

enrichment_file <- readLines(path_enrichment_file)

logo_line <- unlist(strsplit(grep("FUBP3", enrichment_file, value=TRUE),
                             split="\t"))

motifs <- unlist(strsplit(grep("_", logo_line, value=TRUE), split="_"))

motif_matrix <- matrix(motifs, nrow=length(motifs)/3, ncol=3, byrow=TRUE)
print(motif_matrix)
rownames(motif_matrix) <- motif_matrix[, 1]
motif_matrix <- motif_matrix[, -1]

motif_matrix <- data.frame(motif_matrix, stringsAsFactors=FALSE)
colnames(motif_matrix) <- c("Motif #", "R-1")
print(motif_matrix)

motif_matrix <- motif_matrix[order(motif_matrix[, 2]), ]

path_motif_matrix_SM <- file.path("/lab/solexa_bartel/mcgeary",
                                   "RBNS_DominguezFreeseAlexis", rbp,
                                   "equilibrium/top_kmers/3.txt")

motif_matrix_SM <- read.table(path_motif_matrix_SM)

rownames(motif_matrix_SM) <- gsub("T", "U", rownames(motif_matrix_SM))


print(motif_matrix)
print(motif_matrix_SM)

all_motifs <- union(rownames(motif_matrix_SM), rownames(motif_matrix))

print(all_motifs)

matrix_compare <- matrix(0, nrow=length(all_motifs), ncol=2,
                         dimnames=list(all_motifs, c("SM", "PF")))

matrix_compare[rownames(motif_matrix_SM), 1] <- as.numeric(motif_matrix_SM[rownames(motif_matrix_SM), ])

matrix_compare[rownames(motif_matrix), 2] <- as.numeric(motif_matrix[rownames(motif_matrix), 2])


print(matrix_compare)

plot(matrix_compare[, 1], matrix_compare[, 2])

abline(0, 1, lty=2)
break

print(motif_matrix)

print(logo_line)

print(motifs)



break








pwm_path <- "/lab/solexa_bartel/mcgeary/RBNS_DominguezFreeseAlexis/FUBP3/equilibrium/pfm/3_1.txt"


pwm <- read.table(pwm_path, sep="\t", header=FALSE, row.names=1)

row.names(pwm)[4] <- "U"
print(pwm)

dev.new(xpos=20, ypos=20, height=2.8, width=2.8)

p1 <- ggseqlogo(as.matrix(pwm))

print(p1)


pwm_path <- "/lab/solexa_bartel/mcgeary/RBNS_DominguezFreese/FUBP3/equilibrium/pfm/3_2.txt"


pwm <- read.table(pwm_path, sep="\t", header=FALSE, row.names=1)

row.names(pwm)[4] <- "U"
print(pwm)

dev.new(xpos=520, ypos=20, height=2.8, width=2.8)

p1 <- ggseqlogo(as.matrix(pwm))

print(p1)


pwm_path <- "/lab/solexa_bartel/mcgeary/RBNS_DominguezFreese/FUBP3/equilibrium/pfm/3_3.txt"


pwm <- read.table(pwm_path, sep="\t", header=FALSE, row.names=1)

row.names(pwm)[4] <- "U"
print(pwm)

dev.new(xpos=1020, ypos=20, height=2.8, width=2.8)

p1 <- ggseqlogo(as.matrix(pwm))

print(p1)










dev.new(xpos=20, ypos=20, height=5, width=5)

plot(data_R[rownames(data_R), 1], test_data[rownames(data_R), 1])







