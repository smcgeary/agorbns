source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")
print("Out of general")
source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/ModelingFunctions.R")
print("out of modeling functions")

args        <- commandArgs(trailingOnly=TRUE)
mirna       <- args[1]
mirna.start <- as.integer(args[2])
if (mirna == "miR-1") {
  data_i <- 1
  str.nocombI <- "nocombInput_"
} else {
  data_i <- 2
  str.nocombI <- ""
}
if (mirna == "miR-7") {
  experiment <- "equilibrium2_nb"
} else {
  experiment <- "equilibrium"
}
n_constant <- 5
pars_temp_path <- sprintf("16merfits/%s_%s_%s_%s/temp_global.txt", mirna,
                              experiment, n_constant, str.nocombI)
print(pars_temp_path)
n_x <- 4^8 + 1
ind_splits <- c(0, 1)
names(ind_splits) <- c("left", "right")
    # This accounts for the fact tha the order is 1_left, 1_right, 2_left, etc.
kd_l_start_ind  <-  2*(mirna.start - 1)     *n_x
kd_r_start_ind <- (2*(mirna.start - 1) + 1)*n_x
kds_vec_l  <- NaN
kds_vec_r <- NaN
while(class(kds_vec_l) != "data.frame") {
  kds_vec_l <- data.frame(fread(pars_temp_path, skip=kd_l_start_ind,
                                   nrow=n_x-1, header=FALSE), row.names=1)
}
while(class(kds_vec_r) != "data.frame") {
  kds_vec_r <- data.frame(fread(pars_temp_path, skip=kd_r_start_ind,
                                    nrow=n_x-1, header=FALSE), row.names=1)
}
kmers_l <- substr(rownames(kds_vec_l), 1, 12)
kmers_r <- substr(rownames(kds_vec_r), 1, 12)
base_kmers_l <- sort(unique(substr(kmers_l, 5, 12)))
base_kmers_r <- sort(unique(substr(kmers_r, 1, 8)))
base_kmers <- base_kmers_l
kds_16mers <- c()
tick <- 0
output.matrix <- matrix(NaN, nrow=256, ncol=514)
rownames(output.matrix) <- base_kmers
print(kmers_r[1:256])
kmers_4nt <- substr(kmers_r[1:256], 9, 12)
print(kmers_4nt)
print(paste0(kmers_4nt, "_left"))
print(paste0(kmers_4nt, "_right"))
colnames(output.matrix) <- c("mean_l", "mean_r",
                             paste0(kmers_4nt, "_l"),
                             paste0(kmers_4nt, "_r"))

# colnames(output.matrix) <- c("mean_left", "mean_right", paste0())
print(head(output.matrix))
sapply(1:length(base_kmers), function(i_k) {
  base_kmer <- base_kmers[i_k]
  kds_l <- kds_vec_l[grep(sprintf("%s$", base_kmer), kmers_l), , drop=FALSE]
  kds_r <- kds_vec_r[grep(sprintf("^%s", base_kmer), kmers_r), , drop=FALSE]
  temp_names_l <- rownames(kds_l)
  temp_names_r <- rownames(kds_r)
  kds_l <- kds_l[, 1]
  kds_r <- kds_r[, 1]
  names(kds_l) <- paste0(substr(temp_names_l, 1, 4), "_l")
  names(kds_r) <- paste0(substr(temp_names_r, 9, 12), "_r")
  # kmers_16 <- c(sapply(names(kds_l), function(kmer_l) {
  #   c(sapply(names(kds_r), function(kmer_r) {
  #     paste0(kmer_l, base_kmer, kmer_r)
  #   }))
  # }))
  del_kds_l <- kds_l - mean(kds_l)
  del_kds_r <- kds_r - mean(kds_r)
  # print(c("Mean_left"))
  # checknames <- cbind(colnames(output.matrix),
  #                     c("mean_l", "mean_r", names(del_kds_l), names(del_kds_r)))
  # apply(checknames, 1, function(row) {
  #   if (row[1] != row[2]) {
  #     print("something is not right!")
  #     break
  #   }
  # })
  output.matrix[i_k, ] <<- c(mean(kds_l), mean(kds_r),
                                       del_kds_l, del_kds_r)
  # kds_16 <- log(exp(kds_l)%o%exp(kds_r)) - mean(c(kds_l, kds_r))
  # kds_16_vec <- c(t(kds_16))
  # names(kds_16_vec) <- kmers_16
  # kds_16mers <<- c(kds_16mers, kds_16_vec)
  tick <<- tick + 1
  if (tick%%16 == 0) {
    print(tick/256)
  }
})

print(head(output.matrix))
output.matrix <- cbind(rownames(output.matrix), output.matrix)
colnames(output.matrix)[1] <- "base_kmer"

outputDir <- sprintf("%s%s/%s/16mer_kds", kSolexaDir, mirna,
                         experiment)
print(outputDir)
if (!file.exists(outputDir)) {
dir.create(outputDir)
}
file_name <- sprintf("%s/nt%s-%s_table.txt", outputDir, mirna.start,
                         as.integer(mirna.start) + 3)

print(file_name)

write.table(file=file_name, output.matrix,
            row.names=FALSE, sep="\t", quote=FALSE, col.names=TRUE)
