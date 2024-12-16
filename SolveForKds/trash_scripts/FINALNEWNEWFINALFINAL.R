 
library(data.table)
library(colorspace)
source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")



# data_new <- read.table(paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
#                      experiment, "/kds_PAPER/", start, "-", stop, "_", 
#                      sitelist, "_singlebg_multinomialbootstrap_PAPER.txt"),header=TRUE)
# data_new <- data_new[, -ncol(data_new)]

# plot(colMeans(data_new),colMeans(data_new))
# points(colMeans(data_new),apply(data_new,2, function(column) {
#   sort(column)[5]
# }))
# points(colMeans(data_new),apply(data_new,2, function(column) {
#   sort(column)[95]
# }))
# break
# sitesXcounts <- GetSitesXCounts(mirna,
#                                 experiment,
#                                 start,
#                                 stop,
#                                 sitelist,
#                                 mirna.start = mirna.start,
#                                 mirna.stop = mirna.stop)
# # Separate site sequences from data file.
# if (sitelist %in% kmer_list_names) {
#   seqs <- rownames(sitesXcounts)
# } else {
#   seqs <- sitesXcounts[,1]
#   sitesXcounts <- sitesXcounts[,-1]
# }

# # Remove any rows with no reads whoatsoever.
# sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
# # Remove any rows for which there are no input reads:
# sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]

# data <- sitesXcounts[,2:6]
# data <- data.matrix(data)
# par(mfrow=c(1, 2))
# PlotMeanVersusSD <- function(vector, bins) {
#   bins_vector <- seq(0, max(vector), length = bins + 1)
#   lefts <- bins_vector[1 : bins]
#   rights <- bins_vector[2 : (bins + 1)]
#   print(rbind (lefts, rights))
#   XandY <- apply(cbind(lefts, rights)[1:20,], 1, function(row) {
#     # print(row)
#     data <- vector[vector > row[1] & vector <= row[2]]
#     print(data)
#     return(c(mean(data), mean((data - mean(data))^2), length(data), row[1], row[2]))
#   })
#   print(XandY)
#   plot(XandY[1, ], XandY[2, ],type="o")
#   plot(XandY[1, ], XandY[3, ], type="o")

# }


# PlotMeanVersusSD(sitesXcounts[-nrow(sitesXcounts),1], 70)
# break
# sitesXcounts <- GetSitesXCounts(mirna,
#                                 experiment,
#                                 start,
#                                 stop,
#                                 sitelist,
#                                 mirna.start = mirna.start,
#                                 mirna.stop = mirna.stop)
# # Separate site sequences from data file.
# if (sitelist %in% kmer_list_names) {
#   seqs <- rownames(sitesXcounts)
# } else {
#   seqs <- sitesXcounts[,1]
#   sitesXcounts <- sitesXcounts[,-1]
# }





# # Remove any rows with no reads whoatsoever.
# sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
# # Remove any rows for which there are no input reads:
# sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]
# # plot(rowMeans(sitesXcounts[-nrow(sitesXcounts),c(1, 7)]),
# #   ((sitesXcounts[,1] - rowMeans(sitesXcounts[, c(1, 7)]))^2 +
# #    (sitesXcounts[,7] - rowMeans(sitesXcounts[, c(1, 7)]))^2)[-nrow(sitesXcounts)]/2)
# # break
# data <- sitesXcounts[,2:6]
# data <- data.matrix(data)

# rownames(data)[nrow(data)] <- "None"
# # data <- rbind(data_new,data[nrow(data),,drop=FALSE])

# num.kds <- nrow(data)-1
# num.bgs <- 1

# # Define the vector of total site type concentrations:
# c.I.tots <- Norm(sitesXcounts[, 1]) * k.c.lib
# c.totals <- matrix(c.I.tots, nrow=length(c.I.tots), ncol=5, byrow=FALSE)

# colnames(c.totals) <- colnames(data)
# rownames(c.totals) <- rownames(data)

# # Define the vector of AGO concentrations:
# c.agos <- sapply(colnames(data), function(x) {
#   as.numeric(x) / 100
# })

par(mfrow=c(5, 4))
# ModelFunction <- function(pars) {
#   # Split up the parameters into the kd and background parameters.
#   kds  <- c(Logistic(pars[1 : num.kds], 10),1)
#   names(kds) <- rownames(data)
#   bgs  <- rep(10^pars[(num.kds + 1)], 5)
#   stock.ago <- 10^pars[num.kds + 2]
#   # Get the bound counts in each experiment:
#   c.bounds <- as.matrix(
#     sapply(c.agos, function(ago.percent) {
#       return(GetBoundRNA(kds, c.I.tots, ago.percent * stock.ago))
#     }
#   ))

#   # Get the amount of background binding by subtracting the bound from the
#   # total sites in each experiment, normalizing. Must transpose to multiply
#   # each column.
#   c.frees <- c.totals - c.bounds
#   c.bgs <- t(t(c.frees) * bgs / colSums(c.frees))
#   c.all <- c.bounds + c.bgs
#   c.final <- data.frame(t(t(c.all) / colSums(c.all) * colSums(data)))
#   print(dim(c.final))
#   plot(unlist(c.final),c(data), col = rep(c(rgb(0,0,0,alpha =0.4),rgb(1,0,0,alpha=0.4),rgb(0,1,0,alpha=0.4),
#                                       rgb(1,1,0,alpha=0.4), rgb(1,0,1, alpha = 0.4)),each=length(unlist(c.final))/5),log='xy')

#   return(c.final)
# }
for (mirna in c("let-7a", "miR-1", "miR-155", "miR-124", "lsy-6")) {
  for (mirna.start in seq(2, 5)) {

  in_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                       experiment, "/kds_PAPER/", start, "-", stop, "_", 
                       sitelist, "_", mirna.start, "-", mirna.start+3,
                       "_singlebg_poisson_PAPER_last.txt")

pars_raw <- fread(in_file, sep = "\t", header=TRUE, colClasses = "numeric",stringsAsFactors = FALSE)
pars_table <- matrix(unlist(pars_raw), nrow=nrow(pars_raw), ncol = ncol(pars_raw), byrow=FALSE)
pars_logresidual <- unlist(pars_table[nrow(pars_table),])
names(pars_logresidual) <- colnames(pars_raw)

# c.final <- ModelFunction(pars_logresidual)

  in_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                       experiment, "/kds_PAPER/", start, "-", stop, "_", 
                       sitelist, "_", mirna.start, "-", mirna.start+3,
                       "_singlebg_multinomial_PAPER_last.txt")
pars_raw <- fread(in_file, sep = "\t", header=TRUE, colClasses = "numeric",stringsAsFactors = FALSE)
pars_table <- matrix(unlist(pars_raw), nrow=nrow(pars_raw), ncol = ncol(pars_raw), byrow=FALSE)
pars_multi <- unlist(pars_table[nrow(pars_table),])
names(pars_multi) <- colnames(pars_raw)
# ModelFunction(pars_multi)
colors_all <- c(rainbow_hcl(n=length(pars_multi)-3), "black","black")
plot(pars_logresidual[-length(pars_logresidual)],pars_multi[-length(pars_multi)],xlim=c(-10, 4), ylim = c(-10, 4),col=colors_all)
segments(-10, -10, 4, 4, lty = 2)
}}

# plot(rowSums(data)[1:num.kds],pars_logresidual[1:num.kds] - pars_multi[1:num.kds],log='x',col=colors_all,ylim = c(-4, 4))