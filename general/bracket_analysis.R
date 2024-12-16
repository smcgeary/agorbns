# Load Data
source("general/general.R")
library("gplots")
graphics.off()

data <- read.xls("bracket.xlsx")

data_bool <- as.matrix(data[, 4:ncol(data)])

competitions <- sprintf("%s: %s-%s", data[, 1], data[, 2], data[, 3])
rownames(data_bool) <- competitions
# Conditional asking if both 0 and 1 exists in a row.
rows_use <- unlist(lapply(apply(data_bool, 1, unique), function(ind) {
  0 %in% ind & 1 %in% ind
}))

totals <- apply(data_bool, 2, function(col) {
  sum(!(is.na(col)))  
})

data_bool <- data_bool[rows_use, ]
# data_bool <- data_bool[grep("Round 1", rownames(data_bool)), ]
# Standardize the data
# data_bool <- data_bool/apply(data_bool, 1, sd, na.rm=TRUE)
# data_bool <- data_bool - rowMeans(data_bool, na.rm=TRUE)




dev.new(xpos=20, ypos=20, height=4, width=4)
image(t(data_bool))

cor_matrix <- cor(data_bool, use="pairwise.complete.obs")

dev.new(xpos=20, ypos=420, height=4, width=4)
image(cor_matrix)

dev.new(xpos=420, ypos=20, height=7, width=7)

heatmap.2(data_bool, trace="none", col="topo.colors", lwid=c(0.1,1), lhei=c(0.1, 1),
          key=FALSE, margins=c(3, 20))
