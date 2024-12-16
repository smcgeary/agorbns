
temp_all <- read.table("temp_Jarrett_all.txt", row.names=NULL, header=FALSE)

dev.new(xpos=20, ypos=20, height=5, width=5)
plot(temp_all[, 1], temp_all[, 2], type="l")

trans_strings <- c("ENSG00000198727.2", "ENSG00000198938.2", "ENSG00000198695.2",
                   "ENSG00000198712.1", "ENSG00000198840.2", "ENSG00000198804.2",
                   "ENSG00000198886.2", "ENSG00000198786.2", "ENSG00000198888.2",
                   "ENSG00000198899.2", "ENSG00000212907.2", "ENSG00000198763.3",
                   "ENSG00000228253.1")
for (trans_string_i in trans_strings) {
  temp_all <- read.table(sprintf("temp_Jarrett_%s.txt", trans_string_i), row.names=NULL, header=FALSE)
  lines(temp_all[, 1], temp_all[, 2], type="l", col="gray")
}