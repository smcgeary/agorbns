     # 1. Get data table:
mirnas <- c("let-7a", "miR-1", "miR-155", "miR-124", "lsy-6")


mirnas <- c("miR-1")

library(multicore)

library(data.table)
par(mfcol=c(4,5))
GetData <- function(mirna, mirna.start, mirna.stop){
sitesXcounts <- GetSitesXCounts(mirna,
                                        experiment,
                                        start,
                                        stop,
                                        sitelist,
                                        mirna.start = mirna.start,
                                        mirna.stop = mirna.stop)
# Separate site sequences from data file.
if (sitelist == "12mers") {
  seqs <- rownames(sitesXcounts)
} else {
  seqs <- sitesXcounts[,1]
  sitesXcounts <- sitesXcounts[,-1]
}
# if (sitelist == "12mers") {
#   sitesXcounts_new <- t(sapply(1:((nrow(sitesXcounts)-1)/256),function(x){colSums(sitesXcounts[1:256+(x-1)*256,])}))
#   sitesXcounts_new <- rbind(sitesXcounts_new,sitesXcounts[nrow(sitesXcounts),,drop=FALSE])
#   sitesXcounts <- sitesXcounts_new
# }


print("23")
# Remove any rows with no reads whoatsoever.
sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
# Remove any rows for which there are no input reads:
sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]

data <- sitesXcounts[,2:6]
data <- data.matrix(data)
return(data)
}
print(dim(data))

rownames(data)[nrow(data)] <- "None"

num.kds <- nrow(data)
num.bgs <- 1
print(46)

# Define the vector of total site type concentrations:
c.I.tots <- Norm(sitesXcounts[, 1]) * k.c.lib
colnames(c.totals) <- colnames(data)
rownames(c.totals) <- rownames(data)
names(c.I.tots) <- rownames(data)
# Define the vector of AGO concentrations:
c.agos <- sapply(colnames(data), function(x) {
   as.numeric(x) / 100
  })
print(c.agos)




GetPars <- function(mirna,mirna.start, mirna.stop){
in_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                   experiment, "/kds_PAPER/", start, "-", stop, "_", 
                   sitelist, "_", mirna.start, "-", mirna.stop,
                   "_singlebg_combinedinput_PAPER_final_last.txt")




pars_table <- try(fread(in_file, sep = "\t", header=TRUE, colClasses = "numeric",stringsAsFactors = FALSE))

while(class(pars_table) == "try-error") {
  pars_table <- try(fread(in_file, sep = "\t", header=TRUE, colClasses = "numeric",stringsAsFactors = FALSE))

}
table_names <- names(pars_table[1])
pars_table_new <- matrix(unlist(pars_table), nrow=nrow(pars_table), ncol = ncol(pars_table), byrow=FALSE)
pars <- unlist(pars_table_new[nrow(pars_table_new),])
names(pars) <- table_names
pars_old <- pars
pars <- pars_old[-length(pars_old)]

print(length(pars))
print(sort(pars)[1:20])
return(pars)
}
 
GetModel <- function(mirna, pars){
    kds  <- c(Logistic(pars[1 : (num.kds-1)], 10),1)
    kds_first <- kds
    bgs  <- rep(10^pars[(num.kds)], 5)
  B <- bgs[1]
    stock.ago <- 10^pars[num.kds + 1]
  names(kds) <- rownames(data)
  print("hi")
  f.mat <- matrix(
      rep(sapply(c.agos, function(ago.percent) {
           return(GetFreeAgo(kds, c.I.tots, ago.percent * stock.ago))
         }
         ), nrow(data)), nrow=nrow(data), ncol=ncol(data), byrow=TRUE)


  A.mat <- matrix(c.agos*stock.ago, nrow=nrow(data), ncol = ncol(data), byrow=TRUE)




  K.mat <- matrix(kds, nrow=nrow(data), ncol = ncol(data), byrow=FALSE)


  l.mat <- matrix(c.I.tots, nrow=nrow(data), ncol = ncol(data), byrow=FALSE)


  R.mat <- matrix(colSums(data), nrow = nrow(data), ncol = ncol(data), byrow=TRUE)



  L <- 100

  c.final <- (-R.mat*l.mat * (f.mat^2 + (L - A.mat)*f.mat + K.mat*B) *
                  (f.mat^3 - (2*A.mat + B - L - K.mat)*f.mat^2 +
                    (A.mat^2 + (B - L)*A.mat - L*B - 
                     (2*A.mat + B - L)*K.mat)*f.mat +
                    (A.mat^2 + (B - L)*A.mat - L*B)*K.mat)^-1)

  colnames(c.final) <- colnames(data)
  plot(c(c.final),c(data), col = rep(c(rgb(0,0,0,alpha =0.4),rgb(1,0,0,alpha=0.4),rgb(0,1,0,alpha=0.4),
                                      rgb(1,1,0,alpha=0.4), rgb(1,0,1, alpha = 0.4)),each=length(c(c.final))/5),log='xy')
  title(main = mirna)

  sumoflogsquares <<- sum((c.final - data)^2)
  print(sumoflogsquares)
}


sitesXcounts <- GetSitesXCounts(mirna,
                                        experiment,
                                        start,
                                        stop,
                                        sitelist,
                                        mirna.start = 3,
                                        mirna.stop = 6)
# Separate site sequences from data file.
if (sitelist == "12mers") {
  seqs <- rownames(sitesXcounts)
} else {
  seqs <- sitesXcounts[,1]
  sitesXcounts <- sitesXcounts[,-1]
}
# if (sitelist == "12mers") {
#   sitesXcounts_new <- t(sapply(1:((nrow(sitesXcounts)-1)/256),function(x){colSums(sitesXcounts[1:256+(x-1)*256,])}))
#   sitesXcounts_new <- rbind(sitesXcounts_new,sitesXcounts[nrow(sitesXcounts),,drop=FALSE])
#   sitesXcounts <- sitesXcounts_new
# }


print("23")
# Remove any rows with no reads whoatsoever.
sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
# Remove any rows for which there are no input reads:
sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]

data <- sitesXcounts[,2:6]
data <- data.matrix(data)

print(dim(data))

rownames(data)[nrow(data)] <- "None"

num.kds <- nrow(data)
num.bgs <- 1
print(46)

# Define the vector of total site type concentrations:
c.I.tots <- Norm(sitesXcounts[, 1]) * k.c.lib
c.totals <- matrix(c.I.tots, nrow=length(c.I.tots), ncol=5, byrow=FALSE)
colnames(c.totals) <- colnames(data)
rownames(c.totals) <- rownames(data)
names(c.I.tots) <- rownames(data)
# Define the vector of AGO concentrations:
c.agos <- sapply(colnames(data), function(x) {
   as.numeric(x) / 100
  })
print(c.agos)




# plot(ecdf(log(c.I.tots)))
print(sort(c.I.tots,decreasing=TRUE)[1:10])
print(sort(c.I.tots)[1:10])

in_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                   experiment, "/kds_PAPER/", start, "-", stop, "_", 
                   sitelist, "_", 3, "-", 6,
                   "_singlebg_combinedinput_PAPER_final_last.txt")




pars_table <- try(fread(in_file, sep = "\t", header=TRUE, colClasses = "numeric",stringsAsFactors = FALSE))

while(class(pars_table) == "try-error") {
  pars_table <- try(fread(in_file, sep = "\t", header=TRUE, colClasses = "numeric",stringsAsFactors = FALSE))

}
table_names <- names(pars_table[1])
pars_table_new <- matrix(unlist(pars_table), nrow=nrow(pars_table), ncol = ncol(pars_table), byrow=FALSE)
pars <- unlist(pars_table_new[nrow(pars_table_new),])
names(pars) <- table_names
pars_old <- pars
pars <- pars_old[-length(pars_old)]

print(length(pars))
print(sort(pars)[1:20])

  time_prior <- proc.time()

  # Split up the parameters into the kd and background parameters.
    kds  <- c(Logistic(pars[1 : (num.kds-1)], 10),1)
    kds_first <- kds
    bgs  <- rep(10^pars[(num.kds)], 5)
  B <- bgs[1]
    stock.ago <- 10^pars[num.kds + 1]
  names(kds) <- rownames(data)
  print("hi")
  f.mat <- matrix(
      rep(sapply(c.agos, function(ago.percent) {
           return(GetFreeAgo(kds, c.I.tots, ago.percent * stock.ago))
         }
         ), nrow(data)), nrow=nrow(data), ncol=ncol(data), byrow=TRUE)


  A.mat <- matrix(c.agos*stock.ago, nrow=nrow(data), ncol = ncol(data), byrow=TRUE)




  K.mat <- matrix(kds, nrow=nrow(data), ncol = ncol(data), byrow=FALSE)


  l.mat <- matrix(c.I.tots, nrow=nrow(data), ncol = ncol(data), byrow=FALSE)


  R.mat <- matrix(colSums(data), nrow = nrow(data), ncol = ncol(data), byrow=TRUE)



  L <- 100

  c.final <- (-R.mat*l.mat * (f.mat^2 + (L - A.mat)*f.mat + K.mat*B) *
                  (f.mat^3 - (2*A.mat + B - L - K.mat)*f.mat^2 +
                    (A.mat^2 + (B - L)*A.mat - L*B - 
                     (2*A.mat + B - L)*K.mat)*f.mat +
                    (A.mat^2 + (B - L)*A.mat - L*B)*K.mat)^-1)

  colnames(c.final) <- colnames(data)
  plot(c(c.final),c(data), col = rep(c(rgb(0,0,0,alpha =0.4),rgb(1,0,0,alpha=0.4),rgb(0,1,0,alpha=0.4),
                                      rgb(1,1,0,alpha=0.4), rgb(1,0,1, alpha = 0.4)),each=length(c(c.final))/5),log='xy')
  title(main = mirna)

  sumoflogsquares <<- sum((c.final - data)^2)
  print(sumoflogsquares)


















sitesXcounts <- GetSitesXCounts(mirna,
                                        experiment,
                                        start,
                                        stop,
                                        sitelist,
                                        mirna.start = 2,
                                        mirna.stop = 5)
# Separate site sequences from data file.
if (sitelist == "12mers") {
  seqs <- rownames(sitesXcounts)
} else {
  seqs <- sitesXcounts[,1]
  sitesXcounts <- sitesXcounts[,-1]
}
# if (sitelist == "12mers") {
#   sitesXcounts_new <- t(sapply(1:((nrow(sitesXcounts)-1)/256),function(x){colSums(sitesXcounts[1:256+(x-1)*256,])}))
#   sitesXcounts_new <- rbind(sitesXcounts_new,sitesXcounts[nrow(sitesXcounts),,drop=FALSE])
#   sitesXcounts <- sitesXcounts_new
# }


print("23")
# Remove any rows with no reads whoatsoever.
sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
# Remove any rows for which there are no input reads:
sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]

data <- sitesXcounts[,2:6]
data <- data.matrix(data)

print(dim(data))

rownames(data)[nrow(data)] <- "None"

num.kds <- nrow(data)
num.bgs <- 1
print(46)

# Define the vector of total site type concentrations:
c.I.tots <- Norm(sitesXcounts[, 1]) * k.c.lib
c.totals <- matrix(c.I.tots, nrow=length(c.I.tots), ncol=5, byrow=FALSE)
colnames(c.totals) <- colnames(data)
rownames(c.totals) <- rownames(data)
names(c.I.tots) <- rownames(data)
# Define the vector of AGO concentrations:
c.agos <- sapply(colnames(data), function(x) {
   as.numeric(x) / 100
  })
print(c.agos)




# plot(ecdf(log(c.I.tots)))
print(sort(c.I.tots,decreasing=TRUE)[1:10])
print(sort(c.I.tots)[1:10])

in_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                   experiment, "/kds_PAPER/", start, "-", stop, "_", 
                   sitelist, "_", 2, "-", 5,
                   "_singlebg_logresiduals_combinedinput_PAPER_final_last.txt")




pars_table <- try(fread(in_file, sep = "\t", header=TRUE, colClasses = "numeric",stringsAsFactors = FALSE))

while(class(pars_table) == "try-error") {
  pars_table <- try(fread(in_file, sep = "\t", header=TRUE, colClasses = "numeric",stringsAsFactors = FALSE))

}
table_names <- names(pars_table[1])
pars_table_new <- matrix(unlist(pars_table), nrow=nrow(pars_table), ncol = ncol(pars_table), byrow=FALSE)
pars <- unlist(pars_table_new[nrow(pars_table_new),])
names(pars) <- table_names
pars_old <- pars
pars <- pars_old[-length(pars_old)]

print(length(pars))
print(sort(pars)[1:20])

  time_prior <- proc.time()

  # Split up the parameters into the kd and background parameters.
    kds  <- c(Logistic(pars[1 : (num.kds-1)], 10),1)
    kds_first <- kds
    bgs  <- rep(10^pars[(num.kds)], 5)
  B <- bgs[1]
    stock.ago <- 10^pars[num.kds + 1]
  names(kds) <- rownames(data)
  print("hi")
  f.mat <- matrix(
      rep(sapply(c.agos, function(ago.percent) {
           return(GetFreeAgo(kds, c.I.tots, ago.percent * stock.ago))
         }
         ), nrow(data)), nrow=nrow(data), ncol=ncol(data), byrow=TRUE)


  A.mat <- matrix(c.agos*stock.ago, nrow=nrow(data), ncol = ncol(data), byrow=TRUE)




  K.mat <- matrix(kds, nrow=nrow(data), ncol = ncol(data), byrow=FALSE)


  l.mat <- matrix(c.I.tots, nrow=nrow(data), ncol = ncol(data), byrow=FALSE)


  R.mat <- matrix(colSums(data), nrow = nrow(data), ncol = ncol(data), byrow=TRUE)



  L <- 100

  c.final <- (-R.mat*l.mat * (f.mat^2 + (L - A.mat)*f.mat + K.mat*B) *
                  (f.mat^3 - (2*A.mat + B - L - K.mat)*f.mat^2 +
                    (A.mat^2 + (B - L)*A.mat - L*B - 
                     (2*A.mat + B - L)*K.mat)*f.mat +
                    (A.mat^2 + (B - L)*A.mat - L*B)*K.mat)^-1)

  colnames(c.final) <- colnames(data)
  plot(c(c.final),c(data), col = rep(c(rgb(0,0,0,alpha =0.4),rgb(1,0,0,alpha=0.4),rgb(0,1,0,alpha=0.4),
                                      rgb(1,1,0,alpha=0.4), rgb(1,0,1, alpha = 0.4)),each=length(c(c.final))/5),log='xy')
  title(main = mirna)

  sumoflogsquares <<- sum((c.final - data)^2)
  print(sumoflogsquares)



sitesXcounts <- GetSitesXCounts(mirna,
                                        experiment,
                                        start,
                                        stop,
                                        sitelist,
                                        mirna.start = 3,
                                        mirna.stop = 6)
# Separate site sequences from data file.
if (sitelist == "12mers") {
  seqs <- rownames(sitesXcounts)
} else {
  seqs <- sitesXcounts[,1]
  sitesXcounts <- sitesXcounts[,-1]
}
# if (sitelist == "12mers") {
#   sitesXcounts_new <- t(sapply(1:((nrow(sitesXcounts)-1)/256),function(x){colSums(sitesXcounts[1:256+(x-1)*256,])}))
#   sitesXcounts_new <- rbind(sitesXcounts_new,sitesXcounts[nrow(sitesXcounts),,drop=FALSE])
#   sitesXcounts <- sitesXcounts_new
# }


print("23")
# Remove any rows with no reads whoatsoever.
sitesXcounts <- sitesXcounts[rowSums(sitesXcounts)>0,]
# Remove any rows for which there are no input reads:
sitesXcounts <- sitesXcounts[sitesXcounts[,1]>0,]

data <- sitesXcounts[,2:6]
data <- data.matrix(data)

print(dim(data))

rownames(data)[nrow(data)] <- "None"

num.kds <- nrow(data)
num.bgs <- 1
print(46)

# Define the vector of total site type concentrations:
c.I.tots <- Norm(sitesXcounts[, 1]) * k.c.lib
c.totals <- matrix(c.I.tots, nrow=length(c.I.tots), ncol=5, byrow=FALSE)
colnames(c.totals) <- colnames(data)
rownames(c.totals) <- rownames(data)
names(c.I.tots) <- rownames(data)
# Define the vector of AGO concentrations:
c.agos <- sapply(colnames(data), function(x) {
   as.numeric(x) / 100
  })
print(c.agos)




# plot(ecdf(log(c.I.tots)))
print(sort(c.I.tots,decreasing=TRUE)[1:10])
print(sort(c.I.tots)[1:10])

in_file <- paste0("/lab/solexa_bartel/mcgeary/AgoRBNS/", mirna, "/",
                   experiment, "/kds_PAPER/", start, "-", stop, "_", 
                   sitelist, "_", 3, "-", 6,
                   "_singlebg_logresiduals_combinedinput_PAPER_final_last.txt")




pars_table <- try(fread(in_file, sep = "\t", header=TRUE, colClasses = "numeric",stringsAsFactors = FALSE))

while(class(pars_table) == "try-error") {
  pars_table <- try(fread(in_file, sep = "\t", header=TRUE, colClasses = "numeric",stringsAsFactors = FALSE))

}
table_names <- names(pars_table[1])
pars_table_new <- matrix(unlist(pars_table), nrow=nrow(pars_table), ncol = ncol(pars_table), byrow=FALSE)
pars <- unlist(pars_table_new[nrow(pars_table_new),])
names(pars) <- table_names
pars_old <- pars
pars <- pars_old[-length(pars_old)]

print(length(pars))
print(sort(pars)[1:20])

  time_prior <- proc.time()

  # Split up the parameters into the kd and background parameters.
    kds  <- c(Logistic(pars[1 : (num.kds-1)], 10),1)
    kds_first <- kds
    bgs  <- rep(10^pars[(num.kds)], 5)
  B <- bgs[1]
    stock.ago <- 10^pars[num.kds + 1]
  names(kds) <- rownames(data)
  print("hi")
  f.mat <- matrix(
      rep(sapply(c.agos, function(ago.percent) {
           return(GetFreeAgo(kds, c.I.tots, ago.percent * stock.ago))
         }
         ), nrow(data)), nrow=nrow(data), ncol=ncol(data), byrow=TRUE)


  A.mat <- matrix(c.agos*stock.ago, nrow=nrow(data), ncol = ncol(data), byrow=TRUE)




  K.mat <- matrix(kds, nrow=nrow(data), ncol = ncol(data), byrow=FALSE)


  l.mat <- matrix(c.I.tots, nrow=nrow(data), ncol = ncol(data), byrow=FALSE)


  R.mat <- matrix(colSums(data), nrow = nrow(data), ncol = ncol(data), byrow=TRUE)



  L <- 100

  c.final <- (-R.mat*l.mat * (f.mat^2 + (L - A.mat)*f.mat + K.mat*B) *
                  (f.mat^3 - (2*A.mat + B - L - K.mat)*f.mat^2 +
                    (A.mat^2 + (B - L)*A.mat - L*B - 
                     (2*A.mat + B - L)*K.mat)*f.mat +
                    (A.mat^2 + (B - L)*A.mat - L*B)*K.mat)^-1)

  colnames(c.final) <- colnames(data)
  plot(c(c.final),c(data), col = rep(c(rgb(0,0,0,alpha =0.4),rgb(1,0,0,alpha=0.4),rgb(0,1,0,alpha=0.4),
                                      rgb(1,1,0,alpha=0.4), rgb(1,0,1, alpha = 0.4)),each=length(c(c.final))/5),log='xy')
  title(main = mirna)

  sumoflogsquares <<- sum((c.final - data)^2)
  print(sumoflogsquares)






}




Gradient <- function(pars) {
  time_prior <- proc.time()
  # Split up the parameters into the kd and background parameters.
  kds  <- c(Logistic(pars[1 : num.kds], 10),1)
  names(kds) <- rownames(data)
  bgs  <- rep(10^pars[(num.kds + 1)], 5)
  B <- bgs[1]
  stock.ago <- 10^pars[num.kds + 2]

   f.mat <- matrix(
      rep(sapply(c.agos, function(ago.percent) {
           return(GetFreeAgo(kds, c.I.tots, ago.percent * stock.ago))
         }
         ), nrow(data)), nrow=nrow(data), ncol=ncol(data), byrow=TRUE) 
  f.jvec <- sapply(c.agos, function(ago.percent) {
           return(GetFreeAgo(kds, c.I.tots, ago.percent * stock.ago))})

  A.mat <- matrix(c.agos*stock.ago, nrow=nrow(data), ncol = ncol(data), byrow=TRUE)
  A.jvec <- c.agos * stock.ago

  K.mat <- matrix(kds, nrow=nrow(data), ncol = ncol(data), byrow=FALSE)
  K.ivec <- kds

  l.mat <- matrix(c.I.tots, nrow=nrow(data), ncol = ncol(data), byrow=FALSE)
  l.ivec <- c.I.tots

  R.mat <- matrix(colSums(data), nrow = nrow(data), ncol = ncol(data), byrow=TRUE)
  R.jvec <- colSums(data)

  L <- 100
  # Get the amount of background binding by subtracting the bound from the
  # total sites in each experiment, normalizing. Must transpose to multiply
  # each column.

  time.init <- proc.time()

  c.vec_num <- -(
                 (l.ivec %*% t(R.jvec * f.jvec * (f.jvec + L - A.jvec)) +
                 (l.ivec * K.ivec * B)  %*% t(R.jvec))
                ) 


  C1.jvec <- L - B - 2 * A.jvec
  C2.jvec <- A.jvec^2 + (B - L) * A.jvec - L * B


  c.vec_dem <- t(
                (f.jvec^3 + f.jvec^2 * C1.jvec + f.jvec * C2.jvec) +
                t(K.ivec %*% t(f.jvec^2 + f.jvec * C1.jvec + C2.jvec))
                )^-1




c.final <- c.vec_num * c.vec_dem




  time.init <- proc.time()
  dF.base <- (1 + colSums(l.mat * K.mat * (f.mat + K.mat)^(-2)))^(-1)
  time.new <- proc.time()
  dF.base.jvec <- (1 + colSums(l.ivec * K.ivec * (t(f.jvec^2 + t(2 * K.ivec %*% t(f.jvec) + K.ivec^2)))^-1))^-1
  time.new2 <- proc.time()


  dF.dK.mat <- sweep(f.mat * l.mat * (f.mat + K.mat)^(-2),MARGIN=2,dF.base, "*")


  C1.mat <- L - B - 2*A.mat


  C2.mat <- A.mat^2 + (B - L)*A.mat - L*B


  # The d (each model point) d Free derivative:)
  dc.ai.dF <- R.mat * l.mat * (
    (
      f.mat^4
    ) + (
      2 * (L - A.mat) * f.mat^3
    ) + (
      ((L - A.mat)^2 + K.mat * (4 * B + A.mat)) * f.mat^2
    ) + (
      2 * K.mat * (A.mat * (L - B) + B * (K.mat - 2 * L) - (A.mat + B)^2) * f.mat
    ) + (
      K.mat * (A.mat^3 + 2 * A.mat^2 * (B - L) - B * (B - L) * (K.mat + L) +
               A.mat * (B^2 - 2 * B * K.mat - 3 * B * L + L^2))
    )
  ) * ((K.mat + f.mat)^2 * (C2.mat + C1.mat * f.mat + f.mat^2)^2)^(-1)

  time.init <- proc.time()

  # dc.ai.dF_vec_num <- R.mat * l.mat * (
  #   (
  #     f.mat^4
  #   ) + (
  #     2 * (L - A.mat) * f.mat^3
  #   ) + (
  #     ((L - A.mat)^2 + K.mat * (4 * B + A.mat)) * f.mat^2
  #   ) + (
  #     2 * K.mat * (A.mat * (L - B) + B * (K.mat - 2 * L) - (A.mat + B)^2) * f.mat
  #   ) + (
  #     K.mat * (A.mat^3 + 2 * A.mat^2 * (B - L) - B * (B - L) * (K.mat + L) +
  #              A.mat * (B^2 - 2 * B * K.mat - 3 * B * L + L^2))
  #   )
  # )

  # time.new <- proc.time()
  # print(dc.ai.dF_vec_num[1:2, 1:2])
  # print(time.new - time.init)

  # time.init <- time.new

  # dc.ai.dF.vec <- (
  #                   (l.ivec %*% (R.jvec * (f.jvec^4 + 2 * (L - A.jvec) * f.jvec^4) + (L - A.mat)^2)) +

  #                   ((l.jvec * ))


  dc.ai.dKj_specific <- R.mat * l.mat * (
    (
      f.mat^4
    ) + (
      (2 * C1.mat + A.mat) * f.mat^3
    ) + (
      (C2.mat + (C1.mat + A.mat) * C1.mat) * f.mat^2
    ) + (
      (C2.mat * (C1.mat + A.mat)) * f.mat
    )
  ) * ((K.mat + f.mat)^2 * (C2.mat + C1.mat * f.mat + f.mat^2)^2)^(-1)

  dc.ai.dA_specific <- R.mat * l.mat * (
    (
      -f.mat^3
    ) + (
      (2 * (A.mat - L)) * f.mat^2
    ) + (
      ((A.mat - L) * C1.mat - 2 * B * K.mat + C2.mat) * f.mat
    ) + (
      - C1.mat * B * K.mat
    )
  ) * ((K.mat + f.mat) * (C2.mat + C1.mat * f.mat + f.mat^2)^2)^(-1)



  dc.ai.dB <- R.mat * l.mat * (
    (
      -f.mat^3
    ) + (
      (2 * (A.mat - L) - K.mat) * f.mat^2
    ) + (
      ((L - A.mat) * (A.mat - L) - K.mat * (B + C1.mat)) * f.mat
    ) + (
      + A.mat * B * K.mat - B * K.mat * L - K.mat * C2.mat
    )
  ) * ((K.mat + f.mat) * (C2.mat + C1.mat * f.mat + f.mat^2)^2)^(-1)


  residuals <- (log(c.final+1) - log(data+1))*(c.final+1)^(-1)
  

  grad_derivs <- 10*exp(pars[1:num.kds]) * (exp(pars[1:num.kds]) + 1)^(-2) * 2* (colSums(residuals*dc.ai.dF) %*% t(dF.dK.mat[-nrow(dF.dK.mat),]) + rowSums(dc.ai.dKj_specific*residuals)[-nrow(dF.dK.mat)])
  

  dc.ai.dA <- sweep(sweep(dc.ai.dF,MARGIN=2,dF.base,"*") + dc.ai.dA_specific,
                    MARGIN=2,c.agos, "*")

  colnames(c.final) <- colnames(data)


  gradient_with_Ago <- log(10)*stock.ago*2*sum(residuals * dc.ai.dA)


  gradient_with_bg <- log(10)*B*2*sum(residuals * dc.ai.dB)

  return(c(grad_derivs, gradient_with_bg, gradient_with_Ago))
}

CalculateDerivative <- function(index,del,pars){
  start <- ModelLikelihood(pars)
  pars[index] <- pars[index]+del
  out <- ModelLikelihood(pars)
  return((out - start)/del)
}
