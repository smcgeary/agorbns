source("general/general.R")

     # 1. Get data table:
mirnas <- c("let-7a", "miR-1", "miR-155", "miR-124", "lsy-6")


# mirnas <- c("lsy-6")
sitelist <- "12mers"
library(multicore)

library(data.table)
par(mfcol=c(3, 5))
for (mirna in mirnas) {
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
pars_old_new <- pars
pars <- pars_old_new[-length(pars_old_new)]

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
  pars_x <- log10(Logistic(pars_old[-length(pars_old)],10))
  pars_y <- log10(Logistic(pars_old_new[-length(pars_old_new)],10))
  # plot(pars_x[grep("^CGATC",names(pars_x))], pars_y[grep("^CGATC",names(pars_y))],
  #   xlim=c(min(c(pars_x, pars_y)),max(c(pars_x, pars_y))),
  #   ylim=c(min(c(pars_x, pars_y)),max(c(pars_x, pars_y))))

  # segments(min(c(pars_x, pars_y)),
  #          min(c(pars_x, pars_y)),
  #          max(c(pars_x, pars_y)),
  #          max(c(pars_x, pars_y)),lty=2)


  plot(pars_x, pars_y,
    xlim=c(min(c(pars_x, pars_y)),max(c(pars_x, pars_y))),
    ylim=c(min(c(pars_x, pars_y)),max(c(pars_x, pars_y))))

  segments(min(c(pars_x, pars_y)),
           min(c(pars_x, pars_y)),
           max(c(pars_x, pars_y)),
           max(c(pars_x, pars_y)),lty=2)


# plot(pars_x[1:(num.kds-1)]-pars_y[1:(num.kds-1)],c.I.tots[1:(num.kds-1)],log='y')

}

PlotKds <- function(){
    plot(pars_x, pars_y,
    xlim=c(min(c(pars_x, pars_y)),max(c(pars_x, pars_y))),
    ylim=c(min(c(pars_x, pars_y)),max(c(pars_x, pars_y))))

  segments(min(c(pars_x, pars_y)),
           min(c(pars_x, pars_y)),
           max(c(pars_x, pars_y)),
           max(c(pars_x, pars_y)),lty=2)

}



Gradient <- function(pars) {
  time_prior <- proc.time()
  # Split up the parameters into the kd and background parameters.
  kds  <- c(Logistic(pars[1 : num.kds-1], 10),1)
  names(kds) <- rownames(data)
  bgs  <- rep(10^pars[(num.kds)], 5)
  B <- bgs[1]
  stock.ago <- 10^pars[num.kds + 1]

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
  # Get the amount of background binding by subtracting the bound from the
  # total sites in each experiment, normalizing. Must transpose to multiply
  # each column.
  k_8 <- kds["8mer"]
  c.final <- (-R.mat*l.mat * (f.mat^2 + (L - A.mat)*f.mat + K.mat*B) *
                  (f.mat^3 - (2*A.mat + B - L - K.mat)*f.mat^2 +
                    (A.mat^2 + (B - L)*A.mat - L*B - 
                     (2*A.mat + B - L)*K.mat)*f.mat +
                    (A.mat^2 + (B - L)*A.mat - L*B)*K.mat)^-1)
  # print("up to c.final")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new


  dF.base <- (1 + colSums(l.mat * K.mat * (f.mat + K.mat)^(-2)))^(-1)

  # print("dF.base")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new



  dF.dK.mat <- sweep(f.mat * l.mat * (f.mat + K.mat)^(-2),MARGIN=2,dF.base, "*")

  # print("dF.dK.mat")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new


  C1.mat <- L - B - 2*A.mat

  # print("C1.mat")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new


  C2.mat <- A.mat^2 + (B - L)*A.mat - L*B

  # print("C2.mat")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new


  # The d (each model point) d Free derivative:)
  # print("dc.ai.dF")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new

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

  # print("dc.ai.dKj_specific")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new



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

  # print("dc.ai.dA_specific")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new


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

  # print("dc.ai.dB")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new



  # grad_derivs <- sapply(seq(length(kds)), function(index) {
  #   base <- sweep(dc.ai.dF,MARGIN=2,c(dF.dK.mat[index,]),"*")
  #   base[index,] <- base[index,] + dc.ai.dKj_specific[index,]
  #   return(10*exp(pars[index]) * (exp(pars[index]) + 1)^(-2) * 2 * sum((log(c.final) - log(data+1)*base*c.final^(-1))))
  #   })
  lc.final <- log(c.final)
  ldata <- log(data + 1)

  time_new <- proc.time()
  time_prior <- time_new


  # grad_derivs <- sapply(seq(length(kds)), function(index) {
  #   base <- sweep(dc.ai.dF,MARGIN=2,c(dF.dK.mat[index,]),"*")
  #   base[index,] <- base[index,] + dc.ai.dKj_specific[index,]
  #   return(10*exp(pars[index]) * (exp(pars[index]) + 1)^(-2) * 2 * sum((lc.final - ldata)*base*c.final^(-1)))
  #   })

  # time_new <- proc.time()
  # print('no parallel')
  # print(time_new - time_prior)
  # time_prior <- time_new
  residuals <- lc.final - ldata
  residuals_norm <- residuals * c.final^(-1)

  raw_deriv_base_temp <- sweep(dc.ai.dF,MARGIN=2,c(dF.dK.mat[2992,]),"*")
      raw_deriv_base_temp[2992,] <- raw_deriv_base_temp[2992,] + dc.ai.dKj_specific[2992,]
    raw_deriv <<- 10*exp(pars[2992]) * (exp(pars[2992]) + 1)^(-2) * 2 * residuals_norm*raw_deriv_base_temp


  grad_derivs <- unlist(mclapply(seq(num.kds-1), function(index) {
    base <- sweep(dc.ai.dF,MARGIN=2,c(dF.dK.mat[index,]),"*")
    base[index,] <- base[index,] + dc.ai.dKj_specific[index,]
    return(10*exp(pars[index]) * (exp(pars[index]) + 1)^(-2) * 2 * sum(residuals_norm*base))

    }, mc.cores=16))

  time_new <- proc.time()
  print("parallel")
  print(time_new - time_prior)

    # print("grad deriv")

    # time_new <- proc.time()
    # print(time_new - time_prior)
    # time_prior <- time_new


  dc.ai.dA <- sweep(sweep(dc.ai.dF,MARGIN=2,dF.base,"*") + dc.ai.dA_specific,
                    MARGIN=2,c.agos, "*")



  # print("AGO deriv")
  # time_new <- proc.time()
  # print(time_new - time_prior)
  # time_prior <- time_new



  colnames(c.final) <- colnames(data)
  use <- which(c.final > 0 & data > 0)


  gradient_with_Ago <- log(10)*stock.ago*2*sum((log(c.final) - log(data + 1))*dc.ai.dA * c.final^(-1))
   Ago_check <- log(10)*stock.ago*2*colSums((log(c.final) - log(data + 1))*dc.ai.dA * c.final^(-1))


  gradient_with_bg <- log(10)*B*2*sum((log(c.final) - log(data + 1))*dc.ai.dB * c.final^(-1))
  bg_check <- log(10)*B*2*colSums((log(c.final) - log(data + 1))*dc.ai.dB * c.final^(-1))

  return(c(grad_derivs, gradient_with_bg, gradient_with_Ago))
}

CalculateDerivative <- function(index,del,pars){
  start <- ModelLikelihood(pars)
  pars[index] <- pars[index]+del
  out <- ModelLikelihood(pars)
  return((out - start)/del)
}
