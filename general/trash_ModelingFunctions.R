GradientC <- function(pars, sXc, L=100, combined=TRUE, plot_=FALSE,
                      plotname=NULL) {
  data <- GetDataEquil(sXc)
  l <- GetInputEquil(sXc=sXc, combined=combined)
  dil <- sapply(colnames(data), as.numeric)/100
  kds <- 10^pars[1:length(l)]
  As <- 10^pars["AGO"]*dil
  bg <- 10^pars["bg"]
  as <- sapply(As, FreeAgoC, kds=kds, l=l)
  n_x <- nrow(sXc)
  n_col <- ncol(data)
  .C("GradientEquil", as.double(kds), as.double(l), as.double(L),
      as.double(bg), as.double(dil), as.double(As), as.double(as),
      as.double(as.numeric(as.matrix(data))), as.integer(n_x),
      as.integer(n_col),
      gradient=double(length(pars)))[["gradient"]]
}

CostC <- function(pars, sXc, L=100, combined=TRUE) {
  data <- GetDataEquil(sXc)
  l <- GetInputEquil(sXc=sXc, combined=combined)
  dil <- sapply(colnames(data), as.numeric)/100
  kds <- exp(pars[1:length(l)])
  # print("kds old")
  # print(kds)
  # kds_old <<- kds
  As <- exp(pars["AGO"])*dil
  # As_old <<- As
  bg <- exp(pars["bg"])
  # print("bg old")
  # print(bg)
  as <- sapply(As, FreeAgoC, kds=kds, l=l)
  # as_old <<- as
  # print(FreeAgo(kds=kds, l=l, A=As[1]))
  # print("As[1]")
  # print(As[1])
  # print("l")
  # print(l)
  # l_old <<- l
  # x1_old <<- sapply(as, function(as_) l*as_/(as_ + kds))
  n_x <- nrow(sXc)
  n_col <- ncol(data)
  # print(n_col)
  # print(data[, 1])
  tick <<- tick + 1
  out <- .C("CostEquil", as.double(kds), as.double(l), as.double(L),
             as.double(bg), as.double(dil), as.double(As), as.double(as),
             as.double(as.numeric(as.matrix(data))), as.integer(n_x),
             as.integer(n_col), cost = double(1))[["cost"]]
  # if (tick%%100 == 0) {
  #   print(out)
  # }
  out
}

sXa_of_pars.l.dil <- function(pars, l, dil, cols=NULL) {
  # Sites–by–[ago-bound] (nM)
  if (is.null(cols))
    cols <- 1:length(dil)
  kds     <- 10^pars[grep("_Kd", names(pars))]
  # kds.new <<- kds
  A.stock <- 10^pars["AGO"]
  # print("new")
  # print(pars["AGO"])
  # print(A.stock)

  # stock.ago.new <<- A.stock
  A.vec <- A.stock*dil
  # c.agos.new <<- A.vec
  # l.new <<- l
  return(sapply(A.vec, BoundRNA, kds=kds, l=l)[, cols])
}
# 1.b
Ds.a_Dlogp <- function(pars, l, dil, sXa, column=1) {
  bg <- 10^pars["bg"]
  kds     <- 10^pars[grep("_Kd", names(pars))]
  A <- 10^pars["AGO"]*dil[column]
  a <- A - sum(sXa[, column])
  # Partial derivative of all x with respect to free ago:
  print("ds.a_dda")
  ds.a_da <- l*kds/(a + kds)^2
  print("ds.a_dlogp")
  # Partial derivatives of each species with respect to the parameters:
  ds.a_dlogp <- structure(cbind(diag(-a*l/(a + kds)^2),
                                 rep(0, nrow(sXa)),
                                 dil[column]*ds.a_da),
                           dimnames=list(rownames(sXa), names(pars)))
  # Final term:
  print("final term")
  ds.a_dlogp - ds.a_da%o%colSums(ds.a_dlogp)/(1 + sum(ds.a_da))
}
# 2.a Function giving the total recovered site type concentrations:
sXm_of_pars.sXa <- function(pars, sXa, l, cols=NULL) {
  # Sites-by-[model] (nM)
  if (is.null(cols))
    cols <- 1:ncol(sXa)
  bg <- 10^pars["bg"]
  # bgs.new <<- bg
  # c.frees.new <<- l - sXa
  # c.bgs.new <<- bg*MatNorm(l - sXa)
  (sXa + bg*MatNorm(l - sXa))[, cols]
}
# 2.b.i
Ds.m_Dlogp <- function(pars, sXa, l, column=1) {
  bg <- 10^pars["bg"]
  out <- matrix(0, nrow=nrow(sXa), ncol=length(pars),
                dimnames=list(rownames(sXa), names(pars)))
  out[, "bg"] <- MatNorm(l - sXa)[, column]
  out
}
# 2.b.ii
Ds.m_Ds.a <- function(pars, sXa, l, column=1) {
  bg <- 10^pars["bg"]
  s.a <- sXa[, column]
  denom <- sum(l - s.a)
  row.mat <- matrix(bg*(l - s.a)/denom^2, nrow=nrow(sXa), ncol=nrow(sXa),
                dimnames=list(rownames(sXa), rownames(sXa)))
  out <- row.mat + diag(rep(1 - bg/denom, nrow(sXa)))
  # out <- diag(rep(1 - bg/denom, nrow(sXa)))
  out
}
Ds.m_Ds.a <- function(pars, sXa, l, column=1) {
  bg <- 10^pars["bg"]
  s.a <- sXa[, column]
  denom <- sum(l - s.a)
  row.mat <- matrix(bg*(l - s.a)/denom^2, nrow=nrow(sXa), ncol=nrow(sXa),
                dimnames=list(rownames(sXa), rownames(sXa)))
  out <- row.mat + diag(rep(1 - bg/denom, nrow(sXa)))
  # out <- diag(rep(1 - bg/denom, nrow(sXa)))
  out
}

# 3.a Cost function that normalizes the output of 2 and fits a multinomial model:
L_of_sXm <- function(sXm, data, cols.final=NULL) {
  if (is.null(cols.final))
    cols.final <- 1:ncol(sXm)
  # data <- data/1000
  -sum((data*log(MatNorm(sXm)))[, cols.final])
}
# 3.b
DL_DsXm <- function(sXm, data, cols=NULL) {
  # s.m : A single column of sXm, sites-by-model
  # s.c : A signle column of data, sites-by-counts
  if (is.null(cols))
    cols <- 1:ncol(sXm)
  # data <- data/1000
  sapply(1:ncol(sXm), function(j) {
    s.m <- sXm[, j]
    s.c <- data[, j]
    sum(s.c)/sum(s.m) - s.c/s.m
  })
}

DL_DsXm <- function(sXm, data, cols=NULL) {
  # s.m : A single column of sXm, sites-by-model
  # s.c : A signle column of data, sites-by-counts
  if (is.null(cols))
    cols <- 1:ncol(sXm)
  # data <- data/1000
  sapply(1:ncol(sXm), function(j) {
    s.m <- sXm[, j]
    s.c <- data[, j]
    sum(s.c)/sum(s.m) - s.c/s.m
  })
}



DL_DsXa <- function(sXm, data, cols=NULL) {
  # s.m : A single column of sXm, sites-by-model
  # s.c : A signle column of data, sites-by-counts
  if (is.null(cols))
    cols <- 1:ncol(sXm)
  # data <- data/1000
  sapply(1:ncol(sXm), function(j) {
    s.m <- sXm[, j]
    s.c <- data[, j]
    sum(s.c)/sum(s.m) - s.c/s.m
  })
}

FreeAgo <- function(kds, l, A, tol=3.000214e-13) { 
  # COMPATIBLE WITH MULTISITE
  # Args:
  # kds: Either a vector (1-site model), or a list of vectors (multisite) of 
  #      Ago-target site Kd values.
  # l: A vector of [total site type] for each single/multi site target RNA in
  #    the binding reaction
  # A: The concentration of free Ago in the binding reaction
  # tol: The tolerance for the root solve or the optimization.
  if (A > 0) {
    Residual <- function(a, exponent=1) {
      # Next two lines are generic functions that allow for single or multisite.
      oc <- sapply(kds, function(kd.set) sum(SiteOcc(a, kd.set)))
      l <- l*sapply(kds, length)
      return((A - a - sum(oc*l))^exponent)
    }
    a <- NaN
    try(a <- uniroot(Residual, c(0, A), tol=tol)$root)
    if (is.na(a)) {
      a <- optimize(Residual, c(0, A), exponent=2, tol=tol)$minimum
    }
    return(a)
  } else {
    return(0)
  }
}

