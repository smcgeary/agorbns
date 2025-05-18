#23456789!123456789@123456789#123456789$123456789%123456789^123456789&123456789*
# Required to compile the C code:
system("R CMD SHLIB general/AGO_RBNS_functions.c")

# If this does not work, this is an alternative approach that is less robust.
## UNCOMMENT IF above system command doesn't work #############################
# Get the R library path from the environment variable
# conda_path <- system("echo $CONDA_PREFIX", intern=TRUE)
# Compile the C file
# system(paste0("gcc -std=gnu99 -I", conda_path, "/lib/R/include -DNDEBUG -fpic -O3 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 -c general/AGO_RBNS_functions.c -o general/AGO_RBNS_functions.o"))
# Create the shared library
# system(paste0("gcc -std=gnu99 -shared -L", conda_path, "/lib/R/lib -Wl,-Bsymbolic-functions -Wl,-z,relro -o general/AGO_RBNS_functions.so general/AGO_RBNS_functions.o -L", conda_path, "/lib/R/lib -lR"))
## END SECTION TO UNCOMMENT ####################################################


dyn.load("general/AGO_RBNS_functions.so")
options(width=100)
# 1.CORE BIOCHEMICAL FUNCTIONS FOR EQUILIBRIUM:_________________________________
SiteOcc <- function(a, kds) { # AGNOSTIC TO MULTISITE
  # Args:
  # a: [free AGO] in the binding reaction
  # kds: a list of all kds corresponding to the individual site types including
  # Returns:
  return(a/(a + kds))
}


# BetaDist <- function(b_mean, b_var, n_bin, exp=2.5, plot.=FALSE) {
#   p <- seq(0, 1, length.out=n_bin + 1)

#   p <- (p[1:n_bin] + p[2:(n_bin + 1)])/2
#   C <- b_mean*(1 - b_mean)/b_var - 1
#   shape1 <- b_mean*C
#   shape2 <- (1 - b_mean)*C
#   if (length(shape1) == 1) {
#     if (plot.) {
#       p_ecdf <- pbeta(p, shape1=shape1, shape2=shape2)
#       plot(p, p_ecdf, type="l")
#       p_norm <- pnorm(p, mean=b_mean, sd=sqrt(b_var))
#       lines(p, p_norm, col="red")    
#     }
#     p_f <- dbeta(p, shape1=shape1, shape2=shape2)  
#   } else {
#     if (plot.) {
#       plot(p, rep(0, length(p)), type="l", xlim=c(0, 1), ylim=c(0, 1))
#     }
#     p_f <- t(sapply(1:length(shape1), function(ind) {
#       if (plot.) {
#         p_ecdf <- pbeta(p, shape1=shape1[ind], shape2=shape2[ind])
#         lines(p, p_ecdf)            
#       }
#       dbeta(p, shape1=shape1[ind], shape2=shape2[ind])  
#     }))
#   }
#   p_f <- p_f^exp
#   return(list(p, p_f))
# }




##3456789!123456789@123456789#123456789$123456789%123456789^123456789&123456789*
FreeAgo <- function(kds, l, A) {
  # COMPATIBLE WITH MULTISITE
  # Args:
  # kds: Either a vector (1-site model), or a list of vectors (multisite) of 
  #      Ago-target site Kd values.
  # l: A vector of [total site type] for each single/multi site target RNA in
  #    the binding reaction
  # A: The concentration of free Ago in the binding reaction
  # tol: The tolerance for the root solve or the optimization.
  # A <- 1
  # print(kds)
  # print(l)
  kds_check <<- kds
  l_check <<- l
  if (A > 0) {
        C_function <-.C("FreeAgoR", as.double(A), as.double(l), as.double(kds),
                        as.integer(length(kds)), a_f = double(1), n_steps = integer(1),
                        xbs = double(20))
        # print(C_function[["n_steps"]])
        # print(C_function[["xbs"]])
        return(C_function[["a_f"]])

  } else {
    return(0)
  }
}

# FreeAgoR <- function(kds, l, A) {
#   # COMPATIBLE WITH MULTISITE
#   # Args:
#   # kds: Either a vector (1-site model), or a list of vectors (multisite) of 
#   #      Ago-target site Kd values.
#   # l: A vector of [total site type] for each single/multi site target RNA in
#   #    the binding reaction
#   # A: The concentration of free Ago in the binding reaction
#   # tol: The tolerance for the root solve or the optimization.
#   cost_func <- function(a, A_j) {
#     res <- A_j - a - sum(l*a/(a + kds))
#     res
#   }
#   sapply(A, function(A_j) {
#     if (A_j > 0) {
#       return(uniroot(cost_func, A_j=A_j, interval=c(0, A_j), tol=2.220446e-16)$root)
#     } else {
#       return(0)
#     }
#   })
# }

# FreeAgoWithStructureR <- function(kds, l, A, b_mean, b_var, n_bin) { 
#   # COMPATIBLE WITH MULTISITE
#   # Args:
#   # kds: Either a vector (1-site model), or a list of vectors (multisite) of 
#   #      Ago-target site Kd values.
#   # l: A vector of [total site type] for each single/multi site target RNA in
#   #    the binding reaction
#   # A: The concentration of total Ago in the binding reaction
#   # tol: The tolerance for the root solve or the optimization.
#   p_list <- SubfunctionCall(BetaDist)
#   p <- p_list[[1]]
#   p_f <- p_list[[2]]
#   if (length(b_mean) == 1) {
#     l_matrix <- l %*% t(p_f)/length(p)    
#   } else {
#     l_matrix <- l * p_f/length(p)
#   }
#   kds_matrix <- kds %*% t(1/p)
#   if (length(b_mean) == 1) {
#     l_matrix_1 <<- l_matrix
#     kd_matrix_1 <<- kds_matrix
#   } else {
#     l_matrix_2 <<- l_matrix
#     kd_matrix_2 <<- kds_matrix
#   }
#   cost_func <- function(a, A_j) {
#     res <- A_j - a - sum(l_matrix*a/(a + kds_matrix))
#     res
#   }
#   sapply(A, function(A_j) {
#     if (A_j > 0) {
#       return(uniroot(cost_func, A_j=A_j, interval=c(0, A_j), tol=2.220446e-16)$root)
#     } else {
#       return(0)
#     }
#   })
# }





# TotalBoundAgo <- function(kds, l, A, tol=.Machine$double.eps) { 
#   # COMPATIBLE WITH MULTISITE
#   # Args:
#   # kds: Either a vector (1-site model), or a list of vectors (multisite) of 
#   #      Ago-target site Kd values.
#   # l: A vector of [total site type] for each single/multi site target RNA in
#   #    the binding reaction
#   # A: The concentration of free Ago in the binding reaction
#   # tol: The tolerance for the root solve or the optimization.
#   if (A > 0) {
#     Residual <- function(X, exponent=1) {
#       # Next two lines are generic functions that allow for single or multisite.
#       oc <- sapply(kds, function(kd.set) sum(SiteOcc(A - X, kd.set)))
#       l <- l*sapply(kds, length)
#       return((X - sum(oc*l))^exponent)
#     }
#     X <- NaN
#     try(X <- uniroot(Residual, c(0, A), tol=tol)$root)
#     if (is.na(X)) {
#       X <- optimize(Residual, c(0, A), exponent=2, tol=tol)$minimum
#     }
#     return(X)
#   } else {
#     return(0)
#   }
# }

BoundRNA <- function(kds, l, A, total=FALSE) {
  # COMPATIBLE WITH MULTISITE
  # Args:
  # kds: Either a vector (1-site model), or a list of vectors (multisite) of 
  #      Ago-target site Kd values.
  # l: A vector of [total site type] for each single/multi site target RNA in
  #    the binding reaction
  # A: The concentration of free Ago in the binding reaction
  if (total) a <- A - TotalBoundAgo(kds, l, A)
  else       a <- FreeAgo(kds, l, A)
  # Determine the percent of each target site that is entirely unoccupied:
  unoc <- sapply(kds, function(kd.set) prod(1 - SiteOcc(a, kd.set)))
  return(l*(1 - unoc))
}




# BoundRNAWithStrucR <- function(kds, l, A, b_mean, b_var, n_bin) {
#   # COMPATIBLE WITH MULTISITE
#   # Args:
#   # kds: Either a vector (1-site model), or a list of vectors (multisite) of 
#   #      Ago-target site Kd values.
#   # l: A vector of [total site type] for each single/multi site target RNA in
#   #    the binding reaction
#   # A: The concentration of free Ago in the binding reaction
#   a <- FreeAgoWithStructureR(kds, l, A, b_mean, b_var, n_bin)

#   p_list <- SubfunctionCall(BetaDist)
#   p <- p_list[[1]]
#   p_f <- p_list[[2]]
#   if (length(b_mean) == 1) {
#     l_matrix <- l %*% t(p_f)/length(p)  
#   } else {
#     l_matrix <- l * p_f/length(p)
#   }
#   kds_matrix <- kds %*% t(1/p)
#   # Determine the percent of each target site that is entirely unoccupied:
#   occs <- SiteOcc(a, kds_matrix)
#   return(rowSums(l_matrix*occs))
# }




################################################################################
# EQUILIBRIUM COST FUNCTION
# 1.a Function giving the marix of all the ago-bound site type concentrations:

GetInputEquil <- function(sXc, combined=TRUE, conc=100) {
  if (combined) ind <- which(colnames(sXc) == "I_combined")
  else          ind <- 1
  s.c <- sXc[, ind]
  s.c[s.c==0] <- 1
  l <- Norm(s.c)*conc
  names(l) <- rownames(sXc)
  l
}

GetDataEquil <- function(sXc, combine_reps=FALSE) {
  cols <- setdiff(colnames(sXc),
                  c("I", "I_combined", "I.1", "I,1", "I,2", "0", "0,1", "0,2"))
  data <- data.matrix(sXc[, cols])
  if (combine_reps) {
    # Extract the sample names from the column names (40,1 becomes 40)
    col_samples <- gsub("^(.*)\\,(.*)$", colnames(data), replace="\\1",
                        perl=TRUE)
    # Get the unique samples from this list.
    data_cols <- unique(col_samples)
    # Pre-allocate the new data matrix.
    data_new <- data.matrix(matrix(NA, nrow=nrow(data), ncol=length(data_cols),
                                   dimnames=list(rownames(data), data_cols)))
    # Loop over the unique column names, using the indeces for those samples in
    # `col_samples` to add the appropriate colums.
    sapply(data_cols, function(data_col) {
      col_inds <- which(col_samples == data_col)
      data_new[, data_col] <<- rowSums(data[, col_inds, drop=FALSE])
    })
    # Assign the new data matrix to the original.
    data <- data_new
  }
  data
}

EquilEnrichments <- function(mat, l) {
  apply(mat, 2, function(x) Norm(x)/Norm(l))
}








CostC <- function(pars, data, dil, l, L, Y, n_i, n_j, n_mir, zero_grad, n_z,
                  fixed, norm_constant, upper_, lower_, plot_=FALSE,
                  plotname=NULL, tempname=NULL) {
  # This is the function that is used to solve for Kds in the FitSiteKds.R
  # script that is the workhorse function used in the equilibrium Kd paper.
  # message("in Cost C.")
  # print(length(which(is.na(pars))))
  # print(length(which(is.na(data))))
  # print(length(which(is.na(dil))))
  # print(dil)
  # print(length(which(is.na(l))))
  # print(n_i)
  # print(n_j)
  # print(n_mir)
  # print(length(zero_grad))
  # print(n_z)
  # print(fixed)
  # print(norm_constant)
  # print(length(upper_))
  # print(length(lower_))
  # print(sum(which(pars > upper_)))
  # print(sum(which(pars < lower_)))
  out <- .C("CostEquilMultiMirnas",
            as.double(pars),
            as.double(data),
            as.double(dil),
            as.double(l),
            as.double(L),
            as.double(Y),
            as.integer(n_i),
            as.integer(n_j),
            as.integer(n_mir),
            as.double(norm_constant),
            cost=double(1),
            cost_each=double(sum(n_i*n_j)),
            A_j_check=double(sum(n_j)),
            a_j_check=double(sum(n_j)),
            k_i_check=double(n_i))
  cost <- out[["cost"]]
  # if (t_global%%10 == 0) {
  #   print(cost)
  # #   time_new <- proc.time()[3]
  # #   print(time_new - time_interval)
  # #   time_interval <<- time_new

  # #   # print(zero_grad)
  # #   # print(n_z)
  # #   # print(exp(pars[1:16]))
  # #   # print(exp(pars[n_i + 1]))
  # #   # print(exp(pars[n_i + 2]))
  # #   out_mat <- matrix(out[["cost_each"]], nrow=n_i, byrow=FALSE)
  # #   out_mat <<- out_mat
  # #   data_mat <- matrix(data, nrow=n_i, byrow=FALSE)
  # #   data_mat <<- data_mat
  # #   data_mat_norm <- as.numeric(t(t(data_mat)/colSums(data_mat)))
  # #   data_mat_norm <<- data_mat_norm
  # #   l <<- l
  # #   # # message("data in C function")
  # #   # # print(head(out_mat))
  # #   # # R_mat <- t(t(data_mat)/colSums(data_mat))/Norm(l)
  # #   # # print(head(R_mat))
  # #   # # print(head(data))
  # #   # # dev.set(2)
  # #   # # plot(l, data_mat[, 1], log="xy")
  # #   # # dev.set(3)
  # #   # # plot(l, out_mat[, 1], log="xy")
  # #   # # dev.set(4)
  # #   cols <- rgb(0, 0, 0, alpha=0.2)
  # #   plot(data_mat_norm, out[["cost_each"]], col=cols, log="xy", pch=20)
  # }
  # t_global <<- t_global + 1
  cost
}

# # This function is used specifically for the three-prime paper; it allows the
# # deviation of the Kds from the collapsed or suppcomp value to be penalized so
# # as to prevent the very high Kd values that happen sometimes.
# CostCPrior <- function(
#   pars, data, dil, l, L, Y, n_i, n_j, n_mir, zero_grad, n_z, fixed,
#   norm_constant, lambda, lambda0, reg_mean, reg_inds, reg_vals, n_reg, upper_,
#   lower_, plot_=FALSE, plotname=NULL, tempname=NULL
# ) {
#   # This is the function that is used to solve for Kds in the FitSiteKds.R
#   # script that is the workhorse function used in the equilibrium Kd paper.
#   out <- .C("CostEquilMultiMirnasPrior",
#             as.double(pars),
#             as.double(data),
#             as.double(dil),
#             as.double(l),
#             as.double(L),
#             as.double(Y),
#             as.integer(n_i),
#             as.integer(n_j),
#             as.integer(n_mir),
#             as.double(norm_constant),
#             as.double(lambda),
#             as.double(lambda0),
#             as.double(reg_mean),
#             as.integer(reg_inds),
#             as.double(reg_vals),
#             as.integer(n_reg),
#             cost=double(1),
#             cost_each=double(sum(n_i*n_j)),
#             A_j_check=double(sum(n_j)),
#             a_j_check=double(sum(n_j)),
#             k_i_check=double(n_i),
#             residual_check=double(n_reg),
#             reg_ind_check=integer(n_reg))
#   cost <- out[["cost"]]
#   cost
# }


# CostNew <- function(pars, sXc, combined=TRUE, cols.final=NULL) {
#   # This is an R version of the C function that gives the same answer as the 
#   # C function. It's useful in terms of seeing a matrix, non-loop version of
#   # the specific formalism that is used in the C script, which differs from that
#   # of the methods in the paper because the methods in the paper are written
#   # better.
# 	l <- SubfunctionCall(GetInputEquil)
#   print(l)
# 	L <- sum(l)
#   data <- SubfunctionCall(GetDataEquil)
#   print(head(data))
#   dil <- sapply(colnames(data), as.numeric)/100
#   print(dil)
#   bg <- 10^pars[grep("bg", names(pars))]
#   kds <- 10^pars[grep("_Kd", names(pars))]
#   if (is.null(cols.final)) cols.final <- 1:ncol(data)
#   out <- sum(sapply(cols.final, function(column) {
#     A <- 10^pars[grep("AGO", names(pars))]*dil[column]
#     a <- FreeAgo(kds, l, A)
#     c <- data[, column]
#     x <- l*a/(a + kds)
#     f <- l - x
#     m <- x + bg*f/sum(f)
#     -sum(c*log(m/sum(m)))
#   }))
#   out
# }


# CostMethodsCheck <- function(theta, l, y) {
#   # Here the function takes in theta, l, and y, exactly as it is written in the
#   # methods of the paper.
#   print("start")
#   m <- nrow(y)
#   n <- ncol(y)
#   kds <- exp(theta[1:m])
#   a <- exp(theta[m + 2])
#   b <- exp(theta[m + 1])
#   L <- sum(l)
#   x <- sapply(colnames(x), function(col_j) {
#     a_j <- a*as.numeric(col_j)/100
#     a_j_f <- FreeAgo(kds, l, a_j)
#     x_ij <- l*(a_j_f/(a_j_f + kds)*(1 - b/(L - sum(l*a_j_f/(a_j_f + kds)))) + b/(L - sum(l*a_j_f/(a_j_f + kds))))
#     return(x_ij)
#   })
#   x <<- x
#   X <- colSums(x)
#   Y <- colSums(y)
#   return(sum(Y*log(X)) - sum(y*log(x)))
# }

# CostMethodsCheckNew <- function(theta, l, y) {
#   # Here the function takes in theta, l, and y, exactly as it is written in the
#   # methods of the paper.
#   print("start")
#   m <- nrow(y)
#   n <- ncol(y)
#   kds <- exp(theta[1:m])
#   a <- exp(theta[m + 2])
#   b <- exp(theta[m + 1])
#   L <- sum(l)
#   x <- sapply(colnames(x), function(col_j) {
#     a_j <- a*as.numeric(col_j)/100
#     a_j_f <- FreeAgo(kds, l, a_j)
#     x_ij <- l*(a_j_f/(a_j_f + kds)*(1 - b/(L - sum(l*a_j_f/(a_j_f + kds)))) + b/(L - sum(l*a_j_f/(a_j_f + kds))))
#     return(x_ij)
#   })
#   x <<- x
#   X <- colSums(x)
#   Y <- colSums(y)
#   print(sum(Y*log(X)))
#   print(-sum(y*log(x)))
#   return(sum(Y*log(X)) - sum(y*log(x)) + sum(lgamma(y+1)) - sum(lgamma(Y+1)))
# }



GradC <- function(pars, data, dil, l, L, Y, n_i, n_j, n_mir, zero_grad, n_z,
                  fixed, norm_constant, upper_, lower_, plot_=FALSE,
                  plotname=NULL, tempname=NULL) {
  out <- .C("GradEquilMultiMirna", as.double(pars), as.double(data),
            as.double(dil), as.double(l), as.double(L), as.double(Y),
            as.integer(n_i), as.integer(n_j), as.integer(n_mir),
            as.integer(zero_grad), as.integer(n_z), as.integer(fixed),
            as.double(norm_constant), grad=double(n_i + 2*n_mir))
  return(out[["grad"]])
}

# GradCPrior <- function(
#   pars, data, dil, l, L, Y, n_i, n_j, n_mir, zero_grad, n_z, fixed,
#   norm_constant, lambda, lambda0, reg_mean, reg_inds, reg_vals, n_reg, upper_, lower_,
#   plot_=FALSE, plotname=NULL, tempname=NULL
# ) {
#   out <- .C("GradEquilMultiMirnaPrior", as.double(pars), as.double(data),
#             as.double(dil), as.double(l), as.double(L), as.double(Y),
#             as.integer(n_i), as.integer(n_j), as.integer(n_mir),
#             as.integer(zero_grad), as.integer(n_z), as.integer(fixed),
#             as.double(norm_constant), as.double(lambda), as.double(lambda0),
#             as.double(reg_mean), as.integer(reg_inds), as.double(reg_vals),
#             as.integer(n_reg), grad=double(n_i + 2*n_mir))
#   return(out[["grad"]])
# }



# GradientNew <- function(pars, sXc, combined=TRUE, cols.final=NULL) {
#   print("start")
#    l <- SubfunctionCall(GetInputEquil)
#    L <- sum(l)
#   data <- SubfunctionCall(GetDataEquil)
#   dil <- sapply(colnames(data), as.numeric)/100
#   bg <- 10^pars["bg"]
#   kds <- 10^pars[grep("_Kd", names(pars))]
#   if (is.null(cols.final))
#     cols.final <- 1:ncol(data)
#   out <- rowSums(sapply(cols.final, function(column) {
#     A <- 10^pars["AGO"]*dil[column]
#     a <- SubfunctionCall(FreeAgo)

#     x <- l*a/(a + kds)
#     f <- l - x
#     F <- sum(f)
#     c <- data[, column]
#     delxdela <- kds*l/(a + kds)^2
#     delXdela <- sum(delxdela)
#     delxdelkds <- -a*l/(a + kds)^2

#     m <- x + bg*f/F
#     dmdb <- f/F
#     dmdx <- bg*f/F^2
#     dmdxij <- (F - bg)/F
#     dxdA <- delxdela/(1 + sum(delxdela))
#     dXdkds <- delxdelkds/(1 + sum(delxdela))  # Right

#     dLdm <- sum(c)/sum(m) - c/m
#     dLdbg <- sum(dLdm*dmdb)  # Right
#     dLdkds <- (delxdelkds*(sum(dLdm*dmdx) + dLdm*dmdxij) -
#                dXdkds*sum(dLdm*(dmdx*sum(delxdela) + dmdxij*delxdela)))
#     dLdA  <- dil[column] * sum(dLdm*(dmdx*sum(dxdA) + dmdxij*dxdA)) # Right
#   log(10)*(10^pars)*c(dLdkds, dLdbg, dLdA)
#   }))
#   print(out)
#   out["None_Kd"] <- 0
#   print(out)
#   out
# }




# GradientMethodsCheck <- function(theta, l, y) {
#   print("start")
#   m <- nrow(sXc)
#   n <- ncol(sXc)
#   kds <- exp(theta[1:m])
#   a <- exp(theta[m + 2])
#   b <- exp(theta[m + 1])
#   L <- sum(l)
#   out <- rep(0, m + 2)

#   c <- sapply(colnames(x), function(col_j) {
#     a_j <- a*as.numeric(col_j)/100
#     a_j_f <- FreeAgo(kds, l, a_j)
#     c_ij <- l*(a_j_f/(a_j_f + kds))
#     return(c_ij)
#   })
#   x <- sapply(colnames(x), function(col_j) {
#     a_j <- a*as.numeric(col_j)/100
#     a_j_f <- FreeAgo(kds, l, a_j)
#     x_ij <- l*(a_j_f/(a_j_f + kds)*(1 - b/(L - sum(l*a_j_f/(a_j_f + kds)))) + b/(L - sum(l*a_j_f/(a_j_f + kds))))
#     return(x_ij)
#   })
#   C <- colSums(c)
#   X <- colSums(x)
#   Y <- colSums(y)
#   del_f_del_x <- sapply(1:ncol(x), function(j) {
#     Y[j]/X[j] - y[, j]/x[, j]
#   })
#   print(del_f_del_x)
#   f_1_over_L_minus_C <- 1/(L - C)
#   Xi <- sapply(colnames(x), function(col_j) {
#     a_j <- a*as.numeric(col_j)/100
#     a_j_f <- FreeAgo(kds, l, a_j)
#     Xi_j <- 1/(1 - sum(l*kds/(a_j_f + kds)^2))
#     return(Xi_j)
#   })
#   out[m+1] <- sum(sapply(1:ncol(x), function(j) {
#     b*f_1_over_L_minus_C[j]*sum(del_f_del_x[, j]*(l - c[, j]))
#   }))
#   out
# }




# CostC_8mer <- function(pars, data, dil, l, l_adj, L, Y, n_i, n_j, n_mir, zero_grad,
#                        n_z, i_6m8, fixed, upper_, lower_, plot_=FALSE,
#                        plotname = NULL, tempname=NULL) {
#   out <- .C("CostEquilMultiMirnas_8mer",
#             as.double(pars),
#             as.double(data),
#             as.double(dil),
#             as.double(l),
#             as.double(l_adj),
#             as.double(L),
#             as.double(Y),
#             as.integer(n_i),
#             as.integer(n_j),
#             as.integer(n_mir),
#             as.integer(i_6m8),
#             cost=double(1),
#             cost_each=double(sum(n_i*n_j)),
#             A_j_check=double(sum(n_j)))
#   cost <- out[["cost"]]
#   tick <<- tick + 1
#   if (tick%%10000 == 0) {
#     print(cost)
#   }
#   cost
# }

# GradC_8mer <- function(pars, data, dil, l, l_adj, L, Y, n_i, n_j, n_mir, zero_grad,
#                        n_z, i_6m8, fixed, upper_, lower_, plot_=FALSE,
#                        plotname = NULL, tempname=NULL) {
#   out <- grad(CostC_8mer, pars, data=data, dil=dil, l=l, l_adj=l_adj, L=L, Y=Y,
#        n_i=n_i, n_j=n_j, n_mir=n_mir, i_6m8=i_6m8, fixed=fixed, method="simple")
#   out[zero_grad] <- 0
#   out

# }





# CostC12mers <- function(pars, data, dil, l_all, L, Y, n_i, n_j, n_mir, n_pos,
#                         zero_grad, n_z, upper_, lower_, plot_=FALSE, plotname = NULL,
#                   tempname=NULL) {
#   # print("in modeling function")
#   # print(n_j)
#   # message("dil:")
#   # print(dil)
#   # print(pars)
#   out <- .C("CostEquilMultiMirnas12mers",
#             as.double(pars),
#             as.double(data),
#             as.double(dil),
#             as.double(l_all),
#             as.double(L),
#             as.double(Y),
#             as.integer(n_i),
#             as.integer(n_j),
#             as.integer(n_mir),
#             as.integer(n_pos),
#             cost=double(1),
#             cost_each=double(sum(n_j)*n_i))
#   # print(out[["A_j_check"]])
#   # print(out[["j_check"]])

#   # print(head(cbind(data, out[["model_check"]])))
#   # break
#   cost <- out[["cost"]]
#   tick <<- tick + 1
#   # print(tick)
#   # print(cost)
#   cost_each_global <<- out[["cost_each"]]
#   if (tick%%20 == 0) {
#   #   kds_1 <- pars[1:(n_i)]
#   #   kds_2 <- pars[(n_i + 1):(2*n_i)]
#   #   kds_3 <- pars[(2*n_i + 1):(3*n_i)]
#   #   kds_4 <- pars[(3*n_i + 1):(4*n_i)]
#   #   kds_5 <- pars[(4*n_i + 1):(5*n_i)]
#   #   write.table(file="temp_1-4.txt", kds_1, sep="\t", quote=FALSE)
#   #   write.table(file="temp_2-5.txt", kds_2, sep="\t", quote=FALSE)
#   #   write.table(file="temp_3-6.txt", kds_3, sep="\t", quote=FALSE)
#   #   write.table(file="temp_4-7.txt", kds_4, sep="\t", quote=FALSE)
#   #   write.table(file="temp_5-8.txt", kds_5, sep="\t", quote=FALSE)
#     print(cost)
#   #   # print(length(data))
#   #   random_rows <- sample(length(data), 1e5)
#   #   plot(data[random_rows], out[["model_check"]][random_rows], log='xy',
#   #        xlim=c(1, 1e6), ylim=c(1, 1e6))
#   }
#   cost
# }

# GradC12mers <- function(pars, data, dil, l_all, L, Y, n_i, n_j, n_mir, n_pos,
#                         zero_grad, n_z, upper_, lower_, plot_=FALSE, plotname = NULL,
#                         tempname=NULL) {
#   out <- .C("GradEquilMultiMirnas12mers",
#             as.double(pars),
#             as.double(data),
#             as.double(dil),
#             as.double(l_all),
#             as.double(L),
#             as.double(Y),
#             as.integer(n_i),
#             as.integer(n_j),
#             as.integer(n_mir),
#             as.integer(n_pos),
#             as.integer(zero_grad),
#             as.integer(n_z),
#             grad=double(n_pos*n_i + 2*n_mir))
#   return(out[["grad"]])
# }



# CostRDoubleSite <- function(pars, data_1s, data_2s, dil, l_1s, l_2s, l, L,
#                                  Y, n_i, n_j, upper_, lower_,
#                                  plot_=FALSE, plotname = NULL, tempname=NULL) {

#   kds <- exp(pars[1:n_i])
#   b <- exp(pars[n_i + 1])
#   A <- exp(pars[n_i + 2])
#   a <- FreeAgoR(kds, l, A*c(40, 12.65, 4, 1.265, 0.4)/100)
#   x_1s <- matrix(sapply(a, function(a_j) {
#     a_j*l_1s/(a_j + kds)
#   }), nrow=n_i, byrow=FALSE)
#   kds_d <- kds[-length(kds)]
#   x_2s <- matrix(sapply(a, function(a_j) {
#     l_2s*(1 - c(c(kds_d/(a_j + kds_d))%o%c(kds_d/(a_j + kds_d))))
#   }), nrow=(n_i - 1)^2, byrow=FALSE)
#   x <- rbind(x_1s, x_2s)
#   X <- colSums(x)
#   F <- 100 - X
#   l_all <- c(l_1s, l_2s)
#   xb <- t(t(l_all - x)*b/F)
#   m <- x + xb
#   M <- colSums(m)
#   m_norm <- t(t(m)/M)
#   y <- rbind(matrix(data_1s, nrow=n_i), matrix(data_2s, nrow=(n_i - 1)^2))
#   cost <- -sum(y*log(m_norm))
#   cost
# }

# GradRDoubleSite <- function(pars, data_1s, data_2s, dil, l_1s, l_2s, l, L,
#                                  Y, n_i, n_j, upper_, lower_,
#                                  plot_=FALSE, plotname = NULL, tempname=NULL) {
#   kds <- exp(pars[1:n_i])
#   b <- exp(pars[n_i + 1])
#   A <- exp(pars[n_i + 2])
#   a <- FreeAgoR(kds, l, A*c(40, 12.65, 4, 1.265, 0.4)/100)
#   x_1s <- matrix(sapply(a, function(a_j) {
#     a_j*l_1s/(a_j + kds)
#   }), nrow=n_i, byrow=FALSE)
#   kds_d <- kds[-length(kds)]
#   x_2s <- matrix(sapply(a, function(a_j) {
#     l_2s*(1 - c(c(kds_d/(a_j + kds_d))%o%c(kds_d/(a_j + kds_d))))
#   }), nrow=(n_i - 1)^2, byrow=FALSE)
#   x <- rbind(x_1s, x_2s)
#   X <- colSums(x)
#   F <- 100 - X
#   l_all <- c(l_1s, l_2s)
#   f <- l_all - x
#   xb <- b*t(t(f)/F)
#   m <- x + xb
#   M <- colSums(m)
#   m_norm <- t(t(m)/M)
#   y <- rbind(matrix(data_1s, nrow=n_i), matrix(data_2s, nrow=(n_i - 1)^2))
#   cost <- -sum(y*log(m_norm))
#   dcost_dx <- t(Y/M - t(y/m))
#   C0 <- colSums(sapply(a, function(a_j) {
#     l*kds/(a_j + kds)^2
#   }))

#   kds_d_1 <- rep(kds_d, length(kds_d))
#   kds_d_2 <- rep(kds_d, each=length(kds_d))

#   dc1_da <- sapply(a, function(a_j) {
#     l_1s*kds/(a_j + kds)^2
#   })

#   dc1_dkds <- sapply(a, function(a_j) {
#     -l_1s*a_j/(a_j + kds)^2
#   })

#   dc2_dkds <- sapply(1:length(kds_d), function(ind) {
#     r_inds <- (ind - 1)*length(kds_d) + 1:length(kds_d)

#     l_inds <- c(seq(length(kds_d)) - 1)*length(kds_d) + ind

#     lr_ind <- intersect(l_inds, r_inds)
#     l_inds <- setdiff(l_inds, lr_ind)
#     r_inds <- setdiff(r_inds, lr_ind)
#     print("l")
#     print(sXc_new[[2]][l_inds, ])
#     print("r")
#     print(sXc_new[[2]][r_inds, ])
#     print("lr")
#     print(sXc_new[[2]][lr_ind, , drop=FALSE])

#     dc_l <- sapply(a, function(a_j) {
#       (-l_2s*(a_j)/(a_j + kds_d_1)^2*kds_d_2/(a_j + kds_d_2))[l_inds]
#       })

#     dc_r <- sapply(a, function(a_j) {
#       (-l_2s*(a_j)/(a_j + kds_d_1)*kds_d_1/(a_j + kds_d_2)^2)[r_inds]
#       })

#     dc_lr <- sapply(a, function(a_j) {
#       (-l_2s*(a_j)*kds_d_1/(a_j + kds_d_1)^2)[lr_ind]
#       })


#     print(dc_l)
#     print(dc_r)
#     print(dc_lr)


#     # break
#   })

#   print(length(kds_d))
#   print(nrow(x_2s))
#   # break

#   dc2_da <- sapply(a, function(a_j) {
#     l_2s*kds_d_1*kds_d_2*(kds_d_1 + kds_d_2 + 2*a_j)/(a_j + kds_d_1)^2/(a_j + kds_d_2)^2  
#   })

#   # dcost_db   <- t(Y/M - t(y/m))*xb

#   dc_da <- rbind(dc1_da, dc2_da)
#   dC_da <- colSums(dc_da)

#   da_dA <- A*dil/(1 + C0)

#   dcost_db <- dcost_dx*xb
#   dcost_dA <- t(t(dcost_dx)/F*(t(xb)*dC_da + (F - b)*t(dc_da))*da_dA)
#   dcost_dKd <- t(t(dcost_dx)/F*(t(xb)*dC_da + (F - b)*t(dc_da))*da_dA)


#   grad_out <- rep(0, length(pars))
#   grad_out[n_i + 2] <- sum(dcost_dA)
#   grad_out[n_i + 1] <- sum(dcost_db)
#   print(cost)
#   grad_out
# }

# FreeAgoRDoubleSite <- function(pars, data_1s, data_2s, dil, l_1s, l_2s, l, L,
#                                  Y, n_i, n_j, upper_, lower_,
#                                  plot_=FALSE, plotname = NULL, tempname=NULL) {
#   kds <- exp(pars[1:n_i])
#   b <- exp(pars[n_i + 1])
#   A <- exp(pars[n_i + 2])
#   a <- FreeAgoR(kds, l, A*c(40, 12.65, 4, 1.265, 0.4)/100)
#   return(a)
# }

# X_RDoubleSite <- function(a, pars, data_1s, data_2s, dil, l_1s, l_2s, l, L,
#                                  Y, n_i, n_j, upper_, lower_,
#                                  plot_=FALSE, plotname = NULL, tempname=NULL) {
#   kds <- exp(pars[1:n_i])
#   b <- exp(pars[n_i + 1])
#   A <- exp(pars[n_i + 2])
#   x_1s <- matrix(sapply(a, function(a_j) {
#     a_j*l_1s/(a_j + kds)
#   }), nrow=n_i, byrow=FALSE)
#   kds_d <- kds[-length(kds)]
#   x_2s <- matrix(sapply(a, function(a_j) {
#     l_2s*(1 - c(c(kds_d/(a_j + kds_d))%o%c(kds_d/(a_j + kds_d))))
#   }), nrow=(n_i - 1)^2, byrow=FALSE)
#   x <- rbind(x_1s, x_2s)
#   X <- colSums(x)
#   X
# }

# X1_RDoubleSite <- function(a, pars, data_1s, data_2s, dil, l_1s, l_2s, l, L,
#                                  Y, n_i, n_j, upper_, lower_,
#                                  plot_=FALSE, plotname = NULL, tempname=NULL) {
#   kds <- exp(pars[1:n_i])
#   b <- exp(pars[n_i + 1])
#   A <- exp(pars[n_i + 2])
#   x_1s <- matrix(sapply(a, function(a_j) {
#     a_j*l_1s/(a_j + kds)
#   }), nrow=n_i, byrow=FALSE)
#   kds_d <- kds[-length(kds)]
#   x_2s <- matrix(sapply(a, function(a_j) {
#     l_2s*(1 - c(c(kds_d/(a_j + kds_d))%o%c(kds_d/(a_j + kds_d))))
#   }), nrow=(n_i - 1)^2, byrow=FALSE)
#   X <- colSums(x_1s)
#   X
# }

# X2_RDoubleSite <- function(a, pars, data_1s, data_2s, dil, l_1s, l_2s, l, L,
#                                  Y, n_i, n_j, upper_, lower_,
#                                  plot_=FALSE, plotname = NULL, tempname=NULL) {
#   kds <- exp(pars[1:n_i])
#   b <- exp(pars[n_i + 1])
#   A <- exp(pars[n_i + 2])
#   x_1s <- matrix(sapply(a, function(a_j) {
#     a_j*l_1s/(a_j + kds)
#   }), nrow=n_i, byrow=FALSE)
#   kds_d <- kds[-length(kds)]
#   x_2s <- matrix(sapply(a, function(a_j) {
#     l_2s*(1 - c(c(kds_d/(a_j + kds_d))%o%c(kds_d/(a_j + kds_d))))
#   }), nrow=(n_i - 1)^2, byrow=FALSE)
#   X <- colSums(x_2s)
#   X
# }




# CostCFinalDoubleSite <- function(pars, data_1s, data_2s, dil, l_1s, l_2s, l, L,
#                                  Y, n_i, n_j, upper_, lower_, addNone=FALSE,
#                                  plot_=FALSE, plotname = NULL, tempname=NULL) {
#   if (addNone) {
#     pars <- c(pars[1:(length(pars) - 2)], 0, pars[(length(pars) - 1):length(pars)])
#   }
#   out <- .C("CostEquilDoubleSite", as.double(pars), as.double(data_1s),
#             as.double(data_2s), as.double(dil), as.double(l_1s),
#             as.double(l_2s), as.double(l), as.double(L), as.double(Y),
#             as.integer(n_i), as.integer(n_j), cost=double(1),
#             a_j_check=double(n_j), x_1s_check=double(length(data_1s)),
#             x_2s_check=double(length(data_2s)), F_j_check=double(n_j),
#             b_1s_check=double(length(data_1s)),
#             b_2s_check=double(length(data_2s)),
#             cost_1s_check=double(length(data_1s)),
#             cost_2s_check=double(length(data_2s)),
#             m_1s_norm_check=double(length(data_1s)),
#             m_2s_norm_check=double(length(data_2s)),
#             X_j_check=double(n_j))


#   if (tick%%1000 == 0 && plot_ == TRUE) {
#     # kds_original <- EquilPars("miR-155", experiment="equilibrium", n_constant=5,
#     #                       sitelist="paper")
#     model_1s <- matrix(out[["m_1s_norm_check"]], nrow=n_i)
#     model_2s <- matrix(out[["m_2s_norm_check"]], nrow=(n_i - 1)^2)
#     model_totals <- colSums(model_1s) + colSums(model_2s)
#     model_1s_norm <- t(t(model_1s)/model_totals)
#     model_2s_norm <- t(t(model_2s)/model_totals)
#     model_1s_fitting <<- model_1s
#     model_2s_fitting <<- model_2s
#     model <- c(c(model_1s_norm), c(model_2s_norm))
#     print(tick)
#     print(out[["cost"]])
#     data_norm <- c(data_1s, data_2s)/c(rep(Y, each=n_i), rep(Y, each=(n_i - 1)^2))
#     data_1s_norm <- c(data_1s)/c(rep(Y, each=n_i))
#     # print(exp(pars[n_i + 1]))
#     # print(kds_original["bg_miR-155", ])
#     # print(exp(pars[n_i + 2]))
#     # print(kds_original["AGO_miR-155", ])
#     # model <- 
#     xmin <- 1e-6
#     xmax <- 1
#     ymin <- xmin
#     ymax <- xmax
#     par(kPlotParameters)
#     BlankPlot(log='xy')
#     AddLogAxis(1, "model")
#     AddLogAxis(2, "data")

#     points(c(model_1s_norm), c(data_1s_norm), col=ConvertRColortoRGB("black", alpha=0.4))
#     # abline(0, 1, lty=2)
#   }
#   if (plot_ == TRUE) {
#     tick <<- tick + 1   
#   }

#   out[["cost"]]
# }


# Model_C_DoubleSite <- function(pars, data_1s, data_2s, dil, l_1s, l_2s, l, L,
#                                  Y, n_i, n_j, upper_, lower_,
#                                  plot_=FALSE, plotname = NULL, tempname=NULL) {
#   out <- .C("CostEquilDoubleSite", as.double(pars), as.double(data_1s),
#             as.double(data_2s), as.double(dil), as.double(l_1s),
#             as.double(l_2s), as.double(l), as.double(L), as.double(Y),
#             as.integer(n_i), as.integer(n_j), cost=double(1),
#             a_j_check=double(n_j), x_1s_check=double(length(data_1s)),
#             x_2s_check=double(length(data_2s)), F_j_check=double(n_j),
#             b_1s_check=double(length(data_1s)),
#             b_2s_check=double(length(data_2s)),
#             cost_1s_check=double(length(data_1s)),
#             cost_2s_check=double(length(data_2s)),
#             m_1s_vec=double(length(data_1s)),
#             m_2s_vec=double(length(data_2s)),
#             X_j_check=double(n_j))

#   return(list(matrix(out[["m_1s_vec"]], nrow=n_i),
#               matrix(out[["m_2s_vec"]], nrow=(n_i - 1)^2)))
# }


# GradCNumericalDoubleSite <- function(pars, data_1s, data_2s, dil, l_1s, l_2s, l, L,
#                                  Y, n_i, n_j, upper_, lower_, addNone=FALSE,
#                                  plot_=FALSE, plotname = NULL, tempname=NULL) {
  
#   out <- grad(CostCFinalDoubleSite, pars, data_1s=data_1s, data_2s=data_2s, dil=dil,
#        l_1s=l_1s, l_2s=l_2s, l=l, L=L, Y=Y, n_i=n_i, n_j=n_j, upper_=upper_,
#        lower_=lower_, addNone=addNone, plot_=FALSE, plotname=plotname, tempname=tempname, method="simple")
#   out[n_i] <- 0
#   out
# }



# # GradCFinal <- function(pars, data, dil, l, L, Y, n_i, n_j, plot_=FALSE,
# #                     plotname = NULL, tempname=NULL) {
# #   .C("GradEquilOneMirna", as.double(pars), as.double(data), as.double(dil),
# #      as.double(l), as.double(L), as.double(Y), as.integer(n_i), as.integer(n_j),
# #      grad=double(n_i + 2))[["grad"]]
# # }





# GradCFinalDoubleSite <- function(pars, data_1s, data_2s, dil, l_1s, l_2s, l, L,
#                                  Y, n_i, n_j, upper_, lower_, plot_=FALSE,
#                                  plotname = NULL, tempname=NULL) {
#   out <- .C("GradEquilDoubleSite", as.double(pars), as.double(data_1s), 
#             as.double(data_2s), as.double(dil), as.double(l), as.double(l_1s),
#             as.double(l_2s), as.double(L), as.double(Y), as.integer(n_i),
#             as.integer(n_j), as.integer(n_mir), as.integer(fixed),
#             grad=double(n_i + 2*n_mir))
#   return(out[["grad"]])
# }



# CostC_Global_16mers <- function(pars, sXc, dils, n_x, data_ls, 
#                                 data_rs, l_is, num_exps, L=100, plot_=FALSE,
#                                 plotname=NULL, tempname=NULL) {
#   out <- .C("CostEquilGlobalSingleMirna", as.double(pars),
#             as.double(as.numeric(as.matrix(sXc))), as.double(L),
#             as.double(dils), as.integer(data_ls), as.integer(data_rs),
#             as.integer(l_is),
#             as.integer(n_x), as.integer(num_exps),
#             cost=double(1), kd_check=double(num_exps),
#             A_check=double(num_exps),
#             b_check=double(num_exps), data_check=double(num_exps), l_check=double(num_exps))
#   if (tick%%100 == 0 && plot_ == TRUE) {
#     print(tick)
#     print(tempname)
#     print(out[["cost"]])
#     write.table(pars, file=paste0(tempname, "_global.txt"), sep="\t",
#                 row.names=names(pars), col.names=FALSE, quote=FALSE)
#     print(head(exp(pars[((9*n_x)+1):(10*n_x)])))
#     print(exp(pars[(num_exps*n_x+1):(num_exps*n_x+2)]))

#   }
#   if (plot_ == TRUE) {
#     tick <<- tick + 1   

#   }
#   out[["cost"]]
# }

# CostEquilGlobal16mersMir7_C <- function(pars, sXc, dils, data_ls, 
#                                 data_rs, l_is, n_x, n_exp, L=100, plot_=FALSE,
#                                 plotname=NULL, tempname=NULL) {
#   out <- .C("CostEquilGlobal16mersMir7", as.double(pars),
#             as.double(sXc), as.double(L),
#             as.double(dils), as.integer(data_ls), as.integer(data_rs),
#             as.integer(l_is), as.integer(n_x), as.integer(n_exp),
#             cost=double(3*n_exp), kd_check=double(n_exp),
#             A_check=double(3*n_exp),
#             b_check=double(3*n_exp), data_check=double(3*n_exp), l_check=double(3*n_exp))
#   kds_check_new <<- out[["kd_check"]]
#   data_check_new <<- out[["data_check"]]
#   l_check_new <<- out[["l_check"]]
#   A_check_new <<- out[["A_check"]]
#   b_check_new <<- out[["b_check"]]
#   if (tick%%100 == 0) {
#     print(tick)
#     print(tempname)
#     print(sum(out[["cost"]]))
#     write.table(pars, file=paste0(tempname, "_global.txt"), sep="\t",
#                 row.names=names(pars), col.names=FALSE, quote=FALSE)
#   }
#   tick <<- tick + 1
#   sum(out[["cost"]])
# }



# GradC_Global_16mers <- function(pars, sXc, dils, n_x, data_ls, 
#                                 data_rs, l_is, num_exps, L=100, plot_=FALSE,
#                                 plotname=NULL, tempname=NULL) {
#   num_exps <- 10
#   n_x <- 4^8 + 1
#   out <- .C("GradEquilGlobalSingleMirna", as.double(pars),
#             as.double(as.numeric(as.matrix(sXc))), as.double(L),
#             as.double(dils), as.integer(data_ls), as.integer(data_rs),
#             as.integer(l_is),
#             as.integer(n_x), as.integer(num_exps),
#             dgdp=double(num_exps*n_x + 2), kd_check=double(num_exps),
#             A_check=double(num_exps),
#             b_check=double(num_exps), data_check=double(num_exps), l_check=double(num_exps))[["dgdp"]]
#   if (tick%%10 == 0) {
#     print(tail(out))
#   }
#   out
# }


# GradEquilGlobal16mersMir7_C <- function(pars, sXc, dils, data_ls, 
#                                 data_rs, l_is, n_x, n_exp, L=100, plot_=FALSE,
#                                 plotname=NULL, tempname=NULL) {
#   out <- .C("GradEquilGlobal16mersMir7", as.double(pars), as.double(sXc),
#             as.double(L), as.double(dils), as.integer(data_ls),
#             as.integer(data_rs), as.integer(l_is), as.integer(n_x),
#             as.integer(n_exp),
#             dgdp=double(n_exp*n_x + 6), kd_check=double(n_exp),
#             A_check=double(n_exp),
#             b_check=double(n_exp), data_check=double(n_exp), l_check=double(n_exp))[["dgdp"]]
#   out
# }


# GradC_Global_16mers_Parallel <- function(pars, sXc, dils, n_x, data_ls, 
#                                 data_rs, l_is, num_exps, L=100, plot_=FALSE,
#                                 plotname=NULL, tempname=NULL) {
#   time.start <- proc.time()[3]
#   out <- .C("GradEquilGlobalSingleMirnaParallel", as.double(pars),
#             as.double(as.numeric(as.matrix(sXc))), as.double(L),
#             as.double(dils), as.integer(data_ls), as.integer(data_rs),
#             as.integer(l_is),
#             as.integer(n_x), as.integer(num_exps),
#             dgdp=double(num_exps*n_x + 2), kd_check=double(num_exps),
#             A_check=double(num_exps),
#             b_check=double(num_exps), data_check=double(num_exps), l_check=double(num_exps))
#   print(paste0("\t", proc.time()[3] - time.start))
#   # tick <<- tick + 1
#   # out <- out[["dgdp"]]
#   # names(out) <- names(pars)
#   out
# }



# CostCMir7AllTemp <- function(pars, sXc, dils, n_x, L=100,
#                          data_i_start=c(3, 10, 17, 24, 31, 39, 47, 55, 63, 71, 79, 87),
#                          data_i_end  =c(6, 13, 20, 27, 35, 43, 51, 59, 67, 75, 83, 91),
#                          l_i_p=c(2,9,16,23), plot_=FALSE,
#                          plotname = NULL, tempname=NULL) {
#   # print(sXc[, input_i])
#   mirnas <- c("miR-7-23nt", "miR-7-24nt", "miR-7-25nt")
#   poss <- c("nt2-5", "nt3-6", "nt4-7", "nt5-8")
#   sXc_vec <- unlist(sapply(sXc, function(x) {
#     sapply(x, function(y) as.double(as.numeric(as.matrix(y))))
#   }))
#   print("entering C function:")
#   out <- .C("CostEquilMir7All",
#             as.double(pars), as.double(sXc_vec), as.double(L), as.double(dils),
#             as.integer(data_i_start), as.integer(data_i_end),
#             as.integer(l_i_p), as.integer(n_x),
#             cost=double(3*4), kd_check=double(length(dils)),
#             A_check = double(length(dils)), b_check = double(length(dils)),
#             dils_check = double(length(dils)), data_check=double(length(dils)),
#             l_check = double(length(dils)))
#   print(out[["cost"]])
#   par(kPlotParameters)

#   plot(out[["cost"]], costs_new_all, xlim=c(0, 1.5*max(out[["cost"]])),
#        ylim=c(0, 1.5*max(out[["cost"]])))
#   abline(0, 1)
#   tick <<- tick + 1
#   sum(out[["cost"]])
# }

# CostCMir7All <- function(pars, sXc, dils, n_x, L=100,
#                          #                           23 nt|              24 nt|              25 nt|#
#                          #               1...2...3...4...5|..1...2...3...4...5|..1...2...3...4...5|#
#                          data_i_start=c( 3, 10, 17, 24, 31,  38, 46, 54, 62, 70, 78, 86, 94,102,110),
#                          data_i_end  =c( 6, 13, 20, 27, 34,  42, 50, 58, 66, 74, 82, 90, 98,106,114),
#                          l_i_p=c(2, 9, 16, 23, 30), l_i_p_single = c(1, 8, 15, 22, 29), plot_=FALSE,
#                          plotname = NULL, tempname=NULL, combined=TRUE) {
#   if (!combined) {
#     l_i_p <- l_i_p_single
#   }
#   # print(l_i_p)
#   columns_per_exp <- data_i_end - data_i_start + 1
#   print(columns_per_exp)
#   all_columns <- sum(columns_per_exp)
#   # print(all_columns)
#   out <- .C("CostEquilMir7All", as.double(pars), as.double(sXc), as.double(L),
#             as.double(dils), as.integer(data_i_start), as.integer(data_i_end),
#             as.integer(l_i_p), as.integer(n_x), cost=as.double(0),
#             # cost2=as.double(0),
#             # l_out=double(n_x*5),
#             # data_out=double(n_x*sum(data_i_end - data_i_start + 1)),
#             # data_out_2=double(n_x*sum(data_i_end - data_i_start + 1)),
#             # model_out=double(n_x*sum(data_i_end - data_i_start + 1)),
#             # columns_all=integer(all_columns),
#             # inds_all=integer(n_x*sum(data_i_end - data_i_start + 1)),
#             # i_all=integer(n_x*sum(data_i_end - data_i_start + 1)),
#             # j_all=integer(n_x*sum(data_i_end - data_i_start + 1)),
#             # order_all=integer(n_x*sum(data_i_end - data_i_start + 1)),
#             # Aj_all=double(sum(data_i_end - data_i_start + 1)),
#             # A_m_all=double(sum(data_i_end - data_i_start + 1)),
#             # bg_m_all=double(sum(data_i_end - data_i_start + 1)),
#             cost_each=double(n_x*sum(data_i_end - data_i_start + 1)))
#             # prob_each=double(n_x*sum(data_i_end - data_i_start + 1)),
#             # A_j_sum_alt=double(sum(data_i_end - data_i_start + 1)))
#             # )

#   # l_out <<- matrix(out[["l_out"]], nrow=n_x, ncol=5)
#   # pars_out <<- pars

#   cost_each_2_global <<- out[["cost_each"]]
#   # message("None Kds")
#   # print(pars_out[seq(5)*(4^8+1)])
#   # sXc_vec_2 <<- sXc
#   # model_global <<- out[["model_out"]]
#   # data_global <<- out[["data_out"]]
#   # ind_global <<- out[["inds_all"]]
#   # i_global <<- out[["i_all"]]
#   # j_global <<- out[["j_all"]]
#   # order_global <<- out[["order_all"]]
#   # plot(out[["model_out"]], out[["data_out"]], log='xy')
#   # print(out[["cost"]])
#   # print(out[["cost2"]])
#   # print(out[["Aj_all"]])
#   # print(out[["A_m_all"]])
#   # print(out[["bg_m_all"]])
#   # print(out[["cost_each"]])
#   # prob_each_global <<- out[["prob_each"]]
#   # print(out[["Aj_all"]])
#   # print(out[["A_j_sum_alt"]])
#   # plot(out[["Aj_all"]], out[["A_j_sum_alt"]], log='xy')
#   # abline(0, 1, lty=2)
#   # break

#   # if (tick%%10 == 0) {
#   #   print(tick)
#   #   # row_sample <- c(sample(1:nrow(sXc), 10000, replace=FALSE), n_i)

#   #   plot(out[["model_out"]], out[["data_out"]], log='xy')
#   #   print(out[["cost"]])
#   #   print(out[["cost2"]])
#   #   kds_1_4 <- pars[       1 :   n_x ]
#   #   kds_2_5 <- pars[(  n_x+1):(2*n_x)]
#   #   kds_3_6 <- pars[(2*n_x+1):(3*n_x)]
#   #   kds_4_7 <- pars[(3*n_x+1):(4*n_x)]
#   #   kds_5_8 <- pars[(4*n_x+1):(5*n_x)]
#   #   print(exp(pars[(5*n_x+1):(5*n_x+6)]))
#   #   print(paste0(tempname, "_1-4.txt"))
#   #   write.table(file=paste0(tempname, "_1-4.txt"), kds_1_4)    
#   #   write.table(file=paste0(tempname, "_2-5.txt"), kds_2_5)    
#   #   write.table(file=paste0(tempname, "_3-6.txt"), kds_3_6)    
#   #   write.table(file=paste0(tempname, "_4-7.txt"), kds_4_7)    
#   #   write.table(file=paste0(tempname, "_5-8.txt"), kds_5_8)
#   # }
#   tick <<- tick + 1
#   out[["cost"]]
# }



# GradCMir7All <- function(pars, sXc, dils, n_x, L=100,
#                          #                           23 nt|              24 nt|              25 nt|#
#                          #               1...2...3...4...5|..1...2...3...4...5|..1...2...3...4...5|#
#                          data_i_start=c( 3, 10, 17, 24, 31,  38, 46, 54, 62, 70, 78, 86, 94,102,110),
#                          data_i_end  =c( 6, 13, 20, 27, 34,  42, 50, 58, 66, 74, 82, 90, 98,106,114),
#                          l_i_p=c(2, 9, 16, 23, 30), l_i_p_single = c(1, 8, 15, 22, 29), plot_=FALSE,
#                          plotname = NULL, tempname=NULL, combined=TRUE) {
#   if (!combined) {
#     l_i_p <- l_i_p_single
#   }
#   out <- .C("GradientEquilMir7All", as.double(pars), as.double(sXc),
#             as.double(L), as.double(dils), as.integer(data_i_start),
#             as.integer(data_i_end), as.integer(l_i_p), as.integer(n_x),
#             dgdp=double(5*n_x + 6))
#   out[["dgdp"]]
# }


# CostCNewInteractive <- function(pars, sXc, dil, l, L, Y, n_i, n_j,
#                                 plot_=FALSE, plotname=NULL, tempname=NULL) {
#                     #             pars, y, dil, n_x, L=100, data_i=3:7, input_i=2, plot_=FALSE,
#                     # plotname = NULL) {
#   # print(sXc[, input_i])
#   y <- as.numeric(as.matrix(sXc))
#   out <- .C("CostEquilOneMirna",
#             as.double(pars),
#             as.double(y),
#             as.double(dil),
#             as.double(l),
#             as.double(L),
#             as.double(Y),
#             as.integer(n_i),
#             as.integer(n_j),
#             cost=double(1))[["cost"]]
#   # if (tick%%100 == 0) {
#   #   print(out)
#   # }
#   if (tick%%100 == 0 & plot_) {
#     print(out)
#     print(tick)
#     dev.set(2)
#     model <- ModelC(pars/log(10), sXc)
#     row_sample <- c(sample(1:nrow(sXc), 10000, replace=FALSE), n_i)
#     # row_sample <- 1:nrow(sXc)
#     plot(unlist(sXc[row_sample, seq(n_j) + 2]), c(model[row_sample, ]), log='xy',
#          xlim=c(0.1, 3e7), ylim=c(0.1, 3e7),
#          col=rep(
#                  c("blue", "green", "orange", "red", "purple"),
#                  each=length(row_sample)))
#     legend("bottomright", legend=colnames(sXc)[seq(n_j) + 2], lwd=2,
#            col=c("blue", "green", "orange", "red", "purple"))
#     if (class(tempname) == "character") {
#       print(tempname)
#       write.table(file=paste0(tempname), pars, quote=FALSE, col.names=FALSE)
#       print(head(fread(file=paste0(tempname))))
#     }

#     print(exp(pars["AGO"]))
#     print(exp(pars["bg"]))
#     print(proc.time()[3] - time_start)
#   }
#   tick <<- tick + 1
#   out
# }

# GradCNew <- function(pars, sXc, dil, l, L, Y, n_i, n_j, plot_=FALSE,
#                     plotname=NULL, tempname=NULL) {
#   y <- as.numeric(as.matrix(sXc))
#   .C("GradEquilOneMirna",
#      as.double(pars),
#      as.double(y),
#      as.double(dil),
#      as.double(l),
#      as.double(L),
#      as.double(Y),
#      as.integer(n_i),
#      as.integer(n_j),
#      grad=double(n_i + 2))[["grad"]]
# }




# ModelC <- function(pars, sXc, L=100, combined=TRUE) {
#   data <- GetDataEquil(sXc)
#   l <- GetInputEquil(sXc=sXc, combined=combined)
#   dil <- sapply(colnames(data), as.numeric)/100
#   kds <- 10^pars[1:length(l)]
#   As <- 10^pars["AGO"]*dil
#   bg <- 10^pars["bg"]
#   as <- sapply(As, FreeAgo, kds=kds, l=l)
#   n_x <- nrow(sXc)
#   n_col <- ncol(data)
#   matrix(.C("ModelEquil", as.double(kds), as.double(l), as.double(L),
#              as.double(bg), as.double(dil), as.double(As), as.double(as),
#              as.double(colSums(data)), as.integer(n_x), as.integer(n_col),
#              model=double(n_x*n_col))[["model"]],
#         ncol=n_col, nrow=n_x, byrow=FALSE, dimnames=list(rownames(data),
#                                                          colnames(data)))
# }


# # CostC <- function(pars, sXc, L=100, combined=TRUE, plot_=FALSE, plotname=NULL) {
# #   data <- GetDataEquil(sXc)
# #   l <- GetInputEquil(sXc=sXc, combined=combined)
# #   dil <- sapply(colnames(data), as.numeric)/100
# #   kds <- 10^pars[1:length(l)]
# #   As <- 10^pars["AGO"]*dil
# #   bg <- 10^pars["bg"]
# #   as <- sapply(As, FreeAgoC, kds=kds, l=l)
# #   n_x <- nrow(sXc)
# #   n_col <- ncol(data)
# #   tick <<- tick + 1
# #   out <- .C("CostEquil", as.double(kds), as.double(l), as.double(L),
# #              as.double(bg), as.double(dil), as.double(As), as.double(as),
# #              as.double(as.numeric(as.matrix(data))), as.integer(n_x),
# #              as.integer(n_col), cost = double(1))[["cost"]]
# #   if (tick%%100 == 0) {
# #     print(out)
# #   }  
# #   if (tick%%1000 == 0 & plot_) {
# #     print(out)
# #     print(tick)
# #     pdf(file=plotname)
# #     model <- ModelC(pars, sXc)
# #     plot(data, model, log='xy', xlim=c(0.1, 1e7), ylim=c(0.1, 1e7))
# #     dev.off()
# #     print(10^pars["AGO"])
# #     print(10^pars["bg"])
# #   }
# #   out
# # }

# CostCInteractive <- function(pars, sXc, L=100) {
#   time_start <- proc.time()[3]
#   data <- GetDataEquil(sXc)
#   as <- sapply(As, FreeAgoC, kds=kds, l=l)
#   n_x <- nrow(sXc)
#   out <- SubfuntionCall(CostC)
#   if (tick%%100 == 0) {
#     l <- GetInputEquil(sXc=sXc)
#     dil <- sapply(colnames(data), as.numeric)/100
#     kds <- 10^pars[1:length(l)]
#     As <- 10^pars["AGO"]*dil
#     bg <- 10^pars["bg"]
#     dev.set(2)
#     m <- sapply(1:length(as), function(i) {
#       as[i]*l/(as[i] + kds) + bg*kds*l/(as[i] + kds)/(L - As[i] + as[i])
#       })
#     sample_points <- sample(n_x, 10000)
#     plot(c(MatNorm(m)[sample_points,]),
#           c(MatNorm(data)[sample_points,]), xlim=c(1e-8, 1), ylim=c(1e-8, 1),
#            log='xy', col=ConvertRColortoRGB("black", alpha=0.2))
#   }
#   print(out[["cost"]])

#    tick <<- tick + 1
#   out[["cost"]]
# }

# GetInputEquil <- function(sXc, combined=TRUE, conc=100) {
#   if (combined) ind <- 2
#   else          ind <- 1
#   s.c <- sXc[, ind]
#   s.c[s.c==0] <- 1
#   l <- Norm(s.c)*conc
#   names(l) <- rownames(sXc)
#   l
# }

# # GetDataEquil <- function(sXc) {
# #   data.matrix(sXc[, 3:(ncol(sXc) - 1)])
# # }


# # Overall cost function:
# CostOld <- function(pars, sXc, combined=TRUE,
#                                     cols.final=NULL, zero.vector=NULL,
#                                     print.=TRUE) {
#   # print("getting l")
#   l <- SubfunctionCall(GetInputEquil)
#   # print("getting data")
#   data <- SubfunctionCall(GetDataEquil)
#   if (is.null(cols.final))
#     cols.final <- 1:ncol(data)
#   dil <- sapply(colnames(data), as.numeric)/100
#   sXa <- SubfunctionCall(sXa_of_pars.l.dil)
#   sXm <- SubfunctionCall(sXm_of_pars.sXa)
#   # if (print.) {
#     # R.data <- EquilEnrichments(data, l)
#     # R.model <- EquilEnrichments(sXm, l)
#     # pdf(file="current_optimization.pdf")
#     # print("opened pdf")
#     # plot(MatNorm(data), MatNorm(sXm), log='xy')
#     # dev.off()
#     # par(kPlotParameters)
#     # plot(c(1), type = "n", xlim=c(0.09,60), ylim=c(0.3, 300), log='xy')
#     # sapply(1:nrow(R.data), function(i) {
#     #   x <- colnames(R.data)
#     #   x[length(x)] <- 0.1
#     #   points(x, R.data[i, ], col=kSiteColors[rownames(R.data)][i])
#     #   lines(x, R.model[i, ], col=kSiteColors[rownames(R.data)][i])
#     # })
#     # # print(pars)
#     # grad.n <- grad(CostFunctionEquilibrium, pars, sXc=sXc, print=FALSE)
#     # grad.s <- grad(CostFunctionEquilibrium, pars, sXc=sXc, print=FALSE, method="simple")
#     # grad.a <- GradientEquilibrium(pars, sXc)
#     # dev.set(4)
#     # kSiteColors["bg"] <- "brown"
#     # plot(c(grad.a, grad.a), c(grad.n, grad.s), pch=rep(c(19, 1), each=length(pars)), col=kSiteColors[c(rownames(R.data), "bg", "AGO")])
#     # abline(0, 1)
#     # print(SubfunctionCall(L_of_sXm))
#   # }
#   tick <<- tick + 1
#   if (print.)
#     print(SubfunctionCall(L_of_sXm))
#   # print(SubfunctionCall(L_of_sXm))
#   return(SubfunctionCall(L_of_sXm))
# }



# GradientOld <- function(pars, sXc, combined=TRUE, print.=TRUE,
#                                 zero.vector=1, cols.final=NULL) {
#   # print("in gradient")
#   # print("getting l")
#   l <- SubfunctionCall(GetInputEquil)
#   # print("getting data")
#   data <- SubfunctionCall(GetDataEquil)
#   if (is.null(cols.final))
#     cols.final <- 1:ncol(data)
#   dil <- sapply(colnames(data), as.numeric)/100
#   # print(dil)
#   sXa <- SubfunctionCall(sXa_of_pars.l.dil)
#   sXm <- SubfunctionCall(sXm_of_pars.sXa)
#   dL_dsXm <- SubfunctionCall(DL_DsXm)
#   dL_dlogp <- rowSums(sapply(cols.final, function(column) {
#     ds.a_dlogp <- SubfunctionCall(Ds.a_Dlogp)
#     ds.m_ds.a  <- SubfunctionCall(Ds.m_Ds.a)
#     ds.m_dlogp <- SubfunctionCall(Ds.m_Dlogp)
#     #    vector          %matrix
#     dL_dsXm[, column]%*%(ds.m_dlogp + ds.m_ds.a%*%ds.a_dlogp)
#   }))
#   out <- log(10)*10^pars*dL_dlogp
#   out*zero.vector
# }


BgRNA <- function(x, l, bg) {
  bg*Norm(l - x)  
}


# EquilBoundTotal <- function(pars, sXc, A.dil=NULL, combined=TRUE,
#                                  addbg=TRUE, flowthrough=FALSE) {
#   # Assign the kds and stock Ago concentration:
#   sXc <- sXc[,colnames(sXc)!="seq"]
#   names(pars) <- gsub("(.*)_Kd", names(pars), replace="\\1")
#   names(pars) <- gsub("(.*)_Kd", names(pars), replace="\\1")
#   kds     <- 10^pars[rownames(sXc)]
#   A.stock <- 10^pars["AGO"]
#   # Assign background term based on flag:
#   # Assign the input RNA concentration based on flag:
#   l <- 100*Norm(sXc[, 1 + combined])
#   # Assign Ago concentrations:
#   if (length(A.dil) == 0) A.dil <- as.numeric(colnames(sXc[, 3:7]))/100
#   else                    A.dil <- A.dil/100
#   A.vec <- A.stock*A.dil
#   # Assign the bound RNA, x:
#   x <- sapply(A.vec, BoundRNA, kds=kds, l=l)
#   x_bound <<- x
#   # Add RNA:
#   if (addbg) {
#     bg <- 10^pars["bg"]
#     x.bg <- apply(x, 2, BgRNA, l=l, bg=bg)
#     x <- x + x.bg
#   }
#   x
# }



EquilSingleSiteModelFreq <- function(pars, sXc, A.dil=NULL, combined=TRUE,
                                 addbg=TRUE, flowthrough=FALSE) {
  # Assign the kds and stock Ago concentration:
  # print("combined:")
  # print(combined)
  # print(pars)
  # print(sXc)
  # break
  names(pars) <- gsub("(.*)_Kd", names(pars), replace="\\1")
  names(pars) <- gsub("(.*)_Kd", names(pars), replace="\\1")
  kds     <- 10^pars[rownames(sXc)]
  A.stock <- 10^pars["AGO"]
  # Assign background term based on flag:
  # Assign the input RNA concentration based on flag:
  l <- 100*Norm(sXc[, 1 + combined])
  # If this the figure scripts never break, get rid of the `l_global` variable.
  # l_global <<- l
  # Assign Ago concentrations:
  if (length(A.dil) == 0) A.dil <- as.numeric(colnames(sXc[, 3:7]))/100
  else                    A.dil <- A.dil/100
  A.vec <- A.stock*A.dil
  A.vec_global <<- A.vec
  # Assign the bound RNA, x:
  x <- sapply(A.vec, BoundRNA, kds=kds, l=l)
  x_bound <<- x
  # Add RNA:
  if (addbg) {
    bg <- 10^pars["bg"]
    x.bg <- apply(x, 2, BgRNA, l=l, bg=bg)
    x <- x + x.bg
  }
  x_global <<- x
  if (flowthrough) {
    x <- apply(x, 2, function(x_i) {c(l) - x_i})
  }
  apply(x, 2, Norm)      
}

# EquilSingleSiteModelPerSite <- function(pars, sXc, A.dil=NULL, combined=TRUE,
#                                  addbg=TRUE) {
#   # Assign the kds and stock Ago concentration:
#   sXc <- sXc[,colnames(sXc)!="seq"]
#   names(pars) <- gsub("(.*)_Kd", names(pars), replace="\\1")
#   names(pars) <- gsub("(.*)_Kd", names(pars), replace="\\1")

#   kds     <- 10^pars[rownames(sXc)]
#   A.stock <- 10^pars["AGO"]
#   # Assign background term based on flag:
#   # Assign the input RNA concentration based on flag:
#   l <- 100*Norm(sXc[, 1 + combined])
#   l_global <<- l
#   # Assign Ago concentrations:
#   if (length(A.dil) == 0) A.dil <- as.numeric(colnames(sXc[, 3:7]))/100
#   else                    A.dil <- A.dil/100
#   A.vec <- A.stock*A.dil
#   # Assign the bound RNA, x:
#   x <- sapply(A.vec, BoundRNA, kds=kds, l=l)
#   # Add RNA:
#   if (addbg) {
#     bg <- 10^pars["bg"]
#     x.bg <- apply(x, 2, BgRNA, l=l, bg=bg)
#     x <- x + x.bg
#   }
#   apply(x, 2, function(x_i) {
#     x_i / l
#   })      
# }


# AddPseudoCount <- function(model.freq, sXc) {
#   sXc <- sXc[,colnames(sXc)!="seq"]
#   pcounts <- 1/colSums(sXc[, 3:7])
#   model.ps <- t(t(model.freq) + pcounts)
#   apply(model.ps, 2, Norm)
# }

# CostMulti <- function(model.freq, sXc) {
#   data <- sXc[,colnames(sXc)!="seq"][,3:7]
#   -sum(data*log(model.freq))
# }




# DCostMultiDModel <- function(model.p, data) {
#   -data/model.p
# }

# DerivPseudo <- function(model.freq, sXc) {
#   sXc <- sXc[, colnames(sXc)!="seq"]
#   pcounts <- 1/colSums(sXc[, 3:7])
   
# }





# FullOptFunctionEquil <- function(pars, sXc, printcost=FALSE) {
#   model.freq <- SubfunctionCall(EquilSingleSiteModelFreq)
#   cost <- SubfunctionCall(CostMulti)
#   # if (printcost) print(cost)
#   data <- sXc[,colnames(sXc)!="seq"][,3:7]
#   data <- apply(data, 2, Norm)
#   if (tick %% 10 == 0) {
#     plot(c(as.numeric(model.freq)), as.numeric(unlist(data)), log ='xy')
#     print(pars)
#   }

#   tick <<- tick + 1

#   cost
# }


# MultinomialCost2 <- function(model.p, data) {
#   -sum(data*log(model.p)) + sum(data*log(apply(data, 2, Norm)))
# }




# ModelFunctionNew <- function(pars, num.kds, num.bgs) {
#   pars.current <<- pars
#   kds  <- 10^pars[1 : num.kds]
#   bgs  <- rep(10^pars[(num.kds + 1)], ncol(data))
#   stock.ago <- 10^pars[num.kds + 2]

#   kds.multi <- lapply(rownames(data.m), function(name) {
#   names <- unlist(strsplit(name, split = ","))
#   kd <- kds[names]
#   return(kd)
#   })
#   print(kds.multi)
#   # Get the bound counts in each exp:
#   c.agos <- sapply(colnames(data), function(x) {
#     as.numeric(x) / 100
#   })
#   c.bounds <- as.matrix(
#     sapply(c.agos, function(percent) {
#       return(GetBoundRNA(kds, c.I.tots, percent * stock.ago))
#     }
#   ))

#     c.boundsmulti <- as.matrix(
#     sapply(c.agos, function(percent) {
#       return(GetBoundRNAMulti(kds.multi, c.I.tots.multi, percent * stock.ago))
#     }
#   ))
#   c.bounds <<- c.bounds
#   c.boundsmulti <<- c.boundsmulti
#   # Get the amount of background binding by subtracting the bound from the
#   # total sites in each exp, normalizing. Must transpose to multiply
#   # each column.
#   c.frees <- c.totals - c.boundsmulti
#   c.bgs <- t(t(c.frees) * bgs / colSums(c.frees))
#   c.all <- c.boundsmulti + c.bgs
#   c.final <- data.frame(t(t(c.all) / colSums(c.all) * colSums(data.m)))
#   return(c.final)
# }


# ## Save, for working out multisite binding:
# # kds.single <- c(0.5, 0.5, 0.5)
# # kds.multi <- list(c(0.5, 0.5), c(0.5))

# # l.single <- c(1, 1, 1)
# # l.multi  <- c(0.5, 1)

# # a.single <- FreeAgo(kds.single, l.single, 2)
# # a.multi <- FreeAgo(kds.multi, l.multi, 2)
# # a.multi2 <- FreeAgo(kds.multi, l.multi, 2, tol=.Machine$double.eps^0.25)
# # a.multi3 <- FreeAgo(kds.multi, l.multi, 2, tol=1e-1*.Machine$double.eps^0.25)
# # a.multi4 <- FreeAgo(kds.multi, l.multi, 2, tol=1e-2*.Machine$double.eps^0.25)
# # a.multi5 <- FreeAgo(kds.multi, l.multi, 2, tol=.Machine$double.eps)

# # ocs <- SiteOcc(a.single, kds.single)


# # x.single <- BoundRNA(kds.single, l.single, 2)
# # x.multi <- BoundRNA(kds.multi, l.multi, 2)


# #23456789!123456789@123456789#123456789$123456789%123456789^123456789&123456789*
# GetContamResidual <- function(b, kd.b, l, A, a, B, exponent=1){
#   # a: [free AGO] in the binding reaction
#   # b: [free contaminant] in the binding reaction
#   # A: [total AGO] in the binding reaction
#   # B: [total contaminant] in the binding reaction
#   # kd.b: The Kd value for the contaminant for all RNA.
#   res <- (b^2 + (sum(l) - A - B + a + kd.b)*b - B*kd.b)^exponent
#   return(res)
# }

# #23456789!123456789@123456789#123456789$123456789%123456789^123456789&123456789*
# GetFreeContam <- function(kds.a, kd.b, l, A, a, B) {
#   solution <- NaN
#   try(solution <- uniroot(GetContamResidual, c(0, B), kd.b=kd.b, l=l, A=A, B=B,
#                           a=a, tol=1e-6*.Machine$double.eps^0.25)$root,
#       silent=TRUE)
#   if (is.na(solution)) {
#     solution <- optimize(GetContamResidual, c(0, B), kd.b=kd.b, l=l, A=A, B=B,
#                          a=a, exponent=2, tol=1e-6*.Machine$double.eps^0.25)$minimum
#   }
#   return(solution)
# }

# #23456789!123456789@123456789#123456789$123456789%123456789^123456789&123456789*
# FreeAgoAndContamCombinedRoot <- function(a, kds.a, kd.b, l, A, B){
#   # a: The concentration of free AGO in the binding reaction
#   # b: The concentration of free contaminant in the binding reaction
#   # A: The concentration of free AGO in the binding reaction
#   # B: The concentration of total contaminant in the binding reaction
#   # kds.a: The Kd values specific to Ago binding each target site.
#   # kd.b: The Kd value for the contaminant for all RNA.
#   b <- GetFreeContam(kds.a, kd.b, l, A, a, B)
#   ocs <- GetOccupancy(a, kds.a*(b/kd.b + 1))
#   root <- A - a - sum(ocs*l)
#   return(root)
# }

# #23456789!123456789@123456789#123456789$123456789%123456789^123456789&123456789*
# GetFreeAgoAndContam <- function(kds.a, kd.b, l, A, B) {
#   # a: The concentration of free AGO in the binding reaction
#   # b: The concentration of free contaminant in the binding reaction
#   # A: The concentration of free AGO in the binding reaction
#   # B: The concentration of total contaminant in the binding reaction
#   # kds.a: The Kd values specific to Ago binding each target site.
#   # kd.b: The Kd value for the contaminant for all RNA.
#   a <- uniroot(FreeAgoAndContamCombinedRoot, c(0, A), kds.a=kds.a, kd.b=kd.b,
#                l=l, A=A, B=B, tol=1e-6*.Machine$double.eps^0.25)$root
#   b <- GetFreeContam(kds.a, kd.b, l, A, a, B)
#   return(c(a, b))
# }

# #23456789!123456789@123456789#123456789$123456789%123456789^123456789&123456789*
# GetOccupanciesContaminant <- function(a, b, kds.a, kd.b) {
#   # Calculate common denomenator:
#   ocs.den <- kds.a*kd.b + kds.a*b + kd.b*a
#   return(list(a*kd.b/ocs.den, b*kds.a/ocs.den))
# }

# #23456789!123456789@123456789#123456789$123456789%123456789^123456789&123456789*




