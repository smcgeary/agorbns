FitOffsetAndPairingModel <- function(
  mirna, experiment, n_constant=0, offsetmin=0, offsetmax=16,
  sitelist="programmed_suppcomp", collapsemm=FALSE, supp_base=FALSE
) {
  # Make a dataframe concatenating all the offset matrices into one.
  df_all <- do.call("rbind", lapply(offsetmin:offsetmax,
                                MakeNucleotideContributionDf,
                                mirna=mirna,
                                experiment=experiment,
                                n_constant=n_constant,
                                collapsemm=collapsemm,
                                supp_base=supp_base))
  print("made df")
  # Remove all rows for which the kd is NA
  # Fit the linear model on the resuting data set.
  model <- lm(logkd ~ pos_5p3p + offset - 1, data=df_all)
  model <<- model
  coefs <- model$coefficients
  # Get the offset coefficients:
  coefs_offset <- coefs[sprintf("offset%s", (offsetmin + 1):offsetmax)]
  coefs_pairing <- coefs[grep("pos_5p3p", names(coefs), value=TRUE)]
  coefs_base <- coefs[1]
  names(coefs_pairing) <- gsub("pos_5p3p", replacement="", names(coefs_pairing))
  len_mir <- nchar(kMirnaSeqs[mirna])
  # Define the possible 5-prime starting nucleotides and possible 3-prime
  # starting nucleotides, for overall matrix.
  # if (downto2mers) {
  #   len_k <- 2
  # } else if (downto3mers) {
  #   len_k <- 3
  # } else if (downto4mers) {
  len_k <- 4
  # } else {
  #   len_k <- 5
  # }
  nucs_5p <- 9:(len_mir - len_k + 1)
  nucs_3p <- (9 + len_k - 1):len_mir

  # Define the output matrix.
  out_matrix <- matrix(NA, nrow=length(nucs_3p), ncol=length(nucs_5p),
                       dimnames=list(nucs_3p, nucs_5p))
  for (i in 1:length(coefs_pairing)) {
    coef_pairing <- coefs_pairing[i]
    name_pairing <- names(coefs_pairing)[i]
    name_split <- unlist(strsplit(name_pairing, split="\\|"))
    rowname_val <- name_split[2]
    colname_val <- name_split[1]
    out_matrix[rowname_val, colname_val] <- coef_pairing
  }
  # out_matrix["14", "10"] <- 0
  return(list(`offsets`=coefs_offset, `pairing`=out_matrix))
}

DetermineKdCorrectionParameters <- function(kds_p, kds_r) {
  kas_p_fit <- 1/kds_p[1:24, 2]
  kas_r_fit <- 1/kds_r[rownames(kds_p)[1:24], 2]
  min_val <- min(1/kds_p[1:(nrow(kds_p) - 2), 2])
  print(min_val)
  # print(kas_p_fit)
  # print(kas_r_fit)
  pairs <- cbind(kas_p_fit, kas_r_fit)
  rownames(pairs) <- rownames(kds_p)[1:24]
  # print(pairs)
  # dev.new(xpos=20, ypos=20, height=5, width=5)
  # plot(kas_r_fit, kas_p_fit)
  # 3.B Define the cost function.
  # dev.new(xpos=20, ypos=20, height=5, width=5)
  # plot(kas_r_fit, kas_p_fit, log="xy", xlim=c(0.1, 1000), ylim=c(0.1, 1000))
  CostFunction <- function(pars) {
    # Use the exponential transormation because this doesn't cause warning
    # messages from R, owing to negative numbers, which makes the log
    # transformation fundamentally impossible.
    m <- exp(pars[1])
    b <- pars[2]
    # b <- 0
    # This is the new equation derived on 200707, that 1.) with increasing Ka,
    # it approaches Ka * m, but causes Ka of the worst site in the random
    # library to set to 1. The function is log transformed for the comparison.
    # model <- log(m*(kas_r_fit - min(kas_r_fit)) + 1)
    # model <- log(m*(kas_r_fit - 1) + 1)
    inds <- which(kas_r_fit > b)
    model <- kas_r_fit
    model[inds] <- m*(kas_r_fit[inds] - b) + b
    # model <- log(m*(kas_r_fit - b) + b)
    # lines(kas_r_fit, exp(model), lty=2)
    # Calculate the squared residuals.
    res <- (log(kas_p_fit) - log(model))^2
    # Return the sum of the squared residuals.
    return(sum(res))
  }
  # Perform the optimization to get parameter m and b.
  pars_init <- c(0, (1 + min(kas_r_fit))/2)
  pars_lower <- c(log(1/100), 1)
  pars_upper <- c(log(100), min(kas_r_fit))
  # pars_init <- c(0)
  # pars_lower <- c(log(1/100))
  # pars_upper <- c(log(100))
  print(pars_init)
  print(pars_lower)
  print(pars_upper)
  pars <- optim(pars_init, CostFunction, method="L-BFGS-B",
                lower=pars_lower, upper=pars_upper)$par
  m <- exp(pars[1])
  b <- pars[2]
  # b <- 0
  x <- 10^seq(0, 4, length.out=50)
  inds <- which(x > b)
  y <- x
  y[inds] <- m*(y[inds] - b) + b


  # lines(x, y, lty=2)

  return(list(`m`=m, `b`=b))
}

FitTriplePairingAndOffsetModel <- function(
  n_constant=3, offsetmin=-4, offsetmax=16,
  sitelist="programmed_suppcomp", collapsemm=FALSE, supp_base=FALSE
) {
  # Make a dataframe concatenating all the offset matrices into one.
  df_all_l7 <- do.call(
    "rbind",
    lapply(offsetmin:offsetmax, MakeNucleotideContributionDf, mirna="let-7a-21nt",
           experiment="equil_c2_nb", n_constant=n_constant, collapsemm=collapsemm,
           supp_base=supp_base, nucletters=TRUE)
  )

  df_all_m1 <- do.call(
    "rbind",
    lapply(offsetmin:offsetmax, MakeNucleotideContributionDf, mirna="miR-1",
           experiment="equil_c_nb", n_constant=n_constant, collapsemm=collapsemm,
           supp_base=supp_base, nucletters=TRUE)
  )

  df_all_m155 <- do.call(
    "rbind",
    lapply(offsetmin:offsetmax, MakeNucleotideContributionDf, mirna="miR-155",
           experiment="equil_sc_nb", n_constant=n_constant, collapsemm=collapsemm,
           supp_base=supp_base, nucletters=TRUE)
  )


  # df_all_l7$pos_5p3p <- paste0("l7|", df_all_l7$pos_5p3p)
  # df_all_m1$pos_5p3p <- paste0("m1|", df_all_m1$pos_5p3p)
  # df_all_m155$pos_5p3p <- paste0("m155|", df_all_m155$pos_5p3p)

  df_all <- do.call("rbind", list(df_all_l7, df_all_m1, df_all_m155))


  df_all$mir <- rep(c("l7", "m1", "m155"),
                    times=c(nrow(df_all_l7),
                            nrow(df_all_m1),
                            nrow(df_all_m155)))


  print(head(df_all))
  print(tail(df_all))
  break
  # Remove any levels that are not within the dataframe    
  df_all <- drop.levels(df_all)
  n_pars_5p3p <- length(unique(df_all$pos_5p3p))

  pars_5p3pl7 <- unique(grep("^l7\\|", df_all$pos_5p3p, perl=TRUE, value=TRUE))
  pars_5p3pm1 <- unique(grep("^m1\\|", df_all$pos_5p3p, perl=TRUE, value=TRUE))
  pars_5p3pm155 <- unique(grep("^m155\\|", df_all$pos_5p3p, perl=TRUE, value=TRUE))

  n_pars_5p3pl7   <- 2*length(pars_5p3pl7) 
  n_pars_5p3pm1   <- 2*length(pars_5p3pm1)
  n_pars_5p3pm155 <- 3*length(pars_5p3pm155)

  n_pars_offset <- 2*length(unique(df_all$offset))

  print(n_pars_5p3pl7)
  print(n_pars_5p3pm1)
  print(n_pars_5p3pm155)
  print(n_pars_offset)
  n_pars_5p3p <- n_pars_5p3pl7 + n_pars_5p3pm1 + n_pars_5p3pm155

  pars.init <- rep(0, n_pars_5p3p + n_pars_offset + 1)
  pars.init[(n_pars_5p3p + 1):(n_pars_5p3p + 5)] <- 5
  names(pars.init) <- c(paste0(pars_5p3pl7, "_1"),
                        paste0(pars_5p3pl7, "_2"),
                        paste0(pars_5p3pm1, "_1"),
                        paste0(pars_5p3pm1, "_2"),
                        paste0(pars_5p3pm155, "_1"),
                        paste0(pars_5p3pm155, "_2"),
                        paste0(pars_5p3pm155, "_3"),
                        paste0(unique(as.character(df_all$offset)), "_1"),
                        paste0(unique(as.character(df_all$offset)), "_2"),
                        "base")

  coefs_1 <- -log(FitPairingAndOffsetModel("let-7a-21nt", "equil_c2_nb")[[2]])
  coefs_2 <- -log(FitPairingAndOffsetModel("miR-1", "equil_c_nb")[[2]])

  pairing_l7 <- FitPairingAndOffsetModel("let-7a-21nt", "equil_c2_nb")[[1]]
  print(pairing_l7)
  print(as.numeric(pairing_l7))
  names_pairing_l7 <- paste0(rep(colnames(pairing_l7), each=nrow(pairing_l7)),
                             rep(rownames(pairing_l7), times=ncol(pairing_l7)), sep="|")
  print(names_pairing_l7)
  print(cbind(names(pairing_l7), as.numeric(pairing_l7)))
  break

  pars.init[(n_pars_5p3p + 1): (n_pars_5p3p + n_pars_offset/2)] <- coefs_1
  pars.init[(n_pars_5p3p + n_pars_offset/2 + 1):(n_pars_5p3p + n_pars_offset)] <- coefs_2


  print(pars.init)
  # break

  # n_pars_5p3p <- length(unique(grep("^m1\\|", df_all$pos_5p3p, perl=TRUE, value=TRUE)))

  # n_pars_offset <- length(unique(df_all$offset))
  # Define the loss and gradient functions.
  graphics.off()
  dev.new(xpos= 20, ypos= 20, height=3, width=3)
  dev.new(xpos=320, ypos= 20, height=3, width=3)
  dev.new(xpos=620, ypos= 20, height=3, width=3)
  dev.new(xpos= 20, ypos=320, height=3, width=3)
  dev.new(xpos=320, ypos=320, height=3, width=3)
  dev.new(xpos=620, ypos=320, height=3, width=3)

  LossFunction <- function(pars) {
    # Initialize the parameters.
    pars.5p3pl7 <- pars[1:n_pars_5p3pl7]
    pars.5p3pm1 <- pars[(n_pars_5p3pl7 + 1):(n_pars_5p3pl7 + n_pars_5p3pm1)]
    pars.5p3pm155 <- pars[(n_pars_5p3pl7 + n_pars_5p3pm1 + 1):(n_pars_5p3pl7 + n_pars_5p3pm1 + n_pars_5p3pm155)]
    pars.offset <- pars[(n_pars_5p3p + 1):(n_pars_5p3p + n_pars_offset)]
    pars.base <- exp(pars[length(pars)])
    # print(pars.5p3pl7)
    # print(pars.5p3pm1)
    # print(pars.5p3pm155)
    # print(pars.offset)
    # print(pars.base)
    n_l7 <- length(pars.5p3pl7)
    n_m1 <- length(pars.5p3pm1)
    n_m155 <- length(pars.5p3pm155)
    #Binding model 1
    pars.5p3p1 <- c(pars.5p3pl7[1:(n_l7/2)], pars.5p3pm1[1:(n_m1/2)], pars.5p3pm155[1:(n_m155/3)])
    pars.5p3p2 <- c(pars.5p3pl7[(n_l7/2 + 1):(n_l7)], pars.5p3pm1[(n_m1/2 + 1):(n_m1)], pars.5p3pm155[(n_m155/3 + 1):(n_m155*2/3)])
    pars.5p3p3 <- c(rep(0, n_l7/2), rep(0, n_m1/2), pars.5p3pm155[(n_m155*2/3 + 1):(n_m155)])

    pars.offset1 <- pars.offset[1:(length(pars.offset)/2)]
    pars.offset2 <- pars.offset[(length(pars.offset)/2 + 1):length(pars.offset)]
    names(pars.offset1) <- gsub("^(.*)_.*$", replacement="\\1", names(pars.offset1))
    names(pars.offset2) <- names(pars.offset1)
    names(pars.5p3p1) <- gsub("^(.*)_.*$", replacement="\\1", names(pars.5p3p1),
                              perl=TRUE)
    names(pars.5p3p2) <- names(pars.5p3p1)
    names(pars.5p3p3) <- names(pars.5p3p1)
    # Calculate the dG and offset values for each datapoint.
    dG_5p3p1 <- exp(pars.5p3p1[as.character(df_all$pos_5p3p)])
    dG_5p3p2 <- exp(pars.5p3p2[as.character(df_all$pos_5p3p)])
    dG_5p3p3 <- exp(pars.5p3p3[as.character(df_all$pos_5p3p)])
    dG_5p3p3[1:(nrow(df_all_l7) + nrow(df_all_m1))] <- 0

    sigmoids1 <- 1/(exp(pars.offset1[as.character(df_all$offset)]) + 1)
    sigmoids2 <- 1/(exp(pars.offset2[as.character(df_all$offset)]) + 1)
    # Calculate the final model.
    model <- dG_5p3p1*sigmoids1 + dG_5p3p2*sigmoids2 + dG_5p3p3 + pars.base
    if (tick%%10000 == 0) {
      message("pars.offset1:")
      print(pars.offset1)
      message("pars.offset2:")
      print(pars.offset2)
      message("pars.offset")
      print(pars.offset)
      print(pars[(n_pars_5p3p + 1):(n_pars_5p3p + n_pars_offset)])
      cols <- rep(c("red", "blue", "green"), times=c(nrow(df_all_l7),
                                                   nrow(df_all_m1), 
                                                   nrow(df_all_m155)))
      dev.set(2)
      plot(df_all$logkd, model, col=cols)
      dev.set(3)
      print(length(sigmoids1))
      plot(-4:16, 1/(exp(pars.offset1) + 1), type="o", col="red", ylim=c(0, 0.5))
      lines(-4:16, 1/(exp(pars.offset2) + 1), type="o", col="blue")
      dev.set(4)
      plot(exp(pars.5p3pl7[1:(n_l7/2)]), exp(pars.5p3pl7[(n_l7/2 + 1):(n_l7)]))
      dev.set(5)
      plot(exp(pars.5p3pm1[1:(n_m1/2)]), exp(pars.5p3pm1[(n_m1/2 + 1):(n_m1)]))
      dev.set(6)
      plot(exp(pars.5p3pm155[1:(n_m155/3)]), exp(pars.5p3pm155[(n_m155/3 + 1):(n_m155*2/3)]))
      dev.set(7)
      plot(exp(pars.5p3pm155[1:(n_m155/3)]), exp(pars.5p3pm155[(n_m155*2/3 + 1):(n_m155)]))
    }
    tick <<- tick + 1
    loss <- sum((df_all$logkd - model)^2)
  }
  Gradient_numerical <- function(pars) {
    return(grad(LossFunction, pars, method="simple"))
  }
  # Gradient <- function(pars) {
  #   # Split up the parameters.
  #   pars.5p3p <- pars[1:n_pars_5p3p]
  #   pars.offset <- pars[(n_pars_5p3p + 1):(n_pars_5p3p + n_pars_offset)]
  #   pars_base <- exp(pars[length(pars)])
  #   # Assign the dG and sigmoids for each row of data.
  #   dG_5p3p <- exp(pars.5p3p[as.character(df_all$pos_5p3p)])
  #   sigmoids <- 1/(exp(pars.offset[as.character(df_all$offset)]) + 1)
  #   # Calculate the final model.
  #   model <- dG_5p3p*sigmoids + pars_base
  #   # Calculate the residuals.
  #   residuals <- 2*(model - df_all$logkd)
  #   # Calculate the derivative of each of the pairing coefficients.
  #   p_start <- 1
  #   p_stop <- n_pars_5p3p
  #   grad.5p3p <- sapply(p_start:p_stop, function(ind) {
  #     inds_sum <- which(as.character(df_all$pos_5p3p) == names(pars)[ind])
  #     sum((residuals*dG_5p3p*sigmoids)[inds_sum])
  #   })
  #   # Calculate the derivative of each of the offset coefficients.
  #   p_start <- p_stop + 1
  #   p_stop <- p_stop + n_pars_offset
  #   grad.offset <- sapply(p_start:p_stop, function(ind) {
  #     inds_sum <- which(as.character(df_all$offset) == names(pars)[ind])
  #     -sum((residuals*dG_5p3p*sigmoids*(1 - sigmoids))[inds_sum])
  #   })
  #   # Calculate the derivative of the base coefficient.
  #   grad.base <- sum(residuals)*pars_base
  #   # Concatenate the three components to make the entire gradient vector.
  #   grad_all <- c(grad.5p3p, grad.offset, grad.base)
  #   names(grad_all) <- names(pars)
  #   return(grad_all)
  # }
  # Initialize the parameters.
  # Solve for the parameters using optim.
  tick <- 0
  coefs <- optim(pars.init, LossFunction, gr=Gradient_numerical, method="L-BFGS-B",
                 control=list(maxit=1e7))$par
  coefs <<- coefs
  break
  # Extract the pairing coefficients.
  coefs_pairing <- exp(coefs[grep("\\|", names(coefs), value=TRUE, perl=TRUE)])
  # Remove the pairing coefficients.
  coefs <- coefs[grep("\\|", names(coefs), perl=TRUE, invert=TRUE)]
  # Extract the base coefficient.
  coefs_base <- exp(coefs[length(coefs)])
  # Extract the offset coefficients.
  coefs_offset <- 1/(exp(coefs[-length(coefs)]) + 1)
  # Re-order the offset coefficients according to -4:16.
  coefs_offset <- coefs_offset[order(as.integer(names(coefs_offset)))]
  len_mir <- nchar(kMirnaSeqs[mirna])
  # Define the possible 5-prime starting nucleotides and possible 3-prime
  # starting nucleotides, for overall matrix.
  len_k <- 4
  nucs_5p <- 9:(len_mir - len_k + 1)
  nucs_3p <- (9 + len_k - 1):len_mir

  # Define the pairing matrix.
  pairing_matrix <- matrix(NA, nrow=length(nucs_3p), ncol=length(nucs_5p),
                       dimnames=list(nucs_3p, nucs_5p))
  for (i in 1:length(coefs_pairing)) {
    coef_pairing <- coefs_pairing[i]
    name_pairing <- names(coefs_pairing)[i]
    name_split <- unlist(strsplit(name_pairing, split="\\|"))
    rowname_val <- name_split[2]
    colname_val <- name_split[1]
    pairing_matrix[rowname_val, colname_val] <- coef_pairing
  }
  coefs_pairing <- pairing_matrix
  # Normalize for the maximum offvalue, such that the maximum ofset is `1`.
  coefs_pairing <- coefs_pairing*max(coefs_offset)
  coefs_offset <- coefs_offset/max(coefs_offset)

  return(list(`pairing`=coefs_pairing, `offsets`=coefs_offset,
              `base`=coefs_base, `data`=df_all))
}



FitPairingOffsetAndMismatchModel <- function(
  mirna, experiment, n_constant=3, offsetmin=-4, offsetmax=16, len_min=4,
  len_max=11, pos_3p_max=23, pos_3p_min=9, sitelist="programmed",
  corrected_kds=TRUE, collapsemm=FALSE, supp_base=FALSE, kd_fc=TRUE
) {
  # Make a dataframe concatenating all the offset matrices into one.
  print("FitPairingOffsetAndMismatchModel")
  print("corrected_kds")
  print(corrected_kds)
  df_all <- SubfunctionCall(MakeFullNucleotideAndMisMatchContributionDf)
  print(head(df_all))
  print(tail(df_all))
  lens <- as.integer(df_all$pos_3p) - as.integer(df_all$pos_5p) + 1
  df_all <- df_all[which(lens >= len_min &
                         lens <= len_max &
                         as.integer(df_all$pos_3p) <= pos_3p_max &
                         as.integer(df_all$pos_3p) >= pos_3p_min), ]
  df_all <- drop.levels(df_all)
  print(nrow(df_all))
  break
  n_pars_5p3p <- length(unique(df_all$pos_5p3p))
  n_pars_offset <- length(unique(df_all$offset))
  n_pars_mm <- length(unique(df_all$mm))
  t_global <<- 0
  LossFunction <- function(pars, printout=TRUE) {
      # Initialize the parameters.
    pars.5p3p <- pars[1:n_pars_5p3p]
    pars.offset <- pars[(n_pars_5p3p + 1):(n_pars_5p3p + n_pars_offset)]
    pars.mm <- pars[(n_pars_5p3p + n_pars_offset + 1):(n_pars_5p3p + n_pars_offset + n_pars_mm)]
    pars_base <- exp(pars[length(pars)])

    dG_5p3p <- exp(pars.5p3p[as.character(df_all$pos_5p3p)])
    dG_mm <- exp(pars.mm[as.character(df_all$mm)])

    sigmoids <- 1/(exp(pars.offset[as.character(df_all$offset)]) + 1)
    model_vals <- dG_5p3p*dG_mm*sigmoids
    model_vals <<- model_vals
    residuals_global <<- (df_all$logkd - model_vals)^2
    sum((df_all$logkd - model_vals)^2)
  }
  Gradient_numerical <- function(pars) {
    return(grad(LossFunction, pars, printout=FALSE))
  }
  Gradient <- function(pars, printout=FALSE) {
    pars.5p3p <- pars[1:n_pars_5p3p]
    pars.offset <- pars[(n_pars_5p3p + 1):(n_pars_5p3p + n_pars_offset)]
    pars.mm <- pars[(n_pars_5p3p + n_pars_offset + 1):(n_pars_5p3p + n_pars_offset + n_pars_mm)]
    pars_base <- exp(pars[length(pars)])

    dG_5p3p <- exp(pars.5p3p[as.character(df_all$pos_5p3p)])
    dG_mm <- exp(pars.mm[as.character(df_all$mm)])
    sigmoids <- 1/(exp(pars.offset[as.character(df_all$offset)]) + 1)

    model <- dG_5p3p*dG_mm*sigmoids + pars_base
    model <- dG_5p3p*dG_mm*sigmoids
    residuals <- 2*(model - df_all$logkd)
    p_start <- 1
    p_stop <- n_pars_5p3p
    grad.5p3p <- sapply(p_start:p_stop, function(ind) {
      inds_sum <- which(as.character(df_all$pos_5p3p) == names(pars)[ind])
      sum((residuals*dG_5p3p*dG_mm*sigmoids)[inds_sum])
    })
    p_start <- p_stop + 1
    p_stop <- p_stop + n_pars_offset
    grad.offset <- sapply(p_start:p_stop, function(ind) {
      inds_sum <- which(as.character(df_all$offset) == names(pars)[ind])
      -sum((residuals*(dG_5p3p*dG_mm)*sigmoids*(1 - sigmoids))[inds_sum])
    })
    p_start <- p_stop + 1
    p_stop <- p_stop + n_pars_mm
    grad.mm <- sapply(p_start:p_stop, function(ind) {
      inds_sum <- which(as.character(df_all$mm) == names(pars)[ind])
      sum((residuals*dG_5p3p*dG_mm*sigmoids)[inds_sum])
    })
    grad.base <- sum(residuals)*pars_base
    grad.base <- 0
    grad_all <- c(grad.5p3p, grad.offset, grad.mm, grad.base)
    names(grad_all) <- names(pars)
    return(grad_all)
  }
  # Initialize the parameters
  pars.init <- rep(0, n_pars_5p3p + n_pars_offset + n_pars_mm + 1)
  names(pars.init) <- c(unique(as.character(df_all$pos_5p3p)),
                        unique(as.character(df_all$offset)),
                        unique(as.character(df_all$mm)),
                        "base")
  # Get the coefficients through the optimization routine.
  coefs <- optim(pars.init, LossFunction, gr=Gradient, method="L-BFGS-B",
                 control=list(maxit=1e7))$par
  # Extract the pairing coefficients.
  coefs_pairing <- exp(coefs[grep("\\|", names(coefs), value=TRUE, perl=TRUE)])
  # Remove the pairing coefficients from the list.
  coefs <- coefs[grep("\\|", names(coefs), perl=TRUE, invert=TRUE)]
  # Extract the mm coefficients.
  coefs_mm <- exp(coefs[grep("8mer", names(coefs), value=TRUE)])
  # Remove the mm coefficients.
  coefs <- coefs[grep("8mer", names(coefs), value=TRUE, invert=TRUE)]
  # Define the base coefficient.
  coefs_base <- exp(coefs[length(coefs)])
  coefs_base <- 0

  # Extract the offset coefficients.
  coefs_offset <- 1/(exp(coefs[-length(coefs)]) + 1)
  # Reorder the offset coefficients from -4 to +16.
  coefs_offset <- coefs_offset[order(as.integer(names(coefs_offset)))]
  len_mir <- nchar(kMirnaSeqs[mirna])
  # Define the possible 5-prime starting nucleotides and possible 3-prime
  # starting nucleotides, for overall matrix.
  len_k <- 4
  nucs_5p <- 9:(len_mir - len_k + 1)
  nucs_3p <- (9 + len_k - 1):len_mir
  # Define the output matrix.
  out_matrix <- matrix(NA, nrow=length(nucs_3p), ncol=length(nucs_5p),
                       dimnames=list(nucs_3p, nucs_5p))
  for (i in 1:length(coefs_pairing)) {
    coef_pairing <- coefs_pairing[i]
    name_pairing <- names(coefs_pairing)[i]
    name_split <- unlist(strsplit(name_pairing, split="\\|"))
    rowname_val <- name_split[2]
    colname_val <- name_split[1]
    out_matrix[rowname_val, colname_val] <- coef_pairing
  }
  coefs_pairing <- out_matrix
  # Normalize for the maximum offset coefficient value, such that the maximum
  # offset is `1`.
  coefs_pairing <- coefs_pairing*max(coefs_offset)
  coefs_offset <- coefs_offset/max(coefs_offset)
  # Normalize for the mean mm coefficient value, such that the mean mm
  # coefficient is `1`.
  coefs_pairing <- coefs_pairing*mean(coefs_mm)
  coefs_mm <- coefs_mm/mean(coefs_mm)

  return(list(`pairing`=coefs_pairing, `offsets`=coefs_offset, `mm`=coefs_mm,
              `base`=coefs_base, `data`=df_all))
}

GetAverageCenterOfPairing <- function(
  n_constant=3, offsetmin=-4, offsetmax=16,
  sitelist="programmed_suppcomp", collapsemm=FALSE, supp_base=FALSE
) {
  # let-7a
  df_all_l7 <- do.call(
    "rbind",
    lapply(offsetmin:offsetmax, MakeNucleotideContributionDf, mirna="let-7a-21nt",
           experiment="equil_c2_nb", n_constant=n_constant, collapsemm=collapsemm,
           supp_base=supp_base)
  )
  # miR-1
  df_all_m1 <- do.call(
    "rbind",
    lapply(offsetmin:offsetmax, MakeNucleotideContributionDf, mirna="miR-1",
           experiment="equil_c_nb", n_constant=n_constant, collapsemm=collapsemm,
           supp_base=supp_base)
  )
  # miR-155
  df_all_m155 <- do.call(
    "rbind",
    lapply(offsetmin:offsetmax, MakeNucleotideContributionDf, mirna="miR-155",
           experiment="equil_sc_nb", n_constant=n_constant, collapsemm=collapsemm,
           supp_base=supp_base)
  )
  # Now apply function to each of the three datasets taking the average center
  # of pairing over all three.
  xpos_move <<- 0
  len_k <- 7
  graphics.off()
  lapply(list(df_all_l7, df_all_m1, df_all_m155), function(df_all) {
    averages <- rep(0, 23 - 8 - len_k + 1)
    names(averages) <- as.character(9:(23 - len_k + 1))
    weights <- rep(0, 23 - 8 - len_k + 1)
    names(weights) <- names(averages)
    # center_average <- 0
    # weight <- 0
    pos_5p <- as.integer(as.character(df_all[, 2]))
    pos_3p <- as.integer(as.character(df_all[, 3]))
    df_all$lengths <- pos_3p - pos_5p + 1
    df_all <- df_all[which(df_all$lengths == len_k), ]
    apply(df_all, 1, function(row) {
      # print(row)
      logkd <- as.numeric(row[1])

      pos_5p <- as.character(row[2])
      averages[[pos_5p]] <<- (averages[[pos_5p]]*weights[[pos_5p]] + logkd)/(weights[[pos_5p]] + 1)
      weights[[pos_5p]] <<- weights[[pos_5p]] + 1
    })
    dev.new(xpos=20 + xpos_move, ypos=20, height=3, width=3)
    plot(as.integer(names(averages)), averages, type="o")
    xpos_move <<- xpos_move + 300
  })
}

MakePositionAndRegisterKdMatrix <- function(
  mirna, experiment, len_k, n_constant=3, sitelist="programmed", kdrel=FALSE
) {
  kds_collapsed <- SubfunctionCall(EquilPars, sitelist="programmed_collapsed")
  len_mir <- nchar(kMirnaSeqs[mirna])
  num_sites <- len_mir - 8 - len_k + 1
  starts <- 9:(9 + num_sites - 1)
  stops <- starts + len_k - 1
  sites <- sprintf("%smer-m%s.%s", len_k, starts, stops)
    vals <- t(sapply(sites, function(site) {
      kds <- SubfunctionCall(GetPositionalProgKds)
      num_pos <- 25 - len_k + 1
      mir_pos <- c("Seed", as.character(9:(9 + num_pos - 1)))

      kds_out <- matrix(0, nrow=18, ncol=1 + 25 - len_k + 1,
                        dimnames=list(mm_8mer_sites, mir_pos))
      for (mm_site in mm_8mer_sites) {
        for (pos in mir_pos) {
          if (pos == "Seed") {
            site_name <- sprintf("%s_Kd", mm_site)

            val <- kds_collapsed[site_name, 2]
          } else {
            site_name <- sprintf("%s|%s|%s_Kd", site, pos, mm_site)
            val <- kds[site_name, 2]
          }
          kds_out[mm_site, as.character(pos)] <- val
        }
      }
      if (kdrel) {
        kds_out <- (t(t(kds_out/kds_out[, 1])))
        kds_out <- kds_out[, 2:ncol(kds_out)]
      }
      return(exp(colMeans(log(kds_out), na.rm=TRUE)))
  }))
}


CollapseSeedSitesTable <- function(kd_data) {
  # Get all the unique sites by grabbing the last position of the rownames.
  mm_sites <- unique(gsub("^(.*)\\|(.*)\\|(.*)_Kd", rownames(kd_data),
                          replace="\\3"))
  thrp_sites <- unique(gsub("^(.*)\\|(.*)\\|(.*)_Kd", rownames(kd_data),
                          replace="\\1"))
  # Remove the AGO and bg rownames (do this by explicitly looking for them,
  # rather than assuming they are always the last two indeces).
  mm_sites <- grep("bg", grep("AGO", mm_sites, invert=TRUE, value=TRUE),
           invert=TRUE, value=TRUE) 
  all_sites <- c(kSeedSites, mm_sites)
  average_site_kds <- matrix(
    NaN, nrow=length(all_sites), ncol=1,
    dimnames=list(sprintf("%s|NA|NA_Kd", all_sites), "x")
  )
  for (site in all_sites) {
    grep_target <- sprintf("^%s\\|", site)
    kd_rows <- grep(grep_target, rownames(kd_data), perl=TRUE, value=TRUE)
    average <- mean(kd_data[kd_rows, 1])
    average_site_kds[sprintf("%s|NA|NA_Kd", site), 1] <- average
    kd_data <- kd_data[setdiff(rownames(kd_data), kd_rows), , drop=FALSE]
  }
  average_site_kds <- average_site_kds[!is.na(average_site_kds[, 1]), 1, drop=FALSE]
  kd_data <- rbind(kd_data[1:18, 1, drop=FALSE], average_site_kds, kd_data[19:nrow(kd_data), 1, drop=FALSE])
  kd_data
}

GetSeedKdMatrix <- function(kd_data, site) {
  # Define the relevant regular expression. This finds the site of interest,
  # And all spacings and mismatchtypes in programmed_region.
  grep_query <- sprintf("^%s\\|.*\\|8mer-mm[ACGT]._Kd$", site)
  # Get the row indeces using the regular expression.
  inds_kd <- grep(grep_query, rownames(kd_data), perl=TRUE)
  # Sub-select the Kds of relevance.
  kds_use <- kd_data[inds_kd, , drop=FALSE]
  # Get the positions of the site within the random library, and also get the
  # set of 18 8mer-mm sites within the programmed region of the library, to be
  # used as column and row names, respectively.
  lib_pos <- GetRandomPositionOfKds(kds_use)
  mm_sites <-GetProgrammedSiteOfKds(kds_use)
  # Convert the n x 1 matrix of Kd values into an 18 x (n/18) matrix, where the
  # rows and columns correspond to the 8mer-mm sites in the programmed and the
  # position of the non-programmed site within the random portion of library,
  # respectively.
  out_matrix <- matrix(kds_use[, 1], nrow=18, ncol=length(inds_kd)/18,
                       byrow=TRUE, dimnames=list(mm_sites, lib_pos))
  out_matrix
}


MakeMismatchAndPositionKdMatrix <- function(
  mirna, experiment, site, n_constant=3, sitelist="programmed", kdrel=FALSE,
  prog_canon=FALSE
) {
  kds_collapsed <- SubfunctionCall(EquilPars, sitelist="programmed_collapsed")

  mm_8mer_sites <- GetAll8merMmSites(mirna)
  if (prog_canon) {
    prog_sites <- c("8mer", "7mer-m8", "7mer-A1", "6mer", "6mer-m8", "6mer-A1",
                    mm_8mer_sites)
  } else {
    prog_sites <- mm_8mer_sites
  }
  len_k <- as.integer(unlist(strsplit(site, split="mer"))[[1]])
  if (prog_canon) {
    new <- TRUE
  } else {
    new <- FALSE
  }
  kds <- SubfunctionCall(GetPositionalProgKds)
  num_pos <- 25 - len_k + 1
  mir_pos <- c("Seed", as.character(9:(9 + num_pos - 1)))

  out <- matrix(0, nrow=length(prog_sites), ncol=1 + 25 - len_k + 1,
                dimnames=list(prog_sites, mir_pos))
  for (prog_site in prog_sites) {
    for (pos in mir_pos) {
      if (pos == "Seed") {
        site_name <- sprintf("%s_Kd", prog_site)

        val <- kds_collapsed[site_name, 2]
      } else {
        site_name <- sprintf("%s|%s|%s_Kd", site, pos, prog_site)
        val <- kds[site_name, 2]
      }
      out[prog_site, as.character(pos)] <- val
    }
  }
  # If `kdrel` is true, divide each positional kd by the kd corresponding to the
  # site in the programmed 
  if (kdrel) {
    out <- (t(t(out/out[, 1])))
    out <- out[, -1]
  }
  out
}

MakeNucleotideContributionDf <- function(
  mirna, experiment, offset, n_constant=3, sitelist="programmed_suppcomp",
  corrected_kds=FALSE, len_min=4, len_max=11, pos_3p_min=9, pos_3p_max=12,
  collapsemm=FALSE, supp_base=FALSE, offset_base=FALSE, site_base=NULL,
  nucletters=FALSE, loop=FALSE, cutoff=FALSE, kd_fc=TRUE
) {
  logkd_mat <- SubfunctionCall(MakePairingMatrix))
  logkd_df <- c(logkd_mat)
  pos_3p <- rep(rownames(logkd_mat), ncol(logkd_mat))
  pos_5p <- rep(colnames(logkd_mat), each=nrow(logkd_mat))
  offset_df <- rep(offset, ncol(logkd_mat)*nrow(logkd_mat))
  if (nucletters) {
    pos_mat <- matrix(0, nrow=length(logkd_df), ncol=15)
    colnames(pos_mat) <- paste0("p", 9:23)
    for (i in 1:length(logkd_df)) {
      minstr <- as.character(pos_3p[i])
      maxstr <- as.character(pos_5p[i])
      cols_str <- paste0("p", as.integer(minstr):as.integer(maxstr))
      inds_pos <- sapply(cols_str, function(col_str_i) {
        which(colnames(pos_mat) == col_str_i)
      })
      pos_mat[i, inds_pos] <- 1
    }
    output <- data.frame(`logkd`=logkd_df, pos_mat, `offset`=as.factor(offset_df))
  } else {
    output <- data.frame(`logkd`=logkd_df, `pos_5p`=pos_5p, `pos_3p`=pos_3p,
                         `pairing`=sprintf("%s|%s", pos_5p, pos_3p),
                         `offset`=as.factor(offset_df))      
  }
  output <- output[!(is.na(output[, 1])), ]
  len <- as.integer(as.character(output$pos_3p)) - as.integer(as.character(output$pos_5p)) + 1
  output <- output[which(len >= len_min & len <= len_max & as.integer(as.character(output$pos_3p)) >= pos_3p_min & as.integer(as.character(output$pos_3p)) <= pos_3p_max), ]
  print(dim(output))
  # This is a new conditional that imposes a cutoff of less than 20.
  if (cutoff & (nrow(output) < 25)) {
    output <- data.frame(`logkd`=c(), `pos_5p`=c(), `pos_3p`=c(), `pairing`=c(),
                         `offset`=c())
  }
  return(output)
}



GetBestSuppCompKd <- function(
  mirna, experiment, n_constant=0, sitelist="programmed_suppcomp",
  collapsemm=FALSE, supp_base=FALSE
) {
  print("in MakeNucleotideContributionDf")
  kds <- SubfunctionCall(EquilPars)
  grep_string <- "^[0-9]{1,2}mer-m[0-9]{1,2}\\.[0-9]{1,2}\\|.*\\|Comp_Kd$"
  print(grep_string)
  inds <- grep(grep_string, rownames(kds), perl=TRUE)
  kds_use <- kds[inds, ]
  ind_min <- which.min(kds_use[, 2])
  print(rownames(kds_use)[ind_min])
}


# Used for Figures showing the pairing
MakePairingMatrix <- function(
  mirna, experiment, offset, n_constant=3, sitelist="progthrp_suppcomp",
  len_lim=c(4, 11), corrected_kds=TRUE, sumseed=FALSE, supp_base=FALSE,
  offset_base=FALSE, site_base=NULL, error=FALSE, loop=FALSE, kd_fc=TRUE
) {
  print("old function")
  # Generate the site names.
  mm8mer_sites <- GetAll8merMmSites(mirna)
  # canon_sites <- c("8mer", "7mer-m8", "7mer-A1", "6mer")
  if (corrected_kds & !grepl("equilibrium", experiment)) {
    kds_thrp <- SubfunctionCall(ApplyKdCorrection, rand_n_constant=n_constant,
                                prog_n_constant=n_constant,
                                prog_sitelist=sitelist)
  } else {
    if (mirna == "miR-1" & experiment == "equilibrium") {
      buffer <- TRUE
      combined <- FALSE
    } else if (mirna == "miR-7-23nt" & experiment == "equilibrium2_nb") {
      combined <- FALSE
    }
    kds_thrp <- SubfunctionCall(EquilPars)
  }
  if (!is.null(site_base) & !grepl("suppcomp", sitelist)) { # For one base site.
    kd_base <- kds_thrp[sprintf("%s_Kd", site_base), 2]
    base_str <- site_base   
  } else if (supp_base) {
    if (sumseed) {
      kd_base <- kds_thrp["Supp_Kd", 2]
    } else {
      kd_base <- GeoMean(kds_thrp[sprintf("%s_Kd", kCanonicalSites), 2])
    }
    base_str <- "Supp"
  } else if (offset_base) {
    if (sumseed) {
      kd_base <- kds_thrp["Offset_Kd", 2]
    } else {
      kd_base <- GeoMean(kds_thrp[sprintf("%s_Kd", c("6mer-m8", "6mer-A1")), 2])
    }    
    base_str <- "Offset" 
  } else {
    if (sumseed) {
      kd_base <- kds_thrp["Comp_Kd", 2]
    } else {
      kd_base <- GeoMean(kds_thrp[sprintf("%s_Kd", mm8mer_sites), 2])
    }
    base_str <- "Comp"    
  }
  if (!kd_fc) {
    kd_base <- 1
  }
  len_mir <- nchar(kMirnaSeqs[mirna])
  # Define the possible 5-prime starting nucleotides and possible 3-prime
  # starting nucleotides, for overall matrix.
  nucs_5p <- 9:(len_mir - len_lim[1] + 1)
  nucs_3p <- (9 + len_lim[1] - 1):len_mir
  # Define the output matrix.
  out_matrix <- matrix(NA, nrow=length(nucs_3p), ncol=length(nucs_5p),
                       dimnames=list(nucs_3p, nucs_5p))
  for (kmer_len in len_lim[1]:len_lim[2]) {
    starts_i <- 9:(len_mir - kmer_len + 1)
    stops_i <- starts_i + kmer_len - 1
    sites <- sprintf("%smer-m%s.%s", kmer_len, starts_i, stops_i)
    if (loop) {
      dist_str <- 9 + offset
    } else {
      dist_str <- starts_i + offset
    }
    sites_use <- sprintf("%s|%s|%s_Kd", sites, dist_str, base_str)
    if (error) {
      kds <- kds_thrp[sites_use, 5]/kds_thrp[sites_use, 3]
    } else {
      kds <- kds_thrp[sites_use, 2]/kd_base
    }
    for (i in 1:length(starts_i)) {
      out_matrix[as.character(stops_i[i]), as.character(starts_i[i])] <- kds[i]
    }
  }
  return(log10(1/out_matrix))
}

MakeMismatchByOffsetMatrix <- function(
  mirna, experiment, site_len, start_pos, n_constant=3, sitelist="progthrp",
  corrected_kds=TRUE, kd_fc=TRUE, rand_base=FALSE, supp_base=FALSE,
  offset_lim=c(-4, 16), add_comp_sites=FALSE, collapsemm_addcomp=FALSE
) {
  print("in MakeMismatchByOffsetMatrix")
  # 1. Load either the corrected Kd Matrix or the uncorrected Kd matrix.
  if (corrected_kds & !grepl("equilibrium", sitelist)) {
    kds <- SubfunctionCall(ApplyKdCorrection, prog_sitelist=sitelist,
                           rand_sitelist="randthrp")
  } else {
    if (mirna == "miR-1" & experiment == "equilibrium") {
      combined <- FALSE
      buffer <- TRUE      
    } else if (mirna == "miR-7-23nt" & experiment == "equilibrium2_nb") {
      combined <- FALSE
    }
    kds <- SubfunctionCall(EquilPars)
    if (add_comp_sites) {
      kds_comp <- SubfunctionCall(
        EquilPars, sitelist=sprintf("%s_suppcomp", sitelist), 
        collapsemm=collapsemm_addcomp
      )
    }

  }
  print(kd_fc)
  # Make the three-prie site string.
  thrp_site <- sprintf("%smer-m%s\\.%s", site_len, start_pos,
                      start_pos + site_len - 1)
  # Make the overall grep string that uses `thrp_site`.
  if (supp_base) {
    print("supp base")
    if (add_comp_sites) {
      target_comp <- sprintf("%s\\|.*\\|Comp_Kd", thrp_site)
    }
    target <- sprintf("%s\\|.*\\|(8mer|7mer-A1|7mer-m8|6mer|6mer-A1|6mer-m8)_Kd", thrp_site)
  } else {
    target <- sprintf("%s\\|.*\\|8mer-mm[ATCG][2-7]_Kd", thrp_site)
  }
  # Get the rownames corresponding to that threeprime site.
  inds_use <- grep(target, rownames(kds), perl=TRUE, value=TRUE)
  # Get the strings associated with the mismatch site and the offset position.
  mm_sites <- gsub("^(.*)\\|(.*)\\|(.*)_Kd", inds_use, replacement="\\3",
                   perl=TRUE)
  offsets <- as.integer(
    gsub("^(.*)\\|(.*)\\|(.*)_Kd", inds_use, replacement="\\2", perl=TRUE)
  ) - start_pos
  # Pre-allocate the output matrix, with rows corresponding to the mismatch
  # sites and the columns corresponding to the offset positions.
  if (supp_base) {
    if (add_comp_sites) {
      nrow_mat <- 7
      row_names <- c(kSeedSites, "8mer-mm2-7")
    } else {
      nrow_mat <- 6
      row_names <- kSeedSites
    }
  } else {
    nrow_mat <- 18
    row_names <- GetAll8merMmSites(mirna)
  }
  mat_out <- matrix(NA, nrow=nrow_mat,
                    ncol=length(-4:16),
                    dimnames=list(row_names, -4:16))
  # Boolean allowing the mismatch Kds from the random libraries (i.e., when
  # `sitelist` == "equilibrium") rather than those calcualted from the
  # programmed libraries.
  if (rand_base) {
    if (mirna %in% c("let-7a-21nt", "let-7a_plus1", "let-7a_minus1",
                     "let-7a_miR-155")) {
      mirna_rand <- "let-7a"
    } else if (mirna == "miR-155_let-7a") {
      mirna_rand <- mirna
    } else {
      mirna_rand <- mirna
    }
    kds_base <- EquilPars(
      mirna, "equilibrium", n_constant, "bipartite_random"
    )[sprintf("%s_Kd", row_names), 2]
  } else {
    kds_base <- kds[sprintf("%s_Kd", row_names), 2]
    if (add_comp_sites) {
      # Get the indeces for the 8mer-mm sites in the compensatory kd data set.
      if (collapsemm_addcomp) {
        kd_base_comp <- kds_comp["Comp_Kd", 2]
      } else {
        inds_mm_comp <- grep("^8mer-mm[ACTG][2-7]_Kd", rownames(kds_comp),
                             perl=TRUE, value=TRUE)
        if (length(inds_mm_comp) != 18) {
          print("problem with number of mismatch kds.")
          break
        }
        # Get the geometric average of the sites.
        kd_base_comp <- GeoMean(kds_comp[inds_mm_comp, 2])
      }
      kds_base[length(kds_base)] <- kd_base_comp
    }
  }
  # Assign the actual miRNA sites to be used for the columns that aren't the
  # the base Kds.
  kds_use <- kds[inds_use, 2]
  # Make the matrix to be thread to input values into the matrix.
  df_out <- cbind(mm_sites, offsets, kds_use)
  # Make sure that the values used are only those that conform to the offset
  # limits.
  df_out <- df_out[which(offsets >= offset_lim[1] & offsets <= offset_lim[2]), ]
  if (length(dim(df_out)) != 0) {
    apply(df_out, 1, function(row) {
      if (row[1] != "") {
        mat_out[row[1], row[2]] <<- as.numeric(row[3])
      }
    })    
  }
  print(mat_out)
  if (add_comp_sites) {
    inds_use <- grep(target_comp, rownames(kds_comp), perl=TRUE, value=TRUE)
    print(inds_use)
    # Get the strings associated with the mismatch site and the offset position.
    offsets <- as.integer(
      gsub("^(.*)\\|(.*)\\|(.*)_Kd", inds_use, replacement="\\2", perl=TRUE)
    ) - start_pos
    print(offsets)
      kds_use <- kds_comp[inds_use, 2]
  # Make the matrix to be thread to input values into the matrix.
    df_out <- cbind(as.character(offsets), kds_use)
  # Make sure that the values used are only those that conform to the offset
  # limits.
    df_out <- df_out[which(offsets >= offset_lim[1] & offsets <= offset_lim[2]), ]
    if (length(dim(df_out)) != 0) {
      apply(df_out, 1, function(row) {
        if (row[1] != "") {
          mat_out["8mer-mm2-7", row[1]] <<- as.numeric(row[2])
        }
      })    
    }
    rownames(mat_out)[nrow(mat_out)] <- "Seed mm"
  }
  # Populate the matrix as appropriate.
  # Add the left-hand-most column to the matrix, and name it "None."
  mat_out <- cbind(kds_base, mat_out)
  colnames(mat_out)[1] <- "None"
  # Boolean allowing the matrix to be a fold-change matrix, rather than the
  # absolute values of the kds, made by dividing the entire matrix by the first.
  if (kd_fc) {
    mat_out <- mat_out/(mat_out[, 1])
  }
  if (supp_base) {
    mat_out <- mat_out[rev(1:nrow(mat_out)), ]
  }
  # Return the matrix
  return(mat_out)
}



# FitTriplePairingAndOffsetModel()
# break


