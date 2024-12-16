################################################################################
#GenerateSiteKds_PAPER.py
################################################################################

source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")
print("Out of general")

# Initial parameters and constants.
# args       <- commandArgs(trailingOnly=TRUE)
# mirna      <- args[1]
# experiment <- args[2]
# n_constant <- args[3]
# sitelist   <- args[4]

# if ("-reps" %in% args) {
#   reps <- as.integer(args[which(args == "-reps") + 1])
# } else {
#   reps <- 40
# }

# if (sitelist == "12mers") {
#   mirna.start = as.integer(args[which(args == "-mir_start") + 1])
#   str.mir_start = sprintf("_%s-%s", mirna.start, mirna.start + 3)
# } else {
#   mirna.start = FALSE
#   str.mir_start = ""
# }

# if ("-uniq" %in% args) {
#   uniq <- TRUE
#   str.unique <- "_uniq"
# } else {
#   uniq <- FALSE
#   str.unique <- ""
# }

# if ("-buffer" %in% args) {
#   buffer <- TRUE
#   str.buffer <- "_buffer3p"
# } else {
#   buffer <- FALSE
#   str.buffer <- ""
# }

# if ("-nocombI" %in% args) {
  combined <- FALSE
  str.combined <- "_nocombInput"
# } else {
#   combined <- TRUE
#   str.combined <- ""
# }

# if ("-fixed" %in% args & mirna %in% c("miR-7-23nt", "miR-7-24nt", "miR-7-25nt")
#     & experiment == "equilibrium2_nb") {
#   fixed <- TRUE
#   altfixed <- FALSE
#   str.fixed <- "_fixed"
# } else if ("-altfixed" %in% args & mirna %in% c("miR-7-23nt", "miR-7-24nt", "miR-7-25nt")
#     & experiment == "equilibrium2_nb") {
#   fixed <- TRUE
#   altfixed <- TRUE
#   str.fixed <- "_altfixed"
# } else {
  altfixed <- FALSE
  fixed <- FALSE
#   str.fixed <- ""
# }

# if ("-L" %in% args) {
#   L <- as.integer(args[which(args == "-L") + 1])
#   str.L <- sprintf("_L%s", L)
# } else {
  L <- 100
#   str.L <- ""
# }

# if ("-single" %in% args) {
#   single <- TRUE
#   str.single <- "_singleonly"
# } else {
#   single <- FALSE
#   str.single <- ""
# }


sXc_path <- file.path("/lab/solexa_bartel/nbisaria/analysis/miR-1/equilibrium/",
                      "seandata/site_counts_streamline3p_combined/Counts.txt")

kds_path <- file.path("/lab/solexa_bartel/nbisaria/analysis/miR-1/equilibrium",
                      "seandata/kds_site_counts_streamline3p_combined",
                      "kds_singlebg_multinomial_gradient_Kds.txt")

sXc_i <- read.table(sXc_path)

colnames(sXc_i) <- gsub("X", "", colnames(sXc_i))

kds_nb <- read.table(kds_path)

print(head(sXc_i))
print(head(kds_nb))

plot(kds_nb[, 1], sXc_i[rownames(kds_nb), 2], log='xy')







# MAIN #########################################################################
# 1. Get data table:
mirna <- "miR-1"
mir_list <- unlist(strsplit(mirna, ","))
# exp_list  <- unlist(strsplit(experiment, ","))
# mir_exp_list <- data.frame(mirna=rep(mir_list, length(exp_list)),
#                            experiment=rep(exp_list,each=length(mir_list)))

print("Mirna experiment list:")
# print(mir_exp_list)
# sXc <- apply(mir_exp_list, 1, function(row) {
#   sXc_i <- SitesXCounts(row[1], experiment=row[2], n_constant=n_constant,
#                         sitelist=sitelist, mirna.start=mirna.start,
#                         uniq=uniq, buffer=buffer)
#   if (single) {
#     msXc <- SitesXCounts(row[1], experiment=row[2], n_constant=n_constant,
#                          sitelist=sitelist, mirna.start=mirna.start,
#                          uniq=uniq, buffer=buffer, multisite=TRUE)
#     regex_split         <- ",\\(\\d*\\),"
#     regex_split_capture <- ",\\((\\d*)\\),"
#     regex_site <- "[^,]*"
#     regex_site_capture <- "([^,]*)"

#     single_sites = grep(regex_split, rownames(msXc), invert=TRUE)
#     double_sites = grep(sprintf("^%s%s%s$", regex_site, regex_split,
#                                 regex_site), rownames(msXc), perl=TRUE)
#     nonoverlapping_sites = grep("\\|", rownames(msXc), invert=TRUE, perl=TRUE)
#     single_nonoverlapping_sites = intersect(single_sites, nonoverlapping_sites)
#     sXc_out <- msXc[single_nonoverlapping_sites, ][rownames(sXc_i), ]
#   } else {
#     sXc_out <- sXc_i
#   }
#   sXc_out
# })

sXc <- list(sXc_i)


if (length(sXc) > 1) {
  str.global <- "_global"
} else {
  str.global <- ""
}
# str_condition <- sprintf("%s_%s%s%s%s%s%s%s%s%s_PAPER", n_constant, sitelist,
#                          str.unique, str.buffer, str.mir_start, str.combined,
#                          str.global, str.fixed, str.L, str.single)

# kOutputFileMean <- GetAnalysisPath(mir_list[[1]], experiment,
#                                    condition=sprintf("%s_%s%s%s%s_PAPER_logmean",
#                                                      n_constant, sitelist,
#                                                      str.mir_start,
#                                                      str.combined, str.global),
#                                    analysis_type="kds_PAPER")

kOutputFileMean <- "temp_nb_kds_miR-1.txt"
pars_sd <- 0.1
plot_ <- FALSE
plotname <- NULL


InitializeEquilSitePars <- function(sXc, combined=TRUE, fixed=FALSE) {
  n_mir = length(sXc)
  print(head(sXc[[1]]))
  l <- SubfunctionCall(GetInputEquil, sXc=sXc[[1]])
  data <- SubfunctionCall(GetDataEquil, sXc=sXc[[1]])
  kds <- log(Norm(l)/Norm(rowSums(data)))
  kds <- kds - kds[length(kds)]
  kds <- kds/2


  names(kds) <- paste0(rownames(sXc[[1]]), "_Kd")
  
  if (fixed) {
    if (altfixed) {
      combined_fixed <- TRUE
    } else {
      combined_fixed <- FALSE
    }
    pars <- EquilPars(mirna, experiment=experiment, n_constant=n_constant,
                      sitelist=sitelist, uniq=uniq, combined=combined_fixed, global=TRUE)
    bgs <- log(pars[sprintf("bg_%s", mirna),]$Mean)
    As <- log(pars[sprintf("AGO_%s", mirna),]$Mean)
  } else {
    bgs <- rep(log(0.1), n_mir)
    As <- rep(log(10), n_mir)    
  }
  names(bgs) <- sprintf("bg_%s", mir_list)
  names(As) <- sprintf("AGO_%s", mir_list)
  pars <- c(kds, bgs, As)
  random_pars <- rnorm(length(pars), 0, pars_sd)
  if (fixed) {
    random_pars[length(pars) - 1] = 0
    random_pars[length(pars)    ] = 0
  }
  names(random_pars) <- names(pars)
  # pars <- pars + random_pars
  pars["None_Kd"] <- 0
  return(pars)
}

combined <- FALSE
fixed <- FALSE
tick <- 0
OptimizeEquilSitePars <- function(sXc, pars=NULL, fixed=fixed, AGOfixed=FALSE,
                                  plotname_=plotname) {
  print(combined)
  print(fixed)
  time_start <- proc.time()[3]
  tick <- 0
  if (is.null(pars) == TRUE) {
    initial.pars <- InitializeEquilSitePars(sXc, combined=combined, fixed=fixed)
  } else {
    initial.pars <- pars
  }
    L_ <- L
  if (length(sXc) > 1) {
    data <- do.call(cbind, sapply(sXc, function(sXc_i) {
      as.matrix(sXc_i[, which(!(colnames(sXc_i) %in% c("I", "I_combined", "0")))])
    }))
  } else {
    data <- sXc[[1]][, which(!(colnames(sXc[[1]]) %in% c("I", "I_combined", "0")))]
  }
  n_j <- sapply(sXc, function(sXc_i) {
    ncol(sXc_i[, which(!(colnames(sXc_i) %in% c("I", "I_combined", "0")))])
  })
  l <- Norm(sXc[[1]][, "I"] + 1)*L_
  Y <- colSums(data)
  print(head(data))
  print(head(data))
  dil <- as.numeric(colnames(data))/100.0
  n_i <- nrow(data)
  n_mir <- length(sXc)
  data_vec <- as.numeric(as.matrix(data))
  lower <- log(rep(c(1e-6, 0.001, 0.1), c(n_i, n_mir, n_mir)))
  upper <- log(rep(c(1e4, 10, 10), c(n_i, n_mir, n_mir)))

  if (fixed) {
    zero_grad <- n_i + c(0, 1, 2)
  } else if (AGOfixed) {
    zero_grad <- n_i + c(0, 1, 2)
  } else {
    zero_grad <- c(n_i)
  }
  # print(grep("None", rownames(sXc[[1]])))
  # zero_grad <- grep("None", rownames(sXc[[1]]))
  # print(dil)
  # # print(l)
  # print(L_)
  # print(n_i)
  # print(n_j)
  # print(n_mir)
  # print(zero_grad)
  # print(initial.pars[zero_grad])
  # break
  # print(fixed)
  # print(upper)
  # print(lower)
  n_z <- length(zero_grad)
  solution <- optim(initial.pars,
                    CostC,
                    gr = GradC,
                    data = data_vec,
                    dil  = dil,
                    l = l,
                    L = L_,
                    Y = Y,
                    n_i = n_i,
                    n_j = n_j,
                    n_mir = n_mir,
                    zero_grad = zero_grad,
                    n_z = n_z,
                    fixed = fixed,
                    upper_ = upper,
                    lower_ = lower,
                    plot_ = plot_,
                    plotname = plotname_,
                    method = "L-BFGS-B",
                    lower=lower,
                    upper=upper,
                    control = list(maxit=10000000, factr=10000, fnscale=1))
  output.pars <- solution$par
  time_elapsed <- time_start - proc.time()[3]
  # print(time_elapsed)
  output.pars/log(10)
}


print(head(sXc[[1]]))

pars.init <- InitializeEquilSitePars(sXc)

graphics.off()

pars.MLE <- OptimizeEquilSitePars(sXc, fixed=fixed)
# pars.MLE_2 <- OptimizeEquilSitePars(sXc, fixed=fixed)

dev.new(xpos=20, ypos=320, height=5, width=5)


names(pars.MLE) <- gsub("_Kd", "", names(pars.MLE))
plot(10^pars.MLE[rownames(kds_nb)], kds_nb[, 1],
     col = c(rep("black", nrow(kds_nb) - 1), "blue"), log='xy',
     xlim=c(1e-5, 1e4), ylim=c(1e-5, 1e4))
abline(0, 1, lty=2)
print(10^mean(log10(kds_nb[, 1]) - pars.MLE[rownames(kds_nb)]))
break
# print(sprintf("%.16f", 10^pars.MLE))
resampled.pars <- matrix(NaN, nrow=length(pars.MLE), ncol=5*reps,
                         dimnames=list(names(pars.MLE), 1:(5*reps)))
print("HERE")

for (j in 0:(ncol(sXc[[1]]) - 4)) {
  message("column")
  message(j)
  for (i in 1:reps) {
    i_full <- j*reps + i
    tick <- 0
    # Resample the data:
    sXc.resample <- lapply(sXc, function(sXc_i) {
      sXc_new <- apply(sXc_i, 2, function(col) {rmultinom(1, sum(col), col)})
      rownames(sXc_new) <-rownames(sXc_i)
      return(sXc_new)
    })
    # If not 12mers, remove a column:
    if (sitelist != "12mers") {
      # if (mirna == "miR-7-23nt" & experiment == "equilibrium2_nb") {
      #   min_col <- 4
      # } else {
      #   min_col <- 3
      # }
      sXc.resample.withhold <- lapply(sXc.resample, function(sXc_i) {
        leave_out <- j + 3
        sXc_i[, -leave_out]
      })
    } else {
      sXc.resample.withhold <- sXc.resample
      plotname <- paste0("12merfits/", paste(mirna, experiment, n_constant, sitelist, mir_start,
                     sep="_"), "_rep", j*reps + i, ".pdf")
    }
    if (length(sXc.resample.withhold) != 1) {
      sXc.resample.withhold <- sXc.resample
    }
    if (i == 1) {
      print(colnames(sXc.resample.withhold[[1]]))    
    }
    if (mirna == "miR-7-23nt" & experiment == "equilibrium2_nb") {
      AGOfixed <- TRUE
    } else {
      AGOfixed <- FALSE
    }
    print(sXc.resample.withhold)
    resampled.pars[, i_full] <- OptimizeEquilSitePars(sXc.resample.withhold,
                                                      pars = log(10)*pars.MLE,
                                                      fixed=fixed, AGOfixed=AGOfixed,
                                                      plotname_=plotname)
    resampled.pars.sort <- t(apply(resampled.pars[, 1:i_full, drop=FALSE], 1,
                                   sort))
    output <- 10^cbind(pars.MLE,
                       rowMeans(resampled.pars, na.rm=TRUE),
                       resampled.pars.sort[, ceiling(0.025*i_full)],
                       resampled.pars.sort[, ceiling(0.5*i_full)],
                       resampled.pars.sort[, ceiling(0.975*i_full)])
    colnames(output) <- c("Full","Mean", "Lower_CI", "Median", "Upper_CI")
    sapply(mir_list, function(mirna) {
      kOutputFile <- GetAnalysisPath(mirna, experiment, condition=str_condition,
                               analysis_type="kds_PAPER")
      print(kOutputFile)
      write.table(file=kOutputFile, output, sep="\t", quote = FALSE)
    })
    if (sitelist == "12mers") {
      kds.mean <- log(10^rowMeans(resampled.pars, na.rm=TRUE))
      # print(kds.mean[1:5])
      write.table(file=kOutputFileMean, kds.mean, sep="\t", quote=FALSE,
                  row.names=FALSE)
    }
  }
}

print(ceiling(0.025*i_full))
print(ceiling(0.975*i_full))
write.table(file=sprintf("%s_temp_sampleout.txt", mirna), cbind(pars.MLE, resampled.pars))
print(sXc)
print(output)
warnings()





