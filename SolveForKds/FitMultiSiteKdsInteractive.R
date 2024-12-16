################################################################################
#GenerateSiteKds_PAPER.py
################################################################################

# THIS IS THE OFFICIAL SCRIPT FOR FINALIZING THE MANUSCRIPT.
# THIS SCRIPT HAS THESE DECISIONS MADE:
# 1.) PROTEIN IS FIT AS A PARAMETER
# 2.) A SINGLE BACKGROUND PARAMETER IS GIVEN TO ALL FIVE CONCENTRATIONS.
# 3.) LIST, START, AND STOP POSITION ARE MANDATORY PARAMETERS.
# 4.) ALL INPUTS ARE COMBINED ACROSS FOR THESE ANALYSES.

#### UNCOMMENT HERE (191211)

library(colorspace)
# library(parallel)
library(data.table)
library(numDeriv)
graphics.off()
source("/lab/bartel1_ata/mcgeary/computation/AgoRBNS/general/general.R")
# Initial parameters and constants.
# args       <- commandArgs(trailingOnly=TRUE)
# mirna      <- args[1]
# experiment <- args[2]
# n_constant <- args[3]
# sitelist   <- args[4]


mirna <- "let-7a"
experiment <- "equilibrium_mmseed_nb"
n_constant <- "5"
sitelist <- "programmedversusrandom"


# if (length(args) == 5) {
#   reps = as.integer(args[5])
# } else {
#   reps = 200
# }

# if ("-nocombI" %in% args) {
#   combined <- FALSE
#   str.combined <- "_nocombInput"
# } else {
#   combined <- TRUE
#   str.combined <- ""
# }

# mirna <- "miR-155"
# experiment <- "equilibrium"
# n_constant <- 5
# sitelist <- "mismatch_and_threeprime"
reps <- 200
combined <- FALSE
str.combined <- ""
# Loads general functions used in AgoRBNS analysis.

# TODO make a new table of miRNA concentrations or just convert within the
# script.
stockago <- read.table(paste0("/lab/bartel1_ata/mcgeary/computation/AgoRBNS",
                       "/SolveForKds/k_c_stockago.txt"), row.names=1,
                       header=TRUE, sep="\t")
k.c.stockago <- stockago[mirna,experiment]
k.c.lib <- 100

# MAIN #########################################################################
# 1. Get data table:
sXc <- SitesXCounts(mirna, experiment, n_constant, sitelist)
print(sXc["None", ])
msXc <- SitesXCounts(mirna, experiment, n_constant, sitelist, multisite=TRUE)
print(msXc["None", ])



# This is equivalent to ",(`any digits`),".
re_split         <- ",\\(\\d+\\),"
# This is the same but captures the digits within the parentheses.
re_split_capture <- ",\\((\\d+)\\),"
# This is equivalent to `stretch of characters without a comma`.
re_site <- "[^,\\|]+"
# This is the same but captures the site.
re_site_capture <- "([^,\\|]+)"
# This is equivalent to "|`any digits`|".
re_ol <- "\\|\\d+\\|"
# This is the same but captures the digits within.

# grab single sites.
ss_inds <- grep(sprintf("^%s$", re_site), rownames(msXc), perl=TRUE)
# Grab double sites.
# "`stretch without a comma`,(`any digits`),`strech without a comma`."
ds_inds <- grep(paste0("^", re_site, re_split, re_site, "$"),
                rownames(msXc), perl=TRUE)
ts_inds <- grep(paste0("^", re_site, re_split, re_site, re_split, re_site, "$"),
                rownames(msXc), perl=TRUE)

dol_inds <- grep(paste0("^", re_site, re_ol, re_site, "$"), rownames(msXc),
                 perl=TRUE)
tol_inds <- grep(paste0("^", re_site, re_ol, re_site, re_ol, re_site, "$"),
                 rownames(msXc), perl=TRUE)

# Gets the indeces of sites that have one site on the left followed by two
# overlapping sites on the right.
r_ol_ts_inds <- grep(paste0("^", re_site, re_ol, re_site, re_split, re_site,
                            "$"), rownames(msXc), perl=TRUE)
l_ol_ts_inds <- grep(paste0("^", re_site, re_split, re_site, re_ol, re_site),
                     rownames(msXc), perl=TRUE)

sXc_1s <- msXc[ss_inds, ][rownames(sXc), ]
sXc_2s_dist <- msXc[ds_inds, ]
sXc_3s_dist <- msXc[ts_inds, ]

sXc_dol <- msXc[dol_inds, ]
sXc_tol <- msXc[tol_inds, ]

sXc_r_ol_ts <- msXc[r_ol_ts_inds, ]
sXc_l_ol_ts <- msXc[l_ol_ts_inds, ]

message("This is the sum of all single site reads, double-site reads, and double")
message("partially overlapping reads:")
print((colSums(sXc_1s) + colSums(sXc_2s_dist) + colSums(sXc_3s_dist) +
       colSums(sXc_dol) + colSums(sXc_tol) + colSums(sXc_r_ol_ts) +
       colSums(sXc_l_ol_ts))/colSums(sXc))


# This is a check making sure that each of the sXc objects within the list are
# mutually orthogonal.
for (i in seq(length(list_sXcs))) {
  for (j in setdiff(seq(length(list_sXcs)), i)) {
    print(i)
    print(j)
    print(intersect(rownames(list_sXcs[[i]]), rownames(list_sXcs[[j]])))
  }
}


# Separate all instances of all of the 3p supplementary sites from the double
# site matrix.
ThrPSupInds <- grep("^.*mer-m[^238m].*8mer-mm.*$", rownames(sXc_2s_dist), perl=TRUE)
sXc_3psup_sites <- sXc_2s_dist[ThrPSupInds, ]
sXc_2s_dist <- sXc_2s_dist[-ThrPSupInds, ]

# Collect all instances of all of the 3p supplementary sites.
ThrPSupInds <- grep("^.*mer-m[^238m].*8mer-mm.*$", rownames(sXc_2s_dist), perl=TRUE)
sXc_3psup_sites <- sXc_2s_dist[ThrPSupInds, ]

# Now make a regex capture group to get the six canonical seed sites and the 18
# different single-mismatch seed sites, in order to parse the positional
# dependence of the sites in the programmed library.
regex_seed_double <- "^(8mer|7mer-m8|7mer-A1|6mer|6mer-A1|6mer-m8|8mer-mm.{2,2}),.*8mer-mm.*$"
DoubleSiteInds <- grep(regex_seed_double, rownames(sXc_2s_dist), perl=TRUE)
sXc_2s_seed_sites <- sXc_2s_dist[DoubleSiteInds, ]
sXc_2s_dist <- sXc_2s_dist[-DoubleSiteInds, ]



# This portion of the script splits up each two-site distance name into the
# left-hand site, the distance between the two sites, and the right-hand site.
site_2s_l <- gsub(sprintf("^%s%s%s$", regex_site_capture,
                         regex_split_capture, regex_site_capture),
                 rownames(sXc_2s_dist), replace="\\1", perl=TRUE)
dist <- as.integer(gsub(sprintf("^%s%s%s$", regex_site_capture,
                         regex_split_capture, regex_site_capture),
                 rownames(sXc_2s_dist), replace="\\2", perl=TRUE))

site_2s_r <- gsub(sprintf("^%s%s%s$", regex_site_capture,
                         regex_split_capture, regex_site_capture),
                 rownames(sXc_2s_dist), replace="\\3", perl=TRUE)




# sXc_norm <- t(t(sXc)/colSums(sXc))
# sXc_1s_norm <- t(t(sXc_1s)/colSums(sXc))


# sXc_R <- sXc_norm[, 3:(ncol(sXc) - 1)]/sXc_norm[, 1 + combined]
# sXc_1s_R <- sXc_1s_norm[, 3:(ncol(sXc) - 1)]/sXc_1s_norm[, 1 + combined]

# x <- c(40, 12.65, 4, 1.265, 0.4)
# dev.new(xpos=20, ypos=20, height=4, width=4)
# par(kPlotParameters)
# xmin <- 0.3
# xmax <- 60
# ymin <- 0.5
# ymax <- 300
# BlankPlot(log="xy")
# AddLogAxis(1, label="Conc")
# AddLogAxis(2, label="R")
# sapply(rownames(sXc), function(site) {
#   points(x, sXc_R[site, ], type="o", col=kSiteColors[site], pch=20)
# })
# dev.new(xpos=420, ypos=20, height=4, width=4)
# par(kPlotParameters)
# BlankPlot(log="xy")
# AddLogAxis(1, label="Conc")
# AddLogAxis(2, label="R")
# sapply(rownames(sXc), function(site) {
#   points(x, sXc_1s_R[site, ], type="o", col=kSiteColors[site], pch=20)
# })




double.site.df <- data.frame(left=site_2s_l, dist=dist, right=site_2s_r)
double.site.df <- cbind(double.site.df, sXc_2s_dist)
print(head(double.site.df))
double.pairs <- expand.grid(setdiff(rownames(sXc), "None"),
                            setdiff(rownames(sXc), "None"),
                            stringsAsFactors=FALSE)
print(head(double.pairs))
double.pairs.strings <- unique(apply(double.pairs, 1, function(pair) {
  paste0(sort(pair), collapse=" ")
}))
print(head(double.pairs.strings))


# Make the triple site strings.
# First remove the intervening ",(N)," portions of each string, sort the
# resulting names and make a string with spaces.
ts_name_vec <- sapply(rownames(sXc_3s_dist), function(name) {
  paste0(sort(unlist(strsplit(name, split=","))[c(1, 3, 5)]), collapse=" ")
})


ts_names_uniq <- unique(ts_name_vec)

print(head(ts_names_uniq))
print(head(ts_name_vec))
time_1 <- proc.time()[3]

tick <- 0
sXc_3s <- t(sapply(ts_names_uniq, function(ts_name_uniq) {
  if (tick%%1000 == 0) {
    print(tick/length(ts_names_uniq))
    print(proc.time()[3] - time_1)
  }
  tick <<- tick + 1
  colSums(sXc_3s_dist[which(ts_name_vec == ts_name_uniq), ])
}))

print(head(sXc_3s))

break

# Make the collapsed double-site containing matrix:
time_1 <- proc.time()[3]
tick <- 0

sXc_2s <- t(sapply(double.pairs.strings, function(col) {
  if (tick%%1000 == 0) {
    print(tick/length(double.pairs.strings))
    print(proc.time()[3] - time_1)
  }
  tick <<- tick + 1
  splits <- unlist(strsplit(col, split=" "))
  sum_1 <- colSums(subset(double.site.df,
                 left==splits[1] & right==splits[2])[, 4:ncol(double.site.df)])
  sum_2 <- colSums(subset(double.site.df,
                 left==splits[2] & right==splits[1])[, 4:ncol(double.site.df)])
  sum_1 + sum_2
}))
print(proc.time()[3] - time_1)

print(head(sXc_2s))

list_sXcs <- list(sXc_1s, sXc_2s, sXc_3s_dist, sXc_dol, sXc_tol,
                  sXc_r_ol_ts, sXc_l_ol_ts)
names(list_sXcs) <- c("single sites", "double sites", "triple sites",
                      "double ol sites", "triple ol sites",
                      "right-ol triple sites", "left_ol, triple sites")


# sXc_final <- list(sXc_1s, double.site.combined)

print(head(sXc))

ind_None <- nrow(sXc_1s)
sXc_final_list <- list(sXc_1s[-ind_None, ], double.site.combined,
                       sXc_3psup_sites, sXc_2s_seed_sites,
                       sXc_1s[ind_None, , drop=FALSE])

# print(dim(sXc_final))
print(dim(sXc))
print(head(sXc_final))
print(tail(sXc_final))

print(colSums(sXc_final) / colSums(sXc))

break

print(sXc)

break

pars_sd <- 0.1
plot_ <- FALSE
plotname <- NULL

InitializeEquilSitePars <- function(sXc, combined=TRUE) {
  L_ <- 100
  l_1s <- ((sXc[[1]][, 1 + combined] + 1)*L_/
           sum(sum(sXc[[1]][, 1 + combined] + 1),
               sum(sXc[[2]][, 1 + combined] + 1)))
  l_2s <- ((sXc[[2]][, 1 + combined] + 2)*L_/
           sum(sum(sXc[[1]][, 1 + combined] + 1),
               sum(sXc[[2]][, 1 + combined] + 1)))
  l <- c(l_1s, l_2s)
  data_2s <- SubfunctionCall(GetDataEquil, sXc=sXc[[1]])
  data_1s <- SubfunctionCall(GetDataEquil, sXc=sXc[[1]])
  kds <- log(l_1s/sum(l)/
             (rowSums(data_1s)/(sum(data_2s) + sum(data_1s))))
  kds <- kds - kds[length(kds)]
  bgs <- log(0.1)
  As <- log(1)
  names(kds) <- paste0(rownames(sXc[[1]]), "_Kd")  
  names(bgs) <- sprintf("bg_%s", mirna)
  names(As) <- sprintf("AGO_%s", mirna)
  pars <- c(kds, bgs, As)
  pars["None_Kd"] <- 0
  return(pars)
}

tick <- 0
OptimizeEquilSitePars <- function(sXc, pars=NULL,
                                  plotname_=plotname) {
  time_start <- proc.time()[3]
  tick <- 0
  if (is.null(pars) == TRUE) {
    initial.pars <- InitializeEquilSitePars(sXc, combined=combined)
  } else {
    initial.pars <- pars
  }
  L_ <- 100

  data_1s <- SubfunctionCall(GetDataEquil, sXc=sXc[[1]])
  data_2s <- SubfunctionCall(GetDataEquil, sXc=sXc[[2]])
  pars_old <- log(EquilPars(mirna, sitelist="mismatch_and_threeprime")$Full)
  names(pars_old) <- names(initial.pars)
  initial.pars <- pars_old
  n_j <- ncol(sXc[[1]]) - 3
  l_1s <- ((sXc[[1]][, 1 + combined] + 1)*L_/
           sum(sum(sXc[[1]][, 1 + combined] + 1),
               sum(sXc[[2]][, 1 + combined] + 1)))
  l_2s <- ((sXc[[2]][, 1 + combined] + 2)*L_/
           sum(sum(sXc[[1]][, 1 + combined] + 1),
               sum(sXc[[2]][, 1 + combined] + 1)))
  names(l_1s) <- rownames(sXc[[1]])
  names(l_2s) <- rownames(sXc[[2]])
  # Make the combined l term, for the root solve:
  # Sum the rows and columns across each of the double sites:
  l_2s_l <- c(colSums(matrix(l_2s, nrow=sqrt(length(l_2s)))), 0)
  l_2s_r <- c(rowSums(matrix(l_2s, nrow=sqrt(length(l_2s)))), 0)
  # Add them to the single site vector to make the total l vector for Free Ago:
  l <- l_1s + l_2s_l + l_2s_r
  Y <- colSums(data_1s) + colSums(data_2s)
  Y_1s <- colSums(data_1s)
  dil <- as.numeric(colnames(data_1s))/100.0
  n_i <- nrow(data_1s)
  data_1s <- as.numeric(as.matrix(data_1s))
  data_2s <- as.numeric(as.matrix(data_2s))
  n_mir <- 1
  lower <- log(rep(c(1e-4, 0.01, 0.1), c(n_i, n_mir, n_mir)))
  upper <- log(rep(c(10, 10, 10), c(n_i, n_mir, n_mir)))

  solution_1s <- optim(initial.pars,
                    CostCFinalTemp,
                    gr = GradCFinalTemp,
                    data = data_1s,
                    dil  = dil,
                    l = l_1s/sum(l_1s)*L_,
                    L = L_,
                    Y = Y_1s,
                    n_i = n_i,
                    n_j = n_j,
                    n_mir = 1,
                    fixed = FALSE,
                    upper_ = upper,
                    lower_ = lower,
                    plot_ = TRUE,
                    plotname = plotname_,
                    method = "L-BFGS-B",
                    lower=lower,
                    upper=upper,
                    control = list(maxit=10000000, factr=10000, fnscale=1))

  pars_1s <- solution_1s$par
  # solution_2s <- optim(pars_1s[-(length(initial.pars) - 2)],
  #                   CostCFinalDoubleSite,
  #                   # gr = GradCNumericalDoubleSite,
  #                   data_1s = data_1s,
  #                   data_2s = data_2s,
  #                   dil  = dil,
  #                   l_1s = l_1s,
  #                   l_2s = l_2s,
  #                   l = l,
  #                   L = L_,
  #                   Y = Y,
  #                   n_i = n_i,
  #                   n_j = n_j,
  #                   upper_ = upper,
  #                   lower_ = lower,
  #                   plot_ = TRUE,
  #                   plotname = plotname_,
  #                   addNone = TRUE,
  #                   method = "Nelder-Mead",
  #                   # lower=lower,
  #                   # upper=upper,
  #                   control = list(maxit=10000000, factr=10000, fnscale=1))
  # pars_2s <- solution_2s$par
  # solution_2s_b <- optim(pars_2s,
  #                   CostCFinalDoubleSite,
  #                   # gr = GradCNumericalDoubleSite,
  #                   data_1s = data_1s,
  #                   data_2s = data_2s,
  #                   dil  = dil,
  #                   l_1s = l_1s,
  #                   l_2s = l_2s,
  #                   l = l,
  #                   L = L_,
  #                   Y = Y,
  #                   n_i = n_i,
  #                   n_j = n_j,
  #                   upper_ = upper,
  #                   lower_ = lower,
  #                   plot_ = TRUE,
  #                   plotname = plotname_,
  #                   addNone = TRUE,
  #                   method = "Nelder-Mead",
  #                   # lower=lower,
  #                   # upper=upper,
  #                   control = list(maxit=10000000, factr=10000, fnscale=1))


  print("done 2")

  # pars_2s_b <- solution_2s_b$par
  print(proc.time()[3] - time_start)
  # print(tick)
  # pars_2s <- c(pars_2s[1:(length(pars_2s) - 2)], 0, pars_2s[(length(pars_2s) - 1):length(pars_2s)])
  # pars_2s_b <- c(pars_2s_b[1:(length(pars_2s_b) - 2)], 0, pars_2s_b[(length(pars_2s_b) - 1):length(pars_2s_b)])
  pars_1s/log(10)
}

GetModel <- function(pars, sXc) {
  pars <- pars*log(10)
  time_start <- proc.time()[3]
  tick <- 0
  L_ <- 100
  data_1s_mat <- SubfunctionCall(GetDataEquil, sXc=sXc[[1]])
  data_2s_mat <- SubfunctionCall(GetDataEquil, sXc=sXc[[2]])
  n_j <- ncol(sXc[[1]]) - 3
  l_1s <- ((sXc[[1]][, 1 + combined] + 1)*L_/
           sum(sum(sXc[[1]][, 1 + combined] + 1),
               sum(sXc[[2]][, 1 + combined] + 1)))
  l_2s <- ((sXc[[2]][, 1 + combined] + 2)*L_/
           sum(sum(sXc[[1]][, 1 + combined] + 1),
               sum(sXc[[2]][, 1 + combined] + 1)))
  names(l_1s) <- rownames(sXc[[1]])
  names(l_2s) <- rownames(sXc[[2]])
  # Make the combined l term, for the root solve:
  # Sum the rows and columns across each of the double sites:
  l_2s_l <- c(colSums(matrix(l_2s, nrow=sqrt(length(l_2s)))), 0)
  l_2s_r <- c(rowSums(matrix(l_2s, nrow=sqrt(length(l_2s)))), 0)
  # Add them to the single site vector to make the total l vector for Free Ago:
  l <- l_1s + l_2s_l + l_2s_r
  Y <- colSums(data_1s_mat) + colSums(data_2s_mat)
  dil <- as.numeric(colnames(data_1s_mat))/100.0
  n_i <- nrow(data_1s_mat)
  data_1s <- as.numeric(as.matrix(data_1s_mat))
  data_2s <- as.numeric(as.matrix(data_2s_mat))
  n_mir <- 1
  lower <- log(rep(c(1e-4, 0.01, 0.1), c(n_i, n_mir, n_mir)))
  upper <- log(rep(c(10, 10, 10), c(n_i, n_mir, n_mir)))

  model <- Model_C_DoubleSite(pars, data_1s, data_2s, dil, l_1s, l_2s, l, L_,
                                 Y, n_i, n_j, upper_, lower_)
  # print(dim(sXc[[1]]))
  # print(dim(sXc[[2]]))
  # print(dim(model[[1]]))
  # print(dim(model[[2]]))

  rownames(model[[1]]) <- rownames(data_1s_mat)
  rownames(model[[2]]) <- rownames(data_2s_mat)
  colnames(model[[1]]) <- colnames(data_1s_mat)
  colnames(model[[2]]) <- colnames(data_2s_mat)
  # print(head(data_1s))
  # print(head(model[[1]]))
  # print(head(data_2s))
  # print(head(model[[2]]))

  # graphics.off()
  # dev.new(xpos=20, ypos=20, height=3, width=3)
  # par(kPlotParameters)
  # xmin <- 1
  # xmax <- 1e3
  # ymin <- xmin
  # ymax <- xmax
  # BlankPlot(log='xy')
  # inds5p <- grep("9mer-m11.19 ", rownames(data_2s))
  # inds3p <- grep(" 9mer-m11.19", rownames(data_2s))

  # points(model[[2]][inds3p,], data_2s[inds3p, ], col=ConvertRColortoRGB("red", alpha=0.5))
  # points(model[[2]][inds5p,], data_2s[inds5p, ], col=ConvertRColortoRGB("blue", alpha=0.5))
  # AddLogAxis(1, label="model")
  # AddLogAxis(2, label="data")

  # dev.new(xpos=20, ypos=330, height=3, width=3)
  # par(kPlotParameters)
  # BlankPlot(log='xy')
  # inds5p <- grep("9mer-m10.18 ", rownames(data_2s))
  # inds3p <- grep(" 9mer-m10.18", rownames(data_2s))

  # points(model[[2]][inds3p,], data_2s[inds3p, ], col=ConvertRColortoRGB("red", alpha=0.5))
  # points(model[[2]][inds5p,], data_2s[inds5p, ], col=ConvertRColortoRGB("blue", alpha=0.5))
  # AddLogAxis(1, label="model")
  # AddLogAxis(2, label="data")
  # dev.new(xpos=20, ypos=640, height=3, width=3)
  # par(kPlotParameters)
  # BlankPlot(log='xy')
  # inds5p <- grep("9mer-m9.17 ", rownames(data_2s))
  # inds3p <- grep(" 9mer-m9.17", rownames(data_2s))

  # points(model[[2]][inds3p,], data_2s[inds3p, ], col=ConvertRColortoRGB("red", alpha=0.5))
  # points(model[[2]][inds5p,], data_2s[inds5p, ], col=ConvertRColortoRGB("blue", alpha=0.5))
  # AddLogAxis(1, label="model")
  # AddLogAxis(2, label="data")

  # dev.new(xpos=320, ypos=20, height=3, width=3)
  # par(kPlotParameters)
  # BlankPlot(log='xy')
  # inds5p <- grep("10mer-m9.18 ", rownames(data_2s))
  # inds3p <- grep(" 10mer-m9.18", rownames(data_2s))

  # points(model[[2]][inds3p,], data_2s[inds3p, ], col=ConvertRColortoRGB("red", alpha=0.5))
  # points(model[[2]][inds5p,], data_2s[inds5p, ], col=ConvertRColortoRGB("blue", alpha=0.5))
  # AddLogAxis(1, label="model")
  # AddLogAxis(2, label="data")

  # dev.new(xpos=320, ypos=330, height=3, width=3)
  # par(kPlotParameters)
  # BlankPlot(log='xy')
  # inds5p <- grep("10mer-m10.19 ", rownames(data_2s))
  # inds3p <- grep(" 10mer-m10.19", rownames(data_2s))
  # points(model[[2]][inds3p,], data_2s[inds3p, ], col=ConvertRColortoRGB("red", alpha=0.5))
  # points(model[[2]][inds5p,], data_2s[inds5p, ], col=ConvertRColortoRGB("blue", alpha=0.5))
  # AddLogAxis(1, label="model")
  # AddLogAxis(2, label="data")

  # dev.new(xpos=320, ypos=640, height=3, width=3)
  # par(kPlotParameters)
  # BlankPlot(log='xy')
  # inds5p <- grep("11mer-m9.19 ", rownames(data_2s))
  # inds3p <- grep(" 11mer-m9.19", rownames(data_2s))
  # print(inds5p)
  # points(model[[2]][inds3p,], data_2s[inds3p, ], col=ConvertRColortoRGB("red", alpha=0.5))
  # points(model[[2]][inds5p,], data_2s[inds5p, ], col=ConvertRColortoRGB("blue", alpha=0.5))
  # AddLogAxis(1, label="model")
  # AddLogAxis(2, label="data")
  model

}

pars.MLE <- OptimizeEquilSitePars(sXc)
print(10^pars.MLE)
dev.new(xpos=420, ypos=20, height=4, width=4)
par(kPlotParameters)
xmin <- 0.001
xmax <- 2
ymin <- 1e-1
ymax <- 1e1
kds_paper <- EquilPars(mirna, sitelist="paper")
sites_paper <- intersect(rownames(kds_paper), names(pars.MLE))
print(sites_paper)
print(kds_paper)
BlankPlot(log="xy", inv="x")
AddLogAxis(1, label="Relative Kd; single site fit")
AddLogAxis(2, label="single-site/paper fit; relative Kd")

kds_paper_vec <- kds_paper[sites_paper, ]$Full
kds_paper_del <- kds_paper[sites_paper, ]$Full/c((10^pars.MLE[sites_paper]))


points(kds_paper[sites_paper,]$Full, kds_paper_del,
       col=kSiteColors[gsub(sprintf("_(Kd|%s)", mirna), "", sites_paper)],)
xy <- GetPlotFractionalCoords(0.05, 0.95, log='xy', inv='x')
text(xy[1], xy[2], labels=paste0(mirna, " MultiSite"), adj=0)


break

# pars.MLE <- OptimizeEquilSitePars(sXc)

model <- GetModel(pars.MLE, sXc)

kds_old <- EquilPars(mirna, sitelist="paper")
kdssingle <- EquilPars(mirna, sitelist="mismatch_and_threeprime")

sites_old <- rownames(kds_old)
xmin <- 0.001
xmax <- 2
ymin <- 1e-1
ymax <- 1e1
BlankPlot(log="xy", inv="x")
AddLogAxis(1, label="Conc")
AddLogAxis(2, label="R")

points(10^pars.MLE[, 1][sites_old], 10^pars.MLE_del[sites_old],
       col=kSiteColors[c(gsub("_Kd", "", sites_old), "bg", "AGO")],)





break

points(c(10^pars.MLE[, 1][sites_old]), kds_old$Full,
     col=kSiteColors[c(gsub("_Kd", "", rownames(kds_old)), "bg", "AGO")], pch=1, lwd=2)

points(c(10^pars.MLE[, 2][sites_old]), kds_old$Full,
     col=kSiteColors[c(gsub("_Kd", "", rownames(kds_old)), "bg", "AGO")])



break
dev.new(xpos=20, ypos=420, height=4, width=4)
par(kPlotParameters)
xmin <- 0.3
xmax <- 60
ymin <- 0.5
ymax <- 300
BlankPlot(log="xy")
AddLogAxis(1, label="Conc")
AddLogAxis(2, label="R")


model_1s <- model[[1]]
model_2s <- model[[2]]

model_totals <- colSums(model_1s) + colSums(model_2s)

model_1s_norm <- t(t(model_1s)/model_totals)
model_1s_R <- model_1s_norm/sXc_1s_norm[, 1 + combined]
print(379)
sapply(rownames(sXc[[1]]), function(site) {
  print(site)
  points(x, sXc_R[site, ], col=kSiteColors[site], pch=20)
  lines(x, model_1s_R[site, ], col=kSiteColors[site], pch=20)

})
break
dev.new(xpos=420, ypos=420, height=4, width=4)
par(kPlotParameters)
BlankPlot(log="xy")
AddLogAxis(1, label="Conc")
AddLogAxis(2, label="R")
sapply(rownames(sXc), function(site) {
  points(x, sXc_1s_R[site, ], type="o", col=kSiteColors[site], pch=20)
})


break
pars_loocv <- matrix(NaN,nrow=length(pars),ncol=reps)
rownames(pars_loocv) <- names(pars)
colnames(pars_loocv) <- seq(reps)
for (i_trial in seq(reps)) {
  print(i_trial)
  tick <- 0
  i_col <- sample(1:ncol(data.all),1)
  data.temp <- data.all[,-i_col]
  data <- apply(data.temp, 2, function(col) {rmultinom(1,size=sum(col), prob=col)})
  rownames(data) <- rownames(data.temp)
  colnames(data) <- colnames(data.temp)
  input_resample <- rmultinom(1,size=sum(sitesXcounts[,2]),prob=sitesXcounts[,2])
  c.I.tots <- Norm(input_resample+1)*k.c.lib
  c.totals <- matrix(
                rep(c.I.tots,ncol(data)), nrow=length(c.I.tots), ncol=ncol(data), byrow=FALSE)
  colnames(c.totals) <- colnames(data)
  rownames(c.totals) <- rownames(data)
  pars <- pars.save
  solution <- optim(pars,
                    ModelLikelihood,
                    gr = Gradient,
                    method = "L-BFGS-B",
                    lower=c(rep(-13, length=num.kds), -2, -1),
                    upper=c(rep(10, length=num.kds), 1, 2),
                    control = list(maxit=100000, factr=1e2))

  pars <- solution$par

  pars_loocv[,i_trial] <- pars

}

pars_loocv_sort <- t(apply(pars_loocv,1,sort))
print(pars_loocv_sort)
print(pars.save)
print(10^pars.save)
print(10^cbind(pars.save,rowMeans(pars_loocv)))
output <- 10^cbind(pars.save,rowMeans(pars_loocv),
                   pars_loocv_sort[,ceiling(0.025*reps)],
                   pars_loocv_sort[,ceiling(0.5*reps)],
                   pars_loocv_sort[,ceiling(0.975*reps)])
colnames(output) <- c("Full","Mean", "Lower_CI", "Median", "Upper_CI")
write.table(file=out_file, output, sep="\t", quote = FALSE)
print(output)
print(out_file)


warnings()





